"""候选变星与已知星表源的交叉证认(一对一最优分配)。

规则
----
* 候选 CSV 含 candidate_id、job_id、ra_deg、dec_deg,不超过 20 行,
  candidate_id 唯一、job_id 不重复;星表 CSV 含 source_id、ra_deg、
  dec_deg,不超过 50 行,source_id 唯一;
* 角距按球面大圆距离(haversine)计算,天然处理赤经跨 0°;
* 只保留角距不超过阈值(角秒)的候选—源边;
* 每颗候选、每个星表源最多使用一次;先使匹配数最多,再使总角距最小;
* 最优性比较一律使用四舍五入到 0.001 角秒的整数角距(毫角秒单位),
  全程整数费用,避免浮点误差;
* 若存在多个同等最优分配,受影响候选列出其在全部最优分配中出现过的
  所有可选源以及未匹配可能,标为待复核,不任意认定其中一种。

实现
----
把分配建模为最小费用最大流:源点 → 候选 → 星表源 → 汇点,另给每条
候选 → 汇点的「未匹配」边赋予远大于任何真实角距费用的惩罚 BIG,
于是最小费用首先最小化未匹配数(=最大化匹配数),其次最小化角距和。
歧义判定采用强制重解:对每颗候选的每个可选源(以及强制未匹配)各
解一次剩余问题,达到全局最优值即说明该选择属于某个最优分配。
规模上限 20×50,重解次数不超过 20×51,图很小,耗时可忽略。
"""
from __future__ import annotations

import heapq
import io
import math

import pandas as pd
from fastapi import HTTPException

ARCSEC_PER_RAD = 206264.80624709636  # 180*3600/pi
MAX_CANDIDATES = 20
MAX_SOURCES = 50
# 全天最大角距 648000 角秒 = 648_000_000 毫角秒,20 条边总和不超过
# 1.296e10;BIG 取 1e12,保证任何一条真实边都比未匹配惩罚便宜。
BIG = 10**12


# --------------------------------------------------------------------------- #
# 输入校验
# --------------------------------------------------------------------------- #

def _read_df(raw: bytes, name: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(io.BytesIO(raw))
    except Exception as exc:
        raise HTTPException(400, f"{name} CSV 解析失败:{exc}")
    df.columns = [str(c).strip() for c in df.columns]
    return df


def _check_ids(df: pd.DataFrame, col: str, label: str) -> list[str]:
    ids = df[col].map(lambda x: str(x).strip()).to_list()
    if any(not x for x in ids):
        raise HTTPException(400, f"{label}不允许为空值")
    if len(set(ids)) != len(ids):
        raise HTTPException(400, f"{label}必须唯一,存在重复")
    return ids


def _check_coords(df: pd.DataFrame) -> tuple[list[float], list[float]]:
    ra = pd.to_numeric(df["ra_deg"], errors="coerce").to_numpy(dtype=float)
    dec = pd.to_numeric(df["dec_deg"], errors="coerce").to_numpy(dtype=float)
    if not (len(ra) > 0 and len(ra) == len(dec)):
        raise HTTPException(400, "CSV 至少需要一行数据")
    if not (pd.Series(ra).notna().all() and pd.Series(dec).notna().all()
            and all(math.isfinite(x) for x in ra)
            and all(math.isfinite(x) for x in dec)):
        raise HTTPException(400, "ra_deg、dec_deg 必须全部为有限数值"
                                 "(不允许空值/NaN/Inf)")
    if any(d < -90.0 or d > 90.0 for d in dec):
        raise HTTPException(400, "dec_deg 必须落在 [-90, 90] 度之间")
    # 赤经归一化到 [0, 360),跨 0° 的配对在距离函数里处理
    ra = [x % 360.0 for x in ra]
    return ra, list(dec)


def parse_candidates(raw: bytes) -> list[dict]:
    """解析并校验候选 CSV,返回规范化记录列表。"""
    df = _read_df(raw, "候选")
    missing = [c for c in ("candidate_id", "job_id", "ra_deg", "dec_deg")
               if c not in df.columns]
    if missing:
        raise HTTPException(400, f"候选 CSV 缺少必需列:{', '.join(missing)}"
                                 "(需要 candidate_id、job_id、ra_deg、dec_deg)")
    if len(df) > MAX_CANDIDATES:
        raise HTTPException(400, f"候选 CSV 不超过 {MAX_CANDIDATES} 行,"
                                 f"当前 {len(df)} 行")
    if len(df) < 1:
        raise HTTPException(400, "候选 CSV 至少需要一行数据")
    ids = _check_ids(df, "candidate_id", "candidate_id")
    job_ids = _check_ids(df, "job_id", "job_id")  # 不重复即唯一
    ra, dec = _check_coords(df)
    return [{"candidate_id": cid, "job_id": jid, "ra_deg": r, "dec_deg": d}
            for cid, jid, r, d in zip(ids, job_ids, ra, dec)]


def parse_sources(raw: bytes) -> list[dict]:
    """解析并校验星表 CSV,返回规范化记录列表。"""
    df = _read_df(raw, "星表")
    missing = [c for c in ("source_id", "ra_deg", "dec_deg")
               if c not in df.columns]
    if missing:
        raise HTTPException(400, f"星表 CSV 缺少必需列:{', '.join(missing)}"
                                 "(需要 source_id、ra_deg、dec_deg)")
    if len(df) > MAX_SOURCES:
        raise HTTPException(400, f"星表 CSV 不超过 {MAX_SOURCES} 行,"
                                 f"当前 {len(df)} 行")
    if len(df) < 1:
        raise HTTPException(400, "星表 CSV 至少需要一行数据")
    ids = _check_ids(df, "source_id", "source_id")
    ra, dec = _check_coords(df)
    return [{"source_id": sid, "ra_deg": r, "dec_deg": d}
            for sid, r, d in zip(ids, ra, dec)]


# --------------------------------------------------------------------------- #
# 球面距离
# --------------------------------------------------------------------------- #

def angular_separation_arcsec(ra1: float, dec1: float,
                              ra2: float, dec2: float) -> float:
    """两点球面大圆距离(角秒)。ΔRA 折回 [-180,180),正确处理跨 0°。"""
    r1, d1 = math.radians(ra1), math.radians(dec1)
    r2, d2 = math.radians(ra2), math.radians(dec2)
    dra = (r2 - r1 + math.pi) % (2.0 * math.pi) - math.pi
    ddec = d2 - d1
    a = (math.sin(ddec / 2.0) ** 2
         + math.cos(d1) * math.cos(d2) * math.sin(dra / 2.0) ** 2)
    a = min(1.0, max(0.0, a))
    return 2.0 * math.asin(math.sqrt(a)) * ARCSEC_PER_RAD


def round_milli_arcsec(distance_arcsec: float) -> int:
    """四舍五入到 0.001 角秒,返回整数毫角秒(半入,非银行家舍入)。"""
    return math.floor(distance_arcsec * 1000.0 + 0.5)


# --------------------------------------------------------------------------- #
# 最小费用最大流(整数费用,Dijkstra + 势函数)
# --------------------------------------------------------------------------- #

class _Edge:
    __slots__ = ("to", "rev", "cap", "cost")

    def __init__(self, to: int, rev: int, cap: int, cost: int) -> None:
        self.to = to
        self.rev = rev
        self.cap = cap
        self.cost = cost


class MinCostFlow:
    def __init__(self, n: int) -> None:
        self.g: list[list[_Edge]] = [[] for _ in range(n)]

    def add_edge(self, fr: int, to: int, cap: int, cost: int) -> _Edge:
        fwd = _Edge(to, len(self.g[to]), cap, cost)
        rev = _Edge(fr, len(self.g[fr]), 0, -cost)
        self.g[fr].append(fwd)
        self.g[to].append(rev)
        return fwd

    def flow(self, s: int, t: int, maxf: int) -> tuple[int, int]:
        n = len(self.g)
        h = [0] * n          # 势函数:初始费用均非负,从零开始
        prev_v = [-1] * n
        prev_e = [-1] * n
        total_flow = total_cost = 0
        INF = 10**30
        while total_flow < maxf:
            dist = [INF] * n
            dist[s] = 0
            pq: list[tuple[int, int]] = [(0, s)]
            while pq:
                d, v = heapq.heappop(pq)
                if d != dist[v]:
                    continue
                for ei, e in enumerate(self.g[v]):
                    if e.cap <= 0:
                        continue
                    nd = d + e.cost + h[v] - h[e.to]
                    if nd < dist[e.to]:
                        dist[e.to] = nd
                        prev_v[e.to] = v
                        prev_e[e.to] = ei
                        heapq.heappush(pq, (nd, e.to))
            if dist[t] == INF:
                break  # 无增广路
            for v in range(n):
                if dist[v] < INF:
                    h[v] += dist[v]
            add = maxf - total_flow
            v = t
            while v != s:
                add = min(add, self.g[prev_v[v]][prev_e[v]].cap)
                v = prev_v[v]
            v = t
            while v != s:
                e = self.g[prev_v[v]][prev_e[v]]
                e.cap -= add
                self.g[v][e.rev].cap += add
                v = prev_v[v]
            total_flow += add
            total_cost += add * h[t]
        return total_flow, total_cost


def _solve(cand_idx: list[int], src_idx: list[int],
           weights: dict[tuple[int, int], int]) -> tuple[int, int]:
    """在指定候选/源子集上求最优分配。

    返回 (匹配数, 匹配边的毫角秒费用和);未匹配惩罚 BIG 不计入角距和。
    """
    n, m = len(cand_idx), len(src_idx)
    node_s = 0
    cand_node = {c: 1 + k for k, c in enumerate(cand_idx)}
    src_node = {s: 1 + n + k for k, s in enumerate(src_idx)}
    node_t = 1 + n + m
    mcf = MinCostFlow(node_t + 1)
    unmatched_edges: list[_Edge] = []
    for c in cand_idx:
        mcf.add_edge(node_s, cand_node[c], 1, 0)
        ue = mcf.add_edge(cand_node[c], node_t, 1, BIG)  # 未匹配
        unmatched_edges.append(ue)
    for s in src_idx:
        mcf.add_edge(src_node[s], node_t, 1, 0)
    for (c, s), w in weights.items():
        if c in cand_node and s in src_node:
            mcf.add_edge(cand_node[c], src_node[s], 1, w)
    _, total_cost = mcf.flow(node_s, node_t, n)
    n_unmatched = sum(1 for e in unmatched_edges if e.cap == 0)
    match_cost = total_cost - BIG * n_unmatched
    return n - n_unmatched, match_cost


# --------------------------------------------------------------------------- #
# 主流程
# --------------------------------------------------------------------------- #

def crossmatch(candidates: list[dict], sources: list[dict],
               max_sep_arcsec: float,
               periods: dict[str, float] | None = None) -> dict:
    """执行交叉证认并组装前端需要的结果结构。

    periods:job_id → 该已完成批次的最高峰周期(天)。
    """
    n, m = len(candidates), len(sources)
    # 阈值本身也按 0.001 角秒取整后比较,与最优性比较口径一致
    limit_milli = round_milli_arcsec(max_sep_arcsec)

    # 全量边表与费用矩阵
    weights: dict[tuple[int, int], int] = {}
    edges_of: list[list[dict]] = [[] for _ in range(n)]
    for i, c in enumerate(candidates):
        for j, s in enumerate(sources):
            d = angular_separation_arcsec(c["ra_deg"], c["dec_deg"],
                                          s["ra_deg"], s["dec_deg"])
            milli = round_milli_arcsec(d)
            if milli <= limit_milli:
                weights[(i, j)] = milli
                edges_of[i].append({
                    "source_index": j,
                    "source_id": s["source_id"],
                    "ra_deg": s["ra_deg"],
                    "dec_deg": s["dec_deg"],
                    "separation_arcsec": milli / 1000.0,
                    "milli": milli,
                })

    all_c, all_s = list(range(n)), list(range(m))
    k0, cost0 = _solve(all_c, all_s, weights)

    # 强制重解:判定每个候选的每种选择是否属于某个全局最优分配
    source_feasible_any = [False] * m
    results = []
    n_review = n_matched = n_unmatched = 0
    for i, c in enumerate(candidates):
        other_c = [x for x in all_c if x != i]
        options = []
        for opt in edges_of[i]:
            j = opt["source_index"]
            kk, cc = _solve(other_c, [x for x in all_s if x != j], weights)
            feasible = (kk + 1 == k0 and cc + opt["milli"] == cost0)
            if feasible:
                source_feasible_any[j] = True
            options.append({
                "kind": "source",
                "source_id": opt["source_id"],
                "ra_deg": opt["ra_deg"],
                "dec_deg": opt["dec_deg"],
                "separation_arcsec": opt["separation_arcsec"],
                "feasible": feasible,
            })
        # 强制该候选未匹配
        kk, cc = _solve(other_c, all_s, weights)
        unmatched_feasible = (kk == k0 and cc == cost0)
        options.append({"kind": "unmatched", "feasible": unmatched_feasible})

        feasible_sources = [o for o in options
                            if o["kind"] == "source" and o["feasible"]]
        if len([o for o in options if o["feasible"]]) > 1:
            status = "review"
            n_review += 1
            match = None
        elif len(feasible_sources) == 1:
            status = "matched"
            n_matched += 1
            o = feasible_sources[0]
            match = {"source_id": o["source_id"],
                     "separation_arcsec": o["separation_arcsec"]}
        else:
            # 唯一可行结局就是未匹配
            status = "unmatched"
            n_unmatched += 1
            match = None

        item = {
            "candidate_id": c["candidate_id"],
            "job_id": c["job_id"],
            "ra_deg": c["ra_deg"],
            "dec_deg": c["dec_deg"],
            "best_period_day": (periods or {}).get(c["job_id"]),
            "status": status,
            "match": match,
        }
        if status == "review":
            item["options"] = options
        results.append(item)

    unmatched_sources = [sources[j]["source_id"] for j in all_s
                         if not source_feasible_any[j]]
    return {
        "max_sep_arcsec": max_sep_arcsec,
        "summary": {
            "n_candidates": n,
            "n_sources": m,
            "optimal_matches": k0,
            "min_total_sep_arcsec": cost0 / 1000.0,
            "n_review": n_review,
            "n_matched": n_matched,
            "n_unmatched": n_unmatched,
        },
        "candidates": results,
        "unmatched_sources": unmatched_sources,
    }
