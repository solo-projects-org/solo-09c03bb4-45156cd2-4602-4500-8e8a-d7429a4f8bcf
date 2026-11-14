"""生成交叉证认演示用的虚构坐标样例(固定数值,可复现)。

输出两个文件:
  candidates_sample.csv —— 5 颗候选(candidate_id/job_id/ra_deg/dec_deg),
  catalog_sample.csv    —— 5 个星表源(source_id/ra_deg/dec_deg)。

候选 CSV 里的 job_id 是占位符 JOB_A / JOB_A2 / JOB_B / JOB_B2 / JOB_C
(每颗候选各占一个已完成测光批次,满足 job_id 不重复);网页
「载入关联样例并演示」会先确保五批样例测光任务完成,再把占位符
替换成真实批号后提交。本地用 curl 演示时同样需要先替换为真实 job_id。

配合阈值 1.5 角秒,样例覆盖以下情形
------------------------------------
A 组(赤经跨 0° + 同等最优的完美匹配二义性):
  候选 C1、C2 分居 0° 子午线西侧与东侧,源 S1、S2 在同一 RA 上南北
  各偏 1 角秒;4 条候选—源边严格等于 1.000 角秒,构成 2×2 完全二义,
  两种完美匹配总角距相同(2.000 角秒),C1、C2 均标为待复核,
  各自列出两个可选源。

B 组(源争用 + 未匹配也是同等最优选项):
  候选 C3、C4 都只能匹配唯一的源 S3,且两条边都恰为 1.000 角秒;
  最优匹配数只能为 1,两种分配(给 C3 或给 C4)同等最优,因此两颗
  候选都待复核:可选源 S3 与「未匹配」同时出现。

C 组(确定性匹配):
  C5 与 S4 相距 0.500 角秒,唯一最优,直接认定。

S5 离所有候选都超过阈值,演示「确定无法匹配的星表源」清单。
"""
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent

# A 组放在高纬,让 RA 方向每角秒对应的度数足够大,跨 0° 的坐标简洁可读。
DEC0 = 89.5


def offset(ra0: float, dec0: float, dra_arcsec: float,
           ddec_arcsec: float) -> tuple[float, float]:
    """以切点平面上的东/北角秒偏移构造近似坐标。"""
    cos_dec = math.cos(math.radians(dec0)) or 1e-12
    return (ra0 + dra_arcsec / 3600.0 / cos_dec,
            dec0 + ddec_arcsec / 3600.0)


# ---------------------------------------------------------------- A 组 ----
# 两颗候选同居 0° 子午线两侧(RA 几乎相同,一负一正,验证跨 0° 配对),
# 两个源在同一 RA 上南北各偏 1";于是 4 条候选—源边严格等于 1.000",
# 构成 2×2 完全二义:两种完美匹配的总角距都是 2.000"。
ra_a, dec_a = 0.0, DEC0
c1 = offset(ra_a, dec_a, -0.001, 0.0)    # C1 在子午线西侧(RA 为负)
c2 = offset(ra_a, dec_a, 0.001, 0.0)     # C2 在子午线东侧
s1 = offset(ra_a, dec_a, 0.0, 1.0)       # S1 偏北 1"
s2 = offset(ra_a, dec_a, 0.0, -1.0)      # S2 偏南 1"

# ---------------------------------------------------------------- B 组 ----
ra_b, dec_b = 120.0, -30.0
s3 = offset(ra_b, dec_b, 0.0, 0.0)
c3 = offset(ra_b, dec_b, 0.0, 1.0)
c4 = offset(ra_b, dec_b, 0.0, -1.0)

# ---------------------------------------------------------------- C 组 ----
ra_c, dec_c = 240.0, 10.0
# 切点平面东 0.3"、北 0.4" → 与 C5 的大圆距离 0.500"
s4 = offset(ra_c, dec_c, 0.3, 0.4)
c5 = (ra_c, dec_c)
s5 = offset(ra_c, dec_c, 10.0, 0.0)    # 离 C5 十角秒,确定无法匹配的源


def fmt(x: float) -> str:
    s = f"{x:.7f}".rstrip("0").rstrip(".")
    return s if s else "0"


candidates = [
    ("C1", "JOB_A", c1),
    ("C2", "JOB_A2", c2),   # 题面要求候选内 job_id 不重复,各占一个占位符
    ("C3", "JOB_B", c3),
    ("C4", "JOB_B2", c4),
    ("C5", "JOB_C", c5),
]
sources = [
    ("S1", s1), ("S2", s2), ("S3", s3), ("S4", s4), ("S5", s5),
]

with (HERE / "candidates_sample.csv").open("w", encoding="utf-8") as fh:
    fh.write("candidate_id,job_id,ra_deg,dec_deg\n")
    for cid, jid, (ra, dec) in candidates:
        fh.write(f"{cid},{jid},{fmt(ra)},{fmt(dec)}\n")

with (HERE / "catalog_sample.csv").open("w", encoding="utf-8") as fh:
    fh.write("source_id,ra_deg,dec_deg\n")
    for sid, (ra, dec) in sources:
        fh.write(f"{sid},{fmt(ra)},{fmt(dec)}\n")

print("已写出 candidates_sample.csv(5 行)与 catalog_sample.csv(5 行)")
print("建议最大角距阈值:1.5 角秒")
