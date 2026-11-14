# 周期变星候选筛选与星表交叉证认(Lomb–Scargle + 置换检验 + 一对一最优分配)

面向天文巡测数据员的网页工具,包含两个环节:

1. **周期筛查**:上传单颗星的测光序列,自动搜索周期信号,输出最高峰周期、
   功率谱、相位折叠散点与置换尾概率,供人工复核;
2. **星表关联(交叉证认)**:上传多颗候选与已知变星星表,按球面角距做
   一对一最优分配,避免把同一个星表源认给两颗候选;存在同等最优分配时
   不任意认定,列出全部可选结局并标为待复核。

## 快速开始

```bash
docker compose up --build
# 打开 http://localhost:8000
```

页面上可直接点击「载入虚构样例并提交」(`sample/sample_star.csv`,
一颗周期约 4.25 天的虚构脉动变星,160 个观测点)。

关联演示点击「载入关联样例并演示」:页面会自动提交并等待 5 批样例测光
任务(受最多两批并行限制会排队),再把候选样例中的 job_id 占位符替换成
真实批号并发起关联,阈值已预填 1.5 角秒。

无 Docker 时本地运行:

```bash
python -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python -m uvicorn app.main:app --host 0.0.0.0 --port 8000
```

## 输入要求

- CSV 含 `time_day`、`flux` 两列(允许其他列存在,只取这两列);
- 至少 10 行;两列均为有限数值(无空值/NaN/Inf);观测时刻互异;
- 周期范围须满足 `0 < 下限 < 上限`。

## 计算流程

1. 周期范围取倒数得到频率范围 `[1/上限, 1/下限]`,均匀取 **10000** 个频率点;
2. 用**带浮动均值**的 Lomb–Scargle(`astropy.timeseries.LombScargle`,
   `fit_mean=True`)计算功率谱,取最高峰周期;
3. **置换检验**:固定时刻,随机重排光通量 **200** 次,每次重算全频段
   最大功率;尾概率 = `(不低于原峰的次数 + 1) / 201`;
4. 展示最高峰周期、峰值功率、尾概率、功率谱图与按峰周期折叠的相位散点图。

## 并发与生命周期语义

- 每批任务使用 **2 个计算子进程**(各承担 100 次置换,worker 0 另算原始
  功率谱);最多 **2 批并行**,超额提交返回 **HTTP 429**;
- 任务可查询、可取消(`POST /api/jobs/{id}/cancel`);
- 无论成功、失败还是取消,都要等本批全部子进程退出、临时文件清理完毕,
  才公布结局并释放名额;另一批不受影响;
- 取消一经受理即终止子进程,**不会发布迟到结果**;重复取消幂等,不会
  重复释放名额;
- 全部任务状态只存内存(重启即清空);计算期临时文件随收尾删除。

## 星表关联(交叉证认)

### 输入要求

- 候选 CSV:含 `candidate_id`、`job_id`、`ra_deg`、`dec_deg`,**不超过 20 行**;
  `candidate_id` 唯一、`job_id` 不重复;
- 星表 CSV:含 `source_id`、`ra_deg`、`dec_deg`,**不超过 50 行**;
  `source_id` 唯一;
- 坐标为十进制度,允许其他列存在;赤经归一化到 `[0,360)`,负赤经(如
  `-0.0000318`)会被折回 `359.9999…`,**赤经跨 0° 的配对在距离函数内处理**;
  赤纬须落在 `[-90, 90]`;
- 每个 `job_id` 必须指向当前内存中**已完成(成功)**的测光批次,否则整次
  请求返回 409;结果中每颗候选都会带上该批的**最高峰周期**;
- 最大角距为**正数**,单位角秒。

### 匹配规则

1. 角距按球面大圆距离(haversine)计算,只保留不超过阈值的候选—源边;
2. 每颗候选、每个星表源**最多使用一次**;优化目标按字典序:
   **先使匹配数最多,再使总角距最小**;
3. 所有角距先**四舍五入到 0.001 角秒**(整数毫角秒),再参与阈值过滤与
   最优性比较,全程整数费用,避免浮点误差;
4. 分配建模为最小费用最大流(源点→候选→星表源→汇点,未匹配边带大惩罚),
   天然实现上述字典序目标;
5. **同等最优分配不任意认定**:对每颗候选的每个可选源(以及强制未匹配)
   各重解一次剩余问题,凡能达到全局最优值的选择,都属于某个最优分配。
   - 候选只有一个可行结局(某源,或未匹配)→ 直接给出匹配源与角距,
     或标明未匹配;
   - 候选存在多个可行结局 → 状态为**待复核**,表格列出全部可选源及角距
     (并注明哪些属于最优分配)以及「未匹配」是否可行;
6. 同时输出确定不会出现在任何最优分配中的星表源清单。

### 虚构坐标样例(`sample/`,阈值 1.5 角秒)

| 组别 | 设计 | 预期结果 |
|---|---|---|
| A 组 C1/C2 vs S1/S2 | 候选分居 0° 子午线两侧(验证跨 0°),4 条边严格 1.000″ | 两种完美匹配等价,两颗候选均待复核,各列两个可选源 |
| B 组 C3/C4 vs S3 | 两颗候选都只能配同一个源,边均为 1.000″ | 源争用;两颗候选待复核,可选源与「未匹配」均可行 |
| C 组 C5 vs S4 | 相距 0.500″ | 唯一最优,直接认定 |
| S5 | 距所有候选 10″ 以上 | 进入「未使用星表源」清单 |

重新生成:`python sample/make_crossmatch_sample.py`(纯标准库,可复现)。
候选样例中的 `JOB_A/JOB_A2/JOB_B/JOB_B2/JOB_C` 是 job_id 占位符,
网页演示按钮会自动替换;用 curl 手动演示时需先提交测光批次并自行替换。

## API 摘要

| 方法 | 路径 | 说明 |
|---|---|---|
| GET | `/` | 网页 |
| POST | `/api/jobs` | 提交批处理(multipart:`file`、`period_min`、`period_max`)→ 201 / 400 / 429 |
| GET | `/api/jobs` | 批处理列表与名额占用 |
| GET | `/api/jobs/{id}` | 单批详情(含指标) |
| POST | `/api/jobs/{id}/cancel` | 取消(幂等) |
| GET | `/api/jobs/{id}/plot/power.png` | 功率谱图(成功后) |
| GET | `/api/jobs/{id}/plot/phase.png` | 相位折叠散点图(成功后) |
| GET | `/api/sample.csv` | 周期筛查虚构样例数据 |
| POST | `/api/crossmatch` | 交叉证认(multipart:`candidate_file`、`catalog_file`、`max_sep_arcsec`)→ 200 / 400 / 409 |
| GET | `/api/sample_candidates.csv` | 关联演示候选样例(含 job_id 占位符) |
| GET | `/api/sample_catalog.csv` | 关联演示星表样例 |

示例:

```bash
curl -F "file=@sample/sample_star.csv" -F "period_min=2" -F "period_max=10" \
     http://localhost:8000/api/jobs

# 交叉证认(先把候选 CSV 里的 job_id 占位符换成真实已完成批号)
curl -F "candidate_file=@sample/candidates_sample.csv" \
     -F "catalog_file=@sample/catalog_sample.csv" \
     -F "max_sep_arcsec=1.5" \
     http://localhost:8000/api/crossmatch
```

## 目录结构

```
app/main.py                FastAPI 主程序:校验、任务队列、子进程编排、收尾出图、关联接口
app/worker.py              计算子进程:Lomb–Scargle 功率谱 + 置换检验
app/matching.py            交叉证认:CSV 校验、球面距离、最小费用最大流、同等最优枚举
app/static/index.html      网页(原生 JS,无外部依赖)
sample/make_sample.py      虚构测光样例生成脚本(固定种子,可复现)
sample/sample_star.csv     已生成的测光样例(160 行,真周期 4.25 天)
sample/make_crossmatch_sample.py  关联坐标样例生成脚本(纯标准库,可复现)
sample/candidates_sample.csv      关联候选样例(5 行,job_id 为占位符)
sample/catalog_sample.csv         关联星表样例(5 个源)
Dockerfile / docker-compose.yml
```

重新生成测光样例:`python sample/make_sample.py`
重新生成关联样例:`python sample/make_crossmatch_sample.py`
