# 岩芯 CT 灰度切片 —— 一维持久同源分析服务（纯后端）

输入一张整数灰度切片（2..96 的方阵/矩形矩阵，取值 0..65535）与最小持续度
门槛，服务自行完成 GF(2) 稀疏边界归约，返回持续度不小于门槛的 **一维 (H1)
持久条码**，并为每条码给出可复算的拓扑证据：**出生边、死亡方格、按固定
归约顺序确定的出生环边集，以及独立路线的复核证书**。

不包含任何前端页面、不依赖在线服务、不调用任何持久同调第三方库
（归约为手写的稀疏 GF(2) 支点算法；Web 框架仅用 FastAPI/uvicorn）。

## 运行

```bash
# 默认宿主机端口 8080
docker compose up --build

# 通过环境变量配置宿主机端口（例如 9090）
API_HOST_PORT=9090 docker compose up --build

# 或使用 .env
cp .env.example .env   # 修改 API_HOST_PORT
docker compose up --build
```

健康检查：

```bash
curl http://localhost:8080/health
# {"status":"ok", ...}
```

镜像内置 `HEALTHCHECK`，Compose 也声明了 `healthcheck`。容器内监听端口由
`PORT` 控制（默认 8000），宿主机映射端口由 `API_HOST_PORT` 控制
（默认 8080）。

## 接口

### `POST /api/v1/analyze`

请求：

```json
{
  "matrix": [[0,0,0],[0,9,0],[0,0,0]],
  "min_persistence": 5
}
```

响应（节选）：

```json
{
  "request": {"rows": 3, "cols": 3, "min_persistence": 5},
  "all_certificates_verified": true,
  "barcodes": [
    {
      "rank": 1,
      "persistence": 9,
      "birth": {
        "threshold": 0,
        "edge": {"kind": "vertical", "row": 1, "col": 2,
                 "endpoints": [[1,2],[2,2]], "filtration_order": 15}
      },
      "death": {
        "threshold": 9,
        "square": {"row": 1, "col": 1, "filtration_order": 24}
      },
      "birth_cycle": {
        "edge_count": 8,
        "ordering": "filtration_total_order_ascending",
        "edges": [ ... 出生环上每条边，按固定归约顺序 ... ]
      },
      "certificate": {
        "even_degree_closed_chain_at_birth": true,
        "birth_edge_is_reduction_pivot": true,
        "endpoints_connected_without_birth_edge": true,
        "absent_at_earlier_threshold": true,
        "death_square_enclosed_by_birth_cycle": true,
        "interior_square_count": 4,
        "interior_unfilled_squares_at_birth": 4,
        "interior_filled_at_death": true,
        "death_square_is_last_filled": true,
        "death_threshold_strictly_after_birth": true
      }
    }
  ]
}
```

非法请求返回 `400/422` 与定位明确的结构化错误，例如：

```json
{
  "error": {
    "code": "INVALID_MATRIX",
    "message": "请求校验失败，共 1 项问题，详见 issues。",
    "issues": [
      {"location": "matrix[0][1]", "reason": "灰度值必须在 0..65535 之间，实际为 70000。",
       "row": 0, "column": 1, "min": 0, "max": 65535}
    ]
  }
}
```

一次请求中的全部问题会同时报出（尺寸、行列定位、非整数/布尔/浮点、
越界、缺字段、JSON 语法错误、空请求体）。

## 算法（可复算说明）

### 1. 单元滤流值

对 m×n 灰度矩阵 `g`：

| 单元 | 记号 | 滤流值 |
|---|---|---|
| 顶点 | V(r,c) | `g[r][c]` |
| 横边 | H(r,c)：(r,c)–(r,c+1) | `max(g[r][c], g[r][c+1])` |
| 竖边 | E(r,c)：(r,c)–(r+1,c) | `max(g[r][c], g[r+1][c])` |
| 方格 | Q(r,c) | 四角灰度最大值 |

这是灰度网格的下星滤流：阈值 t 的子复形由灰度 ≤ t 的顶点及其间已被
“点亮”的边、四角都 ≤ t 的方格组成。

### 2. 全序（固定归约顺序）

滤流值升序；同值时按 **顶点 < 横边 < 竖边 < 方格**，同类再按
**(行, 列) 升序**。该次序保证面严格先于余面（边界矩阵严格上三角）。

### 3. GF(2) 稀疏边界归约

* 边的边界 = 两个端点顶点；方格的边界 = 四条边。
* 因分块结构，H1 的 (出生边, 死亡方格) 配对只需归约 D2 子块
  （行=边，列=方格）。
* 标准从左到右支点（pivot）归约：维护“支点行 → 归约列”，新方格列不断
  异或占用其当前支点的旧列，支点单调下降；列空为正单元，否则支点边与
  该方格配对。
* 死亡方格的归约列就是所报 **出生环边集**：其支点恰为出生边，所有边的
  全序位置不晚于出生边，即在出生阈值已存在。
* 仅保留 **死亡值严格大于出生值** 的环（同值生灭的短命/瞬时配对丢弃），
  再按 **持续度降序 → 出生值升序 → 出生单元全序 → 死亡单元全序** 裁决。

### 4. 每条码的独立复核（`app/certify.py`）

复核不调用归约结果，而是用图搜索独立验证：

1. 出生环每个顶点的关联边数为偶数（GF(2) 下 ∂z=0，偶度闭链）；
2. 出生边是归约列支点；删去出生边后其两端点仍由环上其余边连通
   （出生边一加入就闭合成环）；
3. 在平面方格对偶图上，从图外按“只能穿过尚未出现的原边”洪水泛滥：
   取严格小于出生值的阈值，死亡方格位置仍可被淹没（孔洞在更早阈值
   不存在）；
4. 用出生环边封堵对偶图得到的不可达区域即环内：死亡方格位于环内，
   环内方格按全序最后填入者恰为死亡方格（消解由它触发），出生时环内
   至少一个方格未填入、死亡时全部填入，且死亡值严格大于出生值。

任一复核失败服务返回 `500 CERTIFICATE_FAILED` 而不会发出条码。

### 5. 测试中的第二条独立对照

`tests/test_topology.py` 对 60 个随机小矩阵，在滤流每个前缀上用**独立的
稠密 GF(2) 高斯消元**计算 `dim H1 = |C1| − rank D1 − rank D2`，与归约
配对逐前缀比较，结果必须完全一致。

## 本地开发

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/pip install pytest httpx
.venv/bin/python -m pytest -q
.venv/bin/uvicorn app.main:app --port 8080
```

## 目录

```
app/topology.py    # 滤流、GF(2) 稀疏支点归约、条码与排序
app/certify.py     # 独立对偶图洪水复核与证据
app/validation.py  # 定位明确的结构化请求校验
app/main.py        # FastAPI 路由（/health、/api/v1/analyze）
tests/             # 归约 vs 稠密高斯消元对照、证书、接口与错误用例
Dockerfile         # 非 root 运行 + HEALTHCHECK
docker-compose.yml # API_HOST_PORT / PORT 环境变量配置
```
