"""岩芯 CT 灰度切片一维持久同源分析 —— 纯后端 HTTP 服务。

接口
----
* ``GET  /health``        健康检查（容器编排探活）
* ``GET  /``              服务信息（JSON）
* ``POST /api/v1/analyze`` 提交灰度矩阵与最小持续度，返回 H1 持久条码及
                           出生边、死亡方格、出生环边集与独立复核证据

服务不提供任何页面（docs/redoc 已关闭），不访问任何在线服务。
"""

import json
import logging
import time
from typing import Any, Dict

from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse

from .certify import CertificateError, verify_barcode
from .topology import (
    DIM_HORIZONTAL,
    DIM_VERTICAL,
    build_filtration,
    compute_barcodes,
    select_barcodes,
)
from .validation import ValidationError, validate_payload

logger = logging.getLogger("corect")
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s %(message)s")

SERVICE_NAME = "core-ct-persistence-service"
SERVICE_VERSION = "1.0.0"
STARTED_AT = time.time()

app = FastAPI(
    title=SERVICE_NAME,
    version=SERVICE_VERSION,
    docs_url=None,
    redoc_url=None,
    openapi_url="/openapi.json",
)


def _edge_kind(dim: int) -> str:
    if dim == DIM_HORIZONTAL:
        return "horizontal"
    if dim == DIM_VERTICAL:
        return "vertical"
    return "vertex" if dim == 0 else "square"


def _edge_payload(edge: Any) -> Dict[str, Any]:
    if edge.dim == DIM_HORIZONTAL:
        endpoints = [[edge.row, edge.col], [edge.row, edge.col + 1]]
    else:
        endpoints = [[edge.row, edge.col], [edge.row + 1, edge.col]]
    return {
        "kind": _edge_kind(edge.dim),
        "row": edge.row,
        "col": edge.col,
        "endpoints": endpoints,
        "threshold": edge.value,
        "filtration_order": edge.index,
    }


def _barcode_payload(barcode: Any, cells_by_index: Dict[int, Any], g: Any, cells: Any) -> Dict[str, Any]:
    birth_edge = cells_by_index[barcode.birth_index]
    death_square = cells_by_index[barcode.death_index]

    cycle_edges = [cells_by_index[i] for i in barcode.cycle_indices]
    certificate = verify_barcode(g, cells, barcode)

    return {
        "persistence": barcode.persistence,
        "birth": {
            "threshold": barcode.birth_value,
            "edge": _edge_payload(birth_edge),
        },
        "death": {
            "threshold": barcode.death_value,
            "square": {
                "row": death_square.row,
                "col": death_square.col,
                "vertices": [
                    [death_square.row, death_square.col],
                    [death_square.row, death_square.col + 1],
                    [death_square.row + 1, death_square.col],
                    [death_square.row + 1, death_square.col + 1],
                ],
                "threshold": death_square.value,
                "filtration_order": death_square.index,
            },
        },
        "birth_cycle": {
            "edge_count": len(barcode.cycle_indices),
            # 边集按固定归约顺序（滤流全序）升序给出，可独立复算
            "ordering": "filtration_total_order_ascending",
            "edges": [_edge_payload(e) for e in cycle_edges],
        },
        "certificate": certificate.as_dict(),
    }


@app.get("/health")
async def health() -> Dict[str, Any]:
    return {
        "status": "ok",
        "service": SERVICE_NAME,
        "version": SERVICE_VERSION,
        "uptime_seconds": round(time.time() - STARTED_AT, 3),
    }


@app.get("/")
async def root() -> Dict[str, Any]:
    return {
        "service": SERVICE_NAME,
        "version": SERVICE_VERSION,
        "endpoints": {
            "health": "GET /health",
            "analyze": "POST /api/v1/analyze",
            "openapi": "GET /openapi.json",
        },
    }


@app.exception_handler(ValidationError)
async def validation_error_handler(_request: Request, exc: ValidationError) -> JSONResponse:
    return JSONResponse(
        status_code=422,
        content={
            "error": {
                "code": exc.code,
                "message": exc.message,
                "issues": exc.issues,
            }
        },
    )


@app.post("/api/v1/analyze")
async def analyze(request: Request) -> JSONResponse:
    raw = await request.body()
    if not raw:
        return JSONResponse(
            status_code=400,
            content={
                "error": {
                    "code": "EMPTY_BODY",
                    "message": "请求体为空，需要 JSON 对象。",
                    "issues": [{"location": "$", "reason": "empty request body"}],
                }
            },
        )
    try:
        payload = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        return JSONResponse(
            status_code=400,
            content={
                "error": {
                    "code": "MALFORMED_JSON",
                    "message": "请求体不是合法 JSON。",
                    "issues": [
                        {
                            "location": "$",
                            "reason": str(exc),
                        }
                    ],
                }
            },
        )

    g, min_persistence = validate_payload(payload)
    m, n = len(g), len(g[0])

    started = time.perf_counter()
    cells, columns = build_filtration(g)
    all_barcodes = compute_barcodes(g, cells, columns)
    barcodes = select_barcodes(all_barcodes, min_persistence)
    cells_by_index = {c.index: c for c in cells}

    result_barcodes = []
    try:
        for rank, barcode in enumerate(barcodes, start=1):
            item = _barcode_payload(barcode, cells_by_index, g, cells)
            item["rank"] = rank
            result_barcodes.append(item)
    except CertificateError as exc:
        # 归约结果与独立复核不一致属于服务内部错误，绝不静默返回
        logger.exception("barcode certificate verification failed")
        return JSONResponse(
            status_code=500,
            content={
                "error": {
                    "code": "CERTIFICATE_FAILED",
                    "message": "归约结果未通过独立拓扑复核。",
                    "issues": [{"location": "topology", "reason": str(exc)}],
                }
            },
        )

    elapsed_ms = round((time.perf_counter() - started) * 1000.0, 3)
    positive_pairs = sum(1 for b in all_barcodes if b.death_value > b.birth_value)

    return JSONResponse(
        content={
            "service": SERVICE_NAME,
            "version": SERVICE_VERSION,
            "request": {
                "rows": m,
                "cols": n,
                "min_persistence": min_persistence,
            },
            "filtration": {
                "value_rules": {
                    "vertex(r,c)": "g[r][c]",
                    "horizontal_edge(r,c)": "max(g[r][c], g[r][c+1])",
                    "vertical_edge(r,c)": "max(g[r][c], g[r+1][c])",
                    "square(r,c)": "max(g[r][c], g[r][c+1], g[r+1][c], g[r+1][c+1])",
                },
                "tie_break": [
                    "filtration value ascending",
                    "vertex < horizontal_edge < vertical_edge < square",
                    "(row ascending, col ascending)",
                ],
                "total_cells": len(cells),
                "reduction": "GF(2) sparse left-to-right pivot reduction on D2 block",
            },
            "sort_rule": [
                "persistence descending",
                "birth threshold ascending",
                "birth cell filtration order ascending",
                "death cell filtration order ascending",
            ],
            "stats": {
                "all_birth_death_pairs": len(all_barcodes),
                "positive_persistence_pairs": positive_pairs,
                "zero_persistence_pairs_dropped": len(all_barcodes) - positive_pairs,
                "below_threshold_dropped": positive_pairs - len(barcodes),
                "returned_barcodes": len(barcodes),
                "elapsed_ms": elapsed_ms,
            },
            "all_certificates_verified": True,
            "barcodes": result_barcodes,
        }
    )
