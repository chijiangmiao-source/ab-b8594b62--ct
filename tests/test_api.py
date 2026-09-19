"""HTTP 接口测试：健康检查、成功响应契约、结构化校验错误、排序与证据。"""

import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)


def test_health_ok():
    resp = client.get("/health")
    assert resp.status_code == 200
    body = resp.json()
    assert body["status"] == "ok"
    assert body["service"]


def test_root_is_json_not_html():
    resp = client.get("/")
    assert resp.status_code == 200
    assert resp.headers["content-type"].startswith("application/json")
    assert "analyze" in resp.json()["endpoints"]


def test_analyze_simple_hole_contract():
    g = [
        [0, 0, 0],
        [0, 9, 0],
        [0, 0, 0],
    ]
    resp = client.post("/api/v1/analyze", json={"matrix": g, "min_persistence": 1})
    assert resp.status_code == 200, resp.text
    body = resp.json()

    assert body["request"] == {"rows": 3, "cols": 3, "min_persistence": 1}
    assert body["all_certificates_verified"] is True
    assert body["stats"]["returned_barcodes"] == 1
    assert body["stats"]["zero_persistence_pairs_dropped"] >= 1

    bar = body["barcodes"][0]
    assert bar["rank"] == 1
    assert bar["persistence"] == 9
    assert bar["birth"]["threshold"] == 0
    assert bar["birth"]["edge"]["kind"] == "vertical"
    assert bar["birth"]["edge"]["row"] == 1
    assert bar["birth"]["edge"]["col"] == 2
    assert bar["death"]["threshold"] == 9
    assert bar["death"]["square"]["row"] == 1
    assert bar["death"]["square"]["col"] == 1

    edges = bar["birth_cycle"]["edges"]
    assert bar["birth_cycle"]["edge_count"] == 8
    assert len(edges) == 8
    # 出生环边集按固定归约顺序（滤流全序）给出
    orders = [e["filtration_order"] for e in edges]
    assert orders == sorted(orders)
    # 出生边是出生环中全序最大者（支点）
    assert max(orders) == bar["birth"]["edge"]["filtration_order"]

    cert = bar["certificate"]
    assert cert["even_degree_closed_chain_at_birth"] is True
    assert cert["absent_at_earlier_threshold"] is True
    assert cert["death_square_enclosed_by_birth_cycle"] is True
    assert cert["interior_square_count"] == 4
    assert cert["interior_unfilled_squares_at_birth"] == 4
    assert cert["death_threshold_strictly_after_birth"] is True


def test_min_persistence_inclusive_boundary():
    g = [
        [0, 0, 0],
        [0, 5, 0],
        [0, 0, 0],
    ]
    r1 = client.post("/api/v1/analyze", json={"matrix": g, "min_persistence": 5})
    assert r1.json()["stats"]["returned_barcodes"] == 1
    r2 = client.post("/api/v1/analyze", json={"matrix": g, "min_persistence": 6})
    assert r2.json()["stats"]["returned_barcodes"] == 0


def test_response_sorted_by_persistence_descending():
    g = [[0] * 5 for _ in range(5)]
    g[1][1] = 3
    g[3][3] = 8
    resp = client.post("/api/v1/analyze", json={"matrix": g, "min_persistence": 0})
    persistences = [b["persistence"] for b in resp.json()["barcodes"]]
    assert persistences == sorted(persistences, reverse=True)


def test_empty_body_is_structured_error():
    resp = client.post("/api/v1/analyze", content=b"", headers={"content-type": "application/json"})
    assert resp.status_code == 400
    body = resp.json()
    assert body["error"]["code"] == "EMPTY_BODY"
    assert body["error"]["issues"][0]["location"] == "$"


def test_malformed_json_is_structured_error():
    resp = client.post("/api/v1/analyze", content="{not json", headers={"content-type": "application/json"})
    assert resp.status_code == 400
    assert resp.json()["error"]["code"] == "MALFORMED_JSON"


def test_missing_matrix_field():
    resp = client.post("/api/v1/analyze", json={"min_persistence": 1})
    assert resp.status_code == 422
    body = resp.json()
    assert body["error"]["code"] == "MISSING_FIELD"
    assert body["error"]["issues"][0]["location"] == "matrix"


def test_matrix_not_array():
    resp = client.post("/api/v1/analyze", json={"matrix": "x", "min_persistence": 1})
    assert resp.status_code == 422
    assert resp.json()["error"]["issues"][0]["location"] == "matrix"


def test_row_count_out_of_range_is_located():
    g = [[0, 0]]  # 1 行 < 2
    resp = client.post("/api/v1/analyze", json={"matrix": g, "min_persistence": 0})
    assert resp.status_code == 422
    locations = {i["location"] for i in resp.json()["error"]["issues"]}
    assert "matrix.rows" in locations


def test_97_columns_rejected():
    g = [[0] * 97 for _ in range(2)]
    resp = client.post("/api/v1/analyze", json={"matrix": g, "min_persistence": 0})
    assert resp.status_code == 422
    locations = {i["location"] for i in resp.json()["error"]["issues"]}
    assert "matrix[0].length" in locations


def test_ragged_matrix_pinpoints_row():
    g = [[0, 0, 0], [0, 0], [0, 0, 0]]
    resp = client.post("/api/v1/analyze", json={"matrix": g, "min_persistence": 0})
    assert resp.status_code == 422
    issues = resp.json()["error"]["issues"]
    assert any(i["location"] == "matrix[1].length" and i["expected"] == 3 for i in issues)


def test_value_out_of_range_pinpoints_cell():
    g = [[0, 70000], [0, 0]]
    resp = client.post("/api/v1/analyze", json={"matrix": g, "min_persistence": 0})
    assert resp.status_code == 422
    issues = resp.json()["error"]["issues"]
    hit = [i for i in issues if i["location"] == "matrix[0][1]"]
    assert hit and hit[0]["max"] == 65535


def test_negative_value_pinpoints_cell():
    g = [[-1, 0], [0, 0]]
    resp = client.post("/api/v1/analyze", json={"matrix": g, "min_persistence": 0})
    assert resp.status_code == 422
    assert any(i["location"] == "matrix[0][0]" for i in resp.json()["error"]["issues"])


def test_float_and_bool_rejected_as_non_integer():
    g = [[1.5, 0], [True, 0]]
    resp = client.post("/api/v1/analyze", json={"matrix": g, "min_persistence": 0})
    assert resp.status_code == 422
    locations = {i["location"] for i in resp.json()["error"]["issues"]}
    assert "matrix[0][0]" in locations
    assert "matrix[1][0]" in locations


def test_min_persistence_out_of_range_located():
    g = [[0, 0], [0, 0]]
    resp = client.post("/api/v1/analyze", json={"matrix": g, "min_persistence": 65536})
    assert resp.status_code == 422
    locations = {i["location"] for i in resp.json()["error"]["issues"]}
    assert "min_persistence" in locations


def test_multiple_errors_reported_together():
    # 行长度越界 + 单元格越界 + min_persistence 缺失，一次性全部定位
    g = [[0, 0], [0, 99999]]
    resp = client.post("/api/v1/analyze", json={"matrix": g})
    assert resp.status_code == 422
    locations = {i["location"] for i in resp.json()["error"]["issues"]}
    assert "matrix[1][1]" in locations
    assert "min_persistence" in locations


def test_96x96_extremes_accepted_and_fast():
    # 边界尺寸 + 极端值：96x96 约 36k 单元，必须秒级完成
    g = [[((r * 31 + c * 17) % 65536) for c in range(96)] for r in range(96)]
    resp = client.post("/api/v1/analyze", json={"matrix": g, "min_persistence": 0})
    assert resp.status_code == 200, resp.text
    body = resp.json()
    assert body["request"] == {"rows": 96, "cols": 96, "min_persistence": 0}
    assert body["stats"]["elapsed_ms"] < 10000
