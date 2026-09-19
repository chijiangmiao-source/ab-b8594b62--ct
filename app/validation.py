"""分析接口的请求校验：所有错误都带明确定位，供调用方复算与排障。"""

from typing import Any, List, Tuple

MIN_SIDE = 2
MAX_SIDE = 96
MIN_VALUE = 0
MAX_VALUE = 65535


class ValidationError(Exception):
    """请求不合法。``issues`` 中每项至少含 location/reason。"""

    def __init__(self, code: str, message: str, issues: List[dict]):
        super().__init__(message)
        self.code = code
        self.message = message
        self.issues = issues


def _is_strict_int(x: Any) -> bool:
    """JSON 整数校验。bool 是 int 的子类，须显式排除。"""
    return isinstance(x, int) and not isinstance(x, bool)


def _issue(location: str, reason: str, value: Any = None, **extra: Any) -> dict:
    item: dict = {"location": location, "reason": reason}
    if value is not None:
        item["value"] = value
    item.update(extra)
    return item


def validate_payload(payload: Any) -> Tuple[List[List[int]], int]:
    """校验解析后的 JSON 请求体，返回 (矩阵, 最小持续度)。

    矩阵：高宽各 2..96 的矩形整数数组，元素取值 0..65535；
    最小持续度：0..65535 的整数。
    不合法时抛出 :class:`ValidationError`，issues 逐项定位到
    matrix[row][col] / min_persistence / matrix.rows 等。
    """
    issues: List[dict] = []

    if not isinstance(payload, dict):
        raise ValidationError(
            "INVALID_REQUEST",
            "请求体必须是 JSON 对象，包含 matrix 与 min_persistence 字段。",
            [_issue("$", "expected JSON object")],
        )

    if "matrix" not in payload:
        raise ValidationError(
            "MISSING_FIELD",
            "缺少必填字段 matrix。",
            [_issue("matrix", "field is required")],
        )

    matrix = payload["matrix"]
    if not isinstance(matrix, list):
        raise ValidationError(
            "INVALID_MATRIX",
            "matrix 必须是二维整数数组。",
            [_issue("matrix", "expected array of rows", matrix)],
        )

    rows = len(matrix)
    if rows < MIN_SIDE or rows > MAX_SIDE:
        issues.append(
            _issue(
                "matrix.rows",
                f"矩阵高度必须在 {MIN_SIDE}..{MAX_SIDE} 之间，实际为 {rows}。",
                rows,
                actual=rows,
                min=MIN_SIDE,
                max=MAX_SIDE,
            )
        )

    # 先确定列数（取首行长度），用于定位越界元素；无法确定则跳过逐格校验
    width: int = -1
    if matrix and isinstance(matrix[0], list):
        width = len(matrix[0])
        if width < MIN_SIDE or width > MAX_SIDE:
            issues.append(
                _issue(
                    "matrix[0].length",
                    f"矩阵宽度必须在 {MIN_SIDE}..{MAX_SIDE} 之间，实际为 {width}。",
                    width,
                    actual=width,
                    min=MIN_SIDE,
                    max=MAX_SIDE,
                )
            )

    for r, row in enumerate(matrix):
        row_loc = f"matrix[{r}]"
        if not isinstance(row, list):
            issues.append(_issue(row_loc, "该行必须是整数数组。", row))
            continue
        if width >= 0 and len(row) != width:
            issues.append(
                _issue(
                    f"{row_loc}.length",
                    f"矩阵必须为矩形：第 {r} 行长度为 {len(row)}，"
                    f"与首行宽度 {width} 不一致。",
                    len(row),
                    row=r,
                    expected=width,
                    actual=len(row),
                )
            )
        for c, value in enumerate(row):
            if not _is_strict_int(value):
                issues.append(
                    _issue(
                        f"matrix[{r}][{c}]",
                        "灰度值必须是整数。",
                        value,
                        row=r,
                        column=c,
                    )
                )
            elif value < MIN_VALUE or value > MAX_VALUE:
                issues.append(
                    _issue(
                        f"matrix[{r}][{c}]",
                        f"灰度值必须在 {MIN_VALUE}..{MAX_VALUE} 之间，实际为 {value}。",
                        value,
                        row=r,
                        column=c,
                        min=MIN_VALUE,
                        max=MAX_VALUE,
                    )
                )

    # 宽度无法从首行判定时，仍要给出高宽错误
    if width < 0 and not issues:  # pragma: no cover - 已被 matrix.rows 分支覆盖
        issues.append(_issue("matrix[0]", "首行必须是数组。"))

    if "min_persistence" not in payload:
        issues.append(_issue("min_persistence", "字段是必填项。"))
    else:
        mp = payload["min_persistence"]
        if not _is_strict_int(mp):
            issues.append(
                _issue("min_persistence", "最小持续度必须是整数。", mp)
            )
        elif mp < MIN_VALUE or mp > MAX_VALUE:
            issues.append(
                _issue(
                    "min_persistence",
                    f"最小持续度必须在 {MIN_VALUE}..{MAX_VALUE} 之间，实际为 {mp}。",
                    mp,
                    min=MIN_VALUE,
                    max=MAX_VALUE,
                )
            )

    if issues:
        raise ValidationError(
            "INVALID_MATRIX",
            f"请求校验失败，共 {len(issues)} 项问题，详见 issues。",
            issues,
        )

    return matrix, payload["min_persistence"]
