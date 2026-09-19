"""条码的独立拓扑复核（与边界归约不同的算法路线，供“可复算证据”使用）。

对每条返回的条码独立验证：

1. **偶度闭链**：出生环在每个顶点的关联边数均为偶数（GF(2) 下 ∂z = 0）。
2. **出生时成环**：出生边是归约列在全序中的最大元（支点），环上其余边
   在全序中更早；删去出生边后，其两端点仍由环上其余边连通——加入出生
   边恰好把一条已存在的路径闭合成环。
3. **更早阈值不存在**：在平面方格对偶图上从图外洪水泛滥，只能穿过在该
   阈值“尚未出现”的原边（已出现的边构成堤坝）。取严格小于出生值的阈
   值，死亡方格所在位置必须仍能被图外淹没（孔洞尚未形成）。
4. **死亡方格触发消解**：用出生环边封堵对偶图得到的不可达区域即环内；
   死亡方格必须位于环内，且环内方格的滤流全序最大者恰为死亡方格
   （它是最后填入的方格，消解由它触发）；出生阈值下环内至少还有一个
   方格未填入（孔洞真实存在），死亡阈值下全部填入（环成为这些方格之
   和的边界），并且死亡值严格大于出生值。
"""

from dataclasses import dataclass
from typing import Dict, List, Sequence, Set, Tuple

from .topology import (
    DIM_HORIZONTAL,
    DIM_SQUARE,
    DIM_VERTICAL,
    Barcode,
    Cell,
)


@dataclass(frozen=True)
class Certificate:
    cycle_edge_count: int
    even_degree_vertex_count: int
    interior_square_count: int
    interior_unfilled_at_birth: int

    def as_dict(self) -> dict:
        return {
            "even_degree_closed_chain_at_birth": True,
            "cycle_edge_count": self.cycle_edge_count,
            "even_degree_vertex_count": self.even_degree_vertex_count,
            "birth_edge_is_reduction_pivot": True,
            "endpoints_connected_without_birth_edge": True,
            "absent_at_earlier_threshold": True,
            "death_square_enclosed_by_birth_cycle": True,
            "interior_square_count": self.interior_square_count,
            "interior_unfilled_squares_at_birth": self.interior_unfilled_at_birth,
            "interior_filled_at_death": True,
            "death_square_is_last_filled": True,
            "death_threshold_strictly_after_birth": True,
        }


class CertificateError(AssertionError):
    """复核失败：意味着归约实现或滤流顺序有误。"""


def _edge_endpoints(cell: Cell) -> Tuple[Tuple[int, int], Tuple[int, int]]:
    if cell.dim == DIM_HORIZONTAL:  # 横边：(r,c) → (r,c+1)
        return (cell.row, cell.col), (cell.row, cell.col + 1)
    # 竖边：(r,c) → (r+1,c)
    return (cell.row, cell.col), (cell.row + 1, cell.col)


def _dual_flood(
    m: int,
    n: int,
    blocked_edges: Set[int],
    index_of: Dict[Tuple[int, int, int], int],
) -> Set[Tuple[int, int]]:
    """从图外对偶节点 BFS。

    对偶节点为 (m-1)*(n-1) 个方格；每条原边分隔两个相邻方格（边界边的
    另一侧是图外）。原边在 ``blocked_edges`` 中时视为堤坝不可穿过，否则
    两侧方格连通。返回图外不可达的方格坐标集合（堤坝围出的区域）。
    """
    sq_rows, sq_cols = m - 1, n - 1
    total = sq_rows * sq_cols
    outer = total
    adj: List[List[int]] = [[] for _ in range(total + 1)]

    def node(rr: int, cc: int) -> int:
        if rr < 0 or rr >= sq_rows or cc < 0 or cc >= sq_cols:
            return outer
        return rr * sq_cols + cc

    def link(edge_index: int, a: int, b: int) -> None:
        if a == b or edge_index in blocked_edges:
            return
        adj[a].append(b)
        adj[b].append(a)

    # 横边 H(r,c) 分隔上侧 Q(r-1,c) 与下侧 Q(r,c)
    for rr in range(m):
        for cc in range(n - 1):
            link(index_of[(DIM_HORIZONTAL, rr, cc)], node(rr - 1, cc), node(rr, cc))
    # 竖边 V(r,c) 分隔左侧 Q(r,c-1) 与右侧 Q(r,c)
    for rr in range(m - 1):
        for cc in range(n):
            link(index_of[(DIM_VERTICAL, rr, cc)], node(rr, cc - 1), node(rr, cc))

    seen = {outer}
    stack = [outer]
    while stack:
        cur = stack.pop()
        for nxt in adj[cur]:
            if nxt not in seen:
                seen.add(nxt)
                stack.append(nxt)

    return {
        (rr, cc)
        for rr in range(sq_rows)
        for cc in range(sq_cols)
        if rr * sq_cols + cc not in seen
    }


def verify_barcode(
    g: Sequence[Sequence[int]],
    cells: Sequence[Cell],
    barcode: Barcode,
) -> Certificate:
    """对单条条码执行全部独立检查，失败抛 :class:`CertificateError`。"""
    m = len(g)
    n = len(g[0])

    index_of: Dict[Tuple[int, int, int], int] = {
        (c.dim, c.row, c.col): c.index for c in cells
    }
    by_index: Dict[int, Cell] = {c.index: c for c in cells}

    cycle: List[Cell] = [by_index[i] for i in barcode.cycle_indices]
    cycle_set: Set[int] = set(barcode.cycle_indices)
    problems: List[str] = []

    def fail_if(condition: bool, message: str) -> None:
        if condition:
            problems.append(message)

    # ---- 1. 偶度闭链 -----------------------------------------------------
    degree: Dict[Tuple[int, int], int] = {}
    adjacency: Dict[Tuple[int, int], List[Tuple[Tuple[int, int], int]]] = {}
    for edge in cycle:
        if edge.dim not in (DIM_HORIZONTAL, DIM_VERTICAL):
            fail_if(True, f"出生环含非边单元（dim={edge.dim}）")
            continue
        u, v = _edge_endpoints(edge)
        degree[u] = degree.get(u, 0) + 1
        degree[v] = degree.get(v, 0) + 1
        adjacency.setdefault(u, []).append((v, edge.index))
        adjacency.setdefault(v, []).append((u, edge.index))

    odd = sorted(vtx for vtx, d in degree.items() if d % 2 != 0)
    fail_if(bool(odd), f"存在奇度顶点，不是偶度闭链：{odd[:4]}（共 {len(odd)} 个）")

    # ---- 2. 支点与出生路径 ----------------------------------------------
    birth_cell = by_index[barcode.birth_index]
    death_cell = by_index[barcode.death_index]
    fail_if(birth_cell.dim not in (DIM_HORIZONTAL, DIM_VERTICAL), "出生单元不是边")
    fail_if(death_cell.dim != DIM_SQUARE, "死亡单元不是方格")
    fail_if(max(cycle_set) != barcode.birth_index, "出生边不是归约列支点（全序最大元）")
    fail_if(
        any(by_index[i].value > barcode.birth_value for i in cycle_set),
        "环上存在出生阈值尚未加入的边",
    )

    endpoints_connected = False
    if birth_cell.dim in (DIM_HORIZONTAL, DIM_VERTICAL):
        bu, bv = _edge_endpoints(birth_cell)
        seen = {bu}
        stack = [bu]
        while stack:
            cur = stack.pop()
            for nxt, eidx in adjacency.get(cur, []):
                if eidx == barcode.birth_index or nxt in seen:
                    continue
                seen.add(nxt)
                stack.append(nxt)
        endpoints_connected = bv in seen
    fail_if(not endpoints_connected, "删去出生边后其两端点在环上不连通")

    # ---- 3. 更早阈值（严格小于出生值）下孔洞不存在 ----------------------
    earlier_edges: Set[int] = {
        c.index
        for c in cells
        if c.dim in (DIM_HORIZONTAL, DIM_VERTICAL) and c.value < barcode.birth_value
    }
    enclosed_before = _dual_flood(m, n, earlier_edges, index_of)
    dpos = (death_cell.row, death_cell.col)
    fail_if(dpos in enclosed_before, "孔洞在更早阈值已存在（出生边不是真正的出生单元）")

    # ---- 4. 环内区域与死亡方格 ------------------------------------------
    interior = _dual_flood(m, n, cycle_set, index_of)
    fail_if(not interior, "出生环未围出任何方格区域")
    fail_if(dpos not in interior, "死亡方格不在出生环围出的区域内")

    interior_indices = [index_of[(DIM_SQUARE, rr, cc)] for rr, cc in interior]
    last_filled = by_index[max(interior_indices)]
    fail_if(
        max(interior_indices) != barcode.death_index,
        f"区域内最后填入的方格是 ({last_filled.row},{last_filled.col})，"
        f"而非所报死亡方格 ({death_cell.row},{death_cell.col})",
    )
    unfilled_at_birth = [
        (rr, cc)
        for rr, cc in interior
        if index_of[(DIM_SQUARE, rr, cc)] > barcode.birth_index
    ]
    fail_if(not unfilled_at_birth, "出生阈值下环内已被方格填满，孔洞不存在")
    fail_if(
        barcode.death_value <= barcode.birth_value,
        "死亡值不严格大于出生值",
    )

    if problems:
        raise CertificateError(
            f"条码 birth=({birth_cell.row},{birth_cell.col}) "
            f"death=({death_cell.row},{death_cell.col})：" + "；".join(problems)
        )

    cert = Certificate(
        cycle_edge_count=len(cycle),
        even_degree_vertex_count=len(degree),
        interior_square_count=len(interior),
        interior_unfilled_at_birth=len(unfilled_at_birth),
    )
    return cert
