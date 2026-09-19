"""岩芯 CT 灰度切片的一维持久同源分析（纯后端，无持久同调第三方库）。

本模块在二维规则方格网格上构造下星式滤流（filtration），并自行完成
GF(2) 上的稀疏边界矩阵归约，输出 H1 持久条码、出生边、死亡方格以及
按固定归约顺序确定的出生环边集。

单元与滤流值定义
----------------
设输入灰度矩阵 ``g`` 为 m 行 n 列：

* 顶点      V(r, c)：滤流值 ``g[r][c]``
* 横边      H(r, c)：连接 (r, c)-(r, c+1)，滤流值 ``max(g[r][c], g[r][c+1])``
* 竖边      E(r, c)：连接 (r, c)-(r+1, c)，滤流值 ``max(g[r][c], g[r+1][c])``
* 方格      Q(r, c)：四个顶点灰度的最大值

全序（“固定归约顺序”）：滤流值升序；同值时按
顶点(0) < 横边(1) < 竖边(2) < 方格(3)，同类按 (行, 列) 升序。
面单元的滤流值不大于其余面，且同值时维度次序保证面严格排在余面之前，
因此边界矩阵严格上三角。

归约得到的 (出生边 e, 死亡方格 q) 配对中，q 的归约列就是被 q 消掉的
1-闭链，其支点（最大非零行）恰为 e；该边集即出生时的代表环。
"""

from dataclasses import dataclass
from typing import List, Sequence, Tuple

# 维度编号（同时作为同值平局时的维度次序）
DIM_VERTEX = 0
DIM_HORIZONTAL = 1
DIM_VERTICAL = 2
DIM_SQUARE = 3

CellKey = Tuple[int, int, int]  # (dim, row, col)


@dataclass(frozen=True)
class Cell:
    """滤流全序中的一个单元。"""

    index: int
    value: int
    dim: int
    row: int
    col: int


@dataclass(frozen=True)
class Barcode:
    """一条 H1 持久条码及其拓扑证据。"""

    birth_index: int
    death_index: int
    birth_value: int
    death_value: int
    persistence: int
    birth_dim: int
    birth_row: int
    birth_col: int
    death_row: int
    death_col: int
    # 出生环边集，按滤流全序（即固定归约顺序）升序排列的全局单元下标
    cycle_indices: Tuple[int, ...]


def _cell_keys_and_values(
    g: Sequence[Sequence[int]],
) -> List[Tuple[int, int, int, int]]:
    """枚举全部单元，返回 (value, dim_rank, row, col) 列表。"""
    m = len(g)
    n = len(g[0])
    records: List[Tuple[int, int, int, int]] = []

    for r in range(m):
        for c in range(n):
            records.append((g[r][c], DIM_VERTEX, r, c))

    for r in range(m):
        for c in range(n - 1):
            records.append((max(g[r][c], g[r][c + 1]), DIM_HORIZONTAL, r, c))

    for r in range(m - 1):
        for c in range(n):
            records.append((max(g[r][c], g[r + 1][c]), DIM_VERTICAL, r, c))

    for r in range(m - 1):
        for c in range(n - 1):
            value = max(g[r][c], g[r][c + 1], g[r + 1][c], g[r + 1][c + 1])
            records.append((value, DIM_SQUARE, r, c))

    return records


def _faces(key: CellKey) -> Tuple[CellKey, ...]:
    """返回单元的余维一（codimension-1）面。"""
    dim, r, c = key
    if dim == DIM_VERTEX:
        return ()
    if dim == DIM_HORIZONTAL:  # (r,c)-(r,c+1)
        return ((DIM_VERTEX, r, c), (DIM_VERTEX, r, c + 1))
    if dim == DIM_VERTICAL:  # (r,c)-(r+1,c)
        return ((DIM_VERTEX, r, c), (DIM_VERTEX, r + 1, c))
    # 方格的四条边：上边、下边（横边），左边、右边（竖边）
    return (
        (DIM_HORIZONTAL, r, c),
        (DIM_HORIZONTAL, r + 1, c),
        (DIM_VERTICAL, r, c),
        (DIM_VERTICAL, r, c + 1),
    )


def build_filtration(
    g: Sequence[Sequence[int]],
) -> Tuple[List[Cell], List[set]]:
    """按全序排列单元并构造稀疏边界列（每列只含 2 或 4 个非零项）。

    返回 (按全序排列的单元表, 与之一一对应的边界列)；边界列以全局下标
    集合表示 GF(2) 系数，面严格先于余面。
    """
    records = _cell_keys_and_values(g)
    records.sort(key=lambda x: (x[0], x[1], x[2], x[3]))

    index_of = {}
    cells: List[Cell] = []
    for i, (value, dim, r, c) in enumerate(records):
        index_of[(dim, r, c)] = i
        cells.append(Cell(i, value, dim, r, c))

    columns: List[set] = []
    for cell in cells:
        boundary = {index_of[face] for face in _faces((cell.dim, cell.row, cell.col))}
        # 严格上三角自检：所有面必须在全序中更早出现
        if boundary and max(boundary) >= cell.index:  # pragma: no cover - 算法不变量
            raise RuntimeError("filtration order invariant violated")
        columns.append(boundary)

    return cells, columns


def reduce_boundary_matrix(
    cells: Sequence[Cell], columns: Sequence[set]
) -> List[Tuple[int, int, set]]:
    """GF(2) 稀疏边界矩阵的标准从左到右支点归约（仅求 H1 配对）。

    边界矩阵按维度分块：方格（2-单元）列只含边（1-单元）行，而边列的支点
    只可能落在顶点行，永远不会出现在方格列中。因此 (1,2) 维持久配对只需
    归约 D2 子块（行=边，列=方格），无需归约顶点列与边列。

    维护 ``pivot_owner[row] = column``：归约列的最大非零行（支点）被哪一
    列占用。处理新方格列时，只要其支点已被占用，就把占用列（更早加入
    的方格的归约列）整体异或进来，支点严格单调下降直到：

    * 列非空且支点未被占用 → 该支点边与本方格配对（环在此方格处死亡），
      归约列即为该环从出生到死亡沿用的代表 1-闭链，其最大非零行恰为
      出生边；
    * 列归约为空 → 该方格为正单元（产生 H2；盘状网格中不存在）。

    返回所有 (出生边下标, 死亡方格下标, 归约列) 配对。
    """
    pivot_owner: dict = {}
    reduced: dict = {}
    pairs: List[Tuple[int, int, set]] = []

    for cell in cells:
        if cell.dim != DIM_SQUARE:
            continue
        column = set(columns[cell.index])  # 拷贝，后续原地异或
        while column:
            pivot = max(column)
            owner = pivot_owner.get(pivot)
            if owner is None:
                pivot_owner[pivot] = cell.index
                reduced[cell.index] = column
                pairs.append((pivot, cell.index, column))
                break
            column ^= reduced[owner]
        # 列归约为空：正 2-单元（盘状有限网格中不会出现）

    return pairs


def compute_barcodes(
    g: Sequence[Sequence[int]],
    cells: Sequence[Cell] | None = None,
    columns: Sequence[set] | None = None,
) -> List[Barcode]:
    """计算完整 H1 持久条码（含同值生灭的零长度配对，由调用方过滤）。"""
    if cells is None or columns is None:
        cells, columns = build_filtration(g)
    pairs = reduce_boundary_matrix(cells, columns)

    barcodes: List[Barcode] = []
    for birth_index, death_index, cycle_column in pairs:
        birth = cells[birth_index]
        death = cells[death_index]
        # 边界矩阵按维度分块：方格列只会与边行配对
        if birth.dim not in (DIM_HORIZONTAL, DIM_VERTICAL):  # pragma: no cover
            continue
        if death.dim != DIM_SQUARE:  # pragma: no cover
            continue
        barcodes.append(
            Barcode(
                birth_index=birth_index,
                death_index=death_index,
                birth_value=birth.value,
                death_value=death.value,
                persistence=death.value - birth.value,
                birth_dim=birth.dim,
                birth_row=birth.row,
                birth_col=birth.col,
                death_row=death.row,
                death_col=death.col,
                cycle_indices=tuple(sorted(cycle_column)),
            )
        )
    return barcodes


def select_barcodes(
    barcodes: Sequence[Barcode], min_persistence: int
) -> List[Barcode]:
    """过滤零持续度配对并按裁决顺序排序。

    仅保留死亡值严格大于出生值的环，且持续度不小于门槛；排序为
    持续度降序、出生值升序、出生单元全序升序、死亡单元全序升序。
    """
    kept = [b for b in barcodes if b.death_value > b.birth_value and b.persistence >= min_persistence]
    kept.sort(
        key=lambda b: (
            -b.persistence,
            b.birth_value,
            b.birth_index,
            b.death_index,
        )
    )
    return kept
