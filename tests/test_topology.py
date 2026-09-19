"""拓扑核心测试：手工例子 + 随机矩阵对照独立稠密 GF(2) 秩计算。

独立对照路线（与边界归约完全不同）：对滤流全序的每个前缀子复形，用稠密
GF(2) 高斯消元计算

    dim H1 = |C1| - rank(D1) - rank(D2)

并与“出生 ≤ 前缀 < 死亡”的条码条数逐前缀比较，二者必须处处相等。
"""

import random

from app.certify import verify_barcode
from app.topology import (
    DIM_HORIZONTAL,
    DIM_SQUARE,
    DIM_VERTICAL,
    DIM_VERTEX,
    build_filtration,
    compute_barcodes,
    select_barcodes,
)


def _gf2_rank(rows, ncols):
    """对若干个位向量（int 位掩码）做 GF(2) 高斯消元求秩。"""
    basis = {}
    for vec in rows:
        x = vec
        while x:
            p = x.bit_length() - 1
            if p in basis:
                x ^= basis[p]
            else:
                basis[p] = x
                break
    return len(basis)


def _betti1_prefixes(g):
    """返回每个全序前缀位置 i（已加入 cells[0..i]）的 dim H1 序列。"""
    cells, columns = build_filtration(g)
    n = len(cells)

    # 统一的全局位编号：行/列都用全序下标；构造 D1（边->顶点）与 D2（方格->边）
    betti = []
    for i in range(n):
        present = [c for c in cells if c.index <= i]
        present_idx = {c.index for c in present}
        c1 = [c for c in present if c.dim in (DIM_HORIZONTAL, DIM_VERTICAL)]
        d1_rows = []
        for c in c1:
            mask = 0
            for f in columns[c.index]:
                mask |= 1 << f
            d1_rows.append(mask)
        d2_rows = []
        for c in present:
            if c.dim != DIM_SQUARE:
                continue
            mask = 0
            for f in columns[c.index]:
                if f in present_idx:
                    mask |= 1 << f
            d2_rows.append(mask)
        rank_d1 = _gf2_rank(d1_rows, n)
        rank_d2 = _gf2_rank(d2_rows, n)
        betti.append(len(c1) - rank_d1 - rank_d2)
    return cells, betti


def _alive_counts(cells, pairs):
    """归约配对给出的每个前缀存活 H1 条数。"""
    alive = [0] * len(cells)
    for b, d, _ in pairs:
        for i in range(b, d):
            alive[i] += 1
    return alive


def test_single_square_grid_has_no_positive_barcode():
    # 2x2 只有一个方格，任何灰度分布都不可能产生真孔洞
    for g in (
        [[0, 0], [0, 0]],
        [[0, 65535], [65535, 0]],
        [[5, 9], [1, 2]],
    ):
        barcodes = compute_barcodes(g)
        assert all(b.death_value == b.birth_value for b in barcodes)
        assert select_barcodes(barcodes, 0) == []


def test_handcomputed_ring_with_high_center():
    # 3x3，中心高灰度：外圈 8 条边在阈值 0 闭合成环，
    # 4 个方格都在阈值 9 填入，唯一孔洞持续度 9。
    g = [
        [0, 0, 0],
        [0, 9, 0],
        [0, 0, 0],
    ]
    barcodes = select_barcodes(compute_barcodes(g), 0)
    assert len(barcodes) == 1
    b = barcodes[0]
    assert b.birth_value == 0
    assert b.death_value == 9
    assert b.persistence == 9
    # 全序裁决：竖边 V(1,2) 是外圈环上最后加入的边
    assert (b.birth_dim, b.birth_row, b.birth_col) == (DIM_VERTICAL, 1, 2)
    # 四个方格按 (r,c) 顺序归约，Q(1,1) 触发消解
    assert (b.death_row, b.death_col) == (1, 1)
    # 出生环恰为外圈 8 条边
    assert len(b.cycle_indices) == 8
    cells, _ = build_filtration(g)
    by_index = {c.index: c for c in cells}
    ring = [
        (DIM_HORIZONTAL, 0, 0),
        (DIM_HORIZONTAL, 0, 1),
        (DIM_HORIZONTAL, 2, 0),
        (DIM_HORIZONTAL, 2, 1),
        (DIM_VERTICAL, 0, 0),
        (DIM_VERTICAL, 0, 2),
        (DIM_VERTICAL, 1, 0),
        (DIM_VERTICAL, 1, 2),
    ]
    expected_indices = {
        next(c.index for c in cells if (c.dim, c.row, c.col) == key)
        for key in ring
    }
    assert set(b.cycle_indices) == expected_indices
    # cycle_indices 必须按全序（固定归约顺序）升序
    assert list(b.cycle_indices) == sorted(b.cycle_indices)

    # 独立复核全部通过
    cert = verify_barcode(g, cells, b)
    assert cert.interior_square_count == 4
    assert cert.interior_unfilled_at_birth == 4


def test_flat_matrix_yields_nothing():
    g = [[7] * 5 for _ in range(4)]
    # 同值下顶点→边→方格的次序使所有配对都是零持续度
    assert select_barcodes(compute_barcodes(g), 0) == []


def test_no_barcodes_when_threshold_too_high():
    g = [
        [0, 0, 0],
        [0, 9, 0],
        [0, 0, 0],
    ]
    assert select_barcodes(compute_barcodes(g), 10) == []
    assert len(select_barcodes(compute_barcodes(g), 9)) == 1


def test_tie_break_total_order():
    g = [[1, 2], [3, 4]]
    cells, _ = build_filtration(g)
    keys = [(c.value, c.dim, c.row, c.col) for c in cells]
    assert keys == sorted(keys)
    # 同值时维度次序：顶点 < 横边 < 竖边 < 方格
    # g=[[1,2],[3,4]] 中值 4 的单元：顶点(1,1)、H(1,0)、V(0,1)、Q(0,0)
    group = [(c.dim, c.row, c.col) for c in cells if c.value == 4]
    assert group == sorted(group, key=lambda x: (x[0], x[1], x[2]))
    assert [x[0] for x in group] == [
        DIM_VERTEX,
        DIM_HORIZONTAL,
        DIM_VERTICAL,
        DIM_SQUARE,
    ]


def test_reduction_matches_independent_betti_random():
    rng = random.Random(20260919)
    for trial in range(60):
        m = rng.randint(2, 6)
        n = rng.randint(2, 6)
        # 偏置采样制造孔洞：少量高值点
        g = [
            [rng.choice([0, 0, 0, 1, 5, 9]) for _ in range(n)]
            for _ in range(m)
        ]
        cells, betti = _betti1_prefixes(g)
        raw = compute_barcodes(g)
        pairs = [(b.birth_index, b.death_index, set(b.cycle_indices)) for b in raw]
        alive = _alive_counts(cells, pairs)
        # 前缀 0..len-1 上，归约配对数与独立秩计算的 Betti1 完全一致
        assert alive == betti, f"mismatch on trial {trial}: {g}"
        # 每条正持续度条码都必须通过独立复核
        for b in raw:
            if b.death_value > b.birth_value:
                verify_barcode(g, cells, b)


def test_sort_rule_descending_persistence():
    # 构造两个不同持续度的孔洞：5x5 中两个高值中心
    g = [[0] * 5 for _ in range(5)]
    g[1][1] = 3
    g[3][3] = 8
    barcodes = select_barcodes(compute_barcodes(g), 0)
    persistences = [b.persistence for b in barcodes]
    assert persistences == sorted(persistences, reverse=True)
    assert persistences[0] >= 8 - 0


def test_min_persistence_boundary_is_inclusive():
    g = [
        [0, 0, 0],
        [0, 5, 0],
        [0, 0, 0],
    ]
    assert len(select_barcodes(compute_barcodes(g), 5)) == 1
    assert select_barcodes(compute_barcodes(g), 6) == []


def test_short_lived_noise_loop_filtered_but_true_void_kept():
    # 两个隔离高值中心：值 1 的“细裂隙噪声环”（持续度 1）与值 9 的
    # 真孔洞（持续度 9）。门槛 5 必须只留下真孔洞。
    g = [[0] * 5 for _ in range(5)]
    g[1][1] = 1
    g[3][3] = 9

    all_bars = select_barcodes(compute_barcodes(g), 0)
    assert sorted(b.persistence for b in all_bars) == [1, 9]

    kept = select_barcodes(compute_barcodes(g), 5)
    assert [b.persistence for b in kept] == [9]
    # 噪声环自身也是合法条码：通过同样的独立复核
    noise = next(b for b in all_bars if b.persistence == 1)
    cells, _ = build_filtration(g)
    verify_barcode(g, cells, noise)
