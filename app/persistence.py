"""1-dimensional persistence via GF(2) sparse boundary reduction.

The boundary matrix of the cubical complex is reduced column by column in
filtration order over GF(2) — the classic left-to-right reduction, implemented
here from scratch on sparse sets (no external persistent-homology library):

* a column that reduces to zero is *positive*: it creates a homology class;
* a column whose reduced form is non-empty is *negative*: it kills the class
  created by the row index of its lowest entry (``low``).

A 1-dimensional bar is a pair (edge ``e``, square ``s``) where the reduced
boundary of square ``s`` has ``low = e``.  The reduced boundary itself is
reported as the bar's birth cycle: it is a closed even-degree edge chain whose
youngest edge is exactly ``e`` (so the chain exists from the birth threshold
on and no earlier), and it is annihilated precisely when square ``s`` enters —
it is the boundary of the 2-chain accumulated while reducing column ``s``,
and it bounds no 2-chain below ``s``.

The full grid complex is a solid rectangle, hence contractible, so every
1-class that is born also dies: the barcode is finite and every reported bar
has a concrete birth edge and death square.
"""

from __future__ import annotations

from dataclasses import dataclass

from .filtration import SQUARE, VERTEX, Cell, Filtration


@dataclass(frozen=True)
class Bar:
    """One 1-dimensional persistence bar.

    ``cycle_edges`` are the edges of the reported birth cycle (the reduced
    boundary of the death square), sorted by filtration order — the fixed
    order in which the reduction consumed them.  The last edge is always the
    birth edge itself.
    """

    birth_edge: Cell
    death_square: Cell
    cycle_edges: tuple[Cell, ...]

    @property
    def persistence(self) -> int:
        return self.death_square.value - self.birth_edge.value


def compute_h1_bars(filtration: Filtration) -> list[Bar]:
    """Reduce the boundary matrix and return every 1-dimensional bar.

    Bars with ``death_value == birth_value`` (zero persistence) are included
    here; the service layer applies the persistence filters.
    """
    cells = filtration.cells
    boundaries = filtration.boundaries

    low_to_col: dict[int, int] = {}     # lowest row of a reduced column -> that column
    reduced_cols: dict[int, set[int]] = {}  # negative columns: reduced boundary

    bars: list[Bar] = []

    for j, cell in enumerate(cells):
        if cell.kind == VERTEX:
            # Empty boundary: positive, creates an H0 class; irrelevant for H1.
            continue

        col = set(boundaries[j])
        while col:
            low = max(col)
            killer = low_to_col.get(low)
            if killer is None:
                break
            col ^= reduced_cols[killer]

        if not col:
            # Positive cell: a class is born.  For an edge this is an H1
            # birth; nothing needs to be stored — the bar is reported when
            # the killing square is reduced.
            continue

        low = max(col)
        low_to_col[low] = j
        reduced_cols[j] = col
        if cell.kind == SQUARE:
            # `low` is a positive edge (standard invariant of the reduction);
            # pair it with this square.  The reduced boundary `col` is the
            # birth cycle: a closed edge chain, youngest edge == low, killed
            # exactly by this square.
            bars.append(
                Bar(
                    birth_edge=cells[low],
                    death_square=cell,
                    cycle_edges=tuple(cells[i] for i in sorted(col)),
                )
            )

    return bars
