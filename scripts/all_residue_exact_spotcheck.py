from __future__ import annotations

from fractions import Fraction as F


def c2(i: int | F) -> F:
    return F(i * (i - 1), 2)


def c3(i: int | F) -> F:
    return F(i * (i - 1) * (i - 2), 6)


def cut_coeff(a: int, k: int, i: int) -> F:
    kk = k - 2
    alpha = (2 * (a - 2)) % kk
    if i < k:
        return c2(i) * F(i - 2 - alpha, kk - alpha)
    return c2(i)


def cut_rhs(a: int, k: int) -> F:
    kk = k - 2
    alpha = (2 * (a - 2)) % kk
    return c2(a) * F(2 * (a - 2) - alpha, kk)


def solve_linear(matrix: list[list[F]], rhs: list[F]) -> list[F]:
    n = len(rhs)
    aug = [row[:] + [rhs[i]] for i, row in enumerate(matrix)]
    for col in range(n):
        pivot = next(r for r in range(col, n) if aug[r][col] != 0)
        aug[col], aug[pivot] = aug[pivot], aug[col]
        pv = aug[col][col]
        aug[col] = [x / pv for x in aug[col]]
        for r in range(n):
            if r != col and aug[r][col] != 0:
                factor = aug[r][col]
                aug[r] = [aug[r][j] - factor * aug[col][j] for j in range(n + 1)]
    return [aug[i][-1] for i in range(n)]


def right_shift(rho: int) -> int:
    return 1 if rho <= 4 else 2


def candidate_value(q: int, rho: int, c: int, b: int) -> tuple[F, dict[int, F]]:
    a = (q * c + rho) // 2
    t = (c - q + 2) // (q + 1)
    h = right_shift(rho)
    sizes = [c - t, c - t + 1, c + h, c + h + 1]
    matrix = [
        [F(1) for _ in sizes],
        [c3(i) for i in sizes],
        [cut_coeff(a, c - t + 1, i) for i in sizes],
        [cut_coeff(a, c + h - 1, i) for i in sizes],
    ]
    rhs = [F(b), 2 * c3(a), cut_rhs(a, c - t + 1), cut_rhs(a, c + h - 1)]
    weights = solve_linear(matrix, rhs)
    return sum(F(sizes[i]) * weights[i] for i in range(4)), dict(zip(sizes, weights))


def es_value(a: int, b: int) -> tuple[F, dict[int, F]]:
    sizes = list(range(2, a + 1))
    nv = len(sizes)
    rows = [
        ([F(1)] * nv, F(b)),
        ([c3(i) for i in sizes], 2 * c3(a)),
    ]
    for k in range(3, a + 1):
        rows.append(([cut_coeff(a, k, i) for i in sizes], cut_rhs(a, k)))
    m = len(rows)
    obj = [F(i) for i in sizes]
    ncol = nv + m
    table = [
        rows[r][0][:] + [F(1) if j == r else F(0) for j in range(m)] + [rows[r][1]]
        for r in range(m)
    ]
    cost = obj[:] + [F(0)] * m
    basis = [nv + r for r in range(m)]
    while True:
        enter = next((j for j in range(ncol) if cost[j] > 0), -1)
        if enter == -1:
            break
        leave = -1
        best = None
        for r in range(m):
            if table[r][enter] > 0:
                ratio = table[r][ncol] / table[r][enter]
                if best is None or ratio < best or (ratio == best and basis[r] < basis[leave]):
                    best = ratio
                    leave = r
        if leave == -1:
            raise RuntimeError("unbounded")
        pivot = table[leave][enter]
        table[leave] = [x / pivot for x in table[leave]]
        for r in range(m):
            if r != leave and table[r][enter] != 0:
                factor = table[r][enter]
                table[r] = [table[r][j] - factor * table[leave][j] for j in range(ncol + 1)]
        if cost[enter] != 0:
            factor = cost[enter]
            cost = [cost[j] - factor * table[leave][j] for j in range(ncol)]
        basis[leave] = enter
    support: dict[int, F] = {}
    value = F(0)
    for r in range(m):
        if basis[r] < nv and table[r][ncol] != 0:
            support[sizes[basis[r]]] = table[r][ncol]
            value += obj[basis[r]] * table[r][ncol]
    return value, support


def midpoint_b(q: int, rho: int, c: int) -> int:
    # For rho>=5 the right pair shifts, so the clean closed roots below are
    # no longer the h=1 symbolic roots.  Use a small exact search for an integer
    # b with all candidate weights nonnegative.
    t = (c - q + 2) // (q + 1)
    h = right_shift(rho)
    a = (q * c + rho) // 2
    good = []
    for b in range(1, max(2 * a, 10)):
        _, weights = candidate_value(q, rho, c, b)
        if all(v >= 0 for v in weights.values()):
            good.append(b)
    if not good:
        raise RuntimeError(f"no integer candidate band for q={q}, rho={rho}, c={c}, t={t}, h={h}")
    return good[len(good) // 2]


def main() -> None:
    samples = [(5, 4, 26), (6, 4, 36), (7, 2, 40), (8, 0, 40), (10, 8, 50)]
    for q, rho, c in samples:
        if (q * c + rho) % 2:
            continue
        a = (q * c + rho) // 2
        b = midpoint_b(q, rho, c)
        cand, weights = candidate_value(q, rho, c, b)
        lp, support = es_value(a, b)
        rel = tuple(i - c for i in support)
        print(f"q={q} rho={rho} c={c} a={a} b={b}")
        print(f"  candidate == LP: {cand == lp}; value={lp}")
        print(f"  LP support rel={rel}")
        print(f"  candidate weights={{{', '.join(f'{k}: {v}' for k, v in weights.items() if v)}}}")


if __name__ == "__main__":
    main()
