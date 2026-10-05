"""Verify the r=4 progression certificate used in Theorem thm:r4.

The symbolic checks cover the displayed primal identities and the moving
two-cut dual formulas. Exact rational checks then verify every n=16,...,31
and every primal/dual row at four separated multiples of 30.
"""

from __future__ import annotations

from fractions import Fraction as F

import sympy as sp


def c2(x):
    return x * (x - 1) / 2


def c3(x):
    return x * (x - 1) * (x - 2) / 6


def assert_zero(expression, name: str) -> None:
    value = sp.factor(expression)
    if value != 0:
        raise AssertionError(f"{name}: {value}")


def symbolic_checks() -> None:
    m, n, k, i = sp.symbols("m n k i")
    d = 61 * m**2 + 270 * m + 200
    h = 549 * m**2 + 1130 * m + 800
    a = 242 * m**3 + 605 * m**2 + 450 * m
    e3 = (m + 2) * n - 16 * (m - 1)
    e4 = (4 * m**2 + 30 * m + 50) * n - 125 * m**2 - 375 * m + 500
    sizes = [2 * m / 5 + 1, 2 * m / 5 + 2, m / 2, m / 2 + 1]
    weights = [100 * (m + 5) * e3 / d, 25 * m * e3 / d, 0, -16 * e4 / d]
    f4 = (
        2
        * (m + 2)
        * (9 * m**2 * n + 100 * m**2 + 55 * m * n - 100 * m + 50 * n)
        / d
    )

    assert_zero(sum(weights) - n, "primal edge count")
    assert_zero(sum(weight * c3(size) for weight, size in zip(weights, sizes)) - 2 * c3(m),
                "primal third moment")
    assert_zero(sum(weight * size for weight, size in zip(weights, sizes)) - f4,
                "primal objective")

    alpha_k = 2 * m - 4 * k + 4
    alpha_k1 = 2 * m - 4 * k

    def low_coeff(size, cut, alpha):
        return c2(size) * (size - 2 - alpha) / (cut - 2 - alpha)

    matrix = sp.Matrix(
        [
            [
                1,
                c3(size),
                low_coeff(size, k, alpha_k) if index < 2 else c2(size),
                low_coeff(size, k + 1, alpha_k1) if index < 2 else c2(size),
            ]
            for index, size in enumerate(sizes)
        ]
    )
    rho = [sp.factor(value) for value in matrix.inv() * sp.Matrix(sizes)]
    expected_g = 24 * (11 * m**2 + 170 * m + 200) / (m * (m - 2) * d)
    expected_k = (
        (m - 10) * (5 * k - 2 * m - 6) * (h * k - a)
        / (m * (m - 2) * (2 * m + 5) * d)
    )
    expected_k1 = (
        -(m - 10) * (5 * k - 2 * m - 1) * (h * k - a - h)
        / (m * (m - 2) * (2 * m + 5) * d)
    )
    assert_zero(rho[1] - expected_g, "third-moment multiplier")
    assert_zero(rho[2] - expected_k, "D_k multiplier")
    assert_zero(rho[3] - expected_k1, "D_k+1 multiplier")

    slack_low = sp.factor(
        rho[0]
        + rho[1] * c3(i)
        + rho[2] * low_coeff(i, k, alpha_k)
        + rho[3] * low_coeff(i, k + 1, alpha_k1)
        - i
    )
    expected_low = (
        (2 * m - 5 * i + 5)
        * (2 * m - 5 * i + 10)
        * (18 * m**2 + 29 * m * i + 56 * m + 10 * i + 40)
        / (2 * (2 * m + 5) * d)
    )
    assert_zero(slack_low - expected_low, "low dual slack")

    slack_high = sp.factor(
        rho[0] + rho[1] * c3(i) + (rho[2] + rho[3]) * c2(i) - i
    )
    expected_high = (
        (m - 2 * i)
        * (m - 2 * i + 2)
        * (
            18 * m**3
            + 11 * m**2 * i
            + 74 * m**2
            + 170 * m * i
            - 120 * m
            + 200 * i
            - 200
        )
        / (m * (m - 2) * d)
    )
    assert_zero(slack_high - expected_high, "high dual slack")

    slack_at_k = sp.factor(
        rho[0]
        + rho[1] * c3(k)
        + rho[2] * c2(k)
        + rho[3] * low_coeff(k, k + 1, alpha_k1)
        - k
    )
    expected_at_k = (
        (5 * k - 2 * m - 10)
        * (5 * k - 2 * m - 5)
        * (29 * k * m + 10 * k + 18 * m**2 + 56 * m + 40)
        / (2 * (2 * m + 5) * d)
    )
    assert_zero(slack_at_k - expected_at_k, "intervening dual slack")


def fc2(x: int) -> F:
    return F(x * (x - 1), 2)


def fc3(x: int) -> F:
    return F(x * (x - 1) * (x - 2), 6)


def cut_coefficient(m: int, cut: int, size: int) -> F:
    alpha = (2 * (m - 2)) % (cut - 2)
    if size < cut:
        return fc2(size) * F(size - 2 - alpha, cut - 2 - alpha)
    return fc2(size)


def cut_rhs(m: int, cut: int) -> F:
    alpha = (2 * (m - 2)) % (cut - 2)
    return fc2(m) * F(2 * (m - 2) - alpha, cut - 2)


def solve_linear(matrix: list[list[F]], vector: list[F]) -> list[F]:
    count = len(vector)
    augmented = [row[:] + [vector[index]] for index, row in enumerate(matrix)]
    for column in range(count):
        pivot = next(row for row in range(column, count) if augmented[row][column])
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(count):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor:
                augmented[row] = [
                    augmented[row][entry] - factor * augmented[column][entry]
                    for entry in range(count + 1)
                ]
    return [augmented[index][-1] for index in range(count)]


def exact_full_band_check(m: int) -> None:
    d = 61 * m * m + 270 * m + 200
    h = 549 * m * m + 1130 * m + 800
    a = 242 * m**3 + 605 * m * m + 450 * m
    k = (a + h - 1) // h
    sizes = [2 * m // 5 + 1, 2 * m // 5 + 2, m // 2, m // 2 + 1]
    dual_matrix = [
        [
            F(1),
            fc3(size),
            cut_coefficient(m, k, size),
            cut_coefficient(m, k + 1, size),
        ]
        for size in sizes
    ]
    rho = solve_linear(dual_matrix, [F(size) for size in sizes])
    if any(multiplier < 0 for multiplier in rho[1:]):
        raise AssertionError((m, "negative dual multiplier", rho))
    for size in range(2, m + 1):
        slack = (
            rho[0]
            + rho[1] * fc3(size)
            + rho[2] * cut_coefficient(m, k, size)
            + rho[3] * cut_coefficient(m, k + 1, size)
            - size
        )
        if slack < 0:
            raise AssertionError((m, "negative dual slack", size, slack))

    for n in range(16, 32):
        e3 = (m + 2) * n - 16 * (m - 1)
        e4 = (4 * m * m + 30 * m + 50) * n - 125 * m * m - 375 * m + 500
        weights = [F(100 * (m + 5) * e3, d), F(25 * m * e3, d), F(0), F(-16 * e4, d)]
        if any(weight < 0 for weight in weights):
            raise AssertionError((m, n, "negative primal weight", weights))
        if sum(weights) != n:
            raise AssertionError((m, n, "edge count"))
        if sum(weight * fc3(size) for weight, size in zip(weights, sizes)) != 2 * fc3(m):
            raise AssertionError((m, n, "third moment"))
        for cut in range(3, m + 1):
            lhs = sum(
                weight * cut_coefficient(m, cut, size)
                for weight, size in zip(weights, sizes)
            )
            if lhs > cut_rhs(m, cut):
                raise AssertionError((m, n, "primal cut", cut, lhs - cut_rhs(m, cut)))
        objective = sum(weight * size for weight, size in zip(weights, sizes))
        f4 = F(
            2
            * (m + 2)
            * (9 * m * m * n + 100 * m * m + 55 * m * n - 100 * m + 50 * n),
            d,
        )
        if objective != f4:
            raise AssertionError((m, n, "objective", objective, f4))
        dual_objective = rho[0] * n + rho[1] * 2 * fc3(m)
        dual_objective += rho[2] * cut_rhs(m, k) + rho[3] * cut_rhs(m, k + 1)
        if dual_objective != f4:
            raise AssertionError((m, n, "dual objective", dual_objective, f4))


def main() -> None:
    symbolic_checks()
    print("r=4 symbolic primal-dual identities: PASS")
    for m in (570, 600, 1500, 15000):
        exact_full_band_check(m)
        print(f"m={m}: exact certificates for every n=16..31 PASS")
    print("r=4 subprogression verification: PASS")


if __name__ == "__main__":
    main()
