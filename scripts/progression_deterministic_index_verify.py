"""Verify a deterministic cut index for the five progression bands.

For r=5,...,9 and m=M_r*j, define

    k_det = ceil(Ahat_r(m) / Hhat_r(m)).

The script proves symbolically that Ahat/Hhat is the lower endpoint of the
unit-width cut-index interval obtained from the multiplier-feasible remainder
window in Theorem 3.1. It also proves that D_k and D_{k+1} have quotient r for
all j>=r+1, and checks the finite bridge r+1<=j<=14 by exact arithmetic.
"""

from __future__ import annotations

from fractions import Fraction as F
import math

import sympy as sp


j, alpha = sp.symbols("j alpha", positive=True)


DATA = {
    5: {
        "M": 15,
        "H": (1001, 6810, 12375, 6750),
        "A": (362, 2977, 9060, 13275, 6750),
    },
    6: {
        "M": 21,
        "H": (1651, 13314, 28665, 18522),
        "A": (506, 4923, 17724, 30429, 18522),
    },
    7: {
        "M": 28,
        "H": (2535, 23632, 58800, 43904),
        "A": (674, 7571, 31472, 61936, 43904),
    },
    8: {
        "M": 36,
        "H": (3689, 39024, 110160, 93312),
        "A": (866, 11029, 51984, 115344, 93312),
    },
    9: {
        "M": 45,
        "H": (5149, 60930, 192375, 182250),
        "A": (1082, 15405, 81180, 200475, 182250),
    },
}


def c2(x):
    return x * (x - 1) / 2


def c3(x):
    return x * (x - 1) * (x - 2) / 6


def polynomial(coefficients, x):
    degree = len(coefficients) - 1
    return sum(coefficient * x ** (degree - index) for index, coefficient in enumerate(coefficients))


def positive_on_tail(expression, start):
    poly = sp.Poly(sp.expand(expression), j)
    return poly.eval(start) > 0 and poly.LC() > 0 and poly.count_roots(start, sp.oo) == 0


def multiplier_window(r, m):
    k1m2 = (2 * (m - 2) - alpha) / r
    k2m2 = k1m2 + 1
    den1 = k1m2 - alpha
    den2 = k2m2 - (alpha - r)

    def theta1(x):
        return c2(x) * (x - 2 - alpha) / den1

    def theta2(x):
        return c2(x) * (x - 2 - (alpha - r)) / den2

    lower_support = [r * j + 1, r * j + 2]
    upper_support = [(r + 1) * j + 1, (r + 1) * j + 2]
    y0, ys, y1, y2 = sp.symbols("y0 ys y1 y2")
    equations = [
        sp.Eq(y0 + ys * c3(size) + y1 * theta1(size) + y2 * theta2(size), size)
        for size in lower_support
    ] + [
        sp.Eq(y0 + ys * c3(size) + (y1 + y2) * c2(size), size)
        for size in upper_support
    ]
    solution = sp.solve(equations, [y0, ys, y1, y2], dict=True)[0]
    upper_alpha = [
        sp.simplify(root)
        for root in sp.solve(sp.numer(sp.together(solution[y1])), alpha)
        if sp.degree(sp.numer(sp.together(sp.simplify(root))), j) >= 1
    ][0]
    lower_alpha = [
        sp.simplify(root)
        for root in sp.solve(sp.numer(sp.together(solution[y2])), alpha)
        if sp.degree(sp.numer(sp.together(sp.simplify(root))), j) >= 1
    ][0]
    return lower_alpha, upper_alpha


def fc2(x):
    return F(x * (x - 1), 2)


def fc3(x):
    return F(x * (x - 1) * (x - 2), 6)


def cut_coefficient(k, size, m):
    divisor = k - 2
    remainder = (2 * (m - 2)) % divisor
    if size >= k:
        return fc2(size)
    return fc2(size) * F(size - 2 - remainder, divisor - remainder)


def solve4(matrix, rhs):
    augmented = [matrix[row][:] + [rhs[row]] for row in range(4)]
    for column in range(4):
        pivot = next(row for row in range(column, 4) if augmented[row][column] != 0)
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        pivot_value = augmented[column][column]
        augmented[column] = [value / pivot_value for value in augmented[column]]
        for row in range(4):
            if row != column and augmented[row][column] != 0:
                factor = augmented[row][column]
                augmented[row] = [
                    augmented[row][index] - factor * augmented[column][index]
                    for index in range(5)
                ]
    return [augmented[row][4] for row in range(4)]


def exact_dual_feasible(r, jj, k):
    m = DATA[r]["M"] * jj
    support = [r * jj + 1, r * jj + 2, (r + 1) * jj + 1, (r + 1) * jj + 2]
    matrix = [
        [F(1), fc3(size), cut_coefficient(k, size, m), cut_coefficient(k + 1, size, m)]
        for size in support
    ]
    multipliers = solve4(matrix, [F(size) for size in support])
    if any(value < 0 for value in multipliers[1:]):
        return False
    return all(
        multipliers[0]
        + multipliers[1] * fc3(size)
        + multipliers[2] * cut_coefficient(k, size, m)
        + multipliers[3] * cut_coefficient(k + 1, size, m)
        - size
        >= 0
        for size in range(2, m + 1)
    )


def main():
    for r, data in DATA.items():
        m = data["M"] * j
        hhat = polynomial(data["H"], m)
        ahat = polynomial(data["A"], m)
        ratio = sp.cancel(ahat / hhat)

        lower_alpha, upper_alpha = multiplier_window(r, m)
        lower_k = sp.cancel((2 * (m - 2) - upper_alpha) / r + 2)
        upper_k = sp.cancel((2 * (m - 2) - lower_alpha) / r + 2)
        assert sp.simplify(ratio - lower_k) == 0
        assert sp.simplify(upper_k - ratio - 1) == 0

        quotient_upper = sp.expand((2 * m - 4) * hhat - r * ahat)
        quotient_lower = sp.expand((r + 1) * ahat - (2 * m + 2 * r - 2) * hhat)
        assert positive_on_tail(quotient_upper, r + 1)
        assert positive_on_tail(quotient_lower, r + 1)

        for jj in range(r + 1, 15):
            mm = data["M"] * jj
            hh = int(hhat.subs(j, jj))
            aa = int(ahat.subs(j, jj))
            k = math.ceil(F(aa, hh))
            assert (2 * (mm - 2)) // (k - 2) == r
            assert (2 * (mm - 2)) // (k - 1) == r
            assert exact_dual_feasible(r, jj, k)

        print(
            f"r={r}: k_det=ceil(Ahat/Hhat), unit multiplier window, "
            f"quotient tail j>={r + 1}, exact bridge through j=14: PASS"
        )

    print("deterministic progression cut index: PASS")


if __name__ == "__main__":
    main()
