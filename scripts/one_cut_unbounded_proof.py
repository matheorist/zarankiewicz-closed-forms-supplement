from __future__ import annotations

import sympy as sp


s, i, j = sp.symbols("s i j", integer=True, positive=True)
x, y, z = sp.symbols("x y z", integer=True, nonnegative=True)


def c2(value):
    return value * (value - 1) / 2


def c3(value):
    return value * (value - 1) * (value - 2) / 6


def cut_coeff(q, a, k, size, *, left):
    d = k - 2
    alpha = 2 * (a - 2) - q * d
    if left:
        return c2(size) * (size - 2 - alpha) / (d - alpha)
    return c2(size)


def cut_rhs(q, a, k):
    d = k - 2
    alpha = 2 * (a - 2) - q * d
    return c2(a) * (2 * (a - 2) - alpha) / d


def distribution(q, a, k, sizes):
    matrix = sp.Matrix(
        [
            [1, 1, 1],
            [c3(size) for size in sizes],
            [cut_coeff(q, a, k, sizes[0], left=True), c2(sizes[1]), c2(sizes[2])],
        ]
    )
    rhs = sp.Matrix([a, 2 * c3(a), cut_rhs(q, a, k)])
    weights = [sp.factor(value) for value in matrix.inv() * rhs]
    objective = sp.factor(sum(size * weight for size, weight in zip(sizes, weights)))
    return weights, objective


def local_slack(q, a, u, sizes, weights):
    k = u + j
    lhs = cut_coeff(q, a, k, sizes[0], left=True) * weights[0]
    lhs += c2(sizes[1]) * weights[1] + c2(sizes[2]) * weights[2]
    return sp.factor(cut_rhs(q, a, k) - lhs)


def migrated_data(q, c, a, t):
    capital_c = c2(a)
    u = c - t
    tile_r = t
    minus = {c + 1: 4 * capital_c / c2(c + 1), c + 2: (q - 4) * capital_c / c2(c + 2)}
    plus = {u: tile_r * capital_c / c2(u), u + 1: (q + 1 - tile_r) * capital_c / c2(u + 1)}
    b_minus, b_plus = sum(minus.values()), sum(plus.values())
    z_minus = sum(size * weight for size, weight in minus.items())
    z_plus = sum(size * weight for size, weight in plus.items())
    lam = sp.factor((a - b_minus) / (b_plus - b_minus))
    value = sp.factor(z_minus + lam * (z_plus - z_minus))
    return value, sp.factor(a - b_minus), sp.factor(b_plus - a)


def positive_coeff_certificate(expr, substitutions, variables):
    value = sp.cancel(expr.subs(substitutions))
    numerator, denominator = sp.fraction(value)
    for polynomial in (numerator, denominator):
        coefficients = sp.Poly(sp.expand(polynomial), *variables).coeffs()
        assert coefficients and all(coefficient >= 0 for coefficient in coefficients)
        assert any(coefficient > 0 for coefficient in coefficients)


t = 2 * s + 1
q = 4 * s + 4
c = 8 * s**2 + 20 * s + 7
a = 2 * (s + 1) * c
u = c - t
k0 = u + s + 2
a_t = (q + 1) * t + q - 3
b_t = (q + 1) * t + 2 * q - 2
assert sp.expand(u - a_t) == 0
assert sp.expand(c - a_t - t) == 0
assert sp.expand(b_t - c - (q + 1 - t)) == 0

sizes0 = [u + 1, c + 2, c + 3]
weights0, objective0 = distribution(q, a, k0, sizes0)
sizes_left = [u, c + 1, c + 2]
weights_left, objective_left = distribution(q, a, k0 - 1, sizes_left)
for expression in weights0 + weights_left + [objective_left - objective0]:
    positive_coeff_certificate(expression, {s: x + 5}, [x])

slack0 = local_slack(q, a, u, sizes0, weights0)
slack_left_witness = local_slack(q, a, u, sizes_left, weights_left)
positive_coeff_certificate(
    slack0 / (j - s - 2), {s: x + 5, j: x + z + 7}, [x, z]
)
positive_coeff_certificate(
    slack_left_witness / (s + 1 - j), {s: x + 5, j: z + 1}, [x, z]
)

es, lower_mass_margin, upper_mass_margin = migrated_data(q, c, a, t)
for margin in (lower_mass_margin, upper_mass_margin):
    positive_coeff_certificate(margin, {s: x + 5}, [x])
gap = sp.factor(objective0 - es)
positive_coeff_certificate(gap, {s: x + 5}, [x])
assert sp.limit(gap / t, s, sp.oo) == sp.Rational(5, 18)

# Three-row dual certificate at the central cut.
y0, ys, yk = sp.symbols("y0 ys yk")
d = k0 - 2
alpha = 2 * (a - 2) - q * d
theta = lambda size: c2(size) * (size - 2 - alpha) / (d - alpha)
dual = sp.solve(
    [
        y0 + c3(sizes0[0]) * ys + theta(sizes0[0]) * yk - sizes0[0],
        y0 + c3(sizes0[1]) * ys + c2(sizes0[1]) * yk - sizes0[1],
        y0 + c3(sizes0[2]) * ys + c2(sizes0[2]) * yk - sizes0[2],
    ],
    [y0, ys, yk],
    dict=True,
)[0]
for multiplier in dual.values():
    positive_coeff_certificate(multiplier, {s: x + 5}, [x])
slack_left = sp.factor((y0 + c3(i) * ys + theta(i) * yk - i).subs(dual))
slack_right = sp.factor((y0 + c3(i) * ys + c2(i) * yk - i).subs(dual))
w = sizes0[0]
left_quotient = sp.factor(slack_left / (w - i))
positive_coeff_certificate(left_quotient.subs(i, w - 1), {s: x + 5}, [x])
positive_coeff_certificate(-left_quotient.subs(i, w), {s: x + 5}, [x])
positive_coeff_certificate(-sp.diff(left_quotient, i), {s: x + 5, i: z + 2}, [x, z])
right_quotient = sp.factor(slack_right / ((c + 2 - i) * (c + 3 - i)))
positive_coeff_certificate(
    right_quotient.subs(i, k0 + z), {s: x + 5}, [x, z]
)

# The central point itself covers all cuts outside the transition block.
moment0 = sp.factor(sum(c2(size) * weight for size, weight in zip(sizes0, weights0)))
for expression in [moment0 - q * c2(a), (q + 1) * c2(a) - moment0]:
    positive_coeff_certificate(expression, {s: x + 5}, [x])
for special_k, flags in ((c + 2, (True, False, False)), (c + 3, (True, True, False))):
    lhs = sum(
        cut_coeff(q, a, special_k, size, left=left) * weight
        for size, weight, left in zip(sizes0, weights0, flags)
    )
    positive_coeff_certificate(cut_rhs(q, a, special_k) - lhs, {s: x + 5}, [x])

if __name__ == "__main__":
    print("unbounded diagonal one-cut certificate: PASS")
    print("parameters: t=2x+1, q=4x+4, c=8x^2+20x+7, x>=5")
    print("minimising cut: k=u+x+2")
    print("central support: {u+1,c+2,c+3}")
    print("gap/t limit:", sp.limit(gap / t, s, sp.oo))
