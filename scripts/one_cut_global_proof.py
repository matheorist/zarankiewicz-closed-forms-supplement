from __future__ import annotations

import sympy as sp


r, i, j = sp.symbols("r i j", integer=True, positive=True)
x, y, z = sp.symbols("x y z", integer=True, nonnegative=True)


def c2(value):
    return value * (value - 1) / 2


def c3(value):
    return value * (value - 1) * (value - 2) / 6


def cut_coeff(q, n, k, size, *, left):
    d = k - 2
    alpha = 2 * (n - 2) - q * d
    if left:
        return c2(size) * (size - 2 - alpha) / (d - alpha)
    return c2(size)


def cut_rhs(q, n, k):
    d = k - 2
    alpha = 2 * (n - 2) - q * d
    return c2(n) * (2 * (n - 2) - alpha) / d


def distribution(q, n, k, sizes, left_flags):
    matrix = sp.Matrix(
        [
            [1, 1, 1],
            [c3(size) for size in sizes],
            [cut_coeff(q, n, k, size, left=left) for size, left in zip(sizes, left_flags)],
        ]
    )
    rhs = sp.Matrix([n, 2 * c3(n), cut_rhs(q, n, k)])
    weights = [sp.factor(value) for value in matrix.inv() * rhs]
    objective = sp.factor(sum(size * weight for size, weight in zip(sizes, weights)))
    moment = sp.factor(sum(c2(size) * weight for size, weight in zip(sizes, weights)))
    return weights, objective, moment


def local_slack(q, n, u, sizes, weights, left_flags):
    k = u + j
    lhs = sum(
        cut_coeff(q, n, k, size, left=left) * weight
        for size, weight, left in zip(sizes, weights, left_flags)
    )
    return sp.factor(cut_rhs(q, n, k) - lhs)


def migrated_value(q, c, n, t):
    capital_c = c2(n)
    u = c - t
    tile_r = c - ((q + 1) * t + q - 3)
    minus = {c + 1: 4 * capital_c / c2(c + 1), c + 2: (q - 4) * capital_c / c2(c + 2)}
    plus = {u: tile_r * capital_c / c2(u), u + 1: (q + 1 - tile_r) * capital_c / c2(u + 1)}
    b_minus = sum(minus.values())
    b_plus = sum(plus.values())
    z_minus = sum(size * weight for size, weight in minus.items())
    z_plus = sum(size * weight for size, weight in plus.items())
    lam = sp.factor((n - b_minus) / (b_plus - b_minus))
    return sp.factor(z_minus + lam * (z_plus - z_minus))


def positive_coeff_certificate(expr, substitutions, variables):
    value = sp.cancel(expr.subs(substitutions))
    numerator, denominator = sp.fraction(value)
    for polynomial in (numerator, denominator):
        coefficients = sp.Poly(sp.expand(polynomial), *variables).coeffs()
        assert coefficients and all(coefficient >= 0 for coefficient in coefficients)
        assert any(coefficient > 0 for coefficient in coefficients)


def central_certificate(q, c, n, k, w, lower):
    sizes = [w, c + 1, c + 2]
    weights, objective, _ = distribution(q, n, k, sizes, [True, False, False])
    for weight in weights:
        positive_coeff_certificate(weight, {r: x + lower}, [x])

    y0, ys, yk = sp.symbols("y0 ys yk")
    d = k - 2
    alpha = 2 * (n - 2) - q * d
    theta = lambda size: c2(size) * (size - 2 - alpha) / (d - alpha)
    dual = sp.solve(
        [
            y0 + c3(w) * ys + theta(w) * yk - w,
            y0 + c3(c + 1) * ys + c2(c + 1) * yk - (c + 1),
            y0 + c3(c + 2) * ys + c2(c + 2) * yk - (c + 2),
        ],
        [y0, ys, yk],
        dict=True,
    )[0]
    for multiplier in dual.values():
        positive_coeff_certificate(multiplier, {r: x + lower}, [x])

    slack_left = sp.factor((y0 + c3(i) * ys + theta(i) * yk - i).subs(dual))
    slack_right = sp.factor((y0 + c3(i) * ys + c2(i) * yk - i).subs(dual))
    quotient_left = sp.factor(slack_left / (w - i))
    positive_coeff_certificate(quotient_left.subs(i, w - 1), {r: x + lower}, [x])
    positive_coeff_certificate(-quotient_left.subs(i, w), {r: x + lower}, [x])
    positive_coeff_certificate(-sp.diff(quotient_left, i), {r: x + lower, i: z + 2}, [x, z])
    quotient_right = sp.factor(slack_right / ((c + 1 - i) * (c + 2 - i)))
    positive_coeff_certificate(quotient_right, {r: x + lower, i: k.subs(r, x + lower) + z}, [x, z])
    return objective


def check_even():
    h = 2 * r
    q, c, n, u = 2 * h, 2 * h**2, 2 * h**3, 2 * h**2 - h + 2
    k0, w0 = c - r + 2, u + 1
    objective0 = central_certificate(q, c, n, k0, w0, 2)

    sizes_minus = [u - 1, c + 1, c + 2]
    weights_minus, objective_minus, _ = distribution(
        q, n, k0 - 1, sizes_minus, [True, False, False]
    )
    for weight in weights_minus:
        positive_coeff_certificate(weight, {r: x + 2}, [x])
    positive_coeff_certificate(objective_minus - objective0, {r: x + 2}, [x])

    weights0, _, _ = distribution(q, n, k0, [w0, c + 1, c + 2], [True, False, False])
    slack_minus = local_slack(q, n, u, sizes_minus, weights_minus, [True, False, False])
    slack0 = local_slack(q, n, u, [w0, c + 1, c + 2], weights0, [True, False, False])
    positive_coeff_certificate(
        slack_minus / (r - 1 - j), {j: z + 1, r: z + y + 2}, [y, z]
    )
    positive_coeff_certificate(slack0 / (j - r), {r: x + 2, j: x + z + 2}, [x, z])

    es = migrated_value(q, c, n, h - 2)
    gap = sp.factor(objective0 - es)
    positive_coeff_certificate(gap, {r: x + 2}, [x])
    assert sp.limit(gap, r, sp.oo) == sp.Rational(5, 12)
    return gap


def check_odd():
    h = 2 * r + 1
    q, c, n, u = 2 * h, 2 * h**2, 2 * h**3, 2 * h**2 - h + 2
    k0, w0 = c - r + 1, u
    objective0 = central_certificate(q, c, n, k0, w0, 2)
    weights0, _, _ = distribution(q, n, k0, [w0, c + 1, c + 2], [True, False, False])
    slack0 = local_slack(q, n, u, [w0, c + 1, c + 2], weights0, [True, False, False])
    positive_coeff_certificate(slack0 / (r - j), {j: z + 1, r: z + y + 2}, [y, z])

    # For r=2,...,5 the right witness uses {u+1,u+2,c+2}.
    old_sizes = [u + 1, u + 2, c + 2]
    old_weights, old_objective, _ = distribution(q, n, k0 + 1, old_sizes, [True, True, False])
    old_slack = local_slack(q, n, u, old_sizes, old_weights, [True, True, False])
    for rv in range(2, 6):
        assert all(weight.subs(r, rv) > 0 for weight in old_weights)
        assert (old_objective - objective0).subs(r, rv) > 0
        for jv in range(rv + 1, 2 * rv + 1):
            assert old_slack.subs({r: rv, j: jv}) >= 0

    # From r=6 onward the right witness migrates to {u+2,c+1,c+2}.
    new_sizes = [u + 2, c + 1, c + 2]
    new_weights, new_objective, _ = distribution(q, n, k0 + 1, new_sizes, [True, False, False])
    new_slack = local_slack(q, n, u, new_sizes, new_weights, [True, False, False])
    for weight in new_weights:
        positive_coeff_certificate(weight, {r: x + 6}, [x])
    positive_coeff_certificate(new_objective - objective0, {r: x + 6}, [x])
    positive_coeff_certificate(
        new_slack / (j - r - 1), {r: x + 6, j: x + z + 7}, [x, z]
    )

    es = migrated_value(q, c, n, h - 2)
    gap = sp.factor(objective0 - es)
    positive_coeff_certificate(gap, {r: x + 2}, [x])
    assert sp.limit(gap, r, sp.oo) == sp.Rational(7, 12)
    return gap


def check_roman():
    for h, lower in ((2 * r, 2), (2 * r + 1, 2)):
        q, c, n = 2 * h, 2 * h**2, 2 * h**3
        x0, x1 = sp.symbols("x0 x1")
        solution = sp.solve(
            [x0 + x1 - n, c3(c) * x0 + c3(c + 1) * x1 - 2 * c3(n)],
            [x0, x1],
            dict=True,
        )[0]
        weights = [sp.factor(solution[x0]), sp.factor(solution[x1])]
        moment = sp.factor(c2(c) * weights[0] + c2(c + 1) * weights[1])
        for expression in weights + [moment - q * c2(n), (q + 1) * c2(n) - moment]:
            positive_coeff_certificate(expression, {r: x + lower}, [x])
        transition = sp.factor(
            cut_rhs(q, n, c + 1)
            - cut_coeff(q, n, c + 1, c, left=True) * weights[0]
            - c2(c + 1) * weights[1]
        )
        positive_coeff_certificate(-transition, {r: x + lower}, [x])


if __name__ == "__main__":
    even_gap = check_even()
    odd_gap = check_odd()
    check_roman()
    print("one-cut global certificate: PASS")
    print("even gap limit:", sp.limit(even_gap, r, sp.oo))
    print("odd gap limit:", sp.limit(odd_gap, r, sp.oo))
    print("minimising cut: c-floor((h-3)/2)")
