from __future__ import annotations

from fractions import Fraction as F

import sympy as sp


def c2(x: int) -> F:
    return F(x * (x - 1), 2)


def c3(x: int) -> F:
    return F(x * (x - 1) * (x - 2), 6)


def cut_coeff(a: int, k: int, i: int) -> F:
    d = k - 2
    alpha = (2 * (a - 2)) % d
    if i >= k:
        return c2(i)
    return c2(i) * F(i - 2 - alpha, d - alpha)


def cut_rhs(a: int, k: int) -> F:
    d = k - 2
    alpha = (2 * (a - 2)) % d
    return c2(a) * F(2 * (a - 2) - alpha, d)


def symbolic_checks() -> None:
    k, y = sp.symbols("k y")
    f = lambda z: 2 * (y - z) / (z * (z - 1))
    second_difference = sp.factor(f(k + 1) - 2 * f(k) + f(k - 1))
    expected = 4 * (3 * y - k - 1) / ((k - 2) * (k - 1) * k * (k + 1))
    assert sp.simplify(second_difference - expected) == 0

    # The threshold used to prove 2*y_0 > v+1 is checked in the two
    # residue branches.  This sharpened split is what lets the certificate
    # start at t=1.
    q, t, Q, T = sp.symbols("q t Q T", nonnegative=True)
    margin_low = (
        (q * (t + 1) - 2) * (q * (t + 1) - 3)
        - (q + 1) * (t + 1) * (t + 2)
    )
    shifted_low = sp.Poly(sp.expand(margin_low.subs({q: Q + 5, t: T + 1})), Q, T)
    assert all(coefficient > 0 for coefficient in shifted_low.coeffs())

    margin_high = (
        (q * (t + 1) + 3) * (q * (t + 1) + 2)
        - (q + 1) * (t + 2) * (t + 3)
    )
    shifted_high = sp.Poly(sp.expand(margin_high.subs({q: Q + 6, t: T + 1})), Q, T)
    assert all(coefficient > 0 for coefficient in shifted_high.coeffs())

    i, v = sp.symbols("i v", positive=True)
    threshold_identity = (v + 1 - 2 * i) / (i * (i - 1)) + 1 / v
    expected_threshold = (v - i) * (v - i + 1) / (i * v * (i - 1))
    assert sp.simplify(threshold_identity - expected_threshold) == 0


def exact_certificate(q: int, rho: int, t: int, r: int) -> bool:
    c = (q + 1) * t + q + rho - 3 + r
    x = q * c + rho
    if x % 2:
        return False

    a = x // 2
    h = 1 if rho <= 4 else 2
    u = c - t
    v = c + h
    beta = 4 - rho if h == 1 else q + 4 - rho
    C = c2(a)

    right = {
        v: F(beta) * C / c2(v),
        v + 1: F(q - beta) * C / c2(v + 1),
    }
    left = {
        u: F(r) * C / c2(u),
        u + 1: F(q + 1 - r) * C / c2(u + 1),
    }

    for point, moment in ((right, q), (left, q + 1)):
        assert all(weight >= 0 for weight in point.values())
        assert sum(c2(i) * weight for i, weight in point.items()) == moment * C
        assert sum(c3(i) * weight for i, weight in point.items()) == 2 * c3(a)
        for cut in range(3, a + 1):
            lhs = sum(cut_coeff(a, cut, i) * weight for i, weight in point.items())
            assert lhs <= cut_rhs(a, cut)

    A = (
        F(r, u * (u - 1))
        + F(q + 1 - r, u * (u + 1))
        - F(beta, v * (v - 1))
        - F(q - beta, v * (v + 1))
    )
    B = (
        F(r, u - 1)
        + F(q + 1 - r, u)
        - F(beta, v - 1)
        - F(q - beta, v)
    )
    assert A > 0
    y0 = B / A
    assert 2 * y0 > v + 1

    def f(i: int) -> F:
        return 2 * (y0 - i) / (i * (i - 1))

    ys = -3 * (f(v + 1) - f(v))
    assert ys > 0

    multipliers: dict[int, F] = {}
    for cut in range(u + 1, v + 1):
        assert (x - 4) // (cut - 2) == q
        alpha = q * (c - cut + 2) + rho - 4
        delta = cut - 2 - alpha
        multiplier = F(4 * delta) * (3 * y0 - cut - 1)
        multiplier /= (cut - 2) * (cut - 1) * cut * (cut + 1)
        assert multiplier > 0
        multipliers[cut] = multiplier

    Y = sum(multipliers.values(), F(0))
    assert Y == -f(v + 1) + (v - 1) * (f(v + 1) - f(v))

    for i in range(2, a + 1):
        normalized_slack = f(i) + F(i - 2, 3) * ys + Y
        for cut, multiplier in multipliers.items():
            alpha = (x - 4) % (cut - 2)
            delta = cut - 2 - alpha
            normalized_slack -= multiplier * F(max(cut - i, 0), delta)
        assert normalized_slack >= 0
        if u <= i <= v + 1:
            assert normalized_slack == 0

    lam = F(2, 5)
    primal = {
        i: (1 - lam) * right.get(i, F(0)) + lam * left.get(i, F(0))
        for i in set(right) | set(left)
    }
    mass = sum(primal.values())
    primal_value = sum(i * weight for i, weight in primal.items())
    dual_value = mass * y0 + 2 * c3(a) * ys
    dual_value += sum(cut_rhs(a, cut) * yk for cut, yk in multipliers.items())
    assert primal_value == dual_value
    return True


def exact_quotient_gluing(q: int, rho: int, t: int, r: int) -> bool:
    c = (q + 1) * t + q + rho - 3 + r
    x = q * c + rho
    if x % 2:
        return False

    a = x // 2
    C = c2(a)
    u = c - t
    upper = {
        u: F(r) * C / c2(u),
        u + 1: F(q + 1 - r) * C / c2(u + 1),
    }

    qp = q + 1
    cp, rhop = divmod(x, qp)
    hp = 1 if rhop <= 4 else 2
    betap = 4 - rhop if hp == 1 else qp + 4 - rhop
    vp = cp + hp
    next_lower = {
        vp: F(betap) * C / c2(vp),
        vp + 1: F(qp - betap) * C / c2(vp + 1),
    }

    upper = {i: weight for i, weight in upper.items() if weight}
    next_lower = {i: weight for i, weight in next_lower.items() if weight}
    assert upper == next_lower
    return True


def main() -> None:
    symbolic_checks()
    checked = 0
    glued = 0
    for q in range(5, 13):
        for rho in range(q):
            for t in range(1, 6):
                old_max = q - rho + 1
                r_values = {1, (q + 2) // 2, old_max, q + 1}
                if old_max < q + 1:
                    r_values.add(old_max + 1)
                for r in sorted(r_values):
                    checked += exact_certificate(q, rho, t, r)
                    glued += exact_quotient_gluing(q, rho, t, r)

    print("all-residue migrated-support certificate: PASS")
    print("symbolic mechanism: telescoping local cut block")
    print(f"exact finite certificates checked: {checked}")
    print(f"exact quotient gluings checked: {glued}")
    print("right support shift: h=1 for rho<=4; h=2 for rho>=5")


if __name__ == "__main__":
    main()
