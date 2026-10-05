from __future__ import annotations

import sympy as sp


h = sp.symbols("h", integer=True, positive=True)
x, y = sp.symbols("x y", integer=True, nonnegative=True)


def c2(z):
    return z * (z - 1) / 2


def dgh_bk(a, b, k, alpha):
    d = k - 2
    beta = 2 * (d - alpha) / ((k + 1) * d - 2 * alpha)
    intercept = 2 * (a - 2) / 3 * (1 + beta * (k + 1) / d)
    intercept -= alpha * beta * (k - 1) / d
    return sp.factor(c2(a) * intercept / c2(k) + (k + 1) * (2 - beta) * b / 3)


def migrated_value(q, c, a, t, b):
    capital_c = c2(a)
    u = c - t
    r = c - ((q + 1) * t + q - 3)
    p_minus = {
        c + 1: 4 * capital_c / c2(c + 1),
        c + 2: (q - 4) * capital_c / c2(c + 2),
    }
    p_plus = {
        u: r * capital_c / c2(u),
        u + 1: (q + 1 - r) * capital_c / c2(u + 1),
    }
    b_minus = sum(p_minus.values())
    b_plus = sum(p_plus.values())
    z_minus = sum(size * weight for size, weight in p_minus.items())
    z_plus = sum(size * weight for size, weight in p_plus.items())
    lam = (b - b_minus) / (b_plus - b_minus)
    return sp.factor(z_minus + lam * (z_plus - z_minus))


def positive_after_shift(poly, symbol, shift: int) -> bool:
    shifted = sp.Poly(sp.expand(poly.subs(symbol, y + shift)), y)
    return all(coefficient > 0 for coefficient in shifted.coeffs())


q = 2 * h
c = 2 * h**2
n = 2 * h**3
t = h - 2
d0 = c - 2
alpha_c = 4 * (h - 1)

es = migrated_value(q, c, n, t, n)
best_candidate = dgh_bk(n, n, c, alpha_c)
gain = sp.factor(best_candidate - es)
gain_num, gain_den = sp.fraction(gain)
assert positive_after_shift(gain_num, h, 4)
assert positive_after_shift(gain_den, h, 4)
assert sp.limit(gain / q**2, h, sp.oo) == sp.Rational(1, 4)

# Roman's diagonal envelope is attained by R_c.
roman_left_test = sp.factor(2 * (n - 1) * (n - 2) - (c - 1) * c * (c - 2))
roman_right_test = sp.factor(c * (c - 1) * (c + 1) - 2 * (n - 1) * (n - 2))
assert positive_after_shift(roman_left_test, h, 4)
assert positive_after_shift(roman_right_test, h, 4)
roman_c = 2 * n * (n - 1) * (n - 2) / (3 * c * (c - 1)) + 2 * (c + 1) * n / 3
assert positive_after_shift(sp.fraction(sp.factor(roman_c - best_candidate))[0], h, 4)

# If B_k lies below the Roman envelope at b=n, its two crossings with
# R_k and R_{k-1} force c-h+2 <= k <= c+h.
k = sp.symbols("k", integer=True, positive=True)
left_min = sp.expand(
    (n - 1) * (2 * (k + 1) * (n - 2) - 3 * (k - 3) * (k - 1))
    - k * (k - 2) * (k - 1) * (k + 1)
)
assert sp.simplify(sp.diff(left_min, k, 2) + 4 * (3 * h**3 + 3 * k**2 - 3 * k - 2)) == 0
assert positive_after_shift(sp.factor(left_min.subs(k, 3)), h, 4)
assert positive_after_shift(sp.factor(left_min.subs(k, c - h + 1)), h, 4)

right_max = sp.expand((n - 1) * (2 * n - 10 + 3 * k) - k * (k - 1) * (k - 2))
assert positive_after_shift(sp.factor(-right_max.subs(k, c + h + 1)), h, 4)
right_derivative_at_boundary = sp.factor(-sp.diff(right_max, k).subs(k, c + h + 1))
assert positive_after_shift(right_derivative_at_boundary, h, 4)

# In the remaining interval, d=k-2 has quotient q+1 at d=c-h,
# quotient q on c-h+1 <= d <= c-1, and quotient q-1 on c <= d <= c+h-2.
def crossing_numerators(d_value, alpha_value):
    kval = d_value + 2
    left = sp.factor(
        (n - 1) * (2 * (kval + 1) * (n - 2) - 3 * alpha_value * (kval - 1))
        - kval * (kval - 2) * (kval - 1) * (kval + 1)
    )
    right = sp.factor(
        (n - 1) * (2 * n - 10 + 3 * kval - 3 * alpha_value)
        - kval * (kval - 2) * (kval - 1)
    )
    return left, right


left_singleton, _ = crossing_numerators(c - h, h - 4)
assert positive_after_shift(left_singleton, h, 4)

alpha_q_minus_1 = 2 * h**2 - 4 - (2 * h - 1) * x
_, right_q_minus_1 = crossing_numerators(c + x, alpha_q_minus_1)
right_block = sp.Poly(sp.expand((-right_q_minus_1).subs(h, x + y + 2)), x, y)
assert all(coefficient >= 0 for coefficient in right_block.coeffs())

# On the quotient-q block compare every candidate directly with d0=c-2.
d, ell = sp.symbols("d ell", positive=True)
alpha = 2 * (n - 2) - ell * d
cell = dgh_bk(n, n, d + 2, alpha).subs(ell, q)
left_difference = sp.cancel(cell.subs(d, d0 - x) - cell.subs(d, d0))
left_num, left_den = sp.fraction(left_difference)
left_certificate = sp.Poly(sp.expand(left_num.subs(h, x + y + 3)), x, y)
left_den_certificate = sp.Poly(sp.expand(left_den.subs(h, x + y + 3)), x, y)
assert all(coefficient > 0 for coefficient in left_certificate.coeffs())
assert all(coefficient > 0 for coefficient in left_den_certificate.coeffs())

right_difference = sp.factor(cell.subs(d, d0 + 1) - cell.subs(d, d0))
right_num, right_den = sp.fraction(right_difference)
assert positive_after_shift(right_num, h, 4)
assert positive_after_shift(right_den, h, 4)

print("single-cut minimizer: B_sc(2h^3,2h^3)=B_{2h^2} for every h>=4")
print("strict gain: B_sc-Es > 0 for every h>=4")
print("asymptotic gain: B_sc-Es = h^2+O(h) = (1/4+o(1))q^2")
