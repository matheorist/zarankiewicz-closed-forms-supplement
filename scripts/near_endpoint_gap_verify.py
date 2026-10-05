"""Verify the two explicit near-endpoint Delta=2 obstructions.

The finite systems in equation (one-unit-system) are built with exact integer
coefficients.  SCIP is used only for the finite integer-feasibility decision;
the feasibility tolerance is tightened, and any reported feasible point is
checked again with Python integers before it is accepted.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb, gcd

from ortools.linear_solver import pywraplp


CASES = ((297, 878), (576, 1598))
MAX_DEFECT = 8


def cut_coefficient(a: int, cut: int, degree: int) -> Fraction:
    """Return theta_cut(degree) from the manuscript's cut definition."""
    alpha = (2 * (a - 2)) % (cut - 2)
    if degree < cut:
        return Fraction(
            comb(degree, 2) * (degree - 2 - alpha),
            cut - 2 - alpha,
        )
    return Fraction(comb(degree, 2))


def solve_linear_system(matrix: list[list[Fraction]]) -> list[Fraction]:
    """Solve a small nonsingular rational system by exact elimination."""
    work = [row[:] for row in matrix]
    size = len(work)
    for column in range(size):
        pivot = next(row for row in range(column, size) if work[row][column])
        work[column], work[pivot] = work[pivot], work[column]
        scale = work[column][column]
        work[column] = [value / scale for value in work[column]]
        for row in range(size):
            if row == column or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [
                work[row][index] - scale * work[column][index]
                for index in range(size + 1)
            ]
    return [work[index][size] for index in range(size)]


def endpoint_data(q: int, c: int) -> tuple[int, int, int, list[int]]:
    """Return (a, b0, k, S1) using exact endpoint-dual arithmetic."""
    a = (q * c + 4) // 2
    k = c + 2
    b0 = q * comb(a, 2) // comb(k, 2)
    support = (c - 1, c, c + 1, c + 2)
    dual = solve_linear_system(
        [
            [
                Fraction(1),
                Fraction(comb(i, 3)),
                cut_coefficient(a, c, i),
                cut_coefficient(a, c + 1, i),
                Fraction(i),
            ]
            for i in support
        ]
    )
    slacks = []
    for degree in range(2, a + 1):
        slack = (
            dual[0]
            + Fraction(comb(degree, 3)) * dual[1]
            + cut_coefficient(a, c, degree) * dual[2]
            + cut_coefficient(a, c + 1, degree) * dual[3]
            - degree
        )
        if slack <= 1:
            slacks.append(degree)
    return a, b0, k, slacks


def integer_coefficients(values: list[Fraction]) -> list[int]:
    denominator = 1
    for value in values:
        denominator = denominator * value.denominator // gcd(
            denominator, value.denominator
        )
    return [int(value * denominator) for value in values]


def build_model(q: int, c: int, defect: int) -> tuple[pywraplp.Solver, list[int], list, list]:
    a, b0, k, active_degrees = endpoint_data(q, c)
    triple_budget = comb(k - 1, 2)
    point_budget = 2 * comb(a - 1, 2)
    psi = {degree: comb(degree - 1, 2) for degree in active_degrees}
    penalty = {
        degree: comb(degree, 3)
        - comb(k, 3)
        - comb(k - 1, 2) * (degree - k)
        for degree in active_degrees
    }

    # The edge equation forces the total positive deviation to balance the
    # degrees below k-1.  The resulting bounds are valid consequences of the
    # reduced triple inequality and substantially shorten the finite search.
    positive_unit_cost = min(
        Fraction(penalty[degree], degree - k)
        for degree in active_degrees
        if degree > k
    )
    positive_unit_cost = int(positive_unit_cost)
    upper_bounds = {}
    for degree in active_degrees:
        deviation = degree - k
        cost = penalty[degree]
        if deviation > 0:
            upper_bounds[degree] = min(b0, triple_budget // cost)
        elif deviation < 0:
            upper_bounds[degree] = min(
                b0,
                (triple_budget + positive_unit_cost)
                // (cost + positive_unit_cost * abs(deviation)),
            )
        else:
            upper_bounds[degree] = b0

    solver = pywraplp.Solver.CreateSolver("SCIP")
    if solver is None:
        raise RuntimeError("SCIP is required for the finite integer check")
    solver.SetTimeLimit(60_000)
    solver.SetSolverSpecificParametersAsString("numerics/feastol = 1e-12")

    x = [
        solver.IntVar(0, upper_bounds[degree], f"x_{degree}")
        for degree in active_degrees
    ]
    t = [
        solver.IntVar(
            0,
            min(b0, point_budget // psi[degree]),
            f"t_{degree}",
        )
        for degree in active_degrees
    ]

    solver.Add(solver.Sum(x) == b0)
    solver.Add(
        solver.Sum(degree * value for degree, value in zip(active_degrees, x))
        == k * b0 - 1
    )
    solver.Add(
        solver.Sum(penalty[degree] * value for degree, value in zip(active_degrees, x))
        <= triple_budget
    )
    for cut in (c, c + 1):
        differences = [
            cut_coefficient(a, cut, degree) - cut_coefficient(a, cut, k)
            for degree in active_degrees
        ]
        coefficients = integer_coefficients(differences)
        solver.Add(solver.Sum(coefficient * value for coefficient, value in zip(coefficients, x)) <= 0)

    solver.Add(
        solver.Sum(psi[degree] * value for degree, value in zip(active_degrees, t))
        == point_budget - defect
    )
    for x_value, t_value in zip(x, t):
        solver.Add(t_value <= x_value)

    return solver, active_degrees, x, t


def check_case(q: int, c: int) -> None:
    a, b0, k, active_degrees = endpoint_data(q, c)
    expected_first, expected_last = active_degrees[0], active_degrees[-1]
    expected = {
        (297, 878): (848, 909),
        (576, 1598): (1558, 1639),
    }[(q, c)]
    assert (expected_first, expected_last) == expected

    for defect in range(MAX_DEFECT + 1):
        solver, degrees, x, t = build_model(q, c, defect)
        status = solver.Solve()
        if status in (pywraplp.Solver.FEASIBLE, pywraplp.Solver.OPTIMAL):
            x_values = [round(value.solution_value()) for value in x]
            t_values = [round(value.solution_value()) for value in t]
            point_value = sum(
                comb(degree - 1, 2) * value
                for degree, value in zip(degrees, t_values)
            )
            expected_point = 2 * comb(a - 1, 2) - defect
            assert point_value == expected_point, (
                q,
                c,
                defect,
                point_value,
                expected_point,
            )
            raise AssertionError(
                f"unexpected feasible profile at (q,c,d)=({q},{c},{defect})"
            )
        if status != pywraplp.Solver.INFEASIBLE:
            raise RuntimeError(
                f"SCIP did not prove infeasibility at (q,c,d)=({q},{c},{defect}); "
                f"status={status}"
            )
        print(f"(q,c)=({q},{c}), d={defect}: infeasible")

    print(
        f"PASS: (a,b0,k)=({a},{b0},{k}); "
        f"all d=0,...,{MAX_DEFECT} are infeasible"
    )


def main() -> None:
    for q, c in CASES:
        check_case(q, c)


if __name__ == "__main__":
    main()
