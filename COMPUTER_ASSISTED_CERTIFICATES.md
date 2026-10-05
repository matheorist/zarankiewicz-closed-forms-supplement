# Computer-assisted certificate audit

This audit maps every computer-assisted theorem in
`zarankiewicz_closed_forms_20260624.tex` to an executable certificate.
All theorem-level arithmetic is exact. Floating-point scans are discovery
tools only and are not part of any proof.

## Exceptional progression band: `thm:r4`

Script: `progression_r4_subprogression_verify.py`

Certified formulas and domain:

- `30 | m`, `m >= 570`;
- every integer `16 <= n <= 31`;
- support `{2m/5+1, 2m/5+2, m/2+1}`;
- moving adjacent cuts `D_k,D_(k+1)`, where
  `k=ceil((242m^3+605m^2+450m)/(549m^2+1130m+800))`.

The manuscript proof partitions the primal residuals by
`floor(2(m-2)/(j-2))`. Each residual is affine in `n`, so the two endpoints
`n=16,31` suffice. On each quotient region, cancellation leaves the
consecutive support-root factors or one of the boundary sign factors
`-(j-m)`, `-(3j-2m-2)`, and `-(2j-m-2)`, multiplied by positive factors.
The threshold `m>=570` is the endpoint condition
`E_4(31)=-m^2+555m+2050<0`.

The script independently:

- verifies the edge-count, third-moment, and objective identities symbolically;
- derives the three displayed nonzero dual multipliers from the tight-size
  equations;
- verifies the two generic dual-slack factors and the intervening slack at
  `i=k` symbolically;
- checks every primal cut and every dual reduced cost in exact rational
  arithmetic for all `n=16,...,31` at
  `m=570,600,1500,15000`.

The exact spot checks are independent diagnostics; the infinite domain follows
from the symbolic formulas and sign partition in the manuscript proof.

## Progression upper bound: `thm:r5`

Script: `progression_allbands_proof.py`

Certified domain:

- `r in {5,6,7,8,9}`
- `m = binom(r+1,2) j`
- every integer `j >= 5`
- every `n` in the cubic band `I_r`

Obligations:

- derives the two-cut multiplier window and proves its width is exactly `r`;
- verifies remainder legality and nonnegativity of all dual multipliers;
- factors the generic dual slacks and the exceptional slack at `i=k`;
- verifies that the dual objective is the displayed rational function `F_r`;
- proves the symbolic tail for `j >= 15`;
- checks exact rational dual feasibility for every `5 <= j <= 120`.

The symbolic tail and exact bridge overlap, hence cover all `j >= 5`.
Successful output ends without `certificate verification failed` and reports
all five bands feasible.

### Deterministic cut index

Script: `progression_deterministic_index_verify.py`

The existence proof in `thm:r5` chooses a remainder in a feasible interval.
The following equivalent deterministic choice, recovered from the earlier
certificate derivation, is useful for reproduction:

`k_det = ceil(Ahat_r(m)/Hhat_r(m))`.

The hats distinguish these index polynomials from the remainder-window
functions `A_r(j),B_r(j)` in the manuscript.

| `r` | `Hhat_r(m)` | `Ahat_r(m)` |
|---:|---|---|
| 5 | `1001m^3+6810m^2+12375m+6750` | `362m^4+2977m^3+9060m^2+13275m+6750` |
| 6 | `1651m^3+13314m^2+28665m+18522` | `506m^4+4923m^3+17724m^2+30429m+18522` |
| 7 | `2535m^3+23632m^2+58800m+43904` | `674m^4+7571m^3+31472m^2+61936m+43904` |
| 8 | `3689m^3+39024m^2+110160m+93312` | `866m^4+11029m^3+51984m^2+115344m+93312` |
| 9 | `5149m^3+60930m^2+192375m+182250` | `1082m^4+15405m^3+81180m^2+200475m+182250` |

For `m=M_r j`, the script proves symbolically that `Ahat_r/Hhat_r` is exactly
the lower endpoint of the cut-index interval obtained from
`A_r(j) <= alpha <= B_r(j)`, and that the upper endpoint is one larger.
It also proves both cuts have quotient `r` for `j >= r+1` and checks exact dual
feasibility for `r+1 <= j <= 14`. Combined with the all-window symbolic slack
proof for `j >= 15`, this certifies the deterministic index for every
`j >= r+1`. The general theorem still covers `5 <= j < r+1` by its
remainder-window argument and exact bridge.

## Progression identity: `thm:r5id`

Script: `progression_allbands_lower.py`

Certified domain:

- the same five bands;
- every `n in I_r`;
- `j >= j_r`, where `(j_5,...,j_9)=(22,41,41,213,67)`.

Obligations:

- solves the four-row primal system on
  `{rj+1,rj+2,(r+1)j+1,(r+1)j+2}`;
- verifies its objective is `F_r`;
- derives the exact nonnegativity thresholds `j_r`;
- factors the remaining cut slacks by quotient range;
- performs exact feasibility checks beginning at every threshold.

Successful output reports `obj=F_r:True`, all symbolic cut checks true, and
exact feasibility true for every band.

## Fixed-support companion residues: `thm:companion`

Script: `companion_certificate_verify.py`

Certified domain:

- `q >= 5`, `0 <= rho <= 4`;
- `2a = qc+rho`;
- `3q+rho-1 <= c <= 4q+rho`;
- every integer `b` between the displayed endpoint masses.

Obligations:

- verifies the four primal weights, endpoint masses, and objective;
- verifies signs of the four active-row multipliers;
- factors both generic dual-slack ranges;
- checks the left-boundary case where the `D_{c-1}` multiplier vanishes;
- confirms exact primal-dual equality.

Successful output reports that the companion certificate checks pass.
The same script also verifies the diagnostic boundary example
`q=5, a=30, b=28,...,43` by exact rational LP, confirming that the fixed-basis
value differs from `Es` at every listed point.  This finite check is not used
to extend the theorem beyond its stated range.

## All-residue migrated support: `thm:allres-migration` and `cor:coverage-compatibility`

Script: `all_residue_migration_proof.py`

Certified symbolic domain:

- `q >= 5`, `0 <= rho < q`, `t >= 1`;
- `2a=qc+rho` is even;
- `(q+1)t+q+rho-3 < c <= (q+1)t+2q+rho-2`.

Obligations:

- verifies both endpoint distributions, their second and third moments,
  endpoint masses, and endpoint objectives;
- proves primal feasibility for every cut;
- verifies the telescoping local block `K={u+1,...,v}`;
- proves positivity of the dual multipliers and all interior and exterior
  slacks;
- verifies exact strong duality for an interior convex combination;
- checks quotient gluing with the next `(q+1)` band.

The symbolic mechanism is supplemented by exact checks on
`5 <= q <= 12`, all residues, `1 <= t <= 5`, and representative values of the
tile coordinate `r`. The expected totals are 740 migrated-support checks and
740 quotient-gluing identities.

Script: `all_residue_exact_spotcheck.py`

Independent check: compares the candidate-support LP with the full all-cut LP
at five exact rational sample points. All five values must agree.

## Diagonal single-cut comparisons

Labels: `thm:scgain`, `thm:trueonecut`, and `thm:unboundedonecut`.

Script: `single_cut_gain_proof.py`

- certifies the minimizing published DGH cut on `N=2h^3`, `h>=4`;
- verifies the exact positive gap and its `h^2+O(h)` expansion.

Script: `one_cut_global_proof.py`

- certifies the exact best one-cut LP on `N=2h^3`;
- verifies the global primal witnesses and central dual certificate;
- verifies limiting all-cut gaps `5/12` for even `h` and `7/12` for odd `h`.

Script: `one_cut_unbounded_proof.py`

- certifies the second diagonal family stated in Theorem 5.4;
- verifies the mass interval, global one-cut witnesses, and central dual;
- verifies `gap/t -> 5/18`, hence an unbounded `Theta(N^(1/3))` gap.

## Explicit two-edge near-endpoint gaps

Script: `near_endpoint_gap_verify.py`

This script constructs the finite integer systems in equation (49) from the
endpoint dual slacks.  All cut coefficients are cleared to integers before the
finite feasibility check.  The SCIP backend is run with a tightened feasibility
tolerance, and any feasible profile is independently rechecked with Python
integers.

The script checks every `d=0,...,8` for both parameter pairs
`(q,c)=(297,878)` and `(q,c)=(576,1598)`.  It also reconstructs the active
degree windows `{848,...,909}` and `{1558,...,1639}` from the rational endpoint
dual.  Each of the 18 systems must be reported as infeasible.

## Aggregate numerical reproduction

Script: `verification.ps1`

This reproduces the displayed comparison table, fixed-band dual-objective
identities, migrated diagonal ratios, and Brown construction block counts.
Its floating-point values are illustrations, not proof inputs.

## One-command order

Run the commands in `REPRODUCIBILITY.md` from this directory. A failed
assertion, a nonzero process exit code, or an output total different from the
one stated above is a failed certificate and must not be ignored.
