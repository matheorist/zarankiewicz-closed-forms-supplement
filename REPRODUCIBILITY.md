# Reproducibility guide

This file records the computational components for
`zarankiewicz_closed_forms_20260624.tex`.

For a theorem-by-theorem statement of the certified domains, symbolic
obligations, finite bridges, and expected terminal output, see
`COMPUTER_ASSISTED_CERTIFICATES.md`.

## Environment

- Python 3.12.10
- SymPy 1.14.0
- OR-Tools 9.15.6755 (SCIP backend for the finite integer checks)
- PowerShell 7 or Windows PowerShell 5.1 for `verification.ps1`

The theorem-level identities use exact SymPy algebra or Python
`fractions.Fraction`. Floating-point HiGHS solves are used only for pattern
discovery. The two near-endpoint examples use an integer-coefficient
feasibility model with the SCIP backend and recheck any reported feasible
profile with Python integers.

## Core checks

Run these commands from this directory:

```powershell
python scripts/progression_r4_subprogression_verify.py
python scripts/progression_allbands_proof.py
python scripts/progression_deterministic_index_verify.py
python scripts/progression_allbands_lower.py
python scripts/companion_certificate_verify.py
python scripts/all_residue_migration_proof.py
python scripts/all_residue_exact_spotcheck.py
python scripts/single_cut_gain_proof.py
python scripts/one_cut_global_proof.py
python scripts/one_cut_unbounded_proof.py
python scripts/near_endpoint_gap_verify.py
powershell.exe -NoProfile -ExecutionPolicy Bypass -File verification.ps1
```

Expected high-level results:

- `progression_r4_subprogression_verify.py`: verifies the symbolic primal and
  dual identities for the `thm:r4` subprogression certificate and checks
  every `n=16,...,31`, every primal cut, and every dual reduced cost exactly at
  `m=570,600,1500,15000`.
- `progression_allbands_proof.py`: for `r=5,...,9`, the symbolic checks of
  the multiplier window, remainder legality, multiplier signs, generic and
  exceptional dual slacks, and dual objective pass for all `j >= 15`;
  exact dual feasibility passes for `5 <= j <= 120`. Together these checks
  certify the upper bound for every `j >= 5`.
- `progression_deterministic_index_verify.py`: verifies the explicit choice
  `k_det=ceil(Ahat_r(m)/Hhat_r(m))`. It proves that `Ahat/Hhat` is the lower
  endpoint of the unit-width multiplier-feasible cut-index interval, proves
  both selected cuts have quotient `r` for `j >= r+1`, and checks the exact
  bridge through `j=14`.
- `progression_allbands_lower.py`: the five progression identities and
  thresholds `(22,41,41,213,67)` are reproduced.
- `companion_certificate_verify.py`: the fixed-support residue certificates
  pass on the implemented symbolic and exact checks; the sub-threshold example
  `q=5, a=30, b=28,...,43` is also checked pointwise and satisfies `Es != F_q`
  at every point.
- `all_residue_migration_proof.py`: the all-residue endpoint distributions,
  telescoping local-cut-block identity, multiplier signs, dual slacks, and
  strong-duality equality pass throughout the complete migrated range
  `c >= 2q + rho - 1`; 740 exact finite checks and 740 exact quotient-gluing
  identities are also checked.
- `all_residue_exact_spotcheck.py`: all five exact candidate/full-LP
  comparisons agree.
- `single_cut_gain_proof.py`: on the infinite family `N=2h^3`, `h>=4`, the
  best DGH Theorem 1.3 closed form is attained at `k=2h^2`; its exact gap over
  the all-`v=2` optimum is positive and equals `h^2+O(h)`.
- `one_cut_global_proof.py`: the exact best one-cut LP is attained at
  `k=2h^2-floor((h-3)/2)`; the global lower witnesses and the central dual
  certificate pass, and the all-cut gaps tend to `5/12` for even `h` and
  `7/12` for odd `h`.
- `one_cut_unbounded_proof.py`: on the diagonal family
  `t=2x+1`, `q=4x+4`, `c=8x^2+20x+7`, `N=2(x+1)c`, the best one-cut LP is
  attained at `k=c-t+s+2`; its exact gap over the all-cut optimum is positive
  and satisfies `gap/t -> 5/18`, hence grows as `Theta(N^(1/3))`.
- `near_endpoint_gap_verify.py`: constructs the two finite integer systems in
  equation (49), with exact integer coefficients, and proves infeasibility for
  every `d=0,...,8` at `(q,c)=(297,878)` and `(576,1598)`. The script requires
  OR-Tools with its SCIP backend.
- `verification.ps1`: the displayed comparison table, dual-objective
  identities, migrated diagonal ratios, and Brown block counts are reproduced.

## Pattern-discovery scans

```powershell
python all_residue_shifted_candidate_scan.py --q-min 5 --q-max 12 --max-offset 20 --out all_residue_shifted_candidate_scan_20260711.md
```

Expected result: `455/455` tested midpoints match the candidate value.
This table records pattern discovery. The proof of Theorem 4.8 is supplied by
`all_residue_migration_proof.py` and does not depend on floating-point output.

## LaTeX

The manuscript is compiled from the submission folder with two passes of
`pdflatex`; the source file is not duplicated in this repository.
