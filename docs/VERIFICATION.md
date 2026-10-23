# Release verification

Version 1.0.2 was checked in a clean Windows virtual environment with Python
3.12, NumPy 2.5.4 and SciPy 1.18.1. These are observed verification versions;
the broader dependency ranges are convenience specifications, not a claim that
every combination has been tested.

- All original/adapted source files have SHA-256 provenance entries.
- All Python source files parse successfully.
- The release content was scanned for local drive paths, credential token
  patterns and private-key blocks; none was found.
- Archive consistency checks passed: 44 candidates, six directly reached
  buildings plus one additional neighbourhood arrival; 28 nonzero source–building
  pairs; 117 planned FE positions, 96 qualified and 21 NC; 81 capacity constraints
  (43 finite brackets, 36 right-censored and two left-censored).
- All 531 archived fragility bound rows reproduce to maximum absolute difference
  1.11 × 10^-16. All eight signed sensitivity coefficients reproduce exactly.
- All 100,000 annual-example parameter draws match the archived inputs, and all
  EAL outputs agree within floating-point tolerance. Analytic EAL is
  4974.7423623814275 CNY/year. The original compound-Poisson simulation reproduces
  4989.45 CNY/year, standard error 72.58272545687214 CNY/year. All seven Sobol
  first/total indices match the archived table.
- The three mitigation BCR values reproduce as 0.9929452562, 0.8980168866 and
  0.3196949178; none is rounded up to an economically favourable conclusion.
- Municipal tables reconcile to 1,384,739 footprint records and all 24 counties.
  Longnan debris-flow NC remains 42,497 occupied grid cells.

These checks do not rerun susceptibility fitting, 28.84 million runout paths,
commercial FE solving or raw municipal footprint preparation. Requirements for
those stages are listed in INPUTS.md. The archive is not presented as a new
agent-efficiency benchmark.
