# Gansu geohazard and rural-building risk: research archive

Research code and supporting numerical data for:

**A Multi-Agent Framework for Regional Landslide Susceptibility Mapping and Multi-Hazard Fragility Assessment of Rural Buildings**

Hao Wang, Chen Shi, Yawei Li and Hongrong Ma.
Repository: https://github.com/chunfeng6619-pixel/gansu-geohazard-building-risk

This release accompanies manuscript `1022` and its Supplementary Material.
It combines archived numerical evidence with portable calculation and audit
commands. It is a research snapshot, not a claim that all regional hazards or
actual buildings have calibrated annual losses.

## Quick start

Python 3.10 or later, in a fresh virtual environment:

```bash
python -m venv .venv
# Activate .venv with the command appropriate to your operating system.
python -m pip install -e .
python scripts/audit_archive.py
python scripts/reproduce_annual_example.py --output results/annual
python scripts/reproduce_fragility_bounds.py --output results/fragility
python scripts/check_local_sensitivity.py
```

Run these commands from the repository root. The audit verifies file hashes,
source–building counts, the 117-position FE quality ledger, capacity censoring,
municipal counts and probability closure. Output is written only to `results/`.
The annual benchmark regenerates its original 100,000 parameter draws and
200,000 compound-Poisson years with seed 20260912.

## What can be reproduced

| Evidence | Files | Scope |
|---|---|---|
| Spatial model validation, terrain contributions and ALE | `data/D1_Spatial_model`, `research/spatial` | Archived tables are fully accessible; refitting additionally requires the source grid arrays |
| Conditional source-to-building arrival | `data/D2_Runout`, `research/runout` | 5,768 historical source areas × 5,000 paths; direct and neighbourhood arrivals kept distinct |
| Structural response and quality state | `data/D3_Finite_elements`, `data/FE_input_decks`, `research/finite_elements` | 115 extracted positions; 117 planned, 96 qualified and 21 NC; 50 eligible-case INP decks, not 117 independent INP files |
| Capacity intervals and fragility | `data/D4_Fragility`, `scripts/reproduce_fragility_bounds.py` | Local damage thresholds; interval and censored capacities retained |
| Local signed sensitivity | `data/D5_Local_sensitivity`, `scripts/check_local_sensitivity.py` | Central differences, ±2.5% around the LZ1 medium debris-flow case |
| Independent annual-risk example | `src/gansu_research/annual.py`, `data/D6_Annual_example` | Analytic EAL, parameter uncertainty and mitigation economics; independent of mapped assets |
| Municipal spatial scenarios | `data/D7_Spatial_scenarios`, `research/spatial_scenarios` | Full-municipality 100 m occupied cells; archetype assignments are scenarios, not surveyed attributes |

See [REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md),
[MANUSCRIPT_MAPPING.md](docs/MANUSCRIPT_MAPPING.md) and
[INPUTS.md](docs/INPUTS.md) for exact boundaries and additional inputs.
`provenance/source_files.csv` records original and archived file hashes and the
path-only portability adaptations. The raw source files are retained separately
from newly added release-audit utilities.

## Interpretation

1. Spatial susceptibility estimates membership of the historical inventory under
   the selected terrain features. It is not an annual occurrence probability.
2. Runout arrival is conditional on activation of a historical source under the
   declared routing and motion-angle priors. Velocity is an energy proxy, not a
   calibrated pressure input for all nine FE archetypes.
3. Nine archetypes are engineering scenarios. The three prescribed load families
   are landslide, debris flow and equivalent rockfall impact. Local tensile
   damage D_t and EDI do not directly measure whole-building collapse.
4. Capacity dispersion and local damage thresholds are modelling assumptions.
   NC positions stay NC and are never filled with zero damage.
5. The annual benchmark is independent: baseline analytic EAL is
   4,974.742362 CNY/year and simulated EAL is 4,989.45 CNY/year. None is assigned
   to a real footprint. The three cost–benefit ratios are all below one under
   the benchmark cost, lifetime and discount assumptions.
6. Municipal building classes are explicit scenario assignments. Cells with
   any NC-assigned footprint remain NC; empty cells are NoData.

## Licensing and citation

Code is MIT. Data rights differ by component; see [DATA_LICENSES.md](DATA_LICENSES.md).
Do not drop Overture/OSM attribution from reused spatial derivatives.
Use `CITATION.cff` to cite version 1.0.0. A repository URL is not a DOI; the DOI
will be added only after a public Zenodo record has actually been published.

中文说明：[README_zh.md](README_zh.md).
