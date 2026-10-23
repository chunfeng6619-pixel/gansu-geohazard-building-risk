# Manuscript and supplementary correspondence

The study snapshot is manuscript `1022`, not an earlier paper. Section and result
names below are used to avoid ambiguity if figure numbers change at submission.

| Main-text evidence | Supplement/data | Code or archived model input |
|---|---|---|
| Susceptibility verification and terrain contribution | D1; supplementary validation tables | research/spatial/common.py, compare_cv.py, factor_analysis.py |
| ALE and spatial uncertainty methods | D1 ALE table | research/spatial/factor_analysis.py, bootstrap_space.py |
| 44 candidates and source–building arrival | D2 | research/runout/gansu_risk_v11_mc.py |
| Nine archetypes, load families and quality decisions | D3; supplementary model/load/quality tables | research/finite_elements/model_builder.py; data/FE_input_decks/*.inp |
| 81 capacity constraints and censored cases | D4 | scripts/reproduce_fragility_bounds.py |
| Signed local sensitivity | D5 and original eight perturbations | scripts/check_local_sensitivity.py |
| Annual risk, parameter uncertainty and economics | D6 | src/gansu_research/annual.py; scripts/reproduce_annual_example.py |
| Municipal archetype/conditional-damage maps | D7 | research/spatial_scenarios/prepare_municipal_data.py |

The included `data/README.txt` is the supporting-data description prepared for
the 1022 Supplementary Material. Table columns, units and existing labels are
retained. Missing values, right-censored upper capacities and NC are preserved.
The full supplementary Word document is not published in this software archive.
