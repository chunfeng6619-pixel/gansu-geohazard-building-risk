# Reproducibility levels, version 1.0.1

The release distinguishes numerical replay from rerunning external solvers.

**Runs with included files:** archive checksums and scientific consistency audit;
analytic annual benchmark, original parameter draws and compound-Poisson Monte
Carlo; mitigation present values and BCR; local central sensitivities; conditional
fragility interval values evaluated at every archived intensity and assumed beta.

**Needs original geospatial inputs:** spatial-model retraining, spatial bootstrap
and the 28.84 million runout paths. `research/spatial/common.py` now reads
`GANSU_SPATIAL_INPUT` and writes `GANSU_SPATIAL_OUTPUT`. Install
`requirements-spatial.txt` before using the source research scripts. Their
published statistics are archived, but a complete raw-input reanalysis is not
claimed by the quick-start commands.

**Needs licensed Abaqus and model inputs:** solving the 50 archived eligible-case
INP files and extracting ODB fields. The source builder is retained as the base
model generator; the archived INP decks are authoritative for those individual
cases. The base builder alone does not recreate every supplementary intensity,
recovery run or accepted extraction in the quality ledger. NC cases must not be
replaced by predictions from an unvalidated surrogate.

**Needs Overture geometry and county boundaries:** full spatial-scenario rebuild.
The source assignment algorithm and retained prior ID assignments are provided.
It uses 24 spatial waves, ID-stable mixing and fixed seeds; it does not infer
observed structure material or disaster occurrence frequency.

No new multi-agent efficiency experiment is introduced in this release.
The original analysis scripts implement professional calculation stages; the
portable verification commands are release utilities, not retrospective proof
that an autonomous agent performed or accelerated every calculation.

## Key settings preserved

- Susceptibility: HistGradientBoostingClassifier, learning rate 0.08, 100
  iterations, 15 leaves, min_samples_leaf 200, L2 10, no early stopping,
  seed 20260909; logistic control with train-fold standardisation.
- Runout: seed 20260911; 5,000 paths per source; downhill-neighbour weights
  proportional to slope^4; truncated-normal motion angle, mean 25°, SD 5°,
  bounds 10–40°; residual-energy velocity proxy.
- Local sensitivity: ±2.5%, central difference divided by baseline response.
- Annual example: seed 20260912; independent uniforms, 100,000 samples;
  64-point Gauss–Legendre integration; 200,000 simulated years. Its five
  economic states differ from the four local FE damage states.

Original copied code and data hashes: `provenance/source_files.csv`.
