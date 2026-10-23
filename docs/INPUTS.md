# Inputs required beyond this archive

## Regional model

Set GANSU_SPATIAL_INPUT to a directory containing `models/grid_arrays.npz`.
The arrays expected by common.py are X, y, pid, fold, block, rr, cc, mask,
features and transform. Read the loader for exact array indexing and types.
X columns follow the `features` array: elevation, slope, periodic aspect
components, relief and TPI. y is historical inventory membership; pid preserves
positive polygon identity; the geographic block and fold arrays are fixed
before fitting. The source study used 30 m cells in UTM 48N (EPSG:32648).
Source inventory membership is not a complete absence/presence survey.

The baseline `fit_probability.py` additionally needs `y_any` in the grid NPZ
and the original polygon shapefile at `raw/commondata/arcgis/Landslide.shp`
(with its companion SHX, DBF and PRJ files). It creates baseline predictions,
per-fold estimator/calibrator pairs and the logistic coefficient table in the
input-project directory. `compare_cv.py`, `factor_analysis.py` and subsequent
extended analyses consume those baseline products. Therefore a grid array alone
is insufficient to rerun all stages; use the original geospatial inputs and run
the baseline first. No downloaded pickled model should be loaded from an
untrusted source.

## Runout

Set GANSU_RAW_INPUT to the original project-input root containing:

- tianshui_probability_20260909/gis/elevation_m.tif
- tianshui_probability_20260909/gis/study_mask_30m.tif
- tianshui_probability_20260909/gis/landslides_2013_WGS84.geojson
- three_regions_20260909/tianshui/osm_exposure_clipped_WGS84.geojson

GeoJSON polygons require `polygon_id`; footprint features require `kind` and
`osm_id`. Raw geometries and the source DEM are omitted. Source provider rights
must be checked before redistribution. Output root: GANSU_RUNOUT_OUTPUT.

## Finite elements

Archived eligible-case INP decks can be submitted to a separately licensed
Abaqus installation. The base generator writes to GANSU_FE_OUTPUT. Historical
ODB extraction needs GANSU_FE_NEW/GANSU_FE_OLD with their manifest, mesh JSON,
history/metadata and ODB files. These large solver outputs are not included.
Model and load units are SI (m, N, kg, s); response-table IM units are retained
explicitly as kPa_base, kPa_mean_gross or kN_s. Different IM definitions are
not interchangeable.

## Municipal scenarios

Set GANSU_MUNICIPAL_ROOT to a project directory with `data/region_definition.json`,
`data/city_counties_WGS84.geojson`, the archived archetype interval and retained
assignment CSVs, and `raw/{Lanzhou,Tianshui,Longnan}_municipality_bbox.parquet`.
The original Overture buildings release is 2026-09-23.1, accessed 2026-10-20.
Provider documentation: https://docs.overturemaps.org/guides/buildings/ .
Check availability and upstream release version before downloading. Current
footprints should not silently replace the archived source version.
The municipal metadata includes source hashes and regional bounds.

Third-party source dataset referenced by the manuscript:
https://doi.org/10.5281/zenodo.8174931 . Complete provider datasets, current
downloads and licence terms must be obtained from the original providers.
