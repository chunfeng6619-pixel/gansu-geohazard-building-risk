Supplementary Data for the 1022 manuscript

D1 Spatial model: pooled metrics, fold-specific split audit, terrain ablation, spatial permutation, coefficient stability and ALE responses.
D2 Runout: all 44 candidate buildings, source-to-building pairs and the runout calculation audit. Probabilities are source-conditioned or source-mixture quantities; they do not include event rates.
D3 Finite elements: 115 extracted responses, a 117-position quality ledger, force-resultant checks and model definitions. Two missing LN1 landslide positions are explicitly NC. Unqualified maxima remain in the raw extraction table but are excluded from the 117-position eligible-response column.
D4 Fragility: 81 capacity constraints, conditional exceedance bounds, 324 EDI interval comparisons and summaries of the 10,000-draw assumption sensitivity for fully bracketed combinations. Empty upper capacity denotes right censoring; zero lower capacity denotes left censoring. These are local structural scenarios.
D5 Local sensitivity: eight independently perturbed responses and central-difference coefficients about the same LZ1 medium debris-flow baseline.
D6 Annual example: fixed and prior assumptions in the audit, 100,000 parameter/EAL samples in NPZ format, loss-exceedance curves, Sobol estimates and mitigation economics. Annual event rates and economic inputs are independent scenario assumptions. This example is not attached to mapped footprints.
D7 Spatial scenarios: 24 county footprint counts, three municipal metadata files, all occupied 100 m cell summaries, retained earlier type assignments and 27 archetype probability intervals. Full municipal footprint geometry can be recovered from the public release named below; it is not duplicated here. Grid means weight footprint records equally; any NC-assigned footprint makes the cell NC. Unoccupied cells are NoData.

The bilingual Word files have identical section, table, figure, equation and data numbering. The accompanying cross-reference CSV records the main-text correspondence. In particular, Eq. (S8) is the Poisson time-window conversion cited in the main manuscript.

Units and formats: probabilities in CSV follow each column's own naming or metadata, not a blanket percent conversion; figures may display percentages. Structural intensity is kPa_base, kPa_mean_gross or kN_s as indicated. NPZ is a NumPy archive; CSV uses UTF-8. Occupied-cell summaries are gzip-compressed CSV. Numeric scientific records retain the archived precision; no new FE analysis or spatial probability model was run for this supplement.

Footprint attribution: Overture Maps Foundation, release 2026-09-23.1, incorporating OpenStreetMap and the dataset cited in main-text Fig. 8 (DOI:10.5281/zenodo.8174931). © Overture Maps Foundation | © OpenStreetMap contributors; ODbL. Full source references follow the main manuscript.

中文说明：D1—D7分别提供空间模型、随机传播、有限元、条件易损性、局部敏感性、独立年度算例和市域空间情景的数值依据。请按各列字段及元数据识别概率、百分数和单位；NC不作为零损伤或零概率。完整市域的实际结构属性未知，类型为空间情景配置。独立算例EAL不赋给真实建筑。八张图为补充材料内部支持图；正文1022文件保持不变。
