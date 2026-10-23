"""Check immutable files, scientific counts, NC handling and conditional probabilities."""
from collections import Counter
from pathlib import Path
import csv
import gzip
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data'


def rows(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def main():
    manifest = json.loads((ROOT/'MANIFEST_SHA256.json').read_text(encoding='utf-8'))
    for name, digest in manifest.items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
    for row in rows(ROOT/'provenance/source_files.csv'):
        assert hashlib.sha256((ROOT/row['file']).read_bytes()).hexdigest() == row['archive_sha256']

    candidates = rows(DATA/'D2_Runout/building_arrivals.csv')
    pairs = rows(DATA/'D2_Runout/source_building_pairs.csv')
    runout = json.loads((DATA/'D2_Runout/runout_audit.json').read_text(encoding='utf-8-sig'))
    ns = int(runout['sources'])
    assert len(candidates) == 44 and ns == 5768
    direct = sum(int(r['arrival_samples']) > 0 for r in candidates)
    near = sum(float(r['P_neighbourhood_equal_source']) > 0 for r in candidates)
    assert (direct, near) == (6, 7)
    for r in candidates:
        hits = sum(int(p['hits']) for p in pairs if p['osm_id'] == r['osm_id'])
        assert hits == int(r['arrival_samples'])
        assert math.isclose(hits/(ns*5000), float(r['P_reach_equal_source']), abs_tol=1e-15)
        assert r['annual_rate'] == r['EAL'] == r['whole_building_severe_probability'] == ''
    for p in pairs:
        assert int(p['trials']) == 5000
        assert math.isclose(int(p['hits'])/5000, float(p['P_reach_given_source']), abs_tol=1e-15)

    planned = rows(DATA/'D3_Finite_elements/all_117_planned_positions.csv')
    fe_counts = Counter(r['status'] for r in planned)
    assert len(planned) == 117 and fe_counts['qualified'] == 96 and fe_counts['NC'] == 21
    assert len(rows(DATA/'D3_Finite_elements/extracted_response_positions.csv')) == 115
    capacity = rows(DATA/'D4_Fragility/capacity_constraints.csv')
    cap_counts = Counter(r['status'] for r in capacity)
    assert cap_counts == dict(bracketed=43, right_censored=36, left_censored=2)
    for row in capacity:
        if row['status'] == 'right_censored':
            assert row['capacity_upper'] == ''

    municipal = []
    for city in ['Lanzhou', 'Tianshui', 'Longnan']:
        meta = json.loads((DATA/f'D7_Spatial_scenarios/{city}_metadata.json').read_text(encoding='utf-8-sig'))
        count = footprints = 0
        nc = Counter()
        with gzip.open(DATA/f'D7_Spatial_scenarios/{city}_100m_occupied_cells.csv.gz', 'rt', encoding='utf-8-sig') as f:
            for row in csv.DictReader(f):
                count += 1
                n = int(row['n_footprints'])
                footprints += n
                assert sum(int(row[f'n_type{t}']) for t in [1, 2, 3]) == n and n > 0
                for hazard in ['SLIDE', 'DEBRIS', 'ROCK']:
                    if int(row[f'{hazard}_NC_cell']):
                        nc[hazard] += 1
                        assert row[f'{hazard}_mean_p_lower'] == row[f'{hazard}_mean_p_upper'] == ''
                    else:
                        lo, hi = float(row[f'{hazard}_mean_p_lower']), float(row[f'{hazard}_mean_p_upper'])
                        assert 0 <= lo <= hi <= 1
        assert count == meta['n_occupied_cells'] and footprints == meta['n_buildings']
        assert sum(r['count'] for r in meta['county_counts']) == footprints
        assert meta['type_proportions_are_observed'] is False
        for hazard in ['SLIDE', 'DEBRIS', 'ROCK']:
            assert nc[hazard] == meta['status_counts'][hazard]['NC_cells']
        municipal.append(dict(city=city, footprint_records=footprints,
            occupied_cells=count, counties=meta['n_counties'], NC_cells=dict(nc)))
    assert sum(r['footprint_records'] for r in municipal) == 1384739
    report = dict(passed=True, verified_files=len(manifest), source_building_pairs=len(pairs),
        direct_buildings=direct, additional_neighbourhood_buildings=near-direct,
        FE_positions=dict(fe_counts), capacity_status=dict(cap_counts), municipalities=municipal)
    out = ROOT/'results/archive_audit.json'
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
