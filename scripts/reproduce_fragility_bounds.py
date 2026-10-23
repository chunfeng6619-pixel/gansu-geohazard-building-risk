"""Re-evaluate archived local fragility bounds without inventing censored capacities."""
from pathlib import Path
import argparse
import csv
import json
import math
from scipy.special import ndtr

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'results/fragility')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    data = ROOT/'data/D4_Fragility'
    with (data/'capacity_constraints.csv').open(encoding='utf-8-sig') as f:
        capacity = {(r['prototype'], r['hazard'], float(r['local_Dt_trigger'])): r
                    for r in csv.DictReader(f)}
    with (data/'damage_exceedance_intervals.csv').open(encoding='utf-8-sig') as f:
        rows = list(csv.DictReader(f))
    errors = []
    output = []
    for row in rows:
        key = (row['prototype'], row['hazard'], float(row['local_Dt_trigger']))
        cap = capacity[key]
        lower = float(cap['capacity_lower'])
        upper = float(cap['capacity_upper']) if cap['capacity_upper'] else math.inf
        im = float(row['IM'])
        beta = float(row['assumed_log_capacity_dispersion'])
        assert im > 0 and beta > 0
        p_low = 0. if math.isinf(upper) else float(ndtr(math.log(im/upper)/beta))
        p_high = 1. if lower <= 0 else float(ndtr(math.log(im/lower)/beta))
        assert 0 <= p_low <= p_high <= 1
        errors.extend([abs(p_low-float(row['conditional_probability_lower'])),
                       abs(p_high-float(row['conditional_probability_upper']))])
        output.append(dict(prototype=key[0], hazard=key[1], local_Dt_trigger=key[2],
            IM=im, IM_unit=row['IM_unit'], assumed_beta=beta,
            capacity_status=cap['status'], probability_lower=p_low, probability_upper=p_high))
    assert max(errors) < 1e-12, max(errors)
    with (args.output/'recomputed_bounds.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    audit = dict(capacity_records=len(capacity), probability_rows=len(output),
                 maximum_absolute_difference=max(errors), all_archived_bounds_match=True)
    (args.output/'reproduction_audit.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
