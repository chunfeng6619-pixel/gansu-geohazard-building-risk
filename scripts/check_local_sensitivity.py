"""Verify signed central differences from the original accepted perturbation cases."""
from pathlib import Path
import csv
import json
import math

ROOT = Path(__file__).resolve().parents[1]


def main():
    with (ROOT/'data/6-5_original_valid_cases.csv').open(encoding='utf-8-sig') as f:
        cases = list(csv.DictReader(f))
    with (ROOT/'data/D5_Local_sensitivity/signed_central_sensitivity.csv').open(encoding='utf-8-sig') as f:
        sensitivities = list(csv.DictReader(f))
    base = next(r for r in cases if r['name'] == 'V8_ACC_LZ1')
    errors = []
    for row in sensitivities:
        metric = row['metric']
        selection = [r for r in cases if r['parameter'] == row['parameter']]
        lo = next(r for r in selection if float(r['factor']) < 1)
        hi = next(r for r in selection if float(r['factor']) > 1)
        b = float(base[metric])
        eps = float(row['relative_perturbation'])
        assert lo['accepted'] == hi['accepted'] == 'True'
        assert math.isclose(float(lo['factor']), 1-eps)
        assert math.isclose(float(hi['factor']), 1+eps)
        value = (float(hi[metric])-float(lo[metric]))/(2*eps*b)
        errors.append(abs(value-float(row['signed_central_sensitivity'])))
    assert max(errors) < 1e-12
    print(json.dumps(dict(verified_parameter_metric_pairs=len(errors),
        relative_perturbation=.025, maximum_absolute_difference=max(errors),
        passed=True), indent=2))


if __name__ == '__main__':
    main()
