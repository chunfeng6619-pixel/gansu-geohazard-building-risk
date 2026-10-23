"""Replay the independent annual-loss example using the original RNG sequence."""
from pathlib import Path
import argparse
import csv
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from gansu_research.annual import BASELINE, HIGH, LOW, LR, PARAMS, evaluate


def write_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'results/annual')
    parser.add_argument('--with-sobol', action='store_true')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    archived = ROOT/'data/D6_Annual_example'
    reference = json.loads((archived/'annual_risk_audit.json').read_text(encoding='utf-8-sig'))
    baseline = evaluate(BASELINE, details=True)
    eal = float(baseline['EAL'][0])
    np.testing.assert_allclose(eal, reference['EAL'], rtol=1e-12)
    assert evaluate(BASELINE, arrival_mult=0)[0] == 0
    zero_rate = BASELINE.copy()
    zero_rate[0] = 0
    assert evaluate(zero_rate)[0] == 0

    rng = np.random.default_rng(20260912)
    inputs = LOW + (HIGH-LOW)*rng.random((100000, 7))
    parameter_eal = np.concatenate([evaluate(batch) for batch in np.array_split(inputs, 25)])
    with np.load(archived/'uncertainty_100000_parameters_EAL.npz') as old:
        np.testing.assert_array_equal(inputs, old['inputs'])
        np.testing.assert_allclose(parameter_eal, old['EAL'], rtol=1e-12, atol=1e-9)
    np.savez_compressed(args.output/'parameter_uncertainty.npz', inputs=inputs,
                        EAL=parameter_eal, parameter_names=np.array(PARAMS))

    # Sobol evaluation is deterministic and consumes no state from the main RNG.
    if args.with_sobol:
        from scipy.stats import qmc
        ab = qmc.Sobol(14, scramble=True, seed=20260913).random_base2(13)
        a = LOW+(HIGH-LOW)*ab[:, :7]
        b = LOW+(HIGH-LOW)*ab[:, 7:]
        ya, yb = evaluate(a), evaluate(b)
        variance = np.var(np.r_[ya, yb], ddof=1)
        sobol = []
        for j, name in enumerate(PARAMS):
            mixed = a.copy()
            mixed[:, j] = b[:, j]
            yc = evaluate(mixed)
            sobol.append(dict(parameter=name,
                S1=float(np.mean((yb-yb.mean())*(yc-ya))/variance),
                ST=float(np.mean((ya-yc)**2)/(2*variance))))
        with (archived/'Sobol_indices.csv').open(encoding='utf-8-sig') as f:
            old = list(csv.DictReader(f))
        for row, oldrow in zip(sobol, old):
            assert row['parameter'] == oldrow['parameter']
            np.testing.assert_allclose([row['S1'], row['ST']],
                [float(oldrow['S1']), float(oldrow['ST'])], rtol=1e-10, atol=1e-12)
        write_csv(args.output/'Sobol_indices.csv', sobol)

    ny = 200000
    counts = rng.poisson(BASELINE[0], ny)
    year_index = np.repeat(np.arange(ny), counts)
    state = rng.choice(5, size=len(year_index), p=baseline['state_prob'][0])
    loss = BASELINE[-1]*LR[state]
    annual = np.bincount(year_index, weights=loss, minlength=ny)
    simulated = float(annual.mean())
    se = float(annual.std(ddof=1)/np.sqrt(ny))
    np.testing.assert_allclose(simulated, reference['MC_annual_mean'], rtol=1e-12)
    assert abs(simulated-eal) < 5*se

    pv = sum((1+.03)**(-t) for t in range(1, 31))
    interventions = [
        ('capacity_x1.5', dict(capacity_mult=1.5), 30000),
        ('energy_head_minus_10m', dict(energy_removed=10), 50000),
        ('arrival_x0.1', dict(arrival_mult=.1), 250000)]
    economics = []
    for name, kw, cost in interventions:
        after = float(evaluate(BASELINE, **kw)[0])
        benefit = (eal-after)*pv
        present_cost = cost+.005*cost*pv
        economics.append(dict(scenario=name, before_EAL=eal, after_EAL=after,
            PV_benefit=benefit, PV_cost=present_cost, BCR=benefit/present_cost))
        assert 0 <= after < eal and benefit/present_cost < 1
    with (archived/'mitigation_cost_benefit.csv').open(encoding='utf-8-sig') as f:
        old = list(csv.DictReader(f))
    for row, oldrow in zip(economics, old):
        np.testing.assert_allclose([row['after_EAL'], row['BCR']],
            [float(oldrow['after_EAL']), float(oldrow['BCR'])], rtol=1e-12)
    write_csv(args.output/'mitigation_cost_benefit.csv', economics)
    result = dict(scope='independent_parametric_example', analytic_EAL_CNY_per_year=eal,
        simulated_EAL_CNY_per_year=simulated, MC_standard_error=se,
        parameter_samples=100000, simulated_years=ny,
        uncertainty_quantiles=np.quantile(parameter_eal, [.05, .5, .95]).tolist(),
        archived_parameter_draws_match=True, archived_numerical_results_match=True,
        mitigation_BCR=[r['BCR'] for r in economics])
    (args.output/'reproduction_audit.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
