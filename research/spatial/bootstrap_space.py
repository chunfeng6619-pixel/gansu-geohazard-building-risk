"""500 complete block-bootstrap fits; spatially isolated OOB calibration.
Resampling interval is conditional on the existing inventory, blocks and model.
"""
from common import *
import time
from concurrent.futures import ThreadPoolExecutor,as_completed
N=500
path=E/'arrays/bootstrap_probabilities_500x531538.float32'
resume=path.exists() and (E/'tables/bootstrap_replicates.csv').exists()
out=np.memmap(path,dtype='float32',mode='r+' if resume else 'w+',shape=(N,len(y)))
seeds=np.arange(N)+2026091000
block_index=np.searchsorted(blocks,block)
rows=pd.read_csv(E/'tables/bootstrap_replicates.csv').to_dict('records') if resume else []
done={int(r['replicate'])-1 for r in rows};start=time.time()
def one(j):
        seed=seeds[j]
        rng=np.random.default_rng(int(seed));counts=np.bincount(rng.integers(0,len(blocks),len(blocks)),minlength=len(blocks))
        cal=counts[block_index]==0
        fit=(~cal)&eligible(cal,300)
        wt=counts[block_index][fit].astype(float)
        if y[fit].sum()<100 or y[cal].sum()<100:raise RuntimeError('Insufficient fit/calibration events; do not silently accept invalid replicate')
        model=gb();model.fit(X[fit],y[fit],sample_weight=wt)
        c=calibrate(model,X[cal],y[cal]);pair={'estimator':model,'calibrator':c}
        pp=predict(pair,X).astype('float32')
        info={'replicate':j+1,'seed':int(seed),'unique_inbag_blocks':int((counts>0).sum()),'oob_blocks':int((counts==0).sum()),'fit_cells':int(fit.sum()),'weighted_fit_cells':float(wt.sum()),'fit_positive':int(y[fit].sum()),'calibration_cells':int(cal.sum()),'calibration_positive':int(y[cal].sum()),'map_mean_probability':float(pp.mean()),'calibration_logit_intercept':float(c.intercept_[0]),'calibration_logit_slope':float(c.coef_[0,0])}
        if j%25==0:joblib.dump(pair,E/'models'/f'bootstrap_checkpoint_{j+1:03d}.joblib',compress=3)
        return j,pp,info
with threadpool_limits(limits=2):
    with ThreadPoolExecutor(max_workers=3) as pool:
        for future in as_completed([pool.submit(one,j) for j in range(N) if j not in done]):
            j,pp,info=future.result();out[j]=pp;rows.append(info)
            if len(rows)%10==0 or len(rows)==N:
                out.flush();csv('bootstrap_replicates.csv',pd.DataFrame(rows).sort_values('replicate'))
                dump('bootstrap_progress.json',{'completed':len(rows),'total':N,'elapsed_seconds_current_segment':time.time()-start})
                print(json.dumps({'bootstrap_completed':len(rows),'seconds_current_segment':round(time.time()-start,1)}),flush=True)
out.flush()
mean=np.empty(len(y),dtype='float32');quant=np.empty((3,len(y)),dtype='float32')
for begin in range(0,len(y),10000):
    end=min(begin+10000,len(y));z=np.asarray(out[:,begin:end]);mean[begin:end]=z.mean(axis=0);quant[:,begin:end]=np.percentile(z,[2.5,50,97.5],axis=0)
width=quant[2]-quant[0]
high=(mean>=np.quantile(mean,.8))&(width>=np.quantile(width,.8))
for name,value in [('bootstrap_mean',mean),('bootstrap_p025',quant[0]),('bootstrap_p50',quant[1]),('bootstrap_p975',quant[2]),('bootstrap_interval_width',width)]:raster(name,value)
raster('high_probability_high_uncertainty_top20',high.astype('uint8'),255,'uint8')
np.savez_compressed(E/'arrays/bootstrap_summary.npz',mean=mean,p025=quant[0],p50=quant[1],p975=quant[2],width=width,high_high=high)
csv('bootstrap_summary_by_block.csv',pd.DataFrame({'block':block,'mean':mean,'low':quant[0],'median':quant[1],'high':quant[2],'width':width,'high_high':high}).groupby('block').agg(mean_probability=('mean','mean'),mean_interval_width=('width','mean'),high_high_fraction=('high_high','mean'),cells=('mean','size')).reset_index())
dump('bootstrap_results.json',{'replicates':N,'geographic_blocks':len(blocks),'seconds':time.time()-start,'sampling':'47 original 4 km blocks drawn with replacement 47 times; occurrence weights, no artificial class balancing','calibration':'Entire out-of-bag blocks; in-bag fitting locations removed within 300 m or sharing a positive polygon with OOB calibration','mean_probability':float(mean.mean()),'median_interval_width':float(np.median(width)),'interval_width_q95':float(np.quantile(width,.95)),'high_high_area_km2':float(high.sum()*.0009),'high_high_fraction':float(high.mean()),'high_probability_cutoff':float(np.quantile(mean,.8)),'high_uncertainty_cutoff':float(np.quantile(width,.8)),'interval_interpretation':'Pointwise percentile resampling interval conditional on observed inventory and chosen spatial blocks; not simultaneous confidence band, not future/annual occurrence probability, no inventory completeness or DEM measurement uncertainty','raw_array':str(path),'raw_array_shape':[N,len(y)]})
print('BOOTSTRAP FINISHED',flush=True)
