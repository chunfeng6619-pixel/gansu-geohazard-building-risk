from common import *
from scipy.spatial.distance import cdist
from shapely.geometry import shape,mapping,box
from shapely.ops import transform as geotransform
from pyproj import Transformer
import time
start=time.time();p=np.load(B/'models/predictions.npz')['boost_calibrated']
calrows=[];binrows=[]
for protocol in ['random','spatial_0','spatial_300','spatial_1000','random_size_matched_300']:
    pred=np.load(E/'arrays'/f'cv_{protocol}.npz')
    for model in ['logistic','boost_raw','boost_calibrated']:
        prob=pred[model]
        for k in [10,20,50]:
            bins=np.asarray(pd.qcut(prob,k,labels=False,duplicates='drop'));n=np.bincount(bins);mean=np.bincount(bins,weights=prob)/n;obs=np.bincount(bins,weights=y)/n;weights=n/len(y)
            ece=np.sum(weights*np.abs(obs-mean));mce=np.max(np.abs(obs-mean));rel=np.sum(weights*(obs-mean)**2);res=np.sum(weights*(obs-y.mean())**2);unc=y.mean()*(1-y.mean());binned=np.mean((mean[bins]-y)**2);raw=np.mean((prob-y)**2)
            assert abs(binned-(unc-res+rel))<1e-12
            calrows.append({'protocol':protocol,'model':model,'equal_frequency_bins':k,'ece':ece,'mce':mce,'uncertainty':unc,'resolution':res,'reliability':rel,'brier_grouped_predictions':binned,'brier_raw':raw,'within_bin_correction':raw-binned,'identity_error':raw-(unc-res+rel+(raw-binned))})
            if k==10:
                binrows.extend({'protocol':protocol,'model':model,'bin':j,'cells':int(n[j]),'mean_probability':mean[j],'observed_fraction':obs[j]} for j in range(len(n)))
csv('calibration_extended_metrics.csv',pd.DataFrame(calrows));csv('calibration_10bin_details.csv',pd.DataFrame(binrows))

# Empirical area-budget capture, using mapped positive grid centres as denominator.
rng=np.random.default_rng(4897);bix=np.searchsorted(blocks,block);btotal=np.bincount(bix,weights=y,minlength=len(blocks));bsel=rng.integers(0,len(blocks),size=(2000,len(blocks)));bd=btotal[bsel].sum(axis=1)
eff=[]
pred=np.load(E/'arrays/cv_spatial_300.npz')
for model in ['logistic','boost_raw','boost_calibrated']:
    prob=pred[model];order=np.argsort(-prob,kind='stable');cy=np.concatenate([[0],np.cumsum(y[order])]);cp=np.concatenate([[0],np.cumsum(prob[order])])
    for fraction in np.linspace(0,1,101):
        count=min(int(np.ceil(fraction*len(y))),len(y));idx=order[:count];captured=int(cy[count]);perblock=np.bincount(bix[idx],weights=y[idx],minlength=len(blocks));boot=perblock[bsel].sum(axis=1)/bd
        lo,hi=np.percentile(boot,[2.5,97.5]);actualfraction=count/len(y)
        eff.append({'model':model,'area_budget_fraction':actualfraction,'area_budget_percent_nominal':int(round(fraction*100)),'selected_cells':count,'selected_area_km2':count*.0009,'captured_positive_cells':captured,'captured_fraction':captured/y.sum(),'lift_vs_uniform_area':captured/y.sum()/actualfraction if count else np.nan,'observed_fraction_selected':captured/count if count else np.nan,'predicted_mean_selected':cp[count]/count if count else np.nan,'fixed_selected_region_bootstrap_low':lo,'fixed_selected_region_bootstrap_high':hi,'score_cutoff_not_rainfall_threshold':prob[order[count-1]] if count else np.nan})
csv('investigation_efficiency_curve.csv',pd.DataFrame(eff));csv('investigation_efficiency_key_budgets.csv',pd.DataFrame(eff).query('area_budget_percent_nominal in [5,10,20,30,40,50]'))

# Block mean OOF residual and predeclared grid adjacency.
agg=pd.DataFrame({'block_id':block,'residual':y-p,'observed':y,'probability':p,'easting':xx,'northing':yy}).groupby('block_id').agg(cells=('observed','size'),positive_cells=('observed','sum'),observed_fraction=('observed','mean'),mean_probability=('probability','mean'),mean_residual=('residual','mean'),easting_m=('easting','mean'),northing_m=('northing','mean')).reset_index()
coords=np.column_stack([blocks%100,blocks//100]);dv=np.abs(coords[:,None,:]-coords[None,:,:])
def bh(p):
    n=len(p);order=np.argsort(p);q=np.minimum.accumulate((p[order]*n/np.arange(1,n+1))[::-1])[::-1];out=np.empty(n);out[order]=np.minimum(q,1);return out
results=[];lisas=[];perms=9999;rng=np.random.default_rng(69009)
for scheme in ['rook','queen']:
    adj=(dv.sum(axis=2)==1) if scheme=='rook' else (dv.max(axis=2)==1)
    degree=adj.sum(axis=1);weights=adj/np.maximum(degree[:,None],1)
    v=agg.mean_residual.to_numpy();z=v-v.mean();den=np.dot(z,z);n=len(z);S0=weights.sum();expected=-1/(n-1)
    I=n/S0*np.dot(z,weights@z)/den
    simulated=np.empty(perms)
    for i in range(perms):
        zz=rng.permutation(z);simulated[i]=n/S0*np.dot(zz,weights@zz)/den
    globalp=(1+np.sum(np.abs(simulated-expected)>=abs(I-expected)))/(perms+1)
    local=(n-1)*z*(weights@z)/den;pv=[]
    for i in range(n):
        other=np.delete(z,i);k=degree[i]
        if not k:pv.append(1.);continue
        sel=np.argsort(rng.random((perms,n-1)),axis=1)[:,:k];lag=other[sel].mean(axis=1);sim=(n-1)*z[i]*lag/den;ei=-z[i]**2/den
        pv.append((1+np.sum(np.abs(sim-ei)>=abs(local[i]-ei)))/(perms+1))
    pv=np.array(pv);q=bh(pv);qby=np.minimum(q*np.sum(1/np.arange(1,n+1)),1)
    lag=weights@z;quad=np.where(z>=0,np.where(lag>=0,'HH','HL'),np.where(lag>=0,'LH','LL'))
    for i,b in enumerate(blocks):lisas.append({'adjacency':scheme,'block_id':int(b),'neighbors':int(degree[i]),'global_moran_i':I,'global_permutation_p_two_sided':globalp,'local_moran_i':local[i],'conditional_permutation_p_two_sided':pv[i],'bh_q':q[i],'by_q':qby[i],'quadrant':quad[i],'lisa_bh_005':quad[i] if q[i]<.05 else 'not_significant','mean_residual':v[i],'centered_residual':z[i],'spatial_lag_residual':lag[i]})
    results.append({'adjacency':scheme,'n':n,'islands':int((degree==0).sum()),'moran_i':I,'expected_under_randomization':expected,'permutations':perms,'global_p_two_sided':globalp,'local_unadjusted_below_005':int((pv<.05).sum()),'local_bh_below_005':int((q<.05).sum()),'local_by_below_005':int((qby<.05).sum())})
    np.save(E/'arrays'/f'moran_global_null_{scheme}.npy',simulated)
csv('spatial_residual_block_summary.csv',agg);csv('moran_global_results.csv',pd.DataFrame(results));csv('lisa_conditional_permutation_fdr.csv',pd.DataFrame(lisas))
# GIS polygons of nominal blocks clipped to the actual analysis footprint.
foot=shape(json.loads((B/'gis/analysis_footprint_WGS84.geojson').read_text(encoding='utf-8'))['features'][0]['geometry']);forward=Transformer.from_crs(4326,32648,always_xy=True).transform;back=Transformer.from_crs(32648,4326,always_xy=True).transform;foot=geotransform(forward,foot)
bx0=int(np.floor(xx/4000).min());by0=int(np.floor(yy/4000).min());lisa=pd.DataFrame(lisas).query("adjacency=='rook'").set_index('block_id');features_out=[]
for _,row in agg.iterrows():
    b=int(row.block_id);x0=(b%100+bx0)*4000;y0=(b//100+by0)*4000;geom=box(x0,y0,x0+4000,y0+4000).intersection(foot)
    props={**row.to_dict(),**lisa.loc[b].to_dict()};features_out.append({'type':'Feature','geometry':mapping(geotransform(back,geom)),'properties':props})
dump('gis/residual_lisa_blocks_WGS84.geojson',{'type':'FeatureCollection','features':features_out})
dump('diagnostics_results.json',{'calibration':'ECE/MCE equal-frequency bins at 10,20,50; joint intercept/slope and intercept-only offset diagnostic in CV table','brier_decomposition':'Exact identity for bin-mean predictions: uncertainty-resolution+reliability; raw-probability score includes explicit within-bin correction','capture_definition':'Fraction of the 8654 mapped positive 30m centres captured by ranked area; NOT fraction of 5923 distinct landslides, and no monetary cost/benefit assumed','capture_interval':'2000 geographic block resamples of fixed selected regions; resampling region area can differ from nominal budget; intervals conditional on already-fitted OOF map','moran_residual':'Unweighted block mean(y-p), one value per each of 47 blocks; descriptive exploratory inference, different block areas and cell counts create heteroskedasticity','moran_test':'9999 global permutations; local conditional randomization holds own residual fixed, samples neighbors without replacement; two-sided extremeness about theoretical conditional mean','multiple_testing':'BH q at 0.05 is displayed; BY q also supplied for stronger dependence-robust correction','moran_results':results,'seconds':time.time()-start})
print(json.dumps(results),flush=True)
