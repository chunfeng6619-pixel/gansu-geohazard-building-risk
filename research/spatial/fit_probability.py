"""Spatially separated, natural-prevalence inventory probability models.

Primary: fixed-parameter boosted trees, sigmoid calibration on separate blocks.
Comparison: logistic regression and training-prevalence intercept-only model.
Five outer folds: 4 km blocks; 300 m exclusion around held-out cells; all
positive-centre cells of the same polygon excluded together from fitting.
No coordinates, class balancing, random pixel splits, or test tuning.
"""
from pathlib import Path
import sys, json, os, time
B=Path(os.environ.get('GANSU_SPATIAL_INPUT', 'external_inputs/spatial')).resolve()
for sub in ['tables','models','gis','qa']:(B/sub).mkdir(parents=True,exist_ok=True)
# Dependencies are installed using requirements-spatial.txt.
os.environ['OMP_NUM_THREADS']='4'
os.environ['OPENBLAS_NUM_THREADS']='4'
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import Affine
from rasterio.features import rasterize
import shapefile
from shapely.geometry import shape
from shapely import make_valid
from scipy.ndimage import distance_transform_edt
from scipy.special import logit
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss, log_loss, roc_curve, precision_recall_curve
import sklearn, joblib
from threadpoolctl import threadpool_limits

d=np.load(B/'models/grid_arrays.npz')
X=d['X'];y=d['y'];ya=d['y_any'];mask=d['mask'];rr=d['rr'];cc=d['cc'];pid=d['pid'];fold=d['fold'];block=d['block'];features=d['features'].tolist()
tf=Affine(*d['transform'][:6]);h,w=mask.shape
shape_reader=shapefile.Reader(str(B/'raw/commondata/arcgis/Landslide.shp'))
pid_any_grid=rasterize([(make_valid(shape(s.__geo_interface__)),i+1) for i,s in enumerate(shape_reader.shapes())],out_shape=mask.shape,transform=tf,all_touched=True,dtype='int32')
pid_any=pid_any_grid[mask]
def raster(name,v,fill=-9999,dtype='float32'):
    a=np.full(mask.shape,fill,dtype=dtype);a[mask]=v
    with rasterio.open(B/'gis'/f'{name}.tif','w',driver='GTiff',height=h,width=w,count=1,dtype=dtype,crs='EPSG:32648',transform=tf,nodata=fill,compress='deflate',tiled=True) as dst:
        dst.write(a,1)
        dst.update_tags(target='30m centre mapped in 2013 inventory',interpretation='Retrospective inventory membership conditional on footprint, NOT annual probability',model=name)

def eligible(hold,buffer_m=300,polygon_ids=None):
    held_grid=np.zeros(mask.shape,dtype=bool);held_grid[rr[hold],cc[hold]]=True
    dist=distance_transform_edt(~held_grid)*30
    keep=dist[mask]>buffer_m
    used_ids=pid if polygon_ids is None else polygon_ids
    hp=np.unique(used_ids[hold]);hp=hp[hp>0]
    keep&=~np.isin(used_ids,hp)
    return keep

def metrics(label,p):
    return dict(cells=int(len(label)),positives=int(label.sum()),observed_fraction=float(label.mean()),predicted_mean=float(p.mean()),roc_auc=float(roc_auc_score(label,p)),average_precision=float(average_precision_score(label,p)),brier=float(brier_score_loss(label,p)),log_loss=float(log_loss(label,p)))

names=['training_prevalence','logistic','boost_raw','boost_calibrated']
pred={n:np.full(len(y),np.nan) for n in names}
sens_any=np.full(len(y),np.nan);sens_buffer=np.full(len(y),np.nan)
ensemble=np.zeros(len(y)); disagreement=np.zeros((5,len(y)),dtype='float32');ood=np.zeros(len(y),dtype='uint8')
details=[]; per_fold=[]; coefs=[]; sens=[]
start=time.time()
with threadpool_limits(limits=4):
    for f in range(5):
        test=fold==f
        outer_train=(~test)&eligible(test,300)
        cal=(fold==(f+1)%5)&outer_train
        fit=outer_train&(fold!=(f+1)%5)&eligible(cal,300)
        assert not set(pid[fit][pid[fit]>0])&set(pid[test][pid[test]>0])
        assert not set(pid[fit][pid[fit]>0])&set(pid[cal][pid[cal]>0])
        assert not np.any(test&fit) and not np.any(test&cal) and not np.any(cal&fit)
        split=dict(outer_fold=f,fit_cells=int(fit.sum()),calibration_cells=int(cal.sum()),test_cells=int(test.sum()),logistic_train_cells=int(outer_train.sum()),fit_positive=int(y[fit].sum()),calibration_positive=int(y[cal].sum()),test_positive=int(y[test].sum()))
        details.append(split)
        pred['training_prevalence'][test]=y[outer_train].mean()
        lr=make_pipeline(StandardScaler(),LogisticRegression(C=1.0,max_iter=400,solver='lbfgs'))
        lr.fit(X[outer_train],y[outer_train]);pred['logistic'][test]=lr.predict_proba(X[test])[:,1]
        for name,coef in zip(features,lr[-1].coef_[0]):
            coefs.append(dict(outer_fold=f,feature=name,log_odds_per_training_sd=float(coef)))
        gb=HistGradientBoostingClassifier(learning_rate=.08,max_iter=100,max_leaf_nodes=15,min_samples_leaf=200,l2_regularization=10,early_stopping=False,random_state=20260909)
        gb.fit(X[fit],y[fit])
        calibrator=LogisticRegression(C=1e6,max_iter=200)
        cal_logit=gb.decision_function(X[cal]).reshape(-1,1)
        calibrator.fit(cal_logit,y[cal])
        pred['boost_raw'][test]=gb.predict_proba(X[test])[:,1]
        pred['boost_calibrated'][test]=calibrator.predict_proba(gb.decision_function(X[test]).reshape(-1,1))[:,1]
        allp=calibrator.predict_proba(gb.decision_function(X).reshape(-1,1))[:,1]
        ensemble+=allp/5;disagreement[f]=allp
        ood[test]=np.any((X[test]<X[fit].min(axis=0))|(X[test]>X[fit].max(axis=0)),axis=1).astype('uint8')
        joblib.dump({'estimator':gb,'calibrator':calibrator,'features':features,'split':split},B/'models'/f'boost_pair_fold_{f}.joblib',compress=3)
        for name in names:
            m=metrics(y[test],pred[name][test]);m.update(outer_fold=f,model=name)
            m['brier_skill_vs_train_prevalence']=1-m['brier']/brier_score_loss(y[test],pred['training_prevalence'][test])
            per_fold.append(m)
        # Two prespecified protocol sensitivities; descriptive, not model selection.
        alr=make_pipeline(StandardScaler(),LogisticRegression(C=1,max_iter=400))
        # Any-touch sensitivity uses its own polygon IDs for split exclusion.
        any_train=(~test)&eligible(test,300,polygon_ids=pid_any)
        alr.fit(X[any_train],ya[any_train]);sens_any[test]=alr.predict_proba(X[test])[:,1]
        train1000=(~test)&eligible(test,1000)
        blr=make_pipeline(StandardScaler(),LogisticRegression(C=1,max_iter=400))
        blr.fit(X[train1000],y[train1000]);sens_buffer[test]=blr.predict_proba(X[test])[:,1]
        print(json.dumps({'fold_complete':f,'seconds':round(time.time()-start,1),'primary':per_fold[-1]},ensure_ascii=False),flush=True)

assert all(np.isfinite(p).all() and np.all((p>=0)&(p<=1)) for p in pred.values())
summary={name:metrics(y,p) for name,p in pred.items()}
base=summary['training_prevalence']['brier']
for m in summary.values():m['brier_skill_vs_train_prevalence']=1-m['brier']/base
pd.DataFrame([{'model':n,**m} for n,m in summary.items()]).to_csv(B/'tables/model_metrics_oof.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(per_fold).to_csv(B/'tables/model_metrics_by_fold.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(details).to_csv(B/'tables/spatial_split_audit.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(coefs).to_csv(B/'tables/logistic_coefficients_by_fold.csv',index=False,encoding='utf-8-sig')
sens=[{'scenario':'centre_30m_300m_buffer','model':'logistic',**summary['logistic']},{'scenario':'any_touch_30m_300m_buffer','model':'logistic',**metrics(ya,sens_any)},{'scenario':'centre_30m_1000m_buffer','model':'logistic',**metrics(y,sens_buffer)}]
pd.DataFrame(sens).to_csv(B/'tables/sensitivity_metrics.csv',index=False,encoding='utf-8-sig')
p=pred['boost_calibrated']
for name,pr in pred.items():raster(f'probability_oof_{name}',pr)
raster('probability_ensemble_calibrated',ensemble)
raster('ensemble_model_spread_not_confidence_interval',disagreement.std(axis=0))
raster('oof_covariate_range_extrapolation',ood,255,'uint8')
np.savez_compressed(B/'models/predictions.npz',**pred,ensemble=ensemble,spread=disagreement.std(axis=0),sensitivity_any=sens_any,sensitivity_buffer1000=sens_buffer,ood=ood)
df=pd.DataFrame({'cell_id':rr*w+cc,'inventory_label':y,'block_id':block,'outer_fold':fold,'probability_oof_calibrated':p,'probability_ensemble_calibrated':ensemble,'model_spread_not_ci':disagreement.std(axis=0),'oof_outside_training_feature_range':ood})
df.to_csv(B/'tables/grid_probability_30m.csv.gz',index=False,encoding='utf-8-sig',compression='gzip')

# Equal-area ranking classes; 1 lowest, 5 highest. This is relative priority,
# not an official hazard/risk classification or a rainfall threshold.
q=pd.qcut(p,5,labels=False,duplicates='drop').astype(int)+1
raster('relative_priority_quintile_oof',q,255,'uint8')
df['priority_quintile']=q
tab=df.groupby('priority_quintile').agg(cells=('inventory_label','size'),positive_cells=('inventory_label','sum'),observed_fraction=('inventory_label','mean'),predicted_probability_mean=('probability_oof_calibrated','mean'),predicted_probability_min=('probability_oof_calibrated','min'),predicted_probability_max=('probability_oof_calibrated','max'))
tab['area_km2']=tab.cells*.0009;tab['positive_cell_capture_fraction']=tab.positive_cells/y.sum()
# Block-bootstrap intervals describe variation among current mapped blocks;
# they do not include incomplete inventory, future climate or domain uncertainty.
rng=np.random.default_rng(20260909)
blocks=np.unique(block)
bins=pd.qcut(p,10,labels=False,duplicates='drop').astype(int)
caldf=pd.DataFrame({'bin':bins,'block':block,'y':y,'p':p})
cal=caldf.groupby('bin').agg(cells=('y','size'),observed_fraction=('y','mean'),predicted_mean=('p','mean'))
boot_indices=rng.integers(0,len(blocks),size=(2000,len(blocks)))
for code,groups,table in [('priority_quintile',q,tab),('calibration_bin',bins,cal)]:
    for g in np.unique(groups):
        counts=np.array([np.sum((block==b)&(groups==g)) for b in blocks])
        sums=np.array([y[(block==b)&(groups==g)].sum() for b in blocks])
        denom=counts[boot_indices].sum(axis=1)
        rates=sums[boot_indices].sum(axis=1)/np.maximum(denom,1)
        lo,hi=np.percentile(rates[denom>0],[2.5,97.5])
        table.loc[g,'observed_fraction_block_bootstrap_low']=lo
        table.loc[g,'observed_fraction_block_bootstrap_high']=hi
tab.to_csv(B/'tables/priority_quintile_oof.csv',encoding='utf-8-sig')
cal.to_csv(B/'tables/calibration_bins_oof.csv',encoding='utf-8-sig')
blockstats=[]
for b in blocks:
    keep=block==b
    blockstats.append({'block_id':int(b),'n':int(keep.sum()),'positive':int(y[keep].sum()),'sum_p':float(p[keep].sum()),'brier_sum':float(np.square(p[keep]-y[keep]).sum()),'baseline_brier_sum':float(np.square(pred['training_prevalence'][keep]-y[keep]).sum())})
bs=pd.DataFrame(blockstats);bs.to_csv(B/'tables/validation_block_statistics.csv',index=False,encoding='utf-8-sig')
boot_skill=1-bs.brier_sum.to_numpy()[boot_indices].sum(axis=1)/bs.baseline_brier_sum.to_numpy()[boot_indices].sum(axis=1)
# Empirical score threshold for ranking only: fitted from the OOF distribution,
# not evaluated as an independently selected operational alarm threshold.
sort=np.argsort(-p);top20=sort[:int(np.ceil(len(y)*.2))]
metadata={'primary_model':'boost_calibrated; prespecified before evaluating folds','features':features,'hyperparameters':dict(learning_rate=.08,max_iter=100,max_leaf_nodes=15,min_samples_leaf=200,l2_regularization=10,early_stopping=False),'validation':'5 outer geographic folds, 4 km blocks, 300 m guards, whole positive polygons removed from fitting if present in held-out centre labels; separate spatial calibration fold','natural_prevalence':'All usable 30 m cells, no balanced sampling/class weights','calibration':'Logistic map of raw tree logits; separate spatial blocks; no test labels used in fitting/calibration','versions':{'numpy':np.__version__,'sklearn':sklearn.__version__,'rasterio':rasterio.__version__},'pooled_oof_metrics':summary,'brier_skill_block_bootstrap_95_interval':np.percentile(boot_skill,[2.5,97.5]).tolist(),'n_spatial_blocks':int(len(blocks)),'oof_feature_extrapolation_fraction':float(ood.mean()),'top20_area_capture_positive_centres':float(y[top20].sum()/y.sum()),'top20_observed_fraction':float(y[top20].mean()),'top20_mean_probability':float(p[top20].mean()),'ensemble_mean_probability':float(ensemble.mean()),'oof_probability_quantiles':{str(k):float(np.percentile(p,k)) for k in [0,10,25,50,75,90,95,99,100]},'uncertainty_note':'Bootstrap over mapped 4 km blocks (2000 replicates); not confidence in annual/future event probability. Spread of 5 overlapping fitted models is NOT a confidence interval.','scope':'Modelled inventory membership / event-context spatial susceptibility. No verified negative events, no annual rates, no loss or fatality estimate. Calibration is conditional on the dataset population, not proof of future calibration.','elapsed_seconds':time.time()-start}
(B/'model_results.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(metadata,ensure_ascii=False,indent=2),flush=True)
