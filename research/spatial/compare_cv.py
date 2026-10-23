from common import *
from sklearn.model_selection import StratifiedKFold
import time
random_folds=np.zeros(len(y),dtype='uint8')
for f,(_,idx) in enumerate(StratifiedKFold(5,shuffle=True,random_state=20260909).split(X,y)):random_folds[idx]=f
np.save(E/'arrays/random_folds.npy',random_folds)
protocols=[('random',None),('spatial_0',0),('spatial_300',300),('spatial_1000',1000),('random_size_matched_300',None)]
all_metrics=[];fold_metrics=[];audit=[];start=time.time()
with threadpool_limits(limits=2):
    for protocol,buffer in protocols:
        outputs={n:np.zeros(len(y)) for n in ['logistic','boost_raw','boost_calibrated','baseline']}
        for f in range(5):
            rand=buffer is None
            fit,cal,test,outer=split(f,buffer or 0,random_folds if rand else fold,random=rand)
            if protocol=='random_size_matched_300':
                sf,sc,st,so=split(f,300)
                rng=np.random.default_rng(202609091+f)
                for cur,target in [(fit,sf.sum()),(cal,sc.sum()),(outer,so.sum())]:
                    ix=np.flatnonzero(cur);cur[:]=False;cur[rng.choice(ix,min(target,len(ix)),replace=False)]=True
            model_lr=lr();model_lr.fit(X[outer],y[outer]);outputs['logistic'][test]=model_lr.predict_proba(X[test])[:,1]
            if protocol=='spatial_300':
                pair=joblib.load(B/'models'/f'boost_pair_fold_{f}.joblib')
                assert pair['split']['fit_cells']==int(fit.sum()) and pair['split']['calibration_cells']==int(cal.sum())
            else:
                tree=gb();tree.fit(X[fit],y[fit]);c=calibrate(tree,X[cal],y[cal]);pair={'estimator':tree,'calibrator':c}
            outputs['boost_raw'][test]=pair['estimator'].predict_proba(X[test])[:,1]
            outputs['boost_calibrated'][test]=predict(pair,X[test])
            outputs['baseline'][test]=y[outer].mean()
            if protocol=='spatial_300':joblib.dump(pair,E/'models'/f'base_spatial300_fold_{f}.joblib',compress=3)
            testids=np.unique(pid[test]);testids=testids[testids>0]
            shared=np.intersect1d(testids,np.unique(pid[fit]))
            sharedfraction=np.isin(pid[test&(y==1)],shared).mean()
            g=np.zeros(mask.shape,bool);g[rr[fit],cc[fit]]=True
            min_dist=distance_transform_edt(~g)[mask][test]*30
            audit.append({'protocol':protocol,'fold':f,'fit_n':int(fit.sum()),'cal_n':int(cal.sum()),'test_n':int(test.sum()),'lr_fit_n':int(outer.sum()),'fit_positive':int(y[fit].sum()),'test_positive':int(y[test].sum()),'shared_positive_polygons':len(shared),'test_positive_cells_sharing_fit_polygon_fraction':float(sharedfraction),'test_distance_to_fit_median_m':float(np.median(min_dist)),'test_distance_to_fit_min_m':float(min_dist.min())})
            for name,p in outputs.items():fold_metrics.append({'protocol':protocol,'fold':f,'model':name,**metrics(y[test],p[test],cal=name!='baseline')})
            print(json.dumps({'cv_protocol':protocol,'fold':f,'seconds':round(time.time()-start,1)}),flush=True)
        for name,p in outputs.items():all_metrics.append({'protocol':protocol,'model':name,**metrics(y,p,cal=name!='baseline')})
        np.savez_compressed(E/'arrays'/f'cv_{protocol}.npz',**outputs)
        csv('cv_metrics_pooled.csv',pd.DataFrame(all_metrics));csv('cv_metrics_folds.csv',pd.DataFrame(fold_metrics));csv('cv_split_audit.csv',pd.DataFrame(audit))
original=np.load(B/'models/predictions.npz')['boost_calibrated']
assert np.allclose(np.load(E/'arrays/cv_spatial_300.npz')['boost_calibrated'],original,atol=1e-10)
dump('cv_results.json',{'seconds':time.time()-start,'models_fitted_or_reused':50,'probability_scope':'inventory membership; all protocols preserve natural cell prevalence','random_definition':'Stratified random cell assignment; deliberately no polygon protection or spatial guard','spatial_definition':'Same existing 47 geographic blocks and 5 folds; all held-out positive polygons excluded from fitting; guards 0/300/1000m applied at both fit/calibration and fit/test boundaries','extra_control':'Random size matched to spatial+300 fit and calibration counts, same random test folds','calibration_slope_note':'Post-hoc diagnostic logistic regression of observed labels on OOF logits, not recalibration of test predictions','original_spatial300_replay':'PASS'})
