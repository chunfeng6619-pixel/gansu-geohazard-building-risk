from common import *
from scipy.spatial import cKDTree
import time
start=time.time();base=np.load(B/'models/predictions.npz')['boost_calibrated'];pairs=[joblib.load(B/'models'/f'boost_pair_fold_{f}.joblib') for f in range(5)]
groups={features[j]:[j] for j in range(6)};groups['aspect_joint']=[2,3]
ablations=[];foldrows=[]
with threadpool_limits(limits=2):
    for name,drop in groups.items():
        cols=[j for j in range(6) if j not in drop];p=np.zeros(len(y))
        for f in range(5):
            fit,cal,test,_=split(f,300);model=gb();model.fit(X[fit][:,cols],y[fit]);c=calibrate(model,X[cal][:,cols],y[cal]);pair={'estimator':model,'calibrator':c}
            p[test]=predict(pair,X[test][:,cols]);bm=metrics(y[test],base[test],False);am=metrics(y[test],p[test],False)
            foldrows.append({'removed':name,'fold':f,'delta_auc':bm['auc']-am['auc'],'delta_ap':bm['ap']-am['ap'],'delta_brier':am['brier']-bm['brier'],'delta_log_loss':am['log_loss']-bm['log_loss']})
        bm=metrics(y,base,False);am=metrics(y,p,False)
        ablations.append({'removed':name,**am,'delta_auc':bm['auc']-am['auc'],'delta_ap':bm['ap']-am['ap'],'delta_brier':am['brier']-bm['brier'],'delta_log_loss':am['log_loss']-bm['log_loss']})
        np.save(E/'arrays'/f'lofo_{name}.npy',p.astype('float32'))
        csv('factor_ablation_pooled.csv',pd.DataFrame(ablations));csv('factor_ablation_folds.csv',pd.DataFrame(foldrows))
        print(json.dumps({'ablation_done':name,'delta_ap':ablations[-1]['delta_ap'],'seconds':round(time.time()-start,1)}),flush=True)
    # Exchange entire held-out geographic templates; only positions with a
    # donor at <=45 m local-coordinate distance enter each paired comparison.
    pgroups={'elevation_m':[0],'slope_deg':[1],'aspect_joint':[2,3],'relief_330m_m':[4],'tpi_330m_m':[5]}
    permrows=[]
    for f in range(5):
        idx=np.flatnonzero(fold==f);local=np.column_stack([xx[idx]%4000,yy[idx]%4000]);bb=block[idx];ub=np.unique(bb)
        trees={b:cKDTree(local[bb==b]) for b in ub};members={b:np.flatnonzero(bb==b) for b in ub}
        rng=np.random.default_rng(8900+f)
        for rep in range(20):
            order=rng.permutation(ub);donors=dict(zip(order,np.roll(order,1)));mapping_idx=np.empty(len(idx),int);distances=np.empty(len(idx))
            for target in ub:
                pos=members[target];dist,near=trees[donors[target]].query(local[pos]);mapping_idx[pos]=members[donors[target]][near];distances[pos]=dist
            use=distances<=45;original=X[idx[use]];donated=X[idx[mapping_idx[use]]];labels=y[idx[use]];bm=metrics(labels,base[idx[use]],False)
            for name,cols in pgroups.items():
                z=original.copy();z[:,cols]=donated[:,cols]
                if name=='aspect_joint':assert np.max(np.abs(z[:,2]**2+z[:,3]**2-1))<1e-5
                pp=predict(pairs[f],z);met=metrics(labels,pp,False)
                permrows.append({'factor':name,'fold':f,'replicate':rep+1,'paired_cells':int(use.sum()),'paired_fraction_of_test':float(use.mean()),'donor_distance_median_m':float(np.median(distances[use])),'delta_auc':bm['auc']-met['auc'],'delta_ap':bm['ap']-met['ap'],'delta_brier':met['brier']-bm['brier'],'delta_log_loss':met['log_loss']-bm['log_loss']})
        csv('spatial_permutation_repeats.csv',pd.DataFrame(permrows));print(json.dumps({'permutation_fold_done':f,'seconds':round(time.time()-start,1)}),flush=True)
    perm=pd.DataFrame(permrows);csv('spatial_permutation_summary.csv',perm.groupby('factor').agg(mean_delta_auc=('delta_auc','mean'),mean_delta_ap=('delta_ap','mean'),mean_delta_brier=('delta_brier','mean'),paired_fraction_mean=('paired_fraction_of_test','mean'),replicates=('delta_auc','size')).reset_index())
    # ALE main effects for non-angular features only; no impossible cosine/sine interventions.
    alerows=[]
    for j in [0,1,4,5]:
        edges=np.unique(np.quantile(X[:,j],np.linspace(0,1,21)));bi=np.clip(np.searchsorted(edges,X[:,j],side='right')-1,0,len(edges)-2)
        for f in range(5):
            diffs=[];counts=[]
            for k in range(len(edges)-1):
                ix=np.flatnonzero((fold==f)&(bi==k));z=X[ix].copy();z[:,j]=edges[k];lo=predict(pairs[f],z);z[:,j]=edges[k+1];hi=predict(pairs[f],z)
                diffs.append(float(np.mean(hi-lo)));counts.append(len(ix))
            # Bin midpoint ALE = sum previous increments + half current increment.
            effect=np.cumsum(diffs)-.5*np.array(diffs);effect-=np.average(effect,weights=counts)
            for k,value in enumerate(effect):alerows.append({'feature':features[j],'fold':f,'bin':k,'x_low':edges[k],'x_high':edges[k+1],'x_mid':(edges[k]+edges[k+1])/2,'cells':counts[k],'ale_probability_difference':value})
    csv('ale_main_effects_by_fold.csv',pd.DataFrame(alerows))
    # Supported bivariate conditional responses + local mixed differences.
    # Deliberately NOT labelled pure 2D ALE or causal mechanisms.
    theta=np.mod(np.degrees(np.arctan2(X[:,3],X[:,2])),360)
    descriptors={'slope_deg':(X[:,1],lambda z,v:z.__setitem__((slice(None),1),v)),'elevation_m':(X[:,0],lambda z,v:z.__setitem__((slice(None),0),v)),'relief_330m_m':(X[:,4],lambda z,v:z.__setitem__((slice(None),4),v)),'tpi_330m_m':(X[:,5],lambda z,v:z.__setitem__((slice(None),5),v)),'aspect_degrees':(theta,None)}
    def setval(z,name,val):
        if name=='aspect_degrees':z[:,2]=np.cos(np.radians(val));z[:,3]=np.sin(np.radians(val))
        else:descriptors[name][1](z,val)
    inter=[]
    for first,second in [('slope_deg','aspect_degrees'),('slope_deg','relief_330m_m'),('elevation_m','tpi_330m_m')]:
        v1=descriptors[first][0];v2=descriptors[second][0]
        e1=np.unique(np.quantile(v1,np.linspace(0,1,9)));e2=np.arange(0,361,45) if second=='aspect_degrees' else np.unique(np.quantile(v2,np.linspace(0,1,9)))
        b1=np.clip(np.searchsorted(e1,v1,side='right')-1,0,len(e1)-2);b2=np.clip(np.searchsorted(e2,v2,side='right')-1,0,len(e2)-2)
        for k in range(len(e1)-1):
            for l in range(len(e2)-1):
                sel=(b1==k)&(b2==l);n=int(sel.sum());nb=len(np.unique(block[sel]));mix=[];weights=[]
                if n>=100 and nb>=3:
                    for f in range(5):
                        ix=np.flatnonzero(sel&(fold==f))
                        if not len(ix):continue
                        z=X[ix].copy();setval(z,first,e1[k]);setval(z,second,e2[l]);p00=predict(pairs[f],z)
                        setval(z,first,e1[k+1]);p10=predict(pairs[f],z)
                        setval(z,second,e2[l+1]);p11=predict(pairs[f],z)
                        setval(z,first,e1[k]);p01=predict(pairs[f],z)
                        mix.append(np.mean(p11-p10-p01+p00));weights.append(len(ix))
                inter.append({'pair':first+'__'+second,'first':first,'second':second,'i':k,'j':l,'x_low':e1[k],'x_high':e1[k+1],'y_low':e2[l],'y_high':e2[l+1],'cells':n,'blocks':nb,'supported':n>=100 and nb>=3,'oof_mean_probability':float(base[sel].mean()) if n else np.nan,'observed_fraction':float(y[sel].mean()) if n else np.nan,'local_mixed_probability_difference':float(np.average(mix,weights=weights)) if mix else np.nan})
        csv('bivariate_supported_responses.csv',pd.DataFrame(inter));print(json.dumps({'interaction_done':first+'__'+second,'seconds':round(time.time()-start,1)}),flush=True)
    coef=pd.read_csv(B/'tables/logistic_coefficients_by_fold.csv');csv('coefficient_stability.csv',coef.groupby('feature').agg(mean=('log_odds_per_training_sd','mean'),minimum=('log_odds_per_training_sd','min'),maximum=('log_odds_per_training_sd','max'),fold_sd=('log_odds_per_training_sd','std')).reset_index())
    csv('feature_spearman_correlation.csv',pd.DataFrame(X,columns=features).corr(method='spearman').reset_index(names='feature'))
    dump('factor_results.json',{'ablation_fits':35,'spatial_permutation_predicts':500,'permutation_definition':'Deranged donor 4km blocks within each held-out fold, matching local template positions within 45m; paired baseline is evaluated on exactly the same retained cells','permutation_caution':'Approximate preservation of within-block texture; not exact marginal-preserving or conditional permutation. Correlated terrain predictors can still form unusual combinations. Interpret alongside ablation and ALE, not causally.','aspect_handling':'Cosine and sine moved together in permutation and angular mixed differences; individual component ablations also reported, plus joint aspect ablation','ale_definition':'20 empirical quantile bins for elevation/slope/relief/TPI, accumulated held-out local endpoint differences, weighted centering; 5 fold spread is not confidence interval','two_dimensional_definition':'Conditional OOF response surface plus finite mixed endpoint differences per supported 8x8 bin; not centred second-order ALE, not mechanism identification','support_rule':'At least 100 cells and 3 geographic blocks for bivariate surfaces','seconds':time.time()-start})
print('FACTOR ANALYSIS FINISHED',flush=True)
