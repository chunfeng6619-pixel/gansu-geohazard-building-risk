from pathlib import Path
import os,sys,json
B=Path(os.environ.get('GANSU_SPATIAL_INPUT', 'external_inputs/spatial')).resolve()
E=Path(os.environ.get('GANSU_SPATIAL_OUTPUT', 'results/spatial')).resolve()
# Dependencies are installed using requirements-spatial.txt.
os.environ['OMP_NUM_THREADS']='2';os.environ['OPENBLAS_NUM_THREADS']='2'
import numpy as np,pandas as pd
from scipy.ndimage import distance_transform_edt
from scipy.special import expit,logit
from scipy.optimize import minimize,brentq
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score,average_precision_score,brier_score_loss,log_loss
from threadpoolctl import threadpool_limits
import joblib,rasterio
from rasterio.transform import Affine
for sub in ['tables','arrays','models','gis','figures','qa']:(E/sub).mkdir(parents=True,exist_ok=True)
d=np.load(B/'models/grid_arrays.npz')
X=d['X']; y=d['y'];pid=d['pid'];fold=d['fold'];block=d['block'];rr=d['rr'];cc=d['cc'];mask=d['mask'];features=d['features'].tolist()
tf=Affine(*d['transform'][:6]);h,w=mask.shape;blocks=np.unique(block)
xx=tf.c+(cc+.5)*tf.a; yy=tf.f+(rr+.5)*tf.e
def dump(name,obj):
    (E/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=lambda v:v.item() if hasattr(v,'item') else str(v)),encoding='utf-8')
def csv(name,df):df.to_csv(E/'tables'/name,index=False,encoding='utf-8-sig')
def gb():return HistGradientBoostingClassifier(learning_rate=.08,max_iter=100,max_leaf_nodes=15,min_samples_leaf=200,l2_regularization=10,early_stopping=False,random_state=20260909)
def lr():return make_pipeline(StandardScaler(),LogisticRegression(C=1,max_iter=400))
def eligible(hold,buffer=300):
    g=np.zeros(mask.shape,bool);g[rr[hold],cc[hold]]=True
    if buffer>0:keep=distance_transform_edt(~g)[mask]*30>buffer
    else:keep=~hold
    shared=np.unique(pid[hold]);shared=shared[shared>0]
    return keep&~np.isin(pid,shared)
def split(f,buffer=300,folds=fold,random=False):
    test=folds==f
    outer=(~test) if random else eligible(test,buffer)
    cal=(folds==(f+1)%5)&outer
    fit=outer&(folds!=(f+1)%5)
    if not random:fit&=eligible(cal,buffer)
    return fit,cal,test,outer
def calibrate(model,z,label,weights=None):
    c=LogisticRegression(C=1e6,max_iter=300)
    c.fit(model.decision_function(z).reshape(-1,1),label,sample_weight=weights)
    return c
def predict(pair,z):return pair['calibrator'].predict_proba(pair['estimator'].decision_function(z).reshape(-1,1))[:,1]
def calibration_fit(label,p,weights=None):
    z=logit(np.clip(p,1e-10,1-1e-10));wt=np.ones(len(p)) if weights is None else weights;wt=wt/wt.sum()
    def fun(ab):
        t=ab[0]+ab[1]*z;res=expit(t)-label
        return np.sum(wt*(np.logaddexp(0,t)-label*t)),np.array([np.sum(wt*res),np.sum(wt*res*z)])
    opt=minimize(fun,[0,1],jac=True,method='BFGS',options={'gtol':1e-9,'maxiter':100})
    citl=brentq(lambda a:np.sum(wt*(expit(a+z)-label)),-40,40)
    return float(opt.x[0]),float(opt.x[1]),float(citl)
def metrics(label,p,cal=True):
    out={'n':len(label),'positive':int(label.sum()),'prevalence':float(label.mean()),'mean_probability':float(np.mean(p)),'auc':roc_auc_score(label,p),'ap':average_precision_score(label,p),'brier':brier_score_loss(label,p),'log_loss':log_loss(label,p)}
    if cal:
        out['calibration_intercept'],out['calibration_slope'],out['calibration_in_large']=calibration_fit(label,p)
    return out
def raster(name,v,fill=-9999,dtype='float32'):
    a=np.full(mask.shape,fill,dtype=dtype);a[mask]=v
    with rasterio.open(E/'gis'/f'{name}.tif','w',driver='GTiff',height=h,width=w,count=1,dtype=dtype,crs='EPSG:32648',transform=tf,nodata=fill,compress='deflate',tiled=True) as ds:
        ds.write(a,1);ds.update_tags(scope='2013 inventory membership, NOT annual hazard probability',method=name)
