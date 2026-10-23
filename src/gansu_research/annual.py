"""Independent parametric benchmark; not calibrated annual risk for mapped assets."""
import numpy as np
from scipy.special import ndtr

PARAMS=['lambda_per_year','angle_deg','depth_m','capacity_scale','beta','loss_scale','value_CNY']
LOW=np.array([.02,15,.3,.6,.35,.7,150000]);HIGH=np.array([.2,35,1.5,1.6,.9,1.,450000])
BASELINE=np.array([.1,25,.8,1.,.6,1.,300000])
GQ,GW=np.polynomial.legendre.leggauss(64);U=(GQ+1)/2;W=GW/2
THETA=np.array([100,300,600,1200.]) # hypothetical pressure thresholds kPa
LR=np.array([0,.05,.2,.6,1.])

def evaluate(x,angle_add=0,energy_removed=0,capacity_mult=1,arrival_mult=1,details=False):
    x=np.atleast_2d(x);lam,alpha,h,cscale,beta,lscale,C=x.T
    if not(np.isfinite(x).all() and np.all(lam>=0) and np.all((alpha>0)&(alpha<90)) and np.all(h>=0) and np.all(cscale>0) and np.all(beta>0) and np.all((lscale>=0)&(lscale<=1)) and np.all(C>=0)):
        raise ValueError('invalid physical/probability input')
    if not(0<=arrival_mult<=1 and energy_removed>=0 and capacity_mult>0):raise ValueError('invalid intervention')
    threshold=100*np.tan(np.deg2rad(alpha+angle_add))+energy_removed
    reach=np.clip(1-threshold/80,0,1)*arrival_mult
    energy=np.maximum(80-threshold,0)[:,None]*U
    # Idealized energy proxy on a 100m path; uniform initial drop [0,80]m.
    v2=2*9.81*energy
    pressure=(1800*v2+.5*1800*9.81*h[:,None])/1000 # kPa mean over assumed depth
    q=ndtr((np.log(np.maximum(pressure,1e-300))[:,:,None]-np.log(THETA[None,None,:]*cscale[:,None,None]*capacity_mult))/beta[:,None,None])
    conditional=np.sum(q*W[None,:,None],axis=1)
    eventq=reach[:,None]*conditional
    p=np.column_stack([1-eventq[:,0],eventq[:,:-1]-eventq[:,1:],eventq[:,-1]])
    assert np.min(p)>-1e-12 and np.allclose(p.sum(1),1)
    el=C*lscale*(p@LR);eal=lam*el
    if details:return dict(reach=reach,eventq=eventq,state_prob=p,event_EL=el,EAL=eal,annual_exceed=-np.expm1(-lam[:,None]*eventq))
    return eal
