"""Independent measured-landmark geometry; row-vector LPS coordinates in mm."""
import itertools
import numpy as np
from scipy.spatial import cKDTree
def lattice(points):
 pairs=np.array(list(cKDTree(points).query_pairs(12.4)))
 vec=points[pairs[:,1]]-points[pairs[:,0]]
 vec=vec[np.linalg.norm(vec,axis=1)>9]
 axes=[]
 for j in range(3):
  candidates=vec[np.abs(vec[:,j])/np.linalg.norm(vec,axis=1)>.85]
  candidates*=np.sign(candidates[:,j,None])
  middle=np.median(candidates,axis=0)
  good=candidates[np.linalg.norm(candidates-middle,axis=1)<.65]
  axes.append(np.median(good,axis=0))
 basis=np.array(axes).T
 for _ in range(4):
  q=points@np.linalg.inv(basis).T
  phase=np.angle(np.mean(np.exp(q*2j*np.pi),axis=0))/(2*np.pi)
  ids=np.rint(q-phase).astype(int)
  resid=np.linalg.norm((q-phase-ids)@basis.T,axis=1)
  valid=resid<1.2
  fit=np.linalg.lstsq(np.c_[ids[valid],np.ones(valid.sum())],points[valid],rcond=None)[0]
  basis=fit[:3].T
 ids=ids[valid];pts=points[valid]
 # Translate the arbitrary integer origin to the centre of the detected extent.
 centre=np.rint((np.percentile(ids,1,axis=0)+np.percentile(ids,99,axis=0))/2).astype(int)
 ids-=centre
 return pts,ids,basis,valid

def rigid(source,target):
 a=source.mean(0);b=target.mean(0)
 u,_,vt=np.linalg.svd((source-a).T@(target-b))
 sign=np.eye(3);sign[-1,-1]=np.linalg.det(u@vt)
 R=u@sign@vt;t=b-a@R
 return R,t

def unique_map(ids,pts):
 counts={tuple(i):sum(np.all(ids==i,axis=1)) for i in ids}
 return {tuple(i):p for i,p in zip(ids,pts) if counts[tuple(i)]==1}

def metrics(errors):
    e=np.asarray(errors,float)
    if e.ndim!=1 or not len(e) or not np.all(np.isfinite(e)):raise ValueError('Finite nonempty error vector required')
    return {'n':len(e),'mean_mm':float(e.mean()),'median_mm':float(np.median(e)),
            **{f'p{q}_mm':float(np.percentile(e,q)) for q in [5,25,75,95]},'max_mm':float(e.max())}

def align_indices(ct,mr):
    candidates=[]
    for shift in itertools.product(range(-2,3),repeat=3):
        moved={tuple(np.array(k)+shift):v for k,v in mr.items()};keys=sorted(ct.keys()&moved.keys())
        if len(keys)<30:continue
        a=np.array([moved[k] for k in keys]);b=np.array([ct[k] for k in keys]);R,t=rigid(a,b)
        rms=np.sqrt(np.mean(np.sum((a@R+t-b)**2,axis=1)));candidates.append((len(keys),rms,moved))
    if not candidates:raise ValueError('Insufficient common phantom landmarks')
    candidates.sort(key=lambda x:(-x[0],x[1]));n,rms,result=candidates[0]
    if n<650 or rms>=1.5:raise ValueError('Landmark correspondence outside validated phantom range')
    return result

def evaluate(ct,mr,fit):
    ct=np.asarray(ct,float);mr=np.asarray(mr,float);fit=np.asarray(fit,bool)
    if ct.shape!=mr.shape or ct.ndim!=2 or ct.shape[1]!=3 or fit.shape!=(len(ct),):raise ValueError('Coordinate shapes disagree')
    if fit.sum()<30 or (~fit).sum()<1:raise ValueError('Separate fit and evaluation subsets required')
    R,t=rigid(mr[fit],ct[fit]);vectors=mr@R+t-ct;errors=np.linalg.norm(vectors,axis=1)
    return metrics(errors[~fit]),vectors,R,t

def sphere_coverage(diameter_mm,shift_mm):
    D=np.asarray(diameter_mm,float);d=np.asarray(shift_mm,float)
    if np.any(D<=0) or np.any(d<0):raise ValueError('Positive diameter and nonnegative shift required')
    x=d/(D/2)
    return np.where(x>=2,0.,1-3*x/4+x**3/16)
