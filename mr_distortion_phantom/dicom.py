"""Local-only raw DICOM analysis. No identifying headers are exported."""
import json
from pathlib import Path
import numpy as np
import pydicom
from scipy import ndimage as ndi
from .core import lattice,unique_map,align_indices,evaluate

def read_volume(folder):
    ds=[pydicom.dcmread(f) for f in Path(folder).glob('*.dcm')]
    if len(ds)<3:raise ValueError('A classic single-frame DICOM volume is required')
    if len({str(d.SeriesInstanceUID) for d in ds})!=1:raise ValueError('One series per input directory is required')
    iop=np.array(ds[0].ImageOrientationPatient,float);normal=np.cross(iop[:3],iop[3:]);ds.sort(key=lambda d:np.dot(d.ImagePositionPatient,normal))
    pos=np.array([d.ImagePositionPatient for d in ds],float);step=(pos[-1]-pos[0])/(len(ds)-1)
    if not all(np.allclose(d.ImageOrientationPatient,iop,atol=1e-5) for d in ds):raise ValueError('Nonuniform orientation')
    if np.max(np.linalg.norm(pos-(pos[0]+np.arange(len(ds))[:,None]*step),axis=1))>=.002:raise ValueError('Nonuniform slice positions')
    sp=np.array([np.linalg.norm(step),*map(float,ds[0].PixelSpacing)])
    if np.any(sp<=0):raise ValueError('Invalid voxel spacing')
    basis=np.column_stack([step,iop[3:]*sp[1],iop[:3]*sp[2]])
    vox=np.stack([d.pixel_array.astype(np.float32)*float(d.get('RescaleSlope',1))+float(d.get('RescaleIntercept',0)) for d in ds])
    return dict(v=vox,origin=pos[0],basis=basis,spacing=sp)

def detect(data,modality,ct_crop=None,sigma=1.4):
    a=data;v=a['v'];sp=a['spacing'];offset=np.zeros(3)
    if modality=='CT':
        if ct_crop is None:raise ValueError('Set a CT crop covering the cranial grid (z0 z1 y0 y1 x0 x1)')
        starts=np.array(ct_crop[::2]);ends=np.array(ct_crop[1::2])
        if np.any(starts<0) or np.any(ends>v.shape) or np.any(ends<=starts):raise ValueError('Invalid CT crop')
        offset=starts;v=v[starts[0]:ends[0],starts[1]:ends[1],starts[2]:ends[2]];signal=np.clip(v,-50,400)
    else:signal=-v
    smooth=ndi.gaussian_filter(signal,sigma/sp);response=smooth-ndi.gaussian_filter(signal,3.5/sp)
    maxima=response==ndi.maximum_filter(response,size=tuple(int(round(7/s))|1 for s in sp))
    valid=ndi.binary_erosion(np.ones(v.shape,bool),iterations=5)
    if modality!='CT':valid &= ndi.gaussian_filter(v,4/sp)>np.percentile(v,65)*.5
    coords=np.argwhere(maxima&valid&(response>np.percentile(response,98)));scores=response[tuple(coords.T)];refined=coords.astype(float)
    for axis in range(3):
        lo=coords.copy();hi=coords.copy();lo[:,axis]-=1;hi[:,axis]+=1;vl=response[tuple(lo.T)];vh=response[tuple(hi.T)]
        refined[:,axis]+=.5*(vl-vh)/(vl-2*scores+vh)
    return (refined+offset)@a['basis'].T+a['origin']

def analyze(args):
    maps=[]
    for folder,modality in [(args.ct,'CT'),(args.original,'MR'),(args.corrected,'MR')]:
        pts=detect(read_volume(folder),modality,args.ct_crop);pts,ids,_,_=lattice(pts);maps.append(unique_map(ids,pts))
    ct=maps[0];maps=[ct,*[align_indices(ct,m) for m in maps[1:]]];keys=sorted(set.intersection(*(set(m) for m in maps)))
    center=np.median(np.array(list(ct)),axis=0);radius=np.linalg.norm((np.array(keys)-center)*[11,10.5,10],axis=1);fit=radius<=30
    target=np.array([ct[k] for k in keys]);out=Path(args.output);out.mkdir(parents=True,exist_ok=True);result={}
    for state,m in zip(['original','corrected'],maps[1:]):
        measured=np.array([m[k] for k in keys]);stats,v,R,t=evaluate(target,measured,fit);result[state]=stats
        np.savetxt(out/(state+'_residuals.csv'),np.c_[np.array(keys),fit,radius,v,np.linalg.norm(v,axis=1)],delimiter=',',header='i,j,k,registration,radius_mm,L_mm,P_mm,S_mm,error_mm',comments='')
    (out/'metrics.json').write_text(json.dumps(result,indent=2));print('Phantom analysis complete; identifiers were not exported.')
