"""Prior-conditional stochastic runout, not calibrated hazard probability. SI units."""
from pathlib import Path
import sys,json,csv,math,time,os
import numpy as np
# Dependencies are installed using requirements-spatial.txt.
import rasterio
from rasterio.features import geometry_mask
from rasterio.windows import Window
from shapely.geometry import shape,mapping
from shapely.ops import transform as gt
from pyproj import Transformer
from scipy.stats import truncnorm
from scipy.ndimage import binary_dilation
BASE=Path(os.environ.get('GANSU_RAW_INPUT', 'external_inputs')).resolve()
OUT=Path(os.environ.get('GANSU_RUNOUT_OUTPUT', 'results/runout')).resolve()
for p in ['results','gis','figures','qa','scripts','source_data']:(OUT/p).mkdir(parents=True,exist_ok=True)
def source(rel):
    local=OUT/'source_data/project_inputs'/rel
    return local if local.exists() else BASE/rel
def rows(p):return list(csv.DictReader(open(p,encoding='utf-8-sig')))
def write(p,rs):
    if not rs:return
    with open(p,'w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rs[0]));w.writeheader();w.writerows(rs)
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
def run():
    G=source('tianshui_probability_20260909/gis')
    with rasterio.open(G/'elevation_m.tif') as d:
        z=d.read(1).astype('float64');tr=d.transform;profile=d.profile;good=(z!=d.nodata)&np.isfinite(z)
    with rasterio.open(G/'study_mask_30m.tif') as d:good&=d.read(1)==1
    nr,nc=z.shape;inv=~tr;zf=z.ravel();gf=good.ravel()
    to=Transformer.from_crs(4326,profile['crs'],always_xy=True).transform
    fs=json.loads((G/'landslides_2013_WGS84.geojson').read_text())['features']
    sources=[];starts=[];areas=[]
    for f in fs:
        geom=gt(to,shape(f['geometry']))
        if not geom.is_valid:geom=geom.buffer(0)
        c0,r0=inv*(geom.bounds[0],geom.bounds[3]);c1,r1=inv*(geom.bounds[2],geom.bounds[1])
        c0=max(0,math.floor(c0));r0=max(0,math.floor(r0));c1=min(nc,math.ceil(c1)+1);r1=min(nr,math.ceil(r1)+1)
        if c1<=c0 or r1<=r0:continue
        m=geometry_mask([mapping(geom)],out_shape=(r1-r0,c1-c0),transform=rasterio.windows.transform(Window(c0,r0,c1-c0,r1-r0),tr),invert=True,all_touched=True)&good[r0:r1,c0:c1]
        rr,cc=np.where(m)
        if not len(rr):continue
        ix=np.argsort(z[rr+r0,cc+c0],kind='stable')[-3:][::-1];cells=(rr[ix]+r0)*nc+cc[ix]+c0
        starts.append(np.pad(cells,(0,3-len(cells)),mode='edge'));areas.append(geom.area)
        sources.append(dict(source_index=len(sources),polygon_id=f['properties']['polygon_id'],candidate_cells=';'.join(map(str,cells)),candidate_count=len(cells),area_m2=geom.area))
    starts=np.array(starts);areas=np.array(areas);ns=len(starts)
    write(OUT/'source_data/随机起点候选.csv',sources)
    # Downhill roulette transition. slope^4 is a declared scenario choice, not a fitted Flow-R model.
    destinations=np.zeros((z.size,8),dtype='int32');weights=np.zeros((z.size,8),dtype='float64');distances=[]
    for j,(dr,dc) in enumerate([(a,b) for a in [-1,0,1] for b in [-1,0,1] if a or b]):
        dist=30*math.hypot(dr,dc);distances.append(dist)
        rr,cc=np.indices(z.shape);rt=rr+dr;ct=cc+dc;valid=(rt>=0)&(rt<nr)&(ct>=0)&(ct<nc)
        rt=rt.clip(0,nr-1);ct=ct.clip(0,nc-1);de=rt*nc+ct
        slope=np.maximum((z-z[rt,ct])/dist,0);valid&=good&good[rt,ct]
        weights[:,j]=np.where(valid,slope**4,0).ravel();destinations[:,j]=de.ravel()
    total=weights.sum(1);cdf=np.cumsum(weights/np.where(total[:,None]>0,total[:,None],1),axis=1).astype('float32');del weights
    distances=np.array(distances)
    bfs=[f for f in json.loads(source('three_regions_20260909/tianshui/osm_exposure_clipped_WGS84.geojson').read_text())['features'] if f['properties']['kind']=='building']
    # uint64 bitset retains overlapping footprints: 44 buildings < 64.
    bits=np.zeros(z.size,dtype='uint64');widebits=bits.copy()
    for b,f in enumerate(bfs):
        m=geometry_mask([mapping(gt(to,shape(f['geometry'])))],out_shape=z.shape,transform=tr,all_touched=True,invert=True)&good
        bits[m.ravel()]|=np.uint64(1)<<np.uint64(b)
        widebits[(binary_dilation(m,structure=np.ones((3,3)))&good).ravel()]|=np.uint64(1)<<np.uint64(b)
    N=5000;rng=np.random.default_rng(20260911);nb=len(bfs)
    hits=np.zeros((ns,nb),dtype='int32');widehits=hits.copy();visits=np.zeros(z.size,dtype='int64');hitrecords=[]
    stopcounts=dict(energy=0,sink=0);stepmax=0;t0=time.time()
    # Stratification over sources avoids an almost empty 5000-path regional sample.
    for lo in range(0,ns,64):
        hi=min(lo+64,ns);si=np.repeat(np.arange(lo,hi),N);n=len(si)
        # For sources with fewer than 3 candidates sample only actual candidates.
        choices=(rng.random(n)*np.array([sources[i]['candidate_count'] for i in si])).astype(int)
        k=starts[si,choices];z0=zf[k];s=np.zeros(n);seen=bits[k].copy();wseen=widebits[k].copy()
        alpha=truncnorm.ppf(rng.random(n),-3,3,loc=25,scale=5);mu=np.tan(np.deg2rad(alpha))
        # Save first arrival speed to prevent repeated residence from biasing f(v|arrival).
        initial_nonzero=np.flatnonzero(seen)
        for b in range(nb):
            q=initial_nonzero[(seen[initial_nonzero]&(np.uint64(1)<<np.uint64(b)))!=0]
            if len(q):hitrecords.append(np.column_stack([si[q],np.full(len(q),b),np.zeros(len(q)),alpha[q]]))
        alive=np.arange(n);np.add.at(visits,k,1);steps=0
        while len(alive):
            kk=k[alive];has=total[kk]>0;stopcounts['sink']+=int((~has).sum());alive=alive[has];kk=kk[has]
            if not len(alive):break
            j=np.minimum((rng.random(len(alive))[:,None]>cdf[kk]).sum(1),7)
            dest=destinations[kk,j];sn=s[alive]+distances[j];energy=z0[alive]-zf[dest]-mu[alive]*sn
            ok=energy>=0;stopcounts['energy']+=int((~ok).sum());alive=alive[ok];dest=dest[ok];energy=energy[ok]
            if not len(alive):break
            s[alive]=sn[ok];k[alive]=dest
            new=bits[dest]&~seen[alive];seen[alive]|=bits[dest];wseen[alive]|=widebits[dest]
            nonzero=np.flatnonzero(new)
            if len(nonzero):
                for b in range(nb):
                    q=nonzero[(new[nonzero]&(np.uint64(1)<<np.uint64(b)))!=0]
                    if len(q):
                        a=alive[q];hitrecords.append(np.column_stack([si[a],np.full(len(q),b),np.sqrt(2*9.81*energy[q]),alpha[a]]))
            np.add.at(visits,dest,1);steps+=1
            if steps>z.size:raise RuntimeError('strict descent violated')
        stepmax=max(stepmax,steps)
        for b in range(nb):
            hits[lo:hi,b]=((seen&(np.uint64(1)<<np.uint64(b)))!=0).reshape(hi-lo,N).sum(1)
            widehits[lo:hi,b]=((wseen&(np.uint64(1)<<np.uint64(b)))!=0).reshape(hi-lo,N).sum(1)
        print(f'sources {hi}/{ns}, paths {hi*N}, elapsed {time.time()-t0:.1f}s',flush=True)
    rec=np.concatenate(hitrecords) if hitrecords else np.empty((0,4))
    assert np.array_equal(np.bincount(rec[:,1].astype(int),minlength=nb),hits.sum(0))
    assert np.all(widehits>=hits)
    p=hits/N;wp=widehits/N;uniform=np.full(ns,1/ns);aw=areas/areas.sum()
    out=[];features=[];sourcepairs=[];curves=[]
    for b,f in enumerate(bfs):
        vals=rec[rec[:,1]==b,2];active=np.flatnonzero(hits[:,b]);means=p[:,b]
        pm=float(uniform@means);se=float(np.sqrt(np.sum(uniform**2*means*(1-means)/N)))
        row=dict(osm_id=f['properties']['osm_id'],evidence_scope='uncalibrated_prior_single_source_event',P_reach_equal_source=pm,MC_SE_plugin=se,P_reach_area_weighted=float(aw@means),P_neighbourhood_equal_source=float(uniform@wp[:,b]),contributing_sources=len(active),arrival_samples=len(vals),v_proxy_p05=float(np.quantile(vals,.05)) if len(vals) else None,v_proxy_p50=float(np.median(vals)) if len(vals) else None,v_proxy_p95=float(np.quantile(vals,.95)) if len(vals) else None,annual_rate=None,whole_building_severe_probability=None,EAL=None)
        out.append(row);features.append(dict(type='Feature',geometry=f['geometry'],properties=row))
        for a in active:
            sourcepairs.append(dict(osm_id=row['osm_id'],polygon_id=sources[a]['polygon_id'],hits=int(hits[a,b]),trials=N,P_reach_given_source=float(means[a])))
        for v in np.linspace(0,50,101):
            curves.append(dict(osm_id=row['osm_id'],v_proxy_m_s=v,P_v_gt_given_arrival=float(np.mean(vals>v)) if len(vals) else None,P_arrival_and_v_gt=float(np.sum(vals>v)/(ns*N)),scope='kinematic_proxy_not_building_pressure'))
    write(OUT/'results/天水44建筑条件到达概率.csv',out)
    write(OUT/'results/源区与建筑条件到达概率.csv',sourcepairs)
    write(OUT/'results/建筑到达速度代理超越曲线.csv',curves)
    dump(OUT/'gis/天水建筑先验条件到达概率.geojson',dict(type='FeatureCollection',features=features))
    np.savez_compressed(OUT/'source_data/传播随机试算原始统计.npz',hits=hits,widehits=widehits,first_arrivals=rec,source_areas=areas,starts=starts,visits=visits,dem=np.where(good,z,np.nan),transform=np.array(tuple(tr)))
    prof=profile.copy();prof.update(dtype='float32',nodata=-9999,compress='deflate',count=1)
    with rasterio.open(OUT/'gis/单源事件等权条件路径概率.tif','w',**prof) as d:d.write(np.where(good,visits.reshape(z.shape)/(ns*N),-9999).astype('float32'),1)
    audit=dict(seed=20260911,sources=ns,trials_per_source=N,total_paths=ns*N,initial_candidates='up to 3 highest cells touching each historical polygon; equal chance',angle_prior='truncated normal mean25 sd5 bounds10..40 deg, uncalibrated',routing='random strictly downhill neighbour with probability proportional to slope^4',energy='v^2/(2g)=z0-z-mu*s, mu=tan(alpha), v0=0',source_weights='equal source or affected-polygon area; neither calibrated occurrence weights',first_arrival_records=len(rec),direct_buildings=sum(x['arrival_samples']>0 for x in out),neighbourhood_buildings=sum(x['P_neighbourhood_equal_source']>0 for x in out),maximum_steps=stepmax,stops=stopcounts,seconds=time.time()-t0,limitations=['No volume, depth, entrainment, momentum runup or barrier model','Fixed 30m DEM: terrain uncertainty not sampled','Velocity is residual energy proxy, not measured or validated building impact speed','No new sources outside historical polygons','Current OSM is incomplete','MC standard error is conditional simulation error, not epistemic confidence','Zero hits does not prove zero physical probability','No source activation rates or whole-building fragility attached'])
    dump(OUT/'qa/传播计算审计.json',audit);print(json.dumps(audit,ensure_ascii=False),flush=True)
if __name__=='__main__':run()
