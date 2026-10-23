"""Full-municipality public geometry and explicitly hypothetical damage scenarios."""
from pathlib import Path
import os,sys,json,csv,gzip,hashlib,argparse,gc
GIS=os.environ.get('GIS_EXTRA_PATH')
ARROW=os.environ.get('OVERTURE_EXTRA_PATH')
for p in [ARROW,GIS]:
    if p:sys.path.insert(0,p)
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import shapely
from shapely.geometry import shape
from shapely.ops import unary_union,transform
from pyproj import Transformer
import rasterio
from rasterio.transform import from_origin

ROOT=Path(os.environ.get('GANSU_MUNICIPAL_ROOT', 'external_inputs/municipal')).resolve()
OLD=ROOT/'previous_assignments'
PR=Transformer.from_crs(4326,32648,always_xy=True)
HAZARDS=['SLIDE','DEBRIS','ROCK'];CELL=100.

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()

def old_classes(city):
    snapshot=ROOT/'data/retained_previous_type_assignments.csv'
    if snapshot.exists():
        with snapshot.open(encoding='utf-8-sig') as f:
            return {r['gers_id']:int(r['type_code']) for r in csv.DictReader(f) if r['case']==city}
    names={'Lanzhou':['Lanzhou'],'Tianshui':['Tianshui','Tianshui_urban'],'Longnan':['Longnan']}[city]
    result={}
    for key in names:
        f=OLD/(key+'_WGS84.geojson')
        if not f.exists():continue
        for row in json.loads(f.read_text(encoding='utf-8'))['features']:
            result[row['id']]=row['properties']['hypothetical_type_code']
    return result

def assign_types(ids,xy,seed,old):
    rng=np.random.default_rng(seed)
    directions=rng.normal(size=(24,2));directions/=np.linalg.norm(directions,axis=1)[:,None]
    wave=rng.uniform(1000.,5000.,24);phase=rng.uniform(0,2*np.pi,24);amp=rng.normal(size=24)
    result=np.empty(len(ids),dtype=np.uint8);kept=np.zeros(len(ids),dtype=bool)
    for start in range(0,len(ids),50000):
        stop=min(start+50000,len(ids));coords=xy[start:stop]-np.array([500000.,3800000.])
        field=np.sum(np.sin(coords@directions.T/wave*2*np.pi+phase)*amp,axis=1)/np.sqrt(24)
        # ID-stable illustrative local mixing; not a material or hazard parameter.
        uu=np.array([int.from_bytes(hashlib.blake2b((str(seed)+x).encode(),digest_size=8).digest(),'little')/2**64 for x in ids[start:stop]])
        field+=.42*np.sqrt(12)*(uu-.5)
        code=np.where(field<-.15,1,np.where(field<.80,2,3)).astype(np.uint8)
        for j,x in enumerate(ids[start:stop]):
            if x in old:code[j]=old[x];kept[start+j]=True
        result[start:stop]=code
    return result,kept

def prepare(case):
    definitions=json.loads((ROOT/'data/region_definition.json').read_text(encoding='utf-8'))
    row=next(x for x in definitions['regions'] if x['case']==case);city=row['city']
    boundaries=json.loads((ROOT/'data/city_counties_WGS84.geojson').read_text(encoding='utf-8'))['features']
    counties=[f for f in boundaries if f['properties']['region']==city]
    region=unary_union([shape(f['geometry']) for f in counties]);region_utm=transform(PR.transform,region)
    shapely.prepare(region)
    with (ROOT/'data/archetype_probability_intervals.csv').open(encoding='utf-8-sig') as f:
        intervals={(r['prototype'],r['hazard']):r for r in csv.DictReader(f)}
    source=ROOT/'raw'/(case+'_municipality_bbox.parquet');pf=pq.ParquetFile(source)
    print(case,'reading',pf.metadata.num_rows,'bbox footprints',flush=True)
    table=pf.read()
    all_ids=table['id'].to_pylist();wkb=np.asarray(table['geometry'].to_pylist(),dtype=object)
    geoms=shapely.from_wkb(wkb);centres=shapely.centroid(geoms)
    inside=shapely.covers(region,centres);table=table.filter(pa.array(inside))
    geoms=geoms[inside];centres=centres[inside];wkb=wkb[inside]
    ids=table['id'].to_pylist();assert len(ids)==len(set(ids))
    assert all(isinstance(x,str) and x for x in ids)
    xx,yy=PR.transform(shapely.get_x(centres),shapely.get_y(centres));xy=np.column_stack([xx,yy])
    assert np.isfinite(xy).all()
    n=len(ids);seed={'Lanzhou':2026102001,'Tianshui':2026102002,'Longnan':2026102003}[case]
    old=old_classes(case);codes,kept=assign_types(ids,xy,seed,old)
    prefix={'Lanzhou':'LZ','Tianshui':'TS','Longnan':'LN'}[case]
    county_index=np.full(n,-1,dtype=np.int16);county_counts=[]
    for j,f in enumerate(counties):
        cg=shape(f['geometry']);shapely.prepare(cg)
        mask=(county_index<0)&shapely.covers(cg,centres)
        county_index[mask]=j;county_counts.append({'county':f['properties']['county'],'count':int(mask.sum())})
    assert (county_index>=0).all() and sum(x['count'] for x in county_counts)==n
    invalid=int((~shapely.is_valid(geoms)).sum())
    values={
        'municipality':pa.array([city]*n),
        'county':pa.array([counties[j]['properties']['county'] for j in county_index]),
        'actual_structure_type':pa.array(['未知']*n),
        'scenario_type_code':pa.array(codes),
        'scenario_prototype':pa.array([prefix+str(c) for c in codes]),
        'previous_scenario_assignment_retained':pa.array(kept),
        'centroid_easting_UTM48N':pa.array(xx),'centroid_northing_UTM48N':pa.array(yy),
        'hazard_occurrence_probability':pa.nulls(n,type=pa.float64()),
        'annual_EAL':pa.nulls(n,type=pa.float64())}
    a,b,c,d=region_utm.bounds;ox=np.floor(a/CELL)*CELL;oy=np.floor(b/CELL)*CELL
    nx=int(np.ceil((c-ox)/CELL));ny=int(np.ceil((d-oy)/CELL))
    cols=np.floor((xx-ox)/CELL).astype(np.int32);rows=np.floor((yy-oy)/CELL).astype(np.int32)
    assert (cols>=0).all() and (cols<nx).all() and (rows>=0).all() and (rows<ny).all()
    flat=rows.astype(np.int64)*nx+cols
    counts=np.bincount(flat*3+codes-1,minlength=nx*ny*3).reshape(ny*nx,3).astype(np.uint32)
    total=counts.sum(axis=1);occupied=np.flatnonzero(total);assert total.sum()==n
    arrays={'count':total.reshape(ny,nx).astype(np.uint32)};status_counts={}
    metadata={'case':case,'city':city,'scope':'full municipality; all constituent county polygons',
        'source_bbox':row['bbox'],'region_bounds_UTM48N':list(region_utm.bounds),
        'origin_lower_left_UTM48N':[float(ox),float(oy)],'grid_cell_m':CELL,'grid_shape':[ny,nx],
        'n_buildings':n,'n_occupied_cells':len(occupied),'n_counties':len(counties),'county_counts':county_counts,
        'previous_assignment_matches':int(kept.sum()),'source_invalid_footprints':invalid,
        'mean_weight':'equal count of public footprint records',
        'NC_rule':'any NC-assigned footprint makes the occupied grid cell NC; no numerical mean is drawn',
        'unoccupied_rule':'NoData, not zero damage or absence of buildings',
        'hypothetical_assignment_seed':seed,
        'hypothetical_assignment':'spatial Fourier field, 24 waves 1000–5000 m, ID-stable local mixing; cutoffs -0.15 and 0.80; existing IDs retain prior type',
        'type_proportions_are_observed':False,'raw_bbox_rows':len(all_ids),'release':'2026-09-23.1',
        'access_date':'2026-10-20','source_geometry_modified':False,
        'source_sha256':digest(source),'source_url':'https://docs.overturemaps.org/guides/buildings/',
        'conditional_metric':'mean of existing local P(Dt>=0.5 | prototype, prescribed load, beta=.5) bounds',
        'regional_hazard_occurrence_model':False,'new_FE_runs':0,'threshold_Dt':.5,'beta_assumed':.5}
    for hazard in HAZARDS:
        rr=[intervals[(prefix+str(j),hazard)] for j in [1,2,3]]
        lower=np.array([float(r['probability_lower']) for r in rr]);upper=np.array([float(r['probability_upper']) for r in rr])
        bad=np.array([r['FE_status']!='eligible' for r in rr]);point_nc=bad[codes-1]
        values[hazard+'_probability_lower']=pa.array(lower[codes-1])
        values[hazard+'_probability_upper']=pa.array(upper[codes-1])
        values[hazard+'_map_probability_upper']=pa.array(upper[codes-1],mask=point_nc)
        values[hazard+'_NC']=pa.array(point_nc)
        nc_counts=counts[:,bad].sum(axis=1) if bad.any() else np.zeros(len(total),dtype=np.uint32)
        cell_nc=(total>0)&(nc_counts>0);eligible=(total>0)&~cell_nc
        lo=np.full(len(total),np.nan);hi=np.full(len(total),np.nan)
        lo[eligible]=counts[eligible]@lower/total[eligible];hi[eligible]=counts[eligible]@upper/total[eligible]
        assert np.all((lo[eligible]>=0)&(hi[eligible]<=1)&(lo[eligible]<=hi[eligible]))
        arrays[hazard+'_lower']=lo.reshape(ny,nx);arrays[hazard+'_upper']=hi.reshape(ny,nx)
        arrays[hazard+'_NC']=cell_nc.reshape(ny,nx)
        status_counts[hazard]={'eligible_cells':int(eligible.sum()),'NC_cells':int(cell_nc.sum()),
                               'eligible_buildings':int((~point_nc).sum()),'NC_buildings':int(point_nc.sum())}
        output=ROOT/'data'/(case+'_'+hazard+'_100m_probability_bounds.tif')
        affine=from_origin(ox,oy+ny*CELL,CELL,CELL)
        with rasterio.open(output,'w',driver='GTiff',width=nx,height=ny,count=4,
            dtype='float32',crs='EPSG:32648',transform=affine,nodata=-9999.,compress='deflate',tiled=True) as dst:
            bands=[lo.reshape(ny,nx),hi.reshape(ny,nx),
                   np.where(total.reshape(ny,nx)>0,cell_nc.reshape(ny,nx).astype(float),np.nan),
                   np.where(total.reshape(ny,nx)>0,total.reshape(ny,nx).astype(float),np.nan)]
            for j,band in enumerate(bands,1):
                dst.write(np.flipud(np.where(np.isfinite(band),band,-9999.)).astype(np.float32),j)
            for j,name in enumerate(['mean_conditional_p_lower','mean_conditional_p_upper','NC_cell_flag','footprint_count'],1):dst.set_band_description(j,name)
            dst.update_tags(scope='hypothetical archetype allocations; prescribed loads; conditional local damage',
                            cell_size_m='100',threshold_Dt='0.5',beta_assumed='0.5',NC_rule='any NC -> no probability colour')
    for name,array in values.items():table=table.append_column(name,array)
    meta=dict(table.schema.metadata or {});meta[b'scenario_scope']=json.dumps(metadata,ensure_ascii=False).encode('utf-8')
    table=table.replace_schema_metadata(meta)
    pq.write_table(table,ROOT/'data'/(case+'_municipality_building_scenarios.parquet'),compression='zstd')
    metadata['status_counts']=status_counts
    np.savez_compressed(ROOT/'data'/(case+'_100m_grid.npz'),**arrays)
    with gzip.open(ROOT/'data'/(case+'_100m_occupied_cells.csv.gz'),'wt',encoding='utf-8-sig',newline='') as f:
        header=['grid_id','easting_UTM48N','northing_UTM48N','n_footprints','n_type1','n_type2','n_type3']
        for h in HAZARDS:header.extend([h+'_mean_p_lower',h+'_mean_p_upper',h+'_NC_cell'])
        writer=csv.writer(f);writer.writerow(header)
        for k in occupied:
            rr=int(k//nx);cc=int(k%nx)
            rec=[int(k),float(ox+(cc+.5)*CELL),float(oy+(rr+.5)*CELL),int(total[k]),*map(int,counts[k])]
            for h in HAZARDS:
                lo=arrays[h+'_lower'][rr,cc];hi=arrays[h+'_upper'][rr,cc]
                rec.extend(['' if np.isnan(lo) else float(lo),'' if np.isnan(hi) else float(hi),int(arrays[h+'_NC'][rr,cc])])
            writer.writerow(rec)
    (ROOT/'data'/(case+'_metadata.json')).write_text(json.dumps(metadata,ensure_ascii=False,indent=2),encoding='utf-8')
    print(case,n,'municipality footprints;',len(occupied),'occupied 100 m cells;',int(kept.sum()),'retained prior types',flush=True)
    del table,geoms,centres,wkb,arrays,counts;gc.collect()

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--case',choices=['Lanzhou','Tianshui','Longnan','all'],default='all');args=p.parse_args()
    for case in (['Lanzhou','Tianshui','Longnan'] if args.case=='all' else [args.case]):prepare(case)
