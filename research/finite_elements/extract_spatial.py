from __future__ import print_function
import os,json,math,sys,traceback
import numpy as np
from odbAccess import openOdb
from abaqusConstants import MAX_PRINCIPAL
R=os.environ.get('GANSU_FE_NEW', 'external_inputs/fe_new');OLD=os.environ.get('GANSU_FE_OLD', 'external_inputs/fe_old')
MODE=sys.argv[1] if len(sys.argv)>1 else 'old'
OUT=R+'/spatial_'+MODE
if not os.path.isdir(OUT):os.makedirs(OUT)
manifest=json.load(open(R+('/spatial_manifest.json' if MODE=='old' else '/manifest.json')))
report=[]
def items(field,kind):
    for block in field.bulkDataBlocks:
        data=np.array(block.data,copy=True).tolist();labels=tuple(block.nodeLabels if kind=='node' else block.elementLabels)
        for n,v in zip(labels,data):yield int(n),v
def nodal(fo,key):return {n:v for n,v in items(fo[key],'node')} if key in fo.keys() else {}
def scalar(fo,key,allowed,invariant=None):
    d={}
    if key not in fo.keys():return d
    field=fo[key].getScalarField(invariant=invariant) if invariant is not None else fo[key]
    for e,v in items(field,'element'):
        if e in allowed:d[e]=max(d.get(e,-1e100),float(v[0]))
    return d
for entry in manifest:
    case=str(entry['case'] if MODE=='old' else entry['name']);base=OLD if MODE=='old' else R;odb=None
    try:
        path=base+'/'+case+'.odb';target=OUT+'/'+case+'.json'
        if not os.path.exists(path) or os.path.exists(base+'/'+case+'.lck'):continue
        if os.path.exists(target) and os.path.getmtime(target)>os.path.getmtime(path):
            print(case,'cached');continue
        mesh=json.load(open(base+'/'+case+'_mesh.json'))
        oldcache=json.load(open(base+'/'+case+'_result.json')) if MODE=='old' else None
        meta=oldcache['metadata'] if oldcache else json.load(open(base+'/'+case+'_meta.json'))
        steps=[str(s) for s in entry['steps']] if MODE=='old' else ['DEBRIS']
        odb=openOdb(path=path,readOnly=True);wall=set(mesh['groups']['WALL']);areas={int(k):v for k,v in mesh['areas'].items()};total=sum(areas[e] for e in wall)
        elements={int(k):v for k,v in mesh['elements'].items()};coords={int(k):v for k,v in mesh['nodes'].items()}
        wn=set(n for e in wall for n in elements[e]);beam=set(e for g in ['BEAM','COL','RAFTER','PURLIN','JOIST'] for e in mesh['groups'].get(g,[]));base_nodes=set(mesh['base'])
        model=meta['model'];H=model['H'];nlev=model['n'];levels=[mesh['levels'][str(k)] for k in range(nlev+1)]
        levelwall=[[n for n in wn if abs(coords[n][2]-k*H)<1e-6] for k in range(nlev+1)]
        pairs=[]
        for k in range(nlev):
            bot={tuple(round(v,6) for v in coords[n][:2]):n for n in levelwall[k]}
            pairs.append([(bot[tuple(round(v,6) for v in coords[n][:2])],n) for n in levelwall[k+1] if tuple(round(v,6) for v in coords[n][:2]) in bot])
        face={e:('front' if all(abs(coords[n][1])<1e-6 for n in elements[e]) else 'back' if all(abs(coords[n][1]-model['W'])<1e-6 for n in elements[e]) else 'left' if all(abs(coords[n][0])<1e-6 for n in elements[e]) else 'right' if all(abs(coords[n][0]-model['L'])<1e-6 for n in elements[e]) else 'interior') for e in wall}
        groups={g:[e for e in wall if face[e]==g] for g in ['front','back','left','right','interior']}
        grav=odb.steps['GRAVITY'].frames[-1];gu=nodal(grav.fieldOutputs,'U');gravity={}
        def energy_at(sn,t):
            result={}
            for reg in odb.steps[sn].historyRegions.values():
                for key in ['ALLIE','ALLSE','ALLKE','ALLWK','ALLSD','ALLVD','ALLAE']:
                    if key in reg.historyOutputs.keys():
                        values=reg.historyOutputs[key].data
                        pair=min(values,key=lambda p:abs(p[0]-t))
                        assert abs(pair[0]-t)<2e-5,(case,sn,key,pair[0],t)
                        result[key]=float(pair[1])
            return result
        gravity=energy_at('GRAVITY',grav.frameValue)
        hs=[];onset={};frame_refs=[];seen_fields={}
        for sn in steps:
            for fi,fr in enumerate(odb.steps[sn].frames):
                fo=fr.fieldOutputs;seen_fields[sn]=list(fo.keys());u=nodal(fo,'U');rf=nodal(fo,'RF');aa=nodal(fo,'A');vv=nodal(fo,'V')
                du={n:[u[n][j]-gu[n][j] for j in range(3)] for n in u}
                dt=scalar(fo,'DAMAGET',wall);dc=scalar(fo,'DAMAGEC',wall)
                levelmean=[sum(du[n][1] for n in nodes)/len(nodes) for nodes in levels]
                paired=[max([abs(du[b][1]-du[a][1])/H for a,b in pp] or [0.]) for pp in pairs]
                row=dict(step=sn,frame=fi,time=float(fr.frameValue),dt_max=max(list(dt.values()) or [0.]),dc_max=max(list(dc.values()) or [0.]),wall_u_mm=max(abs(du[n][1]) for n in wn)*1000,dt_mean=sum(areas[e]*dt.get(e,0) for e in wall)/total,area_gt_05_pct=100*sum(areas[e] for e in wall if dt.get(e,0)>.5)/total,rf_y_N=sum(rf.get(n,[0,0,0])[1] for n in base_nodes),rf_z_N=sum(rf.get(n,[0,0,0])[2] for n in base_nodes),level_mean_u_mm=[v*1000 for v in levelmean],level_max_u_mm=[max(abs(du[n][1]) for n in nodes)*1000 for nodes in levels],story_mean_drift=[abs(levelmean[k+1]-levelmean[k])/H for k in range(nlev)],story_paired_max_drift=paired)
                if aa:row.update(level_mean_a=[sum(aa[n][1] for n in ns)/len(ns) for ns in levels],level_max_a=[max(abs(aa[n][1]) for n in ns) for ns in levels])
                if vv:row['level_max_v']=[max(abs(vv[n][1]) for n in ns) for ns in levels]
                row.update(energy_at(sn,fr.frameValue))
                row['faces']={g:dict(max_Dt=max([dt.get(e,0) for e in es] or [0.]),mean_Dt=sum(areas[e]*dt.get(e,0) for e in es)/sum(areas[e] for e in es),area_gt_05_pct=100*sum(areas[e] for e in es if dt.get(e,0)>.5)/sum(areas[e] for e in es)) for g,es in groups.items() if es}
                for e,d in dt.items():
                    if d>=.01 and e not in onset:onset[e]=dict(step=sn,time=fr.frameValue,frame=fi)
                hs.append(row);frame_refs.append((sn,fi))
        imax=max(range(len(hs)),key=lambda i:hs[i]['dt_max']);iu=max(range(len(hs)),key=lambda i:hs[i]['wall_u_mm']);first=next((i for i,h in enumerate(hs) if h['dt_max']>=.01),0)
        indices=set([0,first,imax,iu,len(hs)-1])
        if meta['hazard']=='DEBRIS':indices.update(min(range(len(hs)),key=lambda i:abs(hs[i]['time']-t)) for t in [.05,.1,.25,.35])
        elif meta['hazard']=='ROCK':indices.update(min(range(len(hs)),key=lambda i:abs(hs[i]['time']-t)) for t in [.015,.03,.1])
        else:indices.update(i for i in range(len(hs)) if i==len(hs)-1 or hs[i]['step']!=hs[i+1]['step'])
        snapshots=[]
        for idx in sorted(indices):
            sn,fi=frame_refs[idx];fr=odb.steps[sn].frames[fi];fo=fr.fieldOutputs
            dt=scalar(fo,'DAMAGET',wall);dc=scalar(fo,'DAMAGEC',wall);peq=scalar(fo,'PEEQ',wall);le=scalar(fo,'LE',wall,MAX_PRINCIPAL)
            stress={}
            for e,val in items(fo['S'],'element'):
                if e in beam:stress[e]=max(stress.get(e,0),abs(float(val[0])))
            u=nodal(fo,'U');du={n:u[n][1]-gu[n][1] for n in u}
            maps=[[dict(node=b,x=coords[b][0],y=coords[b][1],theta=abs(du[b]-du[a])/H) for a,b in pp] for pp in pairs]
            snapshots.append(dict(history_index=idx,step=sn,frame=fi,time=float(fr.frameValue),Dt=dt,Dc=dc,PEEQ=peq,LE_max_principal=le,beam_abs_S11_Pa=stress,story_drift_points=maps))
        cache_match=None
        if oldcache:
            existing=[h for h in oldcache['history'] if h['step'] in steps];assert len(existing)==len(hs)
            errors=dict(dt=max(abs(a['dt_max']-b['dt_max']) for a,b in zip(hs,existing)),u_mm=max(abs(a['wall_u_mm']-b['hazard_wall_uy_m']*1000) for a,b in zip(hs,existing)))
            print(case,'cache differences',errors)
            assert errors['dt']<1e-6 and errors['u_mm']<1e-4,errors
            cache_match=errors
        complete='THE ANALYSIS HAS COMPLETED SUCCESSFULLY' in open(base+'/'+case+'.sta').read()
        out=dict(case=case,metadata=meta,complete=complete,history=hs,snapshots=snapshots,first_Dt_001=onset,field_keys=seen_fields,gravity_energy=gravity,peak_damage_history_index=imax,peak_displacement_history_index=iu,face_element_groups=groups,cache_match=cache_match,read_only=True)
        json.dump(out,open(target,'w'));report.append(dict(case=case,frames=len(hs),snapshots=len(snapshots),A_present=all('level_mean_a' in h for h in hs),cache_match=cache_match,passed=True))
        print(case,'PASS',len(hs),'frames',len(snapshots),'snapshots');sys.stdout.flush()
    except Exception as e:
        report.append(dict(case=case,passed=False,error=str(e),traceback=traceback.format_exc()));print(case,'ERROR',str(e));sys.stdout.flush()
    finally:
        if odb is not None:odb.close()
json.dump(report,open(R+'/extraction_'+MODE+'_report.json','w'),indent=2)
print('Finished extraction',MODE,len(report))
