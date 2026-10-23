"""SI units m N kg s. Assumed archetypes, not surveyed assets."""
from pathlib import Path
import json, math, csv, os
ROOT=Path(os.environ.get('GANSU_FE_OUTPUT', 'results/generated_fe')).resolve();ROOT.mkdir(parents=True,exist_ok=True)
MATS={'earth':dict(E=40e6,nu=.30,rho=1900,fc=1e6,ft=30000,u=.020),
      'brick':dict(E=3.2e9,nu=.18,rho=1800,fc=2e6,ft=120000,u=.002),
      'stone':dict(E=2.25e9,nu=.20,rho=2200,fc=1.5e6,ft=80000,u=.002)}
MODELS=[
dict(id='LZ1',label='兰州土木单层原型',wall='earth',t=.40,L=9.,W=6.,H=3.,n=1,frame='wood',col=.20,beam=(.18,.24),roof='tile'),
dict(id='LZ2',label='兰州单层砖混原型',wall='brick',t=.24,L=10.,W=7.,H=3.,n=1,frame='rc',col=.24,beam=(.24,.30),roof='rc'),
dict(id='LZ3',label='兰州二层框架填充墙原型',wall='brick',t=.18,L=9.,W=7.5,H=3.,n=2,frame='rc',col=.30,beam=(.25,.45),roof='rc'),
dict(id='TS1',label='天水单层生土原型',wall='earth',t=.45,L=8.,W=6.,H=3.,n=1,frame='wood',col=.18,beam=(.16,.22),roof='tile'),
dict(id='TS2',label='天水砖木单层原型',wall='brick',t=.37,L=9.,W=6.,H=3.,n=1,frame='wood',col=.20,beam=(.18,.24),roof='tile'),
dict(id='TS3',label='天水二层砖混原型',wall='brick',t=.24,L=9.,W=7.5,H=3.,n=2,frame='rc',col=.24,beam=(.24,.30),roof='rc'),
dict(id='LN1',label='陇南木构架土围护原型',wall='earth',t=.20,L=9.,W=6.,H=3.,n=1,frame='wood',col=.24,beam=(.20,.28),roof='tile'),
dict(id='LN2',label='陇南二层石木原型',wall='stone',t=.45,L=8.,W=6.,H=2.8,n=2,frame='wood',col=.24,beam=(.20,.28),roof='slate'),
dict(id='LN3',label='陇南单层砖混原型',wall='brick',t=.24,L=10.,W=7.,H=3.,n=1,frame='rc',col=.24,beam=(.24,.30),roof='rc')]

def make_mesh(m,mesh=.5):
    nodes={};coords={};elems={};groups={};areas={};front=[];beam_seen={};promoted=[]
    wall_nodes={}
    def node(xyz,kind='frame'):
        xyz=tuple(round(v,8) for v in xyz)
        key=(xyz,'wall' if m['id']=='LN1' and kind=='wall' else 'shared')
        if key not in nodes:
            k=len(nodes)+1;nodes[key]=k;coords[k]=xyz
        if kind=='wall':wall_nodes[xyz]=nodes[key]
        return nodes[key]
    def quad(xyz,g):
        n=[node(v,'wall' if g=='WALL' else 'frame') for v in xyz];eid=len(elems)+1;elems[eid]=n;groups.setdefault(g,[]).append(eid)
        a=[xyz[1][i]-xyz[0][i] for i in range(3)];b=[xyz[3][i]-xyz[0][i] for i in range(3)]
        areas[eid]=math.sqrt(sum(v*v for v in (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])))
        return eid
    def line(a,b,g):
        dist=math.dist(a,b);N=max(1,round(dist/mesh))
        for j in range(N):
            v=[node(tuple(a[k]+(b[k]-a[k])*s/N for k in range(3))) for s in (j,j+1)]
            key=tuple(sorted(v))
            if key in beam_seen:
                old_e,old_g=beam_seen[key]
                if old_g in ('JOIST','PURLIN') and g=='BEAM':
                    groups[old_g].remove(old_e);groups.setdefault(g,[]).append(old_e)
                    beam_seen[key]=(old_e,g);promoted.append({'element':old_e,'from':old_g,'to':g})
                else:
                    assert old_g==g, 'Unexpected coincident structural members: %s / %s'%(old_g,g)
                continue
            eid=len(elems)+1;beam_seen[key]=(eid,g);elems[eid]=v;groups.setdefault(g,[]).append(eid)
    def triangle(xyz,g):
        n=[node(v) for v in xyz];eid=len(elems)+1;elems[eid]=n;groups.setdefault(g,[]).append(eid)
        a=[xyz[1][i]-xyz[0][i] for i in range(3)];b=[xyz[2][i]-xyz[0][i] for i in range(3)]
        areas[eid]=.5*math.sqrt(sum(v*v for v in (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])))
    L,W,H=m['L'],m['W'],m['H']; nx=3*max(2,round(L/3/mesh));ny=max(2,round(W/mesh));nz=max(2,round(H/mesh))
    ny=2*max(1,round(ny/2))
    xs=[round(L*i/nx,8) for i in range(nx+1)];ys=[round(W*i/ny,8) for i in range(ny+1)];zs=[round(H*i/nz,8) for i in range(nz*m['n']+1)]
    for z0,z1 in zip(zs,zs[1:]):
        zmid=(z0+z1)/2;zz=zmid%H
        for y in (0,W):
            for x0,x1 in zip(xs,xs[1:]):
                xx=(x0+x1)/2/L
                opening=(y==0 and ((.43<xx<.57 and ((zmid<H and zz<2.1) or .9<zz<2.1)) or ((.10<xx<.26 or .74<xx<.90) and .9<zz<2.1))) or (y==W and (.10<xx<.26 or .43<xx<.57 or .74<xx<.90) and .9<zz<2.1)
                if opening:continue
                eid=quad([(x0,y,z0),(x1,y,z0),(x1,y,z1),(x0,y,z1)],'WALL')
                if y==0:front.append((eid,(x0+x1)/2,zmid,areas[eid]))
        for x in (0,L/3,2*L/3,L):
            for y0,y1 in zip(ys,ys[1:]):
                yy=(y0+y1)/2/W
                if x not in (0,L) and .40<yy<.60 and zz<2.1:continue
                if x in (0,L) and .35<yy<.60 and .9<zz<2.1:continue
                quad([(x,y0,z0),(x,y1,z0),(x,y1,z1),(x,y0,z1)],'WALL')
    # Traditional roof: actual double pitch, overhang, ridge, rafters and purlins.
    for story in range(1,m['n']+1):
        z=story*H;g='ROOF' if story==m['n'] else 'FLOOR'
        if story==m['n'] and m['roof']!='rc':
            rise=1.2;overhang=.4
            xx=[-overhang]+xs+[round(L+overhang,8)];yy=sorted(set(round(v,8) for v in [-overhang]+ys+[W/2,W+overhang]))
            def rz(y):return round(z+rise*(1-abs(y-W/2)/(W/2)),8)
            for x0,x1 in zip(xx,xx[1:]):
                for y0,y1 in zip(yy,yy[1:]):quad([(x0,y0,rz(y0)),(x1,y0,rz(y0)),(x1,y1,rz(y1)),(x0,y1,rz(y1))],g)
            for x in xs:
                for y0,y1 in zip(yy,yy[1:]):line((x,y0,rz(y0)),(x,y1,rz(y1)),'RAFTER')
            for x in (0,L/3,2*L/3,L):
                line((x,W/2,z),(x,W/2,z+rise),'COL')
                # Tie beam at eaves level transfers ridge-post reactions to supports.
            for y in (0,W/2,W):
                for x0,x1 in zip(xx,xx[1:]):line((x0,y,rz(y)),(x1,y,rz(y)),'PURLIN')
            # Light timber-board gable enclosure; not counted as loadbearing masonry.
            for x in (0,L):
                yyinside=[y for y in yy if 0<=y<=W]
                for y0,y1 in zip(yyinside,yyinside[1:]):
                    poly=[(x,y0,z),(x,y1,z),(x,y1,rz(y1)),(x,y0,rz(y0))]
                    unique=[]
                    for a in poly:
                        if a not in unique:unique.append(a)
                    for k in range(1,len(unique)-1):triangle([unique[0],unique[k],unique[k+1]],'GABLE')
        else:
            for x0,x1 in zip(xs,xs[1:]):
                for y0,y1 in zip(ys,ys[1:]):quad([(x0,y0,z),(x1,y0,z),(x1,y1,z),(x0,y1,z)],g)
            if m['frame']=='wood':
                # Timber joists spanning between the cross-walls, floor boards above.
                for y in ys:
                    for x0,x1 in zip(xs,xs[1:]):line((x0,y,z),(x1,y,z),'JOIST')
        for x in (0,L/3,2*L/3,L):
            # Beam subdivisions match wall mesh exactly.
            for y0,y1 in zip(ys,ys[1:]):line((x,y0,z),(x,y1,z),'BEAM')
        for y in (0,W):
            for x0,x1 in zip(xs,xs[1:]):line((x0,y,z),(x1,y,z),'BEAM')
    for x in (0,L/3,2*L/3,L):
        for y in (0,W):
            for z0,z1 in zip(zs,zs[1:]):line((x,y,z0),(x,y,z1),'COL')
    if m['id']=='LN1':
        for xyz,nwall in wall_nodes.items():
            k=(xyz,'shared')
            if k in nodes and xyz[2]>0:
                for g in (('LINKX','LINKY','LINKZ') if m.get('vertical_wall_frame_link',True) else ('LINKX','LINKY')):
                    eid=len(elems)+1;elems[eid]=[nwall,nodes[k]];groups.setdefault(g,[]).append(eid)
    # Node registry rounds partition coordinates consistently; enforce no unconnected nodes.
    used={n for v in elems.values() for n in v};assert used==set(coords)
    assert all(a>1e-8 for a in areas.values()), 'Degenerate shell'
    assert all(len(set(v))==len(v) for g,es in groups.items() if not g.startswith('LINK') for e in es for v in [elems[e]]), 'Repeated element node'
    return dict(nodes=coords,elements=elems,groups=groups,areas=areas,front=front,promoted_edge_beams=promoted,
                base=[i for i,c in coords.items() if abs(c[2])<1e-7],
                levels={str(j):[i for i,c in coords.items() if abs(c[2]-j*H)<1e-7] for j in range(m['n']+1)})

def refine_mesh(mesh,m):
    # Uniform h-refinement of the SAME geometric boundary and openings.
    coords=dict(mesh['nodes']);registry={tuple(round(v,8) for v in xyz):n for n,xyz in coords.items()}
    def node(xyz):
        key=tuple(round(v,8) for v in xyz)
        if key not in registry:
            n=len(coords)+1;registry[key]=n;coords[n]=key
        return registry[key]
    def mid(a,b):return node([(coords[a][k]+coords[b][k])/2 for k in range(3)])
    elems={};groups={};areas={};front=[];frontids={e for e,_,_,_ in mesh['front']}
    for g,es in mesh['groups'].items():
        assert not g.startswith('LINK'),'Spring refinement requires tributary stiffness rescaling'
        for e in es:
            ns=mesh['elements'][e]
            if len(ns)==2:
                a,b=ns;c=mid(a,b);children=[[a,c],[c,b]]
            elif len(ns)==3:
                a,b,c=ns;ab=mid(a,b);bc=mid(b,c);ca=mid(c,a);children=[[a,ab,ca],[ab,b,bc],[ca,bc,c],[ab,bc,ca]]
            else:
                a,b,c,d=ns;ab=mid(a,b);bc=mid(b,c);cd=mid(c,d);da=mid(d,a);center=node([sum(coords[n][k] for n in ns)/4 for k in range(3)])
                children=[[a,ab,center,da],[ab,b,bc,center],[center,bc,c,cd],[da,center,cd,d]]
            for child in children:
                en=len(elems)+1;elems[en]=child;groups.setdefault(g,[]).append(en)
                if e in mesh['areas']:areas[en]=mesh['areas'][e]/len(children)
                if e in frontids:front.append((en,sum(coords[n][0] for n in child)/len(child),sum(coords[n][2] for n in child)/len(child),areas[en]))
    return dict(nodes=coords,elements=elems,groups=groups,areas=areas,front=front,base=[n for n,c in coords.items() if abs(c[2])<1e-7],levels={str(j):[n for n,c in coords.items() if abs(c[2]-j*m['H'])<1e-7] for j in range(m['n']+1)})

def ids(v):return '\n'.join(', '.join(map(str,v[i:i+14])) for i in range(0,len(v),14))
def material(m,nonlinear):
    a=m.get('wall_material',MATS[m['wall']]);s=['*MATERIAL, NAME=WALLMAT','*DENSITY',str(a['rho']),'*ELASTIC','%g, %g'%(a['E'],a['nu'])]
    if nonlinear:
        fc,ft,u=a['fc'],a['ft'],a['u']
        s+=['*CONCRETE DAMAGED PLASTICITY','25., 0.1, 1.16, 0.667, 0.0001',
            '*CONCRETE COMPRESSION HARDENING','%g, 0.'%(.4*fc),'%g, 0.005'%fc,'%g, 0.04'%(.2*fc),
            '*CONCRETE COMPRESSION DAMAGE','0., 0.','0.1, 0.005','0.7, 0.04',
            '*CONCRETE TENSION STIFFENING, TYPE=DISPLACEMENT','%g, 0.'%ft,'%g, %g'%(.5*ft,u*.25),'%g, %g'%(.1*ft,u),'0., %g'%(2*u),
            '*CONCRETE TENSION DAMAGE, TYPE=DISPLACEMENT','0., 0.','0.5, %g'%(u*.25),'0.9, %g'%u,'0.99, %g'%(2*u)]
    if m['frame']=='wood':E,rho,nu=9e9,500,.3
    else:E,rho,nu=28e9,2500,.2
    s+=['*MATERIAL, NAME=FRAMEMAT','*DENSITY',str(rho),'*ELASTIC','%g, %g'%(E,nu),
        '*MATERIAL, NAME=DECKMAT','*DENSITY',str(2500 if m['roof'] in ('rc','tile') else 4000),
        '*ELASTIC','28000000000., .2' if m['roof']=='rc' else '50000000., .25']
    return s

def loads(mesh,hazard,intensity=1,depth=None):
    f={};meta={}
    if hazard=='SLIDE':
        pmax=15000*intensity;h=1.5 if depth is None else depth;meta=dict(p_base_Pa=pmax,loaded_height_m=h)
        rule=lambda x,z:pmax*max(0,1-z/h) if z<h else 0
    elif hazard=='DEBRIS':
        rho,h,v=1800.,(1.5 if depth is None else depth),3.*math.sqrt(intensity);meta=dict(density_kg_m3=rho,depth_m=h,velocity_m_s=v,mean_gross_pressure_Pa=rho*v*v+.5*rho*9.81*h)
        rule=lambda x,z:rho*v*v+rho*9.81*(h-z) if z<h else 0
    else:
        # Triangle impulse: integral F dt = m v. No contact/penetration equivalence claimed.
        mass,v,tau=300.*intensity,6.,.03;impulse=mass*v;peak=2*impulse/tau
        patch=[(eid,x,z,a) for eid,x,z,a in mesh['front'] if .4<x<1.6 and .1<z<1.05]
        area=sum(a for _,_,_,a in patch);meta=dict(mass_kg=mass,velocity_m_s=v,pulse_s=tau,impulse_Ns=impulse,kinetic_energy_J=.5*mass*v*v,patch_area_m2=area,peak_force_N=peak)
        pe={e for e,_,_,_ in patch};rule=None
    for eid,x,z,area in mesh['front']:
        force=(peak/area if False else 0)
        pressure=(peak/meta['patch_area_m2'] if eid in pe else 0) if hazard=='ROCK' else rule(x,z)
        ns=mesh['elements'][eid]
        if hazard in ('SLIDE','DEBRIS'):
            # Exact consistent nodal forces for a linear pressure clipped at flow depth.
            z0=mesh['nodes'][ns[0]][2];z1=mesh['nodes'][ns[2]][2]
            top=min(z1,h);height=z1-z0
            if top<=z0:continue
            dz=top-z0;width=area/height
            a=(pmax*(1-z0/h) if hazard=='SLIDE' else rho*v*v+rho*9.81*(h-z0))
            b=(-pmax/h if hazard=='SLIDE' else -rho*9.81)
            integral=a*dz+b*dz*dz/2
            first=a*dz*dz/2+b*dz**3/3
            lower=width*(integral-first/height)/2;upper=width*first/height/2
            for n,val in zip(ns,(lower,lower,upper,upper)):f[n]=f.get(n,0)+val
        else:
            for n in ns:f[n]=f.get(n,0)+pressure*area/4
    f={n:v for n,v in f.items() if v};meta['resultant_peak_N']=sum(f.values())
    return f,meta

def write_case(m,hazard='ELASTIC',intensity=1,meshsize=.5,suffix=''):
    mesh=refine_mesh(make_mesh(m,.5),m) if meshsize<.5 else make_mesh(m,meshsize)
    if meshsize<.5:meshsize=.25
    nonlin=hazard!='ELASTIC';name=m['id']+'_'+hazard+suffix
    s=['*HEADING',name+' | ASSUMED ARCHETYPE | SI m N kg s | NOT SITE RISK','*PREPRINT, ECHO=NO, MODEL=NO, HISTORY=NO, CONTACT=NO','*NODE']
    s+=['%d, %.8g, %.8g, %.8g'%((k,)+tuple(v)) for k,v in mesh['nodes'].items()]
    for g,els in mesh['groups'].items():
        et='SPRING2' if g.startswith('LINK') else ('B31' if g in ('BEAM','COL','RAFTER','PURLIN','JOIST') else ('S3' if g=='GABLE' else 'S4R'))
        s+=['*ELEMENT, TYPE=%s, ELSET=%s'%(et,g)]
        s+=[str(e)+', '+', '.join(map(str,mesh['elements'][e])) for e in els]
    s+=['*NSET, NSET=BASE',ids(mesh['base']),'*NSET, NSET=ALLN',ids(list(mesh['nodes'])),'*ELSET, ELSET=ALLE',ids(list(mesh['elements']))]
    s+=['*ELSET, ELSET=GRAVITY_ELEMENTS',ids([e for g,es in mesh['groups'].items() if not g.startswith('LINK') for e in es])]
    for level,n in mesh['levels'].items():s+=['*NSET, NSET=LEVEL'+level,ids(n)]
    s+=material(m,nonlin)
    s+=['*SHELL SECTION, ELSET=WALL, MATERIAL=WALLMAT','%g, 5'%m['t']]
    for g in ('ROOF','FLOOR'):
        if g in mesh['groups']:
            mat='FRAMEMAT' if g=='FLOOR' and m['frame']=='wood' else 'DECKMAT'
            s+=['*SHELL SECTION, ELSET=%s, MATERIAL=%s'%(g,mat),('.12, 5' if m['roof']=='rc' else '.03, 5' if g=='FLOOR' else '.05, 5')]
    if 'GABLE' in mesh['groups']:s+=['*SHELL SECTION, ELSET=GABLE, MATERIAL=FRAMEMAT','.02, 3']
    for g,dim in [('BEAM',m['beam']),('COL',(m['col'],m['col'])),('RAFTER',(.10,.14)),('PURLIN',(.12,.16)),('JOIST',(.10,.18))]:
        if g not in mesh['groups']:continue
        # RECT dimension a follows n1. For horizontal beams n1 is vertical:
        # write height,width, keeping the intended physical height vertical.
        dims=(dim[1],dim[0]) if g in ('BEAM','PURLIN','JOIST') else dim
        s+=['*BEAM SECTION, ELSET=%s, MATERIAL=FRAMEMAT, SECTION=RECT'%g,'%g, %g'%tuple(dims),('0., 0., -1.' if g in ('BEAM','PURLIN','JOIST') else '1., 0., 0.')]
    # A pinned bending limit at physical timber beam-column joints, only in
    # explicitly requested sensitivity cases; never at all mesh subdivisions.
    releases=[]
    if m.get('joint_bending_release',False):
        assert m['frame']=='wood'
        for e in mesh['groups'].get('BEAM',[]):
            for end,n in enumerate(mesh['elements'][e],1):
                x,y,z=mesh['nodes'][n]
                if any(abs(x-m['L']*j/3)<1e-6 for j in range(4)) and min(abs(y),abs(y-m['W']))<1e-6 and any(abs(z-m['H']*j)<1e-6 for j in range(1,m['n']+1)):
                    releases.append('%d, S%d, M1-M2'%(e,end))
        s+=['*RELEASE']+releases
    for g,dof in [('LINKX',1),('LINKY',2),('LINKZ',3)]:
        if g in mesh['groups']:s+=['*SPRING, ELSET='+g,'%d, %d'%(dof,dof),str(m.get('link_k',1e5))]
    s+=['*BOUNDARY','BASE, 1, 6']
    def output():
        out=['*OUTPUT, FIELD, FREQUENCY=1','*NODE OUTPUT','U, RF','*ELEMENT OUTPUT, ELSET=WALL','S, LE'+(', DAMAGET, DAMAGEC, PEEQ' if nonlin else '')]
        for g in ('BEAM','COL','RAFTER','PURLIN','JOIST'):
            if g in mesh['groups']:out+=['*ELEMENT OUTPUT, ELSET='+g,'S, SF, SM']
        return out+['*OUTPUT, HISTORY, FREQUENCY=1','*ENERGY OUTPUT','ALLIE, ALLSE, ALLKE, ALLWK, ALLSD, ALLVD, ALLAE']
    s+=['*STEP, NAME=GRAVITY, NLGEOM=YES, INC=600','*STATIC','0.1, 1., 1e-7, 0.2','*DLOAD','GRAVITY_ELEMENTS, GRAV, 9.81, 0., 0., -1.']+output()+['*END STEP']
    meta=dict(model=m,hazard=hazard,intensity=intensity,mesh_m=meshsize,name=name,nonlinear_walls=nonlin,revision='v4_section_edge_beam_and_zero_tension_tail',wall_material=m.get('wall_material',MATS[m['wall']]),beam_end_bending_releases=len(releases))
    if hazard=='ELASTIC':
        s+=['*STEP, NAME=MODES, PERTURBATION','*FREQUENCY, EIGENSOLVER=LANCZOS','6,','*OUTPUT, FIELD','*NODE OUTPUT','U','*END STEP']
    else:
        f,force_meta=loads(mesh,hazard,intensity,depth=m.get('debris_depth'));meta.update(force_meta)
        if hazard=='SLIDE':
            pressures=m.get('slide_pressures_Pa',[5000.,15000.,30000.])
            meta['slide_step_pressures_Pa']={ 'SLIDE_L%d'%(j+1):p*intensity for j,p in enumerate(pressures)}
            for j,mul in enumerate([p/15000 for p in pressures],1):
                s+=['*STEP, NAME=SLIDE_L%d'%j+', NLGEOM=YES, INC=1200','*STATIC, STABILIZE=0.0002','0.03, 1., 1e-7, .1','*CLOAD, OP=NEW']
                s+=['%d, 2, %.9g'%(n,v*mul) for n,v in f.items()];s+=output()+['*END STEP']
        else:
            # Dynamics starts from converged static gravity; finite pulse duration.
            t=.5 if hazard=='DEBRIS' else .30
            s+=['*AMPLITUDE, NAME=PULSE, DEFINITION=TABULAR, TIME=STEP TIME']
            s+=['0., 0., .05, 1., .25, 1., .35, 0., .5, 0.' if hazard=='DEBRIS' else '0., 0., .015, 1., .03, 0., .30, 0.']
            # Amplitudes are model data, move them before the first step.
            amp=s[-2:];s=s[:-2];pos=next(i for i,x in enumerate(s) if x.startswith('*STEP'));s[pos:pos]=amp
            s+=['*STEP, NAME='+hazard+', NLGEOM=YES, INC=2400','*DYNAMIC, APPLICATION=TRANSIENT FIDELITY','0.001, %g, 1e-7, .003'%t,'*CLOAD, AMPLITUDE=PULSE, OP=NEW']
            s+=['%d, 2, %.9g'%(n,v) for n,v in f.items()];s+=output()+['*END STEP']
    (ROOT/(name+'.inp')).write_text('\n'.join(s)+'\n',encoding='ascii')
    (ROOT/(name+'_mesh.json')).write_text(json.dumps(mesh),encoding='ascii')
    (ROOT/(name+'_meta.json')).write_text(json.dumps(meta,ensure_ascii=True,indent=2),encoding='ascii')
    return meta

if __name__=='__main__':
    cases=[]
    for m in MODELS:
        for hazard in ('ELASTIC','SLIDE','DEBRIS','ROCK'):cases.append(write_case(m,hazard))
    (ROOT/'manifest.json').write_text(json.dumps(cases,ensure_ascii=True,indent=2),encoding='ascii')
    print('Generated %d cases in %s'%(len(cases),ROOT))
