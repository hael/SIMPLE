"""Optional analytic known-start-volume solve3D adapter."""
import json
from resource_planner.workflows.common import SpatialAdapter, GENERATOR

class Solve3DAdapter(SpatialAdapter):
    name='solve3D'
    qualification='provisional synthetic single-stage known-start-volume screening; not full de novo solve3D; same-frame volume correlation >=0.8'

    def cases(self,args):return [args.box]

    def science(self,args):
        # nstages is the final stage index; supplied volumes start at stage 4
        # (after symmetry search), so index 4 executes one stage.
        return {'pgrp':'c1','pgrp_start':'c1','mskdiam':args.box*1.3*0.8,
                'nstates':1,'nstages':4,'balance':'none','nsample':args.particles,'filt_mode':'none',
                'automsk':'no','force_lp_range':'yes','lpstart':15,'lpstop':6}

    def metadata(self,args):
        return {'threads':args.threads,'workload':{'box':args.box,'particles':args.particles,
                'smpd':1.3,'ctf':'no','snr':0.1,'fixture_generator':'simulate_particles',
                'initialization':'known_volume',**self.science(args)}}

    def prepare(self,ctx,size):
        inputs=ctx.root/'fixture'
        ctx.command([ctx.args.fixture_python,str(GENERATOR),'volume',str(inputs),
                     '--box',str(size)],ctx.root,'fixture',timeout=300)
        ctx.command([str(ctx.args.simple_exec),'prg=simulate_particles',f'vol1={inputs}/volume.mrc',
                     f'outstk={inputs}/particles.mrcs',f'outfile={inputs}/orientations.txt',
                     f'nptcls={ctx.args.particles}','smpd=1.3',f'mskdiam={size*1.3*0.8}',
                     'nthr=1','pgrp=c1','snr=0.1','ctf=no','sherr=0','even=yes'],inputs,'simulate',timeout=ctx.args.timeout)
        (inputs/'halves.txt').write_text(''.join(
            f'state=1 eo={i % 2}\n' for i in range(ctx.args.particles)))
        return {'folder':ctx.root,'size':size,'input':inputs,
                'files':[inputs/name for name in ('volume.mrc','particles.mrcs','orientations.txt','halves.txt','truth.json')]}

    def jobs(self,ctx,fixture,directory,nthr,workers):
        project=ctx.project(directory,'reconstruction',fixture['input']/'particles.mrcs',particles=True,ctf_flag='no',
                            oritab=fixture['input']/'halves.txt')
        # Supplied-volume initialization supports fresh particles with balance=none.
        params={**self.science(ctx.args),'nthr':nthr,'nthr_ini3D':nthr,
                'vol1':str(fixture['input']/'volume.mrc'),'memreport':'yes','memreport_interval':1}
        return [([str(ctx.args.simple_exec),'prg=solve3D',f'projfile={project}','mkdir=no',
                  'qsys_name=local',*[f'{k}={v}' for k,v in params.items()]],directory)]

    def validate(self,ctx,fixture,work):
        volumes=[p for p in work.glob('recvol_state*.mrc') if re.fullmatch(r'recvol_state\d+\.mrc',p.name)]
        if len(volumes)!=1:raise ValueError('expected one reconstructed state volume')
        if not re.search(r'ITERATION\s+\d+', (work/'run.log').read_text()):
            raise ValueError('no reconstruction iteration observed')
        ctx.command([str(ctx.args.simple_exec),'prg=print_project_field',
                     f'projfile={work}/reconstruction/reconstruction.simple',
                     'oritype=ptcl3D','mkdir=no'],work,'assignments')
        rows=[]
        for line in (work/'assignments.log').read_text().splitlines():
            values=dict(re.findall(r'(\w+)=\s*([^\s]+)',line))
            if 'state' in values:rows.append(values)
        if len(rows)!=ctx.args.particles:raise ValueError('3D particle assignment count mismatch')
        for row in rows:
            values=[float(row.get(k,'0')) for k in ('e1','e2','e3','state','updatecnt')]
            if not all(math.isfinite(v) for v in values) or values[3]!=1:
                raise ValueError('invalid 3D particle orientations/state')
        result=self.check_with_numpy(ctx,fixture,work,'validate_volume',volumes[0])
        updated=sum(float(row.get('updatecnt','0'))>0 for row in rows)
        result.update(particles_assigned=len(rows),particles_updated=updated)
        result['quality_qualified']=result['quality_qualified'] and updated==len(rows)
        return result
