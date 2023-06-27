program simple_test_refine3D
include 'simple_lib.f08'
use simple_image,               only: image
use simple_commander_project,   only: new_project_commander, import_particles_commander
use simple_cmdline,             only: cmdline
use simple_commander_cluster2D, only: cluster2D_commander_distr
implicit none
type(new_project_commander)      :: new_project_com
type(import_particles_commander) :: import_particles_com
type(cluster2D_commander_distr)  :: cluster2D_com
type(cmdline)                    :: new_project_cline, cline, import_particles_cline, cluster2D_cline
character(len=LONGSTRLEN)        :: cwd=''
if( command_argument_count() < 2 )then
    write(logfhandle,'(a)') 'Usage: simple_test_refine3D nparts=xx nthr=yy'
    stop
else
    call cline%parse_oldschool
endif
call cline%checkvar('nparts', 1)
call cline%checkvar('nthr',   2)

!---------- bgal dataset -----------------
! new project
new_project_cline = cline
call new_project_cline%set('projname', 'bgal')
call new_project_com%execute(new_project_cline)
!call simple_getcwd(cwd)

! import particles
import_particles_cline = cline
call import_particles_cline%delete('smpd')
call import_particles_cline%set('smpd',   1.275)
call import_particles_cline%set('deftab', '/mnt/beegfs/elmlund/testbench/bgal/deftab.txt')
call import_particles_cline%set('stk',    '/mnt/beegfs/elmlund/testbench/bgal/sumstack.mrc')
call import_particles_cline%set('cs',     2.7)
call import_particles_cline%set('fraca',  0.1)
call import_particles_cline%set('kv',     300.)
call import_particles_cline%set('projfile', 'bgal.simple')
call import_particles_com%execute(import_particles_cline)

! cluster2D
cluster2D_cline = cline
call cluster2D_cline%delete('smpd')
call cluster2D_cline%set('smpd',    1.275)
call cluster2D_cline%set('ncls',    90.)
call cluster2D_cline%set('mskdiam', 180.)
call cluster2D_cline%set('nparts',  6.)
call cluster2D_cline%set('nthr',    10.)
call cluster2D_cline%set('objfun',  'euclid')
call cluster2D_cline%set('projfile', 'bgal.simple')
call cluster2D_com%execute(cluster2D_cline)
end program simple_test_refine3D
