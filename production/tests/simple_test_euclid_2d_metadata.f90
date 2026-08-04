program simple_test_euclid_2d_metadata
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_builder,    only: builder
use simple_cmdline,    only: cmdline
use simple_core_module_api
use simple_parameters, only: parameters
use simple_sp_project, only: sp_project
use simple_string,     only: string
implicit none
#include "simple_local_flags.inc"

type(cmdline)   :: args, cline
type(parameters) :: params
type(builder)   :: build
type(sp_project) :: spproj
character(len=512) :: projfile
type(string) :: proj_arg, inpl_mode
real :: smpd, mskdiam, lp
integer :: box_crop
integer :: nptcls, ncls, nrots, iptcl, inpl
integer :: n_active, n_offgrid, n_invalid_inpl, n_nonfinite_e3, ncls_out
integer :: n_offgrid_logged
real :: e3, grid_e3, delta_e3, dang
logical :: expect_continuous

call args%parse_oldschool
if( command_argument_count() < 5 )then
    write(logfhandle,'(a)') &
        'simple_test_euclid_2d_metadata projfile=xx mskdiam=xx lp=xx smpd=xx box_crop=xx inpl_refine=yes|no'
    stop 2
endif
call args%checkvar('projfile', 1)
call args%checkvar('mskdiam', 2)
call args%checkvar('lp', 3)
call args%checkvar('smpd', 4)
call args%checkvar('box_crop', 5)
proj_arg = args%get_carg('projfile')
projfile = proj_arg%to_char()
mskdiam  = args%get_rarg('mskdiam')
lp       = args%get_rarg('lp')
smpd     = args%get_rarg('smpd')
box_crop = args%get_iarg('box_crop')
expect_continuous = .false.
inpl_mode = string('no')
if( args%defined('inpl_refine') )then
    inpl_mode = args%get_carg('inpl_refine')
    select case(trim(inpl_mode%to_char()))
    case('yes')
        expect_continuous = .true.
    case('no')
        expect_continuous = .false.
    case default
        THROW_HARD('inpl_refine must be yes or no')
    end select
endif

call spproj%read(string(trim(projfile)))
nptcls = spproj%os_ptcl2D%get_noris()
ncls_out = spproj%os_out%get_noris()
if( nptcls < 1 ) THROW_HARD('2D metadata test found no particles')

ncls = max(1, ncls_out)
call cline%set('projfile', trim(projfile))
call cline%set('smpd', smpd)
call cline%set('mskdiam', mskdiam)
call cline%set('lp', lp)
call cline%set('box_crop', box_crop)
call cline%set('ncls', ncls)
call cline%set('nptcls', nptcls)
call cline%set('ctf', 'no')
call cline%set('objfun', 'euclid')
call cline%set('inpl_refine', inpl_mode%to_char())
call build%init_params_and_build_strategy2D_tbox(cline, params, wthreads=.false.)
call build%pftc%new(params, nptcls, [1, nptcls], params%kfromto)
nrots = build%pftc%get_nrots()
dang  = build%pftc%get_dang()
if( nrots < 1 .or. dang <= 0. ) THROW_HARD('invalid polar rotation grid')

n_active = 0
n_offgrid = 0
n_invalid_inpl = 0
n_nonfinite_e3 = 0
n_offgrid_logged = 0
do iptcl = 1, nptcls
    if( spproj%os_ptcl2D%get(iptcl, 'state') <= 0. ) cycle
    n_active = n_active + 1
    inpl = nint(spproj%os_ptcl2D%get(iptcl, 'inpl'))
    e3 = spproj%os_ptcl2D%e3get(iptcl)
    if( inpl < 1 .or. inpl > nrots )then
        n_invalid_inpl = n_invalid_inpl + 1
        cycle
    endif
    if( .not. ieee_is_finite(real(e3)) )then
        n_nonfinite_e3 = n_nonfinite_e3 + 1
        cycle
    endif
    grid_e3 = 360. - build%pftc%get_rot(inpl)
    delta_e3 = abs(e3 - grid_e3)
    if( delta_e3 > 180. ) delta_e3 = 360. - delta_e3
    if( delta_e3 > max(1.e-4, 0.01*dang) )then
        n_offgrid = n_offgrid + 1
        if( n_offgrid_logged < 8 )then
            write(logfhandle,'(a,2(i8,1x),a,3(f14.6,1x))') &
                'EUCLID_2D_METADATA OFFGRID_SAMPLE IPTCL/INPL: ', iptcl, inpl, &
                'E3/GRID/DELTA: ', e3, grid_e3, delta_e3
            n_offgrid_logged = n_offgrid_logged + 1
        endif
    endif
enddo

write(logfhandle,'(a,4i8)') 'EUCLID_2D_METADATA ACTIVE/OFFGRID/INVALID/NONFINITE: ', &
    n_active, n_offgrid, n_invalid_inpl, n_nonfinite_e3
write(logfhandle,'(a,2i8)') 'EUCLID_2D_METADATA ROTATIONS/CLASS_AVERAGES: ', nrots, ncls_out
write(logfhandle,'(a,f12.6)') 'EUCLID_2D_METADATA ANGLE_STEP_DEG: ', dang
write(logfhandle,'(a,l1)') 'EUCLID_2D_METADATA CONTINUOUS_EXPECTED: ', expect_continuous
write(logfhandle,'(a,i8)') 'EUCLID_2D_METADATA CONTINUOUS_OFFGRID_COUNT: ', n_offgrid

if( n_active < 1 ) THROW_HARD('2D metadata test found no active particles')
if( n_invalid_inpl /= 0 ) THROW_HARD('2D metadata contains invalid integer inpl values')
if( n_nonfinite_e3 /= 0 ) THROW_HARD('2D metadata contains non-finite e3 values')
if( .not. expect_continuous .and. n_offgrid /= 0 ) &
    THROW_HARD('2D metadata contains unexpected off-grid e3 values')
if( ncls_out < 1 ) THROW_HARD('2D metadata contains no class-average metadata')

call build%kill_strategy2D_tbox
call build%kill_general_tbox
call spproj%kill
call cline%kill
call args%kill
call proj_arg%kill
call inpl_mode%kill
call simple_end('SIMPLE_TEST_EUCLID_2D_METADATA NORMAL STOP')
end program simple_test_euclid_2d_metadata
