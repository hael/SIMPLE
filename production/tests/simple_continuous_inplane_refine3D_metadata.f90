module continuous_inplane_refine3D_metadata_test
implicit none
private
public :: run_metadata_state_contract, run_metadata_project_contract

contains

subroutine run_metadata_state_contract()
use simple_strategy3D_utils, only: resolve_inplane_e3
implicit none

real, parameter :: dang = 1.25
real :: e3

e3 = resolve_inplane_e3(17, 17.4, .true., .false., dang)
if( abs(e3 - 340.) > 1.e-5 ) &
    &error stop 'default-off metadata did not retain the integer grid angle'

e3 = resolve_inplane_e3(17, 17.4, .false., .true., dang)
if( abs(e3 - 340.) > 1.e-5 ) &
    &error stop 'invalid continuous metadata replaced the integer grid angle'

e3 = resolve_inplane_e3(17, 17.4, .true., .true., dang)
if( abs(e3 - 339.5) > 1.e-5 ) &
    &error stop 'accepted continuous coordinate was not converted to e3'

write(*,'(a)') 'REFINE3D_INPL_CONT_METADATA GRID/CONTINUOUS: PASS'
end subroutine run_metadata_state_contract

subroutine run_metadata_project_contract()
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_sp_project, only: sp_project
use simple_string, only: string
implicit none
#include "simple_local_flags.inc"

type(sp_project) :: spproj
character(len=4096) :: argument, key, value, projfile
character(len=16) :: inpl_mode
integer :: iarg, separator, iptcl, nptcls, nactive, ncontinuous
integer :: ninvalid, nnonfinite, nnearest_mismatch, inpl, nearest_inpl, nrots
integer :: min_inpl, max_inpl
real :: dang, e3, grid_e3, delta_e3, rot, shift(2), corr, tolerance
logical :: project_exists, expect_continuous

projfile = ''
inpl_mode = 'no'
dang = 0.
nrots = 0
do iarg = 1, command_argument_count()
    call get_command_argument(iarg, argument)
    separator = index(argument, '=')
    if( separator <= 1 ) cycle
    key   = adjustl(argument(:separator-1))
    value = adjustl(argument(separator+1:))
    select case(trim(key))
    case('projfile')
        projfile = trim(value)
    case('dang')
        read(value,*) dang
    case('nrots')
        read(value,*) nrots
    case('inpl_cont')
        inpl_mode = trim(value)
    case('case')
        continue
    case default
        write(logfhandle,'(a)') 'WARNING: ignored argument '//trim(argument)
    end select
enddo

if( dang <= 0. .and. nrots > 0 ) dang = 360. / real(nrots)
if( nrots <= 0 .and. dang > 0. ) nrots = nint(360. / dang)
if( len_trim(projfile) == 0 .or. dang <= 0. .or. nrots <= 0 )then
    write(logfhandle,'(a)') &
        'metadata_project requires projfile=xx inpl_cont=yes|no and dang=xx or nrots=xx'
    stop 2
endif
select case(trim(inpl_mode))
case('yes')
    expect_continuous = .true.
case('no')
    expect_continuous = .false.
case default
    THROW_HARD('inpl_cont must be yes or no; metadata_project')
end select

inquire(file=trim(projfile), exist=project_exists)
if( .not. project_exists ) THROW_HARD('refine3D metadata project does not exist')
call spproj%read(string(trim(projfile)))
nptcls = spproj%os_ptcl3D%get_noris()
if( nptcls < 1 ) THROW_HARD('refine3D metadata project contains no ptcl3D records')

tolerance = max(1.e-4, 0.01*dang)
nactive = 0
ncontinuous = 0
ninvalid = 0
nnonfinite = 0
nnearest_mismatch = 0
min_inpl = huge(1)
max_inpl = 0
do iptcl = 1, nptcls
    if( spproj%os_ptcl3D%get(iptcl, 'state') <= 0. ) cycle
    nactive = nactive + 1
    inpl = nint(spproj%os_ptcl3D%get(iptcl, 'inpl'))
    min_inpl = min(min_inpl, inpl)
    max_inpl = max(max_inpl, inpl)
    e3 = spproj%os_ptcl3D%e3get(iptcl)
    shift = spproj%os_ptcl3D%get_2Dshift(iptcl)
    corr = spproj%os_ptcl3D%get(iptcl, 'corr')
    if( inpl < 1 .or. inpl > nrots )then
        ninvalid = ninvalid + 1
        cycle
    endif
    if( .not. ieee_is_finite(e3) .or. .not. all(ieee_is_finite(shift)) .or. &
        &.not. ieee_is_finite(corr) )then
        nnonfinite = nnonfinite + 1
        cycle
    endif
    grid_e3 = 360. - real(inpl-1) * dang
    delta_e3 = abs(modulo(e3 - grid_e3 + 180., 360.) - 180.)
    if( delta_e3 > tolerance ) ncontinuous = ncontinuous + 1
    rot = modulo(360. - e3, 360.)
    nearest_inpl = modulo(nint(rot / dang), nrots) + 1
    if( nearest_inpl /= inpl ) nnearest_mismatch = nnearest_mismatch + 1
enddo

write(logfhandle,'(a)') 'REFINE3D_INPL_CONT_METADATA PROJECT: '//trim(projfile)
write(logfhandle,'(a,5(i0,1x))') &
    'REFINE3D_INPL_CONT_METADATA ACTIVE/CONTINUOUS/INVALID/NONFINITE/MISMATCH: ', &
    nactive, ncontinuous, ninvalid, nnonfinite, nnearest_mismatch
write(logfhandle,'(a,f12.6)') 'REFINE3D_INPL_CONT_METADATA ANGLE_STEP_DEG: ', dang
write(logfhandle,'(a,2(i0,1x))') 'REFINE3D_INPL_CONT_METADATA INPL_MIN/MAX: ', min_inpl, max_inpl
write(logfhandle,'(a,l1)') 'REFINE3D_INPL_CONT_METADATA CONTINUOUS_EXPECTED: ', expect_continuous

if( nactive < 1 ) THROW_HARD('refine3D metadata project contains no active particles')
if( ninvalid /= 0 ) THROW_HARD('refine3D metadata contains invalid integer inpl values')
if( nnonfinite /= 0 ) THROW_HARD('refine3D metadata contains non-finite pose values')
if( nnearest_mismatch /= 0 ) &
    &THROW_HARD('refine3D metadata integer and continuous angles are inconsistent')
if( expect_continuous .and. ncontinuous == 0 ) &
    &THROW_HARD('refine3D opt-in metadata contains no continuous e3 values')
if( (.not. expect_continuous) .and. ncontinuous /= 0 ) &
    &THROW_HARD('refine3D default-off metadata contains continuous e3 values')

call spproj%kill
call simple_end('SIMPLE_TEST_CONT_INPL_REFINE3D_METADATA NORMAL STOP')
end subroutine run_metadata_project_contract

end module continuous_inplane_refine3D_metadata_test
