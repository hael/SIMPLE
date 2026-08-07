module continuous_inplane_refine3D_baseline_test
implicit none
private
public :: run_legacy_baseline

contains

subroutine run_legacy_baseline()
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_sp_project, only: sp_project
use simple_string,     only: string
implicit none
#include "simple_local_flags.inc"

type(sp_project) :: spproj
type(ori)        :: orientation
character(len=4096) :: argument, key, value, projfile, snapshot
integer :: iarg, separator, iptcl, nptcls, nactive, ninvalid, nnonfinite
integer :: state, proj, inpl, snapshot_unit
real :: euler(3), shift(2), corr
real(kind=8) :: corr_sum, shift_sq_sum
logical :: project_exists

projfile = ''
snapshot = 'continuous_inplane_refine3D_legacy_baseline.tsv'

do iarg = 1, command_argument_count()
    call get_command_argument(iarg, argument)
    separator = index(argument, '=')
    if( separator <= 1 ) cycle
    key   = adjustl(argument(:separator-1))
    value = adjustl(argument(separator+1:))
    select case(trim(key))
    case('projfile')
        projfile = trim(value)
    case('snapshot')
        snapshot = trim(value)
    case('case')
        continue
    case default
        write(logfhandle,'(a)') 'WARNING: ignored argument '//trim(argument)
    end select
enddo

if( len_trim(projfile) == 0 )then
    write(logfhandle,'(a)') &
        'simple_test_continuous_inplane_refine3D projfile=xx [snapshot=xx.tsv] [case=baseline]'
    stop 2
endif
inquire(file=trim(projfile), exist=project_exists)
if( .not. project_exists ) THROW_HARD('refine3D baseline project does not exist')

call spproj%read(string(trim(projfile)))
nptcls = spproj%os_ptcl3D%get_noris()
if( nptcls < 1 ) THROW_HARD('refine3D baseline project contains no ptcl3D metadata')

open(newunit=snapshot_unit, file=trim(snapshot), status='replace', action='write')
write(snapshot_unit,'(a)') 'iptcl state proj inpl e1 e2 e3 sx sy corr'

nactive      = 0
ninvalid     = 0
nnonfinite   = 0
corr_sum     = 0.d0
shift_sq_sum = 0.d0
do iptcl = 1, nptcls
    state = nint(spproj%os_ptcl3D%get(iptcl, 'state'))
    if( state <= 0 ) cycle
    nactive = nactive + 1
    proj = nint(spproj%os_ptcl3D%get(iptcl, 'proj'))
    inpl = nint(spproj%os_ptcl3D%get(iptcl, 'inpl'))
    corr = spproj%os_ptcl3D%get(iptcl, 'corr')
    call spproj%os_ptcl3D%get_ori(iptcl, orientation)
    euler = orientation%get_euler()
    shift = orientation%get_2Dshift()

    if( state < 1 .or. proj < 1 .or. inpl < 1 ) ninvalid = ninvalid + 1
    if( .not. all(ieee_is_finite(euler)) .or. &
        &.not. all(ieee_is_finite(shift)) .or. .not. ieee_is_finite(corr) )then
        nnonfinite = nnonfinite + 1
    endif
    corr_sum = corr_sum + real(corr, kind=8)
    shift_sq_sum = shift_sq_sum + sum(real(shift, kind=8)**2)
    write(snapshot_unit,'(4(i0,1x),6(es24.16,1x))') &
        iptcl, state, proj, inpl, euler, shift, corr
enddo
close(snapshot_unit)

write(logfhandle,'(a)') 'REFINE3D_LEGACY_BASELINE PROJECT: '//trim(projfile)
write(logfhandle,'(a)') 'REFINE3D_LEGACY_BASELINE SNAPSHOT: '//trim(snapshot)
write(logfhandle,'(a,4(i0,1x))') &
    'REFINE3D_LEGACY_BASELINE TOTAL/ACTIVE/INVALID/NONFINITE: ', &
    nptcls, nactive, ninvalid, nnonfinite
write(logfhandle,'(a,es24.16)') 'REFINE3D_LEGACY_BASELINE CORR_SUM: ', corr_sum
write(logfhandle,'(a,es24.16)') 'REFINE3D_LEGACY_BASELINE SHIFT_SQ_SUM: ', shift_sq_sum

if( nactive < 1 ) THROW_HARD('refine3D baseline project contains no active particles')
if( ninvalid /= 0 ) THROW_HARD('refine3D baseline contains invalid state/proj/inpl indices')
if( nnonfinite /= 0 ) THROW_HARD('refine3D baseline contains non-finite pose metadata')

call orientation%kill
call spproj%kill
call simple_end('SIMPLE_TEST_CONT_INPL_REFINE3D_BASELINE NORMAL STOP')
end subroutine run_legacy_baseline

end module continuous_inplane_refine3D_baseline_test
