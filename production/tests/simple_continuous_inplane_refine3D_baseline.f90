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
character(len=4096) :: argument, key, value, projfile, snapshot, expected_snapshot
character(len=4096) :: expected_header
integer :: iarg, separator, iptcl, nptcls, nactive, ninvalid, nnonfinite
integer :: state, proj, inpl, snapshot_unit, expected_unit, ios, ncompared
integer :: current_fields(4), expected_fields(4)
real :: euler(3), shift(2), corr
real(kind=8) :: corr_sum, shift_sq_sum, snapshot_atol, snapshot_rtol
real(kind=8) :: current_values(6), expected_values(6), value_tolerances(6)
logical :: project_exists, expected_exists, compare_snapshot

projfile = ''
snapshot = 'continuous_inplane_refine3D_legacy_baseline.tsv'
expected_snapshot = ''
snapshot_atol = 1.d-6
snapshot_rtol = 1.d-6

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
    case('expected_snapshot')
        expected_snapshot = trim(value)
    case('snapshot_atol')
        read(value,*) snapshot_atol
    case('snapshot_rtol')
        read(value,*) snapshot_rtol
    case('case')
        continue
    case default
        write(logfhandle,'(a)') 'WARNING: ignored argument '//trim(argument)
    end select
enddo

if( len_trim(projfile) == 0 )then
    write(logfhandle,'(a)') &
        'simple_test_continuous_inplane_refine3D projfile=xx [snapshot=xx.tsv] '// &
        &'[expected_snapshot=xx.tsv] [snapshot_atol=xx] [snapshot_rtol=xx] [case=baseline]'
    stop 2
endif
if( snapshot_atol < 0.d0 .or. snapshot_rtol < 0.d0 ) &
    &THROW_HARD('refine3D baseline snapshot tolerances must be non-negative')
inquire(file=trim(projfile), exist=project_exists)
if( .not. project_exists ) THROW_HARD('refine3D baseline project does not exist')

compare_snapshot = len_trim(expected_snapshot) > 0
if( compare_snapshot )then
    if( trim(snapshot) == trim(expected_snapshot) ) &
        &THROW_HARD('refine3D baseline output and expected snapshot must differ')
    inquire(file=trim(expected_snapshot), exist=expected_exists)
    if( .not. expected_exists ) THROW_HARD('expected refine3D baseline snapshot does not exist')
    open(newunit=expected_unit, file=trim(expected_snapshot), status='old', &
        &action='read', iostat=ios)
    if( ios /= 0 ) THROW_HARD('failed to open expected refine3D baseline snapshot')
    read(expected_unit,'(a)',iostat=ios) expected_header
    if( ios /= 0 .or. trim(expected_header) /= &
        &'iptcl state proj inpl e1 e2 e3 sx sy corr' )then
        THROW_HARD('expected refine3D baseline snapshot has an invalid header')
    endif
endif

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
ncompared    = 0
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
    if( compare_snapshot )then
        read(expected_unit,*,iostat=ios) expected_fields, expected_values
        if( ios /= 0 ) THROW_HARD('expected refine3D baseline snapshot ended before the current project')
        current_fields = [iptcl, state, proj, inpl]
        if( any(current_fields /= expected_fields) ) &
            &THROW_HARD('refine3D baseline integer metadata differs at particle '//int2str(iptcl))
        current_values = [real(euler,kind=8), real(shift,kind=8), real(corr,kind=8)]
        if( .not. all(ieee_is_finite(current_values)) ) &
            &THROW_HARD('current refine3D baseline has non-finite metadata at particle '//int2str(iptcl))
        if( .not. all(ieee_is_finite(expected_values)) ) &
            &THROW_HARD('expected refine3D baseline has non-finite metadata at particle '//int2str(iptcl))
        value_tolerances = snapshot_atol + snapshot_rtol * &
            &max(abs(current_values), abs(expected_values))
        if( any(abs(current_values-expected_values) > value_tolerances) ) &
            &THROW_HARD('refine3D baseline floating metadata differs at particle '//int2str(iptcl))
        ncompared = ncompared + 1
    endif
enddo
close(snapshot_unit)
if( compare_snapshot )then
    read(expected_unit,*,iostat=ios) expected_fields, expected_values
    if( ios == 0 ) THROW_HARD('expected refine3D baseline snapshot contains extra particle rows')
    if( ios > 0 ) THROW_HARD('expected refine3D baseline snapshot has a malformed trailing row')
    close(expected_unit)
    if( ncompared /= nactive ) THROW_HARD('refine3D baseline snapshot comparison count mismatch')
endif

write(logfhandle,'(a)') 'REFINE3D_LEGACY_BASELINE PROJECT: '//trim(projfile)
write(logfhandle,'(a)') 'REFINE3D_LEGACY_BASELINE SNAPSHOT: '//trim(snapshot)
if( compare_snapshot )then
    write(logfhandle,'(a)') 'REFINE3D_LEGACY_BASELINE MODE: COMPARE'
    write(logfhandle,'(a)') 'REFINE3D_LEGACY_BASELINE EXPECTED: '//trim(expected_snapshot)
    write(logfhandle,'(a,2(es12.4,1x))') &
        &'REFINE3D_LEGACY_BASELINE ATOL/RTOL: ', snapshot_atol, snapshot_rtol
else
    write(logfhandle,'(a)') 'REFINE3D_LEGACY_BASELINE MODE: GENERATE_ONLY'
    write(logfhandle,'(a)') &
        &'NOTICE: pass expected_snapshot=xx.tsv to perform a regression comparison.'
endif
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
