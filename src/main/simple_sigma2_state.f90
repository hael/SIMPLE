!@descr: canonical per-particle sigma2 state validation, reduction and transactional consolidation
module simple_sigma2_state
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
use simple_fileio,            only: get_fpath
use simple_oris,              only: oris
use simple_sp_project,        only: sp_project
use simple_string,            only: string
use simple_string_utils,      only: int2str_pad
use simple_sigma2_state_file, only: sigma2_state_header, sigma2_state_read_header, &
    &sigma2_state_read_particles, sigma2_state_read_groups, sigma2_state_write_particles, &
    &sigma2_state_write_groups, sigma2_state_read_local_range, sigma2_state_refresh_integrity, &
    &sigma2_state_validate_file, sigma2_state_publish, sigma2_state_digest_begin, &
    &sigma2_state_digest_text, sigma2_state_digest_integer, SIGMA2_GROUP_GLOBAL, SIGMA2_GROUP_STACK, &
    &SIGMA2_STATE_NEXT_FNAME, SIGMA2_PROV_RESIDUAL, SIGMA2_STATE_CANDIDATE, &
    &SIGMA2_STATE_COMMITTED, sigma2_state_create_candidate
implicit none
private

public :: sigma2_state_layout_digest
public :: sigma2_state_project_layout_digest
public :: sigma2_state_candidate_path, sigma2_state_range_path, sigma2_state_prepare_update
public :: sigma2_state_merge_local_ranges, sigma2_state_reduce_groups
public :: sigma2_state_validate_identity, sigma2_state_validate_science, sigma2_state_commit

integer, parameter :: REDUCE_BLOCK_ROWS = 4096
character(len=*), parameter :: SIGMA2_RANGE_FBODY = 'sigma2_state_range_part'

contains

    function sigma2_state_layout_digest(lineage, stack_refs, stack_ids, stack_indices, nrows) result(digest)
        character(len=*), intent(in) :: lineage
        type(string),     intent(in) :: stack_refs(:)
        integer,          intent(in) :: stack_ids(:), stack_indices(:)
        integer, optional, intent(in) :: nrows
        integer(int64) :: digest
        integer :: i, n
        n = size(stack_ids)
        if( present(nrows) ) n = nrows
        if( n < 0 .or. n > size(stack_ids) .or. n > size(stack_indices) )then
            digest = 0_int64
            return
        endif
        digest = sigma2_state_digest_begin()
        call sigma2_state_digest_text(digest, trim(lineage))
        do i = 1, n
            if( stack_ids(i) < 1 .or. stack_ids(i) > size(stack_refs) )then
                digest = 0_int64
                return
            endif
            call sigma2_state_digest_text(digest, trim(stack_refs(stack_ids(i))%to_char()))
            call sigma2_state_digest_integer(digest, stack_indices(i))
        enddo
        if( digest == 0_int64 ) digest = 1_int64
    end function sigma2_state_layout_digest

    function sigma2_state_project_layout_digest(project, particles, nrows) result(digest)
        type(sp_project), intent(in) :: project
        class(oris),      intent(in) :: particles
        integer, optional, intent(in) :: nrows
        integer(int64) :: digest
        type(string), allocatable :: stack_refs(:)
        integer, allocatable :: stack_ids(:), stack_indices(:)
        type(string) :: lineage, stack_ref
        integer :: i, nptcls, nstks
        digest = 0_int64
        nptcls = particles%get_noris(consider_state=.false.)
        nstks  = project%os_stk%get_noris(consider_state=.false.)
        if( nptcls < 1 .or. nstks < 1 ) return
        if( project%projinfo%get_noris() /= 1 ) return
        if( project%projinfo%isthere(1, 'projname') )then
            lineage = project%projinfo%get_str(1, 'projname')
        else if( project%projinfo%isthere(1, 'projfile') )then
            lineage = project%projinfo%get_str(1, 'projfile')
        else
            return
        endif
        allocate(stack_refs(nstks), stack_ids(nptcls), stack_indices(nptcls))
        do i = 1, nstks
            stack_ref = project%os_stk%get_str(i, 'stk')
            stack_refs(i) = trim(adjustl(stack_ref%to_char()))
        enddo
        do i = 1, nptcls
            stack_ids(i)     = particles%get_int(i, 'stkind')
            stack_indices(i) = particles%get_int(i, 'indstk')
        enddo
        if( present(nrows) )then
            digest = sigma2_state_layout_digest(lineage%to_char(), stack_refs, stack_ids, stack_indices, nrows)
        else
            digest = sigma2_state_layout_digest(lineage%to_char(), stack_refs, stack_ids, stack_indices)
        endif
        call lineage%kill
        call stack_ref%kill
        call stack_refs(:)%kill
        deallocate(stack_refs, stack_ids, stack_indices)
    end function sigma2_state_project_layout_digest

    function sigma2_state_candidate_path(committed_path) result(candidate_path)
        character(len=*), intent(in) :: committed_path
        type(string) :: candidate_path
        type(string) :: parent
        parent = get_fpath(string(committed_path))
        candidate_path = parent//SIGMA2_STATE_NEXT_FNAME
        call parent%kill
    end function sigma2_state_candidate_path

    function sigma2_state_range_path(committed_path, part, numlen) result(range_path)
        character(len=*), intent(in) :: committed_path
        integer,          intent(in) :: part, numlen
        type(string) :: range_path
        type(string) :: parent
        parent = get_fpath(string(committed_path))
        range_path = parent//SIGMA2_RANGE_FBODY//int2str_pad(part,max(1,numlen))//'.bin'
        call parent%kill
    end function sigma2_state_range_path

    subroutine sigma2_state_prepare_update(committed_path, candidate_path, status, message)
        character(len=*), intent(in) :: committed_path, candidate_path
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(sigma2_state_header) :: header
        call sigma2_state_validate_file(committed_path, status, message, deep=.true.)
        if( status /= 0 ) return
        call sigma2_state_read_header(committed_path, header, status, message)
        if( status /= 0 ) return
        if( header%state /= SIGMA2_STATE_COMMITTED )then
            status = 1; message = 'canonical sigma2 update source is not committed'; return
        endif
        header%generation = header%generation + 1_int64
        header%provenance = SIGMA2_PROV_RESIDUAL
        call sigma2_state_create_candidate(candidate_path, header, status, message, source_path=committed_path)
    end subroutine sigma2_state_prepare_update

    subroutine sigma2_state_merge_local_ranges(candidate_path, range_paths, scheduled_rows, status, message)
        character(len=*), intent(in) :: candidate_path
        type(string),     intent(in) :: range_paths(:)
        logical,          intent(in) :: scheduled_rows(:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(sigma2_state_header) :: header
        real(real32), allocatable :: spectra(:,:)
        logical, allocatable :: covered(:)
        integer(int64) :: generation, layout_digest
        integer :: i, first_row, last_row, kfrom, kto
        call sigma2_state_read_header(candidate_path, header, status, message)
        if( status /= 0 ) return
        if( header%state /= SIGMA2_STATE_CANDIDATE )then
            status = 1; message = 'local ranges can only merge into a sigma2 candidate'; return
        endif
        if( size(scheduled_rows) /= header%nptcls )then
            status = 1; message = 'sigma2 scheduled-row mask has the wrong size'; return
        endif
        allocate(covered(header%nptcls), source=.false.)
        do i = 1, size(range_paths)
            call sigma2_state_read_local_range(range_paths(i)%to_char(), generation, layout_digest, &
                &first_row, last_row, kfrom, kto, spectra, status, message)
            if( status /= 0 ) return
            if( generation /= header%generation .or. layout_digest /= header%layout_digest )then
                status = 1; message = 'local sigma2 range belongs to another generation'; return
            endif
            if( kfrom /= header%kfrom .or. kto /= header%kto )then
                status = 1; message = 'local sigma2 range has incompatible shell bounds'; return
            endif
            if( first_row < 1 .or. last_row > header%nptcls )then
                status = 1; message = 'local sigma2 range lies outside the project'; return
            endif
            if( any(covered(first_row:last_row)) )then
                status = 1; message = 'overlapping local sigma2 ranges'; return
            endif
            if( .not. all(scheduled_rows(first_row:last_row)) )then
                status = 1; message = 'local sigma2 range contains unscheduled rows'; return
            endif
            call sigma2_state_write_particles(candidate_path, first_row, spectra, status, message)
            if( status /= 0 ) return
            covered(first_row:last_row) = .true.
            deallocate(spectra)
        enddo
        if( any(covered .neqv. scheduled_rows) )then
            status = 1; message = 'local sigma2 ranges do not exactly cover the schedule'; return
        endif
        call sigma2_state_refresh_integrity(candidate_path, status, message)
    end subroutine sigma2_state_merge_local_ranges

    subroutine sigma2_state_reduce_groups(path, active, eo, group_ids, status, message)
        character(len=*), intent(in) :: path
        logical,          intent(in) :: active(:)
        integer,          intent(in) :: eo(:), group_ids(:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(sigma2_state_header) :: header
        real(real32), allocatable :: spectra(:,:), groups(:,:,:)
        integer(int64), allocatable :: checksums(:), counts(:,:)
        real(real64), allocatable :: sums(:,:,:)
        integer :: first_row, last_row, i, row, half, group, nshell
        call sigma2_state_read_header(path, header, status, message)
        if( status /= 0 ) return
        if( size(active) /= header%nptcls .or. size(eo) /= header%nptcls .or. &
            &size(group_ids) /= header%nptcls )then
            status = 1; message = 'sigma2 grouping metadata has the wrong size'; return
        endif
        nshell = int(header%kto-header%kfrom+1)
        allocate(sums(nshell,2,header%ngroups), source=0.0_real64)
        allocate(counts(2,header%ngroups), source=0_int64)
        do first_row = 1, int(header%nptcls), REDUCE_BLOCK_ROWS
            last_row = min(int(header%nptcls), first_row+REDUCE_BLOCK_ROWS-1)
            call sigma2_state_read_particles(path, first_row, last_row, spectra, status, message, checksums)
            if( status /= 0 ) return
            do i = 1, size(spectra,2)
                row = first_row+i-1
                if( .not. active(row) ) cycle
                if( checksums(i) == 0_int64 )then
                    status = 1; message = 'active particle has no canonical sigma2 record'; return
                endif
                if( any(.not. ieee_is_finite(spectra(:,i))) .or. any(spectra(:,i) <= 0.0_real32) )then
                    status = 1; message = 'active particle has invalid canonical sigma2 values'; return
                endif
                if( eo(row) < 0 .or. eo(row) > 1 )then
                    status = 1; message = 'active particle has invalid even/odd assignment'; return
                endif
                half = eo(row)+1
                select case(header%grouping)
                    case(SIGMA2_GROUP_GLOBAL)
                        group = 1
                    case(SIGMA2_GROUP_STACK)
                        group = group_ids(row)
                    case default
                        status = 1; message = 'unsupported canonical sigma2 grouping policy'; return
                end select
                if( group < 1 .or. group > header%ngroups )then
                    status = 1; message = 'active particle has invalid sigma2 group'; return
                endif
                sums(:,half,group) = sums(:,half,group) + real(spectra(:,i),real64)
                counts(half,group) = counts(half,group) + 1_int64
            enddo
            deallocate(spectra, checksums)
        enddo
        if( any(counts == 0_int64) )then
            status = 1; message = 'canonical sigma2 group has an empty even/odd half'; return
        endif
        allocate(groups(nshell,2,header%ngroups))
        do group = 1, int(header%ngroups)
            do half = 1, 2
                groups(:,half,group) = real(sums(:,half,group)/real(counts(half,group),real64), real32)
            enddo
        enddo
        call sigma2_state_write_groups(path, groups, status, message)
        if( status == 0 ) call sigma2_state_refresh_integrity(path, status, message)
    end subroutine sigma2_state_reduce_groups

    subroutine sigma2_state_validate_identity(path, box, smpd, kfrom, kto, nptcls, layout_digest, &
        &status, message, expected_state, expected_grouping, expected_ngroups)
        character(len=*), intent(in) :: path
        integer,          intent(in) :: box, kfrom, kto, nptcls
        real,             intent(in) :: smpd
        integer(int64),   intent(in) :: layout_digest
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        integer(int32), optional, intent(in) :: expected_state
        integer(int32), optional, intent(in) :: expected_grouping
        integer, optional, intent(in) :: expected_ngroups
        type(sigma2_state_header) :: header
        call sigma2_state_validate_file(path, status, message)
        if( status /= 0 ) return
        call sigma2_state_read_header(path, header, status, message)
        if( status /= 0 ) return
        if( header%box /= box .or. header%kfrom /= kfrom .or. header%kto /= kto .or. &
            &header%nptcls /= nptcls .or. header%layout_digest /= layout_digest .or. &
            &abs(header%smpd-real(smpd,real64)) > 1.0e-9_real64 )then
            status = 1; message = 'canonical sigma2 identity does not match the particle project'; return
        endif
        if( present(expected_state) )then
            if( header%state /= expected_state )then
                status = 1; message = 'canonical sigma2 transaction state does not match'; return
            endif
        endif
        if( present(expected_grouping) )then
            if( header%grouping /= expected_grouping )then
                status = 1; message = 'canonical sigma2 grouping policy does not match'; return
            endif
        endif
        if( present(expected_ngroups) )then
            if( header%ngroups /= expected_ngroups )then
                status = 1; message = 'canonical sigma2 group count does not match'; return
            endif
        endif
    end subroutine sigma2_state_validate_identity

    subroutine sigma2_state_validate_science(path, active, eo, group_ids, status, message)
        character(len=*), intent(in) :: path
        logical,          intent(in) :: active(:)
        integer,          intent(in) :: eo(:), group_ids(:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(sigma2_state_header) :: header
        real(real32), allocatable :: stored(:,:,:), spectra(:,:)
        integer(int64), allocatable :: checksums(:), counts(:,:)
        real(real64), allocatable :: sums(:,:,:)
        integer :: first_row, last_row, i, row, half, group, nshell
        real(real64) :: expected, tolerance
        call sigma2_state_validate_file(path, status, message, deep=.true.)
        if( status /= 0 ) return
        call sigma2_state_read_header(path, header, status, message)
        if( status /= 0 ) return
        if( size(active) /= header%nptcls .or. size(eo) /= header%nptcls .or. &
            &size(group_ids) /= header%nptcls )then
            status = 1; message = 'sigma2 validation metadata has the wrong size'; return
        endif
        call sigma2_state_read_groups(path, stored, status, message)
        if( status /= 0 ) return
        nshell = int(header%kto-header%kfrom+1)
        allocate(sums(nshell,2,header%ngroups), source=0.0_real64)
        allocate(counts(2,header%ngroups), source=0_int64)
        do first_row = 1, int(header%nptcls), REDUCE_BLOCK_ROWS
            last_row = min(int(header%nptcls), first_row+REDUCE_BLOCK_ROWS-1)
            call sigma2_state_read_particles(path, first_row, last_row, spectra, status, message, checksums)
            if( status /= 0 ) return
            do i = 1, size(spectra,2)
                row = first_row+i-1
                if( .not. active(row) ) cycle
                if( checksums(i) == 0_int64 .or. any(.not. ieee_is_finite(spectra(:,i))) .or. &
                    &any(spectra(:,i) <= 0.0_real32) )then
                    status = 1; message = 'active particle has invalid canonical sigma2 state'; return
                endif
                if( eo(row) < 0 .or. eo(row) > 1 )then
                    status = 1; message = 'active particle has invalid even/odd assignment'; return
                endif
                half = eo(row)+1
                if( header%grouping == SIGMA2_GROUP_GLOBAL )then
                    group = 1
                else
                    group = group_ids(row)
                endif
                if( group < 1 .or. group > header%ngroups )then
                    status = 1; message = 'active particle has invalid canonical sigma2 group'; return
                endif
                sums(:,half,group) = sums(:,half,group) + real(spectra(:,i),real64)
                counts(half,group) = counts(half,group) + 1_int64
            enddo
            deallocate(spectra, checksums)
        enddo
        if( any(counts == 0_int64) )then
            status = 1; message = 'canonical sigma2 group has an empty even/odd half'; return
        endif
        do group = 1, int(header%ngroups)
            do half = 1, 2
                do i = 1, nshell
                    expected = sums(i,half,group)/real(counts(half,group),real64)
                    tolerance = 8.0_real64*epsilon(1.0_real32)*max(1.0_real64,abs(expected))
                    if( abs(real(stored(i,half,group),real64)-expected) > tolerance )then
                        status = 1; message = 'canonical grouped sigma2 does not match particle records'; return
                    endif
                enddo
            enddo
        enddo
        status = 0
        message = ''
    end subroutine sigma2_state_validate_science

    subroutine sigma2_state_commit(candidate_path, committed_path, active, eo, group_ids, status, message)
        character(len=*), intent(in) :: candidate_path, committed_path
        logical,          intent(in) :: active(:)
        integer,          intent(in) :: eo(:), group_ids(:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        call sigma2_state_validate_science(candidate_path, active, eo, group_ids, status, message)
        if( status /= 0 ) return
        call sigma2_state_publish(candidate_path, committed_path, status, message)
    end subroutine sigma2_state_commit

end module simple_sigma2_state
