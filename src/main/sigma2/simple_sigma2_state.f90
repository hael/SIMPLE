!@descr: canonical per-particle sigma2 state validation, reduction and transactional consolidation
! An update of the committed store (a matcher iteration) runs as one transaction: the workers read their
! rows of the committed store and each writes a local range file; the master then streams a fresh candidate
! (sigma2_state_assemble_candidate), each row from the range that covers it or else from the committed
! store, summing the active, valid rows per half as they pass, and publishes it (sigma2_state_publish).
! The merges of stores (sigma2_state_concatenate) stream the same way. sigma2_state_reduce is the full
! reduction of a store's records, which the tests hold the streamed global spectra against.
module simple_sigma2_state
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
use simple_defs,              only: logfhandle
use simple_string,            only: string
use simple_string_utils,      only: int2str, int2str_pad
use simple_syslib,            only: del_file
use simple_sigma2_state_file, only: sigma2_state_header, sigma2_state_read_header, &
    &sigma2_state_read_particles, sigma2_state_write_global, sigma2_state_validate_file, &
    &sigma2_state_open_candidate, sigma2_state_write_particles_unit, sigma2_state_close_candidate, &
    &sigma2_state_open_store, sigma2_state_read_particles_unit, sigma2_state_open_local_range, &
    &sigma2_state_read_local_range_rows, SIGMA2_PROV_RESIDUAL, SIGMA2_STATE_COMMITTED
implicit none
private

public :: sigma2_state_candidate_path, sigma2_state_range_path, sigma2_state_next_generation
public :: sigma2_state_prepare_update, sigma2_state_assemble_candidate, sigma2_state_concatenate
public :: sigma2_state_reduce, sigma2_state_validate_identity

integer, parameter :: BLOCK_ROWS = 4096
character(len=*), parameter :: SIGMA2_RANGE_FBODY = 'sigma2_state_range_part'

contains

    !> Transaction-scoped name of the candidate that will commit the given
    !! generation: <committed stem>.g<generation>.next. Candidate and range
    !! names carry the generation they belong to, so a file left by another
    !! transaction (a crashed run, a second run in the same directory) can
    !! never be mistaken for this one.
    function sigma2_state_candidate_path(committed_path, generation) result(candidate_path)
        character(len=*), intent(in) :: committed_path
        integer(int64),   intent(in) :: generation
        type(string) :: candidate_path
        candidate_path = transaction_stem(committed_path)//'.g'//int2str(int(generation))//'.next'
    end function sigma2_state_candidate_path

    !> Transaction-scoped name of one worker's range file for the update that
    !! commits the given generation: <committed stem>.g<generation>.part<NN>.range
    function sigma2_state_range_path(committed_path, generation, part, numlen) result(range_path)
        character(len=*), intent(in) :: committed_path
        integer(int64),   intent(in) :: generation
        integer,          intent(in) :: part, numlen
        type(string) :: range_path
        range_path = transaction_stem(committed_path)//'.g'//int2str(int(generation))//'.part'//&
            &int2str_pad(part,max(1,numlen))//'.range'
    end function sigma2_state_range_path

    !> The generation the next update of the committed state will commit
    subroutine sigma2_state_next_generation(committed_path, generation, status, message)
        character(len=*), intent(in)  :: committed_path
        integer(int64),   intent(out) :: generation
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(sigma2_state_header) :: header
        generation = 0_int64
        call sigma2_state_read_header(committed_path, header, status, message)
        if( status /= 0 ) return
        if( header%state /= SIGMA2_STATE_COMMITTED )then
            status = 1; message = 'canonical sigma2 update source is not committed'; return
        endif
        generation = header%generation + 1_int64
    end subroutine sigma2_state_next_generation

    !> the committed path without its .bin extension
    function transaction_stem(committed_path) result(stem)
        character(len=*), intent(in) :: committed_path
        type(string) :: stem
        integer :: l
        l = len_trim(committed_path)
        if( l > 4 )then
            if( committed_path(l-3:l) == '.bin' )then
                stem = committed_path(1:l-4)
                return
            endif
        endif
        stem = committed_path(1:l)
    end function transaction_stem

    !> An active particle record enters the global spectra only when every
    !! shell is finite and positive. Invalid records are skipped and counted
    !! in one warning.
    pure logical function record_is_valid(spectrum)
        real(real32), intent(in) :: spectrum(:)
        record_is_valid = all(ieee_is_finite(spectrum)) .and. all(spectrum > 0.0_real32)
    end function record_is_valid

    !> adds the valid records of the active particles among block(shells, rows first_row:) to their
    !! half's sums and counts; the active records skipped as invalid are counted in ninvalid
    subroutine add_records(block, first_row, active, eo, sums, counts, ninvalid, status, message)
        real(real32),     intent(in)    :: block(:,:)
        integer,          intent(in)    :: first_row
        logical,          intent(in)    :: active(:)
        integer,          intent(in)    :: eo(:)
        real(real64),     intent(inout) :: sums(:,:)
        integer(int64),   intent(inout) :: counts(2), ninvalid
        integer,          intent(out)   :: status
        character(len=*), intent(out)   :: message
        integer :: i, row, half
        status  = 0
        message = ''
        do i = 1, size(block,2)
            row = first_row+i-1
            if( .not. active(row) ) cycle
            if( .not. record_is_valid(block(:,i)) )then
                ninvalid = ninvalid + 1_int64
                cycle
            endif
            if( eo(row) < 0 .or. eo(row) > 1 )then
                status = 1; message = 'active particle has invalid even/odd assignment'; return
            endif
            half = eo(row)+1
            sums(:,half) = sums(:,half) + real(block(:,i),real64)
            counts(half) = counts(half) + 1_int64
        enddo
    end subroutine add_records

    !> the two global spectra (shells, half): per half, the mean of the records summed; an empty half or
    !! a non-finite mean is refused, and the skipped invalid records are reported in one warning
    subroutine form_global(sums, counts, ninvalid, global, status, message)
        real(real64),              intent(in)  :: sums(:,:)
        integer(int64),            intent(in)  :: counts(2), ninvalid
        real(real32), allocatable, intent(out) :: global(:,:)
        integer,                   intent(out) :: status
        character(len=*),          intent(out) :: message
        integer :: half
        status  = 0
        message = ''
        if( ninvalid > 0_int64 ) write(logfhandle,'(A,I0,A)') &
            &'>>> WARNING: canonical sigma2 reduction skipped ', ninvalid, &
            &' active particle record(s) with non-finite or non-positive values'
        if( any(counts == 0_int64) )then
            status = 1; message = 'canonical sigma2 has an empty even/odd half'; return
        endif
        allocate(global(size(sums,1),2))
        do half = 1, 2
            global(:,half) = real(sums(:,half)/real(counts(half),real64), real32)
        enddo
        if( .not. all(ieee_is_finite(global)) )then
            status = 1; message = 'canonical global sigma2 is not finite'
        endif
    end subroutine form_global

    !> Opens the update that commits the next generation, before the workers run: the committed store is
    !! validated, and a candidate or range files of that generation left by another transaction (a crashed
    !! run) are removed. Nothing is copied: the commit streams a fresh candidate.
    subroutine sigma2_state_prepare_update(committed_path, nparts, numlen, status, message)
        character(len=*), intent(in)  :: committed_path
        integer,          intent(in)  :: nparts, numlen
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        integer(int64) :: generation
        integer :: ipart
        call sigma2_state_validate_file(committed_path, status, message)
        if( status /= 0 ) return
        call sigma2_state_next_generation(committed_path, generation, status, message)
        if( status /= 0 ) return
        call del_file(sigma2_state_candidate_path(committed_path, generation))
        do ipart = 1, nparts
            call del_file(sigma2_state_range_path(committed_path, generation, ipart, numlen))
        enddo
    end subroutine sigma2_state_prepare_update

    !> The streaming commit's candidate of the next generation: a fresh file of the committed store's
    !! identity filled in row blocks through one handle, each row from the local range that covers it or
    !! else from the committed store. The active, valid rows are summed per half as they pass and the two
    !! global spectra are their means; the candidate is synced once. The ranges must belong to this
    !! transaction, lie in the schedule and cover it exactly, one range per row. The range files are left
    !! in place: the caller deletes them after publication, so a failed commit can be rerun.
    subroutine sigma2_state_assemble_candidate(committed_path, candidate_path, range_paths, scheduled_rows, &
        &active, eo, status, message)
        character(len=*), intent(in)  :: committed_path, candidate_path
        type(string),     intent(in)  :: range_paths(:)
        logical,          intent(in)  :: scheduled_rows(:), active(:)
        integer,          intent(in)  :: eo(:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(sigma2_state_header) :: committed, header
        real(real32), allocatable :: block(:,:), global(:,:)
        real(real64), allocatable :: sums(:,:)
        integer,      allocatable :: owner(:), range_first(:)
        integer(int64) :: counts(2), ninvalid, generation, layout_digest
        integer :: nrows, nshell, i, first_row, last_row, kfrom, kto, row, seg_last, src, nblock
        integer :: src_unit, range_unit, dst_unit, open_range, first_row_r, last_row_r
        src_unit   = -1
        range_unit = -1
        dst_unit   = -1
        call sigma2_state_validate_file(committed_path, status, message)
        if( status /= 0 ) return
        call sigma2_state_read_header(committed_path, committed, status, message)
        if( status /= 0 ) return
        if( committed%state /= SIGMA2_STATE_COMMITTED )then
            status = 1; message = 'canonical sigma2 update source is not committed'; return
        endif
        nrows = int(committed%nptcls)
        if( size(scheduled_rows) /= nrows )then
            status = 1; message = 'sigma2 scheduled-row mask has the wrong size'; return
        endif
        if( size(active) /= nrows .or. size(eo) /= nrows )then
            status = 1; message = 'sigma2 even/odd metadata has the wrong size'; return
        endif
        header            = committed
        header%generation = committed%generation + 1_int64
        header%provenance = SIGMA2_PROV_RESIDUAL
        ! the range each row comes from (0: the committed store)
        allocate(owner(nrows), source=0)
        allocate(range_first(size(range_paths)))
        do i = 1, size(range_paths)
            call sigma2_state_open_local_range(range_paths(i)%to_char(), range_unit, generation, layout_digest, &
                &first_row, last_row, kfrom, kto, status, message)
            if( status /= 0 ) return
            close(range_unit)
            range_unit = -1
            if( generation /= header%generation .or. layout_digest /= header%layout_digest )then
                status = 1; message = 'local sigma2 range belongs to another generation'; return
            endif
            if( kfrom /= header%kfrom .or. kto /= header%kto )then
                status = 1; message = 'local sigma2 range has incompatible shell bounds'; return
            endif
            if( first_row < 1 .or. last_row > nrows )then
                status = 1; message = 'local sigma2 range lies outside the project'; return
            endif
            if( any(owner(first_row:last_row) /= 0) )then
                status = 1; message = 'overlapping local sigma2 ranges'; return
            endif
            if( .not. all(scheduled_rows(first_row:last_row)) )then
                status = 1; message = 'local sigma2 range contains unscheduled rows'; return
            endif
            owner(first_row:last_row) = i
            range_first(i)            = first_row
        enddo
        if( any(scheduled_rows .and. owner == 0) )then
            status = 1; message = 'local sigma2 ranges do not exactly cover the schedule'; return
        endif
        if( any(owner == 0) )then
            call sigma2_state_open_store(committed_path, committed, src_unit, status, message)
            if( status /= 0 ) return
        endif
        call sigma2_state_open_candidate(candidate_path, header, dst_unit, status, message)
        if( status /= 0 )then
            call close_units
            return
        endif
        nshell = int(header%kto-header%kfrom+1)
        allocate(sums(nshell,2), source=0.0_real64)
        allocate(block(nshell,min(nrows,BLOCK_ROWS)))
        counts     = 0_int64
        ninvalid   = 0_int64
        open_range = 0
        do first_row = 1, nrows, BLOCK_ROWS
            last_row = min(nrows, first_row+BLOCK_ROWS-1)
            nblock   = last_row-first_row+1
            ! the block in runs of rows from one source
            row = first_row
            do while( row <= last_row )
                src      = owner(row)
                seg_last = row
                do while( seg_last < last_row )
                    if( owner(seg_last+1) /= src ) exit
                    seg_last = seg_last+1
                enddo
                if( src == 0 )then
                    call sigma2_state_read_particles_unit(src_unit, committed, row, &
                        &block(:,row-first_row+1:seg_last-first_row+1), status, message)
                else
                    if( src /= open_range )then
                        if( range_unit /= -1 ) close(range_unit)
                        call sigma2_state_open_local_range(range_paths(src)%to_char(), range_unit, generation, &
                            &layout_digest, first_row_r, last_row_r, kfrom, kto, status, message)
                        open_range = src
                    endif
                    if( status == 0 ) call sigma2_state_read_local_range_rows(range_unit, range_first(src), row, &
                        &block(:,row-first_row+1:seg_last-first_row+1), status, message)
                endif
                if( status /= 0 )then
                    call close_units
                    return
                endif
                row = seg_last+1
            enddo
            call add_records(block(:,:nblock), first_row, active, eo, sums, counts, ninvalid, status, message)
            if( status == 0 ) call sigma2_state_write_particles_unit(dst_unit, header, first_row, &
                &block(:,:nblock), status, message)
            if( status /= 0 )then
                call close_units
                return
            endif
        enddo
        if( src_unit   /= -1 ) close(src_unit)
        if( range_unit /= -1 ) close(range_unit)
        call form_global(sums, counts, ninvalid, global, status, message)
        if( status /= 0 )then
            close(dst_unit)
            return
        endif
        call sigma2_state_close_candidate(candidate_path, dst_unit, header, global, status, message)

      contains

        subroutine close_units
            if( src_unit   /= -1 ) close(src_unit)
            if( range_unit /= -1 ) close(range_unit)
            if( dst_unit   /= -1 ) close(dst_unit)
        end subroutine close_units

    end subroutine sigma2_state_assemble_candidate

    !> A fresh candidate of identity header whose rows are the committed stores source_paths one after
    !! another, copied in row blocks through one handle each. The active, valid rows are summed per half
    !! as they pass and the two global spectra are their means; the candidate is synced once. The
    !! sources must have the candidate's shells and together its row count.
    subroutine sigma2_state_concatenate(source_paths, candidate_path, header, active, eo, status, message)
        type(string),              intent(in)    :: source_paths(:)
        character(len=*),          intent(in)    :: candidate_path
        type(sigma2_state_header), intent(inout) :: header
        logical,                   intent(in)    :: active(:)
        integer,                   intent(in)    :: eo(:)
        integer,                   intent(out)   :: status
        character(len=*),          intent(out)   :: message
        type(sigma2_state_header) :: source
        real(real32), allocatable :: block(:,:), global(:,:)
        real(real64), allocatable :: sums(:,:)
        integer(int64) :: counts(2), ninvalid
        integer :: isrc, offset, first_row, nblock, nshell, src_unit, dst_unit
        src_unit = -1
        if( size(active) /= header%nptcls .or. size(eo) /= header%nptcls )then
            status = 1; message = 'sigma2 even/odd metadata has the wrong size'; return
        endif
        call sigma2_state_open_candidate(candidate_path, header, dst_unit, status, message)
        if( status /= 0 ) return
        nshell = int(header%kto-header%kfrom+1)
        allocate(sums(nshell,2), source=0.0_real64)
        allocate(block(nshell,BLOCK_ROWS))
        counts   = 0_int64
        ninvalid = 0_int64
        offset   = 0
        do isrc = 1, size(source_paths)
            call sigma2_state_open_store(source_paths(isrc)%to_char(), source, src_unit, status, message)
            if( status /= 0 )then
                call close_units
                return
            endif
            if( source%kfrom /= header%kfrom .or. source%kto /= header%kto )then
                status = 1; message = 'concatenated sigma2 stores have incompatible shell bounds'
            else if( offset + int(source%nptcls) > header%nptcls )then
                status = 1; message = 'concatenated sigma2 stores hold more rows than the candidate'
            endif
            if( status /= 0 )then
                call close_units
                return
            endif
            do first_row = 1, int(source%nptcls), BLOCK_ROWS
                nblock = min(int(source%nptcls), first_row+BLOCK_ROWS-1) - first_row + 1
                call sigma2_state_read_particles_unit(src_unit, source, first_row, block(:,:nblock), status, message)
                if( status == 0 ) call add_records(block(:,:nblock), offset+first_row, active, eo, sums, counts, &
                    &ninvalid, status, message)
                if( status == 0 ) call sigma2_state_write_particles_unit(dst_unit, header, offset+first_row, &
                    &block(:,:nblock), status, message)
                if( status /= 0 )then
                    call close_units
                    return
                endif
            enddo
            close(src_unit)
            src_unit = -1
            offset   = offset + int(source%nptcls)
        enddo
        if( offset /= header%nptcls )then
            status = 1; message = 'concatenated sigma2 stores do not cover the candidate rows'
            call close_units
            return
        endif
        call form_global(sums, counts, ninvalid, global, status, message)
        if( status /= 0 )then
            call close_units
            return
        endif
        call sigma2_state_close_candidate(candidate_path, dst_unit, header, global, status, message)

      contains

        subroutine close_units
            if( src_unit /= -1 ) close(src_unit)
            close(dst_unit)
        end subroutine close_units

    end subroutine sigma2_state_concatenate

    !> The full reduction of the candidate at path, its global spectra written: per half, the mean of
    !! the valid records of the active particles of that half, every record re-read from the file
    subroutine sigma2_state_reduce(path, active, eo, status, message)
        character(len=*), intent(in) :: path
        logical,          intent(in) :: active(:)
        integer,          intent(in) :: eo(:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(sigma2_state_header) :: header
        real(real32), allocatable :: spectra(:,:), global(:,:)
        real(real64), allocatable :: sums(:,:)
        integer(int64) :: counts(2), ninvalid
        integer :: first_row, last_row, nshell
        call sigma2_state_read_header(path, header, status, message)
        if( status /= 0 ) return
        if( size(active) /= header%nptcls .or. size(eo) /= header%nptcls )then
            status = 1; message = 'sigma2 even/odd metadata has the wrong size'; return
        endif
        nshell = int(header%kto-header%kfrom+1)
        allocate(sums(nshell,2), source=0.0_real64)
        counts   = 0_int64
        ninvalid = 0_int64
        do first_row = 1, int(header%nptcls), BLOCK_ROWS
            last_row = min(int(header%nptcls), first_row+BLOCK_ROWS-1)
            call sigma2_state_read_particles(path, first_row, last_row, spectra, status, message)
            if( status /= 0 ) return
            call add_records(spectra, first_row, active, eo, sums, counts, ninvalid, status, message)
            if( status /= 0 ) return
            deallocate(spectra)
        enddo
        call form_global(sums, counts, ninvalid, global, status, message)
        if( status /= 0 ) return
        call sigma2_state_write_global(path, global, status, message)
    end subroutine sigma2_state_reduce

    subroutine sigma2_state_validate_identity(path, box, smpd, kfrom, kto, nptcls, layout_digest, &
        &status, message, expected_state)
        character(len=*), intent(in) :: path
        integer,          intent(in) :: box, kfrom, kto, nptcls
        real,             intent(in) :: smpd
        integer(int64),   intent(in) :: layout_digest
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        integer(int32), optional, intent(in) :: expected_state
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
    end subroutine sigma2_state_validate_identity

end module simple_sigma2_state
