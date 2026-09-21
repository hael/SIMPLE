!@descr: flex per-state weight files: identity, science validation, transactions, delivery and loading
!!
!! Domain layer over simple_flex_weights_file (the bytes). Builder-free, like simple_sigma2_state:
!! callers pass the project and the particle field explicitly. The particle-layout identity is the
!! canonical sigma2 digest of the SAME project field (lineage + ordered stack reference and stack
!! index per physical row, `state` excluded), so the weight files and a sigma2 store of one project
!! agree on what "row i" means.
!!
!! One file per delivered state (`flex_weights_state_NNN.bin`), registered in the project's out
!! segment as imgkind `flex_weights` with the state index, beside `vol_flex` state NNN, so the
!! workflow's per-state selection and removal apply to the weights as they do to the maps.
!!
!! Producer path (flex_pca master, once per run): flex_weights_deliver scatters the run's
!! selection-indexed weight table onto the full layout, writes one candidate per state, validates
!! each file and the set, then publishes them. Consumer path: flex_weights_consumable +
!! flex_weights_load_state (one state) or flex_weights_load_all (every registered state).
!! The worker range protocol (write_local_range / merge_local_ranges) has no producer yet:
!! flex_pca is master-computes, workers-read.
module simple_flex_weights_state
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
use simple_oris,              only: oris
use simple_sp_project,        only: sp_project
use simple_string,            only: string
use simple_string_utils,      only: int2str, int2str_pad
use simple_syslib,            only: file_exists, del_file
use simple_sigma2_state,      only: sigma2_state_project_layout_digest
use simple_flex_weights_file, only: flex_weights_header, flex_weights_init_header, flex_weights_read_header, &
    &flex_weights_create_candidate, flex_weights_write_rows, flex_weights_read_rows, &
    &flex_weights_write_scalars, flex_weights_read_scalars, flex_weights_read_local_range, &
    &flex_weights_validate_file, flex_weights_publish, FLEX_WEIGHTS_FBODY, FLEX_WEIGHTS_EXT, &
    &FLEX_WEIGHTS_KIND_PARTITION, FLEX_WEIGHTS_KIND_KERNEL, FLEX_WEIGHTS_CANDIDATE, FLEX_WEIGHTS_COMMITTED, &
    &FLEX_WEIGHTS_SCALAR_MASS, FLEX_WEIGHTS_SCALAR_NEFF, FLEX_WEIGHTS_SCALAR_POP, FLEX_WEIGHTS_SCALAR_BW, &
    &FLEX_WEIGHTS_SCALAR_NFIELDS
implicit none
private

public :: flex_weights_state_fname, flex_weights_candidate_path, flex_weights_range_path
public :: flex_weights_next_generation, flex_weights_prepare_update, flex_weights_merge_local_ranges
public :: flex_weights_validate_identity, flex_weights_validate_science, flex_weights_validate_set
public :: flex_weights_commit, flex_weights_infer_kind, flex_weights_state_scalars
public :: flex_weights_write_store, flex_weights_deliver
public :: flex_weights_consumable, flex_weights_load_state, flex_weights_load_all
public :: FLEX_WEIGHTS_STALE_SCAN

integer,        parameter :: ROW_BLOCK = 262144
!> a PARTITION row sums to one across the files within this (responsibilities are renormalised in
!! double, then cast)
real(real64),   parameter :: ROW_SUM_TOL = 1.0e-3_real64
real(real64),   parameter :: MASS_TOL    = 1.0e-3_real64
!> how many state indices past the delivered count a delivery clears of stale files
integer,        parameter :: FLEX_WEIGHTS_STALE_SCAN = 64

contains

    ! ---- naming ----

    !> flex_weights_state_NNN.bin
    function flex_weights_state_fname(state) result(fname)
        integer, intent(in) :: state
        type(string) :: fname
        fname = FLEX_WEIGHTS_FBODY//int2str_pad(state,3)//FLEX_WEIGHTS_EXT
    end function flex_weights_state_fname

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

    !> <committed stem>.g<generation>.next
    function flex_weights_candidate_path(committed_path, generation) result(candidate_path)
        character(len=*), intent(in) :: committed_path
        integer(int64),   intent(in) :: generation
        type(string) :: candidate_path
        candidate_path = transaction_stem(committed_path)//'.g'//int2str(int(generation))//'.next'
    end function flex_weights_candidate_path

    !> <committed stem>.g<generation>.part<NN>.range
    function flex_weights_range_path(committed_path, generation, part, numlen) result(range_path)
        character(len=*), intent(in) :: committed_path
        integer(int64),   intent(in) :: generation
        integer,          intent(in) :: part, numlen
        type(string) :: range_path
        range_path = transaction_stem(committed_path)//'.g'//int2str(int(generation))//'.part'//&
            &int2str_pad(part,max(1,numlen))//'.range'
    end function flex_weights_range_path

    ! ---- transactions ----

    !> The generation the next update of a committed file will commit; 1 when there is none
    subroutine flex_weights_next_generation(committed_path, generation, status, message)
        character(len=*), intent(in)  :: committed_path
        integer(int64),   intent(out) :: generation
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(flex_weights_header) :: header
        generation = 1_int64
        status = 0
        message = ''
        call flex_weights_read_header(committed_path, header, status, message)
        if( status /= 0 ) return
        if( header%state /= FLEX_WEIGHTS_COMMITTED )then
            status = 1; message = 'flex weights update source is not committed'; return
        endif
        generation = header%generation + 1_int64
    end subroutine flex_weights_next_generation

    !> Candidate for the next generation, seeded with every committed row
    subroutine flex_weights_prepare_update(committed_path, candidate_path, status, message)
        character(len=*), intent(in)  :: committed_path, candidate_path
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(flex_weights_header) :: header
        call flex_weights_validate_file(committed_path, status, message, deep=.true.)
        if( status /= 0 ) return
        call flex_weights_read_header(committed_path, header, status, message)
        if( status /= 0 ) return
        if( header%state /= FLEX_WEIGHTS_COMMITTED )then
            status = 1; message = 'flex weights update source is not committed'; return
        endif
        header%generation = header%generation + 1_int64
        call flex_weights_create_candidate(candidate_path, header, status, message, source_path=committed_path)
    end subroutine flex_weights_prepare_update

    !> Exact-coverage merge of worker range files into one state's candidate
    subroutine flex_weights_merge_local_ranges(candidate_path, range_paths, scheduled_rows, status, message)
        character(len=*), intent(in)  :: candidate_path
        type(string),     intent(in)  :: range_paths(:)
        logical,          intent(in)  :: scheduled_rows(:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(flex_weights_header) :: header
        real(real32),   allocatable :: weights(:)
        integer(int32), allocatable :: flags(:)
        logical,        allocatable :: covered(:)
        integer(int64) :: generation, layout_digest
        integer :: i, state_index, first_row, last_row
        call flex_weights_read_header(candidate_path, header, status, message)
        if( status /= 0 ) return
        if( header%state /= FLEX_WEIGHTS_CANDIDATE )then
            status = 1; message = 'local ranges can only merge into a flex weights candidate'; return
        endif
        if( size(scheduled_rows) /= header%nptcls )then
            status = 1; message = 'flex weights scheduled-row mask has the wrong size'; return
        endif
        allocate(covered(header%nptcls), source=.false.)
        do i = 1, size(range_paths)
            call flex_weights_read_local_range(range_paths(i)%to_char(), generation, layout_digest, &
                &state_index, first_row, last_row, weights, flags, status, message)
            if( status /= 0 ) return
            if( generation /= header%generation .or. layout_digest /= header%layout_digest )then
                status = 1; message = 'local flex weights range belongs to another generation'; return
            endif
            if( state_index /= header%state_index )then
                status = 1; message = 'local flex weights range belongs to another state'; return
            endif
            if( first_row < 1 .or. last_row > header%nptcls )then
                status = 1; message = 'local flex weights range lies outside the project'; return
            endif
            if( any(covered(first_row:last_row)) )then
                status = 1; message = 'overlapping local flex weights ranges'; return
            endif
            if( .not. all(scheduled_rows(first_row:last_row)) )then
                status = 1; message = 'local flex weights range contains unscheduled rows'; return
            endif
            call flex_weights_write_rows(candidate_path, first_row, weights, flags, status, message)
            if( status /= 0 ) return
            covered(first_row:last_row) = .true.
            deallocate(weights, flags)
        enddo
        if( any(covered .neqv. scheduled_rows) )then
            status = 1; message = 'local flex weights ranges do not exactly cover the schedule'; return
        endif
    end subroutine flex_weights_merge_local_ranges

    ! ---- validation ----

    subroutine flex_weights_validate_identity(path, box, smpd, nptcls, layout_digest, status, message, &
        &expected_state, expected_index, expected_nstates)
        character(len=*), intent(in)  :: path
        integer,          intent(in)  :: box, nptcls
        real,             intent(in)  :: smpd
        integer(int64),   intent(in)  :: layout_digest
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        integer(int32), optional, intent(in) :: expected_state
        integer,        optional, intent(in) :: expected_index, expected_nstates
        type(flex_weights_header) :: header
        call flex_weights_validate_file(path, status, message)
        if( status /= 0 ) return
        call flex_weights_read_header(path, header, status, message)
        if( status /= 0 ) return
        if( header%box /= box .or. header%nptcls /= nptcls .or. header%layout_digest /= layout_digest .or. &
            &abs(header%smpd-real(smpd,real64)) > 1.0e-9_real64 )then
            status = 1; message = 'flex weights identity does not match the particle project'; return
        endif
        if( present(expected_state) )then
            if( header%state /= expected_state )then
                status = 1; message = 'flex weights transaction state does not match'; return
            endif
        endif
        if( present(expected_index) )then
            if( header%state_index /= expected_index )then
                status = 1; message = 'flex weights file belongs to another state'; return
            endif
        endif
        if( present(expected_nstates) )then
            if( header%nstates /= expected_nstates )then
                status = 1; message = 'flex weights file belongs to a delivery of another state count'; return
            endif
        endif
    end subroutine flex_weights_validate_identity

    !> Invariants of one state's file (checked blockwise, so any particle count fits):
    !!  - every weight finite and in [0,1]; every flag 0 or 1; a flagged row has a positive weight;
    !!  - an inactive row (active=.false.) is zero with flag 0;
    !!  - the scalars equal the reduction of the rows (mass, neff, population).
    subroutine flex_weights_validate_science(path, active, status, message)
        character(len=*), intent(in)  :: path
        logical,          intent(in)  :: active(:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(flex_weights_header) :: header
        real(real64),   allocatable :: scalars(:)
        real(real32),   allocatable :: weights(:)
        integer(int32), allocatable :: flags(:)
        real(real64) :: mass, sumsq, expected, tolerance
        integer(int64) :: pop
        integer :: first_row, last_row, i, row
        call flex_weights_validate_file(path, status, message, deep=.true.)
        if( status /= 0 ) return
        call flex_weights_read_header(path, header, status, message)
        if( status /= 0 ) return
        if( size(active) /= header%nptcls )then
            status = 1; message = 'flex weights validation mask has the wrong size'; return
        endif
        call flex_weights_read_scalars(path, scalars, status, message)
        if( status /= 0 ) return
        mass = 0.0_real64; sumsq = 0.0_real64; pop = 0_int64
        do first_row = 1, int(header%nptcls), ROW_BLOCK
            last_row = min(int(header%nptcls), first_row+ROW_BLOCK-1)
            call flex_weights_read_rows(path, first_row, last_row, weights, flags, status, message)
            if( status /= 0 ) return
            do i = 1, size(weights)
                row = first_row+i-1
                if( .not. ieee_is_finite(weights(i)) )then
                    status = 1; message = 'flex weights row has a non-finite value'; return
                endif
                if( weights(i) < 0.0_real32 .or. weights(i) > 1.0_real32 )then
                    status = 1; message = 'flex weights row has a value outside [0,1]'; return
                endif
                if( flags(i) < 0 .or. flags(i) > 1 )then
                    status = 1; message = 'flex weights row has an invalid label flag'; return
                endif
                if( .not. active(row) )then
                    if( weights(i) > 0.0_real32 .or. flags(i) /= 0 )then
                        status = 1; message = 'inactive particle carries flex weights'; return
                    endif
                    cycle
                endif
                if( flags(i) == 1 .and. weights(i) <= 0.0_real32 )then
                    status = 1; message = 'flex weights label flag points at a zero weight'; return
                endif
                mass  = mass  + real(weights(i),real64)
                sumsq = sumsq + real(weights(i),real64)**2
                if( flags(i) == 1 ) pop = pop + 1_int64
            enddo
            deallocate(weights, flags)
        enddo
        expected  = mass
        tolerance = MASS_TOL*max(1.0_real64, abs(expected))
        if( abs(scalars(FLEX_WEIGHTS_SCALAR_MASS)-expected) > tolerance )then
            status = 1; message = 'flex weights state mass does not match the rows'; return
        endif
        expected = 0.0_real64
        if( sumsq > 0.0_real64 ) expected = mass**2/sumsq
        tolerance = MASS_TOL*max(1.0_real64, abs(expected))
        if( abs(scalars(FLEX_WEIGHTS_SCALAR_NEFF)-expected) > tolerance )then
            status = 1; message = 'flex weights state effective size does not match the rows'; return
        endif
        if( nint(scalars(FLEX_WEIGHTS_SCALAR_POP)) /= pop )then
            status = 1; message = 'flex weights state population does not match the label flags'; return
        endif
        status = 0
        message = ''
    end subroutine flex_weights_validate_science

    !> Invariants across the files of one delivery: identical identity, generation, kind and state
    !! count, a distinct state index each, every particle flagged in at most one file, and for the
    !! PARTITION kind every flagged/weighted row summing to one across the files.
    subroutine flex_weights_validate_set(paths, active, status, message)
        type(string),     intent(in)  :: paths(:)
        logical,          intent(in)  :: active(:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(flex_weights_header) :: ref, header
        real(real32),   allocatable :: weights(:)
        integer(int32), allocatable :: flags(:)
        real(real64),   allocatable :: rowsum(:)
        integer,        allocatable :: nflag(:)
        logical,        allocatable :: seen(:)
        integer :: nfiles, s, first_row, last_row, i, row
        status = 0
        message = ''
        nfiles = size(paths)
        if( nfiles < 1 )then
            status = 1; message = 'empty flex weights set'; return
        endif
        call flex_weights_read_header(paths(1)%to_char(), ref, status, message)
        if( status /= 0 ) return
        if( ref%nstates /= nfiles )then
            status = 1; message = 'flex weights set does not hold every state of its delivery'; return
        endif
        if( size(active) /= ref%nptcls )then
            status = 1; message = 'flex weights set validation mask has the wrong size'; return
        endif
        allocate(rowsum(ref%nptcls), source=0.0_real64)
        allocate(nflag(ref%nptcls), source=0)
        allocate(seen(nfiles), source=.false.)
        do s = 1, nfiles
            call flex_weights_validate_science(paths(s)%to_char(), active, status, message)
            if( status /= 0 ) return
            call flex_weights_read_header(paths(s)%to_char(), header, status, message)
            if( status /= 0 ) return
            if( header%nptcls /= ref%nptcls .or. header%nstates /= ref%nstates .or. &
                &header%ncomp /= ref%ncomp .or. header%weight_kind /= ref%weight_kind .or. &
                &header%provenance /= ref%provenance .or. header%generation /= ref%generation .or. &
                &header%layout_digest /= ref%layout_digest .or. header%box /= ref%box .or. &
                &header%state /= ref%state .or. abs(header%smpd-ref%smpd) > 1.0e-9_real64 )then
                status = 1; message = 'flex weights set mixes files of different deliveries'; return
            endif
            if( seen(header%state_index) )then
                status = 1; message = 'flex weights set holds one state twice'; return
            endif
            seen(header%state_index) = .true.
            do first_row = 1, int(ref%nptcls), ROW_BLOCK
                last_row = min(int(ref%nptcls), first_row+ROW_BLOCK-1)
                call flex_weights_read_rows(paths(s)%to_char(), first_row, last_row, weights, flags, status, message)
                if( status /= 0 ) return
                do i = 1, size(weights)
                    row = first_row+i-1
                    rowsum(row) = rowsum(row) + real(weights(i),real64)
                    nflag(row)  = nflag(row)  + flags(i)
                enddo
                deallocate(weights, flags)
            enddo
        enddo
        if( .not. all(seen) )then
            status = 1; message = 'flex weights set is missing a state'; return
        endif
        do row = 1, int(ref%nptcls)
            if( nflag(row) > 1 )then
                status = 1; message = 'particle hard-labelled to more than one flex state'; return
            endif
            if( ref%weight_kind == FLEX_WEIGHTS_KIND_PARTITION .and. rowsum(row) > 0.0_real64 )then
                if( abs(rowsum(row)-1.0_real64) > ROW_SUM_TOL )then
                    status = 1; message = 'flex weights partition row does not sum to one across the states'; return
                endif
            endif
        enddo
        deallocate(rowsum, nflag, seen)
    end subroutine flex_weights_validate_set

    subroutine flex_weights_commit(candidate_path, committed_path, active, status, message)
        character(len=*), intent(in)  :: candidate_path, committed_path
        logical,          intent(in)  :: active(:)
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        call flex_weights_validate_science(candidate_path, active, status, message)
        if( status /= 0 ) return
        call flex_weights_publish(candidate_path, committed_path, status, message)
    end subroutine flex_weights_commit

    ! ---- reductions ----

    !> PARTITION when every nonzero row of the (nstates, nrows) table sums to one, else KERNEL
    pure integer(int32) function flex_weights_infer_kind(rows) result(kind)
        real(real32), intent(in) :: rows(:,:)
        real(real64) :: rowsum
        integer :: i
        kind = FLEX_WEIGHTS_KIND_PARTITION
        do i = 1, size(rows,2)
            rowsum = sum(real(rows(:,i),real64))
            if( rowsum <= 0.0_real64 ) cycle
            if( abs(rowsum-1.0_real64) > ROW_SUM_TOL )then
                kind = FLEX_WEIGHTS_KIND_KERNEL
                return
            endif
        enddo
    end function flex_weights_infer_kind

    !> One state's scalar section from its column: mass, effective size, hard population, bandwidth, target
    subroutine flex_weights_state_scalars(weights, flags, target, bandwidth, scalars)
        real(real32), intent(in) :: weights(:)
        integer(int32), intent(in) :: flags(:)
        real,         intent(in) :: target(:), bandwidth
        real(real64), allocatable, intent(out) :: scalars(:)
        real(real64) :: mass, sumsq
        allocate(scalars(FLEX_WEIGHTS_SCALAR_NFIELDS+size(target)), source=0.0_real64)
        mass  = sum(real(weights,real64))
        sumsq = sum(real(weights,real64)**2)
        scalars(FLEX_WEIGHTS_SCALAR_MASS) = mass
        if( sumsq > 0.0_real64 ) scalars(FLEX_WEIGHTS_SCALAR_NEFF) = mass**2/sumsq
        scalars(FLEX_WEIGHTS_SCALAR_POP)  = real(count(flags == 1), real64)
        scalars(FLEX_WEIGHTS_SCALAR_BW)   = real(bandwidth, real64)
        scalars(FLEX_WEIGHTS_SCALAR_NFIELDS+1:) = real(target, real64)
    end subroutine flex_weights_state_scalars

    ! ---- producer ----

    !> Write and publish every state file of a delivery from a selection-indexed table, in the
    !! working directory. `pinds` are the project rows of the selection; `weights_sel(nsel,nstates)`
    !! and `labels_sel(nsel)` are indexed like the flex state stage holds them. Rows outside the
    !! selection are zero, flag 0. Every candidate is validated on its own and as a set before any
    !! file is published; a stale file of a state beyond the delivered count is removed.
    subroutine flex_weights_write_store(box, smpd, box_crop, smpd_crop, layout_digest, nptcls, active, pinds, &
        &weights_sel, labels_sel, targets, bandwidths, provenance, status, message)
        integer,          intent(in)  :: box, box_crop, nptcls
        real,             intent(in)  :: smpd, smpd_crop
        integer(int64),   intent(in)  :: layout_digest
        logical,          intent(in)  :: active(:)         !< (nptcls)
        integer,          intent(in)  :: pinds(:)          !< (nsel) project rows
        real,             intent(in)  :: weights_sel(:,:)  !< (nsel, nstates)
        integer,          intent(in)  :: labels_sel(:)     !< (nsel)
        real,             intent(in)  :: targets(:,:)      !< (ncomp, nstates)
        real,             intent(in)  :: bandwidths(:)     !< (nstates)
        integer(int32),   intent(in)  :: provenance
        integer,          intent(out) :: status
        character(len=*), intent(out) :: message
        type(flex_weights_header) :: header
        type(string),   allocatable :: committed(:), candidates(:)
        real(real32),   allocatable :: rows(:,:), column(:)
        integer(int32), allocatable :: flags(:)
        real(real64),   allocatable :: scalars(:)
        logical,        allocatable :: seen(:)
        integer(int32) :: kind
        integer(int64) :: generation
        integer :: nsel, nstates, ncomp, i, s, first_row, last_row
        status  = 0
        message = ''
        nsel    = size(pinds)
        nstates = size(weights_sel,2)
        ncomp   = size(targets,1)
        if( nptcls < 1 .or. nstates < 1 .or. nsel < 1 .or. size(active) /= nptcls .or. &
            &size(weights_sel,1) /= nsel .or. size(labels_sel) /= nsel .or. &
            &size(targets,2) /= nstates .or. size(bandwidths) /= nstates )then
            status = 1; message = 'invalid flex weights delivery dimensions'; return
        endif
        allocate(rows(nstates,nptcls), source=0.0_real32)
        allocate(flags(nptcls), source=0_int32)
        allocate(seen(nptcls), source=.false.)
        do i = 1, nsel
            if( pinds(i) < 1 .or. pinds(i) > nptcls )then
                status = 1; message = 'flex weights delivery row outside the project'; return
            endif
            if( seen(pinds(i)) )then
                status = 1; message = 'duplicate project row in flex weights delivery'; return
            endif
            if( labels_sel(i) < 0 .or. labels_sel(i) > nstates )then
                status = 1; message = 'flex weights delivery label outside the state range'; return
            endif
            seen(pinds(i))   = .true.
            rows(:,pinds(i)) = real(weights_sel(i,:), real32)
        enddo
        deallocate(seen)
        kind = flex_weights_infer_kind(rows)
        allocate(committed(nstates), candidates(nstates))
        ! the files of one delivery share the generation: the next after state 1's committed file
        committed(1) = flex_weights_state_fname(1)
        call flex_weights_next_generation(committed(1)%to_char(), generation, status, message)
        if( status /= 0 )then
            generation = 1_int64
            status = 0
            message = ''
        endif
        allocate(column(nptcls))
        do s = 1, nstates
            committed(s)  = flex_weights_state_fname(s)
            candidates(s) = flex_weights_candidate_path(committed(s)%to_char(), generation)
            column = rows(s,:)
            flags  = 0_int32
            do i = 1, nsel
                if( labels_sel(i) == s ) flags(pinds(i)) = 1_int32
            enddo
            call flex_weights_init_header(header, nptcls, nstates, s, ncomp, box, smpd, box_crop, smpd_crop, &
                &kind, provenance, generation, layout_digest)
            call flex_weights_create_candidate(candidates(s)%to_char(), header, status, message)
            if( status /= 0 ) return
            call flex_weights_state_scalars(column, flags, targets(:,s), bandwidths(s), scalars)
            call flex_weights_write_scalars(candidates(s)%to_char(), scalars, status, message)
            if( status /= 0 ) return
            do first_row = 1, nptcls, ROW_BLOCK
                last_row = min(nptcls, first_row+ROW_BLOCK-1)
                call flex_weights_write_rows(candidates(s)%to_char(), first_row, column(first_row:last_row), &
                    &flags(first_row:last_row), status, message)
                if( status /= 0 ) return
            enddo
            deallocate(scalars)
        enddo
        deallocate(rows, column, flags)
        ! validate the whole delivery before publishing any of it
        call flex_weights_validate_set(candidates, active, status, message)
        if( status /= 0 ) return
        do s = 1, nstates
            call flex_weights_publish(candidates(s)%to_char(), committed(s)%to_char(), status, message)
            if( status /= 0 ) return
        enddo
        ! a previous delivery with more states leaves files this one must not be read with
        do s = nstates+1, nstates+FLEX_WEIGHTS_STALE_SCAN
            if( file_exists(flex_weights_state_fname(s)) ) call del_file(flex_weights_state_fname(s))
        enddo
        do s = 1, nstates
            call committed(s)%kill
            call candidates(s)%kill
        enddo
        deallocate(committed, candidates)
    end subroutine flex_weights_write_store

    !> The producer entry for a workflow: identity from the project and its particle field
    subroutine flex_weights_deliver(project, os, box, smpd, box_crop, smpd_crop, pinds, weights_sel, &
        &labels_sel, targets, bandwidths, provenance, status, message)
        class(sp_project), intent(inout) :: project
        class(oris),       intent(inout) :: os
        integer,           intent(in)    :: box, box_crop
        real,              intent(in)    :: smpd, smpd_crop
        integer,           intent(in)    :: pinds(:), labels_sel(:)
        real,              intent(in)    :: weights_sel(:,:), targets(:,:), bandwidths(:)
        integer(int32),    intent(in)    :: provenance
        integer,           intent(out)   :: status
        character(len=*),  intent(out)   :: message
        logical, allocatable :: active(:)
        integer(int64) :: layout_digest
        integer :: nptcls, i
        status  = 0
        message = ''
        nptcls  = os%get_noris()
        if( nptcls < 1 )then
            status = 1; message = 'flex weights delivery on an empty particle field'; return
        endif
        layout_digest = sigma2_state_project_layout_digest(project, os)
        if( layout_digest == 0_int64 )then
            status = 1; message = 'flex weights layout digest is undefined for this project'; return
        endif
        allocate(active(nptcls))
        do i = 1, nptcls
            active(i) = os%get_state(i) > 0
        enddo
        call flex_weights_write_store(box, smpd, box_crop, smpd_crop, layout_digest, nptcls, active, pinds, &
            &weights_sel, labels_sel, targets, bandwidths, provenance, status, message)
        deallocate(active)
    end subroutine flex_weights_deliver

    ! ---- consumers ----

    !> One state's file: registered in the out segment, file-valid, committed, same native grid,
    !! same ordered particle layout, and the state index it claims
    logical function flex_weights_consumable(project, os, box, smpd, state, message) result(l_ok)
        class(sp_project), intent(inout) :: project
        class(oris),       intent(inout) :: os
        integer,           intent(in)    :: box, state
        real,              intent(in)    :: smpd
        character(len=*),  intent(out)   :: message
        type(string)   :: path
        integer(int64) :: layout_digest
        integer        :: status
        logical        :: found
        l_ok    = .false.
        message = ''
        call project%get_flex_weights(state, path, found)
        if( .not. found )then
            message = 'flex weights of state '//int2str(state)//' are not registered in the project'
            return
        endif
        call flex_weights_validate_file(path%to_char(), status, message, deep=.true.)
        if( status /= 0 )then
            call path%kill
            return
        endif
        layout_digest = sigma2_state_project_layout_digest(project, os)
        if( layout_digest == 0_int64 )then
            message = 'flex weights layout digest is undefined for this project'
            call path%kill
            return
        endif
        call flex_weights_validate_identity(path%to_char(), box, smpd, os%get_noris(), layout_digest, &
            &status, message, expected_state=FLEX_WEIGHTS_COMMITTED, expected_index=state)
        l_ok = status == 0
        call path%kill
    end function flex_weights_consumable

    !> One state over the full layout: weights(nptcls), flags(nptcls), scalars and header
    subroutine flex_weights_load_state(project, os, box, smpd, state, header, weights, flags, scalars, &
        &status, message)
        class(sp_project), intent(inout) :: project
        class(oris),       intent(inout) :: os
        integer,           intent(in)    :: box, state
        real,              intent(in)    :: smpd
        type(flex_weights_header),   intent(out) :: header
        real(real32),   allocatable, intent(out) :: weights(:)
        integer(int32), allocatable, intent(out) :: flags(:)
        real(real64),   allocatable, intent(out) :: scalars(:)
        integer,           intent(out)   :: status
        character(len=*),  intent(out)   :: message
        type(string) :: path
        logical      :: found
        status = 1
        if( .not. flex_weights_consumable(project, os, box, smpd, state, message) ) return
        call project%get_flex_weights(state, path, found)
        if( .not. found )then
            message = 'flex weights path disappeared after validation'; return
        endif
        call flex_weights_read_header(path%to_char(), header, status, message)
        if( status == 0 ) call flex_weights_read_scalars(path%to_char(), scalars, status, message)
        if( status == 0 ) call flex_weights_read_rows(path%to_char(), 1, int(header%nptcls), weights, flags, &
            &status, message)
        call path%kill
    end subroutine flex_weights_load_state

    !> Every state of the registered delivery: weights(nstates,nptcls), labels(nptcls) from the
    !! flags, scalars(:,nstates); the files are checked as a set
    subroutine flex_weights_load_all(project, os, box, smpd, nstates, weights, labels, scalars, status, message)
        class(sp_project), intent(inout) :: project
        class(oris),       intent(inout) :: os
        integer,           intent(in)    :: box
        real,              intent(in)    :: smpd
        integer,           intent(out)   :: nstates
        real(real32),   allocatable, intent(out) :: weights(:,:)
        integer(int32), allocatable, intent(out) :: labels(:)
        real(real64),   allocatable, intent(out) :: scalars(:,:)
        integer,           intent(out)   :: status
        character(len=*),  intent(out)   :: message
        type(flex_weights_header) :: header
        type(string),   allocatable :: paths(:)
        real(real32),   allocatable :: column(:)
        integer(int32), allocatable :: flags(:)
        real(real64),   allocatable :: sc(:)
        logical,        allocatable :: active(:)
        integer :: s, i, nptcls
        logical :: found
        status  = 1
        nstates = 0
        if( .not. flex_weights_consumable(project, os, box, smpd, 1, message) ) return
        call flex_weights_load_state(project, os, box, smpd, 1, header, column, flags, sc, status, message)
        if( status /= 0 ) return
        nstates = int(header%nstates)
        nptcls  = int(header%nptcls)
        allocate(paths(nstates))
        do s = 1, nstates
            call project%get_flex_weights(s, paths(s), found)
            if( .not. found )then
                status = 1; message = 'flex weights of state '//int2str(s)//' are not registered in the project'
                return
            endif
        enddo
        allocate(active(nptcls))
        do i = 1, nptcls
            active(i) = os%get_state(i) > 0
        enddo
        call flex_weights_validate_set(paths, active, status, message)
        if( status /= 0 ) return
        allocate(weights(nstates,nptcls), labels(nptcls), scalars(size(sc),nstates))
        weights(1,:) = column
        scalars(:,1) = sc
        labels = 0_int32
        where( flags == 1 ) labels = 1_int32
        deallocate(column, flags, sc)
        do s = 2, nstates
            call flex_weights_load_state(project, os, box, smpd, s, header, column, flags, sc, status, message)
            if( status /= 0 ) return
            weights(s,:) = column
            scalars(:,s) = sc
            where( flags == 1 ) labels = int(s, int32)
            deallocate(column, flags, sc)
        enddo
        do s = 1, nstates
            call paths(s)%kill
        enddo
        deallocate(paths, active)
    end subroutine flex_weights_load_all

end module simple_flex_weights_state
