!@descr: project-owned even/odd partition and explicit selected-particle repartitioning, owned by the master
! Default partitions keep micrograph groups intact and assign ungrouped rows by index.
! Explicit selected-particle repartitioning may split those groups and records per-index provenance.
! Stack markers distinguish complete, missing, partial and inconsistent partitions.
module simple_eo_partition
use, intrinsic :: iso_fortran_env, only: int64
use simple_core_module_api
use simple_sp_project, only: sp_project
use simple_oris,       only: oris
implicit none

public :: ensure_field_halves, eo_transaction, eo_partition_status, repartition_selected_halves
public :: EO_COMPLETE, EO_MISSING, EO_PARTIAL, EO_INCONSISTENT, EO_SEED_DEFAULT
private
#include "simple_local_flags.inc"

integer,        parameter :: EO_COMPLETE = 0, EO_MISSING = 1, EO_PARTIAL = 2, EO_INCONSISTENT = 3
integer(int64), parameter :: EO_SEED_DEFAULT = 20261014_int64   !< the transaction's own random stream
integer(int64), parameter :: MINSTD_A = 48271_int64, MINSTD_M = 2147483647_int64
integer,        parameter :: MIN_GROUPS_WARN = 20               !< fewer independent groups than this warns
real,           parameter :: IMBALANCE_WARN  = 0.05             !< halves differing by more of the active warns
character(len=*), parameter :: KEY_GROUP = 'eo_group', KEY_INDEX = 'eo_byindex'

contains

    !> The halves of the field oritype, made by the master before any worker reads them. Particle fields: the
    !! transaction on the project file (its micrograph, stack and both particle segments), the labels and the
    !! provenance then copied to the loaded fields of spproj; without a project file, on spproj itself.
    !! Class-average fields: the per-index routine (no micrograph groups by nature), written to the project.
    !! A distributed worker only checks that its field has halves.
    subroutine ensure_field_halves( spproj, oritype, projfile, seed )
        class(sp_project),        intent(inout) :: spproj
        character(len=*),         intent(in)    :: oritype
        class(string),            intent(in)    :: projfile
        integer(int64), optional, intent(in)    :: seed
        type(sp_project), allocatable :: proj
        class(oris), pointer :: field
        character(len=STDLEN) :: msg
        integer :: status
        logical :: l_changed, l_recorded
        if( l_distr_worker_glob )then
            ! a worker consumes the halves its master made
            call spproj%ptr2oritype(oritype, field)
            if( field%get_nevenodd() == 0 ) THROW_HARD('no even/odd partition: the master makes it; ensure_field_halves')
            nullify(field)
            return
        endif
        select case(trim(oritype))
            case('ptcl2D', 'ptcl3D')
                ! every loaded stack row recorded and the field labelled: the loaded fields are validated in
                ! memory, the project file is read only when something is to be labelled
                if( spproj%os_stk%get_noris() > 0 )then
                    call spproj%ptr2oritype(oritype, field)
                    l_recorded = field%get_nodd() > 0 .and. all_marked(spproj)
                    nullify(field)
                    if( l_recorded )then
                        status = eo_partition_status(spproj, msg)
                        if( status == EO_COMPLETE ) return
                        if( status == EO_INCONSISTENT )then
                            write(logfhandle,'(A)') '>>> EO PARTITION: '//trim(msg)
                            THROW_HARD('inconsistent even/odd partition in the project; ensure_field_halves')
                        endif
                    endif
                endif
                if( file_exists(projfile) .and. fname2format(projfile) == 'O' )then
                    allocate(proj)
                    call proj%read_segment('mic',    projfile)
                    call proj%read_segment('stk',    projfile)
                    call proj%read_segment('ptcl2D', projfile)
                    call proj%read_segment('ptcl3D', projfile)
                    call eo_transaction(proj, l_changed, seed)
                    if( l_changed )then
                        if( proj%os_ptcl2D%get_noris() > 0 ) call proj%write_segment_inside('ptcl2D', projfile)
                        if( proj%os_ptcl3D%get_noris() > 0 ) call proj%write_segment_inside('ptcl3D', projfile)
                        if( proj%os_stk%get_noris()    > 0 ) call proj%write_segment_inside('stk',    projfile)
                    endif
                    call copy_halves(proj, spproj)
                    call proj%kill
                    deallocate(proj)
                else
                    call eo_transaction(spproj, l_changed, seed)
                endif
            case default
                call spproj%ptr2oritype(oritype, field)
                if( field%get_nevenodd() == 0 )then
                    call field%partition_eo
                    call spproj%write_segment_inside(oritype, projfile)
                endif
                nullify(field)
        end select
    end subroutine ensure_field_halves

    !> The transaction on the fields in memory: validates the partition, labels the stack rows without provenance
    !! (never changing a labelled one), the same labels in ptcl2D and ptcl3D, the provenance on the stack rows;
    !! l_changed when anything was labelled. Stops on an inconsistent partition.
    subroutine eo_transaction( spproj, l_changed, seed )
        class(sp_project),        intent(inout) :: spproj
        logical,                  intent(out)   :: l_changed
        integer(int64), optional, intent(in)    :: seed
        character(len=STDLEN) :: msg
        integer(int64) :: rng
        integer :: status
        rng = EO_SEED_DEFAULT
        if( present(seed) ) rng = seed
        call run_transaction(spproj, rng, status, l_changed, msg)
        if( status == EO_INCONSISTENT )then
            write(logfhandle,'(A)') '>>> EO PARTITION: '//trim(msg)
            THROW_HARD('inconsistent even/odd partition in the project; eo_transaction')
        endif
    end subroutine eo_transaction

    subroutine repartition_selected_halves( spproj, groups )
        use simple_type_defs, only: class_sample
        class(sp_project), intent(inout) :: spproj
        type(class_sample), intent(in) :: groups(:)
        logical, allocatable :: touched_stacks(:)
        real(dp) :: direction_scores(0:1), total_scores(0:1)
        integer :: counts(0:1), group, member, half, stack, status
        character(len=STDLEN) :: msg
        if( l_distr_worker_glob ) THROW_HARD('only the master may repartition selected particles')
        status = eo_partition_status(spproj,msg)
        if( status /= EO_COMPLETE ) THROW_HARD('repartitioning requires complete existing half-set provenance')
        allocate(touched_stacks(spproj%get_nstks()), source=.false.)
        counts = 0
        total_scores = 0._dp
        do group = 1,size(groups)
            direction_scores = 0._dp
            do member = 1,groups(group)%nsample,2
                half = 0
                if( member == groups(group)%nsample )then
                    if( counts(0) > counts(1) ) half = 1
                    if( counts(0) == counts(1) .and. total_scores(0) > total_scores(1) ) half = 1
                    call assign_half(groups(group)%pinds(member),half,groups(group)%ccs(member))
                else
                    if( direction_scores(0) > direction_scores(1) ) half = 1
                    if( direction_scores(0) == direction_scores(1) .and. total_scores(0) > total_scores(1) ) half = 1
                    call assign_half(groups(group)%pinds(member),half,groups(group)%ccs(member))
                    call assign_half(groups(group)%pinds(member+1),1-half,groups(group)%ccs(member+1))
                endif
            enddo
        enddo
        do stack = 1,size(touched_stacks)
            if( .not.touched_stacks(stack) ) cycle
            call spproj%os_stk%delete_entry(stack,KEY_GROUP)
            call spproj%os_stk%set(stack,KEY_INDEX,1)
        enddo
        status = eo_partition_status(spproj,msg)
        if( status /= EO_COMPLETE ) THROW_HARD('score-balanced repartition produced inconsistent half-set provenance')
    contains
        subroutine assign_half( particle, label, score )
            integer, intent(in) :: particle, label
            real, intent(in) :: score
            call spproj%os_ptcl2D%set(particle,'eo',label)
            call spproj%os_ptcl3D%set(particle,'eo',label)
            touched_stacks(spproj%os_ptcl3D%get_int(particle,'stkind')) = .true.
            counts(label) = counts(label) + 1
            total_scores(label) = total_scores(label) + real(score,dp)
            direction_scores(label) = direction_scores(label) + real(score,dp)
        end subroutine assign_half
    end subroutine repartition_selected_halves

    !> The state of the project's partition without changing it: EO_COMPLETE, EO_MISSING, EO_PARTIAL or
    !! EO_INCONSISTENT (msg says why)
    integer function eo_partition_status( spproj, msg ) result( status )
        class(sp_project), intent(inout) :: spproj
        character(len=*),  intent(out)   :: msg
        integer, allocatable :: stkind(:)
        logical, allocatable :: active(:)
        integer :: nrows, nstks
        call field_layout(spproj, nrows, nstks, stkind, active, status, msg)
        if( status /= EO_COMPLETE ) return
        call check_marked(spproj, nrows, nstks, stkind, active, status, msg)
    end function eo_partition_status

    ! ---- private ----

    subroutine run_transaction( spproj, rng, status, l_changed, msg )
        class(sp_project), intent(inout) :: spproj
        integer(int64),    intent(inout) :: rng
        integer,           intent(out)   :: status
        logical,           intent(out)   :: l_changed
        character(len=*),  intent(out)   :: msg
        integer, allocatable :: stkind(:), stkpop(:), stkhalf(:), order(:)
        logical, allocatable :: active(:), unmarked(:)
        integer :: nrows, nstks, i, s, k, nunmarked, half, pops(2), ngroups, nmics, nmic_with_ptcls, imic
        logical :: l_groups, l_lab2, l_lab3
        l_changed = .false.
        call field_layout(spproj, nrows, nstks, stkind, active, status, msg)
        if( status == EO_INCONSISTENT ) return
        if( nrows == 0 ) return
        if( nstks == 0 )then
            ! particles without stack rows: per index on both fields, no provenance to record
            if( spproj%os_ptcl2D%get_nodd() == 0 .and. spproj%os_ptcl3D%get_nodd() == 0 )then
                call label_by_index_all(spproj, nrows)
                l_changed = .true.
            endif
            status = EO_COMPLETE
            return
        endif
        ! a project of an earlier release: labels but no provenance; its labels are kept and marked per index
        if( .not. any_marked(spproj, nstks) )then
            l_lab2 = spproj%os_ptcl2D%get_nodd() > 0
            l_lab3 = spproj%os_ptcl3D%get_nodd() > 0
            if( l_lab2 .or. l_lab3 )then
                call adopt_legacy(spproj, nrows, nstks, active, l_lab2, l_lab3, status, msg)
                if( status == EO_INCONSISTENT ) return
                l_changed = .true.
                status    = EO_COMPLETE
                return
            endif
        endif
        call check_marked(spproj, nrows, nstks, stkind, active, status, msg)
        if( status == EO_INCONSISTENT .or. status == EO_COMPLETE ) return
        ! the unmarked stack rows with particles get their labels
        allocate(stkpop(nstks), source=0)
        allocate(stkhalf(nstks), source=-1)
        allocate(unmarked(nstks), source=.false.)
        do s = 1, nstks
            unmarked(s) = .not. (spproj%os_stk%isthere(s, KEY_GROUP) .or. spproj%os_stk%isthere(s, KEY_INDEX))
        end do
        pops = 0
        do i = 1, nrows
            s = stkind(i)
            if( unmarked(s) )then
                if( active(i) ) stkpop(s) = stkpop(s) + 1
            else if( active(i) )then
                half = spproj%os_ptcl3D%get_eo(i)
                if( spproj%os_ptcl3D%get_noris() == 0 ) half = spproj%os_ptcl2D%get_eo(i)
                if( half == 0 .or. half == 1 ) pops(half+1) = pops(half+1) + 1
            endif
        end do
        nunmarked = count(unmarked .and. stack_has_rows(nrows, nstks, stkind))
        ! route: micrograph groups when the stack rows are the micrographs with particles, one to one
        l_groups = .false.
        nmics    = spproj%os_mic%get_noris()
        if( nmics > 0 )then
            nmic_with_ptcls = 0
            do imic = 1, nmics
                if( spproj%os_mic%get_int(imic, 'nptcls') > 0 ) nmic_with_ptcls = nmic_with_ptcls + 1
            end do
            l_groups = nmics == nstks .or. nmic_with_ptcls == nstks
        endif
        if( l_groups )then
            ngroups = count(stkpop > 0)
            if( ngroups + count(.not. unmarked) < 2 )then
                write(logfhandle,'(A,I0)') '>>> EO PARTITION: micrograph groups with active particles: ', ngroups
                THROW_HARD('fewer than two micrograph groups: no two half-sets; run_transaction')
            endif
            ! largest group first, a random order among equal ones, each to the half with fewer active particles
            order = pack([(s, s=1,nstks)], unmarked .and. stack_has_rows(nrows, nstks, stkind))
            call eo_shuffle(order, rng)
            call stable_sort_desc(order, stkpop)
            do k = 1, size(order)
                s = order(k)
                if( pops(1) < pops(2) )then
                    half = 0
                else if( pops(2) < pops(1) )then
                    half = 1
                else
                    half = int(mod(next_random(rng), 2_int64))
                endif
                pops(half+1) = pops(half+1) + stkpop(s)
                stkhalf(s)   = half
                call spproj%os_stk%set(s, KEY_GROUP, real(half))
            end do
            do i = 1, nrows
                if( stkhalf(stkind(i)) >= 0 ) call set_label(spproj, i, stkhalf(stkind(i)))
            end do
            write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> EO PARTITION: MICROGRAPH GROUPS, ', size(order), &
                &' labelled; active particles even/odd ', pops(1), ' / ', pops(2), '; independent groups ', &
                &ngroups + count(.not. unmarked)
            if( ngroups + count(.not. unmarked) < MIN_GROUPS_WARN )then
                THROW_WARN('few independent micrograph groups for the half-sets; run_transaction')
            endif
            if( real(abs(pops(1) - pops(2))) > IMBALANCE_WARN * real(max(1, sum(pops))) )then
                THROW_WARN('half-sets differ by more than 5 % of the active particles; run_transaction')
            endif
        else
            ! no scientific grouping known: row parity, as the per-index routine labels
            do i = 1, nrows
                if( unmarked(stkind(i)) ) call set_label(spproj, i, index_half(i))
            end do
            do s = 1, nstks
                if( unmarked(s) ) call spproj%os_stk%set(s, KEY_INDEX, 1.)
            end do
            write(logfhandle,'(A,I0,A)') '>>> EO PARTITION: PER-INDEX ROUTE (a stack row without a micrograph row), ', &
                &nunmarked, ' stack row(s) labelled by particle index'
        endif
        l_changed = .true.
        status    = EO_COMPLETE
    end subroutine run_transaction

    !> rows, stack rows, each row's stack row and activity (ptcl2D's state, ptcl3D's when there is no ptcl2D)
    subroutine field_layout( spproj, nrows, nstks, stkind, active, status, msg )
        class(sp_project),    intent(inout) :: spproj
        integer,              intent(out)   :: nrows, nstks
        integer, allocatable, intent(out)   :: stkind(:)
        logical, allocatable, intent(out)   :: active(:)
        integer,              intent(out)   :: status
        character(len=*),     intent(out)   :: msg
        integer :: n2, n3, i
        status = EO_COMPLETE
        msg    = ''
        n2     = spproj%os_ptcl2D%get_noris()
        n3     = spproj%os_ptcl3D%get_noris()
        nstks  = spproj%os_stk%get_noris()
        nrows  = max(n2, n3)
        allocate(stkind(nrows), source=0)
        allocate(active(nrows), source=.false.)
        if( n2 > 0 .and. n3 > 0 .and. n2 /= n3 )then
            status = EO_INCONSISTENT
            write(msg,'(A,I0,A,I0,A)') 'ptcl2D has ', n2, ' rows and ptcl3D ', n3, ': the half-sets need congruent fields'
            return
        endif
        do i = 1, nrows
            if( n2 > 0 )then
                active(i) = spproj%os_ptcl2D%get_state(i) > 0
                if( nstks > 0 ) stkind(i) = spproj%os_ptcl2D%get_int(i, 'stkind')
            else
                active(i) = spproj%os_ptcl3D%get_state(i) > 0
                if( nstks > 0 ) stkind(i) = spproj%os_ptcl3D%get_int(i, 'stkind')
            endif
            if( nstks > 0 .and. (stkind(i) < 1 .or. stkind(i) > nstks) )then
                status = EO_INCONSISTENT
                write(msg,'(A,I0,A)') 'particle ', i, ' has no stack row'
                return
            endif
        end do
    end subroutine field_layout

    !> EO_COMPLETE when every stack row with particles is marked and its particles agree, EO_PARTIAL when some
    !! are unmarked, EO_MISSING when none is and no particle is labelled, EO_INCONSISTENT otherwise
    subroutine check_marked( spproj, nrows, nstks, stkind, active, status, msg )
        class(sp_project), intent(inout) :: spproj
        integer,           intent(in)    :: nrows, nstks, stkind(nrows)
        logical,           intent(in)    :: active(nrows)
        integer,           intent(out)   :: status
        character(len=*),  intent(out)   :: msg
        logical, allocatable :: has_rows(:)
        integer :: i, s, label, e2, e3, nmarked, nwith
        status = EO_COMPLETE
        msg    = ''
        if( nstks == 0 ) return
        has_rows = stack_has_rows(nrows, nstks, stkind)
        nmarked  = 0
        nwith    = count(has_rows)
        do s = 1, nstks
            if( .not. has_rows(s) ) cycle
            if( spproj%os_stk%isthere(s, KEY_GROUP) .or. spproj%os_stk%isthere(s, KEY_INDEX) ) nmarked = nmarked + 1
        end do
        do i = 1, nrows
            if( .not. active(i) ) cycle
            s  = stkind(i)
            e2 = -1
            e3 = -1
            if( spproj%os_ptcl2D%get_noris() > 0 ) e2 = spproj%os_ptcl2D%get_eo(i)
            if( spproj%os_ptcl3D%get_noris() > 0 ) e3 = spproj%os_ptcl3D%get_eo(i)
            if( spproj%os_stk%isthere(s, KEY_GROUP) )then
                label = nint(spproj%os_stk%get(s, KEY_GROUP))
                if( (e2 >= 0 .and. e2 /= label) .or. (e3 >= 0 .and. e3 /= label) )then
                    status = EO_INCONSISTENT
                    write(msg,'(A,I0,A,I0)') 'particle ', i, ' disagrees with the half of its micrograph group, stack row ', s
                    return
                endif
            else if( spproj%os_stk%isthere(s, KEY_INDEX) )then
                if( e2 >= 0 .and. e3 >= 0 .and. e2 /= e3 )then
                    status = EO_INCONSISTENT
                    write(msg,'(A,I0,A)') 'particle ', i, ' has different halves in ptcl2D and ptcl3D'
                    return
                endif
            endif
        end do
        if( nmarked == nwith )then
            status = EO_COMPLETE
        else if( nmarked == 0 .and. spproj%os_ptcl2D%get_nodd() == 0 .and. spproj%os_ptcl3D%get_nodd() == 0 )then
            status = EO_MISSING
        else
            status = EO_PARTIAL
        endif
    end subroutine check_marked

    !> labels made before the provenance existed: kept, copied to an unlabelled twin field, marked per index
    subroutine adopt_legacy( spproj, nrows, nstks, active, l_lab2, l_lab3, status, msg )
        class(sp_project), intent(inout) :: spproj
        integer,           intent(in)    :: nrows, nstks
        logical,           intent(in)    :: active(nrows), l_lab2, l_lab3
        integer,           intent(out)   :: status
        character(len=*),  intent(out)   :: msg
        integer :: i, s, n2, n3
        status = EO_COMPLETE
        msg    = ''
        n2 = spproj%os_ptcl2D%get_noris()
        n3 = spproj%os_ptcl3D%get_noris()
        if( n2 > 0 .and. n3 > 0 )then
            if( l_lab2 .and. l_lab3 )then
                do i = 1, nrows
                    if( .not. active(i) ) cycle
                    if( spproj%os_ptcl2D%get_eo(i) /= spproj%os_ptcl3D%get_eo(i) )then
                        status = EO_INCONSISTENT
                        write(msg,'(A,I0,A)') 'particle ', i, ' has different halves in ptcl2D and ptcl3D'
                        return
                    endif
                end do
            else if( l_lab2 )then
                do i = 1, nrows
                    call spproj%os_ptcl3D%set(i, 'eo', real(spproj%os_ptcl2D%get_eo(i)))
                end do
            else
                do i = 1, nrows
                    call spproj%os_ptcl2D%set(i, 'eo', real(spproj%os_ptcl3D%get_eo(i)))
                end do
            endif
        endif
        do s = 1, nstks
            call spproj%os_stk%set(s, KEY_INDEX, 1.)
        end do
        if( spproj%os_mic%get_noris() == 0 )then
            write(logfhandle,'(A)') '>>> EO PARTITION: PER-INDEX ROUTE (a stack row without a micrograph row); '//&
                &'the labels the project holds are kept, recorded as per-index'
        else
            write(logfhandle,'(A)') '>>> EO PARTITION: labels without provenance (an earlier release) kept, recorded as per-index'
        endif
    end subroutine adopt_legacy

    subroutine label_by_index_all( spproj, nrows )
        class(sp_project), intent(inout) :: spproj
        integer,           intent(in)    :: nrows
        integer :: i
        do i = 1, nrows
            call set_label(spproj, i, index_half(i))
        end do
    end subroutine label_by_index_all

    !> the half of row i on both particle fields
    subroutine set_label( spproj, i, half )
        class(sp_project), intent(inout) :: spproj
        integer,           intent(in)    :: i, half
        if( spproj%os_ptcl2D%get_noris() > 0 ) call spproj%os_ptcl2D%set(i, 'eo', real(half))
        if( spproj%os_ptcl3D%get_noris() > 0 ) call spproj%os_ptcl3D%set(i, 'eo', real(half))
    end subroutine set_label

    !> row parity, the rule of the per-index routine (oris%partition_eo): odd rows even, even rows odd
    pure integer function index_half( i )
        integer, intent(in) :: i
        index_half = merge(0, 1, mod(i, 2) == 1)
    end function index_half

    logical function any_marked( spproj, nstks )
        class(sp_project), intent(inout) :: spproj
        integer,           intent(in)    :: nstks
        integer :: s
        any_marked = .false.
        do s = 1, nstks
            if( spproj%os_stk%isthere(s, KEY_GROUP) .or. spproj%os_stk%isthere(s, KEY_INDEX) )then
                any_marked = .true.
                return
            endif
        end do
    end function any_marked

    logical function all_marked( spproj )
        class(sp_project), intent(inout) :: spproj
        integer :: s
        all_marked = .false.
        do s = 1, spproj%os_stk%get_noris()
            if( .not. (spproj%os_stk%isthere(s, KEY_GROUP) .or. spproj%os_stk%isthere(s, KEY_INDEX)) ) return
        end do
        all_marked = .true.
    end function all_marked

    function stack_has_rows( nrows, nstks, stkind ) result( has )
        integer, intent(in) :: nrows, nstks, stkind(nrows)
        logical :: has(nstks)
        integer :: i
        has = .false.
        do i = 1, nrows
            has(stkind(i)) = .true.
        end do
    end function stack_has_rows

    !> the labels and the provenance of the transaction's fields into the fields dst has loaded
    subroutine copy_halves( src, dst )
        class(sp_project), intent(inout) :: src, dst
        integer :: i
        if( rows_match(src%os_stk, dst%os_stk) )then
            do i = 1, src%os_stk%get_noris()
                if( src%os_stk%isthere(i, KEY_GROUP) ) call dst%os_stk%set(i, KEY_GROUP, src%os_stk%get(i, KEY_GROUP))
                if( src%os_stk%isthere(i, KEY_INDEX) ) call dst%os_stk%set(i, KEY_INDEX, src%os_stk%get(i, KEY_INDEX))
            end do
        endif
        if( rows_match(src%os_ptcl2D, dst%os_ptcl2D) )then
            do i = 1, src%os_ptcl2D%get_noris()
                call dst%os_ptcl2D%set(i, 'eo', real(src%os_ptcl2D%get_eo(i)))
            end do
        endif
        if( rows_match(src%os_ptcl3D, dst%os_ptcl3D) )then
            do i = 1, src%os_ptcl3D%get_noris()
                call dst%os_ptcl3D%set(i, 'eo', real(src%os_ptcl3D%get_eo(i)))
            end do
        endif
    end subroutine copy_halves

    !> a loaded field matches the project file's row for row; an unloaded one is left alone
    logical function rows_match( src, dst )
        class(oris), intent(inout) :: src, dst
        rows_match = .false.
        if( dst%get_noris() == 0 ) return
        if( dst%get_noris() /= src%get_noris() )then
            THROW_HARD('a loaded field differs from the project file in its rows; copy_halves')
        endif
        rows_match = .true.
    end function rows_match

    !> MINSTD (Park and Miller, multiplier 48271): the transaction's own generator, never the
    !! process-wide random state. Every product fits in 64 bits: integer overflow is undefined in
    !! Fortran, and an optimising compiler may assume it never happens.
    integer(int64) function next_random( state )
        integer(int64), intent(inout) :: state
        if( state <= 0_int64 .or. state >= MINSTD_M ) state = 1_int64 + modulo(state, MINSTD_M - 1_int64)
        state       = mod(state * MINSTD_A, MINSTD_M)
        next_random = state
    end function next_random

    !> Fisher-Yates with the transaction's generator
    subroutine eo_shuffle( v, state )
        integer,        intent(inout) :: v(:)
        integer(int64), intent(inout) :: state
        integer :: i, j, t
        do i = size(v), 2, -1
            j    = 1 + int(mod(next_random(state), int(i,int64)))
            t    = v(i)
            v(i) = v(j)
            v(j) = t
        end do
    end subroutine eo_shuffle

    !> stable merge sort of the stack rows v by key(v) descending
    subroutine stable_sort_desc( v, key )
        integer, intent(inout) :: v(:)
        integer, intent(in)    :: key(:)
        integer, allocatable :: tmp(:)
        integer :: n, width, lo, mid, hi, i, j, k
        n = size(v)
        allocate(tmp(n))
        width = 1
        do while( width < n )
            lo = 1
            do while( lo <= n )
                mid = min(lo + width - 1, n)
                hi  = min(lo + 2*width - 1, n)
                i = lo
                j = mid + 1
                do k = lo, hi
                    if( i <= mid .and. j <= hi )then
                        ! the left run wins ties: stable
                        if( key(v(i)) >= key(v(j)) )then
                            tmp(k) = v(i)
                            i = i + 1
                        else
                            tmp(k) = v(j)
                            j = j + 1
                        endif
                    else if( i <= mid )then
                        tmp(k) = v(i)
                        i = i + 1
                    else
                        tmp(k) = v(j)
                        j = j + 1
                    endif
                end do
                lo = hi + 1
            end do
            v     = tmp
            width = 2 * width
        end do
    end subroutine stable_sort_desc

end module simple_eo_partition
