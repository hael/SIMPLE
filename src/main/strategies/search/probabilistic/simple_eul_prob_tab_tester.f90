!@descr: unit tests for the 3D probability table: the seed shift of an assigned candidate, the dense part files and ref_assign
! The seed table holds the particle-frame seed every candidate of a particle was scored with; an assigned, unrefined
! candidate records it as it is, whatever its in-plane angle. A refined candidate keeps the shift of its search.
! ref_assign is compared with the assignment it replaced, kept here as the reference (reference_assign): on synthetic
! dense tables of one and of several states, with the same random seed, both give every particle the same reference.
! A particle's state is the state of its closest reference; an exact tie goes to the lowest reference index.
module simple_eul_prob_tab_tester
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_builder,            only: builder
use simple_parameters,         only: parameters
use simple_eul_prob_tab,       only: eul_prob_tab
use simple_eul_prob_tab_utils, only: record_seed_shift, calc_athres, sample_likelihood_dist, sample_likelihood_dists, &
    &prob_candidate, prob_dense_entry, prob_refined_entry
use simple_rnd,                only: seed_rnd_fixed
use simple_type_defs,          only: ptcl_ref
use simple_test_utils
!$ use omp_lib
implicit none
private
public :: run_all_eul_prob_tab_tests

! eight in-plane rotations
integer, parameter :: NROTS = 8
real,    parameter :: SEED(2) = [1.5, -0.5]
! the synthetic dense tables of the assignment comparison
integer, parameter :: NPTCLS_TAB = 150
integer, parameter :: SEED_TAB   = 20261016   !< the distances and the particles' angular distances
integer, parameter :: SEED_DRAW  = 20261017   !< the assignment's draws, the same for both assignments

contains

    subroutine run_all_eul_prob_tab_tests()
        write(*,'(A)') '**** running all 3D probability table tests ****'
        call test_record_seed_shift()
        write(*,'(A)') 'test_dense_part_round_trip'
        call test_dense_part_round_trip()
        write(*,'(A)') 'test_dense_read_threads'
        call test_dense_read_threads()
        write(*,'(A)') 'test_likelihood_samplers'
        call test_likelihood_samplers()
        write(*,'(A)') 'test_ref_assign_against_reference'
        call test_ref_assign_against_reference()
    end subroutine run_all_eul_prob_tab_tests

    !> One state of 200 directions, three states of 200, and two states of one direction (each state's draw
    !! has a single reference). A different draw seed changes the reference assignment of the first two
    !! tables, so their agreement is not that of a draw-free assignment. Every table is random but holds the tie cases: particle 1 is flat (every
    !! distance equal), particle 2 ties its minimum between a reference of state 2 and one of state 3,
    !! particle 3 between state 1 and state 3, particle 4 has only invalid distances, and particles 5 and
    !! 6 have identical distances (ties between particles in every reference's order).
    subroutine test_ref_assign_against_reference()
        call compare_with_reference(nspace=200, nstates=1)
        call compare_with_reference(nspace=200, nstates=3)
        call compare_with_reference(nspace=1,  nstates=2)
    end subroutine test_ref_assign_against_reference

    subroutine compare_with_reference( nspace, nstates )
        integer, intent(in) :: nspace, nstates
        class(parameters), allocatable, target :: params
        class(builder),    allocatable, target :: build
        type(eul_prob_tab) :: tab
        real,    allocatable :: dists(:,:), state_athres(:)
        integer, allocatable :: ref_expected(:), ref_other(:), pinds(:)
        integer :: nrefs, i, iref, istate, nsame
        character(len=64) :: label
        write(label,'(A,I0,A,I0)') 'nspace ', nspace, ', nstates ', nstates
        nrefs = nspace * nstates
        allocate(params, build)
        params%nspace     = nspace
        params%nstates    = nstates
        params%prob_athres = 10.
        params%l_doshift  = .false.
        ! the particle field: states in turn, every particle sampled, angular distances of 0 to 20 degrees
        call seed_rnd_fixed(SEED_TAB)
        allocate(build%spproj_field)
        call build%spproj_field%new(NPTCLS_TAB, is_ptcl=.true.)
        do i = 1,NPTCLS_TAB
            call build%spproj_field%set_state(i, mod(i-1, nstates) + 1)
            call build%spproj_field%set(i, 'sampled', 1.)
            call build%spproj_field%set(i, 'dist',    20. * ran3())
        enddo
        pinds = [(i, i=1,NPTCLS_TAB)]
        ! the distances, with the tie cases
        allocate(dists(nrefs,NPTCLS_TAB))
        do i = 1,NPTCLS_TAB
            do iref = 1,nrefs
                dists(iref,i) = 5. * ran3()
            enddo
        enddo
        dists(:,1) = 1.
        if( nstates >= 3 )then
            dists(:,2) = 3. + dists(:,2)
            dists(nspace + min(5,nspace), 2)   = 0.01
            dists(2*nspace + min(7,nspace), 2) = 0.01
            dists(:,3) = 3. + dists(:,3)
            dists(min(10,nspace), 3)    = 0.01
            dists(2*nspace + 1, 3)      = 0.01
        endif
        dists(:,4) = huge(1.)
        dists(:,6) = dists(:,5)
        ! the new assignment
        call tab%new(params, build, pinds)
        do i = 1,NPTCLS_TAB
            do iref = 1,nrefs
                call set_dist(tab, iref, i, dists(iref,i))
            enddo
        enddo
        call seed_rnd_fixed(SEED_DRAW)
        call tab%ref_assign
        ! the reference assignment on the same table and draws
        allocate(state_athres(nstates))
        do istate = 1,nstates
            state_athres(istate) = calc_athres(build%spproj_field, 'dist', params%prob_athres, state=istate)
        enddo
        call seed_rnd_fixed(SEED_DRAW)
        call reference_assign(dists, nspace, nstates, params%prob_athres, state_athres, ref_expected)
        nsame = 0
        do i = 1,NPTCLS_TAB
            iref = ref_expected(i)
            if( tab%assgn_map(i)%istate == (iref - 1)/nspace + 1 .and. &
                &tab%assgn_map(i)%iproj == modulo(iref - 1, nspace) + 1 ) nsame = nsame + 1
        enddo
        call assert_int(NPTCLS_TAB, nsame, 'ref_assign agrees with the reference assignment ('//trim(label)//')')
        if( nspace > 1 )then
            call seed_rnd_fixed(SEED_DRAW + 1)
            call reference_assign(dists, nspace, nstates, params%prob_athres, state_athres, ref_other)
            call assert_true(any(ref_other /= ref_expected), 'the draws do not change the assignment ('//trim(label)//')')
        endif
        ! the written tie rule: the state of the lowest-index reference at the particle's minimum
        call assert_int(1, tab%assgn_map(1)%istate, 'a flat particle takes the state of reference 1 ('//trim(label)//')')
        call assert_int(1, tab%assgn_map(4)%istate, 'an invalid particle takes the state of reference 1 ('//trim(label)//')')
        if( nstates >= 3 )then
            call assert_int(2, tab%assgn_map(2)%istate, 'a state 2 / state 3 tie goes to state 2 ('//trim(label)//')')
            call assert_int(1, tab%assgn_map(3)%istate, 'a state 1 / state 3 tie goes to state 1 ('//trim(label)//')')
        endif
        call tab%kill
        call build%spproj_field%kill
        deallocate(build%spproj_field)
        deallocate(params, build)
    end subroutine compare_with_reference

    !> Two dense part files written by worker tables in batches, read by a master whose particle order is the reverse:
    !! every row, refined shift and seed comes back exactly, an unrefined entry without shift
    subroutine test_dense_part_round_trip()
        integer, parameter :: NSPACE_RT = 50, NSTATES_RT = 2, NPTCLS_RT = 40, NPARTS_RT = 2
        class(parameters), allocatable, target :: params
        class(builder),    allocatable, target :: build
        type(eul_prob_tab) :: master
        type(prob_dense_entry), allocatable :: rows(:,:)
        real,                   allocatable :: shifts(:,:,:), seeds(:,:)
        logical,                allocatable :: refined(:,:)
        integer,                allocatable :: pinds(:)
        type(prob_candidate) :: cand
        integer :: nrefs, i, ri, i_glob, nwrong
        nrefs = NSPACE_RT * NSTATES_RT
        call make_round_trip_parts(params, build, NSPACE_RT, NSTATES_RT, NPTCLS_RT, NPARTS_RT, rows, shifts, refined, seeds)
        pinds = [(i, i=NPTCLS_RT,1,-1)]
        call master%new(params, build, pinds)
        call master%read_tabs_to_glob(string('prob_tab_tester_part'), NPARTS_RT, 1)
        nwrong = 0
        do i = 1,NPTCLS_RT
            i_glob = NPTCLS_RT - i + 1
            if( any(master%seed_shifts(:,i_glob) /= seeds(:,i)) .or. .not. master%seed_has_sh(i_glob) ) nwrong = nwrong + 1
            do ri = 1,nrefs
                cand = master%dense_candidate(ri, i_glob)
                if( cand%iref /= ri .or. cand%dist /= rows(ri,i)%dist .or. cand%inpl /= rows(ri,i)%inpl )then
                    nwrong = nwrong + 1
                else if( refined(ri,i) )then
                    if( .not. cand%has_sh .or. cand%x /= shifts(1,ri,i) .or. cand%y /= shifts(2,ri,i) ) nwrong = nwrong + 1
                else if( cand%has_sh .or. cand%x /= 0. .or. cand%y /= 0. )then
                    nwrong = nwrong + 1
                endif
            enddo
        enddo
        call assert_int(0, nwrong, 'the dense part files read back exactly')
        call master%kill
        call kill_round_trip_parts(params, build, NPARTS_RT)
    end subroutine test_dense_part_round_trip

    !> The likelihood sampler over an array of distances draws exactly as the one over a distance function, from the same
    !! seed: 200 random sets of 1 to 300 distances with ties, a flat set, sentinel distances and every sample size
    subroutine test_likelihood_samplers()
        real,    allocatable :: dists(:), pvec(:), pvec2(:), cdf(:)
        integer, allocatable :: inds(:), inds2(:)
        real    :: dist1, corr1, dist2, corr2
        integer :: iset, n, nsample, which1, which2, i, nwrong
        call seed_rnd_fixed(SEED_TAB)
        nwrong = 0
        do iset = 1,200
            n = 1 + int(299. * ran3())
            allocate(dists(n), pvec(n), pvec2(n), cdf(n), inds(n), inds2(n))
            do i = 1,n
                dists(i) = 10. * ran3()
            enddo
            if( n > 4 ) dists(n/2:n/2+2) = dists(1)
            if( mod(iset,10) == 0 ) dists = 3.
            if( mod(iset,13) == 0 ) dists(1:(n+1)/2) = huge(1.)
            nsample = 1 + mod(iset-1, n)
            call seed_rnd_fixed(SEED_DRAW + iset)
            call sample_likelihood_dist(n, dist_at, nsample, dist1, corr1, which1, pvec, inds)
            call seed_rnd_fixed(SEED_DRAW + iset)
            call sample_likelihood_dists(n, dists, nsample, dist2, corr2, which2, pvec2, inds2, cdf)
            if( which1 /= which2 .or. dist1 /= dist2 .or. corr1 /= corr2 ) nwrong = nwrong + 1
            deallocate(dists, pvec, pvec2, cdf, inds, inds2)
        enddo
        call assert_int(0, nwrong, 'the array likelihood sampler draws exactly as the function sampler')

    contains

        real function dist_at( i )
            integer, intent(in) :: i
            dist_at = dists(i)
        end function dist_at

    end subroutine test_likelihood_samplers

    !> The parallel read of the dense part files at one and at four threads, over 1, 3 and 7 parts: the same rows,
    !! refined lists and seeds, bit for bit
    subroutine test_dense_read_threads()
        integer, parameter :: NSPACE_TH = 40, NSTATES_TH = 3, NPTCLS_TH = 70
        integer, parameter :: NPARTS_TH(3) = [1, 3, 7]
        class(parameters), allocatable, target :: params
        class(builder),    allocatable, target :: build
        type(eul_prob_tab) :: master1, master4
        type(prob_dense_entry), allocatable :: rows(:,:)
        real,                   allocatable :: shifts(:,:,:), seeds(:,:)
        logical,                allocatable :: refined(:,:)
        integer,                allocatable :: pinds(:)
        integer :: i, k, nthr_prev, nwrong
        character(len=32) :: label
        nthr_prev = 1
        !$ nthr_prev = omp_get_max_threads()
        do k = 1,size(NPARTS_TH)
            write(label,'(I0,A)') NPARTS_TH(k), ' parts'
            call make_round_trip_parts(params, build, NSPACE_TH, NSTATES_TH, NPTCLS_TH, NPARTS_TH(k), rows, shifts,&
                &refined, seeds)
            ! the master's particles in a scrambled order
            pinds = [(mod(17*i, NPTCLS_TH) + 1, i=0,NPTCLS_TH-1)]
            !$ call omp_set_num_threads(1)
            call master1%new(params, build, pinds)
            call master1%read_tabs_to_glob(string('prob_tab_tester_part'), NPARTS_TH(k), 1)
            !$ call omp_set_num_threads(4)
            call master4%new(params, build, pinds)
            call master4%read_tabs_to_glob(string('prob_tab_tester_part'), NPARTS_TH(k), 1)
            !$ call omp_set_num_threads(nthr_prev)
            nwrong = count(master1%rows(:,:)%dist /= master4%rows(:,:)%dist) +&
                &count(master1%rows(:,:)%inpl /= master4%rows(:,:)%inpl) +&
                &count(master1%refined_offsets /= master4%refined_offsets) +&
                &count(master1%seed_shifts /= master4%seed_shifts) + count(master1%seed_has_sh .neqv. master4%seed_has_sh)
            if( size(master1%refined) /= size(master4%refined) )then
                nwrong = nwrong + 1
            else
                nwrong = nwrong + count(master1%refined(:)%i /= master4%refined(:)%i) +&
                    &count(master1%refined(:)%ri /= master4%refined(:)%ri) +&
                    &count(master1%refined(:)%x /= master4%refined(:)%x) + count(master1%refined(:)%y /= master4%refined(:)%y)
            endif
            call assert_int(0, nwrong, 'the read at four threads equals the read at one ('//trim(label)//')')
            call assert_int(36, master4%seed_nrots, 'the rotation grid of the parts is taken ('//trim(label)//')')
            call master1%kill
            call master4%kill
            call kill_round_trip_parts(params, build, NPARTS_TH(k))
        enddo
    end subroutine test_dense_read_threads

    !> NPARTS_RT dense part files prob_tab_tester_part<n>.dat of consecutive particles, written by worker tables in
    !! two batches each: random rows, a refined shift on about one entry in ten, a seed per particle
    subroutine make_round_trip_parts( params, build, nspace, nstates, nptcls, nparts, rows, shifts, refined, seeds )
        class(parameters), allocatable, target, intent(inout) :: params
        class(builder),    allocatable, target, intent(inout) :: build
        integer,                                intent(in)    :: nspace, nstates, nptcls, nparts
        type(prob_dense_entry), allocatable,    intent(out)   :: rows(:,:)
        real,                   allocatable,    intent(out)   :: shifts(:,:,:), seeds(:,:)
        logical,                allocatable,    intent(out)   :: refined(:,:)
        type(eul_prob_tab) :: worker
        type(prob_refined_entry), allocatable :: batch_refined(:)
        integer, allocatable :: parts(:,:)
        integer :: nrefs, i, ri, ipart, i_loc, nloc, ibatch, first, last, n
        nrefs = nspace * nstates
        allocate(params, build)
        params%nspace      = nspace
        params%nstates     = nstates
        params%prob_athres = 10.
        params%l_doshift   = .true.
        call seed_rnd_fixed(SEED_TAB)
        allocate(build%spproj_field)
        call build%spproj_field%new(nptcls, is_ptcl=.true.)
        do i = 1,nptcls
            call build%spproj_field%set_state(i, mod(i-1, nstates) + 1)
        enddo
        allocate(rows(nrefs,nptcls), shifts(2,nrefs,nptcls), refined(nrefs,nptcls), seeds(2,nptcls))
        do i = 1,nptcls
            seeds(:,i) = [ran3() - 0.5, ran3() - 0.5]
            do ri = 1,nrefs
                rows(ri,i)     = prob_dense_entry(5. * ran3(), 1 + int(36. * ran3()))
                refined(ri,i) = ran3() < 0.1
                shifts(:,ri,i) = [ran3() - 0.5, ran3() - 0.5]
            enddo
        enddo
        parts = split_nobjs_even(nptcls, nparts)
        do ipart = 1,nparts
            nloc = parts(ipart,2) - parts(ipart,1) + 1
            call worker%new_worker(params, build, [(i, i=parts(ipart,1),parts(ipart,2))])
            worker%seed_nrots = 36
            worker%seed_shifts = seeds(:,parts(ipart,1):parts(ipart,2))
            worker%seed_has_sh = .true.
            call worker%begin_write(string('prob_tab_tester_part')//int2str(ipart)//'.dat')
            do ibatch = 1,2
                first = merge(1, nloc/2 + 1, ibatch == 1)
                last  = merge(nloc/2, nloc, ibatch == 1)
                allocate(batch_refined(count(refined(:,parts(ipart,1)+first-1:parts(ipart,1)+last-1))))
                n = 0
                do i_loc = first,last
                    i = parts(ipart,1) + i_loc - 1
                    do ri = 1,nrefs
                        if( .not. refined(ri,i) ) cycle
                        n = n + 1
                        batch_refined(n) = prob_refined_entry(i_loc, ri, shifts(1,ri,i), shifts(2,ri,i))
                    enddo
                enddo
                call worker%write_rows(first, rows(:,parts(ipart,1)+first-1:parts(ipart,1)+last-1), batch_refined)
                deallocate(batch_refined)
            enddo
            call worker%write_tab(string('prob_tab_tester_part')//int2str(ipart)//'.dat')
            call worker%kill
        enddo
    end subroutine make_round_trip_parts

    subroutine kill_round_trip_parts( params, build, nparts )
        class(parameters), allocatable, target, intent(inout) :: params
        class(builder),    allocatable, target, intent(inout) :: build
        integer,                                intent(in)    :: nparts
        integer :: ipart
        do ipart = 1,nparts
            call del_file(string('prob_tab_tester_part')//int2str(ipart)//'.dat')
        enddo
        call build%spproj_field%kill
        deallocate(build%spproj_field)
        deallocate(params, build)
    end subroutine kill_round_trip_parts

    subroutine set_dist( tab, iref, i, dist )
        type(eul_prob_tab), intent(inout) :: tab
        integer,            intent(in)    :: iref, i
        real,               intent(in)    :: dist
        tab%rows(iref,i)%dist = dist
        tab%rows(iref,i)%inpl = 0
    end subroutine set_dist

    !> The assignment ref_assign replaced (as of 2026-10-15), on a dense table dists(nrefs,nptcls) whose
    !! states all exist: every reference sorts the particles; one state draws among all heads and
    !! advances every head after each assignment; several states label every particle by a global
    !! greedy frontier (the minimal head first, minloc over the references), then draw within each
    !! state among its active references, advancing every head of the state after each assignment.
    !! iref_assigned(i) is the reference particle i was assigned.
    subroutine reference_assign( dists, nspace, nstates, prob_athres, state_athres, iref_assigned )
        real,                 intent(in)  :: dists(:,:), prob_athres, state_athres(:)
        integer,              intent(in)  :: nspace, nstates
        integer, allocatable, intent(out) :: iref_assigned(:)
        integer, allocatable :: stab_inds(:,:), inds_sorted(:), iref_dist_inds(:), greedy_state(:)
        integer, allocatable :: active_refs(:), active_inds(:)
        real,    allocatable :: score_work(:), iref_dist(:), dists_sorted(:), active_dists(:), active_dists_sorted(:)
        logical, allocatable :: ptcl_avail(:)
        integer :: nrefs, nptcls, i, iref, assigned_iref, assigned_ptcl, istate, active_idx, nactive
        real    :: projs_athres, dist_tmp, corr_tmp
        nrefs  = size(dists,1)
        nptcls = size(dists,2)
        allocate(stab_inds(nptcls,nrefs), score_work(nptcls), inds_sorted(nrefs), iref_dist_inds(nrefs),&
            &greedy_state(nptcls), active_refs(nrefs), active_inds(nrefs), iref_dist(nrefs), dists_sorted(nrefs),&
            &active_dists(nrefs), active_dists_sorted(nrefs), ptcl_avail(nptcls), iref_assigned(nptcls))
        iref_assigned = 0
        do iref = 1,nrefs
            do i = 1,nptcls
                stab_inds(i,iref) = i
                score_work(i)     = ref_score(iref, i)
            enddo
            call hpsort(score_work, stab_inds(:,iref))
        enddo
        projs_athres = maxval(state_athres)
        greedy_state = 0
        if( nstates > 1 )then
            call reset_ref_frontier()
            do while( any(ptcl_avail) )
                assigned_iref = minloc(iref_dist, dim=1)
                assigned_ptcl = stab_inds(iref_dist_inds(assigned_iref), assigned_iref)
                greedy_state(assigned_ptcl) = ref_state(assigned_iref)
                ptcl_avail(assigned_ptcl)   = .false.
                do iref = 1,nrefs
                    call advance_ref_head(iref, 0)
                enddo
            enddo
            do istate = 1,nstates
                call reset_ref_frontier()
                ptcl_avail = greedy_state == istate
                do iref = 1,nrefs
                    if( ref_state(iref) == istate ) call advance_ref_head(iref, istate)
                enddo
                do while( any(ptcl_avail) )
                    nactive = 0
                    do iref = 1,nrefs
                        if( ref_state(iref) /= istate ) cycle
                        if( iref_dist_inds(iref) > nptcls ) cycle
                        nactive = nactive + 1
                        active_refs(nactive)  = iref
                        active_dists(nactive) = iref_dist(iref)
                    enddo
                    if( nactive == 1 )then
                        active_idx = 1
                    else
                        call sample_likelihood_dist(nactive, active_ref_dist, likelihood_nsample(nactive, state_athres(istate)),&
                            &dist_tmp, corr_tmp, active_idx, active_dists_sorted, active_inds)
                    endif
                    assigned_iref = active_refs(active_idx)
                    call assign_current_ref()
                    do iref = 1,nrefs
                        if( ref_state(iref) == istate ) call advance_ref_head(iref, istate)
                    enddo
                enddo
            enddo
        else
            call reset_ref_frontier()
            do while( any(ptcl_avail) )
                call sample_likelihood_dist(nrefs, frontier_ref_dist, likelihood_nsample(nrefs, projs_athres),&
                    &dist_tmp, corr_tmp, assigned_iref, dists_sorted, inds_sorted)
                call assign_current_ref()
                do iref = 1,nrefs
                    call advance_ref_head(iref, 0)
                enddo
            enddo
        endif

    contains

        real function ref_score( iref_loc, iptcl_loc ) result( score )
            integer, intent(in) :: iref_loc, iptcl_loc
            real :: invalid_dist
            invalid_dist = 0.1 * huge(invalid_dist)
            score = dists(iref_loc,iptcl_loc)
            if( .not.(ieee_is_finite(score) .and. score < invalid_dist) ) score = invalid_dist
        end function ref_score

        integer function ref_state( iref_loc )
            integer, intent(in) :: iref_loc
            ref_state = (iref_loc - 1) / nspace + 1
        end function ref_state

        integer function likelihood_nsample( n, athres_ub_in ) result( nsample )
            integer, intent(in) :: n
            real,    intent(in) :: athres_ub_in
            nsample = min(n, max(1, int(min(prob_athres, athres_ub_in) * real(n) / 180.)))
        end function likelihood_nsample

        real function frontier_ref_dist( iref_loc ) result( dist )
            integer, intent(in) :: iref_loc
            dist = iref_dist(iref_loc)
        end function frontier_ref_dist

        real function active_ref_dist( active_loc ) result( dist )
            integer, intent(in) :: active_loc
            dist = active_dists(active_loc)
        end function active_ref_dist

        subroutine reset_ref_frontier()
            integer :: iref_loc
            do iref_loc = 1,nrefs
                iref_dist(iref_loc) = ref_score(iref_loc, stab_inds(1,iref_loc))
            enddo
            iref_dist_inds = 1
            ptcl_avail     = .true.
        end subroutine reset_ref_frontier

        subroutine advance_ref_head( iref_current, state_filter )
            integer, intent(in) :: iref_current, state_filter
            integer :: iptcl_current
            do while( iref_dist_inds(iref_current) <= nptcls )
                iptcl_current = stab_inds(iref_dist_inds(iref_current), iref_current)
                if( ptcl_avail(iptcl_current) .and. &
                    &(state_filter == 0 .or. greedy_state(iptcl_current) == state_filter) )then
                    iref_dist(iref_current) = ref_score(iref_current, iptcl_current)
                    return
                endif
                iref_dist_inds(iref_current) = iref_dist_inds(iref_current) + 1
            enddo
            iref_dist(iref_current) = huge(iref_dist(iref_current))
        end subroutine advance_ref_head

        subroutine assign_current_ref()
            assigned_ptcl = stab_inds(iref_dist_inds(assigned_iref), assigned_iref)
            ptcl_avail(assigned_ptcl)    = .false.
            iref_assigned(assigned_ptcl) = assigned_iref
        end subroutine assign_current_ref

    end subroutine reference_assign

    subroutine test_record_seed_shift()
        type(ptcl_ref) :: assgn
        integer :: inpl
        write(*,'(A)') 'test_record_seed_shift'
        ! the particle-frame seed is recorded whatever the angle
        do inpl = 1, NROTS
            assgn = unrefined(inpl)
            call record_seed_shift(assgn, SEED, .true., .true.)
            call assert_true(assgn%has_sh .and. assgn%x == SEED(1) .and. assgn%y == SEED(2), &
                &'the seed is recorded as it is, at every in-plane angle')
        enddo
        ! a refined candidate keeps its own shift, a particle without seed or without shift search is left alone
        assgn        = unrefined(3)
        assgn%has_sh = .true.
        assgn%x      = 0.25
        assgn%y      = -0.75
        call record_seed_shift(assgn, SEED, .true., .true.)
        call assert_true(assgn%x == 0.25 .and. assgn%y == -0.75, 'a refined candidate keeps the shift of its search')
        assgn = unrefined(3)
        call record_seed_shift(assgn, SEED, .false., .true.)
        call assert_false(assgn%has_sh, 'no seed, nothing recorded')
        call record_seed_shift(assgn, SEED, .true., .false.)
        call assert_false(assgn%has_sh, 'no shift search, nothing recorded')
    end subroutine test_record_seed_shift

    function unrefined( inpl ) result( assgn )
        integer, intent(in) :: inpl
        type(ptcl_ref) :: assgn
        assgn%inpl = inpl
    end function unrefined

end module simple_eul_prob_tab_tester
