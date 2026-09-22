!@descr: unit test routines for the oris class
module simple_oris_tester
use simple_core_module_api
use simple_test_utils    ! for assert_* utilities and counters
implicit none
private
public :: run_all_oris_tests

contains

    subroutine run_all_oris_tests()
        write(*,'(A)') '**** running all oris tests ****'
        call test_constructors_and_basic_props()
        call test_getters_setters()
        call test_extract_and_copy()
        call test_compress_and_masks()
        call test_sampling_and_updatecnt()
        call test_randomization_and_symmetry()
        call test_proj_space_and_remap()
        call test_stats_and_ordering()
        call test_rotations_and_errors()
        call test_misc_flags()
        call test_reseed_classes()
        call test_reallocate()
        call test_write_read_roundtrip()
        call test_rnd_oris_bounds()
        ! call report_summary()
    end subroutine run_all_oris_tests

    !---------------------------------------------------------------
    ! 1. Constructors and basic properties
    !---------------------------------------------------------------
    subroutine test_constructors_and_basic_props()
        type(oris) :: os, os2
        type(ori)  :: o
        integer    :: n
        write(*,'(A)') 'test_constructors_and_basic_props'
        n = 5
        call os%new(n, .true.)
        call assert_int(n, os%get_noris(), 'get_noris after new()')
        call assert_true(os%is_particle(), 'is_particle() true for is_ptcl=.true.')
        ! constructor interface oris(n,is_ptcl)
        os2 = oris(n, .false.)
        call assert_int(n, os2%get_noris(), 'oris constructor: size')
        call assert_true(.not. os2%is_particle(), 'oris constructor: is_particle() false')
        ! new_2: from array of ori
        call o%new_ori(.true.)
        call o%e1set(10.0)
        call o%e2set(20.0)
        call o%e3set(30.0)
        call os2%new([o, o])
        call assert_int(2, os2%get_noris(), 'new_2 size from ori array')
        call assert_real(10.0, os2%e1get(1), 1.0e-6, 'new_2 e1get')
        call os%kill
        call os2%kill
        call o%kill
    end subroutine test_constructors_and_basic_props

    !---------------------------------------------------------------
    ! 2. Getters / setters on scalars
    !---------------------------------------------------------------
    subroutine test_getters_setters()
        type(oris) :: os
        integer    :: n, i
        real       :: euls(3), sh(2)
        write(*,'(A)') 'test_getters_setters'
        n = 4
        call os%new(n, .true.)
        ! Initialize some per-ori data
        do i = 1, n
            euls = [real(i*10), real(i*20), real(i*30)]
            call os%set_euler(i, euls)
            sh   = [real(i), real(-i)]
            call os%set_shift(i, sh)
            call os%set_state(i, i)         ! 1..n
            call os%set_class(i, n+1-i)     ! n..1
            call os%set(i, 'proj', real(i))
            call os%set(i, 'eo',   mod(i,2))
            call os%set(i, 'updatecnt', 0.0)
            call os%set(i, 'sampled',   0.0)
            call os%set(i, 'corr',      real(i)/10.0)
        end do
        ! Basic Euler getters
        call assert_real(20.0, os%e1get(2), 1.0e-4, 'e1get')
        call assert_real(60.0, os%e2get(3), 1.0e-4, 'e2get')
        call assert_real(120.0,os%e3get(4), 1.0e-4, 'e3get')
        euls = os%get_euler(1)
        call assert_real(10.0, euls(1), 1.0e-4, 'get_euler(1)%e1')
        call assert_real(20.0, euls(2), 1.0e-4, 'get_euler(1)%e2')
        call assert_real(30.0, euls(3), 1.0e-4, 'get_euler(1)%e3')
        ! Shifts
        sh = os%get_2Dshift(3)
        call assert_real( 3.0, sh(1), 1.0e-4, 'get_2Dshift x')
        call assert_real(-3.0, sh(2), 1.0e-4, 'get_2Dshift y')
        ! State / class / proj / eo
        call assert_int(2, os%get_state(2), 'get_state')
        call assert_int(4, os%get_class(1), 'get_class')
        call assert_int(3, os%get_proj(3),  'get_proj')
        call assert_int(0, os%get_eo(2),    'get_eo')
        ! Simple aggregated accessors
        call assert_int(4, os%get_n('class'), 'get_n(class)')
        call assert_int(1, os%get_pop(1,'state'), 'get_pop single')
        call os%kill
    end subroutine test_getters_setters

    !---------------------------------------------------------------
    ! 3. extract_subset, copy, append, delete
    !---------------------------------------------------------------
    subroutine test_extract_and_copy()
        type(oris) :: os, sub1, sub2, os_copy, os2
        integer    :: n, i
        integer, allocatable :: inds(:)
        write(*,'(A)') 'test_extract_and_copy'
        n = 6
        call os%new(n, .true.)
        call os%set_all('state', [(real(i), i=1,n)])
        ! extract_subset(range)
        sub1 = os%extract_subset(2,4)
        call assert_int(3, sub1%get_noris(), 'extract_subset(range) size')
        call assert_int(2, sub1%get_state(1), 'extract_subset(range) state(1)')
        ! extract_subset(indices)
        inds = [1,3,6]
        sub2 = os%extract_subset(inds)
        call assert_int(3, sub2%get_noris(), 'extract_subset(indices) size')
        call assert_int(3, sub2%get_state(2), 'extract_subset(indices) state(2)')
        ! copy_2
        call os_copy%copy(os)
        call assert_int(n, os_copy%get_noris(), 'copy_2 size')
        call assert_int(5, os_copy%get_state(5), 'copy_2 content')
        ! append(oris)
        os2 = oris(2, .true.)
        call os2%set_all('state', [10.0, 11.0])
        call os%append(os2)
        call assert_int(n+2, os%get_noris(), 'append_2 size')
        call assert_int(11, os%get_state(n+2), 'append_2 content')
        ! delete (single entry)
        call os%delete(3)
        call assert_int(n+1, os%get_noris(), 'delete size')
        call os%kill
        call os2%kill
        call os_copy%kill
        call sub1%kill
        call sub2%kill
        if (allocated(inds)) deallocate(inds)
    end subroutine test_extract_and_copy

    !---------------------------------------------------------------
    ! 4. compress, masks, get_all, get_all_sampled, included
    !---------------------------------------------------------------
    subroutine test_compress_and_masks()
        type(oris) :: os
        logical, allocatable :: mask(:), incl(:)
        real,    allocatable :: arr(:), arrs(:)
        integer, allocatable :: pinds(:)
        integer :: n
        write(*,'(A)') 'test_compress_and_masks'
        n = 5
        call os%new(n, .true.)
        call os%set_all('state', [1.0, 0.0, 1.0, 0.0, 1.0])
        call os%set_all('corr',  [0.1, 0.2, 0.3, 0.4, 0.5])
        incl = os%included()
        call assert_int(3, count(incl), 'included() count')
        mask = incl
        call os%compress(mask)
        call assert_int(3, os%get_noris(), 'compress() size')
        ! get_all
        arr = os%get_all('corr')
        call assert_int(3, size(arr), 'get_all size after compress')
        ! get_all_sampled: fabricate situation
        call os%set_all('state',   [1.0,1.0,1.0])
        call os%set_all('sampled', [1.0,0.0,1.0])
        arrs = os%get_all_sampled('corr')
        call assert_int(2, size(arrs), 'get_all_sampled default')
        call os%mask_from_state(1, mask, pinds)
        call assert_int(3, size(pinds), 'mask_from_state')
        call os%kill
        if (allocated(mask))  deallocate(mask)
        if (allocated(incl))  deallocate(incl)
        if (allocated(arr))   deallocate(arr)
        if (allocated(arrs))  deallocate(arrs)
        if (allocated(pinds)) deallocate(pinds)
    end subroutine test_compress_and_masks

    !---------------------------------------------------------------
    ! 5. Sampling and updatecnt / sampled related methods
    !---------------------------------------------------------------
    subroutine test_sampling_and_updatecnt()
        type(oris) :: os
        type(class_sample), allocatable :: clssmp(:)
        integer, allocatable :: inds(:)
        integer :: i, n, nsamp
        real    :: frac, updfrac
        write(*,'(A)') 'test_sampling_and_updatecnt'
        n = 10
        call os%new(n, .true.)
        call os%set_all2single('state', 1.0)
        call os%set_all2single('updatecnt', 0.0)
        call os%set_all2single('sampled',   0.0)
        ! sample4update_all
        call os%sample4update_all([1, n], nsamp, inds, .true.)
        call assert_int(n, nsamp, 'sample4update_all nsamp')
        call assert_true(os%has_been_sampled(), 'has_been_sampled after sample4update_all')
        ! sample4update_rnd
        frac = 0.5
        call os%set_all2single('updatecnt', 0.0)
        call os%sample4update_rnd([1, n], frac, nsamp, inds, .true.)
        call assert_true(nsamp <= n, 'sample4update_rnd nsamp<=n')
        ! sample4update_cnt exhausts lower updatecnt tiers and samples
        ! uniformly only within the cutoff tier.
        call os%set_all2single('updatecnt', 1.0)
        do i = 1, 3
            call os%set(i, 'updatecnt', 0.)
        end do
        do i = 9, 10
            call os%set(i, 'updatecnt', 2.)
        end do
        call os%sample4update_cnt([1, n], frac, nsamp, inds, .true.)
        call assert_int(5, nsamp, 'sample4update_cnt nsamp')
        call assert_int(3, count(inds <= 3), 'sample4update_cnt exhausts lowest tier')
        call assert_int(2, count(inds >= 4 .and. inds <= 8), 'sample4update_cnt cutoff tier')
        call assert_int(0, count(inds >= 9), 'sample4update_cnt excludes higher tier')
        ! sample4update_class prefers never-updated particles inside each class quota
        allocate(clssmp(2))
        clssmp(1)%clsind = 1
        clssmp(1)%pop    = 5
        allocate(clssmp(1)%pinds(5), source=[1,2,3,4,5])
        allocate(clssmp(1)%ccs(5),   source=[5.,4.,3.,2.,1.])
        clssmp(2)%clsind = 2
        clssmp(2)%pop    = 5
        allocate(clssmp(2)%pinds(5), source=[6,7,8,9,10])
        allocate(clssmp(2)%ccs(5),   source=[5.,4.,3.,2.,1.])
        call os%set_all2single('updatecnt', 2.0)
        call os%set(4,  'updatecnt', 0.)
        call os%set(5,  'updatecnt', 0.)
        call os%set(9,  'updatecnt', 0.)
        call os%set(10, 'updatecnt', 0.)
        frac = 0.4
        call os%sample4update_class(clssmp, [1, n], frac, nsamp, inds, .true., .false.)
        call assert_int(4, nsamp, 'sample4update_class nsamp')
        if( size(inds) == 4 )then
            call assert_true(all(inds == [4,5,9,10]), 'sample4update_class low updatecnt first')
        else
            call assert_true(.false., 'sample4update_class low updatecnt first')
        endif
        ! sampled_only keeps class-balanced updates inside a persistent cohort
        call os%set_all2single('sampled',   0.0)
        call os%set_all2single('updatecnt', 0.0)
        do i = 1, 10
            if( any(i == [1,4,5,6,9,10]) )then
                call os%set(i, 'sampled',   1.)
                call os%set(i, 'updatecnt', 1.)
            endif
        enddo
        frac = 0.3
        call os%sample4update_class(clssmp, [1, n], frac, nsamp, inds, .true., .false., sampled_only=.true.)
        call assert_int(3, nsamp, 'sample4update_class sampled_only exact nsamp')
        do i = 1, size(inds)
            call assert_true(any(inds(i) == [1,4,5,6,9,10]), 'sample4update_class sampled_only cohort')
        enddo
        do i = 1, size(clssmp)
            if( allocated(clssmp(i)%pinds) ) deallocate(clssmp(i)%pinds)
            if( allocated(clssmp(i)%ccs)   ) deallocate(clssmp(i)%ccs)
        end do
        deallocate(clssmp)
        ! sample4update_missing selects only active never-updated particles
        call os%set_all2single('state', 1.0)
        call os%set_all2single('updatecnt', 1.0)
        call os%set(4,  'updatecnt', 0.)
        call os%set(5,  'updatecnt', 0.)
        call os%set(10, 'updatecnt', 0.)
        call os%set(10, 'state',     0.)
        call os%sample4update_missing([1, n], nsamp, inds, .true.)
        call assert_int(2, nsamp, 'sample4update_missing active missing size')
        if( size(inds) == 2 )then
            call assert_true(all(inds == [4,5]), 'sample4update_missing active missing inds')
        else
            call assert_true(.false., 'sample4update_missing active missing inds')
        endif
        call assert_int(1, os%get_updatecnt(4),  'sample4update_missing updatecnt 4')
        call assert_int(0, os%get_updatecnt(10), 'sample4update_missing inactive skip')
        call os%set_all2single('state', 1.0)
        call os%set_all2single('updatecnt', 1.0)
        call os%sample4update_missing([1, n], nsamp, inds, .true.)
        call assert_int(0, nsamp,     'sample4update_missing no-op size')
        call assert_int(0, size(inds),'sample4update_missing no-op inds')
        ! sample4update_updated: requires some updatecnt>0
        call os%set_all2single('updatecnt', 0.0)
        call os%set_updatecnt(1, [1,2,3])
        call os%sample4update_updated([1, n], nsamp, inds, .true.)
        call assert_int(3, nsamp, 'sample4update_updated size')
        ! update fraction
        updfrac = os%get_update_frac()
        call assert_true(updfrac > 0.0, 'get_update_frac>0')
        call os%kill
        if (allocated(inds)) deallocate(inds)
        call test_sample4update_cnt_large()
    end subroutine test_sampling_and_updatecnt

    !---------------------------------------------------------------
    ! Large-population regression and timing test.
    !---------------------------------------------------------------
    subroutine test_sample4update_cnt_large()
        integer, parameter :: n = 500001
        real,    parameter :: update_frac = 0.005
        real,    parameter :: zero_frac   = 0.9 * update_frac
        type(oris) :: os
        integer, allocatable :: candidate_inds(:), zero_inds(:), inds(:)
        integer :: i, nzero, nsamples, nexpected
        integer :: nzero_sampled, nselected_sampled
        integer(timer_int_kind) :: tstart
        real(timer_int_kind)    :: elapsed
        write(*,'(A)') 'test_sample4update_cnt_large'
        nzero     = nint(zero_frac * real(n))
        nexpected = min(n, max(1, nint(update_frac * real(n))))
        call os%new(n, .true.)
        call os%set_all2single('state',     1)
        call os%set_all2single('updatecnt', 1)
        call os%set_all2single('sampled',   0)
        ! Randomly position an exact 0.9 * update_frac population at updatecnt=0.
        allocate(candidate_inds(n))
        do i = 1, n
            candidate_inds(i) = i
        end do
        call partial_shuffle(candidate_inds, nzero)
        allocate(zero_inds(nzero), source=candidate_inds(:nzero))
        deallocate(candidate_inds)
        do i = 1, nzero
            call os%set(zero_inds(i), 'updatecnt', 0)
        end do
        tstart = tic()
        call os%sample4update_cnt([1, n], update_frac, nsamples, inds, .true.)
        elapsed           = toc(tstart)
        nzero_sampled     = 0
        nselected_sampled = 0
        do i = 1, nzero
            if( os%get_sampled(zero_inds(i)) > 0 ) nzero_sampled = nzero_sampled + 1
        end do
        do i = 1, size(inds)
            if( os%get_sampled(inds(i)) > 0 ) nselected_sampled = nselected_sampled + 1
        end do
        call assert_int(nexpected, nsamples, 'large sample4update_cnt sample count')
        call assert_int(nexpected, size(inds), 'large sample4update_cnt returned index count')
        call assert_int(nexpected, nselected_sampled, 'large sample4update_cnt marks every selected particle sampled')
        call assert_int(nzero, nzero_sampled, 'large sample4update_cnt samples every updatecnt=0 particle')
        write(*,'(A,I0)')       '  N:                         ', n
        write(*,'(A,F7.4)')     '  update_frac:               ', update_frac
        write(*,'(A,I0)')       '  updatecnt=0 particles:     ', nzero
        write(*,'(A,I0,A,I0)')  '  zero-count coverage:       ', nzero_sampled, '/', nzero
        call os%kill
        if( allocated(inds) )      deallocate(inds)
        if( allocated(zero_inds) ) deallocate(zero_inds)
    end subroutine test_sample4update_cnt_large

    !---------------------------------------------------------------
    ! 6. Randomization helpers and symmetry / merge / partition_eo
    !---------------------------------------------------------------
    subroutine test_randomization_and_symmetry()
        type(oris) :: os, os2
        integer    :: n
        write(*,'(A)') 'test_randomization_and_symmetry'
        n = 8
        call os%new(n, .false.)
        call os%rnd_oris()
        call os%rnd_inpls()
        call os%rnd_states(2)
        call os%rnd_lps()
        call os%rnd_corrs()
        call os%partition_eo()
        call assert_int(os%get_noris(), os%get_nevenodd(), 'partition_eo nevenodd==n')
        ! symmetrize
        os2 = os
        call os2%symmetrize(3)
        call assert_int(3*n, os2%get_noris(), 'symmetrize size')
        ! merge
        call os%merge(os2)
        call assert_int(n+3*n, os%get_noris(), 'merge sizes')
        call os%kill
    end subroutine test_randomization_and_symmetry

    !---------------------------------------------------------------
    ! 7. Projection space, remap_projs, proj2class
    !---------------------------------------------------------------
    subroutine test_proj_space_and_remap()
        type(oris) :: os_ptcl, e_space
        integer    :: n, ne, i
        integer, allocatable :: mapped(:)
        write(*,'(A)') 'test_proj_space_and_remap'
        n  = 6
        ne = 3
        call os_ptcl%new(n, .true.)
        call e_space%new(ne, .false.)
        ! Some arbitrary euler setup
        do i=1,ne
            call e_space%set_euler(i, [real(i*30), 45.0, 0.0])
        end do
        do i=1,n
            call os_ptcl%set_euler(i, [real(i*30), 45.0, 0.0])
        end do
        call os_ptcl%set_projs(e_space)
        call os_ptcl%proj2class()
        call assert_true(os_ptcl%isthere('proj'), 'proj exists after set_projs')
        allocate(mapped(n))
        call os_ptcl%remap_projs(e_space, mapped)
        call assert_int(n, size(mapped), 'remap_projs size')
        call os_ptcl%kill
        call e_space%kill
        deallocate(mapped)
    end subroutine test_proj_space_and_remap

    !---------------------------------------------------------------
    ! 8. Stats, ordering, spiral
    ! (only using the simpler stats interface to avoid extra types)
    !---------------------------------------------------------------
    subroutine test_stats_and_ordering()
        type(oris) :: os
        real       :: ave, sdev, var
        logical    :: err
        integer, allocatable :: inds(:)
        integer, allocatable :: pops(:)
        integer :: n, ncls
        write(*,'(A)') 'test_stats_and_ordering'
        n = 5
        call os%new(n, .true.)
        call os%set_all2single('state', 1.0)
        call os%set_all('corr',  [0.2, 0.4, 0.1, 0.5, 0.3])
        call os%set_all('class', [1,1,2,2,3])
        ! spiral: just exercise the call
        call os%spiral()
        ! stats (simple interface)
        call os%stats('corr', ave, sdev, var, err)
        call assert_true(.not. err, 'stats_1 no error')
        call os%stats('corr', ave, sdev, var, err, fromto=[2,4])
        call assert_true(.not. err, 'stats_2 no error with fromto')
        ! order / order_cls
        inds = os%order()
        call assert_int(n, size(inds), 'order size')
        ncls = os%get_n('class')
        inds = os%order_cls(ncls)
        call assert_int(ncls, size(inds), 'order_cls size')
        ! get_pops
        call os%get_pops(pops, 'class')
        call assert_int(ncls, size(pops), 'get_pops size')
        call os%kill
        if (allocated(inds)) deallocate(inds)
        if (allocated(pops)) deallocate(pops)
    end subroutine test_stats_and_ordering

    !---------------------------------------------------------------
    ! 9. Rotations, alignment / CTF error introduction
    !---------------------------------------------------------------
    subroutine test_rotations_and_errors()
        type(oris) :: os
        type(ori)  :: e
        integer    :: n
        write(*,'(A)') 'test_rotations_and_errors'
        n = 4
        call os%new(n, .false.)
        call os%set_all2single('state', 1.0)
        call os%spiral()
        call e%new_ori(.false.)
        call e%set_euler([15.0, 30.0, 45.0])
        call os%rot(e)
        call os%rot_transp(e)
        call os%introd_alig_err(5.0, 2.0)
        call os%introd_ctf_err(500.0)
        call os%kill
        call e%kill
    end subroutine test_rotations_and_errors

    !---------------------------------------------------------------
    ! 10. Misc flags / utility methods
    !---------------------------------------------------------------
    subroutine test_misc_flags()
        type(oris) :: os
        integer    :: n
        write(*,'(A)') 'test_misc_flags'
        n = 3
        call os%new(n, .true.)
        call os%set_all('state', [1.0,0.0,1.0])
        call assert_true(os%any_state_zero(), 'any_state_zero()')
        call assert_int(2, os%count_state_gt_zero(), 'count_state_gt_zero')
        ! zero/zero_* utilities
        call os%set_all('x', [1.0,2.0,3.0])
        call os%zero_shifts()
        call assert_real(0.0, os%get(1,'x'), 1.0e-6, 'zero_shifts() (x field zeroed)')
        call os%zero_inpl()
        call os%zero('corr')
        ! revshsgn / revorisgn
        call os%set_all('x', [1.0,2.0,3.0])
        call os%set_all('y', [4.0,5.0,6.0])
        call os%revshsgn()
        call assert_real(-1.0, os%get(1,'x'), 1.0e-6, 'revshsgn x')
        call os%revorisgn()
        ! delete_2Dclustering/delete_3Dalignment: just exercise calls
        call os%delete_2Dclustering()
        call os%delete_3Dalignment()
        call os%kill
    end subroutine test_misc_flags

    !---------------------------------------------------------------
    ! 11. reseed_classes (abinitio2D cls_init=prev seed partition)
    !---------------------------------------------------------------
    subroutine test_reseed_classes()
        type(oris) :: os
        integer, allocatable :: parent_of_seed(:), seed_pops(:)
        integer :: i, n, ndropped, icls, iseed
        real    :: corr_mean_a, corr_mean_b
        write(*,'(A)') 'test_reseed_classes'
        ! parents: class 1 = 60 ptcls, class 2 = 30, class 3 = 10, class 4 = 4 (rejected),
        ! plus 6 state=0 particles in class 1 that must be ignored
        n = 110
        call os%new(n, .true.)
        call os%set_all2single('state', 1)
        do i = 1, n
            if( i <= 60 )then
                icls = 1
            else if( i <= 90 )then
                icls = 2
            else if( i <= 100 )then
                icls = 3
            else if( i <= 104 )then
                icls = 4
            else
                icls = 1
                call os%set_state(i, 0)
            endif
            call os%set_class(i, icls)
            call os%set(i, 'corr', real(i) / real(n))
        end do
        ! K = M = 3: identity allocation on the accepted parents, class 4 unassigned
        call os%reseed_classes([1,2,3], 3, parent_of_seed, seed_pops, ndropped)
        call assert_int(0, ndropped, 'reseed K=M: nothing dropped')
        call assert_int(3, size(parent_of_seed), 'reseed K=M: size(parent_of_seed)')
        call assert_int(100, sum(seed_pops), 'reseed K=M: all active accepted particles seeded')
        call assert_int(60, seed_pops(1), 'reseed K=M: parent 1 population preserved')
        call assert_int(10, seed_pops(3), 'reseed K=M: parent 3 population preserved')
        call assert_int(0, os%get_class(101), 'reseed K=M: rejected class -> class 0')
        call assert_int(0, os%get_class(105), 'reseed K=M: state=0 particle -> class 0')
        ! K = 10 > M: seeds proportional to population (6,3,1), children balanced
        do i = 1, n
            if( i <= 60 )then
                icls = 1
            else if( i <= 90 )then
                icls = 2
            else if( i <= 100 )then
                icls = 3
            else
                icls = 4
            endif
            call os%set_class(i, icls)
        end do
        call os%reseed_classes([1,2,3], 10, parent_of_seed, seed_pops, ndropped)
        call assert_int(0,  ndropped, 'reseed K>M: nothing dropped')
        call assert_int(10, size(seed_pops), 'reseed K>M: K seed classes')
        call assert_int(6,  count(parent_of_seed == 1), 'reseed K>M: parent 1 gets 6 seeds')
        call assert_int(3,  count(parent_of_seed == 2), 'reseed K>M: parent 2 gets 3 seeds')
        call assert_int(1,  count(parent_of_seed == 3), 'reseed K>M: parent 3 gets 1 seed')
        call assert_int(100, sum(seed_pops), 'reseed K>M: all active accepted particles seeded')
        call assert_int(10, minval(seed_pops), 'reseed K>M: balanced children (min)')
        call assert_int(10, maxval(seed_pops), 'reseed K>M: balanced children (max)')
        ! rank interleaving: the children of parent 1 have the same corr mean (within one rank step)
        corr_mean_a = 0.
        corr_mean_b = 0.
        do i = 1, 60
            iseed = os%get_class(i)
            call assert_true(iseed >= 1 .and. iseed <= 6, 'reseed K>M: parent-1 particle lands in a parent-1 seed')
            if( iseed == 1 ) corr_mean_a = corr_mean_a + os%get(i, 'corr')
            if( iseed == 6 ) corr_mean_b = corr_mean_b + os%get(i, 'corr')
        end do
        call assert_real(corr_mean_a / 10., corr_mean_b / 10., 6. / real(n), 'reseed K>M: interleaved children share the corr distribution')
        ! K = 2 < M: the smallest parent is dropped and its particles unassigned
        do i = 1, n
            if( i <= 60 )then
                icls = 1
            else if( i <= 90 )then
                icls = 2
            else if( i <= 100 )then
                icls = 3
            else
                icls = 4
            endif
            call os%set_class(i, icls)
        end do
        call os%reseed_classes([1,2,3], 2, parent_of_seed, seed_pops, ndropped)
        call assert_int(1,  ndropped, 'reseed K<M: one parent dropped')
        call assert_int(2,  size(seed_pops), 'reseed K<M: K seed classes')
        call assert_int(1,  parent_of_seed(1), 'reseed K<M: largest parent kept first')
        call assert_int(2,  parent_of_seed(2), 'reseed K<M: second parent kept')
        call assert_int(90, sum(seed_pops), 'reseed K<M: dropped parent particles unassigned')
        call assert_int(0,  os%get_class(95), 'reseed K<M: dropped-parent particle -> class 0')
        deallocate(parent_of_seed, seed_pops)
        call os%kill
    end subroutine test_reseed_classes

    !---------------------------------------------------------------
    ! 12. reallocate
    !---------------------------------------------------------------
    subroutine test_reallocate()
        type(oris) :: os
        integer    :: n, i
        write(*,'(A)') 'test_reallocate'
        n = 4
        call os%new(n, .true.)
        call os%set_all('state', [(real(i), i=1,n)])
        call os%set_all('e1',    [(10.0*real(i), i=1,n)])
        call os%reallocate(10)
        call assert_int(10, os%get_noris(),         'reallocate grows to the requested size')
        call assert_true(os%is_particle(),          'reallocate preserves is_ptcl')
        call assert_int(4, os%get_state(4),         'reallocate preserves existing entries (state)')
        call assert_real(40.0, os%e1get(4), 1.0e-5, 'reallocate preserves existing entries (e1)')
        call assert_true(os%exists(10),             'reallocate: new entries exist')
        call assert_int(0, os%get_state(10),        'reallocate: new entries are blank')
        call os%kill
    end subroutine test_reallocate

    !---------------------------------------------------------------
    ! 13. write / read round-trip through a text orientation file
    !---------------------------------------------------------------
    subroutine test_write_read_roundtrip()
        type(oris)   :: os, os2
        type(string) :: fname, s
        real         :: e1(3), e2(3), sh1(2), sh2(2)
        integer      :: n, i, j, nst
        write(*,'(A)') 'test_write_read_roundtrip'
        n     = 5
        fname = string('oris_tester_roundtrip.txt')
        call del_file(fname)
        call os%new(n, .true.)
        call os%rnd_oris(3.0)
        call os%set_all('state', [1.0, 2.0, 1.0, 2.0, 1.0])
        call os%set_all('corr',  [(0.1*real(i), i=1,n)])
        call os%set(2, 'tag', 'ABC')
        call os%write(fname)
        call assert_true(file_exists(fname), 'write creates the orientation file')
        call assert_int(n, nlines(fname),    'write emits one line per orientation')
        call os2%new(n, .true.)
        call os2%read(fname, nst=nst)
        call assert_int(2, nst,              'read reports the number of states')
        do i = 1,n
            e1  = os%get_euler(i)
            e2  = os2%get_euler(i)
            sh1 = os%get_2Dshift(i)
            sh2 = os2%get_2Dshift(i)
            do j = 1,3
                call assert_real(e1(j), e2(j), 1.0e-3, 'write/read round-trip Euler angle')
            end do
            call assert_real(sh1(1), sh2(1), 1.0e-4, 'write/read round-trip shift x')
            call assert_real(sh1(2), sh2(2), 1.0e-4, 'write/read round-trip shift y')
            call assert_int(os%get_state(i), os2%get_state(i), 'write/read round-trip state')
            call assert_real(os%get(i,'corr'), os2%get(i,'corr'), 1.0e-5, 'write/read round-trip corr')
        end do
        s = os2%get_str(2, 'tag')
        call assert_char('ABC', s%to_char(), 'write/read round-trip char key')
        call del_file(fname)
        call os%kill
        call os2%kill
    end subroutine test_write_read_roundtrip

    !---------------------------------------------------------------
    ! 14. rnd_oris stays within the shift and Euler limits
    !---------------------------------------------------------------
    subroutine test_rnd_oris_bounds()
        type(oris) :: os
        real       :: e(3), sh(2), eullims(3,2)
        integer    :: n, i
        logical    :: shifts_ok, eulers_ok, lims_ok
        write(*,'(A)') 'test_rnd_oris_bounds'
        n = 50
        call os%new(n, .false.)
        call os%rnd_oris(3.0)
        shifts_ok = .true.
        eulers_ok = .true.
        do i = 1,n
            sh = os%get_2Dshift(i)
            e  = os%get_euler(i)
            if( any(abs(sh) > 3.0) ) shifts_ok = .false.
            if( e(1) < 0.0 .or. e(1) > 360.0 ) eulers_ok = .false.
            if( e(2) < 0.0 .or. e(2) > 180.0 ) eulers_ok = .false.
            if( e(3) < 0.0 .or. e(3) > 360.0 ) eulers_ok = .false.
        end do
        call assert_true(shifts_ok, 'rnd_oris: shifts within +/- trs')
        call assert_true(eulers_ok, 'rnd_oris: Euler angles within their canonical ranges')
        eullims(:,1) = [ 0.0,  0.0,   0.0]
        eullims(:,2) = [90.0, 45.0, 360.0]
        call os%rnd_oris(0.0, eullims)
        lims_ok = .true.
        do i = 1,n
            e  = os%get_euler(i)
            sh = os%get_2Dshift(i)
            if( e(1) >= 90.0 .or. e(2) >= 45.0 ) lims_ok = .false.
            if( any(abs(sh) > 0.0) )             lims_ok = .false.
        end do
        call assert_true(lims_ok, 'rnd_oris: Euler limits honoured and trs=0 gives zero shifts')
        call os%kill
    end subroutine test_rnd_oris_bounds

end module simple_oris_tester
