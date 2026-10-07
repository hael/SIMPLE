!@descr: unit test routines for the oris class
module simple_oris_tester
use simple_core_module_api
use simple_test_utils    ! for assert_* utilities and counters
use simple_oris,        only: population_blend_weights, class_sample_quotas, class_sample_sweep
use simple_oris_utils,  only: oridist_from_oris
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
        call test_empty_partition_sampling()
        call test_nested_class_quota()
        call test_class_sample_sweep()
        call test_cohort_rescore()
        call test_sample4rec_global_coverage()
        call test_group_update_counts()
        call test_blend_weight_cases()
        call test_blend_mass_scripted()
        call test_randomization_and_symmetry()
        call test_proj_space_and_remap()
        call test_stats_and_ordering()
        call test_rotations_and_errors()
        call test_misc_flags()
        call test_reseed_classes()
        call test_reallocate()
        call test_write_read_roundtrip()
        call test_rnd_oris_bounds()
        call test_assignment()
        call test_oridist_from_oris()
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

    ! A distributed partition may have nothing to update: its range holds no active particle, or its
    ! active particles get none of a class-balanced sample drawn over the whole project. With
    ! allow_empty the samplers return an empty sample and stamp nothing (without it they stop the run)
    subroutine test_empty_partition_sampling()
        type(oris)                      :: os
        type(class_sample), allocatable :: clssmp(:)
        integer,            allocatable :: inds(:)
        integer                         :: i, nsamp
        write(*,'(A)') 'test_empty_partition_sampling'
        call os%new(10, .true.)
        call os%set_all2single('state',     1.0)
        call os%set_all2single('updatecnt', 0.0)
        call os%set_all2single('sampled',   0.0)
        ! rows 1-6 masked, as the frozen rows of an add-on: the partition [1,6] holds no active particle
        do i = 1, 6
            call os%set_state(i, 0)
        end do
        call os%sample4update_cnt([1, 6], 0.5, nsamp, inds, .true., allow_empty=.true.)
        call assert_int(0, nsamp,      'sample4update_cnt: a partition without active particles samples nothing')
        call assert_int(0, size(inds), 'sample4update_cnt: no indices')
        call os%sample4update_fillin([1, 6], 0.5, nsamp, inds, .true., allow_empty=.true.)
        call assert_int(0, nsamp,      'sample4update_fillin: a partition without active particles samples nothing')
        call assert_int(0, size(inds), 'sample4update_fillin: no indices')
        ! class-balanced: one class holds rows 9 and 10, and half of the 4 active particles is 2, so the
        ! sample is rows 9 and 10; the partition [5,8] has active rows 7 and 8 but none of the sample
        allocate(clssmp(1))
        clssmp(1)%clsind = 1
        clssmp(1)%pop    = 2
        allocate(clssmp(1)%pinds(2), source=[9,10])
        allocate(clssmp(1)%ccs(2),   source=[2.,1.])
        call os%sample4update_class(clssmp, [5, 8], 0.5, nsamp, inds, .true., .false., allow_empty=.true.)
        call assert_int(0, nsamp,      'sample4update_class: a partition outside the global sample samples nothing')
        call assert_int(0, size(inds), 'sample4update_class: no indices')
        ! nothing was stamped: no sampled mark and no update count anywhere
        call assert_false(os%has_been_sampled(), 'an empty sample stamps no particle')
        call assert_true(all([(os%get(i, 'updatecnt') == 0., i=1,10)]), 'an empty sample counts no update')
        deallocate(clssmp(1)%pinds, clssmp(1)%ccs)
        deallocate(clssmp)
        call os%kill
    end subroutine test_empty_partition_sampling

    ! units of consecutive particle indices with the given populations and group labels
    subroutine make_units( pops, groups, clssmp )
        integer,                         intent(in)    :: pops(:), groups(:)
        type(class_sample), allocatable, intent(inout) :: clssmp(:)
        integer :: i, j, first
        if( allocated(clssmp) ) call free_units(clssmp)
        allocate(clssmp(size(pops)))
        first = 1
        do i = 1, size(pops)
            clssmp(i)%clsind = i
            clssmp(i)%pop    = pops(i)
            clssmp(i)%group  = groups(i)
            allocate(clssmp(i)%pinds(pops(i)), source=[(j, j=first,first+pops(i)-1)])
            allocate(clssmp(i)%ccs(pops(i)),   source=1.)
            first = first + pops(i)
        end do
    end subroutine make_units

    subroutine free_units( clssmp )
        type(class_sample), allocatable, intent(inout) :: clssmp(:)
        integer :: i
        do i = 1, size(clssmp)
            if( allocated(clssmp(i)%pinds) ) deallocate(clssmp(i)%pinds)
            if( allocated(clssmp(i)%ccs)   ) deallocate(clssmp(i)%ccs)
        end do
        deallocate(clssmp)
    end subroutine free_units

    ! particles of every unit in a sample
    function unit_counts( clssmp, inds ) result( cnts )
        type(class_sample), intent(in) :: clssmp(:)
        integer,            intent(in) :: inds(:)
        integer :: cnts(size(clssmp)), i, j
        cnts = 0
        do i = 1, size(clssmp)
            do j = 1, size(inds)
                if( any(clssmp(i)%pinds == inds(j)) ) cnts(i) = cnts(i) + 1
            end do
        end do
    end function unit_counts

    ! balance=class: every unit the same count capped at its population, lowest updatecnt first.
    ! balance=cavg: every group the same count capped at its population and, inside a group, every unit
    ! the same count capped at its population; the remainder of a group's equal split goes to its
    ! least-updated unit
    subroutine test_nested_class_quota()
        type(oris)                      :: os
        type(class_sample), allocatable :: clssmp(:)
        integer,            allocatable :: inds(:)
        integer :: nsamp, i, cnts3(3), cnts5(5)
        write(*,'(A)') 'test_nested_class_quota'
        ! class: pops 2, 5, 10; 9 of 17 -> the equal rounds stop at 4 per unit, the first capped at 2
        call os%new(17, .true.)
        call os%set_all2single('state',     1.0)
        call os%set_all2single('updatecnt', 0.0)
        call os%set_all2single('sampled',   0.0)
        call os%set(8,  'updatecnt', 3.)   ! unit 3 (rows 8-17): rows 8 and 9 were updated before
        call os%set(9,  'updatecnt', 3.)
        call make_units([2,5,10], [0,0,0], clssmp)
        call os%sample4update_class(clssmp, [1, 17], 9./17., nsamp, inds, .true., .false.)
        cnts3 = unit_counts(clssmp, inds)
        call assert_int(10, nsamp, 'class: equal rounds overshoot by less than the number of units')
        call assert_true(all(cnts3 == [2,4,4]), 'class: every unit the same count capped at its population')
        call assert_true(.not. any(inds == 8) .and. .not. any(inds == 9), 'class: lowest updatecnt first within a unit')
        call free_units(clssmp)
        call os%kill
        ! cavg: group 1 holds units of 2, 6 and 6 particles, group 2 one unit of 3, group 3 one unit of 20.
        ! 19 of 37: group quotas 8, 3 (capped), 8; inside group 1: 2 (capped), 3, 3
        call os%new(37, .true.)
        call os%set_all2single('state',     1.0)
        call os%set_all2single('updatecnt', 0.0)
        call os%set_all2single('sampled',   0.0)
        call make_units([2,6,6,3,20], [1,1,1,2,3], clssmp)
        call os%sample4update_class(clssmp, [1, 37], 19./37., nsamp, inds, .true., .false.)
        cnts5 = unit_counts(clssmp, inds)
        call assert_int(19, nsamp, 'cavg: sample size')
        call assert_true(all(cnts5 == [2,3,3,3,8]), 'cavg: nested equal quota over groups, then units')
        call assert_int(8, sum(cnts5(1:3)), 'cavg: group 1 gets the equal group count')
        call free_units(clssmp)
        call os%kill
        ! remainder: one group of three units of 4 and 4 draws -> 1 each and the remainder to the unit
        ! whose particles have the lowest mean updatecnt (unit 2)
        call os%new(12, .true.)
        call os%set_all2single('state',     1.0)
        call os%set_all2single('updatecnt', 1.0)
        call os%set_all2single('sampled',   0.0)
        do i = 5, 8
            call os%set(i, 'updatecnt', 0.)
        end do
        call make_units([4,4,4], [1,1,1], clssmp)
        call os%sample4update_class(clssmp, [1, 12], 4./12., nsamp, inds, .true., .false.)
        cnts3 = unit_counts(clssmp, inds)
        call assert_true(all(cnts3 == [1,2,1]), 'cavg: the group remainder goes to the least-updated unit')
        call free_units(clssmp)
        call os%kill
        if( allocated(inds) ) deallocate(inds)
    end subroutine test_nested_class_quota

    ! a fixed skewed population (the cavg layout above) is visited completely in sweep iterations
    subroutine test_class_sample_sweep()
        type(oris)                      :: os
        type(class_sample), allocatable :: clssmp(:)
        integer,            allocatable :: inds(:)
        real    :: quotas(5)
        integer :: nsamp, it, sweep, i
        write(*,'(A)') 'test_class_sample_sweep'
        call make_units([2,6,6,3,20], [1,1,1,2,3], clssmp)
        call class_sample_quotas(clssmp, 19, quotas)
        call assert_real(2., quotas(1), 1.e-5, 'quota of the capped unit is its population')
        call assert_real(3., quotas(2), 1.e-5, 'quota of an open unit of group 1')
        call assert_real(8., quotas(5), 1.e-5, 'quota of the single-unit group 3')
        sweep = class_sample_sweep(clssmp, 19)
        call assert_int(3, sweep, 'sweep = max over units of ceil(pop/quota)')
        call os%new(37, .true.)
        call os%set_all2single('state',     1.0)
        call os%set_all2single('updatecnt', 0.0)
        call os%set_all2single('sampled',   0.0)
        do it = 1, sweep
            call os%sample4update_class(clssmp, [1, 37], 19./37., nsamp, inds, .true., .false.)
        end do
        call assert_true(all([(os%get_updatecnt(i) > 0, i=1,37)]), 'every particle is visited within one sweep')
        call free_units(clssmp)
        call os%kill
        if( allocated(inds) ) deallocate(inds)
    end subroutine test_class_sample_sweep

    ! sample4update_rescore returns the latest round on every call and advances sampled and updatecnt on
    ! exactly those rows once per call; a later draw prefers the rows outside the cohort
    subroutine test_cohort_rescore()
        type(oris)                      :: os
        type(class_sample), allocatable :: clssmp(:)
        integer,            allocatable :: inds(:), inds2(:)
        integer :: nsamp, i, ucnt_before(10), marker
        write(*,'(A)') 'test_cohort_rescore'
        call os%new(10, .true.)
        call os%set_all2single('state',     1.0)
        call os%set_all2single('updatecnt', 1.0)
        call os%set_all2single('sampled',   0.0)
        do i = 1, 3
            call os%set(i, 'updatecnt', 0.)
        end do
        ! the cohort: rows 1-3 (lowest updatecnt tier)
        call os%sample4update_cnt([1, 10], 0.3, nsamp, inds, .true.)
        call assert_true(size(inds) == 3 .and. all(inds == [1,2,3]), 'cohort draw takes the lowest tier')
        marker = nint(os%get(1, 'sampled'))
        do i = 1, 10
            ucnt_before(i) = os%get_updatecnt(i)
        end do
        call os%sample4update_rescore([1, 10], nsamp, inds2)
        call assert_true(size(inds2) == 3, 'rescore returns the cohort size')
        if( size(inds2) == 3 ) call assert_true(all(inds2 == inds), 'rescore returns the cohort')
        call assert_true(all([(nint(os%get(i, 'sampled')) == marker + 1, i=1,3)]), 'rescore advances the round marker once')
        call assert_true(all([(os%get_updatecnt(i) == ucnt_before(i) + 1, i=1,3)]), 'rescore counts one update per cohort row')
        call assert_true(all([(os%get_updatecnt(i) == ucnt_before(i), i=4,10)]), 'rescore leaves the other rows alone')
        call os%sample4update_rescore([1, 10], nsamp, inds2)
        call assert_true(size(inds2) == 3, 'repeated rescore returns the cohort size')
        if( size(inds2) == 3 ) call assert_true(all(inds2 == inds), 'repeated rescore returns the identical cohort')
        call assert_true(all([(os%get_updatecnt(i) == ucnt_before(i) + 2, i=1,3)]), 'one update per rescore call')
        ! a range without cohort rows
        call os%sample4update_rescore([5, 10], nsamp, inds2, allow_empty=.true.)
        call assert_int(0, nsamp,       'rescore on a range without cohort rows samples nothing')
        call assert_int(0, size(inds2), 'rescore on a range without cohort rows: no indices')
        call assert_true(all([(os%get_updatecnt(i) == 1, i=4,10)]), 'an empty rescore counts no update')
        ! the next draws avoid the previous cohort while lower-updatecnt rows remain
        call os%sample4update_cnt([1, 10], 0.3, nsamp, inds2, .true.)
        call assert_true(.not. any(inds2 <= 3), 'global draw after a cohort takes no cohort row')
        call make_units([10], [0], clssmp)
        call os%sample4update_class(clssmp, [1, 10], 0.3, nsamp, inds2, .true., .false.)
        call assert_true(.not. any(inds2 <= 3), 'unit draw after a cohort takes no cohort row')
        call free_units(clssmp)
        call os%kill
    end subroutine test_cohort_rescore

    ! sample4rec decides "nothing updated yet" over the whole project: when only
    ! the range [1,5] holds updated rows, the range [6,10] must return none of its
    ! never-updated rows, and a project without updated rows returns every active row
    subroutine test_sample4rec_global_coverage()
        type(oris)           :: os
        integer, allocatable :: inds(:), pops(:)
        integer              :: i, nsamp
        write(*,'(A)') 'test_sample4rec_global_coverage'
        call os%new(10, .true.)
        call os%set_all2single('state',     1.0)
        call os%set_all2single('updatecnt', 0.0)
        call os%set(2, 'updatecnt', 1.)
        call os%set(4, 'updatecnt', 3.)
        call os%set_state(5, 0)
        call os%set(5, 'updatecnt', 2.)   ! inactive: never selected
        call os%set_state(7, 0)
        call os%sample4rec([1, 5], nsamp, inds)
        call assert_int(2, nsamp, 'sample4rec: updated range selects its updated active rows')
        if( size(inds) == 2 )then
            call assert_true(all(inds == [2,4]), 'sample4rec: updated range indices')
        else
            call assert_true(.false., 'sample4rec: updated range indices')
        endif
        call os%sample4rec([6, 10], nsamp, inds)
        call assert_int(0, nsamp,      'sample4rec: a range without updated rows selects nothing when others are updated')
        call assert_int(0, size(inds), 'sample4rec: no never-updated indices')
        ! the population a full reconstruction represents: rows 2 and 4
        call os%get_state_rec_pops(1, pops)
        call assert_int(2, pops(1), 'state rec pops: updated active rows when any is updated')
        ! no active updated row anywhere (row 5 is inactive): every active row of the range
        call os%set(2, 'updatecnt', 0.)
        call os%set(4, 'updatecnt', 0.)
        call os%sample4rec([6, 10], nsamp, inds)
        call assert_int(4, nsamp, 'sample4rec: nothing updated selects every active row')
        call os%get_state_rec_pops(1, pops)
        call assert_int(8, pops(1), 'state rec pops: every active row when nothing is updated')
        if( size(inds) == 4 )then
            call assert_true(all(inds == [6,8,9,10]), 'sample4rec: nothing updated indices')
        else
            call assert_true(.false., 'sample4rec: nothing updated indices')
        endif
        call assert_true(all([(os%get(i, 'updatecnt') == merge(2., 0., i == 5), i=1,10)]), 'sample4rec stamps nothing')
        call os%kill
    end subroutine test_sample4rec_global_coverage

    ! The counts behind the realized fraction f = n/N, shared by 2D (class) and 3D (state).
    ! Fixture: class 1 = rows 1,2,5,6; class 2 = rows 3,4,7,8; row 8 inactive, row 6 never
    ! updated; state 1 = rows 1-4, state 2 = rows 5-8; current marker (2) on rows 1,2,3,8.
    ! Counted by hand: class N = [3, 3], n = [2, 1]; state N = [4, 2], n = [3, 0].
    subroutine test_group_update_counts()
        integer, parameter :: NREP_CLS(2) = [3, 3], NSMP_CLS(2) = [2, 1]
        integer, parameter :: NREP_ST(2)  = [4, 2], NSMP_ST(2)  = [3, 0]
        type(oris)           :: os
        integer, allocatable :: nrep(:), nsmp(:)
        real,    allocatable :: rho(:)
        integer              :: i
        write(*,'(A)') 'test_group_update_counts'
        call os%new(8, .true.)
        call os%set_all2single('state',     1.0)
        call os%set_all2single('updatecnt', 1.0)
        call os%set_all2single('sampled',   1.0)
        do i = 1, 8
            call os%set(i, 'class', merge(1., 2., any(i == [1,2,5,6])))
        end do
        call os%set(6, 'updatecnt', 0.)
        call os%set_state(8, 0)
        do i = 1, 3
            call os%set(i, 'sampled', 2.)
        end do
        call os%set(8, 'sampled', 2.)
        call os%get_group_update_counts('class', 2, nrep, nsmp)
        call assert_true(all(nrep == NREP_CLS) .and. all(nsmp == NSMP_CLS), 'group update counts: classes')
        do i = 5, 8
            call os%set_state(i, 2)
        end do
        call os%set_state(8, 0)
        call os%get_group_update_counts('state', 2, nrep, nsmp)
        call assert_true(all(nrep == NREP_ST) .and. all(nsmp == NSMP_ST), 'group update counts: states')
        call os%get_state_update_fracs(2, rho)
        call assert_true(all(abs(rho - [0.75, 0.]) < 1.e-6), 'state update fractions = n/N')
        call os%kill
    end subroutine test_group_update_counts

    ! Closed-form cases of population_blend_weights, each worked out by hand from
    ! s = u/f, w = (1-u)*N/M, mass = s*n + w*M (note Section 4.1)
    subroutine test_blend_weight_cases()
        ! note Section 3.3 example: stored sums represent 100, 10 resampled + 30 first-time
        real, parameter :: W_NOTE    = 0.9     ! (130-40)/100
        real, parameter :: MASS_OLD  = 40. + (90./130.) * 100.  ! current rule, 109.23 (the note rounds to 109)
        ! ufrac override: N = M = 100, n = 10, u = 0.5 -> s = 5, w = 0.5, current share 50/100
        real, parameter :: S_UFRAC   = 5.0
        real, parameter :: W_UFRAC   = 0.5
        real, parameter :: TOL       = 1.e-5
        real :: s, w, mnew, f
        integer :: n
        write(*,'(A)') 'test_blend_weight_cases'
        call population_blend_weights(130, 40, 100., s, w, mnew)
        call assert_real(1.0,    s,    TOL, 'note example: s = 1 by default')
        call assert_real(W_NOTE, w,    TOL, 'note example: w = (N-n)/M')
        call assert_real(130.,   mnew, TOL, 'note example: mass after blend = N')
        f = 40. / 130.
        call assert_real(MASS_OLD, 40. + (1. - f) * 100., TOL, 'note example: the current rule loses mass')
        ! no population change (M = N): the rule is the current recurrence w = 1 - f
        do n = 0, 100, 10
            call population_blend_weights(100, n, 100., s, w, mnew)
            call assert_real(1. - real(n) / 100., w, TOL, 'no population change: w = 1 - f')
            call assert_real(100., mnew, TOL, 'no population change: mass = N')
        end do
        ! f = 1: replace
        call population_blend_weights(50, 50, 80., s, w, mnew)
        call assert_real(1.,  s,    TOL, 'f = 1: s = 1')
        call assert_real(0.,  w,    TOL, 'f = 1: w = 0')
        call assert_real(50., mnew, TOL, 'f = 1: mass = N')
        ! n = 0 (f = 0): keep the previous sums at mass N, w > 1 when rows returned
        call population_blend_weights(60, 0, 40., s, w, mnew)
        call assert_real(0.,  s,    TOL, 'n = 0: no current contribution')
        call assert_real(1.5, w,    TOL, 'n = 0: w = N/M')
        call assert_real(60., mnew, TOL, 'n = 0: mass = N')
        ! M = 0: nothing to carry
        call population_blend_weights(60, 20, 0., s, w, mnew)
        call assert_real(1.,  s,    TOL, 'M = 0: s = 1')
        call assert_real(0.,  w,    TOL, 'M = 0: w = 0')
        call assert_real(20., mnew, TOL, 'M = 0: mass is what is stored, n')
        ! empty group
        call population_blend_weights(0, 0, 25., s, w, mnew)
        call assert_real(0., s + w + mnew, TOL, 'N = 0: s = w = mass = 0')
        ! ufrac override: current-map coefficient stays u, mass stays N
        call population_blend_weights(100, 10, 100., s, w, mnew, ufrac=0.5)
        call assert_real(S_UFRAC, s,    TOL, 'ufrac: s = u/f')
        call assert_real(W_UFRAC, w,    TOL, 'ufrac: w = (1-u)*N/M')
        call assert_real(100.,    mnew, TOL, 'ufrac: mass = N')
        call assert_real(0.5, s * 10. / mnew, TOL, 'ufrac: current-map coefficient = u')
        ! ufrac override with a changed population: 30 first-time rows, N = 130, M = 100
        call population_blend_weights(130, 40, 100., s, w, mnew, ufrac=0.2)
        call assert_real(130., mnew, TOL, 'ufrac with population change: mass = N')
        call assert_real(0.2, s * 40. / mnew, TOL, 'ufrac with population change: coefficient = u')
        ! u = 1 replaces the stored sums
        call population_blend_weights(100, 25, 100., s, w, mnew, ufrac=1.0)
        call assert_real(0.,   w,    TOL, 'u = 1: w = 0')
        call assert_real(100., mnew, TOL, 'u = 1: mass = N')
    end subroutine test_blend_weight_cases

    ! Scripted rounds on a project with unit mass per particle contribution: each
    ! round samples rows (new 'sampled' marker, updatecnt + 1) and changes the
    ! population: group moves, first-time rows, appends, deactivation and
    ! re-activation. The counts come from oris%get_group_update_counts, as in
    ! production. Under the population rule the carried mass of every group equals
    ! N(g) after every blend; the current rule (w = 1 - f) is run as a control and
    ! loses mass on the first-time round (hand-computed value).
    subroutine test_blend_mass_scripted()
        integer, parameter :: NG = 2
        ! control, group 1 by hand: round 1 mass 6; round 2 N = 5 (row 2 left), n = 1,
        ! mass 1 + (1 - 1/5)*6 = 5.8; round 3 N = 7 (first-time rows 13, 14), n = 3,
        ! mass 3 + (1 - 3/7)*5.8 = 6.314 < 7
        real, parameter :: MASS_OLD_R3 = 3. + (1. - 3./7.) * (1. + (1. - 1./5.) * 6.)
        real, parameter :: TOL         = 1.e-5
        type(oris)           :: os
        integer, allocatable :: nrep(:), nsmp(:)
        real    :: mstored(NG), mold(NG), s, w, mnew, f
        integer :: round, g, i
        write(*,'(A)') 'test_blend_mass_scripted'
        call os%new(12, .true.)
        call os%set_all2single('state',     1.0)
        call os%set_all2single('updatecnt', 0.0)
        call os%set_all2single('sampled',   0.0)
        do i = 1, 12
            call os%set(i, 'class', merge(1., 2., i <= 6))
        end do
        mstored = 0.
        mold    = 0.
        do round = 1, 6
            select case(round)
                case(1) ! full first update
                    call sample([(i, i=1,12)])
                case(2) ! resampling, row 2 moves to group 2
                    call os%set(2, 'class', 2.)
                    call sample([1,2,7])
                case(3) ! append rows 13-16 to group 1, first-time rows 13 and 14
                    call os%reallocate(16)
                    do i = 13, 16
                        call os%set_state(i, 1)
                        call os%set(i, 'class',     1.)
                        call os%set(i, 'updatecnt', 0.)
                        call os%set(i, 'sampled',   0.)
                    end do
                    call sample([3,13,14])
                case(4) ! deactivate rows 4 and 5
                    call os%set_state(4, 0)
                    call os%set_state(5, 0)
                    call sample([6])
                case(5) ! re-activate rows 4 and 5
                    call os%set_state(4, 1)
                    call os%set_state(5, 1)
                    call sample([8])
                case(6) ! nothing sampled in group 1
                    call sample([9,10])
            end select
            call os%get_group_update_counts('class', NG, nrep, nsmp)
            do g = 1, NG
                call population_blend_weights(nrep(g), nsmp(g), mstored(g), s, w, mnew)
                ! blend unit-mass sums: current mass is the number of sampled rows
                mstored(g) = s * real(nsmp(g)) + w * mstored(g)
                call assert_real(real(nrep(g)), mstored(g), TOL, 'population rule: carried mass = N(g)')
                call assert_real(mnew, mstored(g), TOL, 'population rule: recorded M = carried mass')
                ! control: the current rule
                f = 0.
                if( nrep(g) > 0 ) f = real(nsmp(g)) / real(nrep(g))
                mold(g) = real(nsmp(g)) + (1. - f) * mold(g)
            end do
            if( round == 3 ) call assert_real(MASS_OLD_R3, mold(1), TOL, 'control: the current rule loses mass')
        end do
        call os%kill

    contains

        ! stamp the rows with a new sampled marker and count their update
        subroutine sample( rows )
            integer, intent(in) :: rows(:)
            integer :: j
            do j = 1, size(rows)
                call os%set(rows(j), 'sampled',   real(round))
                call os%set(rows(j), 'updatecnt', os%get(rows(j), 'updatecnt') + 1.)
            end do
        end subroutine sample

    end subroutine test_blend_mass_scripted

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
        ! order / order_cls: best correlation first, state ignored
        inds = os%order()
        call assert_int(n, size(inds), 'order size')
        call assert_true(all(inds == [4, 2, 5, 1, 3]), 'order sorts by corr, best first')
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
    ! 11. reseed_classes (solve2D cls_init=prev seed partition)
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

    !---------------------------------------------------------------
    ! 15. assignment (from the retired in-module test_oris)
    !---------------------------------------------------------------
    ! intrinsic assignment of an oris goes through the defined assignment of its ori components:
    ! every field is copied and the copy is independent of the original
    subroutine test_assignment()
        type(oris) :: os, os2
        integer    :: i, n
        logical    :: same
        write(*,'(A)') 'test_assignment'
        n = 5
        call os%new(n, .false.)
        call os%rnd_oris(5.)
        os2  = os
        same = os2%get_noris() == n
        do i = 1,n
            if( any(abs(os2%get_euler(i) - os%get_euler(i)) > 0.) ) same = .false.
            if( any(abs(os2%get_2Dshift(i) - os%get_2Dshift(i)) > 0.) ) same = .false.
        end do
        call assert_true(same, 'assignment copies every orientation and shift')
        call os%set_euler(1, [11., 22., 33.])
        call os%set_shift(1, [4., -4.])
        call assert_true(any(abs(os2%get_euler(1) - [11., 22., 33.]) > 1.e-3), 'the copy keeps its angles when the original changes')
        call assert_true(any(abs(os2%get_2Dshift(1) - [4., -4.]) > 1.e-3),     'the copy keeps its shift when the original changes')
        call os%kill
        call os2%kill
    end subroutine test_assignment

    !> The orientation histogram of a state: its particles only, in the elevation band of their
    !! projection direction, on the grid the histogram's shape gives
    subroutine test_oridist_from_oris()
        type(oris) :: os
        integer    :: hist(72,36), coarse(4,2)
        write(*,'(A)') 'test_oridist_from_oris'
        call os%new(5, is_ptcl=.true.)
        ! two north-pole and one south-pole directions in state 1, two in state 2
        call os%set_euler(1, [0.,   0., 0.])
        call os%set_euler(2, [30.,  0., 0.])
        call os%set_euler(3, [0., 180., 0.])
        call os%set_euler(4, [0.,   0., 0.])
        call os%set_euler(5, [0.,   0., 0.])
        call os%set_state(1, 1)
        call os%set_state(2, 1)
        call os%set_state(3, 1)
        call os%set_state(4, 2)
        call os%set_state(5, 2)
        call oridist_from_oris(os, 1, hist)
        call assert_int(3, sum(hist),        'the particles of state 1 only')
        call assert_int(2, sum(hist(:,36)),  'north-pole directions in the top elevation band')
        call assert_int(1, sum(hist(:,1)),   'south-pole directions in the bottom band')
        call oridist_from_oris(os, 2, hist)
        call assert_int(2, sum(hist),        'the particles of state 2 only')
        call oridist_from_oris(os, 1, coarse)
        call assert_int(2, sum(coarse(:,2)), 'the bins follow the histogram''s shape')
        call assert_int(1, sum(coarse(:,1)), 'in both bands')
        call os%kill
    end subroutine test_oridist_from_oris

end module simple_oris_tester
