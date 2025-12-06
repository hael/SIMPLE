module simple_oris_tester
include 'simple_lib.f08'
use simple_ori
use simple_oris          ! provides type(oris), type(ori) and methods
use simple_sym
use simple_defs          ! for dp, etc.
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
        call test_euler_angles
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
            call os%set(i, 'w',         1.0)
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
        integer, allocatable :: inds(:)
        integer :: n, nsamp
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
        ! sample4update_cnt (prefers low updatecnt)
        call os%set_all2single('updatecnt', 0.0)
        call os%sample4update_cnt([1, n], frac, nsamp, inds, .true.)
        call assert_true(nsamp <= n, 'sample4update_cnt nsamp<=n')
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
    end subroutine test_sampling_and_updatecnt

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
    ! 11. Euler angles
    !---------------------------------------------------------------
    subroutine test_euler_angles
        type(oris) :: oris_obj
        type(sym)  :: pgrpsyms !< symmetry elements object
        type(ori)  :: o1, o2
        real       :: R(3,3),eullims(3,2), euls(3), euls2(3), diff, error, threshold
        integer    :: isample, n
        integer, parameter :: NSAMPLE = 100
        call seed_rnd
        call pgrpsyms%new('c4')
        oris_obj = oris(1, is_ptcl=.false.)
        call oris_obj%set_euler(1, [45., 90., 180.])
        call oris_obj%get_ori(1, o1)
        eullims = pgrpsyms%get_eullims()
        print *, 'lim1: ', eullims(1,1), eullims(1,2)
        print *, 'lim2: ', eullims(2,1), eullims(2,2)
        print *, 'lim3: ', eullims(3,1), eullims(3,2)
        print *, 'RND EULS'
        call o2%new(.true.)
        do isample = 1,NSAMPLE
            call pgrpsyms%rnd_euler(o2)
            euls = o2%get_euler()
            print *, euls(1), euls(2), euls(3)
        end do
        ! consistency between spider and new m2euler routines
        n = 0
        error = 0.
        threshold = 0.01 ! degrees
        print *,'--- random'
        do isample = 1,NSAMPLE
            call pgrpsyms%rnd_euler(o2)
            euls  = o2%get_euler()
            R     = o2%get_mat()
            euls  = m2euler(R); euls2 = m2euler_fast(R)
            diff  = angular_error(euls, euls2)
            error = error + diff
            if( diff > threshold )then
                n = n+1
                print *,n,isample,diff,euls,euls2
            endif
        enddo
        print *,'--- phi=0'
        do isample = 1,NSAMPLE
            call pgrpsyms%rnd_euler(o2)
            euls  = o2%get_euler()
            euls(1) = 0.
            call o2%set_euler(euls)
            R     = o2%get_mat()
            euls  = m2euler(R); euls2 = m2euler_fast(R)
            diff  = angular_error(euls, euls2)
            error = error + diff
            if( diff > threshold )then
                n = n+1
                print *,n,isample,diff,euls,euls2
            endif
        enddo
        print *,'--- phi=90'
        do isample = 1,NSAMPLE
            call pgrpsyms%rnd_euler(o2)
            euls  = o2%get_euler()
            euls(1) = 90.
            call o2%set_euler(euls)
            R     = o2%get_mat()
            euls  = m2euler(R); euls2 = m2euler_fast(R)
            diff  = angular_error(euls, euls2)
            error = error + diff
            if( diff > threshold )then
                n = n+1
                print *,n,isample,diff,euls,euls2
            endif
        enddo
        print *,'--- psi=0'
        do isample = 1,NSAMPLE
            call pgrpsyms%rnd_euler(o2)
            euls  = o2%get_euler()
            euls(3) = 0.
            call o2%set_euler(euls)
            R     = o2%get_mat()
            euls  = m2euler(R); euls2 = m2euler_fast(R)
            diff  = angular_error(euls, euls2)
            error = error + diff
            if( diff > threshold )then
                n = n+1
                print *,n,isample,diff,euls,euls2
            endif
        enddo
        print *,'--- psi=270'
        do isample = 1,NSAMPLE
            call pgrpsyms%rnd_euler(o2)
            euls  = o2%get_euler()
            euls(3) = 270.
            call o2%set_euler(euls)
            R     = o2%get_mat()
            euls  = m2euler(R); euls2 = m2euler_fast(R)
            diff  = angular_error(euls, euls2)
            error = error + diff
            if( diff > threshold )then
                n = n+1
                print *,n,isample,diff,euls,euls2
            endif
        enddo
        print *,'--- theta=0'
        do isample = 1,NSAMPLE
            call pgrpsyms%rnd_euler(o2)
            euls  = o2%get_euler()
            euls(2) = 0.
            call o2%set_euler(euls)
            R     = o2%get_mat()
            euls  = m2euler(R); euls2 = m2euler_fast(R)
            diff  = min(angular_error(euls, euls2), angular_error(euls, [euls2(3),euls2(2),euls2(1)]))
            error = error + diff
            if( diff > threshold )then
                n = n+1
                print *,n,isample,diff,euls,euls2
            endif
        enddo
        if( n > 0 )then
            write(*,'(A,I6,A)')'M2EULER TEST NOT PASSED: ',n, ' ERRORS'
        else
            write(*,*)'M2EULER TEST PASSED'
        endif
        write(*,'(A,F9.6,A)')'AVERAGE DIFFERENCE: ',error/real(6*NSAMPLE),' degrees'

        contains

            real function angular_error( angs1, angs2 )
                real, intent(in) :: angs1(3), angs2(3)
                angular_error = min((angs1(1)-angs2(1))**2, (360.-angs1(1)-angs2(1))**2)
                angular_error = angular_error + (angs1(2)-angs2(2))**2
                angular_error = angular_error + min((angs1(3)-angs2(3))**2, (360.-angs1(3)-angs2(3))**2)
                angular_error = sqrt(angular_error/3.)
            end function angular_error

        end subroutine test_euler_angles

end module simple_oris_tester
