!@descr: unit tests for the PCG order prior's closed-form choice: evidence weight fraction, membrane and soluble fixtures, core criterion, threads
submodule (simple_pcg_order_prior) simple_pcg_order_prior_tester
!$ use omp_lib,      only: omp_get_max_threads, omp_set_num_threads
use simple_test_utils, only: assert_true, set_fixed_seed
implicit none

integer, parameter :: TBOX    = 48
integer, parameter :: TSEED   = 20261108
real,    parameter :: TSMPD   = 2.0
real,    parameter :: TRES    = 6.0     !< the fixture pair's FSC=0.143 (A)
real,    parameter :: SUP_R   = 20.     !< radii in voxels: evaluation domain
real,    parameter :: BLOB_R  = 4.5     !< protein blob centres
real,    parameter :: CORE_R  = 9.5     !< protein: NU evidence and fine labels
real,    parameter :: MIC_R   = 12.     !< micelle shell radius and width
real,    parameter :: MIC_W   = 1.5
real,    parameter :: DENSE_R = 15.5    !< the membrane fixture's dense region: protein and micelle

contains

    module subroutine run_all_pcg_order_prior_tests()
        write(logfhandle,'(A)') '**** running all PCG order prior tests ****'
        call test_evidence_fraction()
        call test_membrane_choice()
        call test_soluble_choice()
        call test_core_criterion()
        call test_band_diagonal()
        call test_team_matches_serial()
    end subroutine run_all_pcg_order_prior_tests

    !> a margin without ties: w_ord >= 1/2 exactly where the margin reaches the (1-f) quantile
    subroutine test_evidence_fraction()
        type(image)            :: evmap, w_order(2)
        real(c_float), pointer :: rp(:,:,:)
        real,    allocatable   :: rnd(:,:,:)
        logical, allocatable   :: lsup(:,:,:)
        logical :: l_valid, l_ok
        real    :: frac
        integer :: ifr, nsup
        write(logfhandle,'(A)') '-- test_evidence_fraction'
        call set_fixed_seed(TSEED)
        lsup = ball(SUP_R)
        nsup = count(lsup)
        allocate(rnd(TBOX,TBOX,TBOX))
        call random_number(rnd)
        call evmap%new([TBOX,TBOX,TBOX], TSMPD)
        call evmap%get_rmat_ptr(rp)
        rp(:TBOX,:TBOX,:TBOX) = 0.
        where( lsup ) rp(:TBOX,:TBOX,:TBOX) = 1. + rnd
        nullify(rp)
        l_ok = .true.
        do ifr = 1, size(ORDER_FRAC_GRID)
            call build_order_weights(evmap, lsup, ORDER_FRAC_GRID(ifr), w_order, l_valid)
            call w_order(1)%get_rmat_ptr(rp)
            frac = real(count(lsup .and. rp(:TBOX,:TBOX,:TBOX) >= 0.5)) / real(nsup)
            nullify(rp)
            ! a few voxels: the quantile index rounds and single-precision draws tie
            l_ok = l_ok .and. l_valid .and. abs(frac - ORDER_FRAC_GRID(ifr)) <= 4. / real(nsup)
        end do
        call assert_true(l_ok, 'the logistic at the (1-f) quantile orders a fraction f of the domain')
        call evmap%kill
        call w_order(1)%kill
        call w_order(2)%kill
    end subroutine test_evidence_fraction

    !> protein with evidence and fine labels, a dense micelle without either, solvent: with a flat diagonal the
    !! strongest ridge shrinks the micelle twentyfold per half while the core, fully ordered, keeps its signal
    subroutine test_membrane_choice()
        type(image)             :: base(2), w_dens(2), evmap, locres, w_order(2)
        type(order_prior_stats) :: st
        logical, allocatable    :: lsup(:,:,:)
        write(logfhandle,'(A)') '-- test_membrane_choice'
        call fixture(.true., base, w_dens, evmap, locres, lsup)
        call select_order_prior(base, w_dens, evmap, locres, TRES, lsup, flat_diag(1.), [1.,1.], -1., w_order, st, .false.)
        write(logfhandle,'(A,F6.3,A,F6.2,A,F7.4,A,F7.4)') 'fraction ', st%frac, ', lambda_rel ', st%lambda_rel, &
            &', removed ', st%removed, ', core kept ', st%core_kept
        call assert_true(.not. st%l_off, 'a dense, unordered micelle gets the order ridge')
        call assert_true(st%core_kept >= CORE_SIGNAL_MIN, 'the chosen setting keeps 95 % of the core signal')
        call assert_true(st%removed >= 0.9, 'the chosen setting removes most of the micelle signal')
        ! a strength set as given: the ordered fraction still spares the core and takes the micelle
        call select_order_prior(base, w_dens, evmap, locres, TRES, lsup, flat_diag(1.), [1.,1.], 20., w_order, st, .false.)
        call assert_true(.not. st%l_off .and. st%lambda_rel == 20., 'a set strength is applied')
        call assert_true(st%removed >= 0.9 .and. st%core_kept >= CORE_SIGNAL_MIN, &
            &'at a set strength the order weight takes the micelle and spares the core')
        call kill_set(base, w_dens, evmap, locres, w_order)
    end subroutine test_membrane_choice

    !> everything dense is ordered: nothing around the core to remove, so the strongest strength that keeps the
    !! core (solvent flattening); without evidence of order no setting keeps it, so no prior
    subroutine test_soluble_choice()
        type(image)             :: base(2), w_dens(2), evmap, locres, w_order(2)
        type(order_prior_stats) :: st
        logical, allocatable    :: lsup(:,:,:)
        write(logfhandle,'(A)') '-- test_soluble_choice'
        call fixture(.false., base, w_dens, evmap, locres, lsup)
        call select_order_prior(base, w_dens, evmap, locres, TRES, lsup, flat_diag(1.), [1.,1.], -1., w_order, st, .false.)
        call assert_true(.not. st%l_off .and. st%core_kept >= CORE_SIGNAL_MIN, 'a soluble fixture gets a prior that keeps its core')
        call assert_true(st%lambda_rel == ORDER_LAMBDA_GRID(size(ORDER_LAMBDA_GRID)), &
            &'with nothing around the core, the strongest strength that keeps it')
        call random_evidence(evmap, lsup)
        call select_order_prior(base, w_dens, evmap, locres, TRES, lsup, flat_diag(1.e-2), [1.,1.], -1., w_order, st, .false.)
        call assert_true(st%l_off .and. st%lambda_rel == 0., 'no prior when no setting keeps the core')
        call kill_set(base, w_dens, evmap, locres, w_order)
    end subroutine test_soluble_choice

    !> with a flat diagonal every band shrinks alike: the selector's core figure is sum_core x_e x_o s_e s_o /
    !! sum_core x_e x_o, s = 1/(1 + lambda (1 - w_order)), recomputed on the shared helpers' regions and weights
    subroutine test_core_criterion()
        type(image)             :: base(2), w_dens(2), evmap, locres, w_order(2)
        type(order_prior_stats) :: st
        real(c_float), pointer  :: rx1(:,:,:), rx2(:,:,:), rw1(:,:,:), rw2(:,:,:)
        logical, allocatable    :: lsup(:,:,:), lcore(:,:,:), lperi(:,:,:)
        real(dp) :: num, den, p, s1, s2
        integer  :: i, j, k
        write(logfhandle,'(A)') '-- test_core_criterion'
        call fixture(.true., base, w_dens, evmap, locres, lsup)
        call select_order_prior(base, w_dens, evmap, locres, TRES, lsup, flat_diag(1.), [1.,1.], -1., w_order, st, .false.)
        call order_regions(w_dens, locres, TRES, lsup, lcore, lperi)
        call base(1)%get_rmat_ptr(rx1)
        call base(2)%get_rmat_ptr(rx2)
        call w_order(1)%get_rmat_ptr(rw1)
        call w_order(2)%get_rmat_ptr(rw2)
        num = 0._dp
        den = 0._dp
        do k = 1, TBOX
            do j = 1, TBOX
                do i = 1, TBOX
                    if( .not. lcore(i,j,k) ) cycle
                    p   = real(rx1(i,j,k),dp) * real(rx2(i,j,k),dp)
                    s1  = 1._dp / (1._dp + st%lambda_rel * (1._dp - rw1(i,j,k)))
                    s2  = 1._dp / (1._dp + st%lambda_rel * (1._dp - rw2(i,j,k)))
                    num = num + p * s1 * s2
                    den = den + p
                end do
            end do
        end do
        nullify(rx1, rx2, rw1, rw2)
        ! the bands sum to the half up to transform round-off
        call assert_true(abs(real(num / den) - st%core_kept) <= 1.e-4, &
            &'the selector''s core figure is the cross-half signal kept under the chosen weights')
        call kill_set(base, w_dens, evmap, locres, w_order)
    end subroutine test_core_criterion

    !> no evidence of order, so the core is penalized at every fraction: at one strength, data weak below shell 4 and
    !! ten times stronger above keep more of it than data weak everywhere and less than data strong everywhere
    subroutine test_band_diagonal()
        type(image)             :: base(2), w_dens(2), evmap, locres, w_order(2)
        type(order_prior_stats) :: st_lo, st_step, st_hi
        logical, allocatable    :: lsup(:,:,:)
        real,    allocatable    :: dsh(:,:)
        write(logfhandle,'(A)') '-- test_band_diagonal'
        call fixture(.false., base, w_dens, evmap, locres, lsup)
        call random_evidence(evmap, lsup)
        call select_order_prior(base, w_dens, evmap, locres, TRES, lsup, flat_diag(0.01), [1.,1.], 1., w_order, st_lo, .false.)
        call select_order_prior(base, w_dens, evmap, locres, TRES, lsup, flat_diag(0.1),  [1.,1.], 1., w_order, st_hi, .false.)
        allocate(dsh(0:TBOX,2), source=0.1)
        dsh(0:3,:) = 0.01
        call select_order_prior(base, w_dens, evmap, locres, TRES, lsup, dsh, [1.,1.], 1., w_order, st_step, .false.)
        write(logfhandle,'(A,3F8.4)') 'core kept, weak / stepped / strong data ', st_lo%core_kept, st_step%core_kept, &
            &st_hi%core_kept
        call assert_true(st_lo%core_kept < st_step%core_kept .and. st_step%core_kept < st_hi%core_kept, &
            &'a band is shrunk by the diagonal of its own shells')
        call kill_set(base, w_dens, evmap, locres, w_order)
    end subroutine test_band_diagonal

    !> the selector's OpenMP reduction gives one thread's choice; its figures differ by summation order only
    subroutine test_team_matches_serial()
        type(image)             :: base(2), w_dens(2), evmap, locres, w_order(2)
        type(order_prior_stats) :: st1, st3
        logical, allocatable    :: lsup(:,:,:)
        integer :: nthr_prev
        write(logfhandle,'(A)') '-- test_team_matches_serial'
        call fixture(.true., base, w_dens, evmap, locres, lsup)
        nthr_prev = 1
        !$ nthr_prev = omp_get_max_threads()
        !$ call omp_set_num_threads(1)
        call select_order_prior(base, w_dens, evmap, locres, TRES, lsup, flat_diag(1.), [1.,1.], -1., w_order, st1, .false.)
        !$ call omp_set_num_threads(3)
        call select_order_prior(base, w_dens, evmap, locres, TRES, lsup, flat_diag(1.), [1.,1.], -1., w_order, st3, .false.)
        !$ call omp_set_num_threads(nthr_prev)
        call assert_true(st1%frac == st3%frac .and. st1%lambda_rel == st3%lambda_rel, &
            &'a team of three chooses the setting one thread chooses')
        call assert_true(abs(st1%removed - st3%removed) <= 1.e-6 .and. abs(st1%core_kept - st3%core_kept) <= 1.e-6, &
            &'a team of three gives one thread''s figures to summation order')
        call kill_set(base, w_dens, evmap, locres, w_order)
    end subroutine test_team_matches_serial

    ! ---------------- fixture ----------------

    !> protein blobs with evidence ~2 and labels at half the FSC inside CORE_R; with l_micelle a smooth shell at
    !! MIC_R, dense but without evidence or fine labels; independent noise per half
    subroutine fixture( l_micelle, base, w_dens, evmap, locres, lsup )
        logical,              intent(in)    :: l_micelle
        type(image),          intent(inout) :: base(2), w_dens(2), evmap, locres
        logical, allocatable, intent(inout) :: lsup(:,:,:)
        real,    allocatable   :: sig(:,:,:), rnd(:,:,:), centres(:,:), dist(:,:,:)
        real(c_float), pointer :: rp(:,:,:)
        real    :: c, r2, r_dense
        integer :: i, j, k, b, ih
        call set_fixed_seed(TSEED)
        lsup = ball(SUP_R)
        dist = radius()
        allocate(centres(3,10))
        call random_number(centres)
        centres = (centres - 0.5) * 2. * BLOB_R / sqrt(3.)
        allocate(sig(TBOX,TBOX,TBOX), source=0.)
        c = real(TBOX/2 + 1)
        do k = 1, TBOX
            do j = 1, TBOX
                do i = 1, TBOX
                    do b = 1, size(centres,2)
                        r2 = (real(i) - c - centres(1,b))**2 + (real(j) - c - centres(2,b))**2 + (real(k) - c - centres(3,b))**2
                        sig(i,j,k) = sig(i,j,k) + exp(-0.5 * r2)
                    end do
                    if( l_micelle ) sig(i,j,k) = sig(i,j,k) + 0.5 * exp(-0.5 * ((dist(i,j,k) - MIC_R) / MIC_W)**2)
                end do
            end do
        end do
        allocate(rnd(TBOX,TBOX,TBOX))
        do ih = 1, 2
            call random_number(rnd)
            call base(ih)%new([TBOX,TBOX,TBOX], TSMPD)
            call base(ih)%set_rmat(sig + 0.1 * sqrt(12.) * (rnd - 0.5), .false.)
        end do
        ! soluble: dense exactly where ordered, so nothing dense lies outside the core
        r_dense = CORE_R - 0.5
        if( l_micelle ) r_dense = DENSE_R
        call w_dens(1)%new([TBOX,TBOX,TBOX], TSMPD)
        call w_dens(1)%set_rmat(1. / (1. + exp((dist - r_dense) / 0.7)), .false.)
        call w_dens(2)%copy(w_dens(1))
        call random_number(rnd)
        call evmap%new([TBOX,TBOX,TBOX], TSMPD)
        call evmap%get_rmat_ptr(rp)
        rp(:TBOX,:TBOX,:TBOX) = 0.
        where( lsup ) rp(:TBOX,:TBOX,:TBOX) = 0.01 * rnd
        where( lsup .and. dist < CORE_R ) rp(:TBOX,:TBOX,:TBOX) = 2. + 0.1 * rnd
        nullify(rp)
        call locres%new([TBOX,TBOX,TBOX], TSMPD)
        call locres%get_rmat_ptr(rp)
        rp(:TBOX,:TBOX,:TBOX) = 0.
        where( lsup ) rp(:TBOX,:TBOX,:TBOX) = 4. * TRES
        where( lsup .and. dist < CORE_R ) rp(:TBOX,:TBOX,:TBOX) = 0.5 * TRES
        nullify(rp)
        deallocate(sig, rnd, centres, dist)
    end subroutine fixture

    function radius() result( d )
        real, allocatable :: d(:,:,:)
        real    :: c
        integer :: i, j, k
        allocate(d(TBOX,TBOX,TBOX))
        c = real(TBOX/2 + 1)
        do k = 1, TBOX
            do j = 1, TBOX
                do i = 1, TBOX
                    d(i,j,k) = sqrt((real(i) - c)**2 + (real(j) - c)**2 + (real(k) - c)**2)
                end do
            end do
        end do
    end function radius

    !> evidence without structure: uniform noise over the domain
    subroutine random_evidence( evmap, lsup )
        type(image), intent(inout) :: evmap
        logical,     intent(in)    :: lsup(:,:,:)
        real(c_float), pointer :: rp(:,:,:)
        real, allocatable      :: rnd(:,:,:)
        allocate(rnd(TBOX,TBOX,TBOX))
        call random_number(rnd)
        call evmap%get_rmat_ptr(rp)
        rp(:TBOX,:TBOX,:TBOX) = 0.
        where( lsup ) rp(:TBOX,:TBOX,:TBOX) = 0.01 * rnd
        nullify(rp)
    end subroutine random_evidence

    !> a shell diagonal of one value for both halves
    function flat_diag( d ) result( dsh )
        real, intent(in)  :: d
        real, allocatable :: dsh(:,:)
        allocate(dsh(0:TBOX,2), source=d)
    end function flat_diag

    function ball( r ) result( l )
        real, intent(in)     :: r
        logical, allocatable :: l(:,:,:)
        real,    allocatable :: d(:,:,:)
        d = radius()
        l = d <= r
    end function ball

    subroutine kill_set( base, w_dens, evmap, locres, w_order )
        type(image), intent(inout) :: base(2), w_dens(2), evmap, locres, w_order(2)
        integer :: ih
        do ih = 1, 2
            call base(ih)%kill
            call w_dens(ih)%kill
            call w_order(ih)%kill
        end do
        call evmap%kill
        call locres%kill
    end subroutine kill_set

end submodule simple_pcg_order_prior_tester
