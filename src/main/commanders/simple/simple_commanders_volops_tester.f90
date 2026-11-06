!@descr: library tests of the volume programs symaxis_search, dock_volpair and centervol on phantoms of known symmetry, rotation and offset
! The phantoms are sums of Gaussian blobs placed analytically, so a rotated or displaced copy is exact, not interpolated.
module simple_commanders_volops_tester
use simple_defs,                only: logfhandle
use simple_string,              only: string
use simple_cmdline,             only: cmdline
use simple_image,               only: image
use simple_ui,                  only: make_ui
use simple_commanders_volops,   only: commander_symaxis_search, commander_dock_volpair, commander_centervol
use simple_test_truth_metrics,  only: axis_operator, group_operators
use simple_test_state_fixture,  only: point_group_check
use simple_test_utils,          only: assert_true, enter_fixture, leave_fixture, tests_failed
implicit none
private
public :: run_all_volops_lib_tests

integer, parameter :: BOX = 64, NSEED = 4, NTHR = 4
real,    parameter :: SMPD = 2.0, MSKDIAM = 100., LP = 8., SIGMA = 3.0                !< A
real,    parameter :: SEED(3,NSEED) = reshape([12., 5., 3.,  -4., 14., 6.,  6., -9., 11.,  -10., -6., -8.], [3,NSEED])
real,    parameter :: AMPS(NSEED)   = [1.0, 0.8, 0.6, 0.9]
real,    parameter :: SEED_I(3,2) = reshape([24., 16., 13.,  -10., 22., -15.], [3,2])  !< 32 and 28 A out: 120 blobs apart
real,    parameter :: CORR_FLOOR = 0.95  !< a map moved by the program back onto its exact target
real,    parameter :: CONTROL_CEIL = 0.8 !< the same comparison before the program ran

contains

    subroutine run_all_volops_lib_tests()
        type(string) :: cwd_saved, root
        integer      :: nfail0
        write(*,'(A)') '**** running all volume program library tests ****'
        nfail0 = tests_failed
        call enter_fixture('volops', cwd_saved, root)
        call make_ui
        call test_symaxis_search_d2()
        call test_symaxis_search_i()
        call test_dock_volpair()
        call test_centervol()
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine run_all_volops_lib_tests

    !> symaxis_search pgrp=d2 on a D2 phantom turned by a generic rotation: the map it aligns to the axis it
    !! finds is invariant under the D2 operators about the box axes, which the turned input is not (the
    !! control); rotvol's interpolation bounds the invariance below 1
    subroutine test_symaxis_search_d2()
        type(commander_symaxis_search)  :: xsymaxis
        class(cmdline), allocatable     :: cl
        real, allocatable :: ops(:,:,:), pos(:,:)
        real    :: rot(3,3), cc_aligned, cc_input, cc_larger
        integer :: g, k, n
        write(*,'(A)') 'test_symaxis_search_d2'
        ops = group_operators('d2')
        rot = axis_operator([1., 2., 3.] / norm2([1., 2., 3.]), 37.)
        allocate(pos(3, NSEED * size(ops, 3)))
        n = 0
        do g = 1, size(ops, 3)
            do k = 1, NSEED
                n = n + 1
                pos(:,n) = matmul(rot, matmul(ops(:,:,g), SEED(:,k)))
            enddo
        enddo
        call write_blobs('d2_turned.mrc', pos, [(AMPS, g = 1, size(ops, 3))])
        allocate(cl)
        call cl%set('prg',     'symaxis_search')
        call cl%set('vol1',    'd2_turned.mrc')
        call cl%set('smpd',    SMPD)
        call cl%set('mskdiam', MSKDIAM)
        call cl%set('lp',      LP)
        call cl%set('pgrp',    'd2')
        call cl%set('center',  'no')
        call cl%set('mkdir',   'no')
        call cl%set('nthr',    NTHR)
        call xsymaxis%execute(cl)
        call cl%kill
        deallocate(cl)
        call point_group_check(string('vol_aligned2_d2axis.mrc'), ops, axis_operator([0., 0., 1.], 90.), cc_aligned, cc_larger)
        call point_group_check(string('d2_turned.mrc'), ops, axis_operator([0., 0., 1.], 90.), cc_input, cc_larger)
        write(logfhandle,'(a,f7.4,a,f7.4)') '>>> symaxis_search D2 self-correlation aligned ', cc_aligned, ', input ', cc_input
        call assert_true(cc_aligned >= CORR_FLOOR, 'symaxis_search: the aligned map is invariant under D2')
        call assert_true(cc_input < CONTROL_CEIL, 'symaxis_search: the turned input is not (the control)')
    end subroutine test_symaxis_search_d2

    !> symaxis_search with pgrp=i aligns a turned icosahedral phantom to the setting the builder's operators
    !! (refinement, reconstruction) use: the aligned map is invariant under them (above the control ceiling; the
    !! search's icosahedral accuracy keeps it under the exact floor), the turned input is not
    subroutine test_symaxis_search_i()
        type(commander_symaxis_search)  :: xsymaxis
        class(cmdline), allocatable     :: cl
        real, allocatable :: ops(:,:,:), pos(:,:)
        real    :: rot(3,3), cc_aligned, cc_input, cc_larger
        integer :: g, k, n
        write(*,'(A)') 'test_symaxis_search_i'
        ops = group_operators('i')
        rot = axis_operator([1., 2., 3.] / norm2([1., 2., 3.]), 37.)
        allocate(pos(3, size(SEED_I,2) * size(ops, 3)))
        n = 0
        do g = 1, size(ops, 3)
            do k = 1, size(SEED_I,2)
                n = n + 1
                pos(:,n) = matmul(rot, matmul(ops(:,:,g), SEED_I(:,k)))
            enddo
        enddo
        call write_blobs('i_turned.mrc', pos, [(AMPS(1:2), g = 1, size(ops, 3))])
        allocate(cl)
        call cl%set('prg',     'symaxis_search')
        call cl%set('vol1',    'i_turned.mrc')
        call cl%set('smpd',    SMPD)
        call cl%set('mskdiam', MSKDIAM)
        call cl%set('lp',      LP)
        call cl%set('pgrp',    'i')
        call cl%set('center',  'no')
        call cl%set('mkdir',   'no')
        call cl%set('nthr',    NTHR)
        call xsymaxis%execute(cl)
        call cl%kill
        deallocate(cl)
        call point_group_check(string('vol_aligned2_iaxis.mrc'), ops, axis_operator([0., 0., 1.], 90.), cc_aligned, cc_larger)
        call point_group_check(string('i_turned.mrc'), ops, axis_operator([0., 0., 1.], 90.), cc_input, cc_larger)
        write(logfhandle,'(a,f7.4,a,f7.4)') '>>> symaxis_search I self-correlation aligned ', cc_aligned, ', input ', cc_input
        call assert_true(cc_aligned >= CONTROL_CEIL, 'symaxis_search: the aligned map is in the builder''s I setting')
        call assert_true(cc_input < CONTROL_CEIL, 'symaxis_search: the turned input is not (the control)')
    end subroutine test_symaxis_search_i

    !> dock_volpair docks a copy of the asymmetric phantom turned by a generic rotation and shifted by
    !! 2-3 pixels back onto the phantom: the docked map correlates with it, the undocked copy does not
    subroutine test_dock_volpair()
        type(commander_dock_volpair) :: xdock
        class(cmdline), allocatable  :: cl
        real    :: rot(3,3), shift(3), pos(3,NSEED), cc_docked, cc_input
        integer :: k
        write(*,'(A)') 'test_dock_volpair'
        rot   = axis_operator([-2., 1., 2.] / 3., 52.)
        shift = [5., -4., 3.]
        do k = 1, NSEED
            pos(:,k) = matmul(rot, SEED(:,k)) + shift
        enddo
        call write_blobs('asym.mrc',       SEED, AMPS)
        call write_blobs('asym_moved.mrc', pos,  AMPS)
        allocate(cl)
        call cl%set('prg',     'dock_volpair')
        call cl%set('vol1',    'asym.mrc')
        call cl%set('vol2',    'asym_moved.mrc')
        call cl%set('outvol',  'asym_docked.mrc')
        call cl%set('outfile', 'dock_report.txt')
        call cl%set('smpd',    SMPD)
        call cl%set('mskdiam', MSKDIAM)
        call cl%set('lp',      LP)
        call cl%set('trs',     6.)
        call cl%set('mkdir',   'no')
        call cl%set('nthr',    NTHR)
        call xdock%execute(cl)
        call cl%kill
        deallocate(cl)
        cc_docked = map_corr('asym.mrc', 'asym_docked.mrc')
        cc_input  = map_corr('asym.mrc', 'asym_moved.mrc')
        write(logfhandle,'(a,f7.4,a,f7.4)') '>>> dock_volpair correlation docked ', cc_docked, ', input ', cc_input
        call assert_true(cc_docked >= CORR_FLOOR, 'dock_volpair: the docked map lies on the reference')
        call assert_true(cc_input < CONTROL_CEIL, 'dock_volpair: the moved input does not (the control)')
    end subroutine test_dock_volpair

    !> centervol (centre of mass) on the asymmetric phantom, its mass centred on the box centre, displaced by
    !! about 4.5 pixels: the written map's centre of mass is back within one voxel of the box centre
    subroutine test_centervol()
        type(commander_centervol)   :: xcenter
        class(cmdline), allocatable :: cl
        real, parameter :: OFFSET(3) = [6., -4., 5.]  !< A
        real    :: com(3), pos(3,NSEED), off_in, off_out
        integer :: k
        write(*,'(A)') 'test_centervol'
        ! equal blobs: the centre of mass is the amplitude-weighted mean of their centres
        com = matmul(SEED, AMPS) / sum(AMPS)
        do k = 1, NSEED
            pos(:,k) = SEED(:,k) - com + OFFSET
        enddo
        call write_blobs('asym_offset.mrc', pos, AMPS)
        allocate(cl)
        call cl%set('prg',     'center')
        call cl%set('vol1',    'asym_offset.mrc')
        call cl%set('smpd',    SMPD)
        call cl%set('mskdiam', MSKDIAM)
        call cl%set('mkdir',   'no')
        call cl%set('nthr',    NTHR)
        call xcenter%execute(cl)
        call cl%kill
        deallocate(cl)
        off_in  = norm2(mass_centre('asym_offset.mrc'))
        off_out = norm2(mass_centre('shifted_vol_state1.mrc'))
        write(logfhandle,'(a,f7.3,a,f7.3)') '>>> centervol offset of the centre of mass (A) input ', off_in, ', output ', off_out
        call assert_true(off_out <= SMPD, 'centervol: the centre of mass is back within one voxel')
        call assert_true(off_in > 3. * SMPD, 'centervol: the input was off centre (the control)')
    end subroutine test_centervol

    !> Gaussian blobs of width SIGMA at pos (A from the box centre) with amplitudes amps, written to fname
    subroutine write_blobs( fname, pos, amps )
        character(len=*), intent(in) :: fname
        real,             intent(in) :: pos(:,:), amps(:)
        type(image) :: vol
        real, allocatable :: rmat(:,:,:)
        real    :: c0(3), p(3)
        integer :: i, j, k, b
        allocate(rmat(BOX,BOX,BOX), source=0.)
        c0 = real(BOX / 2 + 1)
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    p = ([real(i), real(j), real(k)] - c0) * SMPD
                    do b = 1, size(amps)
                        rmat(i,j,k) = rmat(i,j,k) + amps(b) * exp(-sum((p - pos(:,b))**2) / (2. * SIGMA**2))
                    enddo
                enddo
            enddo
        enddo
        call vol%new([BOX,BOX,BOX], SMPD)
        call vol%set_rmat(rmat, .false.)
        call vol%write(string(fname), del_if_exists=.true.)
        call vol%kill
    end subroutine write_blobs

    real function map_corr( fname1, fname2 )
        character(len=*), intent(in) :: fname1, fname2
        type(image) :: v1, v2
        call v1%new([BOX,BOX,BOX], SMPD)
        call v2%new([BOX,BOX,BOX], SMPD)
        call v1%read(string(fname1))
        call v2%read(string(fname2))
        map_corr = v1%real_corr(v2)
        call v1%kill
        call v2%kill
    end function map_corr

    !> the centre of mass of a map's positive density, in A from the box centre
    function mass_centre( fname ) result( com )
        character(len=*), intent(in) :: fname
        real :: com(3)
        type(image) :: vol
        real, allocatable :: rmat(:,:,:)
        real    :: c0(3), m
        integer :: i, j, k
        call vol%new([BOX,BOX,BOX], SMPD)
        call vol%read(string(fname))
        rmat = vol%get_rmat()
        call vol%kill
        c0  = real(BOX / 2 + 1)
        com = 0.
        m   = 0.
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    if( rmat(i,j,k) <= 0. ) cycle
                    com = com + rmat(i,j,k) * ([real(i), real(j), real(k)] - c0)
                    m   = m + rmat(i,j,k)
                enddo
            enddo
        enddo
        com = com / m * SMPD
    end function mass_centre

end module simple_commanders_volops_tester
