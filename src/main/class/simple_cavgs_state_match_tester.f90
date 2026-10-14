!@descr: unit tests of class averages matched against state volumes (simple_cavgs_state_match)
! The reduction of the correlations to a state, reference and margin per class, and the composition of
! a matched direction and in-plane angle into a 3D pose, on fixed numbers; then two small asymmetric
! volumes of Gaussian blobs, projected at known poses on the matching grid and once mirrored, matched
! against both: each projection goes to its own volume with a positive margin, at its own pose, and
! the mirrored one at the opposite view (ori%mirror2d: the image mirrored about y).
module simple_cavgs_state_match_tester
use simple_image,             only: image
use simple_oris,              only: oris
use simple_simple_volinterp,  only: reproject
use simple_imgarr_utils,      only: dealloc_imgarr
use simple_linalg,            only: rad2deg
use simple_ori_utils,         only: euler_dist, euler_inplrotdist, euler_mirror
use simple_cavgs_state_match, only: match_cavgs_to_states, reduce_state_matches, matched_pose
use simple_test_utils
implicit none
private
public :: run_all_cavgs_state_match_tests

contains

    subroutine run_all_cavgs_state_match_tests()
        write(*,'(A)') '**** running all class-average state matching tests ****'
        call test_reduce_state_matches()
        call test_matched_pose()
        call test_two_volumes()
    end subroutine run_all_cavgs_state_match_tests

    !> per class the best state and reference, the best correlation, and the margin over the best
    !! of another state; 0 margin with one state
    subroutine test_reduce_state_matches()
        integer, parameter :: NSTATES = 2, NSPACE = 3, NCLS = 2
        real    :: corr(NSTATES*NSPACE,NCLS), best_corr(NCLS), margin(NCLS)
        integer :: best_state(NCLS), best_ref(NCLS)
        write(*,'(A)') 'test_reduce_state_matches'
        ! class 1 fits state 2's second direction; class 2 state 1's third, narrowly
        corr(:,1) = [0.30, 0.50, 0.40,   0.20, 0.90, 0.60]
        corr(:,2) = [0.10, 0.20, 0.75,   0.70, 0.30, 0.10]
        call reduce_state_matches(corr, NSTATES, NSPACE, best_state, best_ref, best_corr, margin)
        call assert_int(2, best_state(1), 'class 1: state 2')
        call assert_int(5, best_ref(1),   'class 1: state 2''s second direction')
        call assert_real(0.90, best_corr(1), 1.e-6, 'class 1: its best correlation')
        call assert_real(0.40, margin(1),    1.e-6, 'class 1: the margin over state 1''s best')
        call assert_int(1, best_state(2), 'class 2: state 1')
        call assert_int(3, best_ref(2),   'class 2: state 1''s third direction')
        call assert_real(0.05, margin(2),    1.e-6, 'class 2: a narrow margin')
        ! one state of six directions: no margin
        call reduce_state_matches(corr, 1, NSTATES*NSPACE, best_state, best_ref, best_corr, margin)
        call assert_true(all(best_state == 1), 'one state: every class in it')
        call assert_int(5, best_ref(1), 'one state: class 1''s best direction')
        call assert_true(all(margin == 0.), 'one state: no margin')
    end subroutine test_reduce_state_matches

    !> the direction's first two angles, and 360 less the alignment's rotation (refine3D's convention)
    subroutine test_matched_pose()
        real :: euls(3)
        write(*,'(A)') 'test_matched_pose'
        euls = matched_pose([40., 75., 123.], 30.)
        call assert_real(40.,  euls(1), 1.e-4, 'the direction''s first angle')
        call assert_real(75.,  euls(2), 1.e-4, 'its second')
        call assert_real(330., euls(3), 1.e-4, '360 less the rotation, the direction''s own ignored')
        euls = matched_pose([10., 20., 0.], 0.)
        call assert_real(0.,   euls(3), 1.e-4, 'no rotation: 0, not 360')
    end subroutine test_matched_pose

    !> projections of two volumes, each matched to its own at its own pose; a mirrored one to the
    !! opposite view
    subroutine test_two_volumes()
        integer, parameter :: BOX = 40, NSPACE = 500, NTARG = 6
        real,    parameter :: SMPD = 3., LP = 9., TRS = 3., MSKDIAM = 84.
        real,    parameter :: ANG_TOL = 5., MIRR_ANG_TOL = 12., SHIFT_TOL = 1.
        integer, parameter :: IDIRS(5)   = [17, 101, 163, 255, 389]
        real,    parameter :: PSIS(5)    = [40., 200., 300., 90., 10.]
        integer, parameter :: EXPECTED(NTARG) = [1, 1, 1, 2, 2, 1]
        type(image), allocatable :: vols(:), targs(:), projs(:)
        type(oris)           :: eulspace, opose
        integer, allocatable :: states(:)
        real,    allocatable :: eulers(:,:), shifts(:,:), corrs(:), margins(:)
        real    :: euls_true(3,NTARG), euls(3)
        integer :: i
        character(len=96) :: msg
        write(*,'(A)') 'test_two_volumes'
        allocate(vols(2), targs(NTARG))
        call blob_volume(vols(1), BOX, SMPD, reshape([ 6., 0., 0.,  -4., 5., 1.,   0.,-3., 6.,  -2.,-4.,-5.], [3,4]),&
            &[1.0, 0.8, 0.6, 0.5], [4., 6., 3., 4.])
        call blob_volume(vols(2), BOX, SMPD, reshape([-5., 2., 3.,   3., 5.,-4.,   2.,-6.,-1.,   5., 1., 5.], [3,4]),&
            &[1.0, 0.7, 0.9, 0.5], [5., 4., 3., 4.])
        ! the true poses: directions of the matching grid, with an in-plane angle; the sixth target is
        ! the first mirrored, the projection from the opposite side
        call eulspace%new(NSPACE, is_ptcl=.false.)
        call eulspace%spiral()
        do i = 1,5
            euls = eulspace%get_euler(IDIRS(i))
            euls_true(:,i) = [euls(1), euls(2), PSIS(i)]
        enddo
        call euler_mirror(euls_true(:,1), euls_true(:,6))
        call eulspace%kill
        ! the first three projected from volume 1, the next two from volume 2
        call opose%new(3, is_ptcl=.false.)
        do i = 1,3
            call opose%set_euler(i, euls_true(:,i))
        enddo
        projs = reproject(vols(1), opose)
        do i = 1,3
            call targs(i)%copy(projs(i))
        enddo
        call targs(6)%copy(projs(1))
        call targs(6)%mirror('y')
        call dealloc_imgarr(projs)
        call opose%new(2, is_ptcl=.false.)
        do i = 1,2
            call opose%set_euler(i, euls_true(:,3+i))
        enddo
        projs = reproject(vols(2), opose)
        do i = 1,2
            call targs(3+i)%copy(projs(i))
        enddo
        call dealloc_imgarr(projs)
        call opose%kill
        call match_cavgs_to_states(targs, vols, NSPACE, LP, TRS, MSKDIAM, 1, states, eulers, shifts, corrs, margins)
        call assert_true(all(states == EXPECTED), 'every projection matched to the volume it came from')
        call assert_true(all(margins > 0.), 'every projection with a positive margin')
        do i = 1,5
            write(msg,'(A,I0,A,F7.2,A)') 'projection ', i, ': direction within ', ANG_TOL, ' degrees'
            call assert_true(rad2deg(euler_dist(eulers(:,i), euls_true(:,i))) <= ANG_TOL, trim(msg))
            write(msg,'(A,I0,A,F7.2,A)') 'projection ', i, ': in-plane angle within ', ANG_TOL, ' degrees'
            call assert_true(rad2deg(euler_inplrotdist(eulers(:,i), euls_true(:,i))) <= ANG_TOL, trim(msg))
            write(msg,'(A,I0,A)') 'projection ', i, ': no shift'
            call assert_true(sqrt(sum(shifts(:,i)**2)) <= SHIFT_TOL, trim(msg))
        enddo
        call assert_true(rad2deg(euler_dist(eulers(:,6), euls_true(:,6))) <= MIRR_ANG_TOL,&
            &'the mirrored projection: the opposite view, without a mirror search')
        call dealloc_imgarr(targs)
        call dealloc_imgarr(vols)
    end subroutine test_two_volumes

    !> a volume of Gaussian blobs at @p centres (voxels from the centre), amplitudes @p amps and
    !! variances @p vars (voxels squared)
    subroutine blob_volume( vol, box, smpd, centres, amps, vars )
        class(image), intent(inout) :: vol
        integer,      intent(in)    :: box
        real,         intent(in)    :: smpd, centres(:,:), amps(:), vars(:)
        real, allocatable :: rmat(:,:,:)
        real    :: r(3)
        integer :: x, y, z, j
        allocate(rmat(box,box,box), source=0.)
        do z = 1,box
            do y = 1,box
                do x = 1,box
                    r = real([x, y, z] - (box / 2 + 1))
                    do j = 1,size(amps)
                        rmat(x,y,z) = rmat(x,y,z) + amps(j) * exp(-sum((r - centres(:,j))**2) / (2. * vars(j)))
                    enddo
                enddo
            enddo
        enddo
        call vol%new([box,box,box], smpd)
        call vol%set_rmat(rmat, .false.)
    end subroutine blob_volume

end module simple_cavgs_state_match_tester
