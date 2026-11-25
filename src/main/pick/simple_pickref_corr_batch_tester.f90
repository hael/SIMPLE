!@descr: unit tests of simple_pickref_corr_batch: the batched Pearson correlation of reference picking against a brute force, its refused windows and references, and a team of three against one thread
module simple_pickref_corr_batch_tester
use simple_test_utils,         only: assert_true, set_fixed_seed
use simple_defs,               only: dp
use simple_image,              only: image
use simple_rnd,                only: gasdev, ran3
use simple_pickref_corr_batch, only: pickref_corr_batch
implicit none
private
public :: run_all_pickref_corr_batch_tests

integer,  parameter :: T_NX      = 96, T_NY = 80, T_BOX = 12, T_NREFS = 3, T_NFIX = 8, T_NPOS = 40, T_BATCH = 7
integer,  parameter :: T_YDRAW   = 44           !< the drawn windows' last corner row: they stay clear of the flat block
real,     parameter :: T_NOISE   = 0.1          !< the micrograph's noise relative to the planted unit-variance references
real(dp), parameter :: TOL_SCORE = 1.d-4        !< single-precision products of T_BOX**2 terms against the double sums
real,     parameter :: PLANT_MIN = 0.98         !< a planted reference under T_NOISE correlates at 1/sqrt(1 + 0.01)

contains

    subroutine run_all_pickref_corr_batch_tests()
        write(*,'(A)') '**** running all pickref corr batch tests ****'
        call test_scores()
        call test_team_identity()
    end subroutine run_all_pickref_corr_batch_tests

    !> A noise micrograph with reference 1 planted at (8,8) and (60,8), reference 2 (refused) at (34,40) and a flat
    !! block in the corner from (72,56); three zero-mean unit-variance references; positions 1-8 fixed (the plants, an
    !! unrelated window, the flat block, three outside the micrograph), the rest drawn above the flat block
    subroutine make_fixture( mic, refs, valid, pos, micn )
        type(image), intent(inout) :: mic, refs(T_NREFS)
        logical,     intent(out)   :: valid(T_NREFS)
        integer,     intent(out)   :: pos(T_NPOS,2)
        real(dp),    intent(out)   :: micn(T_NX,T_NY)
        real    :: rmat(T_NX,T_NY,1), rref(T_BOX,T_BOX,1,T_NREFS)
        integer :: iref, i, j, ip
        call set_fixed_seed(20261123)
        do iref = 1,T_NREFS
            do j = 1,T_BOX
                do i = 1,T_BOX
                    rref(i,j,1,iref) = gasdev()
                enddo
            enddo
            rref(:,:,1,iref) = rref(:,:,1,iref) - sum(rref(:,:,1,iref)) / real(T_BOX**2)
            rref(:,:,1,iref) = rref(:,:,1,iref) / sqrt(sum(rref(:,:,1,iref)**2) / real(T_BOX**2))
            call refs(iref)%new([T_BOX,T_BOX,1], 1.)
            call refs(iref)%set_rmat(rref(:,:,:,iref), .false.)
        enddo
        valid = [.true., .false., .true.]
        do j = 1,T_NY
            do i = 1,T_NX
                rmat(i,j,1) = T_NOISE * gasdev()
            enddo
        enddo
        rmat( 9:20,  9:20,1) = rmat( 9:20,  9:20,1) + rref(:,:,1,1)
        rmat(61:72,  9:20,1) = rmat(61:72,  9:20,1) + rref(:,:,1,1)
        rmat(35:46, 41:52,1) = rmat(35:46, 41:52,1) + rref(:,:,1,2)
        rmat(73:T_NX, 57:T_NY,1) = 0.5
        call mic%new([T_NX,T_NY,1], 1.)
        call mic%set_rmat(rmat, .false.)
        ! the micrograph as the correlation sees it: mean removed, divided by the largest magnitude
        micn = real(rmat(:,:,1),dp) - sum(real(rmat(:,:,1),dp)) / real(T_NX*T_NY,dp)
        micn = micn / maxval(abs(micn))
        pos(1:T_NFIX,1) = [ 8, 60, 34, 30, 76, -1, 90,  5]
        pos(1:T_NFIX,2) = [ 8,  8, 40, 10, 60,  5,  5, 75]
        do ip = T_NFIX+1,T_NPOS
            pos(ip,1) = min(T_NX - T_BOX, int(ran3() * real(T_NX - T_BOX + 1)))
            pos(ip,2) = min(T_YDRAW, int(ran3() * real(T_YDRAW + 1)))
        enddo
    end subroutine make_fixture

    !> the best Pearson correlation over the accepted references of the window at corner pos, in double precision
    real(dp) function brute_score( micn, refs, valid, pos )
        real(dp),    intent(in) :: micn(T_NX,T_NY)
        type(image), intent(in) :: refs(T_NREFS)
        logical,     intent(in) :: valid(T_NREFS)
        integer,     intent(in) :: pos(2)
        real, allocatable :: r(:,:,:)
        real(dp)          :: w(T_BOX,T_BOX), sd, s
        integer           :: iref
        w  = micn(pos(1)+1:pos(1)+T_BOX, pos(2)+1:pos(2)+T_BOX)
        w  = w - sum(w) / real(T_BOX**2,dp)
        sd = sqrt(sum(w**2) / real(T_BOX**2,dp))
        brute_score = -huge(1._dp)
        do iref = 1,T_NREFS
            if( .not. valid(iref) ) cycle
            r = refs(iref)%get_rmat()
            s = sum(real(r(:,:,1),dp) * w) / (real(T_BOX**2,dp) * sd)
            brute_score = max(brute_score, s)
        enddo
    end function brute_score

    !> every window inside the micrograph scores the brute force's correlation over the accepted references; a planted
    !! accepted reference correlates, an unrelated window and a refused reference's plant do not, and a flat window or
    !! one reaching outside the micrograph scores -1
    subroutine test_scores()
        type(pickref_corr_batch) :: corr
        type(image)              :: mic, refs(T_NREFS)
        logical                  :: valid(T_NREFS)
        integer                  :: pos(T_NPOS,2), ip, iref
        real                     :: scores(T_NPOS)
        real(dp)                 :: micn(T_NX,T_NY), err
        write(*,'(A)') 'test_scores'
        call make_fixture(mic, refs, valid, pos, micn)
        call corr%new(mic, refs, valid, T_BATCH)
        call corr%score_positions(pos, scores)
        err = 0._dp
        do ip = 1,T_NPOS
            if( ip >= 5 .and. ip <= T_NFIX ) cycle
            err = max(err, abs(real(scores(ip),dp) - brute_score(micn, refs, valid, pos(ip,:))))
        enddo
        call assert_true(err < TOL_SCORE, 'a window scores the best Pearson correlation over the accepted references')
        call assert_true(scores(1) > PLANT_MIN .and. scores(2) > PLANT_MIN, 'a planted accepted reference correlates')
        call assert_true(scores(3) < 0.5, 'the plant of a refused reference does not correlate')
        call assert_true(scores(4) < 0.5, 'an unrelated window does not correlate')
        call assert_true(scores(5) == -1., 'a flat window scores -1')
        call assert_true(all(scores(6:8) == -1.), 'a window reaching outside the micrograph scores -1')
        call corr%kill
        call mic%kill
        do iref = 1,T_NREFS
            call refs(iref)%kill
        enddo
    end subroutine test_scores

    !> batches of T_BATCH positions on a team of three give the scores of one thread, bit for bit
    subroutine test_team_identity()
        !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
        type(pickref_corr_batch) :: corr
        type(image)              :: mic, refs(T_NREFS)
        logical                  :: valid(T_NREFS)
        integer                  :: pos(T_NPOS,2), nthr_prev, iref
        real                     :: scores1(T_NPOS), scoresn(T_NPOS)
        real(dp)                 :: micn(T_NX,T_NY)
        write(*,'(A)') 'test_team_identity'
        call make_fixture(mic, refs, valid, pos, micn)
        call corr%new(mic, refs, valid, T_BATCH)
        nthr_prev = 1
        !$ nthr_prev = omp_get_max_threads()
        !$ call omp_set_num_threads(1)
        call corr%score_positions(pos, scores1)
        !$ call omp_set_num_threads(3)
        call corr%score_positions(pos, scoresn)
        !$ call omp_set_num_threads(nthr_prev)
        call assert_true(all(transfer(scoresn, [0]) == transfer(scores1, [0])), &
            &'batches on a team of three give the scores of one thread, bit for bit')
        call corr%kill
        call mic%kill
        do iref = 1,T_NREFS
            call refs(iref)%kill
        enddo
    end subroutine test_team_identity

end module simple_pickref_corr_batch_tester
