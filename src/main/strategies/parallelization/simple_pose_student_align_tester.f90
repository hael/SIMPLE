!@descr: unit tests of the pose student's alignment after prediction (simple_pose_student_align)
! White references on the band, two mirror pairs (1 with 3, 2 with 4); each particle one of them rotated and shifted
! beyond the window, inside the disk (a small phantom written to the run directory for the builder and removed): the
! planted reference recovers the planted pose, a direction whose partner is the planted reference finds it, a list of
! every direction keeps the planted one wherever it holds it, a reference is searched once, and a team of four aligns
! as one thread does, bit for bit.
module simple_pose_student_align_tester
use simple_core_module_api,    only: sp, dp
use simple_rnd,                only: gasdev
use simple_cmdline,            only: cmdline
use simple_parameters,         only: parameters
use simple_string,             only: string
use simple_builder,            only: builder
use simple_image,              only: image
use simple_syslib,             only: del_file
use simple_matcher_smpl_and_lplims, only: set_bp_range3D
use simple_shiftmap_srch,      only: setup_shiftmap_srch
use simple_pose_student_align, only: pose_student_align, pose_align_result
use simple_test_polar_fixture, only: rotate_half_ring, shift_args
use simple_test_utils
!$ use omp_lib,                only: omp_get_thread_num
implicit none
private
public :: run_all_pose_student_align_tests

character(len=*), parameter :: PHANTOM     = 'tmp_pose_student_align_tester_phantom.mrc'
integer,          parameter :: NREFS       = 4              !< the references: two mirror pairs
integer,          parameter :: NPTCLS      = 12             !< the particles: three per reference
integer,          parameter :: NTEAM       = 4              !< the team of the thread check
integer,          parameter :: MIRR(NREFS) = [3, 4, 1, 2]   !< each reference's mirror partner
real,             parameter :: TRS         = 2.             !< the window (pixels): the planted shifts lie beyond it
real,             parameter :: SHSIG       = 3.             !< the prior width (A): a disk of 9 A, 6 px at 1.5 A per pixel
real,             parameter :: SH_TOL      = 0.25           !< the recovered shift's tolerance (pixels)

real, allocatable, target :: sigma2_fix(:,:)   !< the fixture's unit noise spectra: the calculator points at them

contains

    subroutine run_all_pose_student_align_tests()
        write(*,'(A)') '**** running all pose student alignment tests ****'
        call test_planted_reference_recovers_pose()
        call test_partner_of_a_direction_searched()
        call test_cone_keeps_planted_direction()
        call test_reference_searched_once()
        call test_team_as_one_thread()
    end subroutine run_all_pose_student_align_tests

    !> the planted reference alone: its rotation and the planted shift (5 px, beyond the 2 px window) as the increment;
    !! its partner is searched beside it and loses
    subroutine test_planted_reference_recovers_pose()
        type(builder),    target      :: b
        type(parameters), target      :: p
        type(cmdline)                 :: cline
        type(pose_student_align)      :: aln
        type(pose_align_result)       :: res
        integer,          allocatable :: ref_true(:), inpl_true(:)
        real,             allocatable :: sh_true(:,:)
        integer                       :: iptcl
        write(*,'(A)') 'test_planted_reference_recovers_pose'
        call make_fixture(b, p, cline, ref_true, inpl_true, sh_true)
        call aln%new(p, b, NREFS)
        do iptcl = 1, NPTCLS
            call aln%align(iptcl, [ref_true(iptcl)], MIRR, res)
            call assert_int(ref_true(iptcl), res%iref, 'the planted direction beats its partner')
            call assert_int(inpl_true(iptcl), res%irot, 'the planted rotation is recovered')
            call assert_true(all(abs(res%sh - sh_true(:,iptcl)) <= SH_TOL), 'the planted shift beyond the window is recovered')
        enddo
        call aln%kill
        call kill_fixture(b, cline)
    end subroutine test_planted_reference_recovers_pose

    !> a list naming only the planted reference's mirror partner (the cone names one direction of a pair) finds the
    !! planted reference, its rotation and shift
    subroutine test_partner_of_a_direction_searched()
        type(builder),    target      :: b
        type(parameters), target      :: p
        type(cmdline)                 :: cline
        type(pose_student_align)      :: aln
        type(pose_align_result)       :: res
        integer,          allocatable :: ref_true(:), inpl_true(:)
        real,             allocatable :: sh_true(:,:)
        integer                       :: iptcl
        write(*,'(A)') 'test_partner_of_a_direction_searched'
        call make_fixture(b, p, cline, ref_true, inpl_true, sh_true)
        call aln%new(p, b, NREFS)
        do iptcl = 1, NPTCLS
            call aln%align(iptcl, [MIRR(ref_true(iptcl))], MIRR, res)
            call assert_int(ref_true(iptcl), res%iref, 'the partner of the listed direction is searched and wins')
            call assert_int(inpl_true(iptcl), res%irot, 'the partner''s rotation is recovered')
            call assert_true(all(abs(res%sh - sh_true(:,iptcl)) <= SH_TOL), 'the partner''s shift is recovered')
        enddo
        call aln%kill
        call kill_fixture(b, cline)
    end subroutine test_partner_of_a_direction_searched

    !> a list of every pair's direction (the cone) in either order keeps the planted reference; the other pair alone (a
    !! wrong mode) stays in that pair at a lower score
    subroutine test_cone_keeps_planted_direction()
        type(builder),    target      :: b
        type(parameters), target      :: p
        type(cmdline)                 :: cline
        type(pose_student_align)      :: aln
        type(pose_align_result)       :: res, res_other
        integer,          allocatable :: ref_true(:), inpl_true(:)
        real,             allocatable :: sh_true(:,:)
        integer                       :: iptcl, iorder, other
        write(*,'(A)') 'test_cone_keeps_planted_direction'
        call make_fixture(b, p, cline, ref_true, inpl_true, sh_true)
        call aln%new(p, b, NREFS)
        do iptcl = 1, NPTCLS
            do iorder = 1, 2
                if( iorder == 1 )then
                    call aln%align(iptcl, [1, 2], MIRR, res)
                else
                    call aln%align(iptcl, [4, 3], MIRR, res)
                endif
                call assert_int(ref_true(iptcl), res%iref, 'the cone keeps the planted direction wherever it holds it')
                call assert_int(inpl_true(iptcl), res%irot, 'the cone recovers the planted rotation')
                call assert_true(all(abs(res%sh - sh_true(:,iptcl)) <= SH_TOL), 'the cone recovers the planted shift')
            enddo
            ! the other pair: the direction after the planted one's pair
            other = 1 + mod(ref_true(iptcl), 2)
            call aln%align(iptcl, [other], MIRR, res_other)
            call assert_true(res_other%iref == other .or. res_other%iref == MIRR(other), 'a list of one pair cannot leave it')
            call assert_true(res_other%corr < res%corr, 'another pair aligns worse than the planted one')
        enddo
        call aln%kill
        call kill_fixture(b, cline)
    end subroutine test_cone_keeps_planted_direction

    !> a list that names a direction twice and both members of a pair searches each reference once: the result of the
    !! pair named by one direction, bit for bit; a list beyond nmax references aligns to none (the caller's status)
    subroutine test_reference_searched_once()
        type(builder),    target      :: b
        type(parameters), target      :: p
        type(cmdline)                 :: cline
        type(pose_student_align)      :: aln
        type(pose_align_result)       :: res1, res2
        integer,          allocatable :: ref_true(:), inpl_true(:)
        real,             allocatable :: sh_true(:,:)
        integer                       :: iptcl
        write(*,'(A)') 'test_reference_searched_once'
        call make_fixture(b, p, cline, ref_true, inpl_true, sh_true)
        ! two references: a pair, however often named, fits
        call aln%new(p, b, 2)
        do iptcl = 1, NPTCLS
            call aln%align(iptcl, [ref_true(iptcl)], MIRR, res1)
            call aln%align(iptcl, [ref_true(iptcl), MIRR(ref_true(iptcl)), ref_true(iptcl)], MIRR, res2)
            call assert_true(res1%iref == res2%iref .and. res1%irot == res2%irot, 'a pair named thrice aligns as once')
            call assert_true(transfer(res1%corr, 0) == transfer(res2%corr, 0) .and. &
                &all(transfer(res1%sh, 0, 2) == transfer(res2%sh, 0, 2)), 'a pair named thrice scores as once, bit for bit')
            call aln%align(iptcl, [1, 2], MIRR, res2)
            call assert_int(0, res2%iref, 'two pairs exceed a list of two references: no reference, the caller''s status')
        enddo
        call aln%kill
        call kill_fixture(b, cline)
    end subroutine test_reference_searched_once

    !> every particle aligned over the whole cone by a team of NTEAM, one object per thread, against one thread
    subroutine test_team_as_one_thread()
        type(builder),    target      :: b
        type(parameters), target      :: p
        type(cmdline)                 :: cline
        type(pose_student_align)      :: aln(NTEAM)
        type(pose_align_result)       :: res1(NPTCLS), rest(NPTCLS)
        integer,          allocatable :: ref_true(:), inpl_true(:)
        real,             allocatable :: sh_true(:,:)
        integer                       :: iptcl, i, ithr
        write(*,'(A)') 'test_team_as_one_thread'
        call make_fixture(b, p, cline, ref_true, inpl_true, sh_true)
        do i = 1, NTEAM
            call aln(i)%new(p, b, NREFS)
        enddo
        do iptcl = 1, NPTCLS
            call aln(1)%align(iptcl, [1, 2], MIRR, res1(iptcl))
        enddo
        !$omp parallel do num_threads(NTEAM) default(shared) private(iptcl,ithr) schedule(dynamic) proc_bind(close)
        do iptcl = NPTCLS, 1, -1
            ithr = 1
            !$ ithr = omp_get_thread_num() + 1
            call aln(ithr)%align(iptcl, [1, 2], MIRR, rest(iptcl))
        enddo
        !$omp end parallel do
        call assert_true(all(res1%iref == rest%iref) .and. all(res1%irot == rest%irot), &
            &'a team aligns to the references and rotations of one thread')
        call assert_true(all(transfer(res1%corr, 0, NPTCLS) == transfer(rest%corr, 0, NPTCLS)), 'a team scores as one thread, bit for bit')
        do iptcl = 1, NPTCLS
            call assert_true(all(transfer(res1(iptcl)%sh, 0, 2) == transfer(rest(iptcl)%sh, 0, 2)), &
                &'a team shifts as one thread, bit for bit')
            call assert_true(transfer(res1(iptcl)%sdev, 0) == transfer(rest(iptcl)%sdev, 0), 'a team reads the width of one thread')
        enddo
        do i = 1, NTEAM
            call aln(i)%kill
        enddo
        call kill_fixture(b, cline)
    end subroutine test_team_as_one_thread

    ! ===== the fixture =====

    !> NREFS white references on the band and NPTCLS particles, particle i of reference 1 + mod(i-1, NREFS), rotated by
    !! a grid rotation and shifted by a planted 5 px in a direction of its own, every particle not yet centred; the
    !! search reports the planted shift as the increment (the optimizer's shift rotated into the particle's axes)
    subroutine make_fixture( b, p, cline, ref_true, inpl_true, sh_true )
        type(builder),            target, intent(inout)        :: b
        type(parameters),         target, intent(inout)        :: p
        type(cmdline),            intent(inout)                :: cline
        integer,                  allocatable,     intent(out) :: ref_true(:), inpl_true(:)
        real,                     allocatable,     intent(out) :: sh_true(:,:)
        complex(sp), allocatable :: refs(:,:,:), ptcl(:,:)
        real                     :: ang, s(2)
        integer                  :: iref, iptcl, k, pp, pdim(3), kfromto(2), pftsz, nrots
        call set_fixed_seed(20261117)
        call write_phantom
        call cline%set('vol1',        PHANTOM)
        call cline%set('mskdiam',     36.)
        call cline%set('smpd',        1.5)
        call cline%set('lp',          6.)
        call cline%set('nptcls',      NPTCLS)
        call cline%set('nspace',      NREFS)
        call cline%set('ctf',         'no')
        call cline%set('objfun',      'euclid')
        call cline%set('trs',         TRS)
        call cline%set('shift_sigma', SHSIG)
        call cline%set('prob_beta',   1.)
        call cline%set('nthr',        NTEAM)
        call cline%check
        call b%init_params_and_build_strategy3D_tbox(cline, p)
        call set_bp_range3D(p, b, cline)
        call b%pftc%new(p, NREFS, [1,NPTCLS], p%kfromto)
        nrots   = b%pftc%get_nrots()
        pftsz   = b%pftc%get_pftsz()
        kfromto = b%pftc%get_kfromto()
        pdim    = b%pftc%get_pdim_interp()
        allocate(refs(pdim(1),pdim(2):pdim(3),NREFS), ptcl(pdim(1),pdim(2):pdim(3)), source=cmplx(0.,0.))
        do iref = 1, NREFS
            do k = kfromto(1), kfromto(2)
                do pp = 1, pftsz
                    refs(pp,k,iref) = cmplx(gasdev(), gasdev())
                enddo
            enddo
            call b%pftc%set_ref_pft(iref, refs(:,:,iref), iseven=.true.)
        enddo
        call b%pftc%memoize_refs
        ! white noise: unit sigma2 in both halves
        allocate(sigma2_fix(kfromto(1):kfromto(2),2), source=1.)
        call b%pftc%assign_sigma2_noise(sigma2_fix)
        allocate(ref_true(NPTCLS), inpl_true(NPTCLS), sh_true(2,NPTCLS))
        do iptcl = 1, NPTCLS
            ref_true(iptcl)  = 1 + mod(iptcl - 1, NREFS)
            inpl_true(iptcl) = 1 + mod(11 * iptcl, nrots)
            ang = 2. * real(iptcl) + 0.7
            s   = 5. * [cos(ang), sin(ang)]
            call rotate_half_ring(pftsz, kfromto(2)-kfromto(1)+1, refs(:,kfromto(1):kfromto(2),ref_true(iptcl)), &
                &inpl_true(iptcl), ptcl(:,kfromto(1):kfromto(2)))
            call shift_ptcl(s, ptcl(:,kfromto(1):kfromto(2)))
            call b%pftc%set_eo(iptcl, .true.)
            call b%pftc%set_ptcl_pft(iptcl, ptcl)
            call b%spproj_field%set_state(iptcl, 1)
            call b%spproj_field%set(iptcl, 'shsig',  -1.)
            call b%spproj_field%set(iptcl, 'shincr',  0.)
            sh_true(:,iptcl) = s
        enddo
        call b%pftc%memoize_ptcls
        call setup_shiftmap_srch(p, b)

    contains

        !> a half ring shifted by s with the calculator's phase, as a particle at shift s
        subroutine shift_ptcl( s_here, pft )
            real,        intent(in)    :: s_here(2)
            complex(sp), intent(inout) :: pft(pftsz,kfromto(1):kfromto(2))
            real(dp) :: ax(pftsz,kfromto(1):kfromto(2)), ay(pftsz,kfromto(1):kfromto(2))
            call shift_args(pftsz, kfromto, p%box, ax, ay)
            pft = pft * cmplx(cos(ax * s_here(1) + ay * s_here(2)), sin(ax * s_here(1) + ay * s_here(2)), kind=sp)
        end subroutine shift_ptcl

    end subroutine make_fixture

    subroutine kill_fixture( b, cline )
        type(builder), intent(inout) :: b
        type(cmdline), intent(inout) :: cline
        call b%pftc%kill
        call b%kill_strategy3D_tbox
        call b%kill_general_tbox
        call cline%kill
        call del_file(PHANTOM)
        if( allocated(sigma2_fix) ) deallocate(sigma2_fix)
    end subroutine kill_fixture

    !> a 32^3 two-blob phantom: the builder's volume (the alignment reads only the calculator)
    subroutine write_phantom
        integer, parameter   :: BOX = 32
        type(image)          :: vol
        real,    allocatable :: rmat(:,:,:)
        integer              :: i, j, k
        allocate(rmat(BOX,BOX,BOX))
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    rmat(i,j,k) = exp(-real((i-14)**2 + (j-18)**2 + (k-16)**2)/18.) + &
                        &0.5*exp(-real((i-20)**2 + (j-13)**2 + (k-19)**2)/12.)
                end do
            end do
        end do
        call vol%new([BOX,BOX,BOX], 1.5)
        call vol%set_rmat(rmat, .false.)
        call vol%write(string(PHANTOM), del_if_exists=.true.)
        call vol%kill
    end subroutine write_phantom

end module simple_pose_student_align_tester
