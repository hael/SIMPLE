!@descr: unit tests of the class student's alignment after prediction (simple_class_student_align)
! White class averages on the band; each particle one of them rotated and shifted beyond the window, inside the disk
! (a small phantom written to the run directory for the builder and removed): a list of the planted class alone
! recovers the planted pose, a list of every class keeps the planted class from any place in it, and a team of four
! aligns as one thread does, bit for bit.
module simple_class_student_align_tester
use simple_core_module_api,     only: sp, dp
use simple_rnd,                 only: gasdev
use simple_cmdline,             only: cmdline
use simple_parameters,          only: parameters
use simple_string,              only: string
use simple_builder,             only: builder
use simple_image,               only: image
use simple_syslib,              only: del_file
use simple_matcher_smpl_and_lplims, only: set_bp_range3D
use simple_shiftmap_srch,       only: setup_shiftmap_srch
use simple_class_student_align, only: class_student_align
use simple_test_polar_fixture,  only: rotate_half_ring, shift_args
use simple_test_utils
!$ use omp_lib,                 only: omp_get_thread_num
implicit none
private
public :: run_all_class_student_align_tests

character(len=*), parameter :: PHANTOM = 'tmp_class_student_align_tester_phantom.mrc'
integer,          parameter :: NREFS  = 4      !< the class averages
integer,          parameter :: NPTCLS = 12     !< the particles: three per class
integer,          parameter :: NTEAM  = 4      !< the team of the thread check
real,             parameter :: TRS    = 2.     !< the window (pixels): the planted shifts lie beyond it
real,             parameter :: SHSIG  = 3.     !< the prior width (A): a disk of 9 A, 6 px at 1.5 A per pixel
real,             parameter :: SH_TOL = 0.25   !< the recovered shift's tolerance (pixels)

real, allocatable, target :: sigma2_fix(:,:)   !< the fixture's unit noise spectra: the calculator points at them

contains

    subroutine run_all_class_student_align_tests()
        write(*,'(A)') '**** running all class student alignment tests ****'
        call test_one_class_recovers_planted_pose()
        call test_list_keeps_planted_class()
        call test_team_as_one_thread()
        ! call report_summary()
    end subroutine run_all_class_student_align_tests

    !> a list of the planted class alone: its rotation, and the planted shift (5 px, beyond the 2 px window) as the increment
    subroutine test_one_class_recovers_planted_pose()
        type(builder),    target :: b
        type(parameters), target :: p
        type(cmdline)            :: cline
        type(class_student_align) :: aln
        integer, allocatable :: cls_true(:), inpl_true(:)
        real,    allocatable :: sh_true(:,:)
        real    :: corr, sh(2), sdev
        integer :: iptcl, cls, irot
        write(*,'(A)') 'test_one_class_recovers_planted_pose'
        call make_fixture(b, p, cline, cls_true, inpl_true, sh_true)
        call aln%new(p, b, NREFS)
        do iptcl = 1, NPTCLS
            call aln%align(iptcl, [cls_true(iptcl)], cls, irot, corr, sh, sdev)
            call assert_int(cls_true(iptcl), cls, 'a list of one class commits that class')
            call assert_int(inpl_true(iptcl), irot, 'one class recovers the planted rotation')
            call assert_true(all(abs(sh - sh_true(:,iptcl)) <= SH_TOL), 'one class recovers the planted shift beyond the window')
        enddo
        call aln%kill
        call kill_fixture(b, cline)
    end subroutine test_one_class_recovers_planted_pose

    !> lists of every class with the planted class first, in the middle and last: the planted class, rotation and shift;
    !! another class alone keeps that class at a lower score
    subroutine test_list_keeps_planted_class()
        type(builder),    target :: b
        type(parameters), target :: p
        type(cmdline)            :: cline
        type(class_student_align) :: aln
        integer, allocatable :: cls_true(:), inpl_true(:)
        real,    allocatable :: sh_true(:,:)
        real    :: corr, corr1, sh(2), sh1(2), sdev
        integer :: iptcl, cls, irot, cls1, irot1, iorder, refs(NREFS), other
        write(*,'(A)') 'test_list_keeps_planted_class'
        call make_fixture(b, p, cline, cls_true, inpl_true, sh_true)
        call aln%new(p, b, NREFS)
        do iptcl = 1, NPTCLS
            other = modulo(cls_true(iptcl), NREFS) + 1
            do iorder = 1, 3
                refs = list_with(cls_true(iptcl), iorder)
                call aln%align(iptcl, refs, cls, irot, corr, sh, sdev)
                call assert_int(cls_true(iptcl), cls, 'the list keeps the planted class wherever it holds it')
                call assert_int(inpl_true(iptcl), irot, 'the list recovers the planted rotation')
                call assert_true(all(abs(sh - sh_true(:,iptcl)) <= SH_TOL), 'the list recovers the planted shift')
            enddo
            call aln%align(iptcl, [other], cls1, irot1, corr1, sh1, sdev)
            call assert_int(other, cls1, 'a list of one class cannot leave it')
            call assert_true(corr1 < corr, 'another class aligns worse than the planted one')
        enddo
        call aln%kill
        call kill_fixture(b, cline)

    contains

        !> the classes with c first (1), in the middle (2) or last (3)
        function list_with( c, iorder ) result( l )
            integer, intent(in) :: c, iorder
            integer :: l(NREFS), j, n, others(NREFS-1)
            n = 0
            do j = 1, NREFS
                if( j == c ) cycle
                n = n + 1
                others(n) = j
            enddo
            select case(iorder)
                case(1)
                    l = [c, others]
                case(2)
                    l = [others(1:NREFS/2), c, others(NREFS/2+1:)]
                case default
                    l = [others, c]
            end select
        end function list_with

    end subroutine test_list_keeps_planted_class

    !> every particle aligned over the whole list by a team of NTEAM, one object per thread, against one thread
    subroutine test_team_as_one_thread()
        type(builder),    target :: b
        type(parameters), target :: p
        type(cmdline)            :: cline
        type(class_student_align) :: aln(NTEAM)
        integer, allocatable :: cls_true(:), inpl_true(:)
        real,    allocatable :: sh_true(:,:)
        real    :: corr1(NPTCLS), corrt(NPTCLS), sh1(2,NPTCLS), sht(2,NPTCLS), sdev1(NPTCLS), sdevt(NPTCLS)
        integer :: cls1(NPTCLS), clst(NPTCLS), irot1(NPTCLS), irott(NPTCLS), iptcl, i, ithr, refs(NREFS)
        write(*,'(A)') 'test_team_as_one_thread'
        call make_fixture(b, p, cline, cls_true, inpl_true, sh_true)
        do i = 1, NTEAM
            call aln(i)%new(p, b, NREFS)
        enddo
        refs = [(i, i = 1, NREFS)]
        do iptcl = 1, NPTCLS
            call aln(1)%align(iptcl, refs, cls1(iptcl), irot1(iptcl), corr1(iptcl), sh1(:,iptcl), sdev1(iptcl))
        enddo
        !$omp parallel do num_threads(NTEAM) default(shared) private(iptcl,ithr) schedule(dynamic) proc_bind(close)
        do iptcl = NPTCLS, 1, -1
            ithr = 1
            !$ ithr = omp_get_thread_num() + 1
            call aln(ithr)%align(iptcl, refs, clst(iptcl), irott(iptcl), corrt(iptcl), sht(:,iptcl), sdevt(iptcl))
        enddo
        !$omp end parallel do
        call assert_true(all(cls1 == clst) .and. all(irot1 == irott), 'a team aligns to the classes and rotations of one thread')
        call assert_true(all(transfer(corr1, 0, NPTCLS) == transfer(corrt, 0, NPTCLS)), 'a team scores as one thread, bit for bit')
        call assert_true(all(transfer(sh1, 0, 2*NPTCLS) == transfer(sht, 0, 2*NPTCLS)), 'a team shifts as one thread, bit for bit')
        call assert_true(all(transfer(sdev1, 0, NPTCLS) == transfer(sdevt, 0, NPTCLS)), 'a team reads the width of one thread')
        do i = 1, NTEAM
            call aln(i)%kill
        enddo
        call kill_fixture(b, cline)
    end subroutine test_team_as_one_thread

    ! ===== the fixture =====

    !> NREFS white class averages on the band and NPTCLS particles, particle i of class 1 + mod(i-1, NREFS), rotated by
    !! a grid rotation and shifted by a planted 5 px in a direction of its own, every particle not yet centred; the
    !! search reports the planted shift as the increment (the optimizer's shift rotated into the particle's axes)
    subroutine make_fixture( b, p, cline, cls_true, inpl_true, sh_true )
        type(builder),    target, intent(inout) :: b
        type(parameters), target, intent(inout) :: p
        type(cmdline),            intent(inout) :: cline
        integer, allocatable,     intent(out)   :: cls_true(:), inpl_true(:)
        real,    allocatable,     intent(out)   :: sh_true(:,:)
        complex(sp), allocatable :: refs(:,:,:), ptcl(:,:)
        real    :: ang, s(2)
        integer :: icls, iptcl, k, pp, pdim(3), kfromto(2), pftsz, nrots
        call set_fixed_seed(20261116)
        call write_phantom
        call cline%set('vol1',        PHANTOM)
        call cline%set('mskdiam',     36.)
        call cline%set('smpd',        1.5)
        call cline%set('lp',          6.)
        call cline%set('nptcls',      NPTCLS)
        call cline%set('nspace',      2)
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
        do icls = 1, NREFS
            do k = kfromto(1), kfromto(2)
                do pp = 1, pftsz
                    refs(pp,k,icls) = cmplx(gasdev(), gasdev())
                enddo
            enddo
            call b%pftc%set_ref_pft(icls, refs(:,:,icls), iseven=.true.)
        enddo
        call b%pftc%memoize_refs
        ! white noise: unit sigma2 in both halves
        allocate(sigma2_fix(kfromto(1):kfromto(2),2), source=1.)
        call b%pftc%assign_sigma2_noise(sigma2_fix)
        allocate(cls_true(NPTCLS), inpl_true(NPTCLS), sh_true(2,NPTCLS))
        do iptcl = 1, NPTCLS
            cls_true(iptcl)  = 1 + mod(iptcl - 1, NREFS)
            inpl_true(iptcl) = 1 + mod(7 * iptcl, nrots)
            ang = 2. * real(iptcl) + 0.3
            s   = 5. * [cos(ang), sin(ang)]
            call rotate_half_ring(pftsz, kfromto(2)-kfromto(1)+1, refs(:,kfromto(1):kfromto(2),cls_true(iptcl)), &
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
        integer, parameter :: BOX = 32
        type(image)        :: vol
        real, allocatable  :: rmat(:,:,:)
        integer :: i, j, k
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

end module simple_class_student_align_tester
