module simple_volpft_srch_tester
include 'simple_lib.f08'
use simple_builder,    only: builder
use simple_parameters, only: parameters
use simple_ori,        only: ori
use simple_image,      only: image
use simple_cmdline,    only: cmdline
use simple_projector,  only: projector
use simple_volpft_srch ! singleton
implicit none

public :: exec_volpft_srch_test
private
#include "simple_local_flags.inc"

! module global constants
integer, parameter :: NTESTS=10, NPEAKS=3
real,    parameter :: SNR=0.5, LPLIM=20.0, ROERR_LIM=2.0

! module global variables
type(parameters) :: p
type(builder)    :: b
type(projector)  :: vol_ref
type(cmdline)    :: cline_here

contains

    subroutine exec_volpft_srch_test( cline, be_verbose )
        class(cmdline),    intent(inout) :: cline
        logical, optional, intent(in)    :: be_verbose
        call setup_testenv( cline, be_verbose )
        write(*,*) '****volpft_srch, init'
        call test_volpft_srch
        write(*,*) '****volpft_srch, completed'
    end subroutine exec_volpft_srch_test

    subroutine setup_testenv( cline, be_verbose )
        use simple_projector_hlev, only: rotvol
        class(cmdline),    intent(inout) :: cline
        logical, optional, intent(in)    :: be_verbose
        type(ori)   :: ranori
        type(image) :: vol_tmp
        verbose = .false.
        if( present(be_verbose) ) verbose = be_verbose
        ! it is assumed that vol1, smpd, msk are part of the inputted command line
        ! setting the remainder of the command line up in here
        cline_here = cline
        call cline_here%set('snr', SNR  )
        call cline_here%set('lp',  LPLIM)
        ! create parameters and build
        call p%new(cline) ! constants & derived constants produced, mode=2
        call b%build_general_tbox(p, cline)   ! general objects built
        ! generate images
        ! deal with reference
        call vol_ref%new([p%box,p%box,p%box], p%smpd)
        call vol_ref%read(p%vols(1))
        call vol_ref%add_gauran(SNR)
        call vol_ref%mask(p%msk,'soft')
        call vol_ref%write('vol_ref.mrc')
        call vol_ref%fft()
        ! deal with target (randomly rotated version of vols(1))
        call vol_tmp%new([p%box,p%box,p%box], p%smpd)
        call vol_tmp%read(p%vols(1))
        call ranori%new
        ! call ranori%rnd_ori
        call ranori%set_euler([20.,40.,160.])
        call b%vol%copy( rotvol(vol_tmp, ranori))
        call b%vol%add_gauran(SNR)
        call b%vol%mask(p%msk,'soft')
        call b%vol%write('rotated.mrc')
        call b%vol%fft()
        call volpft_srch_init(vol_ref,b%vol,p%hp,p%lp)
        call vol_tmp%kill
    end subroutine setup_testenv

    subroutine test_volpft_srch
        type(ori) :: o_best, ranori
        real      :: corr_best, dist, sumdist, sumcorr, corr
        integer   :: itest
        sumdist = 0.
        sumcorr = 0.


        o_best  = volpft_srch_minimize()
        call o_best%print_ori

        ! call ranori%new
        ! do itest=1,NTESTS
        !     call progress(itest,NTESTS)
        !     call ranori%rnd_ori
        !     o_best  = volpft_srch_minimize_eul()
        !     corr    = o_best%get('corr')
        !     sumcorr = sumcorr + corr
        !     dist    = o_best.euldist.ranori
        !     sumdist = sumdist + dist
        ! end do
        ! dist = sumdist/real(NTESTS)
        ! corr = sumcorr/real(NTESTS)
        ! write(*,'(a,1x,f5.2)') 'ROT ERROR (IN DEGREES): ', dist
        ! write(*,'(a,1x,f7.4)') 'CORR: ', corr
        ! if( .not. test_passed() ) THROW_HARD('****volpft_srch_tester FAILURE volpft_srch :: volpft_6dimsrch')
        !
        ! contains
        !
        !     function test_passed() result( passed )
        !         logical :: passed
        !         passed = .false.
        !         if( dist < ROERR_LIM ) passed = .true.
        !     end function test_passed

    end subroutine test_volpft_srch

end module simple_volpft_srch_tester
