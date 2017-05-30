module simple_volpft_srch_tester
use simple_volpft_corrcalc, only: volpft_corrcalc
use simple_build,           only: build
use simple_params,          only: params
use simple_ori,             only: ori
use simple_image,           only: image
use simple_rnd,             only: ran3
use simple_math,            only: euclid
use simple_cmdline,         only: cmdline
use simple_projector,       only: projector
use simple_volpft_srch      ! singleton
use simple_defs             ! singleton
use simple_gridding         ! singleton
use simple_jiffys           ! singleton
implicit none

public :: exec_volpft_srch_test
private

! module global constants
integer, parameter    :: NTESTS=10, NPEAKS=3
real, parameter       :: TRS=5., SNR=0.1, LPLIM=20.0, SHERR_LIM=0.5, ROERR_LIM=2.0

! module global variables
type(params)          :: p
type(build)           :: b
type(projector)       :: vol_ref
type(volpft_corrcalc) :: vpftcc
type(cmdline)         :: cline_here
logical               :: verbose=.false.

contains

    subroutine exec_volpft_srch_test( cline, be_verbose )
        class(cmdline),    intent(inout) :: cline
        logical, optional, intent(in)    :: be_verbose
        call setup_testenv( cline, be_verbose )
        write(*,*) '****volpft_srch, init'
        ! call test_volpft_srch
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
        call cline_here%set('trs', TRS  )
        call cline_here%set('snr', SNR  )
        call cline_here%set('lp',  LPLIM)
        ! create parameters and build
        p = params(cline,checkdistr=.false.) ! constants & derived constants produced, mode=2
        call b%build_general_tbox(p,cline)   ! general objects built
        ! generate images
        ! deal with reference 
        call vol_ref%new([p%box,p%box,p%box], p%smpd)
        call vol_ref%read(p%vols(1))
        call vol_ref%add_gauran(p%snr)
        call vol_ref%mask(p%msk,'soft')
        call vol_ref%fwd_ft
        ! deal with target (randomly rotated version of vols(1))
        call vol_tmp%new([p%box,p%box,p%box], p%smpd)
        call vol_tmp%read(p%vols(1))
        call ranori%rnd_ori
        b%vol = rotvol(vol_tmp, ranori, p)
        call b%vol%add_gauran(p%snr)
        call b%vol%mask(p%msk,'soft')
        call b%vol%fwd_ft
        ! create correlator object
        call vpftcc%new(vol_ref,b%vol,p%hp,p%lp)
        ! initialise srch object
        call volpft_srch_init(vpftcc, TRS)
        call vol_tmp%kill
    end subroutine setup_testenv

    ! subroutine test_volpft_srch
    !     type(ori)   :: o_best, e, ranori
    !     real        :: corr_best, shvec(3), x, y, z, dist, sumdist, sherr, xf, yf, zf
    !     integer     :: itest
    !     sumdist = 0.
    !     sherr   = 0.
    !     do itest=1,NTESTS
    !         call progress(itest,NTESTS)
    !         call ranori%rnd_ori
    !         x = ran3()*2*TRS-TRS
    !         y = ran3()*2*TRS-TRS
    !         z = ran3()*2*TRS-TRS
    !         call ranori%set('x',x)
    !         call ranori%set('y',y)
    !         call ranori%set('z',z)
    !         call vpftcc%extract_ref(ranori)
    !         call vpftcc%shift_orig_ref([x,y,z])
    !         call volpft_6dimsrch(NPEAKS, corr_best, o_best)
    !         dist    = o_best.euldist.ranori
    !         sumdist = sumdist+dist
    !         xf = o_best%get('x')
    !         yf = o_best%get('y')
    !         zf = o_best%get('z')
    !         sherr = sherr+euclid([x,y,z],[-xf,-yf,-zf])
    !     end do
    !     dist  = sumdist/real(NTESTS)
    !     sherr = sherr/real(NTESTS*3)
    !     if( verbose ) write(*,'(a,1x,f5.2)') 'SHIFT ERROR (IN PIXELS ): ', sherr
    !     if( verbose ) write(*,'(a,1x,f5.2)') 'ROT   ERROR (IN DEGREES): ', dist
    !     if( .not. test_passed() ) stop '****volpft_srch_tester FAILURE volpft_srch :: volpft_6dimsrch'

    !     contains

    !         function test_passed() result( passed )
    !             logical :: passed
    !             passed = .false.
    !             if( sherr < SHERR_LIM .and. dist < ROERR_LIM ) passed = .true.
    !         end function test_passed

    ! end subroutine test_volpft_srch

end module simple_volpft_srch_tester
