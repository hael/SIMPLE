module simple_inpl_srch_tester
use simple_image,            only: image
use simple_build,            only: build
use simple_params,           only: params
use simple_cmdline,          only: cmdline
use simple_simulator,        only: simimg
use simple_ori,              only: ori
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_rnd,              only: ran3
use simple_oris,             only: oris
use simple_ctf,              only: ctf
use simple_jiffys            ! singleton
use simple_hadamard_common   ! singleton
use simple_pftcc_inplsrch    ! singleton
use simple_pftcc_shsrch      ! singleton
use simple_math              ! singleton
use simple_defs              ! singleton
implicit none

public :: exec_inpl_srch_test
private

! module global constants
integer, parameter           :: NPROJS=100
real, parameter              :: KV=300., CS=2.7, FRACA=0.07, TRS=5.0, DFX=2.2, DFY=2.5, ANGAST=30., BFAC=50.
real, parameter              :: SNR=0.2, SNR_PINK=SNR/0.2, SNR_DETECTOR=SNR/0.8, LPLIM=8.
real, parameter              :: ROERR_LIM=5., SHERR_LIM=0.5, SHSHERR_LIM=0.8
character(len=32), parameter :: refsname  = 'inpl_srch_test_refs.mrc'
character(len=32), parameter :: ptclsname = 'inpl_srch_test_ptcls.mrc'
character(len=32), parameter :: orisname  = 'inpl_srch_test_oris.txt'

! module global variables
type(polarft_corrcalc)   :: pftcc
type(params)             :: p
type(build)              :: b
type(oris)               :: o_refs, o_ptcls
type(image), allocatable :: imgs_refs(:), imgs_ptcls(:)
type(ori)                :: orientation, orientation_opt
type(ctf)                :: tfun
real                     :: shlims(2,2), crxy(4), cxy(3)
real                     :: sherr, roerr, sherr_avg, roerr_avg
type(cmdline)            :: cline_here
logical                  :: verbose=.false.

contains

    subroutine exec_inpl_srch_test( cline, be_verbose )
        class(cmdline),    intent(inout) :: cline
        logical, optional, intent(in)    :: be_verbose
        call setup_testenv( cline, be_verbose )
        write(*,*) '****inpl_srch_test, init'
        call test_pftcc_inplsrch
        call test_pftcc_shsrch
        write(*,*) '****inpl_srch_test, completed'
    end subroutine exec_inpl_srch_test

    subroutine setup_testenv( cline, be_verbose )
        class(cmdline),    intent(inout) :: cline
        logical, optional, intent(in)    :: be_verbose
        integer :: iproj
        verbose = .false.
        if( present(be_verbose) ) verbose = be_verbose
        ! it is assumed that vol1, smpd, msk are part of the inputted command line
        ! setting the remainder of the command line up in here
        cline_here = cline
        call cline_here%set('ncls',    real(NPROJS))
        call cline_here%set('nptcls',  real(NPROJS))
        call cline_here%set('lp',      LPLIM       )
        ! generate orientations
        call o_refs%new(NPROJS)
        call o_refs%spiral
        o_ptcls = o_refs
        call o_ptcls%rnd_inpls(TRS)
        ! create parameters and build
        p = params(cline_here)                   ! parameters generated
        call b%build_general_tbox(p, cline_here) ! general objects built
        ! generate images
        call b%vol%read(p%vols(1))
        imgs_refs  = b%proj%projvol(b%vol, o_refs,  p)
        imgs_ptcls = b%proj%projvol(b%vol, o_ptcls, p)
        do iproj=1,NPROJS
            call imgs_ptcls(iproj)%shift(o_ptcls%get(iproj,'x'),o_ptcls%get(iproj,'y')) 
            call imgs_refs(iproj)%write(refsname, iproj)
            call imgs_ptcls(iproj)%write(ptclsname, iproj)
        end do
        ! set resolution range
        p%kfromto(1) = 2
        p%kfromto(2) = b%img%get_find(p%lp)
        ! build corr calculator
        call pftcc%new(1, [1,1], [p%box,p%box,1], p%kfromto, p%ring2, 'yes')
        ! set shift limits
        shlims(:,1) = -TRS
        shlims(:,2) =  TRS
        ! generate CTF object
        tfun = ctf(p%smpd, KV,CS, FRACA)
    end subroutine setup_testenv

    subroutine test_pftcc_inplsrch
        integer :: iproj
        if( verbose ) write(*,*) 'testing pftcc_inplsrch :: minimize'
        call pftcc_inplsrch_init(pftcc, shlims)
        call pftcc_inplsrch_set_indices(1,1)
        call conduct_test
        if( .not. test_passed() ) stop '****inpl_srch_tester FAILURE pftcc_inplsrch :: minimize'

    contains

        function test_passed() result( passed )
            logical :: passed
            passed = .false.
            if( sherr_avg < SHERR_LIM .and. roerr_avg < ROERR_LIM ) passed = .true.
        end function test_passed

        subroutine conduct_test
            real    :: corrs(pftcc%get_nrots())
            integer :: loc(1)
            sherr_avg = 0.
            roerr_avg = 0.
            do iproj=1,NPROJS
                b%img_copy = imgs_refs(iproj)
                call insert_ref(b%img_copy)
                call pftcc%apply_ctf(p%smpd, tfun, DFX, DFY, ANGAST)
                b%img = imgs_ptcls(iproj)
                orientation = o_ptcls%get_ori(iproj) 
                call orientation%set('dfx',       DFX)
                call orientation%set('dfy',       DFY)
                call orientation%set('angast', ANGAST)
                call simimg(b%img, orientation, tfun, 'ctf', SNR, SNR_PINK, SNR_DETECTOR, BFAC)
                call prepimg4align(b, p, orientation)
                call b%proj%img2polarft(1, b%img, pftcc)
                corrs = pftcc%gencorrs(1,1)
                loc   = maxloc(corrs)
                crxy = pftcc_inplsrch_minimize(irot=loc(1))
                call orientation_opt%set('corr', crxy(1))
                call orientation_opt%e3set(-crxy(2))
                call orientation_opt%set('x', crxy(3))
                call orientation_opt%set('y', crxy(4))
                sherr     = euclid(crxy(3:4),[orientation%get('x'),orientation%get('y')])
                sherr_avg = sherr_avg + sherr
                roerr     = rad2deg(orientation.inplrotdist.orientation_opt)
                roerr_avg = roerr_avg + roerr
            end do
            sherr_avg = sherr_avg/real(NPROJS)
            roerr_avg = roerr_avg/real(NPROJS)
            if( verbose ) write(*,'(a,1x,f5.2)') 'SHIFT ERROR (IN PIXELS ): ', sherr_avg
            if( verbose ) write(*,'(a,1x,f5.2)') 'ROT   ERROR (IN DEGREES): ', roerr_avg
        end subroutine conduct_test

    end subroutine test_pftcc_inplsrch

    subroutine test_pftcc_shsrch
        integer :: iproj
        if( verbose ) write(*,*) 'testing pftcc_shsrch :: minimize'
        call pftcc_shsrch_init(pftcc, shlims)
        call conduct_test
        if( .not. test_passed() ) stop '****inpl_srch_tester FAILURE pftcc_shsrch :: minimize'

    contains

        function test_passed() result( passed )
            logical :: passed
            passed = .false.
            if( sherr_avg < SHSHERR_LIM ) passed = .true.
        end function test_passed

        subroutine conduct_test
            real    :: corrs(pftcc%get_nrots())
            integer :: loc(1)
            sherr_avg = 0.
            do iproj=1,NPROJS
                b%img_copy = imgs_refs(iproj)
                call insert_ref(b%img_copy)
                call pftcc%apply_ctf(p%smpd, tfun, DFX, DFY, ANGAST)
                b%img = imgs_ptcls(iproj)
                orientation = o_ptcls%get_ori(iproj) 
                call orientation%set('dfx',       DFX)
                call orientation%set('dfy',       DFY)
                call orientation%set('angast', ANGAST)
                call simimg(b%img, orientation, tfun, 'ctf', SNR, SNR_PINK, SNR_DETECTOR, BFAC)
                call prepimg4align(b, p, orientation)
                call b%proj%img2polarft(1, b%img, pftcc)
                corrs = pftcc%gencorrs(1,1)
                loc   = maxloc(corrs)
                call pftcc_shsrch_set_indices(1,1,loc(1))
                cxy = pftcc_shsrch_minimize()
                call orientation_opt%set('corr', cxy(1))
                call orientation_opt%set('x', cxy(2))
                call orientation_opt%set('y', cxy(3))
                sherr     = euclid(cxy(2:3),[orientation%get('x'),orientation%get('y')])
                sherr_avg = sherr_avg + sherr
            end do
            sherr_avg = sherr_avg/real(NPROJS)
            if( verbose ) write(*,'(a,1x,f5.2)') 'SHIFT ERROR (IN PIXELS ): ', sherr_avg
        end subroutine conduct_test

    end subroutine test_pftcc_shsrch

    subroutine insert_ref( img )
        class(image), intent(inout) :: img
        call img%norm
        ! apply mask
        call img%mask(p%msk, 'soft')
        ! move to Fourier space
        call img%fwd_ft
        ! transfer to polar coordinates
        call b%proj%img2polarft(1, img, pftcc, isptcl=.false.)
    end subroutine insert_ref

end module simple_inpl_srch_tester
