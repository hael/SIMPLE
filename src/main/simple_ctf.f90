! defines the Contrast Transfer Function (CTF) of the electron microscope
! based on a class used in CTFFIND4, developed by Alexis Rohou and Nikolaus
! Grigorieff at Janelia Farm. The below copyright statement therefore
! needs to be included here:
! Copyright 2014 Howard Hughes Medical Institute
! All rights reserved
! Use is subject to Janelia Farm Research Campus Software Copyright 1.1
! license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
module simple_ctf
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_image, only: image
implicit none

public :: ctf, ctf_test
private

type ctf
    private
    real    :: smpd        = 0.    !< sampling distance (input unit: A)
    real    :: kV          = 0.    !< acceleration voltage of the electron microscope (input unit: kV)
    real    :: Cs          = 0.    !< spherical aberration (input unit: mm)
    real    :: wl          = 0.    !< wavelength (input unit: A)
    real    :: amp_contr   = 0.07  !< fraction of amplitude contrast ([0.07,0.15] see Mindell 03)
    real    :: dfx         = 0.    !< underfocus x-axis, underfocus is positive; larger value = more underfocus (input unit: microns)
    real    :: dfy         = 0.    !< underfocus y-axis (input unit: microns)
    real    :: angast      = 0.    !< azimuth of x-axis 0.0 means axis is at 3 o'clock (input unit: degrees)
    real    :: phaseq      = 0.    !< phase constrast weight (derived constant)
    real    :: ampliq      = 0.    !< amplitude contrast weight (derived constant)
  contains
    procedure          :: init
    procedure, private :: eval_1
    procedure, private :: eval_2
    procedure, private :: eval_3
    procedure, private :: eval_4
    generic            :: eval => eval_1, eval_2, eval_3, eval_4
    procedure, private :: evalPhSh
    procedure, private :: eval_df
    procedure          :: ctf2img
    procedure          :: apply
    procedure          :: apply_serial
    procedure          :: apply_and_shift
    procedure          :: kV2wl
end type ctf

interface ctf
    module procedure constructor
end interface

contains

    function constructor( smpd, kV, Cs, amp_contr ) result( self )
        real, intent(in) :: smpd      !< sampling distance
        real, intent(in) :: kV        !< accelleration voltage
        real, intent(in) :: Cs        !< constant
        real, intent(in) :: amp_contr !< amplitude contrast
        type(ctf) :: self
        ! set constants
        self%kV        = kV
        self%wl        = self%kV2wl() / smpd
        self%Cs        = (Cs*1.0e7) / smpd
        self%amp_contr = amp_contr
        self%smpd      = smpd
        ! compute derived constants (phase and amplitude contrast weights)
        self%phaseq    = sqrt(1. - amp_contr * amp_contr)
        self%ampliq    = amp_contr
    end function constructor

    !>  \brief  initialise a CTF object with defocus/astigmatism params
    subroutine init( self, dfx, dfy, angast )
        class(ctf), intent(inout) :: self   !< instance
        real,       intent(in)    :: dfx    !< defocus x-axis (um)
        real,       intent(in)    :: dfy    !< defocus y-axis (um)
        real,       intent(in)    :: angast !< astigmatism (degrees)
        self%dfx    = (dfx*1.0e4)/self%smpd
        self%dfy    = (dfy*1.0e4)/self%smpd
        self%angast = deg2rad(angast)
    end subroutine init

    !>  \brief Returns the CTF, based on CTFFIND3 subroutine (see Mindell 2003)
    !
    ! How to use eval:
    !
    !            lFreqSq = tfun%getLowFreq4Fit**2
    !            hFreqSq = tfun%getHighFreq4Fit**2
    !            do k=-ydim,ydim
    !                kinv = inv_logical_kdim*real(k)
    !                kinvsq = kinv*kinv
    !                do h=-xdim,xdim
    !                   hinv = inv_logical_hdim*real(h)
    !                    hinvsq = hinv*hinv
    !                    spaFreqSq = kinvsq+hinvsq
    !                    if( spaFreqSq .gt. lFreqSq .and. spaFreqSq .le. hFreqSq )then
    !                        if( spaFreqSq .gt. 0. )then
    !                            ang = atan2(k,h)
    !                        else
    !                            ang = 0.
    !                        endif
    function eval_1( self, spaFreqSq, dfx, dfy, angast, ang ) result ( val )
        class(ctf), intent(inout) :: self      !< instance
        real,       intent(in)    :: spaFreqSq !< squared reciprocal pixels
        real,       intent(in)    :: dfx       !< Defocus along first axis (micrometers)
        real,       intent(in)    :: dfy       !< Defocus along second axis (for astigmatic CTF, dfx .ne. dfy) (micrometers)
        real,       intent(in)    :: angast    !< Azimuth of first axis. 0.0 means axis is at 3 o'clock. (radians)
        real,       intent(in)    :: ang       !< Angle at which to compute the CTF (radians)
        real :: val, phshift
        ! initialize the CTF object, using the input parameters
        call self%init(dfx, dfy, angast)
        ! compute phase shift
        phshift = self%evalPhSh(spaFreqSq, ang, 0.)
        ! compute value of CTF, assuming white particles
        val = (self%phaseq * sin(phshift) + self%ampliq * cos(phshift))
    end function eval_1

    function eval_2( self, spaFreqSq, dfx, dfy, angast, ang, add_phshift ) result ( val )
        class(ctf), intent(inout) :: self        !< instance
        real,       intent(in)    :: spaFreqSq   !< squared reciprocal pixels
        real,       intent(in)    :: dfx         !< Defocus along first axis (micrometers)
        real,       intent(in)    :: dfy         !< Defocus along second axis (for astigmatic CTF, dfx .ne. dfy) (micrometers)
        real,       intent(in)    :: angast      !< Azimuth of first axis. 0.0 means axis is at 3 o'clock. (radians)
        real,       intent(in)    :: ang         !< Angle at which to compute the CTF (radians)
        real,       intent(in)    :: add_phshift !< aditional phase shift (radians), for phase plate
        real :: val, phshift
        ! initialize the CTF object, using the input parameters
        call self%init(dfx, dfy, angast)
        ! compute phase shift
        phshift = self%evalPhSh(spaFreqSq, ang, add_phshift)
        ! compute value of CTF, assuming white particles
        val = (self%phaseq * sin(phshift) + self%ampliq * cos(phshift))
    end function eval_2

    !>  \brief Returns the CTF with pre-initialize parameters
    pure elemental real function eval_3( self, spaFreqSq, ang )
        class(ctf), intent(in) :: self        !< instance
        real,       intent(in) :: spaFreqSq   !< squared reciprocal pixels
        real,       intent(in) :: ang         !< Angle at which to compute the CTF (radians)
        real :: df, phshift
        ! compute DOUBLE the effective defocus
        df = (self%dfx + self%dfy + cos(2.0 * (ang - self%angast)) * (self%dfx - self%dfy))
        !! compute the ctf argument
        phshift = 0.5 * pi * self%wl * spaFreqSq * (df - self%wl * self%wl * spaFreqSq * self%Cs)
        ! compute value of CTF, assuming white particles
        eval_3 = (self%phaseq * sin(phshift) + self%ampliq * cos(phshift))
    end function eval_3

    !>  \brief Returns the CTF with pre-initialize parameters
    pure elemental real function eval_4( self, spaFreqSq, ang, add_phshift )
        class(ctf), intent(in) :: self        !< instance
        real,       intent(in) :: spaFreqSq   !< squared reciprocal pixels
        real,       intent(in) :: ang         !< Angle at which to compute the CTF (radians)
        real,       intent(in) :: add_phshift !< aditional phase shift (radians), for phase plate
        real :: df, phshift
        ! compute DOUBLE the effective defocus
        df      = (self%dfx + self%dfy + cos(2.0 * (ang - self%angast)) * (self%dfx - self%dfy))
        ! compute the ctf argument
        phshift = 0.5 * pi * self%wl * spaFreqSq * (df - self%wl * self%wl * spaFreqSq * self%Cs) + add_phshift
        ! compute value of CTF, assuming white particles
        eval_4  = (self%phaseq * sin(phshift) + self%ampliq * cos(phshift))
    end function eval_4

    !>  \brief returns the argument (radians) to the sine and cosine terms of the ctf
    !!
    !!  We follow the convention, like the rest of the cryo-EM/3DEM field, that underfocusing the objective lens
    !!  gives rise to a positive phase shift of scattered electrons, whereas the spherical aberration gives a
    !!  negative phase shift of scattered electrons
    pure elemental real function evalPhSh( self, spaFreqSq, ang, add_phshift ) result( phshift )
        class(ctf), intent(in) :: self        !< instance
        real,       intent(in) :: spaFreqSq   !< square of spatial frequency at which to compute the ctf (1/pixels^2)
        real,       intent(in) :: ang         !< angle at which to compute the ctf (radians)
        real,       intent(in) :: add_phshift !< aditional phase shift (radians), for phase plate
        real :: df                            !< defocus at point at which we're evaluating the ctf
        !! compute the defocus
        df = self%eval_df(ang)
        !! compute the ctf argument
        phshift = pi * self%wl * spaFreqSq * (df - 0.5 * self%wl * self%wl * spaFreqSq * self%Cs) + add_phshift
    end function evalPhSh

    !>  \brief  Return the effective defocus given the pre-set CTF parameters (from CTFFIND4)
    pure elemental real function eval_df( self, ang ) result (df)
        class(ctf), intent(in) :: self !< instance
        real,       intent(in) :: ang  !< angle at which to compute the defocus (radians)
        df = 0.5 * (self%dfx + self%dfy + cos(2.0 * (ang - self%angast)) * (self%dfx - self%dfy))
    end function eval_df

    !>  \brief  is for making a CTF image
    !!          modes: abs, ctf, flip, flipneg, neg, square
    subroutine ctf2img( self, img, dfx, mode, dfy, angast, bfac, add_phshift )
        class(ctf),       intent(inout) :: self        !< instance
        class(image),     intent(inout) :: img         !< image (output)
        real,             intent(in)    :: dfx         !< defocus x-axis
        character(len=*), intent(in)    :: mode        !< abs, ctf, flip, flipneg, neg, square
        real, optional,   intent(in)    :: dfy         !< defocus y-axis
        real, optional,   intent(in)    :: angast      !< angle of astigmatism
        real, optional,   intent(in)    :: bfac        !< bfactor
        real, optional,   intent(in)    :: add_phshift !< aditional phase shift (radians), for phase plate
        integer :: lims(3,2),h,k,phys(3),ldim(3)
        real    :: ang, tval, ddfy, aangast, spaFreqSq, hinv, aadd_phshift, kinv, inv_ldim(3)
        if( img%is_3d() )then
            print *, 'ldim: ', img%get_ldim()
            stop 'Only 4 2D images; ctf2img; simple_ctf'
        endif
        ! set CTF params
        ddfy = dfx
        if( present(dfy) ) ddfy = dfy
        aangast = 0.
        if( present(angast) ) aangast = angast
        aadd_phshift = 0.
        if( present(add_phshift) ) aadd_phshift = add_phshift
        ! init object
        call self%init(dfx, ddfy, aangast)
        ! initialize
        call img%set_ft(.true.)
        img      = 0.
        lims     = img%loop_lims(2)
        ldim     = img%get_ldim()
        inv_ldim = 1./real(ldim)
        !$omp parallel do collapse(2) default(shared) private(h,hinv,k,kinv,spaFreqSq,ang,tval,phys) &
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                hinv      = real(h) * inv_ldim(1)
                kinv      = real(k) * inv_ldim(2)
                spaFreqSq = hinv * hinv + kinv * kinv
                ang       = atan2(real(k),real(h))
                tval      = self%eval(spaFreqSq, ang, aadd_phshift)
                select case(mode)
                    case('abs')
                        tval = abs(tval)
                    case('ctf')
                        ! tval = tval
                    case('flip')
                        tval = sign(1.,tval)
                    case('flipneg')
                        tval = -sign(1.,tval)
                    case('neg')
                        tval = -tval
                    case('square')
                        tval = min(1.,max(tval * tval,0.001))
                    case('sqrt')
                        tval = min(1.,max(tval * tval,0.001))
                        tval = sqrt(tval * tval)
                    case DEFAULT
                        write(*,*) 'mode:', mode
                        stop 'unsupported in ctf2img; simple_ctf'
                end select
                phys = img%comp_addr_phys([h,k,0])
                call img%set_fcomp([h,k,0], phys, cmplx(tval,0.))
            end do
        end do
        !$omp end parallel do
        if( present(bfac) ) call img%apply_bfac(bfac)
    end subroutine ctf2img

    !>  \brief  is for applying CTF to an image
    subroutine apply( self, img, dfx, mode, dfy, angast, bfac, add_phshift )
        class(ctf),       intent(inout) :: self        !< instance
        class(image),     intent(inout) :: img         !< image (output)
        real,             intent(in)    :: dfx         !< defocus x-axis
        character(len=*), intent(in)    :: mode        !< abs, ctf, flip, flipneg, neg, square
        real, optional,   intent(in)    :: dfy         !< defocus y-axis
        real, optional,   intent(in)    :: angast      !< angle of astigmatism
        real, optional,   intent(in)    :: bfac        !< bfactor
        real, optional,   intent(in)    :: add_phshift !< aditional phase shift (radians), for phase plate
        integer     :: ldim(3), ldim_pd(3)
        real        :: aadd_phshift
        type(image) :: ctfimg, img_pd
        ldim = img%get_ldim()
        if( img%is_3d() )then
            print *, 'ldim: ', ldim
            stop 'Only 4 2D images; apply; simple_ctf'
        endif
        if( img%is_ft() )then
            call ctfimg%new(ldim, self%smpd)
            call self%ctf2img(ctfimg, dfx, mode, dfy, angast, bfac, add_phshift)
            call img%mul(ctfimg)
        else
            ldim_pd(1:2) = 2*ldim(1:2)
            ldim_pd(3)   = 1
            call ctfimg%new(ldim_pd, self%smpd)
            call img_pd%new(ldim_pd, self%smpd)
            call self%ctf2img(ctfimg, dfx, mode, dfy, angast, bfac, add_phshift)
            call img%pad_mirr(img_pd)
            call img_pd%fft()
            call img_pd%mul(ctfimg)
            call img_pd%ifft()
            call img_pd%clip(img)
            call img_pd%kill()
        endif
        call ctfimg%kill
    end subroutine apply

    !>  \brief  is for optimised serial application of CTF
    !!          modes: abs, ctf, flip, flipneg, neg, square
    subroutine apply_serial( self, img, dfx, mode, dfy, angast, add_phshift )
        class(ctf),       intent(inout) :: self        !< instance
        class(image),     intent(inout) :: img         !< image (output)
        real,             intent(in)    :: dfx         !< defocus x-axis
        character(len=*), intent(in)    :: mode        !< abs, ctf, flip, flipneg, neg, square
        real,             intent(in)    :: dfy         !< defocus y-axis
        real,             intent(in)    :: angast      !< angle of astigmatism
        real,             intent(in)    :: add_phshift !< aditional phase shift (radians), for phase plate
        integer :: lims(3,2),h,k,phys(3),ldim(3)
        real    :: ang,tval,spaFreqSq,hinv,kinv,inv_ldim(3)
        ! init object
        call self%init(dfx, dfy, angast)
        ! initialize
        lims     = img%loop_lims(2)
        ldim     = img%get_ldim()
        inv_ldim = 1./real(ldim)
        select case(mode)
            case('abs')
                do h=lims(1,1),lims(1,2)
                    do k=lims(2,1),lims(2,2)
                        hinv      = real(h) * inv_ldim(1)
                        kinv      = real(k) * inv_ldim(2)
                        spaFreqSq = hinv * hinv + kinv * kinv
                        ang       = atan2(real(k),real(h))
                        tval      = self%eval(spaFreqSq, ang, add_phshift)
                        phys      = img%comp_addr_phys([h,k,0])
                        call img%mul_cmat_at(phys(1),phys(2),phys(3), abs(tval))
                    end do
                end do
            case('ctf')
                do h=lims(1,1),lims(1,2)
                    do k=lims(2,1),lims(2,2)
                        hinv      = real(h) * inv_ldim(1)
                        kinv      = real(k) * inv_ldim(2)
                        spaFreqSq = hinv * hinv + kinv * kinv
                        ang       = atan2(real(k),real(h))
                        tval      = self%eval(spaFreqSq, ang, add_phshift)
                        phys      = img%comp_addr_phys([h,k,0])
                        call img%mul_cmat_at(phys(1),phys(2),phys(3), tval)
                    end do
                end do
            case('flip')
                do h=lims(1,1),lims(1,2)
                    do k=lims(2,1),lims(2,2)
                        hinv      = real(h) * inv_ldim(1)
                        kinv      = real(k) * inv_ldim(2)
                        spaFreqSq = hinv * hinv + kinv * kinv
                        ang       = atan2(real(k),real(h))
                        tval      = self%eval(spaFreqSq, ang, add_phshift)
                        phys      = img%comp_addr_phys([h,k,0])
                        call img%mul_cmat_at(phys(1),phys(2),phys(3), sign(1.,tval))
                    end do
                end do
            case('flipneg')
                do h=lims(1,1),lims(1,2)
                    do k=lims(2,1),lims(2,2)
                        hinv      = real(h) * inv_ldim(1)
                        kinv      = real(k) * inv_ldim(2)
                        spaFreqSq = hinv * hinv + kinv * kinv
                        ang       = atan2(real(k),real(h))
                        tval      = self%eval(spaFreqSq, ang, add_phshift)
                        phys      = img%comp_addr_phys([h,k,0])
                        call img%mul_cmat_at(phys(1),phys(2),phys(3), -sign(1.,tval))
                    end do
                end do
            case('neg')
                do h=lims(1,1),lims(1,2)
                    do k=lims(2,1),lims(2,2)
                        hinv      = real(h) * inv_ldim(1)
                        kinv      = real(k) * inv_ldim(2)
                        spaFreqSq = hinv * hinv + kinv * kinv
                        ang       = atan2(real(k),real(h))
                        tval      = self%eval(spaFreqSq, ang, add_phshift)
                        phys      = img%comp_addr_phys([h,k,0])
                        call img%mul_cmat_at( phys(1),phys(2),phys(3), -tval)
                    end do
                end do
            case('square')
                do h=lims(1,1),lims(1,2)
                    do k=lims(2,1),lims(2,2)
                        hinv      = real(h) * inv_ldim(1)
                        kinv      = real(k) * inv_ldim(2)
                        spaFreqSq = hinv * hinv + kinv * kinv
                        ang       = atan2(real(k),real(h))
                        tval      = self%eval(spaFreqSq, ang, add_phshift)
                        phys      = img%comp_addr_phys([h,k,0])
                        call img%mul_cmat_at( phys(1),phys(2),phys(3), min(1.,max(tval**2.,0.001)))
                    end do
                end do
            case DEFAULT
                write(*,*) 'mode:', mode
                stop 'unsupported in ctf2img; simple_ctf'
        end select
    end subroutine apply_serial

    !>  \brief  is for applying CTF to an image and shifting it (used in classaverager)
    !!          KEEP THIS ROUTINE SERIAL
    subroutine apply_and_shift( self, img, imode, lims, rho, x, y, dfx, dfy, angast, add_phshift)
        class(ctf),     intent(inout) :: self        !< instance
        class(image),   intent(inout) :: img         !< modified image (output)
        integer,        intent(in)    :: imode       !< 1=abs 2=ctf 3=no
        integer,        intent(in)    :: lims(3,2)   !< loop limits
        real,           intent(out)   :: rho(lims(1,1):lims(1,2),lims(2,1):lims(2,2))
        real,           intent(in)    :: x, y        !< rotational origin shift
        real,           intent(in)    :: dfx         !< defocus x-axis
        real,           intent(in)    :: dfy         !< defocus y-axis
        real,           intent(in)    :: angast      !< angle of astigmatism
        real,           intent(in)    :: add_phshift !< aditional phase shift (radians), for phase plate
        integer :: ldim(3),logi(3),h,k,phys(3)
        real    :: ang,tval,spaFreqSq,hinv,kinv,inv_ldim(3)
        complex :: comp
        ! initialize
        call self%init(dfx, dfy, angast)
        ldim     = img%get_ldim()
        inv_ldim = 1./real(ldim)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                ! calculate CTF and CTF**2.0 values
                hinv      = real(h) * inv_ldim(1)
                kinv      = real(k) * inv_ldim(2)
                spaFreqSq = hinv * hinv + kinv * kinv
                ang       = atan2(real(k),real(h))
                tval      = 1.0
                if( imode <  3 ) tval = self%eval(spaFreqSq, ang, add_phshift)
                ! rho
                rho(h,k)  = tval * tval
                if( imode == 1 ) tval = abs(tval)
                ! multiply image with tval & weight
                logi = [h,k,0]
                phys = img%comp_addr_phys(logi)
                call img%mul_cmat_at(phys, tval)
                ! shift image
                call img%mul_cmat_at(phys(1),phys(2),phys(3), img%oshift(logi,[x,y,0.]))
            end do
        end do
    end subroutine apply_and_shift

    pure elemental real function kV2wl( self ) result (wavelength)
        class(ctf), intent(in) :: self
        wavelength = 12.26 / sqrt((1000.0 * self%kV) + (0.9784*(1000.0 * self%kV)**2)/(10.0**6.0))
    end function kV2wl

    !>  \brief This test compare ctf2img implementations
    subroutine ctf_test

        type(ctf)               :: tfun
        type(image)             :: img1, img2
        real                    :: dfx_ran, dfy_ran, angast_ran, phshift_ran, err, minmax(2)
        integer, parameter      :: BOX=512, NTST=10
        real,    parameter      :: SMPD=1.26, KV=300., CS=2.7, AC=0.1
        real,       allocatable :: tvals_old(:,:), tvals_new(:,:)
        integer                 :: itst
        integer(timer_int_kind) :: tctf
         ! test:            1 err    :    0.00000000
         ! test:            2 err    :    0.00000000
         ! test:            3 err    :    0.00000000
         ! test:            4 err    :    0.00000000
         ! test:            5 err    :    0.00000000
         ! test:            6 err    :    0.00000000
         ! test:            7 err    :    0.00000000
         ! test:            8 err    :    0.00000000
         ! test:            9 err    :    0.00000000
         ! test:           10 err    :    0.00000000
         ! err    :    0.00000000
         ! time(s):   0.33581284400000000
        tctf   = tic()
        err = 0.
        do itst=1,NTST
            dfx_ran     = 0.5 + ran3() * 4.5
            dfy_ran     = 0.5 + ran3() * 4.5
            angast_ran  = ran3() * 2. * pi
            phshift_ran = ran3()
            ! previous implementation
            call img1%new([BOX,BOX,1], SMPD)
            call img1%set_ft(.true.)
            call ctf2img_old(img1, dfx_ran, 'ctf', dfy_ran, angast_ran, add_phshift=phshift_ran)
            ! current
            call img2%new([BOX,BOX,1], SMPD)
            call img2%set_ft(.true.)
            tfun = ctf(SMPD, KV, CS, AC)
            call tfun%ctf2img(img2, dfx_ran, 'ctf', dfy_ran, angast_ran, add_phshift=phshift_ran)
            !
            img1   = img1.lone.img2
            minmax = img1%minmax()
            err = err + minmax(2)
            print *, 'test: ', itst, '; err    : '  , minmax(2)
        enddo
        print *, 'err    : '  , err
        print *, 'time(s): ', toc(tctf)


        contains

            !>  \brief  is for making a CTF image, previous implementation
            subroutine ctf2img_old(img, dfx, mode, dfy, angast, bfac, add_phshift )
                class(image),     intent(inout) :: img         !< image (output)
                real,             intent(in)    :: dfx         !< defocus x-axis
                character(len=*), intent(in)    :: mode        !< abs, ctf, flip, flipneg, neg, square
                real, optional,   intent(in)    :: dfy         !< defocus y-axis
                real, optional,   intent(in)    :: angast      !< angle of astigmatism
                real, optional,   intent(in)    :: bfac        !< bfactor
                real, optional,   intent(in)    :: add_phshift !< aditional phase shift (radians), for phase plate
                type(ctf) :: tfun
                integer   :: lims(3,2),h,k,phys(3),ldim(3)
                real      :: ang, tval, ddfy, aangast, spaFreqSq, hinv, aadd_phshift, kinv, inv_ldim(3)
                if( img%is_3d() )then
                    print *, 'ldim: ', img%get_ldim()
                    stop 'Only 4 2D images; ctf2img; simple_ctf'
                endif
                ! set CTF params
                ddfy = dfx
                if( present(dfy) ) ddfy = dfy
                aangast = 0.
                if( present(angast) ) aangast = angast
                aadd_phshift = 0.
                if( present(add_phshift) ) aadd_phshift = add_phshift
                ! init object
                tfun = ctf(SMPD, KV, CS, AC)
                call tfun%init(dfx, ddfy, aangast)
                ! initialize
                img      = cmplx(0.,0.)
                lims     = img%loop_lims(2)
                ldim     = img%get_ldim()
                inv_ldim = 1./real(ldim)
                !$omp parallel do collapse(2) default(shared) private(h,hinv,k,kinv,spaFreqSq,ang,tval,phys) &
                !$omp schedule(static) proc_bind(close)
                do h=lims(1,1),lims(1,2)
                    do k=lims(2,1),lims(2,2)
                        hinv      = real(h) * inv_ldim(1)
                        kinv      = real(k) * inv_ldim(2)
                        spaFreqSq = hinv * hinv + kinv * kinv
                        ang       = atan2(real(k),real(h))
                        tval      = tfun%eval(spaFreqSq, dfx, ddfy, aangast, ang, aadd_phshift)
                        select case(mode)
                            case('abs')
                                tval = abs(tval)
                            case('ctf')
                                ! tval = tval
                            case('flip')
                                tval = sign(1.,tval)
                            case('flipneg')
                                tval = -sign(1.,tval)
                            case('neg')
                                tval = -tval
                            case('square')
                                tval = min(1.,max(tval**2.,0.001))
                            case DEFAULT
                                write(*,*) 'mode:', mode
                                stop 'unsupported in ctf2img; simple_ctf'
                        end select
                        phys = img%comp_addr_phys([h,k,0])
                        call img%set_fcomp([h,k,0], phys, cmplx(tval,0.))
                    end do
                end do
                !$omp end parallel do
                if( present(bfac) ) call img%apply_bfac(bfac)
            end subroutine ctf2img_old

    end subroutine ctf_test

end module simple_ctf
