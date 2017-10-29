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
#include "simple_lib.f08"
implicit none

public :: ctf
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
    procedure, private :: init
    procedure, private :: eval_1
    procedure, private :: eval_2
    generic            :: eval => eval_1, eval_2
    procedure, private :: evalPhSh
    procedure, private :: eval_df
    procedure          :: ctf2img
    procedure          :: ctf2pspecimg
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
        use simple_math, only: deg2rad
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
        use simple_image, only: image
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
                tval      = self%eval(spaFreqSq, dfx, ddfy, aangast, ang, aadd_phshift)
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
    end subroutine ctf2img

    !>  \brief  is for making a CTF power-spec image
    subroutine ctf2pspecimg( self, img, dfx, dfy, angast, add_phshift )
        use simple_image, only: image
        class(ctf),     intent(inout) :: self        !< instance
        class(image),   intent(inout) :: img         !< image (output)
        real,           intent(in)    :: dfx         !< defocus x-axis
        real,           intent(in)    :: dfy         !< defocus y-axis
        real,           intent(in)    :: angast      !< angle of astigmatism
        real, optional, intent(in)    :: add_phshift !< aditional phase shift (radians), for phase plate
        integer :: lims(3,2),h,mh,k,mk,phys(3),ldim(3),inds(3)
        real    :: ang, tval, spaFreqSq, hinv, aadd_phshift, kinv, inv_ldim(3)
        ! initialize
        aadd_phshift = 0.
        if( present(add_phshift) ) aadd_phshift = add_phshift
        call self%init(dfx, dfy, angast)
        img      = 0.
        lims     = img%loop_lims(3)
        mh       = maxval(lims(1,:))
        mk       = maxval(lims(2,:))
        inds     = 1
        ldim     = img%get_ldim()
        inv_ldim = 1./real(ldim)
        !$omp parallel do collapse(2) default(shared) private(h,hinv,k,kinv,inds,spaFreqSq,ang,tval,phys) &
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                inds(1)   = min(max(1,h+mh+1),ldim(1))
                inds(2)   = min(max(1,k+mk+1),ldim(2))
                inds(3)   = 1
                hinv      = real(h) * inv_ldim(1)
                kinv      = real(k) * inv_ldim(2)
                spaFreqSq = hinv * hinv + kinv * kinv
                ang       = atan2(real(k),real(h))
                tval      = self%eval(spaFreqSq, dfx, dfy, angast, ang, aadd_phshift)
                tval      = sqrt(min(1.,max(tval**2.,0.001)))
                call img%set(inds, tval)
            end do
        end do
        !$omp end parallel do
    end subroutine ctf2pspecimg

    !>  \brief  is for applying CTF to an image
    subroutine apply( self, img, dfx, mode, dfy, angast, bfac, add_phshift )
        use simple_image, only: image
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
            call img_pd%fwd_ft
            call img_pd%mul(ctfimg)
            call img_pd%bwd_ft
            call img_pd%clip(img)
            call img_pd%kill
        endif
        call ctfimg%kill
    end subroutine apply

    !>  \brief  is for optimised serial application of CTF
    !!          modes: abs, ctf, flip, flipneg, neg, square
    subroutine apply_serial( self, img, dfx, mode, dfy, angast, add_phshift )
        use simple_image, only: image
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
                        tval      = self%eval(spaFreqSq, dfx, dfy, angast, ang, add_phshift)
                        phys = img%comp_addr_phys([h,k,0])
                        call img%mul_cmat_at(abs(tval), phys)
                    end do
                end do
            case('ctf')
                do h=lims(1,1),lims(1,2)
                    do k=lims(2,1),lims(2,2)
                        hinv      = real(h) * inv_ldim(1)
                        kinv      = real(k) * inv_ldim(2)
                        spaFreqSq = hinv * hinv + kinv * kinv
                        ang       = atan2(real(k),real(h))
                        tval      = self%eval(spaFreqSq, dfx, dfy, angast, ang, add_phshift)
                        phys = img%comp_addr_phys([h,k,0])
                        call img%mul_cmat_at(tval, phys)
                    end do
                end do
            case('flip')
                do h=lims(1,1),lims(1,2)
                    do k=lims(2,1),lims(2,2)
                        hinv      = real(h) * inv_ldim(1)
                        kinv      = real(k) * inv_ldim(2)
                        spaFreqSq = hinv * hinv + kinv * kinv
                        ang       = atan2(real(k),real(h))
                        tval      = self%eval(spaFreqSq, dfx, dfy, angast, ang, add_phshift)
                        phys = img%comp_addr_phys([h,k,0])
                        call img%mul_cmat_at(sign(1.,tval), phys)
                    end do
                end do
            case('flipneg')
                do h=lims(1,1),lims(1,2)
                    do k=lims(2,1),lims(2,2)
                        hinv      = real(h) * inv_ldim(1)
                        kinv      = real(k) * inv_ldim(2)
                        spaFreqSq = hinv * hinv + kinv * kinv
                        ang       = atan2(real(k),real(h))
                        tval      = self%eval(spaFreqSq, dfx, dfy, angast, ang, add_phshift)
                        phys = img%comp_addr_phys([h,k,0])
                        call img%mul_cmat_at(-sign(1.,tval), phys)
                    end do
                end do
            case('neg')
                do h=lims(1,1),lims(1,2)
                    do k=lims(2,1),lims(2,2)
                        hinv      = real(h) * inv_ldim(1)
                        kinv      = real(k) * inv_ldim(2)
                        spaFreqSq = hinv * hinv + kinv * kinv
                        ang       = atan2(real(k),real(h))
                        tval      = self%eval(spaFreqSq, dfx, dfy, angast, ang, add_phshift)
                        phys = img%comp_addr_phys([h,k,0])
                        call img%mul_cmat_at(-tval, phys)
                    end do
                end do
            case('square')
                do h=lims(1,1),lims(1,2)
                    do k=lims(2,1),lims(2,2)
                        hinv      = real(h) * inv_ldim(1)
                        kinv      = real(k) * inv_ldim(2)
                        spaFreqSq = hinv * hinv + kinv * kinv
                        ang       = atan2(real(k),real(h))
                        tval      = self%eval(spaFreqSq, dfx, dfy, angast, ang, add_phshift)
                        phys = img%comp_addr_phys([h,k,0])
                        call img%mul_cmat_at(min(1.,max(tval**2.,0.001)), phys)
                    end do
                end do
            case DEFAULT
                write(*,*) 'mode:', mode
                stop 'unsupported in ctf2img; simple_ctf'
        end select
    end subroutine apply_serial

    !>  \brief  is for applying CTF to an image and shifting it
    subroutine apply_and_shift( self, img, imgctfsq, x, y, dfx, mode, dfy, angast, add_phshift )
        use simple_image, only: image
        class(ctf),       intent(inout) :: self        !< instance
        class(image),     intent(inout) :: img         !< modified image (output)
        class(image),     intent(inout) :: imgctfsq    !< CTF**2 image (output)
        real,             intent(in)    :: x, y        !< rotational origin shift
        real,             intent(in)    :: dfx         !< defocus x-axis
        character(len=*), intent(in)    :: mode        !< abs, ctf, flip, flipneg, neg, square
        real, optional,   intent(in)    :: dfy         !< defocus y-axis
        real, optional,   intent(in)    :: angast      !< angle of astigmatism
        real, optional,   intent(in)    :: add_phshift !< aditional phase shift (radians), for phase plate
        integer :: ldim(3),logi(3),imode,lims(3,2),h,k,phys(3)
        real    :: ang,tval,ddfy,aangast,spaFreqSq,hinv,kinv,inv_ldim(3),tvalsq,aadd_phshift
        complex :: comp
        ldim = img%get_ldim()
        ! check that image is 2D
        if( img%is_3d() )then
            print *, 'ldim: ', ldim
            stop 'Only 4 2D images; apply; simple_ctf'
        endif
        ! check that image is FTed
        if( .not. img%is_ft() )then
            stop 'expecting FT input; simple_ctf :: apply_and_shift'
        endif
        ! set CTF params
        ddfy = dfx
        if( present(dfy) ) ddfy = dfy
        aangast = 0.
        if( present(angast) ) aangast = angast
        aadd_phshift = 0.
        if( present(add_phshift) ) aadd_phshift = add_phshift
        select case(mode)
            case('abs')
                imode = 1
            case('ctf')
                imode = 2
            case DEFAULT
                imode = 3
        end select
        ! initialize
        call self%init(dfx, ddfy, aangast)
        lims     = img%loop_lims(2)
        ldim     = img%get_ldim()
        inv_ldim = 1./real(ldim)
        !$omp parallel do collapse(2) default(shared) proc_bind(close) &
        !$omp private(h,hinv,k,kinv,spaFreqSq,ang,tval,tvalsq,logi,phys,comp) schedule(static)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                ! calculate CTF and CTF**2.0 values
                hinv      = real(h) * inv_ldim(1)
                kinv      = real(k) * inv_ldim(2)
                spaFreqSq = hinv * hinv + kinv * kinv
                ang       = atan2(real(k),real(h))
                tval      = 1.0
                if( imode <  3 ) tval = self%eval(spaFreqSq, dfx, ddfy, aangast, ang, aadd_phshift)
                tvalsq = tval * tval
                if( imode == 1 ) tval = abs(tval)
                ! multiply image
                logi = [h,k,0]
                phys = img%comp_addr_phys(logi)
                comp = img%get_fcomp(logi, phys)
                comp = comp * tval
                ! shift image
                comp = comp * img%oshift(logi, [x,y,0.])
                ! set outputs
                call img%set_fcomp(logi, phys, comp)
                call imgctfsq%set_fcomp(logi, phys, cmplx(tvalsq,0.))
            end do
        end do
        !$omp end parallel do
    end subroutine apply_and_shift

    pure elemental real function kV2wl( self ) result (wavelength)
        class(ctf), intent(in) :: self
        wavelength = 12.26 / sqrt((1000.0 * self%kV) + (0.9784*(1000.0 * self%kV)**2)/(10.0**6.0))
    end function kV2wl

end module simple_ctf
