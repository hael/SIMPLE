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
implicit none
private
public :: ctf
#include "simple_local_flags.inc"

real :: CTF_FIRST_LIM = PI !< Phase shift defining limit for intact CTF: CTF_FIRST_LIM==pi/2=>peak; CTF_FIRST_LIM==pi=>zero)

type ctf
    private
    real    :: smpd            = 0.    !< sampling distance (input unit: A)
    real    :: kV              = 0.    !< acceleration voltage of the electron microscope (input unit: kV)
    real    :: Cs              = 0.    !< spherical aberration (input unit: mm)
    real    :: wl              = 0.    !< wavelength (input unit: A)
    real    :: amp_contr       = 0.07  !< fraction of amplitude contrast ([0.07,0.15] see Mindell 03)
    real    :: amp_contr_const = 0.    !< Amplitude contrast derived term
    real    :: dfx             = 0.    !< underfocus x-axis, underfocus is positive; larger value = more underfocus (input unit: microns)
    real    :: dfy             = 0.    !< underfocus y-axis (input unit: microns)
    real    :: angast          = 0.    !< azimuth of x-axis 0.0 means axis is at 3 o'clock (input unit: degrees)
    contains
    procedure          :: init
    procedure, private :: evalPhSh
    procedure, private :: eval_1, eval_2, eval_3, eval_4, eval_5
    generic            :: eval => eval_1, eval_2, eval_3, eval_4, eval_5
    procedure, private :: eval_sign
    procedure, private :: eval_df
    procedure          :: nextrema
    procedure          :: spafreqsqatnthzero
    procedure          :: apply
    procedure          :: ctf2img
    procedure          :: ctf_1stzero2img
    procedure          :: apply_serial
    procedure          :: wienerlike_restoration
    procedure          :: phaseflip_and_shift_serial
    procedure          :: eval_and_apply, eval_and_apply_before_first_zero
    procedure, private :: kV2wl
    procedure          :: apply_convention
    procedure          :: calc_ice_frac
end type ctf

interface ctf
    module procedure constructor
end interface

contains

    elemental function constructor( smpd, kV, Cs, amp_contr ) result( self )
        real, intent(in) :: smpd      !< sampling distance
        real, intent(in) :: kV        !< accelleration voltage
        real, intent(in) :: Cs        !< constant
        real, intent(in) :: amp_contr !< amplitude contrast
        type(ctf) :: self
        real      :: phaseq
        ! set constants
        self%kV        = kV
        self%wl        = self%kV2wl() / smpd
        self%Cs        = (Cs*1.0e7) / smpd
        self%amp_contr = amp_contr
        self%smpd      = smpd
        ! compute derived constant
        phaseq               = sqrt(1. - self%amp_contr*self%amp_contr)
        self%amp_contr_const = atan(self%amp_contr / phaseq)
    end function constructor

    !>  \brief  initialise a CTF object with defocus/astigmatism params
    elemental subroutine init( self, dfx, dfy, angast )
        class(ctf), intent(inout) :: self   !< instance
        real,       intent(in)    :: dfx    !< defocus x-axis (um)
        real,       intent(in)    :: dfy    !< defocus y-axis (um)
        real,       intent(in)    :: angast !< astigmatism (degrees)
        self%dfx    = (dfx*1.0e4)/self%smpd
        self%dfy    = (dfy*1.0e4)/self%smpd
        self%angast = deg2rad(angast)
    end subroutine init

    !>  \brief returns the argument (radians) to the ctf
    !!
    !!  We follow the convention, like the rest of the cryo-EM/3DEM field, that underfocusing the objective lens
    !!  gives rise to a positive phase shift of scattered electrons, whereas the spherical aberration gives a
    !!  negative phase shift of scattered electrons
    pure elemental real function evalPhSh( self, spaFreqSq, ang, add_phshift )
        class(ctf), intent(in) :: self        !< instance
        real,       intent(in) :: spaFreqSq   !< square of spatial frequency at which to compute the ctf (1/pixels^2)
        real,       intent(in) :: ang         !< angle at which to compute the ctf (radians)
        real,       intent(in) :: add_phshift !< aditional phase shift (radians), for phase plate
        real :: df                            !< defocus at point at which we're evaluating the ctf
        !! compute the defocus
        df = self%eval_df(ang)
        !! compute the ctf argument
        evalPhSh = PI * self%wl * spaFreqSq * (df - 0.5 * self%wl*self%wl * spaFreqSq * self%Cs) + add_phshift
    end function evalPhSh

    !>  \brief Returns the CTF, based on CTFFIND4 subroutine (Rohou & Grigorieff (2015))
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
    real function eval_1( self, spaFreqSq, dfx, dfy, angast, ang )
        class(ctf), intent(inout) :: self      !< instance
        real,       intent(in)    :: spaFreqSq !< squared reciprocal pixels
        real,       intent(in)    :: dfx       !< Defocus along first axis (micrometers)
        real,       intent(in)    :: dfy       !< Defocus along second axis (for astigmatic CTF, dfx .ne. dfy) (micrometers)
        real,       intent(in)    :: angast    !< Azimuth of first axis. 0.0 means axis is at 3 o'clock. (radians)
        real,       intent(in)    :: ang       !< Angle at which to compute the CTF (radians)
        ! initialize the CTF object, using the input parameters
        call self%init(dfx, dfy, angast)
        ! compute phase shift + amplitude constrast term & compute value of CTF, assuming white particles
        eval_1 = sin( self%evalPhSh(spaFreqSq, ang, 0.) + self%amp_contr_const )
    end function eval_1

    real function eval_2( self, spaFreqSq, dfx, dfy, angast, ang, add_phshift )
        class(ctf), intent(inout) :: self        !< instance
        real,       intent(in)    :: spaFreqSq   !< squared reciprocal pixels
        real,       intent(in)    :: dfx         !< Defocus along first axis (micrometers)
        real,       intent(in)    :: dfy         !< Defocus along second axis (for astigmatic CTF, dfx .ne. dfy) (micrometers)
        real,       intent(in)    :: angast      !< Azimuth of first axis. 0.0 means axis is at 3 o'clock. (radians)
        real,       intent(in)    :: ang         !< Angle at which to compute the CTF (radians)
        real,       intent(in)    :: add_phshift !< aditional phase shift (radians), for phase plate
        ! initialize the CTF object, using the input parameters
        call self%init(dfx, dfy, angast)
        ! compute phase shift + amplitude constrast term & compute value of CTF, assuming white particles
        eval_2 = sin( self%evalPhSh(spaFreqSq, ang, add_phshift) + self%amp_contr_const )
    end function eval_2

    !>  \brief Returns the CTF with pre-initialize parameters
    pure elemental real function eval_3( self, spaFreqSq, ang )
        class(ctf), intent(in) :: self        !< instance
        real,       intent(in) :: spaFreqSq   !< squared reciprocal pixels
        real,       intent(in) :: ang         !< Angle at which to compute the CTF (radians)
        ! compute phase shift + amplitude constrast term & compute value of CTF, assuming white particles
        eval_3 = sin( self%evalPhSh(spaFreqSq, ang, 0.) + self%amp_contr_const )
    end function eval_3

    !>  \brief Returns the CTF with pre-initialize parameters
    pure elemental real function eval_4( self, spaFreqSq, ang, add_phshift )
        class(ctf), intent(in) :: self        !< instance
        real,       intent(in) :: spaFreqSq   !< squared reciprocal pixels
        real,       intent(in) :: ang         !< Angle at which to compute the CTF (radians)
        real,       intent(in) :: add_phshift !< aditional phase shift (radians), for phase plate
        ! compute phase shift + amplitude constrast term & compute value of CTF, assuming white particles
        eval_4 = sin( self%evalPhSh(spaFreqSq, ang, add_phshift) + self%amp_contr_const )
    end function eval_4

    !>  \brief Returns the CTF with pre-initialize parameters
    pure elemental real function eval_5( self, spaFreqSq, ang, add_phshift, before1stzero )
        class(ctf), intent(in) :: self        !< instance
        real,       intent(in) :: spaFreqSq   !< squared reciprocal pixels
        real,       intent(in) :: ang         !< Angle at which to compute the CTF (radians)
        real,       intent(in) :: add_phshift !< aditional phase shift (radians), for phase plate
        logical,    intent(in) :: before1stzero !< whether of not the CTF value is calculated before the first zero
        real :: totalphaseshift
        ! compute phase shift + amplitude constrast term & compute value of CTF, assuming white particles
        ! totalphaseshift is > 0 and not modulated by 2pi!
        totalphaseshift = self%evalPhSh(spaFreqSq, ang, add_phshift) + self%amp_contr_const
        if( before1stzero )then
            eval_5 = sin( totalphaseshift )
        else
            if( totalphaseshift > CTF_FIRST_LIM )then
                eval_5 = sin( totalphaseshift )
            else
                eval_5 = 1.0
            endif
        endif
    end function eval_5

    !>  \brief Returns the sign of the CTF with pre-initialize parameters
    pure elemental integer function eval_sign( self, spaFreqSq, ang, add_phshift )
        class(ctf), intent(in) :: self        !< instance
        real,       intent(in) :: spaFreqSq   !< squared reciprocal pixels
        real,       intent(in) :: ang         !< Angle at which to compute the CTF (radians)
        real,       intent(in) :: add_phshift !< aditional phase shift (radians), for phase plate
        real :: angle
        ! compute phase shift + amplitude constrast term, no need to evaluate the sine to workout the sign
        angle = self%evalPhSh(spaFreqSq, ang, add_phshift) + self%amp_contr_const
        do while( angle > TWOPI )
            angle = angle - TWOPI
        enddo
        do while( angle < 0. )
            angle = angle + TWOPI
        enddo
        eval_sign = 1
        if( angle > PI ) eval_sign = -1
    end function eval_sign

    !>  \brief  Return the effective defocus given the pre-set CTF parameters (from CTFFIND4)
    pure elemental real function eval_df( self, ang )
        class(ctf), intent(in) :: self !< instance
        real,       intent(in) :: ang  !< angle at which to compute the defocus (radians)
        eval_df = 0.5 * (self%dfx + self%dfy + cos(2.0 * (ang - self%angast)) * (self%dfx - self%dfy))
    end function eval_df

    !>  \brief  Edited Eq 11 of Rohou & Grigorieff (2015)
    integer function nextrema( self, spaFreqSq, ang, phshift )
        class(ctf), intent(in) :: self
        real,       intent(in) :: spaFreqSq   !< squared reciprocal pixels
        real,       intent(in) :: ang         !< Angle at which to compute the CTF (radians) )
        real,       intent(in) :: phshift     !< aditional phase shift (radians), for phase plate
        real :: nextrema_before_chi_extremum, chi_extremem_spafreqsq, df, rn
        ! spatial frequency of the phase aberration extremum
        if( self%cs < 0.0001 )then
            chi_extremem_spafreqsq = 9999.9999;
        else
            df = self%eval_df(ang)
            chi_extremem_spafreqsq = df / (self%wl * self%wl * self%cs)
        endif
        ! extrema
        if( spaFreqSq <= chi_extremem_spafreqsq )then
            rn = abs(floor( 1. / PI * (self%evalPhSh(spaFreqSq, ang, phshift)+self%amp_contr_const) + 0.5))
        else
            nextrema_before_chi_extremum = floor( 1. / PI * (self%evalPhSh(chi_extremem_spafreqsq, ang, phshift)+self%amp_contr_const) + 0.5)
            rn = floor( 1. / PI * (self%evalPhSh(spaFreqSq, ang, phshift)+self%amp_contr_const) + 0.5)
            rn = nextrema_before_chi_extremum + abs(rn-nextrema_before_chi_extremum)
        endif
        nextrema = nint(abs(rn))
    end function nextrema

    !>  \brief  returns spatial frequency at some zero
    real function SpaFreqSqAtNthZero( self, nzero, add_phshift, ang )
        class(ctf), intent(in)  :: self
        integer,    intent(in)  :: nzero
        real,       intent(in)  :: add_phshift
        real,       intent(in)  :: ang
        real :: phshift, A, B, C, determinant, one,two
        phshift = real(nzero) * PI
        A = -0.5 * PI * self%wl**3. * self%cs
        B = PI * self%wl * self%eval_df(ang)
        C = self%amp_contr_const + add_phshift
        determinant = B**2. - 4.*A*(C-phshift)
        if( abs(self%cs) < 0.0001 )then
            SpaFreqSqAtNthZero = (phshift-C) / B
        else
            if( determinant < 0. )then
                SpaFreqSqAtNthZero = 0.
            else
                one = (-B + sqrt(determinant)) / (2.0 * A)
                two = (-B - sqrt(determinant)) / (2.0 * A)
                if( one > 0. .and. two > 0. )then
                    SpaFreqSqAtNthZero = min(one,two)
                else if( one > 0. )then
                    SpaFreqSqAtNthZero = one
                else if( two > 0. )then
                    SpaFreqSqAtNthZero = two
                else
                    SpaFreqSqAtNthZero = 0.
                endif
            endif
        endif
    end function SpaFreqSqAtNthZero

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
        integer         :: ldim(3), ldim_pd(3)
        type(ctfparams) :: ctfvars
        type(image)     :: img_pd
        ldim = img%get_ldim()
        if( img%is_3d() )then
            write(logfhandle,*) 'ldim: ', ldim
            THROW_HARD('Only 4 2D images; apply')
        endif
        ctfvars%smpd = self%smpd
        ctfvars%dfx  = dfx
        ! defaults
        ctfvars%dfy     = ctfvars%dfx
        ctfvars%angast  = 0.
        ctfvars%phshift = 0.
        ! optionals
        if(present(dfy))         ctfvars%dfy     = dfy
        if(present(angast))      ctfvars%angast  = angast
        if(present(add_phshift)) ctfvars%phshift = add_phshift
        if( img%is_ft() )then
            call self%apply_serial(img, mode, ctfvars)
        else
            ldim_pd(1:2) = 2*ldim(1:2)
            ldim_pd(3)   = 1
            call img_pd%new(ldim_pd, self%smpd)
            call img%pad_mirr(img_pd)
            call img_pd%fft()
            call self%apply_serial(img_pd, mode, ctfvars)
            call img_pd%ifft()
            call img_pd%clip(img)
            call img_pd%kill()
        endif
        if( present(bfac) ) call img%apply_bfac(bfac)
    end subroutine apply

    !>  \brief  is for generating an image of CTF
    subroutine ctf2img( self, img, dfx, dfy, angast, phshift )
        use simple_image, only: image
        class(ctf),       intent(inout) :: self
        class(image),     intent(inout) :: img
        real,             intent(in)    :: dfx, dfy, angast
        real,   optional, intent(in)    :: phshift
        integer :: lims(3,2),h,k,phys(3),ldim(3)
        real    :: ang,tval,spaFreqSq,hinv,kinv,inv_ldim(3),pphshift
        pphshift = 0.
        if( present(phshift) ) pphshift = phshift
        ! flag image as FT
        call img%set_ft(.true.)
        ! init object
        call self%init(dfx, dfy, angast)
        ! initialize
        lims     = img%loop_lims(2)
        ldim     = img%get_ldim()
        inv_ldim = 1./real(ldim)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                hinv      = real(h) * inv_ldim(1)
                kinv      = real(k) * inv_ldim(2)
                spaFreqSq = hinv * hinv + kinv * kinv
                ang       = atan2(real(k),real(h))
                tval      = self%eval(spaFreqSq, ang, pphshift)
                phys      = img%comp_addr_phys([h,k,0])
                call img%set_cmat_at(phys(1),phys(2),phys(3), cmplx(tval,0.))
            end do
        end do
    end subroutine ctf2img

    subroutine ctf_1stzero2img( self, img, dfx, dfy, angast, phshift )
        use simple_image, only: image
        class(ctf),     intent(inout) :: self
        class(image),   intent(inout) :: img
        real,           intent(in)    :: dfx, dfy, angast
        real, optional, intent(in)    :: phshift
        integer :: lims(3,2),h,k,phys(2),ldim(3)
        real    :: ang,tval,spaFreqSq,hinv,kinv,inv_ldim(3),pphshift
        pphshift = 0.
        if( present(phshift) ) pphshift = phshift
        ! init object
        call self%init(dfx, dfy, angast)
        ! initialize
        lims     = img%loop_lims(2)
        ldim     = img%get_ldim()
        inv_ldim = 1./real(ldim)
        ! initialize image and flag as FT
        img = cmplx(0.,0.)
        ! generate image
        do h = lims(1,1),lims(1,2)
            hinv = real(h) * inv_ldim(1)
            do k=lims(2,1),lims(2,2)
                kinv      = real(k) * inv_ldim(2)
                spaFreqSq = hinv * hinv + kinv * kinv
                ang       = atan2(real(k),real(h))
                tval      = self%eval(spaFreqSq, ang, pphshift, before1stzero=.false.)
                phys      = img%comp_addr_phys(h,k)
                call img%set_cmat_at(phys(1),phys(2),1, cmplx(tval,0.0))
            end do
        end do
    end subroutine ctf_1stzero2img

    !>  \brief  is for optimised serial application of CTF
    !!          modes: abs, ctf, flip, flipneg, neg, square
    subroutine apply_serial( self, img, mode, ctfparms )
        use simple_image, only: image
        class(ctf),       intent(inout) :: self        !< instance
        class(image),     intent(inout) :: img         !< image (output)
        character(len=*), intent(in)    :: mode        !< abs, ctf, flip, flipneg, neg, square
        type(ctfparams),  intent(in)    :: ctfparms    !< CTF parameters
        integer :: lims(3,2),h,k,phys(3),ldim(3)
        real    :: ang,tval,spaFreqSq,hinv,kinv,inv_ldim(3)
        ! init object
        call self%init(ctfparms%dfx, ctfparms%dfy, ctfparms%angast)
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
                        tval      = self%eval(spaFreqSq, ang, ctfparms%phshift)
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
                        tval      = self%eval(spaFreqSq, ang, ctfparms%phshift)
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
                        tval      = self%eval(spaFreqSq, ang, ctfparms%phshift)
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
                        tval      = self%eval(spaFreqSq, ang, ctfparms%phshift)
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
                        tval      = self%eval(spaFreqSq, ang, ctfparms%phshift)
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
                        tval      = self%eval(spaFreqSq, ang, ctfparms%phshift)
                        phys      = img%comp_addr_phys([h,k,0])
                        call img%mul_cmat_at( phys(1),phys(2),phys(3), min(1.,max(tval**2.,0.001)))
                    end do
                end do
            case DEFAULT
                write(logfhandle,*) 'mode:', mode
                THROW_HARD('unsupported mode in apply_serial')
        end select
    end subroutine apply_serial

    !>  \brief  is for applycation of wiener-like filter with SNR exponential assumption
    subroutine wienerlike_restoration(self, img, ctfparms, b_fac)
        use simple_image, only: image
        class(ctf),        intent(inout) :: self     !< instance
        class(image),      intent(inout) :: img      !< image (output)
        type(ctfparams),   intent(in)    :: ctfparms !< CTF parameters
        integer, optional, intent(in)    :: b_fac    !< b_fact value
        integer :: lims(3,2),h,k,phys(3),ldim(3), bb_fac
        real    :: ang,tval,spaFreqSq,hinv,kinv,inv_ldim(3)
        logical :: is_in_time
        integer :: nyq
        real    :: w, snr, r
        bb_fac = -50 !default value
        if(present(b_fac)) bb_fac = b_fac
        !fft controls
        if(.not. img%is_ft()) call img%fft(); is_in_time = .true.
        ! init object
        call self%init(ctfparms%dfx, ctfparms%dfy, ctfparms%angast)
        ! initialize
        lims     = img%loop_lims(2)
        ldim     = img%get_ldim()
        inv_ldim = 1./real(ldim)
        nyq = img%get_nyq()
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                hinv      = real(h) * inv_ldim(1)
                kinv      = real(k) * inv_ldim(2)
                spaFreqSq = hinv * hinv + kinv * kinv
                r         = sqrt(spaFreqSq)
                ang       = atan2(real(k),real(h))
                tval      = self%eval(spaFreqSq, ang, ctfparms%phshift)
                snr       = exp(bb_fac*r/ctfparms%smpd)
                if( snr < 1.e-10 )then
                    w = 0.
                else
                    w = tval/(tval*tval + 1./snr)
                endif
                phys = img%comp_addr_phys([h,k,0])
                call img%mul_cmat_at(phys(1),phys(2),phys(3), w)
            end do
        end do
        if( is_in_time ) call img%ifft()
    end subroutine wienerlike_restoration

    !>  \brief  is for optimised serial phase-flipping & shifting for use in prepimg4align
    subroutine phaseflip_and_shift_serial( self, img, x, y, ctfparms )
        use simple_image, only: image
        class(ctf),       intent(inout) :: self        !< instance
        class(image),     intent(inout) :: img         !< image (output)
        real,             intent(in)    :: x, y        !< defocus x/y-axis
        type(ctfparams),  intent(in)    :: ctfparms
        integer :: lims(3,2),h,k,phys(3),ldim(3), logi(3)
        real    :: ang,tval,spaFreqSq,hinv,kinv,inv_ldim(3)
        ! init object
        call self%init(ctfparms%dfx, ctfparms%dfy, ctfparms%angast)
        ! initialize
        lims     = img%loop_lims(2)
        ldim     = img%get_ldim()
        inv_ldim = 1./real(ldim)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                hinv      = real(h) * inv_ldim(1)
                kinv      = real(k) * inv_ldim(2)
                spaFreqSq = hinv * hinv + kinv * kinv
                ang       = atan2(real(k),real(h))
                tval      = self%eval(spaFreqSq, ang, ctfparms%phshift)
                logi      = [h, k, 0]
                phys      = img%comp_addr_phys(logi)
                call img%mul_cmat_at(phys, sign(1.,tval))
            end do
        end do
        ! shift image
        call img%shift2Dserial([x,y])
    end subroutine phaseflip_and_shift_serial

    ! apply CTF to image, CTF values are also returned
    subroutine eval_and_apply( self, img, imode, logi_lims, tvalsdims, tvals, dfx, dfy, angast, add_phshift, before1stzero )
        use simple_image, only: image
        class(ctf),     intent(inout) :: self        !< instance
        class(image),   intent(inout) :: img         !< modified image (output)
        integer,        intent(in)    :: imode       !< CTFFLAG_FLIP=abs CTFFLAG_YES=ctf CTFFLAG_NO=no
        integer,        intent(in)    :: logi_lims(3,2) !< logical limits
        integer,        intent(in)    :: tvalsdims(2)   !< tvals dimensions
        real,           intent(out)   :: tvals(1:tvalsdims(1),1:tvalsdims(2))
        real,           intent(in)    :: dfx         !< defocus x-axis
        real,           intent(in)    :: dfy         !< defocus y-axis
        real,           intent(in)    :: angast      !< angle of astigmatism
        real,           intent(in)    :: add_phshift !< aditional phase shift (radians), for phase plate
        logical,        intent(in)    :: before1stzero
        integer :: ldim(3),h,k,phys(2)
        real    :: ang,tval,spaFreqSq,hinv,hinvsq,kinv,inv_ldim(3)
        real    :: rh,rk
        if( imode == CTFFLAG_NO )then
            tvals = 1.0
            return
        endif
        ! initialize
        call self%init(dfx, dfy, angast)
        ldim     = img%get_ldim()
        inv_ldim = 1./real(ldim)
        do h=logi_lims(1,1),logi_lims(1,2)
            rh     = real(h)
            hinv   = rh * inv_ldim(1)
            hinvsq = hinv*hinv
            do k=logi_lims(2,1),logi_lims(2,2)
                rk = real(k)
                ! calculate CTF
                kinv      = rk * inv_ldim(2)
                spaFreqSq = hinvsq + kinv*kinv
                ang       = atan2(rk,rh)
                tval      = self%eval(spaFreqSq, ang, add_phshift, before1stzero)
                if( imode == CTFFLAG_FLIP ) tval = abs(tval)
                ! store tval and multiply image with tval
                phys = img%comp_addr_phys(h,k)
                tvals(phys(1),phys(2)) = tval
                call img%mul_cmat_at(phys(1),phys(2),1, tval)
            end do
        end do
    end subroutine eval_and_apply

    ! apply CTF to image, CTF values are also returned only before 1st zero
    subroutine eval_and_apply_before_first_zero( self, img, imode, logi_lims, tvalsdims, tvals, dfx, dfy, angast, add_phshift, maxspafreqsq )
        use simple_image, only: image
        class(ctf),     intent(inout) :: self        !< instance
        class(image),   intent(inout) :: img         !< modified image (output)
        integer,        intent(in)    :: imode       !< CTFFLAG_FLIP=abs CTFFLAG_YES=ctf CTFFLAG_NO=no
        integer,        intent(in)    :: logi_lims(3,2) !< logical limits
        integer,        intent(in)    :: tvalsdims(2)   !< tvals dimensions
        real,           intent(out)   :: tvals(1:tvalsdims(1),1:tvalsdims(2))
        real,           intent(in)    :: dfx         !< defocus x-axis
        real,           intent(in)    :: dfy         !< defocus y-axis
        real,           intent(in)    :: angast      !< angle of astigmatism
        real,           intent(in)    :: add_phshift !< aditional phase shift (radians), for phase plate
        real,           intent(out)   :: maxspafreqsq
        integer :: ldim(3),h,k,phys(2)
        real    :: ang,tval,spaFreqSq,hinv,hinvsq,kinv,rh,rk,total_phshift
        if( imode == CTFFLAG_NO )then
            tvals = 1.0
            return
        endif
        ! initialize
        maxspafreqsq = 0.0
        call self%init(dfx, dfy, angast)
        ldim = img%get_ldim()
        do h = logi_lims(1,1),logi_lims(1,2)
            rh     = real(h)
            hinv   = rh / real(ldim(1))
            hinvsq = hinv*hinv
            do k = logi_lims(2,1),logi_lims(2,2)
                rk = real(k)
                ! calculate CTF, reverse logic from eval_5
                kinv      = rk / real(ldim(2))
                spaFreqSq = hinvsq + kinv*kinv
                ang       = atan2(rk,rh)
                total_phshift = self%evalPhSh(spaFreqSq, ang, add_phshift) + self%amp_contr_const
                if( total_phshift <= CTF_FIRST_LIM )then
                    maxspafreqsq = max(maxspafreqsq, spaFreqSq)
                    phys = img%comp_addr_phys(h,k)
                    tval = sin(total_phshift)
                    if( imode == CTFFLAG_FLIP ) tval = abs(tval)
                    call img%mul_cmat_at(phys(1),phys(2),1, tval)
                    tvals(phys(1),phys(2)) = tval
                endif
            end do
        end do
    end subroutine eval_and_apply_before_first_zero

    pure elemental real function kV2wl( self ) result (wavelength)
        class(ctf), intent(in) :: self
        wavelength = 12.26 / sqrt((1000.0 * self%kV) + (0.9784*(1000.0 * self%kV)**2)/(10.0**6.0))
    end function kV2wl

    subroutine apply_convention(self, dfx, dfy, angast)
        class(ctf), intent(in)    :: self
        real,       intent(inout) :: dfx,dfy,angast
        real :: tmp
        if( angast < 0. ) angast = angast + 180.
        if( dfx < dfy )then
            tmp = dfy
            dfy = dfx
            dfx = tmp
            angast = angast + 90.
        endif
        if( angast > 180. ) angast = angast - 180.
    end subroutine apply_convention

    ! Calculate Joe's Ice Fraction Score for images, not intended for micrographs
    subroutine calc_ice_frac( self, img, ctfparms, score )
        use simple_image, only: image
        class(ctf),       intent(inout) :: self
        class(image),     intent(in)    :: img
        class(ctfparams), intent(in)    :: ctfparms
        real,             intent(out)   :: score
        real, allocatable :: res(:)
        real    :: phshift, start_freq, end_freq, ice_avg, band_avg, mag_max, mag, g
        integer :: lims(3,2), box, sh, h,k, start_find, end_find, peak_maxind, ice_maxind, hmax, kmax, cnt
        score = 0.
        if( abs(self%smpd-ctfparms%smpd) > 1.d-4) THROW_HARD('Inconsistent SMPD; calc_ice_frac')
        if( self%smpd > (ICE_BAND1/2.) ) return
        if( .not.img%is_ft() ) THROW_HARD('Image input must be in the Fourier domain!; calc_ice_frac')
        box = img%get_box()
        res = get_resarr(box, self%smpd)
        lims = img%loop_lims(2)
        call self%init(ctfparms%dfx, ctfparms%dfy, ctfparms%angast)
        phshift    = merge(ctfparms%phshift, 0. ,ctfparms%l_phaseplate)
        start_freq = sqrt(self%SpaFreqSqAtNthZero(1, phshift, deg2rad(ctfparms%angast)))
        end_freq   = sqrt(self%SpaFreqSqAtNthZero(2, phshift, deg2rad(ctfparms%angast)))
        call get_find_at_crit(size(res), res, ICE_BAND1, ice_maxind)
        start_find = max(1,     ice_maxind - 3)
        end_find   = min(box/2, ice_maxind + 3)
        hmax     = -1
        kmax     = -1
        mag_max  = -1.
        band_avg = 0.
        cnt      = 0
        do k = lims(2,1),lims(2,2)
            do h = lims(1,1),lims(1,2)
                sh = nint(hyp(h,k))
                g  = real(sh) / real(box)
                if( g > start_freq .and. g < end_freq )then
                    band_avg = band_avg + csq_fast(img%get_fcomp2D(h,k))
                    cnt      = cnt + 1
                else if( sh >= start_find .and. sh < end_find )then
                    mag = csq_fast(img%get_fcomp2D(h,k))
                    if( mag > mag_max )then
                        hmax = h
                        kmax = k
                        mag_max = mag
                    endif
                endif
            end do
        end do
        if( cnt < 1 ) return
        band_avg = band_avg / real(cnt)
        ice_avg  = 0.
        do k = kmax-1,kmax+1
            do h = hmax-1,hmax+1
                ice_avg = ice_avg + csq_fast(img%get_fcomp2D(h,k))
            enddo
        enddo
        ice_avg = ice_avg / 9.
        score   = ice_avg / band_avg
        if( score > 0.5 ) score = (ice_avg - 0.5*band_avg) / band_avg
    end subroutine calc_ice_frac

end module simple_ctf