!@descr: defines the Contrast Transfer Function (CTF) of the electron microscope
! defines the Contrast Transfer Function (CTF) of the electron microscope
! based on a class used in CTFFIND4, developed by Alexis Rohou and Nikolaus
! Grigorieff at Janelia Farm. The below copyright statement therefore
! needs to be included here:
! Copyright 2014 Howard Hughes Medical Institute
! All rights reserved
! Use is subject to Janelia Farm Research Campus Software Copyright 1.1
! license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
module simple_ctf
use simple_core_module_api
implicit none

public :: ctf
private
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
    procedure          :: get_ctfvars
    procedure, private :: evalPhSh
    procedure, private :: eval_1, eval_2, eval_3, eval_4
    generic            :: eval => eval_1, eval_2, eval_3, eval_4
    procedure, private :: eval_sign
    procedure, private :: eval_df
    procedure          :: nextrema
    procedure          :: spafreqsqatnthzero
    procedure, private :: kV2wl
    procedure          :: apply_convention
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

    pure function get_ctfvars( self ) result(ctfvals)
        class(ctf), intent(in) :: self
        type(ctfvars) :: ctfvals
        ctfvals%smpd            = self%smpd
        ctfvals%kv              = self%kv
        ctfvals%cs              = self%cs
        ctfvals%wl              = self%wl
        ctfvals%amp_contr_const = self%amp_contr_const
        ctfvals%dfx             = self%dfx
        ctfvals%dfy             = self%dfy
        ctfvals%angast          = self%angast
    end function get_ctfvars

    !>  \brief returns the argument (radians) to the ctf
    !!
    !!  We follow the convention, like the rest of the cryo-EM/3DEM field, that underfocusing the objective lens
    !!  gives rise to a positive phase shift of scattered electrons, whereas the spherical aberration gives a
    !!  negative phase shift of scattered electrons
    elemental real function evalPhSh( self, spaFreqSq, ang, add_phshift )
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
    elemental real function eval_3( self, spaFreqSq, ang )
        class(ctf), intent(in) :: self        !< instance
        real,       intent(in) :: spaFreqSq   !< squared reciprocal pixels
        real,       intent(in) :: ang         !< Angle at which to compute the CTF (radians)
        ! compute phase shift + amplitude constrast term & compute value of CTF, assuming white particles
        eval_3 = sin( self%evalPhSh(spaFreqSq, ang, 0.) + self%amp_contr_const )
    end function eval_3

    !>  \brief Returns the CTF with pre-initialize parameters
    elemental real function eval_4( self, spaFreqSq, ang, add_phshift )
        class(ctf), intent(in) :: self        !< instance
        real,       intent(in) :: spaFreqSq   !< squared reciprocal pixels
        real,       intent(in) :: ang         !< Angle at which to compute the CTF (radians)
        real,       intent(in) :: add_phshift !< aditional phase shift (radians), for phase plate
        ! compute phase shift + amplitude constrast term & compute value of CTF, assuming white particles
        eval_4 = sin( self%evalPhSh(spaFreqSq, ang, add_phshift) + self%amp_contr_const )
    end function eval_4

    !>  \brief Returns the sign of the CTF with pre-initialize parameters
    elemental integer function eval_sign( self, spaFreqSq, ang, add_phshift )
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
    elemental real function eval_df( self, ang )
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

    elemental real function kV2wl( self ) result( wavelength )
        class(ctf), intent(in) :: self
        wavelength = 12.26 / sqrt((1000.0 * self%kV) + (0.9784*(1000.0 * self%kV)**2)/(10.0**6.0))
    end function kV2wl

    subroutine apply_convention( self, dfx, dfy, angast )
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

end module simple_ctf
