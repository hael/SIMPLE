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
use simple_image, only: image
use simple_memoize_ft_maps
implicit none

public :: ctf, memoize4ctf_apply, unmemoize4ctf_apply
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
    procedure, private :: evalPhSh
    procedure, private :: eval_1, eval_2, eval_3, eval_4
    generic            :: eval => eval_1, eval_2, eval_3, eval_4
    procedure, private :: eval_sign
    procedure, private :: eval_df
    procedure          :: nextrema
    procedure          :: spafreqsqatnthzero
    procedure          :: apply
    procedure          :: ctf2img
    procedure          :: apply_serial
    procedure          :: eval_and_apply
    procedure, private :: kV2wl
    procedure          :: apply_convention
    procedure          :: calc_ice_frac
end type ctf

interface ctf
    module procedure constructor
end interface

type mem_ctf_apply
    real    :: spaFreqSq, ang
    integer :: ph1, ph2
end type mem_ctf_apply

type(mem_ctf_apply), allocatable :: mem_ctf_apply_mat(:,:)
integer                          :: lims_memo(3,2)

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

    !>  \brief  is for applying CTF to an image
    subroutine apply( self, img, dfx, mode, dfy, angast, bfac )
        class(ctf),       intent(inout) :: self   !< instance
        class(image),     intent(inout) :: img    !< image (output)
        real,             intent(in)    :: dfx    !< defocus x-axis
        character(len=*), intent(in)    :: mode   !< abs, ctf, flip, flipneg, neg, square
        real, optional,   intent(in)    :: dfy    !< defocus y-axis
        real, optional,   intent(in)    :: angast !< angle of astigmatism
        real, optional,   intent(in)    :: bfac   !< bfactor
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
        ! optionals
        if(present(dfy))    ctfvars%dfy     = dfy
        if(present(angast)) ctfvars%angast  = angast
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

    subroutine memoize4ctf_apply( img )
        class(image), intent(inout) :: img   !< image (output)
        integer :: h,k,phys(2),ldim(3)
        real    :: rh,hinv,hinvsq,rk,kinv,inv_ldim(3)
        if( OMP_IN_PARALLEL() )then
            THROW_HARD('No memoization inside OpenMP regions')
        endif
        ! initialize
        ldim = img%get_ldim()
        if( any(ldim == 0) .or. ldim(3) /= 1 ) THROW_HARD('Image incompatible with CTF application')
        lims_memo = img%loop_lims(2)
        if( allocated(mem_ctf_apply_mat) ) deallocate(mem_ctf_apply_mat)
        allocate(mem_ctf_apply_mat(lims_memo(1,1):lims_memo(1,2),lims_memo(2,1):lims_memo(2,2)))
        inv_ldim = 1./real(ldim)
        ! pre-calulate
        do h=lims_memo(1,1),lims_memo(1,2)
            rh     = real(h)
            hinv   = rh * inv_ldim(1)
            hinvsq = hinv * hinv
            do k=lims_memo(2,1),lims_memo(2,2)
                rk   = real(k)
                kinv = rk * inv_ldim(2)
                phys = img%comp_addr_phys(h,k)
                ! stash
                mem_ctf_apply_mat(h,k)%spaFreqSq = hinvsq + kinv * kinv
                mem_ctf_apply_mat(h,k)%ang       = atan2(rk,rh)
                mem_ctf_apply_mat(h,k)%ph1       = phys(1)
                mem_ctf_apply_mat(h,k)%ph2       = phys(2)
            end do
        end do
    end subroutine memoize4ctf_apply

    subroutine unmemoize4ctf_apply
        if( allocated(mem_ctf_apply_mat) ) deallocate(mem_ctf_apply_mat)
    end subroutine unmemoize4ctf_apply

    ! local fast CTF kernel (numerically equivalent refactor)
    pure elemental real function fast_ctf_kernel( h, k, dfx, dfy, angast, amp_contr_const, wl, half_wl2_cs ) result( tval )
        integer, intent(in) :: h, k
        real,    intent(in) :: dfx, dfy, angast, amp_contr_const, wl, half_wl2_cs
        real :: df, evalPhSh, sum_df, diff_df, cterm, pi_wl_s2
        ! Defocus term: exactly the same expression, just factored
        sum_df  = dfx + dfy
        diff_df = dfx - dfy
        cterm   = cos( 2.0 * (mem_ctf_apply_mat(h,k)%ang - angast) )
        df      = 0.5 * ( sum_df + cterm * diff_df )
        ! Phase shift term: preserve the same evaluation structure
        pi_wl_s2 = PI * wl * mem_ctf_apply_mat(h,k)%spaFreqSq
        evalPhSh = pi_wl_s2 * (df - half_wl2_cs * mem_ctf_apply_mat(h,k)%spaFreqSq)
        tval     = sin( evalPhSh + amp_contr_const )
    end function fast_ctf_kernel

    ! local fast CTF kernel (numerically equivalent refactor)
    pure elemental real function ft_map_ctf_kernel( h, k, sum_df, diff_df, angast, amp_contr_const, wl, half_wl2_cs ) result( tval )
        integer, intent(in) :: h, k
        real,    intent(in) :: sum_df, diff_df, angast, amp_contr_const, wl, half_wl2_cs
        real :: df, PhSh, cterm, pi_wl_s2
        ! Defocus term: exactly the same expression as eval_df,
        ! but factored the sum & difference terms are memoized
        cterm   = cos( 2.0 * (ft_map_astigang(h,k) - angast) )
        df      = 0.5 * ( sum_df + cterm * diff_df )
        ! Phase shift term: preserve the same evaluation structure
        pi_wl_s2 = PI * wl * ft_map_spaFreqSq(h,k)
        PhSh     = pi_wl_s2 * (df - half_wl2_cs * ft_map_spaFreqSq(h,k))
        tval     = sin( PhSh + amp_contr_const )
    end function ft_map_ctf_kernel

    !>  \brief  is for generating an image of CTF
    subroutine ctf2img( self, img, dfx_in, dfy_in, angast_in )
        class(ctf),       intent(inout) :: self
        class(image),     intent(inout) :: img
        real,             intent(in)    :: dfx_in, dfy_in, angast_in
        integer :: h,k
        real    :: dfx, dfy, angast, amp_contr_const, wl, half_wl2_cs, tval
        ! flag image as FT
        call img%set_ft(.true.)
        ! init
        call self%init(dfx_in, dfy_in, angast_in) ! conversions
        wl              = self%wl
        half_wl2_cs     = 0.5 * wl * wl * self%cs
        dfx             = self%dfx
        dfy             = self%dfy
        angast          = self%angast
        amp_contr_const = self%amp_contr_const
        do h=lims_memo(1,1),lims_memo(1,2)
            do k=lims_memo(2,1),lims_memo(2,2)
                ! calculate CTF
                tval = fast_ctf_kernel(h, k, dfx, dfy, angast, amp_contr_const, wl, half_wl2_cs)
                ! set cmat
                call img%set_cmat_at(mem_ctf_apply_mat(h,k)%ph1,mem_ctf_apply_mat(h,k)%ph2,1, cmplx(tval,0.))
            end do
        end do
    end subroutine ctf2img

    !>  \brief  is for optimised serial application of CTF
    !!          modes: abs, ctf, flip, flipneg, neg, square
    subroutine apply_serial( self, img, mode, ctfparms )
        class(ctf),       intent(inout) :: self     !< instance
        class(image),     intent(inout) :: img      !< image (output)
        character(len=*), intent(in)    :: mode     !< abs, ctf, flip, flipneg, neg, square
        type(ctfparams),  intent(in)    :: ctfparms !< CTF parameters
        real, parameter :: ZERO = 0., ONE = 1.0, MIN_SQUARE = 0.001
        integer :: h,k
        real    :: dfx, dfy, angast, amp_contr_const, wl, half_wl2_cs, tval, t
        logical :: is_abs,is_ctf,is_flip,is_flipneg,is_neg,is_square
        ! Convert mode string to logical flags (SIMD-compatible)
        is_abs     = (mode == 'abs')
        is_ctf     = (mode == 'ctf')
        is_flip    = (mode == 'flip')
        is_flipneg = (mode == 'flipneg')
        is_neg     = (mode == 'neg')
        is_square  = (mode == 'square')
        ! initialize
        call self%init(ctfparms%dfx, ctfparms%dfy, ctfparms%angast) ! conversions
        wl              = self%wl
        half_wl2_cs     = 0.5 * wl * wl * self%cs
        dfx             = self%dfx
        dfy             = self%dfy
        angast          = self%angast
        amp_contr_const = self%amp_contr_const
        do h = lims_memo(1,1), lims_memo(1,2)
            do k = lims_memo(2,1), lims_memo(2,2)
                ! calculate CTF
                tval = fast_ctf_kernel(h, k, dfx, dfy, angast, amp_contr_const, wl, half_wl2_cs)
                ! SIMD-compatible mode handling using masked operations
                t = merge(       abs(tval),                     ZERO, is_abs)     + &
                    merge(           tval,                      ZERO, is_ctf)     + &
                    merge( sign(ONE, tval),                     ZERO, is_flip)    + &
                    merge(-sign(ONE, tval),                     ZERO, is_flipneg) + &
                    merge(          -tval,                      ZERO, is_neg)     + &
                    merge(min(ONE, max(tval*tval, MIN_SQUARE)), ZERO, is_square)
                ! Multiply (this is the key operation)
                call img%mul_cmat_at(mem_ctf_apply_mat(h,k)%ph1,mem_ctf_apply_mat(h,k)%ph2,1, t)
            end do
        end do
    end subroutine apply_serial

    ! apply CTF to image, CTF values are also returned
    subroutine eval_and_apply( self, img, imode, tvalsdims, tvals, dfx_in, dfy_in, angast_in )
        use simple_image, only: image
        class(ctf),     intent(inout) :: self         !< instance
        class(image),   intent(inout) :: img          !< modified image (output)
        integer,        intent(in)    :: imode        !< CTFFLAG_FLIP=abs CTFFLAG_YES=ctf CTFFLAG_NO=no
        integer,        intent(in)    :: tvalsdims(2) !< tvals dimensions
        real,           intent(out)   :: tvals(1:tvalsdims(1),1:tvalsdims(2))
        real,           intent(in)    :: dfx_in       !< defocus x-axis
        real,           intent(in)    :: dfy_in       !< defocus y-axis
        real,           intent(in)    :: angast_in    !< angle of astigmatism
        complex(kind=c_float_complex), pointer :: pcmat(:,:,:)
        real    :: angast, amp_contr_const, wl, half_wl2_cs
        real    :: sum_df, diff_df,tval
        integer :: h,k, physh,physk
        logical :: l_flip
        if( imode == CTFFLAG_NO )then
            tvals = 1.0
            return
        endif
        ! initialize
        call self%init(dfx_in, dfy_in, angast_in) ! conversions
        wl              = self%wl
        half_wl2_cs     = 0.5 * wl * wl * self%cs
        angast          = self%angast
        amp_contr_const = self%amp_contr_const    
        l_flip          = imode == CTFFLAG_FLIP
        sum_df          = self%dfx + self%dfy
        diff_df         = self%dfx - self%dfy
        call img%get_cmat_ptr(pcmat)
        do h = ft_map_lims(1,1), ft_map_lims(1,2)
            do k = ft_map_lims(2,1), ft_map_lims(2,2)
                ! calculate CTF
                tval = merge(abs(ft_map_ctf_kernel(h, k, sum_df, diff_df, angast, amp_contr_const, wl, half_wl2_cs)),&
                                &ft_map_ctf_kernel(h, k, sum_df, diff_df, angast, amp_contr_const, wl, half_wl2_cs), l_flip)
                ! store tval and multiply image with tval
                physh = ft_map_phys_addrh(h,k)
                physk = ft_map_phys_addrk(h,k)
                tvals(physh, physk)    = tval
                pcmat(physh, physk, 1) = tval * pcmat(physh, physk, 1)
            end do
        end do
        nullify(pcmat)
    end subroutine eval_and_apply

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

    ! Calculate Joe's Ice Fraction Score for images, not intended for micrographs
    subroutine calc_ice_frac( self, img, ctfparms, score )
        class(ctf),       intent(inout) :: self
        class(image),     intent(in)    :: img
        class(ctfparams), intent(in)    :: ctfparms
        real,             intent(out)   :: score
        real, allocatable :: res(:)
        real    :: phshift, start_freq, end_freq, ice_avg, band_avg, mag_max, mag, g
        integer :: lims(3,2), box, sh, h,k, start_find, end_find, ice_maxind, hmax, kmax, cnt
        score = 0.
        if( abs(self%smpd-ctfparms%smpd) > 1.d-4) THROW_HARD('Inconsistent SMPD; calc_ice_frac')
        if( self%smpd > (ICE_BAND1/2.) ) return
        if( .not.img%is_ft() ) THROW_HARD('Image input must be in the Fourier domain!; calc_ice_frac')
        box = img%get_box()
        res = get_resarr(box, self%smpd)
        lims = img%loop_lims(2)
        call self%init(ctfparms%dfx, ctfparms%dfy, ctfparms%angast)
        phshift    = 0.
        start_freq = sqrt(self%SpaFreqSqAtNthZero(1, phshift, deg2rad(ctfparms%angast)))
        end_freq   = sqrt(self%SpaFreqSqAtNthZero(2, phshift, deg2rad(ctfparms%angast)))
        ice_maxind = get_find_at_res(res, ICE_BAND1)
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
