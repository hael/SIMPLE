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
!use simple_params_l, only: CTFMODETYPE
implicit none

public :: ctf, test_ctf
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
    ! INITIALISER
    procedure, private :: init
    ! CTF EVALUATION
    procedure          :: eval
    procedure, private :: evalPhSh
    procedure, private :: eval_df
    ! CTF APPLICATION
    procedure          :: ctf2img
    procedure          :: ctf2spec
    procedure          :: apply
    procedure          :: apply_and_shift
    ! CALCULATORS
    procedure          :: kV2wl
    procedure          :: freqOfAZero
    procedure, private :: solve4PhSh
    procedure, private :: sqFreq4PhSh
    procedure, private :: sqSf4PhSh
end type ctf

interface ctf
    module procedure constructor
end interface

contains

    !>  \brief  is a constructor
    function constructor( smpd, kV, Cs, amp_contr ) result( self )
        real, intent(in) :: smpd !< sampling distance
        real, intent(in) :: kV   !< accelleration voltage
        real, intent(in) :: Cs   !< constant
        real, intent(in) :: amp_contr !< amplitude contrast
        type(ctf) :: self
        ! set constants
        self%kV        = kV
        self%wl        = self%kV2wl(kV) / smpd
        self%Cs        = (Cs*1.0e7)/smpd
        self%amp_contr = amp_contr
        self%smpd      = smpd
        ! compute derived constants (phase and amplitude contrast weights)
        self%phaseq = sqrt(1.-amp_contr**2.)
        self%ampliq = amp_contr
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

    ! CTF EVALUATORS

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
    function eval( self, spaFreqSq, dfx, dfy, angast, ang ) result ( val )
        class(ctf), intent(inout) :: self      !< instance
        real,       intent(in)    :: spaFreqSq !< squared reciprocal pixels
        real,       intent(in)    :: dfx       !< Defocus along first axis (micrometers)
        real,       intent(in)    :: dfy       !< Defocus along second axis (for astigmatic CTF, dfx .ne. dfy) (micrometers)
        real,       intent(in)    :: angast    !< Azimuth of first axis. 0.0 means axis is at 3 o'clock. (radians)
        real,       intent(in)    :: ang       !< Angle at which to compute the CTF (radians)
        real :: val, phase_shift
        ! initialize the CTF object, using the input parameters
        call self%init(dfx, dfy, angast)
        ! compute phase shift
        phase_shift = self%evalPhSh(spaFreqSq, ang)
        ! compute value of CTF, assuming white particles
        val = (self%phaseq * sin(phase_shift) + self%ampliq * cos(phase_shift))
    end function eval

    !>  \brief returns the argument (radians) to the sine and cosine terms of the ctf
    !!
    !!  We follow the convention, like the rest of the cryo-EM/3DEM field, that underfocusing the objective lens
    !!  gives rise to a positive phase shift of scattered electrons, whereas the spherical aberration gives a
    !!  negative phase shift of scattered electrons
    pure elemental real function evalPhSh( self, spaFreqSq, ang ) result( phase_shift )
        class(ctf), intent(in) :: self      !< instance
        real,       intent(in) :: spaFreqSq !< square of spatial frequency at which to compute the ctf (1/pixels^2)
        real,       intent(in) :: ang       !< angle at which to compute the ctf (radians)
        real :: df                          !< defocus at point at which we're evaluating the ctf
        !! compute the defocus
        df = self%eval_df(ang)
        !! compute the ctf argument
        phase_shift = pi * self%wl * spaFreqSq * (df - 0.5 * self%wl**2 * spaFreqSq * self%Cs)
    end function evalPhSh

    !>  \brief  Return the effective defocus given the pre-set CTF parameters (from CTFFIND4)
    pure elemental real function eval_df( self, ang ) result (df)
        class(ctf), intent(in) :: self !< instance
        real,       intent(in) :: ang  !< angle at which to compute the defocus (radians)
        df = 0.5 * (self%dfx + self%dfy + cos(2.0 * (ang - self%angast)) * (self%dfx - self%dfy))
    end function eval_df

    ! CTF APPLICATION ROUTINES

    !>  \brief  is for making a CTF image
    !!          modes: abs, ctf, flip, flipneg, neg, square
    subroutine ctf2img( self, img, dfx, mode, dfy, angast, bfac )
        use simple_image, only: image
        class(ctf),                 intent(inout) :: self        !< instance
        class(image),               intent(inout) :: img         !< image (output)
        real,                       intent(in)    :: dfx         !< defocus x-axis
        character(len=*),           intent(in)    :: mode        !< abs, ctf, flip, flipneg, neg, square
        real,             optional, intent(in)    :: dfy         !< defocus y-axis
        real,             optional, intent(in)    :: angast      !< angle of astigmatism
        real,             optional, intent(in)    :: bfac        !< bfactor
        integer :: lims(3,2),h,k,phys(3),ldim(3)
        real    :: ang, tval, ddfy, aangast, spaFreqSq, hinv
        real    :: kinv, inv_ldim(3)
        if( img%is_3d() )then
            print *, 'ldim: ', img%get_ldim()
            stop 'Only 4 2D images; ctf2img; simple_ctf'
        endif
        ! set CTF params
        ddfy = dfx
        if( present(dfy) ) ddfy = dfy
        aangast = 0.
        if( present(angast) ) aangast = angast
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
                tval      = self%eval(spaFreqSq, dfx, ddfy, aangast, ang)
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

    !>  \brief  is for making a CTF spectrum (1d)
    !!          modes: abs, ctf, flip, flipneg, neg, square
    function ctf2spec( self, kfromto, box, dfx, mode, bfac ) result(spec)
        use simple_image,  only: image
        class(ctf),       intent(inout) :: self       !< instance
        integer,          intent(in)    :: kfromto(2) !< Fourier index range
        integer,          intent(in)    :: box        !< box size
        real,             intent(in)    :: dfx        !< defocus
        character(len=*), intent(in)    :: mode       !< abs, ctf, flip, flipneg, neg, square
        real, optional,   intent(in)    :: bfac       !< bfactor
        real, allocatable :: spec(:)
        integer           :: k
        real              :: tval,res,wght,spaFreqSq,kinv
        ! init object
        call self%init(dfx, dfx, 0.)
        ! allocate spectrum
        allocate( spec(kfromto(1):kfromto(2)), stat=alloc_stat )
        if(alloc_stat /= 0) allocchk('In: ctf2spec; simple_ctf')
        ! do the work
        do k=kfromto(1),kfromto(2) ! loop over resolution range
            kinv = real(k)/real(box)
            spaFreqSq = kinv*kinv
            tval = self%eval(spaFreqSq, dfx, dfx, 0., 0.) ! CFT-value at k
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
            spec(k) = tval
        end do
        if( present(bfac) )then
            do k=kfromto(1),kfromto(2)                 ! loop over resolution range
                res = real(k)/(real(box)*self%smpd)    ! assuming square dimensions
                wght = max(0.,exp(-(bfac/4.)*res*res)) ! b-factor weight
                spec(k) = spec(k)*wght                 ! apply it
            end do
        endif
    end function ctf2spec

    !>  \brief  is for applying CTF to an image
    subroutine apply( self, img, dfx, mode, dfy, angast, bfac)
        use simple_image, only: image
        class(ctf),                 intent(inout) :: self        !< instance
        class(image),               intent(inout) :: img         !< image (output)
        real,                       intent(in)    :: dfx         !< defocus x-axis
        character(len=*),           intent(in)    :: mode        !< abs, ctf, flip, flipneg, neg, square
        real,             optional, intent(in)    :: dfy         !< defocus y-axis
        real,             optional, intent(in)    :: angast      !< angle of astigmatism
        real,             optional, intent(in)    :: bfac        !< bfactor
        integer     :: ldim(3), ldim_pd(3)
        type(image) :: ctfimg, img_pd
        ldim = img%get_ldim()
        if( img%is_3d() )then
            print *, 'ldim: ', ldim
            stop 'Only 4 2D images; apply; simple_ctf'
        endif
        if( img%is_ft() )then
            call ctfimg%new(ldim, self%smpd)
            call self%ctf2img(ctfimg, dfx, mode, dfy, angast, bfac)
            call img%mul(ctfimg)
        else
            ldim_pd(1:2) = 2*ldim(1:2)
            ldim_pd(3)   = 1
            call ctfimg%new(ldim_pd, self%smpd)
            call img_pd%new(ldim_pd, self%smpd)
            call self%ctf2img(ctfimg, dfx, mode, dfy, angast, bfac)
            call img%pad_mirr(img_pd)
            call img_pd%fwd_ft
            call img_pd%mul(ctfimg)
            call img_pd%bwd_ft
            call img_pd%clip(img)
            call img_pd%kill
        endif
        call ctfimg%kill
    end subroutine apply

    !>  \brief  is for applying CTF to an image
    subroutine apply_and_shift( self, img, imgctfsq, x, y, dfx, mode, dfy, angast )
        use simple_image, only: image
        class(ctf),                 intent(inout) :: self     !< instance
        class(image),               intent(inout) :: img      !< modified image (output)
        class(image),               intent(inout) :: imgctfsq !< CTF**2 image (output)
        real,                       intent(in)    :: x, y     !< rotational origin shift
        real,                       intent(in)    :: dfx      !< defocus x-axis
        character(len=*),           intent(in)    :: mode     !< abs, ctf, flip, flipneg, neg, square
        real,             optional, intent(in)    :: dfy      !< defocus y-axis
        real,             optional, intent(in)    :: angast   !< angle of astigmatism
        integer :: ldim(3),logi(3),imode,lims(3,2),h,k,phys(3)
        real    :: ang,tval,ddfy,aangast,spaFreqSq,hinv,kinv,inv_ldim(3),tvalsq
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
                if( imode <  3 ) tval = self%eval(spaFreqSq, dfx, ddfy, aangast, ang)
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

    ! CALCULATORS

    !> \brief   Convert acceleration voltage in kV into electron wavelength in Angstroms
    !! \details  DeBroglie showed that  \f$ \lambda = \frac{h}{mv} \f$, where
    !! h is Planck’s constant  \f$ \SI{6.626e-34}{\joule\second} \f$
    !! mv is momentum (mass times velocity), which is  \f$ \sqrt{2 m_e q V} \f$
    !! mass of an electron is \f$ \SI{9.1e-31}{\kilo\gram} \f$ and charge of an electron is  \f$ \SI{1.6e-19}{\coulomb} \f$
    !! Equation for converting acceleration voltage into electron wavelength in Angstroms (\f$ \si{\angstrom} \f$) :
    !! \f[ \lambda (\si{\metre}) =  \frac{12.26 \times 10^{-10}}{\sqrt{V}} \f]
    !! Due to relativistic effects Lorentz law applies:
    !!  \f[ \lambda (\si{\metre}) =  \frac{12.26 \times 10^{-10}}{\sqrt{V}} \times \frac{1}{\sqrt{1+\frac{e V}{2 m_e c^2}}} \f]
    !! where c is the speed of light,\f$ \sim \SI{3e8}{\metre\per\second} \f$
    !!  \f[ \lambda (\si{\angstrom}) =  \frac{12.26}{\sqrt{\left(10^3 kV + 0.9784 \times (10^3 kV)^2 \right) / 10^6 }} 
    !! \f]
    pure elemental real function kV2wl( self, accelkv ) result (wavelength)
        class(ctf), intent(in) :: self      !< instance
        real,       intent(in) :: accelkv        !< acceleration voltage 
        wavelength = 12.26 / sqrt((1000.0*accelkv) + (0.9784*(1000.0*accelkv)**2)/(10.0**6.0))
    end function kV2wl

    !>  \brief  Return the spatial frequency (reciprocal pixels) of the n-th zero, where n is given as which_zero
    real function freqOfAZero( self, dfx, dfy, angast, ang, which_zero )
        class(ctf), intent(inout) :: self       !< instance
        real,       intent(in)    :: dfx        !< Defocus along first axis (micrometers)
        real,       intent(in)    :: dfy        !< Defocus along second axis (for astigmatic CTF, dfx .ne. dfy) (micrometers)
        real,       intent(in)    :: angast     !< Azimuth of first axis. 0.0 means axis is at 3 o'clock. (radians)
        real,       intent(in)    :: ang        !< Radians. Angle at which to evaluate the CTF
        integer,    intent(in)    :: which_zero !< Number (starting at 1) of the zero you are interested in
        real, allocatable         :: phase_shifts_of_zeroes(:)
        real, allocatable         :: sq_sfs_of_zeroes(:)
        integer                   :: num_freqs
        ! initialize the CTF object, using the input parameters
        call self%init(dfx, dfy, angast)
        ! allocate
        allocate(phase_shifts_of_zeroes(which_zero*8))
        allocate(sq_sfs_of_zeroes(2*size(phase_shifts_of_zeroes)))
        ! Solve the CTF for the phase shift - we will get lots of phase shifts, some of which will not be possible with current CTF
        call self%solve4PhSh(phase_shifts_of_zeroes,0.)
        call self%sqFreq4PhSh(phase_shifts_of_zeroes,ang,sq_sfs_of_zeroes,num_freqs)
        freqOfAZero = sqrt(sq_sfs_of_zeroes(which_zero))
    end function freqOfAZero

    !>  \brief  Find the ctf argument (phase shift) values at which CTF=ctf_value
    !!
    !!  According to Wolfram Alpha, the solutions to \f$a*\cos(t) + b*\sin(t) = c \f$ are:
    !!  1)  \f$ t = 2 (\atan((b-\sqrt(a^2+b^2-c^2))/(a+c))+\pi*n) \f$,   with n an integer
    !!  2)  \f$ t = 2 (\atan((b+\sqrt(a^2+b^2-c^2))/(a+c))+\pi*n) \f$,   with n an integer
    !!
    !!  The above two solutions only hold if:
    !!                \f$  a+c != 0 \f$
    !!                \f$ -b*\sqrt(a^2+b^2-c^2)+a^2+ac+b^2 != 0 \f$  [for solution 1]
    !!                \f$  b*\sqrt(a^2+b^2-c^2)+a^2+ac+b^2 != 0 \f$  [for solution 2]
    !!
    !!  In our case, a = - ampl_const
    !!           and b = \f$ - \sqrt(1-\mathrm{ampl_const}^2) \f$
    !!               c = ctf_value
    !!  and t is the "argument" to the ctf (i.e. the phase shift), which can be "converted" to a spatial frequency
    subroutine solve4PhSh( self, sols, ctf_value )
        class(ctf),        intent(in)    :: self      !< instance
        real, allocatable, intent(inout) :: sols(:)   !< solutions
        real, optional,    intent(in)    :: ctf_value !< value of the ctf at the solutions
        integer :: i,j
        real    :: a,b,c,sqrt_arg,sqrt_val,sol_first_term_1,sol_first_term_2
        if( .not. allocated(sols) ) stop 'array of solutions is not allocated; solve4PhSh; simple_ctf'
        c = 0.
        if (present(ctf_value)) c = ctf_value
        a = -self%amp_contr
        b = -sqrt(1-self%amp_contr**2)
        ! Note that since 0.0 <= ampl_cont <= 1.0:
        ! a**2+b**2 = 1.0
        ! sqrt_arg = a**2+b**2-c**2
        sqrt_arg = max(1.0-c**2,0.0)
        sqrt_val = sqrt(sqrt_arg)
        sol_first_term_1 = atan((b-sqrt_val)/(a+c))
        sol_first_term_2 = atan((b+sqrt_val)/(a+c))
        ! Loop over the zeroes
        do i=1,size(sols)/2
            j = (i-1)*2+1
            sols(j) = 2.*(sol_first_term_1+pi*real(i-1-size(sols)/4))
            j = (i-1)*2+2
            sols(j) = 2.*(sol_first_term_2+pi*real(i-size(sols)/4))
        enddo
    end subroutine solve4PhSh

    !>  \brief  sqFreq4PhSh Compute set of squared spatial frequencies at which the given phase shift 
    !!          is obtained by the CTF
    subroutine sqFreq4PhSh( self, phase_shifts, ang, spaFreqSq, nsols )
        use simple_math, only: hpsort
        class(ctf), intent(in)    :: self
        real,       intent(in)    :: phase_shifts(:) !< CTF phase shifts
        real,       intent(in)    :: ang             !< angle (radians)
        real,       intent(inout) :: spaFreqSq(:)    !< assumed to be 2x size of phase_shift array
        integer,    intent(out)   :: nsols           !< num of solutions
        integer :: i
        real    :: curr_sols(2)
        integer :: curr_nsols
        nsols = 0
        do i=1,size(phase_shifts)
            call self%sqSf4PhSh(phase_shifts(i), ang, curr_sols, curr_nsols)
            if( curr_nsols .gt. 0 )then
                spaFreqSq(nsols+1:nsols+curr_nsols) = curr_sols(1:curr_nsols)
                nsols = nsols+curr_nsols
            endif
        enddo
        call hpsort(nsols,spaFreqSq(1:nsols))
    end subroutine sqFreq4PhSh

    !>  \brief  Given the argument to the ctf (phase shift, in radians), return the squared spatial frequency
    subroutine sqSf4PhSh( self, phase_shift, ang, sq_sf, nsols )
        class(ctf), intent(in)  :: self        !< instance
        real,       intent(in)  :: phase_shift !< phase shift (radians)
        real,       intent(in)  :: ang         !< angle at which to compute the ctf (radians)
        real,       intent(out) :: sq_sf(2)    !< squared spatial frequency (pixels^-2)
        integer,    intent(out) :: nsols       !< there may be 0, 1 or 2 solutions
        real :: df
        real :: a,b,c
        real :: det
        ! compute the defocus
        df = self%eval_df(ang)
        ! solve
        a = pi*self%wl*df
        b = 0.5*pi*self%wl**3.*self%Cs
        c = phase_shift
        det = a**2.-4.*b*c
        if( det .ge. 0.0 )then
            sq_sf(1) = (a+sqrt(det))/(-2.*b)
            sq_sf(2) = (a-sqrt(det))/(-2.*b)
            if( det .eq. 0.0 )then
                nsols = 1
            else
                nsols = 2
            endif
        else
            ! no solutions
            sq_sf = 0.
            nsols = 0
        endif
        ! squared spatial frequencies must be >=0.; let's remove bad solutions
        select case( nsols )
            case( 1 )
                if( sq_sf(1) .lt. 0. )then
                    nsols = 0
                endif
            case (2)
                if( sq_sf(2) .lt. 0. .and. sq_sf(1) .ge. 0. )then
                    nsols = 1
                else if( sq_sf(1) .lt. 0. .and. sq_sf(2) .ge. 0. )then
                    nsols = 1
                    sq_sf(1) = sq_sf(2)
                else if( sq_sf(1) .lt. 0. .and. sq_sf(2) .lt. 0. )then
                    nsols = 0
                endif
        end select
    end subroutine sqSf4PhSh

    !>  \brief  ctf class unit test
    subroutine test_ctf( nthr )
        use simple_image,  only: image
        use simple_rnd,    only: ran3
        integer, intent(in) :: nthr
        integer, parameter :: BOX_SIZES(4) = [64,128,256,512]
        real,    parameter :: SMPDS(4)     = [3.7,2.2,1.07,0.6]
        real,    parameter :: MAX_DF       = 10.
        real,    parameter :: CS           = 2.7
        real,    parameter :: FRACA        = 0.07
        real,    parameter :: KVS(3)       = [120.,200.,300.]
        integer, parameter :: NTESTS       = 10000
        real               :: dfx, dfy, angast
        integer            :: ibox, itst, ldim(3), ikv, ntot, cnt
        type(image)        :: img
        type(ctf)          :: tfun
        write(*,'(a)') '**info(simple_ctf_unit_test): testing the ctf2img routine'
!$      call omp_set_num_threads(nthr)
        ntot = size(KVS)*size(BOX_SIZES)*NTESTS
        cnt  = 0
        do ikv=1,size(KVS)
            do ibox=1,size(BOX_SIZES)
                ldim = [BOX_SIZES(ibox),BOX_SIZES(ibox),1]
                call img%new(ldim,SMPDS(ibox))
                tfun = ctf(SMPDS(ibox), KVS(ikv), CS, FRACA)
                do itst=1,NTESTS
                    cnt    = cnt + 1
                    call progress(cnt,ntot)
                    dfx    = ran3()*MAX_DF+0.1
                    dfy    = ran3()*MAX_DF+0.1
                    angast = ran3()*360.
                    call tfun%ctf2img(img, dfx, 'ctf', dfy, angast)
                end do
            end do
        end do
        write(*,'(a)') 'SIMPLE_CTF_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
    end subroutine test_ctf

end module simple_ctf
