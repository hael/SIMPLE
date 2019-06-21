module simple_ctf_estimate_cost
include 'simple_lib.f08'
use simple_image,    only: image
use simple_ctf,      only: ctf
use simple_opt_spec, only: opt_spec
use simple_opt_de,   only: opt_de
implicit none

public :: ctf_estimate_cost1D, ctf_estimate_cost2D
private
#include "simple_local_flags.inc"

type ctf_estimate_cost1D
    private
    class(image),     pointer :: pspec              !< reference SQRT( power spectrum )
    class(ctfparams), pointer :: parms              !< For microscope characteristics
    real,             pointer :: spec1d(:)
    type(ctf)                 :: tfun               !< Transfer function
    integer                   :: ldim(3)      = 0
    integer                   :: reslims1d(2) = 0
    contains
        procedure, private :: init1D
        procedure, private :: cost1D
        procedure, private :: kill1D
        generic            :: init => init1D
        generic            :: cost => cost1D
        generic            :: kill => kill1D
end type ctf_estimate_cost1D

type ctf_estimate_cost2D
    private
    class(image),     pointer :: pspec              !< reference SQRT( power spectrum )
    class(ctfparams), pointer :: parms              !< For microscope characteristics
    integer,          pointer :: inds(:,:)          !< Mask indices
    type(opt_spec)            :: ospec              !< optimizer specification
    type(opt_de)              :: diffevol           !< DE optimizer
    type(ctf)                 :: tfun               !< Transfer function
    real                      :: astigtol     = 0.05  !< Astigmatic tolerance
    integer                   :: ndim         = 0     !< # of dimensions to optimize
    integer                   :: ldim(3)      = 0
    integer                   :: flims(3,2)   = 0
    integer                   :: ninds        = 0
    contains
        procedure, private :: init2D
        procedure          :: minimize
        procedure, private :: cost2D
        procedure, private :: calc_cost2D
        procedure, private :: kill2D
        generic            :: init => init2D
        generic            :: kill => kill2D
end type ctf_estimate_cost2D

contains

    !>  Constructors
    subroutine init1D( self, spec_img, lims, reslims, spec,  parms )
        class(ctf_estimate_cost1D), intent(inout) :: self
        class(image),             intent(in)    :: spec_img
        integer,                  intent(in)    :: lims(2), reslims(2)
        real,             target, intent(inout) :: spec(lims(1):lims(2))
        class(ctfparams), target, intent(inout) :: parms
        call self%kill
        self%spec1d => spec
        self%parms  => parms
        self%reslims1d = reslims
        self%ldim      = spec_img%get_ldim()
        self%tfun      = ctf(self%parms%smpd, self%parms%kV, self%parms%Cs, self%parms%fraca)
    end subroutine init1D

    subroutine init2D( self, pspec, parms, inds, ndim, limits, astigtol, tol )
        class(ctf_estimate_cost2D), intent(inout) :: self
        class(image),     target, intent(inout) :: pspec
        class(ctfparams), target, intent(inout) :: parms
        integer,          target, intent(inout) :: inds(:,:)
        integer,                  intent(in)    :: ndim
        real,                     intent(in)    :: limits(ndim,2)
        real,                     intent(in)    :: astigtol, tol
        call self%kill
        self%pspec => pspec
        self%parms => parms
        self%inds  => inds
        self%ldim  = self%pspec%get_ldim()
        self%flims = self%pspec%loop_lims(3)
        self%ninds = size(self%inds,dim=2)
        self%ndim  = ndim
        self%astigtol = astigtol
        self%tfun = ctf(self%parms%smpd, self%parms%kV, self%parms%Cs, self%parms%fraca)
        call self%ospec%specify('de',self%ndim,limits=limits, maxits=400, ftol=tol)
        call self%ospec%set_costfun(cost2D_wrapper)
    end subroutine init2D

    ! 1D ROUTINE

    real function cost1D( self, df )
        class(ctf_estimate_cost1D), intent(inout) :: self
        real,                       intent(in)    :: df         !< average defocus
        real          :: ang, spaFreqSq, hinv, inv_ldim, phshift
        real(dp)      :: ctf_sqsum,dotproduct,tvalsq,tval,corr,ctf_sum
        integer       :: h, n
        ! assumes that the 1d spectrum has zero mean and
        ! unit variance over the resolution range
        phshift = 0
        if( self%parms%l_phaseplate )phshift = self%parms%phshift
        call self%tfun%init(df, df, 0.)
        n          = self%reslims1d(2)-self%reslims1d(1)+1
        inv_ldim   = 1./real(self%ldim(1))
        ctf_sqsum  = 0.d0
        dotproduct = 0.d0
        ctf_sum    = 0.d0
        do h = self%reslims1d(1),self%reslims1d(2)
            ! |CTF|
            hinv      = real(h) * inv_ldim
            spaFreqSq = hinv * hinv
            ang       = atan2(0.,real(h))
            tval      = real(self%tfun%eval(spaFreqSq, 0., phshift),dp)
            tvalsq    = min(1.d0,max(tval*tval,0.000001d0))
            tval      = dsqrt(tvalsq)
            ! correlation sums
            ctf_sum    = ctf_sum    + tval
            ctf_sqsum  = ctf_sqsum  + tvalsq
            dotproduct = dotproduct + tval*real(self%spec1d(h),dp)
        end do
        ! one pass correlation
        corr = dotproduct / dsqrt(real(n,dp)*ctf_sqsum-ctf_sum*ctf_sum)
        ! cost
        cost1D = -real(corr,sp)
    end function cost1D

    ! 2D ROUTINES

    !>  Optimization routine
    subroutine minimize( self, parms, cc )
        class(ctf_estimate_cost2D), intent(inout) :: self
        class(ctfparams),           intent(inout) :: parms
        real,                       intent(out)   :: cc
        real :: cost_here
        self%parms%dfx     = parms%dfx
        self%parms%dfy     = parms%dfy
        self%parms%angast  = parms%angast
        self%parms%phshift = parms%phshift
        self%parms%l_phaseplate = parms%l_phaseplate
        select case(self%ndim)
            case(2)
                ! defocus
                self%ospec%x = [parms%dfx,parms%dfy]
            case(3)
                ! defocus & astigmatism
                self%ospec%x = [parms%dfx,parms%dfy,parms%angast]
            case(4)
                ! defocus, astigmatism & phase shift
                if( .not.parms%l_phaseplate ) THROW_HARD('Inconsistency; minimize')
                self%ospec%x = [parms%dfx,parms%dfy,parms%angast,parms%phshift]
        end select
        call self%diffevol%new(self%ospec)
        call self%diffevol%minimize(self%ospec, self, cost_here)
        ! report
        cc = -cost_here
        parms%dfx = self%ospec%x(1)
        parms%dfy = self%ospec%x(2)
        select case(self%ndim)
            case(3)
                parms%angast = self%ospec%x(3)
            case(4)
                parms%phshift = self%ospec%x(4)
        end select
    end subroutine minimize

    !>  Wrapper for optimizer
    real function cost2D_wrapper( self, vec, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        cost2D_wrapper = 1.
        select type(self)
        class is(ctf_estimate_cost2D)
                cost2D_wrapper = self%cost2D(D, vec)
            class DEFAULT
                THROW_HARD('Forbidden type; cost_wrapper')
        end select
    end function cost2D_wrapper

    !>  Cost evaluation fork
    real function cost2D( self, D, vec )
        class(ctf_estimate_cost2D), intent(inout) :: self
        integer,                    intent(in)    :: D
        real,                       intent(in)    :: vec(D)
        cost2D = 1.
        if( abs(vec(1) - vec(2)) > self%astigtol )then
            return
        else
            select case(D)
                case(1)
                    ! 1D defocus
                    cost2D = self%calc_cost2D(vec(1),vec(1),0.,0.)
                case(2)
                    ! 2D defocus
                    if( self%parms%l_phaseplate )then
                        cost2D = self%calc_cost2D(vec(1),vec(2),self%parms%angast,add_phshift=self%parms%phshift)
                    else
                        cost2D = self%calc_cost2D(vec(1),vec(2),self%parms%angast)
                    endif
                case(3)
                    ! defocus & astigmatism
                    if( self%parms%l_phaseplate )then
                        cost2D = self%calc_cost2D(vec(1),vec(2),vec(3),add_phshift=self%parms%phshift)
                    else
                        cost2D = self%calc_cost2D(vec(1),vec(2),vec(3))
                    endif
                case(4)
                    ! defocus,astigmatism & phase-shift
                    cost2D = self%calc_cost2D(vec(1),vec(2),vec(3),vec(4))
            end select
        endif
    end function cost2D

    ! cost function is real-space correlation within resolution mask between the CTF
    ! powerspectrum (the model) and the pre-processed micrograph powerspectrum (the data)
    real function calc_cost2D( self, dfx, dfy, angast, add_phshift )
        class(ctf_estimate_cost2D), intent(inout) :: self
        real,                       intent(in)    :: dfx         !< defocus x-axis
        real,                       intent(in)    :: dfy         !< defocus y-axis
        real,                       intent(in)    :: angast      !< angle of astigmatism
        real,             optional, intent(in)    :: add_phshift !< aditional phase shift (radians), for phase plate
        real, pointer :: prmat(:,:,:)
        real          :: ang, spaFreqSq, hinv, aadd_phshift, kinv, inv_ldim(3)
        real(dp)      :: ctf_sqsum,dotproduct,tvalsq,tval,corr,ctf_sum
        integer       :: h,mh,k,mk, i,j, l
        ! initialize, it assumes that the reference(img) ins zero mean and
        ! unit variance over the resolution range
        aadd_phshift = 0.
        if( present(add_phshift) ) aadd_phshift = add_phshift
        call self%tfun%init(dfx, dfy, angast)
        call self%pspec%get_rmat_ptr(prmat)
        mh = abs(self%flims(1,1))
        mk = abs(self%flims(2,1))
        inv_ldim   = 1./real(self%ldim)
        ctf_sqsum  = 0.d0
        dotproduct = 0.d0
        ctf_sum    = 0.d0
        !$omp parallel do default(shared) private(h,hinv,k,kinv,i,j,spaFreqSq,ang,tval,tvalsq,l) &
        !$omp schedule(static) proc_bind(close) reduction(+:ctf_sqsum,dotproduct,ctf_sum)
        do l = 1,self%ninds
            h = self%inds(1,l)
            k = self%inds(2,l)
            i = min(max(1,h+mh+1),self%ldim(1))
            j = min(max(1,k+mk+1),self%ldim(2))
            ! |CTF|
            hinv      = real(h) * inv_ldim(1)
            kinv      = real(k) * inv_ldim(2)
            spaFreqSq = hinv * hinv + kinv * kinv
            ang       = atan2(real(k),real(h))
            tval      = real(self%tfun%eval(spaFreqSq, ang, aadd_phshift),dp)
            tvalsq    = min(1.d0,max(tval*tval,0.000001d0))
            tval      = dsqrt(tvalsq)
            ! correlation sums
            ctf_sum    = ctf_sum    + tval
            ctf_sqsum  = ctf_sqsum  + tvalsq
            dotproduct = dotproduct + tval*real(prmat(i,j,1),dp)
        end do
        !$omp end parallel do
        ! one pass correlation
        ! because reference is normalized (zero-mean, unit variance)
        corr = dotproduct / dsqrt(real(self%ninds,dp)*ctf_sqsum-ctf_sum*ctf_sum)
        ! cost
        calc_cost2D = -real(corr,sp)
    end function calc_cost2D

    !>  DESTRUCTORS

    subroutine kill1D( self )
        class(ctf_estimate_cost1D), intent(inout) :: self
        self%pspec => null()
        self%parms => null()
    end subroutine kill1D

    subroutine kill2D( self )
        class(ctf_estimate_cost2D), intent(inout) :: self
        self%pspec => null()
        self%parms => null()
        self%inds  => null()
        call self%ospec%kill
        call self%diffevol%kill
    end subroutine kill2D

end module simple_ctf_estimate_cost
