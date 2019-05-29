module simple_ctf_estimate_cost
include 'simple_lib.f08'
use simple_image,    only: image
use simple_ctf,      only: ctf
use simple_opt_spec, only: opt_spec
use simple_opt_de,   only: opt_de
implicit none

public :: ctf_estimate_cost
private
#include "simple_local_flags.inc"

real, parameter :: FTOL_REFINE = 1.e-4

type ctf_estimate_cost
    private
    class(image),     pointer :: pspec              !< reference SQRT( power spectrum )
    class(ctfparams), pointer :: parms              !< For microscope characteristics
    integer,          pointer :: inds(:,:)          !< Mask indices
    type(opt_spec)            :: ospec              !< optimizer specification
    type(opt_de)              :: diffevol           !< DE optimizer
    type(ctf)                 :: tfun               !< Transfer function
    real                      :: astigtol   = 0.05  !< Astigmatic tolerance
    integer                   :: ndim       = 0     !< # of dimensions to optimize
    integer                   :: ldim(3)    = 0
    integer                   :: flims(3,2) = 0
    integer                   :: ninds      = 0
    contains
        procedure          :: init
        procedure          :: eval
        procedure          :: minimize
        procedure, private :: cost
        procedure, private :: calc_cost
        procedure          :: kill
end type ctf_estimate_cost

contains

    !>  Constructor
    subroutine init( self, pspec, parms, inds, ndim, limits, astigtol )
        class(ctf_estimate_cost), intent(inout) :: self
        class(image),     target, intent(inout) :: pspec
        class(ctfparams), target, intent(inout) :: parms
        integer,          target, intent(inout) :: inds(:,:)
        integer,                  intent(in)    :: ndim
        real,                     intent(in)    :: limits(ndim,2)
        real,                     intent(in)    :: astigtol
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
        if( self%ndim > 1 )then
            call self%ospec%specify('de',self%ndim,limits=limits, maxits=400, ftol=FTOL_REFINE)
            call self%ospec%set_costfun(cost_wrapper)
        endif
    end subroutine init

    !>  Single evaluation function for convenience and testing
    real function eval( self, dfx, dfy, angast, phshift )
        class(ctf_estimate_cost), intent(inout) :: self
        real,           intent(in) :: dfx,dfy,angast
        real, optional, intent(in) :: phshift
        eval = self%calc_cost(dfx, dfy, angast, phshift )
    end function eval

    !>  Optimization routin
    subroutine minimize( self, parms, cc )
        class(ctf_estimate_cost), intent(inout) :: self
        class(ctfparams),         intent(inout) :: parms
        real,                     intent(out)   :: cc
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
    real function cost_wrapper( self, vec, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        cost_wrapper = 1.
        select type(self)
            class is(ctf_estimate_cost)
                cost_wrapper = self%cost(D, vec)
            class DEFAULT
                THROW_HARD('Forbidden type; cost_wrapper')
        end select
    end function cost_wrapper

    !>  Cost evaluation fork
    real function cost( self, D, vec )
        class(ctf_estimate_cost), intent(inout) :: self
        integer,                  intent(in)    :: D
        real,                     intent(in)    :: vec(D)
        cost = 1.
        if( abs(vec(1) - vec(2)) > self%astigtol )then
            return
        else
            select case(D)
                case(1)
                    ! 1D defocus
                    cost = self%calc_cost(vec(1),vec(1),0.,cost)
                case(2)
                    ! 2D defocus
                    if( self%parms%l_phaseplate )then
                        cost = self%calc_cost(vec(1),vec(2),self%parms%angast,add_phshift=self%parms%phshift)
                    else
                        cost = self%calc_cost(vec(1),vec(2),self%parms%angast)
                    endif
                case(3)
                    ! defocus & astigmatism
                    if( self%parms%l_phaseplate )then
                        cost = self%calc_cost(vec(1),vec(2),vec(3),add_phshift=self%parms%phshift)
                    else
                        cost = self%calc_cost(vec(1),vec(2),vec(3))
                    endif
                case(4)
                    ! defocus,astigmatism & phase-shift
                    cost = self%calc_cost(vec(1),vec(2),vec(3),vec(4))
            end select
        endif
    end function cost

    ! cost function is real-space correlation within resolution mask between the CTF
    ! powerspectrum (the model) and the pre-processed micrograph powerspectrum (the data)
    real function calc_cost( self, dfx, dfy, angast, add_phshift )
        class(ctf_estimate_cost), intent(inout) :: self
        real,           intent(in)    :: dfx         !< defocus x-axis
        real,           intent(in)    :: dfy         !< defocus y-axis
        real,           intent(in)    :: angast      !< angle of astigmatism
        real, optional, intent(in)    :: add_phshift !< aditional phase shift (radians), for phase plate
        real, pointer :: prmat(:,:,:)
        real          :: ang, spaFreqSq, hinv, aadd_phshift, kinv, inv_ldim(3)
        real(dp)      :: ctf_sqsum,dotproduct,tvalsq,tval,corr,ctf_sum
        integer       :: h,mh,k,mk,ldim(3), i,j, l
        ! initialize, it assumes that the reference(img) ins zero mean and
        ! unit variance over the resolution range
        aadd_phshift = 0.
        if( present(add_phshift) ) aadd_phshift = add_phshift
        call self%tfun%init(dfx, dfy, angast)
        call self%pspec%get_rmat_ptr(prmat)
        mh         = maxval(self%flims(1,:))
        mk         = maxval(self%flims(2,:))
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
            tvalsq    = min(1.d0,max(tval*tval,0.0001d0))
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
        calc_cost = -real(corr,sp)
    end function calc_cost

    !>  Destructor
    subroutine kill( self )
        class(ctf_estimate_cost), intent(inout) :: self
        self%pspec => null()
        self%parms => null()
        self%inds  => null()
        call self%ospec%kill
        call self%diffevol%kill
    end subroutine kill

end module simple_ctf_estimate_cost
