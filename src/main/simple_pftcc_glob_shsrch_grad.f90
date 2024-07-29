module simple_pftcc_glob_shsrch_grad
include 'simple_lib.f08'
use simple_opt_spec,         only: opt_spec
use simple_polarft_corrcalc, only: pftcc_glob
use simple_optimizer,        only: optimizer
use simple_parameters,       only: params_glob
implicit none

public :: pftcc_shsrch_grad
private
#include "simple_local_flags.inc"

type :: pftcc_shsrch_grad
    private
    type(opt_spec)            :: ospec                   !< optimizer specification object
    class(optimizer), pointer :: opt_obj      => null()  !< optimizer object
    integer                   :: nrots        =  0       !< # rotations
    integer                   :: maxits       =  100     !< max # iterations
    integer                   :: cur_inpl_idx =  0       !< index of inplane angle for shift search
    integer                   :: max_evals    =  5       !< max # inplrot/shsrch cycles
    logical                   :: opt_angle    =  .true.  !< optimise in-plane angle with callback flag
contains
    procedure :: new         => grad_shsrch_new
    procedure :: minimize    => grad_shsrch_minimize
    procedure :: kill        => grad_shsrch_kill
end type pftcc_shsrch_grad

contains

    subroutine grad_shsrch_new( self, lims, lims_init, maxits, opt_angle )
        use simple_opt_factory, only: opt_factory
        class(pftcc_shsrch_grad),   intent(inout) :: self           !< instance
        real,                       intent(in)    :: lims(:,:)      !< limits for barrier constraint
        real,             optional, intent(in)    :: lims_init(:,:) !< limits for simplex initialisation by randomised bounds
        integer,          optional, intent(in)    :: maxits         !< maximum iterations
        logical,          optional, intent(in)    :: opt_angle      !< optimise in-plane angle with callback flag
        type(opt_factory) :: opt_fact
        call self%kill
        self%maxits = 100
        if( present(maxits) ) self%maxits = maxits
        self%opt_angle = .true.
        if( present(opt_angle) ) self%opt_angle = opt_angle
        ! make optimizer spec
        call self%ospec%specify('lbfgsb', 2, factr=1.0d+7, pgtol=1.0d-5, limits=lims,&
            max_step=0.01, limits_init=lims_init, maxits=self%maxits)
        ! generate the optimizer object
        call opt_fact%new(self%ospec, self%opt_obj)
        ! get # rotations
        self%nrots = pftcc_glob%get_nrots()
        ! set costfun pointers
        self%ospec%costfun_8    => grad_shsrch_costfun
        self%ospec%gcostfun_8   => grad_shsrch_gcostfun
        self%ospec%fdfcostfun_8 => grad_shsrch_fdfcostfun
        if( self%opt_angle ) self%ospec%opt_callback => grad_shsrch_optimize_angle_wrapper
    end subroutine grad_shsrch_new

    function grad_shsrch_costfun( self, vec, D ) result( cost )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(in)    :: vec(D)
        real(dp)                :: cost
        select type(self)
            class is (pftcc_shsrch_grad)
                ! cost = - pftcc_glob%gencorr_for_rot_8(self%reference, self%particle, vec, self%cur_inpl_idx)
            class default
                THROW_HARD('error in grad_shsrch_costfun: unknown type; grad_shsrch_costfun')
        end select
    end function grad_shsrch_costfun

    subroutine grad_shsrch_gcostfun( self, vec, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: grad(D)
        real(dp)                :: corrs_grad(2)
        grad = 0.
        select type(self)
            class is (pftcc_shsrch_grad)
                ! call pftcc_glob%gencorr_grad_only_for_rot_8(self%reference, self%particle, vec, self%cur_inpl_idx, corrs_grad)
                grad = - corrs_grad
            class default
                THROW_HARD('error in grad_shsrch_gcostfun: unknown type; grad_shsrch_gcostfun')
        end select
    end subroutine grad_shsrch_gcostfun

    subroutine grad_shsrch_fdfcostfun( self, vec, f, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: f, grad(D)
        real(dp)                :: corrs
        real(dp)                :: corrs_grad(2)
        f    = 0.
        grad = 0.
        select type(self)
            class is (pftcc_shsrch_grad)
                ! call pftcc_glob%gencorr_grad_for_rot_8(self%reference, self%particle, vec, self%cur_inpl_idx, corrs, corrs_grad)
                f    = - corrs
                grad = - corrs_grad
            class default
                THROW_HARD('error in grad_shsrch_fdfcostfun: unknown type; grad_shsrch_fdfcostfun')
        end select
    end subroutine grad_shsrch_fdfcostfun

    subroutine grad_shsrch_optimize_angle( self )
        class(pftcc_shsrch_grad), intent(inout) :: self
        real                                    :: corrs(self%nrots)
        ! call pftcc_glob%gencorrs(self%reference, self%particle, self%ospec%x, corrs, kweight=params_glob%l_kweight_rot)
        self%cur_inpl_idx = maxloc(corrs, dim=1)
    end subroutine grad_shsrch_optimize_angle

    subroutine grad_shsrch_optimize_angle_wrapper( self )
        class(*), intent(inout) :: self
        select type(self)
            class is (pftcc_shsrch_grad)
                call grad_shsrch_optimize_angle(self)
            class DEFAULT
                THROW_HARD('error in grad_shsrch_optimize_angle_wrapper: unknown type; simple_pftcc_shsrch_grad')
        end select
    end subroutine grad_shsrch_optimize_angle_wrapper

    !> minimisation routine
    function grad_shsrch_minimize( self, irot, xy, prev_sh, prob ) result( cxy )
        class(pftcc_shsrch_grad), intent(inout) :: self
        integer,                  intent(inout) :: irot
        real, optional,           intent(in)    :: xy(2)
        real, optional,           intent(in)    :: prev_sh(2)
        real, optional,           intent(in)    :: prob
        real     :: corrs(self%nrots), rotmat(2,2), cxy(3), lowest_shift(2), lowest_cost, prob_rnd
        real(dp) :: init_xy(2), lowest_cost_overall, coarse_cost, initial_cost
        integer  :: loc, i, lowest_rot, init_rot
        logical  :: found_better
        found_better = .false.
        if( present(xy) )then
            self%ospec%x   = xy
            self%ospec%x_8 = dble(xy)
        else
            self%ospec%x   = [0.,0.]
            self%ospec%x_8 = [0.d0,0.d0]
        endif
        if( self%opt_angle )then
            ! call pftcc_glob%gencorrs(self%reference, self%particle, self%ospec%x, corrs, kweight=params_glob%l_kweight_rot)
            self%cur_inpl_idx   = maxloc(corrs,dim=1)
            lowest_cost_overall = -corrs(self%cur_inpl_idx)
            initial_cost        = lowest_cost_overall
            ! using lowest_cost_overall to probabilistically random start the shift search and stick with the resulted shifts
            prob_rnd = -lowest_cost_overall
            if( present(prob) ) prob_rnd = prob
            if( (trim(params_glob%sh_rnd).eq.'yes') .and. (ran3() > prob_rnd) )then
                init_xy(1)     = 2.*(ran3()-0.5) * params_glob%sh_sig
                init_xy(2)     = 2.*(ran3()-0.5) * params_glob%sh_sig
                if( present(prev_sh) ) init_xy = init_xy - prev_sh
                self%ospec%x_8 = init_xy
                self%ospec%x   = real(init_xy)
                ! shift search / in-plane rot update
                do i = 1,self%max_evals
                    call self%opt_obj%minimize(self%ospec, self, lowest_cost)
                    ! call pftcc_glob%gencorrs(self%reference, self%particle, self%ospec%x, corrs, kweight=params_glob%l_kweight_rot)
                    loc = maxloc(corrs,dim=1)
                    if( loc == self%cur_inpl_idx ) exit
                    self%cur_inpl_idx = loc
                end do
                irot    =   self%cur_inpl_idx
                cxy(1)  = - real(lowest_cost)  ! correlation
                cxy(2:) =   self%ospec%x       ! shift
                ! rotate the shift vector to the frame of reference
                call rotmat2d(pftcc_glob%get_rot(irot), rotmat)
                cxy(2:) = matmul(cxy(2:), rotmat)
            else
                ! shift search / in-plane rot update
                do i = 1,self%max_evals
                    call self%opt_obj%minimize(self%ospec, self, lowest_cost)
                    ! call pftcc_glob%gencorrs(self%reference, self%particle, self%ospec%x, corrs, kweight=params_glob%l_kweight_rot)
                    loc = maxloc(corrs,dim=1)
                    if( loc == self%cur_inpl_idx ) exit
                    self%cur_inpl_idx = loc
                end do
                ! update best
                lowest_cost = -corrs(self%cur_inpl_idx)
                if( lowest_cost < lowest_cost_overall )then
                    found_better        = .true.
                    lowest_cost_overall = lowest_cost
                    lowest_rot          = self%cur_inpl_idx
                    lowest_shift        = self%ospec%x
                endif
                if( found_better )then
                    irot    =   lowest_rot                 ! in-plane index
                    cxy(1)  = - real(lowest_cost_overall)  ! correlation
                    cxy(2:) =   real(lowest_shift)         ! shift
                    ! rotate the shift vector to the frame of reference
                    call rotmat2d(pftcc_glob%get_rot(irot), rotmat)
                    cxy(2:) = matmul(cxy(2:), rotmat)
                else
                    irot = 0 ! to communicate that a better solution was not found
                endif
            endif
        else
            self%cur_inpl_idx   = irot
            ! lowest_cost_overall = -pftcc_glob%gencorr_for_rot_8(self%reference, self%particle, self%ospec%x_8, self%cur_inpl_idx)
            initial_cost        = lowest_cost_overall
            ! using lowest_cost_overall to probabilistically random start the shift search and stick with the resulted shifts
            prob_rnd = -lowest_cost_overall
            if( present(prob) ) prob_rnd = prob
            if( (trim(params_glob%sh_rnd).eq.'yes') .and. (ran3() > prob_rnd) )then
                init_xy(1)     = 2.*(ran3()-0.5) * params_glob%sh_sig
                init_xy(2)     = 2.*(ran3()-0.5) * params_glob%sh_sig
                if( present(prev_sh) ) init_xy = init_xy - prev_sh
                self%ospec%x_8 = init_xy
                self%ospec%x   = real(init_xy)
                call self%opt_obj%minimize(self%ospec, self, lowest_cost)
                cxy(1)  = - real(lowest_cost)  ! correlation
                cxy(2:) =   self%ospec%x       ! shift
                ! rotate the shift vector to the frame of reference
                call rotmat2d(pftcc_glob%get_rot(irot), rotmat)
                cxy(2:) = matmul(cxy(2:), rotmat)
            else
                ! shift search
                call self%opt_obj%minimize(self%ospec, self, lowest_cost)
                if( lowest_cost < lowest_cost_overall )then
                    found_better        = .true.
                    lowest_cost_overall = lowest_cost
                    lowest_shift        = self%ospec%x
                endif
                if( found_better )then
                    cxy(1)  = - real(lowest_cost_overall)  ! correlation
                    cxy(2:) =   lowest_shift               ! shift
                    ! rotate the shift vector to the frame of reference
                    call rotmat2d(pftcc_glob%get_rot(irot), rotmat)
                    cxy(2:) = matmul(cxy(2:), rotmat)
                else
                    irot = 0 ! to communicate that a better solution was not found
                endif
            endif
        end if
    end function grad_shsrch_minimize

    subroutine grad_shsrch_kill( self )
        class(pftcc_shsrch_grad), intent(inout) :: self
        if( associated(self%opt_obj) )then
            call self%ospec%kill
            call self%opt_obj%kill
            nullify(self%opt_obj)
        end if
    end subroutine grad_shsrch_kill

end module simple_pftcc_glob_shsrch_grad
