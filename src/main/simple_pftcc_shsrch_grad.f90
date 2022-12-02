! rotational origin shift alignment of band-pass limited polar projections in the Fourier domain, gradient based minimizer
module simple_pftcc_shsrch_grad
include 'simple_lib.f08'
use simple_opt_spec,         only: opt_spec
use simple_polarft_corrcalc, only: pftcc_glob
use simple_optimizer,        only: optimizer
implicit none

public :: pftcc_shsrch_grad
private
#include "simple_local_flags.inc"

logical,  parameter :: perform_rndstart = .false.
real(dp), parameter :: init_range       = 2.0_dp  ! range for random initialization (negative to positive)
integer,  parameter :: coarse_num_steps = 5       ! no. of coarse search steps in x AND y (hence real no. is its square)

type :: pftcc_shsrch_grad
    private
    type(opt_spec)            :: ospec                  !< optimizer specification object
    class(optimizer), pointer :: nlopt        =>null()  !< optimizer object
    integer                   :: reference    = 0       !< reference pft
    integer                   :: particle     = 0       !< particle pft
    integer                   :: nrots        = 0       !< # rotations
    integer                   :: maxits       = 100     !< max # iterations
    logical                   :: shbarr       = .true.  !< shift barrier constraint or not
    integer                   :: cur_inpl_idx = 0       !< index of inplane angle for shift search
    integer                   :: max_evals    = 5       !< max # inplrot/shsrch cycles
    real                      :: max_shift    = 0.      !< maximal shift
    logical                   :: opt_angle    = .true.  !< optimise in-plane angle with callback flag
    logical                   :: coarse_init  = .false. !< whether to perform an intial coarse search over the range
contains
    procedure :: new         => grad_shsrch_new
    procedure :: set_indices => grad_shsrch_set_indices
    procedure :: minimize    => grad_shsrch_minimize
    procedure :: minimize_exhaustive
    procedure :: kill        => grad_shsrch_kill
    procedure :: does_opt_angle
    procedure :: coarse_search
    procedure :: coarse_search_opt_angle
end type pftcc_shsrch_grad

contains

    subroutine grad_shsrch_new( self, lims, lims_init, shbarrier, maxits, opt_angle, coarse_init )
        use simple_opt_factory, only: opt_factory
        class(pftcc_shsrch_grad),   intent(inout) :: self           !< instance
        real,                       intent(in)    :: lims(:,:)      !< limits for barrier constraint
        real,             optional, intent(in)    :: lims_init(:,:) !< limits for simplex initialisation by randomised bounds
        character(len=*), optional, intent(in)    :: shbarrier      !< shift barrier constraint or not
        integer,          optional, intent(in)    :: maxits         !< maximum iterations
        logical,          optional, intent(in)    :: opt_angle      !< optimise in-plane angle with callback flag
        logical,          optional, intent(in)    :: coarse_init    !< coarse inital search
        type(opt_factory) :: opt_fact
        call self%kill
        ! flag the barrier constraint
        self%shbarr = .true.
        if( present(shbarrier) )then
            if( shbarrier .eq. 'no' ) self%shbarr = .false.
        endif
        self%maxits = 100
        if( present(maxits) ) self%maxits = maxits
        self%opt_angle = .true.
        if( present(opt_angle) ) self%opt_angle = opt_angle
        self%coarse_init = .false.
        if( present(coarse_init) ) self%coarse_init = coarse_init
        ! make optimizer spec
        call self%ospec%specify('lbfgsb', 2, factr=1.0d+7, pgtol=1.0d-5, limits=lims,&
            max_step=0.01, limits_init=lims_init, maxits=self%maxits)
        ! generate the optimizer object
        call opt_fact%new(self%ospec, self%nlopt)
        ! get # rotations
        self%nrots = pftcc_glob%get_nrots()
        ! set costfun pointers
        self%ospec%costfun_8    => grad_shsrch_costfun
        self%ospec%gcostfun_8   => grad_shsrch_gcostfun
        self%ospec%fdfcostfun_8 => grad_shsrch_fdfcostfun
        if( self%opt_angle ) self%ospec%opt_callback => grad_shsrch_optimize_angle_wrapper
    end subroutine grad_shsrch_new

    pure logical function does_opt_angle( self )
        class(pftcc_shsrch_grad), intent(in) :: self
        does_opt_angle = self%opt_angle
    end function does_opt_angle

    function grad_shsrch_costfun( self, vec, D ) result( cost )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(in)    :: vec(D)
        real(dp)                :: cost
        select type(self)
            class is (pftcc_shsrch_grad)
                cost = - pftcc_glob%gencorr_for_rot_8(self%reference, self%particle, vec, self%cur_inpl_idx)
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
                call pftcc_glob%gencorr_grad_only_for_rot_8(self%reference, self%particle, vec, self%cur_inpl_idx, corrs_grad)
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
                call pftcc_glob%gencorr_grad_for_rot_8(self%reference, self%particle, vec, self%cur_inpl_idx, corrs, corrs_grad)
                f    = - corrs
                grad = - corrs_grad
            class default
                THROW_HARD('error in grad_shsrch_fdfcostfun: unknown type; grad_shsrch_fdfcostfun')
        end select
    end subroutine grad_shsrch_fdfcostfun

    subroutine grad_shsrch_optimize_angle( self )
        class(pftcc_shsrch_grad), intent(inout) :: self
        real                                    :: corrs(self%nrots)
        call pftcc_glob%gencorrs(self%reference, self%particle, self%ospec%x, corrs)
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

    !> set indicies for shift search
    subroutine grad_shsrch_set_indices( self, ref, ptcl )
        class(pftcc_shsrch_grad), intent(inout) :: self
        integer,                  intent(in)    :: ref, ptcl
        self%reference = ref
        self%particle  = ptcl
    end subroutine grad_shsrch_set_indices

    !> exhaustive in-plane minimisation routine
    function minimize_exhaustive( self, irot, halfrot ) result( cxy )
        class(pftcc_shsrch_grad), intent(inout) :: self
        integer,                  intent(out)   :: irot
        real, optional,           intent(in)    :: halfrot ! within +/- x degrees of irot
        logical  :: rotmask(self%nrots)
        real     :: corrs(self%nrots), rotmat(2,2), cxy(3), shift(2), coarse_shift(2)
        real     :: cc, optcc, dang, lowest_cost, lowest_cost_overall, sx,sy, stepx,stepy
        integer  :: loc, i, j, nstepx, nstepy
        if( irot <= 0 ) irot = pftcc_glob%get_roind(0.)
        ! rotation mask
        rotmask = .true.
        if( present(halfrot) )then
            if(halfrot >= 180.0 .or. halfrot < 0. .or. halfrot < abs(pftcc_glob%get_rot(2)-pftcc_glob%get_rot(1)) )then
                THROW_HARD('Invalid rotation boundary')
            endif
            do i = 1,self%nrots
                dang = abs(pftcc_glob%get_rot(i)-pftcc_glob%get_rot(irot))
                if( dang > 180. ) dang = abs(dang-360.)
                rotmask(i) = dang <= halfrot
            enddo
        endif
        ! Coarse search
        optcc = -huge(optcc)
        nstepx = max(1,floor((self%ospec%limits(1,2)-self%ospec%limits(1,1))/4.))
        nstepy = max(1,floor((self%ospec%limits(2,2)-self%ospec%limits(2,1))/4.))
        stepx  = (self%ospec%limits(1,2)-self%ospec%limits(1,1))/real(2*nstepx+1)
        stepy  = (self%ospec%limits(2,2)-self%ospec%limits(2,1))/real(2*nstepy+1)
        do i = -nstepx,nstepx
            sx = real(i)*stepx
            do j = -nstepy,nstepy
                sy = real(j)*stepy
                call pftcc_glob%gencorrs(self%reference, self%particle, [sx,sy], corrs)
                loc = maxloc(corrs,dim=1,mask=rotmask)
                cc  = corrs(loc)
                if (cc > optcc) then
                    irot         = loc
                    optcc        = cc
                    coarse_shift = [sx,sy]
                end if
            enddo
        enddo
        shift = coarse_shift
        do i = -1,1
            sx = coarse_shift(1) + real(i)*stepx/2.
            do j = -1,1
                if( i==0 .and. j==0 )cycle
                sy = coarse_shift(2) + real(j)*stepy/2.
                call pftcc_glob%gencorrs(self%reference, self%particle, [sx,sy], corrs)
                loc = maxloc(corrs,dim=1,mask=rotmask)
                cc  = corrs(loc)
                if (cc > optcc) then
                    irot  = loc
                    optcc = cc
                    shift = [sx,sy]
                end if
            enddo
        enddo
        coarse_shift = shift
        do i = -1,1
            sx = coarse_shift(1) + real(i)*stepx/4.
            do j = -1,1
                if( i==0 .and. j==0 )cycle
                sy = coarse_shift(2) + real(j)*stepy/4.
                call pftcc_glob%gencorrs(self%reference, self%particle, [sx,sy], corrs)
                loc = maxloc(corrs,dim=1,mask=rotmask)
                cc  = corrs(loc)
                if (cc > optcc) then
                    irot  = loc
                    optcc = cc
                    shift = [sx,sy]
                end if
            enddo
        enddo
        coarse_shift      = shift
        self%cur_inpl_idx = irot
        self%ospec%x      = coarse_shift
        self%ospec%x_8    = real(self%ospec%x,dp)
        ! refinement
        lowest_cost_overall = -pftcc_glob%gencorr_for_rot(self%reference, self%particle, self%ospec%x, self%cur_inpl_idx)
        call self%nlopt%minimize(self%ospec, self, lowest_cost)
        if( lowest_cost < lowest_cost_overall )then
            ! refinement solution
            lowest_cost_overall = lowest_cost
            cxy(2:3) = self%ospec%x
        else
            ! coarse solution
            cxy(2:3) = coarse_shift
        endif
        cxy(1) = -lowest_cost_overall
        ! rotate the shift vector to the frame of reference
        call rotmat2d(pftcc_glob%get_rot(irot), rotmat)
        cxy(2:3) = matmul(cxy(2:3), rotmat)
    end function minimize_exhaustive

    !> minimisation routine
    function grad_shsrch_minimize( self, irot ) result( cxy )
        class(pftcc_shsrch_grad), intent(inout) :: self
        integer,                  intent(inout) :: irot
        real     :: corrs(self%nrots), rotmat(2,2), cxy(3), lowest_shift(2), lowest_cost
        real(dp) :: init_xy(2), lowest_cost_overall, coarse_cost
        integer  :: loc, i, lowest_rot, init_rot
        logical  :: found_better
        found_better      = .false.
        if( self%opt_angle )then
            self%ospec%x   = [0.,0.]
            self%ospec%x_8 = [0.d0,0.d0]
            call pftcc_glob%gencorrs(self%reference, self%particle, self%ospec%x, corrs)
            self%cur_inpl_idx   = maxloc(corrs,dim=1)
            lowest_cost_overall = -corrs(self%cur_inpl_idx)
            if( self%coarse_init )then
                call self%coarse_search_opt_angle(init_xy, init_rot)
                if( init_rot /= 0 )then
                    self%ospec%x_8    = init_xy
                    self%ospec%x      = real(init_xy)
                    self%cur_inpl_idx = init_rot
                endif
            end if
            ! shift search / in-plane rot update
            do i = 1,self%max_evals
                call self%nlopt%minimize(self%ospec, self, lowest_cost)
                call pftcc_glob%gencorrs(self%reference, self%particle, self%ospec%x, corrs)
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
        else
            self%cur_inpl_idx = irot
            self%ospec%x      = [0.,0.]
            self%ospec%x_8    = [0.d0,0.d0]
            lowest_cost_overall = -pftcc_glob%gencorr_for_rot_8(self%reference, self%particle, self%ospec%x_8, self%cur_inpl_idx)
            if( perform_rndstart )then
                init_xy(1)     = 2.*(ran3()-0.5) * init_range
                init_xy(2)     = 2.*(ran3()-0.5) * init_range
                self%ospec%x_8 = init_xy
                self%ospec%x   = real(init_xy)
            endif
            if( self%coarse_init )then
                call self%coarse_search(coarse_cost, init_xy)
                if( coarse_cost < lowest_cost_overall )then
                    lowest_cost_overall = coarse_cost
                    self%ospec%x_8      = init_xy
                    self%ospec%x        = real(init_xy)
                endif
            end if
            ! shift search
            call self%nlopt%minimize(self%ospec, self, lowest_cost)
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
        end if
    end function grad_shsrch_minimize

    subroutine coarse_search(self, lowest_cost, init_xy)
        class(pftcc_shsrch_grad), intent(inout) :: self
        real(dp),                 intent(out)   :: lowest_cost, init_xy(2)
        real(dp) :: x, y, cost, stepx, stepy
        integer  :: ix, iy
        lowest_cost = huge(lowest_cost)
        init_xy     = 0.d0
        if (coarse_num_steps .le. 1) return
        stepx = real(self%ospec%limits(1,2)-self%ospec%limits(1,1),dp)/real(coarse_num_steps,dp)
        stepy = real(self%ospec%limits(2,2)-self%ospec%limits(2,1),dp)/real(coarse_num_steps,dp)
        do ix = 1,coarse_num_steps
            x = self%ospec%limits(1,1)+stepx/2. + real(ix-1,dp)*stepx
            do iy = 1,coarse_num_steps
                y = self%ospec%limits(2,1)+stepy/2. + real(iy-1,dp)*stepy
                cost = -pftcc_glob%gencorr_for_rot_8(self%reference, self%particle, [x,y], self%cur_inpl_idx)
                if (cost < lowest_cost) then
                    lowest_cost = cost
                    init_xy     = [x,y]
                end if
            enddo
        enddo
    end subroutine coarse_search

    subroutine coarse_search_opt_angle(self, init_xy, irot)
        class(pftcc_shsrch_grad), intent(inout) :: self
        real(dp),                 intent(out)   :: init_xy(2)
        integer,                  intent(out)   :: irot
        real(dp) :: x, y, lowest_cost, cost, stepx,stepy
        real     :: corrs(self%nrots)
        integer  :: loc, ix,iy
        init_xy = 0.d0
        irot    = 0
        if (coarse_num_steps .le. 1) return
        lowest_cost = huge(lowest_cost)
        stepx = real(self%ospec%limits(1,2)-self%ospec%limits(1,1),dp)/real(coarse_num_steps,dp)
        stepy = real(self%ospec%limits(2,2)-self%ospec%limits(2,1),dp)/real(coarse_num_steps,dp)
        do ix = 1,coarse_num_steps
            x = self%ospec%limits(1,1)+stepx/2. + real(ix-1,dp)*stepx
            do iy = 1,coarse_num_steps
                y = self%ospec%limits(2,1)+stepy/2. + real(iy-1,dp)*stepy
                call pftcc_glob%gencorrs(self%reference, self%particle, real([x,y]), corrs)
                loc  = maxloc(corrs,dim=1)
                cost = - corrs(loc)
                if (cost < lowest_cost) then
                    lowest_cost = cost
                    irot        = loc
                    init_xy(1)  = x
                    init_xy(2)  = y
                end if
            end do
        end do
    end subroutine coarse_search_opt_angle

    subroutine grad_shsrch_kill( self )
        class(pftcc_shsrch_grad), intent(inout) :: self
        if( associated(self%nlopt) )then
            call self%ospec%kill
            call self%nlopt%kill
            nullify(self%nlopt)
        end if
    end subroutine grad_shsrch_kill

    function grad_shsrch_get_nevals( self ) result( nevals )
        class(pftcc_shsrch_grad), intent(inout) :: self
        integer :: nevals
        nevals = self%ospec%nevals
    end function grad_shsrch_get_nevals

    subroutine grad_shsrch_get_peaks( self, peaks )
        class(pftcc_shsrch_grad), intent(inout) :: self
        real, allocatable,        intent(out)   :: peaks(:,:) !< output peak matrix
        allocate(peaks(1,2))
        peaks(1,:) = self%ospec%x
    end subroutine grad_shsrch_get_peaks

end module simple_pftcc_shsrch_grad
