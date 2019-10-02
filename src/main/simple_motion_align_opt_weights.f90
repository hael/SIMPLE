! numerically find optimal weights for frames
module simple_motion_align_opt_weights
    !$ use omp_lib
    !$ use omp_lib_kinds
    include 'simple_lib.f08'
    use simple_error
    use simple_image,        only: image
    use simple_parameters,   only: params_glob
    use simple_ft_expanded,  only: ft_expanded
    use simple_ftexp_shsrch, only: ftexp_shsrch
    use simple_opt_factory,  only: opt_factory
    use simple_opt_spec,     only: opt_spec
    implicit none
    public :: motion_align_opt_weights
    private
#include "simple_local_flags.inc"

    integer, parameter :: opt_weights_Nrestarts = 3

    type :: motion_align_opt_weights
        type(image),           pointer :: frames(:) => null()
        type(ft_expanded), allocatable :: frames_ftexp(:)
        type(ft_expanded)              :: R, Rhat                 !< references for *_Ref procedures
        real,              allocatable :: weights(:)
        real,              allocatable :: Dmat(:,:)
        real                           :: hp=-1., lp=-1.
        real                           :: bfactor = -1.
        real                           :: smpd = 0.
        real                           :: scale_factor = 1.
        integer                        :: ldim(3)
        integer                        :: nframes = 0
        logical                        :: existence = .false.
    contains
        procedure :: new
        procedure :: create_ftexp_objs
        procedure :: calc_opt_weights
        procedure :: get_weights
        procedure :: dealloc_ftexp_objs
        procedure :: set_hp_lp
        procedure :: kill
        procedure :: calc_Dmat                                    !< calculate D-matrix (correlations between pairs of frames)
        procedure :: motion_align_opt_weights_cost_Dmat           !< calculate cost using the D-matrix
        procedure :: motion_align_opt_weights_gcost_Dmat          !< calculate gradient using the D-matrix
        procedure :: motion_align_opt_weights_fdf_Dmat            !< calculate cost&gradient using the D-matrix
        procedure :: motion_align_opt_weights_cost_Ref            !< calculate cost using reference (no D-matrix required)
        procedure :: motion_align_opt_weights_gcost_Ref           !< calculate gradient using reference (no D-matrix required)
        procedure :: motion_align_opt_weights_fdf_Ref             !< calculate cost&gradient using reference (no D-matrix required)
    end type motion_align_opt_weights

contains

    subroutine new( self, frames_ptr )
        class(motion_align_opt_weights),  intent(inout) :: self
        type(image), allocatable, target, intent(in)    :: frames_ptr(:)
        call self%kill
        self%frames    => frames_ptr
        self%smpd      = self%frames(1)%get_smpd()
        self%ldim      = self%frames(1)%get_ldim()
        self%hp        = min((real(minval(self%ldim(1:2))) * self%smpd)/4.,2000.)
        self%lp        = params_glob%lpstop
        self%nframes   = size(frames_ptr, 1)
        allocate(self%Dmat(self%nframes,self%nframes))
        self%existence = .true.
    end subroutine new

    subroutine calc_Dmat( self )
        class(motion_align_opt_weights), intent(inout) :: self
        integer :: i,j,k
        integer, allocatable :: pairs(:,:)
        ! generate list of pairs of frames\
        allocate(pairs(2,self%nframes * (self%nframes-1)/2))
        k = 1
        do i = 1, self%nframes
            do j = i+1, self%nframes
                pairs(1,k) = i
                pairs(2,k) = j
                k = k + 1
            end do
        end do
        do i = 1, self%nframes
            self%Dmat(i,i) = self%frames_ftexp(i)%corr_unnorm(self%frames_ftexp(i))
        end do
        !$omp parallel do default(shared) private(i,j,k) schedule(static) proc_bind(close)
        do k = 1, size(pairs,2)
            i = pairs(1,k)
            j = pairs(2,k)
            self%Dmat(i,j) = self%frames_ftexp(i)%corr_unnorm(self%frames_ftexp(j))
            self%Dmat(j,i) = self%Dmat(i,j)
        end do
        !$omp end parallel do
        deallocate(pairs)
    end subroutine calc_Dmat

    subroutine create_ftexp_objs( self )
        class(motion_align_opt_weights), intent(inout) :: self
        logical :: do_alloc_frames_ftexp
        integer :: i
        if( .not. self%existence ) then
            THROW_HARD('not instantiated; simple_motion_align_opt_weights: create_ftexp_objs')
        end if
        do_alloc_frames_ftexp = .true.
        if( allocated(self%frames_ftexp) ) then
            if( size(self%frames_ftexp) == self%nframes ) then
                do_alloc_frames_ftexp = .false.
            else
                call self%dealloc_ftexp_objs
            end if
        end if
        if( do_alloc_frames_ftexp ) allocate(self%frames_ftexp(self%nframes))
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i=1,self%nframes
            call self%frames_ftexp(i)%new(self%frames(i), self%hp, self%lp, .true., bfac=self%bfactor)
        end do
        !$omp end parallel do
        call self%R   %new(self%frames(1), self%hp, self%lp, .false., bfac=self%bfactor)
        call self%Rhat%new(self%frames(1), self%hp, self%lp, .false., bfac=self%bfactor)
        ! normalize objects
        do i=1,self%nframes
            call self%frames_ftexp(i)%normalize_mat
        end do
    end subroutine create_ftexp_objs

    subroutine dealloc_ftexp_objs( self )
        class(motion_align_opt_weights), intent(inout) :: self
        integer :: i
        if ( allocated(self%frames_ftexp) ) then
            do i = 1, size(self%frames_ftexp)
                call self%frames_ftexp(i)%kill
            end do
            deallocate(self%frames_ftexp)
        end if
        call self%R   %kill
        call self%Rhat%kill
    end subroutine dealloc_ftexp_objs

    subroutine set_hp_lp( self, hp, lp )
        class(motion_align_opt_weights), intent(inout) :: self
        real,                            intent(in)    :: hp, lp
        if( .not. self%existence ) then
            THROW_HARD('not instantiated; simple_motion_align_opt_weights: set_hp_lp')
        end if
        self%hp = hp
        self%lp = lp
        if( .not. self%nframes > 0 ) then
            THROW_HARD('nframes < 1; simple_motion_align__opt_weights: set_hp_lp')
        end if
    end subroutine set_hp_lp

    subroutine calc_opt_weights( self, w_init )
        use simple_optimizer, only: optimizer
        class(motion_align_opt_weights), intent(inout) :: self
        real, optional            :: w_init(:)
        real                      :: ww(self%nframes)
        type(opt_factory)         :: ofac
        type(opt_spec)            :: ospec
        class(optimizer), pointer :: nlopt
        real(dp) :: x_dp(self%nframes)
        real     :: x   (self%nframes)
        real     :: opt_lims(self%nframes,2)
        real     :: lowest_cost
        real     :: lowest_val
        real     :: lowest_vec(self%nframes)
        integer  :: nrun, i
        if ( .not. self%existence ) then
            THROW_HARD('not instantiated; simple_motion_align_opt_weights: calc_opt_weights')
        end if
        opt_lims(:,1) = 0.
        opt_lims(:,2) = 1.
        nlopt => null()
        if( .not. self%existence )then
            THROW_HARD('not instantiated; simple_motion_align_opt_weights: align')
        end if
        if( allocated(self%weights) ) deallocate(self%weights)
        allocate(self%weights(self%nframes))
        call ospec%specify('lbfgsb', self%nframes, &
            ftol   = 1e-7,   gtol  = 1e-7, &
            factr  = 1.0d+5, pgtol = 1.0d-7, &
            limits = opt_lims) !ftol=self%ftol, gtol=self%gtol,  default values for now
        call ofac%new(ospec, nlopt)
        call ospec%set_costfun_8(motion_align_opt_weights_cost_wrapper)
        call ospec%set_gcostfun_8(motion_align_opt_weights_gcost_wrapper)
        call ospec%set_fdfcostfun_8(motion_align_opt_weights_fdf_wrapper)
        ! minimize
        if (present(w_init)) then
            ! if w_init present, then only one minimization
            ospec%x   = w_init(1:self%nframes) / sum(w_init(1:self%nframes))
            ospec%x_8 = real(ospec%x, dp)
            call nlopt%minimize(ospec, self, lowest_cost)
            lowest_vec = ospec%x / sum(ospec%x)
        else
            ! if w_init not present, then Nrestarts minimizations with random initial conditions
            lowest_val = HUGE(lowest_val)
            do nrun = 1, opt_weights_Nrestarts
                call ran3arr(x)
                x_dp = real(x, dp)
                ospec%x   = x
                ospec%x_8 = x_dp
                call nlopt%minimize(ospec, self, lowest_cost)
                if (lowest_cost < lowest_val) then
                    lowest_val = lowest_cost
                    lowest_vec = ospec%x
                end if
            end do
        end if
        self%weights = lowest_vec / sum(ospec%x)
        ! cleanup
        call nlopt%kill
        deallocate(nlopt)
    end subroutine calc_opt_weights

    function get_weights( self ) result( res )
        class(motion_align_opt_weights), intent(inout) :: self
        real, allocatable                              :: res(:)
        if ( .not. self%existence ) then
            THROW_HARD('not instantiated; simple_motion_align_opt_weights: get_weights')
        end if
        if ( size(self%weights) .ne. self%nframes ) then
            THROW_HARD('weights not allocated properly; simple_motion_align_opt_weights: get_weights')
        end if
        allocate(res(self%nframes), source=self%weights)
    end function get_weights

    subroutine kill( self )
        class(motion_align_opt_weights), intent(inout) :: self
        call self%dealloc_ftexp_objs
        if( allocated(self%weights) ) deallocate(self%weights)
        if( allocated(self%Dmat)    ) deallocate(self%Dmat   )
        nullify(self%frames)
        self%existence = .false.
    end subroutine kill

    function motion_align_opt_weights_cost_wrapper( self, vec, D ) result( cost )
        class(*),     intent(inout) :: self
        integer,      intent(in)    :: D
        real(dp),     intent(in)    :: vec(D)
        real(dp) :: cost
        select type(self)
        class is (motion_align_opt_weights)
            cost = self%motion_align_opt_weights_cost_Dmat( vec )
            class DEFAULT
            THROW_HARD('unknown type; simple_motion_align_opt_weights: cost_wrapper')
        end select
    end function motion_align_opt_weights_cost_wrapper

    function motion_align_opt_weights_cost_sp_wrapper( self, vec, D ) result( cost )
        class(*),     intent(inout) :: self
        integer,      intent(in)    :: D
        real,         intent(in)    :: vec(D)
        real :: cost
        select type(self)
        class is (motion_align_opt_weights)
            cost = real(self%motion_align_opt_weights_cost_Dmat( real(vec, dp) ))
            class DEFAULT
            THROW_HARD('unknown type; simple_motion_align_opt_weights: cost_only_wrapper')
        end select
    end function motion_align_opt_weights_cost_sp_wrapper

    subroutine motion_align_opt_weights_gcost( self, vec, grad )
        class(motion_align_opt_weights), intent(inout) :: self
        real(dp),                        intent(in)    :: vec(:)
        real(dp),                        intent(out)   :: grad(:)
        real(dp) :: f
        call self%motion_align_opt_weights_fdf_Dmat( vec, f, grad )
    end subroutine motion_align_opt_weights_gcost

    subroutine motion_align_opt_weights_gcost_wrapper( self,  vec, grad, D )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real(dp), intent(inout) :: vec(D)
        real(dp), intent(out)   :: grad(D)
        grad = 0.d0
        select type(self)
        class is (motion_align_opt_weights)
            call self%motion_align_opt_weights_gcost_Dmat( vec, grad )
            class DEFAULT
            THROW_HARD('unknown type; simple_motion_align_opt_weights_gcost: gcost_wrapper')
        end select
    end subroutine motion_align_opt_weights_gcost_wrapper

    function motion_align_opt_weights_cost_Dmat( self, vec ) result( cost )
        class(motion_align_opt_weights), intent(inout) :: self
        real(dp),                        intent(in)    :: vec(:)
        real(dp) :: cost
        real(dp) :: Num, tmp1, tmp2
        real(dp) :: Den
        real(dp) :: ANum
        integer  :: i, j, k
        cost = 0._dp
        do i = 1, self%nframes
            Num = 0._dp
            Den = 0._dp
            do j = 1, self%nframes
                if (i == j) cycle
                Num = Num + vec(j) * self%Dmat(i,j)
                do k = 1, self%nframes
                    if (i == k) cycle
                    Den = Den + vec(j) * vec(k) * self%Dmat(j,k)
                end do
            end do
            cost = cost + Num  / sqrt(Den)
        end do
        cost = - cost
    end function motion_align_opt_weights_cost_Dmat

    subroutine motion_align_opt_weights_gcost_Dmat( self, vec, grad )
        class(motion_align_opt_weights), intent(inout) :: self
        real(dp),                        intent(in)    :: vec(:)
        real(dp),                        intent(out)   :: grad(:)
        real(dp) :: f
        call self%motion_align_opt_weights_fdf_Dmat( vec, f, grad )
    end subroutine motion_align_opt_weights_gcost_Dmat

    subroutine motion_align_opt_weights_fdf_Dmat( self, vec, f, grad )
        class(motion_align_opt_weights), intent(inout) :: self
        real(dp),                        intent(in)    :: vec(:)
        real(dp),                        intent(out)   :: f
        real(dp),                        intent(out)   :: grad(:)
        real(dp) :: Num, tmp1, tmp2
        real(dp) :: Den
        real(dp) :: ANum
        real(dp) :: BNum1, BNum2
        integer  :: i, j, k, n
        f    = 0._dp
        grad = 0._dp
        do i = 1, self%nframes
            Num = 0._dp
            Den = 0._dp
            do j = 1, self%nframes
                if (i == j) cycle
                Num = Num + vec(j) * self%Dmat(i,j)
                do k = 1, self%nframes
                    if (i == k) cycle
                    Den = Den + vec(j) * vec(k) * self%Dmat(j,k)
                end do
            end do
            f = f + Num  / sqrt(Den)
        end do
        do n = 1, self%nframes            ! loop over elements of grad
            do i = 1, self%nframes
                if (i == n) cycle         ! i neq n
                ANum = self%Dmat(i,n)
                Den = 0._dp
                do j = 1, self%nframes
                    if (j == i) cycle     ! j neq n
                    do k = 1, self%nframes
                        if (k == i) cycle ! k neq n
                        Den = Den + vec(j) * vec(k) * self%Dmat(j,k)
                    end do
                end do
                grad(n) = grad(n) + ANum / sqrt(Den)
                BNum1 = 0._dp
                do j = 1, self%nframes
                    if (j == i) cycle     ! j neq i
                    BNum1 = BNum1 + vec(j) * self%Dmat(i,j)
                end do
                BNum2 = 0._dp
                do k = 1, self%nframes
                    if (k == i) cycle     ! k neq i
                    BNum2 = BNum2 + vec(k) * self%Dmat(n,k)
                end do
                grad(n) = grad(n) - (BNum1 * BNum2) / sqrt(Den**3)
            end do
        end do
        f    = - f
        grad = - grad
    end subroutine motion_align_opt_weights_fdf_Dmat

    function motion_align_opt_weights_cost_Ref( self, vec ) result( cost )
        class(motion_align_opt_weights), intent(inout) :: self
        real(dp),                        intent(in)    :: vec(:)
        real(dp) :: cost
        integer  :: i, n
        real(dp) :: val_RR, T1, T2
        real(dp) :: val_e(self%nframes)
        real(dp) :: val_D(self%nframes)
        ! calculate reference
        call self%R%zero
        do i = 1, self%nframes
            call self%R%add_uncond(self%frames_ftexp(i), real(vec(i)))
        end do
        ! calculate e(i)
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i = 1, self%nframes
            val_e(i) = self%R%corr_unnorm(self%frames_ftexp(i))
        end do
        !$omp end parallel do
        val_RR = 0._dp
        do i = 1, self%nframes
            val_RR = val_RR + vec(i) * val_e(i)
        end do
        ! calculate denominator terms
        do i = 1, self%nframes
            val_D (i) = sqrt(val_RR - 2._dp * vec(i) * val_e(i) + vec(i)**2)
        end do
        ! calculate function value
        cost = 0._dp
        do i = 1, self%nframes
            cost = cost + (val_e(i) - vec(i)) / val_D(i)
        end do
        cost = - cost
    end function motion_align_opt_weights_cost_Ref

    subroutine motion_align_opt_weights_gcost_Ref( self, vec, grad )
        class(motion_align_opt_weights), intent(inout) :: self
        real(dp),                        intent(in)    :: vec(:)
        real(dp),                        intent(out)   :: grad(:)
        real(dp) :: f
        call self%motion_align_opt_weights_fdf_Ref( vec, f, grad )
    end subroutine motion_align_opt_weights_gcost_Ref

    subroutine motion_align_opt_weights_fdf_Ref( self, vec, f, grad )
        class(motion_align_opt_weights), intent(inout) :: self
        real(dp),                        intent(in)    :: vec(:)
        real(dp),                        intent(out)   :: f
        real(dp),                        intent(out)   :: grad(:)
        integer  :: i, n
        real(dp) :: val_RR, T1, T2
        real(dp) :: val_e(self%nframes)
        real(dp) :: val_D(self%nframes)
        real(dp) :: val_EE(self%nframes)
        real(dp) :: kk(self%nframes)
        f    = 0._dp
        grad = 0._dp
        ! calculate reference
        call self%R%zero
        do i = 1, self%nframes
            call self%R%add_uncond(self%frames_ftexp(i), real(vec(i)))
        end do
        ! calculate e(i)
        !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
        do i = 1, self%nframes
            val_e(i) = self%R%corr_unnorm(self%frames_ftexp(i))
        end do
        !$omp end parallel do
        val_RR = 0._dp
        do i = 1, self%nframes
            val_RR = val_RR + vec(i) * val_e(i)
        end do
        ! calculate denominator terms
        do i = 1, self%nframes
            val_D (i) = sqrt(val_RR - 2._dp * vec(i) * val_e(i) + vec(i)**2)
            val_EE(i) = val_D(i)**3
        end do
        ! calculate function value
        f = 0._dp
        do i = 1, self%nframes
            f = f + (val_e(i) - vec(i)) / val_D(i)
        end do
        ! calculate kk(i) (coefficients for Rhat)
        do i = 1, self%nframes
            kk(i) = ( vec(i)**2 - vec(i)*val_e(i) ) / val_EE(i)
        end do
        ! calculate Rhat
        call self%Rhat%zero
        do i = 1, self%nframes
            call self%Rhat%add_uncond(self%frames_ftexp(i), real( 1._dp/val_D(i) - kk(i)) )
        end do
        !$omp parallel do default(shared) private(n,i,T1,T2) schedule(static) proc_bind(close)
        do n = 1, self%nframes
            T1 = 0._dp          ! term 1
            do i = 1, self%nframes
                if (i==n) cycle
                T1 = T1 + val_e(n) * (val_e(i) - vec(i)) / val_EE(i)
            end do
            T2 = self%Rhat%corr_unnorm(self%frames_ftexp(n)) - 1._dp/val_D(n) + kk(n) ! term 2
            grad(n) =  T2 - T1
        end do
        !$omp end parallel do
        f    = - f
        grad = - grad
        write (*,*) 'f=', f, 'vec: ', vec / sum(vec)
    end subroutine motion_align_opt_weights_fdf_Ref

     subroutine motion_align_opt_weights_fdf_wrapper( self, vec, f, grad, D )
        class(*),     intent(inout) :: self
        integer,      intent(in)    :: D
        real(kind=8), intent(inout) :: vec(D)
        real(kind=8), intent(out)   :: f, grad(D)
        f    = 0.d0
        grad = 0.d0
        select type(self)
        class is (motion_align_opt_weights)
            call self%motion_align_opt_weights_fdf_Ref( vec, f, grad )
        class DEFAULT
            THROW_HARD('unknown type; simple_motion_align_opt_weights: fdf_wrapper')
        end select
    end subroutine motion_align_opt_weights_fdf_wrapper

end module simple_motion_align_opt_weights
