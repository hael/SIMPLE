! optimiser specification
module simple_opt_spec
include 'simple_lib.f08'
implicit none

public :: opt_spec, costfun
private

!> struct for all opt specifications
type :: opt_spec
    procedure(costfun),      pointer, nopass :: costfun      =>null() !< defines cost function
    procedure(gcostfun),     pointer, nopass :: gcostfun     =>null() !< defines the gradient of the cost function
    procedure(fdfcostfun),   pointer, nopass :: fdfcostfun   =>null() !< defines the gradient of the cost function
    procedure(costfun_8),    pointer, nopass :: costfun_8    =>null() !< defines cost function, double precision
    procedure(gcostfun_8),   pointer, nopass :: gcostfun_8   =>null() !< defines the gradient of the cost function, double precision
    procedure(fdfcostfun_8), pointer, nopass :: fdfcostfun_8 =>null() !< defines the gradient of the cost function, double precision
    procedure(opt_callback), pointer, nopass :: opt_callback =>null() !< callback function for in-between iteration steps
    character(len=STDLEN)     :: str_opt=''                           !< string descriptor (of optimization routine to be used)
    character(len=STDLEN)     :: str_mode=''                          !< mode string descriptor
    integer                   :: ndim=0, ldim=0                       !< dimension parameters
    real                      :: ftol=1e-5, gtol=1e-5                 !< fractional convergence tolerance to be achieved in costfun/gradient
    real                      :: eps=0.5                              !< learning rate
    real                      :: cfac=0.1                             !< convergence factor bfgs
    real                      :: yb                                   !< best cost obtained so far
    real                      :: max_step=0.01                        !< maximum step size
    integer                   :: maxits=100                           !< maximum number of iterations
    integer                   :: nbest=0                              !< number of best solutions used to re-estimate the model in CE
    integer                   :: niter=0                              !< actual number of iterations
    integer                   :: nrestarts=1                          !< number of restarts
    integer                   :: nevals=0                             !< total number of cost function evaluations
    integer                   :: ngevals=0                            !< total number of gradient function evaluations
    integer                   :: nstates=1                            !< number of states
    integer                   :: npeaks=0                             !< number of peaks (local optima to identify)
    integer                   :: npop=0                               !< population size (4 evolutionary optimizers)
    integer                   :: nnn=0                                !< number of nearest neighbors
    integer                   :: peakcnt=0                            !< peak counter
    integer                   :: lbfgsb_m                             !< m-parameter for lbfgsb algorithm (default: 15)
    logical, allocatable      :: cyclic(:)                            !< to indicate which variables are cyclic (Euler angles)
    real, allocatable         :: limits(:,:)                          !< variable bounds
    real, allocatable         :: limits_init(:,:)                     !< variable bounds for initialisation (randomized bounds)
    real, allocatable         :: stepsz(:)                            !< step sizes for brute force search
    real, allocatable         :: x(:)                                 !< best/final solution
    real(dp), allocatable     :: x_8(:)                               !< double precision value for best/final solution
    real, allocatable         :: xi(:)                                !< search direction used in linmin
    real, allocatable         :: xt(:)                                !< test point, used in linmin
    real, allocatable         :: inipopulation(:,:)                   !< input population for the evolutionary approaches
    real, allocatable         :: population(:,:)                      !< output solution population from the evolutionary approaches
    real, allocatable         :: peaks(:,:)                           !< output peaks (local optimal solutions)
    real, allocatable         :: grad_4_tmp(:)                        !< temporary storage for single precision gradient
    real, allocatable         :: x_4_tmp(:)                           !< temporary storage for single precision function value
    logical                   :: debug     = .false.
    logical                   :: warn      = .false.
    logical                   :: verbose   = .false.
    logical                   :: converged = .false.                  !< converged status
    logical                   :: exists    = .false.                  !< to indicate existence
  contains
    procedure :: specify
    procedure :: change_opt
    procedure :: set_limits
    procedure :: set_limits_init
    procedure :: set_costfun
    procedure :: set_gcostfun
    procedure :: set_fdfcostfun
    procedure :: set_costfun_8
    procedure :: set_gcostfun_8
    procedure :: set_fdfcostfun_8
    procedure :: set_inipop
    procedure :: eval_f_4
    procedure :: eval_df_4
    procedure :: eval_fdf_4
    procedure :: eval_f_8
    procedure :: eval_df_8
    procedure :: eval_fdf_8
    generic   :: eval_f  => eval_f_4, eval_f_8
    generic   :: eval_df => eval_df_4, eval_df_8
    generic   :: eval_fdf => eval_fdf_4, eval_fdf_8
    procedure :: kill
end type opt_spec

!>  \brief  defines the cost function interface
abstract interface
    function costfun( fun_self, vec, D ) result( cost )
        class(*), intent(inout) :: fun_self
        integer, intent(in)     :: D
        real,    intent(in)     :: vec(D)
        real                    :: cost
    end function
end interface

!>  \brief  defines the cost function gradient interface
abstract interface
    subroutine gcostfun( fun_self, vec, grad, D )
        class(*), intent(inout) :: fun_self
        integer, intent(in)     :: D
        real,    intent(inout)  :: vec(D)
        real,    intent(out)    :: grad(D)
    end subroutine
end interface

!>  \brief  defines the cost function with simultaneous gradient interface
abstract interface
    subroutine fdfcostfun( fun_self, vec, f, grad, D )
        class(*), intent(inout) :: fun_self
        integer, intent(in)     :: D
        real,    intent(inout)  :: vec(D)
        real,    intent(out)    :: f, grad(D)
    end subroutine
end interface

!>  \brief  defines the cost function interface, double precision
abstract interface
    function costfun_8( fun_self, vec, D ) result( cost )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: D
        real(kind=8), intent(in)    :: vec(D)
        real(kind=8)                :: cost
    end function
end interface

!>  \brief  defines the cost function gradient interface, double precision
abstract interface
    subroutine gcostfun_8( fun_self, vec, grad, D )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: D
        real(kind=8), intent(inout) :: vec(D)
        real(kind=8), intent(out)   :: grad(D)
    end subroutine
end interface

!>  \brief  defines the cost function with simultaneous gradient interface, double precision
abstract interface
    subroutine fdfcostfun_8( fun_self, vec, f, grad, D )
        class(*),     intent(inout) :: fun_self
        integer,      intent(in)    :: D
        real(kind=8), intent(inout) :: vec(D)
        real(kind=8), intent(out)   :: f, grad(D)
    end subroutine
end interface

!>  \brief  defines a callback routine that can be invoked in-between iteration steps
abstract interface
    subroutine opt_callback( fun_self )
        class(*), intent(inout) :: fun_self
    end subroutine opt_callback
end interface

contains

    !>  \brief  specifies the optimization parameters
    subroutine specify( self, str_opt, ndim, mode, ldim, ftol, gtol, maxits,&
    nrestarts, limits, limits_init, cyclic, verbose_arg, debug_arg, npop,&
    cfac, warn_arg, nbest, nstates, stepsz, npeaks, nnn, max_step )
        class(opt_spec),            intent(inout) :: self                !< instance
        character(len=*),           intent(in)    :: str_opt             !< string descriptor (of optimization routine to be used)
        integer,                    intent(in)    :: ndim                !< problem dimensionality
        character(len=*), optional, intent(in)    :: mode                !< mode string descriptor
        integer,          optional, intent(in)    :: ldim                !< second problem dimensionality
        real,             optional, intent(in)    :: ftol                !< fractional convergence tolerance to be achieved in costfun
        real,             optional, intent(in)    :: gtol                !< fractional convergence tolerance to be achieved in gradient
        integer,          optional, intent(in)    :: maxits              !< maximum number of iterations
        integer,          optional, intent(in)    :: nbest               !< nr of best solutions used to update the CE model
        integer,          optional, intent(in)    :: nrestarts           !< number of restarts
        integer,          optional, intent(in)    :: nstates             !< nr of states
        real,             optional, intent(in)    :: limits(ndim,2)      !< variable bounds
        real,             optional, intent(in)    :: limits_init(ndim,2) !< variable bounds for initialisation (randomized bounds)
        logical,          optional, intent(in)    :: cyclic(ndim)        !< to indicate which variables are cyclic (Euler angles)
        logical,          optional, intent(in)    :: verbose_arg         !< verbose output of optimizer
        logical,          optional, intent(in)    :: debug_arg           !< debugging mode
        logical,          optional, intent(in)    :: warn_arg            !< warning mode
        integer,          optional, intent(in)    :: npop                !< population size (4 evolutionary optimizers)
        integer,          optional, intent(in)    :: npeaks              !< number of peaks (local optima to identify)
        integer,          optional, intent(in)    :: nnn                 !< number of nearest neighbors
        real,             optional, intent(in)    :: cfac                !< convergence factor (bfgs)
        real,             optional, intent(in)    :: stepsz(ndim)        !< step sizes for brute force search
        real,             optional, intent(in)    :: max_step            !< initial step size
        integer :: alloc_stat
        call self%kill
        ! take care of optimizer specifics
        select case(str_opt)
            case('bfgs')
                ! do nothing
            case('fr_cg')
                ! do nothing
            case('pr_cg')
                ! do nothing
            case('bfgs2')
                ! do nothing
            case('lfbgsb')
                ! do nothing
            case('stde')
                ! do nothing
            case('lbfgsb')
                ! do nothing
            case('powell')
                self%maxits   = ndim*1000
            case('simplex')
                if( .not. present(limits) ) stop 'need limits (variable bounds) for simplex opt; specify; simple_opt_spec'
                self%maxits   = ndim*1000
            case('pso')
                self%npop     = 150
                self%maxits   = ndim*10000
            case('de')
                self%str_mode = 'multimodal'
                self%maxits   = ndim*1000
            case('linmin')
                ! do nothing
            case('bforce')
                if( .not. present(limits) ) stop 'need limits (variable bounds) for brute force search; specify; simple_opt_spec'
                if( .not. present(stepsz) ) stop 'need step sizes (stepsz) for brute force search; simple_opt_spec'
            case DEFAULT
                write(*,*) 'Unsupported optimizer string descriptor:', str_opt
                write(*,*) 'specify; simple_opt_spec'
        end select
        ! parse input data
        self%str_opt = str_opt
        self%ndim = ndim
        if(present(mode))        self%str_mode  = mode
        if(present(npop))        self%npop      = npop
        if(present(ldim))        self%ldim      = ldim
        if(present(ftol))        self%ftol      = ftol
        if(present(gtol))        self%gtol      = gtol
        if(present(maxits))      self%maxits    = maxits
        if(present(nbest))       self%nbest     = nbest
        if(present(npeaks))      self%npeaks    = npeaks
        if(present(nrestarts))   self%nrestarts = nrestarts
        if(present(nstates))     self%nstates   = nstates
        if(present(nnn))         self%nnn       = nnn
        if(present(max_step))    self%max_step  = max_step
        if(present(limits))      call self%set_limits(limits)
        if(present(limits_init)) call self%set_limits_init(limits_init)
        allocate( self%cyclic(self%ndim), stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk('In: specify; simple_opt_spec, cyclic',alloc_stat)
        self%cyclic = .false.
        if( present(cyclic) )then
            self%cyclic = cyclic
        endif
        if(present(stepsz))then
            allocate( self%stepsz(ndim), stat=alloc_stat )
            if(alloc_stat.ne.0)call allocchk('In: specify; simple_opt_spec, stepsz',alloc_stat)
            self%stepsz = stepsz
        endif
        if(present(verbose_arg)) self%verbose = verbose_arg
        if(present(npop))        self%npop    = npop
        if(present(cfac))        self%cfac    = cfac
        if(present(debug_arg))   self%debug   = debug_arg
        if(present(warn_arg))    self%warn    = warn_arg
        ! allocate
        if( self%npeaks > 0 )then
            allocate( self%peaks(self%npeaks,self%ndim+1), source=1., stat=alloc_stat )
            if(alloc_stat.ne.0)call allocchk('In: specify; simple_opt_spec, peaks',alloc_stat)
        endif
        select case(str_opt)
            case('linmin')
                allocate( self%x(self%ndim), self%x_8(self%ndim),self%xi(self%ndim), self%xt(self%ndim), stat=alloc_stat )
                if(alloc_stat.ne.0)call allocchk('In: specify; simple_opt_spec, linmin',alloc_stat)
                self%xi = 0.
                self%xt = 0.
            case DEFAULT
                allocate( self%x(self%ndim), self%grad_4_tmp(self%ndim), &
                    & self%x_4_tmp(self%ndim), self%x_8(self%ndim),stat=alloc_stat )
                if(alloc_stat.ne.0)call allocchk('In: specify; simple_opt_spec, DEFAULT',alloc_stat)
        end select
        self%x        = 0.
        self%lbfgsb_m = 15
        self%exists   = .true.
    end subroutine specify

    !>  \brief  to change optimizer
    subroutine change_opt( self, str_opt )
        class(opt_spec), intent(inout) :: self    !< instance
        character(len=*), intent(in)   :: str_opt !< string descriptor (of optimization routine to be used)
        self%str_opt = str_opt
    end subroutine change_opt

    !>  \brief  sets the optimizer limits
    subroutine set_limits( self, lims )
        class(opt_spec), intent(inout) :: self              !< instance
        real,               intent(in) :: lims(self%ndim,2) !< new limits
        integer  :: i
        do i=1,self%ndim
            if(lims(i,2) >= lims(i,1)) then
            else
                write(*,*) 'Bad limits for the constrained optimization!'
                write(*,*) 'set_limits; simple_opt_spec'
                write(*,*) 'Lim(',i,1,'):',lims(i,1),'Lim(',i,2,'):',lims(i,2)
                stop
            endif
        end do
        if( allocated(self%limits) ) deallocate(self%limits)
        allocate( self%limits(self%ndim,2), source=lims, stat=alloc_stat )
        if(alloc_stat.ne.0)call allocchk('In: specify; simple_opt_spec, limits',alloc_stat)
    end subroutine set_limits

    !>  \brief  sets the limits for initialisation (randomized bounds)
    subroutine set_limits_init( self, lims_init )
        class(opt_spec), intent(inout) :: self                   !< instance
        real,            intent(in)    :: lims_init(self%ndim,2) !< new limits
        integer :: i, alloc_stat
        do i=1,self%ndim
            if(lims_init(i,2) >= lims_init(i,1)) then
            else
                write(*,*) 'Bad limits for the constrained optimization!'
                write(*,*) 'set_limits_init; simple_opt_spec'
                write(*,*) 'Lim(',i,1,'):',lims_init(i,1),'Lim(',i,2,'):',lims_init(i,2)
                stop
            endif
        end do
        if( allocated(self%limits_init) ) deallocate(self%limits_init)
        allocate( self%limits_init(self%ndim,2), source=lims_init, stat=alloc_stat )
        if(alloc_stat .ne. 0)call allocchk('In: set_limits_init; simple_opt_spec', alloc_stat)
    end subroutine set_limits_init

    !>  \brief  sets the initialization population for de
    subroutine set_inipop( self, pop )
        class(opt_spec), intent(inout) :: self              !< instance
        real,            intent(in)    :: pop(:,:)
        if(self%str_opt.ne.'de')stop 'Only for use with DE; simple_opt_spec%set_inipop'
        if(size(pop,dim=1) > self%npop)stop 'non-congruent initial population 1; simple_opt_spec%set_inipop'
        if(size(pop,dim=2) .ne. self%ndim)stop 'non-congruent initial population 2; simple_opt_spec%set_inipop'
        if(allocated(self%inipopulation))deallocate(self%inipopulation)
        self%inipopulation = pop
    end subroutine set_inipop

    !> \brief  evaluate the cost function
    function eval_f_4( self, fun_self, x ) result(f)
        class(opt_spec), intent(inout) :: self
        class(*),        intent(inout) :: fun_self
        real,            intent(inout) :: x(:)
        real                           :: f
        if (.not. associated(self%costfun)) then
            write (*,*) 'error : simple_opt_spec: costfun not associated'
            return !  ‘f’ may be used uninitialized in this function
        end if
        f = self%costfun(fun_self, x, self%ndim)
        self%nevals = self%nevals  + 1
    end function eval_f_4

    !> \brief  evaluate the gradient
    subroutine eval_df_4( self, fun_self, x, grad )
        class(opt_spec), intent(inout) :: self
        class(*),        intent(inout) :: fun_self
        real,            intent(inout) :: x(:)
        real,            intent(out)   :: grad(:)
        if (.not. associated(self%gcostfun)) then
            write (*,*) 'error : simple_opt_spec: gcostfun not associated'
            return
        end if
        call self%gcostfun(fun_self, x, grad, self%ndim)
        self%ngevals = self%ngevals  + 1
    end subroutine eval_df_4

    !> \brief  evaluate the cost function and gradient simultaneously, if associated, or subsequently otherwise
    subroutine eval_fdf_4( self, fun_self, x, f, gradient )
        class(opt_spec), intent(inout) :: self
        class(*),        intent(inout) :: fun_self
        real,            intent(inout) :: x(:)
        real,            intent(out)   :: f, gradient(:)
        if ((.not. associated(self%costfun)).or.(.not. associated(self%gcostfun))) then
            write (*,*) 'error : simple_opt_spec: costfun or gcostfun not associated'
            return
        end if
        if (associated(self%fdfcostfun)) then
            call self%fdfcostfun(fun_self, x, f, gradient, self%ndim)
        else
            f        = self%costfun(fun_self, x, self%ndim)
            call self%gcostfun(fun_self, x, gradient, self%ndim)
        end if
        self%nevals  = self%nevals  + 1
        self%ngevals = self%ngevals + 1
    end subroutine eval_fdf_4

    !> \brief  evaluate the cost function, double precision
    function eval_f_8( self, fun_self, x ) result(f)
        class(opt_spec), intent(inout) :: self
        class(*),        intent(inout) :: fun_self
        real(dp),        intent(inout) :: x(:)
        real(dp)                       :: f
        if ((.not. associated(self%costfun)).and.(.not. associated(self%costfun_8))) then
            write (*,*) 'error : simple_opt_spec: no costfun associated'
            return !  ‘f’ may be used uninitialized in this function
        end if
        if (associated(self%costfun_8)) then
            f = self%costfun_8(fun_self, x, self%ndim)
        else
            self%x_4_tmp = real(x, kind=4)
            f   = self%costfun(fun_self, self%x_4_tmp, self%ndim)
        end if
        self%nevals = self%nevals  + 1
    end function eval_f_8

    !> \brief  evaluate the gradient, double precision
    subroutine eval_df_8( self, fun_self, x, grad )
        class(opt_spec), intent(inout) :: self
        class(*),        intent(inout) :: fun_self
        real(dp),        intent(inout) :: x(:)
        real(dp),        intent(out)   :: grad(:)
        if ((.not. associated(self%gcostfun)).and.(.not. associated(self%gcostfun_8))) then
            write (*,*) 'error : simple_opt_spec: no gcostfun associated'
            return
        end if
        if (associated(self%gcostfun_8)) then
            call self%gcostfun_8(fun_self, x, grad, self%ndim)
        else
            self%x_4_tmp = real(x, kind=4)
            call self%gcostfun(fun_self, self%x_4_tmp, self%grad_4_tmp, self%ndim)
            grad = self%grad_4_tmp
        end if
        self%ngevals = self%ngevals  + 1
    end subroutine eval_df_8

    !> \brief  evaluate the cost function and gradient simultaneously, if associated, or subsequently otherwise, double precision
    subroutine eval_fdf_8( self, fun_self, x, f, gradient )
        class(opt_spec), intent(inout) :: self
        class(*),        intent(inout) :: fun_self
        real(dp),        intent(inout) :: x(:)
        real(dp),        intent(out)   :: f, gradient(:)
        real                           :: f_4
        if   ( ((.not. associated(self%costfun  )) .or. (.not. associated(self%gcostfun))) .and.  &
            &  ((.not. associated(self%costfun_8)) .or. (.not. associated(self%gcostfun_8)))         ) then
            write (*,*) 'error : simple_opt_spec: no costfun or gcostfun associated'
            return
        end if

        if (associated(self%fdfcostfun_8)) then
            call self%fdfcostfun_8(fun_self, x, f, gradient, self%ndim)
        else if ( (associated(self%costfun_8)) .and. ( associated(self%gcostfun_8)) ) then
            f = self%costfun_8(fun_self, x, self%ndim)
            call self%gcostfun_8(fun_self, x, gradient, self%ndim)
        else
            self%x_4_tmp = real(x, kind=4)
            if (associated(self%fdfcostfun)) then
                call self%fdfcostfun(fun_self, self%x_4_tmp, f_4, self%grad_4_tmp, self%ndim)
                f        = f_4
                gradient = self%grad_4_tmp
            else
                f        = self%costfun(fun_self, self%x_4_tmp, self%ndim)
                call self%gcostfun(fun_self, self%x_4_tmp, self%grad_4_tmp, self%ndim)
                gradient = self%grad_4_tmp
            endif
        end if
        self%nevals  = self%nevals  + 1
        self%ngevals = self%ngevals + 1
    end subroutine eval_fdf_8

    !>  \brief  sets the cost function in the spec object
    subroutine set_costfun( self, fun )
        class(opt_spec), intent(inout) :: self !< instance
        interface
            !< defines cost function interface
            function fun( fun_self, vec, D ) result(cost)
                class(*), intent(inout) :: fun_self
                integer,  intent(in)    :: D
                real,     intent(in)    :: vec(D)
                real                    :: cost
            end function
        end interface
        self%costfun => fun
    end subroutine set_costfun

    !>  \brief  sets the gradient function in the spec object
    subroutine set_gcostfun( self, fun )
        class(opt_spec), intent(inout) :: self !< instance
        interface
            !< defines cost function gradient interface
            subroutine fun( fun_self, vec, grad, D )
                class(*), intent(inout) :: fun_self
                integer,  intent(in)    :: D
                real,     intent(inout) :: vec( D )
                real,     intent(out)   :: grad( D )
            end subroutine fun
        end interface
        self%gcostfun => fun
    end subroutine set_gcostfun

    !>  \brief  sets the cost and simultaneous gradient function in the spec object
    subroutine set_fdfcostfun( self, fun )
        class(opt_spec), intent(inout) :: self !< instance
        interface
            !< defines cost function gradient interface
            subroutine fun( fun_self, vec, f, grad, D )
                class(*), intent(inout) :: fun_self
                integer,  intent(in)    :: D
                real,     intent(inout) :: vec(D)
                real,     intent(out)   :: f, grad(D)
            end subroutine fun
        end interface
        self%fdfcostfun => fun
    end subroutine set_fdfcostfun


    !>  \brief  sets the cost function in the spec object
    subroutine set_costfun_8( self, fun )
        class(opt_spec), intent(inout) :: self !< instance
        interface
            !< defines cost function interface
            function fun( fun_self, vec, D ) result(cost)
                class(*),     intent(inout) :: fun_self
                integer,      intent(in)    :: D
                real(kind=8), intent(in)    :: vec(D)
                real(kind=8)                :: cost
            end function
        end interface
        self%costfun_8 => fun
    end subroutine set_costfun_8

    !>  \brief  sets the gradient function in the spec object
    subroutine set_gcostfun_8( self, fun )
        class(opt_spec), intent(inout) :: self !< instance
        interface
            !< defines cost function gradient interface
            subroutine fun( fun_self, vec, grad, D )
                class(*),     intent(inout) :: fun_self
                integer,      intent(in)    :: D
                real(kind=8), intent(inout) :: vec( D )
                real(kind=8), intent(out)   :: grad( D )
            end subroutine fun
        end interface
        self%gcostfun_8 => fun
    end subroutine set_gcostfun_8

    !>  \brief  sets the cost and simultaneous gradient function in the spec object
    subroutine set_fdfcostfun_8( self, fun )
        class(opt_spec), intent(inout) :: self !< instance
        interface
            !< defines cost function gradient interface
            subroutine fun( fun_self, vec, f, grad, D )
                class(*),     intent(inout) :: fun_self
                integer,      intent(in)    :: D
                real(kind=8), intent(inout) :: vec(D)
                real(kind=8), intent(out)   :: f, grad(D)
            end subroutine fun
        end interface
        self%fdfcostfun_8 => fun
    end subroutine set_fdfcostfun_8

    !>  \brief  is a destructor
    subroutine kill( self )
        class(opt_spec), intent(inout) :: self !< instance
        if( self%exists )then
            if( allocated(self%cyclic)        ) deallocate(self%cyclic)
            if( allocated(self%limits)        ) deallocate(self%limits)
            if( allocated(self%limits_init)   ) deallocate(self%limits_init)
            if( allocated(self%stepsz)        ) deallocate(self%stepsz)
            if( allocated(self%x)             ) deallocate(self%x)
            if( allocated(self%x_8)           ) deallocate(self%x_8)
            if( allocated(self%xi)            ) deallocate(self%xi)
            if( allocated(self%xt)            ) deallocate(self%xt)
            if( allocated(self%population)    ) deallocate(self%population)
            if( allocated(self%inipopulation) ) deallocate(self%inipopulation)
            if( allocated(self%peaks)         ) deallocate(self%peaks)
            if( allocated(self%grad_4_tmp)    ) deallocate(self%grad_4_tmp)
            if( allocated(self%x_4_tmp)       ) deallocate(self%x_4_tmp)
            self%costfun      => null()
            self%gcostfun     => null()
            self%fdfcostfun   => null()
            self%costfun_8    => null()
            self%gcostfun_8   => null()
            self%fdfcostfun_8 => null()
            self%opt_callback => null()
            self%exists = .false.
        endif
    end subroutine kill

end module simple_opt_spec
