!------------------------------------------------------------------------------!
! SIMPLE v2.5         Elmlund & Elmlund Lab          simplecryoem.com          !
!------------------------------------------------------------------------------!
!> Simple optimisation module: core specification module
module simple_opt_spec
use simple_defs
implicit none

public :: opt_spec
private

type :: opt_spec
    procedure(costfun),  pointer, nopass :: costfun  =>null() !< defines cost function
    procedure(gcostfun), pointer, nopass :: gcostfun =>null() !< defines the gradient of the cost function
    character(len=STDLEN) :: str_opt=''                       !< string descriptor (of optimization routine to be used)
    character(len=STDLEN) :: str_mode=''                      !< mode string descriptor
    integer               :: ndim=0, ldim=0                   !< dimension parameters
    real                  :: ftol=1e-5, gtol=1e-5             !< fractional convergence tolerance to be achieved in costfun/gradient
    real                  :: eps=0.5                          !< learning rate
    real                  :: cfac=0.1                         !< convergence factor bfgs
    real                  :: yb                               !< best cost obtained so far
    integer               :: maxits=100                       !< maximum number of iterations
    integer               :: nbest=0                          !< number of best solutions used to re-estimate the model in CE
    integer               :: niter=0                          !< actual number of iterations
    integer               :: nrestarts=1                      !< number of restarts
    integer               :: nsample=0                        !< number of samples per OASIS iteration
    integer               :: nevals=0                         !< total number of cost function evaluations
    integer               :: nstates=1                        !< number of states
    integer               :: npeaks=0                         !< number of peaks (local optima to identify)
    integer               :: npop=0                           !< population size (4 evolutionary optimizers)
    integer               :: nnn=0                            !< number of nearest neighbors
    integer               :: peakcnt=0                        !< peak counter
    logical, allocatable  :: cyclic(:)                        !< to indicate which variables are cyclic (Euler angles)
    real, allocatable     :: limits(:,:)                      !< variable bounds
    real, allocatable     :: stepsz(:)                        !< step sizes for brute force search
    real, allocatable     :: x(:)                             !< best/final solution
    real, allocatable     :: xi(:)                            !< search direction used in linmin
    real, allocatable     :: xt(:)                            !< test point, used in linmin
    real, allocatable     :: sdevs(:)                         !< standard deviations of the variables, used in oasis
    real, allocatable     :: population(:,:)                  !< output solution population from the evolutionary approaches
    real, allocatable     :: peaks(:,:)                       !< output peaks (local optimal solutions)
#include "simple_local_flags.inc"
!    logical               :: verbose   = .false.              !< verbose output of optimizer on/off
!    logical               :: opt_spec_debug = .false.         !< debugging mode on/off unique to opt_spec
!    logical               :: warn      = .false.              !< warning mode on/off
    logical               :: converged = .false.              !< converged status
    logical               :: exists    = .false.              !< to indicate existence
  contains
    procedure :: specify
    procedure :: change_opt
    procedure :: set_limits
    procedure :: set_costfun
    procedure :: set_gcostfun
    procedure :: kill
end type opt_spec

!>  \brief  defines the cost function interface
abstract interface
    function costfun( vec, D ) result( cost )
        integer, intent(in) :: D
        real,    intent(in) :: vec(D)
        real                :: cost
    end function
end interface

!>  \brief  defines the cost function gradient interface
abstract interface
    function gcostfun( vec, D ) result( grad )
        integer, intent(in)    :: D
        real,    intent(inout) :: vec(D)
        real                   :: grad(D)
    end function
end interface

contains

    !>  \brief  specifies the optimization parameters
    subroutine specify( self, str_opt, ndim, mode, ldim, ftol, gtol, maxits,&
    nrestarts, nsample, limits, cyclic, verbose_arg, debug_arg, npop,&
    eps, cfac, warn_arg, nbest, nstates, stepsz, npeaks, nnn )
        use simple_jiffys, only: alloc_err
        class(opt_spec),            intent(inout) :: self           !< instance
        character(len=*),           intent(in)    :: str_opt        !< string descriptor (of optimization routine to be used)
        integer,                    intent(in)    :: ndim           !< problem dimensionality
        character(len=*), optional, intent(in)    :: mode           !< mode string descriptor
        integer,          optional, intent(in)    :: ldim           !< second problem dimensionality
        real,             optional, intent(in)    :: ftol           !< fractional convergence tolerance to be achieved in costfun
        real,             optional, intent(in)    :: gtol           !< fractional convergence tolerance to be achieved in gradient
        integer,          optional, intent(in)    :: maxits         !< maximum number of iterations
        integer,          optional, intent(in)    :: nbest          !< nr of best solutions used to update the CE model
        integer,          optional, intent(in)    :: nrestarts      !< number of restarts
        integer,          optional, intent(in)    :: nsample        !< number of samples per OASIS iteration
        integer,          optional, intent(in)    :: nstates        !< nr of states
        real,             optional, intent(in)    :: limits(ndim,2) !< variable bounds
        logical,          optional, intent(in)    :: cyclic(ndim)   !< to indicate which variables are cyclic (Euler angles)
        logical,          optional, intent(in)    :: verbose_arg    !< verbose output of optimizer
        logical,          optional, intent(in)    :: debug_arg      !< debugging mode
        logical,          optional, intent(in)    :: warn_arg       !< warning mode
        integer,          optional, intent(in)    :: npop           !< population size (4 evolutionary optimizers)
        integer,          optional, intent(in)    :: npeaks         !< number of peaks (local optima to identify)
        integer,          optional, intent(in)    :: nnn            !< number of nearest neighbors
        real,             optional, intent(in)    :: cfac           !< convergence factor (bfgs)
        real,             optional, intent(in)    :: eps            !< learning rate (OASIS)
        real,             optional, intent(in)    :: stepsz(ndim)   !< step sizes for brute force search
        integer :: alloc_stat, i
        call self%kill
        ! take care of optimizer specifics
        select case(str_opt)
            case('bfgs')
               ! do nothing
               ! self%maxits = ndim*1000
            case('powell')
                self%maxits = ndim*1000
            case('simplex')
                if( .not. present(limits) ) stop 'need limits (variable bounds) for simplex opt; specify; simple_opt_spec'
                self%maxits = ndim*1000
            case('oasis')
                if( .not. present(limits) ) stop 'need limits (variable bounds) for oasis opt; specify; simple_opt_spec'
                self%maxits  = ndim*1000
                self%nsample = max(5,ndim)
            case('pso')
                self%npop    = 150
                self%maxits  = ndim*10000
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
        if(present(mode))      self%str_mode  = mode
        if(present(npop))      self%npop      = npop
        if(present(ldim))      self%ldim      = ldim
        if(present(ftol))      self%ftol      = ftol
        if(present(gtol))      self%gtol      = gtol
        if(present(maxits))    self%maxits    = maxits
        if(present(nbest))     self%nbest     = nbest
        if(present(nsample))   self%nsample   = nsample
        if(present(npeaks))    self%npeaks    = npeaks
        if(present(nrestarts)) self%nrestarts = nrestarts
        if(present(nstates))   self%nstates   = nstates
        if(present(nnn))       self%nnn       = nnn
        if(present(limits))    call self%set_limits( limits )
        allocate( self%cyclic(self%ndim), stat=alloc_stat )
        call alloc_err('In: specify; simple_opt_spec, cyclic', alloc_stat)
        self%cyclic = .false.
        if( present(cyclic) )then
            self%cyclic = cyclic
        endif
        if(present(stepsz))then
            allocate( self%stepsz(ndim), stat=alloc_stat )
            call alloc_err('In: specify; simple_opt_spec, stepsz', alloc_stat)
            self%stepsz = stepsz
        endif
        if(present(verbose_arg)) self%verbose = verbose_arg
        if(present(npop))    self%npop    = npop
        if(present(cfac))    self%cfac    = cfac
        if(present(eps))     self%eps     = eps
        if(present(debug_arg))   self%debug   = debug_arg
        if(present(warn_arg))    self%warn    = warn_arg
        ! allocate
        if( self%npeaks > 0 )then
            allocate( self%peaks(self%npeaks,self%ndim+1), stat=alloc_stat )
            call alloc_err('In: specify; simple_opt_spec, peaks', alloc_stat)
        endif
        select case(str_opt)
            case('oasis')
                allocate( self%x(self%ndim), self%sdevs(self%ndim), stat=alloc_stat )
                call alloc_err('In: specify; simple_opt_spec, oasis', alloc_stat)
                self%sdevs  = 0.
            case('linmin')
                allocate( self%x(self%ndim), self%xi(self%ndim), self%xt(self%ndim), stat=alloc_stat )
                call alloc_err('In: specify; simple_opt_spec, linmin', alloc_stat)
                self%xi = 0.
                self%xt = 0.
            case DEFAULT
                allocate( self%x(self%ndim), stat=alloc_stat )
                call alloc_err('In: specify; simple_opt_spec, DEFAULT', alloc_stat)
        end select
        self%x  = 0.
        self%exists = .true.
    end subroutine specify

    !>  \brief  to change optimizer
    subroutine change_opt( self, str_opt )
        class(opt_spec), intent(inout) :: self    !< instance
        character(len=*), intent(in)   :: str_opt !< string descriptor (of optimization routine to be used)
        self%str_opt = str_opt
    end subroutine change_opt

    !>  \brief  sets the optimizer limits
    subroutine set_limits( self, lims )
        use simple_jiffys, only: alloc_err
        class(opt_spec), intent(inout) :: self              !< instance
        real,               intent(in) :: lims(self%ndim,2) !< new limits
        integer  :: i,alloc_stat
        do i=1,self%ndim
            if(lims(i,2) >= lims(i,1)) then
            else
                write(*,*) 'Bad limits for the constrained optimization!'
                write(*,*) 'set_limits; simple_opt_spec'
                write(*,*) 'Lim(',i,1,'):',lims(i,1),'Lim(',i,2,'):',lims(i,2)
                stop
            endif
        end do
        if( .not.allocated(self%limits))then
            allocate( self%limits(self%ndim,2), stat=alloc_stat )
            call alloc_err('In: specify; simple_opt_spec, limits', alloc_stat)
        endif
        self%limits = lims
    end subroutine set_limits

    !>  \brief  sets the cost function in the spec object
    subroutine set_costfun( self, fun )
        class(opt_spec), intent(inout) :: self !< instance
#if defined (PGI)
        ! GNU COMPILER DOES NOT COPE W EXTERNAL
        real, external :: fun !< defines cost function interface
#else
        ! PGI COMPILER DOES NOT COPE W INTERFACE
        interface
            !< defines cost function interface
            function fun( vec, D ) result(cost)
                integer, intent(in) :: D
                real,    intent(in) :: vec(D)
                real :: cost
            end function
        end interface
#endif
        self%costfun => fun
    end subroutine set_costfun

    !>  \brief  sets the gradient function in the spec object
    subroutine set_gcostfun( self, fun )
        class(opt_spec), intent(inout) :: self !< instance
#if defined (PGI)
        ! GNU COMPILER DOES NOT COPE W EXTERNAL
        real, external  :: fun !< defines cost function gradient interface
#else
        ! PGI COMPILER DOES NOT COPE W INTERFACE
        interface
            !< defines cost function gradient interface
            function fun( vec, D ) result( grad )
                integer, intent(in)    :: D
                real,    intent(inout) :: vec( D )
                real :: grad( D )
            end function
        end interface
#endif
        self%gcostfun => fun
    end subroutine set_gcostfun

    !>  \brief  is a destructor
    subroutine kill( self )
        class(opt_spec), intent(inout) :: self !< instance
        if( self%exists )then
            if( allocated(self%cyclic) )      deallocate(self%cyclic)
            if( allocated(self%limits) )      deallocate(self%limits)
            if( allocated(self%stepsz) )      deallocate(self%stepsz)
            if( allocated(self%x) )           deallocate(self%x)
            if( allocated(self%xi) )          deallocate(self%xi)
            if( allocated(self%xt) )          deallocate(self%xt)
            if( allocated(self%sdevs) )       deallocate(self%sdevs)
            if( allocated(self%population) )  deallocate(self%population)
            if( allocated(self%peaks) )       deallocate(self%peaks)
            self%costfun  => null()
            self%gcostfun => null()
            self%exists = .false.
        endif
    end subroutine kill

end module simple_opt_spec
