module simple_prime_srch
use simple_defs    ! singleton
use simple_jiffys  ! singleton
use simple_params, only: params
implicit none

type prime_srch
    private
    class(params), pointer :: pp => null()     !< pointer to parameters singleton
    integer                :: nrefs = 0        !< total number of references (nstates*nprojs)
    integer                :: nrots = 0        !< number of in-plane rotations in polar representation    
    integer, allocatable   :: inplmat(:,:)     !< in-plane rotation index matrix for GPU-based search
    real,    allocatable   :: angtab(:)        !< list of in-plane angles
    real,    allocatable   :: corrmat2d(:,:)   !< 2D correlation matrix for GPU-based search
    logical                :: exists = .false. !< 2 indicate existence
  contains
    ! CONSTRUCTOR 
    procedure :: new
    ! GETTERS
    procedure :: roind
    procedure :: rot
    procedure :: corr
    ! procedure :: get_corrmat
    procedure :: inpl
    ! procedure :: get_inplmat
    ! GPU CORR CALCULATOR
    procedure :: calc_corrs
    ! DESTRUCTOR
    procedure :: kill
end type prime_srch

interface prime_srch
    module procedure constructor
end interface prime_srch

contains
    
    !>  \brief  is a constructor
    function constructor( p, nrefs, nrots ) result( self )
        class(params), target, intent(in) :: p     !< parameters
        integer,               intent(in) :: nrefs !< number of references
        integer,               intent(in) :: nrots !< number of inplane rotations
        type(prime_srch) :: self
        call self%new(p, nrefs, nrots)
    end function constructor
    
    !>  \brief  is a constructor
    subroutine new( self, p, nrefs, nrots )
        use simple_math, only: round2even, rad2deg
        class(prime_srch),     intent(inout) :: self  !< instance
        class(params), target, intent(in)    :: p     !< parameters
        integer,               intent(in)    :: nrefs !< number of references
        integer,               intent(in)    :: nrots !< number of inplane rotations
        integer :: i, alloc_stat
        real    :: dang
        call self%kill
        self%pp    => p
        self%nrefs = nrefs
        self%nrots = nrots
        ! generate the array of in-plane angles
        allocate(self%angtab(self%nrots), stat=alloc_stat)
        call alloc_err('In: new; simple_prime_srch', alloc_stat)
        dang = twopi/real(self%nrots)
        forall( i=1:self%nrots ) self%angtab(i) = rad2deg( real(i-1)*dang )
        self%exists = .true.
    end subroutine new
    
    !>  \brief calculates the in-plane rotational index for the rot in-plane angle
    function roind( self, rot ) result( ind )
        class(prime_srch), intent(in) :: self
        real,              intent(in) :: rot
        integer :: ind, loc(1)
        loc = minloc((self%angtab-rot)**2 ) ! minimum of squared angular distances
        ind = loc(1)
    end function roind
    
    !>  \brief calculates the in-plane angle of in-plane rotational index rot
    function rot( self, ind ) result( r )
        class(prime_srch), intent(in) :: self
        integer,           intent(in) :: ind
        real :: r
        if( ind < 1 .or. ind > self%nrots )then
            write(*,*) 'rotational index is: ', ind, ' which is out of range; calc_rot; simple_prime3D_srch'
            stop
        endif
        r = self%angtab(ind)
        if( r == 360. ) r = 0.
    end function rot
    
    !>  \brief correlation getter
    real function corr( self, cnt_glob, iref )
        class(prime_srch), intent(in) :: self
        integer,           intent(in) :: cnt_glob, iref
        if( allocated(self%corrmat2d) )then
            if( cnt_glob < 1 .or. cnt_glob > self%nrefs )&
            stop 'cnt_glob out of range; simple_prime_srch :: corr'
            if( iref  < 1 .or. iref > self%nrefs )&
            stop 'iref out of range; simple_prime_srch :: corr'
            corr = self%corrmat2d(cnt_glob,iref)
        else
            stop 'corrmat2d not allocated (and hence not calculated); simple_srch :: corr'
        endif
    end function corr
    
    !>  \brief in-plane index getter
    integer function inpl( self, cnt_glob, iref )
        class(prime_srch), intent(in) :: self
        integer,           intent(in) :: cnt_glob, iref
        if( allocated(self%inplmat) )then
            if( cnt_glob < 1 .or. cnt_glob > self%nrefs )&
            stop 'cnt_glob out of range; simple_prime_srch :: inpl'
            if( iref < 1 .or. iref > self%nrefs )&
            stop 'iref out of range; simple_prime_srch :: inpl'
            inpl = self%inplmat(cnt_glob,iref)
        else
            stop 'inplmat not allocated (and hence not calculated); simple_srch :: inpl'
        endif   
    end function inpl
    
    !>  \brief  prepares the matrices for PRIME2D/3D search
    subroutine calc_corrs( self, pftcc, refine, mode, prev_corrs )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        use simple_oris,             only: oris
        use simple_corrmat,          only: project_corrmat3D_greedy, project_corrmat3D_shc
        class(prime_srch),           intent(inout) :: self          !< instance
        class(polarft_corrcalc),     intent(inout) :: pftcc         !< polarft corr calculator 
        character(len=*),            intent(in)    :: refine        !< refine flag
        character(len=*), optional,  intent(in)    :: mode          !< mode of execution
        real,             optional,  intent(in)    :: prev_corrs(:) !< previous corrs (needed for shc-based projection of the 3D corrmat)
        real, allocatable  :: corrmat3d(:,:,:) ! 3D correlation matrix for GPU-based search
        integer            :: alloc_stat
        character(len=32)  :: mmode
        mmode = 'gpu'
        if( present(mode) ) mmode = mode
        ! the resulting corrmat2d and inplmat 2D matrices are indexed (cnt_glob,iref)
        if( allocated(self%corrmat2d) ) deallocate(self%corrmat2d)
        if( allocated(self%inplmat)   ) deallocate(self%inplmat)
        allocate(   self%corrmat2d(self%nrefs,self%nrefs),&
                    self%inplmat(self%nrefs,self%nrefs), stat=alloc_stat)
        call alloc_err("simple_prime_srch :: calc_corrs, 1", alloc_stat)
        self%corrmat2d = 0.
        self%inplmat   = 0
        ! calculate all correlations needed for search
        select case(trim(mmode))
            case('gpu','cpu')
                allocate(corrmat3d(self%nrefs,self%nrefs,self%nrots), stat=alloc_stat)
                call alloc_err("simple_prime_srch :: calc_corrs, 2", alloc_stat)
                corrmat3d = 0.
                ! we need to expand the dimensions in pftcc before calculating correlations using these modes
                call pftcc%expand_dim
                if( trim(mmode) .eq. 'gpu' )then
                    stop 'simple_polarft_corrcalc :: gencorrs_all_gpu taken out at the moment'
                else
                    call pftcc%gencorrs_all_cpu(corrmat3d)
                endif
                ! process corrmat3d to accomplish the different modes of stochastic search
                select case( refine )
                    case('no')
                        call project_corrmat3D_greedy(self%nrefs, self%nrots, corrmat3d, self%corrmat2d, self%inplmat)
                    case('shc')
                        if( present(prev_corrs) )then
                            call project_corrmat3D_shc(self%nrefs, self%nrots, corrmat3d, prev_corrs, self%corrmat2d, self%inplmat)
                        else
                            stop 'need optional input prev_corrs for this mode of search; simple_prime_srch :: calc_corrs'
                        endif
                    case DEFAULT
                        write(*,*) 'Refinement mode: ', trim(refine)
                        stop 'This refinement mode is not currently supported on GPU'
                end select
                deallocate(corrmat3d)
            case('bench')
                call pftcc%gencorrs_all_tester(self%corrmat2d, self%inplmat)
            case DEFAULT
                write(*,*) 'Unsupported flag (testflag): ', trim(mmode)
                stop 'simple_prime_srch :: calc_corrs'
        end select
    end subroutine calc_corrs
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(prime_srch), intent(inout) :: self !< instance
        if( self%exists )then
            deallocate(self%angtab)
            if( allocated(self%corrmat2d) ) deallocate(self%corrmat2d)
            if( allocated(self%inplmat)   ) deallocate(self%inplmat)
            self%pp     => null()
            self%nrefs  =  0
            self%nrots  =  0
            self%exists = .false.
        endif
    end subroutine kill

end module simple_prime_srch
