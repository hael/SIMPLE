! PCA using standard SVD
module simple_pca
use simple_defs ! singleton
implicit none

public :: pca
private

type pca
    private
    integer           :: N           !< nr of data vectors
    integer           :: D           !< nr of components in each data vec
    integer           :: Q           !< nr of components in each latent vec
    real, allocatable :: E_zn(:,:,:) !< expectations (feature vecs)
    logical           :: existence=.false.
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! GETTERS
    procedure          :: get_N
    procedure          :: get_D
    procedure          :: get_Q
    procedure, private :: generate_1, generate_2
    generic            :: generate => generate_1, generate_2
    ! CALCULATORS
    procedure          :: master
    ! DESTRUCTOR
    procedure          :: kill
end type

contains

    ! CONSTRUCTORS

     !>  \brief  is a constructor
    function constructor( N, D, Q ) result( self )
        integer, intent(in) :: N, D, Q
        type(pca) :: self
        call self%new( N, D, Q )
    end function constructor

    !>  \brief  is a constructor
    subroutine new( self, N, D, Q )
        class(pca), intent(inout) :: self
        integer,    intent(in)    :: N, D, Q
        call self%kill
        self%N = N
        self%D = D
        self%Q = Q
        ! allocate principal subspace and feature vectors
        allocate( self%E_zn(1,self%Q,self%N), source=0.)
        self%existence = .true.
    end subroutine new

    ! GETTERS

    pure integer function get_N( self )
        class(pca), intent(in) :: self
        get_N = self%N
    end function get_N

    pure integer function get_D( self )
        class(pca), intent(in) :: self
        get_D = self%D
    end function get_D

    pure integer function get_Q( self )
        class(pca), intent(in) :: self
        get_Q = self%Q
    end function get_Q

    !>  \brief  is for sampling the generative model at a given image index
    subroutine generate_1( self, i, avg, dat )
        class(pca), intent(inout) :: self
        integer,    intent(in)    :: i
        real,       intent(in)    :: avg(self%D)
        real,       intent(inout) :: dat(self%D)

    end subroutine generate_1

    !>  \brief  is for sampling the generative model at a given image index
    subroutine generate_2( self, i, avg, dat, var )
        class(pca), intent(inout) :: self
        integer,    intent(in)    :: i
        real,       intent(in)    :: avg(self%D)
        real,       intent(inout) :: dat(self%D), var
        
    end subroutine generate_2

    ! CALCULATORS

    subroutine master( self, pcavecs )
        class(pca), intent(inout) :: self
        real,       intent(in)    :: pcavecs(self%D,self%N)
        
    end subroutine master

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill( self )
        class(pca), intent(inout) :: self
        if( self%existence )then
            deallocate( self%E_zn )
            self%existence = .false.
        endif
    end subroutine kill

end module simple_pca