! PCA using standard SVD
module simple_pca
include 'simple_lib.f08'
use simple_defs ! singleton
implicit none

public :: pca
private

type pca
    private
    integer           :: N          !< nr of data vectors
    integer           :: D          !< nr of components in each data vec
    integer           :: Q          !< nr of components in each latent vec
    real, allocatable :: E_zn(:,:)  !< expectations (feature vecs)
    real, allocatable :: data(:,:)  !< projected data on feature vecs
    logical           :: existence=.false.
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! GETTERS
    procedure          :: get_N
    procedure          :: get_D
    procedure          :: get_Q
    procedure          :: get_feat
    procedure, private :: generate_1, generate_2
    generic            :: generate => generate_1, generate_2
    ! CALCULATORS
    procedure          :: master
    ! DESTRUCTOR
    procedure          :: kill
end type

logical :: L_PRINT = .true.

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
        allocate( self%E_zn(self%Q,self%N), self%data(self%D,self%N), source=0.)
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

    !>  \brief  is for getting a feature vector
    function get_feat( self, i ) result( feat )
        class(pca), intent(inout) :: self
        integer,    intent(in)    :: i
        real,       allocatable   :: feat(:)
        allocate(feat(self%Q), source=self%E_zn(:,i))
    end function get_feat

    !>  \brief  is for sampling the generative model at a given image index
    subroutine generate_1( self, i, dat )
        class(pca), intent(inout) :: self
        integer,    intent(in)    :: i
        real,       intent(inout) :: dat(self%D)
        dat = self%data(:,i)
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
        real    :: avg(self%D), eig_vecs(self%D,self%N), eig_vals(self%N), tmp(self%N,self%N), pcavecs_cen(self%D,self%N)
        integer :: i
        avg         = sum(pcavecs, dim=2) / real(self%N)
        pcavecs_cen = pcavecs
        do i = 1, self%D
            pcavecs_cen(i,:) = pcavecs_cen(i,:) - avg(i)
        enddo
        eig_vecs = pcavecs_cen
        call svdcmp(eig_vecs, eig_vals, tmp)
        eig_vals  = eig_vals**2 / real(self%D)
        self%E_zn = transpose(matmul(transpose(pcavecs_cen), eig_vecs(:,1:self%Q)))
        if( L_PRINT )then
            print *, 'eigenvalues:'
            print *, eig_vals
            print *, 'eigenvectors:'
            do i = 1, self%Q
                print *, eig_vecs(:,i)
            enddo
        endif
        self%data = matmul(eig_vecs(:,1:self%Q), self%E_zn)
        do i = 1, self%D
            self%data(i,:) = self%data(i,:) + avg(i)
        enddo
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