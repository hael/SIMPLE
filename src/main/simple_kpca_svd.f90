! kPCA using 'Learning to Find Pre-Images', using svd for eigvals/eigvecs
module simple_kpca_svd
include 'simple_lib.f08'
use simple_defs
use simple_pca, only: pca
implicit none

public :: kpca_svd
private

type, extends(pca) :: kpca_svd
    private
    real, allocatable :: E_zn(:,:)  !< expectations (feature vecs)
    real, allocatable :: data(:,:)  !< projected data on feature vecs
    logical           :: existence=.false.
    contains
    ! CONSTRUCTOR
    procedure :: new      => new_kpca
    ! GETTERS
    procedure :: get_feat => get_feat_kpca
    procedure :: generate => generate_kpca
    ! CALCULATORS
    procedure :: master   => master_kpca
    ! DESTRUCTOR
    procedure :: kill     => kill_kpca
    procedure :: kernel_center
end type

logical :: L_PRINT = .false.

contains

    ! CONSTRUCTORS

    !>  \brief  is a constructor
    subroutine new_kpca( self, N, D, Q )
        class(kpca_svd), intent(inout) :: self
        integer,         intent(in)    :: N, D, Q
        call self%kill
        self%N = N
        self%D = D
        self%Q = Q
        ! allocate principal subspace and feature vectors
        allocate( self%E_zn(self%Q,self%N), self%data(self%D,self%N), source=0.)
        self%existence = .true.
    end subroutine new_kpca

    ! GETTERS

    pure integer function get_N( self )
        class(kpca_svd), intent(in) :: self
        get_N = self%N
    end function get_N

    pure integer function get_D( self )
        class(kpca_svd), intent(in) :: self
        get_D = self%D
    end function get_D

    pure integer function get_Q( self )
        class(kpca_svd), intent(in) :: self
        get_Q = self%Q
    end function get_Q

    !>  \brief  is for getting a feature vector
    function get_feat_kpca( self, i ) result( feat )
        class(kpca_svd), intent(inout) :: self
        integer,    intent(in)    :: i
        real,       allocatable   :: feat(:)
        allocate(feat(self%Q), source=self%E_zn(:,i))
    end function get_feat_kpca

    !>  \brief  is for sampling the generative model at a given image index
    subroutine generate_kpca( self, i, avg, dat )
        class(kpca_svd), intent(inout) :: self
        integer,         intent(in)    :: i
        real,            intent(in)    :: avg(self%D)
        real,            intent(inout) :: dat(self%D)
        dat = avg + self%data(:,i)
    end subroutine generate_kpca

    ! CALCULATORS

    subroutine master_kpca( self, pcavecs, maxpcaits )
        class(kpca_svd),   intent(inout) :: self
        real,              intent(in)    :: pcavecs(self%D,self%N)
        integer, optional, intent(in)    :: maxpcaits ! redundant
    end subroutine master_kpca

    subroutine kernel_center( self, ker )
        class(kpca_svd), intent(inout) :: self
        real,            intent(inout) :: ker(self%N,self%N)
        real :: ones(self%N,self%N)
        ones = 1.
        ! Appendix D.2.2 Centering in Feature Space from Schoelkopf, Bernhard, Support vector learning, 1997
        ker  = ker - matmul(ones, ker) - matmul(ker, ones) + matmul(ones, matmul(ones, ker))
    end subroutine kernel_center

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill_kpca( self )
        class(kpca_svd), intent(inout) :: self
        if( self%existence )then
            deallocate( self%E_zn, self%data )
            self%existence = .false.
        endif
    end subroutine kill_kpca

end module simple_kpca_svd