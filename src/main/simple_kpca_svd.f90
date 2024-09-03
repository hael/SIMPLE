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
    procedure :: rbf_kernel
    procedure :: test_all, test_kernel_center, test_rbf_kernel
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

    subroutine kernel_center( self, dim, ker )
        class(kpca_svd), intent(inout) :: self
        integer,         intent(in)    :: dim
        real,            intent(inout) :: ker(dim,dim)
        real :: ones(dim,dim)
        ones = 1. / dim
        ! Appendix D.2.2 Centering in Feature Space from Schoelkopf, Bernhard, Support vector learning, 1997
        ker  = ker - matmul(ones, ker) - matmul(ker, ones) + matmul(matmul(ones, ker), ones)
    end subroutine kernel_center

    subroutine rbf_kernel( self, dim, mat, ncomp, c, ker )
        class(kpca_svd), intent(inout) :: self
        integer,         intent(in)    :: dim
        real,            intent(in)    :: mat(dim,dim)
        integer,         intent(in)    :: ncomp
        real,            intent(in)    :: c
        real,            intent(out)   :: ker(dim,dim)
        integer :: i, j
        ! squared euclidean distance between pairs of rows
        do i = 1,dim
            do j = 1,dim
                ker(i,j) = euclid(mat(i,:), mat(j,:))**2
            enddo
        enddo
        ker = exp(-ker/real(ncomp)/c)
        call self%kernel_center(dim, ker)
    end subroutine rbf_kernel

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill_kpca( self )
        class(kpca_svd), intent(inout) :: self
        if( self%existence )then
            deallocate( self%E_zn, self%data )
            self%existence = .false.
        endif
    end subroutine kill_kpca

    ! UNIT TESTS
    subroutine test_all( self )
        class(kpca_svd), intent(inout) :: self
        call self%test_kernel_center
        call self%test_rbf_kernel
    end subroutine test_all

    subroutine test_kernel_center( self )
        class(kpca_svd), intent(inout) :: self
        real :: ker(2,2), truth(2,2)
        ker(1,:) = [ 1., 2.]
        ker(2,:) = [ 3., 1.]
        call self%kernel_center(2, ker)
        truth(1,:) = [ -0.75,  0.75 ]
        truth(2,:) = [  0.75, -0.75 ]
        if( maxval(abs(ker - truth)) < TINY )then
            print *, 'Unit test of kernel_center: PASSED'
        else
            print *, 'Unit test of kernel_center: FAILED'
            print *, 'ker = ', ker
        endif
    end subroutine test_kernel_center

    subroutine test_rbf_kernel( self )
        class(kpca_svd), intent(inout) :: self
        real :: mat(2,2), truth(2,2), ker(2,2)
        mat(1,:) = [ 0., 2.]
        mat(2,:) = [ 1., 1.]
        call self%rbf_kernel(2, mat, 2, 1., ker)
        truth(1,:) = [  0.31606028, -0.31606028 ]
        truth(2,:) = [ -0.31606028,  0.31606028 ]
        if( maxval(abs(ker - truth)) < TINY )then
            print *, 'Unit test of rbf_kernel: PASSED'
        else
            print *, 'Unit test of rbf_kernel: FAILED'
            print *, 'ker = ', ker
        endif
    end subroutine test_rbf_kernel

end module simple_kpca_svd