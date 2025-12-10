! PCA using standard SVD
module simple_pca_svd
include 'simple_lib.f08'
use simple_defs
use simple_pca, only: pca
implicit none

public :: pca_svd
private

type, extends(pca) :: pca_svd
    private
    real, allocatable :: E_zn(:,:)  !< expectations (feature vecs)
    real, allocatable :: data(:,:)  !< projected data on feature vecs
    logical           :: existence=.false.
  contains
    ! CONSTRUCTOR
    procedure :: new      => new_svd
    ! GETTERS
    procedure :: get_feat => get_feat_svd
    procedure :: generate => generate_svd
    ! CALCULATORS
    procedure :: master   => master_svd
    procedure, private :: master_ori, master_T
    ! DESTRUCTOR
    procedure :: kill     => kill_svd
end type

logical :: L_PRINT = .false.

contains

    ! CONSTRUCTORS

    !>  \brief  is a constructor
    subroutine new_svd( self, N, D, Q )
        class(pca_svd), intent(inout) :: self
        integer,        intent(in)    :: N, D, Q
        call self%kill
        self%N = N
        self%D = D
        self%Q = Q
        ! allocate principal subspace and feature vectors
        allocate( self%E_zn(self%Q,self%N), self%data(self%D,self%N), source=0.)
        self%existence = .true.
    end subroutine new_svd

    ! GETTERS

    pure integer function get_N( self )
        class(pca_svd), intent(in) :: self
        get_N = self%N
    end function get_N

    pure integer function get_D( self )
        class(pca_svd), intent(in) :: self
        get_D = self%D
    end function get_D

    pure integer function get_Q( self )
        class(pca_svd), intent(in) :: self
        get_Q = self%Q
    end function get_Q

    !>  \brief  is for getting a feature vector
    function get_feat_svd( self, i ) result( feat )
        class(pca_svd), intent(inout) :: self
        integer,        intent(in)    :: i
        real,           allocatable   :: feat(:)
        allocate(feat(self%Q), source=self%E_zn(:,i))
    end function get_feat_svd

    !>  \brief  is for sampling the generative model at a given image index
    subroutine generate_svd( self, i, avg, dat )
        class(pca_svd), intent(inout) :: self
        integer,        intent(in)    :: i
        real,           intent(in)    :: avg(self%D)
        real,           intent(inout) :: dat(self%D)
        dat = avg + self%data(:,i)
    end subroutine generate_svd

    ! CALCULATORS

    subroutine master_svd( self, pcavecs, maxpcaits )
        class(pca_svd),    intent(inout) :: self
        real,              intent(in)    :: pcavecs(self%D,self%N)
        integer, optional, intent(in)    :: maxpcaits ! redundant for the svd approach
        if( self%D >= self%N )then
            call self%master_ori(pcavecs)
        else
            call self%master_T(pcavecs)
        endif
    end subroutine master_svd

    subroutine master_ori( self, pcavecs )
        class(pca_svd), intent(inout) :: self
        real,           intent(in)    :: pcavecs(self%D,self%N)
        real    :: eig_vecs(self%D,self%N), eig_vals(self%N), tmp(self%N,self%N)
        integer :: i, inds(self%N)
        eig_vecs = pcavecs
        call svdcmp(eig_vecs, eig_vals, tmp)
        eig_vals  = eig_vals**2 / real(self%D)
        inds      = (/(i,i=1,self%N)/)
        call hpsort(eig_vals, inds)
        call reverse(eig_vals)
        call reverse(inds)
        eig_vecs  = eig_vecs(:, inds)
        self%E_zn = transpose(matmul(transpose(pcavecs), eig_vecs(:,1:self%Q)))
        if( L_PRINT )then
            print *, 'eigenvalues:'
            print *, eig_vals
            print *, 'eigenvectors:'
            do i = 1, self%Q
                print *, eig_vecs(:,i)
            enddo
        endif
        self%data = matmul(eig_vecs(:,1:self%Q), self%E_zn)
    end subroutine master_ori

    subroutine master_T( self, pcavecs )
        class(pca_svd), intent(inout) :: self
        real,           intent(in)    :: pcavecs(self%D,self%N)
        real    :: eig_vecs(  self%D, self%N), eig_vals(  self%N), pcavecs_T(self%N, self%D), &
                   eig_vecs_T(self%N,self%D), eig_vals_T(self%D), tmp_T(self%D, self%D)
        integer :: i, inds(self%D), min_ND
        pcavecs_T  = transpose(pcavecs)
        eig_vecs_T = pcavecs_T
        inds       = (/(i,i=1,self%D)/)
        call svdcmp(eig_vecs_T, eig_vals_T, tmp_T)
        call hpsort(eig_vals_T, inds)
        call reverse(eig_vals_T)
        call reverse(inds)
        eig_vecs_T = eig_vecs_T(:, inds)
        tmp_T      =      tmp_T(:, inds)
        eig_vals   = 0.
        eig_vecs   = 0.
        min_ND     = min(self%N, self%D)
        eig_vals(  1:min_ND) = eig_vals_T(  1:min_ND)**2 / real(self%D)
        eig_vecs(:,1:min_ND) =      tmp_T(:,1:min_ND)
        self%E_zn = transpose(matmul(transpose(pcavecs), eig_vecs(:,1:self%Q)))
        if( L_PRINT )then
            print *, 'eigenvalues:'
            print *, eig_vals
            print *, 'eigenvectors:'
            do i = 1, self%Q
                print *, eig_vecs(:,i)
            enddo
        endif
        self%data = matmul(eig_vecs(:,1:self%Q), self%E_zn)
    end subroutine master_T

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill_svd( self )
        class(pca_svd), intent(inout) :: self
        if( self%existence )then
            deallocate( self%E_zn, self%data )
            self%existence = .false.
        endif
    end subroutine kill_svd

end module simple_pca_svd