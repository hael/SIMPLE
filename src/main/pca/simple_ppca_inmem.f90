!==Class simple_ppca_inmem
!
! simple_ppca_inmem is the SIMPLE class for probabilistic principal component analysis.
! This code should be able to deal with many millions of particle images.
! The code is distributed with the hope that it will be useful, but _WITHOUT_
! _ANY_ _WARRANTY_. Redistribution or modification is regulated by the
! GNU General Public License. *Author:* Hans Elmlund, 2011-09-03.
!
!==Changes are documented below
!
module simple_ppca_inmem
use simple_defs ! singleton
use simple_pca, only: pca
implicit none

public :: ppca_inmem
private

type, extends(pca) :: ppca_inmem
    private
    real, allocatable :: W(:,:)      !< principal subspace, defined by the W columns
    real, allocatable :: E_zn(:,:,:) !< expectations (feature vecs)
    real, allocatable :: W_1(:,:), W_2(:,:), W_3(:,:), Wt(:,:)
    real, allocatable :: M(:,:), Minv(:,:), MinvWt(:,:), X(:,:), E_znzn(:,:)
    logical           :: existence=.false.
  contains
    ! CONSTRUCTOR
    procedure          :: new      => new_inmem
    ! GETTERS
    procedure          :: get_feat => get_feat_inmem
    procedure          :: generate => generate_inmem
    ! CALCULATORS
    procedure          :: master   => master_inmem
    procedure, private :: init
    procedure, private :: em_opt
    ! DESTRUCTOR
    procedure          :: kill     => kill_inmem
end type

logical :: L_PRINT = .false.

contains

    ! CONSTRUCTORS

    !>  \brief  is a constructor
    subroutine new_inmem( self, N, D, Q )
        class(ppca_inmem), intent(inout) :: self
        integer,           intent(in)    :: N, D, Q
        call self%kill
        self%N = N
        self%D = D
        self%Q = Q
        ! allocate principal subspace and feature vectors
        allocate( self%W(self%D,self%Q), self%E_zn(1,self%Q,self%N), self%W_1(self%D,self%Q),&
        self%W_2(self%Q,self%Q), self%W_3(self%Q,self%Q), self%Wt(self%Q,self%D),&
        self%M(self%Q,self%Q), self%Minv(self%Q,self%Q), self%MinvWt(self%Q,self%D),&
        self%X(self%D,1), self%E_znzn(self%Q,self%Q), source=0.)
        self%existence = .true.
    end subroutine new_inmem

    ! GETTERS

    !>  \brief  is for getting a feature vector
    function get_feat_inmem( self, i ) result( feat )
        class(ppca_inmem), intent(inout) :: self
        integer,           intent(in)    :: i
        real,              allocatable   :: feat(:)
        allocate(feat(self%Q), source=self%E_zn(1,:,i))
    end function get_feat_inmem

    !>  \brief  is for getting a pointer to the features (that become exposed)
    function get_feats_ptr( self ) result( ptr )
        class(ppca_inmem), intent(inout), target :: self
        real, pointer :: ptr(:,:)
        ptr => self%E_zn(1,:,:)
    end function get_feats_ptr

    !>  \brief  is for sampling the generative model at a given image index
    subroutine generate_inmem( self, i, avg, dat )
        class(ppca_inmem), intent(inout) :: self
        integer,           intent(in)    :: i
        real,              intent(in)    :: avg(self%D)
        real,              intent(inout) :: dat(self%D)
        dat = avg + [matmul(self%W,transpose(self%E_zn(:,:,i)))]
    end subroutine generate_inmem

    ! CALCULATORS

    subroutine master_inmem( self, pcavecs, maxpcaits )
        class(ppca_inmem), intent(inout) :: self
        real,              intent(in)    :: pcavecs(self%D,self%N)
        integer, optional, intent(in)    :: maxpcaits
        integer :: k, err
        real    :: p, p_prev
        p = 0.
        k = 0
        call self%init
        do
            k      = k + 1
            p_prev = p
            call self%em_opt( pcavecs, p, err )
            if( err == -1 )then
                write(logfhandle,'(A)') 'ERROR, in matrix inversion, iteration:', k
                write(logfhandle,'(A)') 'RESTARTING'
                call self%init
                k = 0
                cycle
            endif
            if( k == 1 .or. mod(k,5) == 0 )then
                if( L_PRINT ) write(logfhandle,"(1X,A,1X,I3,1X,A,1X,F10.0)") 'Iteration:', k, 'Squared error:', p
            endif
            if( (abs(p-p_prev) < 0.1) .or. k == maxpcaits ) exit
        end do
    end subroutine master_inmem

    subroutine init( self )
        use simple_rnd, only: mnorm_smp, ran3
        class(ppca_inmem), intent(inout) :: self
        integer :: i, j
        real    :: meanv(self%Q), Imat(self%Q,self%Q)
        ! make identity matrices
        Imat=0.; do i=1,self%Q ; Imat(i,i)=1. ; end do
        meanv = 0.
        ! initialize latent variables by zero mean gaussians with unit variance
        do i=1,self%N
            self%E_zn(1,:,i) = mnorm_smp(Imat, self%Q, meanv)
        end do
        ! initialize weight matrix with uniform random nrs, normalize over j
        do i=1,self%D
            do j=1,self%Q
                self%W(i,j) = ran3()
            end do
            self%W(i,:) = self%W(i,:)/sum(self%W(i,:))
        end do
        ! transpose W
        self%Wt  = transpose(self%W)
        ! set W_2 to WtW
        self%W_2 = matmul(self%Wt,self%W)
    end subroutine init

    !>  \brief  EM algorithm
    subroutine em_opt( self, pcavecs, p, err )
        use simple_math, only: matinv
        class(ppca_inmem), intent(inout) :: self
        real,              intent(in)    :: pcavecs(self%D,self%N)
        integer,           intent(out)   :: err
        real,              intent(out)   :: p
        integer :: i
        real    :: tmp(self%D,1)
        ! E-STEP
        self%M = matmul(self%Wt,self%W)
        call matinv(self%M, self%Minv, self%Q, err)
        if( err == -1 ) return
        self%MinvWt = matmul(self%Minv,self%Wt)
        self%W_1 = 0.
        self%W_2 = 0.
        do i=1,self%N
            ! set data vec
            self%X(:,1) = pcavecs(:,i)
            ! Expectation step (calculate expectations using the old W)
            self%E_zn(:,:,i) = matmul(self%MinvWt,self%X)
            self%E_znzn = matmul(transpose(self%E_zn(:,:,i)),self%E_zn(:,:,i))
            ! Prepare for update of W (M-step)
            self%W_1 = self%W_1+matmul(self%X,self%E_zn(:,:,i))
            self%W_2 = self%W_2+self%E_znzn
        end do
        ! M-STEP
        call matinv(self%W_2, self%W_3, self%Q, err)
        if( err == -1 ) return
        ! update W
        self%W   = matmul(self%W_1,self%W_3)
        ! update Wt
        self%Wt  = transpose(self%W)
        ! set W_2 to WtW
        self%W_2 = matmul(self%Wt,self%W)
        ! EVAL REC ERR
        p = 0.
        do i=1,self%N
            ! set data vec
            self%X(:,1) = pcavecs(:,i)
            tmp = matmul(self%W,transpose(self%E_zn(:,:,i)))
            p   = p + sqrt(sum((self%X(:,1)-tmp(:,1))**2.))
        end do
    end subroutine em_opt

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill_inmem( self )
        class(ppca_inmem), intent(inout) :: self
        if( self%existence )then
            deallocate( self%W, self%E_zn, self%W_1,&
            self%W_2, self%W_3, self%Wt, self%M, self%Minv,&
            self%MinvWt, self%X, self%E_znzn )
            self%existence = .false.
        endif
    end subroutine kill_inmem

end module simple_ppca_inmem