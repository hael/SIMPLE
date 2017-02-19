!==Class simple_ppca
!
! simple_ppca is the SIMPLE class for probabilistic principal component analysis.
! This code should be able to deal with many millions of particle images. 
! The code is distributed with the hope that it will be useful, but _WITHOUT_ 
! _ANY_ _WARRANTY_. Redistribution or modification is regulated by the 
! GNU General Public License. *Author:* Hans Elmlund, 2011-09-03.
!
!==Changes are documented below
!
module simple_ppca
use simple_defs
use simple_jiffys, only: alloc_err
implicit none

public :: ppca
private

type ppca
    private
    integer              :: N           !< nr of data vectors
    integer              :: D           !< nr of components in each data vec
    integer              :: Q           !< nr of components in each latent vec
    integer              :: funit       !< file handle data stack
    real, pointer        :: ptr_Dmat(:,:) => null()
    real, allocatable    :: evals(:)
    real, allocatable    :: W(:,:)      !< principal subspace, defined by the W columns
    real, allocatable    :: E_zn(:,:,:) !< expectations (feature vecs)
    real, allocatable    :: W_1(:,:), W_2(:,:), W_3(:,:), Wt(:,:)
    real, allocatable    :: M(:,:), Minv(:,:), MinvWt(:,:), X(:,:), E_znzn(:,:)
    logical              :: inplace     = .false. !< controls whether PPCA performed in memory
    logical              :: dorot       = .false. !< is to perform final orthogonal basis rotation
    logical              :: existence   = .false.
    logical              :: doprint     = .false.
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! GETTERS
    procedure          :: get_feat
    procedure          :: get_feats_ptr
    procedure          :: get_loadings_ptr
    procedure          :: get_evals
    procedure          :: generate_cumul
    procedure          :: generate_featavg
    procedure, private :: rotation
    procedure, private :: generate_1
    procedure, private :: generate_2
    generic            :: generate => generate_1, generate_2
    ! CALCULATORS
    procedure          :: master
    procedure          :: init
    procedure, private :: em_opt
    ! DESTRUCTOR
    procedure :: kill
end type

contains

    ! CONSTRUCTORS

     !>  \brief  is a constructor
    function constructor( N, D, Q, Dmat ) result( self )
        real, intent(in), dimension(:,:), optional :: Dmat
        integer, intent(in)           :: N, D, Q
        type(ppca)                    :: self
        call self%new( N, D, Q, Dmat )
    end function constructor

    !>  \brief  is a constructor
    subroutine new( self, N, D, Q, Dmat, dorot, doprint )
        class(ppca), intent(inout)         :: self
        real, intent(in), dimension(:,:), optional, target :: Dmat
        integer, intent(in)                :: N, D, Q
        logical, intent(in), optional      :: dorot, doprint
        integer                            :: alloc_stat
        call self%kill
        if( present(doprint) )self%doprint = doprint
        if( present(dorot) )  self%dorot   = dorot
        if( present(Dmat) )then
            self%inplace  = .true.
            self%ptr_Dmat => Dmat
        endif
        self%N = N
        self%D = D
        self%Q = Q
        ! allocate principal subspace and feature vectors
        allocate( self%W(self%D,self%Q), self%E_zn(self%N,self%Q,1), self%W_1(self%D,self%Q),&
        self%W_2(self%Q,self%Q), self%W_3(self%Q,self%Q), self%Wt(self%Q,self%D),&
        self%M(self%Q,self%Q), self%Minv(self%Q,self%Q), self%MinvWt(self%Q,self%D),&
        self%X(self%D,1), self%E_znzn(self%Q,self%Q), stat=alloc_stat )
        call alloc_err( 'new; simple_ppca 1', alloc_stat )
        if( self%dorot )then
            allocate( self%evals(self%Q), stat=alloc_stat )
            call alloc_err( 'new; simple_ppca 2', alloc_stat )
        endif
        self%W      = 0.
        self%E_zn   = 0.
        self%W_1    = 0.
        self%W_2    = 0.
        self%W_3    = 0.
        self%Wt     = 0.
        self%M      = 0.
        self%Minv   = 0.
        self%MinvWt = 0.
        self%X      = 0.
        self%E_znzn = 0.
        self%existence = .true.
    end subroutine new
    
    ! GETTERS
    
    !>  \brief  is for getting the eigenvalues after rotation to the orthogonal basis
    function get_evals( self )result( e )
        class(ppca), intent(inout) :: self
        real :: e(self%Q)
        e = self%evals
    end function get_evals        

    !>  \brief  is for getting a feature vector
    function get_feat( self, i ) result( feat )
        class(ppca), intent(inout) :: self
        integer, intent(in)        :: i
        real, allocatable          :: feat(:)
        allocate(feat(self%Q), source=self%E_zn(i,:,1))
    end function get_feat
    
    !>  \brief  is for getting a pointer to the features (that become exposed)
    function get_feats_ptr( self ) result( ptr )
        class(ppca), intent(inout), target :: self
        real, pointer :: ptr(:,:)
        ptr => self%E_zn(:,:,1)
    end function get_feats_ptr

    !>  \brief  is for getting a pointer to the loadings (that become exposed)
    function get_loadings_ptr( self ) result( ptr )
        class(ppca), intent(inout), target :: self
        real, pointer :: ptr(:,:)
        ptr => self%W(:,:)
    end function get_loadings_ptr

    !>  \brief  is for rotation the PPCA subspace onto an orthogonal basis
    !>  Only performed when data accessed from memory
    subroutine rotation( self )
        use simple_math, only: jacobi,eigsrt,svdcmp
        class(ppca), intent(inout) :: self
        real, allocatable :: TUt(:,:),evecs(:,:)
        real    :: tmp(self%D,1)
        integer :: i,j,alloc_stat, nrot
        if( .not.self%inplace )stop 'can only perform cumulative sampling from memory; simple_ppca; generate_cumul'
        allocate( TUt(self%Q,self%N),evecs(self%Q,self%Q), &
            stat=alloc_stat )
        call alloc_err( 'generate_cumul; simple_ppca', alloc_stat )
        self%W_1 = self%W
        call svdcmp( self%W_1, self%D, self%Q, self%D, self%Q, self%evals, evecs ) ! SVD; now W_1 holds basis vectors (eg U in A=UWVt)
        TUt = transpose( matmul( self%ptr_Dmat,self%W_1 ) )                        ! Data Rotation
        do i=1,self%Q
           TUt(i,:) = TUt(i,:)-sum(TUt(i,:))/real(self%N)               ! Centering
        enddo
        self%M = matmul( TUt, transpose(TUt) )/real(self%N)             ! Var-covar
        call jacobi( self%M, self%Q, self%Q, self%evals, evecs, nrot )  ! Eigen decomposition
        call eigsrt( self%evals,evecs,self%Q,self%Q )                   ! sorting
        self%W           = matmul( self%W_1,evecs )                     ! new loadings
        self%E_zn(:,:,1) = matmul( self%ptr_Dmat,self%W )               ! new expectations
        deallocate( TUt,evecs )
    end subroutine rotation

    !>  \brief  is for sampling the generative model at a given image index for a subset of features
    !>  should only be used if rotation onto orthogonal basis has been performed
    function generate_cumul( self, i, feats, AVG ) result( dat ) 
        class(ppca),    intent(inout) :: self
        integer,        intent(in)    :: i, feats(2)
        real, optional, intent(in)    :: AVG(self%D)
        real, allocatable :: dat(:)
        real :: tmp(self%D,1)
        tmp = matmul(self%W(:,feats(1):feats(2)),self%E_zn(i,feats(1):feats(2),:))
        allocate(dat(self%D), source=tmp(:,1))
        if( present(AVG) ) dat = dat+AVG
    end function generate_cumul

    !>  \brief  is for sampling the generative model at a given image index
    function generate_1( self, i, AVG ) result( dat ) 
        class(ppca),    intent(inout) :: self
        integer,        intent(in)    :: i
        real, optional, intent(in)    :: AVG(self%D)
        real, allocatable :: dat(:)
        real :: tmp(self%D,1)
        tmp = matmul(self%W,self%E_zn(i,:,:))
        allocate(dat(self%D), source=tmp(:,1))
        if( present(AVG) ) dat = dat+AVG
    end function generate_1

    !>  \brief  produces the feature space average
    function generate_featavg( self,AVG )result( dat )
        class(ppca), intent(inout) :: self
        real,        intent(in)    :: AVG(self%D)
        real, allocatable :: dat(:),feat(:)
        integer :: i
        allocate( dat(self%D),feat(self%Q) )
        feat = 0.
        do i=1,self%N
            feat = feat+self%E_zn(i,:,1)/real(self%N)
        enddo
        dat = self%generate_2( feat,AVG )
        deallocate( feat )
    end function generate_featavg

    !>  \brief  is for sampling the generative model at arbitrary feature
    !!          useful if doing averaging in feature space
    function generate_2( self, feat, AVG ) result( dat )
        class(ppca),    intent(in) :: self
        real,           intent(in) :: feat(self%Q)
        real, optional, intent(in) :: AVG(self%D)
        real, allocatable :: dat(:)
        real :: tmp1(self%Q,1)
        real :: tmp2(self%D,1)
        tmp1(:,1) = feat
        tmp2 = matmul(self%W,tmp1)
        allocate( dat(self%D) )
        if( present(AVG) )then
            dat = tmp2(:,1)+AVG
        else
            dat = tmp2(:,1)
        endif
    end function generate_2
    
    ! CALCULATORS
    
    !>  \brief  doing it all
    subroutine master( self, datastk, recsz, featstk, maxpcaits, feats_txt )
        use simple_filehandling, only: get_fileunit, fopen_err
        class(ppca),                intent(inout) :: self
        character(len=*),           intent(in)    :: datastk, featstk
        integer,                    intent(in)    :: recsz, maxpcaits
        character(len=*), optional, intent(in)    :: feats_txt
        integer :: k, file_stat, funit2, recsz2, err, fhandle_txt
        real    :: p, p_prev
        self%doprint = .true.
        write(*,'(A)') '>>> GENERATIVE ITERATIVE PCA'
        if( .not.self%inplace )then
            self%funit = get_fileunit()
            open(unit=self%funit, status='old', action='read', file=datastk,&
            access='direct', form='unformatted', recl=recsz, iostat=file_stat)
            call fopen_err('master; simple_ppca, 1', file_stat)
            inquire( iolength=recsz2 ) self%E_zn(1,:,1)
            funit2 = get_fileunit()
            open(unit=funit2, status='replace', action='write', file=featstk,&
            access='direct', form='unformatted', recl=recsz2, iostat=file_stat)
            call fopen_err('master; simple_ppca, 2', file_stat)
            fhandle_txt = get_fileunit()
            if( present(feats_txt) )then
                open(unit=fhandle_txt, status='replace', action='write', file=feats_txt, iostat=file_stat)
                call fopen_err('master; simple_ppca, 3', file_stat)
            endif
        endif
        p = 0.
        k = 0
        call self%init
        do
            k = k+1
            p_prev = p
            call self%em_opt( p, err )
            if( err == -1 )then
                write(*,'(A)') 'ERROR, in matrix inversion, iteration:', k
                write(*,'(A)') 'RESTARTING'
                call self%init
                k = 0 
                cycle
            endif
            if( self%doprint.and.(k == 1 .or. mod(k,5) == 0) )then
                write(*,"(1X,A,1X,I3,1X,A,1X,F10.0)") 'Iteration:', k, 'Squared error:', p
            endif
            if( (abs(p-p_prev) < 0.1) .or. k == maxpcaits ) exit
        end do
        ! write features to disk
        if( .not.self%inplace )then
            do k=1,self%N
                write(funit2, rec=k) self%E_zn(k,:,1)
                if( present(feats_txt) ) write(fhandle_txt,*) self%E_zn(k,:,1)
            end do
            close(unit=self%funit)
            close(unit=funit2)
            if( present(feats_txt) ) close(fhandle_txt)
        endif
        ! subspace rotation for unconstrained algorithm
        if( self%dorot ) call self%rotation
    end subroutine master

    subroutine init( self )
        use simple_rnd, only: mnorm_smp, ran3
        class(ppca), intent(inout) :: self
        integer :: i, j
        real    :: meanv(self%Q), Imat(self%Q,self%Q)
        ! make identity matrices
        Imat=0.; do i=1,self%Q ; Imat(i,i)=1. ; end do
        meanv = 0.
        ! initialize latent variables by zero mean gaussians with unit variance
        do i=1,self%N
            self%E_zn(i,:,1) = mnorm_smp(Imat, self%Q, meanv)
        end do
        ! initialize weight matrix with uniform random nrs, normalize over j
        do i=1,self%D
            do j=1,self%Q
                self%W(i,j) = ran3()
            end do
            self%W(i,:) = self%W(i,:)/sum(self%W(i,:))
        end do
        ! transpose W
        self%Wt = transpose(self%W)
    end subroutine init

    !>  \brief  EM algorithm
    subroutine em_opt( self, p, err )
        use ieee_arithmetic
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_math, only: matinv
        class(ppca), intent(inout) :: self
        integer,     intent(out)   :: err
        real,        intent(out)   :: p
        integer                    :: i
        real(sp)                   :: tmp(self%D,1)
        real(dp)                   :: mat_tmp(self%Q,self%Q)
        ! E-STEP
        self%M = matmul(self%Wt,self%W)
        call matinv(self%M, self%Minv, self%Q, err)
        self%Minv = real(mat_tmp)
        if( err == -1 ) return
        self%MinvWt = matmul(self%Minv,self%Wt)
        self%W_1    = 0.
        self%W_2    = 0.
        do i=1,self%N
           if( .not.self%inplace )then
              ! read/copy data vec
              read(self%funit, rec=i) self%X(:,1) ! read from file
           else
              self%X(:,1) = self%ptr_Dmat(i,:)
           endif
           ! Expectation step (calculate expectations using the old W)
           self%E_zn(i,:,:) = matmul(self%MinvWt,self%X)
           self%E_znzn = matmul(self%E_zn(i,:,:),transpose(self%E_zn(i,:,:)))
           ! Prepare for update of W (M-step)
           self%W_1 = self%W_1+matmul(self%X,transpose(self%E_zn(i,:,:)))
           self%W_2 = self%W_2+self%E_znzn
        end do
        ! M-STEP
        call matinv(self%W_2, self%W_3, self%Q, err)
        self%W_3 = real(mat_tmp)
        if( err == -1 ) return
        ! update W
        self%W = matmul(self%W_1,self%W_3)
        ! update Wt
        self%Wt = transpose(self%W)
        ! EVAL REC ERR
        p = 0.
        if( self%inplace )then
            !$omp parallel default(shared), reduction(+:p), private(i,tmp)
            !$omp do
            do i=1,self%N
                tmp = matmul(self%W,self%E_zn(i,:,:))
                p   = p+sqrt(sum((self%ptr_Dmat(i,:)-tmp(:,1))**2.))
            end do
            !$omp end do nowait
            !$omp end parallel
        else
            do i=1,self%N
                ! read data vec
                read(self%funit, rec=i) self%X(:,1)
                !$omp parallel default(shared)
                !$omp workshare
                tmp = matmul(self%W,self%E_zn(i,:,:))
                !$omp end workshare nowait
                !$omp end parallel
                p = p+sqrt(sum((self%X(:,1)-tmp(:,1))**2.))
            end do
        endif
        if( ieee_is_nan(p) )err=-1
    end subroutine em_opt

    ! DESTRUCTOR
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(ppca), intent(inout) :: self
        if( self%existence )then
            deallocate( self%W, self%E_zn, self%W_1,&
            self%W_2, self%W_3, self%Wt, self%M, self%Minv,&
            self%MinvWt, self%X, self%E_znzn )
            if( allocated(self%evals) ) deallocate(self%evals)
            self%existence = .false.
            self%inplace   = .false.
            self%ptr_Dmat  => null()
        endif
    end subroutine kill
    
end module simple_ppca
