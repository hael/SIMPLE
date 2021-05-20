module simple_ppca
include 'simple_lib.f08'
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
    real(dp), allocatable, public :: evals(:)
    real(dp), allocatable, public :: W(:,:)      !< principal subspace, defined by the W columns
    real(dp), allocatable :: E_zn(:,:,:) !< expectations (feature vecs)
    real(dp), allocatable :: W_1(:,:), W_2(:,:), W_3(:,:), Wt(:,:)
    real(dp), allocatable :: M(:,:), Minv(:,:), MinvWt(:,:), E_znzn(:,:)
    real(sp), allocatable :: X(:,:)
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
    procedure          :: calc_pc_corrs
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
        self%E_znzn(self%Q,self%Q), source = 0.d0, stat=alloc_stat )
        if(alloc_stat/=0)call allocchk('new; simple_ppca 1')
        allocate(self%X(self%D,1), source = 0., stat=alloc_stat )
        if(alloc_stat/=0)call allocchk('new; simple_ppca 2')
        if( self%dorot )then
            allocate( self%evals(self%Q), stat=alloc_stat )
            if(alloc_stat/=0)call allocchk('new; simple_ppca 3')
        endif
        self%existence = .true.
    end subroutine new

    ! GETTERS

    !>  \brief  is for getting the eigenvalues after rotation to the orthogonal basis
    function get_evals( self )result( e )
        class(ppca), intent(inout) :: self
        real :: e(self%Q)
        e = real(self%evals,sp)
    end function get_evals

    !>  \brief  is for getting a feature vector
    function get_feat( self, i ) result( feat )
        class(ppca), intent(inout) :: self
        integer, intent(in)        :: i
        real, allocatable          :: feat(:)
        allocate(feat(self%Q), source=real(self%E_zn(i,:,1),sp))
    end function get_feat

    !>  \brief  is for getting a pointer to the features (that become exposed)
    function get_feats_ptr( self ) result( ptr )
        class(ppca), intent(inout), target :: self
        real(dp), pointer :: ptr(:,:)
        ptr => self%E_zn(:,:,1)
    end function get_feats_ptr

    !>  \brief  is for getting a pointer to the loadings (that become exposed)
    function get_loadings_ptr( self ) result( ptr )
        class(ppca), intent(inout), target :: self
        real(dp), pointer :: ptr(:,:)
        ptr => self%W(:,:)
    end function get_loadings_ptr

    !>  \brief  is for rotation the PPCA subspace onto an orthogonal basis
    !>  Only performed when data accessed from memory
    subroutine rotation( self )
        class(ppca), intent(inout) :: self
        real(dp) :: Tut(self%Q,self%N), evecs(self%Q,self%Q)
        integer  :: i,nrot
        if( .not.self%inplace )stop 'can only perform cumulative sampling from memory; simple_ppca; generate_cumul'
        !$omp parallel proc_bind(close) private(i)
        !$omp workshare
        evecs = 0.d0
        self%W_1 = self%W
        !$omp end workshare nowait
        !$omp single
        call svdcmp(self%W_1, self%evals, evecs)                        ! SVD; now W_1 holds basis vectors (eg U in A=UWVt)
        !$omp end single
        !$omp workshare
        TUt = transpose( matmul( self%ptr_Dmat,self%W_1 ) )             ! Data Rotation
        !$omp end workshare
        !$omp do
        do i=1,self%Q
            TUt(i,:) = TUt(i,:)-sum(TUt(i,:))/real(self%N,dp)           ! Centering
        enddo
        !$omp end do
        !$omp workshare
        self%M = matmul( TUt, transpose(TUt) )/real(self%N,dp)          ! Var-covar
        !$omp end workshare
        !$omp single
        call jacobi( self%M, self%Q, self%Q, self%evals, evecs, nrot )  ! Eigen decomposition
        call eigsrt( self%evals,evecs,self%Q,self%Q )                   ! sorting
        !$omp end single
        !$omp workshare
        self%W           = matmul( self%W_1,evecs )                     ! new loadings
        self%E_zn(:,:,1) = matmul( real(self%ptr_Dmat,dp),self%W )      ! new expectations
        !$omp end workshare
        !$omp end parallel
    end subroutine rotation

    !>  \brief  is for sampling the generative model at a given image index for a subset of features
    !>  should only be used if rotation onto orthogonal basis has been performed
    function generate_cumul( self, i, feats, AVG ) result( dat )
        class(ppca),    intent(inout) :: self
        integer,        intent(in)    :: i, feats(2)
        real, optional, intent(in)    :: AVG(self%D)
        real, allocatable :: dat(:)
        real :: tmp(self%D,1)
        !$omp parallel workshare proc_bind(close)
        tmp = real(matmul(self%W(:,feats(1):feats(2)),self%E_zn(i,feats(1):feats(2),:)),sp)
        !$omp end parallel workshare
        allocate(dat(self%D), source=tmp(:,1))
        if( present(AVG) ) dat = dat+AVG
    end function generate_cumul


    subroutine calc_pc_corrs( self, avg, corrs )
        class(ppca),       intent(inout) :: self
        real,              intent(in)    :: AVG(self%D)
        real, allocatable, intent(out)   :: corrs(:)
        real(dp) :: ccs(self%Q), x(self%D), y(self%D), normsq
        integer  :: i,j
        if(allocated(corrs)) deallocate(corrs)
        allocate(corrs(self%Q),source=0.)
        ccs = 0.d0
        !$omp parallel do proc_bind(close) private(i,j,x,y,normsq) reduction(+:ccs)
        do j = 1,self%N
            x      = real(self%ptr_Dmat(j,:),dp)
            y      = real(avg,dp)
            normsq = sum(x*x)
            do i = 1,self%Q
                y      = y + self%E_zn(j,i,1) * self%W(:,i)
                ccs(i) = ccs(i) + sum(x*y) / sqrt(sum(y*y)*normsq)
            enddo
        enddo
        !$omp end parallel do
        corrs = real(ccs / real(self%N,dp))
    end subroutine calc_pc_corrs

    !>  \brief  is for sampling the generative model at a given image index
    function generate_1( self, i, AVG ) result( dat )
        class(ppca),    intent(inout) :: self
        integer,        intent(in)    :: i
        real, optional, intent(in)    :: AVG(self%D)
        real, allocatable :: dat(:)
        real :: tmp(self%D,1)
        tmp = real(matmul(self%W,self%E_zn(i,:,:)),sp)
        allocate(dat(self%D), source=tmp(:,1))
        if( present(AVG) ) dat = dat+AVG
    end function generate_1

    !>  \brief  produces the feature space average
    function generate_featavg( self,AVG )result( dat )
        class(ppca), intent(inout) :: self
        real,        intent(in)    :: AVG(self%D)
        real(dp), allocatable :: feat(:)
        real(sp), allocatable :: dat(:)
        integer :: i
        allocate( dat(self%D),feat(self%Q) )
        feat = 0.d0
        do i=1,self%N
            feat = feat + self%E_zn(i,:,1)/real(self%N,dp)
        enddo
        dat = self%generate_2( real(feat,sp),AVG )
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
        tmp2 = real(matmul(self%W,real(tmp1,dp)),sp)
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
        class(ppca),                intent(inout) :: self
        character(len=*),           intent(in)    :: datastk, featstk
        integer,                    intent(in)    :: recsz, maxpcaits
        character(len=*), optional, intent(in)    :: feats_txt
        integer  :: k, file_stat, funit2, recsz2, err, fhandle_txt
        real(dp) :: p, p_prev
        self%doprint = .true.
        write(*,'(A)') '>>> GENERATIVE ITERATIVE PCA'
        if( .not.self%inplace )then
            open(newunit=self%funit, status='old', action='read', file=datastk,&
            access='direct', form='unformatted', recl=recsz, iostat=file_stat)
            inquire( iolength=recsz2 ) self%E_zn(1,:,1)
            open(newunit=funit2, status='replace', action='write', file=featstk,&
            access='direct', form='unformatted', recl=recsz2, iostat=file_stat)
        endif
        p = 0.d0
        k = 0.d0
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
            write(*,"(A,1X,I3,1X,A,1X,F10.1)") '>>> Iteration:', k, 'Squared error:', p
            if( (abs(p-p_prev) < 0.1d0) .or. k == maxpcaits ) exit
        end do
        ! write features to disk
        if( .not.self%inplace )then
            do k=1,self%N
                write(funit2, rec=k) real(self%E_zn(k,:,1),sp)
                if( present(feats_txt) ) write(fhandle_txt,*) real(self%E_zn(k,:,1),sp)
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
        integer  :: i, j
        real     :: meanv(self%Q), Imat(self%Q,self%Q)
        ! make identity matrices
        Imat=0.; do i=1,self%Q ; Imat(i,i)=1. ; end do
        meanv = 0.
        !$omp parallel default(shared) private(i,j) proc_bind(close)
        !$omp do schedule(static)
        ! initialize latent variables by zero mean gaussians with unit variance
        do i=1,self%N
            self%E_zn(i,:,1) = real(mnorm_smp(Imat, self%Q, meanv),dp)
        end do
        !$omp end do nowait
        ! initialize weight matrix with uniform random nrs, normalize over j
        !$omp do schedule(static)
        do i=1,self%D
            do j=1,self%Q
                self%W(i,j) = real(ran3(),dp)
            end do
            self%W(i,:) = self%W(i,:)/sum(self%W(i,:))
        end do
        !$omp end do
        !$omp workshare
        ! transpose W
        self%Wt = transpose(self%W)
        !$omp end workshare
        !$omp end parallel
    end subroutine init

    !>  \brief  EM algorithm
    subroutine em_opt( self, p, err )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_math, only: matinv, is_a_number
        class(ppca), intent(inout) :: self
        integer,     intent(out)   :: err
        real(dp),    intent(out)   :: p
        integer                    :: i
        real(dp)                   :: tmp(self%D,1), x(self%D,1), tE_zn(1,self%Q), E_zn(self%Q,1)
        ! E-STEP
        !$omp parallel proc_bind(close)
        !$omp workshare
        self%M = matmul(self%Wt,self%W)
        !$omp end workshare
        !$omp single
        call matinv(self%M, self%Minv, self%Q, err)
        !$omp end single
        !$omp workshare
        self%MinvWt = matmul(self%Minv,self%Wt)
        self%W_1    = 0.d0
        self%W_2    = 0.d0
        !$omp end workshare
        !$omp end parallel
        if( err == -1 ) return
        if( self%inplace )then
            !$omp parallel proc_bind(close)
            do i=1,self%N
                !$omp workshare
                x(:,1) = real(self%ptr_Dmat(i,:),dp)
                ! Expectation step (calculate expectations using the old W)
                E_zn = matmul(self%MinvWt, x)
                self%E_zn(i,:,:) = E_zn
                tE_zn            = transpose(E_zn)
                ! Prepare for update of W (M-step)
                self%W_1 = self%W_1 + matmul(   x, tE_zn)
                self%W_2 = self%W_2 + matmul(E_zn, tE_zn)
                !$omp end workshare
            end do
            !$omp end parallel
        else
            do i=1,self%N
                read(self%funit, rec=i) self%X(:,1) ! read from file
                !!$omp parallel workshare
                ! Expectation step (calculate expectations using the old W)
                self%E_zn(i,:,:) = matmul(self%MinvWt,real(self%X,dp))
                self%E_znzn = matmul(self%E_zn(i,:,:),transpose(self%E_zn(i,:,:)))
                ! Prepare for update of W (M-step)
                self%W_1 = self%W_1+matmul(real(self%X,dp),transpose(self%E_zn(i,:,:)))
                self%W_2 = self%W_2+self%E_znzn
                !!$omp end parallel workshare
            end do
        endif
        ! M-STEP
        call matinv(self%W_2, self%W_3, self%Q, err)
        if( err == -1 ) return
        !$omp parallel workshare proc_bind(close)
        ! update W
        self%W = matmul(self%W_1,self%W_3)
        ! update Wt
        self%Wt = transpose(self%W)
        !$omp end parallel workshare
        ! EVAL REC ERR
        p = 0.d0
        if( self%inplace )then
            !$omp parallel proc_bind(close)
            do i=1,self%N
                !$omp workshare
                tmp = matmul(self%W,self%E_zn(i,:,:))
                p   = p+sqrt(sum((real(self%ptr_Dmat(i,:),dp)-tmp(:,1))**2.d0))
                !$omp end workshare
            end do
            !$omp end parallel
        else
            do i=1,self%N
                ! read data vec
                read(self%funit, rec=i) self%X(:,1)
                !!$omp parallel workshare default(shared) does not compile on all platforms for some reason
                tmp = matmul(self%W,self%E_zn(i,:,:))
                !!$omp end parallel workshare
                p = p+sqrt(sum((real(self%X(:,1),dp)-tmp(:,1))**2.d0))
            end do
        endif
        if( .not. is_a_number(p) )err=-1
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
