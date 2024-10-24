! Implements two spectral clustering algorithms:
! 1. Unnormalized spectral clustering + kmeans
!   https://link.springer.com/article/10.1007/s11222-007-9033-z
! 2. Unnormalized spectral clustering + CPQR
!   https://doi.org/10.1093/imaiai/iay008
module simple_spectral_clustering
include 'simple_lib.f08'
use simple_ran_tabu
implicit none

public :: spec_clust, test_spec_clust
private
#include "simple_local_flags.inc"

real,    parameter :: GAMMA_DEFAULT  = 1.0
integer, parameter :: KMEANS_REPEATS = 100 
integer, parameter :: KMEANS_ITERS   = 200 

type spec_clust
    private
    type(ran_tabu)        :: rt                          ! for initial labels
    integer               :: N                           ! similarity matrix size
    integer               :: K                           ! # of clusters
    integer               :: neigs = 0                   ! # of selected eigenvectors
    real,         pointer :: S(:,:)                      ! pointer to similarity matrix
    real,     allocatable :: centroids(:,:)              ! cluster centers
    real,     allocatable :: features(:,:)               ! eigen vectors selected
    integer,  allocatable :: labels(:)                   ! clustering output
    character(len=STDLEN) :: algorithm = 'cpqr'
    real                  :: gamma     = GAMMA_DEFAULT   ! kernel coefficient
    integer               :: niters    = KMEANS_ITERS
    integer               :: nrepeats  = KMEANS_REPEATS
    logical               :: exists    = .false.
  contains
    procedure          :: new
    procedure, private :: spectral_embedding
    procedure, private :: kmeans
    procedure, private :: cpqr
    procedure, private :: calc_centroids
    procedure          :: cluster
    procedure          :: get_labels
    procedure          :: get_centers
    procedure          :: DunnIndex
    procedure          :: kill
end type spec_clust

contains

    ! CONSTRUCTOR

    !>  \brief  is a constructor
    subroutine new( self, N, K, S, algorithm, gamma )
        class(spec_clust),          intent(inout) :: self
        integer,                    intent(in)    :: N          ! similarity matrix size
        integer,                    intent(in)    :: K          ! number of desired clusters
        real,    target,            intent(inout) :: S(N,N)     ! similarity (=pairwise distance) matrix
        character(len=*), optional, intent(in)    :: algorithm
        real,             optional, intent(in)    :: gamma
        call self%kill
        self%S     => S
        self%N     = N
        self%K     = K
        self%neigs = K
        if( present(gamma) )then
            self%gamma = gamma
        else
            self%gamma = GAMMA_DEFAULT
        endif
        if( present(algorithm) )then
            select case(trim(algorithm))
            case('kmeans','cpqr')
                self%algorithm = trim(algorithm)
            case DEFAULT
                THROW_HARD('Unsupported algorithm: '//trim(algorithm))
            end select
        else
            self%algorithm = 'cpqr'
        endif
        self%nrepeats = KMEANS_REPEATS
        self%niters   = KMEANS_ITERS
        self%rt       = ran_tabu(self%N)
        allocate(self%labels(self%N),self%centroids(self%neigs,self%K),&
            &self%features(self%neigs,self%N))
        self%exists = .true.
    end subroutine new

    subroutine spectral_embedding( self )
        class(spec_clust), intent(inout) :: self
        real, allocatable :: Lap(:,:),eigvecs(:,:),eigvals(:),w(:)
        real    :: degree
        integer :: i
        allocate(Lap(self%N,self%N),eigvals(self%neigs),&
        &eigvecs(self%N,self%neigs),w(self%N))
        ! Laplacian
        !$omp parallel do private(i,degree) default(shared) proc_bind(close)
        do i = 1,self%N
            ! negative of Fully connected similarity graph
            Lap(:,i) = -exp(-self%gamma*self%S(:,i)**2)
            ! degree matrix diagonal element
            degree   = -sum(Lap(:,i))
            ! Laplacian
            Lap(i,i) = Lap(i,i) + degree
        enddo
        !$omp end parallel do
        ! Degree normalization
        do i = 1,self%N
            Lap(i,i) = 0.0
            w(i)     = sqrt(-sum(Lap(:,i)))
        enddo
        do i = 1,self%N
            Lap(:,i) = Lap(:,i) / w(i)
            Lap(i,:) = Lap(i,:) / w(i)
        enddo
        forall(i = 1:self%N) Lap(i,i) = 1.0
        ! Eigen decomposition
        call eigh_sp(self%N, Lap, self%neigs, eigvals, eigvecs, smallest=.true.)
        self%features = transpose(eigvecs)
        deallocate(Lap,eigvals,eigvecs,w)
    end subroutine spectral_embedding

    subroutine kmeans( self, K, F )
        class(spec_clust), intent(inout) :: self
        integer,           intent(in)    :: K      ! dynamic number of clusters
        integer,           intent(in)    :: F      ! dynamic number of features
        integer, allocatable :: labels(:), pops(:)
        real    :: dsq, sumdsq, best_sumdsq, dsqmin, prev_sumdsq
        integer :: repeat,iter, i,l
        allocate(labels(self%N),pops(K))
        ! repeats
        sumdsq      = 0.
        best_sumdsq = huge(dsq)
        do repeat = 1,self%nrepeats
            ! initial random labeling
            l = 0
            do i = 1,self%N
                l = l+1
                if( l > K ) l = 1
                labels(i) = l
            enddo
            call self%rt%shuffle(labels)
            ! one run of k-means
            prev_sumdsq = huge(dsq)
            do iter = 1,self%niters
                ! centroids
                self%centroids = 0.
                pops = 0
                do i = 1,self%N
                    l = labels(i)
                    self%centroids(:F,l) = self%centroids(:F,l) + self%features(:F,i)
                    pops(l) = pops(l)+1
                enddo
                do l =1,K
                    if( pops(l) > 1 )then
                        self%centroids(:F,l) = self%centroids(:F,l) / real(pops(l))
                    endif
                enddo
                ! greedy assignment
                sumdsq = 0.
                do i = 1,self%N
                    dsqmin = huge(dsq)
                    do l = 1,K
                        if( pops(l) == 0 ) cycle
                        ! dsq = sum((self%features(:F,i)-self%centroids(:F,l))**2)
                        dsq = sum(abs(self%features(:F,i)-self%centroids(:F,l)))
                        if( dsq < dsqmin )then
                            labels(i) = l
                            dsqmin    = dsq
                        endif
                    enddo
                    sumdsq = sumdsq + dsqmin
                enddo
                ! convergence
                if( prev_sumdsq - sumdsq < 1.e-6 ) exit
                prev_sumdsq = sumdsq
            enddo
            ! lowest wins
            if( sumdsq < best_sumdsq )then
                best_sumdsq = sumdsq
                self%labels = labels
            endif
        enddo
        deallocate(labels, pops)
    end subroutine kmeans

    subroutine CPQR( self )
        class(spec_clust), intent(inout) :: self
        real,    allocatable :: A(:,:), V(:,:),vec(:)
        integer, allocatable :: jpvt(:)
        integer :: i
        ! transpose of selected eigenvectors matrix
        allocate(A(self%neigs,self%N), source=self%features)
        ! QR-factorization
        allocate(jpvt(self%N))
        call qr(self%neigs, self%N, A, self%neigs, jpvt, vec)
        A = transpose(self%features(:,jpvt(1:self%neigs)))
        deallocate(vec,jpvt)
        ! Eigen problem
        allocate(vec(self%neigs), V(self%neigs,self%neigs))
        call svdcmp(A, vec, V)
        ! Labelling
        A = abs(matmul(transpose(self%features),matmul(A,V)))
        forall(i = 1:self%N) self%labels(i) = maxloc(A(i,:),dim=1)
        deallocate(A,V,vec)
    end subroutine CPQR

    subroutine calc_centroids( self )
        class(spec_clust), intent(inout) :: self
        integer :: pops(self%K),i,l
        self%centroids = 0.
        pops = 0
        do i = 1,self%N
            l = self%labels(i)
            self%centroids(:,l) = self%centroids(:,l) + self%features(:,i)
            pops(l) = pops(l) + 1
        enddo
        do l = 1,self%K
            if( pops(l) > 1 )then
                self%centroids(:,l) = self%centroids(:,l) / real(pops(l))
            endif
        enddo
    end subroutine calc_centroids

    subroutine cluster( self )
        class(spec_clust), intent(inout) :: self
        call self%spectral_embedding
        select case(trim(self%algorithm))
        case('kmeans')
            call self%kmeans(self%K, self%neigs)
        case('cpqr')
            call self%cpqr
        end select
        call self%calc_centroids
    end subroutine cluster

    subroutine get_labels(self, l)
        class(spec_clust),    intent(in)    :: self
        integer, allocatable, intent(inout) :: l(:)
        if( allocated(l) ) deallocate(l)
        allocate(l,source=self%labels)
    end subroutine get_labels

    subroutine get_centers(self, c)
        class(spec_clust), intent(in)    :: self
        real, allocatable, intent(inout) :: c(:,:)
        if( allocated(c) ) deallocate(c)
        allocate(c,source=self%centroids)
    end subroutine get_centers

    real function DunnIndex( self )
        class(spec_clust), intent(in) :: self
        real    :: distmat(self%K,self%K), mininter, maxintra
        integer :: distpops(self%K,self%K)
        integer :: i,j,li,lj
        ! intra/inter cluster mean distance
        distpops = 0
        distmat  = 0.
        do i = 1,self%N
            li = self%labels(i)
            do j = i+1,self%N
                lj = self%labels(j)
                distmat(li,lj)  = distmat(li,lj) + self%S(i,j)
                distpops(li,lj) = distpops(li,lj) + 1
            enddo
        enddo
        where( distpops > 1 ) distmat = distmat / real(distpops)
        mininter = huge(mininter)
        maxintra = -1.
        do i=1,self%K
            do j=i,self%K
                if( i==j )then
                    maxintra = max(distmat(i,j), maxintra)
                else
                    mininter = min(distmat(i,j), mininter)
                endif
            enddo
        enddo
        DunnIndex = mininter / maxintra
    end function DunnIndex

    ! UNIT TEST

    !>  \brief  is the spec_clust unit test
    subroutine test_spec_clust
        write(logfhandle,'(a)') 'TESTING SPECTRAL CLUSTERING/KMEANS'
        call test('kmeans')
        write(logfhandle,'(a)') 'TESTING SPECTRAL CLUSTERING/CPQR'
        call test('cpqr')
      contains

        subroutine test( algorithm )
            character(len=*), intent(in) :: algorithm
            type(spec_clust)     :: spc
            integer, allocatable :: labels(:)
            real                 :: datavecs(90,5), simmat(90,90)
            integer              :: i, j, ncls, nerr
            write(logfhandle,'(a)') '**info(simple_spec_clust_unit_test): testing all functionality'
            ! make data
            do i=1,30
                datavecs(i,:) = 1.
            end do
            do i=31,60
                datavecs(i,:) = 5.
            end do
            do i=61,90
                datavecs(i,:) = 10.
            end do
            do i=1,90
                simmat(i,i) = 0.
                do j=i+1,90
                    simmat(i,j) = euclid(datavecs(i,:),datavecs(j,:))**2
                    simmat(j,i) = simmat(i,j)
                end do
            end do
            call spc%new(90, 3, simmat, algorithm=algorithm)
            call spc%cluster
            call spc%get_labels(labels)
            ncls = maxval(labels)
            nerr = 0
            do i=1,29
                do j=i+1,30
                    if( labels(i) /= labels(j) ) nerr = nerr+1
                end do
            end do
            do i=31,59
                do j=i+1,60
                    if( labels(i) /= labels(j) ) nerr = nerr+1
                end do
            end do
            do i=61,89
                do j=i+1,90
                    if( labels(i) /= labels(j) ) nerr = nerr+1
                end do
            end do
            write(logfhandle,*) 'NR OF CLUSTERS FOUND:', ncls
            write(logfhandle,*) 'NR OF ASSIGNMENT ERRORS:', nerr
            if( ncls == 3 .and. nerr == 0 )then
                write(logfhandle,'(a)') 'SPECTRAL_CLUSTERING UNIT TEST COMPLETED ;-)'
            else
                write(logfhandle,'(a)') 'SPECTRAL_CLUSTERING UNIT TEST FAILED!'
            endif
        end subroutine test

    end subroutine test_spec_clust

    ! DESTRUCTOR

    subroutine kill( self )
        class(spec_clust), intent(inout) :: self
        if( self%exists )then
            nullify(self%S)
            deallocate(self%centroids,self%features,self%labels)
            call self%rt%kill
            self%N         = 0
            self%K         = 0
            self%neigs     = 0
            self%algorithm = 'cpqr'
            self%gamma     = GAMMA_DEFAULT
            self%niters    = KMEANS_ITERS
            self%nrepeats  = KMEANS_REPEATS
            self%exists    = .false.
        endif
    end subroutine kill

end module simple_spectral_clustering
