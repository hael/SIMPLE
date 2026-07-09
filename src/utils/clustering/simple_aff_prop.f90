!@descr: clustering of a similartity matrix using affinity propagation
module simple_aff_prop
use simple_core_module_api
implicit none

public :: aff_prop, test_aff_prop
private
#include "simple_local_flags.inc"

type aff_prop
    private
    integer              :: N                    !< nr of data entries
    real,    allocatable :: A(:,:), R(:,:)       !< A is affinities & R responsibilities
    real,    allocatable :: Aold(:,:), Rold(:,:)
    real,    allocatable :: Rp(:,:), tmp(:)
    real,    allocatable :: AS(:,:), dA(:)
    real,    allocatable :: S(:,:)               !< private similarity matrix
    real,    allocatable :: Y(:), Y2(:)          !< maxvals
    integer, allocatable :: I(:), I2(:)          !< index arrays
    integer              :: maxits=500           !< maximum number of iterations
    real                 :: ftol=1e-9            !< fractional convergence tolerance score
    real                 :: lam=0.5              !< dampening factor
    real                 :: Smin, Smax           !< similarity min/max
    logical              :: exists=.false.       !< to indicate existence
  contains
    procedure :: new
    procedure :: propagate
    procedure :: kill
end type aff_prop

contains

    ! CONSTRUCTOR

    !>  \brief  is a constructor
    subroutine new( self, N, S, ftol, lam, pref, maxits )
        class(aff_prop),   intent(inout) :: self
        integer,           intent(in)    :: N               ! # data entries
        real,              intent(in)    :: S(N,N)          ! similarity matrix
        real,    optional, intent(in)    :: ftol, lam, pref ! tolerance, dampening factor, preference thres
        integer, optional, intent(in)    :: maxits
        integer :: i, j
        real    :: ppref, diff
        call self%kill
        ! set constants
        self%N = N
        if( present(ftol)    ) self%ftol    = ftol
        if( present(lam)     ) self%lam     = lam
        if( present(maxits)  ) self%maxits  = maxits
        ! calculate similarity characteristics
        allocate(self%S(N,N), source=S)
        call analyze_smat(self%S, .false., self%Smin, self%Smax)
        ! low pref [0.1,1] leads to small numbers of clusters and
        ! high pref leads to large numbers of clusters. Frey suggests Smin - (Smax - Smin) as a threshold for fewer cluster
        ppref = self%Smin
        if( present(pref) ) ppref = pref
        ! remove degeneracies reproducibly; volcluster/AP restarts must be deterministic
        do i=1,N
            self%S(i,i) = ppref
        end do
        diff = (self%Smax - self%Smin) * self%ftol
        do i=1,self%N
            do j=1,self%N
                self%S(i,j) = self%S(i,j) + deterministic_tiebreak(i, j) * diff
            end do
        end do
        ! allocate
        allocate( self%A(N,N), self%R(N,N), self%Aold(N,N), self%Rp(N,N),&
        self%Rold(N,N), self%AS(N,N), self%Y(N), self%Y2(N), self%tmp(N),&
        self%I(N), self%I2(N), self%dA(N) )
        self%exists = .true.
    end subroutine new

    ! PROPAGATOR

    !>  \brief  is the message passing algorithm
    subroutine propagate( self, centers, labels, simsum )
        class(aff_prop),      intent(inout) :: self
        integer, allocatable, intent(inout) :: centers(:) !< cluster centres
        integer, allocatable, intent(inout) :: labels(:)  !< cluster labels
        real,                 intent(inout) :: simsum     !< similarity sum
        real, allocatable :: similarities(:)
        real              :: x, realmax, maxdiff, rowmax, val, best, second, exemplar_tol, diag_absmax
        integer           :: i, j, k, ncls, nborder, fallback_center
        logical           :: converged, use_fallback_center
        ! initialize
        realmax   = huge(x)
        self%A    = 0.
        self%R    = 0.
        self%Aold = 0.
        self%Rp   = 0.
        self%Rold = 0.
        self%AS   = 0.
        self%Y    = 0.
        self%Y2   = 0.
        self%tmp  = 0.
        self%I    = 0
        self%I2   = 0
        self%dA   = 0.
        ! iterate
        do i=1,self%maxits
            !------------------------
            ! RESPONSIBILITIES
            !------------------------
            !$omp parallel do default(shared) private(j,k) schedule(static)
            do k = 1, self%N
                do j = 1, self%N
                    self%Rold(j,k) = self%R(j,k)
                    self%AS(j,k)   = self%A(j,k) + self%S(j,k)
                end do
            end do
            !$omp end parallel do
            !$omp parallel do default(shared) private(j,k,best,second,val) schedule(static)
            do j = 1, self%N
                self%I(j)  = 1
                self%I2(j) = 1
                best       = self%AS(j,1)
                second     = -realmax
                do k = 2, self%N
                    val = self%AS(j,k)
                    if( val > best )then
                        second     = best
                        self%I2(j) = self%I(j)
                        best       = val
                        self%I(j)  = k
                    else if( val > second )then
                        second     = val
                        self%I2(j) = k
                    endif
                end do
                self%Y(j)  = best
                self%Y2(j) = second
            end do
            !$omp end parallel do
            !$omp parallel do default(shared) private(j,k) schedule(static)
            do j = 1, self%N
                do k = 1, self%N
                    self%R(j,k) = self%S(j,k) - self%Y(k)
                end do
            enddo
            !$omp end parallel do
            !$omp parallel do default(shared) private(j) schedule(static)
            do j=1,self%N
                self%R(j,self%I(j)) = self%S(j,self%I(j)) - self%Y2(j)
            end do
            !$omp end parallel do
            ! update responsibilities (in a dampened fashion)
            !$omp parallel do default(shared) private(j,k) schedule(static)
            do k = 1, self%N
                do j = 1, self%N
                    self%R(j,k) = (1. - self%lam) * self%R(j,k) + self%lam * self%Rold(j,k)
                end do
            end do
            !$omp end parallel do
            !------------------------
            ! AVAILABILITIES
            !------------------------
            !$omp parallel do default(shared) private(j,k) schedule(static)
            do k = 1, self%N
                do j = 1, self%N
                    self%Aold(j,k) = self%A(j,k)
                    if( self%R(j,k) > 0. )then
                        self%Rp(j,k) = self%R(j,k)
                    else
                        self%Rp(j,k) = 0.
                    endif
                end do
            end do
            !$omp end parallel do
            !$omp parallel do default(shared) private(k) schedule(static)
            do k = 1, self%N
                self%Rp(k,k) = self%R(k,k)
            end do
            !$omp end parallel do
            !$omp parallel do default(shared) private(j,k,val) schedule(static)
            do k = 1, self%N
                val = 0.
                do j = 1, self%N
                    val = val + self%Rp(j,k)
                end do
                self%tmp(k) = val
            end do
            !$omp end parallel do
            !$omp parallel do default(shared) private(j,k) schedule(static)
            do k = 1, self%N
                do j = 1, self%N
                    self%A(j,k) = -self%Rp(j,k) + self%tmp(k)
                end do
            end do
            !$omp end parallel do
            !$omp parallel do default(shared) private(k) schedule(static)
            do k = 1, self%N
                self%dA(k) = self%A(k,k)
            end do
            !$omp end parallel do
            !$omp parallel do default(shared) private(j,k) schedule(static)
            do k = 1, self%N
                do j = 1, self%N
                    if( self%A(j,k) > 0. ) self%A(j,k) = 0.
                end do
            end do
            !$omp end parallel do
            !$omp parallel do default(shared) private(k) schedule(static)
            do k = 1, self%N
                self%A(k,k) = self%dA(k)
            end do
            !$omp end parallel do
            ! update availabilities (in a dampened fashion)
            !$omp parallel do default(shared) private(j,k) schedule(static)
            do k = 1, self%N
                do j = 1, self%N
                    self%A(j,k) = (1. - self%lam) * self%A(j,k) + self%lam * self%Aold(j,k)
                end do
            end do
            !$omp end parallel do
            !========================
            ! Convergence check each 5th iter
            !========================
            if (mod(i,5) == 0) then
                maxdiff = 0.0
                !$omp parallel do default(shared) private(j,k,rowmax,val) schedule(static)
                do j = 1, self%N
                    rowmax = 0.
                    do k = 1, self%N
                        val = max( abs(self%A(j,k)-self%Aold(j,k)), abs(self%R(j,k)-self%Rold(j,k)) )
                        if (val > rowmax) rowmax = val
                    end do
                    self%tmp(j) = rowmax
                end do
                !$omp end parallel do
                do j = 1, self%N
                    if( self%tmp(j) > maxdiff ) maxdiff = self%tmp(j)
                end do
                if (maxdiff < self%ftol) then
                    write(logfhandle,'(a,i6,1x,es12.4)') 'aff_prop converged at iter=', i, maxdiff
                    converged = .true.
                    exit
                end if
            endif
        end do
        !$omp parallel do default(shared) private(j,k) schedule(static)
        do k = 1, self%N
            do j = 1, self%N
                self%R(j,k) = self%R(j,k) + self%A(j,k) ! pseudomarginals
            end do
        end do
        !$omp end parallel do
        diag_absmax = 0.
        do j=1,self%N
            diag_absmax = max(diag_absmax, abs(self%R(j,j)))
        end do
        exemplar_tol = max(10. * epsilon(exemplar_tol), self%ftol) * max(1., diag_absmax)
        nborder = 0
        fallback_center = 1
        ! count the number of clusters
        ncls = 0
        do j=1,self%N
            if( self%R(j,j) > self%R(fallback_center,fallback_center) ) fallback_center = j
            if( self%R(j,j) > exemplar_tol )then
                ncls = ncls + 1
            else if( abs(self%R(j,j)) <= exemplar_tol )then
                nborder = nborder + 1
            endif
        end do
        use_fallback_center = ncls == 0
        if( use_fallback_center ) ncls = 1
        if( nborder > 0 ) write(logfhandle,'(a,i0,1x,es12.4)') 'aff_prop near-zero exemplar candidates ignored: ', nborder, exemplar_tol
        if( allocated(centers) ) deallocate(centers)
        if( allocated(labels) )  deallocate(labels)
        allocate( centers(ncls), similarities(ncls), labels(self%N) )
        ! set the cluster centers
        if( use_fallback_center )then
            centers(1) = fallback_center
        else
            ncls = 0
            do j=1,self%N
                if( self%R(j,j) > exemplar_tol )then
                    ncls = ncls + 1
                    centers(ncls) = j
                endif
            end do
        endif
        ! report back the labeling
        simsum = 0.
        do j=1,self%N
            do k=1,ncls
                if( j .ne. centers(k) )then
                    similarities(k) = self%S(centers(k),j)
                else
                    similarities(k) = realmax
                endif
            end do
            labels(j) = 1
            do k=2,ncls
                if( similarities(k) > similarities(labels(j)) ) labels(j) = k
            end do
            if( j .ne. centers(labels(j)) ) simsum = simsum + similarities(labels(j))
        end do
        simsum = simsum / real(self%N)
        deallocate(similarities)
    end subroutine propagate

    pure real function deterministic_tiebreak( i, j ) result( val )
        integer, intent(in) :: i, j
        integer, parameter :: HASH_MOD = 8191, HASH_A = 37, HASH_B = 1009
        integer :: ii, jj, h
        ii  = mod(i - 1, HASH_MOD)
        jj  = mod(j - 1, HASH_MOD)
        h   = mod(HASH_A * ii + HASH_B * jj + mod(ii * jj, HASH_MOD), HASH_MOD)
        val = real(h) / real(HASH_MOD)
    end function deterministic_tiebreak

    ! UNIT TEST

    !>  \brief  is the aff_prop unit test
    subroutine test_aff_prop
        real,    allocatable :: datavecs(:,:)
        type(aff_prop)       :: apcls
        real,    allocatable :: simmat(:,:), simmat_ref(:,:)
        real                 :: simsum, simsum2
        integer, allocatable :: centers(:), centers2(:), labels(:), labels2(:)
        integer              :: i, j, ncls, nerr, ndet, nper, ntot
        write(logfhandle,'(a)') '**info(simple_aff_prop_unit_test): testing all functionality'
#if defined(_WIN32)
        nper = 40
#else
        nper = 300
#endif
        ntot = 3 * nper
        allocate(datavecs(ntot,5), simmat(ntot,ntot), simmat_ref(ntot,ntot))
        ! make data
        do i=1,nper
            datavecs(i,:) = 1.
        end do
        do i=nper+1,2*nper
            datavecs(i,:) = 5.
        end do
        do i=2*nper+1,ntot
            datavecs(i,:) = 10.
        end do
        do i=1,ntot-1
            simmat(i,i) = 0.
            do j=i+1,ntot
                simmat(i,j) = -euclid(datavecs(i,:),datavecs(j,:))
                simmat(j,i) = simmat(i,j)
            end do
        end do
        simmat(ntot,ntot) = 0.
        simmat_ref = simmat
        call apcls%new(ntot, simmat)
        call apcls%propagate(centers, labels, simsum)
        call apcls%new(ntot, simmat_ref)
        call apcls%propagate(centers2, labels2, simsum2)
        ncls = size(centers)
        nerr = 0
        ndet = 0
        if( size(centers2) /= size(centers) )then
            ndet = ndet + 1
        else
            if( any(centers2 /= centers) ) ndet = ndet + 1
        endif
        if( size(labels2) /= size(labels) )then
            ndet = ndet + 1
        else
            if( any(labels2 /= labels) ) ndet = ndet + 1
        endif
        do i=1,nper-1
            do j=i+1,nper
                if( labels(i) /= labels(j) ) nerr = nerr+1
            end do
        end do
        do i=nper+1,2*nper-1
            do j=i+1,2*nper
                if( labels(i) /= labels(j) ) nerr = nerr+1
            end do
        end do
        do i=2*nper+1,ntot-1
            do j=i+1,ntot
                if( labels(i) /= labels(j) ) nerr = nerr+1
            end do
        end do
        write(logfhandle,*) 'NR OF CLUSTERS FOUND:', ncls
        write(logfhandle,*) 'NR OF ASSIGNMENT ERRORS:', nerr
        write(logfhandle,*) 'NR OF DETERMINISM ERRORS:', ndet
        write(logfhandle,*) 'CENTERS'
        do i=1,size(centers)
            write(logfhandle,*) datavecs(centers(i),:)
        end do
        if( ncls == 3 .and. nerr == 0 .and. ndet == 0 )then
            write(logfhandle,'(a)') 'SIMPLE_AFF_PROP_UNIT_TEST COMPLETED ;-)'
        else
            write(logfhandle,'(a)') 'SIMPLE_AFF_PROP_UNIT_TEST FAILED!'
        endif
        call apcls%kill
    end subroutine test_aff_prop

    ! DESTRUCTOR

    subroutine kill( self )
        class(aff_prop), intent(inout) :: self
        if( self%exists )then
            deallocate( self%S, self%A, self%R, self%Aold, self%Rp,&
            self%Rold, self%AS, self%Y, self%Y2, self%tmp, self%I,&
            self%I2, self%dA )
            self%exists = .false.
        endif
    end subroutine kill

end module simple_aff_prop
