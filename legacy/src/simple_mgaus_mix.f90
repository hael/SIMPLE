module simple_mgaus_mix
use simple_mgaus
use simple_jiffys
use simple_rnd
use simple_params
implicit none

! to generate a sample from a mixture of Gaussians, we first chose on of the components at random with probability given
! by the mixing coefficients and then generate a sample vector x from the corresponding Gaussian component. This process
! is repeated N times to genrate a data set of N independent samples

! Usage:
! (1) make singleton object (make_mgaus_mix)
! (2) set all multivariate Gaussians and mixing coefficients (set_mgaus_mix)
! (3) expectation maximization (em_mgaus_mix)
! (4) use the mgaus_mix_gammas to form weighted averages

public :: make_mgaus_mix, kill_mgaus_mix, set_mgaus_mix, em_mgaus_mix, mgaus_mix_gammas
private

type(mgaus), allocatable :: mgs(:)                ! Multivariate Gaussians
real, allocatable        :: mixc(:)               ! mixing coefficients
real, allocatable        :: mgaus_mix_gammas(:,:) ! responsibilities
real, pointer            :: x(:,:)                ! data
real                     :: p                     ! log likelihood
integer                  :: N                     ! nr of contributing data vectors 
integer                  :: K                     ! nr of Gaussian components (=nr of clusters)
integer                  :: D                     ! dimension of contributing data vectors

contains

    subroutine make_mgaus_mix( x_in, N_in, D_in, K_in )
        integer, intent(in)      :: N_in, D_in, K_in
        real, intent(in), target :: x_in(N_in,D_in)
        integer                  :: alloc_stat, i
        N = N_in
        D = D_in
        K = K_in
        allocate( mgs(K), mixc(K), mgaus_mix_gammas(N,K), stat=alloc_stat )
        call alloc_err('In: new_mgaus_mix, module: simple_mgaus_mix', alloc_stat)
        do i=1,K
            mgs(i) = new_mgaus(D)
        end do
        x => x_in
    end subroutine make_mgaus_mix
    
    subroutine kill_mgaus_mix
        integer :: i
        do i=1,K
            call kill_mgaus(mgs(i))
        end do
        deallocate(mgs, mixc)
    end subroutine kill_mgaus_mix
    
    subroutine set_mgaus_mix( i, avg, cov, fracpop )
        integer, intent(in) :: i
        real, intent(in)    :: avg(D), cov(D,D), fracpop
        call set_mgaus( mgs(i), avg, cov, epsilon )
        mixc(i) = fracpop
    end subroutine set_mgaus_mix
    
    subroutine em_mgaus_mix( eps )
    ! if problems:
    ! (1) detect collapse of a Gaussian (var < some value)
    ! (2) reset its mean to some (random???) value
    ! (3) reset the covariance to some large value
        real, intent(in) :: eps
        real             :: prev, s
        integer          :: i, it
        prev = -0.
        p = -100.
        it = 0
        write(*,'(A)') '>>> EM OPTIMIZATION OF GAUSSIAN MIXTURE'
        do while( abs(p-prev) > 0.5*eps )
            it = it+1
            write(*,*) 'eval gammas'
            call eval_gammas
            write(*,*) 're-est params'
            call re_est_params
            prev = p
            write(*,*) 'eval lhood'
            call eval_lhood
            write(*,*) 'Iteration:', it, 'log-likelihood:', p
        end do
    end subroutine em_mgaus_mix
    
    subroutine eval_gammas
    ! E-step
        real    :: s
        integer :: i, j
        do j=1,N
            s = 0.
            do i=1,K
                mgaus_mix_gammas(j,i) = mixc(i)*sample_mgaus(mgs(i),x(j,:))
                s = s+mgaus_mix_gammas(j,i)
            end do
            do i=1,K
                mgaus_mix_gammas(j,i) = mgaus_mix_gammas(j,i)/s
            end do
        end do
    end subroutine eval_gammas
        
    subroutine eval_lhood
    ! evaluate the log likelihood
        integer :: i, j
        real    :: lsum
        p = 0.
        do j=1,N
            lsum = 0.
            do i=1,K
                lsum = lsum+mixc(i)*sample_mgaus(mgs(i),x(j,:))
            end do
            p = p+log(lsum)
        end do
    end subroutine eval_lhood
    
    subroutine re_est_params
    ! M-step
        integer :: i, j
        real    :: nk, avg(D), cov(D,D)
        real    :: diffvec(D,1), nsum, epsilon
        nsum = sum(mgaus_mix_gammas)
        do i=1,K
            ! re-estimate average
            nk = 0.
            avg = 0.
            do j=1,N
                avg = avg+mgaus_mix_gammas(j,i)*x(j,:)
                nk = nk+mgaus_mix_gammas(j,i)
            end do
            avg = avg/nk
            ! re-estimate covariance matrix
            cov = 0.
            do j=1,N
                diffvec(:,1) = x(j,:)-avg
                cov = cov+mgaus_mix_gammas(j,i)*matmul(diffvec,transpose(diffvec))
            end do
            cov = cov/nk
            ! re-estimate mixing coefficients
            mixc(i) = nk/nsum
            ! set the responsible object
            call set_mgaus( mgs(i), avg, cov, epsilon )
        end do
    end subroutine re_est_params
    
end module simple_mgaus_mix