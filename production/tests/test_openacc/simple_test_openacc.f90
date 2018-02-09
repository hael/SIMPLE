program simple_test_openacc
use openacc
implicit none
integer, parameter :: N=10, M=500, P=30
integer :: i, j, k
real    :: particles(N,M,P), total

particles = 1.0

!$acc kernels
do i=1,N
    do j=1,M
        do k=1,P
            particles(i,j,k) = particles(i,j,k) * 2.0
        end do
    end do
end do
!$acc end kernels

print *, 'sum should be: ', N*M*P*2, 'sum is: ', sum(particles)

!$acc parallel loop
do i=1,N
    do j=1,M
        do k=1,P
            particles(i,j,k) = particles(i,j,k) * 2.0
        end do
    end do
end do
!$acc end parallel loop

print *, 'sum should be: ', N*M*P*4, 'sum is: ', sum(particles)

!$acc parallel loop num_gangs(32) vector_length(128)
do i=1,N
    do j=1,M
        do k=1,P
            particles(i,j,k) = particles(i,j,k) * 2.0
        end do
    end do
end do
!$acc end parallel loop

print *, 'sum should be: ', N*M*P*8, 'sum is: ', sum(particles)

!$acc parallel loop collapse(3)
do i=1,N
    do j=1,M
        do k=1,P
            particles(i,j,k) = particles(i,j,k) * 2.0
        end do
    end do
end do
!$acc end parallel loop

print *, 'sum should be: ', N*M*P*16, 'sum is: ', sum(particles)

particles = 1.0

total = 0.
!$acc parallel loop collapse(3) reduction(+:total)
do i=1,N
    do j=1,M
        do k=1,P
            particles(i,j,k) = particles(i,j,k) * 2.0
            total = total + particles(i,j,k)
        end do
    end do
end do
!$acc end parallel loop

print *, 'sum should be: ', N*M*P*2, 'sum is: ', total

end program simple_test_openacc
