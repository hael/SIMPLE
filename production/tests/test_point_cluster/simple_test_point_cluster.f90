program simple_test_point_cluster
include 'simple_lib.f08'
implicit none
integer, parameter :: N_CLUSTER = 5, N_POINTS = 50
real,    parameter :: EPS = 0.2
integer :: i, j, iter, truth, cur_clusN, prev_clusN, inds(N_POINTS), cur_cluster(N_POINTS), prev_cluster(N_POINTS)
real    :: points(N_POINTS), avg_points(N_POINTS), costs(N_POINTS)
logical :: taken(N_POINTS)
call seed_rnd
do i = 1, N_POINTS
    truth     = floor(ran3() * real(N_CLUSTER)) + 1
    points(i) = truth + ran3() * EPS
enddo
print *, 'Simulated points:'
print *, points
! recursive clustering
avg_points = points
cur_clusN  = N_POINTS
do iter = 1, 10
    prev_clusN   = cur_clusN
    points       = avg_points(1:prev_clusN)
    cur_clusN    = 0
    taken        = .false.
    avg_points   = 0.
    prev_cluster = (/(j,j=1,N_POINTS)/)
    do i = 1, prev_clusN
        if( taken(i) ) cycle
        taken(i)  = .true.
        cur_clusN = cur_clusN + 1
        costs     = huge(points(1))
        do j = 1, prev_clusN
            if( taken(j) .or. j == i ) cycle
            costs(j) = abs(points(j) - points(i))
        enddo
        inds = (/(j,j=1,prev_clusN)/)
        call hpsort(costs, inds)
        ! check if the best cost is within EPS
        if( costs(1) < EPS )then
            taken(inds(1)) = .true.
            avg_points(cur_clusN) = (points(i) + points(inds(1))) / 2.
        else
            avg_points(cur_clusN) = points(i)
        endif
    enddo
    print *, 'cur_clusN  = ', cur_clusN
    print *, 'avg_points = '
    print *, avg_points(1:cur_clusN)
    if( prev_clusN == cur_clusN ) exit
enddo

end program simple_test_point_cluster
