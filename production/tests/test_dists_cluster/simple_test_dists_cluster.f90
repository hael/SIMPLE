program simple_test_dists_cluster
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_parameters,         only: parameters
use simple_nanoparticle_utils, only: read_pdb2matrix
implicit none
real,    allocatable :: dists(:), mat(:,:), dists_cen(:)
integer, allocatable :: cnts(:)
type(parameters)     :: p
type(cmdline)        :: cline
integer              :: Natoms, i, j, cnt, nclus, sum_cnt
real                 :: min_dist, eps
! reading pdb file
if( command_argument_count() < 2 )then
    write(logfhandle,'(a)') 'Usage: simple_test_dists_cluster nthr=yy pdbfile=zz'
    stop
endif
call cline%parse_oldschool
call cline%checkvar('nthr',     1)
call cline%checkvar('pdbfile' , 2)
call cline%check
call p%new(cline)
call read_pdb2matrix( p%pdbfile, mat )
Natoms = size(mat,2)
allocate(dists(Natoms**2), source=0.)
cnt      = 0
min_dist = huge(min_dist)
do i = 1, Natoms
    do j = 1, Natoms
        cnt        = cnt + 1
        dists(cnt) = sqrt(sum((mat(:,i) - mat(:,j))**2))
        if( dists(cnt) < min_dist .and. dists(cnt) > TINY ) min_dist = dists(cnt)
    enddo
enddo
eps = min_dist / 2.
call sort_cluster(dists, nclus, dists_cen, cnts)
sum_cnt = 0
do i = 1, nclus
    print *, 'cluster ', i, ' with center = ', dists_cen(i), ' has ', cnts(i), ' points'
    sum_cnt = sum_cnt + cnts(i)
enddo
print *, 'Natoms  = ', Natoms
print *, 'sum_cnt = ', sum_cnt, ' while Natoms**2 = ', Natoms**2

contains

    subroutine sort_cluster( points_in, cur_clusN, out_points, out_cnts )
        real,                 intent(in)    :: points_in(:)
        integer,              intent(inout) :: cur_clusN
        real,    allocatable, intent(inout) :: out_points(:)
        integer, allocatable, intent(inout) :: out_cnts(:)
        real,    allocatable :: points(:), avg_points(:)
        integer, allocatable :: cnts(:)
        integer :: cur_ind, n_points
        real    :: cur_avg
        n_points = size(points_in)
        allocate(points(n_points), source=points_in)
        allocate(avg_points(n_points))
        allocate(cnts(n_points))
        call hpsort(points)
        cur_clusN = 0
        cur_ind   = 1
        cur_avg   = 0.
        do i = 1, n_points - 1
            cur_avg = cur_avg+points(i)
            if( abs(points(i+1) - cur_avg/real(i-cur_ind+1)) > eps )then
                cur_clusN = cur_clusN + 1
                avg_points(cur_clusN) = cur_avg/real(i-cur_ind+1)
                cnts(cur_clusN)       = i-cur_ind+1
                cur_ind = i + 1
                cur_avg = 0.
            endif
        enddo
        ! last cluster
        cur_clusN = cur_clusN + 1
        avg_points(cur_clusN) = sum(points(cur_ind:i))/real(i-cur_ind+1)
        cnts(cur_clusN)       = i-cur_ind+1
        allocate(out_points(1:cur_clusN), source=avg_points(1:cur_clusN))
        allocate(out_cnts(1:cur_clusN),   source=cnts(1:cur_clusN))
    end subroutine sort_cluster

    subroutine point_cluster( points_in, cur_clusN, out_points )
        real,                 intent(in)    :: points_in(:)
        integer,              intent(inout) :: cur_clusN
        real,    allocatable, intent(inout) :: out_points(:)
        real,    allocatable :: avg_points(:), costs(:), points(:)
        logical, allocatable :: taken(:)
        integer :: i, j, iter, prev_clusN, n_points
        n_points = size(points_in)
        allocate(points(n_points), source=points_in)
        allocate(avg_points(n_points), costs(n_points), taken(n_points))
        ! recursive clustering
        avg_points = points
        cur_clusN  = n_points
        do iter = 1, 30
            prev_clusN = cur_clusN
            points     = avg_points(1:prev_clusN)
            cur_clusN  = 0
            taken      = .false.
            avg_points = 0.
            do i = 1, prev_clusN
                if( taken(i) ) cycle
                taken(i)  = .true.
                cur_clusN = cur_clusN + 1
                costs     = huge(points(1))
                do j = 1, prev_clusN
                    if( taken(j) .or. j == i ) cycle
                    costs(j) = abs(points(j) - points(i))
                enddo
                j = minloc(costs, dim=1)
                ! check if the best cost is within eps
                if( costs(j) < eps )then
                    taken(j) = .true.
                    avg_points(cur_clusN) = (points(i) + points(j)) / 2.
                else
                    avg_points(cur_clusN) = points(i)
                endif
            enddo
            ! print *, 'cur_clusN  = ', cur_clusN
            ! print *, 'avg_points = '
            ! print *, avg_points(1:cur_clusN)
            if( prev_clusN == cur_clusN ) exit
        enddo
        allocate(out_points(cur_clusN), source=avg_points(1:cur_clusN))
    end subroutine point_cluster

end program simple_test_dists_cluster
