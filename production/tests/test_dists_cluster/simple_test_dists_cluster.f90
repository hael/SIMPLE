program simple_test_dists_cluster
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_parameters,         only: parameters
use simple_nanoparticle_utils, only: read_pdb2matrix
implicit none
character(len=LONGSTRLEN), allocatable :: pdbfnames(:)
real,    allocatable :: dists(:), mat(:,:), dists_cen(:), vars(:), centers(:), ref_stats(:), cur_mat(:,:), cur_stats(:), probs(:)
integer, allocatable :: cnts(:), cen_cnts(:)
type(parameters)     :: p
type(cmdline)        :: cline
integer              :: Natoms, i, j, cnt, nclus, sum_cnt, cen_cnt, npdbs, ipdb, ithres, stab_cnt
real                 :: min_dist, eps, dist, d, prob, dist_thres, prev_prob, stab_prob
! reading pdb file
if( command_argument_count() < 3 )then
    write(logfhandle,'(a)') 'Usage: simple_test_dists_cluster nthr=yy pdbfile=zz pdbfiles=tt'
    stop
endif
call cline%parse_oldschool
call cline%checkvar('nthr',     1)
call cline%checkvar('pdbfile' , 2)
call cline%checkvar('pdbfiles', 3)
call cline%check
call p%new(cline)
dist_thres = 2.
do ithres=1,8
    call read_pdb2matrix( p%pdbfile, mat )
    Natoms = size(mat,2)
    if( allocated(dists) )deallocate(dists)
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
    eps = min_dist / 10.
    call sort_cluster(dists, nclus, dists_cen, cnts, vars)
    dist_thres = dist_thres + 1.
    sum_cnt = 0
    cen_cnt = 0
    do i = 1, nclus
        if( dists_cen(i) > dist_thres )cycle
        sum_cnt = sum_cnt + cnts(i)
        cen_cnt = cen_cnt + 1
    enddo
    if( allocated(centers) )  deallocate(centers)
    if( allocated(cen_cnts) ) deallocate(cen_cnts)
    allocate(centers(cen_cnt),source=dists_cen(1:cen_cnt))
    allocate(cen_cnts(cen_cnt), source=0)
    do i = 1, cen_cnt
        if( i < cen_cnt )then
            eps = (centers(i+1) - centers(i))/2.
        else
            eps = (centers(i)   - centers(i-1))/2.
        endif
        do j = 1, cnt
            dist = abs(dists(j) - centers(i))
            if( dist < eps ) cen_cnts(i) = cen_cnts(i) + 1
        enddo
    enddo
    sum_cnt = sum(cen_cnts)
    if( allocated(ref_stats) ) deallocate(ref_stats)
    allocate(ref_stats(cen_cnt))
    do i = 1, cen_cnt
        ref_stats(i) = real(cen_cnts(i)) * 100. / real(sum_cnt)
    enddo
    ! reading pdbfiles and compute Kolmogorov-Smirnov test
    call read_filetable(p%pdbfiles, pdbfnames)
    npdbs = size(pdbfnames)
    if( allocated(cur_stats) ) deallocate(cur_stats)
    allocate(cur_stats(cen_cnt))
    if( allocated(probs) )deallocate(probs)
    allocate(probs(npdbs))
    do ipdb = 1, npdbs
        call read_pdb2matrix( trim(pdbfnames(ipdb)), cur_mat )
        Natoms = size(cur_mat,2)
        if( allocated(dists) )deallocate(dists)
        allocate(dists(Natoms**2), source=0.)
        cnt      = 0
        min_dist = huge(min_dist)
        do i = 1, Natoms
            do j = 1, Natoms
                cnt        = cnt + 1
                dists(cnt) = sqrt(sum((cur_mat(:,i) - cur_mat(:,j))**2))
                if( dists(cnt) < min_dist .and. dists(cnt) > TINY ) min_dist = dists(cnt)
            enddo
        enddo
        cen_cnts = 0
        do i = 1, cen_cnt
            if( i < cen_cnt )then
                eps = (centers(i+1) - centers(i))/2.
            else
                eps = (centers(i)   - centers(i-1))/2.
            endif
            do j = 1, cnt
                dist = abs(dists(j) - centers(i))
                if( dist < eps ) cen_cnts(i) = cen_cnts(i) + 1
            enddo
        enddo
        sum_cnt   = sum(cen_cnts)
        cur_stats = 0.
        do i = 1, cen_cnt
            cur_stats(i) = real(cen_cnts(i)) * 100. / real(sum_cnt)
        enddo
        call kstwo( ref_stats, cen_cnt, cur_stats, cen_cnt, d, prob )
        probs(ipdb) = prob
        ! print *, 'd = ', d, ', prob = ', prob
        ! print *, '----'
        ! print *, 'ref_stats = ', ref_stats
        ! print *, 'cur_stats = ', cur_stats
    enddo
    ! print *, '----'
    ! print *, 'probs = ', probs
    call hpsort(probs)
    call reverse(probs)
    prev_prob = 0.
    stab_cnt  = 0
    do ipdb = 1, npdbs
        if( ipdb > 1 .and. abs(prev_prob - probs(ipdb)) > TINY ) exit
        stab_cnt  = stab_cnt + 1
        prev_prob = probs(ipdb)
    enddo
    print *, 'At thres = ', dist_thres, '; ', stab_cnt, ' out of ', npdbs, ' are stable at prob = ', prev_prob
enddo

contains

    subroutine sort_cluster( points_in, cur_clusN, out_points, out_cnts, out_vars )
        real,                 intent(in)    :: points_in(:)
        integer,              intent(inout) :: cur_clusN
        real,    allocatable, intent(inout) :: out_points(:)
        integer, allocatable, intent(inout) :: out_cnts(:)
        real,    allocatable, intent(inout) :: out_vars(:)
        real,    allocatable :: points(:), avg_points(:), vars(:)
        integer, allocatable :: cnts(:)
        integer :: cur_ind, n_points, i, j
        real    :: cur_avg
        n_points = size(points_in)
        allocate(points(n_points), source=points_in)
        allocate(avg_points(n_points),cnts(n_points),vars(n_points))
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
                vars(cur_clusN)       = 0.
                do j = cur_ind, i
                    vars(cur_clusN) = vars(cur_clusN) + (points(i) - avg_points(cur_clusN))**2
                enddo
                vars(cur_clusN) = vars(cur_clusN)/real(i-cur_ind+1)
                cur_ind = i + 1
                cur_avg = 0.
            endif
        enddo
        ! last cluster
        cur_clusN = cur_clusN + 1
        avg_points(cur_clusN) = sum(points(cur_ind:i))/real(i-cur_ind+1)
        cnts(cur_clusN)       = i-cur_ind+1
        vars(cur_clusN)       = 0.
        do j = cur_ind, i
            vars(cur_clusN) = vars(cur_clusN) + (points(i) - avg_points(cur_clusN))**2
        enddo
        vars(cur_clusN) = vars(cur_clusN)/real(i-cur_ind+1)
        if( allocated(out_points) )deallocate(out_points)
        if( allocated(out_cnts)   )deallocate(out_cnts)
        if( allocated(out_vars)   )deallocate(out_vars)
        allocate(out_points(1:cur_clusN), source=avg_points(1:cur_clusN))
        allocate(out_cnts(1:cur_clusN),   source=cnts(1:cur_clusN))
        allocate(out_vars(1:cur_clusN),   source=vars(1:cur_clusN))
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
