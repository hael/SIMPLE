program simple_test_dists_cluster
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_parameters,         only: parameters
use simple_nanoparticle_utils, only: read_pdb2matrix, write_matrix2pdb
use simple_atoms
implicit none
character(len=LONGSTRLEN), allocatable :: pdbfnames(:)
real,    parameter   :: PROB_THRES = 0.7
real,    parameter   :: MIN_THRES = 10.    ! automatically figured from the first part of the codes
real,    allocatable :: dists(:), mat(:,:), dists_cen(:), vars(:), ref_stats(:), cur_mat(:,:), cur_stats(:), probs(:), out_mat(:,:)
integer, allocatable :: cnts(:)
logical, allocatable :: atom_msk(:), max_msk(:), thres_msk(:), taken(:,:)
type(atoms)          :: molecule
type(parameters)     :: p
type(cmdline)        :: cline
integer              :: Natoms, i, j, k, cnt, nclus, sum_cnt, cen_cnt, npdbs, ipdb, ithres, stab_cnt
logical              :: l_found
real                 :: min_dist, eps, dist, d, dist_thres, rad_thres, stab_prob, max_prob, mid(3), cur_mid(3)
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
    do i = 1, Natoms-1
        do j = i+1, Natoms
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
    if( allocated(ref_stats) )deallocate(ref_stats)
    allocate(ref_stats(cen_cnt), source=0.)
    do i = 1, cen_cnt
        if( i < cen_cnt )then
            eps = (dists_cen(i+1) - dists_cen(i))/2.
        else
            eps = (dists_cen(i)   - dists_cen(i-1))/2.
        endif
        do j = 1, cnt
            dist = abs(dists(j) - dists_cen(i))
            if( dist < eps ) ref_stats(i) = ref_stats(i) + 1.
        enddo
    enddo
    ref_stats = ref_stats * 100. / sum(ref_stats)
    ! reading pdbfiles and compute Kolmogorov-Smirnov test
    call read_filetable(p%pdbfiles, pdbfnames)
    npdbs = size(pdbfnames)
    if( allocated(cur_stats) ) deallocate(cur_stats)
    allocate(cur_stats(cen_cnt))
    if( allocated(probs) )deallocate(probs)
    allocate(probs(npdbs))
    do ipdb = 1, npdbs
        call read_pdb2matrix( trim(pdbfnames(ipdb)), cur_mat )
        call compute_stats(cur_mat, dists_cen(1:cen_cnt), cur_stats)
        call kstwo( ref_stats, cen_cnt, cur_stats, cen_cnt, d, probs(ipdb) )
        ! print *, 'd = ', d, ', prob = ', probs(ipdb)
        ! print *, '----'
        ! print *, 'ref_stats = ', ref_stats
        ! print *, 'cur_stats = ', cur_stats
    enddo
    ! print *, '----'
    ! print *, 'probs = ', probs
    call hpsort(probs)
    call reverse(probs)
    stab_cnt  = 0
    do ipdb = 1, npdbs
        ! if( ipdb > 1 .and. abs(prev_prob - probs(ipdb)) > TINY ) exit
        if( ipdb > 1 .and. probs(ipdb) < PROB_THRES ) exit
        stab_cnt  = stab_cnt + 1
    enddo
    print *, 'At thres = ', dist_thres, '; ', stab_cnt, ' out of ', npdbs, ' are stable at prob > ', PROB_THRES
enddo
rad_thres = 2.
l_found   = .false.
! finding the radius corresponding to the current thres
call read_pdb2matrix( p%pdbfile, mat )
Natoms = size(mat,2)
allocate(taken(Natoms,Natoms), source=.false.)
allocate(atom_msk(Natoms),  source=.false.)
do ithres=1,8
    if( l_found ) exit
    rad_thres = rad_thres + 1.
    do i = 1, Natoms
        if( l_found ) exit
        ! find the neighborhood
        atom_msk = .false.
        do j = 1, Natoms
            if( sqrt(sum((mat(:,i) - mat(:,j))**2)) < rad_thres )atom_msk(j) = .true.
        enddo
        taken = .false.
        do j = 1, Natoms
            do k = 1, Natoms
                if( taken(j,k) .or. .not.atom_msk(k) .or. .not.atom_msk(j) )cycle
                taken(j,k) = .true.
                taken(k,j) = .true.
                if( sqrt(sum((mat(:,j) - mat(:,k))**2)) > MIN_THRES )then
                    l_found = .true.
                    exit
                endif
            enddo
        enddo
    enddo
enddo
print *, 'Radius thres ', rad_thres, ' is corresponding to distribution thres ', MIN_THRES
! testing finding the core by maximizing the prob computed above for each pdb file
do ipdb = 1, npdbs
    call molecule%new(trim(pdbfnames(ipdb)))
    call read_pdb2matrix( trim(pdbfnames(ipdb)), cur_mat )
    Natoms = size(cur_mat,2)
    if( allocated(atom_msk) )deallocate(atom_msk)
    if( allocated(max_msk)  )deallocate(max_msk)
    if( allocated(thres_msk))deallocate(thres_msk)
    allocate(atom_msk(Natoms),  source=.false.)
    allocate(max_msk(Natoms),   source=.false.)
    allocate(thres_msk(Natoms), source=.false.)
    do i = 1, Natoms
        ! find the neighborhood
        atom_msk = .false.
        do j = 1, Natoms
            if( sqrt(sum((cur_mat(:,i) - cur_mat(:,j))**2)) < rad_thres )atom_msk(j) = .true.
        enddo
        ! computing the prob and maximize it
        call compute_stats(cur_mat, dists_cen(1:cen_cnt), cur_stats, atom_msk)
        call kstwo( ref_stats, cen_cnt, cur_stats, cen_cnt, d, probs(ipdb) )
        if( probs(ipdb) > PROB_THRES )thres_msk(i) = .true.
    enddo
    ! center of mass
    mid = 0.
    do i = 1, Natoms
        mid = mid + cur_mat(:,i)
    enddo
    mid      = mid / real(Natoms)
    max_prob = huge(max_prob)
    do i = 1, Natoms
        if( .not. thres_msk(i) )cycle
        ! find the neighborhood
        atom_msk = .false.
        cur_mid  = 0.
        do j = 1, Natoms
            if( sqrt(sum((cur_mat(:,i) - cur_mat(:,j))**2)) < rad_thres )then
                atom_msk(j) = .true.
                cur_mid     = cur_mid + cur_mat(:,j)
            endif
        enddo
        cur_mid = cur_mid / real(count(atom_msk .eqv. .true.))
        if( sqrt(sum((cur_mid - mid)**2)) < max_prob )then
            max_prob = sqrt(sum((cur_mid - mid)**2))
            max_msk  = atom_msk
        endif
    enddo
    print *, 'ipdb = ', ipdb, ' count = ', count(max_msk .eqv. .true.), ' out of ', Natoms
    if( allocated(out_mat) )deallocate(out_mat)
    allocate(out_mat(3,count(max_msk .eqv. .true.)))
    cnt = 0
    do i = 1, Natoms
        if( .not. max_msk(i) )cycle
        cnt            = cnt + 1
        out_mat(:,cnt) = cur_mat(:,i)
        call molecule%set_beta(i, 1.0)
    enddo
    call write_matrix2pdb( 'Pt', out_mat, get_fpath(trim(pdbfnames(ipdb)))//'test_core_ATMS.pdb' )
enddo

contains

    ! compute the distribution of the dists around the reference distance points
    subroutine compute_stats(pos_mat, ref_dists, out_stats, msk)
        real,    allocatable, intent(in)    :: pos_mat(:,:)
        real,                 intent(in)    :: ref_dists(:)
        real,    allocatable, intent(inout) :: out_stats(:)
        logical, optional,    intent(in)    :: msk(:)
        integer :: Na, Nr, i, j, k
        real    :: tmp_eps
        Na        = size(pos_mat,2)
        Nr        = size(ref_dists)
        out_stats = 0.
        do i = 1, Nr
            if( i < Nr )then
                tmp_eps = (ref_dists(i+1) - ref_dists(i))/2.
            else
                tmp_eps = (ref_dists(i)   - ref_dists(i-1))/2.
            endif
            do j = 1, Na-1
                do k = j+1, Na
                    if( present(msk) )then
                        if( .not.msk(j) .or. .not.msk(k) )cycle
                    endif
                    dist = abs(sqrt(sum((pos_mat(:,j) - pos_mat(:,k))**2)) - ref_dists(i))
                    if( dist < tmp_eps ) out_stats(i) = out_stats(i) + 1.
                enddo
            enddo
        enddo
        out_stats = out_stats * 100. / sum(out_stats)
    end subroutine compute_stats

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
    
