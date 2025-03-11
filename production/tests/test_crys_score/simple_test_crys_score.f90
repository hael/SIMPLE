program simple_test_crys_score
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_parameters,         only: parameters
use simple_nanoparticle_utils, only: read_pdb2matrix, write_matrix2pdb
use simple_atoms
implicit none
character(len=LONGSTRLEN), allocatable :: pdbfnames(:), cmd, corenames(:)
real,                      allocatable :: dists(:), mat(:,:), dists_cen(:), vars(:), ref_stats(:), cur_mat(:,:), cur_stats(:)
integer,                   allocatable :: cnts(:)
logical,                   allocatable :: l_dists(:)
character(len=LONGSTRLEN) :: t_name
type(parameters)          :: p
type(cmdline)             :: cline
integer                   :: npdbs, ipdb, icore, ncores, Natoms, cnt, i, j, nclus
real                      :: min_dist, dist, prob, d
! reading pdb file
if( command_argument_count() < 2 )then
    write(logfhandle,'(a)') 'Usage: simple_test_crys_score nthr=yy pdbfiles=tt'
    stop
endif
call cline%parse_oldschool
call cline%checkvar('nthr',     1)
call cline%checkvar('pdbfiles', 2)
call cline%check
call p%new(cline)
call read_filetable(p%pdbfiles, pdbfnames)
npdbs = size(pdbfnames)
do ipdb = 1, npdbs
    ! generating stats for the reference lattice
    call read_pdb2matrix( trim(pdbfnames(ipdb)) // 'startvol_ATMS.pdb', mat )
    Natoms = size(mat,2)
    if( allocated(dists)   )deallocate(dists)
    if( allocated(l_dists) )deallocate(l_dists)
    allocate(  dists(Natoms**2), source=0.)
    allocate(l_dists(Natoms**2), source=.false.)
    cnt      = 0
    min_dist = huge(min_dist)
    do i = 1, Natoms-1
        do j = i+1, Natoms
            cnt          = cnt + 1
            dists(cnt)   = sqrt(sum((mat(:,i) - mat(:,j))**2))
            l_dists(cnt) = .true.
            if( dists(cnt) < min_dist .and. dists(cnt) > TINY ) min_dist = dists(cnt)
        enddo
    enddo
    dists = pack(dists, mask=l_dists)
    call sort_cluster(dists, nclus, dists_cen, cnts, vars, delta=min_dist / 10.)
    call hpsort(dists_cen)
    if( allocated(ref_stats) )deallocate(ref_stats)
    if( allocated(cur_stats) )deallocate(cur_stats)
    allocate(ref_stats(nclus), cur_stats(nclus), source=0.)
    call compute_stats(mat, dists_cen(1:nclus), ref_stats)
    ! generating stats for the core
    t_name = trim(pdbfnames(ipdb)) // 'segvols/1_tseries_core_atoms_analysis/'
    print *, trim(t_name)
    cmd = 'ls '// trim(t_name) // '*_core.pdb | xargs -n 1 basename > output_'//int2str(ipdb)//'.log'
    call execute_command_line(cmd)
    call read_filetable('output_'//int2str(ipdb)//'.log', corenames)
    ncores = size(corenames)
    do icore = 1, ncores
        call read_pdb2matrix( trim(t_name) // trim(corenames(icore)), cur_mat )
        call compute_stats(cur_mat, dists_cen(1:nclus), cur_stats)
        call kstwo( ref_stats, nclus, cur_stats, nclus, d, prob )
        print *, trim(corenames(icore)), prob
    enddo
    call execute_command_line('rm output_'//int2str(ipdb)//'.log')
enddo

contains

    ! compute the distribution of the dists around the reference distance points
    subroutine compute_stats(pos_mat, ref_dists, out_stats, msk)
        real,    allocatable, intent(in)    :: pos_mat(:,:)
        real,                 intent(in)    :: ref_dists(:)
        real,    allocatable, intent(inout) :: out_stats(:)
        logical, optional,    intent(in)    :: msk(:)
        integer :: Na, Nr, i, j, k
        real    :: eps_g, eps_l
        Na        = size(pos_mat,2)
        Nr        = size(ref_dists)
        out_stats = 0.
        do i = 1, Nr
            if( i > 1 .and. i < Nr )then
                eps_g = (ref_dists(i+1) - ref_dists(i))/2.
                eps_l = (ref_dists(i)   - ref_dists(i-1))/2.
            elseif( i == 1 )then
                eps_g = (ref_dists(i+1) - ref_dists(i))/2.
                eps_l =  huge(eps_l)
            elseif( i == Nr )then
                eps_g =  huge(eps_g)
                eps_l = (ref_dists(i)   - ref_dists(i-1))/2.
            endif
            do j = 1, Na-1
                do k = j+1, Na
                    if( present(msk) )then
                        if( .not.msk(j) .or. .not.msk(k) )cycle
                    endif
                    dist = sqrt(sum((pos_mat(:,j) - pos_mat(:,k))**2))
                    if( (dist <= ref_dists(i) .and. abs(dist - ref_dists(i)) < eps_l ) .or. &
                        (dist >= ref_dists(i) .and. abs(dist - ref_dists(i)) < eps_g ) )then
                        out_stats(i) = out_stats(i) + 1.
                    endif
                enddo
            enddo
        enddo
        out_stats = out_stats * 100. / sum(out_stats)
    end subroutine compute_stats

    subroutine sort_cluster( points_in, cur_clusN, out_points, out_cnts, out_vars, delta )
        real,                 intent(in)    :: points_in(:)
        integer,              intent(inout) :: cur_clusN
        real,    allocatable, intent(inout) :: out_points(:)
        integer, allocatable, intent(inout) :: out_cnts(:)
        real,    allocatable, intent(inout) :: out_vars(:)
        real,                 intent(in)    :: delta
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
            if( abs(points(i+1) - cur_avg/real(i-cur_ind+1)) > delta )then
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

end program simple_test_crys_score