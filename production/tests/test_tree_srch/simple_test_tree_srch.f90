program simple_test_tree_srch
! Test to see if can find arbitrary reprojection in tree structure
use simple_core_module_api
use simple_image,             only: image
use simple_cmdline,           only: cmdline
use simple_parameters,        only: parameters
use simple_oris
use simple_projector,         only: projector
use simple_commanders_atoms,  only: commander_pdb2mrc
use simple_tree 
use simple_simple_volinterp
use simple_corrmat
use simple_timer
use simple_aff_prop
use simple_stat
implicit none 
type(parameters)         :: p1
type(oris)               :: spiral
type(commander_pdb2mrc)  :: xpdb2mrc
type(cmdline)            :: cline_pdb2mrc, cline
type(multi_dendro)       :: test_tree
real                     :: objs(2), t1, t2, simsum
integer                  :: indxs(2), rc, ifoo, NPROJ = 10, nthr_max = 10, i, j
type(image)              :: vol
logical                  :: done = .false. 
real(timer_int_kind)     :: rt_cc, rt_ap, rt_tr 
integer(timer_int_kind)  ::  t_cc, t_ap, t_tr 
type(aff_prop)           :: affprop    
character(len=:), allocatable   :: cmd  
real, allocatable               :: dist_mat_cc(:,:), sub_distmats(:,:,:)
type(image), allocatable        :: proj_arr(:)
integer, allocatable            :: centers(:), labels(:), clus_pops(:)
type(multi_dendro), allocatable :: roots(:)
type(s2_node), pointer   :: p

! Load Volume 
write(*, *) 'Downloading the example dataset...'
cmd = 'curl -s -O --cacert /etc/ssl/certs/ca-bundle.crt https://files.rcsb.org/download/1JYX.pdb'

call execute_command_line(cmd, exitstat=rc)
write(*, *) 'Converting .pdb to .mrc...'
call cline_pdb2mrc%set('smpd',                            3.)
call cline_pdb2mrc%set('pdbfile',                 '1JYX.pdb')
call cline_pdb2mrc%checkvar('smpd',                        1)
call cline_pdb2mrc%checkvar('pdbfile',                     2)
call cline_pdb2mrc%check()
call xpdb2mrc%execute_safe(cline_pdb2mrc)
call cline_pdb2mrc%kill()
cmd = 'rm 1JYX.pdb'
call execute_command_line(cmd, exitstat=rc)

! Reproject Volume 
call cline%set('smpd'   , 3.)
call cline%set('vol1'   , '1JYX.mrc')
call cline%set('mskdiam', 180.)
call cline%set('hp', 20.)
call cline%set('lp'   ,   6.)
call cline%set('trs', 5.)
call cline%set('ctf', 'no')
call cline%set('nthr', 8)
call cline%set('objfun', 'cc')
call cline%set('mkdir', 'no')


! Spiral Reprojections of Volume 
call p1%new(cline)
call find_ldim_nptcls(p1%vols(1), p1%ldim, ifoo)
call vol%new(p1%ldim, p1%smpd)
call vol%read(p1%vols(1))
call spiral%new(NPROJ, is_ptcl=.false.)
call spiral%spiral

allocate(proj_arr(NPROJ))
proj_arr = reproject(vol, spiral, NPROJ)
params_glob%ldim = proj_arr(1)%get_ldim()
allocate(dist_mat_cc(NPROJ,NPROJ))
t_cc = tic()
dist_mat_cc = calc_inpl_invariant_cc_nomirr(p1%hp, p1%lp, p1%trs, proj_arr)
rt_cc = toc(t_cc)


! deallocate(dist_mat_cc)

! do i = 1, NPROJ 
!     do j = 1, nthr_max  
!         allocate(dist_mat_cc(i, i))
!         params_glob%nthr = j
!         t_cc = tic()
!         dist_mat_cc = calc_inpl_invariant_cc_nomirr(p1%hp, p1%lp, p1%trs, proj_arr(1:i))
!         rt_cc = toc(t_cc)
!         print *, 'nimgs', i, 'nthr', j, 'time_cc', rt_cc
!         deallocate(dist_mat_cc)
!     end do 
! end do 

! print *, 'N_PROJS:', NPROJ, 'cc time:', rt_cc

call affprop%new(NPROJ, dist_mat_cc, pref=0.)
t_ap = tic()
call affprop%propagate(centers, labels, simsum)
rt_ap = toc(t_ap)

print *, 'N_PROJS:', NPROJ, 'ap time:', rt_ap, 'cc time:', rt_cc

! print *, simsum, size(centers), labels
call test_tree%set_distmat(dist_mat_cc)
call test_tree%set_medoids(centers)
call test_tree%set_clus_pops(labels)
call test_tree%set_subsets(labels)
print *, 'dim', shape(test_tree%subsets), 'n_trees', test_tree%n_trees, 'clus_pops', test_tree%clus_pops
print *, 'matrix', test_tree%subsets
call test_tree%build_multi_dendro()
print *, test_tree%root_array(1)%subset
print *, test_tree%root_array(2)%subset
print *, test_tree%root_array(3)%subset

! should be 1, 4, 9
! split into sub_dmats to make the trees from 
! exemplars will be tree roots


! print *, test_dist_mat
! test_dist_mat = reshape( [0., 0.8, 0., 0.8], [2,2])
! call test_tree%set_distmat(test_dist_mat)
! call test_tree%set_npnts(10)
! call cpu_time(t1)
! call test_tree%build_dendro
! call cpu_time(t2)
! print *, 'tree build time', t2 - t1

! ! point to root to preserve tree structure
! p => test_tree%root 
! ! initialize 
! objs(1) = real(p%left%ref_idx)**2
! objs(2) = real(p%right%ref_idx)**2
! do 
!     print *, objs
!     call walk_dendro(p, indxs, objs, done)
!     if(done) exit
!     print *, indxs
!     objs(1) = real(indxs(1))**2
!     objs(2) = real(indxs(2))**2
! end do 

! call print_dendro(test_tree%root)

end program simple_test_tree_srch
    
