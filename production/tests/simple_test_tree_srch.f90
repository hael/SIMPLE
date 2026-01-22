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
real(timer_int_kind)     :: rt_cc, rt_ap, rt_tr_build, rt_tr_search, rt_tr_tot 
integer(timer_int_kind)  ::  t_cc, t_ap, t_tr_build, t_tr_search
type(aff_prop)           :: affprop    
character(len=:), allocatable   :: cmd  
real, allocatable               :: dist_mat_cc(:,:), sub_distmats(:,:,:)
type(image), allocatable        :: proj_arr(:)
integer, allocatable            :: centers(:), labels(:), cls_pops(:)
type(multi_dendro), allocatable :: roots(:)
type(s2_node), pointer   :: rootp, found

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

print *, dist_mat_cc
call normalize_minmax(dist_mat_cc)
call affprop%new(NPROJ, dist_mat_cc, pref=-2.)
t_ap = tic()
call affprop%propagate(centers, labels, simsum)
print *, 'labels', labels
rt_ap = toc(t_ap)
! centers = [5]
! labels = [1,1,1,1,1,1,1,1,1,1]
print *, 'N_PROJS:', NPROJ, 'ap time:', rt_ap, 'cc time:', rt_cc
! print *, simsum, size(centers), labels
dist_mat_cc = 1. - dist_mat_cc

call test_tree%set_distmat(dist_mat_cc)
call test_tree%set_cls_pops(labels)
call test_tree%set_subsets(labels)
t_tr_build = tic()
call test_tree%build_multi_dendro(2)
rt_tr_build = toc(t_tr_build)
print *, 'tree build time', rt_tr_build

! 2n - 1 th node is root 
rootp => test_tree%node_store(1)%nodes(test_tree%node_store(1)%root_idx)
call print_s2_tree(rootp, show_subset = .true.)
t_tr_search = tic()
found => search_tree4ref(rootp, 5)
rt_tr_search = toc(t_tr_search)
print *, 'tree search time', rt_tr_search

if (.not. associated(found)) then
    print *, "found is NULL"
  else
    print *, "found%ref_idx = ", found%ref_idx
end if
! print *, 'tree build time', rt_tr
! print *, 'tree_build_time', rt_tr
! print *, 'pool', test_tree%root_array(1)%subset
! print *, 'right', test_tree%root_array(1)%right%subset
! print *, 'left', test_tree%root_array(1)%left%subset
! print *, 'right => left', test_tree%root_array(1)%right%left%subset
! print *, 'right => right', test_tree%root_array(1)%right%right%subset

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
!     call walk_from_node(p, indxs, objs, done)
!     if(done) exit
!     print *, indxs
!     objs(1) = real(indxs(1))**2
!     objs(2) = real(indxs(2))**2
! end do 

! call print_dendro(test_tree%root)

end program simple_test_tree_srch
    
