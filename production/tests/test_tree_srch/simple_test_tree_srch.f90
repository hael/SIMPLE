program simple_test_tree_srch
! Test to see if can find arbitrary reprojection in tree structure
use simple_core_module_api
use simple_image,             only: image
use simple_cmdline,           only: cmdline
use simple_parameters,        only: parameters
use simple_oris
use simple_projector,         only: projector
use simple_commanders_sim,    only: commander_simulate_particles
use simple_commanders_atoms,  only: commander_pdb2mrc
use simple_tree 
implicit none 
type(commander_pdb2mrc)       :: xpdb2mrc 
type(image)                   :: vol
type(parameters)              :: p1
type(cmdline)                 :: cline, cline_pdb2mrc
type(oris)                    :: spiral
type(ori)                     :: o1  
character(len=:), allocatable :: cmd
type(image), allocatable      :: reproj(:)
type(dendro)        :: test_tree
real                :: objs(2), t1, t2
integer             :: indxs(2), i, ifoo, rc, NPLANES = 500, ORI_IND = 10
real, allocatable   :: test_dist_mat(:,:), smat_ref(:,:), dmat_ref(:,:)
type(s2_node), pointer :: p
logical             :: done = .false. 

! Load Volume 
! write(*, *) 'Downloading the example dataset...'
! cmd = 'curl -s -o 1JYX.pdb https://files.rcsb.org/download/1JYX.pdb'
! call execute_command_line(cmd, exitstat=rc)
! write(*, *) 'Converting .pdb to .mrc...'
! call cline_pdb2mrc%set('smpd',                            1.)
! call cline_pdb2mrc%set('pdbfile',                 '1JYX.pdb')
! call cline_pdb2mrc%checkvar('smpd',                        1)
! call cline_pdb2mrc%checkvar('pdbfile',                     2)
! call cline_pdb2mrc%check()
! call xpdb2mrc%execute_safe(cline_pdb2mrc)
! call cline_pdb2mrc%kill()
! cmd = 'rm 1JYX.pdb'
! call execute_command_line(cmd, exitstat=rc)

! Reproject Volume 
! call cline%set('smpd'   , 1.)
! call cline%set('nthr'   , 16.)
! call cline%set('vol1'   , '1JYX.mrc')
! call cline%set('mskdiam', 180.)
! call cline%set('lp'   ,   3.)

! Spiral Reprojections of Volume 
! call p1%new(cline)
! call find_ldim_nptcls(p1%vols(1), p1%ldim, ifoo)
! call vol%new(p1%ldim, p1%smpd)
! call vol%read(p1%vols(1))
! call spiral%new(NPLANES, is_ptcl=.false.)
! call spiral%spiral
! do i = 1, NPLANES
!     call spiral%get_ori(o1, i)
!     call vol%fproject(o1, reproj(i))
!     call reproj(i)%ifft()
! end do 

! allocate(smat_ref(nspace, nspace))

! do ispace = 1, NPLANES - 1
!     do jspace = ispace + 1, NPLANES 
!         smat_ref(ispace, jspace) = 
!         smat_ref(jspace, ispace) = smat_ref(ispace, jspace)
!     end do
!     smat_ref(ispace, ispace) = 1.  
! end do 

! Calc FM distance matrix no ctf 

! D = 1 - S 

allocate(test_dist_mat(10,10))

! test_dist_mat = reshape( [0., 0.1, 0.8, 0.9, 0.1, 0., 0.6, 0.3, 0.8, 0.6, 0., 0.7, 0.9, 0.3, 0.7, 0.], [4,4] )

test_dist_mat = reshape([ &
! Cluster 1 (1,2)
    0, 1, 10,10,20,20,30,30,40,40, &
    1, 0, 10,10,20,20,30,30,40,40, &

! Cluster 2 (3,4)
    10,10, 0, 1,10,10,20,20,30,30, &
    10,10, 1, 0,10,10,20,20,30,30, &

! Cluster 3 (5,6)
    20,20,10,10, 0, 1,10,10,20,20, &
    20,20,10,10, 1, 0,10,10,20,20, &

! Cluster 4 (7,8)
    30,30,20,20,10,10, 0, 1,10,10, &
    30,30,20,20,10,10, 1, 0,10,10, &

! Cluster 5 (9,10)
    40,40,30,30,20,20,10,10, 0, 1, &
    40,40,30,30,20,20,10,10, 1, 0  &
], [10,10])

! print *, test_dist_mat
! test_dist_mat = reshape( [0., 0.8, 0., 0.8], [2,2])
call test_tree%set_distmat(test_dist_mat)
call test_tree%set_npnts(10)
call cpu_time(t1)
call test_tree%build_dendro
call cpu_time(t2)
print *, 'tree build time', t2 - t1

! point to root to preserve tree structure
p => test_tree%root 
! initialize 
objs(1) = real(p%left%ref_idx)**2
objs(2) = real(p%right%ref_idx)**2
do 
    print *, objs
    call walk_dendro(p, indxs, objs, done)
    if(done) exit
    print *, indxs
    objs(1) = real(indxs(1))**2
    objs(2) = real(indxs(2))**2
end do 

call print_dendro(test_tree%root)

end program simple_test_tree_srch
    
