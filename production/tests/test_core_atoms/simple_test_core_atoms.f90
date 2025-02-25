program simple_test_core_atoms
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_parameters,         only: parameters
use simple_nanoparticle_utils, only: atoms_register, Kabsch_algo, read_pdb2matrix, write_matrix2pdb
implicit none
real,    allocatable :: core_mat(:,:), mat(:,:), dists(:), dists_avg(:)
integer, allocatable :: inds(:)
integer, parameter   :: NCORE = 50
type(parameters)     :: p
type(cmdline)        :: cline
integer              :: Natoms, i, j, cnt
real                 :: centroid(3), min_dist, dist
! reading pdb file
if( command_argument_count() < 2 )then
    write(logfhandle,'(a)') 'Usage: simple_test_core_atoms nthr=yy pdbfile=zz'
    stop
endif
call cline%parse_oldschool
call cline%checkvar('nthr',     1)
call cline%checkvar('pdbfile' , 2)
call cline%check
call p%new(cline)
call read_pdb2matrix( p%pdbfile, mat )
Natoms = size(mat,2)
allocate(dists(Natoms), dists_avg(Natoms), source=0.)
do i = 1, Natoms
    do j = 1, Natoms
        dists(j) = sqrt(sum((mat(:,i) - mat(:,j))**2))
    enddo
    call hpsort(dists)
    dists_avg(i) = sum(dists(1:NCORE)) / real(NCORE)
enddo
inds = (/(j,j=1,Natoms)/)
call hpsort(dists_avg,inds)
allocate(core_mat(3,NCORE),source=mat(:,inds(1:NCORE)))
! finding center of mass of the dists_avg
centroid = 0.
do i = 1, NCORE
    centroid = centroid + core_mat(:,i) * dists_avg(i)
end do
centroid = centroid / sum(dists_avg(1:NCORE))
! compute the index closest to centroid
min_dist = huge(min_dist)
do j = 1, Natoms
    dist = sqrt(sum((mat(:,j) - centroid)**2))
    if( dist < min_dist )then
        min_dist = dist
        i        = j
    endif
enddo
! get the neighbor of i as the core
do j = 1, Natoms
    dists(j) = sqrt(sum((mat(:,i) - mat(:,j))**2))
enddo
inds = (/(j,j=1,Natoms)/)
call hpsort(dists,inds)
core_mat = mat(:,inds(1:NCORE))
if( cline%defined('pdbout') )then
    call write_matrix2pdb( 'Pt', core_mat, p%pdbout )
else
    call write_matrix2pdb( 'Pt', core_mat, 'ATMS_core_atoms.pdb' )
endif
end program simple_test_core_atoms
