program simple_test_pdb_register
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_parameters,         only: parameters
use simple_nanoparticle_utils, only: atoms_register, Kabsch_algo, read_pdb2matrix, write_matrix2pdb
implicit none
character(len=LONGSTRLEN), allocatable :: pdbfnames(:)
type(parameters)     :: p
type(cmdline)        :: cline
integer, parameter   :: NCORE = 100
integer, allocatable :: inds(:)
real,    allocatable :: cur_mat(:,:), mat(:,:), matrix_rot(:,:)
integer :: npdbs, ipdb, j
real    :: core_mat(3,NCORE), core_rot(3,NCORE), cur_core(3,NCORE), rot_mat(3,3), trans_vec(3), scale 
if( command_argument_count() < 3 )then
    write(logfhandle,'(a)') 'Usage: simple_test_atom_posmatch nthr=yy pdbfile=zz pdbfiles=tt'
    stop
endif
call cline%parse_oldschool
call cline%checkvar('nthr',     1)
call cline%checkvar('pdbfile' , 2)
call cline%checkvar('pdbfiles', 3)
call cline%check
call p%new(cline)
if( .not. cline%defined('maxits') ) p%maxits = 1
! reading the reference file and finding its core
call read_pdb2matrix(p%pdbfile, mat)
call find_core(mat, core_mat)
! reading pdb files
call read_filetable(p%pdbfiles, pdbfnames)
npdbs = size(pdbfnames)
do ipdb = 1, npdbs
    call read_pdb2matrix( trim(pdbfnames(ipdb)), cur_mat )
    call find_core(cur_mat, cur_core)
    call atoms_register(cur_core, core_mat, core_rot, maxits=p%maxits, out_mat=rot_mat, out_trans=trans_vec, out_scale=scale)
    if(allocated(matrix_rot)) deallocate(matrix_rot)
    allocate(matrix_rot(3,size(cur_mat,2)))
    do j = 1, size(cur_mat,2)
        matrix_rot(:,j) = scale * matmul(rot_mat, cur_mat(:,j)) + trans_vec
    enddo
    call write_matrix2pdb( 'Pt', cur_core, 'ATMS_cur_core_'//int2str(ipdb)//'.pdb' )
    call write_matrix2pdb( 'Pt', core_mat, 'ATMS_core_mat_'//int2str(ipdb)//'.pdb' )
    if( cline%defined('pdbout') )then
        call write_matrix2pdb( 'Pt', matrix_rot, p%pdbout )
    else
        call write_matrix2pdb( 'Pt', matrix_rot, 'ATMS_rec_'//int2str(ipdb)//'.pdb' )
    endif
enddo

contains

    subroutine find_core(mat_in, core_out)
        real, intent(in)    :: mat_in(:,:)
        real, intent(inout) :: core_out(:,:)
        real, allocatable   :: dists(:)
        real    :: centroid(3)
        integer :: i, N
        N        = size(mat_in,2)
        centroid = 0.
        do i = 1, N
            centroid = centroid + mat_in(:,i)
        enddo
        centroid = centroid / real(N)
        allocate(dists(N))
        do i = 1, N
            dists(i) = sqrt(sum((mat_in(:,i) - centroid)**2))
        enddo
        inds = (/(i,i=1,N)/)
        call hpsort(dists, inds)
        core_out = mat_in(:,inds(1:NCORE))
    end subroutine find_core

end program simple_test_pdb_register
