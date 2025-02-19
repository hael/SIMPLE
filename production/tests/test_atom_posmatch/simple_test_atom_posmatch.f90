program simple_test_atom_posmatch
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_parameters,         only: parameters
use simple_nanoparticle_utils, only: atoms_register, Kabsch_algo, read_pdb2matrix, write_matrix2pdb
implicit none
type(parameters)     :: p
type(cmdline)        :: cline
integer, parameter   :: N1 = 20, N2 = 5, MAX_POS = 10
integer, allocatable :: inds(:), select_inds(:)
real,    allocatable :: matrix1(:,:), matrix2(:,:), matrix_rot(:,:)
real    :: atom1_pos(3,N1), atom2_pos(3,N2), rot_mat(3,3), scale, translation(3), atom2_rnd(3,N2)
real    :: rec_mat(3,3), rec_trans(3), rec_scale, rec_atom2(3,N2), atom1_rot(3,N1), atom2_rot(3,N2)
integer :: i, j, mapping(N2)
call seed_rnd
do i = 1, N1
    do j = 1, 3
        atom1_pos(j,i) = floor(ran3() * real(MAX_POS))
    enddo
    print *, 'atom1 pos ', i, ' = ', atom1_pos(:,i)
enddo
rot_mat(1,:) = [-0.4297,  0.3140, -0.8466]
rot_mat(2,:) = [-0.7475, -0.6496,  0.1385]
rot_mat(3,:) = [-0.5065,  0.6924,  0.5139]
translation  = [0.2, -0.5, 1.]
scale        = 0.5
select_inds  = rnd_inds(N1)
mapping      = select_inds(1:N2)
call hpsort(mapping)
print *, '------------'
print *, 'mapping = ', mapping
atom2_pos    = scale * matmul(rot_mat, atom1_pos(:,mapping))
do i = 1, N2
    atom2_pos(:,i) = atom2_pos(:,i) + translation
enddo
print *, '------ TRUTH POSITIONS ------'
call Kabsch_algo(atom1_pos(:,mapping(1:3)), atom2_pos(:,1:3), rec_mat, rec_trans, rec_scale)
rec_atom2 = rec_scale * matmul(rec_mat, atom1_pos(:,mapping))
print *, 'atom 1 pos : ',      atom1_pos(:,mapping(4))
print *, 'atom 2 pos : ',      atom2_pos(:,4)
print *, 'rec scale = ',       rec_scale
print *, 'rec translation = ', rec_trans
print *, 'rotated atom 2 : ',  rec_atom2(:,4) + rec_trans
print *, 'ret_mat 1 = ', rec_mat(1,:)
print *, 'ret_mat 2 = ', rec_mat(2,:)
print *, 'ret_mat 3 = ', rec_mat(3,:)
! randomize atom2_pos
inds = rnd_inds(N2)
do i = 1, N2
    atom2_rnd(:,i)   = atom2_pos(:,inds(i))
    mapping(inds(i)) = i
enddo
print *, '------------'
print *, 'randomized inds = ', inds
! generating all permutations of atom2_rnd positions
print *, '------ ROTATED ATOMS 1 POSITIONS (SOME WILL MATCH ATOMS 2 POSITIONS BELOW) ------'
call atoms_register(atom1_pos, atom2_rnd, atom1_rot)
do i = 1, N1
    print *, 'rotated atom 1 pos : ', atom1_rot(:,i)
enddo
print *, '------ ATOMS 2 POSITIONS ------'
do i = 1, N2
    print *, 'atom 2 pos : ', atom2_rnd(:,mapping(i))
enddo
print *, '------ ROTATED ATOMS 2 POSITIONS (SOME WILL MATCH ATOMS 1 POSITIONS BELOW) ------'
call atoms_register(atom2_rnd, atom1_pos, atom2_rot)
do i = 1, N2
    print *, 'rotated atom 2 pos : ', atom2_rot(:,mapping(i))
enddo
print *, '------ ATOMS 1 POSITIONS ------'
do i = 1, N1
    print *, 'atom 1 pos : ', atom1_pos(:,i)
enddo
! reading pdb file
if( command_argument_count() < 3 )then
    write(logfhandle,'(a)') 'Usage: simple_test_opt_lp nthr=yy pdbfile=zz pdbfile2=tt'
else
    call cline%parse_oldschool
    call cline%checkvar('nthr',     1)
    call cline%checkvar('pdbfile' , 2)
    call cline%checkvar('pdbfile2', 3)
    call cline%check
    call p%new(cline)
    call read_pdb2matrix( p%pdbfile,  matrix1 )
    call read_pdb2matrix( p%pdbfile2, matrix2 )
    allocate(matrix_rot(3,size(matrix1,2)))
    call atoms_register(matrix1, matrix2, matrix_rot)
    if( cline%defined('pdbout') )then
        call write_matrix2pdb( 'Pt', matrix_rot, p%pdbout )
    else
        call write_matrix2pdb( 'Pt', matrix_rot, 'ATMS_rec.pdb' )
    endif
endif
end program simple_test_atom_posmatch
