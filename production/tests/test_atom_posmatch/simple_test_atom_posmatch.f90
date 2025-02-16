program simple_test_atom_posmatch
include 'simple_lib.f08'
implicit none
integer, parameter :: N1 = 5, N2 = 5, MAX_POS = 10
real    :: atom1_pos(3,N1), atom2_pos(3,N2), rot_mat(3,3), scale, translation(3), rec_mat(3,3), rec_trans(3), rec_scale, rec_atom2(3,N1)
integer :: i, j
call seed_rnd
do i = 1, N1
    do j = 1, 3
        atom1_pos(j,i) = floor(ran3() * real(MAX_POS))
    enddo
enddo
rot_mat(1,:) = [-0.4297,  0.3140, -0.8466]
rot_mat(2,:) = [-0.7475, -0.6496,  0.1385]
rot_mat(3,:) = [-0.5065,  0.6924,  0.5139]
translation  = [0.2, -0.5, 1.]
scale        = 0.5
atom2_pos    = scale * matmul(rot_mat, atom1_pos)
do i = 1, N2
    atom2_pos(:,i) = atom2_pos(:,i) + translation
enddo
call atom_posmatch(atom1_pos, atom2_pos, rec_mat, rec_trans, rec_scale)
rec_atom2 = rec_scale * matmul(rec_mat, atom1_pos)
print *, 'atom 1 pos : ',     atom1_pos(:,1)
print *, 'atom 2 pos : ',     atom2_pos(:,1)
print *, 'rec scale = ', rec_scale
print *, 'rec translation = ', rec_trans
print *, 'rotated atom 2 : ', rec_atom2(:,1) + rec_trans

contains

    subroutine atom_posmatch( pos1, pos2, ret_mat, ret_trans, ret_scale)
        real, intent(in)    :: pos1(:,:), pos2(:,:)
        real, intent(inout) :: ret_mat(3,3)
        real, intent(inout) :: ret_trans(3)
        real, intent(inout) :: ret_scale
        real, allocatable   :: var1(:,:), var2(:,:)
        real    :: mean1(1,3), mean2(1,3), mat(3,3), eig_vals(3), eig_vecs(3,3)
        integer :: i, N
        N          = size(pos1, 2)
        mean1(1,:) = sum(pos1, dim=2) / real(N)
        mean2(1,:) = sum(pos2, dim=2) / real(N)
        allocate(var1(3,N), var2(3,N))
        do i = 1, N
            var1(:,i) = pos1(:,i) - mean1(1,:)
            var2(:,i) = pos2(:,i) - mean2(1,:)
        enddo
        mat = matmul(var1, transpose(var2))
        call svdcmp(mat, eig_vals, eig_vecs)
        ret_mat   = matmul(eig_vecs, transpose(mat))
        ret_scale = sqrt(sum(var2**2) / sum(var1**2))
        mean2     = - ret_scale * matmul(mean1, transpose(ret_mat)) + mean2
        ret_trans = mean2(1,:)
    end subroutine atom_posmatch
end program simple_test_atom_posmatch
