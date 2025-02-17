program simple_test_atom_posmatch
include 'simple_lib.f08'
implicit none
integer, parameter   :: N1 = 5, N2 = 5, MAX_POS = 10
integer, allocatable :: inds(:), perm(:,:)
real,    allocatable :: costs(:)
real    :: atom1_pos(3,N1), atom2_pos(3,N2), rot_mat(3,3), scale, translation(3), atom2_rnd(3,N2)
real    :: rec_mat(3,3), rec_trans(3), rec_scale, rec_atom2(3,N2), tmp_atom2(3,N2)
integer :: i, j, k, l, array(N2 - 1), cnt, N, counter
call seed_rnd
do i = 1, N1
    do j = 1, 3
        atom1_pos(j,i) = floor(ran3() * real(MAX_POS))
    enddo
    print *, 'atom1 pos ', i, ' = ', atom1_pos(:,i)
enddo
! atom1_pos(:,1) = [4., 8., 8.]
! atom1_pos(:,2) = [8., 4., 4.]
! atom1_pos(:,3) = [8., 5., 3.]
! atom1_pos(:,4) = [7., 1., 7.]
! atom1_pos(:,5) = [6., 6., 6.]
rot_mat(1,:) = [-0.4297,  0.3140, -0.8466]
rot_mat(2,:) = [-0.7475, -0.6496,  0.1385]
rot_mat(3,:) = [-0.5065,  0.6924,  0.5139]
translation  = [0.2, -0.5, 1.]
scale        = 0.5
atom2_pos    = scale * matmul(rot_mat, atom1_pos)
do i = 1, N2
    atom2_pos(:,i) = atom2_pos(:,i) + translation
enddo
print *, '------ TRUTH POSITIONS ------'
call Kabsch_algo(atom1_pos(:,1:3), atom2_pos(:,1:3), rec_mat, rec_trans, rec_scale)
rec_atom2 = rec_scale * matmul(rec_mat, atom1_pos)
print *, 'atom 1 pos : ',      atom1_pos(:,4)
print *, 'atom 2 pos : ',      atom2_pos(:,4)
print *, 'rec scale = ',       rec_scale
print *, 'rec translation = ', rec_trans
print *, 'rotated atom 2 : ',  rec_atom2(:,4) + rec_trans
print *, 'det = ', determinant(rec_mat)
print *, 'ret_mat 1 = ', rec_mat(1,:)
print *, 'ret_mat 2 = ', rec_mat(2,:)
print *, 'ret_mat 3 = ', rec_mat(3,:)
! randomize atom2_pos
inds = scramble(N2)
do i = 1, N2
    atom2_rnd(:,i) = atom2_pos(:,inds(i))
enddo
print *, '------------'
print *, 'inds = ', inds
print *, '------ WRONG POSITIONS ------'
call Kabsch_algo(atom1_pos(:,1:3), atom2_rnd(:,1:3), rec_mat, rec_trans, rec_scale)
rec_atom2 = rec_scale * matmul(rec_mat, atom1_pos)
print *, 'atom 1 pos : ',      atom1_pos(:,1)
print *, 'atom 2 pos : ',      atom2_rnd(:,1)
print *, 'rec scale = ',       rec_scale
print *, 'rec translation = ', rec_trans
print *, 'rotated atom 2 : ',  rec_atom2(:,1) + rec_trans
! generating all permutations of atom2_rnd positions
print *, '------ RIGHT POSITIONS ------'
counter = 1
allocate(costs(N2*(N2-1)*(N2-2)),perm(3,N2*(N2-1)*(N2-2)))
costs = 0.
do i = 1, N2
    do j = 1, N2
        do k = 1, N2
            if( k == j .or. i == j .or. i == k )cycle
            perm(:,counter) = [i,j,k]
            ! print *, 'perm = ', perm(:,counter)
            call Kabsch_algo(atom1_pos(:,1:3), atom2_rnd(:,perm(:,counter)), rec_mat, rec_trans, rec_scale)
            do l = 1, N2
                if( l == perm(1,counter) .or. l == perm(2,counter) .or. l == perm(3,counter) )cycle
                rec_atom2(:,l) = rec_scale * matmul(rec_mat, atom1_pos(:,4)) + rec_trans
                tmp_atom2(:,l) = rec_scale * matmul(rec_mat, atom1_pos(:,5)) + rec_trans
                if( sum((tmp_atom2(:,l) - atom2_rnd(:,l))**2) < sum((rec_atom2(:,l) - atom2_rnd(:,l))**2) )then
                    costs(counter) = costs(counter) + sum((tmp_atom2(:,l) - atom2_rnd(:,l))**2)
                else
                    costs(counter) = costs(counter) + sum((rec_atom2(:,l) - atom2_rnd(:,l))**2)
                endif
            enddo
            ! print *, 'counter = ', counter, '; cost = ', costs(counter)
            counter = counter + 1
        enddo
    enddo
enddo
i = minloc(costs, dim=1)
print *, 'i = ', i
print *, 'CORRECT PERM = ', perm(:,i)
call Kabsch_algo(atom1_pos(:,1:3), atom2_rnd(:,perm(:,i)), rec_mat, rec_trans, rec_scale)
do l = 1, N2
    if( l == perm(1,i) .or. l == perm(2,i) .or. l == perm(3,i) )cycle
    rec_atom2(:,l) = rec_scale * matmul(rec_mat, atom1_pos(:,4)) + rec_trans
    tmp_atom2(:,l) = rec_scale * matmul(rec_mat, atom1_pos(:,5)) + rec_trans
    print *, '------'
    if( sum((tmp_atom2(:,l) - atom2_rnd(:,l))**2) < sum((rec_atom2(:,l) - atom2_rnd(:,l))**2) )then
        print *, 'atom 1 pos : ',      atom1_pos(:,4)
        print *, 'atom 2 pos : ',      atom2_rnd(:,l)
        print *, 'rec scale = ',       rec_scale
        print *, 'rec translation = ', rec_trans
        print *, 'rotated atom 2 : ',  tmp_atom2(:,l)
    else
        print *, 'atom 1 pos : ',      atom1_pos(:,5)
        print *, 'atom 2 pos : ',      atom2_rnd(:,l)
        print *, 'rec scale = ',       rec_scale
        print *, 'rec translation = ', rec_trans
        print *, 'rotated atom 2 : ',  rec_atom2(:,l)
    endif
enddo
contains

    real function determinant(A)
        real :: A(3,3)
        determinant =   A(1,1)*A(2,2)*A(3,3)  &
                      - A(1,1)*A(2,3)*A(3,2)  &
                      - A(1,2)*A(2,1)*A(3,3)  &
                      + A(1,2)*A(2,3)*A(3,1)  &
                      + A(1,3)*A(2,1)*A(3,2)  &
                      - A(1,3)*A(2,2)*A(3,1)
        return
    end function determinant

    subroutine Kabsch_algo( pos1, pos2, ret_mat, ret_trans, ret_scale)
        real,    intent(in)    :: pos1(:,:), pos2(:,:)
        real,    intent(inout) :: ret_mat(3,3)
        real,    intent(inout) :: ret_trans(3)
        real,    intent(inout) :: ret_scale
        real,    allocatable   :: var1(:,:), var2(:,:)
        integer, allocatable   :: inds(:)
        real    :: mean1(1,3), mean2(1,3), mat(3,3), eig_vals(3), eig_vecs(3,3), det, eye(3,3)
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
        ! print *, 'mat 1 = ', mat(1,:)
        ! print *, 'mat 2 = ', mat(2,:)
        ! print *, 'mat 3 = ', mat(3,:)
        call svdcmp(mat, eig_vals, eig_vecs)
        inds = (/(i,i=1,3)/)
        call hpsort(eig_vals, inds)
        call reverse(eig_vals)
        call reverse(inds)
        eig_vecs = eig_vecs(:,inds)
        mat      = mat(:,inds)
        ret_mat  = matmul(eig_vecs, transpose(mat))
        ! mat(:,2) = -mat(:,2)
        ! eig_vecs(:,2) = - eig_vecs(:,2)
        ! eig_vecs(:,3) = - eig_vecs(:,3)
        ! print *, 'eig vals = ', eig_vals
        ! print *, 'eig_vecs 1 = ', eig_vecs(1,:)
        ! print *, 'eig_vecs 2 = ', eig_vecs(2,:)
        ! print *, 'eig_vecs 3 = ', eig_vecs(3,:)
        ! print *, 'mat 1 = ', mat(1,:)
        ! print *, 'mat 2 = ', mat(2,:)
        ! print *, 'mat 3 = ', mat(3,:)
        ! print *, 'det = ', determinant(ret_mat)
        ! reflection special case
        if( determinant(ret_mat) < 0. )then
            eye      = 0.
            eye(1,1) =  1.
            eye(2,2) =  1.
            eye(3,3) = -1.
            ret_mat  = matmul(matmul(eig_vecs, eye), transpose(mat))
        endif
        ret_scale = sqrt(sum(var2**2) / sum(var1**2))
        mean2     = - ret_scale * matmul(mean1, transpose(ret_mat)) + mean2
        ret_trans = mean2(1,:)
    end subroutine Kabsch_algo

    function scramble( number_of_values ) result(array)
        !@(#) M_random::scramble(3f): return integer array of random values 1 to N.
        integer,intent(in)    :: number_of_values
        integer,allocatable   :: array(:)
        integer               :: i, j, k, m, n
        integer               :: temp
        real                  :: u
        array=[(i,i=1,number_of_values)]
        ! The intrinsic RANDOM_NUMBER(3f) returns a real number (or an array
        ! of such) from the uniform distribution over the interval [0,1). (ie.
        ! it includes 0 but not 1.).
        !
        ! To have a discrete uniform distribution on
        ! the integers {n, n+1, ..., m-1, m} carve the continuous distribution
        ! up into m+1-n equal sized chunks, mapping each chunk to an integer.
        !
        ! One way is:
        !   call random_number(u)
        !   j = n + FLOOR((m+1-n)*u)  ! choose one from m-n+1 integers
        n=1
        m=number_of_values
        do k=1,2
            do i=1,m
            call random_number(u)
            j = n + FLOOR((m+1-n)*u)
            ! switch values
            temp=array(j)
            array(j)=array(i)
            array(i)=temp
            enddo
        enddo 
    end function scramble

end program simple_test_atom_posmatch
