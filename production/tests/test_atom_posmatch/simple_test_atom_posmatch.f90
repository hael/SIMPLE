program simple_test_atom_posmatch
include 'simple_lib.f08'
use simple_nanoparticle_utils, only: atoms_register, Kabsch_algo
implicit none
integer, parameter   :: N1 = 10, N2 = 5, MAX_POS = 10
integer, allocatable :: inds(:), select_inds(:)
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
select_inds  = scramble(N1)
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
inds = scramble(N2)
do i = 1, N2
    atom2_rnd(:,i)   = atom2_pos(:,inds(i))
    mapping(inds(i)) = i
enddo
print *, '------------'
print *, 'randomized inds = ', inds
! generating all permutations of atom2_rnd positions
print *, '------ ROTATED ATOMS 1 POSITIONS (SOME WILL MATCH ATOMS 2 POSITIONS BELOW) ------'
call atoms_register(atom1_pos, atom2_rnd, atom1_rot, verbose=.false.)
do i = 1, N1
    print *, 'rotated atom 1 pos : ', atom1_rot(:,i)
enddo
print *, '------ ATOMS 2 POSITIONS ------'
do i = 1, N2
    print *, 'atom 2 pos : ', atom2_rnd(:,mapping(i))
enddo
print *, '------ ROTATED ATOMS 2 POSITIONS (SOME WILL MATCH ATOMS 1 POSITIONS BELOW) ------'
call atoms_register(atom2_rnd, atom1_pos, atom2_rot, verbose=.false.)
do i = 1, N2
    print *, 'rotated atom 2 pos : ', atom2_rot(:,mapping(i))
enddo
print *, '------ ATOMS 1 POSITIONS ------'
do i = 1, N1
    print *, 'atom 1 pos : ', atom1_pos(:,i)
enddo

contains

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
