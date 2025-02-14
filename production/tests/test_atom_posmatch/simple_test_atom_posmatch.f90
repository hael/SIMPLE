program simple_test_atom_posmatch
include 'simple_lib.f08'
implicit none
real, parameter :: THETA = PI/6
real    :: u(3), v(3), rot_mat(3,3), u_rec(3), mat_inv(3,3), axis(3), u_rot(3)
u    = [1., 2., 3.]
v    = [1., 3., 5.]
axis = cross(u,v)
print *, 'original vector: ', u
call rotation_matrix_2(axis, theta, rot_mat)
u_rot = matmul(rot_mat, u)
print *, 'rotated vector: ', u_rot
call rotation_matrix_1(u_rot, u, mat_inv)
u_rec = matmul(mat_inv, u_rot/norm(u_rot)) * norm(u)
print *, 'reconstructed vector: ', u_rec

contains

    function norm(u) result(val)
        real, intent(in) :: u(3)
        real :: val
        val = sqrt(sum(u**2))
    end function norm

    subroutine rotation_matrix_1( u, v, rot_mat)
        real, intent(in)    :: u(3), v(3)
        real, intent(inout) :: rot_mat(3,3)
        real :: theta, axis(3)
        ! Calculate the rotation axis
        axis  = cross(u, v)
        ! Calculate the rotation angle
        theta = acos(dot_product(u, v) / norm(u) / norm(v))
        call rotation_matrix_2(axis, theta, rot_mat)
    end subroutine rotation_matrix_1

    subroutine rotation_matrix_2( axis, theta, rot_mat)
        real, intent(in)    :: axis(3)
        real, intent(in)    :: theta
        real, intent(inout) :: rot_mat(3,3)
        real :: u(3)
        u = axis / norm(axis)
        ! Generate rotation matrix using Rodrigues' formula
        rot_mat(1,1) = cos(theta) + u(1)**2 * (1 - cos(theta))
        rot_mat(1,2) = u(1) * u(2) * (1 - cos(theta)) - u(3) * sin(theta)
        rot_mat(1,3) = u(1) * u(3) * (1 - cos(theta)) + u(2) * sin(theta)
        rot_mat(2,1) = u(2) * u(1) * (1 - cos(theta)) + u(3) * sin(theta)
        rot_mat(2,2) = cos(theta) + u(2)**2 * (1 - cos(theta))
        rot_mat(2,3) = u(2) * u(3) * (1 - cos(theta)) - u(1) * sin(theta)
        rot_mat(3,1) = u(3) * u(1) * (1 - cos(theta)) - u(2) * sin(theta)
        rot_mat(3,2) = u(3) * u(2) * (1 - cos(theta)) + u(1) * sin(theta)
        rot_mat(3,3) = cos(theta) + u(3)**2 * (1 - cos(theta))
    end subroutine rotation_matrix_2

end program simple_test_atom_posmatch
