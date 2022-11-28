! from "Uniform rotations from Gaussians" of https://www.sciencedirect.com/science/article/pii/B9780080507552500361
program simple_test_uniform_rot
include 'simple_lib.f08'
implicit none
integer, parameter  :: N_SAMPLES = 4, N_ARR = 2000
real,    parameter  :: SIGMA = 1.
real    :: unif(N_SAMPLES), u1, u2, mat(3,3), q0, q1, q2, q3
integer :: i, j
call srand(time())
do i = 1, N_ARR
    do j = 1, N_SAMPLES
        call rgauss(SIGMA, u1, u2)
        unif(j) = u1
    enddo
    unif = unif/sqrt(sum(unif**2))
    ! print "(f20.15, f20.15, f20.15, f20.15)", unif
    q0 = unif(1)
    q1 = unif(2)
    q2 = unif(3)
    q3 = unif(4)
    ! First row of the rotation matrix
    mat(1,1) = 2 * (q0 * q0 + q1 * q1) - 1;
    mat(1,2) = 2 * (q1 * q2 - q0 * q3);
    mat(1,3) = 2 * (q1 * q3 + q0 * q2);    
    ! Second row of the rotation matrix
    mat(2,1) = 2 * (q1 * q2 + q0 * q3);
    mat(2,2) = 2 * (q0 * q0 + q2 * q2) - 1;
    mat(2,3) = 2 * (q2 * q3 - q0 * q1);
    ! Third row of the rotation matrix
    mat(3,1) = 2 * (q1 * q3 - q0 * q2);
    mat(3,2) = 2 * (q2 * q3 + q0 * q1);
    mat(3,3) = 2 * (q0 * q0 + q3 * q3) - 1;
    print "(f20.15, f20.15, f20.15)", matmul([1., 0., 0.], mat);
    ! print *, mat(1,:)
    ! print *, mat(2,:)
    ! print *, mat(3,:)
enddo
contains
    ! Box—Muller method
    subroutine rgauss(sig, y1, y2)
        real, intent(in)    :: sig
        real, intent(inout) :: y1, y2
        real :: x1, x2, w
        w = 0.
        do while ( (w .ge. 1.0).or.(w.eq.0) )
            x1 = 2.0 * rand(0) - 1.0
            x2 = 2.0 * rand(0) - 1.0
            w  = x1 * x1 + x2 * x2
        end do
        w = sig*sqrt( (-2.0 * log( w ) ) / w )
        y1 = x1 * w
        y2 = x2 * w
    end subroutine rgauss
end program simple_test_uniform_rot