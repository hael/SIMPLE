module pose_cont_refinement_test_helpers
use simple_defs, only: dp
use simple_cartesian_pose_refiner, only: cartesian_pose_refiner, cartesian_pose_data
use simple_type_defs, only: ctfparams, CTFFLAG_NO
implicit none
private

integer, parameter, public :: TEST_BOX = 24

public :: assert_true
public :: build_test_volume
public :: identity_rotation
public :: prepare_unweighted_particle
public :: rotation_distance

contains

    subroutine assert_true(condition,message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if( .not. condition ) error stop trim(message)
    end subroutine assert_true

    subroutine build_test_volume(volume)
        real, allocatable, intent(out) :: volume(:,:,:)
        real, parameter :: centres(3,4) = reshape([ &
            &-5.,-3., 2., 4., 5.,-3., 0.,-6.,-5., 3.,-2., 6.],[3,4])
        real, parameter :: sigmas(4) = [2.,2.5,1.8,2.2]
        real, parameter :: amplitudes(4) = [1.,0.8,0.6,0.5]
        real :: centre, dx, dy, dz
        integer :: blob, i, j, k

        allocate(volume(TEST_BOX,TEST_BOX,TEST_BOX),source=0.)
        centre = real(TEST_BOX)/2.+0.5
        do k = 1, TEST_BOX
            do j = 1, TEST_BOX
                do i = 1, TEST_BOX
                    do blob = 1, 4
                        dx = real(i)-centre-centres(1,blob)
                        dy = real(j)-centre-centres(2,blob)
                        dz = real(k)-centre-centres(3,blob)
                        volume(i,j,k) = volume(i,j,k)+amplitudes(blob)* &
                            &exp(-(dx*dx+dy*dy+dz*dz)/(2.*sigmas(blob)**2))
                    enddo
                enddo
            enddo
        enddo
    end subroutine build_test_volume

    pure function identity_rotation() result(rotation)
        real(dp) :: rotation(3,3)
        rotation = 0._dp
        rotation(1,1) = 1._dp
        rotation(2,2) = 1._dp
        rotation(3,3) = 1._dp
    end function identity_rotation

    subroutine prepare_unweighted_particle(workspace,observed,data,shell_range)
        type(cartesian_pose_refiner), intent(in) :: workspace
        complex, intent(in) :: observed(-TEST_BOX/2:,-TEST_BOX/2:)
        type(cartesian_pose_data), intent(out) :: data
        integer, intent(in), optional :: shell_range(2)
        type(ctfparams) :: no_ctf
        real :: sigma2(0:TEST_BOX/2)
        integer :: active_range(2)

        no_ctf%ctfflag = CTFFLAG_NO
        sigma2 = 1.
        active_range = [2,TEST_BOX/2]
        if( present(shell_range) ) active_range = shell_range
        call workspace%prepare_particle(observed,no_ctf,sigma2,active_range,data)
    end subroutine prepare_unweighted_particle

    pure function rotation_distance(left,right) result(distance)
        real(dp), intent(in) :: left(3,3), right(3,3)
        real(dp) :: distance, cosine
        cosine = 0.5_dp*(sum(left*right)-1._dp)
        distance = acos(max(-1._dp,min(1._dp,cosine)))
    end function rotation_distance

end module pose_cont_refinement_test_helpers
