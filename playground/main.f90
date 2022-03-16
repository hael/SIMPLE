!> Translating MATLAB codes of TV-Denoising into Fortran
!!
!! @author Cong
program main
    implicit none (type, external)
    external :: sgesv
    real     :: in_img(31, 31)  ! input image/matrix
    real     :: in_img_1D(31*31)
    real     :: id_mat(31*31, 31*31)
    real     :: der_mat(31*31 - 1, 31*31 - 1), adjoint(31*31 - 1, 31*31 - 1), banded_mat(31*31 - 1, 31*31 - 1)
    real     :: der_x(31*31 - 1), der_y(31*31 - 1)
    real     :: max_radius = 10
    real     :: min_radius = 7
    real     :: dist
    real     :: lambda = 0.01
    real     :: cur_1D(31*31)
    integer  :: x,y,i, row_ind, col_ind, Niter = 50

    do x = 1, 31
        do y = 1, 31
            dist = sqrt((x-16.0)**2 + (y-16.0)**2)
            if (dist < max_radius .and. dist > min_radius) then
                in_img(x,y) = 1
            end if
        end do
    end do

    in_img_1D = reshape(in_img, [31*31])

    do i = 1, 31*31
        id_mat(i,i) = 1
    end do

    do row_ind = 1, 31*31 - 1
        do col_ind = 1, 31*31 - 1
            der_mat(row_ind, col_ind) = id_mat(row_ind+1, col_ind) - id_mat(row_ind, col_ind)
        end do
    end do

    adjoint = matmul(der_mat, transpose(der_mat))

    cur_1D = in_img_1D
    der_x  = matmul(der_mat, cur_1D)
    der_y  = matmul(der_mat, cur_1D)
    do i = 1, Niter
        do j = 1, 31*31-1
            banded_mat(j,j) = abs(der_x(j))/lambda
        end do
    end do

    
    !a = reshape([ 2., 3., 1., 1. ], [ 2, 2 ])
    !b = [ 5., 6. ]

    !call sgesv(2, 1, a, 2, pivot, b, 2, rc)
    print *, 'temp = ', adjoint(1,1)

end program main