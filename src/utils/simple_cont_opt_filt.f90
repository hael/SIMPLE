module simple_cont_opt_filt
    use simple_defs
    use simple_image,      only: image
    use simple_math,       only: hyp
    use simple_opt_filter, only: butterworth_filter
    implicit none
    type(image) :: odd_img, even_img
    
    !>  \brief  defines the function interface
    abstract interface
        function fun_interface( fun_self, vec, D ) result( cost )
            class(*), intent(inout) :: fun_self
            integer,  intent(in)    :: D
            real,     intent(in)    :: vec(D)
            real                    :: cost
        end function
    end interface
    
    contains
        ! use even as the target and odd as the object to be convolved with Butterworth Kernel
        function filt_cost( fun_self, x, d ) result( r )
            class(*), intent(inout)  :: fun_self
            integer,  intent(in)     :: d
            real,     intent(in)     :: x(d)
            real                     :: r
            integer                  :: ldim(3), k, l
            type(image)              :: img_ker
            integer,     parameter   :: BW_ORDER = 8
            real,        allocatable :: cur_diff(:,:,:), cur_filt(:), x_mat(:, :)
            ldim = odd_img%get_ldim()
            allocate(cur_diff(ldim(1),ldim(2),1), x_mat(ldim(1),ldim(2)), cur_filt(ldim(1)), source=0.)
            x_mat = reshape(x, [ldim(1), ldim(2)])
            r     = 0.
            do l = ldim(2)/2-2, ldim(2)/2+2
                do k = ldim(1)/2-2, ldim(1)/2+2
                    call img_ker%copy(even_img)
                    call butterworth_filter(cur_filt, BW_ORDER, x_mat(k,l))
                    call img_ker%apply_filter(cur_filt)
                    call img_ker%sqeuclid_matrix(odd_img, cur_diff)
                    r = r + cur_diff(k,l,1)
                enddo
            enddo
        end function

        subroutine filt_gcost( fun_self, x, grad, d )
            class(*), intent(inout)  :: fun_self
            integer,  intent(in)     :: d
            real,     intent(inout)  :: x(d)
            real,     intent(out)    :: grad(d)
            integer                  :: ldim(3), k, l
            type(image)              :: img_plus, img_minus
            integer,     parameter   :: BW_ORDER = 8
            real,        parameter   :: EPS = 0.1
            real,        allocatable :: cur_diff(:,:,:), cur_filt(:), x_mat(:, :), grad_mat(:, :)
            ldim = odd_img%get_ldim()
            allocate(cur_diff(ldim(1),ldim(2),1), grad_mat(ldim(1),ldim(2)), x_mat(ldim(1),ldim(2)), cur_filt(ldim(1)), source=0.)
            x_mat = reshape(x, [ldim(1), ldim(2)])
            do l = ldim(2)/2-2, ldim(2)/2+2
                do k = ldim(1)/2-2, ldim(1)/2+2
                    call img_minus%copy(even_img)
                    call img_plus%copy(even_img)
                    call butterworth_filter(cur_filt, BW_ORDER, x_mat(k,l)-EPS)
                    call img_minus%apply_filter(cur_filt)
                    call img_minus%sqeuclid_matrix(odd_img, cur_diff)
                    grad_mat(k,l) = cur_diff(k,l,1)
                    call butterworth_filter(cur_filt, BW_ORDER, x_mat(k,l)+EPS)
                    call img_plus%apply_filter(cur_filt)
                    call img_plus%sqeuclid_matrix(odd_img, cur_diff)
                    grad_mat(k,l) = (cur_diff(k,l,1)-grad_mat(k,l))/(2*EPS)
                enddo
            enddo
            grad = reshape(grad_mat, [product(ldim)])
        end subroutine
end module