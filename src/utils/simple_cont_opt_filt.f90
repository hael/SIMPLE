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
            real                     :: r, smpd
            integer                  :: ldim(3), k, l, ext, k1, l1, sh
            type(image)              :: img_ker, filt_img, img_pad
            integer,     parameter   :: BW_ORDER = 8
            real,        allocatable :: cur_filt(:), x_mat(:, :)
            real,        pointer     :: img_rmat(:,:,:), filt_rmat(:,:,:)
            ldim = odd_img%get_ldim()
            allocate(x_mat(ldim(1),ldim(2)), cur_filt(ldim(1)), source=0.)
            x_mat = reshape(x, [ldim(1), ldim(2)])
            r     = 0.
            smpd  = 1.
            do l = 1,ldim(2)
                do k = 1,ldim(1)
                    call img_ker%copy(even_img)
                    call butterworth_filter(cur_filt, BW_ORDER, x_mat(k,l))
                    ext = nint(minval([real(x_mat(k,l)) + 10, real(ldim(1)/2)]))
                    call filt_img%new([2*ext, 2*ext, 1], smpd)
                    call filt_img%set_ft(.true.)
                    call filt_img%set_cmat((1., 0.))
                    do l1 = -ext, ext-1
                        do k1 = 0, ext
                            sh = nint(hyp(real(k1),real(l1),0.)*ldim(1)/(2*ext))
                            if( sh == 0 )then 
                                call filt_img%mul([k1,l1,0], maxval(cur_filt))
                            elseif( sh <= ldim(1) )then
                                call filt_img%mul([k1,l1,0], cur_filt(sh))
                            else
                                call filt_img%mul([k1,l1,0], 0.)
                            endif
                        enddo
                    enddo
                    call filt_img%ifft()
                    call filt_img%get_rmat_ptr(filt_rmat)
                    call img_pad%new([ldim(1) + 2*ext, ldim(2) + 2*ext,1], smpd)
                    call img_ker%pad(img_pad)
                    call img_pad%get_rmat_ptr(img_rmat)                    
                    r = r + abs(sum(filt_rmat(1:2*ext, 1:2*ext,1)*img_rmat(k:k+2*ext-1,l:l+2*ext-1,1)) - odd_img%get_rmat_at(k,l,1))
                enddo
            enddo
            write(*, *) 'cost = ', r
        end function

        subroutine filt_gcost( fun_self, x, grad, d )
            class(*), intent(inout)  :: fun_self
            integer,  intent(in)     :: d
            real,     intent(inout)  :: x(d)
            real,     intent(out)    :: grad(d)
            integer                  :: ldim(3), k, l, k1, l1, sh, ext
            real                     :: smpd
            type(image)              :: img_plus, img_minus, filt_img, img_pad
            integer,     parameter   :: BW_ORDER = 8
            real,        parameter   :: EPS = 0.1
            real,        allocatable :: cur_diff(:,:,:), cur_filt(:), x_mat(:, :), grad_mat(:, :)
            real,        pointer     :: img_rmat(:,:,:), filt_rmat(:,:,:)
            ldim = odd_img%get_ldim()
            smpd = 1.
            allocate(cur_diff(ldim(1),ldim(2),1), grad_mat(ldim(1),ldim(2)), x_mat(ldim(1),ldim(2)), cur_filt(ldim(1)), source=0.)
            x_mat = reshape(x, [ldim(1), ldim(2)])
            do l = 1,ldim(2)
                do k = 1,ldim(1)
                    call img_minus%copy(even_img)
                    call img_plus%copy(even_img)
                    call butterworth_filter(cur_filt, BW_ORDER, x_mat(k,l)-EPS)
                    ext = nint(minval([real(x_mat(k,l)) + 10, real(ldim(1)/2)]))
                    call filt_img%new([2*ext, 2*ext, 1], smpd)
                    call filt_img%set_ft(.true.)
                    call filt_img%set_cmat((1., 0.))
                    do l1 = -ext, ext-1
                        do k1 = 0, ext
                            sh = nint(hyp(real(k1),real(l1),0.)*ldim(1)/(2*ext))
                            if( sh == 0 )then 
                                call filt_img%mul([k1,l1,0], maxval(cur_filt))
                            elseif( sh <= ldim(1) )then
                                call filt_img%mul([k1,l1,0], cur_filt(sh))
                            else
                                call filt_img%mul([k1,l1,0], 0.)
                            endif
                        enddo
                    enddo
                    call filt_img%ifft()
                    call filt_img%get_rmat_ptr(filt_rmat)
                    call img_pad%new([ldim(1) + 2*ext, ldim(2) + 2*ext,1], smpd)
                    call img_minus%pad(img_pad)
                    call img_pad%get_rmat_ptr(img_rmat)
                    grad_mat(k,l) = abs(sum(filt_rmat(1:2*ext, 1:2*ext,1)*img_rmat(k:k+2*ext-1,l:l+2*ext-1,1)) - odd_img%get_rmat_at(k,l,1))
                    call butterworth_filter(cur_filt, BW_ORDER, x_mat(k,l)+EPS)
                    ext = nint(minval([real(x_mat(k,l)) + 10, real(ldim(1)/2)]))
                    call filt_img%new([2*ext, 2*ext, 1], smpd)
                    call filt_img%set_ft(.true.)
                    call filt_img%set_cmat((1., 0.))
                    do l1 = -ext, ext-1
                        do k1 = 0, ext
                            sh = nint(hyp(real(k1),real(l1),0.)*ldim(1)/(2*ext))
                            if( sh == 0 )then 
                                call filt_img%mul([k1,l1,0], maxval(cur_filt))
                            elseif( sh <= ldim(1) )then
                                call filt_img%mul([k1,l1,0], cur_filt(sh))
                            else
                                call filt_img%mul([k1,l1,0], 0.)
                            endif
                        enddo
                    enddo
                    call filt_img%ifft()
                    call filt_img%get_rmat_ptr(filt_rmat)
                    call img_pad%new([ldim(1) + 2*ext, ldim(2) + 2*ext,1], smpd)
                    call img_plus%pad(img_pad)
                    call img_pad%get_rmat_ptr(img_rmat)
                    grad_mat(k,l) = (abs(sum(filt_rmat(1:2*ext, 1:2*ext,1)*img_rmat(k:k+2*ext-1,l:l+2*ext-1,1)) - odd_img%get_rmat_at(k,l,1))-grad_mat(k,l))/(2*EPS)
                enddo
            enddo
            grad = reshape(grad_mat, [product(ldim)])
        end subroutine
end module