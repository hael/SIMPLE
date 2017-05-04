module simple_gridding
use simple_image,   only: image
use simple_kbinterpol ! use all in there
implicit none

contains
    
    ! bails when zero images

    !>  \brief  prepare image for gridding interpolation in Fourier space
    subroutine prep4cgrid( img, img4grid, msk )
        class(image), intent(inout)  :: img, img4grid
        real,         intent(in)     :: msk
        real                         :: ave, sdev, var, med
        logical                      :: err
        call img%bwd_ft                                           ! make sure not FTed
        med = img%median_pixel()
        call img%pad(img4grid, backgr=med)                        ! padding in real space                     
        call divide_w_instr(img4grid)                             ! division w instr in real space
        call img4grid%fwd_ft                                      ! return the Fourier transform
    end subroutine prep4cgrid
    
    !> \brief  for dividing a real or complex image with the instrument function
    subroutine divide_w_instr( img )
        use simple_jiffys, only: alloc_err
        class(image), intent(inout) :: img
        real, allocatable :: w1(:), w2(:), w3(:)
        integer :: ldim(3), i, j, k, alloc_stat, lims(3,2)
        real    :: ci, cj, ck, arg
        ! get the limits
        ldim = img%get_ldim()      
        if( img%is_ft() )then
            lims = img%loop_lims(2)
        else
            lims(:,1) = 1
            lims(:,2) = ldim
        endif
        ! make the window
        allocate( w1(lims(1,1):lims(1,2)), w2(lims(2,1):lims(2,2)),&
        w3(lims(3,1):lims(3,2)), stat=alloc_stat )
        call alloc_err("In: divide_w_instr; simple_gridding", alloc_stat)
        ! calculate the values
        ci = -real(ldim(1))/2.
        do i=lims(1,1),lims(1,2)
            if( img%is_ft() )then
                arg = real(i)/real(ldim(1))
            else
                arg = ci/real(ldim(1))
            endif
            w1(i) = kb_instr(arg)
            ci = ci+1.
        end do
        if( img%square_dims() )then
            w2 = w1
            if( img%is_3d() ) w3 = w1
        else
            cj = -real(ldim(2))/2.
            do j=lims(2,1),lims(2,2)
                if( img%is_ft() )then
                    arg = real(j)/real(ldim(2))
                else
                    arg = cj/real(ldim(2))
                endif
                w2(j) = kb_instr(arg)
                cj = cj+1.
            end do
            ck = -real(ldim(3))/2.
            do k=lims(3,1),lims(3,2)
                if( img%is_ft() )then
                    arg = real(k)/real(ldim(3))
                else
                    arg = ck/real(ldim(3))
                endif
                w3(k) = kb_instr(arg)
                ck = ck+1.
            end do
        endif
        if( img%is_2d() ) w3 = 1.
        ! divide the image
        !$omp parallel do schedule(auto) default(shared) private(i,j,k)
        do i=lims(1,1),lims(1,2)
            if( w1(i) == 0. ) cycle
            do j=lims(2,1),lims(2,2)
                if( w2(j) == 0. ) cycle
                do k=lims(3,1),lims(3,2)
                    if( w3(k) == 0. ) cycle
                    call img%div([i,j,k], w1(i)*w2(j)*w3(k))
                end do
            end do
        end do
        !$omp end parallel do
        deallocate(w1,w2,w3)
    end subroutine divide_w_instr
    
end module simple_gridding
