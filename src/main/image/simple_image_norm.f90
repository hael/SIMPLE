submodule (simple_image) simple_image_norm
!$ use omp_lib
!$ use omp_lib_kinds
include  'simple_lib.f08'
#include "simple_local_flags.inc"
implicit none
contains

    module subroutine scale_pixels(self, new_range, ssc, oold_range)
        class(image),   intent(inout) :: self
        real,           intent(in)    :: new_range(2)
        real, optional, intent(out)   :: oold_range(2), ssc
        real :: old_range(2), sc
        old_range(1) = minval(self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)))
        old_range(2) = maxval(self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)))
        sc = (new_range(2) - new_range(1))/(old_range(2) - old_range(1))
        self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) = sc*self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))+new_range(1)-sc*old_range(1)
        if(present(ssc)) ssc = sc
        if(present(oold_range)) oold_range = old_range
    end subroutine scale_pixels

    module subroutine norm( self, a_s )
        class(image),   intent(inout) :: self
        real, optional, intent(in)    :: a_s(2)
        integer :: npix
        real    :: ave, var, ep
        npix   = product(self%ldim)
        ave    = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))) / real(npix)
        self%rmat = self%rmat - ave
        ep     = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
        var    = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))**2.0)
        var    = (var-ep**2./real(npix))/(real(npix)-1.) ! corrected two-pass formula
        if( is_a_number(var) )then
            if( var > 0. ) self%rmat = self%rmat / (sqrt(var))
        endif
        if( present(a_s) )then
            self%rmat = self%rmat * a_s(2)
            self%rmat = self%rmat + a_s(1)
        endif
    end subroutine norm

    module subroutine norm_minmax( self  )
        class(image), intent(inout) :: self
        real    :: smin, smax, delta
        smin  = minval(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
        smax  = maxval(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
        delta = smax - smin
        self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) =&
        &(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) - smin)/delta
    end subroutine norm_minmax

    module subroutine norm4viz( self, brightness, maxmin)
        class(image),      intent(inout) :: self
        real,    optional, intent(in)    :: brightness
        logical, optional, intent(in)    :: maxmin
        real    :: brightness_l
        logical :: maxmin_l
        brightness_l = 128.0
        maxmin_l     = .false.
        if(present(brightness)) brightness_l = brightness
        if(present(maxmin))     maxmin_l     = maxmin
        if(self%is_ft())THROW_HARD('real space only; norm4viz')
        if(maxmin_l) then
            call self%norm_minmax
            brightness_l = brightness_l - 128.0
            self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = brightness_l + 256 *&
            &self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))
        else
            call self%norm
            self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = brightness_l + 10.5 *& ! magic numbers from Joe
            &self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))
        endif
        ! thresholding
        where( self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) > 255. )
            self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = 255.
        else where( self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) < 0. )
            self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = 0.
        end where
    end subroutine norm4viz

    !> \brief norm_ext  is for normalization of an image using inputted average and standard deviation
    !! \param avg Average
    !! \param sdev Standard deviation
    module subroutine norm_ext( self, avg, sdev )
        class(image), intent(inout) :: self
        real,         intent(in)    :: avg, sdev
        if( self%ft )then
            THROW_WARN('cannot normalize FTs; norm_ext')
            return
        endif
        if( abs(avg) > TINY ) self%rmat = self%rmat - avg
        if( sdev     > 0.   ) self%rmat = self%rmat / sdev
    end subroutine norm_ext

    !> \brief normalize whole image to standardize image background (outside of logical mask)
    module subroutine norm_noise( self, lmsk, sdev_noise )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: lmsk(self%ldim(1),self%ldim(2),self%ldim(3)) ! foreground must be true
        real,         intent(inout) :: sdev_noise
        integer :: npix
        real    :: ave, var, ep
        npix = product(self%ldim) - count(lmsk) ! # background pixels
        ave  = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=.not. lmsk) / real(npix) ! background average
        if( abs(ave) > TINY ) self%rmat = self%rmat - ave
        ep         = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)),      mask=.not. lmsk)
        var        = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))**2.0, mask=.not. lmsk)
        var        = (var-ep**2./real(npix))/(real(npix)-1.) ! corrected two-pass formula
        sdev_noise = 0.
        if( is_a_number(var) )then
            sdev_noise = sqrt(var)
            if( var > 0. ) self%rmat = self%rmat / sdev_noise
        endif
    end subroutine norm_noise

    !> \brief normalize whole image to standardize image foreground (within logical mask)
    module subroutine norm_within( self, mask )
        class(image),      intent(inout) :: self
        logical,           intent(in)    :: mask(self%ldim(1),self%ldim(2),self%ldim(3))    ! foreground must be true
        real :: npix, ax, sxx
        npix = real(count(mask))
        ax   = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=mask) / npix
        self%rmat = self%rmat - ax
        sxx  = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))*self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=mask)
        sxx  = sxx/real(npix)
        if( sxx > TINY ) self%rmat = self%rmat / sqrt(sxx)
    end subroutine norm_within

    !>  \brief  cure_outliers for replacing extreme outliers with median of a 13x13 neighbourhood window
    !!          only done on negative values, assuming white ptcls on black bkgr
    !! \param ncured number of corrected points
    !! \param nsigma number of std. dev. to set upper and lower limits
    !! \param deadhot output index of corrected pixels
    !! \param outliers
    module subroutine cure_outliers( self, ncured, nsigma, deadhot, outliers )
        class(image),                   intent(inout) :: self
        integer,                        intent(inout) :: ncured
        real,                           intent(in)    :: nsigma
        integer,                        intent(out)   :: deadhot(2)
        logical, allocatable, optional, intent(out)   :: outliers(:,:)
        real, allocatable :: win(:,:), rmat_pad(:,:)
        real    :: ave, sdev, var, lthresh, uthresh
        integer :: i, j, hwinsz, winsz
        logical :: was_fted, err, present_outliers
        if( self%ldim(3)>1 )THROW_HARD('for 2D images only; cure_outliers')
        was_fted = self%is_ft()
        if( was_fted )THROW_HARD('for real space images only; cure_outliers')
        present_outliers = present(outliers)
        ncured   = 0
        hwinsz   = 6
        deadhot  = 0
        if( present_outliers )then
            if( allocated(outliers) ) deallocate(outliers)
            allocate( outliers(self%ldim(1),self%ldim(2)),source =.false. )
        endif
        call moment( self%rmat(1:self%ldim(1),1:self%ldim(2),1), ave, sdev, var, err )
        if( sdev<TINY )return
        lthresh = ave - real(nsigma) * sdev
        uthresh = ave + real(nsigma) * sdev
        if( any(self%rmat(1:self%ldim(1),1:self%ldim(2),1)<lthresh).or.&
            &any(self%rmat(1:self%ldim(1),1:self%ldim(2),1)>uthresh) )then
            winsz = 2*hwinsz+1
            allocate(rmat_pad(1-hwinsz:self%ldim(1)+hwinsz,1-hwinsz:self%ldim(2)+hwinsz), win(winsz,winsz))
            rmat_pad(:,:) = median( reshape(self%rmat(1:self%ldim(1),1:self%ldim(2),1), (/(self%ldim(1)*self%ldim(2))/)) )
            rmat_pad(1:self%ldim(1), 1:self%ldim(2)) = self%rmat(1:self%ldim(1),1:self%ldim(2),1)
            !$omp parallel do collapse(2) schedule(static) default(shared) private(i,j,win)&
            !$omp reduction(+:ncured,deadhot) proc_bind(close)
            do i=1,self%ldim(1)
                do j=1,self%ldim(2)
                    if( self%rmat(i,j,1)<lthresh .or. self%rmat(i,j,1)>uthresh )then
                        if( present_outliers )then
                            outliers(i,j)=.true.
                            if (self%rmat(i,j,1)<lthresh) deadhot(1) = deadhot(1) + 1
                            if (self%rmat(i,j,1)>uthresh) deadhot(2) = deadhot(2) + 1
                        else
                            win = rmat_pad( i-hwinsz:i+hwinsz, j-hwinsz:j+hwinsz )
                            self%rmat(i,j,1) = median( reshape(win,(/winsz**2/)) )
                            ncured = ncured + 1
                        endif
                    endif
                enddo
            enddo
            !$omp end parallel do
            deallocate( win, rmat_pad )
        endif
    end subroutine cure_outliers

    module subroutine quantize_fwd( self, nquanta, transl_tab, l_msk )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: nquanta
        real,              intent(inout) :: transl_tab(nquanta)
        logical, optional, intent(in)    :: l_msk(self%ldim(1),self%ldim(2),self%ldim(3)) 
        real, allocatable :: pixvals(:)
        integer :: i, j, k, ind
        real    :: dist
        if( present(l_msk) )then
            pixvals = pack(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=l_msk)
        else
            pixvals = pack(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=.true.)
        endif
        if( self%ldim(3) > 1 )then
            call quantize_vec(pixvals, nquanta, transl_tab)
            !$omp parallel do schedule(static) default(shared) private(k,j,i,ind,dist) collapse(3) proc_bind(close)
            do k = 1,self%ldim(3)
                do j = 1,self%ldim(2)
                    do i = 1,self%ldim(1)
                        call find(transl_tab, nquanta, self%rmat(i,j,k), ind, dist)
                        self%rmat(i,j,k) = real(ind - 1)
                    end do
                end do
            end do
            !$omp end parallel do
        else
            call quantize_vec_serial(pixvals, nquanta, transl_tab)
            do j = 1,self%ldim(2)
                do i = 1,self%ldim(1)
                    call find(transl_tab, nquanta, self%rmat(i,j,1), ind, dist)
                    self%rmat(i,j,1) = real(ind - 1)
                end do
            end do
        endif
    end subroutine quantize_fwd

    module subroutine quantize_bwd( self, nquanta, transl_tab )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: nquanta
        real,         intent(in)    :: transl_tab(nquanta)
        real    :: pixvals(nquanta)
        integer :: i, j, k, ind
        real    :: dist
        pixvals = real((/(k,k=1,nquanta)/)) - 1.0
        if( self%ldim(3) > 1 )then
            !$omp parallel do schedule(static) default(shared) private(k,j,i,ind,dist) collapse(3) proc_bind(close)
            do k = 1,self%ldim(3)
                do j = 1,self%ldim(2)
                    do i = 1,self%ldim(1)
                        call find(pixvals, nquanta, self%rmat(i,j,k), ind, dist)
                        self%rmat(i,j,k) = transl_tab(ind)
                    end do
                end do
            end do
            !$omp end parallel do
        else
            do j = 1,self%ldim(2)
                do i = 1,self%ldim(1)
                    call find(pixvals, nquanta, self%rmat(i,j,1), ind, dist)
                    self%rmat(i,j,1) = transl_tab(ind)
                end do
            end do
        endif
    end subroutine quantize_bwd

end submodule simple_image_norm
