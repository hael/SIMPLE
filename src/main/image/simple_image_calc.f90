submodule (simple_image) simple_image_calc
!$ use omp_lib
!$ use omp_lib_kinds
include  'simple_lib.f08'
#include "simple_local_flags.inc"
implicit none
contains

    !===========================
    ! Basic stats / local stats
    !===========================

    module function get_sum_int(self) result(sum_int)
        class(image), intent(in) :: self
        real :: sum_int
        sum_int = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
    end function get_sum_int

    module function minmax( self, radius )result( mm )
        class(image),   intent(in) :: self
        real, optional, intent(in) :: radius
        real    :: mm(2), radsq
        integer :: c(3),i,j,k,dksq,djsq,dsq
        if( present(radius) )then
            radsq = radius**2
            c     = self%ldim/2+1
            if( self%ldim(3)==1 ) c(3)=1
            mm    = [huge(0.),-huge(0.)]
            do k = 1,self%ldim(3)
                dksq = (k-c(3))**2
                do j = 1,self%ldim(2)
                    djsq = dksq + (j-c(2))**2
                    do i = 1,self%ldim(1)
                        dsq = djsq + (i-c(1))**2
                        if( real(dsq) > radsq ) cycle
                        mm(1) = min(mm(1),self%rmat(i,j,k))
                        mm(2) = max(mm(2),self%rmat(i,j,k))
                    enddo
                enddo
            enddo
        else
            mm(1) = minval(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
            mm(2) = maxval(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
        endif
    end function minmax

    module subroutine loc_sdev( self, winsz, sdevimg, asdev )
        class(image),   intent(in)    :: self
        integer,        intent(in)    :: winsz
        class(image),   intent(inout) :: sdevimg
        real, optional, intent(inout) :: asdev
        real                 :: avg
        integer              :: i, j, k, ir(2), jr(2), kr(2), isz, jsz, ksz, npix
        logical              :: isvol
        isvol = .false.
        isvol = self%is_3d()
        if( isvol )then
            call sdevimg%new(self%ldim, self%smpd)
            !$omp parallel do private(i,j,k,ir,jr,kr,isz,jsz,ksz,npix,avg) default(shared) proc_bind(close)
            do i = 1, self%ldim(1)
                ir(1) = max(1,            i - winsz)
                ir(2) = min(self%ldim(1), i + winsz)
                isz   = ir(2) - ir(1) + 1
                do j = 1, self%ldim(2)
                    jr(1) = max(1,            j - winsz)
                    jr(2) = min(self%ldim(2), j + winsz)
                    jsz   = jr(2) - jr(1) + 1
                    do k = 1, self%ldim(3)
                        kr(1)               = max(1,            k - winsz)
                        kr(2)               = min(self%ldim(3), k + winsz)
                        ksz                 = kr(2) - kr(1) + 1
                        npix                = isz * jsz * ksz
                        avg                 = sum(self%rmat(ir(1):ir(2),jr(1):jr(2),kr(1):kr(2))) / real(npix)
                        sdevimg%rmat(i,j,1) = sqrt(sum((self%rmat(ir(1):ir(2),jr(1):jr(2),kr(1):kr(2)) - avg)**2.0) / real(npix - 1)) 
                    enddo    
                enddo
            enddo
            !$omp end parallel do
            if( present(asdev) )then
                asdev = sum(sdevimg%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))) / real(self%ldim(1) * self%ldim(2) * self%ldim(3))
            endif
        else
            call sdevimg%new(self%ldim, self%smpd)
            !$omp parallel do private(i,j,ir,jr,isz,jsz,npix,avg) default(shared) proc_bind(close)
            do i = 1,self%ldim(1)
               ir(1) = max(1,            i - winsz)
               ir(2) = min(self%ldim(1), i + winsz)
               isz   = ir(2) - ir(1) + 1
               do j = 1,self%ldim(2)
                   jr(1)               = max(1,            j - winsz)
                   jr(2)               = min(self%ldim(2), j + winsz)
                   jsz                 = jr(2) - jr(1) + 1
                   npix                = isz * jsz
                   avg                 = sum(self%rmat(ir(1):ir(2),jr(1):jr(2),1)) / real(npix)
                   sdevimg%rmat(i,j,1) = sqrt(sum((self%rmat(ir(1):ir(2),jr(1):jr(2),1) - avg)**2.0) / real(npix - 1)) 
                enddo
            enddo
            !$omp end parallel do
            if( present(asdev) )then
                asdev = sum(sdevimg%rmat(:self%ldim(1),:self%ldim(2),1)) / real(self%ldim(1) * self%ldim(2))
            endif
        endif
    end subroutine loc_sdev

    module function avg_loc_sdev( self, winsz ) result( asdev )
        class(image), intent(in) :: self
        integer,      intent(in) :: winsz
        real(dp)             :: sum_sdevs
        real                 :: avg, asdev
        integer              :: i, j, k, ir(2), jr(2), kr(2), isz, jsz, ksz, npix
        logical              :: isvol
        isvol     = .false.
        isvol     = self%is_3d()
        sum_sdevs = 0.d0
        if( isvol )then
            do i = 1, self%ldim(1)
                ir(1) = max(1,            i - winsz)
                ir(2) = min(self%ldim(1), i + winsz)
                isz   = ir(2) - ir(1) + 1
                do j = 1, self%ldim(2)
                    jr(1)     = max(1,            j - winsz)
                    jr(2)     = min(self%ldim(2), j + winsz)
                    jsz       = jr(2) - jr(1) + 1
                    do k = 1, self%ldim(3)
                        kr(1) = max(1,            k - winsz)
                        kr(2) = min(self%ldim(3), k + winsz)
                        ksz       = kr(2) - kr(1) + 1
                        npix      = isz * jsz * ksz
                        avg       = sum(self%rmat(ir(1):ir(2),jr(1):jr(2),kr(1):kr(2))) / real(npix)
                        sum_sdevs = sum_sdevs + sqrt(real(sum((self%rmat(ir(1):ir(2),jr(1):jr(2),kr(1):kr(2)) - avg)**2.0),dp) / real(npix-1,dp))
                    enddo
                enddo
            enddo
            asdev = real(sum_sdevs / real(self%ldim(1) * self%ldim(2) * self%ldim(3),dp))
        else
            do i = 1, self%ldim(1)
                ir(1) = max(1,            i - winsz)
                ir(2) = min(self%ldim(1), i + winsz)
                isz   = ir(2) - ir(1) + 1
                do j = 1, self%ldim(2)
                    jr(1)     = max(1,            j - winsz)
                    jr(2)     = min(self%ldim(2), j + winsz)
                    jsz       = jr(2) - jr(1) + 1
                    npix      = isz * jsz
                    avg       = sum(self%rmat(ir(1):ir(2),jr(1):jr(2),1)) / real(npix)
                    sum_sdevs = sum_sdevs + sqrt(real(sum((self%rmat(ir(1):ir(2),jr(1):jr(2),1) - avg)**2.0),dp) / real(npix-1,dp))
                enddo
            enddo
            asdev = real(sum_sdevs / real(self%ldim(1) * self%ldim(2),dp))
        endif
    end function avg_loc_sdev

    module subroutine loc_var( self, varimg, avar )
        class(image),   intent(in)    :: self
        class(image),   intent(inout) :: varimg
        real, optional, intent(inout) :: avar
        real    :: avg, ep, val, var
        integer :: i, j, l, nsz, n_4(3,8)
        if( self%ldim(3) /= 1 ) THROW_HARD('not for 3d')
        call varimg%new(self%ldim, self%smpd)
        do i = 1,self%ldim(1)
            do j = 1,self%ldim(2)
                call neigh_8(self%ldim, [i,j,1], n_4, nsz)
                avg = 0.
                do l = 1,nsz
                    avg = avg + self%rmat(n_4(1,l),n_4(2,l),n_4(3,l))
                end do
                avg = avg / real(nsz)
                ep  = 0.
                var = 0.
                do l = 1,nsz
                    val = self%rmat(n_4(1,l),n_4(2,l),n_4(3,l))
                    ep  = ep  + val
                    var = var + val * val
                end do
                var = (var-ep**2./real(nsz))/(real(nsz)-1.) ! corrected two-pass formula
                varimg%rmat(i,j,1) = var
            end do
        end do
        if( present(avar) )then
            avar = sum(varimg%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))) / real(product(self%ldim))
        endif
    end subroutine loc_var

    module subroutine loc_var3D( self, varimg, avar )
        class(image),   intent(in)    :: self
        class(image),   intent(inout) :: varimg
        real, optional, intent(inout) :: avar
        real    :: avg, ep, val, var
        integer :: i, j, k, l, nsz, n_4(3,6)
        if( self%ldim(3) == 1 ) THROW_HARD('not for 2d')
        call varimg%new(self%ldim, self%smpd)
        !$omp parallel do private(i,j,k,l,n_4,nsz,avg,ep,val,var) default(shared) proc_bind(close)
        do i = 1,self%ldim(1)
            do j = 1,self%ldim(2)
                do k = 1,self%ldim(3)
                    call neigh_4_3D(self%ldim, [i,j,k], n_4, nsz)
                    avg = 0.
                    do l = 1,nsz
                        avg = avg + self%rmat(n_4(1,l),n_4(2,l),n_4(3,l))
                    end do
                    avg = avg / real(nsz)
                    ep  = 0.
                    var = 0.
                    do l = 1,nsz
                        val = self%rmat(n_4(1,l),n_4(2,l),n_4(3,l))
                        ep  = ep  + val
                        var = var + val * val
                    end do
                    var = (var-ep**2./real(nsz))/(real(nsz)-1.) ! corrected two-pass formula
                    varimg%rmat(i,j,k) = var 
                end do
            end do
        end do
        !$omp end parallel do
        if( present(avar) )then
            avar = sum(varimg%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))) / real(product(self%ldim))
        endif
    end subroutine loc_var3D

    module subroutine rmsd( self, dev, mean )
        class(image),   intent(inout) :: self
        real,           intent(out)   :: dev
        real, optional, intent(out)   :: mean
        real :: avg
        if( self%ft )then
            dev = 0.
            if( present(mean) ) mean = 0.
        else
            avg = self%mean()
            if(present(mean)) mean = avg
            dev = sum((self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) - avg)**2.0)&
                  &/ real(product(self%ldim))
            if( dev > 0. )then
                dev = sqrt(dev)
            else
                dev = 0.
            endif
        endif
    end subroutine rmsd

    module subroutine stats_1( self, which, ave, sdev, maxv, minv, msk, med, errout )
        class(image),      intent(inout) :: self
        character(len=*),  intent(in)    :: which
        real,              intent(out)   :: ave, sdev, maxv, minv
        real,    optional, intent(in)    :: msk
        real,    optional, intent(out)   :: med
        logical, optional, intent(out)   :: errout
        integer           :: i, j, k, npix, minlen
        real              :: ci, cj, ck, mskrad, e, var
        logical           :: err, didft, background
        real, allocatable :: pixels(:)
        ! FT
        didft = .false.
        if( self%ft )then
            call self%ifft()
            didft = .true.
        endif
        ! 2d/3d
        if( self%ldim(3) > 1 )then
            minlen = minval(self%ldim)
        else
            minlen = minval(self%ldim(1:2))
        endif
        ! mask
        if( present(msk) )then
            mskrad = msk
        else
            mskrad = real(minlen)/2.
        endif
        ! back/foreground
        if( which.eq.'background' )then
            background = .true.
        else if( which.eq.'foreground' )then
            background = .false.
        else
            THROW_HARD('unrecognized parameter: which; stats_1')
        endif
        allocate( pixels(product(self%ldim)) )
        pixels = 0.
        npix   = 0
        if( self%ldim(3) > 1 )then
            ! 3d
            ci = -real(self%ldim(1))/2.
            do i=1,self%ldim(1)
                cj = -real(self%ldim(2))/2.
                do j=1,self%ldim(2)
                    ck = -real(self%ldim(3))/2.
                    do k=1,self%ldim(3)
                        e = hardedge(ci,cj,ck,mskrad)
                        if( background )then
                            if( e < 0.5 )then
                                npix = npix+1
                                pixels(npix) = self%rmat(i,j,k)
                            endif
                        else
                            if( e > 0.5 )then
                                npix = npix+1
                                pixels(npix) = self%rmat(i,j,k)
                            endif
                        endif
                        ck = ck + 1.
                    end do
                    cj = cj + 1.
                end do
                ci = ci + 1.
            end do
        else
            ! 2d
            ci = -real(self%ldim(1))/2.
            do i=1,self%ldim(1)
                cj = -real(self%ldim(2))/2.
                do j=1,self%ldim(2)
                    e = hardedge(ci,cj,mskrad)
                    if( background )then
                        if( e < 0.5 )then
                            npix = npix+1
                            pixels(npix) = self%rmat(i,j,1)
                        endif
                    else
                        if( e > 0.5 )then
                            npix = npix+1
                            pixels(npix) = self%rmat(i,j,1)
                        endif
                    endif
                    cj = cj + 1.
                end do
                ci = ci + 1.
            end do
        endif
        maxv = maxval(pixels(:npix))
        minv = minval(pixels(:npix))
        if(npix>1)then
            call moment( pixels(:npix), ave, sdev, var, err )
            if( present(med) ) med  = median_nocopy(pixels(:npix))
        end if
        deallocate( pixels )
        if( present(errout) )then
            errout = err
        else
            if( err ) THROW_WARN('variance zero; stats_1')
        endif
        if( didft ) call self%fft()
    end subroutine stats_1

    module subroutine stats_2( self, ave, sdev, maxv, minv, mskimg, med, errout )
        class(image),           intent(inout) :: self
        real,                   intent(out)   :: ave, sdev, maxv, minv
        class(image), optional, intent(in)    :: mskimg
        real,         optional, intent(out)   :: med
        logical,      optional, intent(out)   :: errout
        real              :: var
        logical           :: err
        real, allocatable :: pixels(:)
        ! FT
        if( self%ft ) THROW_HARD('not for FTed imgs; stats_2')
        if( present(mskimg) )then
            pixels = pack(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)),&
                &mskimg%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) > 0.95 )
        else
          pixels = pack(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), .true.)
        endif
        maxv = maxval(pixels)
        minv = minval(pixels)
        call moment( pixels, ave, sdev, var, err )
        if( present(med) ) med  = median_nocopy(pixels)
        deallocate( pixels )
        if( present(errout) )then
            errout = err
        else
            if( err ) THROW_WARN('variance zero; stats_2')
        endif
    end subroutine stats_2

    module function variance( self ) result( var )
        class(image), intent(in) :: self
        real    :: ave, var, ep, rmat_subtr_avg(self%ldim(1),self%ldim(2),self%ldim(3))
        integer :: npix
        npix           = product(self%ldim)
        ave            = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))) / real(npix)
        rmat_subtr_avg = self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) - ave
        ep             = sum(rmat_subtr_avg)
        var            = sum(rmat_subtr_avg**2.0)
        var            = (var-ep**2./real(npix))/(real(npix)-1.) ! corrected two-pass formula
    end function variance

    module real function skew( self, mask )
        class(image),      intent(in)  :: self
        logical, optional, intent(in)  :: mask(self%ldim(1),self%ldim(2))
        if( self%is_3d() ) THROW_HARD('2D only!')
        if( self%is_ft() ) THROW_HARD('Real space only!')
        skew = skewness(self%rmat(1:self%ldim(1),1:self%ldim(2),1), mask)
    end function skew

    module real function kurt( self, mask )
        class(image),      intent(in)  :: self
        logical, optional, intent(in)  :: mask(self%ldim(1),self%ldim(2))
        if( self%is_3d() ) THROW_HARD('2D only!')
        if( self%is_ft() ) THROW_HARD('Real space only!')
        kurt = kurtosis(self%rmat(1:self%ldim(1),1:self%ldim(2),1), mask)
    end function kurt

    module function noisesdev( self, msk ) result( sdev )
        use simple_online_var, only: online_var
        class(image), intent(inout) :: self
        real,         intent(in)    :: msk
        type(online_var)            :: ovar
        integer                     :: i, j, k
        real                        :: ci, cj, ck, e, sdev, mv(2)
        logical                     :: didft
        ovar = online_var( )
        didft = .false.
        if( self%ft )then
            call self%ifft()
            didft = .true.
        endif
        ci = -real(self%ldim(1))/2.
        do i=1,self%ldim(1)
            cj = -real(self%ldim(2))/2.
            do j=1,self%ldim(2)
                ck = -real(self%ldim(3))/2.
                do k=1,self%ldim(3)
                    if( self%ldim(3) > 1 )then
                        e = hardedge(ci,cj,ck,msk)
                    else
                        e = hardedge(ci,cj,msk)
                    endif
                    if( e < 0.5 )then
                        call ovar%add(self%rmat(i,j,k))
                    endif
                    ck = ck+1
                end do
                cj = cj+1.
            end do
            ci = ci+1.
        end do
        mv(1) = ovar%get_mean()
        mv(2) = ovar%get_var()
        sdev = 0.
        if( mv(2) > 0. ) sdev = sqrt(mv(2))
        if( didft ) call self%fft()
    end function noisesdev

    module function mean( self ) result( avg )
        class(image), intent(inout) :: self
        real    :: avg
        logical :: didft
        didft = .false.
        if( self%ft )then
            call self%ifft()
            didft = .true.
        endif
        avg = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))/real(product(self%ldim))
        if( didft ) call self%ifft()
    end function mean

    module logical function contains_nans( self )
        class(image), intent(in) :: self
        integer :: i, j, k
        contains_nans = .false.
        do i=1,size(self%rmat,1)
            do j=1,size(self%rmat,2)
                do k=1,size(self%rmat,3)
                    if( .not. is_a_number(self%rmat(i,j,k)) )then
                        contains_nans = .true.
                        return
                    endif
                end do
            end do
        end do
    end function contains_nans

    module subroutine checkimg4nans( self )
        class(image), intent(in) :: self
        if( self%ft )then
            call check4nans3D(self%cmat)
        else
            call check4nans3D(self%rmat)
        endif
    end subroutine checkimg4nans

    module subroutine cure( self, maxv, minv, ave, sdev, n_nans )
        class(image), intent(inout) :: self
        real,         intent(out)   :: maxv, minv, ave, sdev
        integer,      intent(out)   :: n_nans
        integer                     :: i, j, k, npix
        real                        :: var, ep, dev
        if( self%ft )then
            THROW_WARN('cannot cure FTs; cure')
            return
        endif
        npix   = product(self%ldim)
        n_nans = 0
        ave    = 0.
        !$omp parallel do default(shared) private(i,j,k) schedule(static)&
        !$omp collapse(3) proc_bind(close) reduction(+:n_nans,ave)
        do i=1,self%ldim(1)
            do j=1,self%ldim(2)
                do k=1,self%ldim(3)
                    if( .not. is_a_number(self%rmat(i,j,k)) )then
                        n_nans = n_nans + 1
                    else
                        ave = ave + self%rmat(i,j,k)
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        if( n_nans > 0 )then
            write(logfhandle,*) 'found NaNs in simple_image; cure:', n_nans
        endif
        ave       = ave/real(npix)
        maxv      = maxval( self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) )
        minv      = minval( self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) )
        self%rmat = self%rmat - ave
        ! calc sum of devs and sum of devs squared
        ep = 0.
        var = 0.
        !$omp parallel do default(shared) private(i,j,k,dev) schedule(static)&
        !$omp collapse(3) proc_bind(close) reduction(+:ep,var)
        do i=1,self%ldim(1)
            do j=1,self%ldim(2)
                do k=1,self%ldim(3)
                    dev = self%rmat(i,j,k)
                    ep  = ep + dev
                    var = var + dev * dev
                end do
            end do
        end do
        !$omp end parallel do
        var  = (var-ep**2./real(npix))/(real(npix)-1.) ! corrected two-pass formula
        sdev = sqrt(var)
        if( sdev > 0. ) self%rmat = self%rmat/sdev
    end subroutine cure

    module function dead_hot_positions( self, frac ) result( pos )
        class(image), intent(in) :: self
        real, intent(in)         :: frac
        logical, allocatable     :: pos(:,:)
        integer :: ipix, jpix, cnt
        allocate(pos(self%ldim(1),self%ldim(2)))
        pos = .false.
        cnt = 0
        do ipix=1,self%ldim(1)
            do jpix=1,self%ldim(2)
                if( ran3() <= frac )then
                    pos(ipix,jpix) = .true.
                    cnt = cnt+1
                endif
            end do
        end do
    end function dead_hot_positions

    !===========================
    ! Gradients / geometry
    !===========================

    ! This function returns a the gradient matrix of the input image.
    ! It is also possible to have derivates row and column
    ! as output (optional).
    ! It uses masks found in http://www.holoborodko.com/pavel/image-processing/edge-detection/
    ! which is better than Sobel masks because of:
    !                     1) isotropic noise suppression
    !                     2) the estimation of the gradient is still precise
    module subroutine calc_gradient(self, grad, Dc, Dr)
        class(image),   intent(inout) :: self
        real,           intent(out)   :: grad(self%ldim(1), self%ldim(2), self%ldim(3)) ! gradient matrix
        real, optional, intent(out)   :: Dc(self%ldim(1), self%ldim(2), self%ldim(3)), Dr(self%ldim(1), self%ldim(2), self%ldim(3)) ! derivates column and row matrices
        type(image)        :: img_p                         ! padded image
        real, allocatable  :: wc(:,:), wr(:,:)              ! row and column Sobel masks
        integer, parameter :: L1 = 5 , L2 = 3               ! dimension of the masks
        integer            :: ldim(3)                       ! dimension of the image, save just for comfort
        integer            :: i,j,m,n                       ! loop indeces
        real :: Ddc(self%ldim(1),self%ldim(2),self%ldim(3))
        real :: Ddr(self%ldim(1),self%ldim(2),self%ldim(3)) ! column and row derivates
        ldim = self%ldim
        allocate(wc((-(L1-1)/2):((L1-1)/2),(-(L2-1)/2):((L2-1)/2)),&
        &        wr(-(L2-1)/2:(L2-1)/2,-(L1-1)/2:(L1-1)/2), source = 0.)
        wc = (1./32.)*reshape([-1,-2,0,2,1,-2,-4,0,4,2,-1,-2,0,2,1], [L1,L2])
        wr = (1./32.)*reshape([-1,-2,-1,-2,-4,-2,0,0,0,2,4,2,1,2,1], [L2,L1])
        Ddc  = 0. ! initialisation
        Ddr  = 0.
        grad = 0.
        call img_p%new([ldim(1)+L1-1,ldim(2)+L1-1,1],1.) ! pad with the biggest among L1 and L2
        call self%pad(img_p)
        do i = 1, ldim(1)
          do j = 1, ldim(2)
              do m = -(L1-1)/2,(L1-1)/2
                  do n = -(L2-1)/2,(L2-1)/2
                      Ddc(i,j,1) = Ddc(i,j,1)+img_p%rmat(i+m+2,j+n+2,1)*wc(m,n)
                  end do
              end do
          end do
        end do
        do i = 1, ldim(1)
          do j = 1, ldim(2)
              do m = -(L2-1)/2,(L2-1)/2
                  do n = -(L1-1)/2,(L1-1)/2
                      Ddr(i,j,1) = Ddr(i,j,1)+img_p%rmat(i+m+2,j+n+2,1)*wr(m,n)
                  end do
              end do
          end do
        end do
        deallocate(wc, wr)
        grad = sqrt(Ddc**2 + Ddr**2)
        if(present(Dc)) Dc = Ddc
        if(present(Dr)) Dr = Ddr
        call img_p%kill
    end subroutine calc_gradient

    module subroutine gradients_magnitude( self, self_out )
        class(image), intent(in)    :: self
        class(image), intent(inout) :: self_out
        integer :: i,j,ni,nj
        if( self%is_ft() ) THROW_HARD('Image input must be in the spatial domain!')
        if( .not.self%is_2d() ) THROW_HARD('Image input must be in 2D!')
        if( .not.self_out%exists() ) call self_out%copy(self)
        ni = self%ldim(1)
        nj = self%ldim(2)
        !$omp parallel private(i,j) proc_bind(close) default(shared)
        !$omp do
        do i = 1,ni
            if( i == 1 )then
                self_out%rmat(i,:nj,1) = (self%rmat(2,:nj,1) - self%rmat(1,:nj,1))**2
            else if( i == ni )then
                self_out%rmat(i,:nj,1) = (self%rmat(i,:nj,1) - self%rmat(i-1,:nj,1))**2
            else
                self_out%rmat(i,:nj,1) = (self%rmat(i+1,:nj,1) - self%rmat(i-1,:nj,1))**2
            endif
        enddo
        !$omp end do
        !$omp do
        do j = 1,nj
            if( j == 1 )then
                self_out%rmat(:ni,j,1) = self_out%rmat(:ni,j,1) + (self%rmat(:ni,2,1) - self%rmat(:ni,1,1))**2
            else if( j == nj )then
                self_out%rmat(:ni,j,1) = self_out%rmat(:ni,j,1) + (self%rmat(:ni,j,1) - self%rmat(:ni,j-1,1))**2
            else
                self_out%rmat(:ni,j,1) = self_out%rmat(:ni,j,1) + (self%rmat(:ni,j+1,1) - self%rmat(:ni,j-1,1))**2
            endif
        enddo
        !$omp end do
        !$omp workshare
        self_out%rmat(1:ni,1:nj,1) = self_out%rmat(1:ni,1:nj,1)/4.
        self_out%rmat(1:ni,1:nj,1) = merge(sqrt(self_out%rmat(1:ni,1:nj,1)), 0., self_out%rmat(1:ni,1:nj,1)>0.)
        !$omp end workshare
        !$omp end parallel
    end subroutine gradients_magnitude

    ! This function returns the derivates row, column, and z as outputs,
    ! together with the optional gradient matrix/volume.
    ! It uses standard central difference scheme
    module subroutine gradient(self, Dc, Dr, Dz, grad)
        class(image),   intent(inout) :: self
        real, optional, intent(out)   :: Dc(self%ldim(1), self%ldim(2), self%ldim(3)), & ! derivates column matrix
                                         Dr(self%ldim(1), self%ldim(2), self%ldim(3)), & ! derivates row matrix
                                         Dz(self%ldim(1), self%ldim(2), self%ldim(3))    ! derivates z matrix
        real, optional, intent(out)   :: grad(self%ldim(1), self%ldim(2), self%ldim(3))  ! gradient matrix
        integer :: ldim(3), k
        real    :: Ddc(self%ldim(1),self%ldim(2),self%ldim(3))
        real    :: Ddr(self%ldim(1),self%ldim(2),self%ldim(3))
        real    :: Ddz(self%ldim(1),self%ldim(2),self%ldim(3))
        ldim = self%ldim
        Ddc  = 0. ! initialisation
        Ddr  = 0.
        Ddz  = 0.
        do k = 2, ldim(1)-1
            Ddc(k,1:ldim(2),1:ldim(3)) = 0.5*(self%rmat(k+1,1:ldim(2),1:ldim(3)) - self%rmat(k-1,1:ldim(2),1:ldim(3)))
        enddo
        do k = 2, ldim(2)-1
            Ddr(1:ldim(1),k,1:ldim(3)) = 0.5*(self%rmat(1:ldim(1),k+1,1:ldim(3)) - self%rmat(1:ldim(1),k-1,1:ldim(3)))
        enddo
        Ddc(1        ,1:ldim(2),1:ldim(3)) = self%rmat(2        ,1:ldim(2),1:ldim(3)) - self%rmat(1          ,1:ldim(2)  ,1:ldim(3))
        Ddc(  ldim(1),1:ldim(2),1:ldim(3)) = self%rmat(  ldim(1),1:ldim(2),1:ldim(3)) - self%rmat(  ldim(1)-1,1:ldim(2)  ,1:ldim(3))
        Ddr(1:ldim(1),1        ,1:ldim(3)) = self%rmat(1:ldim(1),2        ,1:ldim(3)) - self%rmat(1:ldim(1)  ,1          ,1:ldim(3))
        Ddr(1:ldim(1),  ldim(2),1:ldim(3)) = self%rmat(1:ldim(1),  ldim(2),1:ldim(3)) - self%rmat(1:ldim(1)  ,  ldim(2)-1,1:ldim(3))
        if( ldim(3) > 1 )then
            do k = 2, ldim(3)-1
                Ddz(1:ldim(1),1:ldim(2),k) = 0.5*(self%rmat(1:ldim(1),1:ldim(2),k+1) - self%rmat(1:ldim(1),1:ldim(2),k-1))
            enddo
            Ddz(1:ldim(1),1:ldim(2),1      ) = self%rmat(1:ldim(1),1:ldim(2),2)       - self%rmat(1:ldim(1),1:ldim(2),1)
            Ddz(1:ldim(1),1:ldim(2),ldim(3)) = self%rmat(1:ldim(1),1:ldim(2),ldim(3)) - self%rmat(1:ldim(1),1:ldim(2),ldim(3)-1)
        endif
        if(present(Dc))   Dc   = Ddc
        if(present(Dr))   Dr   = Ddr
        if(present(Dz))   Dz   = Ddz
        if(present(grad)) grad = sqrt(Ddc**2 + Ddr**2 + Ddz**2)
    end subroutine gradient

    module subroutine calc_ice_score( self, score )
        class(image), intent(in)  :: self
        real,         intent(out) :: score
        real,   parameter :: START_FREQ = 15.
        real,   parameter :: END_FREQ   = 6.
        real, allocatable :: res(:), tmp(:)
        real    :: powspec(fdim(self%ldim(1)) - 1)
        real    :: g, gs, ge, mag, mag_max, band_avg, ice_avg
        integer :: lims(3,2), ice_maxind, start_find, end_find
        integer :: nbands, s, e, h, k, hmax, kmax, sh, cnt
        score = 0.
        if( self%smpd > (ICE_BAND1/2.) ) return
        if( .not.self%is_ft() ) THROW_HARD('Image input must be in the Fourier domain!; calc_ice_score')
        lims = self%loop_lims(2)
        res  = get_resarr(self%ldim(1), self%smpd)
        ice_maxind = get_find_at_res(res, ICE_BAND1)
        start_find = get_find_at_res(res, START_FREQ)
        end_find   = get_find_at_res(res, END_FREQ)
        call self%power_spectrum(powspec)
        nbands = end_find-start_find+1
        tmp = powspec(start_find:end_find)
        call hpsort(tmp)
        e = nbands
        s = nint(0.5 *real(nbands))
        band_avg = sum(tmp(s:e)) / real(e-s+1)
        ! location of maximum in ice band
        mag_max = -1.
        gs = real(max(            1,   ice_maxind-3)) / real(self%ldim(1))
        ge = real(min(size(powspec)-1, ice_maxind+3)) / real(self%ldim(1))
        do k = lims(2,1),lims(2,2)
            do h = lims(1,1),lims(1,2)
                sh = nint(hyp(h,k))
                g  = real(sh) / real(self%ldim(1))
                if( g > gs .and. g < ge )then
                    mag = csq_fast(self%get_fcomp2D(h,k))
                    if( mag > mag_max )then
                        hmax = h
                        kmax = k
                        mag_max = mag
                    endif
                endif
            end do
        end do
        ! ice peak
        ice_avg = 0.
        cnt     = 0
        do k = kmax-1,kmax+1
            do h = hmax-1,hmax+1
                sh = nint(hyp(h,k))
                g  = real(sh) / real(self%ldim(1))
                if( g < 0.5 )then
                    ice_avg = ice_avg + csq_fast(self%get_fcomp2D(h,k))
                    cnt     = cnt+1
                endif
            enddo
        enddo
        ice_avg = ice_avg / real(cnt)
        score   = ice_avg / (band_avg + TINY)
    end subroutine calc_ice_score

    module subroutine calc_principal_axes_rotmat( self, radius, R )
        class(image), intent(in)  :: self
        real,         intent(in)  :: radius
        real,         intent(out) :: R(3,3)
        real(dp) :: coord(3), ixx, iyy, izz, ixz, ixy, iyz, m
        real(dp) :: inertia(3,3), eigvals(3), eigvecs(3,3)
        real     :: radiussq
        integer  :: icenter(3),i,j,k
        if( self%is_ft() )      THROW_HARD('Real space only; calc_principal_axes_rotmat')
        if( .not.self%is_3d() ) THROW_HARD('Volumes only; calc_principal_axes_rotmat')
        icenter  = nint(real(self%ldim)/2.)+1
        radiussq = radius**2
        ! Inertia Tensor
        ixx = 0.d0; iyy = 0.d0; izz = 0.d0
        ixy = 0.d0; ixz = 0.d0; iyz = 0.d0
        do k =1,self%ldim(3)
        do j =1,self%ldim(2)
        do i =1,self%ldim(1)
            if( (real(sum(([i,j,k]-icenter)**2)) < radiussq) .and. (self%rmat(i,j,k)>0.0) )then
                coord = real([i,j,k]-icenter, dp)
                m     = real(self%rmat(i,j,k), dp)
                ixx   = ixx + m * (coord(2)**2 + coord(3)**2)
                iyy   = iyy + m * (coord(1)**2 + coord(3)**2)
                izz   = izz + m * (coord(1)**2 + coord(2)**2)
                ixy   = ixy + m * coord(1) * coord(2)
                ixz   = ixz + m * coord(1) * coord(3)
                iyz   = iyz + m * coord(2) * coord(3)
            endif
        enddo
        enddo
        enddo
        inertia(1,:) = [ ixx, -ixy, -ixz]
        inertia(2,:) = [-ixy,  iyy, -iyz]
        inertia(3,:) = [-ixz, -iyz,  izz]
        ! Spectral analysis
        call svdcmp(inertia, eigvals, eigvecs)
        call eigsrt(eigvals, eigvecs, 3, 3)
        ! double checking
        ! identity = matmul(eigvecs, transpose(eigvecs))
        ! inertia  = matmul(eigvecs, matmul(eye(3)*eigvals, transpose(eigvecs)))
        ! Reverse rotation matrix
        R = real(transpose(eigvecs))
    end subroutine calc_principal_axes_rotmat

    !===========================
    ! Physical coords helpers
    !===========================

    module function loop_lims( self, mode, lp_dyn ) result( lims )
        class(image), intent(in)   :: self
        integer, intent(in)        :: mode
        real, intent(in), optional :: lp_dyn
        integer                    :: lims(3,2)
        if( present(lp_dyn) )then
            lims = self%fit%loop_lims(mode, lp_dyn)
        else
            lims = self%fit%loop_lims(mode)
        endif
    end function loop_lims

    module pure function comp_addr_phys1(self,logi) result(phys)
        class(image), intent(in) :: self
        integer,      intent(in) :: logi(3) !<  Logical address
        integer                  :: phys(3) !<  Physical address
        phys = self%fit%comp_addr_phys(logi)
    end function comp_addr_phys1

    module pure function comp_addr_phys2(self,h,k,m) result(phys)
        class(image), intent(in) :: self
        integer,      intent(in) :: h,k,m !<  Logical address
        integer                  :: phys(3) !<  Physical address
        phys = self%fit%comp_addr_phys(h,k,m)
    end function comp_addr_phys2

    module pure function comp_addr_phys3(self,h,k) result(phys)
        class(image), intent(in) :: self
        integer,      intent(in) :: h,k     !<  Logical address
        integer                  :: phys(2) !<  Physical address
        phys = self%fit%comp_addr_phys(h,k)
    end function comp_addr_phys3

    !===========================
    ! Correlation / distances
    !===========================

    module function corr( self1, self2, lp_dyn, hp_dyn ) result( r )
        class(image),   intent(inout) :: self1, self2
        real, optional, intent(in)    :: lp_dyn, hp_dyn
        real    :: r, sumasq, sumbsq, eps
        integer :: h, k, l, phys(3), lims(3,2), sqarg, sqlp, sqhp
        logical :: didft1, didft2
        r = 0.
        sumasq = 0.
        sumbsq = 0.
        eps    = epsilon(sumasq)
        if( self1.eqdims.self2 )then
            didft1 = .false.
            if( .not. self1%ft )then
                call self1%fft()
                didft1 = .true.
            endif
            didft2 = .false.
            if( .not. self2%ft )then
                call self2%fft()
                didft2 = .true.
            endif
            if( present(lp_dyn) )then
                lims = self1%fit%loop_lims(1,lp_dyn)
            else
                lims = self1%fit%loop_lims(2) ! Nyqvist default low-pass limit
            endif
            sqlp = (maxval(lims(:,2)))**2
            if( present(hp_dyn) )then
                sqhp = max(2,self1%get_find(hp_dyn))**2
            else
                sqhp = 2 ! index 2 default high-pass limit
            endif
            !$omp parallel do collapse(3) default(shared) private(h,k,l,sqarg,phys)&
            !$omp reduction(+:r,sumasq,sumbsq) schedule(static) proc_bind(close)
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        sqarg = h*h + k*k + l*l
                        if( sqarg <= sqlp .and. sqarg >= sqhp  )then
                            phys = self1%fit%comp_addr_phys([h,k,l])
                            ! real part of the complex mult btw 1 and 2*
                            r = r + real(self1%cmat(phys(1),phys(2),phys(3))*conjg(self2%cmat(phys(1),phys(2),phys(3))))
                            sumasq = sumasq + csq(self2%cmat(phys(1),phys(2),phys(3)))
                            sumbsq = sumbsq + csq(self1%cmat(phys(1),phys(2),phys(3)))
                        endif
                    end do
                end do
            end do
            !$omp end parallel do
            if( r < eps .and. sumasq < eps .and. sumbsq < eps )then
                r = 1.
            elseif( sqrt(sumasq * sumbsq) < eps )then
                r = 0.
            else
                r = r / sqrt(sumasq * sumbsq)
            endif
            if( didft1 ) call self1%ifft()
            if( didft2 ) call self2%ifft()
        else
            write(logfhandle,*) 'self1%ldim:', self1%ldim
            write(logfhandle,*) 'self2%ldim:', self2%ldim
            THROW_HARD('images to be correlated need to have same dimensions; corr')
        endif
    end function corr

    module function corr_shifted( self_ref, self_ptcl, shvec, lp_dyn, hp_dyn ) result( r )
        class(image),   intent(inout) :: self_ref, self_ptcl
        real,           intent(in)    :: shvec(3)
        real, optional, intent(in)    :: lp_dyn, hp_dyn
        real                          :: r, sumasq, sumbsq
        complex                       :: shcomp
        integer                       :: h, k, l, phys(3), lims(3,2), sqarg, sqlp, sqhp
        ! this is for highly optimised code, so we assume that images are always Fourier transformed beforehand
        if( .not. self_ref%ft  ) THROW_HARD('self_ref not FTed;  corr_shifted')
        if( .not. self_ptcl%ft ) THROW_HARD('self_ptcl not FTed; corr_shifted')
        r = 0.
        sumasq = 0.
        sumbsq = 0.
        if( present(lp_dyn) )then
            lims = self_ref%fit%loop_lims(1,lp_dyn)
        else
            lims = self_ref%fit%loop_lims(2) ! Nyqvist default low-pass limit
        endif
        sqlp = (maxval(lims(:,2)))**2
        if( present(hp_dyn) )then
            sqhp = max(2,self_ref%get_find(hp_dyn))**2
        else
            sqhp = 2 ! index 2 default high-pass limit
        endif
        !$omp parallel do collapse(3) default(shared) private(h,k,l,sqarg,phys,shcomp)&
        !$omp reduction(+:r,sumasq,sumbsq) schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    sqarg = h*h + k*k + l*l
                    if( sqarg <= sqlp .and. sqarg >= sqhp  )then
                        phys = self_ref%fit%comp_addr_phys(h,k,l)
                        ! shift particle
                        shcomp = self_ptcl%cmat(phys(1),phys(2),phys(3))*&
                            &self_ptcl%oshift([h,k,l], shvec)
                        ! real part of the complex mult btw 1 and 2*
                        r = r + real(self_ref%cmat(phys(1),phys(2),phys(3))*conjg(shcomp))
                        sumasq = sumasq + csq(shcomp)
                        sumbsq = sumbsq + csq(self_ref%cmat(phys(1),phys(2),phys(3)))
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        if( sumasq > 0. .and. sumbsq > 0. )then
            r = r / sqrt(sumasq * sumbsq)
        else
            r = 0.
        endif
    end function corr_shifted

    module function real_corr_1( self1, self2 ) result( r )
        class(image), intent(inout) :: self1, self2
        real    :: diff1(self1%ldim(1),self1%ldim(2),self1%ldim(3)), diff1sc
        real    :: diff2(self2%ldim(1),self2%ldim(2),self2%ldim(3)), diff2sc
        real    :: r, ax, ay, sxx, syy, sxy, npix
        integer :: i,j,k
        if( self1%wthreads .and. self2%wthreads )then
            npix = real(product(self1%ldim))
            ax   = self1%mean()
            ay   = self2%mean()
            sxx = 0.
            syy = 0.
            sxy = 0.
            !$omp parallel do default(shared) private(i,j,k,diff1sc,diff2sc) collapse(3) proc_bind(close) schedule(static) reduction(+:sxx,syy,sxy)
            do i=1,self1%ldim(1)
                do j=1,self1%ldim(2)
                    do k=1,self1%ldim(3)
                        diff1sc = self1%rmat(i,j,k) -ax
                        diff2sc = self2%rmat(i,j,k) -ay
                        sxx     = sxx + diff1sc * diff1sc
                        syy     = syy + diff2sc * diff2sc
                        sxy     = sxy + diff1sc * diff2sc
                    end do
                end do
            end do
            !$omp end parallel do
            if( sxx > TINY .and. syy > TINY )then
                r = sxy / sqrt(sxx * syy)
            else
                r = 0.
            endif
        else
            diff1 = 0.
            diff2 = 0.
            npix  = real(product(self1%ldim))
            ax    = sum(self1%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3))) / npix
            ay    = sum(self2%rmat(:self2%ldim(1),:self2%ldim(2),:self2%ldim(3))) / npix
            diff1 = self1%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)) - ax
            diff2 = self2%rmat(:self2%ldim(1),:self2%ldim(2),:self2%ldim(3)) - ay
            sxx   = sum(diff1 * diff1)
            syy   = sum(diff2 * diff2)
            sxy   = sum(diff1 * diff2)
            if( sxx > TINY .and. syy > TINY )then
                r = sxy / sqrt(sxx * syy)
            else
                r = 0.
            endif
        endif
    end function real_corr_1

    module function real_corr_2( self1, self2, mask ) result( r )
        class(image), intent(inout) :: self1, self2
        logical,      intent(in)    :: mask(self1%ldim(1),self1%ldim(2),self1%ldim(3))
        real :: diff1(self1%ldim(1),self1%ldim(2),self1%ldim(3))
        real :: diff2(self2%ldim(1),self2%ldim(2),self2%ldim(3))
        real :: r, sxx, syy, sxy, npix, ax, ay
        diff1 = 0.
        diff2 = 0.
        npix  = real(count(mask))
        ax    = sum(self1%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)), mask=mask) / npix
        ay    = sum(self2%rmat(:self2%ldim(1),:self2%ldim(2),:self2%ldim(3)), mask=mask) / npix
        diff1 = self1%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)) - ax
        diff2 = self2%rmat(:self2%ldim(1),:self2%ldim(2),:self2%ldim(3)) - ay
        sxx   = sum(diff1 * diff1, mask=mask)
        syy   = sum(diff2 * diff2, mask=mask)
        sxy   = sum(diff1 * diff2, mask=mask)
        if( sxx > TINY .and. syy > TINY )then
            r = sxy / sqrt(sxx * syy)
        else
            r = 0.
        endif
    end function real_corr_2

    module function euclid_dist_two_imgs(self1, self2, mask1) result(dist)
        use simple_linalg, only: euclid
        class(image),      intent(inout) :: self1, self2
        logical, optional, intent(in)    :: mask1(self1%ldim(1),self1%ldim(2),self1%ldim(3))
        real              :: diff1(self1%ldim(1),self1%ldim(2),self1%ldim(3))
        real              :: diff2(self2%ldim(1),self2%ldim(2),self2%ldim(3))
        real              :: ax, ay, npix, dist
        real, allocatable :: diff1_flat(:), diff2_flat(:)
        logical           :: mask_here(self1%ldim(1),self1%ldim(2),self1%ldim(3))
        if (present(mask1)) then
            mask_here = mask1
        else 
            mask_here = .true.
        end if
        npix       = real(count(mask_here))
        ax         = sum(self1%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)), mask=mask_here) / npix
        ay         = sum(self2%rmat(:self2%ldim(1),:self2%ldim(2),:self2%ldim(3)), mask=mask_here) / npix
        diff1      = self1%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)) - ax
        diff2      = self2%rmat(:self2%ldim(1),:self2%ldim(2),:self2%ldim(3)) - ay
        diff1_flat = pack(diff1, mask=.true.)
        diff2_flat = pack(diff2, mask=.true.)
        dist       = euclid(diff1_flat, diff2_flat)
    end function euclid_dist_two_imgs

    module subroutine phase_corr(self1, self2, pc, lp )
        class(image),      intent(inout) :: self1, self2, pc
        real,              intent(in)    :: lp
        real, parameter :: width = 3.
        complex     :: c1,c2
        real        :: w,rw,rlplim,rsh,normsq
        real(dp)    :: sqsum1,sqsum2
        integer     :: nrflims(3,2),phys(3)
        integer     :: h,k,l,shsq, lplim, lplimsq, bplplimsq
        if( .not. all([self1%is_ft(),self2%is_ft(),pc%is_ft()]) )then
            THROW_HARD('All inputted images must be FTed')
        endif
        sqsum1    = 0.d0
        sqsum2    = 0.d0
        nrflims   = self1%loop_lims(2)
        lplim     = calc_fourier_index(lp,minval(self1%ldim(1:2)),self1%smpd)
        rlplim    = real(lplim)
        lplimsq   = lplim*lplim
        bplplimsq = min(minval(nrflims(1:2,2)),lplim-nint(WIDTH))**2
        call pc%zero_and_flag_ft
        if( self1%wthreads .or. self2%wthreads )then
            !$omp parallel do default(shared) private(h,k,l,w,rw,shsq,phys,rsh,c1,c2)&
            !$omp proc_bind(close) schedule(static) reduction(+:sqsum1,sqsum2)
            do h = nrflims(1,1),nrflims(1,2)
                rw = merge(1., 2., h==0)
                do k = nrflims(2,1),nrflims(2,2)
                    do l = nrflims(3,1),nrflims(3,2)
                        shsq = h*h+k*k+l*l
                        w    = 1.
                        if( shsq > lplimsq )then
                            cycle
                        else if( shsq == 0 )then
                            cycle
                        else if( shsq > bplplimsq )then
                            rsh = sqrt(real(shsq))
                            w   = 0.5*(1.+cos(PI*(rsh-(rlplim-width))/width))
                        endif
                        phys = pc%comp_addr_phys(h,k,l)
                        c1 = w*self1%cmat(phys(1),phys(2),phys(3))
                        c2 = w*self2%cmat(phys(1),phys(2),phys(3))
                        sqsum1 = sqsum1 + real(rw*csq(c1),dp)
                        sqsum2 = sqsum2 + real(rw*csq(c2),dp)
                        pc%cmat(phys(1),phys(2),phys(3)) = c1 * conjg(c2)
                    enddo
                enddo
            enddo
            !$omp end parallel do
        else
            do h = nrflims(1,1),nrflims(1,2)
                rw = merge(1., 2., h==0)
                do k = nrflims(2,1),nrflims(2,2)
                    do l = nrflims(3,1),nrflims(3,2)
                        shsq = h*h+k*k+l*l
                        w    = 1.
                        if( shsq > lplimsq )then
                            cycle
                        else if( shsq == 0 )then
                            cycle
                        else if( shsq > bplplimsq )then
                            rsh = sqrt(real(shsq))
                            w   = 0.5*(1.+cos(PI*(rsh-(rlplim-width))/width))
                        endif
                        phys = pc%comp_addr_phys(h,k,l)
                        c1 = w*self1%cmat(phys(1),phys(2),phys(3))
                        c2 = w*self2%cmat(phys(1),phys(2),phys(3))
                        sqsum1 = sqsum1 + real(rw*csq(c1),dp)
                        sqsum2 = sqsum2 + real(rw*csq(c2),dp)
                        pc%cmat(phys(1),phys(2),phys(3)) = c1 * conjg(c2)
                    enddo
                enddo
            enddo
        endif
        call pc%ifft()
        normsq = real(sqsum1*sqsum2)
        if( is_a_number(normsq) )then
            if( normsq > 1.0e-12 ) pc%rmat = pc%rmat / sqrt(normsq)
        endif
    end subroutine phase_corr

    module subroutine fcorr_shift( self1, self2, trs, shift, peak_interp )
        class(image),      intent(inout) :: self1, self2
        real,              intent(in)    :: trs
        real,              intent(inout) :: shift(2)
        logical, optional, intent(in)    :: peak_interp
        real    :: alpha, beta, gamma, denom
        integer :: center(2), pos(2), itrs
        logical :: l_interp
        if( self1%is_3d() .or. self2%is_3d() ) THROW_HARD('2d only supported')
        if( .not.(self1%is_ft() .and. self2%is_ft()) ) THROW_HARD('FTed only supported')
        if( .not.(self1.eqdims.self2) ) THROW_HARD('Inconsistent dimensions in fcorr_shift')
        l_interp = .false.
        if(present(peak_interp)) l_interp = peak_interp
        ! dimensions
        center = self1%ldim(1:2)/2+1
        itrs   = min(floor(trs),minval(center)-1)
        ! Correlation image
        self1%cmat = self1%cmat * conjg(self2%cmat)
        call self1%ifft
        ! maximum correlation & offset
        pos   = maxloc(self1%rmat(center(1)-itrs:center(1)+itrs, center(2)-itrs:center(2)+itrs, 1)) -itrs-1
        ! peak interpolation
        if( l_interp )then
            shift = real(pos)
            beta  = self1%rmat(center(1)+pos(1),center(2)+pos(2), 1)
            ! along x
            if( abs(pos(1)) < itrs )then ! within limits
                alpha = self1%rmat(center(1)+pos(1)-1,center(2)+pos(2), 1)
                gamma = self1%rmat(center(1)+pos(1)+1,center(2)+pos(2), 1)
                if( alpha<beta .and. gamma<beta )then
                    denom = alpha + gamma - 2.*beta
                    if( abs(denom) > TINY ) shift(1) = shift(1) + 0.5 * (alpha-gamma) / denom
                endif
            endif
            ! along y
            if( abs(pos(2)) < itrs )then
                alpha = self1%rmat(center(1)+pos(1),center(2)+pos(2)-1, 1)
                gamma = self1%rmat(center(1)+pos(1),center(2)+pos(2)+1, 1)
                if( alpha<beta .and. gamma<beta )then
                    denom = alpha + gamma - 2.*beta
                    if( abs(denom) > TINY ) shift(2) = shift(2) + 0.5 * (alpha-gamma) / denom
                endif
            endif
            ! convention
            shift = -shift
        else
            shift = -real(pos)
        endif
    end subroutine fcorr_shift

    module subroutine fcorr_shift3D( self1, self2, trs, shift, peak_interp )
        class(image),      intent(inout) :: self1, self2
        real,              intent(in)    :: trs
        real,              intent(inout) :: shift(3)
        logical, optional, intent(in)    :: peak_interp
        real    :: alpha, beta, gamma, denom
        integer :: cen(3), pos(3), itrs
        logical :: l_interp
        if( self1%is_2d() .or. self2%is_2d() ) THROW_HARD('3d only supported')
        if( .not.(self1%is_ft() .and. self2%is_ft()) ) THROW_HARD('FTed only supported')
        if( .not.(self1.eqdims.self2) ) THROW_HARD('Inconsistent dimensions in fcorr_shift')
        l_interp = .false.
        if(present(peak_interp)) l_interp = peak_interp
        ! dimensions
        cen  = self1%ldim/2+1
        itrs = min(floor(trs),minval(cen)-1)
        ! Correlation image
        self1%cmat = self1%cmat * conjg(self2%cmat)
        call self1%ifft
        ! maximum correlation & offset
        pos = maxloc(self1%rmat(cen(1)-itrs:cen(1)+itrs, cen(2)-itrs:cen(2)+itrs, cen(3)-itrs:cen(3)+itrs)) -itrs-1
        ! peak interpolation
        if( l_interp )then
            shift = real(pos)
            ! pos   = pos+itrs+1
            beta  = self1%rmat(cen(1)+pos(1),cen(2)+pos(2), cen(3)+pos(3))
            ! along x
            if( abs(pos(1)) < itrs )then ! within limits
                alpha = self1%rmat(cen(1)+pos(1)-1, cen(2)+pos(2), cen(3)+pos(3))
                gamma = self1%rmat(cen(1)+pos(1)+1, cen(2)+pos(2), cen(3)+pos(3))
                if( alpha<beta .and. gamma<beta )then
                    denom = alpha + gamma - 2.*beta
                    if( abs(denom) > TINY ) shift(1) = shift(1) + 0.5 * (alpha-gamma) / denom
                endif
            endif
            ! along y
            if( abs(pos(2)) < itrs )then
                alpha = self1%rmat(cen(1)+pos(1), cen(2)+pos(2)-1, cen(3)+pos(3))
                gamma = self1%rmat(cen(1)+pos(1), cen(2)+pos(2)+1, cen(3)+pos(3))
                if( alpha<beta .and. gamma<beta )then
                    denom = alpha + gamma - 2.*beta
                    if( abs(denom) > TINY ) shift(2) = shift(2) + 0.5 * (alpha-gamma) / denom
                endif
            endif
            ! along z
            if( abs(pos(3)) < itrs )then
                alpha = self1%rmat(cen(1)+pos(1), cen(2)+pos(2), cen(3)+pos(3)-1)
                gamma = self1%rmat(cen(1)+pos(1), cen(2)+pos(2), cen(3)+pos(3)+1)
                if( alpha<beta .and. gamma<beta )then
                    denom = alpha + gamma - 2.*beta
                    if( abs(denom) > TINY ) shift(3) = shift(3) + 0.5 * (alpha-gamma) / denom
                endif
            endif
            ! convention
            shift = -shift
        else
            shift = -real(pos)
        endif
    end subroutine fcorr_shift3D

    module subroutine prenorm4real_corr_1( self, sxx )
        class(image), intent(inout) :: self
        real,         intent(out)   :: sxx
        real :: npix, ax
        npix = real(product(self%ldim))
        ax   = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))) / npix
        self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) - ax
        sxx  = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))*self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
    end subroutine prenorm4real_corr_1

    module subroutine prenorm4real_corr_2( self, sxx, mask )
        class(image), intent(inout) :: self
        real,         intent(out)   :: sxx
        logical,      intent(in)    :: mask(self%ldim(1),self%ldim(2),self%ldim(3))
        real :: npix, ax
        npix = real(count(mask))
        ax   = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=mask) / npix
        self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) - ax
        sxx  = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))*self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=mask)
    end subroutine prenorm4real_corr_2

    module subroutine prenorm4real_corr_3( self, err )
        class(image), intent(inout) :: self
        logical,      intent(out)   :: err
        real :: npix, ave, var
        npix      = real(product(self%ldim))
        ave       = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))) / real(npix)
        self%rmat = self%rmat - ave
        var       = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))*self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))) / real(npix)
        if( var > TINY )then
            self%rmat = self%rmat / sqrt(var)
            err = .false.
        else
            err = .true.
        endif
    end subroutine prenorm4real_corr_3

    module function real_corr_prenorm_1( self_ref, self_ptcl, sxx_ref ) result( r )
        class(image), intent(inout) :: self_ref, self_ptcl
        real,         intent(in)    :: sxx_ref
        real :: diff(self_ptcl%ldim(1),self_ptcl%ldim(2),self_ptcl%ldim(3))
        real :: r, ay, syy, sxy, npix
        npix = real(product(self_ptcl%ldim))
        ay   = sum(self_ptcl%rmat(:self_ptcl%ldim(1),:self_ptcl%ldim(2),:self_ptcl%ldim(3))) / npix
        diff = self_ptcl%rmat(:self_ptcl%ldim(1),:self_ptcl%ldim(2),:self_ptcl%ldim(3)) - ay
        syy  = sum(diff * diff)
        sxy  = sum(self_ref%rmat(:self_ref%ldim(1),:self_ref%ldim(2),:self_ref%ldim(3)) * diff)
        if( sxx_ref > 0. .or. syy > 0. )then
            r = sxy / sqrt(sxx_ref * syy)
        else
            r = 0.
        endif
    end function real_corr_prenorm_1

    module function real_corr_prenorm_2( self_ref, self_ptcl, sxx_ref, mask ) result( r )
        class(image), intent(inout) :: self_ref, self_ptcl
        real,         intent(in)    :: sxx_ref
        logical,      intent(in)    :: mask(self_ptcl%ldim(1),self_ptcl%ldim(2),self_ptcl%ldim(3))
        real :: diff(self_ptcl%ldim(1),self_ptcl%ldim(2),self_ptcl%ldim(3))
        real :: r, ay, syy, sxy, npix
        npix = real(count(mask))
        ay   = sum(self_ptcl%rmat(:self_ptcl%ldim(1),:self_ptcl%ldim(2),:self_ptcl%ldim(3)), mask=mask) / npix
        where( mask ) diff = self_ptcl%rmat(:self_ptcl%ldim(1),:self_ptcl%ldim(2),:self_ptcl%ldim(3)) - ay
        syy  = sum(diff * diff, mask=mask)
        sxy  = sum(self_ref%rmat(:self_ref%ldim(1),:self_ref%ldim(2),:self_ref%ldim(3)) * diff, mask=mask)
        if( sxx_ref > 0. .or. syy > 0. )then
            r = sxy / sqrt(sxx_ref * syy)
        else
            r = 0.
        endif
    end function real_corr_prenorm_2

    module real function real_corr_prenorm_3( self_ref, self_ptcl )
        class(image), intent(inout) :: self_ref, self_ptcl
        real_corr_prenorm_3 = sum(self_ptcl%rmat(:self_ptcl%ldim(1),:self_ptcl%ldim(2),:self_ptcl%ldim(3))&
                                & * self_ref%rmat(:self_ref%ldim(1),:self_ref%ldim(2),:self_ref%ldim(3)))
        real_corr_prenorm_3 = real_corr_prenorm_3 / product(self_ref%ldim)
    end function real_corr_prenorm_3

    module subroutine radial_cc( self1, self2, self_w, smpd, rad_corrs, rad_dists )
        class(image), intent(inout):: self1, self2, self_w
        real,         intent(in)   :: smpd
        real,         intent(out)  :: rad_corrs(int(self1%ldim(1)/2.)), rad_dists(int(self1%ldim(1)/2.))
        real                 :: rad_weights(int(self1%ldim(1)/2.))
        type(image)          :: distimg
        logical, allocatable :: mask(:,:,:), shell_mask(:,:,:)
        real,    parameter   :: shell_size_pix = 1
        integer :: ldim3, n, n_shells
        real    :: dist_lbound, dist_ubound
        if( .not. (self1.eqdims.self2) ) THROW_HARD('Nonconforming dimensions in image; radial_cc')
        call distimg%new(self1%ldim,smpd)
        n_shells    = int(self1%ldim(1) / 2.)
        if( self1%is_3d() )then
            ! 3D
            ldim3 = self1%ldim(3)
        else
            ! 2D
            ldim3 = 1
        endif
        allocate(mask(self1%ldim(1), self1%ldim(2), ldim3),&
        &  shell_mask(self1%ldim(1), self1%ldim(2), ldim3), source=.true.)
        call distimg%cendist
        do n = 0, n_shells-1
            dist_lbound = real(n) * shell_size_pix
            dist_ubound = dist_lbound + shell_size_pix
            where( (distimg%rmat(:distimg%ldim(1),:distimg%ldim(2),:ldim3) > dist_lbound) .and. &
                  &(distimg%rmat(:distimg%ldim(1),:distimg%ldim(2),:ldim3) < dist_ubound) .and. &
                  &(mask(:distimg%ldim(1),:distimg%ldim(2),:ldim3) ) )
                shell_mask = .true.
            else where
                shell_mask = .false.
            end where
            if( count(shell_mask) < 3 )then
                rad_corrs(n+1) = 0.
            else
                rad_corrs(n+1) = self1%real_corr(self2, shell_mask)
            endif
            rad_dists(n+1) = ( ( dist_lbound * smpd + dist_ubound * smpd ) / 2. )
            if( rad_corrs(n+1)   > 0.      ) rad_weights(n+1) = 2 * rad_corrs(n+1) / (rad_corrs(n+1) +1)
            if( rad_weights(n+1) > 0.99999 ) rad_weights(n+1) = 0.99999
            rad_dists(n+1) = ( ( dist_lbound * smpd + dist_ubound * smpd ) / 2. )
            where( shell_mask(:,:,:) .eqv. .true. )
                self_w%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)) = self_w%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)) + rad_weights(n+1)
                self1%rmat (:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)) = self1%rmat (:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)) * rad_weights(n+1)
            end where
        enddo
    end subroutine radial_cc

    module function sqeuclid( self1, self2, mask ) result( r )
        class(image), intent(inout) :: self1, self2
        logical,      intent(in)    :: mask(self1%ldim(1),self1%ldim(2),self1%ldim(3))
        real :: r
        if( self1%wthreads )then
            !$omp parallel workshare
            r = sum((self1%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)) -&
            &self2%rmat(:self2%ldim(1),:self2%ldim(2),:self2%ldim(3)))**2.0, mask=mask)
            !$omp end parallel workshare
        else
            r = sum((self1%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)) -&
            &self2%rmat(:self2%ldim(1),:self2%ldim(2),:self2%ldim(3)))**2.0, mask=mask)
        endif
    end function sqeuclid

    module subroutine sqeuclid_matrix_1( self1, self2, sqdiff )
        class(image), intent(in)    :: self1, self2
        real,         intent(inout) :: sqdiff(self1%ldim(1),self1%ldim(2),self1%ldim(3))
        sqdiff = (self1%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)) -&
        &self2%rmat(:self2%ldim(1),:self2%ldim(2),:self2%ldim(3)))**2.0
    end subroutine sqeuclid_matrix_1

    module subroutine sqeuclid_matrix_2( self1, self2, sqdiff_img )
        class(image), intent(in)    :: self1, self2
        class(image), intent(inout) :: sqdiff_img
        sqdiff_img%rmat = (self1%rmat - self2%rmat)**2.0
    end subroutine sqeuclid_matrix_2

    module function euclid_norm( self1, self2 ) result( r )
        class(image), intent(inout) :: self1, self2
        real :: r
        r = sqrt(sum((self1%rmat(:self1%ldim(1),:self1%ldim(2),:self1%ldim(3)) -&
        &self2%rmat(:self2%ldim(1),:self2%ldim(2),:self2%ldim(3)))**2.0)) / real(product(self1%ldim))
    end function euclid_norm

    !===========================
    ! cost / shift
    !===========================

    module subroutine opt_filter_costfun( even_filt, odd_raw, odd_filt, even_raw, sqdiff_img )
        class(image), intent(in)    :: even_filt, odd_raw, odd_filt, even_raw
        class(image), intent(inout) :: sqdiff_img
        sqdiff_img%rmat = abs(even_filt%rmat - odd_raw%rmat + odd_filt%rmat - even_raw%rmat)
    end subroutine opt_filter_costfun

    module subroutine opt_filter_costfun_workshare( even_filt, odd_raw, odd_filt, even_raw, sqdiff_img )
        class(image), intent(in)    :: even_filt, odd_raw, odd_filt, even_raw
        class(image), intent(inout) :: sqdiff_img
        !$omp parallel workshare
        sqdiff_img%rmat = abs(even_filt%rmat - odd_raw%rmat + odd_filt%rmat - even_raw%rmat)
        !$omp end parallel workshare
    end subroutine opt_filter_costfun_workshare

    !>  \brief  returns the real and imaginary parts of the phase shift at point
    !!          logi in a Fourier transform caused by the origin shift in shvec
    module pure function oshift_1( self, logi, shvec ) result( comp )
        class(image), intent(in) :: self
        real,         intent(in) :: logi(3)
        real,         intent(in) :: shvec(3)
        complex :: comp
        real    :: arg
        integer :: ldim
        if( self%ldim(3) == 1 )then
            ldim = 2
        else
            ldim = 3
        endif
        arg  = sum(logi(:ldim)*shvec(:ldim)*self%shconst(:ldim))
        comp = cmplx(cos(arg),sin(arg))
    end function oshift_1

    !>  \brief  returns the real and imaginary parts of the phase shift at point
    !!          logi in a Fourier transform caused by the origin shift in shvec
    module pure function oshift_2( self, logi, shvec ) result( comp )
        class(image), intent(in) :: self
        integer,      intent(in) :: logi(3)
        real,         intent(in) :: shvec(3)
        complex :: comp
        comp = self%oshift_1(real(logi), shvec)
    end function oshift_2

    !>  \brief  returns the real argument transfer matrix components at point logi in a Fourier transform
    module function gen_argtransf_comp( self, logi, ldim ) result( arg )
        class(image), intent(in)      :: self
        real, intent(in)              :: logi(3)
        integer, intent(in), optional :: ldim
        real                          :: arg(3)
        integer                       :: lstop, i
        lstop = 2
        if( self%ldim(3) > 1 ) lstop = 3
        if( present(ldim) )    lstop = ldim
        arg = 0.
        do i=1,lstop
            if( self%ldim(i) == 1 )then
                cycle
            else
                if( is_even(self%ldim(i)) )then
                    arg = arg+logi(i)*(PI/real(self%ldim(i)/2.))
                else
                    arg = arg+logi(i)*(PI/real((self%ldim(i)-1)/2.))
                endif
            endif
        end do
    end function gen_argtransf_comp

end submodule simple_image_calc
