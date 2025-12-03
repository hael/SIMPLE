submodule (simple_image) simple_image_ops
!$ use omp_lib
!$ use omp_lib_kinds
include  'simple_lib.f08'
#include "simple_local_flags.inc"
implicit none

contains

    !===========================
    ! CTF
    !===========================

    ! KEEP SERIAL
    module subroutine ctf_dens_correct( self_sum, self_rho )
        class(image), intent(inout) :: self_sum
        class(image), intent(inout) :: self_rho
        integer :: h, k, l, lims(3,2), phys(3), nyq, sh
        real    :: denom
        ! set constants
        lims = self_sum%loop_lims(2)
        nyq  = self_sum%get_lfny(1)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    sh    = nint(hyp(h,k,l))
                    phys  = self_sum%comp_addr_phys([h,k,l])
                    denom = real(self_rho%cmat(phys(1),phys(2),phys(3)))
                    if(sh <= nyq .and. abs(denom) > 1.e-10 )then
                        self_sum%cmat(phys(1),phys(2),phys(3)) = self_sum%cmat(phys(1),phys(2),phys(3)) / denom
                    else
                        self_sum%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                    endif
                end do
            end do
        end do
    end subroutine ctf_dens_correct

    module subroutine ctf_dens_correct_wiener( self_sum, self_rho, ssnr )
        class(image), intent(inout) :: self_sum
        class(image), intent(in)    :: self_rho
        real,         intent(in)    :: ssnr(:)
        integer :: h, k, l, lims(3,2), phys(3), nyq, sh
        real    :: denom
        ! set constants
        lims = self_sum%loop_lims(2)
        nyq  = self_sum%get_lfny(1)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    sh   = nint(hyp(h,k,l))
                    phys = self_sum%comp_addr_phys([h,k,l])
                    if(sh > nyq )then
                        self_sum%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                    else
                        if( sh==0 )then
                            denom = real(self_rho%cmat(phys(1),phys(2),phys(3))) + 1.
                            self_sum%cmat(phys(1),phys(2),phys(3)) = self_sum%cmat(phys(1),phys(2),phys(3))/denom
                        else
                            denom = ssnr(sh)*real(self_rho%cmat(phys(1),phys(2),phys(3))) + 1.
                            self_sum%cmat(phys(1),phys(2),phys(3)) = ssnr(sh)*self_sum%cmat(phys(1),phys(2),phys(3))/denom
                        endif
                    endif
                end do
            end do
        end do
    end subroutine ctf_dens_correct_wiener

    !===========================
    ! Insertions
    !===========================

    !> \brief insert  inserts a box*box particle image into a micrograph
    module subroutine insert(self_in, coord, self_out )
        class(image), intent(in)    :: self_in
        integer,      intent(in)    :: coord(2)
        type(image),  intent(inout) :: self_out
        integer :: xllim, xulim, yllim, yulim
        if( self_in%ldim(3) > 1 )       THROW_HARD('only 4 2D images; insert')
        if( self_in%is_ft() )           THROW_HARD('only 4 real images; insert')
        if( .not. self_in%even_dims() ) THROW_HARD('only 4 even particle dims; insert')
        if( self_out%exists() )then
            if( self_out%ldim(3) > 1 )  THROW_HARD('only 4 2D images; insert')
            if( self_out%is_ft() )      THROW_HARD('only 4 real images; insert')
            if( self_out%ldim(1) > self_in%ldim(1) .and. self_out%ldim(2) > self_in%ldim(2) .and. self_out%ldim(3) == 1 )then
                if( (coord(1) < self_in%ldim(1)/2+1 .or. coord(1) > self_out%ldim(1)-self_in%ldim(1)/2-1) .or.&
                    (coord(2) < self_in%ldim(2)/2+1 .or. coord(2) > self_out%ldim(2)-self_in%ldim(2)/2-1) )then
                    THROW_HARD('particle outside micrograph area; insert')
                endif
            else
                THROW_HARD('micrograph needs to have dimensions larger than the particle; insert')
            endif
        else
            THROW_HARD('micrograph (self_out) does not exist; insert')
        endif
        ! set range
        xllim = coord(1)-self_in%ldim(1)/2
        xulim = coord(1)+self_in%ldim(1)/2-1
        yllim = coord(2)-self_in%ldim(2)/2
        yulim = coord(2)+self_in%ldim(2)/2-1
        ! insert particle image matrix into micrograph image matrix
        self_out%rmat(xllim:xulim,yllim:yulim,1) = self_in%rmat(1:self_in%ldim(1),1:self_in%ldim(2),1)
    end subroutine insert

    ! inserts the low-resolution information from one image into another
    module subroutine insert_lowres( self, self2insert, find )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self2insert
        integer,      intent(in)    :: find
        integer :: lims(3,2), phys(3), h, k, l, sh
        complex :: comp
        if( .not. self%ft        ) THROW_HARD('image to be modified assumed to be FTed; insert_lowres')
        if( .not. self2insert%ft ) THROW_HARD('image to insert assumed to be FTed; insert_lowres')
        lims = self%fit%loop_lims(2)
        !$omp parallel do collapse(3) default(shared) private(h,k,l,sh,phys,comp)&
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    ! find shell
                    sh = nint(hyp(h,k,l))
                    if( sh <= find )then
                        ! insert component
                        phys = self%comp_addr_phys([h,k,l])
                        comp = self2insert%get_fcomp([h,k,l],phys)
                        call self%set_fcomp([h,k,l],phys,comp)
                    endif
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine insert_lowres

    ! inserts the low-resolution information from one image into another
    module subroutine insert_lowres_serial( self, self2insert, find )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self2insert
        integer,      intent(in)    :: find
        integer :: lims(3,2), phys(3), h, k, l, sh
        complex :: comp
        if( .not. self%ft        ) THROW_HARD('image to be modified assumed to be FTed; insert_lowres')
        if( .not. self2insert%ft ) THROW_HARD('image to insert assumed to be FTed; insert_lowres')
        lims = self%fit%loop_lims(2)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    ! find shell
                    sh = nint(hyp(h,k,l))
                    if( sh <= find )then
                        ! insert component
                        phys = self%comp_addr_phys([h,k,l])
                        comp = self2insert%get_fcomp([h,k,l],phys)
                        call self%set_fcomp([h,k,l],phys,comp)
                    endif
                end do
            end do
        end do
    end subroutine insert_lowres_serial

    !===========================
    ! Noise
    !===========================

    module subroutine ran( self, b )
        class(image),   intent(inout) :: self
        real, optional, intent(in)    :: b      ! support upper bound, defaults to 1
        integer :: i, j, k
        do k=1,self%ldim(3)
            do j=1,self%ldim(2)
                do i=1,self%ldim(1)
                    self%rmat(i,j,k) = ran3()
                end do
            end do
        end do
        self%ft = .false.
        if( present(b) ) self%rmat = self%rmat * b
    end subroutine ran

    module subroutine gauran( self, mean, sdev )
        class(image), intent(inout) :: self
        real,         intent(in)    :: mean, sdev
        integer :: i, j, k
        do k=1,self%ldim(3)
            do j=1,self%ldim(2)
                do i=1,self%ldim(1)
                    self%rmat(i,j,k) = gasdev( mean, sdev )
                end do
            end do
        end do
        self%ft = .false.
    end subroutine gauran

    module subroutine add_gauran( self, snr )
        class(image), intent(inout) :: self
        real,         intent(in)    :: snr
        real    :: sdev_noise, var
        integer :: i, j, k
        var = self%variance()
        if( var > TINY )then
            sdev_noise = sqrt(var/snr)
        else
            THROW_HARD('variance of image is zero')
        endif
        do i=1,self%ldim(1)
            do j=1,self%ldim(2)
                do k=1,self%ldim(3)
                    self%rmat(i,j,k) = self%rmat(i,j,k) + gasdev(0., sdev_noise)
                end do
            end do
        end do
    end subroutine add_gauran

    module subroutine salt_n_pepper( self, pos )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: pos(:,:)
        integer :: ipix, jpix
        if( .not. self%is_2d() ) THROW_HARD('only for 2D images; salt_n_pepper')
        call self%norm_minmax
        do ipix=1,self%ldim(1)
            do jpix=1,self%ldim(2)
                if( pos(ipix,jpix) )then
                    if( ran3() < 0.5 )then
                        self%rmat(ipix,jpix,1) = 0.
                    else
                        self%rmat(ipix,jpix,1) = 1.
                    endif
                endif
            end do
        end do
    end subroutine salt_n_pepper

    !===========================
    ! Background
    !===========================

    module subroutine div_w_instrfun( self, interpfun, alpha, padded_dim )
        class(image),           intent(inout) :: self
        character(len =*),      intent(in)    :: interpfun
        real,         optional, intent(in)    :: alpha
        integer,      optional, intent(in)    :: padded_dim
        type(kbinterpol)  :: kbwin
        real, allocatable :: w(:)
        real    :: arg
        integer :: center(3), i,j,k, dim, iarg
        if( any(self%ldim==0) .or. self%is_ft() .or. .not.self%square_dims() )then
            THROW_HARD('Erroneous image in div_w_instrfun')
        endif
        center = self%ldim/2+1
        dim    = self%ldim(1)
        if( present(padded_dim) ) dim = padded_dim
        select case(trim(interpfun))
        case('kb')
            ! kaiser-bessel window
            if(.not.present(alpha)) THROW_HARD('alpha must be given for KB interpolator')
            kbwin = kbinterpol(KBWINSZ,alpha)
            allocate(w(self%ldim(1)),source=1.)
            do i = 1,self%ldim(1)
                arg  = real(i-center(1))/real(dim)
                w(i) = kbwin%instr(arg)
            end do
            if( self%is_2d() )then
                !$omp parallel do collapse(2) private(i,j) default(shared) proc_bind(close) schedule(static)
                do i = 1,self%ldim(1)
                    do j = 1,self%ldim(2)
                        self%rmat(i,j,1) = self%rmat(i,j,1) / (w(i)*w(j))
                    enddo
                enddo
                !$omp end parallel do
            else
                !$omp parallel do collapse(3) private(i,j,k) default(shared) proc_bind(close) schedule(static)
                do i = 1,self%ldim(1)
                    do j = 1,self%ldim(2)
                        do k = 1,self%ldim(3)
                            self%rmat(i,j,k) = self%rmat(i,j,k) / (w(i)*w(j)*w(k))
                        enddo
                    enddo
                enddo
                !$omp end parallel do
            endif
        case('linear')
            ! Tri-linear interpolation
            !$omp parallel do collapse(3) private(i,j,k,iarg,arg) default(shared) proc_bind(close) schedule(static)
            do i = 1,self%ldim(1)
                do j = 1,self%ldim(2)
                    do k = 1,self%ldim(3)
                        iarg = sum(([i,j,k]-center)**2)
                        if( iarg == 0 )cycle
                        arg = PI * sqrt(real(iarg)) / real(dim)
                        arg = sin(arg) / arg
                        arg = arg*arg ! normalized sinc^2
                        self%rmat(i,j,k) = self%rmat(i,j,k) / arg
                    enddo
                enddo
            enddo
            !$omp end parallel do
        case('nn')
            ! Nearest-neighbour interpolation
            !$omp parallel do collapse(3) private(i,j,k,iarg,arg) default(shared) proc_bind(close) schedule(static)
            do i = 1,self%ldim(1)
                do j = 1,self%ldim(2)
                    do k = 1,self%ldim(3)
                        iarg = sum(([i,j,k]-center)**2)
                        if( iarg == 0 )cycle
                        arg = PI * sqrt(real(iarg)) / real(dim)
                        arg = sin(arg) / arg ! normalized sinc
                        self%rmat(i,j,k) = self%rmat(i,j,k) / arg
                    enddo
                enddo
            enddo
            !$omp end parallel do
        case DEFAULT
            THROW_HARD('Unsupported interpolation method')
        end select
    end subroutine div_w_instrfun

    !>  Estimates background from gaussian filtered image, for micrographs
    module subroutine estimate_background( self, freq, backgr, mode )
        class(image),     intent(in)    :: self
        real,             intent(in)    :: freq
        class(image),     intent(inout) :: backgr
        character(len=*), intent(in)    :: mode
        real,    parameter :: PADDING = sqrt(2.)
        integer, parameter :: CS_DIM  = 1024
        type(image)     :: img_pad, msk
        real            :: smpd_bin
        integer         :: ldim_pd(3), ldim_bin(3), ldim_cs(3), bin_dim
        integer         :: i,j,is,js,ie,je, binning
        if( self%ft )      THROW_HARD('Real space only!, estimate_background')
        if( self%is_3d() ) THROW_HARD('2D images only!, estimate_background')
        select case(trim(mode))
        case('cryosparc','cs')
            ! produces backround of shape CS_DIM x CS_DIM
            ! bin micrograph
            ldim_cs = [CS_DIM, CS_DIM,1]
            bin_dim = 2**ceiling( log(real(maxval(self%ldim))) / log(2.0) )
            binning = ceiling(real(bin_dim) / real(CS_DIM))
            ldim_bin(1:2) = nint(self%ldim(1:2) / real(binning))
            ldim_bin(3)   = 1
            smpd_bin      = self%smpd * real(binning)
            call backgr%new(ldim_bin, smpd_bin)
            !$omp parallel do private(i,j,is,ie,js,je) proc_bind(close) collapse(2) default(shared)
            do j = 1,ldim_bin(2)
                do i = 1,ldim_bin(1)
                    is = (i-1)*binning + 1
                    ie = i*binning
                    js = (j-1)*binning + 1
                    je = j*binning
                    backgr%rmat(i,j,1) = sum(self%rmat(is:ie,js:je,1))
                enddo
            enddo
            !$omp end parallel do
            ! low-pass padded mic & crop
            ldim_pd = [2*CS_DIM, 2*CS_DIM, 1]
            call img_pad%new(ldim_pd, smpd_bin)
            call backgr%pad(img_pad)
            call img_pad%fft
            call img_pad%bpgau2D(0., 2.*freq)
            call img_pad%ifft
            call img_pad%clip(backgr)
            ! low pass padded mask & crop
            call msk%new(ldim_cs,smpd_bin)
            msk = 1.
            call msk%pad(img_pad)
            call img_pad%fft
            call img_pad%bpgau2D(0., 2.*freq)
            call img_pad%ifft
            call img_pad%clip(msk)
            ! correct for padding
            !$omp parallel workshare proc_bind(close)
            where( msk%rmat(1:ldim_cs(1),1:ldim_cs(2),1) > TINY )
                backgr%rmat(1:ldim_cs(1),1:ldim_cs(2),1) = &
                    &backgr%rmat(1:ldim_cs(1),1:ldim_cs(2),1) / msk%rmat(1:ldim_cs(1),1:ldim_cs(2),1)
            else where
                backgr%rmat(1:ldim_cs(1),1:ldim_cs(2),1) = 0.
            end where
            !$omp end parallel workshare
        case DEFAULT
            ! produces backround of shape self%ldim
            ! padded dimensions
            ldim_pd(1:2) = nint(PADDING*real(self%ldim(1:2)))
            ldim_pd(1:2) = find_larger_magic_box(ldim_pd(1:2))
            ldim_pd(3)   = 1
            ! low-pass padded mic & crop
            call backgr%copy(self)
            call img_pad%new(ldim_pd, self%smpd)
            call backgr%pad(img_pad)
            call img_pad%fft
            call img_pad%bpgau2D(0.,freq)
            call img_pad%ifft
            call img_pad%clip(backgr)
            ! low pass padded mask & crop
            call msk%new(self%ldim, self%smpd)
            msk = 1.
            call msk%pad(img_pad)
            call img_pad%fft
            call img_pad%bpgau2D(0.,freq)
            call img_pad%ifft
            call img_pad%clip(msk)
            ! correct for padding
            !$omp parallel workshare proc_bind(close)
            where( msk%rmat(1:self%ldim(1),1:self%ldim(2),1) > TINY )
                backgr%rmat(1:self%ldim(1),1:self%ldim(2),1) = &
                    &backgr%rmat(1:self%ldim(1),1:self%ldim(2),1) / msk%rmat(1:self%ldim(1),1:self%ldim(2),1)
            else where
                backgr%rmat(1:self%ldim(1),1:self%ldim(2),1) = 0.
            end where
            !$omp end parallel workshare
        end select
        call img_pad%kill
        call msk%kill
    end subroutine estimate_background

    module subroutine subtr_backgr( self, lp )
        class(image), intent(inout) :: self
        real,         intent(in)    :: lp
        type(image) :: tmp
        integer     :: winsz
        call tmp%copy(self)
        winsz = nint(real(self%ldim(1)/2)*self%smpd / lp / sqrt(2.))
        call tmp%real_space_filter(winsz, 'average')
        self%rmat = self%rmat - tmp%rmat
        call tmp%kill()
    end subroutine subtr_backgr

    !>  subtracts background linear ramp including mean
    module subroutine subtr_backgr_ramp( self, lmsk )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: lmsk(self%ldim(1),self%ldim(2),self%ldim(3))
        real, allocatable :: xyz(:,:)
        real    :: A,B,C,D
        integer :: cen(2),i,j,npix
        logical :: err
        if( self%ft )      THROW_HARD('Real space only!, subtr_backr_ramp')
        if( self%is_3d() ) THROW_HARD('2D images only!, subtr_backr_ramp')
        npix = product(self%ldim) - count(lmsk)
        if( npix < 2 ) return
        allocate(xyz(npix,3))
        cen  = self%ldim(1:2)/2 + 1
        npix = 0
        do j = 1,self%ldim(2)
            do i = 1,self%ldim(1)
                if( lmsk(i,j,1) ) cycle
                npix = npix + 1
                xyz(npix,:) = [real(i-cen(1)), real(j-cen(2)), self%rmat(i,j,1)]
            enddo
        enddo
        call fit_lsq_plane(npix, xyz, A,B,C, err)
        if( err ) return
        do j = 1,self%ldim(2)
            D = B*real(j-cen(2)) + C
            do i = 1,self%ldim(1)
                self%rmat(i,j,1) = self%rmat(i,j,1) - (D + A*real(i-cen(1)))
            enddo
        enddo
    end subroutine subtr_backgr_ramp

    !>  Subtracts background, for micrographs
    module subroutine subtract_background( self, freq, mode )
        class(image),               intent(inout) :: self
        real,                       intent(in)    :: freq
        character(len=*), optional, intent(in)    :: mode
        integer,    parameter :: CS_DIM  = 1024
        type(image)           :: tmpimg
        character(len=STDLEN) :: cmode
        if( self%ft )      THROW_HARD('Real space only!, subtract_background')
        if( self%is_3d() ) THROW_HARD('2D images only!, subtract_background')
        cmode = trim(NIL)
        if( present(mode) ) cmode = mode
        call self%estimate_background(freq, tmpimg, cmode)
        select case(trim(mode))
        case('cryosparc','cs')
            call tmpimg%upsample_square_background(self%ldim(1:2))
        case DEFAULT
            ! all done
        end select
        !$omp parallel workshare proc_bind(close)
        self%rmat = self%rmat - tmpimg%rmat
        !$omp end parallel workshare
        call tmpimg%kill
    end subroutine subtract_background

    module subroutine upsample_square_background( self, ldim )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: ldim(2)
        type(image) :: upsampled
        integer :: i,j,fx,fy,fpx,fpy
        real :: scales(2), scale, x,y, dx,dy,b
        if( ldim(1) /= ldim(2) ) THROW_HARD('Unsupported dimensions!')
        scales = real(self%ldim(1:2)) / real(ldim(1:2))
        scale  = product(scales)
        call upsampled%new([ldim(1), ldim(2), 1], self%smpd/scales(1))
        !$omp parallel do default(shared) collapse(2) proc_bind(close)&
        !$omp private(i,j,x,y,fx,fy,fpx,fpy,dx,dy,b)
        do j = 1,ldim(1)
            do i = 1,ldim(2)
                x   = (real(i-1) * scales(1)) + 1.
                y   = (real(j-1) * scales(2)) + 1.
                fx  = floor(x)
                fy  = floor(y)
                if( fx == self%ldim(1) )then
                    fpx = fx-1
                else
                    fpx = fx+1
                endif
                if( fy == self%ldim(2) )then
                    fpy = fy-1
                else
                    fpy = fy+1
                endif
                dx = x - real(fx)
                dy = y - real(fy)
                b =     self%rmat(fx, fy, 1) * (1.-dx) * (1.-dy)
                b = b + self%rmat(fpx,fy, 1) *     dx  * (1.-dy)
                b = b + self%rmat(fx, fpy,1) * (1.-dx) * dy
                b = b + self%rmat(fpx,fpy,1) *     dx  * dy
                upsampled%rmat(i,j,1) = b * scale
            enddo
        enddo
        !$omp end parallel do
        call self%copy(upsampled)
        call upsampled%kill
    end subroutine upsample_square_background

    !===========================
    ! Arithmetics
    !===========================

    module subroutine remove_neg( self )
        class(image), intent(inout) :: self
        real :: minv
        minv = minval(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)))
        if( minv < 0. )self%rmat = self%rmat + abs(minv)
    end subroutine remove_neg

    module subroutine neg( self )
        class(image), intent(inout) :: self
        logical :: didft
        didft = .false.
        if( self%ft )then
        else
            call self%fft()
            didft = .true.
        endif
        call self%mul(-1.)
        if( didft ) call self%ifft()
    end subroutine neg

    module subroutine div_below( self, thres, val )
        class(image), intent(inout) :: self
        real,         intent(in)    :: thres, val
        if( is_equal(val,0.) ) return
        where( self%rmat < thres ) self%rmat = self%rmat / val
    end subroutine div_below

    module subroutine inv( self )
        class(image), intent(inout) :: self
        self%rmat = -1.*self%rmat
    end subroutine inv

    !===========================
    ! Zeroing
    !===========================

    module subroutine zero(self)
        class(image), intent(inout) :: self
        if( self%ft )then
            self%cmat = cmplx(0.,0.)
        else
            self%rmat = 0.
        endif
    end subroutine zero

    module subroutine zero_and_flag_ft(self)
        class(image), intent(inout) :: self
        self%cmat = cmplx(0.,0.)
        self%ft   = .true.
    end subroutine zero_and_flag_ft

    module subroutine zero_and_unflag_ft(self)
        class(image), intent(inout) :: self
        self%rmat = 0.
        self%ft   = .false.
    end subroutine zero_and_unflag_ft

    ! estimates median of background along edges of box and subtracts it to flatten the background
    ! modifies pixels along the edges of the box, which ought to be safe as we are masking
    module subroutine zero_background( self )
        class(image), intent(inout) :: self
        integer :: k
        real    :: med, val1,val2,val3,val4,val5,val6,val7,val8,val9,val10,val11,val12
        k = self%ldim(1)/2
        if( self%ldim(3) == 1 )then
            val1  = selec(k,self%ldim(1),self%rmat( 1           , :self%ldim(2),1))
            val2  = selec(k,self%ldim(1),self%rmat( self%ldim(1), :self%ldim(2),1))
            val3  = selec(k,self%ldim(1),self%rmat(:self%ldim(1),  1,           1))
            val4  = selec(k,self%ldim(1),self%rmat(:self%ldim(1),  self%ldim(2),1))
            med   = (val1+val2+val3+val4) / 4.
        else
            val1  = selec(k,self%ldim(1),self%rmat( 1           ,  1            , :self%ldim(3)))
            val2  = selec(k,self%ldim(1),self%rmat( 1           ,  self%ldim(2) , :self%ldim(3)))
            val3  = selec(k,self%ldim(1),self%rmat( self%ldim(1),  1            , :self%ldim(3)))
            val4  = selec(k,self%ldim(1),self%rmat( self%ldim(1),  self%ldim(2) , :self%ldim(3)))
            val5  = selec(k,self%ldim(1),self%rmat( 1           , :self%ldim(2) ,  1           ))
            val6  = selec(k,self%ldim(1),self%rmat( 1           , :self%ldim(2) ,  self%ldim(3)))
            val7  = selec(k,self%ldim(1),self%rmat( self%ldim(1), :self%ldim(2) ,  1           ))
            val8  = selec(k,self%ldim(1),self%rmat( self%ldim(1), :self%ldim(2) ,  self%ldim(3)))
            val9  = selec(k,self%ldim(1),self%rmat(:self%ldim(1),  1            ,  1           ))
            val10 = selec(k,self%ldim(1),self%rmat(:self%ldim(1),  1            ,  self%ldim(3)))
            val11 = selec(k,self%ldim(1),self%rmat(:self%ldim(1),  self%ldim(2) ,  1           ))
            val12 = selec(k,self%ldim(1),self%rmat(:self%ldim(1),  self%ldim(2) ,  self%ldim(3)))
            med   = (val1+val2+val3+val4+val5+val6+val7+val8+val9+val10+val11+val12) / 12.
        endif
        if(abs(med) > TINY) self%rmat = self%rmat - med
    end subroutine zero_background

    module subroutine zero_below( self, thres )
        class(image), intent(inout) :: self
        real,         intent(in)    :: thres
        where( self%rmat < thres ) self%rmat = 0.
    end subroutine zero_below

    !>  \brief  putting the edge around the image to zero (necessary for avoiding FT artefacts)
    module subroutine zero_edgeavg( self )
        class(image), intent(inout) :: self
        real :: edges_sum, edges_ave
        if( self%ft )           THROW_HARD('not for Fted images; zero_edgeavg')
        if( .not.self%is_2d() ) THROW_HARD('only for 2d images; zero_edgeavg')
        edges_sum = sum(self%rmat(1:self%ldim(1),1,1))
        edges_sum = edges_sum + sum(self%rmat(1:self%ldim(1),self%ldim(2),1))
        edges_sum = edges_sum + sum(self%rmat(1,1:self%ldim(2),1))
        edges_sum = edges_sum + sum(self%rmat(self%ldim(1),1:self%ldim(2),1))
        edges_ave = edges_sum / real( 2*(self%ldim(1)+self%ldim(2)) )
        if( abs(edges_ave) > TINY )self%rmat = self%rmat - edges_ave
    end subroutine zero_edgeavg

    ! substracts median of background defined by an envelope mask (eg 0< <1)
    module subroutine zero_env_background( self, volmsk )
        class(image), intent(inout) :: self, volmsk
        real, allocatable :: vals(:)
        real              :: med
        integer           :: cnt, npix, npix_env, i,j,k
        logical           :: env_mask(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))
        npix = product(self%ldim)
        where( volmsk%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) > 0.0001&
        &.and. volmsk%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) < 0.9999 )
            env_mask = .true.
        else where
            env_mask = .false.
        end where
        npix_env = count(env_mask)
        if( npix_env==0 .or. npix_env==npix )then
            THROW_HARD('Volume mask is not a volume mask; simple_image :: zero_env_background')
        endif
        allocate(vals(npix_env))
        cnt = 0
        do i=1,self%ldim(1)
            do j=1,self%ldim(2)
                do k=1,self%ldim(3)
                    if(env_mask(i,j,k))then
                        cnt = cnt + 1
                        vals(cnt) = self%rmat(i,j,k)
                    endif
                enddo
            enddo
        enddo
        med = median_nocopy(vals)
        if(abs(med) > TINY) self%rmat = self%rmat - med
        deallocate(vals)
    end subroutine zero_env_background

    module subroutine zero_neg( self )
        class(image), intent(inout) :: self
        where( self%rmat < TINY ) self%rmat = 0.
    end subroutine zero_neg

end submodule simple_image_ops
