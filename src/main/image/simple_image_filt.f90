submodule (simple_image) simple_image_filt
!$ use omp_lib
!$ use omp_lib_kinds
include  'simple_lib.f08'
#include "simple_local_flags.inc"
implicit none
contains

    module subroutine bp( self, hplim, lplim, width )
        class(image), intent(inout) :: self
        real, intent(in)            :: hplim, lplim
        real, intent(in), optional  :: width
        integer :: h, k, l, lims(3,2), phys(3)
        logical :: didft, dohp, dolp
        real    :: freq, hplim_freq, lplim_freq, wwidth, w
        wwidth =10.
        if( present(width) ) wwidth = width
        didft = .false.
        if( .not. self%ft )then
            call self%fft()
            didft = .true.
        endif
        dohp = abs(hplim) > TINY
        dolp = abs(lplim) > TINY
        hplim_freq = self%fit%get_find(1,hplim)
        lplim_freq = self%fit%get_find(1,lplim)
        lims = self%fit%loop_lims(2)
        if( self%wthreads )then
            !$omp parallel do private(h,k,l,freq,phys,w) default(shared)&
            !$omp collapse(3) proc_bind(close)
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        freq = hyp(h,k,l)
                        phys = self%comp_addr_phys([h,k,l])
                        if( dohp )then
                            if(freq .lt. hplim_freq) then
                                self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                            else if(freq .le. hplim_freq + wwidth) then
                                w = (1.-cos(((freq-hplim_freq)/wwidth)*pi))/2.
                                self%cmat(phys(1),phys(2),phys(3)) = &
                                    &self%cmat(phys(1),phys(2),phys(3)) * w
                            endif
                        endif
                        if( dolp )then
                            if(freq .gt. lplim_freq)then
                                self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                            else if(freq .ge. lplim_freq - wwidth)then
                                w = (cos(((freq-(lplim_freq-wwidth))/wwidth)*pi)+1.)/2.
                                self%cmat(phys(1),phys(2),phys(3)) = &
                                    &self%cmat(phys(1),phys(2),phys(3)) * w
                            endif
                        endif
                    end do
                end do
            end do
            !$omp end parallel do
        else
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        freq = hyp(h,k,l)
                        phys = self%comp_addr_phys([h,k,l])
                        if( dohp )then
                            if(freq .lt. hplim_freq) then
                                self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                            else if(freq .le. hplim_freq + wwidth) then
                                w = (1.-cos(((freq-hplim_freq)/wwidth)*pi))/2.
                                self%cmat(phys(1),phys(2),phys(3)) = &
                                    &self%cmat(phys(1),phys(2),phys(3)) * w
                            endif
                        endif
                        if( dolp )then
                            if(freq .gt. lplim_freq)then
                                self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                            else if(freq .ge. lplim_freq - wwidth)then
                                w = (cos(((freq-(lplim_freq-wwidth))/wwidth)*pi)+1.)/2.
                                self%cmat(phys(1),phys(2),phys(3)) = &
                                    &self%cmat(phys(1),phys(2),phys(3)) * w
                            endif
                        endif
                    end do
                end do
            end do
        endif
        if( didft ) call self%ifft()
    end subroutine bp

    module subroutine lp( self, find, width )
        class(image),   intent(inout) :: self
        integer,        intent(in)    :: find
        real, optional, intent(in)    :: width
        integer :: h, k, l, lims(3,2), phys(3)
        real    :: freq, wwidth, w
        wwidth =10.
        if( present(width) ) wwidth = width
        lims = self%fit%loop_lims(2)
        !$omp parallel do private(h,k,l,freq,phys,w) default(shared)&
        !$omp collapse(3) proc_bind(close)
        do l=lims(3,1),lims(3,2)
            do k=lims(2,1),lims(2,2)
                do h=lims(1,1),lims(1,2)
                    freq = hyp(h,k,l)
                    phys = self%comp_addr_phys([h,k,l])
                    if(freq .gt. find)then
                        self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                    else if(freq .ge. find - wwidth)then
                        w = (cos(((freq-(find-wwidth))/wwidth)*pi)+1.)/2.
                        self%cmat(phys(1),phys(2),phys(3)) = &
                            &self%cmat(phys(1),phys(2),phys(3)) * w
                    endif
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine lp

    module subroutine lp_background( self, mskvol, lp )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: mskvol
        real,         intent(in)    :: lp
        type(image) :: weights, self_filt
        if( self%is_ft() ) THROW_HARD('only 4 real images; lp_background')
        if( self%ldim(3) == 1 ) THROW_HARD('only 4 volumes; lp_background')
        if( any((self%ldim-mskvol%ldim)/=0) ) THROW_HARD('inconsistent image/msk dimensions; lp_background')
        call self%zero_background
        call self_filt%new(self%ldim,self%smpd)
        call self_filt%copy(self)
        call weights%new(self%ldim,1.)
        call weights%copy(mskvol)
        ! self
        call self%mul(weights)
        ! self low-pass
        call self_filt%bp(0., lp)
        weights%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) =&
            &1. - weights%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))
        call self_filt%mul(weights)
        ! addition
        call self%add(self_filt)
        ! clean
        call weights%kill
        call self_filt%kill
    end subroutine lp_background

    module subroutine bpgau2D( self, hp, lp )
        class(image), intent(inout) :: self
        real,         intent(in)    :: hp, lp
        real    :: hp_fwhm(2), hp_halfinvsigsq(2), hpa
        real    :: lp_fwhm(2), lp_halfinvsigsq(2), lpa
        integer :: phys(2), lims(3,2), h,k
        logical :: l_hp, l_lp
        if(.not.self%ft) THROW_HARD('Input image must be in the reciprocal domain')
        if(.not.self%is_2d()) THROW_HARD('Input image must be two-dimensional')
        lims = self%fit%loop_lims(2)
        l_hp = .false.
        if( hp > TINY )then
            l_hp = .true.
            hp_fwhm         = hp / self%smpd / real(self%ldim(1:2))
            hp_halfinvsigsq = 0.5 * (PI * 2.0 * hp_fwhm / 2.35482)**2
        endif
        l_lp = .false.
        if( lp > TINY )then
            l_lp = .true.
            lp_fwhm         = lp / self%smpd / real(self%ldim(1:2))
            lp_halfinvsigsq = 0.5 * (PI * 2.0 * lp_fwhm / 2.35482)**2
        endif
        !$omp parallel do collapse(2) schedule(static) default(shared) proc_bind(close)&
        !$omp private(h,k,hpa,lpa,phys)
        do h = lims(1,1),lims(1,2)
            do k = lims(2,1),lims(2,2)
                phys = self%comp_addr_phys(h,k)
                if( l_hp )then
                    hpa  = real(h*h) * hp_halfinvsigsq(1) + real(k*k) * hp_halfinvsigsq(2)
                    self%cmat(phys(1),phys(2),1) = self%cmat(phys(1),phys(2),1) * (1.0-exp(-hpa))
                endif
                if( l_lp )then
                    lpa  = real(h*h) * lp_halfinvsigsq(1) + real(k*k) * lp_halfinvsigsq(2)
                    self%cmat(phys(1),phys(2),1) = self%cmat(phys(1),phys(2),1) * exp(-lpa)
                endif
            enddo
        enddo
        !$omp end parallel do
    end subroutine bpgau2D

    module subroutine bpgau3D( self, hp, lp )
        class(image), intent(inout) :: self
        real,         intent(in)    :: hp, lp
        real    :: hp_fwhm(3), hp_halfinvsigsq(3), hpa
        real    :: lp_fwhm(3), lp_halfinvsigsq(3), lpa
        integer :: phys(3), lims(3,2), h,k,l
        logical :: l_hp, l_lp
        if(.not.self%ft) THROW_HARD('Input image must be in the reciprocal domain')
        if(.not.self%is_3d()) THROW_HARD('Input image must be three-dimensional')
        lims = self%fit%loop_lims(2)
        l_hp = .false.
        if( hp > TINY )then
            l_hp = .true.
            hp_fwhm         = hp / self%smpd / real(self%ldim)
            hp_halfinvsigsq = 0.5 * (PI * 2.0 * hp_fwhm / 2.35482)**2
        endif
        l_lp = .false.
        if( lp > TINY )then
            l_lp = .true.
            lp_fwhm         = lp / self%smpd / real(self%ldim)
            lp_halfinvsigsq = 0.5 * (PI * 2.0 * lp_fwhm / 2.35482)**2
        endif
        !$omp parallel do collapse(3) schedule(static) default(shared) proc_bind(close)&
        !$omp private(h,k,l,hpa,lpa,phys)
        do h = lims(1,1),lims(1,2)
            do k = lims(2,1),lims(2,2)
                do l = lims(3,1),lims(3,2)
                    phys = self%comp_addr_phys(h,k,l)
                    if( l_hp )then
                        hpa  = sum(hp_halfinvsigsq * real([h,k,l]**2))
                        self%cmat(phys(1),phys(2),phys(3)) =&
                        &self%cmat(phys(1),phys(2),phys(3)) * (1.0-exp(-hpa))
                    endif
                    if( l_lp )then
                        lpa  = sum(lp_halfinvsigsq * real([h,k,l]**2))
                        self%cmat(phys(1),phys(2),phys(3)) =&
                        &self%cmat(phys(1),phys(2),phys(3)) * exp(-lpa)
                    endif
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine bpgau3D

    module subroutine tophat( self, shell, halfwidth )
        class(image),   intent(inout) :: self
        integer,        intent(in)    :: shell
        real, optional, intent(in)    :: halfwidth
        integer :: h, k, l, lims(3,2), phys(3)
        logical :: didft
        real    :: freq, hplim_freq, lplim_freq, hwdth
        hwdth = 0.5
        if( present(halfwidth) ) hwdth = halfwidth
        didft = .false.
        if( .not. self%ft )then
            call self%fft()
            didft = .true.
        endif
        hplim_freq = max(0., real(shell) - hwdth)
        lplim_freq = real(shell) + hwdth
        lims = self%fit%loop_lims(2)
        if( self%wthreads )then
            !$omp parallel do private(h,k,l,freq,phys) default(shared)&
            !$omp collapse(3) proc_bind(close)
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        freq = hyp(h,k,l)
                        phys = self%comp_addr_phys([h,k,l])
                        if(freq .lt. hplim_freq) then
                            self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                        else if(freq .gt. lplim_freq) then
                            self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                        endif
                    end do
                end do
            end do
            !$omp end parallel do
        else
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        freq = hyp(h,k,l)
                        phys = self%comp_addr_phys([h,k,l])
                        if(freq .lt. hplim_freq) then
                            self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                        else if(freq .gt. lplim_freq) then
                            self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                        endif
                    end do
                end do
            end do
        endif
        if( didft ) call self%ifft()
    end subroutine tophat

    module subroutine phase_rand( self, lp )
        class(image), intent(inout) :: self
        real,         intent(in)    :: lp
        integer                     :: h,k,l,phys(3),lims(3,2)
        logical                     :: didft
        real                        :: freq,lp_freq, amp,phase
        real, parameter             :: errfrac=0.5
        didft = .false.
        if( .not. self%ft )then
            call self%fft()
            didft = .true.
        endif
        lp_freq = real(self%fit%get_find(1,lp)) ! assuming square 4 now
        lims    = self%fit%loop_lims(2)
        !$omp parallel do collapse(3) default(shared) private(h,k,l,freq,phys,amp,phase)&
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    freq = hyp(h,k,l)
                    if(freq .gt. lp_freq)then
                        phys  = self%fit%comp_addr_phys([h,k,l])
                        amp   = mycabs(self%cmat(phys(1),phys(2),phys(3)))
                        phase = ran3() * TWOPI
                        self%cmat(phys(1),phys(2),phys(3)) = amp * cmplx(cos(phase), sin(phase))
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        if( didft ) call self%ifft()
    end subroutine phase_rand

    module subroutine ran_phases_below_noise_power( self_even, self_odd )
        class(image), intent(inout) :: self_even, self_odd
        integer :: h, k, l, lims(3,2), phys(3)
        real    :: noise_pow, even_pow, odd_pow, phase
        complex :: diff
        lims  = self_even%fit%loop_lims(2)
        if( .not.self_even%is_ft() ) THROW_HARD('even image needs to be FTed')
        if( .not.self_odd%is_ft()  ) THROW_HARD('odd  image needs to be FTed')
        !$omp parallel do collapse(3) default(shared) private(h,k,l,phys,diff,noise_pow,even_pow,odd_pow,phase)&
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    phys      = self_even%comp_addr_phys([h,k,l])
                    diff      = self_even%cmat(phys(1),phys(2),phys(3)) -&
                               &self_odd%cmat(phys(1),phys(2),phys(3))
                    noise_pow = csq_fast(diff)
                    even_pow  = csq_fast(self_even%cmat(phys(1),phys(2),phys(3)))
                    odd_pow   = csq_fast(self_odd%cmat(phys(1),phys(2),phys(3)))
                    if( noise_pow > even_pow .or. noise_pow > odd_pow )then
                        phase = ran3() * TWOPI
                        self_even%cmat(phys(1),phys(2),phys(3)) = mycabs(self_even%cmat(phys(1),phys(2),phys(3))) * cmplx(cos(phase), sin(phase))
                        phase = ran3() * TWOPI
                        self_odd%cmat(phys(1),phys(2),phys(3))  = mycabs(self_odd%cmat(phys(1),phys(2),phys(3)))  * cmplx(cos(phase), sin(phase))
                    endif
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine ran_phases_below_noise_power

    module subroutine whiten_noise_power( self_even, self_odd, is_ptcl )
        class(image), intent(inout) :: self_even, self_odd
        logical,      intent(in)    :: is_ptcl
        real(dp) :: counts(fdim(self_even%ldim(1)) - 1)
        real(dp) ::  dspec(fdim(self_even%ldim(1)) - 1)
        complex  :: diff
        integer  :: filtsz, h, k, l, sh, lims(3,2), phys(3)
        filtsz = fdim(self_even%ldim(1)) - 1
        dspec  = 0.d0
        counts = 0.d0
        lims   = self_even%fit%loop_lims(2)
        do l = lims(3,1),lims(3,2)
            do k = lims(2,1),lims(2,2)
                do h = lims(1,1),lims(1,2)
                    phys = self_even%fit%comp_addr_phys([h,k,l])
                    sh   = nint(hyp(h,k,l))
                    if( sh == 0 .or. sh > filtsz ) cycle
                    diff      = self_even%cmat(phys(1),phys(2),phys(3)) -&
                                &self_odd%cmat(phys(1),phys(2),phys(3))
                    dspec(sh)  = dspec(sh)  + csq_fast(dcmplx(diff))
                    counts(sh) = counts(sh) + 1.d0
                end do
            end do
        end do        
        if( is_ptcl )then
            where(counts > DTINY) dspec =         dspec / counts
        else
            where(counts > DTINY) dspec = 0.5d0 * dspec / counts ! 0.5 because of the e/o split
        endif
        do l = lims(3,1),lims(3,2)
            do k = lims(2,1),lims(2,2)
                do h = lims(1,1),lims(1,2)
                    phys = self_even%fit%comp_addr_phys([h,k,l])
                    sh   = nint(hyp(h,k,l))
                    if( sh == 0 .or. sh > filtsz ) cycle
                    if( dspec(sh) > DTINY )then
                                            self_even%cmat(phys(1),phys(2),phys(3)) = &
                                           &self_even%cmat(phys(1),phys(2),phys(3)) / real(dsqrt(dspec(sh)),kind=sp)
                        if( .not. is_ptcl ) self_odd%cmat( phys(1),phys(2),phys(3)) = &
                                           &self_odd%cmat( phys(1),phys(2),phys(3)) / real(dsqrt(dspec(sh)),kind=sp)
                    endif
                end do
            end do
        end do
    end subroutine whiten_noise_power

    module subroutine real_space_filter( self, winsz, which )
        class(image),     intent(inout) :: self
        integer,          intent(in)    :: winsz
        character(len=*), intent(in)    :: which
        real, allocatable     :: pixels(:), wfvals(:)
        integer               :: n, i, j, k, cnt, npix
        real                  :: rn, wfun(-winsz:winsz), norm, avg, sdev
        type(winfuns)         :: fwin
        character(len=STDLEN) :: wstr
        type(image)           :: img_filt
        ! check the number of pixels in window
        if( self%is_3d() )then
            npix = (2*winsz+1)**3
        else
            npix = (2*winsz+1)**2
        endif
        pixels = self%win2arr(1, 1, 1, winsz)
        n = size(pixels)
        rn = real(n)
        allocate(wfvals(n))
        ! make the window function
        wstr = 'bman'
        fwin = winfuns(wstr, real(WINSZ), 1.0)
        ! sample the window function
        do i=-winsz,winsz
            wfun(i) = fwin%eval_apod(real(i))
        end do
        ! memoize wfun vals & normalisation constant
        norm = 0.
        cnt  = 0
        if( self%ldim(3) == 1 )then
            do i=-winsz,winsz
                do j=-winsz,winsz
                    cnt = cnt + 1
                    wfvals(cnt) = wfun(i) * wfun(j)
                    norm = norm + wfvals(cnt)
                end do
            end do
        else
            if (which == 'bman')then
                do i=-winsz,winsz
                    do j=-winsz,winsz
                        do k=-winsz,winsz
                            cnt = cnt + 1
                            wfvals(cnt) = wfun(i) * wfun(j) * wfun(k)
                            norm = norm + wfvals(cnt)
                        end do
                    end do
                end do
            else
                !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static) proc_bind(close)&
                !$omp reduction(+:norm)
                do i=-winsz,winsz
                    do j=-winsz,winsz
                        do k=-winsz,winsz
                            norm = norm + wfun(i) * wfun(j) * wfun(k)
                        end do
                    end do
                end do
                !$omp end parallel do
            endif
        endif
        ! make the output image
        call img_filt%new(self%ldim, self%smpd)
        ! filter
        if( self%ldim(3) == 1 )then
            select case(which)
            case('median')
                !$omp parallel do collapse(2) default(shared) private(i,j,pixels) schedule(static) proc_bind(close)
                do i=1,self%ldim(1)
                    do j=1,self%ldim(2)
                        pixels = self%win2arr(i, j, 1, winsz)
                        img_filt%rmat(i,j,1) = median_nocopy(pixels)
                    end do
                end do
                !$omp end parallel do
            case('average')
                !$omp parallel do collapse(2) default(shared) private(i,j,pixels) schedule(static) proc_bind(close)
                do i=1,self%ldim(1)
                    do j=1,self%ldim(2)
                        pixels = self%win2arr(i, j, 1, winsz)
                        img_filt%rmat(i,j,1) = sum(pixels)/rn
                    end do
                end do
                !$omp end parallel do
            case('stdev')
                !$omp parallel do collapse(2) default(shared) private(i,j,pixels,avg,sdev) schedule(static) proc_bind(close)
                do i=1,self%ldim(1)
                    do j=1,self%ldim(2)
                        pixels = self%win2arr(i, j, 1, winsz)
                        call avg_sdev(pixels, avg, sdev)
                        img_filt%rmat(i,j,1) = sdev
                    end do
                end do
                !$omp end parallel do
            case('bman')
                !$omp parallel do collapse(2) default(shared) private(i,j,pixels) schedule(static) proc_bind(close)
                do i=1,self%ldim(1)
                    do j=1,self%ldim(2)
                        pixels = self%win2arr(i, j, 1, winsz)
                        img_filt%rmat(i,j,1) = sum(pixels * wfvals) / norm
                    end do
                end do
                !$omp end parallel do
            case('NLmean')
                call self%NLmean2D()
                img_filt%rmat = self%rmat
            case DEFAULT
                THROW_HARD('unknown filter type; real_space_filter')
            end select
        else !3D
            select case(which)
            case('median')
                !$omp parallel do collapse(3) default(shared) private(i,j,k,pixels) schedule(static) proc_bind(close)
                do i=1,self%ldim(1)
                    do j=1,self%ldim(2)
                        do k=1,self%ldim(3)
                            pixels = self%win2arr(i, j, k, winsz)
                            img_filt%rmat(i,j,k) = median_nocopy(pixels)
                        end do
                    end do
                end do
                !$omp end parallel do
            case('average')
                !$omp parallel do collapse(3) default(shared) private(i,j,k,pixels) !schedule(static) proc_bind(close)
                do i=1,self%ldim(1)
                    do j=1,self%ldim(2)
                        do k=1,self%ldim(3)
                            pixels = self%win2arr(i, j, k, winsz)
                            img_filt%rmat(i,j,k) = sum(pixels)/rn
                        end do
                    end do
                end do
                !$omp end parallel do
            case('stdev')
                !$omp parallel do collapse(3) default(shared) private(i,j,k,pixels,avg,sdev) schedule(static) proc_bind(close)
                do i=1,self%ldim(1)
                    do j=1,self%ldim(2)
                        do k=1,self%ldim(3)
                            pixels = self%win2arr(i, j, k, winsz)
                            call avg_sdev(pixels, avg, sdev)
                            img_filt%rmat(i,j,k) = sdev
                        end do
                    end do
                end do
                !$omp end parallel do
            case('bman')
                !$omp parallel do collapse(3) default(shared) private(i,j,k,pixels) schedule(static) proc_bind(close)
                do i=1,self%ldim(1)
                    do j=1,self%ldim(2)
                        do k=1,self%ldim(3)
                            pixels = self%win2arr(i, j, k, winsz)
                            img_filt%rmat(i,j,k) = sum(pixels * wfvals) / norm
                        end do
                    end do
                end do
                !$omp end parallel do
            case DEFAULT
                THROW_HARD('unknown filter type; real_space_filter')
            end select
        endif
        call self%copy(img_filt)
        call img_filt%kill()
    end subroutine real_space_filter

    module function hannw( self, oshoot_in ) result( w )
        class(image), intent(inout) :: self
        real, intent(in), optional  :: oshoot_in
        integer                     :: lims(3,2), k, kmax, maxl
        type(winfuns)               :: wfuns
        character(len=STDLEN)       :: wstr
        real, allocatable           :: w(:)
        real                        :: oshoot
        oshoot = 0.3
        if( present(oshoot_in) ) oshoot = oshoot_in
        lims = self%loop_lims(2)
        maxl = maxval(lims)
        kmax = maxl+int(oshoot*real(maxl))
        allocate( w(kmax) )
        wstr = 'hann'
        wfuns = winfuns(wstr, real(kmax), 2.)
        do k=1,kmax
            w(k) = wfuns%eval_apod(real(k))
        end do
    end function hannw

    module subroutine apply_bfac( self, b )
        class(image), intent(inout) :: self
        real,         intent(in)    :: b
        integer                     :: i,j,k,phys(3),lims(3,2)
        real                        :: wght, res
        logical                     :: didft
        didft = .false.
        if( .not. self%ft )then
            call self%fft()
            didft = .true.
        endif
        lims = self%fit%loop_lims(2)
        if( self%wthreads )then
            !$omp parallel do collapse(3) default(shared) proc_bind(close)&
            !$omp private(k,j,i,res,phys,wght) schedule(static)
            do k=lims(3,1),lims(3,2)
                do j=lims(2,1),lims(2,2)
                    do i=lims(1,1),lims(1,2)
                        res = sqrt(real(k*k+j*j+i*i))/(real(self%ldim(1))*self%smpd) ! assuming square dimensions
                        phys = self%fit%comp_addr_phys([i,j,k])
                        wght = max(0.,exp(-(b/4.)*res*res))
                        self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))*wght
                    end do
                end do
            end do
            !$omp end parallel do
        else
            do k=lims(3,1),lims(3,2)
                do j=lims(2,1),lims(2,2)
                    do i=lims(1,1),lims(1,2)
                        res = sqrt(real(k*k+j*j+i*i))/(real(self%ldim(1))*self%smpd) ! assuming square dimensions
                        phys = self%fit%comp_addr_phys([i,j,k])
                        wght = max(0.,exp(-(b/4.)*res*res))
                        self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))*wght
                    end do
                end do
            end do
        endif
        if( didft ) call self%ifft()
    end subroutine apply_bfac

    module subroutine apply_filter_1( self, filter )
        class(image), intent(inout) :: self
        real,         intent(in)    :: filter(:)
        integer :: nyq, sh, h, k, l, lims(3,2)
        logical :: didft
        real    :: fwght, wzero
        nyq = size(filter)
        didft = .false.
        if( .not. self%ft )then
            call self%fft()
            didft = .true.
        endif
        wzero = maxval(filter)
        lims = self%fit%loop_lims(2)
        !$omp parallel do collapse(3) default(shared) private(h,k,l,sh,fwght)&
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    ! find shell
                    sh = nint(hyp(h,k,l))
                    ! set filter weight
                    if( sh > nyq )then
                        fwght = 0.
                    else if( sh == 0 )then
                        fwght = wzero
                    else
                        fwght = filter(sh)
                    endif
                    ! multiply with the weight
                    call self%mul([h,k,l], fwght)
                end do
            end do
        end do
        !$omp end parallel do
        if( didft ) call self%ifft()
    end subroutine apply_filter_1

    module subroutine apply_filter_2( self, filter )
        class(image), intent(inout) :: self
        complex,      intent(in)    :: filter(:)
        integer :: nyq, sh, h, k, l, lims(3,2)
        logical :: didft
        complex :: fwght, wzero, wnyq
        nyq = size(filter)
        didft = .false.
        if( .not. self%ft )then
            call self%fft()
            didft = .true.
        endif
        wzero = filter(1)
        wnyq  = cmplx(0.,0.)
        lims = self%fit%loop_lims(2)
        !$omp parallel do collapse(3) default(shared) private(h,k,l,sh,fwght)&
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    ! find shell
                    sh = nint(hyp(h,k,l))
                    ! set filter weight
                    if( sh > nyq )then
                        fwght = wnyq
                    else if( sh == 0 )then
                        fwght = wzero
                    else
                        fwght = filter(sh)
                    endif
                    ! multiply with the weight
                    call self%mul([h,k,l], fwght)
                end do
            end do
        end do
        !$omp end parallel do
        if( didft ) call self%ifft()
    end subroutine apply_filter_2

    module subroutine apply_filter_serial( self, filter )
        class(image), intent(inout) :: self
        real,         intent(in)    :: filter(:)
        integer :: nyq, sh, h, k, l, lims(3,2)
        real    :: fwght, wzero
        nyq   = size(filter)
        wzero = maxval(filter)
        lims  = self%fit%loop_lims(2)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    ! find shell
                    sh = nint(hyp(h,k,l))
                    ! set filter weight
                    if( sh > nyq )then
                        fwght = 0.
                    else if( sh == 0 )then
                        fwght = wzero
                    else
                        fwght = filter(sh)
                    endif
                    ! multiply with the weight
                    call self%mul([h,k,l], fwght)
                end do
            end do
        end do
    end subroutine apply_filter_serial

    ! don't touch default parameters here. This routine is being actively used
    module subroutine NLmean2D( self, msk, sdev_noise )
        class(image),   intent(inout) :: self
        real, optional, intent(in)    :: msk
        real, optional, intent(in)    :: sdev_noise
        real,  allocatable :: rmat_pad(:,:), rmat_threads(:,:,:,:)
        integer, parameter :: DIM_SW  = 3   ! good if use Euclidean distance
        integer, parameter :: CFR_BOX = 10  ! as suggested in the paper, use a box 21x21
        real    :: exponentials(-CFR_BOX:CFR_BOX,-CFR_BOX:CFR_BOX), sw_px(DIM_SW,DIM_SW)
        real    :: z, sigma, h, h_sq, avg, mmsk
        integer :: i, j, m, n, pad, ithr
        if( self%is_3d() ) THROW_HARD('2D images only; NLmean2D')
        if( self%ft )      THROW_HARD('Real space only;NLmean2D')
        mmsk = real(self%ldim(1)) / 2. - real(DIM_SW)
        if( present(msk) ) mmsk = msk
        if( present(sdev_noise) )then
            if( sdev_noise > SMALL )then
                sigma = sdev_noise
            else
                sigma = self%noisesdev(mmsk) ! estimation of noise
            endif
        else
            sigma = self%noisesdev(mmsk) ! estimation of noise
        endif
        pad   = CFR_BOX + 2
        h     = 4.*sigma
        h_sq  = h**2.
        avg   = sum(self%rmat(:self%ldim(1),:self%ldim(2),1)) / real(product(self%ldim))
        allocate(rmat_threads(nthr_glob,self%ldim(1),self%ldim(2),1), source=0.)
        allocate(rmat_pad(-pad:self%ldim(1)+pad,-pad:self%ldim(2)+pad), source=avg)
        rmat_pad(1:self%ldim(1),1:self%ldim(2)) = self%rmat(1:self%ldim(1),1:self%ldim(2),1)
        !$omp parallel do schedule(static) default(shared) private(m,n,ithr,sw_px,exponentials,i,j,z)&
        !$omp proc_bind(close) firstprivate(rmat_pad) collapse(2)
        do n = 1,self%ldim(2)
            do m = 1,self%ldim(1)
                ithr  = omp_get_thread_num() + 1
                sw_px = rmat_pad(m:m+DIM_SW-1,n:n+DIM_SW-1)
                exponentials = 0.
                do j = -CFR_BOX,CFR_BOX
                    do i = -CFR_BOX,CFR_BOX
                      exponentials(i,j) = &
                      & exp( -sum( (sw_px - rmat_pad(m+i:m+i+DIM_SW-1, n+j:n+j+DIM_SW-1))**2. )/h_sq) ! Euclidean norm
                  enddo
                enddo
                z = sum(exponentials)
                if( z < 0.0000001 ) cycle
                rmat_threads(ithr,m,n,1) = sum(exponentials * rmat_pad(m-CFR_BOX:m+CFR_BOX,n-CFR_BOX:n+CFR_BOX)) / z
            enddo
        enddo
        !$omp end parallel do
        self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = sum(rmat_threads, dim=1)
    end subroutine NLmean2D

    module subroutine NLmean2D_eo( even, odd, avg )
        class(image), intent(inout) :: even, odd, avg
        type(image)        :: noise, noise_var
        real,  allocatable :: rmat_pad(:,:), rmat_threads(:,:,:,:)
        integer, parameter :: DIM_SW  = 1
        integer, parameter :: CFR_BOX = 3
        real    :: exponentials(-CFR_BOX:CFR_BOX,-CFR_BOX:CFR_BOX), sw_px(DIM_SW,DIM_SW)
        real    :: z, pix_avg, sigma2t2
        integer :: i, j, m, n, pad, ithr
        if( even%is_3d() ) THROW_HARD('2D images only; NLmean2D_eo')
        if( even%ft )      THROW_HARD('Real space only;NLmean2D_eo')
        call noise%copy(even)
        call noise%subtr(odd)
        call noise%loc_var(noise_var)
        call noise_var%norm([5.,2.])
        call avg%copy(even)
        call avg%add(odd)
        call avg%mul(0.5)
        pad     = CFR_BOX + 2
        pix_avg = sum(avg%rmat(:avg%ldim(1),:avg%ldim(2),1)) / real(product(avg%ldim))
        allocate(rmat_threads(nthr_glob,avg%ldim(1),avg%ldim(2),1),   source=0.)
        allocate(rmat_pad(-pad:avg%ldim(1)+pad,-pad:avg%ldim(2)+pad), source=pix_avg)
        rmat_pad(1:avg%ldim(1),1:avg%ldim(2)) = avg%rmat(1:avg%ldim(1),1:avg%ldim(2),1)
        do n = 1,avg%ldim(2)
            do m = 1,avg%ldim(1)
                ithr         = omp_get_thread_num() + 1
                sw_px        = rmat_pad(m:m+DIM_SW-1,n:n+DIM_SW-1)
                exponentials = 0.
                sigma2t2     = 2. * noise_var%rmat(m,n,1)
                do j = -CFR_BOX,CFR_BOX
                    do i = -CFR_BOX,CFR_BOX
                      exponentials(i,j) = &
                      & exp( -sum( (sw_px - rmat_pad(m+i:m+i+DIM_SW-1, n+j:n+j+DIM_SW-1))**2. )/sigma2t2) ! Euclidean norm
                  enddo
                enddo
                z = sum(exponentials)
                if( z < 0.0000001 ) cycle
                rmat_threads(ithr,m,n,1) = sum(exponentials * rmat_pad(m-CFR_BOX:m+CFR_BOX,n-CFR_BOX:n+CFR_BOX)) / z
            enddo
        enddo
        avg%rmat(:avg%ldim(1),:avg%ldim(2),:avg%ldim(3)) = sum(rmat_threads, dim=1)
        call noise%kill
        call noise_var%kill
    end subroutine NLmean2D_eo

    module subroutine NLmean3D( self, msk, sdev_noise )
        class(image),   intent(inout) :: self
        real, optional, intent(in)    :: msk
        real, optional, intent(in)    :: sdev_noise
        real,  allocatable :: rmat_pad(:,:,:), rmat_threads(:,:,:,:)
        integer, parameter :: DIM_SW  = 3
        integer, parameter :: CFR_BOX = 10 ! as suggested in the paper, use a box 21x21
        real    :: exponentials(-CFR_BOX:CFR_BOX,-CFR_BOX:CFR_BOX,-CFR_BOX:CFR_BOX), sw_px(DIM_SW,DIM_SW,DIM_SW)
        real    :: z, sigma, h, h_sq, mmsk, avg
        integer :: i, j, k, m, n, o, pad, ithr
        if( self%is_2d() ) THROW_HARD('3D images only; NLmean3D')
        if( self%ft )      THROW_HARD('Real space only; 3DNLmean3D')
        mmsk = real(self%ldim(1)) / 2. - real(DIM_SW)
        if( present(msk) ) mmsk = msk
        if( present(sdev_noise) )then
            if( sdev_noise > SMALL )then
                sigma = sdev_noise
            else
                sigma = self%noisesdev(mmsk) ! estimation of noise
            endif
        else
            sigma = self%noisesdev(mmsk) ! estimation of noise
        endif
        pad  = CFR_BOX + 2
        h    = 4.*sigma
        h_sq = h**2.
        avg  = sum(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3))) / real(product(self%ldim))
        allocate(rmat_threads(nthr_glob,self%ldim(1),self%ldim(2),self%ldim(3)), source=0.)
        allocate(rmat_pad(-pad:self%ldim(1)+pad,-pad:self%ldim(2)+pad,-pad:self%ldim(3)+pad), source=avg)
        rmat_pad(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) = self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))
        !$omp parallel do schedule(static) default(shared) private(m,n,o,ithr,sw_px,exponentials,i,j,k,z)&
        !$omp proc_bind(close) firstprivate(rmat_pad) collapse(3)
        do o = 1, self%ldim(3)
            do n = 1, self%ldim(2)
                do m = 1, self%ldim(1)
                    ithr  = omp_get_thread_num() + 1
                    sw_px = rmat_pad(m:m+DIM_SW-1,n:n+DIM_SW-1,o:o+DIM_SW-1)
                    exponentials = 0.
                    do k = -CFR_BOX,CFR_BOX
                        do j = -CFR_BOX,CFR_BOX
                            do i = -CFR_BOX,CFR_BOX
                                exponentials(i,j,k) = &
                                & exp(-sum((sw_px - rmat_pad(m+i:m+i+DIM_SW-1,n+j:n+j+DIM_SW-1,o+k:o+k+DIM_SW-1))**2.)/h_sq) ! euclidean norm
                            enddo
                        enddo
                    enddo
                    z = sum(exponentials)
                    if( z < 0.0000001 ) cycle
                    rmat_threads(ithr,m,n,o) = sum(exponentials * rmat_pad(m-CFR_BOX:m+CFR_BOX,n-CFR_BOX:n+CFR_BOX,o-CFR_BOX:o+CFR_BOX)) / z
                enddo
            enddo
        enddo
        !$omp end parallel do
        self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) = sum(rmat_threads, dim=1)
    end subroutine NLmean3D

    module subroutine NLmean3D_eo( even, odd, avg )
        class(image), intent(inout) :: even, odd, avg
        real,  allocatable :: rmat_pad(:,:,:), rmat_threads(:,:,:,:)
        integer, parameter :: DIM_SW  = 1
        integer, parameter :: CFR_BOX = 3
        type(image) :: noise, noise_var
        real    :: exponentials(-CFR_BOX:CFR_BOX,-CFR_BOX:CFR_BOX,-CFR_BOX:CFR_BOX), sw_px(DIM_SW,DIM_SW,DIM_SW)
        real    :: z, sigma2t2, pixavg
        integer :: i, j, k, m, n, o, pad, ithr
        if( even%is_2d() ) THROW_HARD('3D images only; NLmean3D_eo')
        if( even%ft )      THROW_HARD('Real space only; 3DNLmean3D_eo')
        call noise%copy(even)
        call noise%subtr(odd)
        call noise%loc_var3D(noise_var)
        call noise_var%norm([5.,2.])
        call avg%copy(even)
        call avg%add(odd)
        call avg%mul(0.5)
        pad     = CFR_BOX + 2
        pixavg  = sum(avg%rmat(:avg%ldim(1),:avg%ldim(2),:avg%ldim(3))) / real(product(avg%ldim))
        allocate(rmat_threads(nthr_glob,avg%ldim(1),avg%ldim(2),avg%ldim(3)), source=0.)
        allocate(rmat_pad(-pad:avg%ldim(1)+pad,-pad:avg%ldim(2)+pad,-pad:avg%ldim(3)+pad), source=pixavg)
        rmat_pad(1:avg%ldim(1),1:avg%ldim(2),1:avg%ldim(3)) = avg%rmat(1:avg%ldim(1),1:avg%ldim(2),1:avg%ldim(3))
        !$omp parallel do schedule(static) default(shared) private(m,n,o,ithr,sw_px,exponentials,sigma2t2,i,j,k,z)&
        !$omp proc_bind(close) firstprivate(rmat_pad) collapse(3)
        do o = 1, avg%ldim(3)
            do n = 1, avg%ldim(2)
                do m = 1, avg%ldim(1)
                    ithr         = omp_get_thread_num() + 1
                    sw_px        = rmat_pad(m:m+DIM_SW-1,n:n+DIM_SW-1,o:o+DIM_SW-1)
                    exponentials = 0.
                    sigma2t2     = 2. * noise_var%rmat(m,n,o)
                    do k = -CFR_BOX,CFR_BOX
                        do j = -CFR_BOX,CFR_BOX
                            do i = -CFR_BOX,CFR_BOX
                                exponentials(i,j,k) = &
                                & exp(-sum((sw_px - rmat_pad(m+i:m+i+DIM_SW-1,n+j:n+j+DIM_SW-1,o+k:o+k+DIM_SW-1))**2.)/sigma2t2) ! euclidean norm
                            enddo
                        enddo
                    enddo
                    z = sum(exponentials)
                    if( z < 0.0000001 ) cycle
                    rmat_threads(ithr,m,n,o) = sum(exponentials * rmat_pad(m-CFR_BOX:m+CFR_BOX,n-CFR_BOX:n+CFR_BOX,o-CFR_BOX:o+CFR_BOX)) / z
                enddo
            enddo
        enddo
        !$omp end parallel do
        avg%rmat(:avg%ldim(1),:avg%ldim(2),:avg%ldim(3)) = sum(rmat_threads, dim=1)
    end subroutine NLmean3D_eo

    module subroutine ICM2D( self, lambda, verbose )
        class(image),      intent(inout) :: self
        real,              intent(in)    :: lambda
        logical, optional, intent(in)    :: verbose
        integer, parameter :: MAXITS    = 3
        integer, parameter :: NQUANTA   = 256
        integer     :: n_8(3,8), nsz, i, j, k, m, n
        real        :: pot_term, pix, min, proba, sigma2t2, x, xmin, transl_tab(NQUANTA), eucl, y, sy, syy, diff, rnsz
        type(image) :: self_prev
        logical     :: l_verbose
        if( self%is_3d() ) THROW_HARD('2D images only; ICM')
        if( self%ft )      THROW_HARD('Real space only; ICM')
        l_verbose = .true.
        if( present(verbose) ) l_verbose = verbose
        call self%quantize_fwd(NQUANTA, transl_tab)
        call self_prev%copy(self)
        sigma2t2 = 10.
        do i = 1, MAXITS
            do m = 1,self%ldim(2)
                do n = 1,self%ldim(1)
                    pix  = self_prev%rmat(n,m,1)
                    call neigh_8(self%ldim, [n,m,1], n_8, nsz)
                    rnsz = real(nsz)
                    sy   = 0.
                    syy  = 0.
                    do j = 1, nsz
                        y   = self_prev%rmat(n_8(1,j),n_8(2,j),1)
                        sy  = sy  + y
                        syy = syy + y*y
                    end do
                    pot_term = syy
                    min      = (pix * pix) / sigma2t2 + lambda * pot_term
                    xmin     = 0.
                    ! Every shade of gray is tested to find the a local minimum of the energy corresponding to a Gibbs distribution
                    do k = 1,NQUANTA - 1
                        x        = real(k)
                        pot_term = syy + rnsz*x*x - 2.0*sy*x
                        diff     = pix - x
                        proba    = (diff * diff) / sigma2t2 + lambda * pot_term
                        if( min > proba )then
                            min  = proba
                            xmin = x
                        endif
                    end do
                    self%rmat(n,m,1) = xmin
                end do
            end do
            if( l_verbose )then
                eucl = self%euclid_norm(self_prev)
                ! write(logfhandle,'(A,I2,A,F8.4)') 'ICM Iteration ', i, ', Euclidean distance ', eucl
            endif
            if( i < MAXITS ) self_prev%rmat = self%rmat
        end do
        call self%quantize_bwd(NQUANTA, transl_tab)
        call self_prev%kill
    end subroutine ICM2D

    module subroutine ICM2D_eo( even, odd, lambda, verbose )
        class(image),      intent(inout) :: even, odd
        real,              intent(in)    :: lambda
        logical, optional, intent(in)    :: verbose
        integer, parameter :: MAXITS    = 3
        integer, parameter :: NQUANTA   = 256
        type(image) :: even_prev, odd_prev, noise, noise_var
        integer     :: n_8(3,8), nsz, i, j, k, m, n
        real        :: transl_tab_even(NQUANTA), transl_tab_odd(NQUANTA)
        real        :: pot_term(2), pix(2), minv(2), proba(2), sigma2t2
        real        :: x, xmin(2), eucl, y(2), sy(2), syy(2), diff(2), rnsz
        logical     :: l_verbose
        if( even%is_3d() ) THROW_HARD('2D images only; ICM2D_eo')
        if( even%ft )      THROW_HARD('Real space only; ICM2D_eo')
        l_verbose = .true.
        if( present(verbose) ) l_verbose = verbose
        call noise%copy(even)
        call noise%subtr(odd)
        call noise%loc_var(noise_var)
        call noise_var%norm([5.,2.])
        call even%quantize_fwd(NQUANTA, transl_tab_even)
        call odd%quantize_fwd(NQUANTA, transl_tab_odd)
        call even_prev%copy(even)
        call odd_prev%copy(odd)
        do i = 1, MAXITS
            do m = 1,even%ldim(2)
                do n = 1,even%ldim(1)
                    pix(1)   = odd_prev%rmat(n,m,1)
                    pix(2)   = even_prev%rmat(n,m,1)
                    sigma2t2 = 2. * noise_var%rmat(n,m,1)
                    call neigh_8(even%ldim, [n,m,1], n_8, nsz)
                    rnsz = real(nsz)
                    sy   = 0.
                    syy  = 0.
                    do j = 1, nsz
                        y(1) = odd_prev%rmat(n_8(1,j),n_8(2,j),1)
                        y(2) = even_prev%rmat(n_8(1,j),n_8(2,j),1)
                        sy   = sy  + y
                        syy  = syy + y*y
                    end do
                    xmin     = 0.
                    pot_term = syy ! x=0.
                    minv     = (pix * pix) / sigma2t2 + lambda * pot_term
                    ! Every shade of gray is tested to find the a local minimum of the energy corresponding to a Gibbs distribution
                    do k = 1,NQUANTA - 1
                        x        = real(k)
                        pot_term = syy + rnsz*x*x - 2.0*sy*x
                        diff     = pix - x
                        proba    = (diff * diff) / sigma2t2 + lambda * pot_term
                        if( minv(1) > proba(1) )then
                            minv(1) = proba(1)
                            xmin(1) = x
                        endif
                        if( minv(2) > proba(2) )then
                            minv(2) = proba(2)
                            xmin(2) = x
                        endif
                    end do
                    odd%rmat(n,m,1)  = xmin(1)
                    even%rmat(n,m,1) = xmin(2)
                end do
            end do
            if( l_verbose )then
                eucl = (even%euclid_norm(even_prev) + odd%euclid_norm(odd_prev)) / 2.
                ! write(logfhandle,'(A,I2,A,F8.4)') 'ICM Iteration ', i, ', Euclidean distance ', eucl
            endif
            if( i < MAXITS)then
                even_prev%rmat = even%rmat
                odd_prev%rmat  = odd%rmat
            endif
        end do
        call even%quantize_bwd(NQUANTA, transl_tab_even)
        call odd%quantize_bwd(NQUANTA, transl_tab_odd)
        call even_prev%kill
        call odd_prev%kill
        call noise%kill
        call noise_var%kill
    end subroutine ICM2D_eo

    module subroutine ICM3D( self, lambda )
        class(image), intent(inout) :: self
        real,         intent(in)    :: lambda
        integer, parameter :: MAXITS    = 3
        integer, parameter :: NQUANTA   = 256
        type(image) :: self_prev
        integer     :: n_4(3,6), nsz, i, j, k, m, n, l
        real        :: pot_term, pix, min, proba, sigma2t2, x, xmin, transl_tab(NQUANTA), eucl,y, sy, syy, diff, rnsz
        if( self%is_2d() ) THROW_HARD('3D images only; ICM')
        if( self%ft )      THROW_HARD('Real space only; ICM')
        call self%quantize_fwd(NQUANTA, transl_tab)
        call self_prev%copy(self)
        sigma2t2 = 10.
        do i = 1, MAXITS
            !$omp parallel do schedule(static) default(shared) private(n,m,l,pix,n_4,nsz,rnsz,pot_term,diff,j,xmin,min,k,x,proba,y,sy,syy)&
            !$omp proc_bind(close) collapse(3)
            do l = 1,self%ldim(3)
                do m = 1,self%ldim(2)
                    do n = 1,self%ldim(1)
                        pix  = self_prev%rmat(n,m,l)
                        call neigh_4_3D(self%ldim, [n,m,l], n_4, nsz)
                        rnsz = real(nsz)
                        sy   = 0.
                        syy  = 0.
                        do j = 1, nsz
                            y   = self_prev%rmat(n_4(1,j),n_4(2,j),n_4(3,j))
                            sy  = sy  + y
                            syy = syy + y*y
                        end do
                        pot_term = syy
                        min      = (pix * pix) / sigma2t2 + lambda * pot_term
                        xmin     = 0.
                        ! Every shade of gray is tested to find the a local minimum of the energy corresponding to a Gibbs distribution
                        do k = 1,NQUANTA - 1
                            x        = real(k)
                            pot_term = syy + rnsz*x*x - 2.0*sy*x
                            diff     = pix - x
                            proba    = (diff * diff) / sigma2t2 + lambda * pot_term
                            if( min > proba )then
                                min  = proba
                                xmin = x
                            endif
                        end do
                        self%rmat(n,m,l) = xmin
                    end do
                end do
            end do
            !$omp end parallel do
            eucl = self%euclid_norm(self_prev)
            ! write(logfhandle,'(A,I2,A,F8.4)') 'ICM Iteration ', i, ', Euclidean distance ', eucl 
            call self_prev%copy(self)
        end do
        call self%quantize_bwd(NQUANTA, transl_tab)
        call self_prev%kill
    end subroutine ICM3D

    module subroutine ICM3D_eo( even, odd, lambda, l_msk )
        class(image),      intent(inout) :: even, odd
        real,              intent(in)    :: lambda
        logical, optional, intent(in)    :: l_msk(even%ldim(1),even%ldim(2),even%ldim(3))
        integer, parameter :: MAXITS    = 3
        integer, parameter :: NQUANTA   = 256
        type(image) :: even_prev, odd_prev, noise, noise_var
        integer     :: n_4(3,6), nsz, i, j, k, m, n, l
        real        :: transl_tab_even(NQUANTA), transl_tab_odd(NQUANTA)
        real        :: sy(2), syy(2), y(2), pot_term(2), pix(2), minv(2)
        real        :: proba(2), sigma2t2, x, xmin(2), eucl, diff(2), rnsz
        if( even%is_2d() ) THROW_HARD('3D images only; ICM3D_eo')
        if( even%ft      ) THROW_HARD('Real space only; ICM3D_eo')
        call noise%copy(even)
        call noise%subtr(odd)
        call noise%loc_var3D(noise_var)
        call noise_var%norm([5.,2.])
        if( present(l_msk) )then
            call even%quantize_fwd(NQUANTA, transl_tab_even, l_msk)
            call odd%quantize_fwd(NQUANTA, transl_tab_odd, l_msk)
        else
            call even%quantize_fwd(NQUANTA, transl_tab_even)
            call odd%quantize_fwd(NQUANTA, transl_tab_odd)
        endif
        call even_prev%copy(even)
        call odd_prev%copy(odd)
        do i = 1, MAXITS
            !$omp parallel do private(n,m,l,pix,sigma2t2,n_4,nsz,rnsz,pot_term,diff,j,xmin,minv,k,x,proba,y,syy,sy)&
            !$omp proc_bind(close) collapse(3) schedule(static) default(shared)
            do l = 1,even%ldim(3)
                do m = 1,even%ldim(2)
                    do n = 1,even%ldim(1)
                        pix(1)   = odd_prev%rmat(n,m,l)
                        pix(2)   = even_prev%rmat(n,m,l)
                        sigma2t2 = 2. * noise_var%rmat(n,m,l)
                        call neigh_4_3D(even%ldim, [n,m,l], n_4, nsz)
                        rnsz = real(nsz)
                        ! x: central pixel/candidate value
                        ! y: nsz neighbours, constants
                        ! pot_term = SUMi((x-yi)**2) = SUMi(yi**2) + nsz.x**2 - 2.x.SUMi(yi)
                        sy  = 0.
                        syy = 0.
                        do j = 1, nsz
                            y(1) = odd_prev%rmat( n_4(1,j),n_4(2,j),n_4(3,j))
                            y(2) = even_prev%rmat(n_4(1,j),n_4(2,j),n_4(3,j))
                            sy   = sy  + y
                            syy  = syy + y * y
                        end do
                        xmin     = 0.
                        pot_term = syy ! x=0
                        minv     = (pix * pix) / sigma2t2 + lambda * pot_term
                        ! Every shade of gray is tested to find the a local minimum of the energy corresponding to a Gibbs distribution
                        do k = 1,NQUANTA - 1
                            x        = real(k)
                            pot_term = syy + rnsz*x*x - 2.0*sy*x
                            diff     = pix - x
                            proba    = (diff * diff) / sigma2t2 + lambda * pot_term
                            if( minv(1) > proba(1) )then
                                minv(1) = proba(1)
                                xmin(1) = x
                            endif
                            if( minv(2) > proba(2) )then
                                minv(2) = proba(2)
                                xmin(2) = x
                            endif
                        end do
                        odd%rmat(n,m,l)  = xmin(1)
                        even%rmat(n,m,l) = xmin(2)
                    end do
                end do
            end do
            !$omp end parallel do
            eucl = (even%euclid_norm(even_prev) + odd%euclid_norm(odd_prev)) / 2.
            ! write(logfhandle,'(A,I2,A,F8.4)') 'ICM Iteration ', i, ', Euclidean distance ', eucl
            call even_prev%copy(even)
            call odd_prev%copy(odd)
        end do
        call even%quantize_bwd(NQUANTA, transl_tab_even)
        call odd%quantize_bwd(NQUANTA, transl_tab_odd)
        call even_prev%kill
        call odd_prev%kill
        call noise%kill
        call noise_var%kill
    end subroutine ICM3D_eo

    module subroutine GLCM( self, nquanta, pmat )
        class(image), intent(in)    :: self
        integer,      intent(in)    :: nquanta
        real,         intent(inout) :: pmat(nquanta,nquanta)
        type(image) :: img_q
        real        :: transl_tab(nquanta)
        integer     :: m, n, j, n_8(3,8), nsz, ipix, jpix
        call img_q%copy(self)
        call img_q%quantize_fwd(nquanta, transl_tab)
        do m = 1,img_q%ldim(2)
            do n = 1,img_q%ldim(1)
                ipix = nint(img_q%rmat(n,m,1) + 1.)
                call neigh_8(img_q%ldim, [n,m,1], n_8, nsz)
                do j = 1, nsz
                    jpix = nint(img_q%rmat(n_8(1,j),n_8(2,j),1) + 1.)
                    if( jpix == ipix ) pmat(ipix,jpix) = pmat(ipix,jpix) + 1.
                end do
            end do
        end do
        call normalize_minmax(pmat)
        call img_q%kill
    end subroutine GLCM

end submodule simple_image_filt
