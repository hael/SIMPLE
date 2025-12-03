submodule (simple_image) simple_image_freq_anal
!$ use omp_lib
!$ use omp_lib_kinds
include  'simple_lib.f08'
#include "simple_local_flags.inc"
implicit none
contains

    module subroutine acf( self )
        class(image), intent(inout) :: self
        if( .not. self%ft )then
            call self%fft()
        endif
        self%cmat = self%cmat*conjg(self%cmat) / sqrt(sum(csq(self%cmat)))
        call self%ifft()
    end subroutine acf

    module function ccf( self1, self2 ) result( cc )
        class(image), intent(inout) :: self1, self2
        type(image) :: cc
        if( .not. self1%ft )then
            call self1%fft()
        endif
        if( .not. self2%ft )then
            call self2%fft()
        endif
        cc      = self1
        cc%cmat = cc%cmat*conjg(self2%cmat)
        call cc%ifft()
    end function ccf

    ! keep serial
    module subroutine spectrum( self, which, spec, norm )
        class(image),      intent(inout) :: self
        character(len=*),  intent(in)    :: which
        real, allocatable, intent(inout) :: spec(:)
        logical, optional, intent(in)    :: norm
        real, allocatable :: counts(:)
        integer :: lims(3,2), phys(3), sh, lfny, h, k, l
        logical :: didft, nnorm
        nnorm = .true.
        if( present(norm) ) nnorm = norm
        didft = .false.
        if( which .ne. 'count' )then
            if( .not. self%ft )then
                call self%fft()
                didft = .true.
            endif
        endif
        lfny = self%get_lfny(1)
        if(allocated(spec))deallocate(spec)
        allocate( spec(lfny), counts(lfny) )
        spec   = 0.
        counts = 0.
        lims   = self%fit%loop_lims(2)
        select case(which)
        case('real')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        sh = nint(hyp(h,k,l))
                        if( sh == 0 .or. sh > lfny ) cycle
                        phys       = self%fit%comp_addr_phys([h,k,l])
                        spec(sh)   = spec(sh) + real(self%cmat(phys(1),phys(2),phys(3)))
                        counts(sh) = counts(sh) + 1.
                    end do
                end do
            end do
        case('power')
            do l=lims(3,1),lims(3,2)
                do k=lims(2,1),lims(2,2)
                    do h=lims(1,1),lims(1,2)
                        sh = nint(hyp(h,k,l))
                        if( sh == 0 .or. sh > lfny ) cycle
                        phys       = self%fit%comp_addr_phys([h,k,l])
                        spec(sh)   = spec(sh) + csq(self%cmat(phys(1),phys(2),phys(3)))
                        counts(sh) = counts(sh) + 1.
                    end do
                end do
            end do
        case('sqrt')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        sh = nint(hyp(h,k,l))
                        if( sh == 0 .or. sh > lfny ) cycle
                        phys       = self%fit%comp_addr_phys([h,k,l])
                        spec(sh)   = spec(sh) + sqrt(csq(self%cmat(phys(1),phys(2),phys(3))))
                        counts(sh) = counts(sh) + 1.
                    end do
                end do
            end do
        case('log')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        sh = nint(hyp(h,k,l))
                        if( sh == 0 .or. sh > lfny ) cycle
                        phys       = self%fit%comp_addr_phys([h,k,l])
                        spec(sh)   = spec(sh) + log(csq(self%cmat(phys(1),phys(2),phys(3))))
                        counts(sh) = counts(sh) + 1.
                    end do
                end do
            end do
        case('absreal')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        sh = nint(hyp(h,k,l))
                        if( sh == 0 .or. sh > lfny ) cycle
                        phys       = self%fit%comp_addr_phys([h,k,l])
                        spec(sh)   = spec(sh) + abs(real(self%cmat(phys(1),phys(2),phys(3))))
                        counts(sh) = counts(sh) + 1.
                    end do
                end do
            end do
        case('absimag')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        sh = nint(hyp(h,k,l))
                        if( sh == 0 .or. sh > lfny ) cycle
                        phys       = self%fit%comp_addr_phys([h,k,l])
                        spec(sh)   = spec(sh) + abs(aimag(self%cmat(phys(1),phys(2),phys(3))))
                        counts(sh) = counts(sh) + 1.
                    end do
                end do
            end do
        case('abs')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        sh = nint(hyp(h,k,l))
                        if( sh == 0 .or. sh > lfny ) cycle
                        phys       = self%fit%comp_addr_phys([h,k,l])
                        spec(sh)   = spec(sh) + cabs(self%cmat(phys(1),phys(2),phys(3)))
                        counts(sh) = counts(sh) + 1.
                    end do
                end do
            end do
        case('phase')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        sh = nint(hyp(h,k,l))
                        if( sh == 0 .or. sh > lfny ) cycle
                        phys       = self%fit%comp_addr_phys([h,k,l])
                        spec(sh)   = spec(sh) + phase_angle(self%cmat(phys(1),phys(2),phys(3)))
                        counts(sh) = counts(sh) + 1.
                    end do
                end do
            end do
        case('count')
            do h=lims(1,1),lims(1,2)
                do k=lims(2,1),lims(2,2)
                    do l=lims(3,1),lims(3,2)
                        sh = nint(hyp(h,k,l))
                        if( sh == 0 .or. sh > lfny ) cycle
                        phys       = self%fit%comp_addr_phys([h,k,l])
                        spec(sh)   = spec(sh) + 1.
                        counts(sh) = counts(sh)+1.
                    end do
                end do
            end do
        case DEFAULT
            write(logfhandle,*) 'Spectrum kind: ', trim(which)
            THROW_HARD('Unsupported spectrum kind; spectrum')
        end select
        if( which .ne. 'count' .and. nnorm )then
            where(counts > 0.)
                spec = spec/counts
            end where
        endif
        if( didft ) call self%ifft()
    end subroutine spectrum

    module subroutine power_spectrum( self, spec )
        class(image), intent(in)    :: self
        real,         intent(inout) :: spec(fdim(self%ldim(1)) - 1)
        integer  :: counts(fdim(self%ldim(1)) - 1)
        real(dp) :: dspec( fdim(self%ldim(1)) - 1)
        integer  :: lims(3,2), phys(3), filtsz, h, k, l, sh
        if( .not.self%is_ft() ) THROW_HARD('Only for FTed images! power_spectrum')
        filtsz = fdim(self%ldim(1)) - 1
        dspec  = 0.d0
        counts = 0
        lims   = self%fit%loop_lims(2)
        if( self%is_2d() )then
            do k=lims(2,1),lims(2,2)
                do h=lims(1,1),lims(1,2)
                    sh = nint(hyp(h,k))
                    if( sh == 0 .or. sh > filtsz ) cycle
                    phys(1:2)  = self%fit%comp_addr_phys(h,k)
                    dspec(sh)  = dspec(sh)  + real(csq_fast(self%cmat(phys(1),phys(2),1)),dp)
                    counts(sh) = counts(sh) + 1
                end do
            end do
        else
            do l=lims(3,1),lims(3,2)
                do k=lims(2,1),lims(2,2)
                    do h=lims(1,1),lims(1,2)
                        sh = nint(hyp(h,k,l))
                        if( sh == 0 .or. sh > filtsz ) cycle
                        phys       = self%fit%comp_addr_phys([h,k,l])
                        dspec(sh)  = dspec(sh)  + real(csq_fast(self%cmat(phys(1),phys(2),phys(3))),dp)
                        counts(sh) = counts(sh) + 1
                    end do
                end do
            end do
        endif
        where(counts > 0)
            dspec = dspec / real(counts,dp)
        end where
        spec = real(dspec,kind=sp)
    end subroutine power_spectrum

    module function guinier_bfac( self, hp, lp ) result( bfac )
        class(image), intent(inout) :: self
        real,         intent(in)    :: hp, lp
        real, allocatable :: plot(:,:)
        integer :: fromk, tok, nk
        real    :: slope, intercept, corr, bfac
        plot  = self%guinier(.false.)
        fromk = self%get_find(hp)
        tok   = self%get_find(lp)
        nk    = tok-fromk+1
        call fit_straight_line(nk, plot(fromk:tok,:), slope, intercept, corr)
        bfac  = 4. * slope
        deallocate(plot)
    end function guinier_bfac

    module function guinier( self, verbose ) result( plot )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: verbose
        real, allocatable :: spec(:), plot(:,:)
        integer           :: lfny, k
        call self%spectrum('absreal',spec=spec)
        lfny = self%get_lfny(1)
        allocate( plot(lfny,2) )
        do k=1,lfny
            plot(k,1) = 1./(self%get_lp(k)**2.)
            plot(k,2) = log(spec(k))
            if( verbose ) write(logfhandle,'(A,1X,F8.4,1X,A,1X,F7.3)') '>>> RECIPROCAL SQUARE RES:', plot(k,1), '>>> LOG(ABS(REAL(F))):', plot(k,2)
        end do
        deallocate( spec )
    end function guinier

    module subroutine fsc( self1, self2, corrs )
        class(image), intent(inout) :: self1, self2
        real,         intent(out)   :: corrs(fdim(self1%ldim(1))-1)
        real(dp)    :: corrs_8(fdim(self1%ldim(1))-1)
        real(dp)    :: sumasq(fdim(self1%ldim(1))-1), sumbsq(fdim(self1%ldim(1))-1)
        complex(dp) :: comp1, comp2
        integer     :: n, lims(3,2), phys(3), sh, h, k, l
        corrs_8 = 0.d0
        sumasq  = 0.d0
        sumbsq  = 0.d0
        lims    = self1%fit%loop_lims(2)
        n       = self1%get_filtsz()
        !$omp parallel do collapse(3) default(shared) private(h,k,l,phys,sh,comp1,comp2)&
        !$omp schedule(static) reduction(+:corrs_8,sumasq,sumbsq) proc_bind(close)
        do k=lims(2,1),lims(2,2)
            do h=lims(1,1),lims(1,2)
                do l=lims(3,1),lims(3,2)
                    ! compute physical address
                    phys = self1%fit%comp_addr_phys(h,k,l)
                    ! find shell
                    sh = nint(hyp(h,k,l))
                    if( sh == 0 .or. sh > n ) cycle
                    ! real part of the complex mult btw self1 and targ*
                    comp1 = self1%cmat(phys(1),phys(2),phys(3))
                    comp2 = self2%cmat(phys(1),phys(2),phys(3))
                    corrs_8(sh) = corrs_8(sh)+ real(comp1 * conjg(comp2), kind=dp)
                    sumasq(sh) = sumasq(sh) + csq(comp1)
                    sumbsq(sh) = sumbsq(sh) + csq(comp2)
                end do
            end do
        end do
        !$omp end parallel do
        ! normalize correlations and compute resolutions
        do k=1,n
            if( sumasq(k) > 0.d0 .and. sumbsq(k) > 0.d0 )then
                corrs_8(k) = corrs_8(k)/sqrt(sumasq(k) * sumbsq(k))
            else
                corrs_8(k) = 0.
            endif
        end do
        ! return single-precision corrs
        corrs = real(corrs_8)
    end subroutine fsc

    module subroutine fsc_scaled( self1, self2, sz, corrs )
        class(image), intent(inout) :: self1, self2
        integer,      intent(in)    :: sz
        real,         intent(out)   :: corrs(sz)
        real(dp)    :: corrs_8(sz), sumasq(sz), sumbsq(sz)
        real        :: scale, rsh
        complex(dp) :: comp1, comp2
        integer     :: lims(3,2), phys(3), sh, sh_sc, h, k, l, n
        scale   = real(2*sz) / real(self1%ldim(1))
        corrs_8 = 0.d0
        sumasq  = 0.d0
        sumbsq  = 0.d0
        lims    = self1%fit%loop_lims(2)
        n       = self1%get_filtsz()
        !$omp parallel do collapse(3) default(shared) private(h,k,l,phys,sh,sh_sc,rsh,comp1,comp2)&
        !$omp schedule(static) reduction(+:corrs_8,sumasq,sumbsq) proc_bind(close)
        do k=lims(2,1),lims(2,2)
            do h=lims(1,1),lims(1,2)
                do l=lims(3,1),lims(3,2)
                    ! shell
                    rsh = sqrt(real(h*h) + real(k*k) + real(l*l))
                    sh  = nint(rsh)
                    if( sh == 0 .or. sh > n ) cycle
                    sh_sc = nint(scale*rsh)
                    if( sh_sc == 0 .or. sh_sc > sz ) cycle
                    ! compute physical address
                    phys = self1%fit%comp_addr_phys(h,k,l)
                    ! real part of the complex mult btw self1 and targ*
                    comp1          = self1%cmat(phys(1),phys(2),phys(3))
                    comp2          = self2%cmat(phys(1),phys(2),phys(3))
                    corrs_8(sh_sc) = corrs_8(sh_sc)+ real(comp1 * conjg(comp2), kind=dp)
                    sumasq(sh_sc)  = sumasq(sh_sc) + csq_fast(comp1)
                    sumbsq(sh_sc)  = sumbsq(sh_sc) + csq_fast(comp2)
                end do
            end do
        end do
        !$omp end parallel do
        ! normalize correlations and compute resolutions
        where( sumasq>DTINY.and. sumbsq>DTINY )
            corrs = real(corrs_8/dsqrt(sumasq * sumbsq))
        else where
            corrs = 0.
        end where
    end subroutine fsc_scaled

    module function get_res( self ) result( res )
        class(image), intent(in) :: self
        real, allocatable        :: res(:)
        integer                  :: n, k
        n = self%get_filtsz()
        allocate( res(n) )
        do k=1,n
            res(k) = self%fit%get_lp(1,k)
        end do
    end function get_res

    module subroutine frc_pspec( self1, self2, corrs )
        class(image), intent(inout) :: self1, self2
        real,         intent(out)   :: corrs(fdim(self1%ldim(1))-1)
        integer     :: k, npix
        type(image) :: maskimg
        logical, allocatable :: l_mask(:,:,:)
        corrs = 0.
        do k = 1, fdim(self1%ldim(1))-3
            call maskimg%ring(self1%ldim, self1%smpd, real(k+2), real(k-2), npix )
            l_mask = bin2logical(maskimg)
            corrs(k) = self1%real_corr(self2, l_mask)
        end do
        call maskimg%kill
        deallocate(l_mask)
    end subroutine frc_pspec

    module subroutine fcomps_below_noise_power_stats( self, noise_vol )
        class(image), intent(inout) :: self, noise_vol
        real, allocatable :: res(:)
        integer :: h, k, l, lims(3,2), phys(3), cnt, counts(self%get_filtsz())
        integer :: sh, counts_all(self%get_filtsz()), filtsz
        real    :: sig_pow, noise_pow
        lims  = self%fit%loop_lims(2)
        if( .not.self%is_ft() ) THROW_HARD('Instance need to be FTed')
        if( .not.self%is_ft() ) THROW_HARD('noise_vol need to be FTed')
        cnt        = 0
        counts     = 0
        counts_all = 0
        filtsz     = self%get_filtsz()
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    sh             = nint(hyp(h,k,l))
                    if( sh == 0 ) cycle
                    phys           = self%comp_addr_phys([h,k,l])
                    sig_pow        = csq(self%cmat(phys(1),phys(2),phys(3)))
                    noise_pow      = csq(noise_vol%cmat(phys(1),phys(2),phys(3)))
                    if( sh <= filtsz ) counts_all(sh) = counts_all(sh) + 1
                    if( noise_pow > sig_pow )then
                        call self%mul([h,k,l], 0.)
                        cnt = cnt + 1
                    else
                        if( sh <= filtsz ) counts(sh) = counts(sh) + 1
                    endif
                enddo
            enddo
        enddo
        res = self%get_res()
        do k=1,size(res)
            write(logfhandle,'(A,1X,F6.2,1X,A,1X,F15.3)') '>>> RESOLUTION:',&
            &res(k), '>>> % pow(FCOMPS) > pow(NOISE):', 100 * (real(counts(k)) / real(counts_all(k)))
        end do
        deallocate(res)
        write(logfhandle,'(a,f4.1)')&
        'fraction (%) of Fourier components with power below noise power zeroed: ',&
        &100 * (real(cnt) / product(self%ldim))
    end subroutine fcomps_below_noise_power_stats

    module subroutine resmsk( self, hplim, lplim )
        class(image), intent(inout) :: self
        real,         intent(in)    :: hplim, lplim
        integer :: h, k, lims(3,2), mh, mk, inds(3)
        real    :: freq, hplim_freq, lplim_freq
        hplim_freq = self%fit%get_find(1,hplim)
        lplim_freq = self%fit%get_find(1,lplim)
        lims = self%loop_lims(3)
        mh = abs(lims(1,1))
        mk = abs(lims(2,1))
        inds = 1
        self%rmat = 0.0
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                freq = hyp(h,k)
                if(freq >= hplim_freq .and. freq <= lplim_freq )then
                    inds(1) = min(max(1,h+mh+1),self%ldim(1))
                    inds(2) = min(max(1,k+mk+1),self%ldim(2))
                    call self%set(inds, 1.)
                endif
            end do
        end do
    end subroutine resmsk

    module subroutine img2spec( self, speckind, lp_backgr_subtr, img_out, postproc )
        class(image),      intent(inout) :: self
        character(len=*),  intent(in)    :: speckind
        real,              intent(in)    :: lp_backgr_subtr
        type(image),       intent(inout) :: img_out
        logical, optional, intent(in)    :: postproc
        type(image) :: tmp, tmp2
        logical     :: didft, l_postproc
        if( self%ldim(3) /= 1 ) THROW_HARD('only for 2D images')
        l_postproc = .true.
        if( present(postproc) ) l_postproc = postproc
        didft = .false.
        if( self%ft )then
            call self%ifft()
            didft = .true.
        endif
        call img_out%new(self%ldim, self%smpd)
        call tmp%copy(self)
        call tmp%norm()
        call tmp%zero_edgeavg
        call tmp%fft()
        call tmp%ft2img(speckind, img_out)
        if( l_postproc )then
            call img_out%dampen_pspec_central_cross
            call img_out%subtr_backgr(lp_backgr_subtr)
        endif
        if( didft ) call self%fft()
        call tmp%kill
        call tmp2%kill
    end subroutine img2spec

    module subroutine mic2spec( self, box, speckind, lp_backgr_subtr, img_out, postproc )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: box
        character(len=*),  intent(in)    :: speckind
        real,              intent(in)    :: lp_backgr_subtr
        type(image),       intent(inout) :: img_out
        logical, optional, intent(in)    :: postproc
        type(image) :: tmp, tmp2
        integer     :: xind, yind, cnt
        logical     :: didft, outside, l_postproc
        if( self%ldim(3) /= 1 ) THROW_HARD('only for 2D images')
        if( self%ldim(1) <= box .or. self%ldim(2) <= box )then
            THROW_HARD('cannot use a box larger than the image')
        endif
        l_postproc = .true.
        if( present(postproc) ) l_postproc = postproc
        didft = .false.
        if( self%ft )then
            call self%ifft()
            didft = .true.
        endif
        call img_out%new([box,box,1], self%smpd)
        call tmp%new([box,box,1], self%smpd)
        call tmp2%new([box,box,1], self%smpd)
        cnt = 0
        do xind=0,self%ldim(1)-box,box/2
            do yind=0,self%ldim(2)-box,box/2
                call self%window_slim([xind,yind],box,tmp,outside)
                call tmp%norm()
                call tmp%zero_edgeavg
                call tmp%fft()
                call tmp%ft2img(speckind, tmp2)
                call img_out%add(tmp2)
                cnt = cnt+1
                call tmp%zero_and_unflag_ft
                call tmp2%zero_and_unflag_ft
            end do
        end do
        call img_out%div(real(cnt))
        if( l_postproc )then
            call img_out%dampen_pspec_central_cross
            call img_out%subtr_backgr(lp_backgr_subtr)
        endif
        if( didft ) call self%fft()
        call tmp%kill
        call tmp2%kill
    end subroutine mic2spec

    module subroutine pspec_graphene_mask( self, ldim, smpd )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: ldim(3)
        real,         intent(in)    :: smpd
        logical, allocatable :: graphene_mask(:)
        type(image) :: tmp
        integer     :: h, k, l, lims(3,2), phys(3), sh, lfny
        call self%new(ldim, smpd)
        self%ft = .true.
        call tmp%new(ldim, smpd)
        graphene_mask = calc_graphene_mask(ldim(1), self%smpd)
        lims = self%fit%loop_lims(2)
        lfny = self%get_lfny(1)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    sh = nint(hyp(h,k,l))
                    if( sh == 0 .or. sh > lfny ) cycle
                    phys = self%fit%comp_addr_phys([h,k,l])
                    if( graphene_mask(sh) )then
                        self%cmat(phys(1),phys(2),phys(3)) = cmplx(1.,0.)
                    else
                        self%cmat(phys(1),phys(2),phys(3)) = cmplx(0.,0.)
                    endif
                end do
            end do
        end do
        call self%ft2img('real', tmp)
        call self%copy(tmp)
        call tmp%kill
    end subroutine pspec_graphene_mask

    !> \brief dampens the central cross of a powerspectrum by mean filtering
    module subroutine dampen_pspec_central_cross( self )
        class(image), intent(inout) :: self
        integer, parameter :: HDAMPWINSZ=2
        integer, parameter :: DAMPWINSZ =5
        real :: medi(self%ldim(1)),medj(self%ldim(2)),vals(DAMPWINSZ**2)
        integer :: h,mh,k,mk,lims(3,2),i,j,l,r,u,d,n
        if( self%ft )          THROW_HARD('not intended for FTs; dampen_pspec_central_cross')
        if( self%ldim(3) > 1 ) THROW_HARD('not intended for 3D imgs; dampen_pspec_central_cross')
        lims = self%loop_lims(3)
        mh   = abs(lims(1,1))
        mk   = abs(lims(2,1))
        n    = DAMPWINSZ*DAMPWINSZ
        ! along h
        h = 0
        i = min(max(1,h+mh+1),self%ldim(1))
        l = min(max(1,h-HDAMPWINSZ+mh+1),self%ldim(1))
        r = min(max(1,h+HDAMPWINSZ+mh+1),self%ldim(1))
        do j=1,self%ldim(2)
            d = max(1,j-HDAMPWINSZ)
            u = min(d+DAMPWINSZ-1,self%ldim(2))
            d = u-DAMPWINSZ+1
            vals = reshape(self%rmat(l:r,d:u,1),(/n/))
            medj(j) = median_nocopy(vals)
        enddo
        ! along k
        k = 0
        j = min(max(1,k+mk+1),self%ldim(2))
        d = min(max(1,k-HDAMPWINSZ+mk+1),self%ldim(2))
        u = min(max(1,k+HDAMPWINSZ+mk+1),self%ldim(2))
        do i=1,self%ldim(1)
            l = max(1,i-HDAMPWINSZ)
            r = min(l+DAMPWINSZ-1,self%ldim(1))
            l = r-DAMPWINSZ+1
            vals = reshape(self%rmat(l:r,d:u,1),(/n/))
            medi(i) = median_nocopy(vals)
        enddo
        ! replace
        h = 0
        i = min(max(1,h+mh+1),self%ldim(1))
        self%rmat(i,1:self%ldim(2),1) = medj
        k = 0
        j = min(max(1,k+mk+1),self%ldim(2))
        self%rmat(1:self%ldim(1),j,1) = medi
    end subroutine dampen_pspec_central_cross

end submodule simple_image_freq_anal