submodule (simple_image) simple_image_seg
!$ use omp_lib
!$ use omp_lib_kinds
include  'simple_lib.f08'
#include "simple_local_flags.inc"
implicit none
contains

    module function nforeground( self ) result( n )
        class(image), intent(in) :: self
        integer :: n
        n = count(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) > 0.5)
    end function nforeground

    module function nbackground( self ) result( n )
        class(image), intent(in) :: self
        integer :: n
        n = product(self%ldim)-self%nforeground()
    end function nbackground

    module subroutine binarize_1( self_in, thres, self_out )
        class(image),           intent(inout) :: self_in
        real,                   intent(in)    :: thres
        class(image), optional, intent(inout) :: self_out
        integer :: n_foreground
        n_foreground = count(self_in%rmat > thres)
        if( n_foreground < 1 ) THROW_HARD('Binarization produces empty image!')
        if( self_in%ft ) THROW_HARD('only for real images; bin_1')
        if( present(self_out) )then
            if( any(self_in%ldim /= self_out%ldim)) THROW_HARD('Images dimensions are not compatible; binarize_1')
            where( self_in%rmat >= thres )
                self_out%rmat  = 1.
            elsewhere
                self_out%rmat  = 0.
            end where
        else
            where( self_in%rmat >= thres )
                self_in%rmat = 1.
            elsewhere
                self_in%rmat = 0.
            end where
        endif
    end subroutine binarize_1

    module subroutine binarize_2( self, npix )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: npix
        real, allocatable :: forsort(:)
        real    :: thres
        integer :: npixtot
        if( self%ft ) THROW_HARD('only for real images')
        npixtot = product(self%ldim)
        forsort = pack( self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), .true.)
        call hpsort(forsort)
        thres = forsort(npixtot-npix-1) ! everyting above this value 1 else 0
        call self%binarize_1( thres )
        deallocate( forsort )
    end subroutine binarize_2

    module subroutine binarize_3( self, thres, mask )
        class(image), intent(inout) :: self
        real,         intent(in)    :: thres
        logical,      intent(inout) :: mask(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))
        if( self%ft ) THROW_HARD('only for real images')
        where( self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) >= thres )
            mask = .true.
        elsewhere
            mask = .false.
        end where
    end subroutine binarize_3

    module subroutine cendist( self, c_point )
        class(image),   intent(inout) :: self
        real, optional, intent(in)    :: c_point(3)
        real    :: centre(3)
        integer :: i
        if( self%ft ) THROW_HARD('real space only; cendist')
        ! Builds square distance image
        self   = 0.
        if( present(c_point) )then
            centre = c_point
        else
            centre = real(self%ldim)/2.+1.
        endif
        if( self%is_2d() )then
            ! 2D
            do i=1,self%ldim(1)
                self%rmat(i,:,1) = self%rmat(i,:,1) + (real(i)-centre(1))**2.
            enddo
            do i=1,self%ldim(2)
                self%rmat(:,i,1) = self%rmat(:,i,1) + (real(i)-centre(2))**2.
            enddo
        else
            ! 3D
            do i=1,self%ldim(1)
                self%rmat(i,:,:) = self%rmat(i,:,:) + (real(i)-centre(1))**2.
            enddo
            do i=1,self%ldim(2)
                self%rmat(:,i,:) = self%rmat(:,i,:) + (real(i)-centre(2))**2.
            enddo
            do i=1,self%ldim(3)
                self%rmat(:,:,i) = self%rmat(:,:,i) + (real(i)-centre(3))**2.
            enddo
        endif
        self%rmat = sqrt(self%rmat)
    end subroutine cendist

    module subroutine bin_inv( self )
        class(image), intent(inout) :: self
        self%rmat = -1.*(self%rmat-1.)
    end subroutine bin_inv

    module subroutine remove_edge( self )
        class(image), intent(inout) :: self
        if( self%ft ) THROW_HARD('only for real binary images (not FTed ones); remove_edge')
        if( any(self%rmat > 1.0001) .or. any(self%rmat < 0. ))&
            THROW_HARD('input to remove edge not binary; remove_edge')
        where( self%rmat < 0.999 ) self%rmat = 0.
    end subroutine remove_edge

    module subroutine zero2one( self )
        class(image), intent(inout) :: self
        integer :: i, j, k
        do i=1,self%ldim(1)
            do j=1,self%ldim(2)
                do k=1,self%ldim(3)
                    if( is_zero(self%rmat(i,j,k)) ) self%rmat(i,j,k) = 1.
                end do
            end do
        end do
    end subroutine zero2one

    module subroutine one_at_edge( self )
        class(image), intent(inout) :: self
        if( self%ft ) THROW_HARD('only for real binary images (not FTed ones); one_at_edge')
        if( any(self%rmat > 1.0001) .or. any(self%rmat < 0. ))&
            THROW_HARD('input to one_at_edge not binary')
        where( self%rmat < 0.999 .and. self%rmat > TINY ) self%rmat = 1.
    end subroutine one_at_edge

    module function bin2logical( self ) result( mask )
        class(image), intent(in)  :: self
        logical,      allocatable :: mask(:,:,:)
        allocate(mask(self%ldim(1),self%ldim(2),self%ldim(3)),source=.false.)
        where( self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) > TINY )
            mask = .true.
        end where
    end function bin2logical

    module subroutine logical2bin( self, mask )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: mask(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3))
        self%rmat = 0.
        where(mask) self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) = 1.
    end subroutine logical2bin

    module subroutine density_inoutside( self, msk, nin, nout, nmsk )
        class(image), intent(in)  :: self
        real,         intent(in)  :: msk
        integer,      intent(out) :: nin, nout, nmsk
        integer :: i
        real    :: centre(3), cendists(self%ldim(1),self%ldim(2)), msksq
        ! calculate distances from the center
        centre = real(self%ldim)/2.+1.
        cendists = 0.
        do i=1,self%ldim(1)
            cendists(i,:) = cendists(i,:) + (real(i)-centre(1))**2.
        enddo
        do i=1,self%ldim(2)
            cendists(:,i) = cendists(:,i) + (real(i)-centre(2))**2.
        enddo
        ! pixels forming the mask
        msksq = msk**2
        nmsk  = count(cendists <= msksq)
        ! pixels inside mask that are foreground
        nin = count(cendists <= msksq .and. self%rmat(:self%ldim(1),:self%ldim(2),1) > 0.5)
        ! pixels outside mask that are foreground
        nout = count(cendists > msksq .and. self%rmat(:self%ldim(1),:self%ldim(2),1) > 0.5)
    end subroutine density_inoutside

    module subroutine calc_bin_thres( self, frac_fg_target, thres )
        class(image), intent(inout) :: self
        real,         intent(in)    :: frac_fg_target
        real,         intent(out)   :: thres 
        integer, parameter   :: NQUANTA = 10
        real,    allocatable :: pixvals(:), means(:), pix_ts(:), frac_fgs(:)
        integer, allocatable :: labels(:)
        integer :: iq, n_fg, npix, loc(1), cnt
        allocate( means(NQUANTA), labels(NQUANTA), frac_fgs(NQUANTA), pix_ts(NQUANTA) )
        means   = 0.
        labels  = 0
        pixvals = pack(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)), mask=.true.)
        call sortmeans(pixvals, NQUANTA, means, labels)
        npix    = product(self%ldim)
        cnt     = 0
        do iq = 1, NQUANTA
            if( count(labels == iq) == 0 )cycle
            cnt         = cnt + 1
            pix_ts(cnt) = minval(pixvals, mask=labels == iq )
            n_fg        = count(self%rmat(:self%ldim(1),:self%ldim(2),:self%ldim(3)) >= pix_ts(cnt))
            if( n_fg == npix )then
                frac_fgs(cnt) = 1.
            else if( n_fg == 0 )then
                frac_fgs(cnt) = 0.
            else
                frac_fgs(cnt) = real(n_fg) / real(npix)
            endif
        end do
        loc      = minloc(abs(frac_fgs(:cnt) - frac_fg_target))
        thres    = pix_ts(loc(1))
    end subroutine calc_bin_thres

    module subroutine mask( self, mskrad, which, inner, width, backgr )
        class(image),     intent(inout) :: self
        real,             intent(in)    :: mskrad
        character(len=*), intent(in)    :: which
        real, optional,   intent(in)    :: inner, width, backgr
        real(dp) :: sumv
        real     :: e, wwidth, d_sq, rad_sq, ave
        real     :: cis(self%ldim(1)), cjs(self%ldim(2)), cks(self%ldim(3))
        integer  :: i, j, k, minlen, ir, jr, kr, npix
        logical  :: didft, doinner, soft, avg_backgr
        ! width
        wwidth = 10.
        if( present(width) ) wwidth = width
        ! inner
        doinner = .false.
        if( present(inner) ) doinner = .true.
        ! FT
        didft = .false.
        if( self%ft )then
            call self%ifft()
            didft = .true.
        endif
        ! minlen
        if( self%is_3d() )then
            minlen = minval(self%ldim)
        else
            minlen = minval(self%ldim(1:2))
        endif
        ! soft mask width limited to +/- COSMSKHALFWIDTH pixels
        minlen = min(nint(2.*(mskrad+COSMSKHALFWIDTH)), minlen)
        ! soft/hard
        soft       = .true.
        avg_backgr = .false.
        select case(trim(which))
            case('soft')
                soft  = .true.
                if(present(backgr))then
                    self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) =&
                        &self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) - backgr
                else
                    call self%zero_background
                endif
            case('softavg')
                ! background is set to its average value
                if( doinner )THROW_HARD('Inner masking not supported with softavg; simple_image :: mask')
                soft       = .true.
                avg_backgr = .true.
                rad_sq     = mskrad*mskrad
                sumv       = 0.d0
                npix       = 0
            case('hard')
                soft  = .false.
                if(present(backgr))THROW_HARD('no backround subtraction with hard masking; simple_image :: mask')
            case DEFAULT
                THROW_HARD('undefined which parameter; simple_image :: mask')
        end select
        ! init center as origin
        forall(i=1:self%ldim(1)) cis(i) = -real(self%ldim(1))/2. + real(i-1)
        forall(i=1:self%ldim(2)) cjs(i) = -real(self%ldim(2))/2. + real(i-1)
        if(self%is_3d())forall(i=1:self%ldim(3)) cks(i) = -real(self%ldim(3))/2. + real(i-1)
        ! MASKING
        if( soft )then
            ! Soft masking
            if( self%is_3d() )then
                if( avg_backgr )then
                    do k=1,self%ldim(3)
                        do j=1,self%ldim(2)
                            do i=1,self%ldim(1)
                                d_sq = cis(i)*cis(i) + cjs(j)*cjs(j) + cks(k)*cks(k)
                                if( d_sq > rad_sq )then
                                    npix = npix + 1
                                    sumv = sumv + real(self%rmat(i,j,k),dp)
                                endif
                            enddo
                        enddo
                    enddo
                    if( npix > 0 )then
                        ave = real(sumv/real(npix,dp))
                        do i=1,self%ldim(1)
                            do j=1,self%ldim(2)
                                do k=1,self%ldim(3)
                                    e = cosedge(cis(i),cjs(j),cks(k),minlen,mskrad)
                                    if( e < 0.0001 )then
                                        self%rmat(i,j,k) = ave
                                    else
                                        if( e < 0.9999 )then
                                            self%rmat(i,j,k) = e*self%rmat(i,j,k) + (1.-e)*ave
                                        endif
                                    endif
                                enddo
                            enddo
                        enddo
                    endif
                else
                    ! 3d
                    do i=1,self%ldim(1)/2
                        ir = self%ldim(1)+1-i
                        do j=1,self%ldim(2)/2
                            jr = self%ldim(2)+1-j
                            do k=1,self%ldim(3)/2
                                e = cosedge(cis(i),cjs(j),cks(k),minlen,mskrad)
                                if( doinner )e = e * cosedge_inner(cis(i),cjs(j),cks(k),wwidth,inner)
                                if(e > 0.9999) cycle
                                kr = self%ldim(3)+1-k
                                self%rmat(i,j,k)    = e * self%rmat(i,j,k)
                                self%rmat(i,j,kr)   = e * self%rmat(i,j,kr)
                                self%rmat(i,jr,k)   = e * self%rmat(i,jr,k)
                                self%rmat(i,jr,kr)  = e * self%rmat(i,jr,kr)
                                self%rmat(ir,j,k)   = e * self%rmat(ir,j,k)
                                self%rmat(ir,j,kr)  = e * self%rmat(ir,j,kr)
                                self%rmat(ir,jr,k)  = e * self%rmat(ir,jr,k)
                                self%rmat(ir,jr,kr) = e * self%rmat(ir,jr,kr)
                            enddo
                        enddo
                    enddo
                endif
            else
                ! 2d
                if( avg_backgr )then
                    do j=1,self%ldim(2)
                        do i=1,self%ldim(1)
                            d_sq = cis(i)*cis(i) + cjs(j)*cjs(j)
                            if( d_sq > rad_sq )then
                                npix = npix + 1
                                sumv = sumv + real(self%rmat(i,j,1),dp)
                            endif
                        enddo
                    enddo
                    if( npix > 0 )then
                        ave  = real(sumv/real(npix,dp))
                        do i=1,self%ldim(1)
                            do j=1,self%ldim(2)
                                e = cosedge(cis(i),cjs(j),minlen,mskrad)
                                if( e < 0.0001 )then
                                    self%rmat(i,j,1) = ave
                                else
                                    if( e < 0.9999 )then
                                        self%rmat(i,j,1) = e*self%rmat(i,j,1) + (1.-e)*ave
                                    endif
                                endif
                            enddo
                        enddo
                    endif
                else
                    do i=1,self%ldim(1)/2
                        ir = self%ldim(1)+1-i
                        do j=1,self%ldim(2)/2
                            e = cosedge(cis(i),cjs(j),minlen,mskrad)
                            if( doinner )e = e * cosedge_inner(cis(i),cjs(j),wwidth,inner)
                            if(e > 0.9999)cycle
                            jr = self%ldim(2)+1-j
                            self%rmat(i,j,1)   = e * self%rmat(i,j,1)
                            self%rmat(i,jr,1)  = e * self%rmat(i,jr,1)
                            self%rmat(ir,j,1)  = e * self%rmat(ir,j,1)
                            self%rmat(ir,jr,1) = e * self%rmat(ir,jr,1)
                        enddo
                    enddo
                endif
            endif
        else
            ! Hard masking
            if( self%is_3d() )then
                ! 3d
                do i=1,self%ldim(1)/2
                    ir = self%ldim(1)+1-i
                    do j=1,self%ldim(2)/2
                        jr = self%ldim(2)+1-j
                        do k=1,self%ldim(3)/2
                            e = hardedge(cis(i),cjs(j),cks(k),mskrad)
                            if( doinner )e = e * hardedge_inner(cis(i),cjs(j),cks(k),inner)
                            kr = self%ldim(3)+1-k
                            self%rmat(i,j,k)    = e * self%rmat(i,j,k)
                            self%rmat(i,j,kr)   = e * self%rmat(i,j,kr)
                            self%rmat(i,jr,k)   = e * self%rmat(i,jr,k)
                            self%rmat(i,jr,kr)  = e * self%rmat(i,jr,kr)
                            self%rmat(ir,j,k)   = e * self%rmat(ir,j,k)
                            self%rmat(ir,j,kr)  = e * self%rmat(ir,j,kr)
                            self%rmat(ir,jr,k)  = e * self%rmat(ir,jr,k)
                            self%rmat(ir,jr,kr) = e * self%rmat(ir,jr,kr)
                        enddo
                    enddo
                enddo
            else
                ! 2d
                do i=1,self%ldim(1)/2
                    ir = self%ldim(1)+1-i
                    do j=1,self%ldim(2)/2
                        jr = self%ldim(2)+1-j
                        e = hardedge(cis(i),cjs(j),mskrad)
                        if( doinner )e = e * hardedge_inner(ir,jr,inner)
                        self%rmat(i,j,1)   = e * self%rmat(i,j,1)
                        self%rmat(i,jr,1)  = e * self%rmat(i,jr,1)
                        self%rmat(ir,j,1)  = e * self%rmat(ir,j,1)
                        self%rmat(ir,jr,1) = e * self%rmat(ir,jr,1)
                    enddo
                enddo
            endif
        endif
        if( didft ) call self%fft()
    end subroutine mask

    !>  \brief  Taper edges of image so that there are no sharp discontinuities in real space
    !!          This is a re-implementation of the MRC program taperedgek.for (Richard Henderson, 1987)
    !!          I stole it from CTFFIND4 (thanks Alexis for the beautiful re-implementation)
    module subroutine taper_edges( self )
        class(image), intent(inout) :: self
        real, allocatable  :: avg_curr_edge_start(:,:)
        real, allocatable  :: avg_curr_edge_stop(:,:)
        real, allocatable  :: avg_curr_edge_avg(:,:)
        real, allocatable  :: smooth_avg_curr_edge_start(:,:)
        real, allocatable  :: smooth_avg_curr_edge_stop(:,:)
        integer            :: curr_dim, ndims
        integer            :: dim2, dim3
        integer            :: i,j,k
        integer            :: j_shift,k_shift
        integer            :: jj,kk
        integer            :: nvals_runnavg
        integer, parameter :: avg_strip_width(3)   = 100
        integer, parameter :: taper_strip_width(3) = 500
        integer, parameter :: smooth_half_width(3) = 1
        ndims = 2
        ! initialise vars
        dim2 = 2;  dim3 = 3; nvals_runnavg = 0
        if (self%is_3d()) ndims = 3
        do curr_dim=1,ndims
            ! take care of dimensions
            select case (curr_dim)
            case (1)
                dim2 = 2
                dim3 = 3
            case (2)
                dim2 = 1
                dim3 = 3
            case (3)
                dim2 = 1
                dim3 = 2
            end select
            ! take care of allocation & initialisation
            if(allocated(avg_curr_edge_start))        deallocate(avg_curr_edge_start)
            if(allocated(avg_curr_edge_stop))         deallocate(avg_curr_edge_stop)
            if(allocated(avg_curr_edge_avg))          deallocate(avg_curr_edge_avg)
            if(allocated(smooth_avg_curr_edge_start)) deallocate(smooth_avg_curr_edge_start)
            if(allocated(smooth_avg_curr_edge_stop))  deallocate(smooth_avg_curr_edge_stop)
            allocate( avg_curr_edge_start(self%ldim(dim2),self%ldim(dim3)),&
                avg_curr_edge_stop (self%ldim(dim2),self%ldim(dim3)),&
                avg_curr_edge_avg(self%ldim(dim2),self%ldim(dim3)),&
                smooth_avg_curr_edge_start(self%ldim(dim2),self%ldim(dim3)),&
                smooth_avg_curr_edge_stop (self%ldim(dim2),self%ldim(dim3)))
            avg_curr_edge_start        = 0.0e0
            avg_curr_edge_stop         = 0.0e0
            avg_curr_edge_avg          = 0.0e0
            smooth_avg_curr_edge_start = 0.0e0
            smooth_avg_curr_edge_stop  = 0.0e0
            ! Deal with X=0 and X=self%ldim(1) edges
            i=1
            do k=1,self%ldim(dim3)
                do j=1,self%ldim(dim2)
                    select case (curr_dim)
                    case (1)
                        avg_curr_edge_start(j,k) =&
                            sum(self%rmat(1:avg_strip_width(curr_dim),j,k))&
                            /avg_strip_width(curr_dim)
                        avg_curr_edge_stop(j,k)  =&
                            sum(self%rmat(self%ldim(curr_dim)-avg_strip_width(1)+1:self%ldim(curr_dim),j,k))&
                            /avg_strip_width(curr_dim)
                    case (2)
                        avg_curr_edge_start(j,k) =&
                            sum(self%rmat(j,1:avg_strip_width(curr_dim),k))&
                            /avg_strip_width(curr_dim)
                        avg_curr_edge_stop(j,k) =&
                            sum(self%rmat(j,self%ldim(curr_dim)-avg_strip_width(1)+1:self%ldim(curr_dim),k))&
                            /avg_strip_width(curr_dim)
                    case (3)
                        avg_curr_edge_start(j,k) =&
                            sum(self%rmat(j,k,1:avg_strip_width(curr_dim)))&
                            /avg_strip_width(curr_dim)
                        avg_curr_edge_stop(j,k) =&
                            sum(self%rmat(j,k,self%ldim(curr_dim)-avg_strip_width(1)+1:self%ldim(curr_dim)))&
                            /avg_strip_width(curr_dim)
                    end select
                enddo
            enddo
            avg_curr_edge_avg   = 0.5e0*(avg_curr_edge_stop + avg_curr_edge_start)
            avg_curr_edge_start = avg_curr_edge_start - avg_curr_edge_avg
            avg_curr_edge_stop  = avg_curr_edge_stop - avg_curr_edge_avg
            ! Apply smoothing parallel to edge in the form of a running average
            do k=1,self%ldim(dim3)
                do j=1,self%ldim(dim2)
                    nvals_runnavg = 0
                    ! Loop over neighbourhood of non-smooth arrays
                    do k_shift=-smooth_half_width(dim3),smooth_half_width(dim3)
                        kk = k+k_shift
                        if (kk .lt. 1 .or. kk .gt. self%ldim(dim3)) cycle
                        do j_shift=-smooth_half_width(dim2),smooth_half_width(dim2)
                            jj = j+j_shift
                            if (jj .lt. 1 .or. jj .gt. self%ldim(dim2)) cycle
                            nvals_runnavg = nvals_runnavg + 1
                            smooth_avg_curr_edge_start (j,k) =&
                                smooth_avg_curr_edge_start(j,k)+avg_curr_edge_start(jj,kk)
                            smooth_avg_curr_edge_stop(j,k)   =&
                                smooth_avg_curr_edge_stop(j,k)+avg_curr_edge_stop(jj,kk)
                        enddo
                    enddo
                    ! Now we can compute the average
                    smooth_avg_curr_edge_start(j,k) = smooth_avg_curr_edge_start(j,k)/nvals_runnavg
                    smooth_avg_curr_edge_stop(j,k)   = smooth_avg_curr_edge_stop(j,k)/nvals_runnavg
                enddo
            enddo
            ! Taper the image
            do i=1,self%ldim(curr_dim)
                if (i .le. taper_strip_width(curr_dim)) then
                    select case (curr_dim)
                    case (1)
                        self%rmat(i,:,:) = self%rmat(i,:,:)&
                            - smooth_avg_curr_edge_start (:,:)&
                            * (taper_strip_width(curr_dim)-i+1)&
                            / taper_strip_width(curr_dim)
                    case (2)
                        self%rmat(1:self%ldim(1),i,:) = self%rmat(1:self%ldim(1),i,:)&
                            - smooth_avg_curr_edge_start(:,:)&
                            * (taper_strip_width(curr_dim)-i+1)&
                            / taper_strip_width(curr_dim)
                    case (3)
                        self%rmat(1:self%ldim(1),:,i) = self%rmat(1:self%ldim(1),:,i)&
                            - smooth_avg_curr_edge_start (:,:)&
                            * (taper_strip_width(curr_dim)-i+1)&
                            / taper_strip_width(curr_dim)
                    end select
                else if (i .ge. self%ldim(curr_dim)-taper_strip_width(curr_dim)+1) then
                    select case (curr_dim)
                    case (1)
                        self%rmat(i,:,:) = self%rmat(i,:,:)&
                            - smooth_avg_curr_edge_stop(:,:)&
                            * (taper_strip_width(curr_dim)+i&
                            - self%ldim(curr_dim))&
                            / taper_strip_width(curr_dim)
                    case (2)
                        self%rmat(1:self%ldim(1),i,:) = self%rmat(1:self%ldim(1),i,:)&
                            - smooth_avg_curr_edge_stop(:,:)&
                            * (taper_strip_width(curr_dim)+i&
                            - self%ldim(curr_dim))&
                            / taper_strip_width(curr_dim)
                    case (3)
                        self%rmat(1:self%ldim(1),:,i) = self%rmat(1:self%ldim(1),:,i)&
                            - smooth_avg_curr_edge_stop(:,:)&
                            * (taper_strip_width(curr_dim)+i&
                            - self%ldim(curr_dim))&
                            / taper_strip_width(curr_dim)
                    end select
                endif
            enddo
        enddo
        if(allocated(avg_curr_edge_start))        deallocate(avg_curr_edge_start)
        if(allocated(avg_curr_edge_stop))         deallocate(avg_curr_edge_stop)
        if(allocated(avg_curr_edge_avg))          deallocate(avg_curr_edge_avg)
        if(allocated(smooth_avg_curr_edge_start)) deallocate(smooth_avg_curr_edge_start)
        if(allocated(smooth_avg_curr_edge_stop))  deallocate(smooth_avg_curr_edge_stop)
    end subroutine taper_edges

    module subroutine taper_edges_hann( self, borders )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: borders(2)
        real    :: w
        integer :: i,j,n
        if( self%is_ft() )THROW_HARD('Real space only!')
        if( self%is_3d() )THROW_HARD('2D images only!')
        n = 0
        do j = 1,borders(2)
            n = n+1
            w = 1. - cos(PIO2* real(j-1)/real(borders(2)) )
            w = min(1.,max(0.,w))
            self%rmat(1:self%ldim(1),j,1) = w * self%rmat(1:self%ldim(1),j,1)
        enddo
        n = 0
        do j = self%ldim(2)-borders(2)+1,self%ldim(2)
            n = n+1
            w = cos(PIO2* real(j-self%ldim(2)+borders(2))/real(borders(2)))
            w = min(1.,max(0.,w))
            self%rmat(1:self%ldim(1),j,1) = w * self%rmat(1:self%ldim(1),j,1)
        enddo
        do i = 1,borders(1)
            w = 1. - cos(PIO2* real(i-1)/real(borders(1)) )
            w = min(1.,max(0.,w))
            self%rmat(i,1:self%ldim(2),1) = w * self%rmat(i,1:self%ldim(2),1)
        enddo
        do i = self%ldim(1)-borders(1)+1,self%ldim(1)
            w = cos(PIO2* real(i-self%ldim(1)+borders(1))/real(borders(1)))
            w = min(1.,max(0.,w))
            self%rmat(i,1:self%ldim(2),1) = w * self%rmat(i,1:self%ldim(2),1)
        enddo
    end subroutine taper_edges_hann

end submodule simple_image_seg
