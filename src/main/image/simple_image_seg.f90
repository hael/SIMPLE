!@descr: image segmentation related stuff to support masking
submodule (simple_image) simple_image_seg
implicit none
#include "simple_local_flags.inc"

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
    
    module subroutine memoize_mask_coords( self, box )
        class(image),      intent(in) :: self
        integer, optional, intent(in) :: box
        integer :: i, box_here
        if( OMP_IN_PARALLEL() )then
            THROW_HARD('No memoization inside OpenMP regions')
        endif
        if( present(box) )then
            ! override such that the image does not need to always exist
            if( is_odd(box) ) THROW_HARD('Even image dimensions only!')
            box_here = box
        else
            if( self%existence )then
                if( self%ldim(1) /= self%ldim(2) )then
                    THROW_HARD('square images assumed')
                endif
                if( (self%ldim(3) == 1) .or. (self%ldim(3) == self%ldim(2)) )then
                    ! all good
                else
                    THROW_HARD('square images assumed')
                endif
                box_here = self%ldim(1)
            else
                THROW_HARD('Image has not been initialized!')
            endif
        endif
        ! Rebuild only if geometry in the first two logical dimensions has changed
        if (mem_msk_box == box_here ) then
            return
        endif
        ! Update cache key
        mem_msk_box = box_here
        ! (Re)allocate
        if( allocated(mem_msk_cs)  ) deallocate(mem_msk_cs)
        if( allocated(mem_msk_cs2) ) deallocate(mem_msk_cs2)
        ! memoize 3D coordinates
        allocate(mem_msk_cs(box_here), mem_msk_cs2(box_here))
        ! Fill (origin at center) + squares
        do i = 1, box_here
            mem_msk_cs(i)  = -0.5*real(box_here) + real(i-1)
            mem_msk_cs2(i) = mem_msk_cs(i) * mem_msk_cs(i)
        end do
    end subroutine memoize_mask_coords

    module subroutine mask2D_soft(self, mskrad, width, backgr)
        class(image),     intent(inout) :: self
        real,             intent(in)    :: mskrad
        real, optional,   intent(in)    :: width, backgr
        real     :: wwidth, rad_sq, ave,  r2, e, cjs2
        integer  :: minlen, npix,  i, j, ir, jr, n1, n2, n3, h1, h2
        if( self%ldim(3) > 1 )             THROW_HARD('not for 3D')
        if( self%ldim(1) /= mem_msk_box  ) THROW_HARD('incongruent mask memoization')
        ! width
        wwidth = 10.0
        if (present(width)) wwidth = width
        ! dims
        n1 = self%ldim(1)
        n2 = self%ldim(2)
        n3 = 1
        ! minlen
        minlen = minval(self%ldim(1:2))
        minlen = min(nint(2.0*(mskrad + COSMSKHALFWIDTH)), minlen)
        ! mode
        if (present(backgr)) then
            self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) = &
                self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) - backgr
        else
            call self%zero_background
        endif
        h1 = n1/2; h2 = n2/2
        do j = 1, h2
            jr = n2 + 1 - j
            cjs2 = mem_msk_cs2(j)
            !$omp simd
            do i = 1, h1
                ir = n1 + 1 - i
                r2 = mem_msk_cs2(i) + cjs2
                e  = cosedge_r2_2d(r2, minlen, mskrad)
                if (e > 0.9999) cycle
                self%rmat(i ,j ,1)  = e * self%rmat(i ,j ,1)
                self%rmat(i ,jr,1)  = e * self%rmat(i ,jr,1)
                self%rmat(ir,j ,1)  = e * self%rmat(ir,j ,1)
                self%rmat(ir,jr,1)  = e * self%rmat(ir,jr,1)
            end do
        end do
    end subroutine mask2D_soft

    module subroutine mask2D_softavg(self, mskrad, width, backgr)
        class(image),     intent(inout) :: self
        real,             intent(in)    :: mskrad
        real, optional,   intent(in)    :: width, backgr
        real(dp) :: sumv, sv
        real     :: wwidth, rad_sq, ave, r2, e, cjs2
        integer  :: minlen, npix, i, j, np, n1l, n2l
        integer  :: n1, n2, n3
        logical  :: soft, avg_backgr
        if( self%ldim(3) > 1 )             THROW_HARD('not for 3D')
        if( self%ldim(1) /= mem_msk_box  ) THROW_HARD('incongruent mask memoization')
        ! width
        wwidth = 10.0
        if (present(width)) wwidth = width
        ! dims
        n1 = self%ldim(1)
        n2 = self%ldim(2)
        n3 = 1
        ! minlen
        minlen = minval(self%ldim(1:2))
        minlen = min(nint(2.0*(mskrad + COSMSKHALFWIDTH)), minlen)
        ! mode
        rad_sq     = mskrad * mskrad
        sumv       = 0.0_dp
        npix       = 0            
        n1l = n1; n2l = n2
        sv = 0.0_dp; np = 0
        ! avg outside radius
        do j = 1, n2l
            cjs2 = mem_msk_cs2(j)
            do i = 1, n1l
                r2 = mem_msk_cs2(i) + cjs2 
                if (r2 > rad_sq) then
                    np = np + 1
                    sv = sv + real(self%rmat(i,j,1), dp)
                endif
            end do
        end do
        if (np <= 0) return
        ave = real(sv / real(np, dp))
        ! apply (j,i with i contiguous)
        do j = 1, n2l
            cjs2 = mem_msk_cs2(j)
            do i = 1, n1l
                r2 = mem_msk_cs2(i) + cjs2 
                e  = cosedge_r2_2d(r2, minlen, mskrad)
                if (e < 0.0001) then
                    self%rmat(i,j,1) = ave
                else if (e < 0.9999) then
                    self%rmat(i,j,1) = e*self%rmat(i,j,1) + (1.0-e)*ave
                endif
            end do
        end do
    end subroutine mask2D_softavg

    module subroutine mask2D_hard(self, mskrad)
        class(image),     intent(inout) :: self
        real,             intent(in)    :: mskrad
        integer :: i, j, ir, jr, h1, h2, n1, n2, n3
        real :: r2, e, cjs2
        if( self%ldim(3) > 1 )             THROW_HARD('not for 3D')
        if( self%ldim(1) /= mem_msk_box  ) THROW_HARD('incongruent mask memoization')
        ! dims
        n1 = self%ldim(1)
        n2 = self%ldim(2)
        n3 = 1
        h1 = n1/2; h2 = n2/2
        do j = 1, h2
            jr   = n2 + 1 - j
            cjs2 = mem_msk_cs2(j)
            !$omp simd
            do i = 1, h1
                ir = n1 + 1 - i
                r2 = mem_msk_cs2(i) + cjs2
                e  = hardedge_r2_2d(r2, mskrad)
                if (e > 0.9999) cycle
                self%rmat(i ,j ,1)  = e * self%rmat(i ,j ,1)
                self%rmat(i ,jr,1)  = e * self%rmat(i ,jr,1)
                self%rmat(ir,j ,1)  = e * self%rmat(ir,j ,1)
                self%rmat(ir,jr,1)  = e * self%rmat(ir,jr,1)
            end do
        end do
    end subroutine mask2D_hard

    module subroutine mask3D_soft(self, mskrad, width, backgr)
        class(image),     intent(inout) :: self
        real,             intent(in)    :: mskrad
        real, optional,   intent(in)    :: width, backgr
        real(dp) :: sumv
        real     :: wwidth, rad_sq, ave, r2, e, cjs2, cks2
        integer  :: minlen, npix, n1, n2, n3, i, j, k, ir, jr, kr, h1, h2, h3
        if( self%ldim(3) == 1 ) THROW_HARD('not for 2D')
        if( self%ldim(1) /= mem_msk_box  )then
            if( OMP_IN_PARALLEL() )then
                THROW_HARD('incongruent mask memoization')
            else
                call memoize_mask_coords(self)
            endif
        endif
        ! width
        wwidth = 10.0
        if (present(width)) wwidth = width
        ! dims
        n1 = self%ldim(1)
        n2 = self%ldim(2)
        n3 = self%ldim(3)
        ! minlen
        minlen = minval(self%ldim)
        minlen = min(nint(2.0*(mskrad + COSMSKHALFWIDTH)), minlen)
        ! mode
        if (present(backgr)) then
            self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) = &
                self%rmat(1:self%ldim(1),1:self%ldim(2),1:self%ldim(3)) - backgr
        else
            call self%zero_background
        endif            
        h1 = n1/2; h2 = n2/2; h3 = n3/2
        do j = 1, h2
            jr   = n2 + 1 - j
            cjs2 = mem_msk_cs2(j)
            do k = 1, h3
                kr   = n3 + 1 - k
                cks2 = mem_msk_cs2(k)
                !$omp simd
                do i = 1, h1
                    ir = n1 + 1 - i
                    r2 = mem_msk_cs2(i) + cjs2 + cks2
                    e  = cosedge_r2_3d(r2, minlen, mskrad)
                    if (e > 0.9999) cycle
                    self%rmat(i ,j ,k )  = e * self%rmat(i ,j ,k )
                    self%rmat(i ,j ,kr)  = e * self%rmat(i ,j ,kr)
                    self%rmat(i ,jr,k )  = e * self%rmat(i ,jr,k )
                    self%rmat(i ,jr,kr)  = e * self%rmat(i ,jr,kr)
                    self%rmat(ir,j ,k )  = e * self%rmat(ir,j ,k )
                    self%rmat(ir,j ,kr)  = e * self%rmat(ir,j ,kr)
                    self%rmat(ir,jr,k )  = e * self%rmat(ir,jr,k )
                    self%rmat(ir,jr,kr)  = e * self%rmat(ir,jr,kr)
                end do
            end do
        end do
    end subroutine mask3D_soft

    module subroutine mask3D_softavg(self, mskrad, width, backgr)
        class(image),     intent(inout) :: self
        real,             intent(in)    :: mskrad
        real, optional,   intent(in)    :: width, backgr
        real(dp) :: sumv, sv
        real     :: wwidth, rad_sq, ave, r2, e, cjs2, cks2
        integer  :: minlen, npix,  i, j, k, n1, n2, n3, np, n1l, n2l, n3l
        if( self%ldim(3) == 1 ) THROW_HARD('not for 2D')
        if( self%ldim(1) /= mem_msk_box  )then
            if( OMP_IN_PARALLEL() )then
                THROW_HARD('incongruent mask memoization')
            else
                call memoize_mask_coords(self)
            endif
        endif
        ! width
        wwidth = 10.0
        if (present(width)) wwidth = width
        ! dims
        n1 = self%ldim(1)
        n2 = self%ldim(2)
        n3 = self%ldim(3)
        ! minlen
        minlen = minval(self%ldim)
        minlen = min(nint(2.0*(mskrad + COSMSKHALFWIDTH)), minlen)
        ! mode
        rad_sq     = mskrad * mskrad
        sumv       = 0.0_dp
        npix       = 0                
        n1l = n1; n2l = n2; n3l = n3
        sv = 0.0_dp; np = 0
        ! avg outside radius (r^2 only)
        do k = 1, n3l
            cks2 = mem_msk_cs2(k)
            do j = 1, n2l
                cjs2 = mem_msk_cs2(j)
                do i = 1, n1l
                    r2 = mem_msk_cs2(i) + cjs2 + cks2 
                    if (r2 > rad_sq) then
                        np = np + 1
                        sv = sv + real(self%rmat(i,j,k), dp)
                    endif
                end do
            end do
        end do
        if (np <= 0) return
        ave = real(sv / real(np, dp))
        ! apply (cache-friendly: k,j,i; i contiguous)
        do k = 1, n3l
            cks2 = mem_msk_cs2(k)
            do j = 1, n2l
                cjs2 = mem_msk_cs2(j)
                do i = 1, n1l
                    r2 = mem_msk_cs2(i) + cjs2 + cks2 
                    e  = cosedge_r2_3d(r2, minlen, mskrad)
                    if (e < 0.0001) then
                        self%rmat(i,j,k) = ave
                    else if (e < 0.9999) then
                        self%rmat(i,j,k) = e*self%rmat(i,j,k) + (1.0-e)*ave
                    endif
                end do
            end do
        end do
    end subroutine mask3D_softavg

    module subroutine mask3D_hard(self, mskrad)
        class(image), intent(inout) :: self
        real,         intent(in)    :: mskrad
        real     :: r2, e, cjs2, cks2
        integer  :: i, j, k, ir, jr, kr, h1, h2, h3, n1, n2, n3
        if( self%ldim(3) == 1 ) THROW_HARD('not for 2D')
        if( self%ldim(1) /= mem_msk_box  )then
            if( OMP_IN_PARALLEL() )then
                THROW_HARD('incongruent mask memoization')
            else
                call memoize_mask_coords(self)
            endif
        endif
        ! dims
        n1 = self%ldim(1)
        n2 = self%ldim(2)
        n3 = self%ldim(3)           
        h1 = n1/2; h2 = n2/2; h3 = n3/2
        do j = 1, h2
            jr   = n2 + 1 - j
            cjs2 = mem_msk_cs2(j)
            do k = 1, h3
                kr   = n3 + 1 - k
                cks2 = mem_msk_cs2(k)
                !$omp simd
                do i = 1, h1
                    ir = n1 + 1 - i
                    r2 = mem_msk_cs2(i) + cjs2 + cks2 
                    e  = hardedge_r2_3d(r2, mskrad)
                    if (e > 0.9999) cycle
                    self%rmat(i ,j ,k )  = e * self%rmat(i ,j ,k )
                    self%rmat(i ,j ,kr)  = e * self%rmat(i ,j ,kr)
                    self%rmat(i ,jr,k )  = e * self%rmat(i ,jr,k )
                    self%rmat(i ,jr,kr)  = e * self%rmat(i ,jr,kr)
                    self%rmat(ir,j ,k )  = e * self%rmat(ir,j ,k )
                    self%rmat(ir,j ,kr)  = e * self%rmat(ir,j ,kr)
                    self%rmat(ir,jr,k )  = e * self%rmat(ir,jr,k )
                    self%rmat(ir,jr,kr)  = e * self%rmat(ir,jr,kr)
                end do
            end do
        end do
    end subroutine mask3D_hard

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
