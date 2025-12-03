submodule (simple_image) simple_image_geom
!$ use omp_lib
!$ use omp_lib_kinds
include  'simple_lib.f08'
#include "simple_local_flags.inc"  
implicit none
contains

    ! windowing

    !>  \brief window extracts a particle image from a box as defined by EMAN 1.9
    module subroutine window( self_in, coord, box, self_out, noutside )
        class(image),      intent(in)    :: self_in
        integer,           intent(in)    :: coord(2), box
        class(image),      intent(inout) :: self_out
        integer, optional, intent(inout) :: noutside
        integer :: i,j, fromc(2), toc(2), xoshoot, yoshoot, xushoot, yushoot, xboxrange(2), yboxrange(2)
        if( self_in%ldim(3) > 1 ) THROW_HARD('only 4 2D images; window')
        if( self_in%is_ft() )     THROW_HARD('only 4 real images; window')
        if( self_out%exists() )then
            if( self_out%is_ft() ) THROW_HARD('only 4 real images; window')
            if( self_out%ldim(1) == box .and. self_out%ldim(2) == box .and. self_out%ldim(3) == 1 )then
                ! go ahead
            else
                call self_out%new([box,box,1], self_in%smpd)
            endif
        else
            call self_out%new([box,box,1], self_in%smpd)
        endif
        fromc = coord+1       ! compensate for the c-range that starts at 0
        toc   = fromc+(box-1) ! the lower left corner is 1,1
        if( any(fromc < 1) .or. toc(1) > self_in%ldim(1) .or. toc(2) > self_in%ldim(2) )then
            if( present(noutside) )then
                noutside = noutside+1
            else
                THROW_WARN('box extends outside micrograph; window')
            endif
        endif
        xoshoot = 0
        yoshoot = 0
        xushoot = 0
        yushoot = 0
        if( toc(1)   > self_in%ldim(1) ) xoshoot =  toc(1)   - self_in%ldim(1)
        if( toc(2)   > self_in%ldim(2) ) yoshoot =  toc(2)   - self_in%ldim(2)
        if( fromc(1) < 1               ) xushoot = -fromc(1) + 1
        if( fromc(2) < 1               ) yushoot = -fromc(2) + 1
        toc(1)        = toc(1)   - xoshoot
        toc(2)        = toc(2)   - yoshoot
        fromc(1)      = fromc(1) + xushoot
        fromc(2)      = fromc(2) + yushoot
        xboxrange(1)  = xushoot  + 1
        xboxrange(2)  = box      - xoshoot
        yboxrange(1)  = yushoot  + 1
        yboxrange(2)  = box      - yoshoot
        self_out%rmat = 0.
        self_out%rmat(xboxrange(1):xboxrange(2),yboxrange(1):yboxrange(2),1) = self_in%rmat(fromc(1):toc(1),fromc(2):toc(2),1)
        ! Rather than stretching the last inside line that creates stripes
        ! it is scrambled to mitigate Gibbs phenomenon and preserve image stats
        if( xboxrange(1) > 1 )then
            do i = 1,xboxrange(1)-1
                do j = 1,self_out%ldim(2)
                    self_out%rmat(i,j,1) = self_out%rmat(xboxrange(1),irnd_uni(self_out%ldim(2)),1)
                enddo
            enddo
        endif
        if( xboxrange(2) < self_out%ldim(1) )then
            do i = xboxrange(2)+1,self_out%ldim(1)
                do j = 1,self_out%ldim(2)
                    self_out%rmat(i,j,1) = self_out%rmat(xboxrange(2),irnd_uni(self_out%ldim(2)),1)
                enddo
            enddo
        endif
        if( yboxrange(1) > 1 )then
            do j = 1,yboxrange(1)-1
                do i = 1,self_out%ldim(1)
                    self_out%rmat(i,j,1) = self_out%rmat(irnd_uni(self_out%ldim(1)),yboxrange(1),1)
                enddo
            enddo
        endif
        if( yboxrange(2) < self_out%ldim(2) )then
            do j = yboxrange(2)+1,self_out%ldim(2)
                do i = 1,self_out%ldim(1)
                    self_out%rmat(i,j,1) = self_out%rmat(irnd_uni(self_out%ldim(1)),yboxrange(2),1)
                enddo
            enddo
        endif
    end subroutine window

    !>  window_slim  extracts a particle image from a box as defined by EMAN 1.9
    module subroutine window_slim( self_in, coord, box, self_out, outside )
        class(image), intent(in)    :: self_in
        integer,      intent(in)    :: coord(:), box
        class(image), intent(inout) :: self_out
        logical,      intent(out)   :: outside
        logical                     :: isvol
        integer, allocatable        :: fromc(:), toc(:)
        isvol         = self_in%is_3d()
        if( isvol )then
            allocate(fromc(3), toc(3))
        else
            allocate(fromc(2), toc(2))
        endif
        fromc         = coord + 1         ! compensate for the c-range that starts at 0
        toc           = fromc + (box - 1) ! the lower left corner is 1,1
        self_out%rmat = 0.
        self_out%ft   = .false.
        outside       = .false.
        if( isvol )then
            if( size(coord) /= 3 ) THROW_HARD("Error! expecting 3D coordinates; window_slim")
            if( fromc(1) < 1 .or. fromc(2) < 1 .or. fromc(3) < 1 .or. toc(1) > self_in%ldim(1) .or. toc(2) > self_in%ldim(2) .or. toc(3) > self_in%ldim(3) )then
                outside = .true.
            else
                self_out%rmat(1:box,1:box,1:box) = self_in%rmat(fromc(1):toc(1),fromc(2):toc(2),fromc(3):toc(3))
            endif
        else
            if( size(coord) /= 2 ) THROW_HARD("Error! expecting 2D coordinates; window_slim")
            if( fromc(1) < 1 .or. fromc(2) < 1 .or. toc(1) > self_in%ldim(1) .or. toc(2) > self_in%ldim(2) )then
                outside = .true.
            else
                self_out%rmat(1:box,1:box,1) = self_in%rmat(fromc(1):toc(1),fromc(2):toc(2),1)
            endif
        endif
    end subroutine window_slim

    module subroutine window_center( self_in, center, rad, self_out, outside )
        class(image), intent(in)    :: self_in
        integer,      intent(in)    :: center(:), rad
        class(image), intent(inout) :: self_out
        logical,      intent(out)   :: outside
        integer :: coord(3), box
        box = rad*2
        if( self_in%is_3d() )then
            coord = center - rad
            if( any(coord < 1) )      THROW_HARD('Error! window is out of the image')
            call window_slim(self_in, coord, box, self_out, outside)
        else
            coord(1:2) = center - rad
            if( any(coord(1:2) < 1) ) THROW_HARD('Error! window is out of the image')
            call window_slim(self_in, coord(1:2), box, self_out, outside)
        endif
    end subroutine window_center

    ! for re-generation of micrograph after convolutional PPCA
    module subroutine add_window( self, imgwin, coord, offset )
        class(image), intent(inout)   :: self
        class(image), intent(in)      :: imgwin
        integer,      intent(in)      :: coord(2)
        integer, optional, intent(in) :: offset  !if offset is present it doesn't sum windows, but take just the inner part
        integer ::  fromc(2), toc(2), ld(3), box
        ld    = imgwin%get_ldim()
        box   = ld(1)
        fromc = coord + 1         ! compensate for the c-range that starts at 0
        toc   = fromc + (box - 1) ! the lower left corner is 1,1
        if( fromc(1) < 1 .or. fromc(2) < 1 .or. toc(1) > self%ldim(1) .or. toc(2) > self%ldim(2) )then
          return
        endif
        if(present(offset)) then  !no sovrapposition
            if(offset == int(box/2)) THROW_HARD("invalid offset choice; add_window") ! box is supposet to be even
            if (offset > int(box/2)) then
                self%rmat(fromc(1)+box-offset-1:offset,fromc(2)+box-offset-1:offset,1) = &
                    &imgwin%rmat(box-offset:offset,box-offset:offset,1)
            else
                self%rmat(fromc(1)+offset-1:box-offset,fromc(2)+offset-1:box-offset,1) = &
                    &imgwin%rmat(offset:box-offset,offset:box-offset,1)
            endif
        else
            self%rmat(fromc(1):toc(1),fromc(2):toc(2),1) = self%rmat(fromc(1):toc(1),fromc(2):toc(2),1) &
                &+ imgwin%rmat(1:box,1:box,1) !add everything
        endif
    end subroutine add_window

    !>  \brief win2arr extracts a small window into an array (circular indexing)
    module function win2arr( self, i, j, k, winsz ) result( pixels )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: i, j, k, winsz
        real, allocatable :: pixels(:)
        integer :: s, ss, t, tt, u, uu, cnt, npix
        if( self%is_ft() ) THROW_HARD('only 4 real images; win2arr')
        if( self%is_3d() )then
            npix = (2*winsz+1)**3
        else
            npix = (2*winsz+1)**2
        endif
        allocate(pixels(npix))
        cnt = 1
        do s=i-winsz,i+winsz
            ss = cyci_1d_static(self%ldim(1), s)
            do t=j-winsz,j+winsz
                tt = cyci_1d_static(self%ldim(2), t)
                if( self%ldim(3) > 1 )then
                    do u=k-winsz,k+winsz
                        uu          = cyci_1d_static(self%ldim(3), u)
                        pixels(cnt) = self%rmat(ss,tt,uu)
                        cnt         = cnt+1
                    end do
                else
                    pixels(cnt) = self%rmat(ss,tt,1)
                    cnt         = cnt+1
                endif
            end do
        end do
    end function win2arr

    !>  \brief win2arr_rad extracts a small window into an array (circular indexing) with radial thres
    module subroutine win2arr_rad( self, i, j, k, winsz, npix_in, maxrad, npix_out, pixels )
        class(image), intent(in)    :: self
        integer,      intent(in)    :: i, j, k, winsz, npix_in
        real,         intent(in)    :: maxrad ! in pixels
        integer,      intent(inout) :: npix_out
        real,         intent(inout) :: pixels(npix_in)
        integer :: s, ss, t, tt, u, uu
        real    :: vec_cen(3), vec_displ(3), maxradsq
        if( self%is_ft() ) THROW_HARD('only 4 real images; win2arr_rad')
        if( self%is_3d() )then
            vec_cen = real([i,j,k])
        else
            vec_cen = real([i,j,1])
        endif
        maxradsq = maxrad**2.
        npix_out = 1
        do s = i - winsz,i + winsz
            ss = cyci_1d_static(self%ldim(1), s)
            do t = j - winsz, j + winsz
                tt = cyci_1d_static(self%ldim(2), t)
                if( self%ldim(3) > 1 )then
                    do u = k - winsz, k + winsz
                        vec_displ = real([s,t,u]) - vec_cen
                        if( sum(vec_displ**2.) <= maxradsq )then
                            uu               = cyci_1d_static(self%ldim(3), u)
                            pixels(npix_out) = self%rmat(ss,tt,uu)
                            npix_out         = npix_out + 1
                        endif
                    end do
                else
                    vec_displ = real([s,t,1]) - vec_cen
                    if( sum(vec_displ(1:2))**2. <= maxradsq )then
                        pixels(npix_out) = self%rmat(ss,tt,1)
                        npix_out         = npix_out + 1
                    endif
                endif
            end do
        end do
    end subroutine win2arr_rad

    ! shapes

    !>  \brief corner extracts a corner of a volume with size box
    module subroutine corner( self_in, box, self_out )
        class(image), intent(in)    :: self_in
        integer,      intent(in)    :: box
        type(image),  intent(inout) :: self_out
        if( self_in%ldim(3) <= 1 ) THROW_HARD('only 4 3D images; corner')
        if( self_in%is_ft() )      THROW_HARD('only 4 real images; corner')
        call self_out%new([box,box,box], self_in%smpd)
        self_out%rmat(:box,:box,:box) = self_in%rmat(:box,:box,:box)
    end subroutine corner

    !>  \brief  just a corner filling fun for testing purposes
    !! \param sqrad half width of square
    module subroutine corners( self, sqrad )
        class(image), intent(inout) :: self
        integer, intent(in)         :: sqrad
        integer :: i, j
        self%rmat = 0.
        do i=self%ldim(1)-sqrad+1,self%ldim(1)
            do j=self%ldim(2)-sqrad+1,self%ldim(2)
                self%rmat(i,j,1) = 1.
            end do
        end do
        do i=1,sqrad
            do j=1,sqrad
                self%rmat(i,j,1) = 1.
            end do
        end do
        do i=self%ldim(1)-sqrad+1,self%ldim(1)
            do j=1,sqrad
                self%rmat(i,j,1) = 1.
            end do
        end do
        do i=1,sqrad
            do j=self%ldim(2)-sqrad+1,self%ldim(2)
                self%rmat(i,j,1) = 1.
            end do
        end do
        self%ft = .false.
    end subroutine corners

    module subroutine gauimg_1( self, wsz , alpha )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: wsz
        real, intent(in), optional  :: alpha
        real    :: x, y, z, xw, yw, zw, a
        integer :: i, j, k
        a = 0.5
        if(present(alpha)) a= alpha
        x = -real(self%ldim(1))/2.
        do i=1,self%ldim(1)
            xw = gauwfun(x, a*real(wsz))
            y = -real(self%ldim(2))/2.
            do j=1,self%ldim(2)
                yw = gauwfun(y, a*real(wsz))
                z = -real(self%ldim(3))/2.
                do k=1,self%ldim(3)
                    if( self%ldim(3) > 1 )then
                        zw = gauwfun(z, a*real(wsz))
                    else
                        zw = 1.
                    endif
                    self%rmat(i,j,k) = xw*yw*zw
                    z = z+1.
                end do
                y = y+1.
            end do
            x = x+1.
        end do
        self%ft = .false.
    end subroutine gauimg_1

    module subroutine gauimg_2( self, wsz, offx,offy)
        class(image), intent(inout) :: self
        integer,      intent(in)    :: wsz, offx, offy
        real    :: x, y, z, xw, yw, zw
        integer :: i, j, k
        x = -real(self%ldim(1))/2. -real(offx)
        do i=1,self%ldim(1)
            xw = gauwfun(x, 0.5*real(wsz))
            y = -real(self%ldim(2))/2. - real(offy)
            do j=1,self%ldim(2)
                yw = gauwfun(y, 0.5*real(wsz))
                z = -real(self%ldim(3))/2.
                do k=1,self%ldim(3)
                    if( self%ldim(3) > 1 )then
                        zw = gauwfun(z, 0.5*real(wsz))
                    else
                        zw = 1.
                    endif
                    self%rmat(i,j,k) =  self%rmat(i,j,k) + xw*yw*zw
                    z = z+1.
                end do
                y = y+1.
            end do
            x = x+1.
        end do
        self%ft = .false.
    end subroutine gauimg_2

    module subroutine gauimg2D( self, xsigma, ysigma, cutoff )
        class(image),   intent(inout) :: self
        real,           intent(in)    :: xsigma, ysigma
        real, optional, intent(in)    :: cutoff
        real    :: xx, x, y, cutoffsq_here, center(2)
        integer :: i, j
        cutoffsq_here = huge(cutoffsq_here)
        if(present(cutoff)) cutoffsq_here = cutoff*cutoff
        ! Center of the img is assumed self%ldim/2 + 1
        center = real(self%ldim(1:2))/2.+1.
        do i=1,self%ldim(1)
            x = real(i)-center(1)
            xx = x*x
            do j=1,self%ldim(2)
                y = real(j)-center(2)
                if(xx+y*y > cutoffsq_here) cycle
                self%rmat(i,j,1) = gaussian2D( [0.,0.], x, y, xsigma, ysigma )
            enddo
        enddo
        self%ft = .false.
    end subroutine gauimg2D

    module subroutine gauimg3D( self, xsigma, ysigma, zsigma, cutoff )
        class(image),   intent(inout) :: self
        real,           intent(in)    :: xsigma, ysigma, zsigma
        real, optional, intent(in)    :: cutoff
        real    :: xx, x, yy, y, z, cutoffsq_here, center(3)
        integer :: i, j, k
        if(self%ldim(3) == 1) then
            call self%gauimg2D(xsigma, ysigma, cutoff)
            return
        endif
        cutoffsq_here = huge(cutoffsq_here)
        if(present(cutoff)) cutoffsq_here = cutoff*cutoff
        ! Center of the img is assumed self%ldim/2 + 1
        center = real(self%ldim(1:3))/2.+1.
        do i=1,self%ldim(1)
            x = real(i)-center(1)
            xx = x*x
            do j=1,self%ldim(2)
                y = real(j)-center(2)
                yy = y*y
                do k=1,self%ldim(3)
                    z = real(k)-center(3)
                    if(xx+yy+z*z > cutoffsq_here) cycle
                    self%rmat(i,j,k) = gaussian3D( [0.,0.,0.], x, y, z, xsigma, ysigma, zsigma )
                enddo
            enddo
        enddo
        self%ft = .false.
    end subroutine gauimg3D

    module subroutine reshape2cube( self, self_out )
        class(image), intent(inout) :: self
        class(image), intent(out)   :: self_out
        logical :: isvol
        integer :: ldim(3), ldim_max
        real    :: smpd
        smpd  = self%get_smpd()
        isvol = self%is_3d()
        if(.not.isvol) THROW_HARD('this is only for volumes; reshape2cube')
        ldim  = self%ldim
        if( ldim(1) == ldim(2) .and. ldim(2) == ldim(3) ) return
        ldim_max      = max(ldim(1),ldim(2),ldim(3))
        call self_out%new([ldim_max,ldim_max,ldim_max], self%smpd, wthreads=self%wthreads)
        self_out%rmat = 0.
        self_out%rmat = self%rmat
        self_out%ft   = .false.
        call self_out%set_smpd(smpd)
    end subroutine reshape2cube

    !>  \brief square just a binary square for testing purposes
    !! \param sqrad half width of square
    module subroutine square( self, sqrad )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: sqrad
        integer :: i, j, k
        self%rmat = 0.
        if( all(self%ldim(1:2) .gt. sqrad) .and. self%ldim(3) == 1 ) then
            do i=self%ldim(1)/2-sqrad+1,self%ldim(1)/2+sqrad
                do j=self%ldim(2)/2-sqrad+1,self%ldim(2)/2+sqrad
                    self%rmat(i,j,1) = 1.
                end do
            end do
        else if( all(self%ldim .gt. sqrad) .and. self%ldim(3) > 1 )then
            do i=self%ldim(1)/2-sqrad+1,self%ldim(1)/2+sqrad
                do j=self%ldim(2)/2-sqrad+1,self%ldim(2)/2+sqrad
                    do k=self%ldim(3)/2-sqrad+1,self%ldim(3)/2+sqrad
                        self%rmat(i,j,k) = 1.
                    end do
                end do
            end do
        else
            THROW_HARD('image is to small to fit the square; square')
        endif
        self%ft = .false.
    end subroutine square

    ! pad/clip/crop

    module subroutine pad( self_in, self_out, backgr, antialiasing )
        class(image),      intent(inout) :: self_in, self_out
        real,    optional, intent(in)    :: backgr
        logical, optional, intent(in)    :: antialiasing
        real, allocatable :: antialw(:)
        real              :: w, ratio
        integer           :: starts(3), stops(3), lims(3,2)
        integer           :: h, k, l, phys_in(3), phys_out(3)
        logical           :: l_antialiasing
        if( self_in.eqdims.self_out )then
            call self_out%copy(self_in)
            return
        endif
        l_antialiasing = .true.
        if( present(antialiasing) ) l_antialiasing = antialiasing
        if( self_out%ldim(1) >= self_in%ldim(1) .and. self_out%ldim(2) >= self_in%ldim(2)&
        .and. self_out%ldim(3) >= self_in%ldim(3) )then
            if( self_in%ft )then
                self_out = cmplx(0.,0.)
                lims = self_in%fit%loop_lims(2)
                if( l_antialiasing )then
                    antialw = self_in%hannw()
                    !$omp parallel do collapse(3) schedule(static) default(shared)&
                    !$omp private(h,k,l,w,phys_out,phys_in) proc_bind(close)
                    do h=lims(1,1),lims(1,2)
                        do k=lims(2,1),lims(2,2)
                            do l=lims(3,1),lims(3,2)
                                w = antialw(max(1,abs(h)))*antialw(max(1,abs(k)))*antialw(max(1,abs(l)))
                                phys_out = self_out%fit%comp_addr_phys(h,k,l)
                                phys_in  = self_in%fit%comp_addr_phys(h,k,l)
                                self_out%cmat(phys_out(1),phys_out(2),phys_out(3))=&
                                self_in%cmat(phys_in(1),phys_in(2),phys_in(3))*w
                            end do
                        end do
                    end do
                    !$omp end parallel do
                    deallocate(antialw)
                else
                    !$omp parallel do collapse(3) schedule(static) default(shared)&
                    !$omp private(h,k,l,phys_out,phys_in) proc_bind(close)
                    do h=lims(1,1),lims(1,2)
                        do k=lims(2,1),lims(2,2)
                            do l=lims(3,1),lims(3,2)
                                phys_out = self_out%fit%comp_addr_phys(h,k,l)
                                phys_in  = self_in%fit%comp_addr_phys(h,k,l)
                                self_out%cmat(phys_out(1),phys_out(2),phys_out(3))= self_in%cmat(phys_in(1),phys_in(2),phys_in(3))
                            end do
                        end do
                    end do
                    !$omp end parallel do
                endif
                ratio = real(self_in%ldim(1))/real(self_out%ldim(1))
                self_out%smpd = self_in%smpd*ratio ! padding Fourier transform, so sampling is finer
                self_out%ft = .true.
            else
                starts = (self_out%ldim-self_in%ldim)/2+1
                stops  = self_out%ldim-starts+1
                if( self_in%ldim(3) == 1 )then
                    starts(3) = 1
                    stops(3)  = 1
                endif
                if( present(backgr) )then
                    self_out%rmat = backgr
                else
                    self_out%rmat = 0.
                endif
                self_out%rmat(starts(1):stops(1),starts(2):stops(2),starts(3):stops(3)) =&
                self_in%rmat(:self_in%ldim(1),:self_in%ldim(2),:self_in%ldim(3))
                self_out%ft = .false.
            endif
        endif
    end subroutine pad

    module subroutine pad_inplace( self, ldim, backgr, antialiasing )
        class(image),   intent(inout) :: self
        integer,           intent(in) :: ldim(3)
        real,    optional, intent(in) :: backgr
        logical, optional, intent(in) :: antialiasing
        type(image) :: tmp
        call tmp%new(ldim, self%smpd, wthreads=self%wthreads)
        call self%pad(tmp, backgr=backgr, antialiasing=antialiasing)
        call self%copy(tmp)
        call tmp%kill()
    end subroutine pad_inplace

    module subroutine pad_mirr_1( self_in, self_out )
        class(image), intent(inout) :: self_in, self_out
        integer :: starts(3), stops(3)
        integer :: i,j, i_in, j_in
        if( self_in.eqdims.self_out )then
            call self_out%copy(self_in)
            return
        endif
        if(self_in%is_3d())THROW_HARD('2D images only; pad_mirr')
        if(self_in%ft)THROW_HARD('real space 2D images only; pad_mirr')
        if( self_out%ldim(1) >= self_in%ldim(1) .and. self_out%ldim(2) >= self_in%ldim(2))then
            self_out%rmat = 0.
            starts  = (self_out%ldim-self_in%ldim)/2+1
            stops   = self_out%ldim-starts+1
            ! actual image
            self_out%rmat(starts(1):stops(1),starts(2):stops(2),1) =&
                &self_in%rmat(:self_in%ldim(1),:self_in%ldim(2),1)
            ! left border
            i_in = 0
            do i = starts(1)-1,1,-1
                i_in = i_in + 1
                if(i_in > self_in%ldim(1))exit
                self_out%rmat(i,starts(2):stops(2),1) = self_in%rmat(i_in,:self_in%ldim(2),1)
            enddo
            ! right border
            i_in = self_in%ldim(1)+1
            do i=stops(1)+1,self_out%ldim(1)
                i_in = i_in - 1
                if(i_in < 1)exit
                self_out%rmat(i,starts(2):stops(2),1) = self_in%rmat(i_in,:self_in%ldim(2),1)
            enddo
            ! upper border & corners
            j_in = starts(2)
            do j = starts(2)-1,1,-1
                j_in = j_in + 1
                if(i_in > self_in%ldim(1))exit
                self_out%rmat(:self_out%ldim(1),j,1) = self_out%rmat(:self_out%ldim(1),j_in,1)
            enddo
            ! lower border & corners
            j_in = stops(2)+1
            do j = stops(2)+1, self_out%ldim(2)
                j_in = j_in - 1
                if(j_in < 1)exit
                self_out%rmat(:self_out%ldim(1),j,1) = self_out%rmat(:self_out%ldim(1),j_in,1)
            enddo
            self_out%ft = .false.
        else
            THROW_HARD('inconsistent dimensions; pad_mirr')
        endif
    end subroutine pad_mirr_1

    module subroutine pad_mirr_2( self, ldim )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: ldim(3)
        type(image) :: tmp
        call tmp%new(ldim , self%smpd)
        call self%pad_mirr_1(tmp)
        call self%copy(tmp)
        call tmp%kill
    end subroutine pad_mirr_2

    ! DO NOT PARALLELISE
    module subroutine clip( self_in, self_out )
        class(image), intent(inout) :: self_in, self_out
        real                        :: ratio
        integer                     :: starts(3), stops(3), lims(3,2)
        integer                     :: phys_out(3), phys_in(3), h,k,l,kpi,kpo,hp
        if( self_out%ldim(1) <= self_in%ldim(1) .and. self_out%ldim(2) <= self_in%ldim(2)&
        .and. self_out%ldim(3) <= self_in%ldim(3) )then
            if( self_in%ft )then
                lims = self_out%fit%loop_lims(2)
                if( self_in%is_2d() )then
                    do k=lims(2,1),lims(2,2)
                        kpi = k + 1 + merge(self_in%ldim(2) ,0,k<0)
                        kpo = k + 1 + merge(self_out%ldim(2),0,k<0)
                        do h=lims(1,1),lims(1,2)
                            hp = h+1
                            self_out%cmat(hp,kpo,1) = self_in%cmat(hp,kpi,1)
                        end do
                    end do
                else
                    do l=lims(3,1),lims(3,2)
                        do k=lims(2,1),lims(2,2)
                            do h=lims(1,1),lims(1,2)
                                phys_out = self_out%fit%comp_addr_phys(h,k,l)
                                phys_in = self_in%fit%comp_addr_phys(h,k,l)
                                self_out%cmat(phys_out(1),phys_out(2),phys_out(3)) =&
                                self_in%cmat(phys_in(1),phys_in(2),phys_in(3))
                            end do
                        end do
                    end do
                endif
                ratio = real(self_in%ldim(1))/real(self_out%ldim(1))
                call self_out%set_smpd(self_in%smpd*ratio) ! clipping Fourier transform, so sampling is coarser
                self_out%ft = .true.
            else
                starts = (self_in%ldim-self_out%ldim)/2+1
                stops  = self_in%ldim-starts+1
                if( self_in%ldim(3) == 1 )then
                    starts(3) = 1
                    stops(3)  = 1
                endif
                self_out%rmat(:self_out%ldim(1),:self_out%ldim(2),:self_out%ldim(3))&
                = self_in%rmat(starts(1):stops(1),starts(2):stops(2),starts(3):stops(3))
                self_out%ft = .false.
            endif
        endif
    end subroutine clip

    module subroutine clip_inplace( self, ldim )
        class(image), intent(inout) :: self
        integer,      intent(in)    :: ldim(3)
        type(image) :: tmp
        call tmp%new(ldim, self%smpd, wthreads=self%wthreads)
        call self%clip(tmp)
        call self%copy(tmp)
        call tmp%kill()
    end subroutine clip_inplace

    module subroutine read_and_crop( self, volfname, smpd, box_crop, smpd_crop )
        class(image),   intent(inout) :: self
        class(string),  intent(in)    :: volfname
        integer,        intent(in)    :: box_crop
        real,           intent(in)    :: smpd, smpd_crop
        integer :: ldim(3), ifoo, box
        call find_ldim_nptcls(volfname, ldim, ifoo)
        ! HE, I would not trust the smpd from the header
        if( ldim(3) /= ldim(1) ) THROW_HARD('Only for volumes')
        box = ldim(1)
        call self%new(ldim, smpd)
        call self%read(volfname)
        if( box < box_crop )then
            ! pad
            call self%fft
            call self%pad_inplace([box_crop,box_crop,box_crop], antialiasing=.false.)
            call self%ifft
        else if( box > box_crop )then
            ! clip
            call self%fft
            call self%clip_inplace([box_crop,box_crop,box_crop])
            call self%ifft
        endif
        call self%set_smpd(smpd_crop) ! safety
    end subroutine read_and_crop

    ! flip/mirror/rotate/shift

    module recursive subroutine flip( self, mode )
        class(image),     intent(inout) :: self
        character(len=*), intent(in)    :: mode
        integer :: i,j,ii,jj
        real    :: val
        if( self%ldim(3) > 1 )THROW_HARD('for 2D images only; flip')
        select case(upperCase(trim(mode)))
        case('X')
            !$omp parallel do private(i,j,ii,val) proc_bind(close) default(shared) schedule(static)
            do i = 1,self%ldim(1)/2
                ii = self%ldim(1) - i - 1
                do j = 1,self%ldim(2)
                    val               = self%rmat(i,j,1)
                    self%rmat(i,j,1)  = self%rmat(ii,j,1)
                    self%rmat(ii,j,1) = val
                enddo
            enddo
            !$omp end parallel do
        case('Y')
            !$omp parallel do private(i,j,jj,val) proc_bind(close) default(shared) schedule(static)
            do j = 1,self%ldim(2)/2
                jj = self%ldim(2) - j - 1
                do i = 1,self%ldim(1)
                    val               = self%rmat(i,j,1)
                    self%rmat(i,j,1)  = self%rmat(i,jj,1)
                    self%rmat(i,jj,1) = val
                enddo
            enddo
            !$omp end parallel do
        case('XY','YX')
            call self%flip('X')
            call self%flip('Y')
        case DEFAULT
            THROW_HARD('Unsupported flip mode')
        end select
    end subroutine flip

    module subroutine mirror( self, md, fourier )
        class(image),      intent(inout) :: self
        character(len=*),  intent(in)    :: md
        logical, optional, intent(in)    :: fourier
        integer :: i, j
        logical :: didft, l_fourier
        didft = .false.
        if( self%ft )then
            call self%ifft()
            didft = .true.
        endif
        l_fourier = .false.
        if( present(fourier) ) l_fourier = fourier
        if( md == 'x' )then
            do i=1,self%ldim(2)
                do j=1,self%ldim(3)
                    if( l_fourier )then
                        call reverse_f(self%rmat(1:self%ldim(1),i,j))
                    else
                        call reverse(self%rmat(1:self%ldim(1),i,j))
                    endif
                end do
            end do
        else if( md == 'y' )then
            do i=1,self%ldim(1)
                do j=1,self%ldim(3)
                    if( l_fourier )then
                        call reverse_f(self%rmat(i,1:self%ldim(2),j))
                    else
                        call reverse(self%rmat(i,1:self%ldim(2),j))
                    endif
                end do
            end do
        else if( md == 'z' )then
            do i=1,self%ldim(1)
                do j=1,self%ldim(2)
                    if( l_fourier )then
                        call reverse_f(self%rmat(i,j,1:self%ldim(3)))
                    else
                        call reverse(self%rmat(i,j,1:self%ldim(3)))
                    endif
                end do
            end do
        else
            write(logfhandle,'(a)') 'Mode needs to be either x, y or z; mirror; simple_image'
        endif
        if( didft ) call self%fft()
    end subroutine mirror

    module function calc_shiftcen( self, lp, msk, hp ) result( xyz )
        class(image),   intent(inout) :: self
        real,           intent(in)    :: lp
        real, optional, intent(in)    :: msk, hp
        type(image) :: tmp
        real        :: xyz(3), rmsk
        if( present(msk) )then
            rmsk = msk
        else
            rmsk = real( self%ldim(1) )/2. - 5. ! 5 pixels outer width
        endif
        call tmp%copy(self)
        if( present(hp) )then
            call tmp%bp(hp, lp)
        else
            call tmp%bp(0., lp)
        endif
        call tmp%ifft()
        call tmp%mask(rmsk, 'hard')
        ! such that norm_minmax will neglect everything < 0. and preserve zero
        where(tmp%rmat < TINY) tmp%rmat=0.
        call tmp%norm_minmax
        call tmp%masscen(xyz)
        call tmp%kill
    end function calc_shiftcen

    module function calc_shiftcen_serial( self, lp, msk, hp ) result( xyz )
        class(image),   intent(inout) :: self
        real,           intent(in)    :: lp
        real,           intent(in)    :: msk
        real, optional, intent(in)    :: hp
        real    :: xyz(3)
        integer :: ithr
        ! get thread index
        ithr = omp_get_thread_num() + 1
        if( all(self%ldim == thread_safe_tmp_imgs(ithr)%ldim) )then
            ! copy rmat
            thread_safe_tmp_imgs(ithr)%rmat = self%rmat
            thread_safe_tmp_imgs(ithr)%ft   = .false.
            call thread_safe_tmp_imgs(ithr)%fft()
            if( present(hp) )then
                call thread_safe_tmp_imgs(ithr)%bp(hp, lp)
            else
                call thread_safe_tmp_imgs(ithr)%bp(0., lp)
            endif
            call thread_safe_tmp_imgs(ithr)%ifft()
            call thread_safe_tmp_imgs(ithr)%mask(msk, 'hard')
            where(thread_safe_tmp_imgs(ithr)%rmat < TINY) thread_safe_tmp_imgs(ithr)%rmat = 0.
            call thread_safe_tmp_imgs(ithr)%norm_minmax
            call thread_safe_tmp_imgs(ithr)%masscen(xyz)
        else
            THROW_HARD('Incompatible dimensions bwetween self and thread_safe_tmp_imgs; calc_shiftcen_serial')
        endif
    end function calc_shiftcen_serial

    module subroutine roavg( self, angstep, avg, ang_stop )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: angstep
        class(image),      intent(inout) :: avg
        integer, optional, intent(in)    :: ang_stop
        real(dp):: avgs_rmat(nthr_glob,self%ldim(1),self%ldim(2),1)
        real    :: rotated(nthr_glob,self%ldim(1),self%ldim(2),1)
        integer :: irot, ithr, aang_stop
        aang_stop = 359
        if( present(ang_stop) ) aang_stop = ang_stop
        avgs_rmat = 0._dp
        rotated   = 0.
        !$omp parallel do schedule(static) default(shared) private(irot,ithr) proc_bind(close)
        do irot = 0 + angstep,aang_stop,angstep
            ! get thread index
            ithr = omp_get_thread_num() + 1
            ! rotate & sum
            call self%rtsq_serial(real(irot), 0., 0., rotated(ithr,:,:,:))
            avgs_rmat(ithr,:,:,:) = avgs_rmat(ithr,:,:,:) + real(rotated(ithr,:,:,:), dp)
        end do
        !$omp end parallel do
        ! add in the zero rotation
        avgs_rmat(1,:,:,:) = avgs_rmat(1,:,:,:) + real(self%rmat(:self%ldim(1),:self%ldim(2),:), dp)
        ! normalize and set output image object
        call avg%new(self%ldim, self%smpd)
        call avg%set_rmat(real(sum(avgs_rmat, dim=1)/real(360/angstep,dp)),.false.)
    end subroutine roavg

    !> \brief rtsq  rotation of image by quadratic interpolation (from spider)
    module subroutine rtsq( self_in, ang, shxi, shyi, self_out )
        class(image),           intent(inout) :: self_in
        real,                   intent(in)    :: ang,shxi,shyi
        class(image), optional, intent(inout) :: self_out
        real    :: shx,shy,ry1,rx1,ry2,rx2,cod,sid,xi,fixcenmshx,fiycenmshy
        real    :: rye2,rye1,rxe2,rxe1,yi,ycod,ysid,yold,xold
        integer :: iycen,ixcen,ix,iy
        real    :: mat_in(self_in%ldim(1),self_in%ldim(2))
        real    :: mat_out(self_in%ldim(1),self_in%ldim(2))
        logical :: didft
        if( self_in%ldim(3) > 1 )         THROW_HARD('only for 2D images; rtsq')
        if( .not. self_in%square_dims() ) THROW_HARD('only for square dims; rtsq;')
        didft = .false.
        if( self_in%ft )then
            call self_in%ifft()
            didft = .true.
        endif
        mat_out = 0. ! this is necessary, because it bugs out if I try to use the 3D matrix
        mat_in = self_in%rmat(:self_in%ldim(1),:self_in%ldim(2),1)
        ! shift within image boundary
        shx = amod(shxi,float(self_in%ldim(1)))
        shy = amod(shyi,float(self_in%ldim(2)))
        ! spider image center
        ixcen = self_in%ldim(1)/2+1
        iycen = self_in%ldim(2)/2+1
        ! image dimensions around origin
        rx1 = -self_in%ldim(1)/2
        rx2 =  self_in%ldim(1)/2
        ry1 = -self_in%ldim(2)/2
        ry2 =  self_in%ldim(2)/2
        rye1 = 0.
        rye2 = 0.
        if(mod(self_in%ldim(1),2) == 0)then
            rx2  =  rx2-1.0
            rxe1 = -self_in%ldim(1)
            rxe2 =  self_in%ldim(1)
        else
            rxe1 = -self_in%ldim(1)-1
            rxe2 =  self_in%ldim(1)+1
        endif
        if(mod(self_in%ldim(2),2) == 0)then
            ry2  =  ry2-1.0
            rye1 = -self_in%ldim(2)
            rye2 =  self_in%ldim(2)
        else
            ry2  = -self_in%ldim(2)-1
            rye2 =  self_in%ldim(2)+1
        endif
        ! create transformation matrix
        cod = cos(deg2rad(ang))
        sid = sin(deg2rad(ang))
        !-(center plus shift)
        fixcenmshx = -ixcen-shx
        fiycenmshy = -iycen-shy
        !$omp parallel do default(shared) private(iy,yi,ycod,ysid,ix,xi,xold,yold)&
        !$omp schedule(static) proc_bind(close)
        do iy=1,self_in%ldim(2)
            yi = iy+fiycenmshy
            if(yi < ry1) yi = min(yi+rye2, ry2)
            if(yi > ry2) yi = max(yi+rye1, ry1)
            ycod =  yi*cod+iycen
            ysid = -yi*sid+ixcen
            do ix=1,self_in%ldim(1)
                xi = ix+fixcenmshx
                if(xi < rx1) xi = min(xi+rxe2, rx2)
                if(xi > rx2) xi = max(xi+rxe1, rx1)
                yold = xi*sid+ycod
                xold = xi*cod+ysid
                mat_out(ix,iy) = quadri(xold,yold,mat_in,self_in%ldim(1),self_in%ldim(2))
            enddo
        enddo
        !$omp end parallel do
        if( present(self_out) )then
            if( (self_in.eqdims.self_out) .and. (self_in.eqsmpd.self_out) )then
                ! no need to update dimensions
            else
                call self_out%copy(self_in)
            endif
            self_out%rmat(:self_in%ldim(1),:self_in%ldim(2),1) = mat_out
            self_out%ft = .false.
        else
            self_in%rmat(:self_in%ldim(1),:self_in%ldim(2),1) = mat_out
            self_in%ft = .false.
        endif
        if( didft )then
            call self_in%fft()
        endif
    end subroutine rtsq

    !> \brief rtsq  rotation of image by quadratic interpolation (from spider)
    module subroutine rtsq_serial( self_in, ang, shxi, shyi, rmat_out )
        class(image), intent(inout) :: self_in
        real,         intent(in)    :: ang,shxi,shyi
        real,         intent(inout) :: rmat_out(self_in%ldim(1),self_in%ldim(2),1)
        real    :: shx,shy,ry1,rx1,ry2,rx2,cod,sid,xi,fixcenmshx,fiycenmshy
        real    :: rye2,rye1,rxe2,rxe1,yi,ycod,ysid,yold,xold
        integer :: iycen,ixcen,ix,iy
        real    :: mat_in(self_in%ldim(1),self_in%ldim(2))
        mat_in = self_in%rmat(:self_in%ldim(1),:self_in%ldim(2),1)
        ! shift within image boundary
        shx = amod(shxi,float(self_in%ldim(1)))
        shy = amod(shyi,float(self_in%ldim(2)))
        ! spider image center
        ixcen = self_in%ldim(1)/2+1
        iycen = self_in%ldim(2)/2+1
        ! image dimensions around origin
        rx1 = -self_in%ldim(1)/2
        rx2 =  self_in%ldim(1)/2
        ry1 = -self_in%ldim(2)/2
        ry2 =  self_in%ldim(2)/2
        rye1 = 0.
        rye2 = 0.
        if(mod(self_in%ldim(1),2) == 0)then
            rx2  =  rx2-1.0
            rxe1 = -self_in%ldim(1)
            rxe2 =  self_in%ldim(1)
        else
            rxe1 = -self_in%ldim(1)-1
            rxe2 =  self_in%ldim(1)+1
        endif
        if(mod(self_in%ldim(2),2) == 0)then
            ry2  =  ry2-1.0
            rye1 = -self_in%ldim(2)
            rye2 =  self_in%ldim(2)
        else
            ry2  = -self_in%ldim(2)-1
            rye2 =  self_in%ldim(2)+1
        endif
        ! create transformation matrix
        cod = cos(deg2rad(ang))
        sid = sin(deg2rad(ang))
        !-(center plus shift)
        fixcenmshx = -ixcen-shx
        fiycenmshy = -iycen-shy
        do iy=1,self_in%ldim(2)
            yi = iy+fiycenmshy
            if(yi < ry1) yi = min(yi+rye2, ry2)
            if(yi > ry2) yi = max(yi+rye1, ry1)
            ycod =  yi*cod+iycen
            ysid = -yi*sid+ixcen
            do ix=1,self_in%ldim(1)
                xi = ix+fixcenmshx
                if(xi < rx1) xi = min(xi+rxe2, rx2)
                if(xi > rx2) xi = max(xi+rxe1, rx1)
                yold = xi*sid+ycod
                xold = xi*cod+ysid
                rmat_out(ix,iy,1) = quadri(xold,yold,mat_in,self_in%ldim(1),self_in%ldim(2))
            enddo
        enddo
    end subroutine rtsq_serial

    !>  \brief  an image shifter to prepare for Fourier transformation
    module subroutine shift_phorig( self )
        class(image), intent(inout) :: self
        complex(sp) :: cswap
        real(sp)    :: rswap
        integer     :: i, j, k, h, ok, ii,jj, kfrom, kto, koffset
        if( .not.self%even_dims() )then
            write(logfhandle,*) 'ldim: ', self%ldim
            THROW_HARD('even dimensions assumed; shift_phorig')
        endif
        if( self%ft )then
            if( self%ldim(3) == 1 )then
                ! serial for now
                koffset = self%array_shape(2) / 2
                do k = 1,koffset
                    ok = k + koffset
                    do h = 1,self%array_shape(1)
                        cswap             = self%cmat(h, k,1)
                        self%cmat(h, k,1) = self%cmat(h,ok,1)
                        self%cmat(h,ok,1) = cswap
                    end do
                end do
            else
                THROW_HARD('3D FT not supported yet; shift_phorig')
            endif
        else
            if( self%ldim(3) == 1 )then
                if( self%wthreads )then
                    !$omp parallel do collapse(2) default(shared) private(rswap,i,j)&
                    !$omp schedule(static) proc_bind(close)
                    do i=1,self%ldim(1)/2
                        do j=1,self%ldim(2)/2
                            !(1)
                            rswap = self%rmat(i,j,1)
                            self%rmat(i,j,1) = self%rmat(self%ldim(1)/2+i,self%ldim(2)/2+j,1)
                            self%rmat(self%ldim(1)/2+i,self%ldim(2)/2+j,1) = rswap
                            !(2)
                            rswap = self%rmat(i,self%ldim(2)/2+j,1)
                            self%rmat(i,self%ldim(2)/2+j,1) = self%rmat(self%ldim(1)/2+i,j,1)
                            self%rmat(self%ldim(1)/2+i,j,1) = rswap
                        end do
                    end do
                    !$omp end parallel do
                else
                    do j=1,self%ldim(2)/2
                        jj = self%ldim(2)/2+j
                        do i=1,self%ldim(1)/2
                            ii = self%ldim(1)/2+i
                            !(1)
                            rswap = self%rmat(i,j,1)
                            self%rmat(i,j,1) = self%rmat(ii,jj,1)
                            self%rmat(ii,jj,1) = rswap
                            !(2)
                            rswap = self%rmat(i,jj,1)
                            self%rmat(i,jj,1) = self%rmat(ii,j,1)
                            self%rmat(ii,j,1) = rswap
                        end do
                    end do
                endif
            else
                kfrom = 1
                kto   = self%ldim(3)/2
                if( self%wthreads )then
                    !$omp parallel do collapse(3) default(shared) private(rswap,i,j,k)&
                    !$omp schedule(static) proc_bind(close)
                    do i=1,self%ldim(1)/2
                        do j=1,self%ldim(2)/2
                            do k=1,kto
                                !(1)
                                rswap = self%rmat(i,j,k)
                                self%rmat(i,j,k) = self%rmat(self%ldim(1)/2+i,self%ldim(2)/2+j,self%ldim(3)/2+k)
                                self%rmat(self%ldim(1)/2+i,self%ldim(2)/2+j,self%ldim(3)/2+k) = rswap
                                !(2)
                                rswap = self%rmat(i,self%ldim(2)/2+j,self%ldim(3)/2+k)
                                self%rmat(i,self%ldim(2)/2+j,self%ldim(3)/2+k) = self%rmat(self%ldim(1)/2+i,j,k)
                                self%rmat(self%ldim(1)/2+i,j,k) = rswap
                                !(3)
                                rswap = self%rmat(self%ldim(1)/2+i,j,self%ldim(3)/2+k)
                                self%rmat(self%ldim(1)/2+i,j,self%ldim(3)/2+k) = self%rmat(i,self%ldim(2)/2+j,k)
                                self%rmat(i,self%ldim(2)/2+j,k) = rswap
                                !(4)
                                rswap = self%rmat(i,j,self%ldim(3)/2+k)
                                self%rmat(i,j,self%ldim(3)/2+k) = self%rmat(self%ldim(1)/2+i,self%ldim(2)/2+j,k)
                                self%rmat(self%ldim(1)/2+i,self%ldim(2)/2+j,k) = rswap
                            end do
                        end do
                    end do
                    !$omp end parallel do
                else
                    do i=1,self%ldim(1)/2
                        do j=1,self%ldim(2)/2
                            do k=1,kto
                                !(1)
                                rswap = self%rmat(i,j,k)
                                self%rmat(i,j,k) = self%rmat(self%ldim(1)/2+i,self%ldim(2)/2+j,self%ldim(3)/2+k)
                                self%rmat(self%ldim(1)/2+i,self%ldim(2)/2+j,self%ldim(3)/2+k) = rswap
                                !(2)
                                rswap = self%rmat(i,self%ldim(2)/2+j,self%ldim(3)/2+k)
                                self%rmat(i,self%ldim(2)/2+j,self%ldim(3)/2+k) = self%rmat(self%ldim(1)/2+i,j,k)
                                self%rmat(self%ldim(1)/2+i,j,k) = rswap
                                !(3)
                                rswap = self%rmat(self%ldim(1)/2+i,j,self%ldim(3)/2+k)
                                self%rmat(self%ldim(1)/2+i,j,self%ldim(3)/2+k) = self%rmat(i,self%ldim(2)/2+j,k)
                                self%rmat(i,self%ldim(2)/2+j,k) = rswap
                                !(4)
                                rswap = self%rmat(i,j,self%ldim(3)/2+k)
                                self%rmat(i,j,self%ldim(3)/2+k) = self%rmat(self%ldim(1)/2+i,self%ldim(2)/2+j,k)
                                self%rmat(self%ldim(1)/2+i,self%ldim(2)/2+j,k) = rswap
                            end do
                        end do
                    end do
                endif
            endif
        endif
    end subroutine shift_phorig

    !> \brief shift  is for origin shifting an image
    !! \param x position in axis 0
    !! \param y position in axis 1
    !! \param z position in axis 2
    module subroutine shift( self, shvec )
        class(image), intent(inout) :: self
        real,         intent(in)    :: shvec(3)
        integer :: h, k, l, lims(3,2), phys(3)
        real    :: shvec_here(3)
        logical :: didft
        shvec_here = shvec
        if( self%ldim(3) == 1 ) shvec_here(3) = 0.0
        didft = .false.
        if( .not. self%ft )then
            call self%fft()
            didft = .true.
        endif
        lims = self%fit%loop_lims(2)
        !$omp parallel do collapse(3) default(shared) private(phys,h,k,l)&
        !$omp schedule(static) proc_bind(close)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    phys = self%fit%comp_addr_phys([h,k,l])
                    self%cmat(phys(1),phys(2),phys(3)) = self%cmat(phys(1),phys(2),phys(3))*&
                    self%oshift([h,k,l], shvec_here)
                end do
            end do
        end do
        !$omp end parallel do
        if( didft ) call self%ifft()
    end subroutine shift

    module subroutine shift2Dserial_1( self, shvec  )
        class(image), intent(inout) :: self
        real,         intent(in)    :: shvec(2)
        real(dp),allocatable :: hcos(:), hsin(:)
        real(dp) :: sh(2), arg, ck, sk
        integer  :: h,k, hphys,kphys, lims(3,2)
        lims = self%fit%loop_lims(2)
        sh   = real(shvec * self%shconst(1:2),dp)
        allocate(hcos(lims(1,1):lims(1,2)),hsin(lims(1,1):lims(1,2)))
        do h=lims(1,1),lims(1,2)
            arg = real(h,dp)*sh(1)
            hcos(h) = dcos(arg)
            hsin(h) = dsin(arg)
        enddo
        do k=lims(2,1),lims(2,2)
            kphys = k + 1 + merge(self%ldim(2),0,k<0)
            arg = real(k,dp)*sh(2)
            ck  = dcos(arg)
            sk  = dsin(arg)
            do h=lims(1,1),lims(1,2)
                hphys = h + 1
                self%cmat(hphys,kphys,1) = self%cmat(hphys,kphys,1)&
                    &* cmplx(ck*hcos(h)-sk*hsin(h), ck*hsin(h)+sk*hcos(h),sp)
            end do
        end do
    end subroutine shift2Dserial_1

    module subroutine shift2Dserial_2( self, shvec, self_out )
        class(image), intent(inout) :: self, self_out
        real,         intent(in)    :: shvec(2)
        real(dp),allocatable :: hcos(:), hsin(:)
        real(dp) :: sh(2), arg, ck, sk
        integer  :: h,k, hphys,kphys, lims(3,2)
        lims = self%fit%loop_lims(2)
        sh   = real(shvec * self%shconst(1:2),dp)
        allocate(hcos(lims(1,1):lims(1,2)),hsin(lims(1,1):lims(1,2)))
        do h=lims(1,1),lims(1,2)
            arg = real(h,dp)*sh(1)
            hcos(h) = dcos(arg)
            hsin(h) = dsin(arg)
        enddo
        do k=lims(2,1),lims(2,2)
            kphys = k + 1 + merge(self%ldim(2),0,k<0)
            arg = real(k,dp)*sh(2)
            ck  = dcos(arg)
            sk  = dsin(arg)
            do h=lims(1,1),lims(1,2)
                hphys = h + 1
                self_out%cmat(hphys,kphys,1) = self%cmat(hphys,kphys,1)&
                    &* cmplx(ck*hcos(h)-sk*hsin(h), ck*hsin(h)+sk*hcos(h),sp)
            end do
        end do
    end subroutine shift2Dserial_2

    module subroutine masscen( self, xyz, mask_in )
        class(image),      intent(inout) :: self
        real        ,      intent(out)   :: xyz(3)
        logical, optional, intent(in)    :: mask_in(:,:,:)
        real                 ::  spix, ci, cj, ck
        integer              :: i, j, k
        logical, allocatable :: mask_here(:,:,:)
        allocate(mask_here(self%ldim(1),self%ldim(2),self%ldim(3)), source=.true.)
        if (present(mask_in)) then 
            if (.not.(size(mask_in, dim=1) .eq. self%ldim(1))) THROW_HARD('mask_in dimension must match dimension of image')
            if (.not.(size(mask_in, dim=2) .eq. self%ldim(2))) THROW_HARD('mask_in dimension must match dimension of image')
            if (.not.(size(mask_in, dim=3) .eq. self%ldim(3))) THROW_HARD('mask_in dimension must match dimension of image')
            mask_here = mask_in   
        end if
        if( self%is_ft() ) THROW_HARD('masscen not implemented for FTs; masscen')
        spix = 0.
        xyz  = 0.
        ci   = -real(self%ldim(1))/2.
        do i=1,self%ldim(1)
            cj = -real(self%ldim(2))/2.
            do j=1,self%ldim(2)
                ck = -real(self%ldim(3))/2.
                do k=1,self%ldim(3)
                    if (mask_here(i,j,k)) then
                        xyz  = xyz  + self%rmat(i,j,k) * [ci, cj, ck]
                        spix = spix + self%rmat(i,j,k) 
                    end if
                    ck = ck + 1.
                end do
                cj = cj + 1.
            end do
            ci = ci + 1.
        end do
        if( is_equal(spix,0.) ) return
        xyz = xyz / spix
        if( self%ldim(3) == 1 ) xyz(3) = 0.
    end subroutine masscen

    module subroutine masscen_adjusted( self, xyz, mask_in )
        class(image),      intent(inout) :: self
        real        ,      intent(out)   :: xyz(3)
        logical, optional, intent(in)    :: mask_in(:,:,:)
        real                 ::  spix
        integer              :: i, j, k
        logical, allocatable :: mask_here(:,:,:)
        allocate(mask_here(self%ldim(1),self%ldim(2),self%ldim(3)), source=.true.)
        if (present(mask_in)) then 
            if (.not.(size(mask_in, dim=1) .eq. self%ldim(1))) THROW_HARD('mask_in dimension must match dimension of image')
            if (.not.(size(mask_in, dim=2) .eq. self%ldim(2))) THROW_HARD('mask_in dimension must match dimension of image')
            if (.not.(size(mask_in, dim=3) .eq. self%ldim(3))) THROW_HARD('mask_in dimension must match dimension of image')
            mask_here = mask_in   
        end if
        if( self%is_ft() ) THROW_HARD('masscen not implemented for FTs; masscen')
        spix = 0.
        xyz  = 0.
        do i=1,self%ldim(1)
            do j=1,self%ldim(2)
                do k=1,self%ldim(3)
                    if (mask_here(i,j,k)) then
                        xyz  = xyz  + self%rmat(i,j,k) * [i, j, k]
                        spix = spix + self%rmat(i,j,k) 
                    end if
                end do
            end do
        end do
        if( is_equal(spix,0.) ) return
        xyz = xyz / spix
        if( self%ldim(3) == 1 ) xyz(3) = 0.
    end subroutine masscen_adjusted

end submodule simple_image_geom
