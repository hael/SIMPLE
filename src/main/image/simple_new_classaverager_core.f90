!@descr: Implementations of core types and routines underpinning cavg restoration
submodule (simple_new_classaverager) simple_new_classaverager_core
implicit none
#include "simple_local_flags.inc"

contains

    ! Type cavgs_set

    module subroutine new_set( self, ldim, ncls )
        class(cavgs_set), intent(inout) :: self
        integer,          intent(in)    :: ldim(2), ncls
        call self%kill_set
        self%ldim = ldim
        self%ncls = ncls
        call self%even%new_stack(  ldim, self%ncls)
        call self%odd%new_stack(   ldim, self%ncls)
        call self%merged%new_stack(ldim, self%ncls)
    end subroutine new_set

    module subroutine zero_set( self, ft )
        class(cavgs_set), intent(inout) :: self
        logical,          intent(in)    :: ft
        call self%even%zero(ft)
        call self%odd%zero(ft)
        call self%merged%zero(ft)
    end subroutine zero_set

    module subroutine copy_fast( self, self2copy, is, ft )
        class(cavgs_set), intent(inout) :: self
        class(cavgs_set), intent(in)    :: self2copy
        integer,         intent(in)    :: is
        logical,         intent(in)    :: ft
        self%even%cmat(:,:,is)    = self2copy%even%cmat(:,:,is)
        self%odd%cmat(:,:,is)     = self2copy%odd%cmat(:,:,is)
        self%merged%cmat(:,:,is)  = self2copy%merged%cmat(:,:,is)
        self%even%ctfsq(:,:,is)   = self2copy%even%ctfsq(:,:,is)
        self%odd%ctfsq(:,:,is)    = self2copy%odd%ctfsq(:,:,is)
        self%merged%ctfsq(:,:,is) = self2copy%merged%ctfsq(:,:,is)
        self%even%slices(is)%ft   = ft
        self%odd%slices(is)%ft    = ft
        self%merged%slices(is)%ft = ft
    end subroutine

    module subroutine kill_set( self )
        class(cavgs_set), intent(inout) :: self
        self%ldim = 0
        self%ncls = 0
        call self%even%kill_stack
        call self%odd%kill_stack
        call self%merged%kill_stack
    end subroutine kill_set

    ! Type stack

    module subroutine new_stack( self, ldim, nslices, alloc_ctfsq )
        class(stack),      intent(inout) :: self
        integer,           intent(in)    :: ldim(2), nslices
        logical, optional, intent(in)    :: alloc_ctfsq
        real             :: center(2),e, r
        integer          :: phys(2), h,i,j,k,sh,nyq,minlen
        logical          :: l_alloc_ctfsq
        call self%kill_stack
        l_alloc_ctfsq = .true.
        if( present(alloc_ctfsq) ) l_alloc_ctfsq = alloc_ctfsq
        self%ldim     = ldim
        self%cshape   = [fdim(self%ldim(1)), self%ldim(2)]  ! FFTW convention
        self%rshape   = [2*self%cshape(1), self%ldim(2)]    ! FFTW convention
        self%nslices  = nslices
        self%fit      = ftiter([self%ldim(1), self%ldim(2), 1], 1.0)
        self%flims    = self%fit%loop_lims(2)
        ! Array allocation
        self%p = fftwf_alloc_complex(int(self%nslices*product(self%cshape),c_size_t))
        call c_f_pointer(self%p, self%cmat, [self%cshape(1), self%cshape(2), self%nslices])
        call c_f_pointer(self%p, self%rmat, [self%rshape(1), self%rshape(2), self%nslices])
        self%rmat = 0.0
        ! Allocation & association of pointers to slices
        allocate(self%slices(nslices))
        do i = 1,self%nslices
            self%slices(i)%c => self%cmat(:,:,i)
            self%slices(i)%r => self%rmat(:,:,i)
        enddo
        self%slices(:)%ft = .false. ! real space
        ! fftw plans
        self%plan_fwd = fftwf_plan_dft_r2c_2d(self%ldim(2), self%ldim(1),&
                        &self%slices(1)%r, self%slices(1)%c, FFTW_ESTIMATE)
        self%plan_bwd = fftwf_plan_dft_c2r_2d(self%ldim(2), self%ldim(1),&
                        &self%slices(1)%c, self%slices(1)%r, FFTW_ESTIMATE)
        ! CTF^2 array
        allocate(self%ctfsq(self%cshape(1),self%cshape(2),self%nslices))
        ! Initialise arrays
        call self%zero(.false.)
        ! Nyqist logical mask
        allocate(self%nyq_mask(self%cshape(1), self%cshape(2)), source=.false.)
        nyq = self%fit%get_lfny(1)
        do k = self%flims(2,1),self%flims(2,2)
            do h = self%flims(1,1),self%flims(1,2)
                sh = nint(hyp(h,k))
                if( sh > nyq ) cycle
                phys = self%fit%comp_addr_phys(h,k)
                self%nyq_mask(phys(1),phys(2)) = .true.
            end do
        end do
        ! real space soft mask, squared image assumed
        allocate(self%soft_mask(self%ldim(1),self%ldim(2)),source=1.0)
        center = self%ldim/2 + 1
        do j = 1,self%ldim(2)
            do i = 1,self%ldim(1)
                r = sqrt((i-center(1))**2 + (j-center(2))**2)
                if( r < (params_glob%msk_crop-COSMSKHALFWIDTH) )then
                    cycle
                else if( r > (params_glob%msk_crop+COSMSKHALFWIDTH) )then
                    e = 0.
                else
                    r = (r-(params_glob%msk_crop-COSMSKHALFWIDTH)) / (2.*COSMSKHALFWIDTH)
                    r = min(1.0, max(0.0, r))
                    e = min(1.0, max(0.0, cos(r * pio2)))
                endif
                self%soft_mask(i, j) = e
            end do
        end do
    end subroutine new_stack

    module subroutine zero( self, ft )
        class(stack), intent(inout) :: self
        logical,      intent(in)    :: ft
        !$omp parallel workshare proc_bind(close)
        self%rmat  = 0.0
        self%ctfsq = 0.0
        self%slices(:)%ft = ft
        !$omp end parallel workshare
    end subroutine zero

    module subroutine zero_slice( self, is, ft )
        class(stack), intent(inout) :: self
        integer,      intent(in)    :: is
        logical,      intent(in)    :: ft
        self%rmat(:,:,is)  = 0.0
        self%ctfsq(:,:,is) = 0.0
        self%slices(is)%ft = ft
    end subroutine zero_slice

    module subroutine write( self, fname, ft )
        class(stack),  intent(inout) :: self
        class(string), intent(in)    :: fname
        logical,       intent(in)    :: ft
        type(imgfile) :: ioimg
        real          :: stats(4)
        integer       :: funit, ierr
        if( ft )then
            ! dump matrix on disk
            call fopen(funit, fname, status='UNKNOWN', action='WRITE', access='STREAM', iostat=ierr)
            call fileiochk( 'write: '//fname%to_char(), ierr )
            write(funit, pos=1, iostat=ierr) self%rmat
            if( ierr .ne. 0 ) call fileiochk('write: '//fname%to_char(), ierr)
            call fclose(funit)
        else
            ! used for proper cavgs output, header updated with stats
            call ioimg%open(fname, [self%ldim(1),self%ldim(2),1], smpd_crop, formatchar='M', readhead=.false., rwaction='write')
            call ioimg%wmrcSlices(1, self%nslices, self%rmat(:self%ldim(1),:,:), [self%ldim(1),self%ldim(2),1], .false.)
            call ioimg%setMode(2)
            call self%calc_cavgs_stats(stats)
            call ioimg%update_MRC_stats(stats)
            call ioimg%close
        endif
    end subroutine write

    module subroutine write_cmat_debug( self, fname )
        class(stack),  intent(inout) :: self
        class(string), intent(in)    :: fname
        type(imgfile) :: ioimg
        real          :: stats(4)
        integer       :: funit, ierr, i
        do i =1,self%nslices
            call self%ifft(i)  
        enddo
        ! used for proper cavgs output, header updated with stats
        call ioimg%open(fname, [self%ldim(1),self%ldim(2),1], smpd_crop, formatchar='M', readhead=.false., rwaction='write')
        call ioimg%wmrcSlices(1, self%nslices, self%rmat(:self%ldim(1),:,:), [self%ldim(1),self%ldim(2),1], .false.)
        call ioimg%setMode(2)
        call self%calc_cavgs_stats(stats)
        call ioimg%update_MRC_stats(stats)
        call ioimg%close
        do i =1,self%nslices
            call self%fft(i)  
        enddo
    end subroutine write_cmat_debug

    module subroutine write_ctfsq_debug( self, fname )
        class(stack),  intent(inout) :: self
        class(string), intent(in)    :: fname
        type(imgfile) :: ioimg
        real          :: stats(4)
        integer       :: funit, ierr, i
        ! used for proper cavgs output, header updated with stats
        call ioimg%open(fname, [self%cshape(1),self%cshape(2),1], smpd_crop, formatchar='M', readhead=.false., rwaction='write')
        call ioimg%wmrcSlices(1, self%nslices, self%ctfsq(:self%cshape(1),:self%cshape(2),:), [self%cshape(1),self%cshape(2),1], .false.)
        call ioimg%setMode(2)
        call ioimg%close
    end subroutine write_ctfsq_debug

    module subroutine write_ctfsq( self, fname )
        class(stack),  intent(inout) :: self
        class(string), intent(in)    :: fname
        integer :: funit, ierr
        call fopen(funit, fname, status='UNKNOWN', action='WRITE', access='STREAM', iostat=ierr)
        call fileiochk( 'write_ctfsq: '//fname%to_char(), ierr)
        write(funit, pos=1, iostat=ierr) self%ctfsq
        if( ierr .ne. 0 ) call fileiochk('write_ctfsq: '//fname%to_char(), ierr)
        call fclose(funit)
    end subroutine write_ctfsq

    module subroutine read_cmat( self, fname )
        class(stack),  intent(inout) :: self
        class(string), intent(in)    :: fname
        integer :: funit, ierr
        call fopen(funit, file=fname, status='OLD', action='READ', access='STREAM', iostat=ierr)
        call fileiochk('read_cmat '//fname%to_char(), ierr)
        read(funit, pos=1, iostat=ierr) self%rmat
        if( ierr .ne. 0 ) call fileiochk('read_cmat '//fname%to_char(), ierr)
        call fclose(funit)
        self%slices(:)%ft = .true.
    end subroutine read_cmat

    module subroutine read_ctfsq( self, fname )
        class(stack),  intent(inout) :: self
        class(string), intent(in)    :: fname
        integer :: funit, ierr
        call fopen(funit, file=fname, status='OLD', action='READ', access='STREAM', iostat=ierr)
        call fileiochk('read_ctfsq '//fname%to_char(), ierr)
        read(funit, pos=1, iostat=ierr) self%ctfsq
        if( ierr .ne. 0 ) call fileiochk('read_ctfsq '//fname%to_char(), ierr)
        call fclose(funit)
    end subroutine read_ctfsq

    ! to calculate stats necessary to the MRC header
    module subroutine calc_cavgs_stats( self, stats )
        class(stack), intent(in)  :: self
        real,         intent(out) :: stats(4)
        real(dp) :: avg,avgv,avgvsq, std,var, A,rn
        integer  :: is
        avg = 0.0d0; std = 0.0d0
        rn  = real(product(self%ldim),dp)
        A   = rn / (rn - 1.d0)
        !$omp parallel proc_bind(close) default(shared)
        !$omp workshare
        stats(1) = minval(self%rmat(:self%ldim(1),:self%ldim(2),:))
        stats(2) = maxval(self%rmat(:self%ldim(1),:self%ldim(2),:))
        !$omp end workshare nowait
        !$omp do schedule(static) private(is,avgv,avgvsq,var) reduction(+:avg,std)
        do is = 1,self%nslices
            avgv   = real(sum(self%rmat(:self%ldim(1),:,is)),   dp) / rn
            avgvsq = real(sum(self%rmat(:self%ldim(1),:,is)**2),dp) / rn
            avg    = avg + avgv
            var    = (avgvsq - avgv**2) * A
            if( var > DTINY ) std = std + sqrt(var)
        enddo
        !$omp end do
        !$omp end parallel
        stats(3) = real(avg) / real(self%nslices)
        stats(4) = real(std) / real(self%nslices)
    end subroutine calc_cavgs_stats

    module subroutine quadrant_swap( self, is )
        class(stack), intent(inout) :: self
        integer,      intent(in)    :: is
        real    :: r
        integer :: i,j,ii,jj
        do j=1,self%ldim(2)/2
            jj = self%ldim(2)/2+j
            do i=1,self%ldim(1)/2
                ii = self%ldim(1)/2+i
                r = self%rmat(i, j, is)
                self%rmat(i,  j, is)   = self%rmat(ii, jj, is)
                self%rmat(ii, jj, is) = r
                r = self%rmat(i, jj, is)
                self%rmat(i,  jj, is) = self%rmat(ii, j, is)
                self%rmat(ii,  j, is) = r
            end do
        end do
    end subroutine quadrant_swap

    ! Forward FT
    module subroutine fft( self, i )
        class(stack), intent(inout) :: self
        integer,      intent(in)    :: i
        if( self%slices(i)%ft ) return
        call self%quadrant_swap(i)
        call fftwf_execute_dft_r2c(self%plan_fwd, self%slices(i)%r, self%slices(i)%c)
        self%slices(i)%c  = self%slices(i)%c / real(product(self%ldim))
        self%slices(i)%ft = .true.
    end subroutine fft

    ! Backwards FT
    module subroutine ifft( self, i )
        class(stack), intent(inout) :: self
        integer,      intent(in)    :: i
        if( self%slices(i)%ft )then
            call fftwf_execute_dft_c2r(self%plan_bwd, self%slices(i)%c, self%slices(i)%r)
            call self%quadrant_swap(i)
            self%slices(i)%ft = .false.
        endif
    end subroutine ifft

    ! To calculate FRC between the e/o slices
    module subroutine frc( self1, self2, is, corrs )
        class(stack), intent(in)  :: self1, self2
        integer,      intent(in)  :: is
        real,         intent(out) :: corrs(fdim(self1%ldim(1))-1)
        real(dp)    :: corrs_8(fdim(self1%ldim(1))-1)
        real(dp)    :: sumasq(fdim(self1%ldim(1))-1), sumbsq(fdim(self1%ldim(1))-1)
        complex(dp) :: comp1, comp2
        integer     :: n, phys(2), sh, h,k
        corrs_8 = 0.d0
        sumasq  = 0.d0
        sumbsq  = 0.d0
        n       = fdim(self1%ldim(1))-1
        do k = self1%flims(2,1),self1%flims(2,2)
            do h = self1%flims(1,1),self1%flims(1,2)
                sh = nint(hyp(h,k))
                if( sh == 0 .or. sh > n ) cycle
                phys  = self1%fit%comp_addr_phys(h,k)
                comp1 = self1%cmat(phys(1), phys(2), is)
                comp2 = self2%cmat(phys(1), phys(2), is)
                corrs_8(sh) = corrs_8(sh) + real(comp1 * conjg(comp2),dp)
                sumasq(sh)  = sumasq(sh)  + real(comp1 * conjg(comp1),dp)
                sumbsq(sh)  = sumbsq(sh)  + real(comp2 * conjg(comp2),dp)
            end do
        end do
        do k = 1,n
            if( sumasq(k) > 0.d0 .and. sumbsq(k) > 0.d0 )then
                corrs(k) = real(corrs_8(k)/sqrt(sumasq(k) * sumbsq(k)))
            else
                corrs(k) = 0.0
            endif
        end do
    end subroutine frc

    ! Class average normalization: sum(CTF.I) / sum(CTF2)
    module subroutine ctf_dens_correct( self, is )
        class(stack), intent(inout) :: self
        integer,      intent(in)    :: is
        real    :: ctfsq
        integer :: i,j
        do j = 1,self%cshape(2)
            do i = 1,self%cshape(1)
                if( self%nyq_mask(i,j) )then
                    ctfsq = abs(self%ctfsq(i,j,is)) 
                    if( ctfsq > TINY ) self%cmat(i,j,is) = self%cmat(i,j,is) / ctfsq
                else
                    self%cmat(i,j,is) = CMPLX_ZERO
                endif
            enddo
        enddo
    end subroutine ctf_dens_correct

    !> Average edge norm subtraction followed by soft masking
    module subroutine softmask( self, is )
        class(stack), intent(inout) :: self
        integer,      intent(in)    :: is
        real(dp) :: edge_sum
        real     :: edge_avg
        edge_sum =            sum(self%rmat(1:self%ldim(1),1,is))   + sum(self%rmat(1:self%ldim(1),self%ldim(2),is))
        edge_sum = edge_sum + sum(self%rmat(1,2:self%ldim(2)-1,is)) + sum(self%rmat(self%ldim(1),2:self%ldim(2)-1,is))
        edge_avg = real(edge_sum / real(4*self%ldim(1)-4,dp))
        self%rmat(:self%ldim(1),:self%ldim(2),is) = self%soft_mask * (self%rmat(:self%ldim(1),:self%ldim(2),is) - edge_avg)
    end subroutine softmask

    ! To insert the low resolutions rings of even+odd into even and odd
    ! and keep these cavgs in register
    module subroutine insert_lowres_serial( self, self2insert, is, find )
        class(stack), intent(inout) :: self
        class(stack), intent(in)    :: self2insert
        integer,      intent(in)    :: is, find
        integer :: physh,physk, h,k, sh
        if( .not. self%slices(is)%ft        )THROW_HARD('image to be modified assumed to be FTed; insert_lowres 1')
        if( .not. self2insert%slices(is)%ft )THROW_HARD('image to be modified assumed to be FTed; insert_lowres 2')
        do h = self%flims(1,1),self%flims(1,2)
            do k = self%flims(2,1),self%flims(2,2)
                sh = nint(hyp(h,k))
                if( sh <= find )then
                    physh = ft_map_phys_addrh(h,k)
                    physk = ft_map_phys_addrk(h,k)
                    self%cmat(physh,physk,is) = self2insert%cmat(physh,physk,is)
                endif
            end do
        end do
    end subroutine insert_lowres_serial

    ! Adds noise term to denominator following ML regularization scheme
    module subroutine add_invnoisepower2rho( self, is, frcsz, frc )
        class(stack), intent(inout) :: self
        integer,      intent(in)    :: is, frcsz
        real,         intent(in)    :: frc(1:frcsz)
        real(dp) :: rsum(0:frcsz)
        real     :: ssnr(0:frcsz), tau2(0:frcsz), sig2(0:frcsz)
        integer  :: cnt(0:frcsz)
        real     :: cc, invtau2, fudge, ctfsq
        integer  :: phys(2), h, k, sh, reslim_ind
        fudge = params_glob%tau
        rsum  = 0.d0
        cnt   = 0
        ssnr  = 0.0; tau2 = 0.0; sig2 = 0.0
        ! SSNR
        do k = 1,frcsz
            cc      = min(0.999,max(0.001,frc(k)))
            ssnr(k) = fudge * cc / (1.-cc)
        enddo
        ! Noise
        do h = self%flims(1,1),self%flims(1,2)
            do k = self%flims(2,1),self%flims(2,2)
                sh = nint(hyp(h,k))
                if( sh > frcsz ) cycle                      ! use logical mask
                phys     = self%fit%comp_addr_phys(h, k)    ! replace with ft_map
                cnt(sh)  = cnt(sh) + 1
                rsum(sh) = rsum(sh) + real(self%ctfsq(phys(1), phys(2), is),dp)
            enddo
        enddo
        where( rsum > 1.d-10 ) sig2 = real(real(cnt,dp) / rsum)
        ! Signal power
        tau2 = ssnr * sig2

        print *, is, frc
        print *, is, ssnr
        print *, is, sig2
        print *, is, tau2

        ! add Tau2 inverse to denominator
        ! because signal assumed infinite at very low resolution there is no addition
        reslim_ind = max(6, calc_fourier_index(params_glob%hp, params_glob%box, params_glob%smpd))
        do h = self%flims(1,1),self%flims(1,2)
            do k = self%flims(2,1),self%flims(2,2)
                sh = nint(hyp(h,k))
                if( (sh < reslim_ind) .or. (sh > frcsz) ) cycle ! use logical mask
                phys  = self%fit%comp_addr_phys(h, k)           ! replace with ft_map
                ctfsq = self%ctfsq(phys(1), phys(2), is)
                if( tau2(sh) > TINY)then
                    invtau2 = 1.0 / (fudge*tau2(sh))
                else
                    invtau2 = min(1.e3, 1.e3*ctfsq)
                endif
                self%ctfsq(phys(1), phys(2), is) = ctfsq + invtau2
            enddo
        enddo
    end subroutine add_invnoisepower2rho

    module subroutine kill_stack( self )
        class(stack), intent(inout) :: self
        integer :: i
        if( allocated(self%slices) )then
            do i = 1,self%nslices
                nullify(self%slices(i)%r, self%slices(i)%c)
            enddo
            deallocate(self%slices,self%ctfsq,self%nyq_mask,self%soft_mask)
        endif
        if( c_associated(self%p) )then
            call fftwf_free(self%p)
            call fftwf_destroy_plan(self%plan_fwd)
            call fftwf_destroy_plan(self%plan_bwd)
            self%p        = c_null_ptr
            self%plan_fwd = c_null_ptr
            self%plan_bwd = c_null_ptr
            nullify(self%rmat,self%cmat)
        endif
        self%ldim    = 0
        self%rshape  = 0
        self%cshape  = 0
        self%flims   = 0
        self%nslices = 0
    end subroutine kill_stack

end submodule simple_new_classaverager_core
