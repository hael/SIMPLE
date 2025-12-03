submodule (simple_image) simple_image_core
!$ use omp_lib
!$ use omp_lib_kinds
include  'simple_lib.f08'
#include "simple_local_flags.inc"
implicit none
contains

    !========================================
    ! Constructors / lifecycle
    !========================================

    module subroutine new( self, ldim, smpd, wthreads )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: ldim(3)
        real,              intent(in)    :: smpd
        logical, optional, intent(in)    :: wthreads
        integer(kind=c_int) :: rc
        integer             :: i
        logical             :: do_allocate
        ! we need to be clever about allocation (because it is costly)
        if( self%existence )then
            if( any(self%ldim /= ldim) )then
                do_allocate = .true.
                call self%kill()
            else
                do_allocate = .false.
                !$omp critical
                call fftwf_destroy_plan(self%plan_fwd)
                call fftwf_destroy_plan(self%plan_bwd)
                !$omp end critical
            endif
        else
            do_allocate = .true.
        endif
        self%wthreads = .true.
        if( present(wthreads) ) self%wthreads = wthreads
        self%wthreads = self%wthreads .and. nthr_glob > 1
        self%ldim = ldim
        self%smpd = smpd
        ! Make Fourier iterator
        call self%fit%new(ldim, smpd)
        ! Work out dimensions of the complex array
        self%array_shape(1)   = fdim(self%ldim(1))
        self%array_shape(2:3) = self%ldim(2:3)
        self%nc = int(product(self%array_shape)) ! # components
        if( do_allocate )then
            ! Letting FFTW do the allocation in C ensures that we will be using aligned memory
            self%p = fftwf_alloc_complex(int(product(self%array_shape),c_size_t))
            ! Set up the complex array which will point at the allocated memory
            call c_f_pointer(self%p,self%cmat,self%array_shape)
            ! Work out the shape of the real array
            self%array_shape(1) = 2*(self%array_shape(1))
            ! Set up the real array
            call c_f_pointer(self%p,self%rmat,self%array_shape)
        endif
        ! put back the shape of the complex array
        self%array_shape(1) = fdim(self%ldim(1))
        ! init
        self%rmat = 0.
        self%ft   = .false.
        !$omp critical
        ! make fftw plans
        if( self%wthreads .and. (any(ldim >= 200) .or. ldim(3) >= 100) )then
            rc = fftwf_init_threads()
            call fftwf_plan_with_nthreads(nthr_glob)
        endif
        if(self%ldim(3) > 1)then
            self%plan_fwd = fftwf_plan_dft_r2c_3d(self%ldim(3), self%ldim(2), self%ldim(1), self%rmat, self%cmat, FFTW_ESTIMATE)
            self%plan_bwd = fftwf_plan_dft_c2r_3d(self%ldim(3), self%ldim(2), self%ldim(1), self%cmat, self%rmat, FFTW_ESTIMATE)
        else
            self%plan_fwd = fftwf_plan_dft_r2c_2d(self%ldim(2), self%ldim(1), self%rmat, self%cmat, FFTW_ESTIMATE)
            self%plan_bwd = fftwf_plan_dft_c2r_2d(self%ldim(2), self%ldim(1), self%cmat, self%rmat, FFTW_ESTIMATE)
        endif
        if( self%wthreads .and. (any(ldim >= 200) .or. ldim(3) >= 100) )then
            ! disable threads for subsequent plans
            call fftwf_plan_with_nthreads(1)
        endif
        !$omp end critical
        ! set shift constant (shconst)
        do i=1,3
            if( self%ldim(i) == 1 )then
                self%shconst(i) = 0.
                cycle
            endif
            if( is_even(self%ldim(i)) )then
                self%shconst(i) = PI/real(self%ldim(i)/2.)
            else
                self%shconst(i) = PI/real((self%ldim(i)-1)/2.)
            endif
        end do
        self%existence = .true.
    end subroutine new

    module subroutine set_wthreads( self, wthreads )
        class(image), intent(inout) :: self
        logical,      intent(in)    :: wthreads
        integer(kind=c_int) :: rc
        if( self%wthreads .eqv. wthreads ) return
        self%wthreads = wthreads
        self%wthreads = self%wthreads .and. nthr_glob > 1
        !$omp critical
        call fftwf_destroy_plan(self%plan_fwd)
        call fftwf_destroy_plan(self%plan_bwd)
        ! make fftw plans
        if( self%wthreads )then
            rc = fftwf_init_threads()
            call fftwf_plan_with_nthreads(nthr_glob)
        endif
        if( self%ldim(3) > 1 )then
            self%plan_fwd = fftwf_plan_dft_r2c_3d(self%ldim(3), self%ldim(2), self%ldim(1), self%rmat, self%cmat, FFTW_ESTIMATE)
            self%plan_bwd = fftwf_plan_dft_c2r_3d(self%ldim(3), self%ldim(2), self%ldim(1), self%cmat, self%rmat, FFTW_ESTIMATE)
        else
            self%plan_fwd = fftwf_plan_dft_r2c_2d(self%ldim(2), self%ldim(1), self%rmat, self%cmat, FFTW_ESTIMATE)
            self%plan_bwd = fftwf_plan_dft_c2r_2d(self%ldim(2), self%ldim(1), self%cmat, self%rmat, FFTW_ESTIMATE)
        endif
        if( self%wthreads )then
            ! disable threads for subsequent plans
            call fftwf_plan_with_nthreads(1)
        endif
        !$omp end critical
    end subroutine set_wthreads

    module subroutine construct_thread_safe_tmp_imgs( self, nthr )
        class(image), intent(in) :: self
        integer,      intent(in) :: nthr
        integer :: i, sz, ldim(3)
        logical :: do_allocate
        if( allocated(thread_safe_tmp_imgs) )then
            ldim = thread_safe_tmp_imgs(1)%get_ldim()
            sz   = size(thread_safe_tmp_imgs)
            if( any(self%ldim /= ldim) .or. sz /= nthr )then
                do i=1,size(thread_safe_tmp_imgs)
                    call thread_safe_tmp_imgs(i)%kill
                end do
                deallocate(thread_safe_tmp_imgs)
                do_allocate = .true.
            else
                do_allocate = .false.
            endif
        else
            do_allocate = .true.
        endif
        if( do_allocate )then
            allocate( thread_safe_tmp_imgs(nthr) )
            do i=1,nthr
                call thread_safe_tmp_imgs(i)%new(self%ldim, self%smpd, .false.)
            end do
        endif
    end subroutine construct_thread_safe_tmp_imgs

    module subroutine kill( self )
        class(image), intent(inout) :: self
        if( self%existence )then
            call fftwf_free(self%p)
            self%rmat=>null()
            self%cmat=>null()
            !$omp critical
            call fftwf_destroy_plan(self%plan_fwd)
            call fftwf_destroy_plan(self%plan_bwd)
            !$omp end critical
            self%existence = .false.
        endif
    end subroutine kill

    module subroutine kill_thread_safe_tmp_imgs( self )
        class(image), intent(in) :: self
        integer :: i
        if( allocated(thread_safe_tmp_imgs) )then
            do i=1,size(thread_safe_tmp_imgs)
                call thread_safe_tmp_imgs(i)%kill
            end do
            deallocate(thread_safe_tmp_imgs)
        endif
    end subroutine kill_thread_safe_tmp_imgs

    !========================================
    ! Copy / clone
    !========================================

    module subroutine copy( self, self_in )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self_in
        call self%new(self_in%ldim, self_in%smpd, self%wthreads)
        self%rmat = self_in%rmat
        self%ft   = self_in%ft
    end subroutine copy

    module subroutine copy_fast( self, self_in )
        class(image), intent(inout) :: self
        class(image), intent(in)    :: self_in
        if( self%wthreads )then
            !$omp parallel workshare
            self%rmat = self_in%rmat
            !$omp end parallel workshare
        else
            self%rmat = self_in%rmat
        endif
        self%ft = self_in%ft
    end subroutine copy_fast

    !========================================
    ! Shape / generator helpers
    !========================================

    module subroutine disc_1( self, ldim, smpd, radius, npix )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: ldim(3)
        real,              intent(in)    :: smpd, radius
        integer, optional, intent(inout) :: npix
        call self%new(ldim, smpd)
        call self%cendist
        where(self%rmat <= radius)
            self%rmat = 1.
        else where
            self%rmat = 0.
        end where
        if( present(npix) ) npix = count(self%rmat>0.5)
    end subroutine disc_1

    module subroutine disc_2( self, ldim, smpd, radius, lmsk, npix )
        class(image),         intent(inout) :: self
        integer,              intent(in)    :: ldim(3)
        real,                 intent(in)    :: smpd, radius
        logical, allocatable, intent(out)   :: lmsk(:,:,:)
        integer, optional,    intent(inout) :: npix
        call self%new(ldim, smpd)
        call self%cendist
        if( allocated(lmsk) ) deallocate(lmsk)
        allocate(lmsk(ldim(1),ldim(2),ldim(3)))
        where(self%rmat(:ldim(1),:ldim(2),:ldim(3)) <= radius)
            self%rmat(:ldim(1),:ldim(2),:ldim(3)) = 1.
            lmsk = .true.
        else where
            self%rmat(:ldim(1),:ldim(2),:ldim(3)) = 0.
            lmsk = .false.
        end where
        if( present(npix) ) npix = count(lmsk)
    end subroutine disc_2

    module subroutine ring( self, ldim, smpd, outer_radius, inner_radius, npix )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: ldim(3)
        real,              intent(in)    :: smpd, outer_radius, inner_radius
        integer, optional, intent(inout) :: npix
        call self%new(ldim, smpd)
        call self%cendist
        where(self%rmat <= outer_radius .and. self%rmat >= inner_radius )
            self%rmat = 1.
        else where
            self%rmat = 0.
        end where
        if( present(npix) )npix = count(self%rmat>0.5)
    end subroutine ring

    module subroutine soft_ring( self, ldim, smpd, radius )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: ldim(3)
        real,              intent(in)    :: smpd
        real,              intent(in)    :: radius ! in pixels
        real, parameter :: K     = 2.0
        real, parameter :: THETA = 2.0
        real    :: d,km1,mode,val,scale
        integer :: c(3),i,j,l
        call self%new(ldim, smpd)
        c     = nint(real(self%ldim)/2.)+1
        km1   = K-1.
        mode  = km1 * THETA
        val   = mode**km1 * exp(-mode/THETA)
        scale = 1./val
        do l=1,self%ldim(3)
        do j=1,self%ldim(2)
        do i=1,self%ldim(1)
            d = hyp(i-c(1),j-c(2),l-c(3))
            if( d > radius )then
                val = 0.
            else
                d   = 20. * (radius - d) / radius
                val = max(0.,min(1.0,scale * (d**km1 * exp(-d/THETA))))
            endif
            self%rmat(i,j,l) = val
        enddo
        enddo
        enddo
    end subroutine soft_ring

    module subroutine disc_sideview( self, ldim, smpd, radius )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: ldim(3)
        real,              intent(in)    :: smpd
        real,              intent(in)    :: radius ! in pixels
        integer, parameter :: B = 3
        real    :: t,d
        integer :: c(3),i,j,h,f
        if( ldim(3) > 1 ) THROW_HARD('2D only; membrane')
        call self%new(ldim, smpd)
        c = nint(real(self%ldim)/2.)+1
        ! central line
        j = c(2)
        do i = 1,self%ldim(1)
            h = i-c(1)
            d = real(abs(h))
            if( d >= radius )then
                if( d > radius+real(B) )then
                    self%rmat(i,j,1) = 0.
                else
                    ! horizontal soft cosine edge
                    d = (d-radius)/real(B)
                    self%rmat(i,j,1) = 2.*sqrt(radius-0.25)*max(0.,cos(d*PIO2))
                endif
            else
                self%rmat(i,j,1) = self%rmat(i,j,1) + 2.*sqrt(radius**2-real(h**2))
            endif
        enddo
        ! thickness taken as ~40Angs, soft edge excluded
        t = 40./self%smpd
        f = floor(t/2.)
        do j = c(2)-f-B,c(2)+f+B
            if( j == c(2) ) cycle
            if( (j >= c(2)-f) .and. (j <= c(2)+f) )then
                self%rmat(:,j,1) = self%rmat(:,c(2),1)
            else
                ! vertical soft edge
                d = (real( abs(j-c(2)) - f)) / real(B)
                self%rmat(:,j,1) = self%rmat(:,c(2),1) * max(0.,cos(d*PIO2))
            endif
        enddo
        self%rmat = self%rmat / maxval(self%rmat(1:self%ldim(1),c(2),1))
    end subroutine disc_sideview

    module subroutine cylinder( self, ldim, smpd, height, radius )
        class(image),      intent(inout) :: self
        integer,           intent(in)    :: ldim(3)
        real,              intent(in)    :: smpd, height, radius
        real    :: rsq,radiussq
        integer :: icenter(3),i,j,k
        call self%new(ldim, smpd)
        icenter  = nint(real(self%ldim/2.))+1
        radiussq = radius**2
        do k = 1,self%ldim(3)
            if( real(abs(k-icenter(3))) > height/2. ) cycle
            do j = 1,self%ldim(2)
            do i = 1,self%ldim(1)
                rsq = (i-icenter(1))**2 + (j-icenter(2))**2
                if( rsq < radiussq ) self%rmat(i,j,k) = 1
            enddo
            enddo
        enddo
    end subroutine cylinder

end submodule simple_image_core
