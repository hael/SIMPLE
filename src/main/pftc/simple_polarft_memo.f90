submodule (simple_polarft_calc) simple_polarft_memo
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
#include "simple_local_flags.inc"
implicit none

contains

    module subroutine setup_npix_per_shell(self)
        class(polarft_calc), intent(inout) :: self
        integer :: h,k,sh
        if( allocated(self%npix_per_shell) ) deallocate(self%npix_per_shell)
        allocate(self%npix_per_shell(self%kfromto(1):self%kfromto(2)),source=0.0)
        do h = 0,self%kfromto(2)
            do k = -self%kfromto(2),self%kfromto(2)
                if( (h==0) .and. (k>0) ) cycle
                sh = nint(sqrt(real(h**2+k**2)))
                if( sh < self%kfromto(1) ) cycle
                if( sh > self%kfromto(2) ) cycle
                self%npix_per_shell(sh) = self%npix_per_shell(sh) + 1.0
            end do
        end do
    end subroutine setup_npix_per_shell

    module subroutine memoize_sqsum_ptcl(self, iptcl)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iptcl
        real(dp) :: sumsqk
        integer  :: i, ik
        logical  :: l_sigma
        i       = self%pinds(iptcl)
        l_sigma = associated(self%sigma2_noise)
        self%sqsums_ptcls(i)  = 0.d0
        self%ksqsums_ptcls(i) = 0.d0
        if( l_sigma ) self%wsqsums_ptcls(i) = 0.d0
        do ik = self%kfromto(1),self%kfromto(2)
            sumsqk                = sum(real(self%pfts_ptcls(:,ik,i)*conjg(self%pfts_ptcls(:,ik,i)),dp))
            self%sqsums_ptcls(i)  = self%sqsums_ptcls(i) + sumsqk
            sumsqk                = real(ik,dp) * sumsqk
            self%ksqsums_ptcls(i) = self%ksqsums_ptcls(i) + sumsqk
            if( l_sigma ) self%wsqsums_ptcls(i) = self%wsqsums_ptcls(i) + sumsqk / real(self%sigma2_noise(ik,iptcl),dp)
        enddo
    end subroutine memoize_sqsum_ptcl

    module subroutine memoize_ptcls(self)
        class(polarft_calc), intent(inout) :: self
        integer :: ithr,i,k
        logical :: l_memoize_absptcl_ctf
        l_memoize_absptcl_ctf = trim(params_glob%sh_inv).eq.'yes'
        !$omp parallel do collapse(2) private(i,k,ithr) default(shared) proc_bind(close) schedule(static)
        do i = 1,self%nptcls
            do k = self%kfromto(1),self%kfromto(2)
                ithr = omp_get_thread_num() + 1
                ! FT(X.CTF)
                if( self%with_ctf )then
                    self%cvec2(ithr)%c(1:self%pftsz) = self%pfts_ptcls(:,k,i) * self%ctfmats(:,k,i)
                else
                    self%cvec2(ithr)%c(1:self%pftsz) = self%pfts_ptcls(:,k,i)
                endif
                self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(self%cvec2(ithr)%c(1:self%pftsz))
                call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
                self%ft_ptcl_ctf(:,k,i) = self%cvec2(ithr)%c(1:self%pftsz+1)
                ! FT(CTF2)
                if( self%with_ctf )then
                    self%rvec1(ithr)%r(1:self%pftsz)            = self%ctfmats(:,k,i)*self%ctfmats(:,k,i)
                    self%rvec1(ithr)%r(self%pftsz+1:self%nrots) = self%rvec1(ithr)%r(1:self%pftsz)
                else
                    self%rvec1(ithr)%r = 1.0
                endif
                call fftwf_execute_dft_r2c(self%plan_mem_r2c, self%rvec1(ithr)%r, self%cvec1(ithr)%c)
                self%ft_ctf2(:,k,i) = self%cvec1(ithr)%c(1:self%pftsz+1)
                if( l_memoize_absptcl_ctf )then
                    if( self%with_ctf )then
                        self%cvec2(ithr)%c(1:self%pftsz) = abs(self%pfts_ptcls(:,k,i)) * self%ctfmats(:,k,i)
                    else
                        self%cvec2(ithr)%c(1:self%pftsz) = abs(self%pfts_ptcls(:,k,i))
                    endif
                    self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(self%cvec2(ithr)%c(1:self%pftsz))
                    call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
                    self%ft_absptcl_ctf(:,k,i) = self%cvec2(ithr)%c(1:self%pftsz+1)
                endif
            enddo
        enddo
        !$omp end parallel do
    end subroutine memoize_ptcls

    module subroutine memoize_refs(self)
        class(polarft_calc), intent(inout) :: self
        integer :: k, ithr, iref
        ! allocations
        call self%allocate_refs_memoization
        ! memoization
        !$omp parallel do collapse(2) private(iref,k,ithr) default(shared) proc_bind(close) schedule(static)
        do iref = 1,self%nrefs
            do k = self%kfromto(1),self%kfromto(2)
                ithr = omp_get_thread_num() + 1
                ! FT(REFeven)*
                self%cvec2(ithr)%c(           1:self%pftsz) = self%pfts_refs_even(:,k,iref)
                self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(self%cvec2(ithr)%c(1:self%pftsz))
                call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
                self%ft_ref_even(:,k,iref) = conjg(self%cvec2(ithr)%c(1:self%pftsz+1))
                ! FT(REFodd)*
                self%cvec2(ithr)%c(           1:self%pftsz) = self%pfts_refs_odd(:,k,iref)
                self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(self%cvec2(ithr)%c(1:self%pftsz))
                call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
                self%ft_ref_odd(:,k,iref) = conjg(self%cvec2(ithr)%c(1:self%pftsz+1))
                ! FT(REF2even)*
                self%rvec1(ithr)%r(           1:self%pftsz) = real(self%pfts_refs_even(:,k,iref)*conjg(self%pfts_refs_even(:,k,iref)))
                self%rvec1(ithr)%r(self%pftsz+1:self%nrots) = self%rvec1(ithr)%r(1:self%pftsz)
                call fftwf_execute_dft_r2c(self%plan_mem_r2c, self%rvec1(ithr)%r, self%cvec1(ithr)%c)
                self%ft_ref2_even(:,k,iref) = conjg(self%cvec1(ithr)%c(1:self%pftsz+1))
                ! FT(REF2odd)*
                self%rvec1(ithr)%r(           1:self%pftsz) = real(self%pfts_refs_odd(:,k,iref)*conjg(self%pfts_refs_odd(:,k,iref)))
                self%rvec1(ithr)%r(self%pftsz+1:self%nrots) = self%rvec1(ithr)%r(1:self%pftsz)
                call fftwf_execute_dft_r2c(self%plan_mem_r2c, self%rvec1(ithr)%r, self%cvec1(ithr)%c)
                self%ft_ref2_odd(:,k,iref) = conjg(self%cvec1(ithr)%c(1:self%pftsz+1))
            enddo
        enddo
        !$omp end parallel do
    end subroutine memoize_refs

    module subroutine allocate_ptcls_memoization(self)
        class(polarft_calc), intent(inout) :: self
        allocate(self%ft_ptcl_ctf(self%pftsz+1,self%kfromto(1):self%kfromto(2),self%nptcls),&
                &self%ft_ctf2(    self%pftsz+1,self%kfromto(1):self%kfromto(2),self%nptcls))
        if( trim(params_glob%sh_inv).eq.'yes' )then
            allocate(self%ft_absptcl_ctf(self%pftsz+1,self%kfromto(1):self%kfromto(2),self%nptcls))
        endif
    end subroutine allocate_ptcls_memoization

    module subroutine allocate_refs_memoization(self)
        class(polarft_calc), intent(inout) :: self
        character(kind=c_char, len=:), allocatable :: fft_wisdoms_fname ! FFTW wisdoms (per part or suffer I/O lag)
        integer(kind=c_int) :: wsdm_ret
        integer             :: ithr
        if( allocated(self%ft_ref_even) ) call self%kill_memoized_refs
        allocate(self%ft_ref_even( self%pftsz+1,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%ft_ref_odd(  self%pftsz+1,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%ft_ref2_even(self%pftsz+1,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%ft_ref2_odd( self%pftsz+1,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%rvec1(nthr_glob), self%cvec1(nthr_glob),self%cvec2(nthr_glob),&
                &self%drvec(nthr_glob))
        ! convenience objects
        do ithr = 1,nthr_glob
            self%cvec1(ithr)%p = fftwf_alloc_complex(int(self%pftsz+1, c_size_t))
            self%cvec2(ithr)%p = fftwf_alloc_complex(int(self%nrots, c_size_t))
            call c_f_pointer(self%cvec1(ithr)%p, self%cvec1(ithr)%c, [self%pftsz+1])
            call c_f_pointer(self%cvec1(ithr)%p, self%rvec1(ithr)%r, [self%nrots+2])
            call c_f_pointer(self%cvec2(ithr)%p, self%cvec2(ithr)%c, [self%nrots])
            self%drvec(ithr)%p = fftw_alloc_real(int(self%nrots, c_size_t))
            call c_f_pointer(self%drvec(ithr)%p, self%drvec(ithr)%r, [self%nrots])
        enddo
        ! plans & FFTW3 wisdoms
        if( params_glob%l_distr_exec )then
            allocate(fft_wisdoms_fname, source='fft_wisdoms_part'//int2str_pad(params_glob%part,params_glob%numlen)//'.dat'//c_null_char)
        else
            allocate(fft_wisdoms_fname, source='fft_wisdoms.dat'//c_null_char)
        endif
        wsdm_ret = fftw_import_wisdom_from_filename(fft_wisdoms_fname)
        self%plan_fwd1    = fftwf_plan_dft_1d(    self%nrots, self%cvec2(1)%c, self%cvec2(1)%c, FFTW_FORWARD, ior(FFTW_PATIENT, FFTW_USE_WISDOM))
        self%plan_bwd1    = fftwf_plan_dft_c2r_1d(self%nrots, self%cvec1(1)%c, self%rvec1(1)%r,               ior(FFTW_PATIENT, FFTW_USE_WISDOM))
        self%plan_mem_r2c = fftwf_plan_dft_r2c_1d(self%nrots, self%rvec1(1)%r, self%cvec1(1)%c,               ior(FFTW_PATIENT, FFTW_USE_WISDOM))
        wsdm_ret = fftw_export_wisdom_to_filename(fft_wisdoms_fname)
        deallocate(fft_wisdoms_fname)
        if (wsdm_ret == 0) then
            write (*, *) 'Error: could not write FFTW3 wisdom file! Check permissions.'
        end if
    end subroutine allocate_refs_memoization

    module subroutine kill_memoized_ptcls(self)
        class(polarft_calc), intent(inout) :: self
        if( allocated(self%ft_ptcl_ctf) )    deallocate(self%ft_ptcl_ctf,self%ft_ctf2)
        if( allocated(self%ft_absptcl_ctf) ) deallocate(self%ft_absptcl_ctf)
    end subroutine kill_memoized_ptcls

    module subroutine kill_memoized_refs(self)
        class(polarft_calc), intent(inout) :: self
        integer :: i
        if( allocated(self%ft_ref_even) )then
            do i = 1,size(self%cvec1,dim=1)
                call fftwf_free(self%cvec1(i)%p)
                call fftwf_free(self%cvec2(i)%p)
                call fftw_free( self%drvec(i)%p)
            enddo
            deallocate(self%ft_ref_even,self%ft_ref_odd,self%ft_ref2_even,self%ft_ref2_odd,&
            &self%rvec1,self%cvec1,self%cvec2,self%drvec)
            call fftwf_destroy_plan(self%plan_fwd1)
            call fftwf_destroy_plan(self%plan_bwd1)
            call fftwf_destroy_plan(self%plan_mem_r2c)
        endif
    end subroutine kill_memoized_refs

end submodule simple_polarft_memo
