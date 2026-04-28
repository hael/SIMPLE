!@descr: polarft class submodule for memoization for performance
submodule (simple_polarft_calc) simple_polarft_memo
#include "simple_local_flags.inc"
implicit none

contains

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
        integer :: ithr, i, k, kk, k0
        if( OMP_IN_PARALLEL() )then
            THROW_HARD('No memoization inside OpenMP regions')
        endif
        k0 = self%kfromto(1)
        !$omp parallel do private(i,k,kk,ithr) default(shared) proc_bind(close) schedule(static)
        do i = 1, self%nptcls
            ithr = omp_get_thread_num() + 1
            ! ========================================================================
            ! Batched computation of FT(X.CTF) for all k shells
            ! ========================================================================
            if( self%with_ctf )then
                do k = self%kfromto(1), self%kfromto(2)
                    kk = k - k0 + 1
                    self%cmat2_many(ithr)%c(1:self%pftsz, kk) = self%pfts_ptcls(:,k,i) * self%ctfmats(:,k,i)
                    self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots, kk) = conjg(self%cmat2_many(ithr)%c(1:self%pftsz, kk))
                end do
            else
                do k = self%kfromto(1), self%kfromto(2)
                    kk = k - k0 + 1
                    self%cmat2_many(ithr)%c(1:self%pftsz, kk) = self%pfts_ptcls(:,k,i)
                    self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots, kk) = conjg(self%cmat2_many(ithr)%c(1:self%pftsz, kk))
                end do
            endif
            ! Execute batched C2C FFT
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
            ! Extract and store FT(X.CTF) results
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%ft_ptcl_ctf(:,k,i) = self%cmat2_many(ithr)%c(1:self%pftsz+1, kk)
            end do
            ! ========================================================================
            ! Batched computation of FT(CTF2) for all k shells
            ! ========================================================================
            if( self%with_ctf )then
                do k = self%kfromto(1), self%kfromto(2)
                    kk = k - k0 + 1
                    self%crmat1_many(ithr)%r(1:self%pftsz, kk)            = self%ctfmats(:,k,i) * self%ctfmats(:,k,i)
                    self%crmat1_many(ithr)%r(self%pftsz+1:self%nrots, kk) = self%crmat1_many(ithr)%r(1:self%pftsz, kk)
                end do
            else
                self%crmat1_many(ithr)%r = 1.0
            endif
            ! Execute batched R2C FFT (requires plan_mem_r2c_many)
            call fftwf_execute_dft_r2c(self%plan_mem_r2c_many, self%crmat1_many(ithr)%r, self%crmat1_many(ithr)%c)
            ! Extract and store FT(CTF2) results
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%ft_ctf2(:,k,i) = self%crmat1_many(ithr)%c(:,kk)
            end do
        enddo
        !$omp end parallel do
    end subroutine memoize_ptcls

    module subroutine memoize_refs(self, eulspace)
        class(polarft_calc),           intent(inout) :: self
        class(oris),         optional, intent(in)    :: eulspace
        real     :: psi
        integer :: k, kk, k0, ithr, iref, mref
        logical :: l_even_inds(self%pftsz+1,1:self%nk)
        logical :: l_mirr_proj, l_rot
        if( OMP_IN_PARALLEL() )then
            THROW_HARD('No memoization inside OpenMP regions')
        endif
        l_mirr_proj = trim(self%p_ptr%mirr_proj)=='yes'
        if( l_mirr_proj ) then
            ! sanity checks
            if( .not.present(eulspace) )then
                THROW_HARD('Reference orientations object is required for projection mirroring')
            endif
            if( .not.eulspace%isthere('mirr') )then
                THROW_HARD('Mirror reprojections indices absent from Reference orientations object')
            endif
            ! filling even/odd index map
            l_even_inds(2:self%pftsz+1:2,:) = .true.
            l_even_inds(1:self%pftsz+1:2,:) = .false.
            ! mirror reprojection indices
            !$omp parallel do private(iref,ithr,mref,psi,l_rot) default(shared) proc_bind(close) schedule(static)
            do iref = 1, self%nrefs/2
                ithr = omp_get_thread_num() + 1
                ! index of mirror reprojection
                mref  = eulspace%get_int(iref,'mirr')
                psi   = eulspace%get(mref,'psi')
                l_rot = (psi>0.1) .and. (psi<359.9)
                ! ========================================================================
                ! Batched computation of FT(REFeven)* for all k shells
                ! ========================================================================
                self%cmat2_many(ithr)%c(           1:self%pftsz, 1:self%nk) = self%pfts_refs_even(:,self%kfromto(1):self%kfromto(2),iref)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots, 1:self%nk) = conjg(self%cmat2_many(ithr)%c(1:self%pftsz, 1:self%nk))
                ! Execute batched C2C FFT for even references
                call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
                ! Extract and store FT(REFeven)* results
                self%ft_ref_even(:, self%kfromto(1):self%kfromto(2), iref) = self%cmat2_many(ithr)%c(1:self%pftsz+1, 1:self%nk)
                ! Store FT(Mirror(REFeven))*
                if( l_rot )then
                    self%ft_ref_even(:,:, mref) = conjg(merge(-self%cmat2_many(ithr)%c(1:self%pftsz+1, 1:self%nk),&
                                                              &self%cmat2_many(ithr)%c(1:self%pftsz+1, 1:self%nk), l_even_inds(:,:)))
                else
                    self%ft_ref_even(:, :, mref) = conjg(self%cmat2_many(ithr)%c(1:self%pftsz+1, 1:self%nk))
                endif
                ! ========================================================================
                ! Batched computation of FT(REF2even)* for all k shells
                ! ========================================================================
                self%crmat1_many(ithr)%r(1:self%pftsz, 1:self%nk) = real(self%pfts_refs_even(:,self%kfromto(1):self%kfromto(2),iref)*&
                                                                        &conjg(self%pfts_refs_even(:,self%kfromto(1):self%kfromto(2),iref)))
                self%crmat1_many(ithr)%r(self%pftsz+1:self%nrots, 1:self%nk) = self%crmat1_many(ithr)%r(1:self%pftsz, 1:self%nk)
                ! Execute batched R2C FFT for even reference squares
                call fftwf_execute_dft_r2c(self%plan_mem_r2c_many, self%crmat1_many(ithr)%r, self%crmat1_many(ithr)%c)
                ! Extract and store FT(REF2even)*
                self%ft_ref2_even(:, self%kfromto(1):self%kfromto(2), iref) = conjg(self%crmat1_many(ithr)%c(:, 1:self%nk))
                ! Extract and store FT(Mirr(REF2even))* = FT(REF2even)
                self%ft_ref2_even(:, self%kfromto(1):self%kfromto(2), mref) =       self%crmat1_many(ithr)%c(:, 1:self%nk)
                ! ========================================================================
                ! Batched computation of FT(REFodd)* for all k shells
                ! ========================================================================
                self%cmat2_many(ithr)%c(           1:self%pftsz, 1:self%nk) = self%pfts_refs_odd(:,self%kfromto(1):self%kfromto(2),iref)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots, 1:self%nk) = conjg(self%cmat2_many(ithr)%c(1:self%pftsz, 1:self%nk))
                ! Execute batched C2C FFT for odd references
                call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
                ! Extract and store FT(REFodd)* results
                self%ft_ref_odd(:, self%kfromto(1):self%kfromto(2), iref) = self%cmat2_many(ithr)%c(1:self%pftsz+1, 1:self%nk)
                ! Store FT(Mirror(REFodd))*
                if( l_rot )then
                    self%ft_ref_odd(:, :, mref) = conjg(merge(-self%cmat2_many(ithr)%c(1:self%pftsz+1, 1:self%nk),&
                                                              &self%cmat2_many(ithr)%c(1:self%pftsz+1, 1:self%nk), l_even_inds(:,:)))
                else
                    self%ft_ref_odd(:, :, mref) = conjg(self%cmat2_many(ithr)%c(1:self%pftsz+1, 1:self%nk))
                endif
                ! ========================================================================
                ! Batched computation of FT(REF2odd)* for all k shells
                ! ========================================================================
                self%crmat1_many(ithr)%r(1:self%pftsz, 1:self%nk) = real(self%pfts_refs_odd(:,self%kfromto(1):self%kfromto(2),iref)*&
                                                                        &conjg(self%pfts_refs_odd(:,self%kfromto(1):self%kfromto(2),iref)))
                self%crmat1_many(ithr)%r(self%pftsz+1:self%nrots, 1:self%nk) = self%crmat1_many(ithr)%r(1:self%pftsz, 1:self%nk)
                ! Execute batched R2C FFT for odd reference squares
                call fftwf_execute_dft_r2c(self%plan_mem_r2c_many, self%crmat1_many(ithr)%r, self%crmat1_many(ithr)%c)
                ! Extract and store FT(REF2odd)*
                self%ft_ref2_odd(:, self%kfromto(1):self%kfromto(2), iref) = conjg(self%crmat1_many(ithr)%c(:, 1:self%nk))
                ! Extract and store FT(Mirr(REF2odd))* = FT(REF2odd)
                self%ft_ref2_odd(:, self%kfromto(1):self%kfromto(2), mref) =       self%crmat1_many(ithr)%c(:, 1:self%nk)
            enddo
            !$omp end parallel do
        else
            k0 = self%kfromto(1)
            !$omp parallel do private(iref,k,kk,ithr) default(shared) proc_bind(close) schedule(static)
            do iref = 1, self%nrefs
                ithr = omp_get_thread_num() + 1
                ! ========================================================================
                ! Batched computation of FT(REFeven)* for all k shells
                ! ========================================================================
                do k = self%kfromto(1), self%kfromto(2)
                    kk = k - k0 + 1
                    self%cmat2_many(ithr)%c(1:self%pftsz, kk) = self%pfts_refs_even(:,k,iref)
                    self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots, kk) = conjg(self%cmat2_many(ithr)%c(1:self%pftsz, kk))
                end do
                ! Execute batched C2C FFT for even references
                call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
                ! Extract and store FT(REFeven)* results
                do k = self%kfromto(1), self%kfromto(2)
                    kk = k - k0 + 1
                    self%ft_ref_even(:,k,iref) = self%cmat2_many(ithr)%c(1:self%pftsz+1, kk)
                end do
                ! ========================================================================
                ! Batched computation of FT(REFodd)* for all k shells
                ! ========================================================================
                do k = self%kfromto(1), self%kfromto(2)
                    kk = k - k0 + 1
                    self%cmat2_many(ithr)%c(1:self%pftsz, kk) = self%pfts_refs_odd(:,k,iref)
                    self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots, kk) = conjg(self%cmat2_many(ithr)%c(1:self%pftsz, kk))
                end do
                ! Execute batched C2C FFT for odd references
                call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
                ! Extract and store FT(REFodd)* results
                do k = self%kfromto(1), self%kfromto(2)
                    kk = k - k0 + 1
                    self%ft_ref_odd(:,k,iref) = self%cmat2_many(ithr)%c(1:self%pftsz+1, kk)
                end do
                ! ========================================================================
                ! Batched computation of FT(REF2even)* for all k shells
                ! ========================================================================
                do k = self%kfromto(1), self%kfromto(2)
                    kk = k - k0 + 1
                    self%crmat1_many(ithr)%r(1:self%pftsz, kk) = real(self%pfts_refs_even(:,k,iref) * conjg(self%pfts_refs_even(:,k,iref)))
                    self%crmat1_many(ithr)%r(self%pftsz+1:self%nrots, kk) = self%crmat1_many(ithr)%r(1:self%pftsz, kk)
                end do
                ! Execute batched R2C FFT for even reference squares
                call fftwf_execute_dft_r2c(self%plan_mem_r2c_many, self%crmat1_many(ithr)%r, self%crmat1_many(ithr)%c)
                ! Extract and store FT(REF2even)* results
                do k = self%kfromto(1), self%kfromto(2)
                    kk = k - k0 + 1
                    self%ft_ref2_even(:,k,iref) = conjg(self%crmat1_many(ithr)%c(:,kk))
                end do
                ! ========================================================================
                ! Batched computation of FT(REF2odd)* for all k shells
                ! ========================================================================
                do k = self%kfromto(1), self%kfromto(2)
                    kk = k - k0 + 1
                    self%crmat1_many(ithr)%r(1:self%pftsz, kk) = real(self%pfts_refs_odd(:,k,iref) * conjg(self%pfts_refs_odd(:,k,iref)))
                    self%crmat1_many(ithr)%r(self%pftsz+1:self%nrots, kk) = self%crmat1_many(ithr)%r(1:self%pftsz, kk)
                end do
                ! Execute batched R2C FFT for odd reference squares
                call fftwf_execute_dft_r2c(self%plan_mem_r2c_many, self%crmat1_many(ithr)%r, self%crmat1_many(ithr)%c)
                ! Extract and store FT(REF2odd)* results
                do k = self%kfromto(1), self%kfromto(2)
                    kk = k - k0 + 1
                    self%ft_ref2_odd(:,k,iref) = conjg(self%crmat1_many(ithr)%c(:,kk))
                end do
            enddo
            !$omp end parallel do
        endif
    end subroutine memoize_refs

    module subroutine alloc_memo_ptcls(self)
        class(polarft_calc), intent(inout) :: self
        call self%kill_memo_ptcls
        allocate(self%ft_ptcl_ctf(self%pftsz+1,self%kfromto(1):self%kfromto(2),self%nptcls),&
                &self%ft_ctf2(    self%pftsz+1,self%kfromto(1):self%kfromto(2),self%nptcls))
    end subroutine alloc_memo_ptcls

    module subroutine alloc_memo_refs(self)
        class(polarft_calc), intent(inout) :: self
        call self%kill_memo_refs
        allocate(self%ft_ref_even( self%pftsz+1,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%ft_ref_odd(  self%pftsz+1,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%ft_ref2_even(self%pftsz+1,self%kfromto(1):self%kfromto(2),self%nrefs),&
                &self%ft_ref2_odd( self%pftsz+1,self%kfromto(1):self%kfromto(2),self%nrefs))
    end subroutine alloc_memo_refs

    module subroutine allocate_memo_workspace(self)
        class(polarft_calc), intent(inout) :: self
        character(kind=c_char, len=:), allocatable :: fft_wisdoms_fname ! FFTW wisdoms (per part or suffer I/O lag)
        integer(kind=c_int) :: wsdm_ret, rank, howmany, n(1),  inembed(1), onembed(1), istride, ostride, idist, odist
        integer             :: ithr
        call self%kill_memo_workspace
        allocate(self%drvec(nthr_glob), self%cmat2_many(nthr_glob), self%crmat1_many(nthr_glob))
        ! convenience objects
        do ithr = 1,nthr_glob
            self%drvec(ithr)%p = fftw_alloc_real(int(self%nrots, c_size_t))
            call c_f_pointer(self%drvec(ithr)%p, self%drvec(ithr)%r, [self%nrots])
            self%cmat2_many(ithr)%p = fftwf_alloc_complex(int(self%nrots * self%nk, c_size_t))
            call c_f_pointer(self%cmat2_many(ithr)%p, self%cmat2_many(ithr)%c, [self%nrots, self%nk])
            ! Allocate complex storage for (pftsz+1) * nk transforms
            self%crmat1_many(ithr)%p = fftwf_alloc_complex(int((self%pftsz+1) * self%nk, c_size_t))
            ! Complex view: (pftsz+1, nk)
            call c_f_pointer(self%crmat1_many(ithr)%p, self%crmat1_many(ithr)%c, [self%pftsz+1, self%nk])
            ! Real in-place view on same memory: (nrots+2, nk)
            call c_f_pointer(self%crmat1_many(ithr)%p, self%crmat1_many(ithr)%r, [self%nrots+2, self%nk])
        enddo
        ! plans & FFTW3 wisdoms
        if( self%p_ptr%l_distr_worker )then
            allocate(fft_wisdoms_fname, source='fft_wisdoms_part'//int2str_pad(self%p_ptr%part,self%p_ptr%numlen)//'.dat'//c_null_char)
        else
            allocate(fft_wisdoms_fname, source='fft_wisdoms.dat'//c_null_char)
        endif
        wsdm_ret = fftw_import_wisdom_from_filename(fft_wisdoms_fname)
        ! plan the many
        rank       = 1_c_int
        n(1)       = int(self%nrots, c_int)
        howmany    = int(self%nk, c_int)
        ! Layout: cmat2_many(:,:) is (nrots, nk) column-major, each column contiguous
        istride    = 1_c_int
        ostride    = 1_c_int
        idist      = int(self%nrots, c_int)
        odist      = int(self%nrots, c_int)
        inembed(1) = n(1)
        onembed(1) = n(1)
        self%plan_fwd1_many = fftwf_plan_many_dft( rank, n, howmany, &
        &self%cmat2_many(1)%c, inembed, istride, idist, &
        &self%cmat2_many(1)%c, onembed, ostride, odist, &
        &FFTW_FORWARD, ior(FFTW_MEASURE, FFTW_USE_WISDOM))
        ! Input complex length is (n/2+1) = pftsz+1
        inembed(1) = int(self%pftsz+1, c_int)
        idist      = int(self%pftsz+1, c_int)
        ! Output real is in-place, but FFTW still wants the logical "n" length.
        ! We use the same in-place padding convention as your existing code (nrots+2),
        ! but for plan_many we pass onembed=n (FFTW uses stride/dist to walk outputs).
        onembed(1) = int(self%nrots, c_int)
        odist      = int(self%nrots+2, c_int)
        self%plan_bwd1_many = fftwf_plan_many_dft_c2r( rank, n, howmany, &
        &self%crmat1_many(1)%c, inembed, istride, idist, &
        &self%crmat1_many(1)%r, onembed, ostride, odist, &
        &ior(FFTW_MEASURE, FFTW_USE_WISDOM) )
        ! Plan the many R2C transforms for CTF2 memoization
        rank       = 1_c_int
        n(1)       = int(self%nrots, c_int)
        howmany    = int(self%nk, c_int)
        ! Input real layout: (nrots+2, nk) with padding
        inembed(1) = int(self%nrots, c_int)
        istride    = 1_c_int
        idist      = int(self%nrots+2, c_int)
        ! Output complex layout: (pftsz+1, nk)
        onembed(1) = int(self%pftsz+1, c_int)
        ostride    = 1_c_int
        odist      = int(self%pftsz+1, c_int)
        self%plan_mem_r2c_many = fftwf_plan_many_dft_r2c( rank, n, howmany, &
        &self%crmat1_many(1)%r, inembed, istride, idist, &
        &self%crmat1_many(1)%c, onembed, ostride, odist, &
        &ior(FFTW_MEASURE, FFTW_USE_WISDOM) )
        ! Retain FFTW wisdom
        wsdm_ret = fftw_export_wisdom_to_filename(fft_wisdoms_fname)
        deallocate(fft_wisdoms_fname)
        if (wsdm_ret == 0) then
            write (*, *) 'Error: could not write FFTW3 wisdom file! Check permissions.'
        end if
    end subroutine allocate_memo_workspace

    module subroutine kill_memo_ptcls(self)
        class(polarft_calc), intent(inout) :: self
        if( allocated(self%ft_ptcl_ctf) )    deallocate(self%ft_ptcl_ctf,self%ft_ctf2)
        if( allocated(self%ft_absptcl_ctf) ) deallocate(self%ft_absptcl_ctf)
    end subroutine kill_memo_ptcls

    module subroutine kill_memo_refs(self)
        class(polarft_calc), intent(inout) :: self
        if( allocated(self%ft_ref_even) )then
            deallocate(self%ft_ref_even,self%ft_ref_odd,self%ft_ref2_even,self%ft_ref2_odd)
        endif
    end subroutine kill_memo_refs

    module subroutine kill_memo_workspace(self)
        class(polarft_calc), intent(inout) :: self
        integer :: ithr
        if( allocated(self%cmat2_many) )then
            do ithr = 1,size(self%drvec)
                call fftw_free( self%drvec(ithr)%p)
                call fftwf_free(self%cmat2_many(ithr)%p)
                call fftwf_free(self%crmat1_many(ithr)%p)
            enddo
            deallocate(self%drvec, self%cmat2_many, self%crmat1_many)
            call fftwf_destroy_plan(self%plan_mem_r2c_many)
            call fftwf_destroy_plan(self%plan_fwd1_many)
            call fftwf_destroy_plan(self%plan_bwd1_many)
        endif
    end subroutine kill_memo_workspace

end submodule simple_polarft_memo
