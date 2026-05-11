!@descr: polarft class submodule for objective function evaluations
submodule (simple_polarft_calc) simple_polarft_corr
use simple_eul_prob_tab_utils, only: sample_bounded_dist, sample_power_dist
#include "simple_local_flags.inc"
implicit none

real, parameter :: SHERRSQ = 0.00001

contains

    ! Benchmark for correlation calculation
    ! Is not FFT-accelerated, does not rely on memoization, for reference only
    module real function calc_corr_rot_shift(self, iref, iptcl, shvec, irot, kweight)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl, irot
        real(sp),            intent(in)    :: shvec(2)
        logical, optional,   intent(in)    :: kweight
        complex(dp), pointer :: pft_ref(:,:), shmat(:,:), pft_rot_ref(:,:)
        real(dp)    :: sqsumref, sqsumptcl, num
        integer     :: i, k, ithr
        logical     :: kw
        kw = .true.
        if( present(kweight) ) kw = kweight
        calc_corr_rot_shift = 0.
        i           = self%pinds(iptcl)
        ithr        = omp_get_thread_num() + 1
        pft_ref     => self%heap_vars(ithr)%pft_ref_8
        pft_rot_ref => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat       => self%heap_vars(ithr)%shmat_8
        pft_ref     = merge(self%pfts_refs_even(:,self%kfromto(1):self%kfromto(2),iref), &
                    self%pfts_refs_odd(:,self%kfromto(1):self%kfromto(2),iref), self%iseven(i))
        call self%gen_shmat4aln_8(ithr, real(shvec,dp),shmat)
        pft_ref = pft_ref * shmat(:,:self%kfromto(2))
        call self%rotate_ref_8(pft_ref, irot, pft_rot_ref)
        pft_rot_ref = pft_rot_ref * self%ctfmats(:,:self%kfromto(2),i)
        select case(self%p_ptr%cc_objfun)
        case(OBJFUN_CC)
            sqsumref  = 0.d0
            sqsumptcl = 0.d0
            num       = 0.d0
            do k = self%kfromto(1),self%kfromto(2)
                if( kw )then
                    sqsumptcl = sqsumptcl + real(k,dp) * real(sum(self%pfts_ptcls(:,k,i) * conjg(self%pfts_ptcls(:,k,i))),dp)
                    sqsumref  = sqsumref  + real(k,dp) * real(sum(pft_rot_ref(:,k)       * conjg(pft_rot_ref(:,k))),dp)
                    num       = num       + real(k,dp) * real(sum(pft_rot_ref(:,k)       * conjg(self%pfts_ptcls(:,k,i))),dp)
                else
                    sqsumptcl = sqsumptcl + real(sum(self%pfts_ptcls(:,k,i)              * conjg(self%pfts_ptcls(:,k,i))),dp)
                    sqsumref  = sqsumref  + real(sum(pft_rot_ref(:,k)                    * conjg(pft_rot_ref(:,k))),dp)
                    num       = num       + real(sum(pft_rot_ref(:,k)                    * conjg(self%pfts_ptcls(:,k,i))),dp)
                endif
            enddo
            calc_corr_rot_shift = real(num/sqrt(sqsumref*sqsumptcl))
        case(OBJFUN_EUCLID)
            pft_rot_ref = pft_rot_ref - self%pfts_ptcls(:,:self%kfromto(2),i)
            sqsumptcl = 0.d0
            num       = 0.d0
            do k = self%kfromto(1),self%kfromto(2)
                if( kw )then
                    sqsumptcl = sqsumptcl + (real(k,dp) / self%sigma2_noise(k,iptcl)) * sum(real(csq_fast(self%pfts_ptcls(:,k,i)),dp))
                    num       = num       + (real(k,dp) / self%sigma2_noise(k,iptcl)) * sum(csq_fast(pft_rot_ref(:,k)))
                else
                    sqsumptcl = sqsumptcl + (1.d0 / self%sigma2_noise(k,iptcl))       * sum(real(csq_fast(self%pfts_ptcls(:,k,i)),dp))
                    num       = num       + (1.d0 / self%sigma2_noise(k,iptcl))       * sum(csq_fast(pft_rot_ref(:,k)))
                endif
            end do
            calc_corr_rot_shift = real(exp( -num / sqsumptcl ))
        end select
    end function calc_corr_rot_shift

    module subroutine calc_frc( self, iref, iptcl, irot, shvec, frc )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl, irot
        real(sp),            intent(in)    :: shvec(2)
        real(sp),            intent(out)   :: frc(self%kfromto(1):self%kfromto(2))
        complex(dp), pointer :: pft_ref(:,:), shmat(:,:), pft_rot_ref(:,:)
        real(dp) :: sumsqref, sumsqptcl, denom, num
        integer  :: k, ithr, i
        i    =  self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref     => self%heap_vars(ithr)%pft_ref_8
        pft_rot_ref => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat       => self%heap_vars(ithr)%shmat_8
        pft_ref = merge(self%pfts_refs_even(:,self%kfromto(1):self%kfromto(2),iref), &
                self%pfts_refs_odd(:,self%kfromto(1):self%kfromto(2),iref), self%iseven(i))
        call self%gen_shmat4aln_8(ithr, real(shvec,dp), shmat)
        pft_ref = pft_ref * shmat(:,:self%kfromto(2))
        call self%rotate_ref_8(pft_ref, irot, pft_rot_ref)
        pft_rot_ref = pft_rot_ref * real(self%ctfmats(:,:self%kfromto(2),i),dp)
        do k = self%kfromto(1),self%kfromto(2)
            num       = real(sum(pft_rot_ref(:,k)       * conjg(self%pfts_ptcls(:,k,i))),dp)
            sumsqptcl = real(sum(self%pfts_ptcls(:,k,i) * conjg(self%pfts_ptcls(:,k,i))),dp)
            sumsqref  = real(sum(pft_rot_ref(:,k)       * conjg(pft_rot_ref(:,k))),dp)
            denom     = sumsqptcl * sumsqref
            if( denom < 1.d-16 )then
                frc(k) = 0.
            else
                frc(k) = real(num / sqrt(denom))
            endif
        end do
    end subroutine calc_frc

    module subroutine gen_objfun_vals( self, iref, iptcl, shift, vals )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        real(sp),            intent(in)    :: shift(2)
        real(sp),            intent(out)   :: vals(self%nrots)
         select case(self%p_ptr%cc_objfun)
            case(OBJFUN_CC)
                call self%gen_corrs(iref, iptcl, shift, vals)
            case(OBJFUN_EUCLID)
                call self%gen_euclids(iref, iptcl, shift, vals)
        end select
    end subroutine gen_objfun_vals

    module subroutine gen_best_objfun_val( self, iref, iptcl, shift, dist, irot )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        real(sp),            intent(in)    :: shift(2)
        real(sp),            intent(out)   :: dist
        integer,             intent(out)   :: irot
        select case(self%p_ptr%cc_objfun)
            case(OBJFUN_CC)
                THROW_HARD('gen_best_objfun_val not implemented for OBJFUN_CC')
            case(OBJFUN_EUCLID)
                call self%gen_best_euclid_val(iref, iptcl, shift, dist, irot)
        end select
    end subroutine gen_best_objfun_val

    module subroutine gen_prob_objfun_val( self, iref, iptcl, shift, athres_ub, prob_athres, dist, irot, pvec_sorted, sorted_inds )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        real(sp),            intent(in)    :: shift(2)
        real(sp),            intent(in)    :: athres_ub, prob_athres
        real(sp),            intent(out)   :: dist
        integer,             intent(out)   :: irot
        real(sp),            intent(inout) :: pvec_sorted(self%nrots)
        integer,             intent(inout) :: sorted_inds(self%nrots)
        select case(self%p_ptr%cc_objfun)
            case(OBJFUN_CC)
                THROW_HARD('gen_prob_objfun_val not implemented for OBJFUN_CC')
            case(OBJFUN_EUCLID)
                call self%gen_prob_euclid_val(iref, iptcl, shift, athres_ub, prob_athres, dist, irot, pvec_sorted, sorted_inds)
        end select
    end subroutine gen_prob_objfun_val

    module subroutine gen_prob_power_objfun_val( self, iref, iptcl, shift, power, nsample, dist, corr, irot,&
        &pvec_sorted, sorted_inds )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl, nsample
        real(sp),            intent(in)    :: shift(2)
        real(sp),            intent(in)    :: power
        real(sp),            intent(out)   :: dist, corr
        integer,             intent(out)   :: irot
        real(sp),            intent(inout) :: pvec_sorted(self%nrots)
        integer,             intent(inout) :: sorted_inds(self%nrots)
        select case(self%p_ptr%cc_objfun)
            case(OBJFUN_CC)
                THROW_HARD('gen_prob_power_objfun_val not implemented for OBJFUN_CC')
            case(OBJFUN_EUCLID)
                call self%gen_prob_power_euclid_val(iref, iptcl, shift, power, nsample, dist, corr, irot, pvec_sorted, sorted_inds)
        end select
    end subroutine gen_prob_power_objfun_val

    module subroutine gen_corrs( self, iref, iptcl, shift, cc )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(sp),                    intent(out)   :: cc(self%nrots)
        complex(sp), pointer :: pft_ref(:,:), shmat(:,:)
        real(sp) :: shift_mag_sq
        integer  :: i, ithr, k, kk, k0
        logical  :: even, needs_shift
        ithr    =  omp_get_thread_num() + 1
        i       =  self%pinds(iptcl)
        k0      =  self%kfromto(1)
        even    =  self%iseven(i)
        shmat   => self%heap_vars(ithr)%shmat
        ! ========================================================================
        ! Prepare FFT(S.REF) if needed
        ! ========================================================================
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        needs_shift  = shift_mag_sq > SHERRSQ
        if (needs_shift) then
            pft_ref => self%heap_vars(ithr)%pft_ref
            ! Generate shift matrix
            call self%gen_shmat4aln(ithr, shift, shmat)
            ! Shift reference
            if (even) then
                pft_ref = shmat(:,:self%kfromto(2)) * self%pfts_refs_even(:,self%kfromto(1):self%kfromto(2),iref)
            else
                pft_ref = shmat(:,:self%kfromto(2)) * self%pfts_refs_odd( :,self%kfromto(1):self%kfromto(2),iref)
            endif
            ! Calculate FFT(S.REF)
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%cmat2_many(ithr)%c(1:self%pftsz,            kk) =       pft_ref(:,k)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots, kk) = conjg(pft_ref(:,k))
            end do
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
        else
            ! no shift, memoized reference is used
            if( even )then
                self%cmat2_many(ithr)%c(1:self%pftsz+1,1:self%nk) = self%ft_ref_even(:,self%kfromto(1):self%kfromto(2),iref)
            else
                self%cmat2_many(ithr)%c(1:self%pftsz+1,1:self%nk) = self%ft_ref_odd( :,self%kfromto(1):self%kfromto(2),iref)
            endif
        endif
        ! ========================================================================
        ! Single IFFT #1: IFFT( sum_k FT(CTF2) x FT(REF2) )
        ! ========================================================================
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        if (even) then
            do k = self%kfromto(1), self%kfromto(2)
                self%crvec1(ithr)%c = self%crvec1(ithr)%c + self%ft_ctf2(:,k,i) * self%ft_ref2_even(:,k,iref)
            end do
        else
            do k = self%kfromto(1), self%kfromto(2)
                self%crvec1(ithr)%c = self%crvec1(ithr)%c + self%ft_ctf2(:,k,i) * self%ft_ref2_odd(:,k,iref)
            end do
        endif
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        self%drvec(ithr)%r = real(self%crvec1(ithr)%r(1:self%nrots), dp)
        ! ========================================================================
        ! Single IFFT #2: IFFT( sum_k FT(X.CTF) x FT(S.REF)* )
        ! ========================================================================
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        do k = self%kfromto(1), self%kfromto(2)
            kk = k - k0 + 1
            self%crvec1(ithr)%c = self%crvec1(ithr)%c + &
                self%ft_ptcl_ctf(:,k,i) * conjg(self%cmat2_many(ithr)%c(1:self%pftsz+1,kk))
        end do
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        self%heap_vars(ithr)%kcorrs = real(self%crvec1(ithr)%r(1:self%nrots), dp)
        ! Final correlation computation
        self%drvec(ithr)%r = self%drvec(ithr)%r * real(self%sqsums_ptcls(i) * real(2*self%nrots), dp)
        cc = real(self%heap_vars(ithr)%kcorrs / dsqrt(self%drvec(ithr)%r))
    end subroutine gen_corrs

    module subroutine gen_euclids( self, iref, iptcl, shift, euclids )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(sp),                    intent(out)   :: euclids(self%nrots)
        complex(sp), pointer :: shmat(:,:)
        complex(sp) :: c
        real(dp)    :: ptcl_sqsum, wk, A
        real(sp)    :: shift_mag_sq
        integer     :: k, i, ithr, kk, k0, p
        logical     :: even
        ithr         =  omp_get_thread_num() + 1
        i            =  self%pinds(iptcl)
        k0           =  self%kfromto(1)
        even         =  self%iseven(i)
        ptcl_sqsum   =  self%wsqsums_ptcls(i)
        shmat        => self%heap_vars(ithr)%shmat
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        if ( shift_mag_sq > SHERRSQ ) then
            ! Reference is shifted first
            call self%gen_shmat4aln(ithr, shift, shmat)
            ! Shift reference
            if (even) then
                do k = self%kfromto(1), self%kfromto(2)
                    kk = k - k0 + 1
                    self%cmat2_many(ithr)%c(1:self%pftsz,            kk) =       shmat(:,k) * self%pfts_refs_even(:,k,iref)
                    self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots, kk) = conjg(shmat(:,k) * self%pfts_refs_even(:,k,iref))
                end do
            else
                do k = self%kfromto(1), self%kfromto(2)
                    kk = k - k0 + 1
                    self%cmat2_many(ithr)%c(1:self%pftsz,            kk) =       shmat(:,k) * self%pfts_refs_odd(:,k,iref)
                    self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots, kk) = conjg(shmat(:,k) * self%pfts_refs_odd(:,k,iref))
                end do
            endif
            ! Calculate FT(S.REF)
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
            ! sum_k w_k * (FT(CTF2) x FT(REF2) - 2*FT(X.CTF) x FT(S.REF)*)
            self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
            if (even) then
                do k = self%kfromto(1), self%kfromto(2)
                    kk = k - k0 + 1
                    wk = real(k, dp) / real(self%sigma2_noise(k,iptcl), dp)
                    do p = 1,self%pftsz
                        c = self%ft_ctf2(p,k,i) * self%ft_ref2_even(p,k,iref)
                        c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%cmat2_many(ithr)%c(p, kk))
                        self%crvec1(ithr)%c(p) = self%crvec1(ithr)%c(p) + wk * c
                    enddo
                end do
            else
                do k = self%kfromto(1), self%kfromto(2)
                    kk = k - k0 + 1
                    wk = real(k, dp) / real(self%sigma2_noise(k,iptcl), dp)
                    do p = 1,self%pftsz
                        c = self%ft_ctf2(p,k,i) * self%ft_ref2_odd(p,k,iref)
                        c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%cmat2_many(ithr)%c(p, kk))
                        self%crvec1(ithr)%c(p) = self%crvec1(ithr)%c(p) + wk * c
                    enddo
                end do
            endif
        else
            ! The reference is not shifted, memoized FT(REF) can be used
            ! sum_k w_k * (FT(CTF2) x FT(REF2) - 2*FT(X.CTF) x FT(REF)*)
            self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
            if (even) then
                do k = self%kfromto(1), self%kfromto(2)
                    wk = real(k, dp) / real(self%sigma2_noise(k,iptcl), dp)
                    do p = 1,self%pftsz
                        c = self%ft_ctf2(p,k,i) * self%ft_ref2_even(p,k,iref)
                        c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%ft_ref_even(p,k,iref))
                        self%crvec1(ithr)%c(p) = self%crvec1(ithr)%c(p) + wk * c
                    enddo
                end do
            else
                do k = self%kfromto(1), self%kfromto(2)
                    wk = real(k, dp) / real(self%sigma2_noise(k,iptcl), dp)
                    do p = 1,self%pftsz
                        c = self%ft_ctf2(p,k,i) * self%ft_ref2_odd(p,k,iref)
                        c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%ft_ref_odd(p,k,iref))
                        self%crvec1(ithr)%c(p) = self%crvec1(ithr)%c(p) + wk * c
                    enddo
                end do
            endif
        endif
        ! Single IFFT: IFFT( sum_k w_k * (FT(CTF2) x FT(REF2) - 2*FT(X.CTF) x FT(S.REF)*) )
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        ! Normalization & exponentiation
        A = ptcl_sqsum * real(2*self%nrots, dp)
        do p = 1,self%nrots
            euclids(p) = exp(-1.d0 - self%crvec1(ithr)%r(p) / A)
        end do
    end subroutine gen_euclids

    subroutine gen_euclid_crvec( self, iref, iptcl, shift, norm, ithr )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(dp),                    intent(out)   :: norm
        integer,                     intent(out)   :: ithr
        complex(sp), pointer :: shmat(:,:)
        complex(sp) :: c
        real(dp)    :: wk
        real(sp)    :: shift_mag_sq
        integer     :: k, i, kk, k0, p
        ithr         = omp_get_thread_num() + 1
        i            = self%pinds(iptcl)
        k0           = self%kfromto(1)
        norm         = self%wsqsums_ptcls(i) * real(2*self%nrots, dp)
        shmat        => self%heap_vars(ithr)%shmat
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        if( shift_mag_sq > SHERRSQ )then
            call self%gen_shmat4aln(ithr, shift, shmat)
            if( self%iseven(i) )then
                do k = self%kfromto(1), self%kfromto(2)
                    kk = k - k0 + 1
                    self%cmat2_many(ithr)%c(1:self%pftsz,            kk) =       shmat(:,k) * self%pfts_refs_even(:,k,iref)
                    self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots, kk) = conjg(shmat(:,k) * self%pfts_refs_even(:,k,iref))
                enddo
            else
                do k = self%kfromto(1), self%kfromto(2)
                    kk = k - k0 + 1
                    self%cmat2_many(ithr)%c(1:self%pftsz,            kk) =       shmat(:,k) * self%pfts_refs_odd(:,k,iref)
                    self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots, kk) = conjg(shmat(:,k) * self%pfts_refs_odd(:,k,iref))
                enddo
            endif
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
            self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
            if( self%iseven(i) )then
                do k = self%kfromto(1), self%kfromto(2)
                    kk = k - k0 + 1
                    wk = real(k, dp) / real(self%sigma2_noise(k,iptcl), dp)
                    do p = 1,self%pftsz
                        c = self%ft_ctf2(p,k,i) * self%ft_ref2_even(p,k,iref)
                        c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%cmat2_many(ithr)%c(p,kk))
                        self%crvec1(ithr)%c(p) = self%crvec1(ithr)%c(p) + wk * c
                    enddo
                enddo
            else
                do k = self%kfromto(1), self%kfromto(2)
                    kk = k - k0 + 1
                    wk = real(k, dp) / real(self%sigma2_noise(k,iptcl), dp)
                    do p = 1,self%pftsz
                        c = self%ft_ctf2(p,k,i) * self%ft_ref2_odd(p,k,iref)
                        c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%cmat2_many(ithr)%c(p,kk))
                        self%crvec1(ithr)%c(p) = self%crvec1(ithr)%c(p) + wk * c
                    enddo
                enddo
            endif
        else
            self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
            if( self%iseven(i) )then
                do k = self%kfromto(1), self%kfromto(2)
                    wk = real(k, dp) / real(self%sigma2_noise(k,iptcl), dp)
                    do p = 1,self%pftsz
                        c = self%ft_ctf2(p,k,i) * self%ft_ref2_even(p,k,iref)
                        c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%ft_ref_even(p,k,iref))
                        self%crvec1(ithr)%c(p) = self%crvec1(ithr)%c(p) + wk * c
                    enddo
                enddo
            else
                do k = self%kfromto(1), self%kfromto(2)
                    wk = real(k, dp) / real(self%sigma2_noise(k,iptcl), dp)
                    do p = 1,self%pftsz
                        c = self%ft_ctf2(p,k,i) * self%ft_ref2_odd(p,k,iref)
                        c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%ft_ref_odd(p,k,iref))
                        self%crvec1(ithr)%c(p) = self%crvec1(ithr)%c(p) + wk * c
                    enddo
                enddo
            endif
        endif
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
    end subroutine gen_euclid_crvec

    real(sp) function euclid_dist_from_crvec( self, ithr, irot, norm ) result(dist)
        class(polarft_calc), intent(in) :: self
        integer,             intent(in) :: ithr, irot
        real(dp),            intent(in) :: norm
        real(dp) :: v
        v = 1.d0 + real(self%crvec1(ithr)%r(irot), dp) / norm
        if( v > -log(real(TINY,dp)) )then
            dist = huge(dist)
        else
            dist = real(v,sp)
        endif
    end function euclid_dist_from_crvec

    module subroutine gen_best_euclid_val( self, iref, iptcl, shift, dist, irot )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(sp),                    intent(out)   :: dist
        integer,                     intent(out)   :: irot
        real(dp) :: norm
        real(sp) :: dist_tmp
        integer  :: ithr, p
        call gen_euclid_crvec(self, iref, iptcl, shift, norm, ithr)
        irot = 1
        dist = euclid_dist_from_crvec(self, ithr, 1, norm)
        do p = 2,self%nrots
            dist_tmp = euclid_dist_from_crvec(self, ithr, p, norm)
            if( dist_tmp < dist )then
                irot = p
                dist = dist_tmp
            endif
        enddo
    end subroutine gen_best_euclid_val

    module subroutine gen_prob_euclid_val( self, iref, iptcl, shift, athres_ub, prob_athres, dist, irot, pvec_sorted, sorted_inds )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(sp),                    intent(in)    :: athres_ub, prob_athres
        real(sp),                    intent(out)   :: dist
        integer,                     intent(out)   :: irot
        real(sp),                    intent(inout) :: pvec_sorted(self%nrots)
        integer,                     intent(inout) :: sorted_inds(self%nrots)
        real(dp)    :: norm
        integer     :: ithr
        call gen_euclid_crvec(self, iref, iptcl, shift, norm, ithr)
        call sample_bounded_dist(self%nrots, euclid_dist_at_rot, athres_ub, prob_athres, dist, irot,&
            &pvec_sorted, sorted_inds)

    contains

        real function euclid_dist_at_rot(p_loc) result(dist_loc)
            integer,  intent(in) :: p_loc
            dist_loc = euclid_dist_from_crvec(self, ithr, p_loc, norm)
        end function euclid_dist_at_rot

    end subroutine gen_prob_euclid_val

    module subroutine gen_prob_power_euclid_val( self, iref, iptcl, shift, power, nsample, dist, corr, irot,&
        &pvec_sorted, sorted_inds )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl, nsample
        real(sp),                    intent(in)    :: shift(2)
        real(sp),                    intent(in)    :: power
        real(sp),                    intent(out)   :: dist, corr
        integer,                     intent(out)   :: irot
        real(sp),                    intent(inout) :: pvec_sorted(self%nrots)
        integer,                     intent(inout) :: sorted_inds(self%nrots)
        real(dp) :: norm
        integer  :: ithr
        call gen_euclid_crvec(self, iref, iptcl, shift, norm, ithr)
        call sample_power_dist(self%nrots, euclid_dist_at_rot, power, nsample, dist, corr, irot,&
            &pvec_sorted, sorted_inds)

    contains

        real function euclid_dist_at_rot(p_loc) result(dist_loc)
            integer, intent(in) :: p_loc
            dist_loc = euclid_dist_from_crvec(self, ithr, p_loc, norm)
        end function euclid_dist_at_rot

    end subroutine gen_prob_power_euclid_val

    module real(dp) function gen_corr_for_rot_8_1( self, iref, iptcl, irot )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl, irot
        complex(dp), pointer :: pft_ref_8(:,:)
        integer :: ithr, i
        i         = self%pinds(iptcl)
        ithr      = omp_get_thread_num() + 1
        pft_ref_8 => self%heap_vars(ithr)%pft_ref_8
        if (self%iseven(i)) then
            pft_ref_8 = self%pfts_refs_even(:,self%kfromto(1):self%kfromto(2),iref)
        else
            pft_ref_8 = self%pfts_refs_odd(:,self%kfromto(1):self%kfromto(2),iref)
        endif
        gen_corr_for_rot_8_1 = 0.d0
        select case(self%p_ptr%cc_objfun)
            case(OBJFUN_CC)
                gen_corr_for_rot_8_1 = self%gen_corr_cc_for_rot_8_1(pft_ref_8, i, irot)
            case(OBJFUN_EUCLID)
                gen_corr_for_rot_8_1 = self%gen_euclid_for_rot_8_1(pft_ref_8, iptcl, irot)
        end select
    end function gen_corr_for_rot_8_1

    module real(dp) function gen_corr_for_rot_8_2( self, iref, iptcl, shvec, irot )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(dp),                    intent(in)    :: shvec(2)
        integer,                     intent(in)    :: irot
        complex(dp), pointer :: pft_ref_8(:,:)
        integer :: ithr, i
        i         = self%pinds(iptcl)
        ithr      = omp_get_thread_num() + 1
        pft_ref_8 => self%heap_vars(ithr)%pft_ref_8
        if (self%iseven(i)) then
            pft_ref_8 = self%pfts_refs_even(:,self%kfromto(1):self%kfromto(2),iref)
        else
            pft_ref_8 = self%pfts_refs_odd(:,self%kfromto(1):self%kfromto(2),iref)
        endif
        gen_corr_for_rot_8_2 = 0.d0
        select case(self%p_ptr%cc_objfun)
            case(OBJFUN_CC)
                gen_corr_for_rot_8_2 = self%gen_corr_cc_for_rot_8_2(pft_ref_8, i, shvec, irot)
            case(OBJFUN_EUCLID)
                gen_corr_for_rot_8_2 = self%gen_euclid_for_rot_8_2(pft_ref_8, iptcl, shvec, irot)
        end select
    end function gen_corr_for_rot_8_2

    module real(dp) function gen_corr_cc_for_rot_8_1( self, pft_ref, i, irot )
        class(polarft_calc),  intent(inout) :: self
        complex(dp), pointer, intent(inout) :: pft_ref(:,:)
        integer,              intent(in)    :: i, irot
        complex(dp) :: crefctf
        real(dp)    :: sqsum_ref, sqsumk, wk, fk
        integer     :: k, p, rp
        sqsum_ref = 0.d0
        gen_corr_cc_for_rot_8_1 = 0.d0
        if( irot <= self%pftsz )then
            do k = self%kfromto(1),self%kfromto(2)
                sqsumk = 0.d0
                wk     = real(k, dp)
                fk     = 0.d0
                do p = 1, irot-1
                    ! rotation index
                    rp = p + self%pftsz - irot + 1
                    ! Rot(REFp).CTFp x PTCLp*
                    crefctf = conjg(pft_ref(rp,k)) * real(self%ctfmats(p,k,i),dp)
                    fk      = fk     + real(crefctf * conjg(self%pfts_ptcls(p,k,i)), dp)
                    ! |Rot(REFp).CTFp|^2
                    sqsumk  = sqsumk + real(crefctf * conjg(crefctf), dp)
                enddo
                do p = irot, self%pftsz
                    rp = p - irot + 1
                    crefctf = pft_ref(rp,k) * real(self%ctfmats(p,k,i),dp)
                    fk      = fk     + real(crefctf * conjg(self%pfts_ptcls(p,k,i)), dp)
                    sqsumk  = sqsumk + real(crefctf * conjg(crefctf), dp)
                enddo
                gen_corr_cc_for_rot_8_1 = gen_corr_cc_for_rot_8_1 + wk * fk
                sqsum_ref = sqsum_ref + wk * sqsumk
            enddo
        else
            do k = self%kfromto(1),self%kfromto(2)
                sqsumk = 0.d0
                wk     = real(k, dp)
                fk     = 0.d0
                do p = 1, irot-self%pftsz-1
                    rp = p + self%nrots - irot + 1
                    crefctf = pft_ref(rp,k) * real(self%ctfmats(p,k,i),dp)
                    fk      = fk     + real(crefctf * conjg(self%pfts_ptcls(p,k,i)), dp)
                    sqsumk  = sqsumk + real(crefctf * conjg(crefctf), dp)
                enddo
                do p = irot-self%pftsz, self%pftsz
                    rp = p - irot + self%pftsz + 1
                    crefctf = conjg(pft_ref(rp,k)) * real(self%ctfmats(p,k,i),dp)
                    fk      = fk     + real(crefctf * conjg(self%pfts_ptcls(p,k,i)), dp)
                    sqsumk  = sqsumk + real(crefctf * conjg(crefctf), dp)
                enddo
                gen_corr_cc_for_rot_8_1 = gen_corr_cc_for_rot_8_1 + wk * fk
                sqsum_ref = sqsum_ref + wk * sqsumk
            enddo
        endif
        gen_corr_cc_for_rot_8_1 = gen_corr_cc_for_rot_8_1 / sqrt(sqsum_ref * self%ksqsums_ptcls(i))
    end function gen_corr_cc_for_rot_8_1

    module real(dp) function gen_corr_cc_for_rot_8_2( self, pft_ref, i, shvec, irot )
        class(polarft_calc), target, intent(inout) :: self
        complex(dp),        pointer, intent(inout) :: pft_ref(:,:)
        integer,                     intent(in)    :: i, irot
        real(dp),                    intent(in)    :: shvec(2)
        complex(dp), pointer :: shmat_8(:,:)
        complex(dp) :: crefctf
        real(dp)    :: sqsum_ref, sqsumk, wk, fk
        integer     :: k, p, rp, ithr
        ithr      = omp_get_thread_num() + 1
        sqsum_ref = 0.d0
        shmat_8 => self%heap_vars(ithr)%shmat_8
        ! shift matrix
        call self%gen_shmat4aln_8(ithr, shvec, shmat_8)
        ! evaluation
        gen_corr_cc_for_rot_8_2 = 0.d0
        if( irot <= self%pftsz )then
            do k = self%kfromto(1),self%kfromto(2)
                sqsumk = 0.d0
                wk     = real(k, dp)
                fk     = 0.d0
                do p = 1, irot-1
                    ! rotation index
                    rp = p + self%pftsz - irot + 1
                    ! Rot(REFp).CTFp x PTCLp*
                    crefctf = conjg(pft_ref(rp,k) * shmat_8(rp,k)) * real(self%ctfmats(p,k,i),dp)
                    fk      = fk     + real(crefctf * conjg(self%pfts_ptcls(p,k,i)), dp)
                    ! |Rot(REFp).CTFp|^2
                    sqsumk  = sqsumk + real(crefctf * conjg(crefctf), dp)
                enddo
                do p = irot, self%pftsz
                    rp = p - irot + 1
                    crefctf = pft_ref(rp,k) * shmat_8(rp,k) * real(self%ctfmats(p,k,i),dp)
                    fk      = fk     + real(crefctf * conjg(self%pfts_ptcls(p,k,i)), dp)
                    sqsumk  = sqsumk + real(crefctf * conjg(crefctf), dp)
                enddo
                gen_corr_cc_for_rot_8_2 = gen_corr_cc_for_rot_8_2 + wk * fk
                sqsum_ref = sqsum_ref + wk * sqsumk
            enddo
        else
            do k = self%kfromto(1),self%kfromto(2)
                sqsumk = 0.d0
                wk     = real(k, dp)
                fk     = 0.d0
                do p = 1, irot-self%pftsz-1
                    rp = p + self%nrots - irot + 1
                    crefctf = pft_ref(rp,k) * shmat_8(rp,k) * real(self%ctfmats(p,k,i),dp)
                    fk      = fk     + real(crefctf * conjg(self%pfts_ptcls(p,k,i)), dp)
                    sqsumk  = sqsumk + real(crefctf * conjg(crefctf), dp)
                enddo
                do p = irot-self%pftsz, self%pftsz
                    rp = p - irot + self%pftsz + 1
                    crefctf = conjg(pft_ref(rp,k) * shmat_8(rp,k)) * real(self%ctfmats(p,k,i),dp)
                    fk      = fk     + real(crefctf * conjg(self%pfts_ptcls(p,k,i)), dp)
                    sqsumk  = sqsumk + real(crefctf * conjg(crefctf), dp)
                enddo
                gen_corr_cc_for_rot_8_2 = gen_corr_cc_for_rot_8_2 + wk * fk
                sqsum_ref = sqsum_ref + wk * sqsumk
            enddo
        endif
        gen_corr_cc_for_rot_8_2 = gen_corr_cc_for_rot_8_2 / sqrt(sqsum_ref * self%ksqsums_ptcls(i))
    end function gen_corr_cc_for_rot_8_2

    module real(dp) function gen_euclid_for_rot_8_1( self, pft_ref, iptcl, irot )
        class(polarft_calc), target, intent(inout) :: self
        complex(dp),        pointer, intent(inout) :: pft_ref(:,:)
        integer,                     intent(in)    :: iptcl, irot
        complex(dp) :: c
        real(dp)    :: sumsq, wk
        integer     :: ithr, p, rp, k, i
        ithr = omp_get_thread_num() + 1
        i = self%pinds(iptcl)
        gen_euclid_for_rot_8_1 = 0.d0
        if( irot <= self%pftsz )then
            do k = self%kfromto(1),self%kfromto(2)
                wk    = real(k, dp) / real(self%sigma2_noise(k,iptcl), dp)
                sumsq = 0.d0
                do p = 1, irot-1
                    ! rotation index
                    rp = p + self%pftsz - irot + 1
                    !  Rot(REFp).CTFp - PTCLp
                    c  = conjg(pft_ref(rp,k)) * real(self%ctfmats(p,k,i),dp) - cmplx(self%pfts_ptcls(p,k,i),kind=dp)
                    ! SUMp | Rot(REFp).CTFp - PTCLp |^2
                    sumsq = sumsq + real(c * conjg(c), dp)
                enddo
                do p = irot, self%pftsz
                    rp = p - irot + 1
                    c  = pft_ref(rp,k) * real(self%ctfmats(p,k,i),dp) - cmplx(self%pfts_ptcls(p,k,i),kind=dp)
                    sumsq = sumsq + real(c * conjg(c), dp)
                enddo
                ! DIST^2 = SUMk weightk SUMp |Rot(REFkp).CTFkp - PTCLkp|^2
                gen_euclid_for_rot_8_1 = gen_euclid_for_rot_8_1 + wk * sumsq
            enddo
        else
            do k = self%kfromto(1),self%kfromto(2)
                wk    = real(k, dp) / real(self%sigma2_noise(k,iptcl), dp)
                sumsq = 0.d0
                do p = 1, irot-self%pftsz-1
                    rp = p + self%nrots - irot + 1
                    c  = pft_ref(rp,k) * real(self%ctfmats(p,k,i),dp) - cmplx(self%pfts_ptcls(p,k,i),kind=dp)
                    sumsq = sumsq + real(c * conjg(c), dp)
                enddo
                do p = irot-self%pftsz, self%pftsz
                    rp = p - irot + self%pftsz + 1
                    c  = conjg(pft_ref(rp,k)) * real(self%ctfmats(p,k,i),dp) - cmplx(self%pfts_ptcls(p,k,i),kind=dp)
                    sumsq = sumsq + real(c * conjg(c), dp)
                enddo
                ! accumulate weighted sum over k
                gen_euclid_for_rot_8_1 = gen_euclid_for_rot_8_1 + wk * sumsq
            enddo
        endif
        ! EUCL = EXP( -DIST^2 / |PTCL|^2 )
        gen_euclid_for_rot_8_1 = exp(-gen_euclid_for_rot_8_1 / self%wsqsums_ptcls(i))
    end function gen_euclid_for_rot_8_1

    module real(dp) function gen_euclid_for_rot_8_2( self, pft_ref, iptcl, shvec, irot )
        class(polarft_calc), target, intent(inout) :: self
        complex(dp),        pointer, intent(inout) :: pft_ref(:,:)
        integer,                     intent(in)    :: iptcl, irot
        real(dp),                    intent(in)    :: shvec(2)
        complex(dp), pointer :: shmat_8(:,:)
        complex(dp) :: c
        real(dp)    :: sumsq, wk
        integer     :: ithr, p, rp, k, i
        ithr    = omp_get_thread_num() + 1
        i       = self%pinds(iptcl)
        shmat_8 => self%heap_vars(ithr)%shmat_8
        ! shift matrix
        call self%gen_shmat4aln_8(ithr, shvec, shmat_8)
        ! splitting both the compute based on irot and the loop avoids branching,
        ! and optimizes for memory access patterns & vectorization
        gen_euclid_for_rot_8_2 = 0.d0
        if( irot <= self%pftsz )then
            do k = self%kfromto(1),self%kfromto(2)
                wk    = real(k, dp) / real(self%sigma2_noise(k,iptcl), dp)
                sumsq = 0.d0
                do p = 1, irot-1
                    ! rotation index
                    rp = p + self%pftsz - irot + 1
                    ! Rot(Shift(REFp))
                    c  = conjg(pft_ref(rp,k) * shmat_8(rp,k))
                    !  Rot(Shift(REFp)).CTFp - PTCLp
                    c  = c * real(self%ctfmats(p,k,i),dp) - cmplx(self%pfts_ptcls(p,k,i),kind=dp)
                    ! SUMp | Rot(Shift(REFp)).CTFp - PTCLp |^2
                    sumsq = sumsq + real(c * conjg(c), dp)
                enddo
                do p = irot, self%pftsz
                    rp = p - irot + 1
                    c  = pft_ref(rp,k) * shmat_8(rp,k)
                    c  = c * real(self%ctfmats(p,k,i),dp) - cmplx(self%pfts_ptcls(p,k,i),kind=dp)
                    sumsq = sumsq + real(c * conjg(c), dp)
                enddo
                ! DIST^2 = SUMk weightk SUMp |Rot(Shift(REFkp)).CTFkp - PTCLkp|^2
                gen_euclid_for_rot_8_2 = gen_euclid_for_rot_8_2 + wk * sumsq
            enddo
        else
            do k = self%kfromto(1),self%kfromto(2)
                wk    = real(k, dp) / real(self%sigma2_noise(k,iptcl), dp)
                sumsq = 0.d0
                do p = 1, irot-self%pftsz-1
                    rp = p + self%nrots - irot + 1
                    c  = pft_ref(rp,k) * shmat_8(rp,k)
                    c  = c * real(self%ctfmats(p,k,i),dp) - cmplx(self%pfts_ptcls(p,k,i),kind=dp)
                    sumsq = sumsq + real(c * conjg(c), dp)
                enddo
                do p = irot-self%pftsz, self%pftsz
                    rp = p - irot + self%pftsz + 1
                    c  = conjg(pft_ref(rp,k) * shmat_8(rp,k))
                    c  = c * real(self%ctfmats(p,k,i),dp) - cmplx(self%pfts_ptcls(p,k,i),kind=dp)
                    sumsq = sumsq + real(c * conjg(c), dp)
                enddo
                ! accumulate weighted sum over k
                gen_euclid_for_rot_8_2 = gen_euclid_for_rot_8_2 + wk * sumsq
            enddo
        endif
        ! EUCL = EXP( -DIST^2 / |PTCL|^2 )
        gen_euclid_for_rot_8_2 = exp(-gen_euclid_for_rot_8_2 / self%wsqsums_ptcls(i))
    end function gen_euclid_for_rot_8_2

    module subroutine gen_corr_grad_for_rot_8( self, iref, iptcl, shvec, irot, f, grad )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(dp),                    intent(in)    :: shvec(2)
        integer,                     intent(in)    :: irot
        real(dp),                    intent(out)   :: f, grad(2)
        complex(dp), pointer :: pft_ref_8(:,:)
        integer :: ithr, i
        i         = self%pinds(iptcl)
        ithr      = omp_get_thread_num() + 1
        pft_ref_8 => self%heap_vars(ithr)%pft_ref_8
        if (self%iseven(i)) then
            pft_ref_8 = self%pfts_refs_even(:,self%kfromto(1):self%kfromto(2),iref)
        else
            pft_ref_8 = self%pfts_refs_odd(:,self%kfromto(1):self%kfromto(2),iref)
        endif
        select case(self%p_ptr%cc_objfun)
            case(OBJFUN_CC)
                call self%gen_corr_cc_grad_for_rot_8(pft_ref_8, i, shvec, irot, f, grad)
            case(OBJFUN_EUCLID)
                call self%gen_euclid_grad_for_rot_8(pft_ref_8, iptcl, shvec, irot, f, grad)
        end select
    end subroutine gen_corr_grad_for_rot_8

    module subroutine gen_corr_cc_grad_for_rot_8( self, pft_ref, i, shvec, irot, f, grad )
        class(polarft_calc), target, intent(inout) :: self
        complex(dp), pointer,        intent(inout) :: pft_ref(:,:)
        integer,                     intent(in)    :: i, irot
        real(dp),                    intent(in)    :: shvec(2)
        real(dp),                    intent(out)   :: f, grad(2)
        real(dp),    pointer :: argtransf(:,:)
        complex(dp), pointer :: shmat_8(:,:)
        complex(dp) :: crefctf, cdiff, cg, cptcl
        real(dp)    :: fk, wk, gkx, gky, sqsumk, sqsum, denom
        integer     :: k, ithr, p, rp
        ithr      = omp_get_thread_num() + 1
        denom     = self%wsqsums_ptcls(i)
        f         = 0.d0
        grad      = 0.d0
        sqsum     = 0.d0
        shmat_8   => self%heap_vars(ithr)%shmat_8
        argtransf => self%argtransf
        ! shift matrix
        call self%gen_shmat4aln_8(ithr, shvec, shmat_8)
        if( irot <= self%pftsz )then
            do k = self%kfromto(1),self%kfromto(2)
                fk     = 0.d0
                sqsumk = 0.d0
                gkx    = 0.d0
                gky    = 0.d0
                do p = 1, irot-1
                    ! rotation index
                    rp      = p + self%pftsz - irot + 1
                    ! Rot(Shift(REF)).CTF.PTCL*
                    crefctf = conjg(pft_ref(rp,k) * shmat_8(rp,k)) * real(self%ctfmats(p,k,i), dp)
                    cptcl   = conjg(cmplx(self%pfts_ptcls(p,k,i), kind=dp))
                    fk      = fk     + real(crefctf * cptcl, dp)
                    ! |Rot(Shift(REF)).CTF|^2
                    sqsumk  = sqsumk + real(crefctf * conjg(crefctf), dp)
                    ! gradient along x
                    cg      = cmplx(0.d0, -argtransf(rp,k), kind=dp) * crefctf
                    gkx     = gkx + real(cg * cptcl, dp)
                    ! gradient along y
                    cg      = cmplx(0.d0, -argtransf(self%pftsz+rp,k), kind=dp) * crefctf
                    gky     = gky + real(cg * cptcl, dp)
                enddo
                do p = irot, self%pftsz
                    rp      = p - irot + 1
                    crefctf = pft_ref(rp,k) * shmat_8(rp,k) * real(self%ctfmats(p,k,i), dp)
                    cptcl   = conjg(cmplx(self%pfts_ptcls(p,k,i), kind=dp))
                    fk      = fk     + real(crefctf * cptcl, dp)
                    sqsumk  = sqsumk + real(crefctf * conjg(crefctf), dp)
                    cg      = cmplx(0.d0, argtransf(rp,k), kind=dp) * crefctf
                    gkx     = gkx + real(cg * cptcl, dp)
                    cg      = cmplx(0.d0, argtransf(self%pftsz+rp,k), kind=dp) * crefctf
                    gky     = gky + real(cg * cptcl, dp)
                enddo
                wk      = real(k, dp)
                f       = f       + wk * fk
                sqsum   = sqsum   + wk * sqsumk
                grad(1) = grad(1) + wk * gkx
                grad(2) = grad(2) + wk * gky
            enddo
        else
            do k = self%kfromto(1),self%kfromto(2)
                fk     = 0.d0
                sqsumk = 0.d0
                gkx    = 0.d0
                gky    = 0.d0
                do p = 1, irot-self%pftsz-1
                    rp      = p + self%nrots - irot + 1
                    crefctf = pft_ref(rp,k) * shmat_8(rp,k) * real(self%ctfmats(p,k,i), dp)
                    cptcl   = conjg(cmplx(self%pfts_ptcls(p,k,i), kind=dp))
                    fk      = fk     + real(crefctf * cptcl, dp)
                    sqsumk  = sqsumk + real(crefctf * conjg(crefctf), dp)
                    cg      = cmplx(0.d0, argtransf(rp,k), kind=dp) * crefctf
                    gkx     = gkx + real(cg * cptcl, dp)
                    cg      = cmplx(0.d0, argtransf(self%pftsz+rp,k), kind=dp) * crefctf
                    gky     = gky + real(cg * cptcl, dp)
                enddo
                do p = irot-self%pftsz, self%pftsz
                    rp     = p - irot + self%pftsz + 1
                    crefctf = conjg(pft_ref(rp,k) * shmat_8(rp,k)) * real(self%ctfmats(p,k,i), dp)
                    cptcl   = conjg(cmplx(self%pfts_ptcls(p,k,i), kind=dp))
                    fk      = fk     + real(crefctf * cptcl, dp)
                    sqsumk  = sqsumk + real(crefctf * conjg(crefctf), dp)
                    cg      = cmplx(0.d0, -argtransf(rp,k), kind=dp) * crefctf
                    gkx     = gkx + real(cg * cptcl, dp)
                    cg      = cmplx(0.d0, -argtransf(self%pftsz+rp,k), kind=dp) * crefctf
                    gky     = gky + real(cg * cptcl, dp)
                enddo
                wk      = real(k, dp)
                f       = f       + wk * fk
                sqsum   = sqsum   + wk * sqsumk
                grad(1) = grad(1) + wk * gkx
                grad(2) = grad(2) + wk * gky
            enddo
        endif
        denom = sqrt(sqsum * self%ksqsums_ptcls(i))
        f     = f    / denom
        grad  = grad / denom
    end subroutine gen_corr_cc_grad_for_rot_8

    module subroutine gen_euclid_grad_for_rot_8( self, pft_ref, iptcl, shvec, irot, f, grad )
        class(polarft_calc),  target, intent(inout) :: self
        complex(dp),         pointer, intent(inout) :: pft_ref(:,:)
        integer,                      intent(in)    :: iptcl, irot
        real(dp),                     intent(in)    :: shvec(2)
        real(dp),                     intent(out)   :: f, grad(2)
        real(dp),    pointer :: argtransf(:,:)
        complex(dp), pointer :: shmat_8(:,:)
        complex(dp) :: crefctf, cdiff, cg
        real(dp)    :: fk, wk, gkx, gky, denom
        integer     :: k, i, ithr, p, rp
        ithr      = omp_get_thread_num() + 1
        i         = self%pinds(iptcl)
        denom     = self%wsqsums_ptcls(i)
        f         = 0.d0
        grad      = 0.d0
        shmat_8   => self%heap_vars(ithr)%shmat_8
        argtransf => self%argtransf
        ! shift matrix
        call self%gen_shmat4aln_8(ithr, shvec, shmat_8)
        ! splitting both the compute based on irot and the loop avoids branching,
        ! and optimizes for memory access patterns & vectorization
        if( irot <= self%pftsz )then
            do k = self%kfromto(1),self%kfromto(2)
                fk    = 0.d0
                gkx   = 0.d0
                gky   = 0.d0
                do p = 1, irot-1
                    rp = p + self%pftsz - irot + 1
                    ! |Rot(Shift(REF)).CTF - PTCL|^2
                    crefctf = conjg(pft_ref(rp,k) * shmat_8(rp,k)) * real(self%ctfmats(p,k,i), dp)
                    cdiff   = crefctf - self%pfts_ptcls(p,k,i)
                    fk      = fk + real(cdiff * conjg(cdiff), dp)
                    ! gradient along x
                    cg  = cmplx(0.d0, -argtransf(rp,k), kind=dp) * crefctf
                    gkx = gkx + real(cg * conjg(cdiff), dp)
                    ! gradient along y
                    cg  = cmplx(0.d0, -argtransf(self%pftsz+rp,k), kind=dp) * crefctf
                    gky = gky + real(cg * conjg(cdiff), dp)
                enddo
                do p = irot, self%pftsz
                    rp = p - irot + 1
                    crefctf = pft_ref(rp,k) * shmat_8(rp,k) * real(self%ctfmats(p,k,i), dp)
                    cdiff   = crefctf - self%pfts_ptcls(p,k,i)
                    fk      = fk + real(cdiff * conjg(cdiff), dp)
                    cg  = cmplx(0.d0, argtransf(rp,k), kind=dp) * crefctf
                    gkx = gkx + real(cg * conjg(cdiff), dp)
                    cg  = cmplx(0.d0, argtransf(self%pftsz+rp,k), kind=dp) * crefctf
                    gky = gky + real(cg * conjg(cdiff), dp)
                enddo
                wk      = real(k, dp) / real(self%sigma2_noise(k,iptcl), dp)
                f       = f       + wk * fk
                grad(1) = grad(1) + wk * gkx
                grad(2) = grad(2) + wk * gky
            enddo
        else
            do k = self%kfromto(1),self%kfromto(2)
                fk    = 0.d0
                gkx   = 0.d0
                gky   = 0.d0
                do p = 1, irot-self%pftsz-1
                    rp = p + self%nrots - irot + 1
                    crefctf = pft_ref(rp,k) * shmat_8(rp,k) * real(self%ctfmats(p,k,i), dp)
                    cdiff   = crefctf - self%pfts_ptcls(p,k,i)
                    fk      = fk + real(cdiff * conjg(cdiff), dp)
                    cg  = cmplx(0.d0, argtransf(rp,k), kind=dp) * crefctf
                    gkx = gkx + real(cg * conjg(cdiff), dp)
                    cg  = cmplx(0.d0, argtransf(self%pftsz+rp,k), kind=dp) * crefctf
                    gky = gky + real(cg * conjg(cdiff), dp)
                enddo
                do p = irot-self%pftsz, self%pftsz
                    rp = p - irot + self%pftsz + 1
                    crefctf = conjg(pft_ref(rp,k) * shmat_8(rp,k)) * real(self%ctfmats(p,k,i), dp)
                    cdiff   = crefctf - self%pfts_ptcls(p,k,i)
                    fk      = fk + real(cdiff * conjg(cdiff), dp)
                    cg  = cmplx(0.d0, -argtransf(rp,k), kind=dp) * crefctf
                    gkx = gkx + real(cg * conjg(cdiff), dp)
                    cg  = cmplx(0.d0, -argtransf(self%pftsz+rp,k), kind=dp) * crefctf
                    gky = gky + real(cg * conjg(cdiff), dp)
                enddo
                wk      = real(k, dp) / real(self%sigma2_noise(k,iptcl), dp)
                f       = f       + wk * fk
                grad(1) = grad(1) + wk * gkx
                grad(2) = grad(2) + wk * gky
            enddo
        endif
        f     = exp(-f / denom)
        grad  = -f * 2.d0 * grad / denom
    end subroutine gen_euclid_grad_for_rot_8

    module subroutine gen_corr_grad_only_for_rot_8( self, iref, iptcl, shvec, irot, grad )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(dp),                    intent(in)    :: shvec(2)
        integer,                     intent(in)    :: irot
        real(dp),                    intent(out)   :: grad(2)
        complex(dp), pointer :: pft_ref_8(:,:)
        real(dp) :: f
        integer  :: ithr, i
        i    = self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
        if (self%iseven(i)) then
            pft_ref_8 = self%pfts_refs_even(:,self%kfromto(1):self%kfromto(2),iref)
        else
            pft_ref_8 = self%pfts_refs_odd(:,self%kfromto(1):self%kfromto(2),iref)
        endif
        select case(self%p_ptr%cc_objfun)
            case(OBJFUN_CC)
                call self%gen_corr_cc_grad_for_rot_8(pft_ref_8, i, shvec, irot, f, grad)
            case(OBJFUN_EUCLID)
                call self%gen_euclid_grad_for_rot_8(pft_ref_8, iptcl, shvec, irot, f, grad)
        end select
    end subroutine gen_corr_grad_only_for_rot_8

    !> Evaluation of sigma2 noise particle contribution
    module subroutine gen_sigma_contrib( self, iref, iptcl, shvec, irot, sigma_contrib)
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shvec(2)
        integer,                     intent(in)    :: irot
        real(sp),  optional,         intent(out)   :: sigma_contrib(self%kfromto(1):self%kfromto(2))
        complex(dp), pointer :: pft_ref_8(:,:), shmat_8(:,:), pft_ref_tmp_8(:,:)
        real(dp)    :: norm_factor
        integer     :: i, ithr
        i    = self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat_8       => self%heap_vars(ithr)%shmat_8
        if (self%iseven(i)) then
            pft_ref_8 = self%pfts_refs_even(:,self%kfromto(1):self%kfromto(2),iref)
        else
            pft_ref_8 = self%pfts_refs_odd(:,self%kfromto(1):self%kfromto(2),iref)
        endif
        ! shift
        call self%gen_shmat4aln_8(ithr, real(shvec,dp), shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8(:,:self%kfromto(2))
        ! rotation
        call self%rotate_ref_8(pft_ref_8, irot, pft_ref_tmp_8)
        ! ctf
        pft_ref_tmp_8 = pft_ref_tmp_8 * real(self%ctfmats(:,:self%kfromto(2),i), dp)
        ! difference
        pft_ref_tmp_8 = pft_ref_tmp_8 - self%pfts_ptcls(:,:self%kfromto(2),i)
        ! sigma2 - pre-compute normalization factor
        norm_factor = 2.d0 * real(self%pftsz, dp)
        if (present(sigma_contrib)) then
            sigma_contrib = real(sum(real(pft_ref_tmp_8 * conjg(pft_ref_tmp_8), dp), dim=1) / norm_factor)
        else
            self%sigma2_noise(self%kfromto(1):self%kfromto(2), iptcl) = &
                real(sum(real(pft_ref_tmp_8 * conjg(pft_ref_tmp_8), dp), dim=1) / norm_factor)
        endif
    end subroutine gen_sigma_contrib

    ! Returns ref vs. ptcl objective function values
    ! and mirrored ref vs. ptcl objective function values
    module subroutine gen_objfun_vals_mirr_vals( self, iref, iptcl, vals, mvals )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        real(sp),            intent(out)   :: vals(self%nrots), mvals(self%nrots)
         select case(self%p_ptr%cc_objfun)
            case(OBJFUN_CC)
                call self%gen_corrs_mirr_corrs(iref, iptcl, vals, mvals)
            case(OBJFUN_EUCLID)
                call self%gen_euclids_mirr_euclids(iref, iptcl, vals, mvals)
        end select
    end subroutine gen_objfun_vals_mirr_vals

    module subroutine gen_euclids_mirr_euclids( self, iref, iptcl, euclids, meuclids )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(out)   :: euclids(self%nrots), meuclids(self%nrots)
        real(dp), pointer :: w_weights(:), sumsq_cache(:)
        real(dp) :: ptcl_sumsq
        integer  :: i, ithr, k, kk, k0
        logical  :: even
        ithr        =  omp_get_thread_num() + 1
        i           =  self%pinds(iptcl)
        k0          =  self%kfromto(1)
        even        =  self%iseven(i)
        w_weights   => self%heap_vars(ithr)%w_weights
        sumsq_cache => self%heap_vars(ithr)%sumsq_cache
        ! Memoized reference
        if( even )then
            self%cmat2_many(ithr)%c(1:self%pftsz+1,1:self%nk) = self%ft_ref_even(:,self%kfromto(1):self%kfromto(2),iref)
        else
            self%cmat2_many(ithr)%c(1:self%pftsz+1,1:self%nk) = self%ft_ref_odd( :,self%kfromto(1):self%kfromto(2),iref)
        endif
        ! Pre-compute weights and particle sums
        do k = self%kfromto(1), self%kfromto(2)
            kk = k - k0 + 1
            ! shell & sigma2 weight
            w_weights(kk)   = real(k, dp) / real(self%sigma2_noise(k,iptcl), dp)
            ! shell & sigma2 weighted particle variance
            sumsq_cache(kk) = sum(real(self%pfts_ptcls(:,k,i) * conjg(self%pfts_ptcls(:,k,i)), dp)) * w_weights(kk)
        end do
        ptcl_sumsq = sum(sumsq_cache)
        ! Single IFFT: IFFT( sum_k w_k * (FT(CTF2) x FT(REF2) - 2*FT(X.CTF) x FT(S.REF)*) )
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        if( even )then
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%crvec1(ithr)%c = self%crvec1(ithr)%c + real(w_weights(kk),c_float) * ( &
                    self%ft_ctf2(:,k,i) * self%ft_ref2_even(:,k,iref) - &
                    2.0 * self%ft_ptcl_ctf(:,k,i) * conjg(self%cmat2_many(ithr)%c(1:self%pftsz+1, kk)) )
            end do
        else
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%crvec1(ithr)%c = self%crvec1(ithr)%c + real(w_weights(kk),c_float) * ( &
                    self%ft_ctf2(:,k,i) * self%ft_ref2_odd(:,k,iref) - &
                    2.0 * self%ft_ptcl_ctf(:,k,i) * conjg(self%cmat2_many(ithr)%c(1:self%pftsz+1, kk)) )
            end do
        endif
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        self%heap_vars(ithr)%kcorrs = real(self%crvec1(ithr)%r(1:self%nrots), dp) / real(2*self%nrots, dp)
        euclids = real(dexp(-(self%heap_vars(ithr)%kcorrs + ptcl_sumsq) / self%wsqsums_ptcls(i)))
        ! Mirrored reference:
        ! 1. Mirror of memoized reference: self%cmat2_many(ithr)%c <- conjg(self%cmat2_many(ithr)%c)
        ! 2. Mirror of memoized reference variance: self%ft_ref2 <- conjg(self%ft_ctf2)
        ! Single IFFT: IFFT( sum_k w_k * (FT(CTF2) x FT(M(REF2)) - 2*FT(X.CTF) x FT(S.M(REF))* ) )
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        if( even )then
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%crvec1(ithr)%c = self%crvec1(ithr)%c + real(w_weights(kk),c_float) * ( &
                    self%ft_ctf2(:,k,i) * conjg(self%ft_ref2_even(:,k,iref)) - &
                    2.0 * self%ft_ptcl_ctf(:,k,i) * self%cmat2_many(ithr)%c(1:self%pftsz+1, kk) )
            end do
        else
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%crvec1(ithr)%c = self%crvec1(ithr)%c + real(w_weights(kk),c_float) * ( &
                    self%ft_ctf2(:,k,i) * conjg(self%ft_ref2_odd(:,k,iref)) - &
                    2.0 * self%ft_ptcl_ctf(:,k,i) * self%cmat2_many(ithr)%c(1:self%pftsz+1, kk) )
            end do
        endif
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        self%heap_vars(ithr)%kcorrs = real(self%crvec1(ithr)%r(1:self%nrots), dp) / real(2*self%nrots, dp)
        meuclids = real(dexp(-self%heap_vars(ithr)%kcorrs / self%wsqsums_ptcls(i)))
    end subroutine gen_euclids_mirr_euclids

    module subroutine gen_corrs_mirr_corrs( self, iref, iptcl, ccs, mccs )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(out)   :: ccs(self%nrots), mccs(self%nrots)
        integer  :: i, ithr, k, kk, k0
        logical  :: even
        ithr    =  omp_get_thread_num() + 1
        i       =  self%pinds(iptcl)
        k0      =  self%kfromto(1)
        even    =  self%iseven(i)
        ! Memoized reference
        if( even )then
            self%cmat2_many(ithr)%c(1:self%pftsz+1,1:self%nk) = self%ft_ref_even(:,self%kfromto(1):self%kfromto(2),iref)
        else
            self%cmat2_many(ithr)%c(1:self%pftsz+1,1:self%nk) = self%ft_ref_odd( :,self%kfromto(1):self%kfromto(2),iref)
        endif
        ! Single IFFT #1: IFFT( sum_k FT(CTF2) x FT(REF2) )
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        if (even) then
            do k = self%kfromto(1), self%kfromto(2)
                self%crvec1(ithr)%c = self%crvec1(ithr)%c + self%ft_ctf2(:,k,i) * self%ft_ref2_even(:,k,iref)
            end do
        else
            do k = self%kfromto(1), self%kfromto(2)
                self%crvec1(ithr)%c = self%crvec1(ithr)%c + self%ft_ctf2(:,k,i) * self%ft_ref2_odd(:,k,iref)
            end do
        endif
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        ! Denominator is mirror invariant and re-used below
        self%drvec(ithr)%r = real(self%crvec1(ithr)%r(1:self%nrots), dp)
        ! Single IFFT #2: IFFT( sum_k FT(X.CTF) x FT(S.REF)* )
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        do k = self%kfromto(1), self%kfromto(2)
            kk = k - k0 + 1
            self%crvec1(ithr)%c = self%crvec1(ithr)%c + &
                self%ft_ptcl_ctf(:,k,i) * conjg(self%cmat2_many(ithr)%c(1:self%pftsz+1, kk))
        end do
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        self%heap_vars(ithr)%kcorrs = real(self%crvec1(ithr)%r(1:self%nrots), dp)
        ! Final correlation computation
        self%drvec(ithr)%r = sqrt(self%drvec(ithr)%r * real(self%sqsums_ptcls(i) * real(2*self%nrots), dp))
        ccs = real(self%heap_vars(ithr)%kcorrs / self%drvec(ithr)%r)
        ! Mirrored reference
        ! 1. Mirror of memoized reference: self%cmat2_many(ithr)%c <- conjg(self%cmat2_many(ithr)%c)
        ! 2. Mirror of memoized reference variance: self%ft_ref2 <- conjg(self%ft_ctf2)
        ! Single IFFT #2: IFFT( sum_k FT(X.CTF) x FT(S.M(REF))* )
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        do k = self%kfromto(1), self%kfromto(2)
            kk = k - k0 + 1
            self%crvec1(ithr)%c = self%crvec1(ithr)%c + &
                self%ft_ptcl_ctf(:,k,i) * self%cmat2_many(ithr)%c(1:self%pftsz+1, kk)
        end do
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        self%heap_vars(ithr)%kcorrs = real(self%crvec1(ithr)%r(1:self%nrots), dp)
        ! Final correlation computation
        mccs = real(self%heap_vars(ithr)%kcorrs / self%drvec(ithr)%r)
    end subroutine gen_corrs_mirr_corrs

    ! the values are stored in self%crmat_many(ithr)%r( 1:nrots, ind )
    ! ind is not the natural reference index (iref) and is left to the user
    ! to access correctly with get_precalc_objfun_vals( ind, ithr, vals),
    ! the vector irefs maps ind to iref.
    ! ithr book-keeping is also left to the user
    ! TODO: add logical mask to handle empty classes
    module subroutine gen_many_objfun_vals( self, nr, irefs, iptcl, shift )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: nr
        integer,             intent(in)    :: irefs(nr)
        integer,             intent(in)    :: iptcl
        real(sp),            intent(in)    :: shift(2)
         select case(self%p_ptr%cc_objfun)
            case(OBJFUN_CC)
                THROW_HARD('Not implemented yet')
            case(OBJFUN_EUCLID)
                call self%gen_many_euclids(nr, irefs, iptcl, shift)
        end select
    end subroutine gen_many_objfun_vals

    module subroutine gen_many_euclids( self, nr, irefs, iptcl, shift )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: nr
        integer,                     intent(in)    :: irefs(nr)
        integer,                     intent(in)    :: iptcl
        real(sp),                    intent(in)    :: shift(2)
        complex(sp), pointer :: shmat(:,:)
        complex(sp) :: c
        real(dp)    :: A, v
        real(sp)    :: wk, shift_mag_sq
        integer     :: k, kk, k0, i, j, ithr, ind, iref, p
        logical     :: even
        ithr  = omp_get_thread_num() + 1
        i     = self%pinds(iptcl)
        even  = self%iseven(i)
        k0    = self%kfromto(1)
        shmat => self%heap_vars(ithr)%shmat
        ! zero output
        self%crmat_many(ithr)%r(:,:) = ZERO
        ! Main branch: shift?
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        if ( shift_mag_sq > SHERRSQ ) then
            ! Obtain shift matrix to apply all references
            call self%gen_shmat4aln(ithr, shift, shmat)
            ! Reference loop
            do j = 1, nr
                iref = irefs(j)
                ! S.REF
                if (even) then
                    do k = self%kfromto(1), self%kfromto(2)
                        kk = k - k0 + 1
                        self%cmat2_many(ithr)%c(1:self%pftsz,            kk) =       shmat(:,k) * self%pfts_refs_even(:,k,iref)
                        self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots, kk) = conjg(shmat(:,k) * self%pfts_refs_even(:,k,iref))
                    end do
                else
                    do k = self%kfromto(1), self%kfromto(2)
                        kk = k - k0 + 1
                        self%cmat2_many(ithr)%c(1:self%pftsz,            kk) =       shmat(:,k) * self%pfts_refs_odd(:,k,iref)
                        self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots, kk) = conjg(shmat(:,k) * self%pfts_refs_odd(:,k,iref))
                    end do
                endif
                ! FT(S.REF)
                call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
                ! sum_k w_k * (FT(CTF2) x FT(REF2) - 2*FT(X.CTF) x FT(S.REF)*)
                if( even )then
                    do k = self%kfromto(1), self%kfromto(2)
                        wk = real(k) / self%sigma2_noise(k,iptcl)
                        kk = k - k0 + 1
                        do p = 1, self%pftsz+1
                            c = self%ft_ctf2(p,k,i) * self%ft_ref2_even(p,k,iref)
                            c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%cmat2_many(ithr)%c(p, kk))
                            self%crmat_many(ithr)%c(p,j) = self%crmat_many(ithr)%c(p,j) + wk * c
                        enddo
                    end do
                else
                    do k = self%kfromto(1), self%kfromto(2)
                        wk = real(k) / self%sigma2_noise(k,iptcl)
                        kk = k - k0 + 1
                        do p = 1, self%pftsz+1
                            c = self%ft_ctf2(p,k,i) * self%ft_ref2_odd(p,k,iref)
                            c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%cmat2_many(ithr)%c(p, kk))
                            self%crmat_many(ithr)%c(p,j) = self%crmat_many(ithr)%c(p,j) + wk * c
                        enddo
                    end do
                endif
            enddo
        else
            ! No shift: use memoize FT(REF)
            ! sum_k w_k * (FT(CTF2) x FT(REF2) - 2*FT(X.CTF) x FT(REF)*)
            if( even )then
                do j = 1, nr
                    iref = irefs(j)
                    do k = self%kfromto(1), self%kfromto(2)
                        wk = real(k) / self%sigma2_noise(k,iptcl)
                        do p = 1, self%pftsz+1
                            c = self%ft_ctf2(p,k,i) * self%ft_ref2_even(p,k,iref)
                            c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%ft_ref_even(p,k,iref))
                            self%crmat_many(ithr)%c(p,j) = self%crmat_many(ithr)%c(p,j) + wk * c
                        enddo
                    end do
                enddo
            else
                do j = 1, nr
                    iref = irefs(j)
                    do k = self%kfromto(1), self%kfromto(2)
                        wk = real(k) / self%sigma2_noise(k,iptcl)
                        do p = 1, self%pftsz+1
                            c = self%ft_ctf2(p,k,i) * self%ft_ref2_odd(p,k,iref)
                            c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%ft_ref_odd(p,k,iref))
                            self%crmat_many(ithr)%c(p,j) = self%crmat_many(ithr)%c(p,j) + wk * c
                        enddo
                    end do
                enddo
            endif
        endif
        ! iFFT for all references: IFFT(sum_k w_k*(FT(CTF2) x FT(REF2) - 2*FT(X.CTF) x FT(S.REF)*))
        call fftwf_execute_dft_c2r(self%plan_bwd_many_refs, &
                                    &self%crmat_many(ithr)%c, self%crmat_many(ithr)%r)
        ! normalization & exponentiation
        A = self%wsqsums_ptcls(i) * real(2*self%nrots, dp)
        do j = 1, nr
            do p = 1, self%nrots
                v = 1.d0 + real(self%crmat_many(ithr)%r(p,j), dp) / A
                self%crmat_many(ithr)%r(p,j) = real(exp(-v),sp)
            enddo
        enddo
    end subroutine gen_many_euclids

    ! Benchmark for comparison of euclid calculation between
    ! openmp offload gpu and cpu host
    module subroutine gen_many_euclids_gpu( self, iptcl )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iptcl
        complex(sp), allocatable :: pft_ptcl(:,:)
        real(sp),    allocatable :: ref_vals(:,:), absctf(:,:)
        complex(sp) :: crefctf, cdiff
        real        :: ftvals(self%nrots), w(self%kfromto(1):self%kfromto(2))
        real        :: sqsumptcl, acc, acc_c
        integer(dp) :: t
        integer     :: i, k, rot, rot_c, ir, jr, ref
#ifdef USE_OPENMP_OFFLOAD
        ! particle index
        i = self%pinds(iptcl)
        ! particle variance and shell weights
        sqsumptcl = 0.0
        do k = self%kfromto(1),self%kfromto(2)
            w(k)      = real(k) / self%sigma2_noise(k,iptcl)
            sqsumptcl = sqsumptcl + w(k) * sum(real(self%pfts_ptcls(:,k,i)*conjg(self%pfts_ptcls(:,k,i))))
        enddo
        ! convenience allocation
        allocate(ref_vals(self%nrots,self%nrefs),source=0.0)
        allocate(pft_ptcl(self%kfromto(1):self%kfromto(2), 1:self%pftsz),&
            &absctf(self%kfromto(1):self%kfromto(2), 1:self%pftsz))
        ! tranposed particle and corresponding CTF prior to device upload
        pft_ptcl(self%kfromto(1):self%kfromto(2), 1:self%pftsz) = transpose(self%pfts_ptcls(:,self%kfromto(1):self%kfromto(2),i))
        absctf(self%kfromto(1):self%kfromto(2),   1:self%pftsz) = transpose(self%ctfmats(:,self%kfromto(1):self%kfromto(2),i))
        t = tic()
        !$omp target teams distribute&
        !$omp& map(to:self%kfromto,self%pftsz,self%nrots,self%nrefs,sqsumptcl,w)&
        !$omp& map(to:self%pfts_refs_even(:,self%kfromto(1):self%kfromto(2),:),pft_ptcl,absctf)&
        !$omp& map(tofrom:ref_vals)&
        !$omp& num_teams(self%nrefs)
        do ref = 1,self%nrefs               ! reference
            !$omp parallel do collapse(2) schedule(static)&
            !$omp& private(ir,jr,rot,rot_c,acc,acc_c,k,crefctf,cdiff)
            do ir = 1,self%pftsz            ! reference rotation
                do jr = 1,self%pftsz        ! particle rotation
                    rot   = modulo(jr - ir, self%nrots) + 1
                    rot_c = merge(rot + self%pftsz, rot - self%pftsz, rot <= self%pftsz)
                    acc   = 0.
                    acc_c = 0.
                    do k = self%kfromto(1), self%kfromto(2)
                        crefctf = self%pfts_refs_even(ir, k, ref) * absctf(k, jr)
                        cdiff   = crefctf - pft_ptcl(k, jr)
                        acc     = acc   + w(k) * real(cdiff * conjg(cdiff))
                        cdiff   = conjg(crefctf) - pft_ptcl(k, jr)
                        acc_c   = acc_c + w(k) * real(cdiff * conjg(cdiff))
                    end do
                    !$omp atomic
                    ref_vals(rot,   ref) = ref_vals(rot,  ref)  + acc
                    !$omp atomic
                    ref_vals(rot_c, ref) = ref_vals(rot_c, ref) + acc_c
                enddo
            enddo
            !$omp parallel do private(rot)
            do rot = 1, self%nrots
                ref_vals(rot, ref) = exp(-ref_vals(rot, ref) / real(sqsumptcl))
            enddo
        enddo
        !$omp end target teams distribute
        print *,'GPU: ',toc(t)
        t = tic()
        do ref = self%nrefs, 1, -1
            call self%gen_euclids(ref, iptcl, [0.,0.], ftvals)
        end do
        print *,'CPU: ',toc(t)
        print *,ref_vals(1:10,1)
        print *,ftvals(1:10)
        print *,ref_vals(self%nrots-10:,1)
        print *,ftvals(self%nrots-10:)
#endif
    end subroutine gen_many_euclids_gpu


    ! Benchmark for comparison of euclid calculation between
    ! openmp offload gpu and cpu host
    module subroutine gen_many_euclids_cufft( self, iptcl )
        use  simple_gpu_utils
        use, intrinsic :: iso_c_binding, only: c_loc, c_f_pointer
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iptcl
        type(c_ptr) :: p_c
        complex(sp) :: c,d
        real(sp)    :: ftvals(self%nrots), ftvals_ref1(self%nrots)
        real(sp)    :: w(self%kfromto(1):self%kfromto(2)), wk
        integer(dp) :: t
        integer     :: i, k, iref, p,j
        integer     :: worst_ref, worst_rot
        real(dp)    :: A, df, max_abs, max_rel, rms, sumsq, f_gpu, f_cpu
        logical     :: even
        integer(c_int) :: plan_cufft, ierr
        integer(c_int) :: n(1), inembed(1), onembed(1)
        integer(c_int) :: howmany, idist, odist
        complex(sp), allocatable, target :: buffer_c(:,:)
        real(sp),                pointer :: buffer(:,:)
#ifdef USE_OPENMP_OFFLOAD
        howmany    = int(self%nrefs, c_int)
        idist      = int(self%pftsz+1, c_int)
        odist      = int(self%nrots+2, c_int)
        allocate(buffer_c(idist, howmany), source=cmplx(0.,0.,sp))
        p_c = c_loc(buffer_c(1,1))
        call c_f_pointer(p_c, buffer, [odist, int(howmany)])
        n(1)       = int(self%nrots, c_int)
        inembed(1) = idist
        onembed(1) = odist
        ierr = cufftPlanMany(plan_cufft, 1_c_int, n, inembed, 1_c_int, idist, onembed, 1_c_int, odist, CUFFT_C2R, howmany)
        if (ierr /= 0) stop 68

        ! particle index
        i    = self%pinds(iptcl)
        even = self%iseven(i)
        ! particle variance and shell weights
        do k = self%kfromto(1),self%kfromto(2)
            w(k) = real(k) / real(self%sigma2_noise(k,iptcl))
        enddo
        A = self%wsqsums_ptcls(i) * real(2*self%nrots, dp)

        t = tic()
        !$omp target data map(from:buffer_c)
        !$omp target teams distribute parallel do collapse(2)&
        !$omp& map(to:self%kfromto,self%pftsz,self%nrefs,w,i)&
        !$omp& map(to:self%ft_ctf2(:,:,i),self%ft_ref2_even)&
        !$omp& map(to:self%ft_ptcl_ctf(:,:,i),self%ft_ref_even)&
        !$omp& private(k,c,d)
        do iref = 1, self%nrefs
            do p = 1, self%pftsz
                d = cmplx(0.,0.,sp)
                do k = self%kfromto(1), self%kfromto(2)
                    c = self%ft_ctf2(p,k,i) * self%ft_ref2_even(p,k,iref)
                    c = c - 2. * self%ft_ptcl_ctf(p,k,i) * conjg(self%ft_ref_even(p,k,iref))
                    d = d + w(k) * c
                end do
                buffer_c(p, iref) = d
            enddo
        enddo
        !$omp end target teams distribute parallel do
        p_c  = omp_get_mapped_ptr(c_loc(buffer_c(1,1)), 0)
        ierr = cufftExecC2R(plan_cufft, p_c, p_c)
        !$omp end target data

        ! apply normalization & exponentiation on CPU (buffer_c is now back on host)
        do iref = 1,self%nrefs
            do p = 1, self%nrots
                buffer(p,iref) = real(exp(-(1.d0 + real(buffer(p,iref), dp) / A)), sp)
            enddo
        enddo
        print *,'GPU: ',toc(t)
        t = tic()
        do j = i,i
            max_abs   = 0.d0
            max_rel   = 0.d0
            sumsq     = 0.d0
            worst_ref = 1
            worst_rot = 1
            ftvals_ref1 = 0.0
            do iref = 1, self%nrefs
                call self%gen_euclids(iref, iptcl, [0.,0.], ftvals)
                if (iref == 1) ftvals_ref1 = ftvals
                do p = 1, self%nrots
                    f_gpu = real(buffer(p,iref), dp)
                    f_cpu = real(ftvals(p), dp)
                    df = abs(f_gpu - f_cpu)
                    sumsq = sumsq + df * df
                    if (df > max_abs) then
                        max_abs = df
                        worst_ref = iref
                        worst_rot = p
                    endif
                    max_rel = max(max_rel, df / max(abs(f_cpu), 1.d-14))
                enddo
            end do
            rms = sqrt(sumsq / real(self%nrefs * self%nrots, dp))
        enddo
        print *,'CPU: ',toc(t)
        print '(a,3es16.7)', 'CUFFT/gen_euclids max_abs/max_rel/rms=', max_abs, max_rel, rms
        print '(a,2i8)', 'CUFFT/gen_euclids worst ref/rot=', worst_ref, worst_rot
        print *,buffer(1:10,1)
        print *,ftvals_ref1(1:10)
        print *,buffer(self%nrots-10:,1)
        print *,ftvals_ref1(self%nrots-10:)
        ierr = cufftDestroy(plan_cufft)
        if (ierr /= 0) stop 71
        deallocate(buffer_c)
#endif
    end subroutine gen_many_euclids_cufft

    ! Benchmark for comparison of euclid calculation between
    ! openmp offload gpu and cpu host
    module subroutine gen_many2many_euclids_cufft( self, pfromto )
        use  simple_gpu_utils
        use, intrinsic :: iso_c_binding, only: c_loc, c_f_pointer
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: pfromto(2)
        type(c_ptr) :: p_c
        complex(sp), allocatable, target :: buffer_c(:,:)
        real(sp),                pointer :: buffer(:,:)
        real(sp),    allocatable         :: ftvals(:,:),wks(:,:)
        complex(sp)    :: c,d
        integer(dp)    :: t
        integer        :: i, k, iref, p,j, iptcl, nptcls, is,ie
        integer        :: worst_ref, worst_rot
        real(dp)       :: A, df, max_abs, max_rel, rms, sumsq, f_gpu, f_cpu
        logical        :: even
        integer(c_int) :: plan_cufft, ierr
        integer(c_int) :: n(1), inembed(1), onembed(1), howmany, idist, odist
        nptcls = pfromto(2) - pfromto(1) + 1
#ifdef USE_OPENMP_OFFLOAD
        howmany    = int(self%nrefs*nptcls, c_int)
        idist      = int(self%pftsz+1, c_int)
        odist      = int(self%nrots+2, c_int)
        allocate(buffer_c(idist, howmany), source=cmplx(0.,0.,sp))
        p_c = c_loc(buffer_c(1,1))
        call c_f_pointer(p_c, buffer, [odist, int(howmany)])
        n(1)       = int(self%nrots, c_int)
        inembed(1) = idist
        onembed(1) = odist
        ierr = cufftPlanMany(plan_cufft, 1_c_int, n, inembed, 1_c_int, idist, onembed,&
                            &1_c_int, odist, CUFFT_C2R, howmany)
        if (ierr /= 0) stop 6
        allocate(wks(self%kfromto(1):self%kfromto(2),pfromto(1):pfromto(2)), source=0.)
        do iptcl = pfromto(1),pfromto(2)
            do k = self%kfromto(1),self%kfromto(2)
                wks(k,iptcl) = real(k) / real(self%sigma2_noise(k,iptcl))
            enddo
        enddo
        is = self%pinds(pfromto(1))
        ie = self%pinds(pfromto(2))
        t = tic()
        !$omp target data map(from:buffer_c)
        !$omp target teams distribute parallel do collapse(2)&
        !$omp& map(to:self%kfromto, self%pftsz, self%nrefs, wks)&
        !$omp& map(to:self%ft_ref2_even, self%ft_ref_even)&
        !$omp& map(to:self%wsqsums_ptcls(is:ie), self%pinds(is:ie))&
        !$omp& map(to:self%ft_ctf2(:,:,is:ie),pfromto)&
        !$omp& map(to:self%ft_ptcl_ctf(:,:,is:ie))&
        !$omp& private(i,k,c,d,iref,iptcl,p,j)
        do iptcl = pfromto(1), pfromto(2)
            do iref = 1, self%nrefs
                i = self%pinds(iptcl)
                j = (iptcl - pfromto(1))*self%nrefs + iref
                do p = 1, self%pftsz
                    d = cmplx(0.,0.,sp)
                    do k = self%kfromto(1), self%kfromto(2)
                        c = self%ft_ctf2(p,k,i) * self%ft_ref2_even(p,k,iref)
                        c = c - 2.0_sp * self%ft_ptcl_ctf(p,k,i) * conjg(self%ft_ref_even(p,k,iref))
                        d = d + wks(k,iptcl) * c
                    end do
                    buffer_c(p, j) = d
                enddo
            enddo
        enddo
        p_c  = omp_get_mapped_ptr(c_loc(buffer_c(1,1)), 0)
        ierr = cufftExecC2R(plan_cufft, p_c, p_c)
        if (ierr /= 0) stop 7
        !$omp end target data
        ! apply normalization & exponentiation on CPU
        do iptcl = pfromto(1), pfromto(2)
            i = self%pinds(iptcl)
            A = self%wsqsums_ptcls(i) * real(2*self%nrots, dp)
            do iref = 1,self%nrefs
                j = (iptcl - pfromto(1))*self%nrefs + iref
                do p = 1, self%nrots
                    buffer(p,j) = real(exp(-(1.d0 + real(buffer(p,j), dp) / A)), sp)
                enddo
            enddo
        enddo
        print *,'GPU: ',toc(t)
        ierr = cufftDestroy(plan_cufft)
        if (ierr /= 0) stop 9
        allocate(ftvals(self%nrots,howmany))
        t = tic()
        do iptcl = pfromto(1), pfromto(2)
            do iref = 1, self%nrefs
                j = (iptcl - pfromto(1))*self%nrefs + iref
                call self%gen_euclids(iref, iptcl, [0.,0.], ftvals(:,j))
            enddo
        enddo
        max_abs   = 0.d0
        max_rel   = 0.d0
        sumsq     = 0.d0
        worst_ref = 1
        worst_rot = 1
        do j = 1, howmany
            do p = 1, self%nrots
                f_gpu = real(buffer(p,j), dp)
                f_cpu = real(ftvals(p,j), dp)
                df = abs(f_gpu - f_cpu)
                sumsq = sumsq + df * df
                if (df > max_abs) then
                    max_abs = df
                    worst_ref = j
                    worst_rot = p
                endif
                max_rel = max(max_rel, df / max(abs(f_cpu), 1.d-14))
            enddo
        enddo
        rms = sqrt(sumsq / real(howmany * self%nrots, dp))
        print *,'CPU: ',toc(t)
        print '(a,3es16.7)', 'CUFFT/gen_euclids max_abs/max_rel/rms=', max_abs, max_rel, rms
        print '(a,2i8)', 'CUFFT/gen_euclids worst slot/rot=', worst_ref, worst_rot
        deallocate(buffer_c, ftvals, wks)
#endif
    end subroutine gen_many2many_euclids_cufft

end submodule simple_polarft_corr
