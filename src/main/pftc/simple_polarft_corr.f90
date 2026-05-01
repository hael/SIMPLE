!@descr: polarft class submodule for objective function evaluations
submodule (simple_polarft_calc) simple_polarft_corr
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
        call self%gen_shmat_8(ithr, real(shvec,dp),shmat)
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
        call self%gen_shmat_8(ithr, real(shvec,dp), shmat)
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
            call self%gen_shmat(ithr, shift, shmat)
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
        complex(sp), pointer :: pft_ref(:,:), shmat(:,:)
        real(dp), pointer    :: w_weights(:)
        real(dp) :: ptcl_sqsum
        real(sp) :: shift_mag_sq
        integer  :: k, i, ithr, kk, k0
        logical  :: even, needs_shift
        ithr         =  omp_get_thread_num() + 1
        i            =  self%pinds(iptcl)
        k0           =  self%kfromto(1)
        even         =  self%iseven(i)
        ptcl_sqsum   =  self%wsqsums_ptcls(i)
        shmat        => self%heap_vars(ithr)%shmat
        w_weights    => self%heap_vars(ithr)%w_weights
        ! ========================================================================
        ! Prepare FFT(S.REF) if needed
        ! ========================================================================
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        needs_shift  = shift_mag_sq > SHERRSQ
        if (needs_shift) then
            pft_ref => self%heap_vars(ithr)%pft_ref
            ! Generate shift matrix
            call self%gen_shmat(ithr, shift, shmat)
            ! Shift reference
            if (even) then
                pft_ref = shmat(:,:self%kfromto(2)) * self%pfts_refs_even(:,self%kfromto(1):self%kfromto(2),iref)
            else
                pft_ref = shmat(:,:self%kfromto(2)) * self%pfts_refs_odd(:,self%kfromto(1):self%kfromto(2),iref)
            endif
            ! Calculate FFT(S.REF)
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%cmat2_many(ithr)%c(1:self%pftsz,            kk) =       pft_ref(:,k)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots, kk) = conjg(pft_ref(:,k))
            end do
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
        else
            ! no shift, memoized reference can be used
            if( even )then
                self%cmat2_many(ithr)%c(1:self%pftsz+1,1:self%nk) = self%ft_ref_even(:,self%kfromto(1):self%kfromto(2),iref)
            else
                self%cmat2_many(ithr)%c(1:self%pftsz+1,1:self%nk) = self%ft_ref_odd( :,self%kfromto(1):self%kfromto(2),iref)
            endif
        endif
        ! ========================================================================
        ! Pre-compute weights and particle sums ONCE
        ! This eliminates redundant computation in the second k-loop
        ! ========================================================================
        do k = self%kfromto(1), self%kfromto(2)
            kk = k - k0 + 1
            w_weights(kk)   = real(k, dp) / real(self%sigma2_noise(k,iptcl), dp)
        end do
        ! ========================================================================
        ! Single IFFT: IFFT( sum_k w_k * (FT(CTF2) x FT(REF2) - 2*FT(X.CTF) x FT(S.REF)*) )
        ! ========================================================================
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        if (even) then
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%crvec1(ithr)%c = self%crvec1(ithr)%c + real(w_weights(kk),c_float) * ( &
                    self%ft_ctf2(:,k,i) * self%ft_ref2_even(:,k,iref) - &
                    2.0 * self%ft_ptcl_ctf(:,k,i) * &
                    conjg(self%cmat2_many(ithr)%c(1:self%pftsz+1, kk)) )
            end do
        else
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%crvec1(ithr)%c = self%crvec1(ithr)%c + real(w_weights(kk),c_float) * ( &
                    self%ft_ctf2(:,k,i) * self%ft_ref2_odd(:,k,iref) - &
                    2.0 * self%ft_ptcl_ctf(:,k,i) * &
                    conjg(self%cmat2_many(ithr)%c(1:self%pftsz+1, kk)) )
            end do
        endif
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        euclids = real(exp( -1.d0 -&
            real(self%crvec1(ithr)%r(1:self%nrots),dp) / (ptcl_sqsum * real(2*self%nrots,dp)) ),dp)
    end subroutine gen_euclids

    module real(dp) function gen_corr_for_rot_8_1( self, iref, iptcl, irot )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl, irot
        complex(dp), pointer :: pft_ref_8(:,:), pft_ref_tmp_8(:,:)
        integer :: ithr, i
        i    = self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
        if (self%iseven(i)) then
            pft_ref_8 = self%pfts_refs_even(:,self%kfromto(1):self%kfromto(2),iref)
        else
            pft_ref_8 = self%pfts_refs_odd(:,self%kfromto(1):self%kfromto(2),iref)
        endif
        ! rotation
        call self%rotate_ref_8(pft_ref_8, irot, pft_ref_tmp_8)
        ! ctf
        pft_ref_tmp_8 = pft_ref_tmp_8 * self%ctfmats(:,:self%kfromto(2),i)
        gen_corr_for_rot_8_1 = 0.d0
        select case(self%p_ptr%cc_objfun)
            case(OBJFUN_CC)
                gen_corr_for_rot_8_1 = self%gen_corr_cc_for_rot_8(pft_ref_tmp_8, i)
            case(OBJFUN_EUCLID)
                gen_corr_for_rot_8_1 = self%gen_euclid_for_rot_8(pft_ref_tmp_8, iptcl)
        end select
    end function gen_corr_for_rot_8_1

    module real(dp) function gen_corr_for_rot_8_2( self, iref, iptcl, shvec, irot )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(dp),                    intent(in)    :: shvec(2)
        integer,                     intent(in)    :: irot
        complex(dp), pointer :: pft_ref_8(:,:), pft_ref_tmp_8(:,:), shmat_8(:,:)
        integer :: ithr, i
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
        call self%gen_shmat_8(ithr, shvec, shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8(:,:self%kfromto(2))
        ! rotation
        call self%rotate_ref_8(pft_ref_8, irot, pft_ref_tmp_8)
        ! ctf
        pft_ref_tmp_8 = pft_ref_tmp_8 * self%ctfmats(:,:self%kfromto(2),i)
        gen_corr_for_rot_8_2 = 0.d0
        select case(self%p_ptr%cc_objfun)
            case(OBJFUN_CC)
                gen_corr_for_rot_8_2 = self%gen_corr_cc_for_rot_8(pft_ref_tmp_8, i)
            case(OBJFUN_EUCLID)
                gen_corr_for_rot_8_2 = self%gen_euclid_for_rot_8(pft_ref_tmp_8, iptcl)
        end select
    end function gen_corr_for_rot_8_2

    module real(dp) function gen_corr_cc_for_rot_8( self, pft_ref, i )
        class(polarft_calc),  intent(inout) :: self
        complex(dp), pointer, intent(inout) :: pft_ref(:,:)
        integer,              intent(in)    :: i
        real(dp) :: sqsum_ref, k_weight
        integer  :: k
        sqsum_ref             = 0.d0
        gen_corr_cc_for_rot_8 = 0.d0
        ! Single k-loop with cached weight computation
        do k = self%kfromto(1), self%kfromto(2)
            k_weight = real(k, kind=dp)
            sqsum_ref             = sqsum_ref + k_weight * &
                sum(real(pft_ref(:,k) * conjg(pft_ref(:,k)), dp))
            gen_corr_cc_for_rot_8 = gen_corr_cc_for_rot_8 + k_weight * &
                sum(real(pft_ref(:,k) * conjg(self%pfts_ptcls(:,k,i)), dp))
        end do
        gen_corr_cc_for_rot_8 = gen_corr_cc_for_rot_8 / dsqrt(sqsum_ref * self%ksqsums_ptcls(i))
    end function gen_corr_cc_for_rot_8

    module real(dp) function gen_euclid_for_rot_8( self, pft_ref, iptcl )
        class(polarft_calc),  intent(inout) :: self
        complex(dp), pointer, intent(inout) :: pft_ref(:,:)
        integer,              intent(in)    :: iptcl
        real(dp) :: w, sum_term
        integer  :: i, k
        i       = self%pinds(iptcl)
        pft_ref = pft_ref - self%pfts_ptcls(:,:self%kfromto(2),i)
        gen_euclid_for_rot_8 = 0.d0
        do k = self%kfromto(1), self%kfromto(2)
            w = real(k, dp) / real(self%sigma2_noise(k,iptcl), dp)
            sum_term = sum(real(pft_ref(:,k) * conjg(pft_ref(:,k)), dp))
            gen_euclid_for_rot_8 = gen_euclid_for_rot_8 + w * sum_term
        end do
        gen_euclid_for_rot_8 = dexp(-gen_euclid_for_rot_8 / self%wsqsums_ptcls(i))
    end function gen_euclid_for_rot_8

    module subroutine gen_corr_grad_for_rot_8( self, iref, iptcl, shvec, irot, f, grad )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(dp),                    intent(in)    :: shvec(2)
        integer,                     intent(in)    :: irot
        real(dp),                    intent(out)   :: f, grad(2)
        complex(dp), pointer :: pft_ref_8(:,:), shmat_8(:,:), pft_ref_tmp_8(:,:)
        integer :: ithr, i
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
        call self%gen_shmat_8(ithr, shvec, shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8(:,:self%kfromto(2))
        select case(self%p_ptr%cc_objfun)
            case(OBJFUN_CC)
                call self%gen_corr_cc_grad_for_rot_8(pft_ref_8, pft_ref_tmp_8, iptcl, irot, f, grad)
            case(OBJFUN_EUCLID)
                call self%gen_euclid_grad_for_rot_8(pft_ref_8, pft_ref_tmp_8, iptcl, irot, f, grad)
        end select
    end subroutine gen_corr_grad_for_rot_8

    module subroutine gen_corr_cc_grad_for_rot_8( self, pft_ref, pft_ref_tmp, iptcl, irot, f, grad )
        class(polarft_calc),  intent(inout) :: self
        complex(dp), pointer, intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,              intent(in)    :: iptcl, irot
        real(dp),             intent(out)   :: f, grad(2)
        real(dp) :: sqsum_ref, sqsum_ptcl, denom, k_weight
        integer  :: k, i
        i          = self%pinds(iptcl)
        sqsum_ref  = 0.d0
        f          = 0.d0
        grad       = 0.d0
        sqsum_ptcl = self%ksqsums_ptcls(i)
        ! First rotation and k-loop: compute sqsum_ref and f
        call self%rotate_ref_8(pft_ref, irot, pft_ref_tmp)
        do k = self%kfromto(1), self%kfromto(2)
            k_weight  = real(k, kind=dp)
            sqsum_ref = sqsum_ref + k_weight * sum(real(self%ctfmats(:,k,i) * self%ctfmats(:,k,i) * &
                                                        pft_ref_tmp(:,k) * conjg(pft_ref_tmp(:,k)), dp))
            f         = f + k_weight * sum(real(self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * &
                                                conjg(self%pfts_ptcls(:,k,i)), dp))
        enddo
        ! Second rotation: gradient x-component
        call self%rotate_ref_8(pft_ref * dcmplx(0.d0, self%argtransf(:self%pftsz,:self%kfromto(2))), irot, pft_ref_tmp)
        do k = self%kfromto(1), self%kfromto(2)
            k_weight = real(k, kind=dp)
            grad(1)  = grad(1) + k_weight * sum(real(self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * &
                                                    conjg(self%pfts_ptcls(:,k,i)), dp))
        enddo
        ! Third rotation: gradient y-component
        call self%rotate_ref_8(pft_ref * dcmplx(0.d0, self%argtransf(self%pftsz+1:,:self%kfromto(2))), irot, pft_ref_tmp)
        do k = self%kfromto(1), self%kfromto(2)
            k_weight = real(k, kind=dp)
            grad(2)  = grad(2) + k_weight * sum(real(self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * &
                                                    conjg(self%pfts_ptcls(:,k,i)), dp))
        end do
        ! Final normalization
        denom = dsqrt(sqsum_ref * sqsum_ptcl)
        f     = f / denom
        grad  = grad / denom
    end subroutine gen_corr_cc_grad_for_rot_8

    module subroutine gen_euclid_grad_for_rot_8( self, pft_ref, pft_ref_tmp, iptcl, irot, f, grad )
        class(polarft_calc),  target, intent(inout) :: self
        complex(dp), pointer,         intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,                      intent(in)    :: iptcl, irot
        real(dp),                     intent(out)   :: f, grad(2)
        complex(dp), pointer :: pft_diff(:,:)
        real(dp),    pointer :: w_weights(:)
        real(dp) :: denom, w, sum_diff
        integer  :: k, i, ithr, kk, k0
        ithr     = omp_get_thread_num() + 1
        i        = self%pinds(iptcl)
        k0       = self%kfromto(1)
        f        = 0.d0
        grad     = 0.d0
        denom    = self%wsqsums_ptcls(i)
        pft_diff => self%heap_vars(ithr)%pft_ref_tmp2_8
        w_weights => self%heap_vars(ithr)%w_weights
        do k = self%kfromto(1), self%kfromto(2)
            kk = k - k0 + 1
            w_weights(kk) = real(k, dp) / real(self%sigma2_noise(k,iptcl), dp)
        end do
        ! First rotation: compute difference and f
        call self%rotate_ref_8(pft_ref, irot, pft_ref_tmp)
        pft_ref_tmp = pft_ref_tmp * self%ctfmats(:,:self%kfromto(2),i)
        pft_diff    = pft_ref_tmp - self%pfts_ptcls(:,:self%kfromto(2),i)  ! Ref(shift + rotation + CTF) - Ptcl
        do k = self%kfromto(1), self%kfromto(2)
            kk       = k - k0 + 1
            w        = w_weights(kk)
            sum_diff = sum(real(pft_diff(:,k) * conjg(pft_diff(:,k)), dp))
            f        = f + w * sum_diff
        end do
        ! Second rotation: x-gradient
        call self%rotate_ref_8(pft_ref * dcmplx(0.d0, self%argtransf(:self%pftsz,:self%kfromto(2))), irot, pft_ref_tmp)
        pft_ref_tmp = pft_ref_tmp * self%ctfmats(:,:self%kfromto(2),i)
        do k = self%kfromto(1), self%kfromto(2)
            kk      = k - k0 + 1
            w       = w_weights(kk)
            grad(1) = grad(1) + w * real(sum(pft_ref_tmp(:,k) * conjg(pft_diff(:,k))), dp)
        end do
        ! Third rotation: y-gradient
        call self%rotate_ref_8(pft_ref * dcmplx(0.d0, self%argtransf(self%pftsz+1:,:self%kfromto(2))), irot, pft_ref_tmp)
        pft_ref_tmp = pft_ref_tmp * self%ctfmats(:,:self%kfromto(2),i)
        do k = self%kfromto(1), self%kfromto(2)
            kk      = k - k0 + 1
            w       = w_weights(kk)
            grad(2) = grad(2) + w * real(sum(pft_ref_tmp(:,k) * conjg(pft_diff(:,k))), dp)
        end do
        ! Final computation
        f    = dexp(-f / denom)
        grad = -f * 2.d0 * grad / denom
    end subroutine gen_euclid_grad_for_rot_8

    module subroutine gen_corr_grad_only_for_rot_8( self, iref, iptcl, shvec, irot, grad )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(dp),                    intent(in)    :: shvec(2)
        integer,                     intent(in)    :: irot
        real(dp),                    intent(out)   :: grad(2)
        complex(dp), pointer :: pft_ref_8(:,:), shmat_8(:,:), pft_ref_tmp_8(:,:)
        real(dp) :: f
        integer  :: ithr, i
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
        call self%gen_shmat_8(ithr, shvec, shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8(:,:self%kfromto(2))
        select case(self%p_ptr%cc_objfun)
            case(OBJFUN_CC)
                call self%gen_corr_cc_grad_only_for_rot_8(pft_ref_8, pft_ref_tmp_8, i, irot, grad)
            case(OBJFUN_EUCLID)
                call self%gen_euclid_grad_for_rot_8(pft_ref_8, pft_ref_tmp_8, iptcl, irot, f, grad)
        end select
    end subroutine gen_corr_grad_only_for_rot_8

    module subroutine gen_corr_cc_grad_only_for_rot_8( self, pft_ref, pft_ref_tmp, i, irot, grad )
        class(polarft_calc),  intent(inout) :: self
        complex(dp), pointer, intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,              intent(in)    :: i, irot
        real(dp),             intent(out)   :: grad(2)
        real(dp) :: sqsum_ref, sqsum_ptcl, k_weight
        integer  :: k
        sqsum_ref  = 0.d0
        grad       = 0.d0
        sqsum_ptcl = self%ksqsums_ptcls(i)
        ! First rotation: compute sqsum_ref
        call self%rotate_ref_8(pft_ref, irot, pft_ref_tmp)
        do k = self%kfromto(1), self%kfromto(2)
            k_weight  = real(k, kind=dp)
            sqsum_ref = sqsum_ref + k_weight * sum(real(self%ctfmats(:,k,i) * self%ctfmats(:,k,i) * &
                                                        pft_ref_tmp(:,k) * conjg(pft_ref_tmp(:,k)), dp))
        enddo
        ! Second rotation: x-gradient
        call self%rotate_ref_8(pft_ref * dcmplx(0.d0, self%argtransf(:self%pftsz,:self%kfromto(2))), irot, pft_ref_tmp)
        do k = self%kfromto(1), self%kfromto(2)
            k_weight = real(k, kind=dp)
            grad(1)  = grad(1) + k_weight * sum(real(self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * &
                                                    conjg(self%pfts_ptcls(:,k,i)), dp))
        enddo
        ! Third rotation: y-gradient
        call self%rotate_ref_8(pft_ref * dcmplx(0.d0, self%argtransf(self%pftsz+1:,:self%kfromto(2))), irot, pft_ref_tmp)
        do k = self%kfromto(1), self%kfromto(2)
            k_weight = real(k, kind=dp)
            grad(2)  = grad(2) + k_weight * sum(real(self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * &
                                                    conjg(self%pfts_ptcls(:,k,i)), dp))
        end do
        grad = grad / dsqrt(sqsum_ref * sqsum_ptcl)
    end subroutine gen_corr_cc_grad_only_for_rot_8

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
        call self%gen_shmat_8(ithr, real(shvec,dp), shmat_8)
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

    ! In developpment

    ! Benchmark for comparison of euclid calculation between
    ! openmp offload gpu and cpu host
    module subroutine gen_all_euclids_gpu( self, iptcl )
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
    end subroutine gen_all_euclids_gpu

    module subroutine gen_all_euclids( self, nr, ref_range, iptcl, euclids )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: nr, ref_range(2), iptcl
        real(sp),                    intent(out)   :: euclids(self%nrots,nr)
        character(kind=c_char, len=:), allocatable :: fft_wisdoms_fname
        type(fftw_crmat)    :: cr_array_many
        type(c_ptr)         :: plan_bwd_many
        real(dp), pointer   :: w_weights(:)
        integer(kind=c_int) :: wsdm_ret, howmany, n(1), inembed(1), onembed(1), istride, ostride, idist, odist
        complex(sp) :: c
        real(dp)    :: ptcl_sqsum
        real(sp)    :: wk
        integer(dp) :: t
        integer     :: k, i, ithr, kk, k0, ind, iref, p
        logical     :: even
        t = tic()
        ithr         =  omp_get_thread_num() + 1
        i            =  self%pinds(iptcl)
        k0           =  self%kfromto(1)
        even         =  self%iseven(i)
        ptcl_sqsum   =  self%wsqsums_ptcls(i)
        w_weights    => self%heap_vars(ithr)%w_weights

        ! Allocate complex storage with fftw
        cr_array_many%p = fftwf_alloc_complex(int((self%pftsz+1) * nr, c_size_t))
        call c_f_pointer(cr_array_many%p, cr_array_many%c, [self%pftsz+1, nr])
        call c_f_pointer(cr_array_many%p, cr_array_many%r, [self%nrots+2, nr])

        ! Planning
        allocate(fft_wisdoms_fname, source='fft_wisdoms.dat'//c_null_char)
        wsdm_ret      = fftw_import_wisdom_from_filename(fft_wisdoms_fname)
        n(1)          = int(self%nrots, c_int)
        howmany       = int(nr, c_int)
        onembed(1)    = n(1)
        odist         = int(self%nrots+2, c_int)
        idist         = int(self%pftsz+1, c_int)
        inembed(1)    = idist
        plan_bwd_many = fftwf_plan_many_dft_c2r( 1_c_int, n, howmany, &
                        &cr_array_many%c, inembed, 1_c_int, idist, &
                        &cr_array_many%r, onembed, 1_c_int, odist, &
                        &ior(FFTW_MEASURE, FFTW_USE_WISDOM) )
        wsdm_ret = fftw_export_wisdom_to_filename(fft_wisdoms_fname)
        deallocate(fft_wisdoms_fname)

        print *, 'time init alloc planning: ', toc(t)
        t = tic()

        ! Pre-compute weights
        do k = self%kfromto(1), self%kfromto(2)
            kk = k - k0 + 1
            w_weights(kk) = real(k, dp) / real(self%sigma2_noise(k,iptcl), dp)
        end do

        cr_array_many%c = CMPLX_ZERO
        if( even )then
            ! Computational note: the triple loop optimizes both for
            ! cache and memory texture and seems the best option
            ind = 0
            do iref = ref_range(1), ref_range(2)
                ind = ind + 1
                do k = self%kfromto(1), self%kfromto(2)
                    kk = k - k0 + 1
                    wk = real(w_weights(kk),sp)
                    do p = 1, self%pftsz+1
                        c = self%ft_ctf2(p,k,i) * self%ft_ref2_even(p,k,iref)
                        c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%ft_ref_even(p,k,iref))
                        cr_array_many%c(p,ind) = cr_array_many%c(p,ind) + wk * c
                    enddo
                end do
            enddo
        else
            ind = 0
            do iref = ref_range(1), ref_range(2)
                ind = ind + 1
                do k = self%kfromto(1), self%kfromto(2)
                    kk = k - k0 + 1
                    wk = real(w_weights(kk),sp)
                    do p = 1, self%pftsz+1
                        c = self%ft_ctf2(p,k,i) * self%ft_ref2_odd(p,k,iref)
                        c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%ft_ref_odd(p,k,iref))
                        cr_array_many%c(p,ind) = cr_array_many%c(p,ind) + wk * c
                    enddo
                end do
            enddo
        endif
        call fftwf_execute_dft_c2r(plan_bwd_many, cr_array_many%c, cr_array_many%r)
        euclids = real(exp(-1.d0 -&
            real(cr_array_many%r(:self%nrots,:), dp) / (ptcl_sqsum * real(2*self%nrots, dp))), sp)

        print *, 'time gen_all_euclids: ', toc(t)
        call fftwf_destroy_plan(plan_bwd_many)
        call fftwf_free(cr_array_many%p)
    end subroutine gen_all_euclids

end submodule simple_polarft_corr
