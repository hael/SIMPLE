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
            ! The reference is not shifted, memoized FT(S.REF) can be used
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

end submodule simple_polarft_corr
