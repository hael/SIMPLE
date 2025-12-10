submodule (simple_polarft_calc) simple_polarft_corr
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
#include "simple_local_flags.inc"
implicit none

real, parameter :: SHERR = 0.001

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
        i    = self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref     => self%heap_vars(ithr)%pft_ref_8
        pft_rot_ref => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat       => self%heap_vars(ithr)%shmat_8
        if( self%iseven(i) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        call self%gen_shmat_8(ithr, real(shvec,dp),shmat)
        pft_ref = pft_ref * shmat
        call self%rotate_pft(pft_ref, irot, pft_rot_ref)
        pft_rot_ref = pft_rot_ref * self%ctfmats(:,:,i)
        select case(params_glob%cc_objfun)
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
            pft_rot_ref = pft_rot_ref - self%pfts_ptcls(:,:,i)
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
        if( self%iseven(i) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        call self%gen_shmat_8(ithr, real(shvec,dp), shmat)
        pft_ref = pft_ref * shmat
        call self%rotate_pft(pft_ref, irot, pft_rot_ref)
        pft_rot_ref = pft_rot_ref * real(self%ctfmats(:,:,i),dp)
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
         select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                call self%gen_corrs(iref, iptcl, shift, vals)
            case(OBJFUN_EUCLID)
                call self%gen_euclids(iref, iptcl, shift, vals)
        end select
    end subroutine gen_objfun_vals

    module subroutine gen_corrs( self, iref, iptcl, shift, cc )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        real(sp),            intent(in)    :: shift(2)
        real(sp),            intent(out)   :: cc(self%nrots)
        complex(sp), pointer :: pft_ref(:,:), shmat(:,:)
        integer :: i, ithr, k
        ithr    = omp_get_thread_num() + 1
        i       = self%pinds(iptcl)
        shmat   => self%heap_vars(ithr)%shmat
        pft_ref => self%heap_vars(ithr)%pft_ref
        pft_ref = merge(self%pfts_refs_even(:,:,iref), self%pfts_refs_odd(:,:,iref), self%iseven(i))
        if( arg(shift) > SHERR )then
            call self%gen_shmat(ithr, shift, shmat)
            pft_ref = shmat * pft_ref
        endif
        self%heap_vars(ithr)%kcorrs = 0.d0
        self%drvec(ithr)%r          = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            ! FT(CTF2) x FT(REF2)), REF2 is shift invariant
            if( self%iseven(i) )then
                self%cvec1(ithr)%c = self%ft_ctf2(:,k,i) * self%ft_ref2_even(:,k,iref)
            else
                self%cvec1(ithr)%c = self%ft_ctf2(:,k,i) * self%ft_ref2_odd(:,k,iref)
            endif
            ! IFFT(FT(CTF2) x FT(REF2))
            call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
            self%drvec(ithr)%r = self%drvec(ithr)%r + real(self%rvec1(ithr)%r(1:self%nrots),dp)
            ! FT(S.REF), shifted reference
            self%cvec2(ithr)%c(1:self%pftsz)            = pft_ref(:,k)
            self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(pft_ref(:,k))
            call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
            ! FT(X.CTF) x FT(S.REF)*
            self%cvec1(ithr)%c = self%ft_ptcl_ctf(:,k,i) * conjg(self%cvec2(ithr)%c(1:self%pftsz+1))
            ! IFFT(FT(X.CTF) x FT(S.REF)*)
            call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
            self%heap_vars(ithr)%kcorrs = self%heap_vars(ithr)%kcorrs + real(self%rvec1(ithr)%r(1:self%nrots),dp)
        end do
        self%drvec(ithr)%r = self%drvec(ithr)%r * real(self%sqsums_ptcls(i) * real(2*self%nrots),dp)
        cc = real(self%heap_vars(ithr)%kcorrs / dsqrt(self%drvec(ithr)%r))
    end subroutine gen_corrs

    module subroutine gen_euclids( self, iref, iptcl, shift, euclids )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        real(sp),            intent(in)    :: shift(2)
        real(sp),            intent(out)   :: euclids(self%nrots)
        complex(sp), pointer :: pft_ref(:,:), shmat(:,:)
        real(dp) :: w, sumsqptcl
        integer  :: k, i, ithr
        ithr    = omp_get_thread_num() + 1
        i       = self%pinds(iptcl)
        shmat   => self%heap_vars(ithr)%shmat
        pft_ref => self%heap_vars(ithr)%pft_ref
        pft_ref = merge(self%pfts_refs_even(:,:,iref), self%pfts_refs_odd(:,:,iref), self%iseven(i))
        if( arg(shift) > SHERR )then
            call self%gen_shmat(ithr, shift, shmat)
            pft_ref = shmat * pft_ref
        endif
        self%heap_vars(ithr)%kcorrs = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            w         = real(k,dp) / real(self%sigma2_noise(k,iptcl),dp)
            sumsqptcl = sum(real(self%pfts_ptcls(:,k,i)*conjg(self%pfts_ptcls(:,k,i)),dp))
            ! FT(CTF2) x FT(REF2)*
            if( self%iseven(i) )then
                self%cvec1(ithr)%c = self%ft_ctf2(:,k,i) * self%ft_ref2_even(:,k,iref)
            else
                self%cvec1(ithr)%c = self%ft_ctf2(:,k,i) * self%ft_ref2_odd(:,k,iref)
            endif
            ! FT(S.REF), shifted reference
            self%cvec2(ithr)%c(1:self%pftsz)            =       pft_ref(:,k)
            self%cvec2(ithr)%c(self%pftsz+1:self%nrots) = conjg(pft_ref(:,k))
            call fftwf_execute_dft(self%plan_fwd1, self%cvec2(ithr)%c, self%cvec2(ithr)%c)
            ! FT(CTF2) x FT(REF2)* - 2 * FT(X.CTF) x FT(REF)*
            self%cvec1(ithr)%c = self%cvec1(ithr)%c - 2.0 * self%ft_ptcl_ctf(:,k,i) * conjg(self%cvec2(ithr)%c(1:self%pftsz+1))
            ! IFFT( FT(CTF2) x FT(REF2)* - 2 * FT(X.CTF) x FT(REF)* )
            call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
            ! k/sig2 x ( |CTF.REF|2 - 2X.CTF.REF ), fftw normalized
            self%drvec(ithr)%r = (w / real(2*self%nrots,dp)) * real(self%rvec1(ithr)%r(1:self%nrots),dp)
            ! k/sig2 x ( |X|2 + |CTF.REF|2 - 2X.CTF.REF )
            self%heap_vars(ithr)%kcorrs = self%heap_vars(ithr)%kcorrs + w * sumsqptcl + self%drvec(ithr)%r
        end do
        euclids = real( dexp( -self%heap_vars(ithr)%kcorrs / self%wsqsums_ptcls(i) ) )
    end subroutine gen_euclids

    module real(dp) function gen_corr_for_rot_8_1( self, iref, iptcl, irot )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl, irot
        complex(dp), pointer :: pft_ref_8(:,:), pft_ref_tmp_8(:,:)
        integer :: ithr, i
        i    =  self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
        if( self%iseven(i) )then
            pft_ref_8 = self%pfts_refs_even(:,:,iref)
        else
            pft_ref_8 = self%pfts_refs_odd(:,:,iref)
        endif
        ! rotation
        call self%rotate_pft(pft_ref_8, irot, pft_ref_tmp_8)
        ! ctf
        pft_ref_tmp_8 = pft_ref_tmp_8 * self%ctfmats(:,:,i)
        gen_corr_for_rot_8_1 = 0.d0
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                gen_corr_for_rot_8_1 = self%gen_corr_cc_for_rot_8(pft_ref_tmp_8, i)
            case(OBJFUN_EUCLID)
                gen_corr_for_rot_8_1 = self%gen_euclid_for_rot_8(pft_ref_tmp_8, iptcl)
        end select
    end function gen_corr_for_rot_8_1

    module real(dp) function gen_corr_for_rot_8_2( self, iref, iptcl, shvec, irot )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        real(dp),            intent(in)    :: shvec(2)
        integer,             intent(in)    :: irot
        complex(dp), pointer :: pft_ref_8(:,:), pft_ref_tmp_8(:,:), shmat_8(:,:)
        integer :: ithr, i
        i    =  self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat_8       => self%heap_vars(ithr)%shmat_8
        if( self%iseven(i) )then
            pft_ref_8 = self%pfts_refs_even(:,:,iref)
        else
            pft_ref_8 = self%pfts_refs_odd(:,:,iref)
        endif
        ! shift
        call self%gen_shmat_8(ithr, shvec, shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8
        ! rotation
        call self%rotate_pft(pft_ref_8, irot, pft_ref_tmp_8)
        ! ctf
        pft_ref_tmp_8 = pft_ref_tmp_8 * self%ctfmats(:,:,i)
        gen_corr_for_rot_8_2 = 0.d0
        select case(params_glob%cc_objfun)
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
        real(dp) :: sqsum_ref
        integer  :: k
        sqsum_ref            = 0.d0
        gen_corr_cc_for_rot_8 = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            sqsum_ref            = sqsum_ref +            real(k,kind=dp) * sum(real(pft_ref(:,k) * conjg(pft_ref(:,k)),dp))
            gen_corr_cc_for_rot_8 = gen_corr_cc_for_rot_8 + real(k,kind=dp) * sum(real(pft_ref(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
        end do
        gen_corr_cc_for_rot_8 = gen_corr_cc_for_rot_8 / dsqrt(sqsum_ref * self%ksqsums_ptcls(i))
    end function gen_corr_cc_for_rot_8

    module real(dp) function gen_euclid_for_rot_8( self, pft_ref, iptcl )
        class(polarft_calc),  intent(inout) :: self
        complex(dp), pointer, intent(inout) :: pft_ref(:,:)
        integer,              intent(in)    :: iptcl
        integer :: i,k
        i       = self%pinds(iptcl)
        pft_ref = pft_ref - self%pfts_ptcls(:,:,i)
        gen_euclid_for_rot_8 = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            gen_euclid_for_rot_8 = gen_euclid_for_rot_8 +&
                &(real(k,dp) / self%sigma2_noise(k,iptcl)) * sum(real(pft_ref(:,k)*conjg(pft_ref(:,k)),dp))
        end do
        gen_euclid_for_rot_8 = dexp( -gen_euclid_for_rot_8 / self%wsqsums_ptcls(i))
    end function gen_euclid_for_rot_8

    module subroutine gen_corr_grad_for_rot_8( self, iref, iptcl, shvec, irot, f, grad )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        real(dp),            intent(in)    :: shvec(2)
        integer,             intent(in)    :: irot
        real(dp),            intent(out)   :: f, grad(2)
        complex(dp), pointer :: pft_ref_8(:,:), shmat_8(:,:), pft_ref_tmp_8(:,:)
        integer :: ithr, i
        i    = self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat_8       => self%heap_vars(ithr)%shmat_8
        if( self%iseven(i) )then
            pft_ref_8 = self%pfts_refs_even(:,:,iref)
        else
            pft_ref_8 = self%pfts_refs_odd(:,:,iref)
        endif
        call self%gen_shmat_8(ithr, shvec, shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8
        select case(params_glob%cc_objfun)
            case(OBJFUN_CC)
                call self%gen_corr_cc_grad_for_rot_8(    pft_ref_8, pft_ref_tmp_8, iptcl, irot, f, grad)
            case(OBJFUN_EUCLID)
                call self%gen_euclid_grad_for_rot_8(pft_ref_8, pft_ref_tmp_8, iptcl, irot, f, grad)
        end select
    end subroutine gen_corr_grad_for_rot_8

    module subroutine gen_corr_cc_grad_for_rot_8( self, pft_ref, pft_ref_tmp, iptcl, irot, f, grad )
        class(polarft_calc),  intent(inout) :: self
        complex(dp), pointer, intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,              intent(in)    :: iptcl, irot
        real(dp),             intent(out)   :: f, grad(2)
        real(dp) :: sqsum_ref, sqsum_ptcl, denom
        integer  :: k, i
        i         = self%pinds(iptcl)
        sqsum_ref = 0.d0
        f         = 0.d0
        grad      = 0.d0
        sqsum_ptcl = self%ksqsums_ptcls(i)
        call self%rotate_pft(pft_ref, irot, pft_ref_tmp)
        do k = self%kfromto(1),self%kfromto(2)
            sqsum_ref = sqsum_ref + real(k,kind=dp) * sum(real(self%ctfmats(:,k,i)*self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * conjg(pft_ref_tmp(:,k)),dp))
            f         = f         + real(k,kind=dp) * sum(real(self%ctfmats(:,k,i)                     * pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
        enddo
        call self%rotate_pft(pft_ref * dcmplx(0.d0,self%argtransf(:self%pftsz,:)), irot, pft_ref_tmp)
        do k = self%kfromto(1),self%kfromto(2)
            grad(1) = grad(1) + real(k,kind=dp) * sum(real(self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
        enddo
        call self%rotate_pft(pft_ref * dcmplx(0.d0,self%argtransf(self%pftsz+1:,:)), irot, pft_ref_tmp)
        do k = self%kfromto(1),self%kfromto(2)
            grad(2) = grad(2) + real(k,kind=dp) * sum(real(self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
        end do
        denom = dsqrt(sqsum_ref * sqsum_ptcl)
        f     = f    / denom
        grad  = grad / denom
    end subroutine gen_corr_cc_grad_for_rot_8

    module subroutine gen_euclid_grad_for_rot_8( self, pft_ref, pft_ref_tmp, iptcl, irot, f, grad )
        class(polarft_calc),  intent(inout) :: self
        complex(dp), pointer, intent(inout) :: pft_ref(:,:), pft_ref_tmp(:,:)
        integer,              intent(in)    :: iptcl, irot
        real(dp),             intent(out)   :: f, grad(2)
        complex(dp), pointer :: pft_diff(:,:)
        real(dp) :: denom, w
        integer  :: k, i, ithr
        ithr     = omp_get_thread_num() + 1
        i        = self%pinds(iptcl)
        f        = 0.d0
        grad     = 0.d0
        denom    = self%wsqsums_ptcls(i)
        pft_diff => self%heap_vars(ithr)%shmat_8
        call self%rotate_pft(pft_ref, irot, pft_ref_tmp)
        pft_ref_tmp = pft_ref_tmp * self%ctfmats(:,:,i)
        pft_diff = pft_ref_tmp - self%pfts_ptcls(:,:,i) ! Ref(shift + rotation + CTF) - Ptcl
        call self%rotate_pft(pft_ref * dcmplx(0.d0,self%argtransf(:self%pftsz,:)), irot, pft_ref_tmp)
        pft_ref_tmp = pft_ref_tmp * self%ctfmats(:,:,i)
        do k = self%kfromto(1),self%kfromto(2)
            w       = real(k,dp) / real(self%sigma2_noise(k,iptcl))
            f       = f + w * sum(real(pft_diff(:,k)*conjg(pft_diff(:,k)),dp))
            grad(1) = grad(1) + w * real(sum(pft_ref_tmp(:,k) * conjg(pft_diff(:,k))),dp)
        end do
        call self%rotate_pft(pft_ref * dcmplx(0.d0,self%argtransf(self%pftsz+1:,:)), irot, pft_ref_tmp)
        pft_ref_tmp = pft_ref_tmp * self%ctfmats(:,:,i)
        do k = self%kfromto(1),self%kfromto(2)
            w       = real(k,dp) / real(self%sigma2_noise(k,iptcl))
            grad(2) = grad(2) + w * real(sum(pft_ref_tmp(:,k) * conjg(pft_diff(:,k))),dp)
        end do
        f    = dexp( -f / denom )
        grad = -f * 2.d0 * grad / denom
    end subroutine gen_euclid_grad_for_rot_8

    module subroutine gen_corr_grad_only_for_rot_8( self, iref, iptcl, shvec, irot, grad )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        real(dp),            intent(in)    :: shvec(2)
        integer,             intent(in)    :: irot
        real(dp),            intent(out)   :: grad(2)
        complex(dp), pointer :: pft_ref_8(:,:), shmat_8(:,:), pft_ref_tmp_8(:,:)
        real(dp) :: f
        integer  :: ithr, i
        i    = self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat_8       => self%heap_vars(ithr)%shmat_8
        if( self%iseven(i) )then
            pft_ref_8 = self%pfts_refs_even(:,:,iref)
        else
            pft_ref_8 = self%pfts_refs_odd(:,:,iref)
        endif
        call self%gen_shmat_8(ithr, shvec, shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8
        select case(params_glob%cc_objfun)
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
        real(dp) :: sqsum_ref, sqsum_ptcl
        integer  :: k
        sqsum_ref = 0.d0
        grad      = 0.d0
        sqsum_ptcl = self%ksqsums_ptcls(i)
        call self%rotate_pft(pft_ref, irot, pft_ref_tmp)
        do k = self%kfromto(1),self%kfromto(2)
            sqsum_ref = sqsum_ref + real(k,kind=dp) * sum(real(self%ctfmats(:,k,i)*self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * conjg(pft_ref_tmp(:,k)),dp))
        enddo
        call self%rotate_pft(pft_ref * dcmplx(0.d0,self%argtransf(:self%pftsz,:)), irot, pft_ref_tmp)
        do k = self%kfromto(1),self%kfromto(2)
            grad(1) = grad(1) + real(k,kind=dp) * sum(real(self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
        enddo
        call self%rotate_pft(pft_ref * dcmplx(0.d0,self%argtransf(self%pftsz+1:,:)), irot, pft_ref_tmp)
        do k = self%kfromto(1),self%kfromto(2)
            grad(2) = grad(2) + real(k,kind=dp) * sum(real(self%ctfmats(:,k,i) * pft_ref_tmp(:,k) * conjg(self%pfts_ptcls(:,k,i)),dp))
        end do
        grad = grad / dsqrt(sqsum_ref * sqsum_ptcl)
    end subroutine gen_corr_cc_grad_only_for_rot_8

    module subroutine gen_sigma_contrib( self, iref, iptcl, shvec, irot, sigma_contrib)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        real(sp),            intent(in)    :: shvec(2)
        integer,             intent(in)    :: irot
        real(sp),  optional, intent(out)   :: sigma_contrib(self%kfromto(1):self%kfromto(2))
        complex(dp), pointer :: pft_ref_8(:,:), shmat_8(:,:), pft_ref_tmp_8(:,:)
        integer  :: i, ithr
        i    = self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ref_8     => self%heap_vars(ithr)%pft_ref_8
        pft_ref_tmp_8 => self%heap_vars(ithr)%pft_ref_tmp_8
        shmat_8       => self%heap_vars(ithr)%shmat_8
        ! e/o
        if( self%iseven(i) )then
            pft_ref_8 = self%pfts_refs_even(:,:,iref)
        else
            pft_ref_8 = self%pfts_refs_odd(:,:,iref)
        endif
        ! shift
        call self%gen_shmat_8(ithr, real(shvec,dp), shmat_8)
        pft_ref_8 = pft_ref_8 * shmat_8
        ! rotation
        call self%rotate_pft(pft_ref_8, irot, pft_ref_tmp_8)
        ! ctf
        pft_ref_tmp_8 = pft_ref_tmp_8 * real(self%ctfmats(:,:,i),dp)
        ! difference
        pft_ref_tmp_8 = pft_ref_tmp_8 - self%pfts_ptcls(:,:,i)
        ! sigma2
        if( present(sigma_contrib) )then
            sigma_contrib = real(sum(real(pft_ref_tmp_8 * conjg(pft_ref_tmp_8),dp), dim=1) / (2.d0*real(self%pftsz,dp)))
        else
            self%sigma2_noise(self%kfromto(1):self%kfromto(2), iptcl) = real(sum(real(pft_ref_tmp_8 * conjg(pft_ref_tmp_8),dp), dim=1) / (2.d0*real(self%pftsz,dp)))
        endif
    end subroutine gen_sigma_contrib

end submodule simple_polarft_corr
