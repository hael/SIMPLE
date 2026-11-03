# Polar kernel candidates: source archive

Supporting material for the [x86 kernel report](../polar_kernel_codegen_x86_report.md).
These are the source variants measured in the 2026-10-22/23 experiment, retained verbatim
from the report's appendix. They are historical candidates, not production modules.


Every form is HEAD's `src/main/pftc/simple_polarft_corr.f90` (ee5b68be1) with the routines listed under it in place of
HEAD's routines of the same name; a routine shown as "deleted" is removed. The other routines of the file, the
distance conversion included, are HEAD's. The diagnostic d1 is the opposite: f3ccb5fbd's file with HEAD's six
routines of c2 (and nothing else of HEAD). A routine whose source equals one printed earlier is referred to.
The repository now holds c14 with two lines added to the doc comments of `euclid_sums` and `cc_sums` (why the
shell loop is split into helper calls); the code is the same.

### c2: HEAD (ee5b68be1): explicit-shape helpers, called with the memo sections

`gen_corrs`:

```fortran
    module subroutine gen_corrs( self, iref, iptcl, shift, cc )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(sp),                    intent(out)   :: cc(self%nrots)
        complex(sp), pointer :: shmat(:,:)
        real(sp) :: shift_mag_sq
        integer  :: i, ithr, k, kk, k0, ieo
        ithr    =  omp_get_thread_num() + 1
        i       =  self%pinds(iptcl)
        k0      =  self%kfromto(1)
        ieo     =  merge(REF_EVEN, REF_ODD, self%iseven(i))
        shmat   => self%heap_vars(ithr)%shmat
        ! IFFT( sum_k FT(CTF2) x FT(REF2) )
        call cc_sums(self%pftsz+1, self%pftsz+1, self%nk, self%ft_ctf2(:,:,i), self%ft_ref2(:,:,iref,ieo), .false.,&
            &self%crvec1(ithr)%c)
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots)
        ! IFFT( sum_k FT(X.CTF) x FT(S.REF)* ), the memoized FT(REF) when there is no shift
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        if( shift_mag_sq > SHERRSQ )then
            call self%gen_shmat4aln(ithr, shift, shmat)
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%cmat2_many(ithr)%c(1:self%pftsz,kk) = shmat(:,k) * self%pfts_refs(:,k,iref,ieo)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots,kk) = &
                    conjg(self%cmat2_many(ithr)%c(1:self%pftsz,kk))
            enddo
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
            call cc_sums(self%pftsz+1, self%nrots, self%nk, self%ft_ptcl_ctf(:,:,i), self%cmat2_many(ithr)%c, .true.,&
                &self%crvec1(ithr)%c)
        else
            call cc_sums(self%pftsz+1, self%pftsz+1, self%nk, self%ft_ptcl_ctf(:,:,i), self%ft_ref(:,:,iref,ieo), .true.,&
                &self%crvec1(ithr)%c)
        endif
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        ! Final correlation computation
        cc = self%crvec1(ithr)%r(1:self%nrots) / &
            &sqrt(cc * real(self%sqsums_ptcls(i), sp) * real(2*self%nrots, sp))
    end subroutine gen_corrs
```

`gen_euclid_crvec`:

```fortran
    subroutine gen_euclid_crvec( self, iref, iptcl, shift, norm, ithr )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(dp),                    intent(out)   :: norm
        integer,                     intent(out)   :: ithr
        complex(sp), pointer :: shmat(:,:)
        real(sp)    :: shift_mag_sq
        integer     :: k, i, kk, k0, ieo, ih
        ithr         = omp_get_thread_num() + 1
        i            = self%pinds(iptcl)
        ih = merge(1, 2, self%iseven(i))
        k0           = self%kfromto(1)
        ieo          = merge(REF_EVEN, REF_ODD, self%iseven(i))
        norm         = self%wsqsum_ptcl_now(iptcl) * real(2*self%nrots, dp)
        shmat        => self%heap_vars(ithr)%shmat
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        if( shift_mag_sq > SHERRSQ )then
            call self%gen_shmat4aln(ithr, shift, shmat)
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%cmat2_many(ithr)%c(1:self%pftsz,kk) = shmat(:,k) * self%pfts_refs(:,k,iref,ieo)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots,kk) = &
                    conjg(self%cmat2_many(ithr)%c(1:self%pftsz,kk))
            enddo
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
            call euclid_sums(self%pftsz, self%pftsz+1, self%nrots, self%nk, self%ft_ctf2(:,:,i),&
                &self%ft_ref2(:,:,iref,ieo), self%ft_ptcl_ctf(:,:,i), self%cmat2_many(ithr)%c, self%sigma2_wk8(:,ih),&
                &self%crvec1(ithr)%c)
        else
            call euclid_sums(self%pftsz, self%pftsz+1, self%pftsz+1, self%nk, self%ft_ctf2(:,:,i),&
                &self%ft_ref2(:,:,iref,ieo), self%ft_ptcl_ctf(:,:,i), self%ft_ref(:,:,iref,ieo), self%sigma2_wk8(:,ih),&
                &self%crvec1(ithr)%c)
        endif
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
    end subroutine gen_euclid_crvec
```

`euclid_sums`:

```fortran
    !> The euclid sums of one pair over angular bins 1..np, sum_k w_k (FT(CTF2) x FT(REF2) - 2 FT(X.CTF) x FT(REF)*),
    !! into acc. Plain array dummies: the compiler then knows unit stride and no overlap and vectorizes the bins.
    subroutine euclid_sums( np, n1, ldr, nk, ctf2, ref2, xctf, ref, wk, acc )
        integer,     intent(in)    :: np, n1, ldr, nk
        complex(sp), intent(in)    :: ctf2(n1,nk), ref2(n1,nk), xctf(n1,nk), ref(ldr,nk)
        real(dp),    intent(in)    :: wk(nk)
        complex(sp), intent(inout) :: acc(n1)
        complex(sp) :: c
        integer     :: kk, p
        acc = cmplx(0.,0.,kind=sp)
        do kk = 1,nk
            do p = 1,np
                c = ctf2(p,kk) * ref2(p,kk)
                c = c - 2.0 * xctf(p,kk) * conjg(ref(p,kk))
                acc(p) = cmplx(acc(p) + wk(kk) * c, kind=sp)
            enddo
        enddo
    end subroutine euclid_sums
```

`gen_euclid_crvec_bank`:

```fortran
    subroutine gen_euclid_crvec_bank( self, iref, iptcl, bank, norm, ithr )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        complex(sp),         intent(in)    :: bank(:,:)
        real(dp),            intent(out)   :: norm
        integer,             intent(out)   :: ithr
        integer :: i, ieo, ih
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        ih   = merge(1, 2, self%iseven(i))
        ieo  = merge(REF_EVEN, REF_ODD, self%iseven(i))
        norm = self%wsqsum_ptcl_now(iptcl) * real(2*self%nrots, dp)
        call euclid_sums(self%pftsz+1, self%pftsz+1, self%pftsz+1, self%nk, self%ft_ctf2(:,:,i),&
            &self%ft_ref2(:,:,iref,ieo), bank, self%ft_ref(:,:,iref,ieo), self%sigma2_wk8(:,ih), self%crvec1(ithr)%c)
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
    end subroutine gen_euclid_crvec_bank
```

`gen_corrs_bank`:

```fortran
    subroutine gen_corrs_bank( self, iref, iptcl, bank, cc )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        complex(sp),         intent(in)    :: bank(:,:)
        real(sp),            intent(out)   :: cc(self%nrots)
        integer :: i, ithr, ieo
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        ieo  = merge(REF_EVEN, REF_ODD, self%iseven(i))
        ! IFFT( sum_k FT(CTF2) x FT(REF2) )
        call cc_sums(self%pftsz+1, self%pftsz+1, self%nk, self%ft_ctf2(:,:,i), self%ft_ref2(:,:,iref,ieo), .false.,&
            &self%crvec1(ithr)%c)
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots)
        ! IFFT( sum_k FT(S.X.CTF) x FT(REF)* )
        call cc_sums(self%pftsz+1, self%pftsz+1, self%nk, bank, self%ft_ref(:,:,iref,ieo), .true., self%crvec1(ithr)%c)
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots) / &
            &sqrt(cc * real(self%sqsums_ptcls(i), sp) * real(2*self%nrots, sp))
    end subroutine gen_corrs_bank
```

`cc_sums`:

```fortran
    !> sum_k a_k x b_k (l_conj: a_k x b_k*) over the n1 angular bins into acc, plain array dummies as in euclid_sums
    subroutine cc_sums( n1, ldb, nk, a, b, l_conj, acc )
        integer,     intent(in)    :: n1, ldb, nk
        complex(sp), intent(in)    :: a(n1,nk), b(ldb,nk)
        logical,     intent(in)    :: l_conj
        complex(sp), intent(inout) :: acc(n1)
        integer :: kk
        acc = cmplx(0.,0.,kind=sp)
        if( l_conj )then
            do kk = 1,nk
                acc = acc + a(:,kk) * conjg(b(1:n1,kk))
            enddo
        else
            do kk = 1,nk
                acc = acc + a(:,kk) * b(1:n1,kk)
            enddo
        endif
    end subroutine cc_sums
```

### c1: the old inline loops of f3ccb5fbd (with HEAD's distance conversion); helpers deleted

`gen_corrs`:

```fortran
    module subroutine gen_corrs( self, iref, iptcl, shift, cc )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(sp),                    intent(out)   :: cc(self%nrots)
        complex(sp), pointer :: shmat(:,:)
        real(sp) :: shift_mag_sq
        integer  :: i, ithr, k, kk, k0, ieo
        logical  :: needs_shift
        ithr    =  omp_get_thread_num() + 1
        i       =  self%pinds(iptcl)
        k0      =  self%kfromto(1)
        ieo     =  merge(REF_EVEN, REF_ODD, self%iseven(i))
        shmat   => self%heap_vars(ithr)%shmat
        ! ========================================================================
        ! Prepare FFT(S.REF) if needed
        ! ========================================================================
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        needs_shift  = shift_mag_sq > SHERRSQ
        if (needs_shift) then
            ! Generate shift matrix
            call self%gen_shmat4aln(ithr, shift, shmat)
            ! Shift reference directly into cmat2_many
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%cmat2_many(ithr)%c(1:self%pftsz,kk) = shmat(:,k) * self%pfts_refs(:,k,iref,ieo)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots,kk) = &
                    conjg(self%cmat2_many(ithr)%c(1:self%pftsz,kk))
            enddo
            ! Calculate FFT(S.REF)
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
        else
            ! no shift, memoized reference is used
            self%cmat2_many(ithr)%c(1:self%pftsz+1,1:self%nk) = &
                self%ft_ref(:,self%kfromto(1):self%kfromto(2),iref,ieo)
        endif
        ! ========================================================================
        ! Single IFFT #1: IFFT( sum_k FT(CTF2) x FT(REF2) )
        ! ========================================================================
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        do k = self%kfromto(1), self%kfromto(2)
            self%crvec1(ithr)%c = self%crvec1(ithr)%c + self%ft_ctf2(:,k,i) * self%ft_ref2(:,k,iref,ieo)
        enddo
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots)
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
        ! Final correlation computation
        cc = self%crvec1(ithr)%r(1:self%nrots) / &
            &sqrt(cc * real(self%sqsums_ptcls(i), sp) * real(2*self%nrots, sp))
    end subroutine gen_corrs
```

`gen_euclid_crvec`:

```fortran
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
        integer     :: k, i, kk, k0, p, ieo, ih
        ithr         = omp_get_thread_num() + 1
        i            = self%pinds(iptcl)
        ih = merge(1, 2, self%iseven(i))
        k0           = self%kfromto(1)
        ieo          = merge(REF_EVEN, REF_ODD, self%iseven(i))
        norm         = self%wsqsum_ptcl_now(iptcl) * real(2*self%nrots, dp)
        shmat        => self%heap_vars(ithr)%shmat
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        if( shift_mag_sq > SHERRSQ )then
            call self%gen_shmat4aln(ithr, shift, shmat)
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%cmat2_many(ithr)%c(1:self%pftsz,kk) = shmat(:,k) * self%pfts_refs(:,k,iref,ieo)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots,kk) = &
                    conjg(self%cmat2_many(ithr)%c(1:self%pftsz,kk))
            enddo
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
            self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                wk = self%sigma2_wk8(k,ih)
                do p = 1,self%pftsz
                    c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                    c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%cmat2_many(ithr)%c(p,kk))
                    self%crvec1(ithr)%c(p) = cmplx(self%crvec1(ithr)%c(p) + wk * c, &
                        &kind=kind(self%crvec1(ithr)%c))
                enddo
            enddo
        else
            self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
            do k = self%kfromto(1), self%kfromto(2)
                wk = self%sigma2_wk8(k,ih)
                do p = 1,self%pftsz
                    c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                    c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%ft_ref(p,k,iref,ieo))
                    self%crvec1(ithr)%c(p) = cmplx(self%crvec1(ithr)%c(p) + wk * c, &
                        &kind=kind(self%crvec1(ithr)%c))
                enddo
            enddo
        endif
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
    end subroutine gen_euclid_crvec
```

- `euclid_sums`: deleted

`gen_euclid_crvec_bank`:

```fortran
    subroutine gen_euclid_crvec_bank( self, iref, iptcl, bank, norm, ithr )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        complex(sp),         intent(in)    :: bank(:,:)
        real(dp),            intent(out)   :: norm
        integer,             intent(out)   :: ithr
        complex(sp) :: c
        real(dp)    :: wk
        integer     :: k, i, kk, k0, p, ieo, ih
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        ih   = merge(1, 2, self%iseven(i))
        k0   = self%kfromto(1)
        ieo  = merge(REF_EVEN, REF_ODD, self%iseven(i))
        norm = self%wsqsum_ptcl_now(iptcl) * real(2*self%nrots, dp)
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        do k = self%kfromto(1), self%kfromto(2)
            kk = k - k0 + 1
            wk = self%sigma2_wk8(k,ih)
            do p = 1,self%pftsz+1
                c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                c = c - 2.0 * bank(p,kk) * conjg(self%ft_ref(p,k,iref,ieo))
                self%crvec1(ithr)%c(p) = cmplx(self%crvec1(ithr)%c(p) + wk * c, &
                    &kind=kind(self%crvec1(ithr)%c))
            enddo
        enddo
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
    end subroutine gen_euclid_crvec_bank
```

`gen_corrs_bank`:

```fortran
    subroutine gen_corrs_bank( self, iref, iptcl, bank, cc )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        complex(sp),         intent(in)    :: bank(:,:)
        real(sp),            intent(out)   :: cc(self%nrots)
        integer :: i, ithr, k, kk, k0, ieo
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        k0   = self%kfromto(1)
        ieo  = merge(REF_EVEN, REF_ODD, self%iseven(i))
        ! IFFT( sum_k FT(CTF2) x FT(REF2) )
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        do k = self%kfromto(1), self%kfromto(2)
            self%crvec1(ithr)%c = self%crvec1(ithr)%c + self%ft_ctf2(:,k,i) * self%ft_ref2(:,k,iref,ieo)
        enddo
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots)
        ! IFFT( sum_k FT(S.X.CTF) x FT(REF)* )
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        do k = self%kfromto(1), self%kfromto(2)
            kk = k - k0 + 1
            self%crvec1(ithr)%c = self%crvec1(ithr)%c + bank(:,kk) * conjg(self%ft_ref(:,k,iref,ieo))
        end do
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots) / &
            &sqrt(cc * real(self%sqsums_ptcls(i), sp) * real(2*self%nrots, sp))
    end subroutine gen_corrs_bank
```

- `cc_sums`: deleted

### c3: old inline loops with a local `pointer, contiguous` alias of the FFT buffers; FT(REF) read directly without a shift

`gen_corrs`:

```fortran
    module subroutine gen_corrs( self, iref, iptcl, shift, cc )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(sp),                    intent(out)   :: cc(self%nrots)
        complex(sp), pointer             :: shmat(:,:)
        complex(sp), pointer, contiguous :: acc(:), sref(:,:)
        real(sp) :: shift_mag_sq
        integer  :: i, ithr, k, kk, k0, ieo
        ithr    =  omp_get_thread_num() + 1
        i       =  self%pinds(iptcl)
        k0      =  self%kfromto(1)
        ieo     =  merge(REF_EVEN, REF_ODD, self%iseven(i))
        shmat   => self%heap_vars(ithr)%shmat
        ! contiguous aliases of the per-thread FFT buffers: unit stride without a call
        acc     => self%crvec1(ithr)%c
        ! IFFT( sum_k FT(CTF2) x FT(REF2) )
        acc = cmplx(0.,0.,kind=sp)
        do k = self%kfromto(1), self%kfromto(2)
            acc = acc + self%ft_ctf2(:,k,i) * self%ft_ref2(:,k,iref,ieo)
        enddo
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots)
        ! IFFT( sum_k FT(X.CTF) x FT(S.REF)* ), the memoized FT(REF) when there is no shift
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        acc = cmplx(0.,0.,kind=sp)
        if( shift_mag_sq > SHERRSQ )then
            call self%gen_shmat4aln(ithr, shift, shmat)
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%cmat2_many(ithr)%c(1:self%pftsz,kk) = shmat(:,k) * self%pfts_refs(:,k,iref,ieo)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots,kk) = &
                    conjg(self%cmat2_many(ithr)%c(1:self%pftsz,kk))
            enddo
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
            sref => self%cmat2_many(ithr)%c
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                acc = acc + self%ft_ptcl_ctf(:,k,i) * conjg(sref(1:self%pftsz+1,kk))
            enddo
        else
            do k = self%kfromto(1), self%kfromto(2)
                acc = acc + self%ft_ptcl_ctf(:,k,i) * conjg(self%ft_ref(:,k,iref,ieo))
            enddo
        endif
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        ! Final correlation computation
        cc = self%crvec1(ithr)%r(1:self%nrots) / &
            &sqrt(cc * real(self%sqsums_ptcls(i), sp) * real(2*self%nrots, sp))
    end subroutine gen_corrs
```

`gen_euclid_crvec`:

```fortran
    subroutine gen_euclid_crvec( self, iref, iptcl, shift, norm, ithr )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(dp),                    intent(out)   :: norm
        integer,                     intent(out)   :: ithr
        complex(sp), pointer             :: shmat(:,:)
        complex(sp), pointer, contiguous :: acc(:), sref(:,:)
        complex(sp) :: c
        real(dp)    :: wk
        real(sp)    :: shift_mag_sq
        integer     :: k, i, kk, k0, p, ieo, ih
        ithr         = omp_get_thread_num() + 1
        i            = self%pinds(iptcl)
        ih = merge(1, 2, self%iseven(i))
        k0           = self%kfromto(1)
        ieo          = merge(REF_EVEN, REF_ODD, self%iseven(i))
        norm         = self%wsqsum_ptcl_now(iptcl) * real(2*self%nrots, dp)
        shmat        => self%heap_vars(ithr)%shmat
        ! contiguous aliases of the per-thread FFT buffers: unit stride without a call
        acc          => self%crvec1(ithr)%c
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        if( shift_mag_sq > SHERRSQ )then
            call self%gen_shmat4aln(ithr, shift, shmat)
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%cmat2_many(ithr)%c(1:self%pftsz,kk) = shmat(:,k) * self%pfts_refs(:,k,iref,ieo)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots,kk) = &
                    conjg(self%cmat2_many(ithr)%c(1:self%pftsz,kk))
            enddo
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
            sref => self%cmat2_many(ithr)%c
            acc  = cmplx(0.,0.,kind=sp)
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                wk = self%sigma2_wk8(k,ih)
                do p = 1,self%pftsz
                    c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                    c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(sref(p,kk))
                    acc(p) = cmplx(acc(p) + wk * c, kind=sp)
                enddo
            enddo
        else
            acc = cmplx(0.,0.,kind=sp)
            do k = self%kfromto(1), self%kfromto(2)
                wk = self%sigma2_wk8(k,ih)
                do p = 1,self%pftsz
                    c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                    c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%ft_ref(p,k,iref,ieo))
                    acc(p) = cmplx(acc(p) + wk * c, kind=sp)
                enddo
            enddo
        endif
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
    end subroutine gen_euclid_crvec
```

- `euclid_sums`: deleted

`gen_euclid_crvec_bank`:

```fortran
    subroutine gen_euclid_crvec_bank( self, iref, iptcl, bank, norm, ithr )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        complex(sp),         intent(in)    :: bank(:,:)
        real(dp),            intent(out)   :: norm
        integer,             intent(out)   :: ithr
        complex(sp), pointer, contiguous :: acc(:)
        complex(sp) :: c
        real(dp)    :: wk
        integer     :: k, i, kk, k0, p, ieo, ih
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        ih   = merge(1, 2, self%iseven(i))
        k0   = self%kfromto(1)
        ieo  = merge(REF_EVEN, REF_ODD, self%iseven(i))
        norm = self%wsqsum_ptcl_now(iptcl) * real(2*self%nrots, dp)
        acc  => self%crvec1(ithr)%c
        acc  = cmplx(0.,0.,kind=sp)
        do k = self%kfromto(1), self%kfromto(2)
            kk = k - k0 + 1
            wk = self%sigma2_wk8(k,ih)
            do p = 1,self%pftsz+1
                c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                c = c - 2.0 * bank(p,kk) * conjg(self%ft_ref(p,k,iref,ieo))
                acc(p) = cmplx(acc(p) + wk * c, kind=sp)
            enddo
        enddo
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
    end subroutine gen_euclid_crvec_bank
```

`gen_corrs_bank`:

```fortran
    subroutine gen_corrs_bank( self, iref, iptcl, bank, cc )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        complex(sp),         intent(in)    :: bank(:,:)
        real(sp),            intent(out)   :: cc(self%nrots)
        complex(sp), pointer, contiguous :: acc(:)
        integer :: i, ithr, k, kk, k0, ieo
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        k0   = self%kfromto(1)
        ieo  = merge(REF_EVEN, REF_ODD, self%iseven(i))
        acc  => self%crvec1(ithr)%c
        ! IFFT( sum_k FT(CTF2) x FT(REF2) )
        acc = cmplx(0.,0.,kind=sp)
        do k = self%kfromto(1), self%kfromto(2)
            acc = acc + self%ft_ctf2(:,k,i) * self%ft_ref2(:,k,iref,ieo)
        enddo
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots)
        ! IFFT( sum_k FT(S.X.CTF) x FT(REF)* )
        acc = cmplx(0.,0.,kind=sp)
        do k = self%kfromto(1), self%kfromto(2)
            kk = k - k0 + 1
            acc = acc + bank(:,kk) * conjg(self%ft_ref(:,k,iref,ieo))
        end do
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots) / &
            &sqrt(cc * real(self%sqsums_ptcls(i), sp) * real(2*self%nrots, sp))
    end subroutine gen_corrs_bank
```

- `cc_sums`: deleted

### c4: HEAD's helpers with assumed-shape `contiguous` dummies

`gen_corrs`:

```fortran
    module subroutine gen_corrs( self, iref, iptcl, shift, cc )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(sp),                    intent(out)   :: cc(self%nrots)
        complex(sp), pointer :: shmat(:,:)
        real(sp) :: shift_mag_sq
        integer  :: i, ithr, k, kk, k0, ieo
        ithr    =  omp_get_thread_num() + 1
        i       =  self%pinds(iptcl)
        k0      =  self%kfromto(1)
        ieo     =  merge(REF_EVEN, REF_ODD, self%iseven(i))
        shmat   => self%heap_vars(ithr)%shmat
        ! IFFT( sum_k FT(CTF2) x FT(REF2) )
        call cc_sums(self%ft_ctf2(:,:,i), self%ft_ref2(:,:,iref,ieo), .false.,&
            &self%crvec1(ithr)%c)
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots)
        ! IFFT( sum_k FT(X.CTF) x FT(S.REF)* ), the memoized FT(REF) when there is no shift
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        if( shift_mag_sq > SHERRSQ )then
            call self%gen_shmat4aln(ithr, shift, shmat)
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%cmat2_many(ithr)%c(1:self%pftsz,kk) = shmat(:,k) * self%pfts_refs(:,k,iref,ieo)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots,kk) = &
                    conjg(self%cmat2_many(ithr)%c(1:self%pftsz,kk))
            enddo
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
            call cc_sums(self%ft_ptcl_ctf(:,:,i), self%cmat2_many(ithr)%c, .true.,&
                &self%crvec1(ithr)%c)
        else
            call cc_sums(self%ft_ptcl_ctf(:,:,i), self%ft_ref(:,:,iref,ieo), .true.,&
                &self%crvec1(ithr)%c)
        endif
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        ! Final correlation computation
        cc = self%crvec1(ithr)%r(1:self%nrots) / &
            &sqrt(cc * real(self%sqsums_ptcls(i), sp) * real(2*self%nrots, sp))
    end subroutine gen_corrs
```

`gen_euclid_crvec`:

```fortran
    subroutine gen_euclid_crvec( self, iref, iptcl, shift, norm, ithr )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(dp),                    intent(out)   :: norm
        integer,                     intent(out)   :: ithr
        complex(sp), pointer :: shmat(:,:)
        real(sp)    :: shift_mag_sq
        integer     :: k, i, kk, k0, ieo, ih
        ithr         = omp_get_thread_num() + 1
        i            = self%pinds(iptcl)
        ih = merge(1, 2, self%iseven(i))
        k0           = self%kfromto(1)
        ieo          = merge(REF_EVEN, REF_ODD, self%iseven(i))
        norm         = self%wsqsum_ptcl_now(iptcl) * real(2*self%nrots, dp)
        shmat        => self%heap_vars(ithr)%shmat
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        if( shift_mag_sq > SHERRSQ )then
            call self%gen_shmat4aln(ithr, shift, shmat)
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%cmat2_many(ithr)%c(1:self%pftsz,kk) = shmat(:,k) * self%pfts_refs(:,k,iref,ieo)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots,kk) = &
                    conjg(self%cmat2_many(ithr)%c(1:self%pftsz,kk))
            enddo
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
            call euclid_sums(self%pftsz, self%ft_ctf2(:,:,i),&
                &self%ft_ref2(:,:,iref,ieo), self%ft_ptcl_ctf(:,:,i), self%cmat2_many(ithr)%c, self%sigma2_wk8(:,ih),&
                &self%crvec1(ithr)%c)
        else
            call euclid_sums(self%pftsz, self%ft_ctf2(:,:,i),&
                &self%ft_ref2(:,:,iref,ieo), self%ft_ptcl_ctf(:,:,i), self%ft_ref(:,:,iref,ieo), self%sigma2_wk8(:,ih),&
                &self%crvec1(ithr)%c)
        endif
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
    end subroutine gen_euclid_crvec
```

`euclid_sums`:

```fortran
    !> The euclid sums of one pair over angular bins 1..np, sum_k w_k (FT(CTF2) x FT(REF2) - 2 FT(X.CTF) x FT(REF)*),
    !! into acc. Contiguous dummies: the compiler then knows unit stride and no overlap and vectorizes the bins.
    subroutine euclid_sums( np, ctf2, ref2, xctf, ref, wk, acc )
        integer,                  intent(in)    :: np
        complex(sp), contiguous,  intent(in)    :: ctf2(:,:), ref2(:,:), xctf(:,:), ref(:,:)
        real(dp),    contiguous,  intent(in)    :: wk(:)
        complex(sp), contiguous,  intent(inout) :: acc(:)
        complex(sp) :: c
        integer     :: kk, p
        acc = cmplx(0.,0.,kind=sp)
        do kk = 1,size(ctf2,2)
            do p = 1,np
                c = ctf2(p,kk) * ref2(p,kk)
                c = c - 2.0 * xctf(p,kk) * conjg(ref(p,kk))
                acc(p) = cmplx(acc(p) + wk(kk) * c, kind=sp)
            enddo
        enddo
    end subroutine euclid_sums
```

`gen_euclid_crvec_bank`:

```fortran
    subroutine gen_euclid_crvec_bank( self, iref, iptcl, bank, norm, ithr )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        complex(sp),         intent(in)    :: bank(:,:)
        real(dp),            intent(out)   :: norm
        integer,             intent(out)   :: ithr
        integer :: i, ieo, ih
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        ih   = merge(1, 2, self%iseven(i))
        ieo  = merge(REF_EVEN, REF_ODD, self%iseven(i))
        norm = self%wsqsum_ptcl_now(iptcl) * real(2*self%nrots, dp)
        call euclid_sums(self%pftsz+1, self%ft_ctf2(:,:,i),&
            &self%ft_ref2(:,:,iref,ieo), bank, self%ft_ref(:,:,iref,ieo), self%sigma2_wk8(:,ih), self%crvec1(ithr)%c)
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
    end subroutine gen_euclid_crvec_bank
```

`gen_corrs_bank`:

```fortran
    subroutine gen_corrs_bank( self, iref, iptcl, bank, cc )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        complex(sp),         intent(in)    :: bank(:,:)
        real(sp),            intent(out)   :: cc(self%nrots)
        integer :: i, ithr, ieo
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        ieo  = merge(REF_EVEN, REF_ODD, self%iseven(i))
        ! IFFT( sum_k FT(CTF2) x FT(REF2) )
        call cc_sums(self%ft_ctf2(:,:,i), self%ft_ref2(:,:,iref,ieo), .false.,&
            &self%crvec1(ithr)%c)
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots)
        ! IFFT( sum_k FT(S.X.CTF) x FT(REF)* )
        call cc_sums(bank, self%ft_ref(:,:,iref,ieo), .true., self%crvec1(ithr)%c)
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots) / &
            &sqrt(cc * real(self%sqsums_ptcls(i), sp) * real(2*self%nrots, sp))
    end subroutine gen_corrs_bank
```

`cc_sums`:

```fortran
    !> sum_k a_k x b_k (l_conj: a_k x b_k*) over the angular bins of acc, contiguous dummies as in euclid_sums
    subroutine cc_sums( a, b, l_conj, acc )
        complex(sp), contiguous, intent(in)    :: a(:,:), b(:,:)
        logical,                 intent(in)    :: l_conj
        complex(sp), contiguous, intent(inout) :: acc(:)
        integer :: kk, n1
        n1  = size(acc)
        acc = cmplx(0.,0.,kind=sp)
        if( l_conj )then
            do kk = 1,size(a,2)
                acc = acc + a(:,kk) * conjg(b(1:n1,kk))
            enddo
        else
            do kk = 1,size(a,2)
                acc = acc + a(:,kk) * b(1:n1,kk)
            enddo
        endif
    end subroutine cc_sums
```

### c5: HEAD's helpers accumulating into a local work array, copied once

`euclid_sums`:

```fortran
    !> The euclid sums of one pair over angular bins 1..np, sum_k w_k (FT(CTF2) x FT(REF2) - 2 FT(X.CTF) x FT(REF)*),
    !! into acc. Accumulated in a local array, copied to acc once.
    subroutine euclid_sums( np, n1, ldr, nk, ctf2, ref2, xctf, ref, wk, acc )
        integer,     intent(in)    :: np, n1, ldr, nk
        complex(sp), intent(in)    :: ctf2(n1,nk), ref2(n1,nk), xctf(n1,nk), ref(ldr,nk)
        real(dp),    intent(in)    :: wk(nk)
        complex(sp), intent(inout) :: acc(n1)
        complex(sp) :: c, w(n1)
        integer     :: kk, p
        w = cmplx(0.,0.,kind=sp)
        do kk = 1,nk
            do p = 1,np
                c = ctf2(p,kk) * ref2(p,kk)
                c = c - 2.0 * xctf(p,kk) * conjg(ref(p,kk))
                w(p) = cmplx(w(p) + wk(kk) * c, kind=sp)
            enddo
        enddo
        acc = w
    end subroutine euclid_sums
```

`cc_sums`:

```fortran
    !> sum_k a_k x b_k (l_conj: a_k x b_k*) over the n1 angular bins into acc, accumulated in a local array
    subroutine cc_sums( n1, ldb, nk, a, b, l_conj, acc )
        integer,     intent(in)    :: n1, ldb, nk
        complex(sp), intent(in)    :: a(n1,nk), b(ldb,nk)
        logical,     intent(in)    :: l_conj
        complex(sp), intent(inout) :: acc(n1)
        complex(sp) :: w(n1)
        integer :: kk
        w = cmplx(0.,0.,kind=sp)
        if( l_conj )then
            do kk = 1,nk
                w = w + a(:,kk) * conjg(b(1:n1,kk))
            enddo
        else
            do kk = 1,nk
                w = w + a(:,kk) * b(1:n1,kk)
            enddo
        endif
        acc = w
    end subroutine cc_sums
```

- `gen_corrs`, `gen_euclid_crvec`, `gen_euclid_crvec_bank`, `gen_corrs_bank`: as c2 (HEAD)

### c6: HEAD's helpers, two shells per pass written out

`euclid_sums`:

```fortran
    !> The euclid sums of one pair over angular bins 1..np, sum_k w_k (FT(CTF2) x FT(REF2) - 2 FT(X.CTF) x FT(REF)*),
    !! into acc. Plain array dummies: the compiler then knows unit stride and no overlap and vectorizes the bins.
    !! Two shells per pass over acc, each still rounded to single precision in shell order.
    subroutine euclid_sums( np, n1, ldr, nk, ctf2, ref2, xctf, ref, wk, acc )
        integer,     intent(in)    :: np, n1, ldr, nk
        complex(sp), intent(in)    :: ctf2(n1,nk), ref2(n1,nk), xctf(n1,nk), ref(ldr,nk)
        real(dp),    intent(in)    :: wk(nk)
        complex(sp), intent(inout) :: acc(n1)
        complex(sp) :: c1, c2
        integer     :: kk, p
        acc = cmplx(0.,0.,kind=sp)
        do kk = 1,nk-1,2
            do p = 1,np
                c1 = ctf2(p,kk) * ref2(p,kk)
                c1 = c1 - 2.0 * xctf(p,kk) * conjg(ref(p,kk))
                c2 = ctf2(p,kk+1) * ref2(p,kk+1)
                c2 = c2 - 2.0 * xctf(p,kk+1) * conjg(ref(p,kk+1))
                acc(p) = cmplx(cmplx(acc(p) + wk(kk) * c1, kind=sp) + wk(kk+1) * c2, kind=sp)
            enddo
        enddo
        if( mod(nk,2) == 1 )then
            do p = 1,np
                c1 = ctf2(p,nk) * ref2(p,nk)
                c1 = c1 - 2.0 * xctf(p,nk) * conjg(ref(p,nk))
                acc(p) = cmplx(acc(p) + wk(nk) * c1, kind=sp)
            enddo
        endif
    end subroutine euclid_sums
```

`cc_sums`:

```fortran
    !> sum_k a_k x b_k (l_conj: a_k x b_k*) over the n1 angular bins into acc, plain array dummies as in euclid_sums,
    !! two shells per pass, summed in shell order
    subroutine cc_sums( n1, ldb, nk, a, b, l_conj, acc )
        integer,     intent(in)    :: n1, ldb, nk
        complex(sp), intent(in)    :: a(n1,nk), b(ldb,nk)
        logical,     intent(in)    :: l_conj
        complex(sp), intent(inout) :: acc(n1)
        integer :: kk
        acc = cmplx(0.,0.,kind=sp)
        if( l_conj )then
            do kk = 1,nk-1,2
                acc = (acc + a(:,kk) * conjg(b(1:n1,kk))) + a(:,kk+1) * conjg(b(1:n1,kk+1))
            enddo
            if( mod(nk,2) == 1 ) acc = acc + a(:,nk) * conjg(b(1:n1,nk))
        else
            do kk = 1,nk-1,2
                acc = (acc + a(:,kk) * b(1:n1,kk)) + a(:,kk+1) * b(1:n1,kk+1)
            enddo
            if( mod(nk,2) == 1 ) acc = acc + a(:,nk) * b(1:n1,nk)
        endif
    end subroutine cc_sums
```

- `gen_corrs`, `gen_euclid_crvec`, `gen_euclid_crvec_bank`, `gen_corrs_bank`: as c2 (HEAD)

### c7: old inline loops, two shells per pass written out

`gen_corrs`:

```fortran
    module subroutine gen_corrs( self, iref, iptcl, shift, cc )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(sp),                    intent(out)   :: cc(self%nrots)
        complex(sp), pointer :: shmat(:,:)
        real(sp) :: shift_mag_sq
        integer  :: i, ithr, k, kk, k0, ieo
        logical  :: needs_shift
        ithr    =  omp_get_thread_num() + 1
        i       =  self%pinds(iptcl)
        k0      =  self%kfromto(1)
        ieo     =  merge(REF_EVEN, REF_ODD, self%iseven(i))
        shmat   => self%heap_vars(ithr)%shmat
        ! ========================================================================
        ! Prepare FFT(S.REF) if needed
        ! ========================================================================
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        needs_shift  = shift_mag_sq > SHERRSQ
        if (needs_shift) then
            ! Generate shift matrix
            call self%gen_shmat4aln(ithr, shift, shmat)
            ! Shift reference directly into cmat2_many
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%cmat2_many(ithr)%c(1:self%pftsz,kk) = shmat(:,k) * self%pfts_refs(:,k,iref,ieo)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots,kk) = &
                    conjg(self%cmat2_many(ithr)%c(1:self%pftsz,kk))
            enddo
            ! Calculate FFT(S.REF)
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
        else
            ! no shift, memoized reference is used
            self%cmat2_many(ithr)%c(1:self%pftsz+1,1:self%nk) = &
                self%ft_ref(:,self%kfromto(1):self%kfromto(2),iref,ieo)
        endif
        ! ========================================================================
        ! Single IFFT #1: IFFT( sum_k FT(CTF2) x FT(REF2) )
        ! ========================================================================
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        do k = self%kfromto(1), self%kfromto(2)-1, 2
            self%crvec1(ithr)%c = (self%crvec1(ithr)%c + self%ft_ctf2(:,k,i) * self%ft_ref2(:,k,iref,ieo)) +&
                &self%ft_ctf2(:,k+1,i) * self%ft_ref2(:,k+1,iref,ieo)
        enddo
        if( mod(self%nk,2) == 1 ) self%crvec1(ithr)%c = self%crvec1(ithr)%c +&
            &self%ft_ctf2(:,self%kfromto(2),i) * self%ft_ref2(:,self%kfromto(2),iref,ieo)
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots)
        ! ========================================================================
        ! Single IFFT #2: IFFT( sum_k FT(X.CTF) x FT(S.REF)* )
        ! ========================================================================
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        do k = self%kfromto(1), self%kfromto(2)-1, 2
            kk = k - k0 + 1
            self%crvec1(ithr)%c = (self%crvec1(ithr)%c + &
                self%ft_ptcl_ctf(:,k,i) * conjg(self%cmat2_many(ithr)%c(1:self%pftsz+1,kk))) +&
                &self%ft_ptcl_ctf(:,k+1,i) * conjg(self%cmat2_many(ithr)%c(1:self%pftsz+1,kk+1))
        end do
        if( mod(self%nk,2) == 1 ) self%crvec1(ithr)%c = self%crvec1(ithr)%c + &
            &self%ft_ptcl_ctf(:,self%kfromto(2),i) * conjg(self%cmat2_many(ithr)%c(1:self%pftsz+1,self%nk))
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        ! Final correlation computation
        cc = self%crvec1(ithr)%r(1:self%nrots) / &
            &sqrt(cc * real(self%sqsums_ptcls(i), sp) * real(2*self%nrots, sp))
    end subroutine gen_corrs
```

`gen_euclid_crvec`:

```fortran
    subroutine gen_euclid_crvec( self, iref, iptcl, shift, norm, ithr )
        class(polarft_calc), target, intent(inout) :: self
        integer,                     intent(in)    :: iref, iptcl
        real(sp),                    intent(in)    :: shift(2)
        real(dp),                    intent(out)   :: norm
        integer,                     intent(out)   :: ithr
        complex(sp), pointer :: shmat(:,:)
        complex(sp) :: c, c2
        real(dp)    :: wk, wk2
        real(sp)    :: shift_mag_sq
        integer     :: k, i, kk, k0, p, ieo, ih
        ithr         = omp_get_thread_num() + 1
        i            = self%pinds(iptcl)
        ih = merge(1, 2, self%iseven(i))
        k0           = self%kfromto(1)
        ieo          = merge(REF_EVEN, REF_ODD, self%iseven(i))
        norm         = self%wsqsum_ptcl_now(iptcl) * real(2*self%nrots, dp)
        shmat        => self%heap_vars(ithr)%shmat
        shift_mag_sq = shift(1)*shift(1) + shift(2)*shift(2)
        if( shift_mag_sq > SHERRSQ )then
            call self%gen_shmat4aln(ithr, shift, shmat)
            do k = self%kfromto(1), self%kfromto(2)
                kk = k - k0 + 1
                self%cmat2_many(ithr)%c(1:self%pftsz,kk) = shmat(:,k) * self%pfts_refs(:,k,iref,ieo)
                self%cmat2_many(ithr)%c(self%pftsz+1:self%nrots,kk) = &
                    conjg(self%cmat2_many(ithr)%c(1:self%pftsz,kk))
            enddo
            call fftwf_execute_dft(self%plan_fwd1_many, self%cmat2_many(ithr)%c, self%cmat2_many(ithr)%c)
            self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
            ! two shells per pass, each rounded to single precision in shell order
            do k = self%kfromto(1), self%kfromto(2)-1, 2
                kk = k - k0 + 1
                wk = self%sigma2_wk8(k,ih)
                wk2 = self%sigma2_wk8(k+1,ih)
                do p = 1,self%pftsz
                    c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                    c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%cmat2_many(ithr)%c(p,kk))
                    c2 = self%ft_ctf2(p,k+1,i) * self%ft_ref2(p,k+1,iref,ieo)
                    c2 = c2 - 2.0 * self%ft_ptcl_ctf(p,k+1,i) * conjg(self%cmat2_many(ithr)%c(p,kk+1))
                    self%crvec1(ithr)%c(p) = cmplx(cmplx(self%crvec1(ithr)%c(p) + wk * c, kind=sp) + wk2 * c2, &
                        &kind=kind(self%crvec1(ithr)%c))
                enddo
            enddo
            if( mod(self%nk,2) == 1 )then
                k  = self%kfromto(2)
                kk = k - k0 + 1
                wk = self%sigma2_wk8(k,ih)
                do p = 1,self%pftsz
                    c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                    c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%cmat2_many(ithr)%c(p,kk))
                    self%crvec1(ithr)%c(p) = cmplx(self%crvec1(ithr)%c(p) + wk * c, &
                        &kind=kind(self%crvec1(ithr)%c))
                enddo
            endif
        else
            self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
            ! two shells per pass, each rounded to single precision in shell order
            do k = self%kfromto(1), self%kfromto(2)-1, 2
                wk  = self%sigma2_wk8(k,ih)
                wk2 = self%sigma2_wk8(k+1,ih)
                do p = 1,self%pftsz
                    c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                    c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%ft_ref(p,k,iref,ieo))
                    c2 = self%ft_ctf2(p,k+1,i) * self%ft_ref2(p,k+1,iref,ieo)
                    c2 = c2 - 2.0 * self%ft_ptcl_ctf(p,k+1,i) * conjg(self%ft_ref(p,k+1,iref,ieo))
                    self%crvec1(ithr)%c(p) = cmplx(cmplx(self%crvec1(ithr)%c(p) + wk * c, kind=sp) + wk2 * c2, &
                        &kind=kind(self%crvec1(ithr)%c))
                enddo
            enddo
            if( mod(self%nk,2) == 1 )then
                k  = self%kfromto(2)
                wk = self%sigma2_wk8(k,ih)
                do p = 1,self%pftsz
                    c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                    c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%ft_ref(p,k,iref,ieo))
                    self%crvec1(ithr)%c(p) = cmplx(self%crvec1(ithr)%c(p) + wk * c, &
                        &kind=kind(self%crvec1(ithr)%c))
                enddo
            endif
        endif
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
    end subroutine gen_euclid_crvec
```

- `euclid_sums`: deleted

`gen_euclid_crvec_bank`:

```fortran
    subroutine gen_euclid_crvec_bank( self, iref, iptcl, bank, norm, ithr )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        complex(sp),         intent(in)    :: bank(:,:)
        real(dp),            intent(out)   :: norm
        integer,             intent(out)   :: ithr
        complex(sp) :: c, c2
        real(dp)    :: wk, wk2
        integer     :: k, i, kk, k0, p, ieo, ih
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        ih   = merge(1, 2, self%iseven(i))
        k0   = self%kfromto(1)
        ieo  = merge(REF_EVEN, REF_ODD, self%iseven(i))
        norm = self%wsqsum_ptcl_now(iptcl) * real(2*self%nrots, dp)
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        ! two shells per pass, each rounded to single precision in shell order
        do k = self%kfromto(1), self%kfromto(2)-1, 2
            kk  = k - k0 + 1
            wk  = self%sigma2_wk8(k,ih)
            wk2 = self%sigma2_wk8(k+1,ih)
            do p = 1,self%pftsz+1
                c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                c = c - 2.0 * bank(p,kk) * conjg(self%ft_ref(p,k,iref,ieo))
                c2 = self%ft_ctf2(p,k+1,i) * self%ft_ref2(p,k+1,iref,ieo)
                c2 = c2 - 2.0 * bank(p,kk+1) * conjg(self%ft_ref(p,k+1,iref,ieo))
                self%crvec1(ithr)%c(p) = cmplx(cmplx(self%crvec1(ithr)%c(p) + wk * c, kind=sp) + wk2 * c2, &
                    &kind=kind(self%crvec1(ithr)%c))
            enddo
        enddo
        if( mod(self%nk,2) == 1 )then
            k  = self%kfromto(2)
            kk = k - k0 + 1
            wk = self%sigma2_wk8(k,ih)
            do p = 1,self%pftsz+1
                c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                c = c - 2.0 * bank(p,kk) * conjg(self%ft_ref(p,k,iref,ieo))
                self%crvec1(ithr)%c(p) = cmplx(self%crvec1(ithr)%c(p) + wk * c, &
                    &kind=kind(self%crvec1(ithr)%c))
            enddo
        endif
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
    end subroutine gen_euclid_crvec_bank
```

`gen_corrs_bank`:

```fortran
    subroutine gen_corrs_bank( self, iref, iptcl, bank, cc )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        complex(sp),         intent(in)    :: bank(:,:)
        real(sp),            intent(out)   :: cc(self%nrots)
        integer :: i, ithr, k, kk, k0, ieo
        ithr = omp_get_thread_num() + 1
        i    = self%pinds(iptcl)
        k0   = self%kfromto(1)
        ieo  = merge(REF_EVEN, REF_ODD, self%iseven(i))
        ! IFFT( sum_k FT(CTF2) x FT(REF2) )
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        do k = self%kfromto(1), self%kfromto(2)-1, 2
            self%crvec1(ithr)%c = (self%crvec1(ithr)%c + self%ft_ctf2(:,k,i) * self%ft_ref2(:,k,iref,ieo)) +&
                &self%ft_ctf2(:,k+1,i) * self%ft_ref2(:,k+1,iref,ieo)
        enddo
        if( mod(self%nk,2) == 1 ) self%crvec1(ithr)%c = self%crvec1(ithr)%c +&
            &self%ft_ctf2(:,self%kfromto(2),i) * self%ft_ref2(:,self%kfromto(2),iref,ieo)
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots)
        ! IFFT( sum_k FT(S.X.CTF) x FT(REF)* )
        self%crvec1(ithr)%c = cmplx(0.,0.,kind=c_float_complex)
        do k = self%kfromto(1), self%kfromto(2)-1, 2
            kk = k - k0 + 1
            self%crvec1(ithr)%c = (self%crvec1(ithr)%c + bank(:,kk) * conjg(self%ft_ref(:,k,iref,ieo))) +&
                &bank(:,kk+1) * conjg(self%ft_ref(:,k+1,iref,ieo))
        end do
        if( mod(self%nk,2) == 1 ) self%crvec1(ithr)%c = self%crvec1(ithr)%c +&
            &bank(:,self%nk) * conjg(self%ft_ref(:,self%kfromto(2),iref,ieo))
        call fftwf_execute_dft_c2r(self%plan_bwd1_single, self%crvec1(ithr)%c, self%crvec1(ithr)%r)
        cc = self%crvec1(ithr)%r(1:self%nrots) / &
            &sqrt(cc * real(self%sqsums_ptcls(i), sp) * real(2*self%nrots, sp))
    end subroutine gen_corrs_bank
```

- `cc_sums`: deleted

### c8: c6 with the cc sums as explicit bin loops (compiled only: cc still jammed)

- `euclid_sums`: as c6

`cc_sums`:

```fortran
    !> sum_k a_k x b_k (l_conj: a_k x b_k*) over the n1 angular bins into acc, plain array dummies as in euclid_sums,
    !! two shells per pass, summed in shell order
    subroutine cc_sums( n1, ldb, nk, a, b, l_conj, acc )
        integer,     intent(in)    :: n1, ldb, nk
        complex(sp), intent(in)    :: a(n1,nk), b(ldb,nk)
        logical,     intent(in)    :: l_conj
        complex(sp), intent(inout) :: acc(n1)
        integer :: kk, p
        acc = cmplx(0.,0.,kind=sp)
        if( l_conj )then
            do kk = 1,nk-1,2
                do p = 1,n1
                    acc(p) = (acc(p) + a(p,kk) * conjg(b(p,kk))) + a(p,kk+1) * conjg(b(p,kk+1))
                enddo
            enddo
            if( mod(nk,2) == 1 )then
                do p = 1,n1
                    acc(p) = acc(p) + a(p,nk) * conjg(b(p,nk))
                enddo
            endif
        else
            do kk = 1,nk-1,2
                do p = 1,n1
                    acc(p) = (acc(p) + a(p,kk) * b(p,kk)) + a(p,kk+1) * b(p,kk+1)
                enddo
            enddo
            if( mod(nk,2) == 1 )then
                do p = 1,n1
                    acc(p) = acc(p) + a(p,nk) * b(p,nk)
                enddo
            endif
        endif
    end subroutine cc_sums
```

- `gen_corrs`, `gen_euclid_crvec`, `gen_euclid_crvec_bank`, `gen_corrs_bank`: as c2 (HEAD)

### c9: one helper call per shell, euclid and cc

`euclid_sums` (with the helpers it calls):

```fortran
    !> The euclid sums of one pair over angular bins 1..np, sum_k w_k (FT(CTF2) x FT(REF2) - 2 FT(X.CTF) x FT(REF)*),
    !! into acc, one shell at a time (euclid_shell_sum)
    subroutine euclid_sums( np, n1, ldr, nk, ctf2, ref2, xctf, ref, wk, acc )
        integer,     intent(in)    :: np, n1, ldr, nk
        complex(sp), intent(in)    :: ctf2(n1,nk), ref2(n1,nk), xctf(n1,nk), ref(ldr,nk)
        real(dp),    intent(in)    :: wk(nk)
        complex(sp), intent(inout) :: acc(n1)
        integer :: kk
        acc = cmplx(0.,0.,kind=sp)
        do kk = 1,nk
            call euclid_shell_sum(np, ctf2(:,kk), ref2(:,kk), xctf(:,kk), ref(:,kk), wk(kk), acc)
        enddo
    end subroutine euclid_sums

    !> One shell of euclid_sums: plain array dummies, unit stride and no overlap
    subroutine euclid_shell_sum( np, ctf2, ref2, xctf, ref, wk, acc )
        integer,     intent(in)    :: np
        complex(sp), intent(in)    :: ctf2(np), ref2(np), xctf(np), ref(np)
        real(dp),    intent(in)    :: wk
        complex(sp), intent(inout) :: acc(np)
        complex(sp) :: c
        integer     :: p
        do p = 1,np
            c = ctf2(p) * ref2(p)
            c = c - 2.0 * xctf(p) * conjg(ref(p))
            acc(p) = cmplx(acc(p) + wk * c, kind=sp)
        enddo
    end subroutine euclid_shell_sum
```

`cc_sums` (with the helpers it calls):

```fortran
    !> sum_k a_k x b_k (l_conj: a_k x b_k*) over the n1 angular bins into acc, one shell at a time (cc_shell_sum)
    subroutine cc_sums( n1, ldb, nk, a, b, l_conj, acc )
        integer,     intent(in)    :: n1, ldb, nk
        complex(sp), intent(in)    :: a(n1,nk), b(ldb,nk)
        logical,     intent(in)    :: l_conj
        complex(sp), intent(inout) :: acc(n1)
        integer :: kk
        acc = cmplx(0.,0.,kind=sp)
        do kk = 1,nk
            call cc_shell_sum(n1, a(:,kk), b(:,kk), l_conj, acc)
        enddo
    end subroutine cc_sums

    !> One shell of cc_sums: plain array dummies, unit stride and no overlap
    subroutine cc_shell_sum( n1, a, b, l_conj, acc )
        integer,     intent(in)    :: n1
        complex(sp), intent(in)    :: a(n1), b(n1)
        logical,     intent(in)    :: l_conj
        complex(sp), intent(inout) :: acc(n1)
        if( l_conj )then
            acc = acc + a * conjg(b)
        else
            acc = acc + a * b
        endif
    end subroutine cc_shell_sum
```

- `gen_corrs`, `gen_euclid_crvec`, `gen_euclid_crvec_bank`, `gen_corrs_bank`: as c2 (HEAD)

### c10: HEAD with `!$omp simd` on the bin loops (compiled only: jam not prevented)

`euclid_sums`:

```fortran
    !> The euclid sums of one pair over angular bins 1..np, sum_k w_k (FT(CTF2) x FT(REF2) - 2 FT(X.CTF) x FT(REF)*),
    !! into acc. Plain array dummies: the compiler then knows unit stride and no overlap and vectorizes the bins.
    subroutine euclid_sums( np, n1, ldr, nk, ctf2, ref2, xctf, ref, wk, acc )
        integer,     intent(in)    :: np, n1, ldr, nk
        complex(sp), intent(in)    :: ctf2(n1,nk), ref2(n1,nk), xctf(n1,nk), ref(ldr,nk)
        real(dp),    intent(in)    :: wk(nk)
        complex(sp), intent(inout) :: acc(n1)
        complex(sp) :: c
        integer     :: kk, p
        acc = cmplx(0.,0.,kind=sp)
        do kk = 1,nk
            !$omp simd private(c)
            do p = 1,np
                c = ctf2(p,kk) * ref2(p,kk)
                c = c - 2.0 * xctf(p,kk) * conjg(ref(p,kk))
                acc(p) = cmplx(acc(p) + wk(kk) * c, kind=sp)
            enddo
        enddo
    end subroutine euclid_sums
```

`cc_sums`:

```fortran
    !> sum_k a_k x b_k (l_conj: a_k x b_k*) over the n1 angular bins into acc, plain array dummies as in euclid_sums
    subroutine cc_sums( n1, ldb, nk, a, b, l_conj, acc )
        integer,     intent(in)    :: n1, ldb, nk
        complex(sp), intent(in)    :: a(n1,nk), b(ldb,nk)
        logical,     intent(in)    :: l_conj
        complex(sp), intent(inout) :: acc(n1)
        integer :: kk, p
        acc = cmplx(0.,0.,kind=sp)
        if( l_conj )then
            do kk = 1,nk
                !$omp simd
                do p = 1,n1
                    acc(p) = acc(p) + a(p,kk) * conjg(b(p,kk))
                enddo
            enddo
        else
            do kk = 1,nk
                !$omp simd
                do p = 1,n1
                    acc(p) = acc(p) + a(p,kk) * b(p,kk)
                enddo
            enddo
        endif
    end subroutine cc_sums
```

- `gen_corrs`, `gen_euclid_crvec`, `gen_euclid_crvec_bank`, `gen_corrs_bank`: as c2 (HEAD)

### c11: HEAD with `!GCC$ unroll 1` on the shell loops (compiled only: jam not prevented)

`euclid_sums`:

```fortran
    !> The euclid sums of one pair over angular bins 1..np, sum_k w_k (FT(CTF2) x FT(REF2) - 2 FT(X.CTF) x FT(REF)*),
    !! into acc. Plain array dummies: the compiler then knows unit stride and no overlap and vectorizes the bins.
    subroutine euclid_sums( np, n1, ldr, nk, ctf2, ref2, xctf, ref, wk, acc )
        integer,     intent(in)    :: np, n1, ldr, nk
        complex(sp), intent(in)    :: ctf2(n1,nk), ref2(n1,nk), xctf(n1,nk), ref(ldr,nk)
        real(dp),    intent(in)    :: wk(nk)
        complex(sp), intent(inout) :: acc(n1)
        complex(sp) :: c
        integer     :: kk, p
        acc = cmplx(0.,0.,kind=sp)
        !GCC$ unroll 1
        do kk = 1,nk
            do p = 1,np
                c = ctf2(p,kk) * ref2(p,kk)
                c = c - 2.0 * xctf(p,kk) * conjg(ref(p,kk))
                acc(p) = cmplx(acc(p) + wk(kk) * c, kind=sp)
            enddo
        enddo
    end subroutine euclid_sums
```

`cc_sums`:

```fortran
    !> sum_k a_k x b_k (l_conj: a_k x b_k*) over the n1 angular bins into acc, plain array dummies as in euclid_sums
    subroutine cc_sums( n1, ldb, nk, a, b, l_conj, acc )
        integer,     intent(in)    :: n1, ldb, nk
        complex(sp), intent(in)    :: a(n1,nk), b(ldb,nk)
        logical,     intent(in)    :: l_conj
        complex(sp), intent(inout) :: acc(n1)
        integer :: kk
        acc = cmplx(0.,0.,kind=sp)
        if( l_conj )then
            !GCC$ unroll 1
            do kk = 1,nk
                acc = acc + a(:,kk) * conjg(b(1:n1,kk))
            enddo
        else
            !GCC$ unroll 1
            do kk = 1,nk
                acc = acc + a(:,kk) * b(1:n1,kk)
            enddo
        endif
    end subroutine cc_sums
```

- `gen_corrs`, `gen_euclid_crvec`, `gen_euclid_crvec_bank`, `gen_corrs_bank`: as c2 (HEAD)

### c12: euclid as c6, cc as c9

- `euclid_sums`: as c6

- `cc_sums`: as c9

- `gen_corrs`, `gen_euclid_crvec`, `gen_euclid_crvec_bank`, `gen_corrs_bank`: as c2 (HEAD)

### c13: euclid as c6, cc two shells per helper call

- `euclid_sums`: as c6

`cc_sums` (with the helpers it calls):

```fortran
    !> sum_k a_k x b_k (l_conj: a_k x b_k*) over the n1 angular bins into acc, two shells per call of cc_shells_sum
    subroutine cc_sums( n1, ldb, nk, a, b, l_conj, acc )
        integer,     intent(in)    :: n1, ldb, nk
        complex(sp), intent(in)    :: a(n1,nk), b(ldb,nk)
        logical,     intent(in)    :: l_conj
        complex(sp), intent(inout) :: acc(n1)
        integer :: kk
        acc = cmplx(0.,0.,kind=sp)
        do kk = 1,nk-1,2
            call cc_shells_sum(n1, ldb, 2, a(:,kk:kk+1), b(:,kk:kk+1), l_conj, acc)
        enddo
        if( mod(nk,2) == 1 ) call cc_shells_sum(n1, ldb, 1, a(:,nk:nk), b(:,nk:nk), l_conj, acc)
    end subroutine cc_sums

    !> One or two shells of cc_sums in one pass over acc, summed in shell order; plain array dummies
    subroutine cc_shells_sum( n1, ldb, ns, a, b, l_conj, acc )
        integer,     intent(in)    :: n1, ldb, ns
        complex(sp), intent(in)    :: a(n1,ns), b(ldb,ns)
        logical,     intent(in)    :: l_conj
        complex(sp), intent(inout) :: acc(n1)
        integer :: p
        if( ns == 2 )then
            if( l_conj )then
                do p = 1,n1
                    acc(p) = (acc(p) + a(p,1) * conjg(b(p,1))) + a(p,2) * conjg(b(p,2))
                enddo
            else
                do p = 1,n1
                    acc(p) = (acc(p) + a(p,1) * b(p,1)) + a(p,2) * b(p,2)
                enddo
            endif
        else
            if( l_conj )then
                acc = acc + a(:,1) * conjg(b(1:n1,1))
            else
                acc = acc + a(:,1) * b(1:n1,1)
            endif
        endif
    end subroutine cc_shells_sum
```

- `gen_corrs`, `gen_euclid_crvec`, `gen_euclid_crvec_bank`, `gen_corrs_bank`: as c2 (HEAD)

### c14: euclid as c9, cc as c13: the form kept

- `euclid_sums`: as c9

- `cc_sums`: as c13

- `gen_corrs`, `gen_euclid_crvec`, `gen_euclid_crvec_bank`, `gen_corrs_bank`: as c2 (HEAD)
