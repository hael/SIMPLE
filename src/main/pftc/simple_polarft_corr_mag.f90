submodule (simple_polarft_calc) simple_polarft_corr_mag
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
#include "simple_local_flags.inc"
implicit none

contains

    module real function calc_magcorr_rot(self, iref, iptcl, irot, kweight)
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl, irot
        logical, optional,   intent(in)    :: kweight
        complex(dp), pointer :: pft_ref(:,:), pft_rot_ref(:,:)
        real(dp),    pointer :: mag_rot_ref(:,:)
        real(dp) :: sqsumref, sqsumptcl, num
        integer  :: i, k, ithr
        logical  :: kw
        kw = .true.
        if( present(kweight) ) kw = kweight
        i    = self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        calc_magcorr_rot = 0.
        pft_ref     => self%heap_vars(ithr)%pft_ref_8
        pft_rot_ref => self%heap_vars(ithr)%pft_ref_tmp_8
        mag_rot_ref => self%heap_vars(ithr)%pft_r1_8
        if( self%iseven(i) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        call self%rotate_pft(pft_ref, irot, pft_rot_ref)
        pft_rot_ref = pft_rot_ref * self%ctfmats(:,:,i)
        select case(params_glob%cc_objfun)
        case(OBJFUN_CC)
            mag_rot_ref = abs(pft_rot_ref)
            sqsumref    = 0.d0
            sqsumptcl   = 0.d0
            num         = 0.d0
            do k = self%kfromto(1),self%kfromto(2)
                if( kw )then
                    sqsumptcl = sqsumptcl + real(k,dp) * sum(real(csq_fast(self%pfts_ptcls(:,k,i)),dp))
                    sqsumref  = sqsumref  + real(k,dp) * sum(mag_rot_ref(:,k)**2)
                    num       = num       + real(k,dp) * sum(mag_rot_ref(:,k) * real(abs(self%pfts_ptcls(:,k,i)),dp))
                else
                    sqsumptcl = sqsumptcl + sum(real(csq_fast(self%pfts_ptcls(:,k,i)),dp))
                    sqsumref  = sqsumref  + sum(mag_rot_ref(:,k)**2)
                    num       = num       + sum(mag_rot_ref(:,k) * real(abs(self%pfts_ptcls(:,k,i)),dp))
                endif
            enddo
            calc_magcorr_rot = real(num/sqrt(sqsumref*sqsumptcl))
        case(OBJFUN_EUCLID)
            ! not implemented
        end select
    end function calc_magcorr_rot

    module subroutine gen_corrs_mag_cc( self, iref, iptcl, ccs, kweight )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        real,                intent(inout) :: ccs(self%pftsz)
        logical,   optional, intent(in)    :: kweight
        complex(dp), pointer :: pft_ref(:,:)
        real(dp),    pointer :: pft_mag_ptcl(:,:)
        real(dp)    :: sumsqptcl
        integer     :: i, k, ithr
        logical     :: kw
        kw = .true.
        if( present(kweight) ) kw = kweight
        i            = self%pinds(iptcl)
        ithr         = omp_get_thread_num() + 1
        pft_ref      => self%heap_vars(ithr)%pft_ref_8
        pft_mag_ptcl => self%heap_vars(ithr)%pft_r1_8
        if( self%iseven(i) )then
            pft_ref = self%pfts_refs_even(:,:,iref)
        else
            pft_ref = self%pfts_refs_odd(:,:,iref)
        endif
        pft_mag_ptcl(:,:)           = real(abs(self%pfts_ptcls(:,:,i)),dp)
        sumsqptcl                   = 0.d0
        self%heap_vars(ithr)%kcorrs = 0.d0
        self%drvec(ithr)%r          = 0.d0
        do k = self%kfromto(1),self%kfromto(2)
            ! |X|2
            if( kw )then
                sumsqptcl = sumsqptcl + real(k,dp) * sum(pft_mag_ptcl(:,k)**2)
            else
                sumsqptcl = sumsqptcl + sum(pft_mag_ptcl(:,k)**2)
            endif
            ! FT(CTF2)
            self%rvec1(ithr)%r(1:self%pftsz) = self%ctfmats(:,k,i)**2
            self%rvec1(ithr)%r(self%pftsz+1:self%nrots) = self%rvec1(ithr)%r(1:self%pftsz)
            call fftwf_execute_dft_r2c(self%plan_mem_r2c, self%rvec1(ithr)%r, self%cvec1(ithr)%c)
            self%cvec2(ithr)%c(1:self%pftsz+1) = self%cvec1(ithr)%c
            ! FT(|REF|2)
            self%rvec1(ithr)%r(1:self%pftsz)            = real(pft_ref(:,k)*conjg(pft_ref(:,k)),sp)
            self%rvec1(ithr)%r(self%pftsz+1:self%nrots) = self%rvec1(ithr)%r(1:self%pftsz)
            call fftwf_execute_dft_r2c(self%plan_mem_r2c, self%rvec1(ithr)%r, self%cvec1(ithr)%c)
            ! FT(CTF2) x FT(|REF|2)*
            self%cvec1(ithr)%c = self%cvec2(ithr)%c(1:self%pftsz+1) * conjg(self%cvec1(ithr)%c)
            ! IFFT( FT(CTF2) x FT(|REF|2)* )
            call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
            if( kw )then
                self%drvec(ithr)%r(1:self%pftsz) = self%drvec(ithr)%r(1:self%pftsz) + real(k,dp) * real(self%rvec1(ithr)%r(1:self%pftsz),dp)
            else
                self%drvec(ithr)%r(1:self%pftsz) = self%drvec(ithr)%r(1:self%pftsz) + real(self%rvec1(ithr)%r(1:self%pftsz),dp)
            endif
            ! FT(|REF|)*
            self%rvec1(ithr)%r(1:self%pftsz) = real(abs(pft_ref(:,k)), kind=sp)
            self%rvec1(ithr)%r(self%pftsz+1:self%nrots) = self%rvec1(ithr)%r(1:self%pftsz)
            call fftwf_execute_dft_r2c(self%plan_mem_r2c, self%rvec1(ithr)%r, self%cvec1(ithr)%c)
            ! FT(|X|.CTF) x FT(|REF|)*
            self%cvec1(ithr)%c = self%ft_absptcl_ctf(:,k,i) * conjg(self%cvec1(ithr)%c)
            ! IFFT( FT(|X|.CTF) x FT(|REF|)* )
            call fftwf_execute_dft_c2r(self%plan_bwd1, self%cvec1(ithr)%c, self%rvec1(ithr)%r)
            if( kw )then
                self%heap_vars(ithr)%kcorrs(1:self%pftsz) = self%heap_vars(ithr)%kcorrs(1:self%pftsz) +&
                    &real(k,dp) * real(self%rvec1(ithr)%r(1:self%pftsz),dp)
            else
                self%heap_vars(ithr)%kcorrs(1:self%pftsz) = self%heap_vars(ithr)%kcorrs(1:self%pftsz) +&
                    &real(self%rvec1(ithr)%r(1:self%pftsz),dp)
            endif
        end do
        self%drvec(ithr)%r(1:self%pftsz) = self%drvec(ithr)%r(1:self%pftsz) * (sumsqptcl * real(2*self%nrots,dp))
        ccs = real(self%heap_vars(ithr)%kcorrs(1:self%pftsz) / dsqrt(self%drvec(ithr)%r(1:self%pftsz)))
    end subroutine gen_corrs_mag_cc

    module subroutine gen_corrs_mag( self, iref, iptcl, ccs, kweight )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl
        real,                intent(inout) :: ccs(self%pftsz)
        logical,   optional, intent(in)    :: kweight
        select case(params_glob%cc_objfun)
        case(OBJFUN_CC)
            call self%gen_corrs_mag_cc( iref, iptcl, ccs, kweight )
        case(OBJFUN_EUCLID)
            ! unsupported
            ccs = 0.
        end select
    end subroutine gen_corrs_mag

    module subroutine bidirectional_shift_search( self, iref, iptcl, irot, hn, shifts, grid1, grid2 )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: iref, iptcl, irot, hn
        real,                intent(in)    :: shifts(-hn:hn)
        real,                intent(out)   :: grid1(-hn:hn,-hn:hn), grid2(-hn:hn,-hn:hn)
        complex(dp), pointer :: pft_ref_8(:,:), diff_8(:,:), pft_shref_8(:,:)
        complex(sp), pointer :: pft_ptcl(:,:)
        real(sp),    pointer :: rctf(:,:)
        real(dp)             :: shvec(2), score, sqsum_ref, denom
        integer              :: ithr, i, ix, iy, k, prot
        i    =  self%pinds(iptcl)
        ithr = omp_get_thread_num() + 1
        pft_ptcl    => self%heap_vars(ithr)%pft_ref
        pft_ref_8   => self%heap_vars(ithr)%pft_ref_8
        diff_8      => self%heap_vars(ithr)%pft_ref_tmp_8
        pft_shref_8 => self%heap_vars(ithr)%shmat_8
        rctf        => self%heap_vars(ithr)%pft_r
        if( self%iseven(i) )then
            pft_ref_8 = self%pfts_refs_even(:,:,iref)
        else
            pft_ref_8 = self%pfts_refs_odd(:,:,iref)
        endif
        ! Rotate particle
        prot = self%nrots-irot+2
        if( prot > self%nrots ) prot = prot-self%nrots
        call self%rotate_pft(self%pfts_ptcls(:,:,i), prot, pft_ptcl(:,:))
        ! Reference CTF modulation
        call self%rotate_pft(self%ctfmats(:,:,i), prot, rctf)
        pft_ref_8 = pft_ref_8 * real(rctf,dp)
        select case(params_glob%cc_objfun)
        case(OBJFUN_CC)
            do ix = -hn,hn
                do iy = -hn,hn
                    shvec = real([shifts(ix), shifts(iy)],dp)
                    call self%gen_shmat_8(ithr, shvec, pft_shref_8) ! shift matrix
                    pft_shref_8 = pft_ref_8 * pft_shref_8           ! shifted reference
                    ! first orientation
                    sqsum_ref = 0.d0
                    score     = 0.d0
                    do k = self%kfromto(1),self%kfromto(2)
                        sqsum_ref = sqsum_ref + real(k,kind=dp) * sum(real(pft_shref_8(:,k) * conjg(pft_shref_8(:,k)),dp))
                        score     = score     + real(k,kind=dp) * sum(real(pft_shref_8(:,k) * conjg(pft_ptcl(:,k)),dp))
                    end do
                    denom        = dsqrt(sqsum_ref * self%ksqsums_ptcls(i))
                    grid1(ix,iy) = real(score / denom)
                    ! second orientation (first+pi)
                    score = 0.d0
                    do k = self%kfromto(1),self%kfromto(2)
                        score = score + real(k,kind=dp) * sum(real(pft_shref_8(:,k) * pft_ptcl(:,k),dp))
                    end do
                    grid2(ix,iy) = real(score / denom)
                enddo
            enddo
        case(OBJFUN_EUCLID)
            do ix = -hn,hn
                do iy = -hn,hn
                    shvec = real([shifts(ix), shifts(iy)],dp)
                    call self%gen_shmat_8(ithr, shvec, pft_shref_8) ! shift matrix
                    pft_shref_8 = pft_ref_8 * pft_shref_8           ! shifted reference
                    ! first orientation
                    diff_8 = pft_shref_8 - dcmplx(pft_ptcl)
                    score = 0.d0
                    do k = self%kfromto(1),self%kfromto(2)
                        score = score + (real(k,dp) / self%sigma2_noise(k,iptcl)) * sum(real(diff_8(:,k)*conjg(diff_8(:,k)),dp))
                    end do
                    grid1(ix,iy) = real(exp( -score / self%wsqsums_ptcls(i) ))
                    ! second orientation (first+pi)
                    diff_8 = pft_shref_8 - dcmplx(conjg(pft_ptcl))
                    score = 0.d0
                    do k = self%kfromto(1),self%kfromto(2)
                        score = score + (real(k,dp) / self%sigma2_noise(k,iptcl)) * sum(real(diff_8(:,k)*conjg(diff_8(:,k)),dp))
                    end do
                    grid2(ix,iy) = real(exp( -score / self%wsqsums_ptcls(i) ))
                enddo
            enddo
        end select
    end subroutine bidirectional_shift_search

end submodule simple_polarft_corr_mag
