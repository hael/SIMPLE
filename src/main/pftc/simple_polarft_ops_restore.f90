!@descr: submodule for class average restoration in the polar Fourier domain
submodule (simple_polarft_calc) simple_polarft_ops_restore
implicit none
#include "simple_local_flags.inc"
contains               

    ! 2D restoration routines

    !>  \brief  Restores 2D class-averages only
    !>  Calculates and writes FRCs and uses regularized wiener restoration
    module subroutine polar_cavger_merge_eos_and_norm2D( self, clsfrcs, fname )
        class(polarft_calc), intent(inout) :: self
        class(class_frcs),   intent(inout) :: clsfrcs
        class(string),       intent(in)    :: fname
        complex(dp) :: pft(self%pftsz,self%kfromto(1):self%interpklim)
        complex(dp) :: even(self%pftsz,self%kfromto(1):self%interpklim)
        complex(dp) :: odd(self%pftsz,self%kfromto(1):self%interpklim)
        real(dp)    :: ctf2(self%pftsz,self%kfromto(1):self%interpklim)
        real        :: frc(fdim(self%p_ptr%box_crop)-1)
        integer     :: icls, eo_pop(2), pop, find, filtsz
        select case(trim(self%p_ptr%ref_type))
            case('polar_cavg')
                ! all good
            case DEFAULT
                THROW_HARD('polar_cavger_merge_eos_and_norm2D only for 2D cavgs restoration')
        end select
        ! In case nspace/self%ncls has changed OR volume/frcs were downsampled
        filtsz = fdim(self%p_ptr%box_crop)-1
        if( (clsfrcs%get_ncls() /= self%ncls) .or. (clsfrcs%get_filtsz() /= filtsz) )then
            call clsfrcs%new(self%ncls, self%p_ptr%box_crop, self%p_ptr%smpd_crop, self%p_ptr%nstates)
        endif
        !$omp parallel do default(shared) schedule(static) proc_bind(close)&
        !$omp private(icls,eo_pop,pop,pft,ctf2,even,odd,frc,find)
        do icls = 1,self%ncls
            self%pfts_merg(:,:,icls) = DCMPLX_ZERO
            frc    = 0.0
            eo_pop = self%prev_eo_pops(:,icls) + self%eo_pops(:,icls)
            pop    = sum(eo_pop)
            if(pop == 0)then
                ! Empty class
                self%pfts_even(:,:,icls) = DCMPLX_ZERO
                self%pfts_odd(:,:,icls)  = DCMPLX_ZERO
                self%ctf2_even(:,:,icls) = 0.d0
                self%ctf2_odd(:,:,icls)  = 0.d0
            else
                ! e/o Normalization
                if( eo_pop(1) > 1) call self%safe_norm(self%pfts_even(:,:,icls), self%ctf2_even(:,:,icls), even)
                if( eo_pop(2) > 1) call self%safe_norm(self%pfts_odd(:,:,icls),  self%ctf2_odd(:,:,icls),  odd)
                ! FRC
                call self%polar_cavger_calc_frc(even, odd, filtsz, frc)
                ! Regularization
                if( self%p_ptr%l_ml_reg )then                    
                    if( pop > 1 )then
                        call add_invtausq2rho(frc(self%kfromto(1):self%interpklim), self%ctf2_even(:,:,icls), self%ctf2_odd(:,:,icls))
                    endif
                    ! e/o re-normalization
                    if( eo_pop(1) > 1) call self%safe_norm(self%pfts_even(:,:,icls), self%ctf2_even(:,:,icls), even)
                    if( eo_pop(2) > 1) call self%safe_norm(self%pfts_odd(:,:,icls),  self%ctf2_odd(:,:,icls),  odd)
                endif
                ! merged class
                if(pop > 1)then
                    pft  = self%pfts_even(:,:,icls) + self%pfts_odd(:,:,icls)
                    ctf2 = self%ctf2_even(:,:,icls) + self%ctf2_odd(:,:,icls)
                    call self%safe_norm(pft, ctf2, self%pfts_merg(:,:,icls))
                endif
                ! average low-resolution info between eo pairs to keep things in register
                find = min(self%interpklim, clsfrcs%estimate_find_for_eoavg(icls, 1))
                if( find >= self%kfromto(1) )then
                    even(:,self%kfromto(1):find) = self%pfts_merg(:,self%kfromto(1):find,icls)
                    odd(:,self%kfromto(1):find)  = self%pfts_merg(:,self%kfromto(1):find,icls)
                endif
                ! updates arrays
                self%pfts_even(:,:,icls) = even
                self%pfts_odd(:,:,icls)  = odd
            endif
            ! store FRC
            call clsfrcs%set_frc(icls, frc, 1)
        end do
        !$omp end parallel do
        ! Write FRCs
        call clsfrcs%write(fname)
        contains

            subroutine add_invtausq2rho( frc, ctf2e, ctf2o )
                real,     intent(in)    :: frc(self%kfromto(1):self%interpklim)
                real(dp), intent(inout) :: ctf2e(self%pftsz,self%kfromto(1):self%interpklim)
                real(dp), intent(inout) :: ctf2o(self%pftsz,self%kfromto(1):self%interpklim)
                real(dp) :: sig2e(self%kfromto(1):self%interpklim), sig2o(self%kfromto(1):self%interpklim)
                real(dp) :: ssnr(self%kfromto(1):self%interpklim), tau2e(self%kfromto(1):self%interpklim)
                real(dp) :: tau2o(self%kfromto(1):self%interpklim)
                real(dp) :: cc, fudge
                integer  :: k, kstart, p
                sig2e = 0.d0
                sig2o = 0.d0
                do k = self%kfromto(1),self%interpklim
                    sig2e(k) = sig2e(k) + sum(ctf2e(:,k))
                    sig2o(k) = sig2o(k) + sum(ctf2o(:,k))
                enddo
                ! SSNR
                fudge = real(self%p_ptr%tau,dp)
                do k = self%kfromto(1),self%interpklim
                    cc      = max(0.001d0, min(0.999d0, frc(k)))
                    ssnr(k) = fudge * cc / (1.d0 - cc)
                enddo
                where( sig2e > DTINY )
                    sig2e = real(self%ncls*self%pftsz,dp) / sig2e
                elsewhere
                    sig2e = 0.d0
                end where
                tau2e = ssnr * sig2e
                where( sig2o > DTINY )
                    sig2o = real(self%ncls*self%pftsz,dp) / sig2o
                elsewhere
                    sig2o = 0.d0
                end where
                tau2o = ssnr * sig2o
                ! Add Tau2 inverse to denominators
                ! because signal assumed infinite at very low resolution there is no addition
                kstart = max(6, calc_fourier_index(self%p_ptr%hp, self%p_ptr%box_crop, self%p_ptr%smpd_crop))
                do k = kstart,self%interpklim
                    if( tau2e(k) > DTINY )then
                        ! CTF2 <- CTF2 + avgCTF2/(tau*SSNR)
                        ctf2e(:,k) = ctf2e(:,k) + 1.d0 / (fudge * tau2e(k))
                    else
                        do p = 1,self%pftsz
                            ctf2e(p,k) = ctf2e(p,k) + min(1.d3, 1.d3 * ctf2e(p,k))
                        enddo
                    endif
                    if( tau2o(k) > DTINY )then
                        ctf2o(:,k) = ctf2o(:,k) + 1.d0 / (fudge * tau2o(k))
                    else
                        do p = 1,self%pftsz
                            ctf2o(p,k) = ctf2o(p,k) + min(1.d3, 1.d3 * ctf2o(p,k))
                        enddo
                    endif
                enddo
            end subroutine add_invtausq2rho

    end subroutine polar_cavger_merge_eos_and_norm2D

    ! 3D restoration routines

    module subroutine polar_cavger_merge_eos_and_norm( self, reforis, symop, cline, update_frac )
        class(polarft_calc), intent(inout) :: self
        type(oris),          intent(in)    :: reforis
        type(sym),           intent(in)    :: symop
        type(cmdline),       intent(in)    :: cline
        real,                intent(in)    :: update_frac
        type(class_frcs)         :: frcs
        type(string)             :: fname
        complex(dp), allocatable :: prev_even(:,:,:), prev_odd(:,:,:)
        complex(dp) :: pfts_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        complex(dp) :: pfts_odd(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp)    :: ctf2_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp)    :: ctf2_odd(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp)    :: fsc(self%kfromto(1):self%interpklim), ufrac_trec
        real        :: fsc_boxcrop(1:fdim(self%p_ptr%box_crop)-1)
        integer     :: tmp_pftsz, tmp_kfromto(2), tmp_nrefs
        integer     :: find4eoavg, i
        select case(trim(self%p_ptr%ref_type))
        case('comlin')
        case DEFAULT
            THROW_HARD('Invalid REF_TYPE='//trim(self%p_ptr%ref_type)//' in polar_cavger_merge_eos_and_norm')
        end select
        ! Mirror Fourier & CTF2 slices
        call mirror_slices( reforis )
        ! Common-lines conribution
        call self%calc_comlin_contrib(reforis, symop, pfts_even, pfts_odd, ctf2_even, ctf2_odd)
        ! e/o trailing reconstruction part 1
        if( self%p_ptr%l_trail_rec )then
            ufrac_trec = real(merge(self%p_ptr%ufrac_trec ,update_frac , cline%defined('ufrac_trec')),dp)
            ! read previous polar references to be used first for FSC calculation
            call self%read_pft_array(string(POLAR_REFS_FBODY)//'_even'//BIN_EXT, prev_even)
            call self%read_pft_array(string(POLAR_REFS_FBODY)//'_odd'//BIN_EXT,  prev_odd)
        endif
        ! Calculate and write global FSC, FRCs
        call calc_fsc
        fsc_boxcrop(                 :self%kfromto(1)) = 1.0
        fsc_boxcrop(self%kfromto(1)  :self%interpklim) = real(fsc(self%kfromto(1):self%interpklim))
        if( self%interpklim < size(fsc_boxcrop) )then
            fsc_boxcrop(self%interpklim+1:)            = 0.0
        endif
        call arr2file(fsc_boxcrop, string(FSC_FBODY//int2str_pad(1,2)//BIN_EXT))
        call frcs%new(self%ncls, self%p_ptr%box_crop, self%p_ptr%smpd_crop, 1)
        do i = 1,self%ncls
            ! FRCs are set to the FSC. to check if we are using those
            call frcs%set_frc(i, fsc_boxcrop, 1)
        enddo
        call frcs%write(string(FRCS_FILE))
        call frcs%kill
        ! ML regularization
        if( self%p_ptr%l_ml_reg ) call add_invtausq2rho
        ! Wiener normalization
        call restore_references
        ! Keeping e/o in register
        find4eoavg = max(K4EOAVGLB,  calc_fourier_index(FREQ4EOAVG3D, self%p_ptr%box, self%p_ptr%smpd))
        find4eoavg = min(find4eoavg, get_find_at_crit(fsc_boxcrop, FSC4EOAVG3D))
        find4eoavg = min(self%kfromto(2), find4eoavg)
        if( find4eoavg >= self%kfromto(1) )then
            !$omp parallel workshare proc_bind(close)
            self%pfts_even(:,self%kfromto(1):find4eoavg,:) = self%pfts_merg(:,self%kfromto(1):find4eoavg,:)
            self%pfts_odd(:, self%kfromto(1):find4eoavg,:) = self%pfts_merg(:,self%kfromto(1):find4eoavg,:)
            !$omp end parallel workshare
        endif
        ! e/o trailing reconstruction part 2
        if( self%p_ptr%l_trail_rec )then
            ! adding weighted previous refs after restoration of current refs
            !$omp parallel workshare proc_bind(close)
            self%pfts_even = ufrac_trec * self%pfts_even + (1.d0-ufrac_trec) * prev_even
            self%pfts_odd  = ufrac_trec * self%pfts_odd  + (1.d0-ufrac_trec) * prev_odd
            !$omp end parallel workshare
            fname = string(POLAR_REFS_FBODY)//BIN_EXT
            call self%get_pft_array_dims(fname, tmp_pftsz, tmp_kfromto, tmp_nrefs)
            if( tmp_nrefs /= self%ncls )then
                !$omp parallel workshare proc_bind(close)
                self%pfts_merg = ufrac_trec * self%pfts_merg + (1.d0-ufrac_trec) * 0.5d0 * (prev_even + prev_odd)
                !$omp end parallel workshare
                deallocate(prev_even,prev_odd)
            else
                deallocate(prev_even,prev_odd)
                call self%read_pft_array(fname, prev_even)
                !$omp parallel workshare proc_bind(close)
                self%pfts_merg = ufrac_trec * self%pfts_merg + (1.d0-ufrac_trec) * prev_even
                !$omp end parallel workshare
                deallocate(prev_even)
            endif
        endif
    contains

        ! Calculate global FSC within [self%kfromto(1);self%interpklim]
        subroutine calc_fsc
            complex(dp) :: even(self%pftsz,self%kfromto(1):self%interpklim)
            complex(dp) :: odd(self%pftsz,self%kfromto(1):self%interpklim)
            complex(dp) :: pft(self%pftsz,self%kfromto(1):self%interpklim)
            real(dp)    :: vare(self%kfromto(1):self%interpklim)
            real(dp)    :: varo(self%kfromto(1):self%interpklim)
            real(dp)    :: ctf2(self%pftsz,self%kfromto(1):self%interpklim)
            integer     :: icls, k
            ! Calculates per slice contribution to global FSC/variance
            ! no need to loop over mirrored slices
            fsc  = 0.d0; vare = 0.d0; varo = 0.d0
            !$omp parallel do default(shared) schedule(static) proc_bind(close)&
            !$omp private(icls,even,odd,k,pft,ctf2) reduction(+:fsc,vare,varo)
            do icls = 1,self%ncls/2
                ! e/o restoration
                pft  = pfts_even(:,:,icls)
                ctf2 = ctf2_even(:,:,icls)
                call self%safe_norm(pft, ctf2, even)
                pft  = pfts_odd(:,:,icls)
                ctf2 = ctf2_odd(:,:,icls)
                call self%safe_norm(pft, ctf2, odd)
                if( self%p_ptr%l_trail_rec )then
                    ! adding previous reference
                    even = ufrac_trec * even + (1.d0-ufrac_trec) * prev_even(:,:,icls)
                    odd  = ufrac_trec * odd  + (1.d0-ufrac_trec) * prev_odd(:,:,icls)
                endif
                ! FSC contribution
                do k = self%kfromto(1),self%interpklim
                    fsc(k)   = fsc(k)  + sum(real(even(:,k) * conjg(odd(:,k)), dp))
                    vare(k)  = vare(k) + sum(real(even(:,k) * conjg(even(:,k)),dp))
                    varo(k)  = varo(k) + sum(real(odd(:,k)  * conjg(odd(:,k)), dp))
                enddo
            enddo
            !$omp end parallel do
            ! Variance
            vare = vare * varo
            ! FSC
            where( vare > DTINY )
                fsc = fsc / sqrt(vare)
            elsewhere
                fsc = 0.d0
            end where
        end subroutine calc_fsc

        subroutine add_invtausq2rho
            real(dp) :: sig2e(self%kfromto(1):self%interpklim), sig2o(self%kfromto(1):self%interpklim)
            real(dp) :: ssnr(self%kfromto(1):self%interpklim), tau2(self%kfromto(1):self%interpklim)
            real(dp) :: cc, fudge, invtau2
            integer  :: icls, k, kstart, p, nhalf
            nhalf = self%ncls/2
            ! Radial CTF2 sum
            sig2e = 0.d0; sig2o = 0.d0
            !$omp parallel do default(shared) schedule(static) proc_bind(close)&
            !$omp private(k) reduction(+:sig2e,sig2o)
            do k = self%kfromto(1),self%interpklim
                sig2e(k) = sig2e(k) + sum(ctf2_even(:,k,1:nhalf))
                sig2o(k) = sig2o(k) + sum(ctf2_odd(:,k,1:nhalf))
            enddo
            !$omp end parallel do
            ! SSNR
            fudge = real(self%p_ptr%tau,dp)
            do k = self%kfromto(1),self%interpklim
                cc = max(0.001d0,min(0.999d0,fsc(k)))
                ssnr(k) = fudge * cc / (1.d0 - cc)
            enddo
            ! Add Tau2 inverse to denominators
            ! because signal assumed infinite at very low resolution there is no addition
            kstart = max(6, calc_fourier_index(self%p_ptr%hp, self%p_ptr%box_crop, self%p_ptr%smpd_crop))
            ! Even
            where( sig2e > DTINY )
                sig2e = real(self%ncls*self%pftsz,dp) / sig2e
            elsewhere
                sig2e = 0.d0
            end where
            tau2 = ssnr * sig2e
            !$omp parallel do default(shared) schedule(static) proc_bind(close) private(k,icls,p,invtau2)
            do k = kstart,self%interpklim
                if( tau2(k) > DTINY )then
                    ! CTF2 <- CTF2 + avgCTF2/(tau*SSNR)
                    invtau2 = 1.d0 / (fudge * tau2(k))
                    ctf2_even(:,k,:) = ctf2_even(:,k,:) + invtau2
                else
                    do icls = 1,self%ncls
                        do p = 1,self%pftsz
                            invtau2 = min(1.d3, 1.d3 * ctf2_even(p,k,icls))
                            ctf2_even(p,k,icls) = ctf2_even(p,k,icls) + invtau2
                        enddo
                    enddo
                endif
            enddo
            !$omp end parallel do
            ! Odd
            where( sig2o > DTINY )
                sig2o = real(self%ncls*self%pftsz,dp) / sig2o
            elsewhere
                sig2o = 0.d0
            end where
            tau2 = ssnr * sig2o
            !$omp parallel do default(shared) schedule(static) proc_bind(close) private(k,icls,p,invtau2)
            do k = kstart,self%interpklim
                if( tau2(k) > DTINY )then
                    invtau2 = 1.d0 / (fudge * tau2(k))
                    ctf2_odd(:,k,:) = ctf2_odd(:,k,:) + invtau2
                else
                    do icls = 1,self%ncls
                        do p = 1,self%pftsz
                            invtau2 = min(1.d3, 1.d3 * ctf2_odd(p,k,icls))
                            ctf2_odd(p,k,icls) = ctf2_odd(p,k,icls) + invtau2
                        enddo
                    enddo
                endif
            enddo
            !$omp end parallel do
        end subroutine add_invtausq2rho

        ! Deals with summing slices and their mirror
        subroutine mirror_slices( ref_space )
            type(oris), intent(in) :: ref_space
            complex(dp) :: pft(self%pftsz,self%kfromto(1):self%interpklim)
            real(dp)    :: ctf2(self%pftsz,self%kfromto(1):self%interpklim)
            real        :: psi
            integer     :: iref, m
            logical     :: l_rotm
            if( .not.ref_space%isthere('mirr') )then
                THROW_HARD('Mirror index missing in reference search space')
            endif
            !$omp parallel do default(shared) proc_bind(close) private(iref,m,pft,ctf2,psi,l_rotm)
            do iref = 1,self%ncls/2
                m      = ref_space%get_int(iref,'mirr')
                psi    = abs(ref_space%get(m, 'psi'))
                l_rotm = (psi > 0.1) .and. (psi < 359.9)
                ! Fourier components
                if( l_rotm )then
                    call self%mirror_pft(self%pfts_even(:,:,m), pft)
                    self%pfts_even(:,:,iref) = self%pfts_even(:,:,iref) + conjg(pft)
                    call self%mirror_pft(self%pfts_odd(:,:,m), pft)
                    self%pfts_odd(:,:,iref)  = self%pfts_odd(:,:,iref)  + conjg(pft)
                    call self%mirror_pft(conjg(self%pfts_even(:,:,iref)), self%pfts_even(:,:,m))
                    call self%mirror_pft(conjg(self%pfts_odd(:,:,iref)),  self%pfts_odd(:,:,m))
                else
                    call self%mirror_pft(self%pfts_even(:,:,m), pft)
                    self%pfts_even(:,:,iref) = self%pfts_even(:,:,iref) + pft
                    call self%mirror_pft(self%pfts_odd(:,:,m), pft)
                    self%pfts_odd(:,:,iref)  = self%pfts_odd(:,:,iref)  + pft
                    call self%mirror_pft(self%pfts_even(:,:,iref), self%pfts_even(:,:,m))
                    call self%mirror_pft(self%pfts_odd(:,:,iref),  self%pfts_odd(:,:,m))
                endif
                ! CTF
                call self%mirror_ctf2(self%ctf2_even(:,:,m), ctf2)
                self%ctf2_even(:,:,iref) = self%ctf2_even(:,:,iref) + ctf2
                call self%mirror_ctf2(self%ctf2_odd(:,:,m),  ctf2)
                self%ctf2_odd(:,:,iref)  = self%ctf2_odd(:,:,iref)  + ctf2
                call self%mirror_ctf2(self%ctf2_even(:,:,iref), self%ctf2_even(:,:,m))
                call self%mirror_ctf2(self%ctf2_odd(:,:,iref),  self%ctf2_odd(:,:,m))
            enddo
            !$omp end parallel do
        end subroutine mirror_slices

        ! Restores slices
        subroutine restore_references
            complex(dp) :: pft(self%pftsz,self%kfromto(1):self%interpklim)
            complex(dp) :: pfte(self%pftsz,self%kfromto(1):self%interpklim)
            complex(dp) :: pfto(self%pftsz,self%kfromto(1):self%interpklim)
            real(dp)    :: ctf2(self%pftsz,self%kfromto(1):self%interpklim)
            real(dp)    :: ctf2e(self%pftsz,self%kfromto(1):self%interpklim)
            real(dp)    :: ctf2o(self%pftsz,self%kfromto(1):self%interpklim)
            real        :: psi
            integer     :: icls, m
            logical     :: l_rotm
            ! Restoration
            !$omp parallel do default(shared) schedule(static) proc_bind(close)&
            !$omp private(icls,m,pft,ctf2,pfte,pfto,ctf2e,ctf2o,psi,l_rotm)
            do icls = 1,self%ncls/2
                ! merged then e/o
                pfte  = pfts_even(:,:,icls)
                pfto  = pfts_odd(:,:,icls)
                ctf2e = ctf2_even(:,:,icls)
                ctf2o = ctf2_odd(:,:,icls)
                pft   = pfte  + pfto
                ctf2  = ctf2e + ctf2o
                call self%safe_norm(pft,  ctf2,  self%pfts_merg(:,:,icls))
                call self%safe_norm(pfte, ctf2e, self%pfts_even(:,:,icls))
                call self%safe_norm(pfto, ctf2o, self%pfts_odd(:,:,icls))
                ! mirror the restored images
                m = reforis%get_int(icls,'mirr')
                call self%mirror_pft(self%pfts_merg(:,:,icls), self%pfts_merg(:,:,m))
                call self%mirror_pft(self%pfts_even(:,:,icls), self%pfts_even(:,:,m))
                call self%mirror_pft(self%pfts_odd(:,:,icls),  self%pfts_odd(:,:,m))
                psi    = abs(reforis%get(m, 'psi'))
                l_rotm = (psi > 0.1) .and. (psi < 359.9)
                if( l_rotm )then
                    self%pfts_merg(:,:,m) = conjg(self%pfts_merg(:,:,m))
                    self%pfts_even(:,:,m) = conjg(self%pfts_even(:,:,m))
                    self%pfts_odd(:,:,m)  = conjg(self%pfts_odd(:,:,m))
                endif
            enddo
            !$omp end parallel do
        end subroutine restore_references

    end subroutine polar_cavger_merge_eos_and_norm

    ! Alignment and filtering routines

    !>  \brief  Filters one polar cluster centre for alignment
    module subroutine polar_prep2Dref( self, clsfrcs, icls, gaufilt )
        class(polarft_calc), intent(inout) :: self
        class(class_frcs),   intent(inout) :: clsfrcs
        integer,             intent(in) :: icls
        logical,             intent(in) :: gaufilt
        real    :: frc(clsfrcs%get_filtsz()), filter(clsfrcs%get_filtsz())
        if( self%p_ptr%l_ml_reg )then
            ! no filtering
        else
            ! FRC-based optimal filter
            call clsfrcs%frc_getter(icls, frc)
            if( any(frc > 0.143) )then
                call fsc2optlp_sub(size(frc), frc, filter, merged=self%p_ptr%l_lpset)
                call self%polar_filterrefs(icls, filter)
            endif
        endif
        ! Optional gaussian filter
        if( gaufilt )then
            call gaussian_filter(self%p_ptr%gaufreq, self%p_ptr%smpd, self%p_ptr%box, filter)
            call self%polar_filterrefs(icls, filter)
        endif
    end subroutine polar_prep2Dref

    !>  \brief prepares a 2D class document with class index, resolution,
    !!         population, average correlation and weight
    module subroutine polar_cavger_gen2Dclassdoc( self, spproj, clsfrcs )
        class(polarft_calc),       intent(in) :: self
        class(sp_project), target, intent(inout) :: spproj
        class(class_frcs),         intent(inout) :: clsfrcs
        class(oris), pointer :: ptcl_field, cls_field
        integer  :: pops(self%ncls)
        real(dp) :: corrs(self%ncls), ws(self%ncls)
        real     :: frc05, frc0143, rstate, w
        integer  :: iptcl, icls, pop, nptcls
        logical  :: l_3D
        l_3D = .false.
        select case(trim(self%p_ptr%oritype))
        case('ptcl2D')
            ptcl_field => spproj%os_ptcl2D
            cls_field  => spproj%os_cls2D
        case('ptcl3D')
            ptcl_field => spproj%os_ptcl3D
            cls_field  => spproj%os_cls3D
            l_3D       = .true.
        case DEFAULT
            THROW_HARD('Unsupported ORITYPE: '//trim(self%p_ptr%oritype))
        end select
        nptcls = ptcl_field%get_noris()
        pops   = 0
        corrs  = 0.d0
        ws     = 0.d0
        !$omp parallel do default(shared) private(iptcl,rstate,icls,w) schedule(static)&
        !$omp proc_bind(close) reduction(+:pops,corrs,ws)
        do iptcl=1,nptcls
            rstate = ptcl_field%get(iptcl,'state')
            if( rstate < 0.5 ) cycle
            w = ptcl_field%get(iptcl,'w')
            if( w < SMALL ) cycle
            if( l_3D )then
                icls = ptcl_field%get_proj(iptcl)
            else
                icls = ptcl_field%get_class(iptcl)
            endif
            if( icls<1 .or. icls>self%ncls )cycle
            pops(icls)  = pops(icls)  + 1
            corrs(icls) = corrs(icls) + real(ptcl_field%get(iptcl,'corr'),dp)
            ws(icls)    = ws(icls)    + real(w,dp)
        enddo
        !$omp end parallel do
        where(pops>1)
            corrs = corrs / real(pops)
            ws    = ws / real(pops)
        elsewhere
            corrs = -1.
            ws    = 0.
        end where
        call cls_field%new(self%ncls, is_ptcl=.false.)
        do icls=1,self%ncls
            pop = pops(icls)
            call clsfrcs%estimate_res(icls, frc05, frc0143)
            if( l_3D )then
                call cls_field%set(icls, 'proj',  icls)
            else
                call cls_field%set(icls, 'class', icls)
            endif
            call cls_field%set(icls, 'pop',       pop)
            call cls_field%set(icls, 'res',       frc0143)
            call cls_field%set(icls, 'corr',      corrs(icls))
            call cls_field%set(icls, 'w',         ws(icls))
            call cls_field%set_state(icls, 1) ! needs to be default val if no selection has been done
            if( pop == 0 )call cls_field%set_state(icls, 0)
        end do
    end subroutine polar_cavger_gen2Dclassdoc

    !>  \brief  Filter references
    !! keep serial
    module subroutine polar_filterrefs( self, icls, filter )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in)    :: icls
        real,                intent(in)    :: filter(:)
        integer :: k
        real    :: fk
        do k = self%kfromto(1),self%interpklim
            fk = filter(k)
            self%pfts_merg(:,k,icls) = fk * self%pfts_merg(:,k,icls)
            self%pfts_even(:,k,icls) = fk * self%pfts_even(:,k,icls)
            self%pfts_odd(:,k,icls)  = fk * self%pfts_odd(:,k,icls)
        enddo
    end subroutine polar_filterrefs

    ! PRIVATE UTILITIES

    ! Calculate common-lines contributions from all the slices
    module subroutine calc_comlin_contrib( self, ref_space, symop, pfts_cl_even, pfts_cl_odd, ctf2_cl_even, ctf2_cl_odd )
        class(polarft_calc), intent(in)    :: self
        type(oris),          intent(in)    :: ref_space
        type(sym),           intent(in)    :: symop
        complex(kind=dp),    intent(inout) :: pfts_cl_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        complex(kind=dp),    intent(inout) :: pfts_cl_odd(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(kind=dp),       intent(inout) :: ctf2_cl_even(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(kind=dp),       intent(inout) :: ctf2_cl_odd(self%pftsz,self%kfromto(1):self%interpklim,self%ncls)
        real(dp), allocatable :: Rsym(:,:,:)
        complex(dp) :: cl_l(self%kfromto(1):self%interpklim), cl_r(self%kfromto(1):self%interpklim)
        complex(dp) :: cl_e(self%kfromto(1):self%interpklim), cl_o(self%kfromto(1):self%interpklim)
        real(dp)    :: rl_l(self%kfromto(1):self%interpklim), rl_r(self%kfromto(1):self%interpklim)
        real(dp)    :: rl_e(self%kfromto(1):self%interpklim), rl_o(self%kfromto(1):self%interpklim)
        real(dp)    :: wl(self%kfromto(1):self%interpklim), wr(self%kfromto(1):self%interpklim)
        real(dp)    :: R(3,3,self%ncls), Rj(3,3), tRi(3,3), eulers(3), psi
        real        :: Rtmp(3,3)
        integer     :: rotl, rotr, iref, jref, m, isym, nsym
        logical     :: l_rotm
        if( .not.ref_space%isthere('mirr') )then
            THROW_HARD('Mirror index missing in reference search space')
        endif
        ! Symmetry rotation matrices
        nsym = symop%get_nsym()
        allocate(Rsym(3,3,nsym),source=0.d0)
        !$omp parallel default(shared) proc_bind(close)&
        !$omp private(iref,jref,m,isym,tRi,Rj,Rtmp,eulers,wl,wr,psi,l_rotm)&
        !$omp& private(cl_l,cl_r,cl_e,cl_o,rl_l,rl_r,rl_e,rl_o,rotl,rotr)
        ! Init
        !$omp workshare
        pfts_cl_even = DCMPLX_ZERO
        pfts_cl_odd  = DCMPLX_ZERO
        ctf2_cl_even = 0.d0
        ctf2_cl_odd  = 0.d0
        !$omp end workshare
        ! Caching rotation matrices
        !$omp do schedule(static)
        do iref = 1,self%ncls
            R(:,:,iref) = real(ref_space%get_mat(iref),dp)
        enddo
        !$omp end do nowait
        !$omp do schedule(static)
        do isym = 1,nsym
            call symop%get_sym_rmat(isym, Rtmp)
            Rsym(:,:,isym) = real(Rtmp,dp)
        end do
        !$omp end do
        ! Common lines contribution
        !$omp do schedule(static)
        do iref = 1,self%ncls/2
            tRi = transpose(R(:,:,iref))
            m   = ref_space%get_int(iref,'mirr')
            do jref = 1,self%ncls/2
                do isym = 1,nsym
                    if( isym == 1 )then
                        if( jref == iref ) cycle    ! self   exclusion
                        if( jref == m )    cycle    ! mirror exclusion
                        ! Rotation of both planes by transpose of Ri (tRixRi=I & Rsym=I)
                        Rj = matmul(R(:,:,jref), tRi)
                    else
                        ! Symmetry operator
                        Rj = matmul(R(:,:,jref), Rsym(:,:,isym))
                        ! Rotation of both planes by transpose of Ri (tRixRi -> I)
                        Rj = matmul(Rj, tRi)
                    endif
                    ! Euler angles identification
                    eulers = dm2euler(Rj)
                    ! in plane rotation angle of jref slice intersecting iref
                    psi = 360.d0 - eulers(3)
                    ! get the weights, rotation indices and compute the interpolated common line
                    call self%gen_clin_weights(psi, rotl, rotr, wl, wr)
                    ! even
                    call self%get_line(jref, rotl, .true., cl_l, rl_l)
                    call self%get_line(jref, rotr, .true., cl_r, rl_r)
                    cl_e = wl*cl_l + wr*cl_r
                    rl_e = wl*rl_l + wr*rl_r
                    ! odd
                    call self%get_line(jref, rotl, .false., cl_l, rl_l)
                    call self%get_line(jref, rotr, .false., cl_r, rl_r)
                    cl_o = wl*cl_l + wr*cl_r
                    rl_o = wl*rl_l + wr*rl_r
                    ! in plane rotation angle of iref slice
                    psi = eulers(1)
                    ! get the weights, rotation indices and extrapolate the common line
                    call self%gen_clin_weights(psi, rotl, rotr, wl, wr)
                    ! leftmost line
                    call extrapolate_line(iref, rotl, wl, cl_e, cl_o, rl_e, rl_o)
                    ! rightmost line
                    call extrapolate_line(iref, rotr, wr, cl_e, cl_o, rl_e, rl_o)
                enddo
            enddo
        enddo
        !$omp end do
        ! Mirroring contributions
        !$omp do schedule(static)
        do iref = 1,self%ncls/2
            m      = ref_space%get_int(iref,'mirr')
            psi    = abs(ref_space%get(m, 'psi'))
            l_rotm = (psi>0.1) .and. (psi<359.9)
            call self%mirror_pft(pfts_cl_even(:,:,iref), pfts_cl_even(:,:,m))
            call self%mirror_pft(pfts_cl_odd(:,:,iref),  pfts_cl_odd(:,:,m))
            if( l_rotm )then
                pfts_cl_even(:,:,m) = conjg(pfts_cl_even(:,:,m))
                pfts_cl_odd(:,:,m) = conjg(pfts_cl_odd(:,:,m))
            endif
            call self%mirror_ctf2(ctf2_cl_even(:,:,iref), ctf2_cl_even(:,:,m))
            call self%mirror_ctf2(ctf2_cl_odd(:,:,iref),  ctf2_cl_odd(:,:,m))
        enddo
        !$omp end do
        !$omp end parallel

    contains

        ! Extrapolate cline, rline to pfts_clin and ctf2_clin
        subroutine extrapolate_line(ref, rot, weight, cle, clo, rle, rlo)
            integer,     intent(in) :: ref, rot
            real(dp),    intent(in) :: weight(self%kfromto(1):self%interpklim)
            complex(dp), intent(in) :: cle(self%kfromto(1):self%interpklim), clo(self%kfromto(1):self%interpklim)
            real(dp),    intent(in) :: rle(self%kfromto(1):self%interpklim), rlo(self%kfromto(1):self%interpklim)
            integer :: irot
            irot = rot
            if( irot < 1 )then
                irot = irot + self%pftsz
                pfts_cl_even(irot,:,ref) = pfts_cl_even(irot,:,ref) + weight * conjg(cle)
                pfts_cl_odd( irot,:,ref) = pfts_cl_odd( irot,:,ref) + weight * conjg(clo)
            elseif( irot > self%pftsz )then
                irot = irot - self%pftsz
                pfts_cl_even(irot,:,ref) = pfts_cl_even(irot,:,ref) + weight * conjg(cle)
                pfts_cl_odd( irot,:,ref) = pfts_cl_odd( irot,:,ref) + weight * conjg(clo)
            else
                pfts_cl_even(irot,:,ref) = pfts_cl_even(irot,:,ref) + weight * cle
                pfts_cl_odd( irot,:,ref) = pfts_cl_odd( irot,:,ref) + weight * clo
            endif
            ctf2_cl_even(irot,:,ref) = ctf2_cl_even(irot,:,ref) + weight * rle
            ctf2_cl_odd( irot,:,ref) = ctf2_cl_odd( irot,:,ref) + weight * rlo
        end subroutine extrapolate_line

    end subroutine calc_comlin_contrib

    ! produces y-mirror of real (reciprocal) matrix 
    module pure subroutine mirror_ctf2( self, ctf2in, ctf2out )
        class(polarft_calc), intent(in)    :: self
        real(dp),            intent(in)    :: ctf2in(self%pftsz,self%kfromto(1):self%interpklim)
        real(dp),            intent(inout) :: ctf2out(self%pftsz,self%kfromto(1):self%interpklim)
        integer :: i,j
        ctf2out(1,:) = ctf2in(1,:)
        do i = 2,self%pftsz/2
            j = self%pftsz-i+2
            ctf2out(i,:) = ctf2in(j,:)
            ctf2out(j,:) = ctf2in(i,:)
        enddo
        i = self%pftsz/2 + 1
        ctf2out(i,:) = ctf2in(i,:)
    end subroutine mirror_ctf2

    ! produces y-mirror of complex (reciprocal) matrix, for particles only
    module pure subroutine mirror_pft( self, pftin, pftout )
        class(polarft_calc), intent(in)    :: self
        complex(dp),         intent(in)    :: pftin(self%pftsz,self%kfromto(1):self%interpklim)
        complex(dp),         intent(inout) :: pftout(self%pftsz,self%kfromto(1):self%interpklim)
        integer :: i,j
        pftout(1,:) = conjg(pftin(1,:))
        do i = 2,self%pftsz/2
            j = self%pftsz-i+2
            pftout(i,:) = pftin(j,:)
            pftout(j,:) = pftin(i,:)
        enddo
        i = self%pftsz/2 + 1
        pftout(i,:) = pftin(i,:)
    end subroutine mirror_pft

    ! Private utility
    module pure subroutine safe_norm( self, Mnum, Mdenom, Mout )
        class(polarft_calc), intent(in)    :: self
        complex(dp),         intent(in)    :: Mnum(self%pftsz,self%kfromto(1):self%interpklim)
        real(dp),            intent(inout) :: Mdenom(self%pftsz,self%kfromto(1):self%interpklim)
        complex(dp),         intent(inout) :: Mout(self%pftsz,self%kfromto(1):self%interpklim)
        logical  :: msk(self%pftsz)
        real(dp) :: avg, t
        integer  :: k, n
        do k = self%kfromto(1),self%interpklim
            msk = Mdenom(:,k) > DSMALL
            n   = count(msk)
            if( n == 0 ) cycle
            avg = sum(Mdenom(:,k),mask=msk) / real(n,dp)
            t   = avg/50.d0
            where((Mdenom(:,k) < t).and.msk) Mdenom(:,k) = Mdenom(:,k) + t
        enddo
        where( Mdenom > DSMALL)
            Mout = Mnum / Mdenom
        elsewhere
            Mout = DCMPLX_ZERO
        end where
    end subroutine safe_norm

    ! Returns complex and ctf2 polar lines given ref and rotational indices
    module pure subroutine get_line( self, ref, rot, even, pftline, ctf2line )
        class(polarft_calc), intent(in)  :: self
        integer,             intent(in)  :: ref, rot
        logical,             intent(in)  :: even
        complex(dp),         intent(out) :: pftline(self%kfromto(1):self%interpklim)
        real(dp),            intent(out) :: ctf2line(self%kfromto(1):self%interpklim)
        integer :: irot
        if( rot >  self%nrots )then
            irot = rot - self%nrots
        else
            irot = rot
        endif
        if( even )then
            if( irot < 1 )then
                irot    = irot + self%pftsz
                pftline = conjg(self%pfts_even(irot,:,ref))
            elseif( irot > self%pftsz )then
                irot    = irot - self%pftsz
                pftline = conjg(self%pfts_even(irot,:,ref))
            else
                pftline = self%pfts_even(irot,:,ref)
            endif
            ctf2line = self%ctf2_even(irot,:,ref)
        else
            if( irot < 1 )then
                irot    = irot + self%pftsz
                pftline = conjg(self%pfts_odd(irot,:,ref))
            elseif( irot > self%pftsz )then
                irot    = irot - self%pftsz
                pftline = conjg(self%pfts_odd(irot,:,ref))
            else
                pftline = self%pfts_odd(irot,:,ref)
            endif
            ctf2line = self%ctf2_odd(irot,:,ref)
        endif
    end subroutine get_line

    module subroutine polar_cavger_calc_frc( self, pft1, pft2, n, frc )
        class(polarft_calc), intent(in)    :: self
        complex(dp),         intent(in)    :: pft1(self%pftsz,self%kfromto(1):self%interpklim), pft2(self%pftsz,self%kfromto(1):self%interpklim)
        integer,             intent(in)    :: n
        real(sp),            intent(inout) :: frc(1:n)
        real(dp) :: var1, var2, denom
        integer  :: k
        frc(1:self%kfromto(1)-1) = 0.999
        do k = self%kfromto(1), self%interpklim
            var1  = sum(csq_fast(pft1(:,k)))
            var2  = sum(csq_fast(pft2(:,k)))
            if( (var1>DTINY) .and. (var2>DTINY) )then
                denom  = sqrt(var1*var2)
                frc(k) = real(sum(pft1(:,k)*conjg(pft2(:,k))) / denom, sp)
            else
                frc(k) = 0.0
            endif
        enddo
        if( self%interpklim < n ) frc(self%interpklim+1:) = 0.0
    end subroutine polar_cavger_calc_frc

end submodule simple_polarft_ops_restore
