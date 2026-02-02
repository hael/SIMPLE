!@descr: submodule for class average restoration in the polar Fourier domain
submodule (simple_polarft_calc) simple_polarft_ops_restore
use simple_class_frcs, only: class_frcs
use simple_builder,    only: build_glob
use simple_cmdline,    only: cmdline
implicit none
#include "simple_local_flags.inc"
contains               

    !>  \brief  Restores 2D class-averages only
    module subroutine polar_cavger_merge_eos_and_norm2D( self )
        class(polarft_calc), intent(inout) :: self
        complex(dp) :: pft(self%pftsz,self%kfromto(1):self%kfromto(2))
        real(dp)    :: ctf2(self%pftsz,self%kfromto(1):self%kfromto(2))
        integer     :: icls, eo_pop(2), pop
        select case(trim(params_glob%ref_type))
        case('polar_cavg')
            ! all good
        case DEFAULT
            THROW_HARD('polar_cavger_merge_eos_and_norm2D only for 2D cavgs restoration')
        end select
        self%pfts_merg = DCMPLX_ZERO
        !$omp parallel do default(shared) schedule(static) proc_bind(close)&
        !$omp private(icls,eo_pop,pop,pft,ctf2)
        do icls = 1,self%ncls
            eo_pop = self%prev_eo_pops(:,icls) + self%eo_pops(:,icls)
            pop    = sum(eo_pop)
            if(pop == 0)then
                self%pfts_even(:,:,icls) = DCMPLX_ZERO
                self%pfts_odd(:,:,icls)  = DCMPLX_ZERO
                self%ctf2_even(:,:,icls) = 0.d0
                self%ctf2_odd(:,:,icls)  = 0.d0
            else
                if(pop > 1)then
                    pft  = self%pfts_even(:,:,icls) + self%pfts_odd(:,:,icls)
                    ctf2 = self%ctf2_even(:,:,icls) + self%ctf2_odd(:,:,icls)
                    call self%safe_norm(pft, ctf2, self%pfts_merg(:,:,icls))
                endif
                if(eo_pop(1) > 1)then
                    pft = self%pfts_even(:,:,icls)
                    call self%safe_norm(pft, self%ctf2_even(:,:,icls), self%pfts_even(:,:,icls))
                endif
                if(eo_pop(2) > 1)then
                    pft = self%pfts_odd(:,:,icls)
                    call self%safe_norm(pft, self%ctf2_odd(:,:,icls), self%pfts_odd(:,:,icls))
                endif
            endif
        end do
        !$omp end parallel do
    end subroutine polar_cavger_merge_eos_and_norm2D

    !>  \brief  Restores 3D slices
    module subroutine polar_cavger_merge_eos_and_norm( self, reforis, cl_weight )
        class(polarft_calc),  intent(inout) :: self
        type(oris),           intent(in)    :: reforis
        real,       optional, intent(in)    :: cl_weight
        type(class_frcs)   :: cavg2clin_frcs
        real,  allocatable :: cavg_clin_frcs(:,:,:)
        complex(dp) :: pfts_clin_even(self%pftsz,self%kfromto(1):self%kfromto(2),self%ncls)
        complex(dp) :: pfts_clin_odd(self%pftsz,self%kfromto(1):self%kfromto(2),self%ncls)
        real(dp)    :: ctf2_clin_even(self%pftsz,self%kfromto(1):self%kfromto(2),self%ncls)
        real(dp)    :: ctf2_clin_odd(self%pftsz,self%kfromto(1):self%kfromto(2),self%ncls)
        real(dp)    :: fsc(self%kfromto(1):self%kfromto(2)), clw
        clw = 1.d0
        select case(trim(params_glob%ref_type))
            case('comlin_hybrid')
                ! 2.5D: cavgs + variable CLs contribution
                ! Common-line contribution weight
                if( present(cl_weight) ) clw = min(max(0.d0,real(cl_weight,dp)),1.d0)
                if( clw > 1.d-6 )then
                    ! Mirroring slices
                    call mirror_slices(reforis, build_glob%pgrpsyms)
                    ! Common-lines conribution
                    call self%calc_comlin_contrib(reforis, build_glob%pgrpsyms,&
                    &pfts_clin_even, pfts_clin_odd, ctf2_clin_even, ctf2_clin_odd)
                endif
            case('comlin_noself', 'comlin')
                ! Mirroring slices
                call mirror_slices(reforis, build_glob%pgrpsyms)
                ! Common-lines conribution
                call self%calc_comlin_contrib(reforis, build_glob%pgrpsyms,&
                &pfts_clin_even, pfts_clin_odd, ctf2_clin_even, ctf2_clin_odd)
            case DEFAULT
                THROW_HARD('Invalid REF_TYPE='//trim(params_glob%ref_type)//' in polar_cavger_merge_eos_and_norm')
        end select
        ! ML-regularization
        if( params_glob%l_ml_reg ) call add_invtausq2rho
        ! Restoration of references
        self%pfts_merg = DCMPLX_ZERO
        select case(trim(params_glob%ref_type))
            case('comlin_hybrid')
                if( clw > 1.d-6 )then
                    call calc_cavg_comlin_frcs(cavg2clin_frcs)
                    call restore_cavgs_comlins(clw)
                else
                    call mirr_and_calc_cavg_comlin_frcs(cavg2clin_frcs)
                    call self%polar_cavger_merge_eos_and_norm2D
                endif
                call cavg2clin_frcs%write(string('frcs_cavg2clin'//BIN_EXT))
            case('comlin_noself')
                call restore_comlins
            case('comlin')
                call calc_cavg_comlin_frcs( cavg2clin_frcs )
                call cavg2clin_frcs%write(string('frcs_cavg2clin'//BIN_EXT))
                call restore_cavgs_comlins( 1.d0 )
        end select
        ! cleanup
        call cavg2clin_frcs%kill
        if( allocated(cavg_clin_frcs) ) deallocate(cavg_clin_frcs)

    contains

        subroutine add_invtausq2rho
            logical, parameter :: DEBUG = .false.
            complex(dp) :: even(self%pftsz,self%kfromto(1):self%kfromto(2)), odd(self%pftsz,self%kfromto(1):self%kfromto(2))
            complex(dp) :: pft(self%pftsz,self%kfromto(1):self%kfromto(2))
            real(dp)    :: vare(self%kfromto(1):self%kfromto(2)), varo(self%kfromto(1):self%kfromto(2))
            real(dp)    :: sig2e(self%kfromto(1):self%kfromto(2)), sig2o(self%kfromto(1):self%kfromto(2))
            real(dp)    :: ssnr(self%kfromto(1):self%kfromto(2)), tau2(self%kfromto(1):self%kfromto(2))
            real(dp)    :: ctf2(self%pftsz,self%kfromto(1):self%kfromto(2)), cc, fudge, invtau2, s, a,b
            integer     :: icls, k, kstart, p
            ! FSC per slice contribution
            fsc  = 0.d0
            vare = 0.d0; varo = 0.d0
            sig2e = 0.d0; sig2o = 0.d0
            select case(trim(params_glob%ref_type))
            case('comlin_hybrid')
                THROW_HARD('Not supported yet')
            case('comlin_noself')
                !$omp parallel do default(shared) schedule(static) proc_bind(close)&
                !$omp private(icls,even,odd,k,pft,ctf2) reduction(+:fsc,vare,varo,sig2e,sig2o)
                do icls = 1,self%ncls/2
                    ! e/o restoration
                    pft  = pfts_clin_even(:,:,icls)
                    ctf2 = ctf2_clin_even(:,:,icls)
                    call self%safe_norm(pft, ctf2, even)
                    pft  = pfts_clin_odd(:,:,icls)
                    ctf2 = ctf2_clin_odd(:,:,icls)
                    call self%safe_norm(pft, ctf2, odd)
                    ! FSC contribution
                    do k = self%kfromto(1),self%kfromto(2)
                        fsc(k)   = fsc(k)   + sum(real(even(:,k) * conjg(odd(:,k)), dp))
                        vare(k)  = vare(k)  + sum(real(even(:,k) * conjg(even(:,k)),dp))
                        varo(k)  = varo(k)  + sum(real(odd(:,k)  * conjg(odd(:,k)), dp))
                        sig2e(k) = sig2e(k) + sum(ctf2_clin_even(:,k,icls))
                        sig2o(k) = sig2o(k) + sum(ctf2_clin_odd(:,k,icls))
                    enddo
                enddo
                !$omp end parallel do
            case('comlin')
                !$omp parallel do default(shared) schedule(static) proc_bind(close)&
                !$omp private(icls,even,odd,k,pft,ctf2) reduction(+:fsc,vare,varo,sig2e,sig2o)
                do icls = 1,self%ncls/2
                    ! e/o restoration
                    pft  = self%pfts_even(:,:,icls) + pfts_clin_even(:,:,icls)
                    ctf2 = self%ctf2_even(:,:,icls) + ctf2_clin_even(:,:,icls)
                    call self%safe_norm(pft, ctf2, even)
                    pft  = self%pfts_odd(:,:,icls) + pfts_clin_odd(:,:,icls)
                    ctf2 = self%ctf2_odd(:,:,icls) + ctf2_clin_odd(:,:,icls)
                    call self%safe_norm(pft, ctf2, odd)
                    ! FSC contribution
                    do k = self%kfromto(1),self%kfromto(2)
                        fsc(k)   = fsc(k)   + real(sum(even(:,k) * conjg(odd(:,k))) ,dp)
                        vare(k)  = vare(k)  + real(sum(even(:,k) * conjg(even(:,k))),dp)
                        varo(k)  = varo(k)  + real(sum(odd(:,k)  * conjg(odd(:,k))) ,dp)
                        sig2e(k) = sig2e(k) + sum(self%ctf2_even(:,k,icls) + ctf2_clin_even(:,k,icls))
                        sig2o(k) = sig2o(k) + sum(self%ctf2_odd(:,k,icls)  + ctf2_clin_odd(:,k,icls))
                    enddo
                enddo
                !$omp end parallel do
            end select
            ! Variances
            vare = vare * varo
            ! FSC
            where( vare > DTINY )
                fsc = fsc / sqrt(vare)
            elsewhere
                fsc = 0.d0
            end where
            ! SSNR
            fudge = real(params_glob%tau,dp)
            do k = self%kfromto(1),self%kfromto(2)
                cc = max(0.001d0,min(0.999d0,fsc(k)))
                ssnr(k) = fudge * cc / (1.d0 - cc)
            enddo
            ! Add Tau2 inverse to denominators
            ! because signal assumed infinite at very low resolution there is no addition
            kstart = max(6, calc_fourier_index(params_glob%hp, params_glob%box_crop, params_glob%smpd_crop))
            ! Even
            where( sig2e > DTINY )
                sig2e = real(self%ncls*self%pftsz,dp) / sig2e
            elsewhere
                sig2e = 0.d0
            end where
            tau2 = ssnr * sig2e
            if( DEBUG )then
                do k = self%kfromto(1),self%kfromto(2)
                    cc = max(0.001d0,min(0.999d0,fsc(k)))
                    s = cc / (1.d0 - cc)
                    print *, k, real(s), real(fsc(k)),&
                    &real(fudge**2*s/(fudge**2*s+1.d0)), tau2(k), 1.d0 / (fudge * tau2(k))
                enddo
                do k = self%kfromto(1),self%kfromto(2)
                    if( tau2(k) > DTINY )then
                        invtau2 = 1.d0 / (fudge * tau2(k))
                        a = sum(sqrt(real(pfts_clin_even(:,k,:)*conjg(pfts_clin_even(:,k,:))))) / real(self%pftsz*self%ncls,dp)
                        b =      sum(ctf2_clin_even(:,k,:))                                     / real(self%pftsz*self%ncls,dp)
                        print *,k,(a/(b+invtau2)) / (a/b)
                    endif
                enddo
            endif
            !$omp parallel do default(shared) schedule(static) proc_bind(close) private(k,icls,p,invtau2)
            do k = kstart,self%kfromto(2)
                if( tau2(k) > DTINY )then
                    ! CTF2 <- CTF2 + avgCTF2/(tau*SSNR)
                    invtau2 = 1.d0 / (fudge * tau2(k))
                    ctf2_clin_even(:,k,:) = ctf2_clin_even(:,k,:) + invtau2
                else
                    do icls = 1,self%ncls
                        do p = 1,self%pftsz
                            invtau2 = min(1.d3, 1.d3 * ctf2_clin_even(p,k,icls))
                            ctf2_clin_even(p,k,icls) = ctf2_clin_even(p,k,icls) + invtau2
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
            do k = kstart,self%kfromto(2)
                if( tau2(k) > DTINY )then
                    invtau2 = 1.d0 / (fudge * tau2(k))
                    ctf2_clin_odd(:,k,:) = ctf2_clin_odd(:,k,:) + invtau2
                else
                    do icls = 1,self%ncls
                        do p = 1,self%pftsz
                            invtau2 = min(1.d3, 1.d3 * ctf2_clin_odd(p,k,icls))
                            ctf2_clin_odd(p,k,icls) = ctf2_clin_odd(p,k,icls) + invtau2
                        enddo
                    enddo
                endif
            enddo
            !$omp end parallel do
        end subroutine add_invtausq2rho

        ! Mirror 2D classes, calculate CL contribution & FRCS
        ! Module arrays are untouched on exit
        subroutine mirr_and_calc_cavg_comlin_frcs( frcs )
            class(class_frcs), intent(inout) :: frcs
            complex(dp) :: pfte_backup(self%pftsz,self%kfromto(1):self%kfromto(2),self%ncls)
            complex(dp) :: pfto_backup(self%pftsz,self%kfromto(1):self%kfromto(2),self%ncls)
            real(dp)    :: ctf2e_backup(self%pftsz,self%kfromto(1):self%kfromto(2),self%ncls)
            real(dp)    :: ctf2o_backup(self%pftsz,self%kfromto(1):self%kfromto(2),self%ncls)
            ! backup classes
            pfte_backup(:,:,:)  = self%pfts_even(:,:,:); pfto_backup(:,:,:)  = self%pfts_odd(:,:,:)
            ctf2e_backup(:,:,:) = self%ctf2_even(:,:,:); ctf2o_backup(:,:,:) = self%ctf2_odd(:,:,:)
            ! mirror classes belonging to the same slice
            call mirror_slices(reforis, build_glob%pgrpsyms)
            ! calculate the per slice CLs
            call self%calc_comlin_contrib(reforis, build_glob%pgrpsyms,&
            &pfts_clin_even, pfts_clin_odd, ctf2_clin_even, ctf2_clin_odd)
            ! CLs vs. CLS FRCs
            call calc_cavg_comlin_frcs(frcs)
            ! restores classes
            self%pfts_even(:,:,:) = pfte_backup(:,:,:);  self%pfts_odd(:,:,:) = pfto_backup(:,:,:)
            self%ctf2_even(:,:,:) = ctf2e_backup(:,:,:); self%ctf2_odd(:,:,:) = ctf2o_backup(:,:,:)
        end subroutine mirr_and_calc_cavg_comlin_frcs

        ! Calculate CL contribution & FRCS from mirrored 2D classes
        ! Module arrays are untouched on exit
        subroutine calc_cavg_comlin_frcs( frcs )
            class(class_frcs), intent(inout) :: frcs
            real, allocatable :: frc(:)
            complex(dp)       :: cavg(self%pftsz,self%kfromto(1):self%kfromto(2)), clin(self%pftsz,self%kfromto(1):self%kfromto(2))
            complex(dp)       :: pft(self%pftsz,self%kfromto(1):self%kfromto(2))
            real(dp)          :: ctf2(self%pftsz,self%kfromto(1):self%kfromto(2))
            integer           :: icls, pop, nk
            call frcs%new(self%ncls, params_glob%box, params_glob%smpd, 1)
            nk = frcs%get_filtsz()
            allocate(frc(1:nk),source=0.)
            !$omp parallel do default(shared) schedule(static) proc_bind(close)&
            !$omp private(icls,pop,pft,ctf2,cavg,clin,frc)
            do icls = 1,self%ncls
                pop = sum(self%prev_eo_pops(:,icls) + self%eo_pops(:,icls))
                if( pop > 1 )then
                    ! cavg
                    pft  = self%pfts_even(:,:,icls) + self%pfts_odd(:,:,icls)
                    ctf2 = self%ctf2_even(:,:,icls) + self%ctf2_odd(:,:,icls)
                    call self%safe_norm(pft, ctf2, cavg)
                    ! comlin
                    pft  = pfts_clin_even(:,:,icls) + pfts_clin_odd(:,:,icls)
                    ctf2 = ctf2_clin_even(:,:,icls) + ctf2_clin_odd(:,:,icls)
                    call self%safe_norm(pft, ctf2, clin)
                    ! FRC
                    call self%polar_cavger_calc_frc(cavg, clin, nk, frc)
                else
                    frc = 0.0
                endif
                call frcs%set_frc(icls, frc)
            end do
            !$omp end parallel do
            deallocate(frc)
        end subroutine calc_cavg_comlin_frcs

        ! Deals with summing slices and their mirror
        subroutine mirror_slices( ref_space, symop )
            type(oris), intent(in) :: ref_space
            type(sym),  intent(in) :: symop
            complex(dp) :: pft(self%pftsz,self%kfromto(1):self%kfromto(2))
            real(dp)    :: ctf2(self%pftsz,self%kfromto(1):self%kfromto(2))
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

        ! Restore common lines contributions only
        subroutine restore_comlins
            complex(dp) :: pft(self%pftsz,self%kfromto(1):self%kfromto(2)),pfte(self%pftsz,self%kfromto(1):self%kfromto(2)),pfto(self%pftsz,self%kfromto(1):self%kfromto(2))
            real(dp)    :: ctf2(self%pftsz,self%kfromto(1):self%kfromto(2)),ctf2e(self%pftsz,self%kfromto(1):self%kfromto(2)),ctf2o(self%pftsz,self%kfromto(1):self%kfromto(2))
            real        :: psi
            integer     :: icls, m
            logical     :: l_rotm
            !$omp parallel do default(shared) schedule(static) proc_bind(close)&
            !$omp private(icls,m,pft,ctf2,pfte,pfto,ctf2e,ctf2o,psi,l_rotm)
            do icls = 1,self%ncls/2
                ! already mirrored common-line contribution
                pfte  = pfts_clin_even(:,:,icls)
                pfto  = pfts_clin_odd(:,:,icls)
                ctf2e = ctf2_clin_even(:,:,icls)
                ctf2o = ctf2_clin_odd(:,:,icls)
                ! merged then e/o
                pft   = pfte  + pfto
                ctf2  = ctf2e + ctf2o
                call self%safe_norm(pft,  ctf2,  self%pfts_merg(:,:,icls))
                call self%safe_norm(pfte, ctf2e, self%pfts_even(:,:,icls))
                call self%safe_norm(pfto, ctf2o, self%pfts_odd(:,:,icls))
                ! mirroring the restored images
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
        end subroutine restore_comlins

        ! Restores slices (=cavgs + comlin)
        subroutine restore_cavgs_comlins( clw )
            real(dp), intent(in) :: clw
            complex(dp) :: pft(self%pftsz,self%kfromto(1):self%kfromto(2)),pfte(self%pftsz,self%kfromto(1):self%kfromto(2))
            complex(dp) :: pfto(self%pftsz,self%kfromto(1):self%kfromto(2))
            real(dp)    :: ctf2(self%pftsz,self%kfromto(1):self%kfromto(2)),ctf2e(self%pftsz,self%kfromto(1):self%kfromto(2))
            real(dp)    :: ctf2o(self%pftsz,self%kfromto(1):self%kfromto(2))
            real        :: psi
            integer     :: icls, m
            logical     :: l_rotm
            ! Restoration
            !$omp parallel do default(shared) schedule(static) proc_bind(close)&
            !$omp private(icls,m,pft,ctf2,pfte,pfto,ctf2e,ctf2o,psi,l_rotm)
            do icls = 1,self%ncls/2
                ! calculate sum of already mirrored class + common-line contribution
                pfte  = self%pfts_even(:,:,icls) + clw * pfts_clin_even(:,:,icls)
                pfto  = self%pfts_odd(:,:,icls)  + clw * pfts_clin_odd(:,:,icls)
                ctf2e = self%ctf2_even(:,:,icls) + clw * ctf2_clin_even(:,:,icls)
                ctf2o = self%ctf2_odd(:,:,icls)  + clw * ctf2_clin_odd(:,:,icls)
                ! merged then e/o
                pft   = pfte  + pfto
                ctf2  = ctf2e + ctf2o
                call self%safe_norm(pft,  ctf2,  self%pfts_merg(:,:,icls))
                call self%safe_norm(pfte, ctf2e, self%pfts_even(:,:,icls))
                call self%safe_norm(pfto, ctf2o, self%pfts_odd(:,:,icls))
                ! mirroring the restored images
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
        end subroutine restore_cavgs_comlins

    end subroutine polar_cavger_merge_eos_and_norm

    !>  \brief  calculates Fourier ring correlations
    module subroutine polar_cavger_calc_and_write_frcs_and_eoavg( self, fname, cline )
        class(polarft_calc), intent(inout) :: self
        class(string),       intent(in)    :: fname
        type(cmdline),       intent(in)    :: cline
        complex(dp), allocatable :: prev_prefs(:,:,:)
        real,        allocatable :: frc(:)
        real(dp) :: ufrac_trec
        integer  :: icls, find, pop, filtsz
        filtsz = fdim(params_glob%box_crop) - 1
        allocate(frc(filtsz),source=0.)
        ! In case nspace/self%ncls has changed OR volume/frcs were downsampled
        if( (build_glob%clsfrcs%get_ncls() /= self%ncls) .or. (build_glob%clsfrcs%get_filtsz() /= filtsz) )then
            call build_glob%clsfrcs%new(self%ncls, params_glob%box_crop, params_glob%smpd_crop, params_glob%nstates)
        endif
        !$omp parallel do default(shared) private(icls,frc,find,pop) schedule(static) proc_bind(close)
        do icls = 1,self%ncls
            if( self%l_comlin )then
                ! calculate FRC (pseudo-cavgs are never empty)
                call self%polar_cavger_calc_frc(self%pfts_even(:,:,icls), self%pfts_odd(:,:,icls), filtsz, frc)
                call build_glob%clsfrcs%set_frc(icls, frc, 1)
                ! average low-resolution info between eo pairs to keep things in register
                find = min(self%kfromto(2), build_glob%clsfrcs%estimate_find_for_eoavg(icls, 1))
                if( find >= self%kfromto(1) )then
                    self%pfts_even(:,self%kfromto(1):find,icls) = self%pfts_merg(:,self%kfromto(1):find,icls)
                    self%pfts_odd(:,self%kfromto(1):find,icls)  = self%pfts_merg(:,self%kfromto(1):find,icls)
                endif
            else
                pop = sum(self%prev_eo_pops(:,icls) + self%eo_pops(:,icls))
                if( pop == 0 )then
                    frc = 0.
                    call build_glob%clsfrcs%set_frc(icls, frc, 1)
                else
                    ! calculate FRC
                    call self%polar_cavger_calc_frc(self%pfts_even(:,:,icls), self%pfts_odd(:,:,icls), filtsz, frc)
                    call build_glob%clsfrcs%set_frc(icls, frc, 1)
                    ! average low-resolution info between eo pairs to keep things in register
                    find = min(self%kfromto(2), build_glob%clsfrcs%estimate_find_for_eoavg(icls, 1))
                    if( find >= self%kfromto(1) )then
                        self%pfts_even(:,self%kfromto(1):find,icls) = self%pfts_merg(:,self%kfromto(1):find,icls)
                        self%pfts_odd(:,self%kfromto(1):find,icls)  = self%pfts_merg(:,self%kfromto(1):find,icls)
                    endif
                endif
            endif
        end do
        !$omp end parallel do
        ! write FRCs
        call build_glob%clsfrcs%write(fname)
        ! e/o Trailing reconstruction
        if( params_glob%l_trail_rec )then
            if( cline%defined('ufrac_trec') )then
                ufrac_trec = real(params_glob%ufrac_trec,dp)
            else
                ufrac_trec = real(build_glob%spproj_field%get_update_frac(),dp)
            endif
            call self%read_pft_array(string(POLAR_REFS_FBODY)//'_even'//BIN_EXT, prev_prefs)
            !$omp parallel workshare proc_bind(close)
            self%pfts_even = ufrac_trec * self%pfts_even + (1.d0-ufrac_trec) * prev_prefs
            !$omp end parallel workshare
            call self%read_pft_array(string(POLAR_REFS_FBODY)//'_odd'//BIN_EXT, prev_prefs)
            !$omp parallel workshare proc_bind(close)
            self%pfts_odd  = ufrac_trec * self%pfts_odd   + (1.d0-ufrac_trec) * prev_prefs
            self%pfts_merg = 0.5d0 * (self%pfts_even + self%pfts_odd)
            !$omp end parallel workshare
            deallocate(prev_prefs)
        endif
    end subroutine polar_cavger_calc_and_write_frcs_and_eoavg

    !>  \brief  Filters one polar cluster centre for alignment
    module subroutine polar_prep2Dref( self, icls, gaufilt )
        class(polarft_calc), intent(inout) :: self
        integer,             intent(in) :: icls
        logical,             intent(in) :: gaufilt
        real    :: frc(build_glob%clsfrcs%get_filtsz())
        real    :: filter(build_glob%clsfrcs%get_filtsz())
        real    :: gaufilter(build_glob%clsfrcs%get_filtsz())
        integer :: k
        if( params_glob%l_ml_reg )then
            ! no filtering, not supported yet in 2D
        else
            ! FRC-based optimal filter
            call build_glob%clsfrcs%frc_getter(icls, frc)
            if( any(frc > 0.143) )then
                call fsc2optlp_sub(size(frc), frc, filter, merged=params_glob%l_lpset)
            else
                filter = 1.0
            endif
            ! gaussian filter
            if( gaufilt )then
                call gaussian_filter(params_glob%gaufreq, params_glob%smpd, params_glob%box, gaufilter)
                ! minimum of FRC-based & gaussian filters
                forall(k = 1:size(filter)) filter(k) = min(filter(k), gaufilter(k))
            endif
            call self%polar_filterrefs(icls, filter)
        endif
    end subroutine polar_prep2Dref

    !>  \brief prepares a 2D class document with class index, resolution,
    !!         population, average correlation and weight
    module subroutine polar_cavger_gen2Dclassdoc( self, spproj )
        class(polarft_calc),       intent(in) :: self
        class(sp_project), target, intent(inout) :: spproj
        class(oris), pointer :: ptcl_field, cls_field
        integer  :: pops(self%ncls)
        real(dp) :: corrs(self%ncls), ws(self%ncls)
        real     :: frc05, frc0143, rstate, w
        integer  :: iptcl, icls, pop, nptcls
        logical  :: l_3D
        l_3D = .false.
        select case(trim(params_glob%oritype))
        case('ptcl2D')
            ptcl_field => spproj%os_ptcl2D
            cls_field  => spproj%os_cls2D
        case('ptcl3D')
            ptcl_field => spproj%os_ptcl3D
            cls_field  => spproj%os_cls3D
            l_3D       = .true.
        case DEFAULT
            THROW_HARD('Unsupported ORITYPE: '//trim(params_glob%oritype))
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
            call build_glob%clsfrcs%estimate_res(icls, frc05, frc0143)
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
        do k = self%kfromto(1),self%kfromto(2)
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
        complex(kind=dp),    intent(inout) :: pfts_cl_even(self%pftsz,self%kfromto(1):self%kfromto(2),self%ncls)
        complex(kind=dp),    intent(inout) :: pfts_cl_odd(self%pftsz,self%kfromto(1):self%kfromto(2),self%ncls)
        real(kind=dp),       intent(inout) :: ctf2_cl_even(self%pftsz,self%kfromto(1):self%kfromto(2),self%ncls)
        real(kind=dp),       intent(inout) :: ctf2_cl_odd(self%pftsz,self%kfromto(1):self%kfromto(2),self%ncls)
        real, allocatable :: Rsym(:,:,:)
        complex(dp) :: cl_l(self%kfromto(1):self%kfromto(2)), cl_r(self%kfromto(1):self%kfromto(2)), cl_e(self%kfromto(1):self%kfromto(2))
        complex(dp) :: cl_o(self%kfromto(1):self%kfromto(2)), pft(self%pftsz,self%kfromto(1):self%kfromto(2))
        real(dp)    :: rl_l(self%kfromto(1):self%kfromto(2)), rl_r(self%kfromto(1):self%kfromto(2)), rl_e(self%kfromto(1):self%kfromto(2))
        real(dp)    :: rl_o(self%kfromto(1):self%kfromto(2)), ctf2(self%pftsz,self%kfromto(1):self%kfromto(2))
        real(dp)    :: w(self%kfromto(1):self%kfromto(2)), wl(self%kfromto(1):self%kfromto(2)), wr(self%kfromto(1):self%kfromto(2)), sumw(self%kfromto(1):self%kfromto(2))
        real        :: eulers(3),R(3,3,self%ncls),Rj(3,3),tRi(3,3),psi,drot,d
        integer     :: rotl, rotr, iref, jref, m, self_irot, targ_irot, isym, nsym
        logical     :: l_rotm
        if( .not.ref_space%isthere('mirr') )then
            THROW_HARD('Mirror index missing in reference search space')
        endif
        drot         = self%get_dang()
        pfts_cl_even = DCMPLX_ZERO; pfts_cl_odd = DCMPLX_ZERO
        ctf2_cl_even = 0.d0; ctf2_cl_odd  = 0.d0
        ! Symmetry rotation matrices
        nsym = symop%get_nsym()
        allocate(Rsym(3,3,nsym),source=0.)
        !$omp parallel default(shared) proc_bind(close)&
        !$omp private(iref,jref,m,tRi,Rj,eulers,targ_irot,self_irot)&
        !$omp& private(d,w,wl,wr,sumw,psi,l_rotm,cl_l,cl_r,cl_e,cl_o)&
        !$omp& private(rl_l,rl_r,rl_e,rl_o,pft,ctf2,rotl,rotr,isym)
        ! Caching rotation matrices
        !$omp do schedule(static)
        do iref = 1,self%ncls
            R(:,:,iref) = ref_space%get_mat(iref)
        enddo
        !$omp end do nowait
        !$omp do schedule(static)
        do isym = 1,nsym
            call symop%get_sym_rmat(isym, Rsym(:,:,isym))
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
                    eulers = m2euler_fast(Rj)
                    ! in plane rotation angle of jref slice intersecting iref
                    psi = 360.0 - eulers(3)
                    ! get the weights, rotation indeces and compute the interpolated common line
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
                    ! get the weights, rotation indeces and extrapolate the common line
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
            if( l_rotm )then
                call self%mirror_pft(conjg(pfts_cl_even(:,:,iref)), pfts_cl_even(:,:,m))
                call self%mirror_pft(conjg(pfts_cl_odd(:,:,iref)),  pfts_cl_odd(:,:,m))
            else
                call self%mirror_pft(pfts_cl_even(:,:,iref), pfts_cl_even(:,:,m))
                call self%mirror_pft(pfts_cl_odd(:,:,iref),  pfts_cl_odd(:,:,m))
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
            real(dp),    intent(in) :: weight(self%kfromto(1):self%kfromto(2))
            complex(dp), intent(in) :: cle(self%kfromto(1):self%kfromto(2)), clo(self%kfromto(1):self%kfromto(2))
            real(dp),    intent(in) :: rle(self%kfromto(1):self%kfromto(2)), rlo(self%kfromto(1):self%kfromto(2))
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
        real(dp),            intent(in)    :: ctf2in(self%pftsz,self%kfromto(1):self%kfromto(2))
        real(dp),            intent(inout) :: ctf2out(self%pftsz,self%kfromto(1):self%kfromto(2))
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

    ! produces y-mirror of complex (reciprocal) matrix 
    module pure subroutine mirror_pft( self, pftin, pftout )
        class(polarft_calc), intent(in)    :: self
        complex(dp),         intent(in)    :: pftin(self%pftsz,self%kfromto(1):self%kfromto(2))
        complex(dp),         intent(inout) :: pftout(self%pftsz,self%kfromto(1):self%kfromto(2))
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
        complex(dp),         intent(in)    :: Mnum(self%pftsz,self%kfromto(1):self%kfromto(2))
        real(dp),            intent(inout) :: Mdenom(self%pftsz,self%kfromto(1):self%kfromto(2))
        complex(dp),         intent(inout) :: Mout(self%pftsz,self%kfromto(1):self%kfromto(2))
        logical  :: msk(self%pftsz)
        real(dp) :: avg, t
        integer  :: k
        do k = self%kfromto(1),self%kfromto(2)
            msk = Mdenom(:,k) > DSMALL
            avg = sum(Mdenom(:,k),mask=msk) / real(count(msk),dp)
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
        complex(dp),         intent(out) :: pftline(self%kfromto(1):self%kfromto(2))
        real(dp),            intent(out) :: ctf2line(self%kfromto(1):self%kfromto(2))
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
        complex(dp),         intent(in)    :: pft1(self%pftsz,self%kfromto(1):self%kfromto(2)), pft2(self%pftsz,self%kfromto(1):self%kfromto(2))
        integer,             intent(in)    :: n
        real(sp),            intent(inout) :: frc(1:n)
        real(dp) :: var1, var2, denom
        integer  :: k
        frc(1:self%kfromto(1)-1) = 0.999
        do k = self%kfromto(1), self%kfromto(2)
            var1  = sum(csq_fast(pft1(:,k)))
            var2  = sum(csq_fast(pft2(:,k)))
            if( (var1>DTINY) .and. (var2>DTINY) )then
                denom  = sqrt(var1*var2)
                frc(k) = real(sum(pft1(:,k)*conjg(pft2(:,k))) / denom, sp)
            else
                frc(k) = 0.0
            endif
        enddo
        if( self%kfromto(2) < n ) frc(self%kfromto(2)+1:) = 0.0
    end subroutine polar_cavger_calc_frc

end submodule simple_polarft_ops_restore
