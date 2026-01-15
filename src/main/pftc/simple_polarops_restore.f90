submodule (simple_polarops) simple_polarops_restore
implicit none
#include "simple_local_flags.inc"
contains

    !>  \brief  Restores 2D class-averages only
    module subroutine polar_cavger_merge_eos_and_norm2D
        complex(dp) :: pft(pftsz,kfromto(1):kfromto(2))
        real(dp)    :: ctf2(pftsz,kfromto(1):kfromto(2))
        integer     :: icls, eo_pop(2), pop
        select case(trim(params_glob%ref_type))
        case('polar_cavg')
            ! all good
        case DEFAULT
            THROW_HARD('polar_cavger_merge_eos_and_norm2D only for 2D cavgs restoration')
        end select
        pfts_merg = DCMPLX_ZERO
        !$omp parallel do default(shared) schedule(static) proc_bind(close)&
        !$omp private(icls,eo_pop,pop,pft,ctf2)
        do icls = 1,ncls
            eo_pop = prev_eo_pops(:,icls) + eo_pops(:,icls)
            pop    = sum(eo_pop)
            if(pop == 0)then
                pfts_even(:,:,icls) = DCMPLX_ZERO
                pfts_odd(:,:,icls)  = DCMPLX_ZERO
                ctf2_even(:,:,icls) = 0.d0
                ctf2_odd(:,:,icls)  = 0.d0
            else
                if(pop > 1)then
                    pft  = pfts_even(:,:,icls) + pfts_odd(:,:,icls)
                    ctf2 = ctf2_even(:,:,icls) + ctf2_odd(:,:,icls)
                    call safe_norm(pft, ctf2, pfts_merg(:,:,icls))
                endif
                if(eo_pop(1) > 1)then
                    pft = pfts_even(:,:,icls)
                    call safe_norm(pft, ctf2_even(:,:,icls), pfts_even(:,:,icls))
                endif
                if(eo_pop(2) > 1)then
                    pft = pfts_odd(:,:,icls)
                    call safe_norm(pft, ctf2_odd(:,:,icls), pfts_odd(:,:,icls))
                endif
            endif
        end do
        !$omp end parallel do
    end subroutine polar_cavger_merge_eos_and_norm2D

    !>  \brief  Restores 3D slices
    module subroutine polar_cavger_merge_eos_and_norm( reforis, cl_weight )
        use simple_class_frcs
        type(oris),           intent(in) :: reforis
        real,       optional, intent(in) :: cl_weight
        type(class_frcs)   :: cavg2clin_frcs
        real,  allocatable :: cavg_clin_frcs(:,:,:)
        complex(dp) :: pfts_clin_even(pftsz,kfromto(1):kfromto(2),ncls)
        complex(dp) :: pfts_clin_odd(pftsz,kfromto(1):kfromto(2),ncls)
        real(dp)    :: ctf2_clin_even(pftsz,kfromto(1):kfromto(2),ncls)
        real(dp)    :: ctf2_clin_odd(pftsz,kfromto(1):kfromto(2),ncls)
        real(dp)    :: fsc(kfromto(1):kfromto(2)), clw
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
                call calc_comlin_contrib(reforis, build_glob%pgrpsyms,&
                &pfts_clin_even, pfts_clin_odd, ctf2_clin_even, ctf2_clin_odd)
            endif
        case('comlin_noself', 'comlin')
            ! Mirroring slices
            call mirror_slices(reforis, build_glob%pgrpsyms)
            ! Common-lines conribution
            call calc_comlin_contrib(reforis, build_glob%pgrpsyms,&
            &pfts_clin_even, pfts_clin_odd, ctf2_clin_even, ctf2_clin_odd)
        case DEFAULT
            THROW_HARD('Invalid REF_TYPE='//trim(params_glob%ref_type)//' in polar_cavger_merge_eos_and_norm')
        end select
        ! ML-regularization
        if( params_glob%l_ml_reg ) call add_invtausq2rho
        ! Restoration of references
        pfts_merg = DCMPLX_ZERO
        select case(trim(params_glob%ref_type))
            case('comlin_hybrid')
                if( clw > 1.d-6 )then
                    call calc_cavg_comlin_frcs(cavg2clin_frcs)
                    call restore_cavgs_comlins(clw)
                else
                    call mirr_and_calc_cavg_comlin_frcs(cavg2clin_frcs)
                    call polar_cavger_merge_eos_and_norm2D
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
            complex(dp) :: even(pftsz,kfromto(1):kfromto(2)), odd(pftsz,kfromto(1):kfromto(2))
            complex(dp) :: pft(pftsz,kfromto(1):kfromto(2))
            real(dp)    :: vare(kfromto(1):kfromto(2)), varo(kfromto(1):kfromto(2))
            real(dp)    :: sig2e(kfromto(1):kfromto(2)), sig2o(kfromto(1):kfromto(2))
            real(dp)    :: ssnr(kfromto(1):kfromto(2)), tau2(kfromto(1):kfromto(2))
            real(dp)    :: ctf2(pftsz,kfromto(1):kfromto(2)), cc, fudge, invtau2, s, a,b
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
                do icls = 1,ncls/2
                    ! e/o restoration
                    pft  = pfts_clin_even(:,:,icls)
                    ctf2 = ctf2_clin_even(:,:,icls)
                    call safe_norm(pft, ctf2, even)
                    pft  = pfts_clin_odd(:,:,icls)
                    ctf2 = ctf2_clin_odd(:,:,icls)
                    call safe_norm(pft, ctf2, odd)
                    ! FSC contribution
                    do k = kfromto(1),kfromto(2)
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
                do icls = 1,ncls/2
                    ! e/o restoration
                    pft  = pfts_even(:,:,icls) + pfts_clin_even(:,:,icls)
                    ctf2 = ctf2_even(:,:,icls) + ctf2_clin_even(:,:,icls)
                    call safe_norm(pft, ctf2, even)
                    pft  = pfts_odd(:,:,icls) + pfts_clin_odd(:,:,icls)
                    ctf2 = ctf2_odd(:,:,icls) + ctf2_clin_odd(:,:,icls)
                    call safe_norm(pft, ctf2, odd)
                    ! FSC contribution
                    do k = kfromto(1),kfromto(2)
                        fsc(k)   = fsc(k)   + real(sum(even(:,k) * conjg(odd(:,k))) ,dp)
                        vare(k)  = vare(k)  + real(sum(even(:,k) * conjg(even(:,k))),dp)
                        varo(k)  = varo(k)  + real(sum(odd(:,k)  * conjg(odd(:,k))) ,dp)
                        sig2e(k) = sig2e(k) + sum(ctf2_even(:,k,icls) + ctf2_clin_even(:,k,icls))
                        sig2o(k) = sig2o(k) + sum(ctf2_odd(:,k,icls)  + ctf2_clin_odd(:,k,icls))
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
            do k = kfromto(1),kfromto(2)
                cc = max(0.001d0,min(0.999d0,fsc(k)))
                ssnr(k) = fudge * cc / (1.d0 - cc)
            enddo
            ! Add Tau2 inverse to denominators
            ! because signal assumed infinite at very low resolution there is no addition
            kstart = max(6, calc_fourier_index(params_glob%hp, params_glob%box_crop, params_glob%smpd_crop))
            ! Even
            where( sig2e > DTINY )
                sig2e = real(ncls*pftsz,dp) / sig2e
            elsewhere
                sig2e = 0.d0
            end where
            tau2 = ssnr * sig2e
            if( DEBUG )then
                do k = kfromto(1),kfromto(2)
                    cc = max(0.001d0,min(0.999d0,fsc(k)))
                    s = cc / (1.d0 - cc)
                    print *, k, real(s), real(fsc(k)),&
                    &real(fudge**2*s/(fudge**2*s+1.d0)), tau2(k), 1.d0 / (fudge * tau2(k))
                enddo
                do k = kfromto(1),kfromto(2)
                    if( tau2(k) > DTINY )then
                        invtau2 = 1.d0 / (fudge * tau2(k))
                        a = sum(sqrt(real(pfts_clin_even(:,k,:)*conjg(pfts_clin_even(:,k,:))))) / real(pftsz*ncls,dp)
                        b =      sum(ctf2_clin_even(:,k,:))                                     / real(pftsz*ncls,dp)
                        print *,k,(a/(b+invtau2)) / (a/b)
                    endif
                enddo
            endif
            !$omp parallel do default(shared) schedule(static) proc_bind(close) private(k,icls,p,invtau2)
            do k = kstart,kfromto(2)
                if( tau2(k) > DTINY )then
                    ! CTF2 <- CTF2 + avgCTF2/(tau*SSNR)
                    invtau2 = 1.d0 / (fudge * tau2(k))
                    ctf2_clin_even(:,k,:) = ctf2_clin_even(:,k,:) + invtau2
                else
                    do icls = 1,ncls
                        do p = 1,pftsz
                            invtau2 = min(1.d3, 1.d3 * ctf2_clin_even(p,k,icls))
                            ctf2_clin_even(p,k,icls) = ctf2_clin_even(p,k,icls) + invtau2
                        enddo
                    enddo
                endif
            enddo
            !$omp end parallel do
            ! Odd
            where( sig2o > DTINY )
                sig2o = real(ncls*pftsz,dp) / sig2o
            elsewhere
                sig2o = 0.d0
            end where
            tau2 = ssnr * sig2o
            !$omp parallel do default(shared) schedule(static) proc_bind(close) private(k,icls,p,invtau2)
            do k = kstart,kfromto(2)
                if( tau2(k) > DTINY )then
                    invtau2 = 1.d0 / (fudge * tau2(k))
                    ctf2_clin_odd(:,k,:) = ctf2_clin_odd(:,k,:) + invtau2
                else
                    do icls = 1,ncls
                        do p = 1,pftsz
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
            complex(dp) :: pfte_backup(pftsz,kfromto(1):kfromto(2),ncls)
            complex(dp) :: pfto_backup(pftsz,kfromto(1):kfromto(2),ncls)
            real(dp)    :: ctf2e_backup(pftsz,kfromto(1):kfromto(2),ncls)
            real(dp)    :: ctf2o_backup(pftsz,kfromto(1):kfromto(2),ncls)
            ! backup classes
            pfte_backup(:,:,:)  = pfts_even(:,:,:); pfto_backup(:,:,:)  = pfts_odd(:,:,:)
            ctf2e_backup(:,:,:) = ctf2_even(:,:,:); ctf2o_backup(:,:,:) = ctf2_odd(:,:,:)
            ! mirror classes belonging to the same slice
            call mirror_slices(reforis, build_glob%pgrpsyms)
            ! calculate the per slice CLs
            call calc_comlin_contrib(reforis, build_glob%pgrpsyms,&
            &pfts_clin_even, pfts_clin_odd, ctf2_clin_even, ctf2_clin_odd)
            ! CLs vs. CLS FRCs
            call calc_cavg_comlin_frcs(frcs)
            ! restores classes
            pfts_even(:,:,:) = pfte_backup(:,:,:);  pfts_odd(:,:,:) = pfto_backup(:,:,:)
            ctf2_even(:,:,:) = ctf2e_backup(:,:,:); ctf2_odd(:,:,:) = ctf2o_backup(:,:,:)
        end subroutine mirr_and_calc_cavg_comlin_frcs

        ! Calculate CL contribution & FRCS from mirrored 2D classes
        ! Module arrays are untouched on exit
        subroutine calc_cavg_comlin_frcs( frcs )
            class(class_frcs), intent(inout) :: frcs
            real, allocatable :: frc(:)
            complex(dp)       :: cavg(pftsz,kfromto(1):kfromto(2)), clin(pftsz,kfromto(1):kfromto(2))
            complex(dp)       :: pft(pftsz,kfromto(1):kfromto(2))
            real(dp)          :: ctf2(pftsz,kfromto(1):kfromto(2))
            integer           :: icls, pop, nk
            call frcs%new(ncls, params_glob%box, params_glob%smpd, 1)
            nk = frcs%get_filtsz()
            allocate(frc(1:nk),source=0.)
            !$omp parallel do default(shared) schedule(static) proc_bind(close)&
            !$omp private(icls,pop,pft,ctf2,cavg,clin,frc)
            do icls = 1,ncls
                pop = sum(prev_eo_pops(:,icls) + eo_pops(:,icls))
                if( pop > 1 )then
                    ! cavg
                    pft  = pfts_even(:,:,icls) + pfts_odd(:,:,icls)
                    ctf2 = ctf2_even(:,:,icls) + ctf2_odd(:,:,icls)
                    call safe_norm(pft, ctf2, cavg)
                    ! comlin
                    pft  = pfts_clin_even(:,:,icls) + pfts_clin_odd(:,:,icls)
                    ctf2 = ctf2_clin_even(:,:,icls) + ctf2_clin_odd(:,:,icls)
                    call safe_norm(pft, ctf2, clin)
                    ! FRC
                    call calc_frc(cavg, clin, nk, frc)
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
            complex(dp) :: pft(pftsz,kfromto(1):kfromto(2))
            real(dp)    :: ctf2(pftsz,kfromto(1):kfromto(2))
            real        :: psi
            integer     :: iref, m
            logical     :: l_rotm
            if( .not.ref_space%isthere('mirr') )then
                THROW_HARD('Mirror index missing in reference search space')
            endif
            !$omp parallel do default(shared) proc_bind(close) private(iref,m,pft,ctf2,psi,l_rotm)
            do iref = 1,ncls/2
                m      = ref_space%get_int(iref,'mirr')
                psi    = abs(ref_space%get(m, 'psi'))
                l_rotm = (psi > 0.1) .and. (psi < 359.9)
                ! Fourier components
                if( l_rotm )then
                    call mirror_pft(pfts_even(:,:,m), pft)
                    pfts_even(:,:,iref) = pfts_even(:,:,iref) + conjg(pft)
                    call mirror_pft(pfts_odd(:,:,m), pft)
                    pfts_odd(:,:,iref)  = pfts_odd(:,:,iref)  + conjg(pft)
                    call mirror_pft(conjg(pfts_even(:,:,iref)), pfts_even(:,:,m))
                    call mirror_pft(conjg(pfts_odd(:,:,iref)),  pfts_odd(:,:,m))
                else
                    call mirror_pft(pfts_even(:,:,m), pft)
                    pfts_even(:,:,iref) = pfts_even(:,:,iref) + pft
                    call mirror_pft(pfts_odd(:,:,m), pft)
                    pfts_odd(:,:,iref)  = pfts_odd(:,:,iref)  + pft
                    call mirror_pft(pfts_even(:,:,iref), pfts_even(:,:,m))
                    call mirror_pft(pfts_odd(:,:,iref),  pfts_odd(:,:,m))
                endif
                ! CTF
                call mirror_ctf2(ctf2_even(:,:,m), ctf2)
                ctf2_even(:,:,iref) = ctf2_even(:,:,iref) + ctf2
                call mirror_ctf2(ctf2_odd(:,:,m),  ctf2)
                ctf2_odd(:,:,iref)  = ctf2_odd(:,:,iref)  + ctf2
                call mirror_ctf2(ctf2_even(:,:,iref), ctf2_even(:,:,m))
                call mirror_ctf2(ctf2_odd(:,:,iref),  ctf2_odd(:,:,m))
            enddo
            !$omp end parallel do
        end subroutine mirror_slices

        ! Restore common lines contributions only
        subroutine restore_comlins
            complex(dp) :: pft(pftsz,kfromto(1):kfromto(2)),pfte(pftsz,kfromto(1):kfromto(2)),pfto(pftsz,kfromto(1):kfromto(2))
            real(dp)    :: ctf2(pftsz,kfromto(1):kfromto(2)),ctf2e(pftsz,kfromto(1):kfromto(2)),ctf2o(pftsz,kfromto(1):kfromto(2))
            real        :: psi
            integer     :: icls, m
            logical     :: l_rotm
            !$omp parallel do default(shared) schedule(static) proc_bind(close)&
            !$omp private(icls,m,pft,ctf2,pfte,pfto,ctf2e,ctf2o,psi,l_rotm)
            do icls = 1,ncls/2
                ! already mirrored common-line contribution
                pfte  = pfts_clin_even(:,:,icls)
                pfto  = pfts_clin_odd(:,:,icls)
                ctf2e = ctf2_clin_even(:,:,icls)
                ctf2o = ctf2_clin_odd(:,:,icls)
                ! merged then e/o
                pft   = pfte  + pfto
                ctf2  = ctf2e + ctf2o
                call safe_norm(pft,  ctf2,  pfts_merg(:,:,icls))
                call safe_norm(pfte, ctf2e, pfts_even(:,:,icls))
                call safe_norm(pfto, ctf2o, pfts_odd(:,:,icls))
                ! mirroring the restored images
                m = reforis%get_int(icls,'mirr')
                call mirror_pft(pfts_merg(:,:,icls), pfts_merg(:,:,m))
                call mirror_pft(pfts_even(:,:,icls), pfts_even(:,:,m))
                call mirror_pft(pfts_odd(:,:,icls),  pfts_odd(:,:,m))
                psi    = abs(reforis%get(m, 'psi'))
                l_rotm = (psi > 0.1) .and. (psi < 359.9)
                if( l_rotm )then
                    pfts_merg(:,:,m) = conjg(pfts_merg(:,:,m))
                    pfts_even(:,:,m) = conjg(pfts_even(:,:,m))
                    pfts_odd(:,:,m)  = conjg(pfts_odd(:,:,m))
                endif
            enddo
            !$omp end parallel do
        end subroutine restore_comlins

        ! Restores slices (=cavgs + comlin)
        subroutine restore_cavgs_comlins( clw )
            real(dp), intent(in) :: clw
            complex(dp) :: pft(pftsz,kfromto(1):kfromto(2)),pfte(pftsz,kfromto(1):kfromto(2))
            complex(dp) :: pfto(pftsz,kfromto(1):kfromto(2))
            real(dp)    :: ctf2(pftsz,kfromto(1):kfromto(2)),ctf2e(pftsz,kfromto(1):kfromto(2))
            real(dp)    :: ctf2o(pftsz,kfromto(1):kfromto(2))
            real        :: psi
            integer     :: icls, m
            logical     :: l_rotm
            ! Restoration
            !$omp parallel do default(shared) schedule(static) proc_bind(close)&
            !$omp private(icls,m,pft,ctf2,pfte,pfto,ctf2e,ctf2o,psi,l_rotm)
            do icls = 1,ncls/2
                ! calculate sum of already mirrored class + common-line contribution
                pfte  = pfts_even(:,:,icls) + clw * pfts_clin_even(:,:,icls)
                pfto  = pfts_odd(:,:,icls)  + clw * pfts_clin_odd(:,:,icls)
                ctf2e = ctf2_even(:,:,icls) + clw * ctf2_clin_even(:,:,icls)
                ctf2o = ctf2_odd(:,:,icls)  + clw * ctf2_clin_odd(:,:,icls)
                ! merged then e/o
                pft   = pfte  + pfto
                ctf2  = ctf2e + ctf2o
                call safe_norm(pft,  ctf2,  pfts_merg(:,:,icls))
                call safe_norm(pfte, ctf2e, pfts_even(:,:,icls))
                call safe_norm(pfto, ctf2o, pfts_odd(:,:,icls))
                ! mirroring the restored images
                m = reforis%get_int(icls,'mirr')
                call mirror_pft(pfts_merg(:,:,icls), pfts_merg(:,:,m))
                call mirror_pft(pfts_even(:,:,icls), pfts_even(:,:,m))
                call mirror_pft(pfts_odd(:,:,icls),  pfts_odd(:,:,m))
                psi    = abs(reforis%get(m, 'psi'))
                l_rotm = (psi > 0.1) .and. (psi < 359.9)
                if( l_rotm )then
                    pfts_merg(:,:,m) = conjg(pfts_merg(:,:,m))
                    pfts_even(:,:,m) = conjg(pfts_even(:,:,m))
                    pfts_odd(:,:,m)  = conjg(pfts_odd(:,:,m))
                endif
            enddo
            !$omp end parallel do
        end subroutine restore_cavgs_comlins

    end subroutine polar_cavger_merge_eos_and_norm

    !>  \brief  calculates Fourier ring correlations
    module subroutine polar_cavger_calc_and_write_frcs_and_eoavg( fname, cline )
        use simple_cmdline, only: cmdline
        class(string), intent(in) :: fname
        type(cmdline), intent(in) :: cline
        complex(dp), allocatable :: prev_prefs(:,:,:)
        real,        allocatable :: frc(:)
        real(dp) :: ufrac_trec
        integer  :: icls, find, pop, filtsz
        filtsz = fdim(params_glob%box_crop) - 1
        allocate(frc(filtsz),source=0.)
        ! In case nspace/ncls has changed
        if( build_glob%clsfrcs%get_ncls() /= ncls )then
            call build_glob%clsfrcs%new(ncls, params_glob%box_crop, params_glob%smpd_crop, params_glob%nstates)
        endif
        !$omp parallel do default(shared) private(icls,frc,find,pop) schedule(static) proc_bind(close)
        do icls = 1,ncls
            if( l_comlin )then
                ! calculate FRC (pseudo-cavgs are never empty)
                call calc_frc(pfts_even(:,:,icls), pfts_odd(:,:,icls), filtsz, frc)
                call build_glob%clsfrcs%set_frc(icls, frc, 1)
                ! average low-resolution info between eo pairs to keep things in register
                find = min(kfromto(2), build_glob%clsfrcs%estimate_find_for_eoavg(icls, 1))
                if( find >= kfromto(1) )then
                    pfts_even(:,kfromto(1):find,icls) = pfts_merg(:,kfromto(1):find,icls)
                    pfts_odd(:,kfromto(1):find,icls)  = pfts_merg(:,kfromto(1):find,icls)
                endif
            else
                pop = sum(prev_eo_pops(:,icls) + eo_pops(:,icls))
                if( pop == 0 )then
                    frc = 0.
                    call build_glob%clsfrcs%set_frc(icls, frc, 1)
                else
                    ! calculate FRC
                    call calc_frc(pfts_even(:,:,icls), pfts_odd(:,:,icls), filtsz, frc)
                    call build_glob%clsfrcs%set_frc(icls, frc, 1)
                    ! average low-resolution info between eo pairs to keep things in register
                    find = min(kfromto(2), build_glob%clsfrcs%estimate_find_for_eoavg(icls, 1))
                    if( find >= kfromto(1) )then
                        pfts_even(:,kfromto(1):find,icls) = pfts_merg(:,kfromto(1):find,icls)
                        pfts_odd(:,kfromto(1):find,icls)  = pfts_merg(:,kfromto(1):find,icls)
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
            call read_pft_array(string(POLAR_REFS_FBODY)//'_even'//BIN_EXT, prev_prefs)
            !$omp parallel workshare proc_bind(close)
            pfts_even = ufrac_trec * pfts_even + (1.d0-ufrac_trec) * prev_prefs
            !$omp end parallel workshare
            call read_pft_array(string(POLAR_REFS_FBODY)//'_odd'//BIN_EXT, prev_prefs)
            !$omp parallel workshare proc_bind(close)
            pfts_odd  = ufrac_trec * pfts_odd   + (1.d0-ufrac_trec) * prev_prefs
            pfts_merg = 0.5d0 * (pfts_even + pfts_odd)
            !$omp end parallel workshare
            deallocate(prev_prefs)
        endif
    end subroutine polar_cavger_calc_and_write_frcs_and_eoavg

    !>  \brief  prepares one polar cluster centre image for alignment
    module subroutine polar_prep2Dref( icls, cavg, center, xyz )
        integer,                intent(in)    :: icls
        class(image), optional, intent(inout) :: cavg
        logical,      optional, intent(in)    :: center
        real,         optional, intent(out)   :: xyz(3)
        real, allocatable :: frc(:), filter(:), gaufilter(:)
        real    :: xy_cavg(2), crop_factor
        integer :: filtsz, k
        logical :: do_center
        ! Centering
        if( present(cavg).and.present(center).and.present(xyz) )then
            ! The three optional arguments are required
            do_center = .false.
            if( present(center) )then
                xyz       = 0.
                do_center = center.and.(params_glob%center .eq. 'yes')
            endif
            if( do_center )then
                crop_factor = real(params_glob%box_crop) / real(params_glob%box)
                select case(trim(params_glob%center_type))
                case('params')
                    ! offset from document in original pixel unit
                    call build_glob%spproj_field%calc_avg_offset2D(icls, xy_cavg)
                    if( arg(xy_cavg) < CENTHRESH )then
                        xyz = 0.
                    else if( arg(xy_cavg) > MAXCENTHRESH2D )then
                        xyz = [xy_cavg(1), xy_cavg(2), 0.]
                    else
                        xyz = cavg%calc_shiftcen_serial(params_glob%cenlp, params_glob%msk_crop)
                        xyz = xyz / crop_factor         ! scaled pixel unit
                        if( arg(xyz(1:2) - xy_cavg) > MAXCENTHRESH2D ) xyz = 0.
                    endif
                case('seg')
                    call calc_cavg_offset(cavg, params_glob%cenlp, params_glob%msk_crop, xy_cavg)
                    xyz(1:2) = xy_cavg / crop_factor    ! scaled pixel unit
                case('mass')
                    xyz = cavg%calc_shiftcen_serial(params_glob%cenlp, params_glob%msk_crop)
                    xyz = xyz / crop_factor             ! scaled pixel unit
                end select
            endif
        endif
        ! Filtering
        if( params_glob%l_ml_reg )then
            ! no filtering, not supported yet
        else
            ! FRC-based optimal filter
            filtsz = build_glob%clsfrcs%get_filtsz()
            allocate(frc(filtsz),filter(filtsz),source=0.)
            call build_glob%clsfrcs%frc_getter(icls, frc)
            if( any(frc > 0.143) )then
                if( params_glob%beta > 0.001 )then
                    call fsc2boostfilter(params_glob%beta, filtsz, frc, filter, merged=params_glob%l_lpset)
                else
                    call fsc2optlp_sub(filtsz, frc, filter, merged=params_glob%l_lpset)
                endif
            else
                filter = 1.0
            endif
            ! gaussian filter
            if(trim(params_glob%gauref).eq.'yes')then
                allocate(gaufilter(filtsz),source=0.)
                call gaussian_filter(params_glob%gaufreq, params_glob%smpd, params_glob%box, gaufilter)
                ! take the minimum of FRC-based & gaussian filters
                forall(k = 1:filtsz) filter(k) = min(filter(k), gaufilter(k))
                deallocate(gaufilter)
            endif
            call filterrefs(icls, filter)
            deallocate(frc,filter)
        endif
    end subroutine polar_prep2Dref

    !>  \brief  prepares one polar cluster centre image for alignment
    module subroutine polar_prep3Dref( icls )
        integer, intent(in) :: icls
        real, allocatable   :: gaufilter(:)
        integer :: filtsz
        ! Gaussian filter if needed
        if(trim(params_glob%gauref).eq.'yes')then
            filtsz = build_glob%clsfrcs%get_filtsz()
            allocate(gaufilter(filtsz),source=0.)
            call gaussian_filter(params_glob%gaufreq, params_glob%smpd, params_glob%box, gaufilter)
            call filterrefs(icls, gaufilter)
            deallocate(gaufilter)
        endif
    end subroutine polar_prep3Dref

    !>  \brief prepares a 2D class document with class index, resolution,
    !!         population, average correlation and weight
    module subroutine polar_cavger_gen2Dclassdoc( spproj )
        use simple_sp_project, only: sp_project
        class(sp_project), target, intent(inout) :: spproj
        class(oris), pointer :: ptcl_field, cls_field
        integer  :: pops(ncls)
        real(dp) :: corrs(ncls), ws(ncls)
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
            if( icls<1 .or. icls>params_glob%ncls )cycle
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
        call cls_field%new(ncls, is_ptcl=.false.)
        do icls=1,ncls
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

    ! PRIVATE UTILITIES

    ! Calculate common-lines contributions from all the slices
    subroutine calc_comlin_contrib( ref_space, symop, pfts_cl_even, pfts_cl_odd, ctf2_cl_even, ctf2_cl_odd )
        type(oris),       intent(in)    :: ref_space
        type(sym),        intent(in)    :: symop
        complex(kind=dp), intent(inout) :: pfts_cl_even(pftsz,kfromto(1):kfromto(2),ncls)
        complex(kind=dp), intent(inout) :: pfts_cl_odd(pftsz,kfromto(1):kfromto(2),ncls)
        real(kind=dp),    intent(inout) :: ctf2_cl_even(pftsz,kfromto(1):kfromto(2),ncls)
        real(kind=dp),    intent(inout) :: ctf2_cl_odd(pftsz,kfromto(1):kfromto(2),ncls)
        real, allocatable :: Rsym(:,:,:)
        complex(dp) :: cl_l(kfromto(1):kfromto(2)), cl_r(kfromto(1):kfromto(2)), cl_e(kfromto(1):kfromto(2))
        complex(dp) :: cl_o(kfromto(1):kfromto(2)), pft(pftsz,kfromto(1):kfromto(2))
        real(dp)    :: rl_l(kfromto(1):kfromto(2)), rl_r(kfromto(1):kfromto(2)), rl_e(kfromto(1):kfromto(2))
        real(dp)    :: rl_o(kfromto(1):kfromto(2)), ctf2(pftsz,kfromto(1):kfromto(2))
        real(dp)    :: w(kfromto(1):kfromto(2)), wl(kfromto(1):kfromto(2)), wr(kfromto(1):kfromto(2)), sumw(kfromto(1):kfromto(2))
        real        :: eulers(3),R(3,3,ncls),Rj(3,3),tRi(3,3),psi,drot,d
        integer     :: rotl, rotr, iref, jref, m, self_irot, targ_irot, isym, nsym
        logical     :: l_rotm
        if( .not.ref_space%isthere('mirr') )then
            THROW_HARD('Mirror index missing in reference search space')
        endif
        drot = pftc_glob%get_dang()
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
        do iref = 1,ncls
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
        do iref = 1,ncls/2
            tRi = transpose(R(:,:,iref))
            m   = ref_space%get_int(iref,'mirr')
            do jref = 1,ncls/2
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
                    call pftc_glob%gen_clin_weights(psi, rotl, rotr, wl, wr)
                    ! even
                    call get_line(jref, rotl, .true., cl_l, rl_l)
                    call get_line(jref, rotr, .true., cl_r, rl_r)
                    cl_e = wl*cl_l + wr*cl_r
                    rl_e = wl*rl_l + wr*rl_r
                    ! odd
                    call get_line(jref, rotl, .false., cl_l, rl_l)
                    call get_line(jref, rotr, .false., cl_r, rl_r)
                    cl_o = wl*cl_l + wr*cl_r
                    rl_o = wl*rl_l + wr*rl_r
                    ! in plane rotation angle of iref slice
                    psi = eulers(1)
                    ! get the weights, rotation indeces and extrapolate the common line
                    call pftc_glob%gen_clin_weights(psi, rotl, rotr, wl, wr)
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
        do iref = 1,ncls/2
            m      = ref_space%get_int(iref,'mirr')
            psi    = abs(ref_space%get(m, 'psi'))
            l_rotm = (psi>0.1) .and. (psi<359.9)
            if( l_rotm )then
                call mirror_pft(conjg(pfts_cl_even(:,:,iref)), pfts_cl_even(:,:,m))
                call mirror_pft(conjg(pfts_cl_odd(:,:,iref)),  pfts_cl_odd(:,:,m))
            else
                call mirror_pft(pfts_cl_even(:,:,iref), pfts_cl_even(:,:,m))
                call mirror_pft(pfts_cl_odd(:,:,iref),  pfts_cl_odd(:,:,m))
            endif
            call mirror_ctf2(ctf2_cl_even(:,:,iref), ctf2_cl_even(:,:,m))
            call mirror_ctf2(ctf2_cl_odd(:,:,iref),  ctf2_cl_odd(:,:,m))
        enddo
        !$omp end do
        !$omp end parallel

    contains

        ! Extrapolate cline, rline to pfts_clin and ctf2_clin
        subroutine extrapolate_line(ref, rot, weight, cle, clo, rle, rlo)
            integer,     intent(in) :: ref, rot
            real(dp),    intent(in) :: weight(kfromto(1):kfromto(2))
            complex(dp), intent(in) :: cle(kfromto(1):kfromto(2)), clo(kfromto(1):kfromto(2))
            real(dp),    intent(in) :: rle(kfromto(1):kfromto(2)), rlo(kfromto(1):kfromto(2))
            integer :: irot
            irot = rot
            if( irot < 1 )then
                irot = irot + pftsz
                pfts_cl_even(irot,:,ref) = pfts_cl_even(irot,:,ref) + weight * conjg(cle)
                pfts_cl_odd( irot,:,ref) = pfts_cl_odd( irot,:,ref) + weight * conjg(clo)
            elseif( irot > pftsz )then
                irot = irot - pftsz
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
    pure subroutine mirror_ctf2( ctf2in, ctf2out )
        real(dp), intent(in)    :: ctf2in(pftsz,kfromto(1):kfromto(2))
        real(dp), intent(inout) :: ctf2out(pftsz,kfromto(1):kfromto(2))
        integer :: i,j
        ctf2out(1,:) = ctf2in(1,:)
        do i = 2,pftsz/2
            j = pftsz-i+2
            ctf2out(i,:) = ctf2in(j,:)
            ctf2out(j,:) = ctf2in(i,:)
        enddo
        i = pftsz/2 + 1
        ctf2out(i,:) = ctf2in(i,:)
    end subroutine mirror_ctf2

    ! produces y-mirror of complex (reciprocal) matrix 
    pure subroutine mirror_pft( pftin, pftout )
        complex(dp), intent(in)    :: pftin(pftsz,kfromto(1):kfromto(2))
        complex(dp), intent(inout) :: pftout(pftsz,kfromto(1):kfromto(2))
        integer :: i,j
        pftout(1,:) = conjg(pftin(1,:))
        do i = 2,pftsz/2
            j = pftsz-i+2
            pftout(i,:) = pftin(j,:)
            pftout(j,:) = pftin(i,:)
        enddo
        i = pftsz/2 + 1
        pftout(i,:) = pftin(i,:)
    end subroutine mirror_pft

    ! Private utility
    pure subroutine safe_norm( Mnum, Mdenom, Mout )
        complex(dp), intent(in)    :: Mnum(pftsz,kfromto(1):kfromto(2))
        real(dp),    intent(inout) :: Mdenom(pftsz,kfromto(1):kfromto(2))
        complex(dp), intent(inout) :: Mout(pftsz,kfromto(1):kfromto(2))
        logical  :: msk(pftsz)
        real(dp) :: avg, t
        integer  :: k
        do k = kfromto(1),kfromto(2)
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
    pure subroutine get_line( ref, rot, even, pftline, ctf2line )
        integer,     intent(in)    :: ref, rot
        logical,     intent(in)    :: even
        complex(dp), intent(out)   :: pftline(kfromto(1):kfromto(2))
        real(dp),    intent(out)   :: ctf2line(kfromto(1):kfromto(2))
        integer :: irot
        if( rot >  nrots )then
            irot = rot - nrots
        else
            irot = rot
        endif
        if( even )then
            if( irot < 1 )then
                irot    = irot + pftsz
                pftline = conjg(pfts_even(irot,:,ref))
            elseif( irot > pftsz )then
                irot    = irot - pftsz
                pftline = conjg(pfts_even(irot,:,ref))
            else
                pftline = pfts_even(irot,:,ref)
            endif
            ctf2line = ctf2_even(irot,:,ref)
        else
            if( irot < 1 )then
                irot    = irot + pftsz
                pftline = conjg(pfts_odd(irot,:,ref))
            elseif( irot > pftsz )then
                irot    = irot - pftsz
                pftline = conjg(pfts_odd(irot,:,ref))
            else
                pftline = pfts_odd(irot,:,ref)
            endif
            ctf2line = ctf2_odd(irot,:,ref)
        endif
    end subroutine get_line

    !>  \brief  Filter references
    !! keep serial
    subroutine filterrefs( icls, filter )
        integer, intent(in) :: icls
        real,    intent(in) :: filter(:)
        integer :: n, k
        n = size(filter)
        if( n < kfromto(2) )then
            THROW_HARD('Incompatible filter size!; polar_cavger_filterref')
        endif
        do k = kfromto(1),kfromto(2)
            pfts_merg(:,k,icls) = filter(k) * pfts_merg(:,k,icls)
            pfts_even(:,k,icls) = filter(k) * pfts_even(:,k,icls)
            pfts_odd(:,k,icls)  = filter(k) * pfts_odd(:,k,icls)
        enddo
    end subroutine filterrefs

    subroutine calc_frc( pft1, pft2, n, frc )
        complex(dp), intent(in)    :: pft1(pftsz,kfromto(1):kfromto(2)), pft2(pftsz,kfromto(1):kfromto(2))
        integer,     intent(in)    :: n
        real(sp),    intent(inout) :: frc(1:n)
        real(dp) :: var1, var2, denom
        integer  :: k
        frc(1:kfromto(1)-1) = 0.999
        do k = kfromto(1), kfromto(2)
            var1  = sum(csq_fast(pft1(:,k)))
            var2  = sum(csq_fast(pft2(:,k)))
            if( (var1>DTINY) .and. (var2>DTINY) )then
                denom  = sqrt(var1*var2)
                frc(k) = real(sum(pft1(:,k)*conjg(pft2(:,k))) / denom, sp)
            else
                frc(k) = 0.0
            endif
        enddo
        if( kfromto(2) < n ) frc(kfromto(2)+1:) = 0.0
    end subroutine calc_frc

end submodule simple_polarops_restore
