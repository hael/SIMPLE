!@descr: flex_pca: the cross-fit-FSC driver context: setup, per-iteration ridge, paired records, teardown
submodule (simple_flex_probe_fit) simple_flex_probe_fit_crossfsc
use simple_core_module_api, only: dtiny, fdim, logfhandle
use simple_image,              only: image
use simple_flex_pca_crossfsc,  only: crossfsc_record, crossfsc_to_invtau2, COV_XFSC_FNAME
use simple_flex_pca_basis,     only: load_probe_state
implicit none
#include "simple_local_flags.inc"

contains

    !> Arm the cross-fit FSC ridge (arm 1, hard-wired; no switch selects another arm), reload the
    !! artifact when the ridge or the paired writer is live, and fail fast on contract violations.
    !! Master-only: workers pass l_master=.false. and stay inert.
    module subroutine xfsc_setup( ctx, params, kfr_ann, l_paired, l_master )
        type(xfsc_ctx_t),  intent(inout) :: ctx
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: kfr_ann(2)
        logical,           intent(in)    :: l_paired, l_master
        call ctx%new
        ctx%v_reg      = 1     ! the cross-fit FSC ridge (the internal e/o arm was the scaffolding)
        ctx%l_paired   = l_paired
        ! the paired master writes a record EVERY iteration; the single-fit engine writes none
        ctx%l_writer   = l_paired .and. l_master
        ctx%l_any      = ctx%l_writer .or. ctx%v_reg > 0
        ctx%l_loaded   = .false.
        ctx%filtsz     = max(1, fdim(params%box_crop) - 1)
        ! crossfsc low-resolution exemption index: the reslim_ind analog
        ctx%klo        = max(6, kfr_ann(1))
        if( .not. l_master )then
            ctx%l_any    = .false.
            ctx%l_writer = .false.
            return
        endif
        if( ctx%l_any )then
            ! restart-complete series reload (the load_probe_state idiom, applied to the artifact)
            call ctx%xf%load(ctx%l_loaded)
            if( ctx%l_loaded )then
                if( ctx%xf%box_crop /= params%box_crop .or. ctx%xf%filtsz /= ctx%filtsz ) &
                    &THROW_HARD('existing '//COV_XFSC_FNAME//' was written on a different lattice; delete it to restart the series')
                write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA XFSC artifact reloaded: ',ctx%xf%nrec, &
                    &' record(s)'
                call flush(logfhandle)
            endif
            if( ctx%v_reg > 0 .and. .not. ctx%l_writer .and. .not. ctx%l_loaded )then
                write(logfhandle,'(A)') '>>> FLEX_PCA XFSC ridge requested but no artifact exists &
                    &and no writer is active; it stays inert (the paired engine writes the records)'
                call flush(logfhandle)
            endif
        endif
    end subroutine xfsc_setup

    !> Build one fit's ridge for THIS iteration (fit%diag%xf_invtau2, applied by fit_iter_finish
    !! immediately before the coupled solves), from the latest record stamped <= t-1 -- the
    !! independence timing rule. Iteration 1 (no previous record) degrades to arm 0 for that
    !! iteration. Also sets l_xf_harvest when the writer needs the payloads.
    module subroutine xfsc_prep_iter( ctx, params, fit, it_eff, tag )
        type(xfsc_ctx_t),     intent(inout) :: ctx
        class(parameters),    intent(in)    :: params
        type(flex_probe_fit), intent(inout) :: fit
        integer,              intent(in)    :: it_eff
        character(len=*),     intent(in)    :: tag   !< '' single-fit; '  fit=A'/'  fit=B' paired
        logical :: l_fitb
        integer :: irec
        fit%diag%l_xf_harvest = ctx%l_writer
        if( allocated(fit%diag%xf_invtau2) ) deallocate(fit%diag%xf_invtau2)
        if( ctx%v_reg <= 0 ) return
        irec = ctx%xf%latest_upto(it_eff - 1)
        if( irec < 1 )then
            write(logfhandle,'(A,I0,A,A)') '>>> FLEX_PCA XFSC REG it=',it_eff, &
                &'  arm degraded to 0: no record stamped <= t-1',tag
        else
            ! which side of the record is THIS fit: the paired engine's fit%spec%id names it
            l_fitb = fit%spec%id == 2
            allocate(fit%diag%xf_invtau2(fit%model%ncomp,ctx%filtsz), source=0.0)
            call xfsc_build_invtau2(ctx, params, ctx%xf%recs(irec), l_fitb, it_eff, fit%diag%xf_invtau2)
            ! fit-level invtau2, added once per halfset by fit_iter_finish (the gold-standard
            ! precedent adds the shared curve into both halves' rho every iteration)
            write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> FLEX_PCA XFSC REG it=',it_eff, &
                &'  arm ',ctx%v_reg,' ridge added from record it=',ctx%xf%recs(irec)%it_eff,tag
        endif
        call flush(logfhandle)
    end subroutine xfsc_prep_iter

    !> Per-component, per-shell invtau2 from one paired record, in this fit's component indices:
    !! the matched cross-fit FSC with this fit's own H. Unmatched components take the kill branch;
    !! components beyond the record's rank get no ridge.
    subroutine xfsc_build_invtau2( ctx, params, rec, l_fitb, it_eff, invtau2 )
        type(xfsc_ctx_t),      intent(in)  :: ctx
        class(parameters),     intent(in)  :: params
        type(crossfsc_record), intent(in)  :: rec
        integer,               intent(in)  :: it_eff
        logical,               intent(in)  :: l_fitb
        real,                  intent(out) :: invtau2(:,:)   !< (ncomp, filtsz)
        real    :: fcur(ctx%filtsz), hcur(ctx%filtsz)
        integer :: qc, kc, kmatch_q, ncomp_fit, nun
        invtau2 = 0.0
        nun     = 0
        do qc = 1, size(invtau2,1)
            ncomp_fit = rec%ncomp_a
            if( l_fitb ) ncomp_fit = rec%ncomp_b
            if( qc > ncomp_fit )then
                ! component born after t-1: no curve, no H -> no ridge this iteration
                write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA XFSC REG: component ',qc, &
                    &' exceeds the previous record''s rank; no ridge applied to it'
                cycle
            endif
            if( l_fitb )then
                hcur = rec%h_b(:,qc)
            else
                hcur = rec%h_a(:,qc)
            endif
            ! locate this fit's component in the matched pairing
            kmatch_q = 0
            do kc = 1, rec%kmatch
                if( .not. l_fitb .and. rec%match_a(kc) == qc ) kmatch_q = kc
                if(       l_fitb .and. rec%match_b(kc) == qc ) kmatch_q = kc
            end do
            if( kmatch_q > 0 )then
                fcur = rec%fsc_cross(:,kmatch_q)
            else
                nun = nun + 1
                ! F<0 across the band -> kill-branch invtau2 everywhere in band
                fcur = -1.0
                write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA XFSC UNMATCHED component ',qc, &
                    &': no partner cleared the match floor; arm-1 ridge SUPPRESSES it (kill branch)'
            endif
            call crossfsc_to_invtau2(fcur, hcur, params%tau, ctx%klo, invtau2(qc,:))
        end do
        if( nun > 0 )then
            write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA XFSC REG it=',it_eff,'  ',nun, &
                &' unmatched component(s) this iteration (see lines above)'
        endif
        call flush(logfhandle)
    end subroutine xfsc_build_invtau2

    !> crossfsc state teardown (the artifact on disk is the persistent series); the per-fit
    !! payloads are freed by kill_probe_fit / scope exit.
    module subroutine xfsc_teardown( ctx )
        type(xfsc_ctx_t), intent(inout) :: ctx
        call ctx%kill
    end subroutine xfsc_teardown

    !> Append one crossfsc record per iteration, after both fits' tails: each fit's own e+o H, plus
    !! greedy signed |cos| matching on the delivered bases (prev_real) and each matched pair's
    !! cross-fit FSC with the sign resolved. Greedy, not varimax+Hungarian.
    module subroutine xfsc_paired_record( ctx, params, fits, it_eff )
        type(xfsc_ctx_t),     intent(inout) :: ctx
        class(parameters),    intent(in)    :: params
        type(flex_probe_fit), intent(inout) :: fits(2)
        integer,              intent(in)    :: it_eff
        type(crossfsc_record) :: rec
        type(image) :: imga, imgb
        real,     pointer     :: pra(:,:,:), prb(:,:,:)
        real(dp), allocatable :: cosm(:,:)
        real,     allocatable :: corrs(:)
        logical,  allocatable :: used_a(:), used_b(:)
        integer  :: nc_a, nc_b, km, q, p, k, best_q, best_p
        real(dp) :: nrm_a, nrm_b, c, best_c, cos_lead
        if( .not. (allocated(fits(1)%diag%xf_h_e) .and. allocated(fits(2)%diag%xf_h_e) .and. &
            &allocated(fits(1)%diag%xf_h_o) .and. allocated(fits(2)%diag%xf_h_o)) ) &
            &THROW_HARD('xfsc_paired_record: writer payloads missing (l_xf_harvest not set?)')
        if( .not. (allocated(fits(1)%history%prev_real) .and. allocated(fits(2)%history%prev_real)) ) &
            &THROW_HARD('xfsc_paired_record: realized bases (prev_real) missing')
        ! delivered ranks; the stashes were taken at entry rank >= delivered rank, truncated here
        nc_a = min(fits(1)%model%ncomp, size(fits(1)%diag%xf_h_e,2))
        nc_b = min(fits(2)%model%ncomp, size(fits(2)%diag%xf_h_e,2))
        km   = min(nc_a, nc_b)
        if( .not. ctx%l_loaded .and. ctx%xf%nrec == 0 )then
            ! file header, frozen at first record
            ctx%xf%box_crop = params%box_crop
            ctx%xf%filtsz   = ctx%filtsz
        endif
        rec%it_eff  = it_eff
        rec%ncomp_a = nc_a
        rec%ncomp_b = nc_b
        rec%kmatch  = km
        allocate(rec%match_a(km), rec%match_b(km), rec%fsc_cross(ctx%filtsz,km))
        allocate(rec%h_a(ctx%filtsz,nc_a), rec%h_b(ctx%filtsz,nc_b))
        ! a paired-engine fit stores its OWN e+o sampling sum: H is the per-shell
        ! mean of (rho_e + rho_o) diagonals = mean_e + mean_o (same voxel counts)
        rec%h_a = fits(1)%diag%xf_h_e(:,1:nc_a) + fits(1)%diag%xf_h_o(:,1:nc_a)
        rec%h_b = fits(2)%diag%xf_h_e(:,1:nc_b) + fits(2)%diag%xf_h_o(:,1:nc_b)
        ! ---- greedy signed |cos| matching on the realized bases ----
        allocate(cosm(nc_a,nc_b), used_a(nc_a), used_b(nc_b))
        used_a = .false.; used_b = .false.
        do p = 1, nc_b
            call fits(2)%history%prev_real(p)%get_rmat_ptr(prb)
            do q = 1, nc_a
                call fits(1)%history%prev_real(q)%get_rmat_ptr(pra)
                nrm_a = sqrt(sum(real(pra,dp)**2))
                nrm_b = sqrt(sum(real(prb,dp)**2))
                cosm(q,p) = sum(real(pra,dp)*real(prb,dp)) / max(nrm_a*nrm_b, DTINY)
            end do
        end do
        allocate(corrs(ctx%filtsz))
        cos_lead = 0.d0
        do k = 1, km
            best_c = -1.d0; best_q = 0; best_p = 0
            do q = 1, nc_a
                if( used_a(q) ) cycle
                do p = 1, nc_b
                    if( used_b(p) ) cycle
                    c = abs(cosm(q,p))
                    if( c > best_c )then
                        best_c = c; best_q = q; best_p = p
                    endif
                end do
            end do
            used_a(best_q) = .true.
            used_b(best_p) = .true.
            rec%match_a(k) = best_q
            rec%match_b(k) = best_p
            if( k == 1 ) cos_lead = best_c
            ! honest cross-fit FSC between A's component and sign * B's matched component
            call imga%copy(fits(1)%history%prev_real(best_q))
            call imgb%copy(fits(2)%history%prev_real(best_p))
            if( cosm(best_q,best_p) < 0.d0 ) call imgb%mul(-1.0)
            call imga%fft
            call imgb%fft
            call imga%fsc(imgb, corrs)
            rec%fsc_cross(:,k) = corrs
            call imga%kill
            call imgb%kill
        end do
        deallocate(cosm, used_a, used_b, corrs)
        write(logfhandle,'(A,I0,A,I0,A,F6.3)') '>>> FLEX_PCA XFSC PAIRED record it=', &
            &it_eff,'  matched pairs=',km,'  |cos| lead=',real(cos_lead)
        call ctx%xf%append(rec)
        call ctx%xf%write
        call rec%kill
    end subroutine xfsc_paired_record

end submodule simple_flex_probe_fit_crossfsc
