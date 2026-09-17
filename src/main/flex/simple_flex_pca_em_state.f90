!@descr: flex_pca EM: probe_fit_t lifecycle (the paired-engine per-fit state hoist)
submodule (simple_flex_pca_em) simple_flex_pca_em_state
use simple_flex_pca_polar, only: polar_grid_kill
implicit none
#include "simple_local_flags.inc"

contains

    !> Construct a fit shell: identity, file namespaces and the fit's particle selection.
    !! Model handles (basis/eigvals/sig2), the stage subsample and all iteration state are
    !! populated by the driver, mirroring the single-fit initialisation order.
    module subroutine new_probe_fit( fit, id, fprefix, meta_fname, pinds, nptcls )
        type(probe_fit_t), intent(inout) :: fit
        integer,           intent(in)    :: id
        character(len=*),  intent(in)    :: fprefix, meta_fname
        integer,           intent(in)    :: pinds(:)
        integer,           intent(in)    :: nptcls
        call kill_probe_fit(fit)
        fit%id         = id
        fit%fprefix    = fprefix
        fit%meta_fname = meta_fname
        allocate(fit%pinds(nptcls), source=pinds(:nptcls))
        fit%nptcls     = nptcls
    end subroutine new_probe_fit

    !> Free everything a fit may hold, in dependency order: per-iteration fields first (in case a
    !! crash or early exit left them allocated), then the polar bank, the MCFA state (freed ONLY
    !! here or by the rank-change resize -- see the type's lifecycle contract), the previous-basis
    !! images, the model handles, and finally the selection. Mirrors the
    !! single-fit teardown at the end of probe_subspace_iteration plus the worker early-exit one.
    module subroutine kill_probe_fit( fit )
        type(probe_fit_t), intent(inout) :: fit
        integer :: q, ithr
        ! ---- per-iteration M-step accumulators and batch scratch ----
        if( allocated(fit%Yeven) )then
            do q = 1, size(fit%Yeven)
                call fit%Yeven(q)%dealloc_rho; call fit%Yeven(q)%kill
            end do
            deallocate(fit%Yeven)
        endif
        if( allocated(fit%Yodd) )then
            do q = 1, size(fit%Yodd)
                call fit%Yodd(q)%dealloc_rho; call fit%Yodd(q)%kill
            end do
            deallocate(fit%Yodd)
        endif
        if( allocated(fit%rho_e)  ) deallocate(fit%rho_e)
        if( allocated(fit%rho_o)  ) deallocate(fit%rho_o)
        if( allocated(fit%prior)  ) deallocate(fit%prior)
        if( allocated(fit%Gth)    ) deallocate(fit%Gth)
        if( allocated(fit%Ath)    ) deallocate(fit%Ath)
        if( allocated(fit%bth)    ) deallocate(fit%bth)
        if( allocated(fit%cth)    ) deallocate(fit%cth)
        if( allocated(fit%zth)    ) deallocate(fit%zth)
        if( allocated(fit%Ainvth) ) deallocate(fit%Ainvth)
        if( allocated(fit%Acpth)  ) deallocate(fit%Acpth)
        if( allocated(fit%hth)    ) deallocate(fit%hth)
        if( allocated(fit%nll_thr)) deallocate(fit%nll_thr)
        if( allocated(fit%gam_dbg)) deallocate(fit%gam_dbg)
        if( allocated(fit%zbatch) ) deallocate(fit%zbatch)
        if( allocated(fit%dens)   ) deallocate(fit%dens)
        if( allocated(fit%valid)  ) deallocate(fit%valid)
        if( allocated(fit%valid_e)) deallocate(fit%valid_e)
        if( allocated(fit%valid_o)) deallocate(fit%valid_o)
        if( allocated(fit%gam_thr)) deallocate(fit%gam_thr)
        if( allocated(fit%gam_acc)) deallocate(fit%gam_acc)
        if( allocated(fit%gam_sum)) deallocate(fit%gam_sum)
        if( allocated(fit%nval_thr)) deallocate(fit%nval_thr)
        if( allocated(fit%mean_fpl) )then
            do ithr = 1, size(fit%mean_fpl)
                call cleanup_plane(fit%mean_fpl(ithr))
            end do
            deallocate(fit%mean_fpl)
        endif
        if( allocated(fit%basis_fpls) )then
            do ithr = 1, size(fit%basis_fpls,2)
                do q = 1, size(fit%basis_fpls,1)
                    call cleanup_plane(fit%basis_fpls(q,ithr))
                end do
            end do
            deallocate(fit%basis_fpls)
        endif
        if( allocated(fit%sec_proj_thr) ) deallocate(fit%sec_proj_thr)
        if( allocated(fit%sec_gram_thr) ) deallocate(fit%sec_gram_thr)
        ! ---- cross-fit-FSC per-fit payloads (the artifact on disk is the persistent series) ----
        fit%l_xf_harvest = .false.
        if( allocated(fit%xf_h_e)    ) deallocate(fit%xf_h_e)
        if( allocated(fit%xf_h_o)    ) deallocate(fit%xf_h_o)
        if( allocated(fit%xf_cnt)    ) deallocate(fit%xf_cnt)
        if( allocated(fit%xf_fscq)   ) deallocate(fit%xf_fscq)
        if( allocated(fit%xf_gam)    ) deallocate(fit%xf_gam)
        if( allocated(fit%xf_invtau2)) deallocate(fit%xf_invtau2)
        ! ---- paired-merge stash ----
        fit%l_mg_stash = .false.
        fit%mg_ncomp   = 0
        fit%mg_npairs  = 0
        if( allocated(fit%mg_ye)  ) deallocate(fit%mg_ye)
        if( allocated(fit%mg_yo)  ) deallocate(fit%mg_yo)
        if( allocated(fit%mg_rhe) ) deallocate(fit%mg_rhe)
        if( allocated(fit%mg_rho) ) deallocate(fit%mg_rho)
        if( allocated(fit%mg_prev) )then
            do q = 1, size(fit%mg_prev)
                call fit%mg_prev(q)%kill
            end do
            deallocate(fit%mg_prev)
        endif
        ! ---- polar E-step bank (grid via polar_grid_kill) ----
        call polar_grid_kill(fit%pg_es)
        if( allocated(fit%UsallE) ) deallocate(fit%UsallE)
        if( allocated(fit%CfE)    ) deallocate(fit%CfE)
        if( allocated(fit%Cm0E)   ) deallocate(fit%Cm0E)
        if( allocated(fit%c00E)   ) deallocate(fit%c00E)
        if( allocated(fit%UbankE) ) deallocate(fit%UbankE)
        if( allocated(fit%CspE)   ) deallocate(fit%CspE)
        if( allocated(fit%xws_es) ) deallocate(fit%xws_es)
        if( allocated(fit%wr_es)  ) deallocate(fit%wr_es)
        if( allocated(fit%wrd_es) ) deallocate(fit%wrd_es)
        if( allocated(fit%Reb_es) ) deallocate(fit%Reb_es)
        if( allocated(fit%rmatb_es) ) deallocate(fit%rmatb_es)
        if( allocated(fit%nrmb_es)  ) deallocate(fit%nrmb_es)
        if( allocated(fit%dir_es) ) deallocate(fit%dir_es)
        if( allocated(fit%cae)    ) deallocate(fit%cae)
        if( allocated(fit%sae)    ) deallocate(fit%sae)
        if( allocated(fit%dused_es) ) deallocate(fit%dused_es)
        if( allocated(fit%hex_es) ) deallocate(fit%hex_es)
        if( allocated(fit%kex_es) ) deallocate(fit%kex_es)
        ! ---- MCFA mixture state (see the type's lifecycle contract) ----
        if( allocated(fit%mix_xi)   ) deallocate(fit%mix_xi)
        if( allocated(fit%mix_Om)   ) deallocate(fit%mix_Om)
        if( allocated(fit%mix_Ominv)) deallocate(fit%mix_Ominv)
        if( allocated(fit%mix_pi)   ) deallocate(fit%mix_pi)
        if( allocated(fit%mix_Omxi) ) deallocate(fit%mix_Omxi)
        if( allocated(fit%mix_xiOx) ) deallocate(fit%mix_xiOx)
        if( allocated(fit%mix_lpi)  ) deallocate(fit%mix_lpi)
        if( allocated(fit%rhs0th)   ) deallocate(fit%rhs0th)
        if( allocated(fit%mkth)     ) deallocate(fit%mkth)
        if( allocated(fit%lwth)     ) deallocate(fit%lwth)
        if( allocated(fit%rkth)     ) deallocate(fit%rkth)
        if( allocated(fit%mxa_sr)   ) deallocate(fit%mxa_sr)
        if( allocated(fit%mxa_sm)   ) deallocate(fit%mxa_sm)
        if( allocated(fit%mxa_smm)  ) deallocate(fit%mxa_smm)
        if( allocated(fit%mxa_sainv)) deallocate(fit%mxa_sainv)
        if( allocated(fit%dm_sr)    ) deallocate(fit%dm_sr)
        if( allocated(fit%dm_sm)    ) deallocate(fit%dm_sm)
        if( allocated(fit%dm_smm)   ) deallocate(fit%dm_smm)
        if( allocated(fit%dm_sai)   ) deallocate(fit%dm_sai)
        if( allocated(fit%dm_z)     ) deallocate(fit%dm_z)
        fit%dm_nz = 0
        ! ---- previous-basis images ----
        if( allocated(fit%prev_real) )then
            do q = 1, size(fit%prev_real)
                call fit%prev_real(q)%kill
            end do
            deallocate(fit%prev_real)
        endif
        ! ---- model handles ----
        if( allocated(fit%basis_recs) )then
            do q = 1, size(fit%basis_recs)
                call fit%basis_recs(q)%dealloc_rho; call fit%basis_recs(q)%kill
            end do
            deallocate(fit%basis_recs)
        endif
        if( allocated(fit%eigvals) ) deallocate(fit%eigvals)
        call fit%mean_rec%dealloc_rho
        call fit%mean_rec%kill
        if( allocated(fit%z) ) deallocate(fit%z)
        ! ---- selection ----
        if( allocated(fit%ppinds) ) deallocate(fit%ppinds)
        if( allocated(fit%pinds)  ) deallocate(fit%pinds)
        call fit%fprefix%kill
        call fit%meta_fname%kill
        ! PCG M-step accumulators (kernel and rhs, per half, plus the pair-merge stash) and the operator
        if( allocated(fit%kacc_e) ) deallocate(fit%kacc_e)
        if( allocated(fit%kacc_o) ) deallocate(fit%kacc_o)
        if( allocated(fit%kpk_e)  ) deallocate(fit%kpk_e)
        if( allocated(fit%kpk_o)  ) deallocate(fit%kpk_o)
        if( allocated(fit%mg_kpe) ) deallocate(fit%mg_kpe)
        if( allocated(fit%mg_kpo) ) deallocate(fit%mg_kpo)
        if( allocated(fit%racc_e) ) deallocate(fit%racc_e)
        if( allocated(fit%racc_o) ) deallocate(fit%racc_o)
        if( allocated(fit%rpk_e)  ) deallocate(fit%rpk_e)
        if( allocated(fit%rpk_o)  ) deallocate(fit%rpk_o)
        if( allocated(fit%mg_rpe) ) deallocate(fit%mg_rpe)
        if( allocated(fit%mg_rpo) ) deallocate(fit%mg_rpo)
        call fit%pcg%kill
    end subroutine kill_probe_fit

end submodule simple_flex_pca_em_state
