!@descr: white-box seam of the polar shift search for its testers: the Levenberg-Marquardt callback's terms
submodule (simple_pftc_shsrch_grad) simple_pftc_shsrch_grad_tester
implicit none

contains

    module subroutine shsrch_lm_terms( srch, vec, f, grad, curv, irot )
        class(pftc_shsrch_grad), intent(inout)              :: srch
        real(dp),                intent(inout)              :: vec(:)
        real(dp),                intent(out)                :: f, grad(:), curv(:,:)
        integer,                 optional,       intent(in) :: irot
        if( srch%search_mode == SHSRCH_SHIFT )then
            if( present(irot) ) srch%cur_inpl_idx = irot
            call srch%b_ptr%pftc%load_ref_aln_8(srch%reference, srch%particle)
            srch%shift_curv_irot = 0
        else
            call srch%b_ptr%pftc%load_ref_joint_curv_8(srch%reference, srch%particle)
        endif
        ! through the spec's pointer, as the optimiser reaches it
        call srch%ospec%fdfcurvfun_8(srch, vec, f, grad, curv, size(vec))
    end subroutine shsrch_lm_terms

    module subroutine shsrch_callback_profile( srch, callbacks, angle_changes, reevals )
        class(pftc_shsrch_grad), intent(in)  :: srch
        integer(int64),          intent(out) :: callbacks, angle_changes, reevals
        callbacks     = srch%profile_callbacks
        angle_changes = srch%profile_angle_changes
        reevals       = srch%profile_reevals
    end subroutine shsrch_callback_profile

end submodule simple_pftc_shsrch_grad_tester
