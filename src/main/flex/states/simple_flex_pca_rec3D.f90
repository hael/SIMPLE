!@descr: flex_pca: the state-reconstruction service -- selects the gridding or PCG backend, runs the state stage, and hands every state's maps to the common 3-D delivery
module simple_flex_pca_rec3D
use simple_core_module_api
use simple_builder,    only: builder
use simple_parameters, only: parameters
use simple_flex_pca_rounds,          only: flex_pca_rounds
use simple_flex_pca_stages,          only: flex_stage_request, PCA_STAGE_STATES
use simple_flex_pca_run_types,       only: flex_run_settings
use simple_flex_pca_state_parts,     only: write_state_weights_round, read_state_weights_round
use simple_flex_pca_states_backend,  only: flex_states_backend, flex_state_maps, flex_rec_box, flex_rec_smpd
use simple_flex_pca_states_gridding, only: flex_states_gridding
use simple_flex_pca_states_pcg,      only: flex_states_pcg
use simple_flex_pca_state_delivery,  only: flex_state_delivery
implicit none

public :: reconstruct_flex_weighted_states
public :: read_state_weights_round
public :: flex_rec_box, flex_rec_smpd
private
#include "simple_local_flags.inc"

contains

    !> Kernel-weighted state reconstruction: each state is a weighted backprojection of all particles.
    !! With outvol_even/outvol_odd present the combined, even and odd maps come from one halfset-split
    !! round. The master ships the weight table and fans the particle range out; every backend then
    !! runs the same way (begin, accumulate or fold the parts, finalize per state, kill) and the
    !! common delivery applies the backend's declared policy.
    subroutine reconstruct_flex_weighted_states( params, cfg, build, pinds, state_weights, nstates, &
        &floor_rho, outvol_even, outvol_odd, split_eo , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters), intent(inout) :: params
        type(flex_run_settings), intent(in)    :: cfg
        class(builder),    intent(inout) :: build
        integer,           intent(in)    :: pinds(:), nstates
        real,              intent(in)    :: state_weights(:,:)
        ! shellwise rho floor before the divide (flex_pca opts in; the trial half maps too)
        logical,                intent(in) :: floor_rho
        type(string), optional, intent(in) :: outvol_even, outvol_odd
        !! worker-side halfset split flag, read from the round-weights table (a worker has no output
        !! names); both routes must set l_fuse identically or master and workers disagree on the halves
        logical,      optional, intent(in) :: split_eo
        class(flex_states_backend), allocatable :: be
        type(flex_state_delivery) :: delivery
        type(flex_state_maps)     :: maps
        logical :: l_floor_rho, l_fuse
        integer :: state, box_rec
        real    :: smpd_rec
        l_floor_rho = floor_rho
        l_fuse = present(outvol_even) .and. present(outvol_odd)
        if( present(split_eo) ) l_fuse = l_fuse .or. split_eo
        if( size(pinds)<1 .or. nstates<1 ) THROW_HARD('invalid flex weighted state reconstruction dimensions')
        if( any(shape(state_weights)/=[size(pinds),nstates]) ) THROW_HARD('flex weighted state table mismatch')
        box_rec  = flex_rec_box(params)
        smpd_rec = flex_rec_smpd(params)
        if( box_rec /= params%box_crop )then
            write(logfhandle,'(A,I0,A,F6.3,A,I0,A,F6.3,A)') '>>> FLEX STATE RECONSTRUCTION decoupled box: rec box=',box_rec, &
                &' smpd=',smpd_rec,' A (covariance box=',params%box_crop,' smpd=',params%smpd_crop,' A)'
        endif
        ! rec_states_backend=pcg: the same weighted least-squares problems on reconstructor_pcg with the
        ! support inside the solve; the weights round and the stage fan-out are shared with the gridding
        ! path. This is deliberately NOT rec_backend: the M-step and the state maps are separate
        ! decisions (doc/refactoring_notes/flex_pca_branch_reconciliation_2026_09_15.md 4.3 -- PCG wins
        ! the basis, gridding wins the state maps), so the default here is gridding even under
        ! rec_backend=pcg.
        if( trim(params%rec_states_backend) == 'pcg' )then
            write(logfhandle,'(A)') '>>> FLEX STATE RECONSTRUCTION: kernel PCG backend (rec_states_backend=pcg)'
            call flush(logfhandle)
            allocate(flex_states_pcg :: be)
        else
            allocate(flex_states_gridding :: be)
        endif
        ! Distributed: the master ships the weight table and fans the particle range out, then folds
        ! the parts; every nonlinear finalisation runs once on the global sums.
        if( rounds%is_master() )then
            call write_state_weights_round(pinds, state_weights, size(pinds), nstates, l_fuse)
            call rounds%run_stage(params, flex_stage_request(stage=PCA_STAGE_STATES, label='state reconstruction'))
        endif
        call be%begin(params, build, rounds, pinds, state_weights, nstates, l_fuse, l_floor_rho, box_rec, smpd_rec)
        if( rounds%is_worker() )then
            call be%accumulate_local_or_write_part(params, build, rounds)
            call be%kill
            deallocate(be)
            return
        endif
        if( rounds%distributed() )then
            call be%fold_parts(params, build, rounds)
        else
            call be%accumulate_local_or_write_part(params, build, rounds)
        endif
        call delivery%new(params, be%delivery_policy(cfg), nstates, l_fuse, box_rec, smpd_rec, &
            &outvol_even=outvol_even, outvol_odd=outvol_odd)
        do state = 1, nstates
            call be%finalize_maps(params, build, rounds, state, maps)
            call delivery%deliver(params, build, state, maps)
            call maps%kill
        end do
        call delivery%finish(params, build)
        call delivery%kill
        call be%kill
        deallocate(be)
    end subroutine reconstruct_flex_weighted_states

end module simple_flex_pca_rec3D
