!@descr: flex_pca probe-fit parent; E-step, update, cross-FSC and engine bodies live in four submodule files
module simple_flex_probe_fit
use simple_core_module_api, only: dp, fplane_type, ori, simple_exception, string
use simple_flex_pca_records,   only: flex_fit_model, flex_selection
use simple_builder,            only: builder
use simple_parameters,         only: parameters
use simple_reconstructor,      only: reconstructor
use simple_ori,                only: ori
use simple_flex_pca_rounds,    only: flex_pca_rounds
use simple_flex_pca_fit_types, only: flex_fit, flex_probe_part, xfsc_ctx_t
use simple_flex_pca_planes,    only: flex_plane_store
use simple_flex_pca_pcg,       only: flex_pcg_environment
implicit none
private
#include "simple_local_flags.inc"

public :: flex_probe_fit, flex_mean_ref, probe_subspace_iteration, fit_engine_iterate
! Narrow seams used by the adjacent tester.  The fixtures live in
! simple_flex_pca_tester so test-only code is not part of these production submodules.
public :: fit_estep_former_polar, polar_ring_gram, polar_ring_selfpower
public :: write_probe_part, reduce_probe_parts, xfsc_paired_record

!> Mean principal-angle cosine vs the previous basis at which a rank-1 fit stops early. Higher ranks
!! run the full n_probe_iters budget (the paired merge needs the last two iterations' frames).
real(dp), parameter :: COV_PROBE_CONV     = 0.999999d0
!> mixture width of the MCFA E-step: the deconvolution picks the macro-clusters downstream, so this
!> is only the E-step's flexibility budget
integer,  parameter :: COV_EM_MIX         = 16
!> consensus resolution shells deflated out of every basis volume each M-step (the background and
!> dilation templates are added on top)
integer,  parameter :: COV_EM_DEFLATE     = 4
integer,  parameter :: PROBE_PART_VERSION = 13
integer,  parameter :: MIX_ZSUB_MAX       = 2000

type, extends(flex_fit) :: flex_probe_fit
  contains
    ! ---- the E-step contract: selection at stage begin, one batch value per batch ----
    procedure, pass(fit) :: estep_begin_stage  => fit_estep_begin_stage   !< formulation and policy, once per stage
    procedure, pass(fit) :: estep_bank_prepare => fit_estep_bank_prepare  !< per-iteration bank (polar), built at the first batch
    procedure, pass(fit) :: estep_batch_begin  => fit_estep_batch_begin   !< reset the batch rows
    procedure, pass(fit) :: estep_particle     => fit_estep_particle      !< the same for one particle (interleaved owners)
    ! ---- the M-step contract ----
    procedure, pass(fit) :: mstep_insert_batch => fit_batch_insert        !< accumulate one batch into Y/rho (and the PCG kernels)
    procedure, pass(fit) :: iter_begin  => fit_iter_begin
    procedure, pass(fit) :: iter_reduce => fit_iter_reduce
    procedure, pass(fit) :: iter_finish => fit_iter_finish
    procedure, pass(fit) :: merge_stash => probe_fit_merge_stash
end type flex_probe_fit

type :: flex_mean_ref
    type(reconstructor), pointer :: p => null()
end type flex_mean_ref

interface

    module subroutine fit_polar_bank_build( build, fit, mean_rec, fpl1, nthr )
        type(builder),        intent(inout) :: build
        type(flex_probe_fit), intent(inout) :: fit
        type(reconstructor),  intent(inout) :: mean_rec
        type(fplane_type),    intent(in)    :: fpl1
        integer,              intent(in)    :: nthr
    end subroutine fit_polar_bank_build

    module subroutine fit_estep_former_polar( fit, mean_rec, o, fpl, row, ithr, a, e_mm, myv )
        type(flex_probe_fit), intent(inout) :: fit
        type(reconstructor),  intent(inout) :: mean_rec
        class(ori),           intent(inout) :: o
        type(fplane_type),    intent(inout) :: fpl
        integer,              intent(in)    :: row, ithr
        real(dp),             intent(out)   :: a, e_mm, myv
    end subroutine fit_estep_former_polar

    module subroutine fit_estep_solve_stats( fit, fpl, i, row, ithr, a )
        type(flex_probe_fit), intent(inout) :: fit
        type(fplane_type),    intent(inout) :: fpl
        integer,              intent(in)    :: i, row, ithr
        real(dp),             intent(in)    :: a
    end subroutine fit_estep_solve_stats

    module subroutine fit_estep_bank_prepare( fit, params, build, mean_rec, fpl1, nthr, it_eff, tag )
        class(flex_probe_fit), intent(inout) :: fit
        class(parameters),     intent(inout) :: params
        type(builder),         intent(inout) :: build
        type(reconstructor),   intent(inout) :: mean_rec
        type(fplane_type),     intent(in)    :: fpl1
        integer,               intent(in)    :: nthr, it_eff
        character(len=*),      intent(in)    :: tag
    end subroutine fit_estep_bank_prepare

    module subroutine fit_estep_batch_begin( fit, batchsz )
        class(flex_probe_fit), intent(inout) :: fit
        integer,               intent(in)    :: batchsz
    end subroutine fit_estep_batch_begin

    module subroutine fit_estep_particle( fit, mean_rec, o, fpl, row, i, ithr )
        class(flex_probe_fit), intent(inout) :: fit
        type(reconstructor),   intent(inout) :: mean_rec
        class(ori),            intent(inout) :: o
        type(fplane_type),     intent(inout) :: fpl
        integer,               intent(in)    :: row, i, ithr
    end subroutine fit_estep_particle

    module subroutine fit_batch_insert( build, fit, orientations, fpls, eo, batchsz )
        type(builder),         intent(inout) :: build
        class(flex_probe_fit), intent(inout) :: fit
        type(ori),             intent(inout) :: orientations(:)
        type(fplane_type),     intent(inout) :: fpls(:)
        integer,               intent(in)    :: eo(:), batchsz
    end subroutine fit_batch_insert

    module subroutine fit_iter_reduce( fit, it_eff, nthr , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        class(flex_probe_fit),  intent(inout) :: fit
        integer,                intent(in)    :: it_eff, nthr
    end subroutine fit_iter_reduce

    module subroutine fit_estep_pass( params, build, plane_store, fits, means, nfits, it_eff, nthr )
        class(parameters),       intent(inout) :: params
        type(builder),           intent(inout) :: build
        class(flex_plane_store), intent(inout) :: plane_store
        integer,                 intent(in)    :: nfits
        type(flex_probe_fit),    intent(inout) :: fits(nfits)
        type(flex_mean_ref),     intent(in)    :: means(nfits)
        integer,                 intent(in)    :: it_eff, nthr
    end subroutine fit_estep_pass

    logical module function cov_polar_enabled()
    end function cov_polar_enabled

    integer module function cov_polar_ndir( nptcls )
        integer, intent(in) :: nptcls
    end function cov_polar_ndir

    module subroutine polar_ring_gram( Us, ldu, ncomp, row0, nrow, Csp, Cout, Mout )
        integer,  intent(in)    :: ldu, ncomp, row0, nrow
        real,     intent(in)    :: Us(ldu,0:ncomp)
        real,     intent(inout) :: Csp(0:ncomp,0:ncomp)      !< caller-owned scratch
        real(dp), intent(out)   :: Cout(ncomp*ncomp), Mout(ncomp)
    end subroutine polar_ring_gram

    real(dp) module function polar_ring_selfpower( Us, ldu, row0, nrow )
        integer, intent(in) :: ldu, row0, nrow
        real,    intent(in) :: Us(ldu,0:*)
    end function polar_ring_selfpower

    module subroutine polar_hybrid_exact_accum( rec0, recs, ncomp, o, fpl, hex, kex, npos, &
            &Gd, bd, cd, e_mm, myv )
        type(reconstructor), intent(in)    :: rec0
        type(reconstructor), intent(in)    :: recs(ncomp)
        integer,             intent(in)    :: ncomp, npos
        class(ori),          intent(inout) :: o
        type(fplane_type),   intent(in)    :: fpl
        integer,             intent(in)    :: hex(npos), kex(npos)
        real(dp),            intent(inout) :: Gd(ncomp,ncomp), bd(ncomp), cd(ncomp)
        real(dp),            intent(inout) :: e_mm, myv
    end subroutine polar_hybrid_exact_accum

    module subroutine subtract_mean_banded( fpl, mean_fpl, a, rec_nyq )
        type(fplane_type), intent(inout) :: fpl
        type(fplane_type), intent(in)    :: mean_fpl
        real,              intent(in)    :: a
        integer,           intent(in)    :: rec_nyq
    end subroutine subtract_mean_banded

    module subroutine fit_iter_finish( params, build, fit, it_eff, nthr )
        class(parameters),     intent(inout) :: params
        type(builder),         intent(inout) :: build
        class(flex_probe_fit), intent(inout) :: fit
        integer,               intent(in)    :: it_eff, nthr
    end subroutine fit_iter_finish

    module subroutine fit_estep_begin_stage( fit, params, nthr )
        class(flex_probe_fit), intent(inout) :: fit
        class(parameters),     intent(inout) :: params
        integer,               intent(in)    :: nthr
    end subroutine fit_estep_begin_stage

    module subroutine fit_iter_begin( params, build, fit, mean_rec, it_eff, niters_eff, nthr )
        class(parameters),     intent(inout) :: params
        type(builder),         intent(inout) :: build
        class(flex_probe_fit), intent(inout) :: fit
        type(reconstructor),   intent(inout) :: mean_rec
        integer,               intent(in)    :: it_eff, niters_eff, nthr
    end subroutine fit_iter_begin

    module subroutine probe_fit_merge_stash( fit )
        class(flex_probe_fit), intent(inout) :: fit
    end subroutine probe_fit_merge_stash

    module subroutine paired_reduce_parts( params, fits, rounds )
        class(flex_pca_rounds), intent(inout) :: rounds
        class(parameters),      intent(in)    :: params
        type(flex_probe_fit),   intent(inout) :: fits(2)
    end subroutine paired_reduce_parts

    module subroutine write_probe_part( fname, part )
        class(string),         intent(in) :: fname
        type(flex_probe_part), intent(in) :: part
    end subroutine write_probe_part

    module subroutine reduce_probe_parts( params, rounds, part )
        class(parameters),     intent(in)    :: params
        class(flex_pca_rounds), intent(inout) :: rounds
        type(flex_probe_part), intent(inout) :: part
    end subroutine reduce_probe_parts

    module subroutine open_probe_part_write( fname, nfits, funit, tmp_fname )
        class(string), intent(in)  :: fname
        integer,       intent(in)  :: nfits
        integer,       intent(out) :: funit
        type(string),  intent(out) :: tmp_fname
    end subroutine open_probe_part_write

    module subroutine write_probe_part_fit( funit, part )
        integer,               intent(in) :: funit
        type(flex_probe_part), intent(in) :: part
    end subroutine write_probe_part_fit

    module subroutine close_probe_part_write( funit, tmp_fname, fname )
        integer,       intent(in)    :: funit
        type(string),  intent(inout) :: tmp_fname
        class(string), intent(in)    :: fname
    end subroutine close_probe_part_write

    module subroutine open_probe_part_read( fname, nfits, funit )
        class(string), intent(in)  :: fname
        integer,       intent(in)  :: nfits
        integer,       intent(out) :: funit
    end subroutine open_probe_part_read

    module subroutine fold_probe_part_fit( funit, part )
        integer,               intent(in)    :: funit
        type(flex_probe_part), intent(inout) :: part
    end subroutine fold_probe_part_fit

    module subroutine close_probe_part_read( funit, fname )
        integer,       intent(in) :: funit
        class(string), intent(in) :: fname
    end subroutine close_probe_part_read

    module subroutine probe_subspace_iteration( params, build, plane_store, pcg_env, model, sel, niters, it_glob, &
        &niters_glob, fprefix, meta_fname, rounds )
        class(flex_pca_rounds),     intent(inout)         :: rounds
        class(parameters),          intent(inout)         :: params
        type(builder),              intent(inout)         :: build
        class(flex_plane_store),    intent(inout)         :: plane_store
        class(flex_pcg_environment), intent(in)           :: pcg_env
        type(flex_fit_model),       intent(inout), target :: model
        type(flex_selection),       intent(in)            :: sel
        integer,                    intent(in)            :: niters
        integer,          optional, intent(in)            :: it_glob, niters_glob
        character(len=*), optional, intent(in)            :: fprefix, meta_fname
    end subroutine probe_subspace_iteration

    module subroutine fit_engine_iterate( params, build, plane_store, fits, means, nfits, niters, &
        &it_glob, niters_glob, l_merge_stash, rounds )
        class(flex_pca_rounds),  intent(inout) :: rounds
        class(parameters),       intent(inout) :: params
        type(builder),           intent(inout) :: build
        class(flex_plane_store), intent(inout) :: plane_store
        integer,                 intent(in)    :: nfits
        type(flex_probe_fit),    intent(inout) :: fits(nfits)
        type(flex_mean_ref),     intent(in)    :: means(nfits)
        integer,                 intent(in)    :: niters, it_glob, niters_glob
        logical,                 intent(in)    :: l_merge_stash
    end subroutine fit_engine_iterate

    module subroutine xfsc_setup( ctx, params, kfr_ann, l_paired, l_master )
        type(xfsc_ctx_t),  intent(inout) :: ctx
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: kfr_ann(2)
        logical,           intent(in)    :: l_paired, l_master
    end subroutine xfsc_setup

    module subroutine xfsc_prep_iter( ctx, params, fit, it_eff, tag )
        type(xfsc_ctx_t),     intent(inout) :: ctx
        class(parameters),    intent(in)    :: params
        type(flex_probe_fit), intent(inout) :: fit
        integer,              intent(in)    :: it_eff
        character(len=*),     intent(in)    :: tag
    end subroutine xfsc_prep_iter

    module subroutine xfsc_teardown( ctx )
        type(xfsc_ctx_t), intent(inout) :: ctx
    end subroutine xfsc_teardown

    module subroutine xfsc_paired_record( ctx, params, fits, it_eff )
        type(xfsc_ctx_t),     intent(inout) :: ctx
        class(parameters),    intent(in)    :: params
        type(flex_probe_fit), intent(inout) :: fits(2)
        integer,              intent(in)    :: it_eff
    end subroutine xfsc_paired_record

end interface

end module simple_flex_probe_fit
