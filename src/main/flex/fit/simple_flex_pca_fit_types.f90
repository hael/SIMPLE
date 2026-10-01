!@descr: flex_pca fit state: the composed per-fit value, its lifecycle, and the distributed payload value
!!
!! What a probe EM fit IS, as data: identity and resolved policy (spec), the model handles
!! (model), the cross-iteration state (history), the E-step bank (estep), the M-step system
!! (mstep), diagnostics (diag) and the per-iteration workspace (iter). Each piece frees itself;
!! the fit's kill runs the pieces in dependency order. The EM engine (simple_flex_pca_em) extends
!! this type with its iteration procedures; nothing here computes.
module simple_flex_pca_fit_types
use simple_core_module_api
use simple_image,              only: image
use simple_reconstructor,      only: reconstructor
use simple_flex_pca_polar,     only: polar_grid_t, polar_grid_kill
use simple_flex_pca_crossfsc,  only: crossfsc_file
use simple_flex_pca_mstep,     only: flex_fit_mstep
use simple_flex_pca_run_types, only: flex_run_settings
use simple_flex_pca_records,   only: flex_fit_model, flex_selection
implicit none
private
#include "simple_local_flags.inc"

public :: flex_fit, flex_fit_spec, flex_fit_model, flex_fit_history, flex_fit_estep, flex_fit_mstep, flex_fit_diag, flex_fit_iter
public :: flex_probe_part, probe_part_borrow, probe_part_restore
public :: xfsc_ctx_t
public :: cleanup_plane

!> Per-fit EM state, owned by the driver so one loop can advance two resident fits.
!! Lifecycle: mix_* and their work arrays (rhs0th, mkth, lwth, rkth, mxa_*) live across
!! iterations; only kill_probe_fit or the rank-change resize in fit_iter_begin may free them.
!> identity, selection, artifact namespaces and the resolved per-fit policy (never mutated by an iteration)
type :: flex_fit_spec
    integer :: id = 0  !< 1 = fit A, 2 = fit B, 0 = single-fit (legacy)
    type(flex_run_settings) :: cfg  !< the run's resolved switches, stamped at construction
    type(string) :: fprefix  !< eigenvolume namespace (default 'flex_pca_pc')
    type(string) :: meta_fname  !< probe-state file (default COV_PROBE_META)
    type(flex_selection) :: sel  !< this fit's half of the master's selection
    integer, allocatable :: ppinds(:)  !< probe-stage subsample; owned per fit (the banked-
    integer :: npp = 0
    real(dp) :: conv_thresh = 0.d0
    logical :: l_mix_req = .false.
    integer :: kmix = 0, n_mix_warm = 3
    integer :: khi_full = 0, kfr_ann(2) = 0
    real :: dstep_ann = 0., lp_it = 0.
    integer :: n_probe_cm = 0, nml_plain = 0
    logical :: l_probe_mls = .false.
    logical :: l_deflate_mean = .false.
    integer :: vdfl = 0
  contains
    procedure :: kill => flex_fit_spec_kill
end type flex_fit_spec

!> cross-iteration state: convergence, the previous basis, and the MCFA mixture with its work arrays and distributed-reduce buffers (freed only here or by the rank-change resize)
type :: flex_fit_history
    real(dp) :: nll_prev = 0.d0
    logical :: l_converged = .false.
    type(image), allocatable :: prev_real(:)  !< previous iteration's orthonormal basis
    logical :: l_mix_active = .false.
    logical :: l_mix_used = .false.
    real(dp) :: ldOm_mix = 0.d0, ldOm_used = 0.d0
    real(dp), allocatable :: mix_xi(:,:), mix_Om(:,:), mix_Ominv(:,:), mix_pi(:)
    real(dp), allocatable :: mix_Omxi(:,:), mix_xiOx(:), mix_lpi(:)
    real(dp), allocatable :: rhs0th(:,:), mkth(:,:,:), lwth(:,:), rkth(:,:)
    real(dp), allocatable :: mxa_sr(:,:), mxa_sm(:,:,:), mxa_smm(:,:,:,:), mxa_sainv(:,:,:)
    real(dp), allocatable :: dm_sr(:), dm_sm(:,:), dm_smm(:,:,:), dm_sai(:,:), dm_z(:,:)
    integer :: dm_nz = 0
  contains
    procedure :: kill => flex_fit_history_kill
end type flex_fit_history

!> the polar E-step bank: grid geometry, pose-fixed direction assignment (per stage) and the per-iteration ring tables and thread scratch
type :: flex_fit_estep
    logical :: l_pol_es = .false., l_pol_grid = .false.
    logical :: l_pol_bank_it = .false., l_pol_hyb = .false., l_rhyb_off = .false.
    integer :: rhyb_req = 0, osamp_pol = 1
    integer :: ndir_es = 0, nsamp_es = 0, nsamp2_es = 0, nk_es = 0
    integer :: ph0_es = 0, pk0_es = 0, hlo_es = 0, hhi_es = 0, klo_es = 0
    integer :: nyqr_es = 0, nyqb_es = 0, rhyb_es = 0, npos_es = 0
    integer,  allocatable :: hex_es(:), kex_es(:)
    type(polar_grid_t) :: pg_es
    real,     allocatable :: rmatb_es(:,:,:), nrmb_es(:,:)
    real,     allocatable :: cae(:), sae(:)
    integer,  allocatable :: dir_es(:)
    logical,  allocatable :: dused_es(:)
    real,     allocatable :: UsallE(:,:,:)
    real(dp), allocatable :: CfE(:,:,:), Cm0E(:,:,:), c00E(:,:)
    complex,  allocatable :: UbankE(:,:,:)
    real,     allocatable :: CspE(:,:,:)
    real,     allocatable :: xws_es(:,:), wr_es(:,:), Reb_es(:,:)
    real(dp), allocatable :: wrd_es(:,:)
    real :: sec_bank = 0.
  contains
    procedure :: kill => flex_fit_estep_kill
end type flex_fit_estep

!> cross-fit-FSC payloads and ridge, and the per-thread timings
type :: flex_fit_diag
    real(dp), allocatable :: gam_dbg(:,:)
    real(dp), allocatable :: sec_proj_thr(:), sec_gram_thr(:)
    integer :: khi_fit = 0  !< per-fit band, written by the BAND/RANK diagnostic;
    ! ---- cross-fit-FSC (crossfsc) per-fit hooks; the driver (xfsc_ctx_t) owns every decision ----
    ! fit_iter_finish harvests the writer payloads when l_xf_harvest is set (H before any ridge touches rho)
    ! and applies xf_invtau2 (record t-1, xfsc_prep_iter) to the rho_e/rho_o diagonals before the solves, once.
    logical :: l_xf_harvest = .false.
    real,     allocatable :: xf_h_e(:,:), xf_h_o(:,:)  !< (filtsz,ncomp) per-shell H, pre-ridge
    integer,  allocatable :: xf_cnt(:)  !< (filtsz) shared per-shell voxel counts
    real,     allocatable :: xf_fscq(:,:)  !< (filtsz,ncomp) internal e/o FSC curves
    real(dp), allocatable :: xf_gam(:)  !< (ncomp) Gamma at the writer site
    real,     allocatable :: xf_invtau2(:,:)  !< (ncomp,filtsz) ridge from record t-1
  contains
    procedure :: kill => flex_fit_diag_kill
end type flex_fit_diag

!> per-iteration workspace: prior, thread scratch, batch rows, Gamma and likelihood accumulators, projection planes
type :: flex_fit_iter
    real(dp), allocatable :: prior(:)
    real(dp), allocatable :: Gth(:,:,:), Ath(:,:,:), bth(:,:), cth(:,:), zth(:,:)
    real(dp), allocatable :: Ainvth(:,:,:), Acpth(:,:,:)
    real(dp), allocatable :: hth(:,:)
    real(dp), allocatable :: nll_thr(:)
    real(dp), allocatable :: zbatch(:,:), dens(:,:,:)
    logical,  allocatable :: valid(:), valid_e(:), valid_o(:)
    real(dp), allocatable :: gam_thr(:,:), gam_acc(:), gam_sum(:)
    integer,  allocatable :: nval_thr(:)
    real(dp) :: nll_tot = 0.d0
    integer :: nval = 0
    type(fplane_type), allocatable :: basis_fpls(:,:), mean_fpl(:)
  contains
    procedure :: kill => flex_fit_iter_kill
end type flex_fit_iter

!> One probe iteration's distributed payload: what a worker ships and a master folds. NOT fit
!! state: the numerators are the part's own copies of the even/odd Y_q lattices; the coupled
!! densities, Gamma, the PCG kernels/right-hand sides and the mixture buffers are BORROWED
!! from the fit around the codec call (probe_part_borrow / probe_part_restore, by move_alloc,
!! never a copy) so the value describes exactly what the codec reads or writes.
type :: flex_probe_part
    integer  :: ncomp = 0, nval = 0
    real(dp) :: nll_sum = 0.d0
    complex,  allocatable :: cmat_e(:,:,:,:), cmat_o(:,:,:,:)   !< even/odd numerators (cmat_exp lattices)
    real,     allocatable :: rho_ex(:,:,:,:), rho_ox(:,:,:,:)   !< even/odd numerator densities
    real,     allocatable :: rho_e(:,:,:,:), rho_o(:,:,:,:)     !< coupled per-voxel densities (packed pairs)
    real(dp), allocatable :: gam_sum(:)                          !< Gamma numerator (summed, divided by the reduced nval)
    real,     allocatable :: kpk_e(:,:), kpk_o(:,:)              !< packed PCG pair kernels (rec_backend=pcg)
    complex,  allocatable :: rpk_e(:,:), rpk_o(:,:)              !< packed PCG right-hand sides
    real(dp), allocatable :: mix_sr(:), mix_sm(:,:), mix_smm(:,:,:), mix_sainv(:,:)  !< MCFA additive statistics
    real(dp), allocatable :: z_sub(:,:)                          !< writer: this part's latent subsample; reader: the pool
    integer  :: nz = 0                                           !< rows of z_sub in use (the pool's fill on the reader)
    logical  :: l_mix_borrowed = .false.                         !< restore returns the mixture buffers to the fit
  contains
    procedure :: kill => probe_part_kill
end type flex_probe_part

!> One probe EM fit (the paired engine's per-fit state, composed). The pieces own their
!! allocations and free them through their own kill; the fit's kill runs them in dependency order.
!! LIFECYCLE CONTRACT (the MCFA free-on-iteration crash): history%mix_* and its work arrays must
!! survive from one iteration's M-step to the next iteration's E-step -- they are freed ONLY by
!! kill or by the rank-change resize at iteration start, never by per-iteration cleanup.
type :: flex_fit
    type(flex_fit_spec)    :: spec
    type(flex_fit_model)   :: model
    type(flex_fit_history) :: history
    type(flex_fit_estep)   :: estep
    type(flex_fit_mstep)   :: mstep
    type(flex_fit_diag)    :: diag
    type(flex_fit_iter)    :: iter
  contains
    procedure, pass(fit) :: new  => flex_fit_new
    procedure, pass(fit) :: kill => flex_fit_kill
end type flex_fit


!> Cross-fit-FSC driver context, one per driver loop (fit_engine_iterate, entered through
!! probe_subspace_iteration or run_flex_pca_paired): the paired master's per-iteration record
!! writer and the SSNR ridge (arm 1, hard-wired in xfsc_setup) built from record t-1. fit_iter_finish
!! stashes the writer payloads in probe_fit_t; the driver owns every consumer decision.
type :: xfsc_ctx_t
    logical  :: l_writer   = .false.   !< this driver writes records (paired master)
    logical  :: l_any      = .false.   !< ridge or writer live -> artifact state resident
    logical  :: l_loaded   = .false.   !< artifact reloaded from disk (restart-complete series)
    logical  :: l_paired   = .false.   !< this driver is the paired engine
    integer  :: v_reg      = 0         !< par.4 arms: 0 internal / 1 cross / 2 blend
    integer  :: reg_active = 0         !< arm ACTIVE this iteration (0 when degraded)
    integer  :: klo        = 6         !< low-resolution exemption index (reslim_ind analog)
    integer  :: filtsz     = 0
    integer  :: pairing_id = 0         !< balanced mod-4 pairing id (paired engine; header field)
    type(crossfsc_file) :: xf
end type xfsc_ctx_t

contains

    !> Construct a fit shell: identity, file namespaces and the fit's particle selection.
    !! Model handles (basis/eigvals/sig2), the stage subsample and all iteration state are
    !! populated by the driver, mirroring the single-fit initialisation order.
    subroutine flex_fit_new( fit, cfg, id, fprefix, meta_fname, sel )
        class(flex_fit), intent(inout) :: fit
        type(flex_run_settings), intent(in) :: cfg
        integer,           intent(in)    :: id
        character(len=*),  intent(in)    :: fprefix, meta_fname
        type(flex_selection), intent(in) :: sel
        call fit%kill
        fit%spec%cfg        = cfg
        fit%spec%id         = id
        fit%spec%fprefix    = fprefix
        fit%spec%meta_fname = meta_fname
        fit%spec%sel        = sel
    end subroutine flex_fit_new


    !> Free everything a fit may hold, in dependency order: the per-iteration workspace first (in
    !! case a crash or early exit left it allocated), the diagnostics, the M-step system with its
    !! paired-merge stash, the polar bank, the cross-iteration history (the MCFA state is freed
    !! ONLY here or by the rank-change resize -- see the lifecycle contract), the model handles,
    !! and finally the identity and selection.
    subroutine flex_fit_kill( fit )
        class(flex_fit), intent(inout) :: fit
        call fit%iter%kill
        call fit%diag%kill
        call fit%mstep%kill
        call fit%estep%kill
        call fit%history%kill
        call fit%model%kill
        call fit%spec%kill
    end subroutine flex_fit_kill


    subroutine flex_fit_spec_kill( self )
        class(flex_fit_spec), intent(inout) :: self
        if( allocated(self%ppinds) ) deallocate(self%ppinds)
        call self%sel%kill
        call self%fprefix%kill
        call self%meta_fname%kill
        call self%cfg%kill
        self%id = 0; self%npp = 0
    end subroutine flex_fit_spec_kill


    subroutine flex_fit_history_kill( self )
        class(flex_fit_history), intent(inout) :: self
        integer :: q
        ! ---- MCFA mixture state (see the fit's lifecycle contract) ----
        if( allocated(self%mix_xi)   ) deallocate(self%mix_xi)
        if( allocated(self%mix_Om)   ) deallocate(self%mix_Om)
        if( allocated(self%mix_Ominv)) deallocate(self%mix_Ominv)
        if( allocated(self%mix_pi)   ) deallocate(self%mix_pi)
        if( allocated(self%mix_Omxi) ) deallocate(self%mix_Omxi)
        if( allocated(self%mix_xiOx) ) deallocate(self%mix_xiOx)
        if( allocated(self%mix_lpi)  ) deallocate(self%mix_lpi)
        if( allocated(self%rhs0th)   ) deallocate(self%rhs0th)
        if( allocated(self%mkth)     ) deallocate(self%mkth)
        if( allocated(self%lwth)     ) deallocate(self%lwth)
        if( allocated(self%rkth)     ) deallocate(self%rkth)
        if( allocated(self%mxa_sr)   ) deallocate(self%mxa_sr)
        if( allocated(self%mxa_sm)   ) deallocate(self%mxa_sm)
        if( allocated(self%mxa_smm)  ) deallocate(self%mxa_smm)
        if( allocated(self%mxa_sainv)) deallocate(self%mxa_sainv)
        if( allocated(self%dm_sr)    ) deallocate(self%dm_sr)
        if( allocated(self%dm_sm)    ) deallocate(self%dm_sm)
        if( allocated(self%dm_smm)   ) deallocate(self%dm_smm)
        if( allocated(self%dm_sai)   ) deallocate(self%dm_sai)
        if( allocated(self%dm_z)     ) deallocate(self%dm_z)
        self%dm_nz = 0
        self%l_mix_active = .false.; self%l_mix_used = .false.
        ! ---- previous-basis images ----
        if( allocated(self%prev_real) )then
            do q = 1, size(self%prev_real)
                call self%prev_real(q)%kill
            end do
            deallocate(self%prev_real)
        endif
        self%l_converged = .false.
    end subroutine flex_fit_history_kill


    subroutine flex_fit_estep_kill( self )
        class(flex_fit_estep), intent(inout) :: self
        ! ---- polar E-step bank (grid via polar_grid_kill) ----
        call polar_grid_kill(self%pg_es)
        if( allocated(self%UsallE) ) deallocate(self%UsallE)
        if( allocated(self%CfE)    ) deallocate(self%CfE)
        if( allocated(self%Cm0E)   ) deallocate(self%Cm0E)
        if( allocated(self%c00E)   ) deallocate(self%c00E)
        if( allocated(self%UbankE) ) deallocate(self%UbankE)
        if( allocated(self%CspE)   ) deallocate(self%CspE)
        if( allocated(self%xws_es) ) deallocate(self%xws_es)
        if( allocated(self%wr_es)  ) deallocate(self%wr_es)
        if( allocated(self%wrd_es) ) deallocate(self%wrd_es)
        if( allocated(self%Reb_es) ) deallocate(self%Reb_es)
        if( allocated(self%rmatb_es) ) deallocate(self%rmatb_es)
        if( allocated(self%nrmb_es)  ) deallocate(self%nrmb_es)
        if( allocated(self%dir_es) ) deallocate(self%dir_es)
        if( allocated(self%cae)    ) deallocate(self%cae)
        if( allocated(self%sae)    ) deallocate(self%sae)
        if( allocated(self%dused_es) ) deallocate(self%dused_es)
        if( allocated(self%hex_es) ) deallocate(self%hex_es)
        if( allocated(self%kex_es) ) deallocate(self%kex_es)
        self%l_pol_grid = .false.; self%l_pol_bank_it = .false.
    end subroutine flex_fit_estep_kill


    subroutine flex_fit_diag_kill( self )
        class(flex_fit_diag), intent(inout) :: self
        ! ---- cross-fit-FSC per-fit payloads (the artifact on disk is the persistent series) ----
        self%l_xf_harvest = .false.
        if( allocated(self%xf_h_e)    ) deallocate(self%xf_h_e)
        if( allocated(self%xf_h_o)    ) deallocate(self%xf_h_o)
        if( allocated(self%xf_cnt)    ) deallocate(self%xf_cnt)
        if( allocated(self%xf_fscq)   ) deallocate(self%xf_fscq)
        if( allocated(self%xf_gam)    ) deallocate(self%xf_gam)
        if( allocated(self%xf_invtau2)) deallocate(self%xf_invtau2)
        if( allocated(self%sec_proj_thr) ) deallocate(self%sec_proj_thr)
        if( allocated(self%sec_gram_thr) ) deallocate(self%sec_gram_thr)
        if( allocated(self%gam_dbg)   ) deallocate(self%gam_dbg)
        self%khi_fit = 0
    end subroutine flex_fit_diag_kill


    subroutine flex_fit_iter_kill( self )
        class(flex_fit_iter), intent(inout) :: self
        integer :: q, ithr
        if( allocated(self%prior)  ) deallocate(self%prior)
        if( allocated(self%Gth)    ) deallocate(self%Gth)
        if( allocated(self%Ath)    ) deallocate(self%Ath)
        if( allocated(self%bth)    ) deallocate(self%bth)
        if( allocated(self%cth)    ) deallocate(self%cth)
        if( allocated(self%zth)    ) deallocate(self%zth)
        if( allocated(self%Ainvth) ) deallocate(self%Ainvth)
        if( allocated(self%Acpth)  ) deallocate(self%Acpth)
        if( allocated(self%hth)    ) deallocate(self%hth)
        if( allocated(self%nll_thr)) deallocate(self%nll_thr)
        if( allocated(self%zbatch) ) deallocate(self%zbatch)
        if( allocated(self%dens)   ) deallocate(self%dens)
        if( allocated(self%valid)  ) deallocate(self%valid)
        if( allocated(self%valid_e)) deallocate(self%valid_e)
        if( allocated(self%valid_o)) deallocate(self%valid_o)
        if( allocated(self%gam_thr)) deallocate(self%gam_thr)
        if( allocated(self%gam_acc)) deallocate(self%gam_acc)
        if( allocated(self%gam_sum)) deallocate(self%gam_sum)
        if( allocated(self%nval_thr)) deallocate(self%nval_thr)
        if( allocated(self%mean_fpl) )then
            do ithr = 1, size(self%mean_fpl)
                call cleanup_plane(self%mean_fpl(ithr))
            end do
            deallocate(self%mean_fpl)
        endif
        if( allocated(self%basis_fpls) )then
            do ithr = 1, size(self%basis_fpls,2)
                do q = 1, size(self%basis_fpls,1)
                    call cleanup_plane(self%basis_fpls(q,ithr))
                end do
            end do
            deallocate(self%basis_fpls)
        endif
        self%nll_tot = 0.d0; self%nval = 0
    end subroutine flex_fit_iter_kill


    !> Hand the fit's reducible accumulators to the payload value without copying: the coupled
    !! densities, Gamma, the PCG kernels and right-hand sides, the scalars, and (on a reducing
    !! master) the mixture reduce buffers. probe_part_restore gives them back.
    subroutine probe_part_borrow( part, fit, l_mix_buffers )
        type(flex_probe_part), intent(inout) :: part
        class(flex_fit),       intent(inout) :: fit
        logical,               intent(in)    :: l_mix_buffers
        part%ncomp   = fit%model%ncomp
        part%nval    = fit%iter%nval
        part%nll_sum = fit%iter%nll_tot
        call move_alloc(fit%mstep%rho_e,   part%rho_e)
        call move_alloc(fit%mstep%rho_o,   part%rho_o)
        call move_alloc(fit%iter%gam_sum,  part%gam_sum)
        call move_alloc(fit%mstep%kpk_e,   part%kpk_e)
        call move_alloc(fit%mstep%kpk_o,   part%kpk_o)
        call move_alloc(fit%mstep%rpk_e,   part%rpk_e)
        call move_alloc(fit%mstep%rpk_o,   part%rpk_o)
        part%l_mix_borrowed = l_mix_buffers
        if( l_mix_buffers )then
            call move_alloc(fit%history%dm_sr,  part%mix_sr)
            call move_alloc(fit%history%dm_sm,  part%mix_sm)
            call move_alloc(fit%history%dm_smm, part%mix_smm)
            call move_alloc(fit%history%dm_sai, part%mix_sainv)
            call move_alloc(fit%history%dm_z,   part%z_sub)
            part%nz = fit%history%dm_nz
        endif
    end subroutine probe_part_borrow


    subroutine probe_part_restore( part, fit )
        type(flex_probe_part), intent(inout) :: part
        class(flex_fit),       intent(inout) :: fit
        fit%iter%nval    = part%nval
        fit%iter%nll_tot = part%nll_sum
        call move_alloc(part%rho_e,   fit%mstep%rho_e)
        call move_alloc(part%rho_o,   fit%mstep%rho_o)
        call move_alloc(part%gam_sum, fit%iter%gam_sum)
        call move_alloc(part%kpk_e,   fit%mstep%kpk_e)
        call move_alloc(part%kpk_o,   fit%mstep%kpk_o)
        call move_alloc(part%rpk_e,   fit%mstep%rpk_e)
        call move_alloc(part%rpk_o,   fit%mstep%rpk_o)
        if( part%l_mix_borrowed )then
            call move_alloc(part%mix_sr,    fit%history%dm_sr)
            call move_alloc(part%mix_sm,    fit%history%dm_sm)
            call move_alloc(part%mix_smm,   fit%history%dm_smm)
            call move_alloc(part%mix_sainv, fit%history%dm_sai)
            call move_alloc(part%z_sub,     fit%history%dm_z)
            fit%history%dm_nz = part%nz
            part%l_mix_borrowed = .false.
        endif
    end subroutine probe_part_restore


    subroutine probe_part_kill( self )
        class(flex_probe_part), intent(inout) :: self
        if( self%l_mix_borrowed ) THROW_HARD('probe_part_kill: borrowed mixture buffers were never restored')
        if( allocated(self%cmat_e) )    deallocate(self%cmat_e)
        if( allocated(self%cmat_o) )    deallocate(self%cmat_o)
        if( allocated(self%rho_ex) )    deallocate(self%rho_ex)
        if( allocated(self%rho_ox) )    deallocate(self%rho_ox)
        if( allocated(self%rho_e) )     deallocate(self%rho_e)
        if( allocated(self%rho_o) )     deallocate(self%rho_o)
        if( allocated(self%gam_sum) )   deallocate(self%gam_sum)
        if( allocated(self%kpk_e) )     deallocate(self%kpk_e)
        if( allocated(self%kpk_o) )     deallocate(self%kpk_o)
        if( allocated(self%rpk_e) )     deallocate(self%rpk_e)
        if( allocated(self%rpk_o) )     deallocate(self%rpk_o)
        if( allocated(self%mix_sr) )    deallocate(self%mix_sr)
        if( allocated(self%mix_sm) )    deallocate(self%mix_sm)
        if( allocated(self%mix_smm) )   deallocate(self%mix_smm)
        if( allocated(self%mix_sainv) ) deallocate(self%mix_sainv)
        if( allocated(self%z_sub) )     deallocate(self%z_sub)
        self%ncomp = 0; self%nval = 0; self%nll_sum = 0.d0; self%nz = 0
    end subroutine probe_part_kill


    subroutine cleanup_plane( fpl )
        type(fplane_type), intent(inout) :: fpl
        if( allocated(fpl%cmplx_plane)    ) deallocate(fpl%cmplx_plane)
        if( allocated(fpl%ctfsq_plane)    ) deallocate(fpl%ctfsq_plane)
        if( allocated(fpl%transfer_plane) ) deallocate(fpl%transfer_plane)
    end subroutine cleanup_plane



end module simple_flex_pca_fit_types
