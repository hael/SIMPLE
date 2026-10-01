!@descr: flex_pca: the M-step backend -- storage and lifecycle of one fit's coupled M-step system, gridding or PCG
module simple_flex_pca_mstep
use simple_core_module_api
use simple_parameters,    only: parameters
use simple_builder,       only: builder
use simple_image,         only: image
use simple_ori,           only: ori
use simple_reconstructor, only: reconstructor
use simple_flex_pca_pcg,  only: flex_pcg_t, flex_pcg_outcome_t, flex_pcg_install_window
use simple_flex_pca_run_types, only: flex_run_settings
use simple_flex_reconstructor_latent_ops, only: insert_planes_oversamp_coupled_batch_scaled, solve_coupled_basis_exp,&
    &add_invtausq2rho_coupled
implicit none

public :: flex_fit_mstep, init_basis_reconstructor
private
#include "simple_local_flags.inc"

!> the M-step system: even/odd numerators and coupled densities, the PCG operator and its packed kernels, and the paired-merge snapshot
type :: flex_fit_mstep
    type(reconstructor), allocatable :: Yeven(:), Yodd(:)
    real,     allocatable :: rho_e(:,:,:,:), rho_o(:,:,:,:)
    integer :: npairs = 0
    integer :: es(3) = 0
    logical :: l_mg_stash = .false.  !< stash present (driver-gated)
    integer :: mg_ncomp = 0, mg_npairs = 0  !< entry rank of the stashed iteration
    complex,  allocatable :: mg_ye(:,:,:,:), mg_yo(:,:,:,:)  !< (es) x ncomp numerators
    real,     allocatable :: mg_rhe(:,:,:,:), mg_rho(:,:,:,:)  !< (npairs, es) packed densities
    real,     allocatable :: mg_kpe(:,:), mg_kpo(:,:)  !< (npairs, npk) packed PCG pair kernels on the band list (rec_backend=pcg)
    complex,  allocatable :: mg_rpe(:,:), mg_rpo(:,:)  !< (ncomp, npk) packed PCG right-hand sides on the band list (rec_backend=pcg)
    logical :: l_pcg = .false.  !< this fit solves its M-step by PCG
    type(flex_pcg_t) :: pcg  !< lattice, envelopes, support, finalized kernels
    real,     allocatable :: kacc_e(:,:), kacc_o(:,:)  !< expanded-lattice accumulators on the band list (inserting processes)
    real,     allocatable :: kpk_e(:,:), kpk_o(:,:)  !< packed pair kernels on the physical band list: transport, reduction and merge form
    complex,  allocatable :: racc_e(:,:), racc_o(:,:)  !< right-hand-side accumulators on the band list
    complex,  allocatable :: rpk_e(:,:), rpk_o(:,:)  !< packed right-hand sides on the physical band list
    type(image), allocatable :: mg_prev(:)  !< entry-frame orthonormal basis images
  contains
    ! ---- the M-step backend contract: storage and lifecycle of the coupled system, gridding
    ! (the existing coupled insertion and solve kernels) or PCG (flex_pcg_t) ----
    procedure :: begin_iteration         => mstep_begin_iteration          !< the system at this rank (and the PCG lattice)
    procedure :: accumulate_batch        => mstep_accumulate_batch         !< Y_q, rho (and the PCG kernels/rhs) += one batch
    procedure :: reduce_local            => mstep_reduce_local             !< fold this process's PCG accumulators into the packed sums
    procedure :: apply_ridge             => mstep_apply_ridge              !< the cross-fit SSNR ridge on rho (and the PCG operator)
    procedure :: solve_halves            => mstep_solve_halves             !< the coupled per-halfset solve
    procedure :: snapshot_for_pair_merge => mstep_snapshot_for_pair_merge  !< raw statistics + entry frame for the paired merge
    procedure :: kill_iteration          => mstep_kill_iteration           !< free the numerators and densities after the solve
    procedure :: kill                    => mstep_kill
end type flex_fit_mstep


contains

    !> A band-limited basis reconstructor on the cropped lattice, expanded accumulators zeroed.
    subroutine init_basis_reconstructor( params, build, rec )
        class(parameters),   intent(inout) :: params
        type(builder),       intent(inout) :: build
        type(reconstructor), intent(inout) :: rec
        call rec%new([params%box_crop,params%box_crop,params%box_crop],params%smpd_crop, wthreads=.true.)
        call rec%alloc_rho(params,build%spproj,expand=.true.)
        call rec%reset
        call rec%reset_exp
    end subroutine init_basis_reconstructor

    !> The system at this rank: even/odd Y_q accumulators (half-set FSC regularization), the
    !! COUPLED latent normal matrix and, under rec_backend=pcg, the solver's lattice and support
    !! with the packed kernel sums zeroed (the full-range accumulators are allocated by the first
    !! batch insert of an inserting process).
    subroutine mstep_begin_iteration( self, params, build, cfg, ncomp )
        class(flex_fit_mstep),   intent(inout) :: self
        class(parameters),       intent(inout) :: params
        type(builder),           intent(inout) :: build
        type(flex_run_settings), intent(in)    :: cfg
        integer,                 intent(in)    :: ncomp
        integer :: q
        ! even/odd Y_q accumulators (half-set FSC regularization) + the COUPLED latent normal matrix.
        ! rho carries one entry per (q,r) pair, not one shared density: the M-step below solves the
        ! components together at every grid point.
        allocate(self%Yeven(ncomp), self%Yodd(ncomp))
        do q = 1, ncomp
            call init_basis_reconstructor(params, build, self%Yeven(q)); call self%Yeven(q)%reset; call self%Yeven(q)%reset_exp
            call init_basis_reconstructor(params, build, self%Yodd(q));  call self%Yodd(q)%reset;  call self%Yodd(q)%reset_exp
        end do
        self%es     = shape(self%Yeven(1)%cmat_exp)
        ! The per-voxel coupled normal matrix is always the FULL ncomp x ncomp SPD system.
        ! A diagonal approximation used to be selectable here; it dropped the cross-component
        ! terms, and on EMPIAR-10028 that cost the 40S rotation outright -- the mode only
        ! appeared once the full solve was restored. It is not a speed/accuracy trade worth
        ! offering, so the option is gone rather than merely defaulted off.
        self%npairs = (ncomp*(ncomp+1))/2
        allocate(self%rho_e(self%npairs,self%es(1),self%es(2),self%es(3)), self%rho_o(self%npairs,self%es(1),self%es(2),self%es(3)), source=0.)
        write(logfhandle,'(A,I0,A,F8.2,A)') '>>> FLEX_PCA PROBE coupled normal matrix rows=',self%npairs, &
            &' (full)  rho even+odd ', 8.d0*real(self%npairs,dp)*real(self%es(1),dp)*real(self%es(2),dp)*real(self%es(3),dp)/1.d9,' GB'
        call flush(logfhandle)
        ! rec_backend=pcg: the solver's lattice and support at this rank, the packed kernel sums zeroed
        ! (the full-range accumulators are allocated by the first batch insert of an inserting process)
        self%l_pcg = trim(params%rec_backend) == 'pcg'
        if( allocated(self%kacc_e) ) deallocate(self%kacc_e)
        if( allocated(self%kacc_o) ) deallocate(self%kacc_o)
        if( allocated(self%kpk_e) ) deallocate(self%kpk_e)
        if( allocated(self%kpk_o) ) deallocate(self%kpk_o)
        if( allocated(self%racc_e) ) deallocate(self%racc_e)
        if( allocated(self%racc_o) ) deallocate(self%racc_o)
        if( allocated(self%rpk_e) ) deallocate(self%rpk_e)
        if( allocated(self%rpk_o) ) deallocate(self%rpk_o)
        if( self%l_pcg )then
            call self%pcg%new(params%box_crop, params%smpd_crop, ncomp)
            call self%pcg%set_verbose(cfg%pcg_verbose)
            if( cfg%l_pcg_lambda_set ) call self%pcg%set_lambda_relative(cfg%pcg_lambda_rel)
            call flex_pcg_install_window(self%pcg, params)
            call self%pcg%alloc_packed(self%kpk_e)
            call self%pcg%alloc_packed(self%kpk_o)
            call self%pcg%alloc_rhs_packed(self%rpk_e)
            call self%pcg%alloc_rhs_packed(self%rpk_o)
            write(logfhandle,'(A,F8.2,A,F8.2,A,F8.2,A,F8.2,A)') '>>> FLEX_PCA PROBE PCG pair kernels: packed even+odd ', &
                &2.d0*self%pcg%bytes_packed()/1.d9, ' GB, accumulators even+odd ', 2.d0*self%pcg%bytes_accum()/1.d9, &
                &' GB; right-hand sides: packed ', 2.d0*self%pcg%bytes_rhs_packed()/1.d9, ' GB, accumulators ', &
                &2.d0*self%pcg%bytes_rhs_accum()/1.d9, ' GB'
            call flush(logfhandle)
        endif
    end subroutine mstep_begin_iteration

    !> One batch into the halfset numerators and the coupled normal matrix (batched KB), and
    !! under rec_backend=pcg the same batch's pair-weighted Gram kernels and right-hand sides.
    subroutine mstep_accumulate_batch( self, build, orientations, fpls, zbatch, dens, valid_e, valid_o, batchsz )
        class(flex_fit_mstep), intent(inout) :: self
        type(builder),         intent(inout) :: build
        type(ori),             intent(inout) :: orientations(:)
        type(fplane_type),     intent(inout) :: fpls(:)
        real(dp),              intent(in)    :: zbatch(:,:), dens(:,:,:)
        logical,               intent(in)    :: valid_e(:), valid_o(:)
        integer,               intent(in)    :: batchsz
        call insert_planes_oversamp_coupled_batch_scaled(self%Yeven, self%rho_e, build%pgrpsyms, &
            &orientations(:batchsz), fpls(:batchsz), zbatch(:,:batchsz), dens(:,:,:batchsz), &
            &valid_e(:batchsz), batchsz)
        call insert_planes_oversamp_coupled_batch_scaled(self%Yodd, self%rho_o, build%pgrpsyms, &
            &orientations(:batchsz), fpls(:batchsz), zbatch(:,:batchsz), dens(:,:,:batchsz), &
            &valid_o(:batchsz), batchsz)
        ! rec_backend=pcg: the pair-weighted Gram kernels of the same batch at doubled coordinates
        ! (accumulators live only on inserting processes; the master reduces the packed sums)
        if( self%l_pcg )then
            if( .not. allocated(self%kacc_e) ) call self%pcg%alloc_accum(self%kacc_e)
            if( .not. allocated(self%kacc_o) ) call self%pcg%alloc_accum(self%kacc_o)
            if( .not. allocated(self%racc_e) ) call self%pcg%alloc_rhs_accum(self%racc_e)
            if( .not. allocated(self%racc_o) ) call self%pcg%alloc_rhs_accum(self%racc_o)
            call self%pcg%accumulate(self%kacc_e, build%pgrpsyms, orientations(:batchsz), fpls(:batchsz), &
                &dens(:,:,:batchsz), valid_e(:batchsz), batchsz)
            call self%pcg%accumulate(self%kacc_o, build%pgrpsyms, orientations(:batchsz), fpls(:batchsz), &
                &dens(:,:,:batchsz), valid_o(:batchsz), batchsz)
            ! the same batch's right-hand sides at doubled coordinates (the exact adjoint for the solve)
            call self%pcg%accumulate_rhs(self%racc_e, build%pgrpsyms, orientations(:batchsz), fpls(:batchsz), &
                &zbatch(:,:batchsz), valid_e(:batchsz), batchsz)
            call self%pcg%accumulate_rhs(self%racc_o, build%pgrpsyms, orientations(:batchsz), fpls(:batchsz), &
                &zbatch(:,:batchsz), valid_o(:batchsz), batchsz)
        endif
    end subroutine mstep_accumulate_batch

    !> Fold this process's PCG kernel and right-hand-side accumulators into the packed sums.
    subroutine mstep_reduce_local( self )
        class(flex_fit_mstep), intent(inout) :: self
        ! rec_backend=pcg: fold this process's kernel accumulators into the packed sums (adds; the
        ! distributed master has no accumulators and receives the parts' packed sums instead)
        if( self%l_pcg )then
            if( allocated(self%kacc_e) ) call self%pcg%fold_accum(self%kacc_e, self%kpk_e)
            if( allocated(self%kacc_o) ) call self%pcg%fold_accum(self%kacc_o, self%kpk_o)
            if( allocated(self%racc_e) ) call self%pcg%fold_rhs(self%racc_e, self%rpk_e)
            if( allocated(self%racc_o) ) call self%pcg%fold_rhs(self%racc_o, self%rpk_o)
        endif
    end subroutine mstep_reduce_local

    !> The cross-fit SSNR shrinkage ridge (the flex analog of add_invtausq2rho) on the diagonal
    !! rows of rho_e and rho_o; the PCG operator carries the same ridge as the preconditioner.
    subroutine mstep_apply_ridge( self, ncomp, invtau2 )
        class(flex_fit_mstep), intent(inout) :: self
        integer,               intent(in)    :: ncomp
        real,                  intent(in)    :: invtau2(:,:)
        call add_invtausq2rho_coupled(self%Yeven, self%rho_e, ncomp, invtau2)
        call add_invtausq2rho_coupled(self%Yodd,  self%rho_o, ncomp, invtau2)
        if( self%l_pcg ) call self%pcg%set_ridge(invtau2)
    end subroutine mstep_apply_ridge

    !> Coupled M-step solve: at every grid point the components share one k x k normal matrix
    !! sum_i |CTF|^2 E[z_i z_i'], so the basis volumes are solved together per halfset. A unit
    !! density is then handed to compress_exp by the caller (the divide has already happened).
    subroutine mstep_solve_halves( self, params, ncomp, it_eff )
        class(flex_fit_mstep), intent(inout) :: self
        class(parameters),     intent(in)    :: params
        integer,               intent(in)    :: ncomp, it_eff
        type(flex_pcg_outcome_t) :: pcg_out
        if( self%l_pcg )then
            ! rec_backend=pcg: the per-voxel solve is the preconditioner of the coupled normal equations
            ! on the pair Gram kernels; each half is solved from it on the spherical support
            call self%pcg%finalize(self%kpk_e)
            call self%pcg%solve(self%Yeven, self%rho_e, self%rpk_e, params%maxits_pcg, params%rtol, pcg_out, &
                &'FLEX_PCA PCG MSTEP even')
            call log_pcg_outcome('even', pcg_out)
            call self%pcg%finalize(self%kpk_o)
            call self%pcg%solve(self%Yodd,  self%rho_o, self%rpk_o, params%maxits_pcg, params%rtol, pcg_out, &
                &'FLEX_PCA PCG MSTEP odd')
            call log_pcg_outcome('odd', pcg_out)
            call self%pcg%clear_ridge
        else
            call solve_coupled_basis_exp(self%Yeven, self%rho_e, ncomp)
            call solve_coupled_basis_exp(self%Yodd,  self%rho_o, ncomp)
        endif

      contains

        subroutine log_pcg_outcome( half, res )
            character(len=*),         intent(in) :: half
            type(flex_pcg_outcome_t), intent(in) :: res
            write(logfhandle,'(A,I0,A,A,A,I0,A,ES10.3,A,ES10.3,A,ES10.3,A,A,A,ES10.3,A,ES10.3,A,F7.4,A,ES10.3,A,F8.1)') &
                &'>>> FLEX_PCA PCG MSTEP it=', it_eff, ' ', trim(half), '  iters=', res%iteration_count, &
                &'  init=', res%initial_rel_residual, '  resid=', res%final_rel_residual, '  update=', &
                &res%final_rel_update, '  stop=', trim(res%stop_reason), '  |b|=', res%rhs_norm, '  |x0|=', &
                &res%start_norm, '  corr(b,Bx0)=', res%start_corr, '  scale=', res%start_scale, '  seconds=', res%seconds
            call flush(logfhandle)
        end subroutine log_pcg_outcome

    end subroutine mstep_solve_halves

    !> Snapshot the LAST-iteration raw M-step sufficient statistics + entry frame for the paired
    !! merge: called after the reductions and BEFORE apply_ridge / solve_halves (which ridge rho,
    !! mutate the numerators in place and free everything). Overwrites every iteration: any
    !! iteration can turn out to be the last. prev_real holds the PREVIOUS iteration's delivered
    !! basis == the frame the statistics' latents live in.
    subroutine mstep_snapshot_for_pair_merge( self, ncomp, prev_real )
        class(flex_fit_mstep),    intent(inout) :: self
        integer,                  intent(in)    :: ncomp
        type(image), allocatable, intent(in)    :: prev_real(:)
        integer :: q, es(3)
        if( .not. (allocated(self%Yeven) .and. allocated(self%rho_e)) ) &
            &THROW_HARD('snapshot_for_pair_merge: no live M-step accumulators to stash')
        es = self%es
        ! (re)size on rank or lattice change
        if( allocated(self%mg_ye) )then
            if( size(self%mg_ye,4) /= ncomp .or. size(self%mg_ye,1) /= es(1) )then
                deallocate(self%mg_ye, self%mg_yo, self%mg_rhe, self%mg_rho)
            endif
        endif
        if( .not. allocated(self%mg_ye) )then
            allocate(self%mg_ye (es(1),es(2),es(3),ncomp), self%mg_yo(es(1),es(2),es(3),ncomp))
            allocate(self%mg_rhe(self%npairs,es(1),es(2),es(3)), self%mg_rho(self%npairs,es(1),es(2),es(3)))
        endif
        do q = 1, ncomp
            self%mg_ye(:,:,:,q) = self%Yeven(q)%cmat_exp
            self%mg_yo(:,:,:,q) = self%Yodd(q)%cmat_exp
        end do
        self%mg_rhe = self%rho_e
        self%mg_rho = self%rho_o
        if( allocated(self%mg_kpe) ) deallocate(self%mg_kpe)
        if( allocated(self%mg_kpo) ) deallocate(self%mg_kpo)
        if( allocated(self%mg_rpe) ) deallocate(self%mg_rpe)
        if( allocated(self%mg_rpo) ) deallocate(self%mg_rpo)
        if( self%l_pcg )then
            allocate(self%mg_kpe, source=self%kpk_e)
            allocate(self%mg_kpo, source=self%kpk_o)
            allocate(self%mg_rpe, source=self%rpk_e)
            allocate(self%mg_rpo, source=self%rpk_o)
        endif
        self%mg_ncomp  = ncomp
        self%mg_npairs = self%npairs
        ! entry-frame basis images (absent only at iteration 1; the driver requires >= 2 iterations)
        if( allocated(self%mg_prev) )then
            do q = 1, size(self%mg_prev)
                call self%mg_prev(q)%kill
            end do
            deallocate(self%mg_prev)
        endif
        if( allocated(prev_real) )then
            if( size(prev_real) /= ncomp ) &
                &THROW_HARD('snapshot_for_pair_merge: entry frame rank /= accumulator rank')
            allocate(self%mg_prev(ncomp))
            do q = 1, ncomp
                call self%mg_prev(q)%copy(prev_real(q))
            end do
        endif
        self%l_mg_stash = .true.
    end subroutine mstep_snapshot_for_pair_merge

    !> Free the numerators and densities after the solve (the packed PCG kernels persist to the
    !! next begin_iteration, which reallocates them at the new rank).
    subroutine mstep_kill_iteration( self )
        class(flex_fit_mstep), intent(inout) :: self
        integer :: q
        do q = 1, size(self%Yeven)
            call self%Yeven(q)%dealloc_rho; call self%Yeven(q)%kill
            call self%Yodd(q)%dealloc_rho;  call self%Yodd(q)%kill
        end do
        deallocate(self%Yeven, self%Yodd, self%rho_e, self%rho_o)
    end subroutine mstep_kill_iteration

    subroutine mstep_kill( self )
        class(flex_fit_mstep), intent(inout) :: self
        integer :: q
        ! ---- M-step accumulators ----
        if( allocated(self%Yeven) )then
            do q = 1, size(self%Yeven)
                call self%Yeven(q)%dealloc_rho; call self%Yeven(q)%kill
            end do
            deallocate(self%Yeven)
        endif
        if( allocated(self%Yodd) )then
            do q = 1, size(self%Yodd)
                call self%Yodd(q)%dealloc_rho; call self%Yodd(q)%kill
            end do
            deallocate(self%Yodd)
        endif
        if( allocated(self%rho_e)  ) deallocate(self%rho_e)
        if( allocated(self%rho_o)  ) deallocate(self%rho_o)
        ! ---- paired-merge stash ----
        self%l_mg_stash = .false.
        self%mg_ncomp   = 0
        self%mg_npairs  = 0
        if( allocated(self%mg_ye)  ) deallocate(self%mg_ye)
        if( allocated(self%mg_yo)  ) deallocate(self%mg_yo)
        if( allocated(self%mg_rhe) ) deallocate(self%mg_rhe)
        if( allocated(self%mg_rho) ) deallocate(self%mg_rho)
        if( allocated(self%mg_prev) )then
            do q = 1, size(self%mg_prev)
                call self%mg_prev(q)%kill
            end do
            deallocate(self%mg_prev)
        endif
        ! PCG M-step accumulators (kernel and rhs, per half, plus the pair-merge stash) and the operator
        if( allocated(self%kacc_e) ) deallocate(self%kacc_e)
        if( allocated(self%kacc_o) ) deallocate(self%kacc_o)
        if( allocated(self%kpk_e)  ) deallocate(self%kpk_e)
        if( allocated(self%kpk_o)  ) deallocate(self%kpk_o)
        if( allocated(self%mg_kpe) ) deallocate(self%mg_kpe)
        if( allocated(self%mg_kpo) ) deallocate(self%mg_kpo)
        if( allocated(self%racc_e) ) deallocate(self%racc_e)
        if( allocated(self%racc_o) ) deallocate(self%racc_o)
        if( allocated(self%rpk_e)  ) deallocate(self%rpk_e)
        if( allocated(self%rpk_o)  ) deallocate(self%rpk_o)
        if( allocated(self%mg_rpe) ) deallocate(self%mg_rpe)
        if( allocated(self%mg_rpo) ) deallocate(self%mg_rpo)
        call self%pcg%kill
        self%l_pcg = .false.; self%npairs = 0; self%es = 0
    end subroutine mstep_kill

end module simple_flex_pca_mstep
