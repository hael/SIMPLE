!@descr: flex_pca state maps on the reconstruct3D PCG backend (rec_states_backend=pcg): the kernel weight of a
!  particle enters through its noise model as sigma2/w, so the right-hand side and the density carry it
!  identically and the preconditioner, Gram kernel, ridge scale and raw artifacts follow without change
!  (doc/implementation_notes/flex_pca_envelope_support.md, section 3.3); one cold solve per (state, half)
module simple_flex_pca_states_pcg
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_builder,           only: builder
use simple_parameters,        only: parameters
use simple_image,             only: image
use simple_reconstructor_pcg, only: reconstructor_pcg, pcg_solver_outcome, PCG_OP_KERNEL, PCG_STOP_INDEFINITE
use simple_flex_pca_pcg,      only: flex_pcg_support_volume, flex_mskfile_set
use simple_matcher_ptcl_io,   only: prepimgbatch, discrete_read_imgbatch, prep_rec_observation
use simple_math_ft,           only: resample_sigma2
use simple_flex_pca_rounds,         only: flex_pca_rounds
use simple_flex_pca_run_types,      only: flex_run_settings
use simple_flex_pca_state_parts,    only: flex_pcg_state_raw_fname, flex_pcg_state_provenance
use simple_flex_pca_states_backend, only: flex_states_backend, flex_state_maps, flex_state_delivery_policy
!$ use omp_lib, only: omp_get_max_active_levels, omp_set_max_active_levels, omp_set_num_threads
implicit none

public :: flex_states_pcg
private
#include "simple_local_flags.inc"

!> kernel weights below this floor drop the particle from the state's selection: 1/w has to stay
!! finite in single precision and a vanishing weight contributes nothing to the solve
real,    parameter :: FLEX_PCG_WEIGHT_FLOOR = 1.0e-3
!> Tikhonov ridge relative to the weighted data scale: with weights in [0,1] the effective particle
!! count of a sparse state is a small fraction of N, so the refinement backend's absolute PCG_LAMBDA
!! would be a materially stronger prior on it than on a populated state
real,    parameter :: FLEX_PCG_LAMBDA_REL = 1.0e-3
!> iteration budget of the cold state solves: the reconstruct3D PCG default (maxits_pcg=2). It is NOT
!! params%maxits_pcg, which is the warm-started basis M-step's budget; a positive rtol still stops earlier
integer, parameter :: FLEX_PCG_STATE_MAXITS = 2
!> thread budget of the paired even/odd solve (PCG_MASTER_NTHR_CAP of the reconstruct3D PCG master)
integer, parameter :: FLEX_PCG_PAIR_NTHR_CAP = 32
character(len=*), parameter :: PCG_STATE_TABLE = 'flex_pca_state_pcg.txt'

!> outcome of one (state, half) solve, recorded inside the concurrent even/odd sections and reported after
type :: half_solve_rec
    type(pcg_solver_outcome) :: outcome
    integer :: niters  = 0
    real    :: seconds = 0.0
    logical :: l_solved = .false.   !< false: no particles above the floor, zero map delivered
end type half_solve_rec


!> The same weighted least-squares problems as the gridding backend on reconstructor_pcg with the
!! support inside the solve. Workers accumulate the raw (B,D) of their particle range and publish
!! one artifact per (state, half, part); the master (or the shared-memory process, which
!! accumulates in place) reduces in part order, finalizes the kernel and solves each half cold on
!! the spherical support -- one operator pair resident at a time, so the raw parts of a state are
!! folded when its maps are finalized. Delivered maps are window*u; nothing masks them afterwards.
type, extends(flex_states_backend) :: flex_states_pcg
    type(reconstructor_pcg) :: pcgops(0:1)
    type(half_solve_rec)    :: sol(0:1)
    character(len=256) :: provenance = ''
    real    :: msk_rec = 0.
    integer :: neo = 0, nthr_half = 1, nthr_pair = 1
    logical :: l_pair = .false.
  contains
    procedure :: begin                          => pcg_begin
    procedure :: accumulate_local_or_write_part => pcg_accumulate
    procedure :: fold_parts                     => pcg_fold_parts
    procedure :: finalize_maps                  => pcg_finalize_maps
    procedure :: delivery_policy                => pcg_delivery_policy
    procedure :: kill                           => pcg_kill
end type flex_states_pcg

contains

    subroutine pcg_begin( self, params, build, rounds, pinds, state_weights, nstates, l_fuse, l_floor_rho, box_rec, smpd_rec )
        class(flex_states_pcg),  intent(inout) :: self
        class(parameters),       intent(inout) :: params
        class(builder),          intent(inout) :: build
        class(flex_pca_rounds),  intent(inout) :: rounds
        integer,                 intent(in)    :: pinds(:), nstates, box_rec
        real,                    intent(in)    :: state_weights(:,:), smpd_rec
        logical,                 intent(in)    :: l_fuse, l_floor_rho
        call self%kill
        call self%set_selection(pinds, state_weights, nstates, l_fuse, l_floor_rho, box_rec, smpd_rec)
        ! the same spherical support as the gridding delivery mask, capped at the box edge
        self%msk_rec    = min(0.5*params%mskdiam/smpd_rec, real(box_rec/2) - COSMSKHALFWIDTH - 1.)
        self%provenance = flex_pcg_state_provenance(params, box_rec, smpd_rec)
        self%neo = 0
        if( l_fuse ) self%neo = 1
        ! even and odd of a state solve concurrently on half the threads each (the reconstruct3D PCG
        ! master's pairing); the accumulation/reduction before them stays serial (shared image batch, disk)
        ! total capped at the reconstruct3D PCG master's PCG_MASTER_NTHR_CAP (32): the FFT scaling of a
        ! single solve flattens well before the boosted master budget, and the cap bounds the thread stacks
        self%nthr_pair = min(nthr_glob, FLEX_PCG_PAIR_NTHR_CAP)
        self%l_pair    = l_fuse .and. self%nthr_pair >= 2
        self%nthr_half = self%nthr_pair
        if( self%l_pair ) self%nthr_half = max(1, self%nthr_pair/2)
    end subroutine pcg_begin

    !> Worker: the raw (B,D) per state and half of this part, one artifact each. Shared memory:
    !! only the image batch is prepared here; each state accumulates in place when its maps are
    !! finalized (one operator pair resident).
    subroutine pcg_accumulate( self, params, build, rounds )
        class(flex_states_pcg),  intent(inout) :: self
        class(parameters),       intent(inout) :: params
        class(builder),          intent(inout) :: build
        class(flex_pca_rounds),  intent(inout) :: rounds
        type(reconstructor_pcg) :: pcgop
        type(string) :: fname
        integer :: state, eo, nsel
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        if( .not. rounds%is_worker() ) return
        do state = 1, self%nstates
            do eo = 0, self%neo
                call accumulate_state_half(self, params, build, rounds, state, eo, pcgop, nsel)
                fname = flex_pcg_state_raw_fname(params, params%part, state, eo)
                ! the part count comes from the command line: rounds%nparts() is the master's plan (1 on a worker)
                call pcgop%write_raw_accum(fname, state, eo, params%part, max(1,params%nparts), nsel, self%provenance)
                write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> FLEX_PCA PCG STATE RAW: part ', params%part, &
                    &' state ', state, ' half ', eo, ' particles ', nsel
                call pcgop%kill
                call fname%kill
            end do
        end do
        call flush(logfhandle)
    end subroutine pcg_accumulate

    !> Distributed master: nothing to do up front -- the raw parts of one state are reduced when
    !! its maps are finalized, so one operator pair is resident rather than nstates.
    subroutine pcg_fold_parts( self, params, build, rounds )
        class(flex_states_pcg),  intent(inout) :: self
        class(parameters),       intent(inout) :: params
        class(builder),          intent(inout) :: build
        class(flex_pca_rounds),  intent(inout) :: rounds
        write(logfhandle,'(A,I0,A)') '>>> FLEX_PCA PCG STATES: ', max(1,params%nparts), &
            &' raw parts per (state, half), reduced per state at its solve'
        call flush(logfhandle)
    end subroutine pcg_fold_parts

    !> Reduce (distributed) or accumulate (shared memory) both halves of one state, solve them
    !! concurrently, report, and hand the maps over: window*u per half, combined = (even+odd)/2.
    subroutine pcg_finalize_maps( self, params, build, rounds, state, maps )
        class(flex_states_pcg),  intent(inout) :: self
        class(parameters),       intent(inout) :: params
        class(builder),          intent(inout) :: build
        class(flex_pca_rounds),  intent(inout) :: rounds
        integer,                 intent(in)    :: state
        type(flex_state_maps),   intent(inout) :: maps
        type(string) :: fname
        integer :: eo, nsel, nsel_part, ipart, nsel_eo(0:1), prev_levels
        call maps%kill
        nsel_eo = 0
        do eo = 0, self%neo
            if( rounds%distributed() )then
                call new_state_operator(self, params, build, rounds, self%pcgops(eo))
                write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> FLEX_PCA PCG STATE ', state, ' half ', eo, &
                    &': reducing ', max(1,params%nparts), ' raw parts and finalizing the kernel'
                call flush(logfhandle)
                call self%pcgops(eo)%begin_reduction
                nsel = 0
                do ipart = 1, max(1,params%nparts)
                    fname = flex_pcg_state_raw_fname(params, ipart, state, eo)
                    if( .not. file_exists(fname) ) THROW_HARD('missing flex PCG raw part: '//fname%to_char())
                    call self%pcgops(eo)%add_raw_accum(fname, state, eo, ipart, max(1,params%nparts), self%provenance, nsel_part)
                    nsel = nsel + nsel_part
                    call del_file(fname)
                    call fname%kill
                end do
            else
                call accumulate_state_half(self, params, build, rounds, state, eo, self%pcgops(eo), nsel)
            endif
            nsel_eo(eo) = nsel
        end do
        if( self%l_pair )then
            prev_levels = 1
            !$ prev_levels = omp_get_max_active_levels()
            !$ call omp_set_max_active_levels(max(2, prev_levels))
            !$omp parallel sections num_threads(2) default(shared)
            !$omp section
            !$ call omp_set_num_threads(self%nthr_half)
            call solve_state_half(self, params, self%pcgops(0), nsel_eo(0), maps%even, self%sol(0))
            !$omp section
            !$ call omp_set_num_threads(self%nthr_half)
            call solve_state_half(self, params, self%pcgops(1), nsel_eo(1), maps%odd, self%sol(1))
            !$omp end parallel sections
            !$ call omp_set_max_active_levels(prev_levels)
            !$ call omp_set_num_threads(nthr_glob)
        else
            call solve_state_half(self, params, self%pcgops(0), nsel_eo(0), maps%even, self%sol(0))
            if( self%l_fuse ) call solve_state_half(self, params, self%pcgops(1), nsel_eo(1), maps%odd, self%sol(1))
        endif
        ! serial: the outcome checks, the log lines and the table, in even/odd order
        do eo = 0, self%neo
            call report_state_half(self, state, eo, nsel_eo(eo), self%sol(eo))
            call self%pcgops(eo)%kill
        end do
        if( self%l_fuse )then
            call maps%combined%copy(maps%even)
            call maps%combined%add(maps%odd)
            call maps%combined%mul(0.5)
            maps%l_halves = .true.
        else
            ! single-set state: one solve over every weighted particle
            call maps%combined%copy(maps%even)
            call maps%even%kill
            maps%l_halves = .false.
        endif
    end subroutine pcg_finalize_maps

    !> The PCG delivery: the same per-state eo-FSC filter as the gridding path (selected by the run
    !! settings), on the windowed solutions as they come out of the solve -- no background
    !! removal, no second mask, no project-FSC fallback (a single-set state is delivered unfiltered).
    function pcg_delivery_policy( self, cfg ) result( policy )
        class(flex_states_pcg),  intent(in) :: self
        type(flex_run_settings), intent(in) :: cfg
        type(flex_state_delivery_policy) :: policy
        policy%l_state_eofilt = cfg%l_state_eofilt
        policy%l_state_filt   = cfg%l_state_filt
        policy%l_mask         = .false.
        policy%l_project_fsc_fallback = .false.
        policy%tag = ' (PCG)'
    end function pcg_delivery_policy

    subroutine pcg_kill( self )
        class(flex_states_pcg), intent(inout) :: self
        integer :: eo
        do eo = 0, 1
            call self%pcgops(eo)%kill
        end do
        self%provenance = ''
        self%msk_rec = 0.; self%neo = 0; self%nthr_half = 1; self%nthr_pair = 1; self%l_pair = .false.
        call self%kill_selection
    end subroutine pcg_kill

    subroutine new_state_operator( self, params, build, rounds, op )
        class(flex_states_pcg),  intent(in)    :: self
        class(parameters),       intent(inout) :: params
        class(builder),          intent(inout) :: build
        class(flex_pca_rounds),  intent(inout) :: rounds
        type(reconstructor_pcg), intent(inout) :: op
        type(image) :: envimg
        if( rounds%is_worker() )then
            call op%new(self%box_rec, self%smpd_rec)
        else
            call op%new(self%box_rec, self%smpd_rec, fft_nthreads=self%nthr_half)
        endif
        call op%set_sym(build%pgrpsyms)
        call op%set_lambda_relative(FLEX_PCG_LAMBDA_REL)
        if( flex_mskfile_set(params) )then
            ! step 2: the envelope resampled to the state box is the solve support
            call flex_pcg_support_volume(params, self%box_rec, self%smpd_rec, envimg, 'state box')
            call op%set_mask_volume(envimg)
            call envimg%kill
        else
            call op%set_mask(self%msk_rec)
        endif
    end subroutine new_state_operator

    subroutine accumulate_state_half( self, params, build, rounds, state_here, eo_here, op, nsel_out )
        class(flex_states_pcg),  intent(inout) :: self
        class(parameters),       intent(inout) :: params
        class(builder),          intent(inout) :: build
        class(flex_pca_rounds),  intent(inout) :: rounds
        integer,                 intent(in)    :: state_here, eo_here
        type(reconstructor_pcg), intent(inout) :: op
        integer,                 intent(out)   :: nsel_out
        type(oris)      :: selection
        type(ori)       :: orientation
        type(ctfparams) :: ctfparms
        type(image)     :: obs
        integer, allocatable :: sel(:)
        real,    allocatable :: sig2(:,:), wsel(:)
        complex, allocatable :: y_batch(:,:,:)
        integer :: lims2(2,2), R, kfromto(2), batchlims(2), batchsz, i, ii, iptcl, ibatch, n
        real    :: shift(2), crop_factor, w
        call new_state_operator(self, params, build, rounds, op)
        allocate(sel(size(self%pinds)), wsel(size(self%pinds)))
        n = 0
        do i = 1, size(self%pinds)
            iptcl = self%pinds(i)
            w     = self%state_weights(i, state_here)
            if( w < FLEX_PCG_WEIGHT_FLOOR ) cycle
            if( build%spproj_field%get_state(iptcl) <= 0 ) cycle
            if( self%l_fuse )then
                if( build%spproj_field%get_eo(iptcl) /= eo_here ) cycle
            endif
            n       = n + 1
            sel(n)  = iptcl
            wsel(n) = w
        end do
        nsel_out = n
        if( n == 0 )then
            deallocate(sel, wsel)
            return
        endif
        lims2 = op%get_lims2()
        R     = lims2(1,2)
        allocate(sig2(0:R,n), source=1.0)
        if( params%cc_objfun == OBJFUN_EUCLID )then
            kfromto = build%esig%get_kfromto()
            do i = 1, n
                call resample_sigma2(kfromto(1), kfromto(2), &
                    &build%esig%sigma2_noise(kfromto(1):kfromto(2),sel(i)), R, 1.0, sig2(0:R,i))
            end do
        endif
        ! the kernel weight as an inverse noise scale: it multiplies the RHS term and |T|^2 alike
        do i = 1, n
            sig2(:,i) = sig2(:,i) / wsel(i)
        end do
        crop_factor = real(self%box_rec) / real(params%box)
        call selection%new(n, .true.)
        call orientation%new(.false.)
        do i = 1, n
            iptcl = sel(i)
            call build%spproj_field%get_ori(iptcl, orientation)
            ctfparms      = build%spproj%get_ctfparams(params%oritype, iptcl)
            ctfparms%smpd = self%smpd_rec
            shift         = build%spproj_field%get_2Dshift(iptcl) * crop_factor
            call orientation%set_ctfvars(ctfparms)
            call orientation%set_shift(shift)
            call selection%set_ori(i, orientation)
        end do
        call op%prep_particles(selection, use_ctf=.true., sig2=sig2)
        allocate(y_batch(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), MAXIMGBATCHSZ))
        call op%begin_accum
        call obs%new([self%box_rec,self%box_rec,1], self%smpd_rec)
        do ibatch = 1, n, MAXIMGBATCHSZ
            batchlims = [ibatch, min(n, ibatch+MAXIMGBATCHSZ-1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call discrete_read_imgbatch(params, build, n, sel, batchlims)
            do ii = 1, batchsz
                ! the backend-neutral observation (normalize, crop, taper) of the rec3D PCG strategy
                call prep_rec_observation(build%imgbatch(ii), build%lmsk, obs, .true.)
                call obs%fft
                y_batch(:,:,ii) = op%extract_native_plane(obs)
            end do
            call op%accumulate_batch(y_batch, batchsz, batchlims(1))
        end do
        call obs%kill
        call selection%kill
        call orientation%kill
        deallocate(y_batch, sig2, sel, wsel)
    end subroutine accumulate_state_half

    subroutine solve_state_half( self, params, op, nsel_here, vol, rec )
        class(flex_states_pcg),  intent(in)    :: self
        class(parameters),       intent(in)    :: params
        type(reconstructor_pcg), intent(inout) :: op
        integer,                 intent(in)    :: nsel_here
        type(image),             intent(inout) :: vol
        type(half_solve_rec),    intent(inout) :: rec
        real, allocatable :: x(:,:,:)
        integer(timer_int_kind) :: t_solve
        rec%l_solved = .false.
        rec%niters   = 0
        call vol%new([self%box_rec,self%box_rec,self%box_rec], self%smpd_rec)
        if( nsel_here == 0 ) return
        t_solve = tic()
        call op%end_accum(.true.)
        call op%set_op_mode(PCG_OP_KERNEL)
        allocate(x(self%box_rec,self%box_rec,self%box_rec), source=0.0)
        call op%solve_accum(x, maxits=FLEX_PCG_STATE_MAXITS, rtol=params%rtol, niters=rec%niters, &
            &outcome=rec%outcome)
        rec%l_solved = all(ieee_is_finite(x))
        if( rec%l_solved ) call vol%set_rmat(x, .false.)
        rec%seconds = real(toc(t_solve))
        deallocate(x)
    end subroutine solve_state_half

    subroutine report_state_half( self, state_here, eo_here, nsel_here, rec )
        class(flex_states_pcg), intent(in) :: self
        integer,                intent(in) :: state_here, eo_here, nsel_here
        type(half_solve_rec), intent(in) :: rec
        character(len=4) :: half
        half = 'all '
        if( self%l_fuse ) half = merge('odd ', 'even', eo_here == 1)
        if( nsel_here == 0 )then
            write(logfhandle,'(A,I0,A,A,A)') '>>> FLEX_PCA PCG STATE ', state_here, ' ', trim(half), &
                &': no particles above the weight floor; zero map delivered'
            return
        endif
        if( trim(rec%outcome%stop_reason) == PCG_STOP_INDEFINITE ) &
            &THROW_HARD('flex PCG state solve lost positive-definiteness (cold start)')
        if( .not. rec%l_solved ) THROW_HARD('flex PCG state solve returned non-finite values')
        write(logfhandle,'(A,I0,A,A,A,I0,A,I0,A,ES10.3,A,ES10.3,A,A,A,F8.1)') '>>> FLEX_PCA PCG STATE ', &
            &state_here, ' ', trim(half), '  nptcls=', nsel_here, '  iters=', rec%niters, &
            &'  resid=', rec%outcome%final_rel_residual, '  update=', rec%outcome%final_rel_update, &
            &'  stop=', trim(rec%outcome%stop_reason), '  seconds=', rec%seconds
        call flush(logfhandle)
        call append_state_table(self, state_here, half, nsel_here, rec%niters, rec%outcome)
    end subroutine report_state_half

    subroutine append_state_table( self, state_here, half, nsel_here, niters_here, res )
        class(flex_states_pcg),   intent(in) :: self
        integer,                  intent(in) :: state_here, nsel_here, niters_here
        character(len=*),         intent(in) :: half
        type(pcg_solver_outcome), intent(in) :: res
        integer :: io_stat, tunit
        logical :: l_exists
        inquire(file=PCG_STATE_TABLE, exist=l_exists)
        if( l_exists )then
            open(newunit=tunit, file=PCG_STATE_TABLE, status='old', position='append', action='write', iostat=io_stat)
        else
            open(newunit=tunit, file=PCG_STATE_TABLE, status='new', action='write', iostat=io_stat)
            if( io_stat == 0 ) write(tunit,'(A)') '# state half nptcls iters init_resid final_resid final_update stop box_rec msk_rec'
        endif
        if( io_stat /= 0 ) return
        write(tunit,'(I4,1X,A,1X,I8,1X,I4,1X,ES11.4,1X,ES11.4,1X,ES11.4,1X,A,1X,I5,1X,F8.2)') state_here, &
            &trim(half), nsel_here, niters_here, res%initial_rel_residual, res%final_rel_residual, &
            &res%final_rel_update, trim(res%stop_reason), self%box_rec, self%msk_rec
        close(tunit)
    end subroutine append_state_table

end module simple_flex_pca_states_pcg
