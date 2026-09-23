!@descr: flex_pca state maps on the reconstruct3D PCG backend (rec_backend=pcg): the kernel weight of a
!  particle enters through its noise model as sigma2/w, so the right-hand side and the density carry it
!  identically and the preconditioner, Gram kernel, ridge scale and raw artifacts follow without change
!  (doc/implementation_notes/flex_pca_envelope_support.md, section 3.3)
module simple_flex_pca_rec3D_pcg
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_builder,           only: builder
use simple_parameters,        only: parameters
use simple_image,             only: image
use simple_reconstructor_pcg, only: reconstructor_pcg, pcg_solver_outcome, PCG_OP_KERNEL, PCG_STOP_INDEFINITE
use simple_flex_pca_pcg,      only: flex_pcg_support_volume, flex_mskfile_set
use simple_matcher_ptcl_io,   only: prepimgbatch, discrete_read_imgbatch, discrete_read_imgbatch_source, &
    &prep_rec_observation
use simple_math_ft,           only: resample_sigma2
use simple_estimate_ssnr,     only: fsc2optlp_sub, get_resolution
use simple_flex_pca_rounds,   only: flex_pca_rounds, flex_pca_part_path
use simple_flex_pca_util,     only: flex_pca_write_state
!$ use omp_lib, only: omp_get_max_active_levels, omp_set_max_active_levels, omp_set_num_threads
implicit none

public :: reconstruct_flex_weighted_states_pcg
public :: flex_pcg_state_raw_fname, flex_pcg_state_provenance
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

contains

    !> One state per weight column, per halfset when l_fuse: workers accumulate the raw (B,D) of their
    !! particle range and publish one artifact per (state, half, part); the master (or the shared-memory
    !! process, which accumulates in place) reduces in part order, finalizes the kernel and solves each
    !! half cold on the spherical support. Delivered maps are window*u; nothing masks them afterwards.
    subroutine reconstruct_flex_weighted_states_pcg( params, build, pinds, state_weights, nstates, l_fuse, &
        &box_rec, smpd_rec, outvol_even, outvol_odd, rounds )
        class(parameters),      intent(inout) :: params
        class(builder),         intent(inout) :: build
        integer,                intent(in)    :: pinds(:), nstates, box_rec
        real,                   intent(in)    :: state_weights(:,:), smpd_rec
        logical,                intent(in)    :: l_fuse
        type(string), optional, intent(in)    :: outvol_even, outvol_odd
        class(flex_pca_rounds), intent(inout) :: rounds
        type(reconstructor_pcg)  :: pcgop, pcgops(0:1)
        type(half_solve_rec)     :: sol(0:1)
        type(image), allocatable :: maps_e(:), maps_o(:)
        type(image)   :: img_e, img_o, img_c, state_img, fsc_e, fsc_o, envimg
        type(string)  :: outvol_bak, fname, state_vol_fname
        character(len=256) :: provenance
        character(len=8)   :: envval
        real,    allocatable :: fsc_eo(:), res_arr(:), filt_half(:), filt_merged(:)
        integer, allocatable :: nsel_e(:), nsel_o(:)
        real    :: msk_rec, fsc05, fsc0143, kc_lp
        integer :: state, eo, neo, nsel, nsel_part, ipart, filtsz, k_lp, iv, tunit
        integer :: envlen, envstat, nthr_half, nthr_pair, prev_levels
        logical :: l_state_eofilt, l_state_filt, l_eo_fsc, l_pair
        if( size(pinds) < 1 .or. nstates < 1 ) THROW_HARD('invalid flex PCG state reconstruction dimensions')
        if( any(shape(state_weights) /= [size(pinds),nstates]) ) THROW_HARD('flex PCG weighted state table mismatch')
        if( l_fuse .and. .not. rounds%is_worker() )then
            if( .not. (present(outvol_even) .and. present(outvol_odd)) ) &
                &THROW_HARD('flex PCG halfset delivery needs outvol_even/outvol_odd')
        endif
        ! the same spherical support as the gridding delivery mask, capped at the box edge
        msk_rec    = min(0.5*params%mskdiam/smpd_rec, real(box_rec/2) - COSMSKHALFWIDTH - 1.)
        provenance = flex_pcg_state_provenance(params, box_rec, smpd_rec)
        neo = 0
        if( l_fuse ) neo = 1
        ! ---- worker: raw (B,D) per state and half of this part ----
        if( rounds%is_worker() )then
            call prepimgbatch(params, build, MAXIMGBATCHSZ)
            do state = 1, nstates
                do eo = 0, neo
                    call accumulate_state_half(state, eo, pcgop, nsel)
                    fname = flex_pcg_state_raw_fname(params, params%part, state, eo)
                    ! the part count comes from the command line: rounds%nparts() is the master's plan (1 on a worker)
                    call pcgop%write_raw_accum(fname, state, eo, params%part, max(1,params%nparts), nsel, provenance)
                    write(logfhandle,'(A,I0,A,I0,A,I0,A,I0)') '>>> FLEX_PCA PCG STATE RAW: part ', params%part, &
                        &' state ', state, ' half ', eo, ' particles ', nsel
                    call pcgop%kill
                    call fname%kill
                end do
            end do
            call flush(logfhandle)
            return
        endif
        ! ---- master (reduce the parts) or shared memory (accumulate in place), then solve ----
        allocate(maps_e(nstates))
        allocate(nsel_e(nstates), source=0)
        if( l_fuse )then
            allocate(maps_o(nstates))
            allocate(nsel_o(nstates), source=0)
        endif
        if( .not. rounds%distributed() ) call prepimgbatch(params, build, MAXIMGBATCHSZ)
        ! ---- delivery policy (applied per state right after its halves are solved): the same
        ! per-state eo-FSC filter as the gridding path, on the windowed solutions as they come out
        ! of the solve (no background removal, no second mask) ----
        l_state_eofilt = .false.
        call get_environment_variable('SIMPLE_COV_STATE_EOFILT', envval, envlen, envstat)
        if( envstat == 0 .and. envlen > 0 )then
            if( trim(adjustl(envval)) == '1' ) l_state_eofilt = .true.
        endif
        l_state_filt = .true.
        call get_environment_variable('SIMPLE_COV_STATE_FILT', envval, envlen, envstat)
        if( envstat == 0 .and. envlen > 0 )then
            if( trim(adjustl(envval)) == '0' ) l_state_filt = .false.
        endif
        filtsz = fdim(box_rec) - 1
        allocate(fsc_eo(filtsz), filt_half(filtsz), filt_merged(filtsz))
        outvol_bak = params%outvol
        ! even and odd of a state solve concurrently on half the threads each (the reconstruct3D PCG
        ! master's pairing); the accumulation/reduction before them stays serial (shared image batch, disk)
        ! total capped at the reconstruct3D PCG master's PCG_MASTER_NTHR_CAP (32): the FFT scaling of a
        ! single solve flattens well before the boosted master budget, and the cap bounds the thread stacks
        nthr_pair = min(nthr_glob, FLEX_PCG_PAIR_NTHR_CAP)
        l_pair    = l_fuse .and. nthr_pair >= 2
        nthr_half = nthr_pair
        if( l_pair ) nthr_half = max(1, nthr_pair/2)
        do state = 1, nstates
            do eo = 0, neo
                if( rounds%distributed() )then
                    call new_state_operator(pcgops(eo))
                    write(logfhandle,'(A,I0,A,I0,A,I0,A)') '>>> FLEX_PCA PCG STATE ', state, ' half ', eo, &
                        &': reducing ', max(1,params%nparts), ' raw parts and finalizing the kernel'
                    call flush(logfhandle)
                    call pcgops(eo)%begin_reduction
                    nsel = 0
                    do ipart = 1, max(1,params%nparts)
                        fname = flex_pcg_state_raw_fname(params, ipart, state, eo)
                        if( .not. file_exists(fname) ) THROW_HARD('missing flex PCG raw part: '//fname%to_char())
                        call pcgops(eo)%add_raw_accum(fname, state, eo, ipart, max(1,params%nparts), provenance, nsel_part)
                        nsel = nsel + nsel_part
                        call del_file(fname)
                        call fname%kill
                    end do
                else
                    call accumulate_state_half(state, eo, pcgops(eo), nsel)
                endif
                if( eo == 0 )then
                    nsel_e(state) = nsel
                else
                    nsel_o(state) = nsel
                endif
            end do
            if( l_pair )then
                prev_levels = 1
                !$ prev_levels = omp_get_max_active_levels()
                !$ call omp_set_max_active_levels(max(2, prev_levels))
                !$omp parallel sections num_threads(2) default(shared)
                !$omp section
                !$ call omp_set_num_threads(nthr_half)
                call solve_state_half(pcgops(0), nsel_e(state), maps_e(state), sol(0))
                !$omp section
                !$ call omp_set_num_threads(nthr_half)
                call solve_state_half(pcgops(1), nsel_o(state), maps_o(state), sol(1))
                !$omp end parallel sections
                !$ call omp_set_max_active_levels(prev_levels)
                !$ call omp_set_num_threads(nthr_glob)
            else
                call solve_state_half(pcgops(0), nsel_e(state), maps_e(state), sol(0))
                if( l_fuse ) call solve_state_half(pcgops(1), nsel_o(state), maps_o(state), sol(1))
            endif
            ! serial: the outcome checks, the log lines and the table, in even/odd order
            do eo = 0, neo
                if( eo == 0 )then
                    call report_state_half(state, eo, nsel_e(state), sol(eo))
                else
                    call report_state_half(state, eo, nsel_o(state), sol(eo))
                endif
                call pcgops(eo)%kill
            end do
            ! deliver this state now: its maps go to disk as soon as both halves are solved,
            ! and its half maps are freed (one pair resident, not nstates)
            call deliver_state(state)
        end do
        call build%spproj%write_segment_inside('out', params%projfile)
        params%outvol = outvol_bak
        call state_vol_fname%kill
        call outvol_bak%kill
        deallocate(fsc_eo, filt_half, filt_merged, maps_e, nsel_e)
        if( l_fuse ) deallocate(maps_o, nsel_o)

    contains

        !> combine even+odd, take the eo-FSC, filter and write the three maps of one solved state
        subroutine deliver_state( state )
            integer, intent(in) :: state
            if( l_fuse )then
                call img_e%copy(maps_e(state))
                call img_o%copy(maps_o(state))
                call img_c%copy(img_e)
                call img_c%add(img_o)
                call img_c%mul(0.5)
                call fsc_e%copy(img_e)
                call fsc_o%copy(img_o)
                call fsc_e%fft
                call fsc_o%fft
                call fsc_e%fsc(fsc_o, fsc_eo)
                call fsc_e%kill
                call fsc_o%kill
                l_eo_fsc = any(fsc_eo > 0.143)
                if( l_eo_fsc )then
                    res_arr = img_e%get_res()
                    call get_resolution(fsc_eo, res_arr, fsc05, fsc0143)
                    if( l_state_eofilt )then
                        call fsc2optlp_sub(filtsz, fsc_eo, filt_half,   merged=.false.)
                        call fsc2optlp_sub(filtsz, fsc_eo, filt_merged, merged=.true.)
                        write(logfhandle,'(A,I3,A,F7.2,A,F7.2,A)') '>>> FLEX STATE (PCG) eo-FSC state=', &
                            &state,'  res(0.143)=',fsc0143,' A  res(0.5)=',fsc05,' A -- per-state optimal filter applied'
                    else
                        kc_lp = real(box_rec) * smpd_rec / fsc0143
                        do k_lp = 1, filtsz
                            filt_merged(k_lp) = 1.0 / (1.0 + (real(k_lp)/max(kc_lp,1.0))**8)
                        end do
                        filt_half = filt_merged
                        write(logfhandle,'(A,I3,A,F7.2,A,F7.2,A)') '>>> FLEX STATE (PCG) eo-FSC state=', &
                            &state,'  res(0.143)=',fsc0143,' A  res(0.5)=',fsc05,' A -- low-pass at the state eo-FSC(0.143) applied'
                    endif
                    deallocate(res_arr)
                else
                    write(logfhandle,'(A,I3,A)') '>>> FLEX STATE (PCG) eo-FSC state=',state, &
                        &'  unmeasurable (no shell above 0.143); delivered unfiltered'
                endif
                call flush(logfhandle)
                do iv = 1, 3
                    select case(iv)
                    case(1); call state_img%copy(img_c); params%outvol = outvol_bak
                    case(2); call state_img%copy(img_e); params%outvol = outvol_even
                    case(3); call state_img%copy(img_o); params%outvol = outvol_odd
                    end select
                    if( l_state_filt .and. l_eo_fsc )then
                        if( iv == 1 )then
                            call state_img%apply_filter(filt_merged)
                        else
                            call state_img%apply_filter(filt_half)
                        endif
                    endif
                    call flex_pca_write_state(params, state_img, state, state_vol_fname)
                    if( iv == 1 )then
                        call build%spproj%add_vol2os_out(state_vol_fname, state_img%get_smpd(), state, 'vol_flex', &
                            &box=state_img%get_box())
                    endif
                    call state_img%kill
                end do
                call img_e%kill
                call img_o%kill
                call img_c%kill
                call maps_o(state)%kill
            else
                ! single-set state: one solve over every weighted particle, delivered unfiltered
                params%outvol = outvol_bak
                call state_img%copy(maps_e(state))
                write(logfhandle,'(A,I3,A)') '>>> FLEX STATE (PCG) state=', state, '  single-set solve delivered unfiltered'
                call flex_pca_write_state(params, state_img, state, state_vol_fname)
                call build%spproj%add_vol2os_out(state_vol_fname, state_img%get_smpd(), state, 'vol_flex', &
                    &box=state_img%get_box())
                call state_img%kill
            endif
            call maps_e(state)%kill
        end subroutine deliver_state

        !> operator of one (state, half) solve at the reconstruction box: relative ridge, symmetry,
        !! spherical support (installed before accumulation so the RHS is projected in end_accum)
        subroutine new_state_operator( op )
            type(reconstructor_pcg), intent(inout) :: op
            if( rounds%is_worker() )then
                call op%new(box_rec, smpd_rec)
            else
                call op%new(box_rec, smpd_rec, fft_nthreads=nthr_half)
            endif
            call op%set_sym(build%pgrpsyms)
            call op%set_lambda_relative(FLEX_PCG_LAMBDA_REL)
            if( flex_mskfile_set(params) )then
                ! step 2: the envelope resampled to the state box is the solve support
                call flex_pcg_support_volume(params, box_rec, smpd_rec, envimg, 'state box')
                call op%set_mask_volume(envimg)
                call envimg%kill
            else
                call op%set_mask(msk_rec)
            endif
        end subroutine new_state_operator

        !> accumulate the raw (B,D) of one (state, half): the selection keeps particles above the
        !! weight floor, in the wanted halfset when l_fuse, and passes sigma2/w as their noise model
        subroutine accumulate_state_half( state_here, eo_here, op, nsel_out )
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
            call new_state_operator(op)
            allocate(sel(size(pinds)), wsel(size(pinds)))
            n = 0
            do i = 1, size(pinds)
                iptcl = pinds(i)
                w     = state_weights(i, state_here)
                if( w < FLEX_PCG_WEIGHT_FLOOR ) cycle
                if( build%spproj_field%get_state(iptcl) <= 0 ) cycle
                if( l_fuse )then
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
            crop_factor = real(box_rec) / real(params%box)
            call selection%new(n, .true.)
            call orientation%new(.false.)
            do i = 1, n
                iptcl = sel(i)
                call build%spproj_field%get_ori(iptcl, orientation)
                ctfparms      = build%spproj%get_ctfparams(params%oritype, iptcl)
                ctfparms%smpd = smpd_rec
                shift         = build%spproj_field%get_2Dshift(iptcl) * crop_factor
                call orientation%set_ctfvars(ctfparms)
                call orientation%set_shift(shift)
                call selection%set_ori(i, orientation)
            end do
            call op%prep_particles(selection, use_ctf=.true., sig2=sig2)
            allocate(y_batch(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), MAXIMGBATCHSZ))
            call op%begin_accum
            call obs%new([box_rec,box_rec,1], smpd_rec)
            do ibatch = 1, n, MAXIMGBATCHSZ
                batchlims = [ibatch, min(n, ibatch+MAXIMGBATCHSZ-1)]
                batchsz   = batchlims(2) - batchlims(1) + 1
                if( params%l_ptcl_src_den )then
                    call discrete_read_imgbatch_source(params, build, 'den', n, sel, batchlims, build%imgbatch(:batchsz))
                else
                    call discrete_read_imgbatch(params, build, n, sel, batchlims)
                endif
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

        !> finalize (kernel + preconditioner) and solve one (state, half) cold; the map is window*u.
        !! Runs inside the concurrent even/odd sections: no I/O, no logging, no THROW here -- the
        !! outcome is recorded and judged by report_state_half afterwards
        subroutine solve_state_half( op, nsel_here, vol, rec )
            type(reconstructor_pcg), intent(inout) :: op
            integer,                 intent(in)    :: nsel_here
            type(image),             intent(inout) :: vol
            type(half_solve_rec),    intent(inout) :: rec
            real, allocatable :: x(:,:,:)
            integer(timer_int_kind) :: t_solve
            rec%l_solved = .false.
            rec%niters   = 0
            call vol%new([box_rec,box_rec,box_rec], smpd_rec)
            if( nsel_here == 0 ) return
            t_solve = tic()
            call op%end_accum(.true.)
            call op%set_op_mode(PCG_OP_KERNEL)
            allocate(x(box_rec,box_rec,box_rec), source=0.0)
            call op%solve_accum(x, maxits=FLEX_PCG_STATE_MAXITS, rtol=params%rtol, niters=rec%niters, &
                &outcome=rec%outcome)
            rec%l_solved = all(ieee_is_finite(x))
            if( rec%l_solved ) call vol%set_rmat(x, .false.)
            rec%seconds = real(toc(t_solve))
            deallocate(x)
        end subroutine solve_state_half

        !> the serial tail of one half solve: the failure checks, the log line and the table row
        subroutine report_state_half( state_here, eo_here, nsel_here, rec )
            integer,              intent(in) :: state_here, eo_here, nsel_here
            type(half_solve_rec), intent(in) :: rec
            character(len=4) :: half
            half = 'all '
            if( l_fuse ) half = merge('odd ', 'even', eo_here == 1)
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
            call append_state_table(state_here, half, nsel_here, rec%niters, rec%outcome)
        end subroutine report_state_half

        !> one line per (state, half) solve; the run's record of what the PCG delivered
        subroutine append_state_table( state_here, half, nsel_here, niters_here, res )
            integer,                  intent(in) :: state_here, nsel_here, niters_here
            character(len=*),         intent(in) :: half
            type(pcg_solver_outcome), intent(in) :: res
            integer :: io_stat
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
                &res%final_rel_update, trim(res%stop_reason), box_rec, msk_rec
            close(tunit)
        end subroutine append_state_table

    end subroutine reconstruct_flex_weighted_states_pcg

    !> raw artifact of one worker's (state, half) accumulation, under the part directory
    function flex_pcg_state_raw_fname( params, part, state, eo ) result( fname )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: part, state, eo
        type(string) :: fname
        character(len=2) :: half
        half = '_e'
        if( eo == 1 ) half = '_o'
        fname = flex_pca_part_path('flex_pca_pcgraw_part'//int2str_pad(part,max(1,params%numlen))// &
            &'_'//int2str_pad(state,2)//half//'.bin')
    end function flex_pcg_state_raw_fname

    !> identity every raw artifact of one run carries; the master refuses a part that disagrees
    function flex_pcg_state_provenance( params, box_rec, smpd_rec ) result( provenance )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: box_rec
        real,              intent(in) :: smpd_rec
        character(len=256) :: provenance
        provenance = 'flexpcg-v1|pgrp='//trim(params%pgrp)//'|objfun='//trim(params%objfun)// &
            &'|ptcl_src='//trim(params%ptcl_src)//'|box='//trim(int2str(params%box))// &
            &'|smpd='//trim(real2str(params%smpd))//'|box_rec='//trim(int2str(box_rec))// &
            &'|smpd_rec='//trim(real2str(smpd_rec))//'|mskdiam='//trim(real2str(params%mskdiam))// &
            &'|ctf='//trim(params%ctf)
    end function flex_pcg_state_provenance

end module simple_flex_pca_rec3D_pcg
