!@descr: experimental CTF/sigma-weighted PCG 3D reconstruction command, see
!  doc/policies/reconstruct3D_pcg_policy.md.
!  reconstruct3D_pcg reads an existing project and reconstructs ONE volume for
!  ONE state from its current project-carried orientation/CTF/shift
!  assignments, using the matrix-free (or, optionally, kernelized) PCG operator
!  in simple_reconstructor_pcg. Fixed-input contract:
!  orientations/shifts are inputs, not optimized; single state; nparts=1, no
!  even/odd split, no distributed execution; point-group symmetry is applied by
!  coordinate replication (c1 is a no-op); writes to a new
!  experimental execution directory; never writes anything back to the project;
!  does not enter through commander_volassemble or reuse its output filenames.
!  Kept in its own file, deliberately separate from production reconstruction.
module simple_commanders_reconstruct3D_pcg
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_commanders_api
use simple_reconstructor_pcg, only: reconstructor_pcg, pcg_solver_outcome, PCG_OP_KERNEL
use simple_matcher_ptcl_io,    only: prepimgbatch, discrete_read_imgbatch, killimgbatch
use simple_sigma2_files,       only: load_sigma2_groups
use simple_math_ft,            only: upsample_sigma2
implicit none

public :: commander_reconstruct3D_pcg
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_reconstruct3D_pcg
  contains
    procedure :: execute => exec_reconstruct3D_pcg
end type commander_reconstruct3D_pcg

contains

    subroutine exec_reconstruct3D_pcg( self, cline )
        class(commander_reconstruct3D_pcg), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        real,    parameter :: LAMBDA         = 1.0e-3
        integer, parameter :: MAXITS_DEFAULT = 30
        type(parameters)         :: params
        type(builder)            :: build
        type(reconstructor_pcg) :: pcgop
        type(pcg_solver_outcome) :: solver_outcome
        type(oris)               :: selection
        type(ori)                :: e
        type(ctfparams)          :: ctfparms
        type(image)              :: outvol
        type(string)             :: opmode
        integer, allocatable     :: pinds(:), tmpinds(:)
        real,    allocatable     :: sig2(:,:), x(:,:,:), rel_res_hist(:)
        ! one batch of observed planes, never the whole selection
        complex, allocatable     :: y_batch(:,:,:)
        integer :: nptcls, i, ii, iptcl, ibatch, batchlims(2), batchsz
        integer :: lims2(2,2), R, niters, funit, cnt, maxits, kfromto(2), nspace_dummy
        real    :: rtol, sdev_noise, edge_mean
        logical :: l_use_ctf, l_kernel, l_sig_loaded, l_norm_noise
        integer(timer_int_kind) :: t0, t1, t_acc
        real(timer_int_kind)    :: rt_tot, rt_setup, rt_solve
        ! cumulative setup checkpoints, differenced into a per-phase profile
        real(timer_int_kind)    :: c_build, c_sel, c_opnew, c_sigma, c_prep, c_accum
        ! accumulation sub-split: I/O read vs per-particle prep vs |T|^2 scatter
        real(timer_int_kind)    :: rt_read, rt_prep, rt_scatter

        t0 = tic()
        ! ---- cline defaults, mirroring exec_rec3D's shape. NOTE: projfile is
        !      deliberately NOT declared/required -- params%new ->
        !      setup_execution_context auto-discovers a unique *.simple in the
        !      cwd because the UI record declares requires_sp_project, exactly
        !      as reconstruct3D does. mkdir=yes is what creates the numbered
        !      execution directory and copies the project into it.
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',   'yes')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        if( .not. cline%defined('state')   ) call cline%set('state',   1.)
        if( .not. cline%defined('objfun')  ) call cline%set('objfun',  'euclid')
        if( .not. cline%defined('maxits')  ) call cline%set('maxits',  MAXITS_DEFAULT)
        ! nspace is only consulted by build_general_tbox's eulspace grid, which
        ! this operator never uses (it works from each particle's own continuous
        ! orientation, not a discretized reference grid). Keep a harmless even
        ! placeholder to satisfy internal setup and avoid build_refspiral's odd
        ! direction-count guard.
        if( .not. cline%defined('nspace') ) call cline%set('nspace', 2.)
        nspace_dummy = nint(cline%get_rarg('nspace'))
        if( is_odd(nspace_dummy) )then
            nspace_dummy = nspace_dummy + 1
            call cline%set('nspace', real(nspace_dummy))
            write(logfhandle,'(a,i0)') '>>> RECONSTRUCT3D_PCG: adjusted odd nspace to even placeholder=', nspace_dummy
        endif
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
        c_build = toc(t0)

        ! ---- fixed-input contract: reject the inputs
        !      this path deliberately does not support, rather than silently
        !      ignoring them and returning a volume that does not match what the
        !      user asked for ----
        if( params%nparts > 1 ) THROW_HARD('reconstruct3D_pcg is single-part; nparts>1 is not supported')

        ! pcgop and rtol are real parameters/type members, not cline-only flags:
        ! simple_args.f90 (the CLI allowlist) is GENERATED from the declarations
        ! in simple_parameters.f90, so a UI record alone is not enough -- an
        ! undeclared key is rejected as "argument is not allowed" before the
        ! commander ever runs.
        maxits = params%maxits
        rtol   = params%rtol
        if( maxits < 1 ) THROW_HARD('maxits must be at least 1; reconstruct3D_pcg')
        if( .not. ieee_is_finite(rtol) ) THROW_HARD('rtol must be finite; reconstruct3D_pcg')
        opmode = string(trim(params%pcgop))
        select case( opmode%to_char() )
            case( 'kernel', 'matrixfree' )
                ! valid
            case default
                THROW_HARD('unknown pcgop='//opmode%to_char()//'; use pcgop=kernel or pcgop=matrixfree')
        end select
        l_kernel = opmode .eq. 'kernel'
        l_norm_noise = params%msk > 0.5
        if( .not. l_norm_noise )then
            write(logfhandle,'(a)') '>>> RECONSTRUCT3D_PCG: WARNING, no mskdiam given, skipping noise normalization'
        endif

        ! ---- particle selection: sample4rec (state>0, updatecnt>0 when
        !      available), then filtered to exactly one state, since the
        !      contract reconstructs a single state ----
        nptcls = 0 ! sample4rec's nsamples argument is intent(inout)
        call build%spproj_field%sample4rec([params%fromp,params%top], nptcls, pinds)
        allocate(tmpinds(max(nptcls,1)))
        cnt = 0
        do i = 1, nptcls
            if( build%spproj_field%get_state(pinds(i)) == params%state )then
                cnt = cnt + 1
                tmpinds(cnt) = pinds(i)
            endif
        end do
        deallocate(pinds)
        ! NOTE: THROW_HARD is a preprocessor macro, so its argument must stay on
        ! one source line -- a Fortran continuation inside it does not compile.
        if( cnt < 1 ) THROW_HARD('no particles for state '//int2str(params%state)//'; reconstruct3D_pcg')
        allocate(pinds(cnt), source=tmpinds(1:cnt))
        nptcls = cnt
        deallocate(tmpinds)
        write(logfhandle,'(a,i0,a,i0)') '>>> RECONSTRUCT3D_PCG: particles selected = ', nptcls, &
            &' for state ', params%state
        c_sel = toc(t0)

        call pcgop%new(params%box, params%smpd, LAMBDA)
        lims2 = pcgop%get_lims2()
        R     = lims2(1,2)
        ! Point-group symmetry by coordinate replication: each plane pixel is
        ! scattered at all M operators R_i.S_g inside the accumulators. c1
        ! (nsym=1, identity) reproduces the asymmetric pass bit-for-bit; see
        ! simple_reconstructor_pcg%set_sym and pcg note section 2.
        call pcgop%set_sym(build%pgrpsyms)
        if( build%pgrpsyms%get_nsym() > 1 )then
            write(logfhandle,'(a,a,a,i0,a)') '>>> RECONSTRUCT3D_PCG: applying point group ', &
                &trim(params%pgrp), ' (', build%pgrpsyms%get_nsym(), ' operators) by coordinate replication'
        endif
        ! Solve for the density INSIDE mskdiam only. This is a constraint on the
        ! normal equations, not a mask applied to the output -- see
        ! reconstructor_pcg%set_mask. Skipped when mskdiam is absent, in which
        ! case the solve runs over the whole box as before.
        if( params%msk > 0.5 )then
            call pcgop%set_mask(params%msk)
            write(logfhandle,'(a,f7.1,a)') '>>> RECONSTRUCT3D_PCG: solving inside a soft sphere of radius ', &
                &params%msk, ' pixels'
        endif
        c_opnew = toc(t0)

        ! ---- sigma2, fetched the way flex_analysis does: euclid_sigma2 over
        !      the group star file, giving per-particle-per-shell spectra.
        !      Gated on objfun exactly as reconstruct3D gates ml_reg: objfun=cc
        !      needs no sigmas, objfun=euclid requires them. Failing loudly
        !      rather than silently falling back matters here -- a missing
        !      sigma2 file under euclid would quietly change what objective is
        !      actually being minimized. ----
        allocate(sig2(0:R,nptcls), source=1.0)
        l_use_ctf = .true.
        if( params%cc_objfun == OBJFUN_EUCLID )then
            ! discovery, carry-over and group loading are shared with
            ! flex_analysis; see simple_sigma2_files
            call load_sigma2_groups(params, build%pftc, build%esig, build%spproj_field, cline, l_sig_loaded)
            if( .not. l_sig_loaded )then
                THROW_HARD('objfun=euclid requires sigma2 files; none found, use objfun=cc instead')
            endif
            kfromto = build%esig%get_kfromto()
            do i = 1, nptcls
                call upsample_sigma2(kfromto(1), kfromto(2), &
                    &build%esig%sigma2_noise(kfromto(1):kfromto(2),pinds(i)), R, sig2(0:R,i))
            end do
            write(logfhandle,'(a)') '>>> RECONSTRUCT3D_PCG: objfun=euclid, per-particle sigma2 loaded'
        else
            write(logfhandle,'(a)') '>>> RECONSTRUCT3D_PCG: objfun=cc, no sigma weighting'
        endif
        c_sigma = toc(t0)

        ! ---- per-particle orientation/CTF/shift onto a selection oris. This is
        !      metadata only and needs no image I/O, so it runs as its own cheap
        !      pass: prep_particles must be complete before the streaming
        !      accumulation below, which reads the cached rotations and CTF. ----
        call selection%new(nptcls, .true.)
        call e%new(.false.)
        do i = 1, nptcls
            iptcl = pinds(i)
            call build%spproj_field%get_ori(iptcl, e)
            ctfparms = build%spproj%get_ctfparams(params%oritype, iptcl)
            call e%set_ctfvars(ctfparms)
            call e%set_shift(build%spproj_field%get_2Dshift(iptcl))
            call selection%set_ori(i, e)
        end do
        call pcgop%prep_particles(selection, use_ctf=l_use_ctf, sig2=sig2)
        c_prep = toc(t0)

        ! ---- streaming accumulation: read a batch, fold it into the two Fourier
        !      accumulators, discard it. The observed planes are needed for one
        !      thing only (forming the RHS) and for one pass only, so nothing
        !      proportional to nptcls is ever resident -- see
        !      doc/policies/reconstruct3D_pcg_policy.md section 1. Batch shape
        !      follows calc_3Drec. ----
        if( l_kernel ) write(logfhandle,'(a)') &
            &'>>> RECONSTRUCT3D_PCG: building kernelized (Toeplitz) normal operator'
        allocate(y_batch(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), MAXIMGBATCHSZ))
        sdev_noise = 0.  ! norm_noise takes it intent(inout)
        edge_mean  = 0.
        call pcgop%begin_accum
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        rt_read = 0.; rt_prep = 0.; rt_scatter = 0.
        do ibatch = 1, nptcls, MAXIMGBATCHSZ
            batchlims = [ibatch, min(nptcls, ibatch+MAXIMGBATCHSZ-1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            t_acc = tic()
            call discrete_read_imgbatch(params, build, nptcls, pinds, batchlims)
            rt_read = rt_read + toc(t_acc)
            t_acc = tic()
            do ii = 1, batchsz
                ! Normalize, TAPER, then transform -- the same three steps, in the
                ! same order, that production fuses into
                ! norm_noise_taper_edge_pad_fft (calc_3Drec's prep).
                !
                ! norm_noise: a bare fft() leaves each particle with its own
                ! arbitrary background level and variance, so the least-squares
                ! fit is dominated by whichever micrographs happen to be noisiest
                ! and every box contributes a large orientation-independent
                ! low-frequency term.
                !
                ! taper_edges_particle: a particle box is a hard-truncated crop
                ! out of a micrograph, so its border is a step discontinuity. The
                ! DFT reads that as a periodic wrap-around jump and answers with
                ! the familiar cross of energy along the axes -- which, scattered
                ! into the volume over thousands of orientations, accumulates as
                ! artefacts at the box edges of the reconstruction. Rolling the
                ! border off to the local edge mean over COSMSKHALFWIDTH pixels
                ! removes the discontinuity. Production has always done this; the
                ! first version of this commander did not.
                !
                ! build%lmsk is the box-sized logical disc built from mskdiam by
                ! build_general_tbox. Without mskdiam params%msk is 0, so lmsk
                ! selects nothing and the background statistics would be
                ! undefined -- fall back to a bare transform in that case.
                if( l_norm_noise )then
                    call build%imgbatch(ii)%norm_noise(build%lmsk, sdev_noise)
                    call build%imgbatch(ii)%taper_edges_particle(nint(COSMSKHALFWIDTH), edge_mean)
                endif
                call build%imgbatch(ii)%fft()
                y_batch(:,:,ii) = pcgop%extract_native_plane(build%imgbatch(ii))
            end do
            rt_prep = rt_prep + toc(t_acc)
            ! one call folds this batch into BOTH accumulators: the RHS (which
            ! needs the planes) and the |T|^2 sampling density (which does not)
            t_acc = tic()
            call pcgop%accumulate_batch(y_batch, batchsz, batchlims(1))
            rt_scatter = rt_scatter + toc(t_acc)
        end do
        call killimgbatch(build)
        deallocate(y_batch)
        c_accum = toc(t0)

        ! ---- close accumulation: folds the RHS once and derives the
        !      preconditioner and, when requested, the Gram kernel. Both come
        !      from the same |T_i|^2 accumulator, so neither costs a second
        !      particle pass. ----
        call pcgop%end_accum(l_kernel)
        if( l_kernel ) call pcgop%set_op_mode(PCG_OP_KERNEL)
        rt_setup = toc(t0)
        write(logfhandle,'(a,f9.2,a)') '>>> RECONSTRUCT3D_PCG: setup time = ', rt_setup, ' s'
        ! per-phase profile: the absT2/rim work the win_wraps and LUT changes
        ! touched lives entirely in the accumulation line -- watch that one
        write(logfhandle,'(a)')      '>>> RECONSTRUCT3D_PCG: setup profile (s)'
        write(logfhandle,'(a,f9.2)') '    build_general_tbox   : ', c_build
        write(logfhandle,'(a,f9.2)') '    particle selection   : ', c_sel   - c_build
        write(logfhandle,'(a,f9.2)') '    operator allocation  : ', c_opnew - c_sel
        write(logfhandle,'(a,f9.2)') '    sigma2 load          : ', c_sigma - c_opnew
        write(logfhandle,'(a,f9.2)') '    prep_particles       : ', c_prep  - c_sigma
        write(logfhandle,'(a,f9.2)') '    accumulation (absT2) : ', c_accum - c_prep
        write(logfhandle,'(a,f9.2)') '      read I/O           : ', rt_read
        write(logfhandle,'(a,f9.2)') '      prep (norm/tap/fft): ', rt_prep
        write(logfhandle,'(a,f9.2)') '      scatter |T|^2      : ', rt_scatter
        write(logfhandle,'(a,f9.2)') '    end_accum (FFT+kern)  : ', rt_setup - c_accum

        ! ---- solve ----
        t1 = tic()
        allocate(x(params%box,params%box,params%box), source=0.0)
        call pcgop%solve_accum(x, maxits=maxits, rtol=rtol, &
            &rel_res_hist=rel_res_hist, niters=niters, outcome=solver_outcome)
        rt_solve = toc(t1)
        rt_tot   = toc(t0)

        ! ---- output: experimental volume + diagnostics table. Deliberately NOT
        !      any simple_refine3D_fnames name, and no spproj%write or anything
        !      else that would mutate the project. ----
        call outvol%new([params%box,params%box,params%box], params%smpd)
        call outvol%set_rmat(x, .false.)
        call outvol%write(string('reconstruct3D_pcg_state'//int2str_pad(params%state,2)//'.mrc'))
        call fopen(funit, string('reconstruct3D_pcg_diagnostics.txt'), status='replace', action='write')
        write(funit,'(a,i0)')    'nptcls=',           nptcls
        write(funit,'(a,i0)')    'state=',            params%state
        write(funit,'(a,i0)')    'box=',              params%box
        write(funit,'(a,a)')     'objfun=',           trim(params%objfun)
        write(funit,'(a,a)')     'operator=',         opmode%to_char()
        write(funit,'(a,i0)')    'nthr=',             params%nthr
        write(funit,'(a,i0)')    'niters=',           niters
        write(funit,'(a,i0)')    'requested_maxits=', solver_outcome%requested_maxits
        write(funit,'(a,a)')     'stop_reason=',      trim(solver_outcome%stop_reason)
        write(funit,'(a,l1)')    'converged=',        solver_outcome%converged
        write(funit,'(a,es14.6)') 'initial_rel_resid_l2=', solver_outcome%initial_rel_residual
        write(funit,'(a,es14.6)') 'final_rel_resid_l2=',   solver_outcome%final_rel_residual
        write(funit,'(a,es14.6)') 'final_rel_update=',     solver_outcome%final_rel_update
        write(funit,'(a,f12.3)') 'setup_time_sec=',   rt_setup
        write(funit,'(a,f12.3)') 'solve_time_sec=',   rt_solve
        write(funit,'(a,f12.3)') 'total_time_sec=',   rt_tot
        if( niters > 0 ) write(funit,'(a,f12.4)') 'sec_per_iteration=', rt_solve/real(niters)
        do i = 1, niters
            ! true relative residual ||r||_2/||b||_2, not the preconditioned norm
            write(funit,'(a,i0,a,es14.6)') 'iter', i, '_rel_resid_l2=', rel_res_hist(i)
        end do
        call fclose(funit)
        write(logfhandle,'(a,i0,a)')   '>>> RECONSTRUCT3D_PCG: PCG ran ', niters, ' iterations'
        write(logfhandle,'(a,f9.2,a)') '>>> RECONSTRUCT3D_PCG: total time = ', rt_tot, ' s'

        call outvol%kill
        call pcgop%kill
        call selection%kill
        call e%kill
        call opmode%kill
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_RECONSTRUCT3D_PCG NORMAL STOP ****')
    end subroutine exec_reconstruct3D_pcg

end module simple_commanders_reconstruct3D_pcg
