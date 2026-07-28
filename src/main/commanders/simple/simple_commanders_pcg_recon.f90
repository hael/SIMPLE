!@descr: experimental CTF/sigma-weighted PCG 3D reconstruction command, see
!  doc/implementation_notes/ctf_sigma_weighted_pcg_reconstruction.md.
!  reconstruct3D_pcg reads an existing project and reconstructs ONE volume for
!  ONE state from its current project-carried orientation/CTF/shift
!  assignments, using the matrix-free (or, optionally, kernelized) PCG operator
!  in simple_pcg_reconstruction. Phase-1 fixed-input contract (note section 2):
!  orientations/shifts are inputs, not optimized; single state; nparts=1, no
!  even/odd split, no symmetry, no distributed execution; writes to a new
!  experimental execution directory; never writes anything back to the project;
!  does not enter through commander_volassemble or reuse its output filenames.
!  Kept in its own file, deliberately separate from production reconstruction.
module simple_commanders_pcg_recon
use simple_commanders_api
use simple_pcg_reconstruction, only: pcg_reconstruction, PCG_OP_KERNEL
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
        real,    parameter :: RTOL_DEFAULT   = 1.0e-3
        type(parameters)         :: params
        type(builder)            :: build
        type(pcg_reconstruction) :: pcgop
        type(oris)               :: selection
        type(ori)                :: e
        type(ctfparams)          :: ctfparms
        type(image)              :: outvol
        type(string)             :: opmode
        integer, allocatable     :: pinds(:), tmpinds(:)
        real,    allocatable     :: sig2(:,:), x(:,:,:), rel_res_hist(:)
        complex, allocatable     :: y_planes(:,:,:)
        integer :: nptcls, i, ii, iptcl, ibatch, batchlims(2), batchsz
        integer :: lims2(2,2), R, niters, funit, cnt, maxits, kfromto(2), nspace_dummy
        real    :: rtol
        logical :: l_use_ctf, l_kernel, l_sig_loaded
        integer(timer_int_kind) :: t0, t1
        real(timer_int_kind)    :: rt_tot, rt_setup, rt_solve

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

        maxits = MAXITS_DEFAULT
        if( cline%defined('maxits') ) maxits = nint(cline%get_rarg('maxits'))
        rtol = RTOL_DEFAULT
        if( cline%defined('rtol')   ) rtol   = cline%get_rarg('rtol')
        opmode = string('matrixfree')
        if( cline%defined('pcgop')  ) opmode = cline%get_carg('pcgop')
        l_kernel = opmode .eq. 'kernel'

        ! ---- particle selection: sample4rec (state>0, updatecnt>0 when
        !      available), then filtered to exactly one state, since the phase-1
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

        call pcgop%new(params%box, params%smpd, LAMBDA)
        lims2 = pcgop%get_lims2()
        R     = lims2(1,2)

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

        ! ---- per-particle orientation/CTF/shift onto a selection oris, plus
        !      the particle planes, read in batches exactly as calc_3Drec does
        !      rather than opening the stack once per particle ----
        call selection%new(nptcls, .true.)
        allocate(y_planes(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), nptcls))
        call e%new(.false.)
        call prepimgbatch(params, build, MAXIMGBATCHSZ)
        do ibatch = 1, nptcls, MAXIMGBATCHSZ
            batchlims = [ibatch, min(nptcls, ibatch+MAXIMGBATCHSZ-1)]
            batchsz   = batchlims(2) - batchlims(1) + 1
            call discrete_read_imgbatch(params, build, nptcls, pinds, batchlims)
            do ii = 1, batchsz
                i     = batchlims(1) + ii - 1
                iptcl = pinds(i)
                call build%spproj_field%get_ori(iptcl, e)
                ctfparms = build%spproj%get_ctfparams(params%oritype, iptcl)
                call e%set_ctfvars(ctfparms)
                call e%set_shift(build%spproj_field%get_2Dshift(iptcl))
                call selection%set_ori(i, e)
                call build%imgbatch(ii)%fft()
                y_planes(:,:,i) = pcgop%extract_native_plane(build%imgbatch(ii))
            end do
        end do
        call killimgbatch(build)

        ! ---- operator setup ----
        call pcgop%prep_particles(selection, use_ctf=l_use_ctf, sig2=sig2)
        call pcgop%build_precond
        if( l_kernel )then
            write(logfhandle,'(a)') '>>> RECONSTRUCT3D_PCG: building kernelized (Toeplitz) normal operator'
            call pcgop%build_kernel
            call pcgop%set_op_mode(PCG_OP_KERNEL)
        endif
        rt_setup = toc(t0)
        write(logfhandle,'(a,f9.2,a)') '>>> RECONSTRUCT3D_PCG: setup time = ', rt_setup, ' s'

        ! ---- solve ----
        t1 = tic()
        allocate(x(params%box,params%box,params%box), source=0.0)
        call pcgop%solve(y_planes, x, maxits=maxits, rtol=rtol, &
            &rel_res_hist=rel_res_hist, niters=niters)
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
        write(funit,'(a,f12.3)') 'setup_time_sec=',   rt_setup
        write(funit,'(a,f12.3)') 'solve_time_sec=',   rt_solve
        write(funit,'(a,f12.3)') 'total_time_sec=',   rt_tot
        if( niters > 0 ) write(funit,'(a,f12.4)') 'sec_per_iteration=', rt_solve/real(niters)
        do i = 1, niters
            write(funit,'(a,i0,a,es14.6)') 'iter', i, '_rel_resid=', rel_res_hist(i)
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

end module simple_commanders_pcg_recon
