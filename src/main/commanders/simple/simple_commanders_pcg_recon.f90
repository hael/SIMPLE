!@descr: milestone 2 of doc/implementation_notes/ctf_sigma_weighted_pcg_reconstruction.md.
!  reconstruct3D_pcg reads an existing project and reconstructs ONE volume for
!  ONE state from its current project-carried orientation/CTF/shift
!  assignments, using the experimental matrix-free PCG operator
!  (simple_pcg_reconstruction.f90). Phase-1 fixed-input contract (note section
!  2): orientations/shifts are inputs, not optimized; single state selected;
!  nparts=1, no even/odd split, no symmetry, no distributed execution; writes
!  to a new experimental execution directory; never writes anything back to
!  the project; does not enter through commander_volassemble or reuse its
!  output filenames. Kept in its own file, deliberately separate from
!  production reconstruction commanders.
module simple_commanders_pcg_recon
use simple_commanders_api
use simple_pcg_reconstruction, only: pcg_reconstruction
use simple_matcher_ptcl_io,    only: prepimgbatch, discrete_read_imgbatch
use simple_sigma2_binfile,     only: sigma2_binfile
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
        real, parameter :: LAMBDA = 1.0e-3
        type(parameters)         :: params
        type(builder)             :: build
        type(pcg_reconstruction) :: pcgop
        type(oris)                :: selection
        type(ori)                 :: e
        type(ctfparams)           :: ctfparms
        type(sigma2_binfile)      :: s2file
        type(string)               :: sigma2_fname
        type(image)                :: outvol
        integer, allocatable      :: pinds(:), state_filtered(:)
        real,    allocatable      :: sig2_per_ptcl(:,:), sig2arr(:), x(:,:,:), rel_res_hist(:)
        complex, allocatable      :: y_planes(:,:,:)
        integer :: nptcls, i, ii, iptcl, ibatch, batchlims(2), batchsz, R, lims2(2,2)
        integer :: niters, funit, cnt, kk, kclamp
        real    :: shift(2)
        logical :: l_have_sigma2
        integer(timer_int_kind) :: t0
        real(timer_int_kind)    :: rt_tot

        t0 = tic()
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',   'yes')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        if( .not. cline%defined('state')   ) call cline%set('state',   1.)
        ! eulspace/nspace is built by init_params_and_build_general_tbox(do3d=.true.)
        ! but never used by this operator (each particle's own continuous
        ! orientation is used directly, not a discretized reference grid) --
        ! set a dummy value so it never becomes a user-visible requirement.
        if( .not. cline%defined('nspace')  ) call cline%set('nspace',  1.)
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)

        ! ---- particle selection: sample4rec (state>0, updatecnt>0), then
        !      filter to exactly one state (phase-1 contract: single state) ----
        nptcls = 0 ! sample4rec's nsamples argument is intent(inout)
        call build%spproj_field%sample4rec([params%fromp,params%top], nptcls, pinds)
        allocate(state_filtered(nptcls))
        cnt = 0
        do i = 1, nptcls
            if( build%spproj_field%get_state(pinds(i)) == params%state )then
                cnt = cnt + 1
                state_filtered(cnt) = pinds(i)
            endif
        end do
        deallocate(pinds)
        allocate(pinds(cnt), source=state_filtered(1:cnt))
        nptcls = cnt
        deallocate(state_filtered)
        if( nptcls < 1 ) THROW_HARD('no particles selected for state='//int2str(params%state)//'; reconstruct3D_pcg')
        write(logfhandle,'(a,i0,a)') '>>> RECONSTRUCT3D_PCG: ', nptcls, ' particles selected'

        ! ---- operator + selection oris (rotation + ctf + shift, per particle,
        !      read-only from the project -- nothing here is ever written back) ----
        call pcgop%new(params%box, params%smpd, LAMBDA)
        lims2 = pcgop%get_lims2()
        R     = lims2(1,2)
        call selection%new(nptcls, .true.)
        allocate(y_planes(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), nptcls))

        ! ---- batched particle image read, matching production's own I/O
        !      pattern (calc_3Drec in simple_matcher_3Drec.f90) instead of
        !      opening/closing the stack file once per particle ----
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
                shift    = build%spproj_field%get_2Dshift(iptcl)
                call e%set_ctfvars(ctfparms)
                call e%set_shift(shift)
                call selection%set_ori(i, e)
                call build%imgbatch(ii)%fft()
                y_planes(:,:,i) = pcgop%extract_native_plane(build%imgbatch(ii))
            end do
        end do

        ! ---- sigma2: gated on objfun, matching reconstruct3D's own convention
        !      (objfun=cc needs no sigmas; objfun=euclid requires them -- not
        !      a silent fallback, since a missing sigma2 file under euclid
        !      would silently change what the objective actually is). Read
        !      via the lightweight, project/pftc-independent sigma2_binfile
        !      and averaged over the selection into one shared profile --
        !      true per-particle sigma2 is a later, explicitly deferred
        !      refinement (same scope decision as milestone 1). ----
        l_have_sigma2 = .false.
        if( trim(params%objfun) == 'euclid' )then
            sigma2_fname = string(SIGMA2_FBODY//int2str_pad(params%part,params%numlen)//'.dat')
            if( .not. file_exists(sigma2_fname) )then
                THROW_HARD('objfun=euclid requires a sigma2 file ('//sigma2_fname%to_char()//'), none found for this partition; reconstruct3D_pcg')
            endif
            call s2file%new_from_file(sigma2_fname)
            call s2file%read(sig2_per_ptcl)
            allocate(sig2arr(0:R), source=0.0)
            do kk = 0, R
                kclamp = min(max(kk, lbound(sig2_per_ptcl,1)), ubound(sig2_per_ptcl,1))
                sig2arr(kk) = sum(sig2_per_ptcl(kclamp, pinds)) / real(nptcls)
            end do
            l_have_sigma2 = .true.
            write(logfhandle,'(a)') '>>> RECONSTRUCT3D_PCG: objfun=euclid, using sigma2 from '//&
                &sigma2_fname%to_char()//' (averaged over the selection)'
        else
            write(logfhandle,'(a)') '>>> RECONSTRUCT3D_PCG: objfun='//trim(params%objfun)//', no sigma weighting'
        endif

        ! ---- solve ----
        allocate(x(params%box,params%box,params%box), source=0.0)
        if( l_have_sigma2 )then
            call pcgop%solve(y_planes, selection, x, use_ctf=.true., sig2arr=sig2arr, &
                &rel_res_hist=rel_res_hist, niters=niters)
        else
            call pcgop%solve(y_planes, selection, x, use_ctf=.true., &
                &rel_res_hist=rel_res_hist, niters=niters)
        endif
        rt_tot = toc(t0)

        ! ---- output: experimental volume + diagnostics table. Deliberately
        !      NOT any simple_refine3D_fnames name, and no call to
        !      spproj%write or anything else that would mutate the project. ----
        call outvol%new([params%box,params%box,params%box], params%smpd)
        call outvol%set_rmat(x, .false.)
        call outvol%write(string('reconstruct3D_pcg_state'//int2str_pad(params%state,2)//'.mrc'))
        open(newunit=funit, file='reconstruct3D_pcg_diagnostics.txt', status='replace', action='write')
        write(funit,'(a,i0)')   'nptcls=', nptcls
        write(funit,'(a,i0)')   'state=', params%state
        write(funit,'(a,i0)')   'niters=', niters
        write(funit,'(a,l1)')   'used_sigma2_file=', l_have_sigma2
        write(funit,'(a,f8.2)') 'wall_time_sec=', rt_tot
        do i = 1, niters
            write(funit,'(a,i0,a,es14.6)') 'iter', i, '_rel_resid=', rel_res_hist(i)
        end do
        close(funit)
        write(logfhandle,'(a,i0,a)')   '>>> RECONSTRUCT3D_PCG: PCG ran ', niters, ' iterations'
        write(logfhandle,'(a,f8.2,a)') '>>> RECONSTRUCT3D_PCG: wall time = ', rt_tot, ' s'

        call build%kill_general_tbox
        call simple_end('**** SIMPLE_RECONSTRUCT3D_PCG NORMAL STOP ****')
    end subroutine exec_reconstruct3D_pcg

end module simple_commanders_pcg_recon
