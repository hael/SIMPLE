!@descr: NU-evidence nonuniform postprocessing (isolated from the standard postprocess path)
!
! postprocess_nu is the commander for the NU-evidence local sharpening
! experiment (nu_evidence_local_sharpening.md): model-free LocScale-style
! local amplitude restoration in which both the confidence field and the
! target spectrum derive from the frozen cross-half NU evidence state -- the
! same state that parameterizes the Q_NU replay prior. It is deliberately
! isolated from the standard postprocess commander (global B-factor + FSC
! filter), which remains untouched: a single isotropic B-factor does not
! serve most specimens, and this path is the recorded alternative.
!
! Inputs are the UNREGULARIZED even/odd half maps (evidence authority is the
! base pair; regularized maps flatten the evidence margin). Products carry
! the _nu_sharp suffix and are display/interpretation maps only -- they must
! never feed FSC correction or resolution claims.
module simple_commanders_postprocess_nu
use simple_commanders_api
use simple_nu_filter, only: setup_nu_dmats, optimize_nu_cutoff_finds, extend_nu_filter_highres_shells, &
    &get_nu_filter_bank_finest_lp, build_nu_evidence_state, cleanup_nu_filter, nu_evidence_state, &
    &assert_nu_evidence_replay_ready, print_nu_evidence_summary, nu_evidence_sharpen_vol, NU_EVIDENCE_SOURCE_BASE
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_postprocess_nu
  contains
    procedure :: execute      => exec_postprocess_nu
end type commander_postprocess_nu

contains

    subroutine exec_postprocess_nu( self, cline )
        class(commander_postprocess_nu), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters)        :: params
        type(nu_evidence_state) :: evstate
        type(image)             :: even, odd, vol_sharp
        type(string)            :: vol_out
        integer                 :: nsteps_ext
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir',     'yes')
        ! finer evidence granularity is affordable post-hoc (cost paid once
        ! per map, not once per CG iteration), so the evidence-gated shell
        ! walk defaults ON here, unlike in abinitio3D
        if( .not. cline%defined('nu_refine') ) call cline%set('nu_refine', 'yes')
        call params%new(cline)
        call odd%new([params%box,params%box,params%box], params%smpd)
        call even%new([params%box,params%box,params%box], params%smpd)
        call odd%read(params%vols(1))
        call even%read(params%vols(2))
        ! frozen evidence from the unregularized half pair, exactly the Q_NU
        ! replay lifecycle: bank -> optional accepted shell walk -> compact
        ! immutable state (one evidence identity, no second NU analysis)
        call setup_nu_dmats(even, odd, params%mskdiam, [real ::], evidence_source=NU_EVIDENCE_SOURCE_BASE)
        call optimize_nu_cutoff_finds()
        if( params%l_nu_refine )then
            call extend_nu_filter_highres_shells(even, odd, nsteps=nsteps_ext)
            if( nsteps_ext > 0 )then
                write(logfhandle,'(A,I0,A,F8.3,A)') '>>> POSTPROCESS_NU: EVIDENCE BANK EXTENDED BY ', &
                    &nsteps_ext, ' ACCEPTED SHELL STEP(S) TO ', get_nu_filter_bank_finest_lp(), ' A'
            endif
        endif
        call build_nu_evidence_state(even, odd, evstate)
        call cleanup_nu_filter()
        call assert_nu_evidence_replay_ready(evstate)
        call print_nu_evidence_summary(evstate)
        ! classical shrink-then-sharpen localized by the evidence; the shipped
        ! product is the single sharpened merged volume
        call nu_evidence_sharpen_vol(evstate, even, odd, vol_sharp)
        if( params%outvol .ne. '' )then
            vol_out = params%outvol
        else
            vol_out = 'vol'//NUSHARP_SUFFIX//params%ext%to_char()
        endif
        call vol_sharp%write(vol_out, del_if_exists=.true.)
        call wait_for_closure(vol_out)
        ! destruct
        call even%kill
        call odd%kill
        call vol_sharp%kill
        call vol_out%kill
        call simple_end('**** SIMPLE_POSTPROCESS_NU NORMAL STOP ****')
    end subroutine exec_postprocess_nu

end module simple_commanders_postprocess_nu
