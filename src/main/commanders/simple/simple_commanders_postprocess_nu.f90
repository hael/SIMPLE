!@descr: NU-evidence nonuniform postprocessing (isolated from the standard postprocess path)
!
! postprocess_nu is the commander for the NU-evidence local sharpening
! experiment (nu_evidence_local_sharpening.md): model-free LocScale-style
! local amplitude restoration in which both the confidence field and the
! target spectrum derive from the frozen cross-half NU evidence state (the
! compact state the NU competition also uses for its envelope). It is deliberately
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
use simple_nu_filter, only: setup_nu_dmats, optimize_nu_cutoff_finds, &
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
        real, allocatable       :: corrs(:), res(:)
        real                    :: fsc05, fsc0143
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir',     'yes')
        call params%new(cline)
        call odd%new([params%box,params%box,params%box], params%smpd)
        call even%new([params%box,params%box,params%box], params%smpd)
        call odd%read(params%vols(1))
        call even%read(params%vols(2))
        ! the same bank as the refinement (2026-09-16, the shell walk
        ! retired): coarse ladder plus fine rungs generated from the box,
        ! bounded by the pair's own FSC=0.143 / NU_BANK_FSC_HEADROOM. No
        ! regularized member: the evidence compaction takes the base pair's
        ! candidates only (its contract), and this program has no
        ! regularized pair input
        allocate(corrs(fdim(params%box)-1), source=0.)
        ! image%fsc reads Fourier coefficients: transform, correlate, and
        ! return the halves to real space for the evidence setup
        call even%fft()
        call odd%fft()
        call even%fsc(odd, corrs)
        call even%ifft()
        call odd%ifft()
        res = get_resarr(params%box, params%smpd)
        call get_resolution(corrs, res, fsc05, fsc0143)
        write(logfhandle,'(A,F8.3,A,F8.3,A)') '>>> POSTPROCESS_NU: HALF-MAP FSC=0.5 ', fsc05, ' A, FSC=0.143 ', fsc0143, ' A'
        ! identical halves (vol1 = vol2, or two copies of one map) give an FSC
        ! of 1 on every shell: no 0.143 crossing, an unbounded bank to Nyquist
        ! and degenerate cross-half evidence that awards the finest cutoff
        ! everywhere; the sharpened product is then noise (2026-09-17)
        if( params%vols(1) == params%vols(2) ) &
            &THROW_HARD('vol1 and vol2 are the same file; postprocess_nu needs the two independent half maps')
        if( fsc0143 <= TINY .or. fsc0143 <= 2.*params%smpd + TINY ) &
            &THROW_HARD('the half-map FSC never falls below 0.143: are vol1/vol2 independent half maps? postprocess_nu')
        ! frozen evidence from the unregularized half pair, the standard
        ! lifecycle: bank -> compact immutable state (one evidence identity,
        ! no second NU analysis)
        call setup_nu_dmats(even, odd, params%mskdiam, [real ::], evidence_source=NU_EVIDENCE_SOURCE_BASE, &
            &fsc_res=fsc0143)
        call optimize_nu_cutoff_finds()
        write(logfhandle,'(A,F8.3,A)') '>>> POSTPROCESS_NU: EVIDENCE BANK FINEST MEMBER ', get_nu_filter_bank_finest_lp(), ' A'
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
        if( allocated(corrs) ) deallocate(corrs)
        if( allocated(res)   ) deallocate(res)
        call simple_end('**** SIMPLE_POSTPROCESS_NU NORMAL STOP ****')
    end subroutine exec_postprocess_nu

end module simple_commanders_postprocess_nu
