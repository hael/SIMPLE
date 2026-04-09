!@descr: polar-reference preparation helpers for matcher workflows
module simple_matcher_refpolar_utils
use simple_pftc_srch_api
use simple_builder,              only: builder
use simple_matcher_refvol_utils, only: report_resolution, estimate_lp_refvols3D
implicit none

public :: prep_pftc_polar_mode
private
#include "simple_local_flags.inc"

contains

    subroutine prep_pftc_polar_mode( params, build, cline, batchsz )
        class(parameters),        intent(inout) :: params
        class(builder),           intent(inout) :: build
        class(cmdline),           intent(in)    :: cline
        integer,                  intent(in)    :: batchsz
        real, allocatable :: gaufilter(:)
        integer           :: iproj, nrefs, filtsz, state
        logical           :: l_filtrefs
        ! Resolution limit estimation
        call report_resolution(params, build, state)
        if( cline%defined('lpstart') .and. cline%defined('lpstop') )then
            call estimate_lp_refvols3D(params, build, cline, params%lpstart, params%lpstop, state)
        endif
        ! Calculator init
        nrefs = params%nspace * params%nstates
        call build%pftc%new(params, nrefs, [1,batchsz], params%kfromto)
        ! Read polar references
        call build%pftc%polar_cavger_new(.true.)
        call build%pftc%polar_cavger_read_all(string(POLAR_REFS_FBODY//BIN_EXT))
        call build%clsfrcs%read(string(FRCS_FILE))
        ! prepare filter
        l_filtrefs = .false.
        if(trim(params%gauref).eq.'yes')then
            l_filtrefs = .true.
            filtsz = build%clsfrcs%get_filtsz()
            allocate(gaufilter(filtsz),source=0.)
            call gaussian_filter(params%gaufreq, params%smpd, params%box, gaufilter)
        endif
        ! PREPARATION OF REFERENCES IN pftc
        !$omp parallel do default(shared) private(iproj)&
        !$omp schedule(static) proc_bind(close)
        do iproj = 1,params%nspace
            if( l_filtrefs ) call build%pftc%polar_filterrefs(iproj, gaufilter)
            ! transfer to pftc
            if( params%l_lpset )then
                ! put the merged class average in both even and odd positions
                call build%pftc%polar_cavger_set_ref_pft(iproj, 'merged')
                call build%pftc%cp_even2odd_ref(iproj)
            else
                ! transfer e/o refs to pftc
                call build%pftc%polar_cavger_set_ref_pft(iproj, 'even')
                call build%pftc%polar_cavger_set_ref_pft(iproj, 'odd')
            endif
        end do
        !$omp end parallel do
        if( allocated(gaufilter) ) deallocate(gaufilter)
    end subroutine prep_pftc_polar_mode

end module simple_matcher_refpolar_utils