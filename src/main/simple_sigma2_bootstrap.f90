!@descr: single owner of the sigma2 bootstrap for work that has alignments (or
!  nothing) but no noise-power estimate yet. One rule everywhere: seed with
!  particle image power (calc_pspec) and let the first euclid pass replace it
!  with residual sigmas; where no euclid pass follows (final reconstructions),
!  run one residual pass (refine=sigma) against the seeded map. Sigma files are
!  never discovered outside the current directory (2026-09-06).
module simple_sigma2_bootstrap
use simple_commanders_api
use simple_commanders_euclid, only: commander_calc_pspec, commander_calc_group_sigmas
use simple_sigma2_files,      only: pick_sigma_group_file_for_iter, canonical_sigma2_consumable
implicit none

public :: sigma2_estimate_available, ensure_sigma2_for_iteration
public :: prepare_residual_sigma2_pass_cline, consolidate_sigma2_groups
private
#include "simple_local_flags.inc"

contains

    !> Can the intended euclid consumer load a sigma2 estimate from the
    !! current directory: the consumer's own selection rule for the legacy
    !! store (the STAR named for the iteration, else the latest STAR), and the
    !! consumer's own identity validation for the canonical store (registered,
    !! file-valid, committed, native grid, ordered layout, grouping)
    logical function sigma2_estimate_available( l_canonical, projfile, iter, box, smpd, l_sigma_glob ) &
            &result( l_available )
        logical,       intent(in) :: l_canonical
        class(string), intent(in) :: projfile
        integer,       intent(in) :: iter, box
        real,          intent(in) :: smpd
        logical,       intent(in) :: l_sigma_glob
        type(sp_project) :: spproj
        type(string)     :: fname
        character(len=STDLEN) :: message
        l_available = .false.
        if( l_canonical )then
            call spproj%read_segment('projinfo', projfile)
            call spproj%read_segment('ptcl3D',   projfile)
            l_available = canonical_sigma2_consumable(spproj, spproj%os_ptcl3D, box, smpd, l_sigma_glob, message)
            if( .not. l_available ) write(logfhandle,'(A)') '>>> SIGMA2 BOOTSTRAP: canonical state not consumable: '//trim(message)
            call spproj%kill
        else
            call pick_sigma_group_file_for_iter(iter, fname, l_available)
            call fname%kill
        endif
    end function sigma2_estimate_available

    !> Guarantee a sigma2 estimate for the given iteration in the current
    !! directory. When one is available nothing happens. Otherwise the
    !! particle power spectra are estimated (calc_pspec) as the grouped STAR
    !! of that iteration plus per-particle files in the template's partition
    !! layout, and, when a consumer command line is given, the handover key
    !! sigma_transition_ready=yes is set on it so the consuming refine3D's
    !! first iteration initializes its workers from that STAR.
    subroutine ensure_sigma2_for_iteration( template_cline, projfile, iter, l_canonical, box, smpd, l_sigma_glob, &
            &label, l_bootstrapped, consumer_cline )
        class(cmdline),           intent(in)    :: template_cline
        class(string),            intent(in)    :: projfile
        integer,                  intent(in)    :: iter, box
        logical,                  intent(in)    :: l_canonical
        real,                     intent(in)    :: smpd
        logical,                  intent(in)    :: l_sigma_glob
        character(len=*),         intent(in)    :: label
        logical,                  intent(out)   :: l_bootstrapped
        class(cmdline), optional, intent(inout) :: consumer_cline
        type(commander_calc_pspec) :: xcalc_pspec
        type(cmdline) :: cline_pspec
        integer       :: state
        l_bootstrapped = .false.
        if( sigma2_estimate_available(l_canonical, projfile, iter, box, smpd, l_sigma_glob) ) return
        cline_pspec = template_cline
        call cline_pspec%set('prg',                    'calc_pspec')
        call cline_pspec%set('mkdir',                          'no')
        call cline_pspec%set('projfile',                   projfile)
        call cline_pspec%set('objfun',                     'euclid')
        call cline_pspec%set('sigma_est',                  'global')
        call cline_pspec%set('cc_emit_sigma',                  'no')
        call cline_pspec%set('sigma_transition_ready',         'no')
        call cline_pspec%set('which_iter',           max(1, iter))
        call cline_pspec%delete('part')
        call cline_pspec%delete('update_frac')
        call cline_pspec%delete('nsample')
        call cline_pspec%delete('fillin')
        call cline_pspec%delete('endit')
        call cline_pspec%delete('startit')
        call cline_pspec%delete('ml_reg')
        call cline_pspec%delete('postprocess')
        call cline_pspec%delete('combine_eo')
        call cline_pspec%delete('rec_backend')
        call cline_pspec%delete('maxits_pcg')
        call cline_pspec%delete('rtol')
        call cline_pspec%delete('trail_seed')
        call cline_pspec%delete('trail_rec')
        call cline_pspec%delete('outfile')
        do state = 1, 99
            if( .not. cline_pspec%defined('vol'//int2str(state)) ) exit
            call cline_pspec%delete('vol'//int2str(state))
        enddo
        write(logfhandle,'(A,I0)') '>>> '//trim(label)// &
            &': no sigma2 estimate in this directory; seeding from particle power spectra at iteration ', max(1, iter)
        call xcalc_pspec%execute(cline_pspec)
        call cline_pspec%kill
        l_bootstrapped = .true.
        if( present(consumer_cline) ) call consumer_cline%set('sigma_transition_ready', 'yes')
    end subroutine ensure_sigma2_for_iteration

    !> refine3D command line for one residual sigma2 pass: no search, every
    !! particle's sigma2 re-estimated from its residual against the given state
    !! volumes at the template's sampling; no volume assembly, no orientation
    !! output. The per-particle files it writes are consolidated by
    !! consolidate_sigma2_groups(iter+1).
    subroutine prepare_residual_sigma2_pass_cline( template_cline, iter, nstates, vols, cline_sigma )
        class(cmdline), intent(in)    :: template_cline
        integer,        intent(in)    :: iter, nstates
        class(string),  intent(in)    :: vols(:)
        type(cmdline),  intent(inout) :: cline_sigma
        integer :: state
        if( size(vols) < nstates ) THROW_HARD('residual sigma2 pass needs one volume per state')
        cline_sigma = template_cline
        call cline_sigma%set('prg',                      'refine3D')
        call cline_sigma%set('mkdir',                          'no')
        call cline_sigma%set('refine',                      'sigma')
        call cline_sigma%set('objfun',                     'euclid')
        call cline_sigma%set('ml_reg',                         'no')
        call cline_sigma%set('sigma_est',                  'global')
        call cline_sigma%set('sigma_transition_ready',         'no')
        call cline_sigma%set('volrec',                         'no')
        call cline_sigma%set('continue',                       'no')
        call cline_sigma%set('center',                         'no')
        call cline_sigma%set('maxits',                            1)
        call cline_sigma%set('minits',                            1)
        call cline_sigma%set('startit',              max(1, iter))
        call cline_sigma%set('which_iter',           max(1, iter))
        call cline_sigma%set('extr_iter',            max(1, iter))
        call cline_sigma%set('nstates',                   nstates)
        call cline_sigma%delete('update_frac')
        call cline_sigma%delete('nsample')
        call cline_sigma%delete('fillin')
        call cline_sigma%delete('endit')
        call cline_sigma%delete('trail_rec')
        call cline_sigma%delete('trail_seed')
        call cline_sigma%delete('ufrac_trec')
        call cline_sigma%delete('objfun_den')
        call cline_sigma%delete('objfun_den_w')
        call cline_sigma%delete('sticky_class_sampling')
        call cline_sigma%delete('postprocess')
        call cline_sigma%delete('combine_eo')
        call cline_sigma%delete('outfile')
        call cline_sigma%delete('vol_even')
        call cline_sigma%delete('vol_odd')
        do state = 1, nstates
            call cline_sigma%set('vol'//int2str(state), vols(state))
        enddo
    end subroutine prepare_residual_sigma2_pass_cline

    !> Consolidate the per-particle sigma2 files a residual pass left in the
    !! current directory into the grouped STAR of the given iteration. Legacy
    !! store only: a canonical refine3D pass owns its complete transaction
    !! (prepare, merge range files, reduce, commit, delete ranges) and returns
    !! with a committed state, so a second consolidation would find nothing
    !! to merge.
    subroutine consolidate_sigma2_groups( template_cline, projfile, iter, l_canonical )
        class(cmdline), intent(in) :: template_cline
        class(string),  intent(in) :: projfile
        integer,        intent(in) :: iter
        logical,        intent(in) :: l_canonical
        type(commander_calc_group_sigmas) :: xcalc_group_sigmas
        type(cmdline) :: cline_groups
        if( l_canonical )then
            write(logfhandle,'(A)') '>>> SIGMA2 BOOTSTRAP: residual groups committed by the canonical refine3D pass'
            return
        endif
        cline_groups = template_cline
        call cline_groups%set('prg',        'calc_group_sigmas')
        call cline_groups%set('mkdir',                     'no')
        call cline_groups%set('projfile',              projfile)
        call cline_groups%set('which_iter',        max(1, iter))
        call cline_groups%delete('part')
        call cline_groups%delete('outfile')
        call xcalc_group_sigmas%execute(cline_groups)
        call cline_groups%kill
    end subroutine consolidate_sigma2_groups

end module simple_sigma2_bootstrap
