!@descr: single owner of the sigma2 bootstrap for work that has alignments (or
!  nothing) but no noise-power estimate yet. One rule everywhere: seed with
!  particle image power (calc_pspec) and let the first euclid pass replace it
!  with residual sigmas; where no euclid pass follows (final reconstructions),
!  run one residual pass (refine=sigma) against the seeded map.
module simple_sigma2_bootstrap
use simple_commanders_api
use simple_commanders_euclid, only: commander_calc_pspec
use simple_sigma2_files,      only: canonical_sigma2_consumable
implicit none

public :: sigma2_estimate_available, ensure_sigma2_for_iteration
public :: prepare_residual_sigma2_pass_cline
private
#include "simple_local_flags.inc"

contains

    !> Can the intended euclid consumer load a registered, file-valid,
    !! committed canonical state with the expected grid, layout, and grouping?
    logical function sigma2_estimate_available( projfile, box, smpd, l_sigma_glob ) &
            &result( l_available )
        class(string), intent(in) :: projfile
        integer,       intent(in) :: box
        real,          intent(in) :: smpd
        logical,       intent(in) :: l_sigma_glob
        type(sp_project) :: spproj
        character(len=STDLEN) :: message
        call spproj%read_segment('projinfo', projfile)
        call spproj%read_segment('ptcl3D',   projfile)
        l_available = canonical_sigma2_consumable(spproj, spproj%os_ptcl3D, box, smpd, l_sigma_glob, message)
        if( .not. l_available ) write(logfhandle,'(A)') '>>> SIGMA2 BOOTSTRAP: canonical state not consumable: '//trim(message)
        call spproj%kill
    end function sigma2_estimate_available

    !> Guarantee a canonical sigma2 estimate for the given particle project.
    !! When one is available nothing happens. Otherwise calc_pspec derives the
    !! particle power spectra and atomically publishes a committed state.
    subroutine ensure_sigma2_for_iteration( template_cline, projfile, iter, box, smpd, l_sigma_glob, &
            &label, l_bootstrapped )
        class(cmdline),           intent(in)    :: template_cline
        class(string),            intent(in)    :: projfile
        integer,                  intent(in)    :: iter, box
        real,                     intent(in)    :: smpd
        logical,                  intent(in)    :: l_sigma_glob
        character(len=*),         intent(in)    :: label
        logical,                  intent(out)   :: l_bootstrapped
        type(commander_calc_pspec) :: xcalc_pspec
        type(cmdline) :: cline_pspec
        integer       :: state
        l_bootstrapped = .false.
        if( sigma2_estimate_available(projfile, box, smpd, l_sigma_glob) ) return
        cline_pspec = template_cline
        call cline_pspec%set('prg',                    'calc_pspec')
        call cline_pspec%set('mkdir',                          'no')
        call cline_pspec%set('projfile',                   projfile)
        call cline_pspec%set('objfun',                     'euclid')
        call cline_pspec%set('sigma_est',                  'global')
        call cline_pspec%set('cc_emit_sigma',                  'no')
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
            &': no compatible canonical sigma2 state; seeding from particle power spectra at iteration ', max(1, iter)
        call xcalc_pspec%execute(cline_pspec)
        call cline_pspec%kill
        l_bootstrapped = .true.
    end subroutine ensure_sigma2_for_iteration

    !> refine3D command line for one residual sigma2 pass: no search, every
    !! particle's sigma2 re-estimated from its residual against the given state
    !! volumes at the template's sampling; no volume assembly, no orientation
    !! output. The pass commits its candidate as the next canonical generation.
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

end module simple_sigma2_bootstrap
