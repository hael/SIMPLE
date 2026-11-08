!@descr: offline order prior of an even/odd pair (program order_prior)
! Outside the refinement loop: the pcg_order closed-form choice on the prior-free pair, checked against full
! surrogate solves of the PCG re-solve; the pair's NU label field is composed from the prior'd pair, the
! references the matcher would see. Design: doc/implementation_notes/planned/pcg_order_prior.md.
module simple_commanders_order_prior
use simple_defs,           only: logfhandle
use simple_defs_fname,     only: NUFILT_SUFFIX, NUEVIDENCE_SUFFIX, NULOCRES_SUFFIX
use simple_commander_base, only: commander_base
use simple_cmdline,        only: cmdline
use simple_parameters,     only: parameters
use simple_image,          only: image
use simple_string,         only: string
use simple_fileio,         only: add2fbody, wait_for_closure
use simple_jiffys,         only: simple_end
implicit none
private

public :: commander_order_prior

type, extends(commander_base) :: commander_order_prior
  contains
    procedure :: execute => exec_order_prior
end type commander_order_prior

contains

    subroutine exec_order_prior( self, cline )
        use simple_order_prior_offline, only: order_prior_offline
        use simple_nu_filter,           only: setup_nu_dmats, optimize_nu_cutoff_finds, nu_filter_vols, cleanup_nu_filter, &
            &get_nu_evidence_map, get_nu_local_resolution_map, print_nu_filtmap_lowpass_stats, NU_ENVMASK_RELATIVE
        class(commander_order_prior), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters)          :: params
        type(order_prior_offline) :: oprior
        type(image)               :: even, odd, evmap, locres, w_dens(2), w_order(2), order_even, order_odd
        type(string)              :: outvol
        real                      :: lambda_fixed
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        if( cline%defined('outvol') )then
            outvol = params%outvol
        else
            outvol = 'order_prior'//params%ext%to_char()
        endif
        lambda_fixed = -1.
        if( .not. params%l_pcg_order_lambda_auto ) lambda_fixed = params%pcg_order_lambda
        ! vol1 odd, vol2 even, as nu_filt3D; the prior-free (_unfil) pair is the intended input
        call odd%new([params%box,params%box,params%box],  params%smpd)
        call even%new([params%box,params%box,params%box], params%smpd)
        call odd%read(params%vols(1))
        call even%read(params%vols(2))
        call oprior%new(even, odd, params%mskdiam)
        call oprior%density_weights(w_dens)
        ! NU on the prior-free pair, capped as in assembly: the evidence margin and the label field; it stays
        ! live until the prior'd pair is composed
        call setup_nu_dmats(even, odd, params%mskdiam, [real ::], fsc_res=oprior%get_res0143())
        call optimize_nu_cutoff_finds()
        call print_nu_filtmap_lowpass_stats()
        call get_nu_evidence_map(evmap, params%amsklp, NU_ENVMASK_RELATIVE)
        call get_nu_local_resolution_map(locres)
        call oprior%run(evmap, locres, w_dens, lambda_fixed, w_order, order_even, order_odd)
        call write_vol(evmap,      NUEVIDENCE_SUFFIX)
        call write_vol(locres,     NULOCRES_SUFFIX)
        call write_vol(w_dens(1),  '_density_weight_even')
        call write_vol(w_dens(2),  '_density_weight_odd')
        call write_vol(w_order(1), '_order_weight')
        call write_pair(even,       odd,       '_base')
        call write_pair(order_even, order_odd, '_order')
        call write_nu_refs(even,       odd,       '_base')
        call write_nu_refs(order_even, order_odd, '_order')
        call cleanup_nu_filter()
        call oprior%kill
        call even%kill
        call odd%kill
        call evmap%kill
        call locres%kill
        call w_dens(1)%kill
        call w_dens(2)%kill
        call w_order(1)%kill
        call w_order(2)%kill
        call order_even%kill
        call order_odd%kill
        call outvol%kill
        call simple_end('**** SIMPLE_ORDER_PRIOR NORMAL STOP ****')

    contains

        subroutine write_vol( vol, suffix )
            class(image),     intent(inout) :: vol
            character(len=*), intent(in)    :: suffix
            type(string) :: fname
            fname = add2fbody(outvol, params%ext, suffix)
            call vol%write(fname, del_if_exists=.true.)
            call wait_for_closure(fname)
            write(logfhandle,'(A)') '>>> ORDER PRIOR: WROTE '//fname%to_char()
            call fname%kill
        end subroutine write_vol

        !> the halves (not of the input pair) and their average
        subroutine write_pair( vol_even, vol_odd, suffix )
            class(image),     intent(inout) :: vol_even, vol_odd
            character(len=*), intent(in)    :: suffix
            type(image) :: avg
            if( suffix /= '_base' )then
                call write_vol(vol_even, suffix//'_even')
                call write_vol(vol_odd,  suffix//'_odd')
            endif
            call avg%copy(vol_even)
            call avg%add(vol_odd)
            call avg%mul(0.5)
            call write_vol(avg, suffix)
            call avg%kill
        end subroutine write_pair

        !> the label field of the prior-free pair composed from this pair, averaged
        subroutine write_nu_refs( vol_even, vol_odd, suffix )
            class(image),     intent(in) :: vol_even, vol_odd
            character(len=*), intent(in) :: suffix
            type(image) :: even_nu, odd_nu
            call nu_filter_vols(even_nu, odd_nu, vol_even, vol_odd)
            call even_nu%add(odd_nu)
            call even_nu%mul(0.5)
            call write_vol(even_nu, suffix//NUFILT_SUFFIX)
            call even_nu%kill
            call odd_nu%kill
        end subroutine write_nu_refs

    end subroutine exec_order_prior

end module simple_commanders_order_prior
