!@descr: assembly-owned nonuniform (NU) filtering of one state's half-map pair
!
!  One routine runs the NU competition for a state exactly as the gridding
!  volassemble always has: a discrete low-pass candidate bank built from the
!  BASE (unregularized) even/odd pair, the ML-regularized pair joining the
!  competition as the finest auxiliary member (ml_reg=yes, nu_refine=no), the
!  high-resolution shell walk extending the bank (nu_refine=yes), the
!  NU-evidence envelope fixing the filter-field background (automsk=yes),
!  synthesis of the filtered even/odd/merged references, the local-resolution
!  map, and the raw finest selected label as the matching low-pass handoff.
!  Both reconstruction backends call it (policy 2026-09-06): the PCG path
!  mirrors gridding and carries no prior of its own beyond the P_tau replay.
module simple_nu_state_filter
use simple_core_module_api
use simple_image,            only: image
use simple_parameters,       only: parameters
use simple_nu_filter,        only: setup_nu_dmats, optimize_nu_cutoff_finds, nu_filter_vols, &
    &cleanup_nu_filter, print_nu_filtmap_lowpass_stats, analyze_filtmap_neighbor_continuity, &
    &NU_DEV_OUTPUT, extend_nu_filter_highres_shell_next, refine_nu_extension_filtmap_ordered_labels, &
    &nu_highres_extension_stats, get_nu_filtmap_finest_selected_lp, get_nu_bank_cap_find, &
    &get_nu_filtmap_highres_shell_depth, write_nu_local_resolution_map, write_nu_evidence_envmask
use simple_vol_pproc_policy, only: vol_pproc_plan, plan_state_postprocess
implicit none

public :: nonuniform_filter_state, nu_state_filter_timings, nu_static_aux_replacement
private
#include "simple_local_flags.inc"

type nu_state_filter_timings
    real(timer_int_kind) :: envmask = 0.
    real(timer_int_kind) :: filter  = 0.
end type nu_state_filter_timings

contains

    !> The ML-regularized pair joins the competition as the auxiliary member
    !! only in static-bank mode; with nu_refine=yes the shell walk owns the
    !! resolution-extension experiment.
    pure logical function nu_static_aux_replacement( params ) result( l_use_aux )
        class(parameters), intent(in) :: params
        l_use_aux = params%l_ml_reg .and. .not. params%l_nu_refine
    end function nu_static_aux_replacement

    !> Run the NU competition for one state and write its derived products.
    !! vol_base_even/odd: the unregularized pair (consumed and killed here).
    !! vol_aux_even/odd:  the ML-regularized pair when l_use_aux (consumed and
    !!                    killed here), ignored otherwise.
    !! res0143:           FSC=0.143 crossing of the base pair, the auxiliary
    !!                    member's effective resolution (clamped by a set lp).
    !! volname/eonames:   the state's merged and even/odd file names; the
    !!                    _nu_filt and _nu_locres products derive from them.
    !! align_lp:          raw finest selected label (0 when none), the
    !!                    matching low-pass handoff for the next iteration.
    subroutine nonuniform_filter_state( params, state, which_iter, vol_base_even, vol_base_odd, &
            &vol_aux_even, vol_aux_odd, l_use_aux, res0143, volname, eonames, align_lp, timings )
        class(parameters),            intent(in)    :: params
        integer,                      intent(in)    :: state, which_iter
        type(image),                  intent(inout) :: vol_base_even, vol_base_odd
        type(image),                  intent(inout) :: vol_aux_even, vol_aux_odd
        logical,                      intent(in)    :: l_use_aux
        real,                         intent(in)    :: res0143
        class(string),                intent(in)    :: volname, eonames(2)
        real,                         intent(out)   :: align_lp
        type(nu_state_filter_timings), optional, intent(inout) :: timings
        type(image), allocatable :: nu_aux_even(:), nu_aux_odd(:)
        type(image)              :: vol_even_nu, vol_odd_nu
        type(vol_pproc_plan)     :: pp_plan
        integer(timer_int_kind)  :: t_filter, t_envmask
        integer :: n_highres_steps
        real    :: aux_resolution
        align_lp = 0.
        if( L_BENCH_GLOB ) t_filter = tic()
        call plan_state_postprocess(params, state, which_iter, pp_plan)
        if( pp_plan%l_nu_envmask_incompatible )then
            write(logfhandle,'(A,1X,A)') &
                &'>>> Existing NU evidence envelope incompatible with current box/sampling, regenerating:', &
                &pp_plan%nu_envmask_file%to_char()
        endif
        ! candidate bank from the base pair, auxiliary member from the ML pair
        n_highres_steps = nu_highres_steps_for_state()
        if( l_use_aux )then
            allocate(nu_aux_even(1), nu_aux_odd(1))
            call nu_aux_even(1)%copy(vol_aux_even)
            call nu_aux_odd(1)%copy(vol_aux_odd)
            aux_resolution = nu_aux_effective_resolution()
            call setup_nu_dmats(vol_base_even, vol_base_odd, params%mskdiam, [aux_resolution], &
                &nu_aux_even, nu_aux_odd, n_highres_steps=n_highres_steps, fsc_res=res0143)
        else
            call setup_nu_dmats(vol_base_even, vol_base_odd, params%mskdiam, [real ::], &
                &n_highres_steps=n_highres_steps, fsc_res=res0143)
        endif
        if( trim(params%automsk).ne.'no' )then
            ! automsk=yes: the filter-field background is the complement of
            ! the NU evidence envelope, derived from the unaries of the setup
            ! that just ran (same pass, no second compute). The objective
            ! domain remains the spherical mskdiam support; nu_refine
            ! independently controls high-resolution extension.
            if( L_BENCH_GLOB ) t_envmask = tic()
            call write_nu_evidence_envmask(params%nu_msk_sig, params%amsklp, &
                &vol_base_even%get_smpd(), state, pp_plan%nu_envmask_file, l_arm_background=.true.)
            if( L_BENCH_GLOB .and. present(timings) ) timings%envmask = timings%envmask + toc(t_envmask)
        endif
        ! the auxiliary inputs are copied into the bank; release them
        call cleanup_nu_aux_images()
        call vol_aux_even%kill
        call vol_aux_odd%kill
        call optimize_nu_cutoff_finds()
        call refine_nonuniform_filter_bank()
        call vol_base_even%kill
        call vol_base_odd%kill
        call nu_filter_vols(vol_even_nu, vol_odd_nu)
        call print_nu_filtmap_lowpass_stats()
        if( NU_DEV_OUTPUT .and. params%part == 1 ) call analyze_filtmap_neighbor_continuity()
        call write_nonuniform_outputs()
        call record_nu_alignment_lowpass_limit()
        call vol_even_nu%kill
        call vol_odd_nu%kill
        call cleanup_nu_filter()
        call pp_plan%nu_envmask_file%kill
        if( L_BENCH_GLOB .and. present(timings) ) timings%filter = timings%filter + toc(t_filter)

    contains

        real function nu_aux_effective_resolution() result(aux_res)
            aux_res = res0143
            if( params%l_lpset .and. params%lp > TINY )then
                if( NU_DEV_OUTPUT .and. params%part == 1 .and. aux_res > params%lp + TINY )then
                    write(logfhandle,'(A,F8.3,A,F8.3,A)') &
                        &'>>> NU auxiliary effective resolution clamped by matching low-pass: FSC ', &
                        &aux_res, ' A; matching LP ', params%lp, ' A'
                endif
                aux_res = min(aux_res, params%lp)
            endif
        end function nu_aux_effective_resolution

        subroutine refine_nonuniform_filter_bank()
            type(nu_highres_extension_stats) :: ext_stats
            integer :: nsteps, n_accepted_this_iteration, cap_find
            if( .not. params%l_nu_refine ) return
            n_accepted_this_iteration = 0
            ! the shell walk cannot pass the FSC-anchored candidate cap
            cap_find = get_nu_bank_cap_find()
            do
                if( cap_find > 0 )then
                    call extend_nu_filter_highres_shell_next(vol_base_even, vol_base_odd, stats=ext_stats, &
                        &max_find=cap_find)
                else
                    call extend_nu_filter_highres_shell_next(vol_base_even, vol_base_odd, stats=ext_stats)
                endif
                if( .not. ext_stats%attempted )then
                    if( NU_DEV_OUTPUT .and. params%part == 1 )then
                        if( ext_stats%n_mask == 0 )then
                            write(logfhandle,'(A)') &
                                &'>>> NU high-resolution extension stopped: empty NU refinement mask'
                        else if( ext_stats%n_tested == 0 )then
                            write(logfhandle,'(A,F8.3,A,I0,A)') &
                                &'>>> NU high-resolution extension stopped: no frontier voxels at current finest label ', &
                                &ext_stats%old_limit, ' A (k=', ext_stats%old_find, ')'
                        else
                            write(logfhandle,'(A,F8.3,A,I0,A)') &
                                &'>>> NU high-resolution extension stopped: no valid next shell after ', &
                                &ext_stats%old_limit, ' A (k=', ext_stats%old_find, ')'
                        endif
                    endif
                    exit
                endif
                if( .not. ext_stats%applied      ) exit
                if( .not. ext_stats%promote_next ) exit
                n_accepted_this_iteration = n_accepted_this_iteration + 1
            end do
            if( n_accepted_this_iteration > 0 )then
                call refine_nu_extension_filtmap_ordered_labels
                nsteps = get_nu_filtmap_highres_shell_depth()
                call write_nu_highres_steps_for_state(nsteps)
                if( NU_DEV_OUTPUT .and. params%part == 1 )then
                    write(logfhandle,'(A,I0,A,I0)') &
                        &'>>> NU high-resolution extension accepted shell steps this iteration: ', &
                        &n_accepted_this_iteration, '; promoted depth for next iteration: ', nsteps
                endif
            endif
        end subroutine refine_nonuniform_filter_bank

        integer function nu_highres_steps_for_state() result(nsteps)
            type(string) :: fname
            integer :: funit, io_stat
            nsteps = 0
            if( .not. params%l_nu_refine ) return
            if( params%startit <= 1 .and. params%which_iter <= params%startit )then
                call write_nu_highres_steps_for_state(0)
                return
            endif
            fname = nu_highres_steps_fname()
            if( .not.file_exists(fname) )then
                call fname%kill
                return
            endif
            call fopen(funit, status='OLD', action='READ', file=fname, iostat=io_stat)
            if( io_stat == 0 )then
                read(funit, *, iostat=io_stat) nsteps
                call fclose(funit)
            endif
            if( io_stat /= 0 ) nsteps = 0
            nsteps = max(0, nsteps)
            call fname%kill
        end function nu_highres_steps_for_state

        subroutine write_nu_highres_steps_for_state( nsteps )
            integer, intent(in) :: nsteps
            type(string) :: fname
            integer :: funit, io_stat
            if( .not. params%l_nu_refine ) return
            fname = nu_highres_steps_fname()
            call fopen(funit, status='REPLACE', action='WRITE', file=fname, iostat=io_stat)
            if( io_stat == 0 )then
                write(funit,'(I0)') max(0, nsteps)
                call fclose(funit)
            else
                write(logfhandle,'(A,1X,A)') '>>> WARNING: failed to write NU high-resolution depth file:', &
                    &fname%to_char()
            endif
            call fname%kill
        end subroutine write_nu_highres_steps_for_state

        function nu_highres_steps_fname() result( fname )
            type(string) :: fname
            fname = 'nu_highres_depth_state'//int2str_pad(state,2)//'.txt'
        end function nu_highres_steps_fname

        subroutine write_nonuniform_outputs()
            type(string) :: eonames_nu(2), volname_nu, locres_name
            eonames_nu(1) = add2fbody(eonames(1), MRC_EXT, NUFILT_SUFFIX)
            eonames_nu(2) = add2fbody(eonames(2), MRC_EXT, NUFILT_SUFFIX)
            volname_nu    = add2fbody(volname,    MRC_EXT, NUFILT_SUFFIX)
            locres_name   = add2fbody(volname,    MRC_EXT, NULOCRES_SUFFIX)
            call vol_even_nu%write(eonames_nu(1), del_if_exists=.true.)
            call vol_odd_nu%write(eonames_nu(2), del_if_exists=.true.)
            call vol_even_nu%add(vol_odd_nu)
            call vol_even_nu%mul(0.5)
            call vol_even_nu%write(volname_nu, del_if_exists=.true.)
            call write_nu_local_resolution_map(locres_name)
            call wait_for_closure(volname_nu)
            call wait_for_closure(locres_name)
            call eonames_nu(1)%kill
            call eonames_nu(2)%kill
            call volname_nu%kill
            call locres_name%kill
        end subroutine write_nonuniform_outputs

        subroutine record_nu_alignment_lowpass_limit()
            real :: selected_lp
            ! raw finest selected label (min_assigned_pct=0): the 5% support
            ! gate introduced 2026-08-30 capped the PfCRT matching band at
            ! 5-6 A against a 4.1 A map and refine3D_auto degraded from there
            selected_lp = get_nu_filtmap_finest_selected_lp(min_assigned_pct=0.)
            if( selected_lp <= TINY ) return
            align_lp = selected_lp
            if( NU_DEV_OUTPUT .and. params%part == 1 )then
                write(logfhandle,'(A,I0,A,F8.3,A)') &
                    &'>>> NU filter state ', state, ' matching low-pass limit for next iteration: ', selected_lp, ' A'
            endif
        end subroutine record_nu_alignment_lowpass_limit

        subroutine cleanup_nu_aux_images()
            integer :: i
            if( allocated(nu_aux_even) )then
                do i = 1, size(nu_aux_even)
                    call nu_aux_even(i)%kill
                enddo
                deallocate(nu_aux_even)
            endif
            if( allocated(nu_aux_odd) )then
                do i = 1, size(nu_aux_odd)
                    call nu_aux_odd(i)%kill
                enddo
                deallocate(nu_aux_odd)
            endif
        end subroutine cleanup_nu_aux_images

    end subroutine nonuniform_filter_state

end module simple_nu_state_filter
