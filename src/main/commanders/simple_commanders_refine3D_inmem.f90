!@descr: in-memory probabilistic alignment helper (not wired into main flow)
module simple_commanders_refine3D_inmem
use simple_core_module_api
use simple_cmdline,             only: cmdline
use simple_builder,             only: builder
use simple_parameters,          only: parameters
use simple_strategy2D3D_common, only: sample_ptcls4fillin, sample_ptcls4update
use simple_eul_prob_tab,        only: eul_prob_tab
use simple_euclid_sigma2,       only: euclid_sigma2
use simple_image,               only: image
use simple_qsys_env,            only: qsys_env
use simple_qsys_funs,           only: qsys_cleanup
implicit none

public :: exec_prob_align_inmem
private
#include "simple_local_flags.inc"

contains

    !> In-memory probabilistic alignment. Caller supplies a prepared builder and parameters.
    subroutine exec_prob_align_inmem( build, params, cline )
        type(builder),    intent(inout) :: build
        type(parameters), intent(inout) :: params
        class(cmdline),   intent(inout) :: cline
        integer, allocatable :: pinds(:)
        type(string)        :: fname
        type(eul_prob_tab)  :: eulprob_obj_glob
        type(cmdline)       :: cline_prob_tab
        type(qsys_env)      :: qenv
        type(chash)         :: job_descr
        integer :: nptcls, ipart
        if( params%startit == 1 ) call build%spproj_field%clean_entry('updatecnt', 'sampled')
        if( params%l_fillin .and. mod(params%startit,5) == 0 )then
            call sample_ptcls4fillin(build, [1,params%nptcls], .true., nptcls, pinds)
        else
            call sample_ptcls4update(params, build, [1,params%nptcls], .true., nptcls, pinds)
        endif
        call build%spproj%write_segment_inside(params%oritype)
        call eulprob_obj_glob%new(params, build, pinds)
        cline_prob_tab = cline
        call cline_prob_tab%set('prg', 'prob_tab' )
        if( .not.cline_prob_tab%defined('nparts') )then
            call exec_prob_tab_inmem(build, params, cline_prob_tab)
        else
            call qenv%new(params, params%nparts, nptcls=params%nptcls)
            call cline_prob_tab%gen_job_descr(job_descr)
            call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        endif
        if( str_has_substr(params%refine, 'prob_state') )then
            do ipart = 1, params%nparts
                fname = string(DIST_FBODY)//int2str_pad(ipart,params%numlen)//'.dat'
                call eulprob_obj_glob%read_state_tab(fname)
            enddo
            call eulprob_obj_glob%state_assign
        else
            do ipart = 1, params%nparts
                fname = string(DIST_FBODY)//int2str_pad(ipart,params%numlen)//'.dat'
                call eulprob_obj_glob%read_tab_to_glob(fname)
            enddo
            call eulprob_obj_glob%ref_assign
        endif
        fname = string(ASSIGNMENT_FBODY)//'.dat'
        call eulprob_obj_glob%write_assignment(fname)
        call eulprob_obj_glob%kill
        call cline_prob_tab%kill
        call qenv%kill
        call job_descr%kill
        call qsys_cleanup(params)
    end subroutine exec_prob_align_inmem

    !> In-memory prob_tab path used by exec_prob_align_inmem when running locally.
    subroutine exec_prob_tab_inmem( build, params, cline )
        use simple_strategy2D3D_common, only: set_bp_range, prepare_refs_sigmas_ptcls, build_batch_particles, killimgbatch
        use simple_pftc_srch_api,       only: polarft_calc
        use simple_imgarr_utils,        only: dealloc_imgarr
        type(builder),    intent(inout) :: build
        type(parameters), intent(inout) :: params
        class(cmdline),   intent(inout) :: cline
        integer, allocatable :: pinds(:)
        type(image), allocatable :: tmp_imgs(:), tmp_imgs_pad(:)
        type(string) :: fname
        type(eul_prob_tab) :: eulprob_obj_part
        integer :: nptcls
        call set_bp_range(params, build, cline)
        if( build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls, pinds)
        else
            THROW_HARD('exec_prob_tab_inmem requires prior particle sampling (in exec_prob_align_inmem)')
        endif
        call prepare_refs_sigmas_ptcls(params, build, cline, tmp_imgs, tmp_imgs_pad, nptcls, params%which_iter,&
                                        do_polar=(params%l_polar .and. (.not.cline%defined('vol1'))) )
        call build_batch_particles(params, build, nptcls, pinds, tmp_imgs, tmp_imgs_pad)
        call eulprob_obj_part%new(params, build, pinds)
        fname = string(DIST_FBODY)//int2str_pad(params%part,params%numlen)//'.dat'
        if( str_has_substr(params%refine, 'prob_state') )then
            call eulprob_obj_part%fill_tab_state_only
            call eulprob_obj_part%write_state_tab(fname)
        else
            call eulprob_obj_part%fill_tab
            call eulprob_obj_part%write_tab(fname)
        endif
        call eulprob_obj_part%kill
        call killimgbatch(build)
        call build%pftc%kill
        call dealloc_imgarr(tmp_imgs)
        call dealloc_imgarr(tmp_imgs_pad)
    end subroutine exec_prob_tab_inmem

end module simple_commanders_refine3D_inmem
