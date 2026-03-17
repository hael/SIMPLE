module simple_commanders_prob
use simple_commanders_api
use simple_pftc_srch_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_prob_tab
  contains
    procedure :: execute      => exec_prob_tab
end type commander_prob_tab

type, extends(commander_base) :: commander_prob_tab_neigh
    contains
        procedure :: execute      => exec_prob_tab_neigh
end type commander_prob_tab_neigh

type, extends(commander_base) :: commander_prob_align
  contains
    procedure :: execute      => exec_prob_align
end type commander_prob_align

type, extends(commander_base) :: commander_prob_align_neigh
    contains
        procedure :: execute      => exec_prob_align_neigh
end type commander_prob_align_neigh

contains

    subroutine exec_prob_tab( self, cline )
        use simple_strategy2D3D_common
        use simple_eul_prob_tab, only: eul_prob_tab
        use simple_imgarr_utils, only: dealloc_imgarr
        class(commander_prob_tab), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        integer,     allocatable :: pinds(:)
        type(image), allocatable :: tmp_imgs(:), tmp_imgs_pad(:)
        type(string)             :: fname
        type(builder)            :: build
        type(parameters)         :: params
        type(eul_prob_tab)       :: eulprob_obj_part
        integer :: nptcls
        call cline%set('mkdir', 'no')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.true.)
        call set_bp_range( params, build, cline )
        ! The policy here ought to be that nothing is done with regards to sampling other than reproducing
        ! what was generated in the driver (prob_align, below). Sampling is delegated to prob_align (below)
        ! and merely reproduced here
        if( build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls, pinds)
        else
            THROW_HARD('exec_prob_tab requires prior particle sampling (in exec_prob_align)')
        endif
        ! PREPARE REFERENCES, SIGMAS, POLAR_CORRCALC, PTCLS
        call prepare_refs_sigmas_ptcls( params, build, cline, tmp_imgs, tmp_imgs_pad, nptcls, params%which_iter,&
                                        do_polar=(params%l_polar .and. (.not.cline%defined('vol1'))) )           
        ! Build polar particle images
        call build_batch_particles(params, build, nptcls, pinds, tmp_imgs, tmp_imgs_pad)
        ! Filling prob table in eul_prob_tab
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
        call build%kill_general_tbox
        call dealloc_imgarr(tmp_imgs)
        call dealloc_imgarr(tmp_imgs_pad)
        call qsys_job_finished(params, string('simple_commanders_refine3D :: exec_prob_tab'))
        call simple_end('**** SIMPLE_PROB_TAB NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_tab

    subroutine exec_prob_tab_neigh( self, cline )
        use simple_strategy2D3D_common
        use simple_eul_prob_tab_neigh, only: eul_prob_tab_neigh
        use simple_imgarr_utils,       only: dealloc_imgarr
        class(commander_prob_tab_neigh), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        integer,     allocatable :: pinds(:)
        type(image), allocatable :: tmp_imgs(:), tmp_imgs_pad(:)
        type(string)             :: fname
        type(builder)            :: build
        type(parameters)         :: params
        type(eul_prob_tab_neigh) :: eulprob_obj_part_neigh
        integer :: nptcls
        call cline%set('mkdir', 'no')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.true.)
        ! exception handling for required neighborhood data structures and incompatible refine policies
        if( str_has_substr(params%refine, 'prob_state') )then
            THROW_HARD('exec_prob_tab_neigh does not support refine=prob_state; use the dense probabilistic state path')
        endif
        if( .not. allocated(build%subspace_inds) )then
            THROW_HARD('exec_prob_tab_neigh requires neighborhood representative projections; enable l_neigh and set nspace_sub')
        endif
        if( size(build%subspace_inds) /= params%nspace_sub )then
            THROW_HARD('exec_prob_tab_neigh: size(subspace_inds) must equal nspace_sub')
        endif
        if( .not. allocated(build%subspace_full2sub_map) )then
            THROW_HARD('exec_prob_tab_neigh requires full-space neighborhood labels (subspace_full2sub_map)')
        endif
        if( size(build%subspace_full2sub_map) /= params%nspace )then
            THROW_HARD('exec_prob_tab_neigh: size(subspace_full2sub_map) must equal nspace')
        endif
        call set_bp_range( params, build, cline )
        ! Sampling policy mirrors exec_prob_tab: only reproduce already sampled particles.
        if( build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls, pinds)
        else
            THROW_HARD('exec_prob_tab_neigh requires prior particle sampling (in exec_prob_align)')
        endif
        ! PREPARE REFERENCES, SIGMAS, POLAR_CORRCALC, PTCLS
        call prepare_refs_sigmas_ptcls( params, build, cline, tmp_imgs, tmp_imgs_pad, nptcls, params%which_iter,&
                                        do_polar=(params%l_polar .and. (.not.cline%defined('vol1'))) )
        ! Build polar particle images
        call build_batch_particles(params, build, nptcls, pinds, tmp_imgs, tmp_imgs_pad)
        call eulprob_obj_part_neigh%new(params, build, pinds)
        call eulprob_obj_part_neigh%fill_tab
        fname = string(DIST_FBODY)//'_neigh_'//int2str_pad(params%part,params%numlen)//'.dat'
        call eulprob_obj_part_neigh%write_tab(fname)
        call eulprob_obj_part_neigh%kill
        call killimgbatch(build)
        call build%pftc%kill
        call build%kill_general_tbox
        call dealloc_imgarr(tmp_imgs)
        call dealloc_imgarr(tmp_imgs_pad)
        call qsys_job_finished(params, string('simple_commanders_refine3D :: exec_prob_tab_neigh'))
        call simple_end('**** SIMPLE_PROB_TAB_NEIGH NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_tab_neigh

    subroutine exec_prob_align( self, cline )
        use simple_eul_prob_tab,        only: eul_prob_tab
        use simple_strategy2D3D_common, only: sample_ptcls4fillin, sample_ptcls4update
        use simple_builder,             only: builder
        class(commander_prob_align), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        integer,     allocatable :: pinds(:)
        type(string)             :: fname
        type(builder)            :: build
        type(parameters)         :: params
        type(commander_prob_tab) :: xprob_tab
        type(eul_prob_tab)       :: eulprob_obj_glob
        type(cmdline)            :: cline_prob_tab
        type(qsys_env)           :: qenv
        type(chash)              :: job_descr
        integer :: nptcls, ipart
        call cline%set('mkdir',  'no')
        call cline%set('stream', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
        if( params%startit == 1 ) call build%spproj_field%clean_entry('updatecnt', 'sampled')
        ! sampled incremented
        if( params%l_fillin .and. mod(params%startit,5) == 0 )then
            call sample_ptcls4fillin(params, build, [1,params%nptcls], .true., nptcls, pinds)
        else
            call sample_ptcls4update(params, build, [1,params%nptcls], .true., nptcls, pinds)
        endif
        ! communicate to project file
        call build%spproj%write_segment_inside(params%oritype)
        ! more prep
        call eulprob_obj_glob%new(params, build, pinds)
        ! generating all corrs on all parts
        cline_prob_tab = cline
        call cline_prob_tab%set('prg', 'prob_tab' ) ! required for distributed call
        ! execution
        if( .not.cline_prob_tab%defined('nparts') )then
            call xprob_tab%execute(cline_prob_tab)
        else
            ! setup the environment for distributed execution
            call qenv%new(params, params%nparts, nptcls=params%nptcls)
            call cline_prob_tab%gen_job_descr(job_descr)
            ! schedule
            call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        endif
        ! reading corrs from all parts
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
        ! write the iptcl->(iref,istate) assignment
        fname = string(ASSIGNMENT_FBODY)//'.dat'
        call eulprob_obj_glob%write_assignment(fname)
        ! cleanup
        call eulprob_obj_glob%kill
        call cline_prob_tab%kill
        call qenv%kill
        call job_descr%kill
        call build%kill_general_tbox
        call qsys_job_finished(params, string('simple_commanders_refine3D :: exec_prob_align'))
        call qsys_cleanup(params)
        call simple_end('**** SIMPLE_PROB_ALIGN NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_align

    subroutine exec_prob_align_neigh( self, cline )
        use simple_eul_prob_tab_neigh,  only: eul_prob_tab_neigh
        use simple_strategy2D3D_common, only: sample_ptcls4fillin, sample_ptcls4update
        use simple_builder,             only: builder
        class(commander_prob_align_neigh), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        integer,           allocatable :: pinds(:)
        type(string)                   :: fname
        type(builder)                  :: build
        type(parameters)               :: params
        type(commander_prob_tab_neigh) :: xprob_tab_neigh
        type(eul_prob_tab_neigh)       :: eulprob_obj_glob_neigh
        type(cmdline)                  :: cline_prob_tab
        type(qsys_env)                 :: qenv
        type(chash)                    :: job_descr
        integer :: nptcls
        call cline%set('mkdir',  'no')
        call cline%set('stream', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
        ! exception handling for required neighborhood setup
        if( .not. params%l_neigh )then
            THROW_HARD('exec_prob_align_neigh requires l_neigh=yes')
        endif
        if( str_has_substr(params%refine, 'prob_state') )then
            THROW_HARD('exec_prob_align_neigh does not support refine=prob_state; use exec_prob_align')
        endif
        if( .not. allocated(build%subspace_inds) )then
            THROW_HARD('exec_prob_align_neigh requires neighborhood representative projections; enable l_neigh and set nspace_sub')
        endif
        if( size(build%subspace_inds) /= params%nspace_sub )then
            THROW_HARD('exec_prob_align_neigh: size(subspace_inds) must equal nspace_sub')
        endif
        if( .not. allocated(build%subspace_full2sub_map) )then
            THROW_HARD('exec_prob_align_neigh requires full-space neighborhood labels (subspace_full2sub_map)')
        endif
        if( size(build%subspace_full2sub_map) /= params%nspace )then
            THROW_HARD('exec_prob_align_neigh: size(subspace_full2sub_map) must equal nspace')
        endif
        if( params%startit == 1 ) call build%spproj_field%clean_entry('updatecnt', 'sampled')
        if( params%l_fillin .and. mod(params%startit,5) == 0 )then
            call sample_ptcls4fillin(params, build, [1,params%nptcls], .true., nptcls, pinds)
        else
            call sample_ptcls4update(params, build, [1,params%nptcls], .true., nptcls, pinds)
        endif
        call build%spproj%write_segment_inside(params%oritype)
        call eulprob_obj_glob_neigh%new(params, build, pinds)
        cline_prob_tab = cline
        call cline_prob_tab%set('prg', 'prob_tab_neigh')
        if( .not. cline_prob_tab%defined('nparts') )then
            call xprob_tab_neigh%execute(cline_prob_tab)
        else
            call qenv%new(params, params%nparts, nptcls=params%nptcls)
            call cline_prob_tab%gen_job_descr(job_descr)
            call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        endif
        call eulprob_obj_glob_neigh%read_tabs_to_glob(string(DIST_FBODY)//'_neigh_', params%nparts, params%numlen)
        call eulprob_obj_glob_neigh%ref_assign
        ! write the iptcl->(iref,istate) assignment
        fname = string(ASSIGNMENT_FBODY)//'.dat'
        call eulprob_obj_glob_neigh%write_assignment(fname)
        ! cleanup
        call eulprob_obj_glob_neigh%kill
        call cline_prob_tab%kill
        call qenv%kill
        call job_descr%kill
        call build%kill_general_tbox
        call qsys_job_finished(params, string('simple_commanders_refine3D :: exec_prob_align_neigh'))
        call qsys_cleanup(params)
        call simple_end('**** SIMPLE_PROB_ALIGN_NEIGH NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_align_neigh

end module simple_commanders_prob