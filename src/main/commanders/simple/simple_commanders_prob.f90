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

type, extends(commander_base) :: commander_prob_tab2D
  contains
    procedure :: execute      => exec_prob_tab2D
end type commander_prob_tab2D

type, extends(commander_base) :: commander_prob_align2D
  contains
    procedure :: execute      => exec_prob_align2D
end type commander_prob_align2D

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
        logical :: do_polar_prepare
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
        do_polar_prepare = (params%l_polar .and. (.not.cline%defined('vol1')))
        ! PREPARE REFERENCES, SIGMAS, POLAR_CORRCALC, PTCLS
        call prepare_refs_sigmas_ptcls( params, build, cline, tmp_imgs, tmp_imgs_pad, nptcls, params%which_iter,&
                                        do_polar=do_polar_prepare )
        if( .not. do_polar_prepare )then
            call read_mask_filter_reproject_refvols(params, build, cline, nptcls)
            call build%vol%kill
            call build%vol_odd%kill
            call build%vol2%kill
        endif
        call build%pftc%memoize_refs
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
        logical :: do_polar_prepare
        call cline%set('mkdir', 'no')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.true.)
        call set_bp_range( params, build, cline )
        ! Sampling policy mirrors exec_prob_tab: only reproduce already sampled particles.
        if( build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls, pinds)
        else
            THROW_HARD('exec_prob_tab_neigh requires prior particle sampling (in exec_prob_align)')
        endif
        do_polar_prepare = (params%l_polar .and. (.not.cline%defined('vol1')))
        ! PREPARE REFERENCES, SIGMAS, POLAR_CORRCALC, PTCLS
        call prepare_refs_sigmas_ptcls( params, build, cline, tmp_imgs, tmp_imgs_pad, nptcls, params%which_iter,&
                                        do_polar=do_polar_prepare )
        if( .not. do_polar_prepare )then
            call read_mask_filter_reproject_refvols(params, build, cline, nptcls)
            call build%vol%kill
            call build%vol_odd%kill
            call build%vol2%kill
        endif
        call build%pftc%memoize_refs
        ! Build polar particle images
        call build_batch_particles(params, build, nptcls, pinds, tmp_imgs, tmp_imgs_pad)
        call eulprob_obj_part_neigh%new_neigh(params, build, pinds)
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
        if( params%startit == 1 ) call build%spproj_field%clean_entry('updatecnt', 'sampled')
        if( params%l_fillin .and. mod(params%startit,5) == 0 )then
            call sample_ptcls4fillin(params, build, [1,params%nptcls], .true., nptcls, pinds)
        else
            call sample_ptcls4update(params, build, [1,params%nptcls], .true., nptcls, pinds)
        endif
        call build%spproj%write_segment_inside(params%oritype)
        ! Global object only needs sampled-set maps before reading partition sparse tables.
        ! The neighborhood scoring itself is performed in each prob_tab_neigh partition job.
        call eulprob_obj_glob_neigh%new_neigh(params, build, pinds)
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

    subroutine exec_prob_tab2D( self, cline )
        use simple_strategy2D3D_common, only: set_bp_range2D
        use simple_strategy2D_matcher,  only: set_b_p_ptrs2D, prep_batch_particles2D, preppftc4align2D, prep_polar_pftc4align2D, &
                                              build_batch_particles2D, clean_batch_particles2D
        use simple_classaverager,       only: cavger_new, cavger_read_all, cavger_kill
        use simple_eul_prob_tab2D,      only: eul_prob_tab2D
        class(commander_prob_tab2D), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        integer,     allocatable :: pinds(:)
        type(string)             :: fname
        type(builder)            :: build
        type(parameters)         :: params
        type(eul_prob_tab2D)     :: eulprob_obj_part
        real    :: frac_srch_space
        integer :: nptcls
        logical :: l_polar, l_use_polar_refs, l_alloc_read_cavgs
        call cline%set('mkdir', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        if( build%spproj_field%get_nevenodd() == 0 )then
            call build%spproj_field%partition_eo
            call build%spproj%write_segment_inside(params%oritype, params%projfile)
        endif
        l_polar          = trim(params%polar).eq.'yes'
        l_use_polar_refs = l_polar .and. (params%which_iter > 1)
        frac_srch_space  = build%spproj_field%get_avg('frac')
        call set_bp_range2D(params, build, cline, params%which_iter, frac_srch_space)
        ! reproduce particle sampling from exec_prob_align2D
        if( build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls, pinds)
        else
            THROW_HARD('exec_prob_tab2D requires prior particle sampling (in exec_prob_align2D)')
        endif
        call set_b_p_ptrs2D(params, build)
        call prep_batch_particles2D(nptcls)
        ! mirror cluster2D_exec reference setup: polar refs only for polar=yes and iter>1
        if( l_use_polar_refs )then
            call prep_polar_pftc4align2D(nptcls, params%which_iter, .false.)
        else
            l_alloc_read_cavgs = l_distr_worker_glob .or. (params%which_iter==1)
            call cavger_new(params, build, alloccavgs=l_alloc_read_cavgs)
            if( .not. cline%defined('refs') ) THROW_HARD('exec_prob_tab2D requires refs on the command line')
            call cavger_read_all
            call preppftc4align2D(nptcls, params%which_iter, .false.)
        endif
        ! build polar particle images
        call build_batch_particles2D(nptcls, pinds)
        ! fill and write the 2D probability table
        call eulprob_obj_part%new(params, build, pinds)
        call eulprob_obj_part%fill_tab
        fname = string(DIST_FBODY)//int2str_pad(params%part,params%numlen)//'.dat'
        call eulprob_obj_part%write_tab(fname)
        call eulprob_obj_part%kill
        call clean_batch_particles2D
        if( (.not.l_use_polar_refs).and.l_distr_worker_glob ) call cavger_kill
        call build%pftc%kill
        call build%kill_general_tbox
        call qsys_job_finished(params, string('simple_commanders_prob :: exec_prob_tab2D'))
        call simple_end('**** SIMPLE_PROB_TAB2D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_tab2D

    subroutine exec_prob_align2D( self, cline )
        use simple_eul_prob_tab2D,       only: eul_prob_tab2D
        use simple_strategy2D_matcher,   only: set_b_p_ptrs2D, sample_ptcls4update2D
        use simple_builder,              only: builder
        class(commander_prob_align2D), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        integer,       allocatable :: pinds(:)
        type(string)               :: fname
        type(builder)              :: build
        type(parameters)           :: params
        type(commander_prob_tab2D) :: xprob_tab2D
        type(eul_prob_tab2D)       :: eulprob_obj_glob
        type(cmdline)              :: cline_prob_tab
        type(qsys_env)             :: qenv
        type(chash)                :: job_descr
        integer :: nptcls, ipart
        call cline%set('mkdir',  'no')
        call cline%set('stream', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.false.)
        call set_b_p_ptrs2D(params, build)
        if( build%spproj_field%get_nevenodd() == 0 )then
            call build%spproj_field%partition_eo
            call build%spproj%write_segment_inside(params%oritype, params%projfile)
        endif
        if( params%startit == 1 ) call build%spproj_field%clean_entry('updatecnt', 'sampled')
        ! sample particles for this iteration
        call sample_ptcls4update2D([params%fromp,params%top], params%l_update_frac, nptcls, pinds)
        ! write sampling to project
        call build%spproj%write_segment_inside(params%oritype)
        ! build the global prob table (nclasses x nptcls)
        call eulprob_obj_glob%new(params, build, pinds)
        ! generate partition-wise dist tables
        cline_prob_tab = cline
        call cline_prob_tab%set('prg', 'prob_tab2D')
        if( .not. cline_prob_tab%defined('nparts') )then
            call xprob_tab2D%execute(cline_prob_tab)
        else
            call qenv%new(params, params%nparts, nptcls=params%nptcls)
            call cline_prob_tab%gen_job_descr(job_descr)
            call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        endif
        ! merge all partition tables into global
        do ipart = 1, params%nparts
            fname = string(DIST_FBODY)//int2str_pad(ipart,params%numlen)//'.dat'
            call eulprob_obj_glob%read_tab_to_glob(fname)
        end do
        ! global probabilistic class assignment
        call eulprob_obj_glob%ref_assign
        ! write assignment to file
        fname = string(ASSIGNMENT_FBODY)//'.dat'
        call eulprob_obj_glob%write_assignment(fname)
        ! cleanup
        call eulprob_obj_glob%kill
        call cline_prob_tab%kill
        call qenv%kill
        call job_descr%kill
        call build%kill_general_tbox
        call qsys_job_finished(params, string('simple_commanders_prob :: exec_prob_align2D'))
        call qsys_cleanup(params)
        call simple_end('**** SIMPLE_PROB_ALIGN2D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_align2D

end module simple_commanders_prob