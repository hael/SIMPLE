!@descr: Strategy pattern for refine_motion_model
module simple_refine_motion_model_strategy
use simple_commanders_api
use simple_parameters, only: parameters
use simple_cmdline,    only: cmdline
use simple_qsys_env,   only: qsys_env
use simple_sp_project, only: sp_project
use simple_builder,    only: builder
implicit none

public :: refine_motion_model_strategy
public :: refine_motion_model_inmem_strategy
public :: refine_motion_model_distr_strategy
public :: create_refine_motion_model_strategy
private
#include "simple_local_flags.inc"

! --------------------------------------------------------------------
! Strategy interface
! --------------------------------------------------------------------
type, abstract :: refine_motion_model_strategy
contains
    procedure :: apply_defaults => strategy_apply_defaults
    procedure(init_interface),           deferred :: initialize
    procedure(exec_interface),           deferred :: execute
    procedure(finalize_interface),       deferred :: finalize_run
    procedure(cleanup_interface),        deferred :: cleanup
    procedure(endmsg_interface),         deferred :: end_message
end type refine_motion_model_strategy

! Worker/shared-memory strategy
type, extends(refine_motion_model_strategy) :: refine_motion_model_inmem_strategy
contains
    procedure :: initialize     => inmem_initialize
    procedure :: execute        => inmem_execute
    procedure :: finalize_run   => inmem_finalize_run
    procedure :: cleanup        => inmem_cleanup
    procedure :: end_message    => inmem_end_message
end type refine_motion_model_inmem_strategy

! Distributed-master strategy
type, extends(refine_motion_model_strategy) :: refine_motion_model_distr_strategy
    type(parameters)                 :: params_snapshot
    type(sp_project)                 :: spproj
    type(qsys_env)                   :: qenv
    type(chash)                      :: job_descr
    type(chash), allocatable         :: part_params(:)
    integer, allocatable             :: parts(:,:)
    integer                          :: nmovies   = 0
    logical                          :: skip_run  = .false.
contains
    procedure :: initialize     => distr_initialize
    procedure :: execute        => distr_execute
    procedure :: finalize_run   => distr_finalize_run
    procedure :: cleanup        => distr_cleanup
    procedure :: end_message    => distr_end_message
end type refine_motion_model_distr_strategy

abstract interface
    subroutine init_interface(self, params, cline)
        import :: refine_motion_model_strategy, parameters, cmdline
        class(refine_motion_model_strategy), intent(inout) :: self
        type(parameters),          intent(inout) :: params
        class(cmdline),            intent(inout) :: cline
    end subroutine init_interface

    subroutine exec_interface(self, params, cline)
        import :: refine_motion_model_strategy, parameters, cmdline
        class(refine_motion_model_strategy), intent(inout) :: self
        type(parameters),          intent(inout) :: params
        class(cmdline),            intent(inout) :: cline
    end subroutine exec_interface

    subroutine finalize_interface(self, params, cline)
        import :: refine_motion_model_strategy, parameters, cmdline
        class(refine_motion_model_strategy), intent(inout) :: self
        type(parameters),          intent(in)    :: params
        class(cmdline),            intent(inout) :: cline
    end subroutine finalize_interface

    subroutine cleanup_interface(self, params, cline)
        import :: refine_motion_model_strategy, parameters, cmdline
        class(refine_motion_model_strategy), intent(inout) :: self
        type(parameters),          intent(in)    :: params
        class(cmdline),            intent(inout) :: cline
    end subroutine cleanup_interface

    function endmsg_interface(self) result(msg)
        import :: refine_motion_model_strategy
        class(refine_motion_model_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
    end function endmsg_interface
end interface

contains

    subroutine strategy_apply_defaults(self, cline)
        class(refine_motion_model_strategy), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        call validate_refine_motion_model_cline(cline)
        call set_refine_motion_model_defaults(cline)
    end subroutine strategy_apply_defaults

    ! --------------------------------------------------------------------
    ! Factory
    ! --------------------------------------------------------------------
    function create_refine_motion_model_strategy(cline) result(strategy)
        class(cmdline), intent(in) :: cline
        class(refine_motion_model_strategy), allocatable :: strategy
        logical :: is_master
        ! master if nparts defined but no explicit worker range/part.
        is_master = cline%defined('nparts') .and. (.not.cline%defined('part')) &
                   .and. (.not.cline%defined('fromp')) .and. (.not.cline%defined('top'))
        if( is_master )then
            allocate(refine_motion_model_distr_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DISTRIBUTED REFINE_MOTION_MODEL (MASTER)'
        else
            allocate(refine_motion_model_inmem_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> REFINE_MOTION_MODEL (WORKER / SHARED-MEMORY)'
        endif
    end function create_refine_motion_model_strategy

    ! --------------------------------------------------------------------
    ! Shared defaults + early validation (private; no common module)
    ! --------------------------------------------------------------------
    subroutine validate_refine_motion_model_cline(cline)
        class(cmdline), intent(inout) :: cline
    end subroutine validate_refine_motion_model_cline

    subroutine set_refine_motion_model_defaults(cline)
        class(cmdline), intent(inout) :: cline
        if( .not. cline%defined('oritype')      ) call cline%set('oritype',     'ptcl3D')
        if( .not. cline%defined('mkdir')        ) call cline%set('mkdir',          'yes')
        if( .not. cline%defined('pcontrast')    ) call cline%set('pcontrast',    'black')
        if( .not. cline%defined('backgr_subtr') ) call cline%set('backgr_subtr',    'no')
        if( .not. cline%defined('wfloat16')     ) call cline%set('wfloat16',        'no')
        if( .not.cline%defined('fromf')         ) call cline%set('fromf',              1)
        if( .not.cline%defined('tof')           ) call cline%set('tof',                0)
        if( .not.cline%defined('stepf')         ) call cline%set('stepf',              5)
    end subroutine set_refine_motion_model_defaults

    subroutine validate_fraction_range( fromf, tof, stepf )
        integer, intent(in) :: fromf, tof, stepf
        if( fromf < 1 ) THROW_HARD('Starting frame must be positive')
        if( tof < fromf ) THROW_HARD('Final frame must not precede starting frame')
        if( stepf < 1 ) THROW_HARD('Incremental frame step size must be positive')
    end subroutine validate_fraction_range

    pure integer function count_fractioned_stacks( fromf, tof, stepf )
        integer, intent(in) :: fromf, tof, stepf
        integer :: nframes
        nframes = tof - fromf + 1
        count_fractioned_stacks = (nframes + stepf - 1) / stepf
    end function count_fractioned_stacks

    ! Full frame blocks; anchor a final partial block at tof, overlapping its predecessor.
    ! If the entire selected range is shorter than stepf, use that range once.
    pure function fraction_frame_range( index, fromf, tof, stepf ) result( frames )
        integer, intent(in) :: index, fromf, tof, stepf
        integer :: frames(2)
        frames(2) = min(fromf - 1 + index * stepf, tof)
        frames(1) = max(fromf, frames(2) - stepf + 1)
    end function fraction_frame_range

    function fraction_output_dir( fromf, tof ) result( dirname )
        integer, intent(in) :: fromf, tof
        type(string)        :: dirname
        dirname = 'frames_'//int2str(fromf)//'-'//int2str(tof)
    end function fraction_output_dir

    ! ====================================================================
    ! WORKER / SHARED-MEMORY STRATEGY
    ! ====================================================================

    subroutine inmem_initialize(self, params, cline)
        class(refine_motion_model_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        class(cmdline),                  intent(inout) :: cline
        call params%new(cline)
        ! mirror original worker routine
        call cline%set('mkdir', 'no')
    end subroutine inmem_initialize

    subroutine inmem_execute(self, params, cline)
        use simple_particle_extractor, only: ptcl_extractor, dealloc_particles_frames
        use simple_motion_model,       only: motion_model
        use simple_matcher_ptcl_io,    only: prepimgbatch, killimgbatch
        class(refine_motion_model_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        class(cmdline),                  intent(inout) :: cline
        type(sp_project), allocatable :: fraction_projects(:)
        type(sp_project)        :: spproj, spproj_in
        type(builder)           :: build
        type(ori)               :: o_mov, o_stk
        type(stack_io)          :: stkio_w
        type(ptcl_extractor)    :: extractor
        type(motion_model)      :: mcmodel
        type(string)            :: mov_name, ext, stack, range, str, frame_dir
        logical,    allocatable :: mov_mask(:), ptcl_mask(:), fraction_mask(:)
        integer,    allocatable :: mic2stk_inds(:), boxcoords(:,:), coords(:,:)
        real    :: prev_shift(2),shift2d(2),shift3d(2),prev_shift_sc(2), translation(2)
        real    :: prev_center_sc(2), stk_min, stk_max, stk_mean, stk_sdev
        integer :: prev_pos(2), new_pos(2), ishift(2), prev_center(2), new_center(2)
        integer :: ffromto(2), i,ind, imov, iptcl, nmovs, prev_box, box_foo, cnt, nmovies_tot, stk_ind
        integer :: fromp, top, istk, nptcls2extract, nptcls, nfractioned_stacks, imov_local
        logical :: l_write16bits, l_from3D
        call validate_input
        params%msk    = RADFRAC_NORM_EXTRACT * real(params%box/2)
        l_write16bits = trim(params%wfloat16).eq.'yes'
        l_from3D      = trim(params%oritype).eq.'ptcl3D'
        ! proceed
        nfractioned_stacks = 0
        if( params%tof > 0 )then
            call validate_fraction_range(params%fromf, params%tof, params%stepf)
            nfractioned_stacks = count_fractioned_stacks(params%fromf, params%tof, params%stepf)
            write(logfhandle,'(A,I0)') '>>> Number of fractioned stacks to extract: ', nfractioned_stacks
            allocate(fraction_projects(nfractioned_stacks))
            allocate(fraction_mask(nmovs), source=.true.)
            do i = 1, nfractioned_stacks
                call fraction_projects(i)%os_stk%new(nmovs, is_ptcl=.false.)
            enddo
        endif
        if( nmovs > 0 )then
            call build%build_general_tbox(params, cline, do3d=.false.)
            ! More input validation
            do imov = 1, nmovies_tot
                if( .not.mov_mask(imov) )cycle
                call spproj_in%os_mic%get_ori(imov, o_mov)
                call o_mov%getter('movie', mov_name)
                if( .not.file_exists(mov_name) )cycle
                stk_ind = mic2stk_inds(imov)
                call spproj_in%os_stk%get_ori(stk_ind, o_stk)
                box_foo = o_stk%get_int('box')
                if( prev_box == 0 ) prev_box = box_foo
                if( prev_box /= box_foo ) THROW_HARD('Inconsistent box size; exec_refine_motion_model')
            enddo
            if( .not.cline%defined('box') ) params%box = prev_box
            if( is_odd(params%box) ) THROW_HARD('Box size must be of even dimension! exec_refine_motion_model')
            ! read particle segments
            call spproj_in%read_segment('ptcl2D', params%projfile)
            call spproj_in%read_segment('ptcl3D', params%projfile)
            allocate(ptcl_mask(spproj_in%os_ptcl2D%get_noris()),source=.false.)
            ! loop over movies within the specified range
            imov_local = 0
            do imov = params%fromp, params%top
                if( .not.mov_mask(imov) ) cycle
                imov_local = imov_local + 1
                call spproj_in%os_mic%get_ori(imov, o_mov)
                call o_mov%getter('movie', mov_name)
                stk_ind  = mic2stk_inds(imov)
                call spproj_in%os_stk%get_ori(stk_ind, o_stk)
                prev_box = o_stk%get_int('box')
                fromp    = o_stk%get_fromp()
                top      = o_stk%get_top()
                ext      = fname2ext(basename(mov_name))
                ! Updates particles positions within stage
                if( allocated(boxcoords) ) deallocate(boxcoords)
                allocate(boxcoords(2,fromp:top),source=0)
                !$omp parallel do default(shared) proc_bind(close) schedule(static)&
                !$omp private(iptcl,prev_pos,prev_shift,prev_center,prev_center_sc,prev_shift_sc)&
                !$omp private(new_center,new_pos,translation,shift2d,shift3d,ishift)
                do iptcl = fromp, top
                    if( spproj_in%os_ptcl2D%get_state(iptcl) == 0 ) cycle
                    if( spproj_in%os_ptcl3D%get_state(iptcl) == 0 ) cycle
                    call spproj_in%get_boxcoords(iptcl, prev_pos)
                    shift2d = spproj_in%os_ptcl2D%get_2Dshift(iptcl)
                    shift3d = spproj_in%os_ptcl3D%get_2Dshift(iptcl)
                    if( l_from3D ) then
                        prev_shift = shift3d
                    else
                        prev_shift = shift2d
                    endif
                    ishift      = nint(prev_shift)
                    new_pos     = prev_pos - ishift
                    translation = -real(ishift)
                    shift2d     = shift2d + translation
                    shift3d     = shift3d + translation
                    if( prev_box /= params%box ) new_pos = new_pos + (prev_box-params%box)/2
                    call spproj_in%set_boxcoords(iptcl, new_pos)
                    call spproj_in%os_ptcl2D%set_shift(iptcl, shift2d)
                    call spproj_in%os_ptcl3D%set_shift(iptcl, shift3d)
                    boxcoords(:,iptcl) = new_pos
                    ptcl_mask(iptcl)   = .true.
                enddo
                !$omp end parallel do
                nptcls2extract = count(ptcl_mask(fromp:top))
                if( nptcls2extract > 0 )then
                    allocate(coords(2,nptcls2extract),source=0)
                    cnt = 0
                    do iptcl = fromp, top
                        if( .not.ptcl_mask(iptcl) ) cycle
                        cnt = cnt + 1
                        call spproj_in%os_ptcl2D%set(iptcl, 'indstk', cnt)
                        call spproj_in%os_ptcl3D%set(iptcl, 'indstk', cnt)
                        coords(:,cnt) = boxcoords(:,iptcl)
                    enddo
                    call prepimgbatch(params, build, nptcls2extract)
                    ! parse movie metatdata and frames
                    call extractor%init_model(o_mov, params)
                    call spproj_in%os_mic%set(imov, 'nptcls', nptcls2extract)
                    stack = string(EXTRACT_STK_FBODY)//get_fbody(basename(mov_name), ext)//STK_EXT
                    call spproj_in%os_stk%set(stk_ind, 'stk', simple_abspath(stack,check_exists=.false.))
                    ! extract all frames particles
                    call extractor%extract_all_particles_frames(coords, [params%fromf, params%tof], 1)
                    deallocate(coords)
                    ! Generate stacks for individual frame blocks
                    do ind = 1, nfractioned_stacks
                        ffromto = fraction_frame_range(ind, params%fromf, params%tof, params%stepf)
                        ! accumulate
                        call extractor%generate_cumulative_particles_stack(nptcls2extract, build%imgbatch,&
                            &ffromto, stk_min, stk_max, stk_mean, stk_sdev)
                        ! write frame-block stack to disk
                        frame_dir = fraction_output_dir(ffromto(1), ffromto(2))
                        call simple_mkdir(frame_dir, verbose=.false.)
                        range = '_'//int2str(ffromto(1))//'_'//int2str(ffromto(2))
                        stack = string(EXTRACT_STK_FBODY)//get_fbody(basename(mov_name), ext)//range//STK_EXT
                        stack = filepath(frame_dir, stack)
                        call stkio_w%open(stack, params%smpd, 'write', box=params%box, wfloat16=l_write16bits)
                        do i = 1, nptcls2extract
                            call stkio_w%write(i, build%imgbatch(i))
                        enddo
                        call stkio_w%close
                        call build%imgbatch(1)%update_header_stats(stack, [stk_min, stk_max, stk_mean, stk_sdev])
                        ! book-keeping
                        call spproj_in%os_stk%set(stk_ind, 'box',        params%box)
                        call spproj_in%os_stk%set(stk_ind, 'nptcls',     nptcls2extract)
                        call spproj_in%os_stk%set(stk_ind, 'nptcls_stk', nptcls2extract)
                        call spproj_in%os_stk%set(stk_ind, 'stk',        simple_abspath(stack,check_exists=.false.))
                        call fraction_projects(ind)%os_stk%transfer_ori(imov_local, spproj_in%os_stk, stk_ind)
                    enddo
                else
                    call spproj_in%os_stk%set(stk_ind,'state',0)
                    call spproj_in%os_mic%set(imov,'state',0)
                    mov_mask(imov) = .false.
                    mic2stk_inds(imov) = 0
                    fraction_mask(imov_local) = .false.
                    do i = 1,nfractioned_stacks
                        call fraction_projects(i)%os_stk%set(imov_local,'state',0)
                    enddo
                endif
            enddo
        endif
        ! Output a binary project containing the stk field for each frame range.
        if( allocated(fraction_projects) )then
            do ind = 1, nfractioned_stacks
                if( size(fraction_mask) > 0 ) call fraction_projects(ind)%os_stk%compress(fraction_mask)
                ffromto = fraction_frame_range(ind, params%fromf, params%tof, params%stepf)
                range = int2str(ffromto(1))//'-'//int2str(ffromto(2))//'_part'//int2str(params%part)
                str   = trim(PTCLS_FRACTIONS_FBODY)//range%to_char()//METADATA_EXT
                call fraction_projects(ind)%write(str)
                call fraction_projects(ind)%kill
            enddo
            deallocate(fraction_projects)
            nmovs = count(fraction_mask)
            deallocate(fraction_mask)
        endif
        call extractor%kill
        call killimgbatch(build)
        call dealloc_particles_frames
        call build%kill_general_tbox
        ! Output project
        call spproj%read_non_data_segments(params%projfile)
        call spproj%projinfo%set(1,'projname', get_fbody(params%outfile,METADATA_EXT,separator=.false.))
        call spproj%projinfo%set(1,'projfile', params%outfile)
        call spproj%os_mic%new(nmovs, is_ptcl=.false.)
        call spproj%os_stk%new(nmovs, is_ptcl=.false.)
        call spproj_in%read_segment('stk', params%projfile)
        cnt = 0
        do imov = params%fromp, params%top
            if( .not.mov_mask(imov) )cycle
            cnt = cnt+1
            call spproj%os_mic%transfer_ori(cnt, spproj_in%os_mic, imov)
            stk_ind = mic2stk_inds(imov)
            call spproj%os_stk%transfer_ori(cnt, spproj_in%os_stk, stk_ind)
        enddo
        if( nmovs > 0 ) then
            nptcls = count(ptcl_mask)
            call spproj%os_ptcl2D%new(nptcls, is_ptcl=.true.)
            call spproj%os_ptcl3D%new(nptcls, is_ptcl=.true.)
            cnt = 0
            do iptcl = 1, size(ptcl_mask)
                if( .not.ptcl_mask(iptcl) )cycle
                cnt = cnt+1
                call spproj%os_ptcl2D%transfer_ori(cnt, spproj_in%os_ptcl2D, iptcl)
                call spproj%os_ptcl3D%transfer_ori(cnt, spproj_in%os_ptcl3D, iptcl)
            enddo
        endif
        call spproj%write(params%outfile)
        ! cleanup
        call spproj_in%kill
        call spproj%kill
        call o_mov%kill
        call o_stk%kill
        contains

            subroutine validate_input()
                call spproj_in%read_segment('mic', params%projfile)
                nmovies_tot = spproj_in%os_mic%get_noris()
                if( spproj_in%get_nmovies() == 0 ) THROW_HARD('No movie to process!')
                if( .not.spproj_in%os_mic%isthere('mcmodel') ) THROW_HARD('No Motion model to process!')
                call spproj_in%read_segment('stk', params%projfile)
                if( spproj_in%get_nstks() == 0 ) THROW_HARD('No particles extracted!')
                box_foo    = 0
                prev_box   = 0
                allocate(mic2stk_inds(nmovies_tot), source=0)
                allocate(mov_mask(nmovies_tot),     source=.false.)
                stk_ind = 0
                do imov = 1, nmovies_tot
                    if( imov > params%top ) exit
                    call spproj_in%os_mic%get_ori(imov, o_mov)
                    if( o_mov%isthere('state') )then
                        if( o_mov%get_state() == 0 )cycle
                    endif
                    if( .not. o_mov%isthere('movie')   ) cycle
                    if( .not. o_mov%isthere('mcmodel') ) cycle
                    do istk = stk_ind, spproj_in%os_stk%get_noris()
                        stk_ind = stk_ind + 1
                        if( spproj_in%os_stk%isthere(stk_ind,'state') )then
                            if( spproj_in%os_stk%get_state(stk_ind) == 1 ) exit
                        else
                            exit
                        endif
                    enddo
                    if( imov>=params%fromp .and. imov<=params%top )then
                        mov_mask(imov) = .true.
                        mic2stk_inds(imov) = stk_ind
                    endif
                enddo
                nmovs = count(mov_mask)
                if( nmovs > 0 ) then
                    write(logfhandle,'(A,I0)') '>>> Number of movies to process: ', nmovs
                    if( params%tof <= 0 ) then
                        do imov = 1, nmovies_tot
                            if( mov_mask(imov) )then
                                call spproj_in%os_mic%get_ori(imov, o_mov)
                                call o_mov%getter('mcmodel', str)
                                call mcmodel%read(str, params)
                                params%tof = mcmodel%nframes
                                call mcmodel%kill
                                exit
                            endif
                        enddo
                    endif
                endif
            end subroutine validate_input

    end subroutine inmem_execute

    subroutine inmem_finalize_run(self, params, cline)
        class(refine_motion_model_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        class(cmdline),                  intent(inout) :: cline
        call qsys_job_finished(params, string('simple_commanders_motion :: exec_refine_motion_model'))
    end subroutine inmem_finalize_run

    subroutine inmem_cleanup(self, params, cline)
        class(refine_motion_model_inmem_strategy), intent(inout) :: self
        type(parameters),                          intent(in)    :: params
        class(cmdline),                            intent(inout) :: cline
        ! No-op
    end subroutine inmem_cleanup

    function inmem_end_message(self) result(msg)
        class(refine_motion_model_inmem_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
        msg = '**** SIMPLE_REFINE_MOTION_MODEL NORMAL STOP ****'
    end function inmem_end_message

    ! ====================================================================
    ! DISTRIBUTED MASTER STRATEGY
    ! ====================================================================

    subroutine distr_initialize(self, params, cline)
        use simple_motion_model, only: motion_model
        class(refine_motion_model_distr_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        class(cmdline),                  intent(inout) :: cline
        type(ori)          :: o_mov
        type(motion_model) :: mcmodel
        type(string)       :: mov_name, model_name
        integer        :: imov, state, ipart
        integer        :: nmovies_valid
        call params%new(cline)
        call cline%set('mkdir', 'no')
        call self%spproj%read(params%projfile)
        if( self%spproj%get_nmovies() == 0 ) THROW_HARD('No movies to process!')
        if( self%spproj%get_nstks()   == 0 ) THROW_HARD('This project file does not contain stacks!')
        if( .not.self%spproj%os_mic%isthere('mcmodel') ) THROW_HARD('No Motion model to process!')
        self%nmovies = self%spproj%os_mic%get_noris()
        if( self%nmovies < params%nparts ) params%nparts = self%nmovies
        if( .not.cline%defined('box') ) then
            params%box = self%spproj%get_box()
            call cline%set('box', params%box)
        endif
        ! validate input
        nmovies_valid = 0
        do imov = 1, self%nmovies
            call self%spproj%os_mic%get_ori(imov, o_mov)
            state = 1
            if( o_mov%isthere('state') ) state = o_mov%get_state()
            if( state == 0 ) cycle
            if( .not. o_mov%isthere('movie')   )cycle
            if( .not. o_mov%isthere('mcmodel') )cycle
            call o_mov%getter('movie', mov_name)
            if( .not.file_exists(mov_name) )cycle
            call o_mov%getter('mcmodel', model_name)
            if( .not.file_exists(model_name) )cycle
            nmovies_valid = nmovies_valid + 1
            if( params%tof <= 0 )then
                call mcmodel%read(model_name, params)
                params%tof = mcmodel%nframes
                call mcmodel%kill
                call cline%set('tof', params%tof)
            endif
        enddo
        call o_mov%kill
        if( nmovies_valid == 0 )then
            THROW_WARN('No particles to extract! exec_refine_motion_model_distr')
            self%skip_run = .true.
            call self%spproj%kill
            call o_mov%kill
        else
            call validate_fraction_range(params%fromf, params%tof, params%stepf)
            ! build explicit fromp/top for each partition
            self%parts = split_nobjs_even(self%nmovies, params%nparts)
            allocate(self%part_params(params%nparts))
            do ipart = 1, params%nparts
                call self%part_params(ipart)%new(2)
                call self%part_params(ipart)%set('fromp', int2str(self%parts(ipart,1)))
                call self%part_params(ipart)%set('top',   int2str(self%parts(ipart,2)))
            end do
            call self%qenv%new(params, params%nparts)
            call cline%gen_job_descr(self%job_descr)
        endif
    end subroutine distr_initialize

    subroutine distr_execute(self, params, cline)
        use simple_oris, only: oris
        class(refine_motion_model_distr_strategy), intent(inout) :: self
        type(parameters),                          intent(inout) :: params
        class(cmdline),                            intent(inout) :: cline
        type(sp_project), allocatable :: spproj_parts(:)
        type(string),     allocatable :: stktab(:)
        type(oris)   :: os_stk
        type(string) :: ptcl_project, mic_project, frame_project, stk_project, range, frame_dir
        integer      :: numlen, ffromto(2), ipart, imic, istk, ind, nfractioned_stacks
        integer      :: nmovs, cnt, nstks, nptcls, i, stkind
        if( self%skip_run ) return
        call validate_fraction_range(params%fromf, params%tof, params%stepf)
        nfractioned_stacks = count_fractioned_stacks(params%fromf, params%tof, params%stepf)
        do ind = 1, nfractioned_stacks
            ffromto = fraction_frame_range(ind, params%fromf, params%tof, params%stepf)
            frame_dir = fraction_output_dir(ffromto(1), ffromto(2))
            call simple_mkdir(frame_dir, verbose=.false.)
        enddo
        ! schedule & run
        call self%qenv%gen_scripts_and_schedule_jobs( self%job_descr, algnfbody=string(ALGN_FBODY), &
            &part_params=self%part_params, array=L_USE_SLURM_ARR, extra_params=params)
        ! ASSEMBLY
        allocate(spproj_parts(params%nparts))
        numlen = len(int2str(params%nparts))
        do ind = 1, nfractioned_stacks
            ffromto = fraction_frame_range(ind, params%fromf, params%tof, params%stepf)
            call self%spproj%os_mic%kill
            call self%spproj%os_stk%kill
            call self%spproj%os_ptcl2D%kill
            call self%spproj%os_ptcl3D%kill
            ! Count micrographs across parts
            nmovs = 0
            do ipart = 1, params%nparts
                mic_project = ALGN_FBODY//int2str_pad(ipart,numlen)//METADATA_EXT
                call spproj_parts(ipart)%read_segment('mic', mic_project)
                nmovs = nmovs + spproj_parts(ipart)%os_mic%get_noris()
            enddo
            if( nmovs > 0 )then
                call self%spproj%os_mic%new(nmovs, is_ptcl=.false.)
                ! transfer micrographs + count stacks
                cnt   = 0
                nstks = 0
                do ipart = 1, params%nparts
                    ! micrographs
                    do imic = 1, spproj_parts(ipart)%os_mic%get_noris()
                        cnt = cnt + 1
                        call self%spproj%os_mic%transfer_ori(cnt, spproj_parts(ipart)%os_mic, imic)
                    enddo
                    call spproj_parts(ipart)%kill
                    ! stacks
                    range = int2str(ffromto(1))//'-'//int2str(ffromto(2))//'_part'//int2str(ipart)
                    stk_project = trim(PTCLS_FRACTIONS_FBODY)//range%to_char()//METADATA_EXT
                    call spproj_parts(ipart)%read_segment('stk', stk_project)
                    nstks = nstks + spproj_parts(ipart)%os_stk%get_noris()
                enddo
                if( nstks /= nmovs ) THROW_HARD('Inconstistent number of stacks in individual projects')
                ! stacks table
                call os_stk%new(nstks, is_ptcl=.false.)
                allocate(stktab(nstks))
                cnt = 0
                do ipart = 1, params%nparts
                    do istk = 1, spproj_parts(ipart)%os_stk%get_noris()
                        cnt = cnt + 1
                        call os_stk%transfer_ori(cnt, spproj_parts(ipart)%os_stk, istk)
                        stktab(cnt) = os_stk%get_str(cnt,'stk')
                    enddo
                    call spproj_parts(ipart)%kill
                enddo
                call self%spproj%add_stktab(stktab, os_stk)
                call os_stk%kill
                if( allocated(stktab) ) deallocate(stktab)
                ! transfer particle params (preserve stkind)
                cnt = 0
                do ipart = 1, params%nparts
                    ptcl_project = ALGN_FBODY//int2str_pad(ipart,numlen)//METADATA_EXT
                    call spproj_parts(ipart)%read_segment('ptcl2D', ptcl_project)
                    call spproj_parts(ipart)%read_segment('ptcl3D', ptcl_project)
                    nptcls = spproj_parts(ipart)%os_ptcl2D%get_noris()
                    if( nptcls /= spproj_parts(ipart)%os_ptcl3D%get_noris() )then
                        THROW_HARD('Inconsistent number of particles')
                    endif
                    do i = 1, nptcls
                        cnt    = cnt + 1
                        ! keep mapping to stacks
                        stkind = self%spproj%os_ptcl2D%get_int(cnt,'stkind')
                        call self%spproj%os_ptcl2D%transfer_ori(cnt, spproj_parts(ipart)%os_ptcl2D, i)
                        call self%spproj%os_ptcl3D%transfer_ori(cnt, spproj_parts(ipart)%os_ptcl3D, i)
                        call self%spproj%os_ptcl2D%set(cnt,'stkind',stkind)
                        call self%spproj%os_ptcl3D%set(cnt,'stkind',stkind)
                    enddo
                    call spproj_parts(ipart)%kill
                enddo
            endif
            ! Write frame range project
            frame_dir     = fraction_output_dir(ffromto(1), ffromto(2))
            frame_project = filepath(frame_dir, basename(params%projfile))
            call self%spproj%projinfo%set(1, 'projname', &
                &get_fbody(basename(params%projfile), METADATA_EXT, separator=.false.))
            call self%spproj%projinfo%set(1, 'projfile', simple_abspath(frame_project, check_exists=.false.))
            call self%spproj%write(frame_project)
            ! Remove per-part fraction projects
            do ipart = 1, params%nparts
                range = int2str(ffromto(1))//'-'//int2str(ffromto(2))//'_part'//int2str(ipart)
                stk_project = trim(PTCLS_FRACTIONS_FBODY)//range%to_char()//METADATA_EXT
                if( file_exists(stk_project) ) call del_file(stk_project)
            enddo
        enddo
        if( allocated(spproj_parts) ) deallocate(spproj_parts)
    end subroutine distr_execute

    subroutine distr_finalize_run(self, params, cline)
        class(refine_motion_model_distr_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        class(cmdline),                  intent(inout) :: cline
        ! No-op
    end subroutine distr_finalize_run

    subroutine distr_cleanup(self, params, cline)
        class(refine_motion_model_distr_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        class(cmdline),                  intent(inout) :: cline
        integer :: i
        if( .not. self%skip_run )then
            call qsys_cleanup(params)
        endif
        call self%spproj%kill
        call self%qenv%kill
        call self%job_descr%kill
        if( allocated(self%part_params) )then
            do i = 1, size(self%part_params)
                call self%part_params(i)%kill
            enddo
            deallocate(self%part_params)
        endif
        if( allocated(self%parts) ) deallocate(self%parts)
    end subroutine distr_cleanup

    function distr_end_message(self) result(msg)
        class(refine_motion_model_distr_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
        msg = '**** SIMPLE_REFINE_MOTION_MODEL_DISTR NORMAL STOP ****'
    end function distr_end_message

    pure logical function box_inside(ildim, coord, box)
        integer, intent(in) :: ildim(3), coord(2), box
        integer             :: fromc(2), toc(2)
        fromc = coord + 1
        toc   = fromc + (box - 1)
        box_inside = .true.
        if( any(fromc < 1) .or. toc(1) > ildim(1) .or. toc(2) > ildim(2) ) box_inside = .false.
    end function box_inside

end module simple_refine_motion_model_strategy
