!@descr: Strategy pattern for refine_motion_model
module simple_refine_motion_model_strategy
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
use simple_core_module_api
use simple_parameters,   only: parameters
use simple_cmdline,      only: cmdline
use simple_qsys_env,     only: qsys_env
use simple_qsys_funs,    only: qsys_cleanup, qsys_declare_part_finished
use simple_sp_project,   only: sp_project
use simple_stack_io,     only: stack_io
use simple_builder,      only: builder
use simple_motion_model, only: motion_model
implicit none

public :: refine_motion_model_strategy
public :: refine_motion_model_inmem_strategy
public :: refine_motion_model_distr_strategy
public :: create_refine_motion_model_strategy
public :: prepare_fraction_dose_manifest, reconstruct_motion_fractions
private
#include "simple_local_flags.inc"

real(dp),         parameter :: DOSE_FRACTION_RELTOL = 1.e-5_dp
real(dp),         parameter :: DOSE_FRACTION_ABSTOL = 1.e-7_dp
real(dp),         parameter :: DOSE_VOLTAGE_TOL     = 1.e-3_dp
integer,          parameter :: DOSE_PARTITION_ROW_WIDTH = 15 !< I12, colon, half digit, semicolon
character(len=*), parameter :: DOSE_MANIFEST_FNAME = 'dose_series.txt'

! --------------------------------------------------------------------
! Strategy interface
! --------------------------------------------------------------------
type, abstract :: refine_motion_model_strategy
contains
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
        use simple_particle_extractor, only: ptcl_extractor
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
        logical :: l_from3D
        call validate_input
        params%msk    = RADFRAC_NORM_EXTRACT * real(params%box/2)
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
                    call extractor%extract_all_particles_frames(coords, [params%fromf, params%tof])
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
                        call stkio_w%open(stack, params%smpd, 'write', box=params%box, wfloat16=.true.)
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
                    mov_mask(imov)            = .false.
                    mic2stk_inds(imov)        = 0
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
        call qsys_declare_part_finished(params, string('simple_commanders_motion :: exec_refine_motion_model'))
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
        use simple_eo_partition, only: ensure_field_halves
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
            call ensure_field_halves(self%spproj, 'ptcl3D', params%projfile)
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
        class(refine_motion_model_distr_strategy), intent(inout) :: self
        type(parameters),                          intent(inout) :: params
        class(cmdline),                            intent(inout) :: cline
        type(string), allocatable :: frame_projects(:)
        integer,      allocatable :: frame_ranges(:,:)
        type(string) :: frame_dir
        integer      :: ffromto(2), ind, nfractioned_stacks
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
        call assemble_motion_fractions(self, params, frame_projects, frame_ranges)
        if( self%spproj%os_ptcl3D%get_noris() > 0 )then
            call reconstruct_motion_fractions(params, cline, frame_projects, frame_ranges, &
                &self%spproj%os_ptcl3D%get_n('state'))
        else
            THROW_WARN('No extracted particles available for interval reconstruction')
        endif
    end subroutine distr_execute

    !> Read extraction part projects; write interval projects and remove consumed stack-table parts.
    subroutine assemble_motion_fractions(self, params, frame_projects, frame_ranges)
        use simple_oris, only: oris
        class(refine_motion_model_distr_strategy), intent(inout) :: self
        type(parameters),                          intent(in)    :: params
        type(string), allocatable,                 intent(out)   :: frame_projects(:)
        integer,      allocatable,                 intent(out)   :: frame_ranges(:,:)
        type(sp_project), allocatable :: spproj_parts(:)
        type(string),     allocatable :: stktab(:)
        type(oris)   :: os_stk
        type(string) :: ptcl_project, mic_project, frame_project, stk_project, range, frame_dir
        integer      :: numlen, ffromto(2), ipart, imic, istk, ind, nfractioned_stacks
        integer      :: nmovs, cnt, nstks, nptcls, i, stkind
        nfractioned_stacks = count_fractioned_stacks(params%fromf, params%tof, params%stepf)
        allocate(frame_projects(nfractioned_stacks), frame_ranges(2,nfractioned_stacks))
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
            frame_projects(ind) = simple_abspath(frame_project)
            frame_ranges(:,ind) = ffromto
            ! Remove per-part fraction projects
            do ipart = 1, params%nparts
                range = int2str(ffromto(1))//'-'//int2str(ffromto(2))//'_part'//int2str(ipart)
                stk_project = trim(PTCLS_FRACTIONS_FBODY)//range%to_char()//METADATA_EXT
                if( file_exists(stk_project) ) call del_file(stk_project)
            enddo
        enddo
        if( allocated(spproj_parts) ) deallocate(spproj_parts)
    end subroutine assemble_motion_fractions

    !> Reconstruct copied interval projects, register maps/FSCs, then export raw FSC observations.
    subroutine reconstruct_motion_fractions(params, cline, frame_projects, frame_ranges, nstates)
        use simple_rec3D_strategy,         only: run_rec3D
        use simple_dose_series_manifest,   only: dose_series_manifest
        use simple_refine3D_fnames,        only: refine3D_state_vol_fname, refine3D_fsc_fname
        !$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads
        type(parameters), intent(in)    :: params
        class(cmdline),   intent(inout) :: cline
        type(string),     intent(in)    :: frame_projects(:)
        integer,          intent(in)    :: frame_ranges(:,:)
        integer,          intent(in)    :: nstates
        type(sp_project)      :: spproj
        type(cmdline)         :: cline_rec
        type(dose_series_manifest) :: dose_manifest
        type(string)          :: cwd_saved, rec_dir, rec_project, volname, fscname
        type(string)          :: manifest_fname, manifest_message
        integer, allocatable  :: state1_rows(:), state1_halves(:), current_state1_rows(:), current_state1_halves(:)
        integer, allocatable  :: pinds(:), selected(:), halves(:), current_halves(:), state_pops(:)
        real,    allocatable  :: fsc_matrix(:,:,:), fsc_curve(:), resolution(:)
        real(dp), allocatable :: reciprocal_frequency(:)
        integer :: i, j, state, nselected, nptcls, old_out_states, nregistered
        integer :: saved_threads, saved_nthr_glob
        logical :: l_manifest
        if( size(frame_ranges,1) /= 2 ) THROW_HARD('Fraction ranges must have two endpoints')
        if( size(frame_ranges,2) /= size(frame_projects) ) THROW_HARD('Fraction project/range count mismatch')
        if( size(frame_projects) == 0 ) return
        if( params%cc_objfun /= OBJFUN_CC ) THROW_HARD('Fraction reconstruction requires objfun=cc')
        if( params%l_ml_reg ) THROW_HARD('Fraction reconstruction requires ml_reg=no')
        if( params%cc_emit_sigma /= 'no' ) THROW_HARD('Fraction reconstruction cannot emit sigma2')
        if( nstates < 1 ) THROW_HARD('Fraction reconstruction requires at least one state')
        resolution = get_resarr(params%box, params%smpd)
        reciprocal_frequency = [(real(j,dp)/(real(params%box,dp)*real(params%smpd,dp)), j=1,size(resolution))]
        ! Missing interval/state curves remain distinct from measured zero correlation.
        allocate(fsc_matrix(size(resolution),size(frame_projects),nstates), source=ieee_value(0.,ieee_quiet_nan))
        call simple_getcwd(cwd_saved)
        manifest_fname = filepath(cwd_saved, DOSE_MANIFEST_FNAME)
        if( file_exists(manifest_fname) ) THROW_HARD('Dose series manifest already exists in the reconstruction directory')
        l_manifest = .false.
        saved_nthr_glob = nthr_glob
        saved_threads   = nthr_glob
        !$ saved_threads = omp_get_max_threads()
        do i = 1, size(frame_projects)
            rec_dir = filepath(cwd_saved, 'rec_frames_'//int2str_pad(frame_ranges(1,i),2)//'-'// &
                &int2str_pad(frame_ranges(2,i),2))
            if( dir_exists(rec_dir) ) THROW_HARD('Reconstruction directory already exists: '//rec_dir%to_char())
            call simple_mkdir(rec_dir)
            rec_project = filepath(rec_dir, basename(frame_projects(i)))
            call simple_copy_file(frame_projects(i), rec_project)
            ! prepare and verify the copied project; the master owns it until the program runs
            call spproj%read(rec_project)
            call spproj%update_projinfo(rec_project)
            if( spproj%get_box() /= params%box ) THROW_HARD('Fraction reconstruction box mismatch')
            if( .not.is_equal(spproj%get_smpd(),params%smpd) ) THROW_HARD('Fraction reconstruction sampling mismatch')
            nptcls  = spproj%os_ptcl3D%get_noris()
            if( spproj%os_ptcl3D%get_n('state') /= nstates ) THROW_HARD('Fraction reconstruction state count changed')
            call spproj%os_ptcl3D%sample4rec([1,nptcls], nselected, pinds)
            call spproj%os_ptcl3D%get_state_rec_pops(nstates, state_pops)
            if( sum(state_pops) /= nselected ) THROW_HARD('Fraction reconstruction state populations disagree')
            if( .not.spproj%os_ptcl3D%isthere('eo') ) THROW_HARD('Fraction particles have no half-set labels')
            allocate(current_halves(nselected))
            do j = 1, nselected
                current_halves(j) = spproj%os_ptcl3D%get_eo(pinds(j))
            enddo
            if( any(current_halves < 0) .or. any(current_halves > 1) ) THROW_HARD('Invalid fraction half-set labels')
            if( allocated(selected) )then
                if( size(selected) /= nselected ) THROW_HARD('Fraction reconstruction particle count changed')
                if( any(selected /= pinds) ) THROW_HARD('Fraction reconstruction particle selection changed')
                if( any(halves /= current_halves) ) THROW_HARD('Fraction reconstruction half-set labels changed')
            else
                selected = pinds
                halves   = current_halves
            endif
            deallocate(current_halves)
            if( i == 1 )then
                l_manifest = prepare_fraction_dose_manifest(params, spproj, pinds, frame_ranges, dose_manifest, &
                    &state1_rows, state1_halves, manifest_message)
            else if( l_manifest )then
                call fraction_state1_rows(spproj, pinds, current_state1_rows, current_state1_halves)
                if( size(current_state1_rows) /= size(state1_rows) )then
                    l_manifest = .false.
                else
                    l_manifest = all(current_state1_rows == state1_rows) .and. all(current_state1_halves == state1_halves)
                endif
                if( .not. l_manifest ) manifest_message = 'State-1 reconstruction particles or half assignments changed between intervals'
            endif
            old_out_states = spproj%os_out%get_n('state')
            do state = 1, max(old_out_states,nstates)
                call spproj%remove_state_artifacts_from_osout(state)
            enddo
            call spproj%write(rec_project)
            call spproj%kill
            write(logfhandle,'(A,2(I0,1X),A,I0)') '>>> RECONSTRUCTING FRACTIONS ', frame_ranges(:,i), &
                &'PARTICLES: ', nselected
            if( nselected > 0 )then
                ! a nested, complete reconstruct3D run on its own parameters and builder
                call cline_rec%set('prg',           'reconstruct3D')
                call cline_rec%set('projfile',      rec_project)
                call cline_rec%set('mkdir',         'no')
                call cline_rec%set('oritype',       'ptcl3D')
                call cline_rec%set('objfun',        'cc')
                call cline_rec%set('ml_reg',        'no')
                call cline_rec%set('unfil_unmasked','yes')
                call cline_rec%set('cc_emit_sigma', 'no')
                call cline_rec%set('rec_backend',   'gridding')
                call cline_rec%set('filt_mode',     'none')
                call cline_rec%set('automsk',       'no')
                call cline_rec%set('postprocess',   'no')
                call cline_rec%set('m_estimator',   'no')
                call cline_rec%set('update_frac',   1.)
                call cline_rec%set('trail_rec',     'no')
                call cline_rec%set('pgrp',    trim(params%pgrp))
                call cline_rec%set('mskdiam', params%mskdiam)
                call cline_rec%set('trs',     5.)
                call cline_rec%set('nthr',    params%nthr)
                call cline_rec%set('nstates', nstates)
                if( params%nparts > 1 ) call cline_rec%set('nparts', params%nparts)
                if( cline%defined('qsys_name') ) call cline_rec%set('qsys_name', trim(params%qsys_name))
                call simple_chdir(rec_dir)
                CWD_GLOB = rec_dir%to_char()
                call run_rec3D(cline_rec)
                ! register artefacts generated by the reconstruction
                call spproj%read_segment('out',rec_project)
                nregistered = 0
                do state = 1, nstates
                    volname = refine3D_state_vol_fname(state)
                    fscname = refine3D_fsc_fname(state)
                    if( .not.file_exists(volname) .or. .not.file_exists(fscname) )then
                        if( state_pops(state) > 0 ) THROW_HARD('Missing fraction map/FSC; fraction folder retained')
                        cycle
                    endif
                    fsc_curve = file2rarr(fscname)
                    if( size(fsc_curve) /= size(resolution) ) THROW_HARD('Fraction FSC shell count mismatch')
                    if( .not.all(ieee_is_finite(fsc_curve)) ) THROW_HARD('Nonfinite fraction FSC values')
                    fsc_matrix(:,i,state) = fsc_curve
                    call spproj%add_vol2os_out(simple_abspath(volname), params%smpd, state, 'vol')
                    call spproj%add_fsc2os_out(simple_abspath(fscname), state, params%box)
                    nregistered = nregistered + 1
                enddo
                if( nregistered > 0 ) call spproj%write_segment_inside('out', rec_project)
                call spproj%kill
                call simple_chdir(cwd_saved)
                CWD_GLOB = cwd_saved%to_char()
                nthr_glob = saved_nthr_glob
                !$ call omp_set_num_threads(saved_threads)
                call remove_fraction_partials(cwd_saved, frame_ranges(:,i), rec_dir, nstates, max(1,params%nparts))
            else
                THROW_WARN('No particles selected for fraction reconstruction')
            endif
            call cline_rec%kill
        enddo
        do state = 1, nstates
            fscname = filepath(cwd_saved, 'fsc_fractions_state'//int2str_pad(state,2))
            call write_fraction_fsc_matrix(fscname//'.tsv', reciprocal_frequency, frame_ranges, fsc_matrix(:,:,state))
            call plot_fraction_fsc(fscname, resolution, reciprocal_frequency, frame_ranges, fsc_matrix(:,:,state), state)
        enddo
        if( l_manifest ) l_manifest = dose_manifest%write(manifest_fname, manifest_message)
        if( l_manifest )then
            write(logfhandle,'(A)') '>>> DOSE SERIES MANIFEST: '//manifest_fname%to_char()
        else
            write(logfhandle,'(A)') '>>> DOSE SERIES MANIFEST NOT GENERATED: '//manifest_message%to_char()
        endif
        call dose_manifest%kill
    end subroutine reconstruct_motion_fractions

    !> Read the loaded interval project and its contributing saved models; prepare metadata without project I/O.
    logical function prepare_fraction_dose_manifest(params, spproj, pinds, frame_ranges, manifest, rows, halves, message) result(ok)
        use simple_dose_series_manifest,   only: dose_series_manifest
        use simple_refine3D_fnames,        only: refine3D_state_halfvol_fname
        class(parameters),   target, intent(in)    :: params
        class(sp_project),           intent(inout) :: spproj
        integer,                     intent(in)    :: pinds(:), frame_ranges(:,:)
        class(dose_series_manifest), intent(inout) :: manifest
        integer,        allocatable, intent(out)   :: rows(:), halves(:)
        type(string),                intent(out)   :: message
        type(motion_model)            :: model
        type(string)                  :: model_name, digest_text, digest, rec_dir, even_name, odd_name
        integer,          allocatable :: mic2stk(:), stk2mic(:)
        logical,          allocatable :: contributing(:)
        character(len=:), allocatable :: partition_rows
        real(dp) :: dose, voltage, current_dose, current_voltage, dose_lo, dose_hi
        integer  :: nmics, nstks, imic, istk, i, neven, nodd, last_fraction, nwithptcls
        ok = .false.
        message = ''
        call manifest%kill
        call fraction_state1_rows(spproj, pinds, rows, halves)
        if( size(rows) == 0 )then
            message = 'No selected particles contribute to state 1'
            return
        endif
        if( any(halves < 0) .or. any(halves > 1) )then
            message = 'State-1 particles have invalid half assignments'
            return
        endif
        neven = count(halves == 0)
        nodd  = count(halves == 1)
        if( neven == 0 .or. nodd == 0 )then
            message = 'State 1 does not have particles in both halves'
            return
        endif
        if( size(frame_ranges,1) /= 2 .or. size(frame_ranges,2) == 0 )then
            message = 'No valid fraction intervals are available for the dose series'
            return
        endif
        last_fraction = maxval(frame_ranges(2,:))
        nmics = spproj%os_mic%get_noris()
        nstks = spproj%os_stk%get_noris()
        if( nmics == 0 .or. nstks == 0 .or. nmics < nstks )then
            message = 'State-1 particle stacks cannot be mapped to saved movie metadata'
            return
        endif
        if( nmics > nstks )then
            nwithptcls = count(spproj%os_mic%get_all('nptcls') > 0.)
            if( nwithptcls /= nstks )then
                message = 'Micrograph particle counts do not define a unique stack mapping'
                return
            endif
        endif
        call spproj%get_mic2stk_inds(mic2stk, stk2mic)
        allocate(contributing(nmics), source=.false.)
        do i = 1, size(rows)
            istk = spproj%os_ptcl3D%get_int(rows(i), 'stkind')
            if( istk < 1 .or. istk > nstks )then
                message = 'A selected state-1 particle has no valid stack index'
                return
            endif
            imic = stk2mic(istk)
            if( imic < 1 .or. imic > nmics )then
                message = 'A selected state-1 stack has no contributing micrograph'
                return
            endif
            contributing(imic) = .true.
        enddo
        dose = 0._dp
        voltage = 0._dp
        do imic = 1, nmics
            if( .not. contributing(imic) ) cycle
            if( .not. spproj%os_mic%isthere(imic, 'mcmodel') )then
                message = 'No saved motion model for contributing micrograph '//int2str(imic)
                return
            endif
            model_name = spproj%os_mic%get_str(imic, 'mcmodel')
            if( .not. file_exists(model_name) )then
                message = 'Saved motion model is missing: '//model_name%to_char()
                return
            endif
            call model%read(model_name, params)
            current_dose = real(model%dose_per_frame,dp)
            current_voltage = real(model%voltage,dp)
            if( .not. ieee_is_finite(current_dose) .or. .not. ieee_is_finite(current_voltage) )then
                message = 'Nonfinite saved dose or voltage for contributing micrograph '//int2str(imic)
            else if( current_dose <= 0._dp )then
                message = 'Saved dose per decoded fraction is not positive for contributing micrograph '//int2str(imic)
            else if( minval(abs(current_voltage - [100._dp,200._dp,300._dp])) > DOSE_VOLTAGE_TOL )then
                message = 'Unsupported saved voltage for contributing micrograph '//int2str(imic)
            else if( model%nframes < last_fraction )then
                message = 'Saved motion model does not cover every requested fraction for micrograph '//int2str(imic)
            else if( dose > 0._dp )then
                if( abs(current_dose-dose) > max(DOSE_FRACTION_ABSTOL, DOSE_FRACTION_RELTOL*max(dose,current_dose)) )then
                    message = 'Saved dose per fraction differs between contributing state-1 micrographs'
                else if( abs(current_voltage-voltage) > DOSE_VOLTAGE_TOL )then
                    message = 'Saved voltage differs between contributing state-1 micrographs'
                endif
            else
                dose = current_dose
                voltage = current_voltage
            endif
            call model%kill
            if( message%strlen_trim() > 0 ) return
        enddo
        allocate(character(len=DOSE_PARTITION_ROW_WIDTH*size(rows)) :: partition_rows)
        do i = 1, size(rows)
            write(partition_rows((i-1)*DOSE_PARTITION_ROW_WIDTH+1:i*DOSE_PARTITION_ROW_WIDTH),'(I12,":",I1,";")') rows(i), halves(i)
        enddo
        digest_text = partition_rows
        digest = digest_text%to_fnv1a_hash64()
        call manifest%new(params%box, real(params%smpd,dp), voltage, dose, size(rows), neven, nodd, digest)
        do i = 1, size(frame_ranges,2)
            rec_dir = 'rec_frames_'//int2str_pad(frame_ranges(1,i),2)//'-'//int2str_pad(frame_ranges(2,i),2)
            even_name = filepath(rec_dir, refine3D_state_halfvol_fname(1, 'even', unfil=.true.))
            odd_name  = filepath(rec_dir, refine3D_state_halfvol_fname(1, 'odd', unfil=.true.))
            dose_lo = real(frame_ranges(1,i)-1,dp) * dose
            dose_hi = real(frame_ranges(2,i),dp) * dose
            call manifest%add_window(frame_ranges(1,i), frame_ranges(2,i), dose_lo, dose_hi, even_name, odd_name)
        enddo
        ok = .true.
    end function prepare_fraction_dose_manifest

    !> Select the actual state-1 reconstruction rows and their existing even/odd labels.
    subroutine fraction_state1_rows(spproj, pinds, rows, halves)
        class(sp_project),    intent(in)  :: spproj
        integer,              intent(in)  :: pinds(:)
        integer, allocatable, intent(out) :: rows(:), halves(:)
        logical, allocatable :: selected(:)
        integer :: i
        allocate(selected(size(pinds)))
        do i = 1, size(pinds)
            selected(i) = spproj%os_ptcl3D%get_state(pinds(i)) == 1
        enddo
        rows = pack(pinds, selected)
        allocate(halves(size(rows)))
        do i = 1, size(rows)
            halves(i) = spproj%os_ptcl3D%get_eo(rows(i))
        enddo
    end subroutine fraction_state1_rows

    !> Remove a validated interval's partial reconstruction volumes; retain its stacks and projects.
    subroutine remove_fraction_partials(cwd, frames, rec_dir, nstates, nparts)
        use simple_syslib, only: simple_rmfile
        use simple_refine3D_fnames, only: refine3D_partial_rec_fname, refine3D_partial_rho_fname
        class(string), intent(in) :: cwd, rec_dir
        integer,       intent(in) :: frames(2), nstates, nparts
        character(len=4), parameter :: HALVES(2)=[character(len=4) :: 'even','odd']
        type(string),   allocatable :: partials(:)
        type(string) :: root, reconstruction_dir, expected_path, resolved
        integer      :: j, state, part, half, numlen, status, nremoved
        if( nstates<1 .or. nparts<1 ) THROW_HARD('Invalid fraction cleanup inputs')
        if( frames(1)<1 .or. frames(2)<frames(1) ) THROW_HARD('Invalid fraction cleanup interval')
        root=simple_abspath(cwd,status=status)
        if( status/=0 ) THROW_HARD('Cannot resolve fraction cleanup root')
        expected_path=filepath(root,'rec_frames_'//int2str_pad(frames(1),2)//'-'//int2str_pad(frames(2),2))
        reconstruction_dir=simple_abspath(rec_dir,status=status)
        if( status/=0 ) THROW_HARD('Cannot resolve reconstruction directory for cleanup')
        if( reconstruction_dir/=expected_path ) THROW_HARD('Unexpected reconstruction directory for cleanup')
        allocate(partials(4*nstates*nparts))
        numlen=len(int2str(nparts))
        j=0
        do state=1,nstates
            do part=1,nparts
                do half=1,2
                    j=j+1
                    partials(j)=filepath(reconstruction_dir,refine3D_partial_rec_fname(state,part,numlen,trim(HALVES(half))))
                    j=j+1
                    partials(j)=filepath(reconstruction_dir,refine3D_partial_rho_fname(state,part,numlen,trim(HALVES(half))))
                enddo
            enddo
        enddo
        do j=1,size(partials)
            if( .not.file_exists(partials(j)) ) cycle
            if( file_exists(partials(j)//'/') ) THROW_HARD('Partial reconstruction cleanup target is a directory')
            resolved=simple_abspath(partials(j),status=status)
            if( status/=0 ) THROW_HARD('Cannot resolve partial reconstruction for cleanup')
            if( resolved/=partials(j) ) THROW_HARD('Refusing redirected partial reconstruction cleanup')
        enddo
        nremoved=0
        do j=1,size(partials)
            if( file_exists(partials(j)) ) nremoved=nremoved+1
            call simple_rmfile(partials(j),status=status)
            if( status/=0 ) THROW_HARD('Failed to remove partial reconstruction')
        enddo
        write(logfhandle,'(A,I0,A)') '>>> REMOVED ', nremoved, ' PARTIAL RECONSTRUCTION FILES FROM '//reconstruction_dir%to_char()
    end subroutine remove_fraction_partials

    !> Write raw FSC columns on their native shell grid.
    subroutine write_fraction_fsc_matrix(fname, reciprocal_frequency, frame_ranges, fsc_matrix)
        class(string), intent(in) :: fname
        real(dp),      intent(in) :: reciprocal_frequency(:)
        real,          intent(in) :: fsc_matrix(:,:)
        integer,       intent(in) :: frame_ranges(:,:)
        character, parameter :: TAB = achar(9)
        integer :: funit, i, j, io_stat
        call fopen(funit, file=fname, status='replace', action='write', form='formatted')
        write(funit,'(A)',advance='no',iostat=io_stat) 'shell'//TAB//'resolution_A'//TAB//'reciprocal_frequency_A^-1'
        call fileiochk('Writing fraction FSC header', io_stat)
        do i = 1, size(fsc_matrix,2)
            write(funit,'(A)',advance='no',iostat=io_stat) TAB//'FSC_fractions_'// &
                &int2str(frame_ranges(1,i))//'-'//int2str(frame_ranges(2,i))
            call fileiochk('Writing fraction FSC intervals', io_stat)
        enddo
        write(funit,'(A)',iostat=io_stat) ''
        call fileiochk('Writing fraction FSC header newline', io_stat)
        do j = 1, size(reciprocal_frequency)
            write(funit,'(I0,*(A,ES24.16E3))',iostat=io_stat) j, TAB, 1._dp/reciprocal_frequency(j), &
                &TAB, reciprocal_frequency(j), (TAB, real(fsc_matrix(j,i),dp), i=1,size(fsc_matrix,2))
            call fileiochk('Writing fraction FSC row', io_stat)
        enddo
        call fclose(funit)
        write(logfhandle,'(A)') '>>> WROTE FRACTION FSC MATRIX: '//fname%to_char()
    end subroutine write_fraction_fsc_matrix

    !> Render raw interval FSC curves and their surface with native gnuplot.
    subroutine plot_fraction_fsc(prefix, resolution, reciprocal_frequency, frame_ranges, fsc_matrix, state)
        class(string), intent(in) :: prefix
        real,          intent(in) :: resolution(:), fsc_matrix(:,:)
        real(dp),      intent(in) :: reciprocal_frequency(:)
        integer,       intent(in) :: frame_ranges(:,:), state
        character, parameter :: TAB = achar(9)
        character(len=:), allocatable :: base, gridfile, gpfile, plotlog, line, interval
        real(dp) :: midpoint, xlo, xhi, ylo, yhi, flo, fhi
        integer :: unit, ios, close_status, i, j, stride, cmdstatus, exitstatus, bytes
        logical :: have_values, have_cells, exists, rendered
        if( size(resolution)==0 .or. size(fsc_matrix,2)==0 ) return
        base = prefix%to_char()
        gridfile = base//'_surface.tsv'
        gpfile = base//'.gp'
        plotlog = base//'_gnuplot.log'
        open(newunit=unit,file=gridfile,status='replace',action='write',iostat=ios)
        if( ios/=0 )then
            THROW_WARN('Cannot create raw FSC surface table; matrix TSV retained')
            return
        endif
        write(unit,'(A)',iostat=ios) '# reciprocal_frequency_A^-1'//TAB//'fraction_midpoint'//TAB//'FSC'//TAB// &
            &'fraction_from'//TAB//'fraction_to'//TAB//'resolution_A'
        intervals: do i=1,size(fsc_matrix,2)
            if( ios/=0 ) exit
            midpoint = (real(frame_ranges(1,i),dp)+real(frame_ranges(2,i),dp))/2._dp
            do j=1,size(resolution)
                write(unit,'(ES24.16E3,A,F12.3,A,ES16.8E3,2(A,I0),A,ES16.8E3)',iostat=ios) &
                    &reciprocal_frequency(j), TAB, midpoint, TAB, fsc_matrix(j,i), TAB, frame_ranges(1,i), &
                    &TAB, frame_ranges(2,i), TAB, resolution(j)
                if( ios/=0 ) exit intervals
            enddo
            write(unit,'(A)',iostat=ios) ''
        enddo intervals
        close(unit,iostat=close_status)
        if( ios/=0 .or. close_status/=0 )then
            THROW_WARN('Cannot finish raw FSC surface table; matrix TSV retained')
            return
        endif
        have_values = any(ieee_is_finite(fsc_matrix))
        have_cells = .false.
        do i=1,size(fsc_matrix,2)-1
            do j=1,size(fsc_matrix,1)-1
                if( all(ieee_is_finite(fsc_matrix(j:j+1,i:i+1))) ) have_cells=.true.
            enddo
        enddo
        flo=-0.2_dp; fhi=1.2_dp
        if( have_values )then
            flo=min(flo,real(minval(fsc_matrix,mask=ieee_is_finite(fsc_matrix)),dp))
            fhi=max(fhi,real(maxval(fsc_matrix,mask=ieee_is_finite(fsc_matrix)),dp))
        endif
        xlo=minval(reciprocal_frequency); xhi=maxval(reciprocal_frequency)
        if( xhi<=xlo )then
            xlo=max(0._dp,xlo-0.001_dp); xhi=xhi+0.001_dp
        endif
        ylo=(real(frame_ranges(1,1),dp)+real(frame_ranges(2,1),dp))/2._dp
        i=size(frame_ranges,2)
        yhi=(real(frame_ranges(1,i),dp)+real(frame_ranges(2,i),dp))/2._dp
        if( yhi<=ylo )then
            ylo=ylo-0.5_dp; yhi=yhi+0.5_dp
        endif
        open(newunit=unit,file=gpfile,status='replace',action='write',iostat=ios)
        if( ios/=0 )then
            THROW_WARN('Cannot create raw FSC plotting commands; TSVs retained')
            return
        endif
        call gp('set encoding utf8')
        call gp('set palette defined (0 "#440154", 0.25 "#3b528b", 0.5 "#21918c", 0.75 "#5ec962", 1 "#fde725")')
        call gp('set terminal pngcairo size 1800,850 enhanced font "Sans,11"')
        call gp('set output '//gnu_quote(base//'_diagnostics.png'))
        call gp('set multiplot layout 1,2 title "Raw FSC observations: state '//int2str(state)//'"')
        call gp('unset pm3d; unset colorbox; unset label; unset zlabel; set ytics autofreq offset 0,0; set xtics offset 0,0')
        call gp('set title "FSC curves of fraction intervals"')
        call gp('set xlabel "Reciprocal frequency (Angstrom^{-1})" offset 0,0; set ylabel "FSC" offset 0,0')
        call gp('set xrange ['//number(xlo)//':'//number(xhi)//']; set yrange ['//number(flo)//':'//number(fhi)//']')
        call gp('set key font ",7" maxcols 2')
        line='plot '
        do i=1,size(fsc_matrix,2)
            if( i>1 ) line=line//', '
            interval=int2str(frame_ranges(1,i))//'-'//int2str(frame_ranges(2,i))
            line=line//gnu_quote(base//'.tsv')//' using 3:'//int2str(i+3)
            if( size(resolution)==1 )then
                line=line//' with points pt 7'
            else
                line=line//' with lines lw 1'
            endif
            line=line//' title '//gnu_quote('Fraction interval: '//interval)
        enddo
        if( have_values )then
            call gp(line)
        else
            call gp('set label 20 "No FSC observations available" at graph 0.1,0.5')
            call gp('plot NaN notitle; unset label 20')
        endif
        call gp('unset key; set colorbox; set view 60,35; set xyplane relative 0; set xtics offset -1,0')
        call gp('set title "FSC by reciprocal frequency and fraction interval"; set zlabel "FSC"')
        call gp('set xlabel "Reciprocal frequency (Angstrom^{-1})" offset 0,-1; set ylabel "Fraction interval" offset 1,-1')
        call gp('set yrange ['//number(ylo)//':'//number(yhi)//']')
        call gp('set zrange ['//number(flo)//':'//number(fhi)//']; set cbrange ['//number(flo)//':'//number(fhi)//']')
        line='set ytics ('
        stride=max(1,(size(frame_ranges,2)+7)/8)
        do i=1,size(frame_ranges,2)
            if( mod(i-1,stride)/=0 .and. i/=size(frame_ranges,2) ) cycle
            if( i>1 .and. i<size(frame_ranges,2) .and. i+stride>size(frame_ranges,2) ) cycle
            if( i>1 ) line=line//', '
            midpoint=(real(frame_ranges(1,i),dp)+real(frame_ranges(2,i),dp))/2._dp
            interval=int2str(frame_ranges(1,i))//'-'//int2str(frame_ranges(2,i))
            line=line//gnu_quote(interval)//' '//number(midpoint)
        enddo
        call gp(line//') offset 3,0')
        line='splot '//gnu_quote(gridfile)//' using 1:2:3'
        if( have_cells )then
            call gp('set pm3d depthorder corners2color mean')
            call gp(line//' with pm3d notitle')
        else if( have_values )then
            call gp('set label 20 "Insufficient grid for a surface: showing points" at graph 0,1 font ",8"')
            call gp(line//' with points pt 7 palette notitle')
        else
            call gp('set label 20 "No FSC observations available" at graph 0.1,0.5')
            call gp('splot NaN notitle')
        endif
        call gp('unset label; unset multiplot; unset output')
        close(unit,iostat=close_status)
        if( ios/=0 .or. close_status/=0 )then
            THROW_WARN('Cannot finish fraction plotting commands; TSVs retained')
            return
        endif
        if( file_exists(base//'_diagnostics.png') ) call del_file(base//'_diagnostics.png')
        exitstatus=1; cmdstatus=0
        call execute_command_line('gnuplot '//shell_quote(gpfile)//' > '//shell_quote(plotlog)//' 2>&1', &
            &wait=.true.,exitstat=exitstatus,cmdstat=cmdstatus)
        rendered=cmdstatus==0 .and. exitstatus==0
        inquire(file=base//'_diagnostics.png',exist=exists,size=bytes,iostat=ios)
        if( ios/=0 .or. .not.exists )then
            rendered=.false.
        else if( bytes<=0 )then
            rendered=.false.
        endif
        if( rendered )then
            write(logfhandle,'(A)') '>>> WROTE FSC DIAGNOSTICS: '//base//'_diagnostics.png'
        else
            THROW_WARN('Fraction diagnostic rendering failed; TSVs retained, see '//plotlog)
        endif
      contains

        subroutine gp(command)
            character(len=*), intent(in) :: command
            if( ios==0 ) write(unit,'(A)',iostat=ios) command
        end subroutine gp

    end subroutine plot_fraction_fsc

    function number(value) result(text)
        real(dp), intent(in) :: value
        character(len=24) :: text
        write(text,'(ES24.16E3)') value
    end function number

    function gnu_quote(value) result(text)
        character(len=*), intent(in) :: value
        character(len=:), allocatable :: text
        integer :: n
        text='"'
        do n=1,len_trim(value)
            select case(iachar(value(n:n)))
            case(34,92,96)
                text=text//achar(92)//value(n:n)
            case(10)
                text=text//achar(92)//'n'
            case(13)
                text=text//achar(92)//'r'
            case(9)
                text=text//achar(92)//'t'
            case default
                text=text//value(n:n)
            end select
        enddo
        text=text//'"'
    end function gnu_quote

    function shell_quote(value) result(text)
        character(len=*), intent(in) :: value
        character(len=:), allocatable :: text
        integer :: n
        text=achar(39)
        do n=1,len_trim(value)
            if( iachar(value(n:n))==39 )then
                text=text//achar(39)//achar(34)//achar(39)//achar(34)//achar(39)
            else
                text=text//value(n:n)
            endif
        enddo
        text=text//achar(39)
    end function shell_quote

    subroutine distr_finalize_run(self, params, cline)
        class(refine_motion_model_distr_strategy), intent(inout) :: self
        type(parameters),                          intent(in)    :: params
        class(cmdline),                            intent(inout) :: cline
        ! No-op
    end subroutine distr_finalize_run

    subroutine distr_cleanup(self, params, cline)
        class(refine_motion_model_distr_strategy), intent(inout) :: self
        type(parameters),                          intent(in)    :: params
        class(cmdline),                            intent(inout) :: cline
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

end module simple_refine_motion_model_strategy
