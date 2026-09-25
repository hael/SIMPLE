!@descr Strategy pattern for refine_motion_model
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
    integer                          :: nmics_tot = 0
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
        if( cline%defined('osmpd') )then
            if( .not.cline%defined('box') ) THROW_HARD('BOX must be defined with OSMPD!')
        endif
    end subroutine validate_refine_motion_model_cline

    subroutine set_refine_motion_model_defaults(cline)
        class(cmdline), intent(inout) :: cline
        call cline%set('oritype',     'ptcl3D')
        if( .not. cline%defined('mkdir')          ) call cline%set('mkdir',          'yes')
        if( .not. cline%defined('pcontrast')      ) call cline%set('pcontrast',    'black')
        if( .not. cline%defined('backgr_subtr')   ) call cline%set('backgr_subtr',    'no')
        if( .not. cline%defined('wfloat16')       ) call cline%set('wfloat16',        'no')
    end subroutine set_refine_motion_model_defaults

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
        use simple_particle_extractor,  only: ptcl_extractor
        class(refine_motion_model_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        class(cmdline),                  intent(inout) :: cline
        type(sp_project)              :: spproj, spproj_in
        type(builder)                 :: build
        type(image)                   :: micrograph
        type(ori)                     :: o_mic, o_stk
        type(ctfparams)               :: ctfparms
        type(stack_io)                :: stkio_w
        type(ptcl_extractor)          :: extractor
        type(string)                  :: mic_name, imgkind, ext
        logical,          allocatable :: mic_mask(:), ptcl_mask(:)
        integer,          allocatable :: mic2stk_inds(:), boxcoords(:,:), ptcl_inds(:), mic_dims(:,:)
        type(string)                  :: stack
        real    :: prev_shift(2),shift2d(2),shift3d(2),prev_shift_sc(2), translation(2), prev_center_sc(2)
        real    :: stk_min, stk_max, stk_mean, stk_sdev
        integer :: prev_pos(2), new_pos(2), ishift(2), ldim(3), prev_center(2), new_center(2)
        integer :: i, nframes, imic, iptcl, nmics, prev_box, box_foo, cnt, nmics_tot, stk_ind
        integer :: fromp, top, istk, nptcls2extract, nptcls
        params%msk = RADFRAC_NORM_EXTRACT * real(params%box/2)
        call spproj_in%read_segment('mic', params%projfile)
        nmics_tot = spproj_in%os_mic%get_noris()
        if( spproj_in%get_nintgs() == 0 ) THROW_HARD('No micrograph to process!')
        if( .not.spproj_in%os_mic%isthere('mcmodel') ) THROW_HARD('No Motion model to process!')
        call spproj_in%read_segment('stk', params%projfile)
        if( spproj_in%get_nstks() == 0 ) THROW_HARD('No particles extracted!')
        box_foo  = 0
        prev_box = 0
        ldim     = 0
        allocate(mic2stk_inds(nmics_tot), source=0)
        allocate(mic_mask(nmics_tot),     source=.false.)
        stk_ind = 0
        do imic = 1, nmics_tot
            if( imic > params%top ) exit
            call spproj_in%os_mic%get_ori(imic, o_mic)
            if( o_mic%isthere('state') )then
                if( o_mic%get_state() == 0 )cycle
            endif
            if( .not. o_mic%isthere('imgkind') )cycle
            if( .not. o_mic%isthere('intg')    )cycle
            call o_mic%getter('imgkind', imgkind)
            if( imgkind.ne.'mic') cycle
            do istk = stk_ind, spproj_in%os_stk%get_noris()
                stk_ind = stk_ind + 1
                if( spproj_in%os_stk%isthere(stk_ind,'state') )then
                    if( spproj_in%os_stk%get_state(stk_ind) == 1 ) exit
                else
                    exit
                endif
            enddo
            if( imic>=params%fromp .and. imic<=params%top )then
                mic_mask(imic) = .true.
                mic2stk_inds(imic) = stk_ind
            endif
        enddo
        nmics = count(mic_mask)
        if( nmics > 0 )then
            call build%build_general_tbox(params, cline, do3d=.false.)
            allocate(mic_dims(3,nmics_tot),source=0)
            do imic = 1, nmics_tot
                if( .not.mic_mask(imic) )cycle
                call spproj_in%os_mic%get_ori(imic, o_mic)
                call o_mic%getter('intg', mic_name)
                if( .not.file_exists(mic_name) )cycle
                call find_ldim_nptcls(mic_name, ldim, nframes )
                if( nframes > 1 ) THROW_HARD('multi-frame extraction not supported; exec_refine_motion_model')
                mic_dims(:,imic) = [ldim(1),ldim(2),1]
                stk_ind = mic2stk_inds(imic)
                call spproj_in%os_stk%get_ori(stk_ind, o_stk)
                box_foo = o_stk%get_int('box')
                if( prev_box == 0 ) prev_box = box_foo
                if( prev_box /= box_foo ) THROW_HARD('Inconsistent box size; exec_refine_motion_model')
            enddo
            if( .not.cline%defined('box') ) params%box = prev_box
            if( is_odd(params%box) ) THROW_HARD('Box size must be of even dimension! exec_refine_motion_model')
            write(logfhandle,'(A)')'>>> EXTRACTING... '
            call spproj_in%read_segment('ptcl2D', params%projfile)
            call spproj_in%read_segment('ptcl3D', params%projfile)
            allocate(ptcl_mask(spproj_in%os_ptcl2D%get_noris()),source=.false.)
            ldim = 0
            do imic = params%fromp, params%top
                if( .not.mic_mask(imic) ) cycle
                call spproj_in%os_mic%get_ori(imic, o_mic)
                call o_mic%getter('intg', mic_name)
                ctfparms = o_mic%get_ctfvars()
                if( any(ldim /= mic_dims(:,imic)) )then
                    ldim = mic_dims(:,imic)
                    call micrograph%new(ldim, ctfparms%smpd)
                endif
                stk_ind = mic2stk_inds(imic)
                call spproj_in%os_stk%get_ori(stk_ind, o_stk)
                prev_box = o_stk%get_int('box')
                fromp    = o_stk%get_fromp()
                top      = o_stk%get_top()
                ext      = fname2ext(basename(mic_name))
                stack    = string(EXTRACT_STK_FBODY)//get_fbody(basename(mic_name), ext)//STK_EXT
                if( allocated(boxcoords) ) deallocate(boxcoords)
                allocate(boxcoords(2,fromp:top),source=0)
                !$omp parallel do default(shared) proc_bind(close) schedule(static)&
                !$omp private(iptcl,prev_pos,prev_shift,prev_center,prev_center_sc,prev_shift_sc,new_center)&
                !$omp private(new_pos,translation,shift2d,shift3d,ishift)
                do iptcl = fromp, top
                    if( spproj_in%os_ptcl2D%get_state(iptcl) == 0 ) cycle
                    if( spproj_in%os_ptcl3D%get_state(iptcl) == 0 ) cycle
                    call spproj_in%get_boxcoords(iptcl, prev_pos)
                    prev_shift = spproj_in%os_ptcl3D%get_2Dshift(iptcl)
                    ishift      = nint(prev_shift)
                    new_pos     = prev_pos - ishift
                    translation = -real(ishift)
                    if( prev_box /= params%box ) new_pos = new_pos + (prev_box-params%box)/2
                    ptcl_mask(iptcl) = box_inside(ldim, new_pos, params%box)
                    if( ptcl_mask(iptcl) )then
                        shift2d = spproj_in%os_ptcl2D%get_2Dshift(iptcl) + translation
                        shift3d = prev_shift                             + translation
                    endif
                    if( ptcl_mask(iptcl) )then
                        call spproj_in%set_boxcoords(iptcl, new_pos)
                        call spproj_in%os_ptcl2D%set_shift(iptcl, shift2d)
                        call spproj_in%os_ptcl3D%set_shift(iptcl, shift3d)
                    else
                        call spproj_in%os_ptcl2D%set_state(iptcl, 0)
                        call spproj_in%os_ptcl3D%set_state(iptcl, 0)
                    endif
                    boxcoords(:,iptcl) = new_pos
                enddo
                !$omp end parallel do
                nptcls2extract = count(ptcl_mask(fromp:top))
                if( nptcls2extract > 0 )then
                    if( allocated(ptcl_inds) ) deallocate(ptcl_inds)
                    allocate(ptcl_inds(nptcls2extract),source=0)
                    cnt = 0
                    do iptcl = fromp, top
                        if( .not.ptcl_mask(iptcl) ) cycle
                        cnt = cnt + 1
                        ptcl_inds(cnt) = iptcl
                        call spproj_in%os_ptcl2D%set(iptcl, 'indstk', cnt)
                        call spproj_in%os_ptcl3D%set(iptcl, 'indstk', cnt)
                    enddo
                    ptcl_inds = ptcl_inds - fromp + 1
                    call prepimgbatch_local(nptcls2extract, .false.)
                    call extractor%init_mov(o_mic, params%box, (params%pcontrast .eq. 'black'))
                    call extractor%extract_particles(ptcl_inds, boxcoords, build%imgbatch, stk_min,stk_max,stk_mean,stk_sdev)
                    call stkio_w%open(stack, params%smpd, 'write', box=params%box, wfloat16=trim(params%wfloat16).eq.'yes')
                    do i = 1, nptcls2extract
                        call stkio_w%write(i, build%imgbatch(i))
                    enddo
                    call stkio_w%close
                    call micrograph%update_header_stats(stack, [stk_min, stk_max, stk_mean, stk_sdev])
                    call spproj_in%os_stk%set(stk_ind,'stk',       simple_abspath(stack,check_exists=.false.))
                    call spproj_in%os_stk%set(stk_ind,'box',       params%box)
                    call spproj_in%os_stk%set(stk_ind,'nptcls',    nptcls2extract)
                    call spproj_in%os_stk%set(stk_ind,'nptcls_stk',nptcls2extract)
                    call spproj_in%os_mic%set(imic,   'nptcls',    nptcls2extract)
                    call spproj_in%os_mic%delete_entry(imic,'boxfile')
                else
                    call spproj_in%os_stk%set(stk_ind,'state',0)
                    call spproj_in%os_mic%set(imic,'state',0)
                    mic_mask(imic) = .false.
                    mic2stk_inds(imic) = 0
                endif
            enddo
        endif
        call micrograph%kill
        call extractor%kill
        call killimgbatch_local
        call build%kill_general_tbox
        ! OUTPUT project
        call spproj%read_non_data_segments(params%projfile)
        call spproj%projinfo%set(1,'projname', get_fbody(params%outfile,METADATA_EXT,separator=.false.))
        call spproj%projinfo%set(1,'projfile', params%outfile)
        nmics = count(mic_mask)
        call spproj%os_mic%new(nmics, is_ptcl=.false.)
        call spproj%os_stk%new(nmics, is_ptcl=.false.)
        nptcls = count(ptcl_mask)
        cnt = 0
        do imic = params%fromp, params%top
            if( .not.mic_mask(imic) )cycle
            cnt = cnt+1
            call spproj%os_mic%transfer_ori(cnt, spproj_in%os_mic, imic)
            stk_ind = mic2stk_inds(imic)
            call spproj%os_stk%transfer_ori(cnt, spproj_in%os_stk, stk_ind)
        enddo
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
        call spproj_in%kill
        call spproj%write(params%outfile)
        write(logfhandle,'(A,I8)')'>>> RE-EXTRACTED  PARTICLES: ', nptcls
        call spproj%kill
        call o_mic%kill
        call o_stk%kill

        contains

            subroutine prepimgbatch_local(batchsz, l_scaled)
                integer, intent(in) :: batchsz
                logical, intent(in) :: l_scaled
                integer :: i
                logical :: doprep
                real    :: smpd_batch
                smpd_batch = params%smpd
                if( l_scaled ) smpd_batch = params%osmpd
                doprep = .false.
                if( .not. allocated(build%imgbatch) )then
                    doprep = .true.
                else
                    if( batchsz > size(build%imgbatch) ) doprep = .true.
                    if( doprep ) call killimgbatch_local
                endif
                if( doprep )then
                    allocate(build%imgbatch(batchsz))
                    !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
                    do i = 1, batchsz
                        call build%imgbatch(i)%new([params%box, params%box, 1], smpd_batch, wthreads=.false.)
                    end do
                    !$omp end parallel do
                endif
            end subroutine prepimgbatch_local

            subroutine killimgbatch_local
                integer :: i
                if( allocated(build%imgbatch) )then
                    do i = 1, size(build%imgbatch)
                        call build%imgbatch(i)%kill
                    enddo
                    deallocate(build%imgbatch)
                endif
            end subroutine killimgbatch_local

    end subroutine inmem_execute

    subroutine inmem_finalize_run(self, params, cline)
        class(refine_motion_model_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        class(cmdline),                  intent(inout) :: cline
        call qsys_job_finished(params, string('simple_commanders_motion :: exec_refine_motion_model'))
    end subroutine inmem_finalize_run

    subroutine inmem_cleanup(self, params, cline)
        class(refine_motion_model_inmem_strategy), intent(inout) :: self
        type(parameters),                intent(in)    :: params
        class(cmdline),                  intent(inout) :: cline
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
        class(refine_motion_model_distr_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        class(cmdline),                  intent(inout) :: cline
        type(ori)      :: o_mic
        type(string)   :: mic_name, imgkind, model_name
        integer        :: imic, nmics_tot, state, ipart
        integer        :: nmics_valid
        call params%new(cline)
        call cline%set('mkdir', 'no')
        call self%spproj%read(params%projfile)
        if( self%spproj%get_nintgs() == 0 ) THROW_HARD('No integrated micrograph to process!')
        if( self%spproj%get_nstks()  == 0 ) THROW_HARD('This project file does not contain stacks!')
        if( .not.self%spproj%os_mic%isthere('mcmodel') ) THROW_HARD('No Motion model to process!')
        nmics_tot = self%spproj%os_mic%get_noris()
        self%nmics_tot = nmics_tot
        if( nmics_tot < params%nparts ) params%nparts = nmics_tot
        ! validate input
        nmics_valid = 0
        do imic = 1, nmics_tot
            call self%spproj%os_mic%get_ori(imic, o_mic)
            state = 1
            if( o_mic%isthere('state') ) state = o_mic%get_state()
            if( state == 0 ) cycle
            if( .not. o_mic%isthere('imgkind') )cycle
            if( .not. o_mic%isthere('movie')   )cycle
            if( .not. o_mic%isthere('mcmodel') )cycle
            call o_mic%getter('imgkind', imgkind)
            if( imgkind.ne.'mov') cycle
            call o_mic%getter('movie', mic_name)
            if( .not.file_exists(mic_name) )cycle
            call o_mic%getter('mcmodel', model_name)
            if( .not.file_exists(model_name) )cycle
            nmics_valid = nmics_valid + 1
        enddo
        call o_mic%kill
        if( nmics_valid == 0 )then
            THROW_WARN('No particles to re-extract! exec_refine_motion_model_distr')
            self%skip_run = .true.
            call self%spproj%kill
            call o_mic%kill
            return
        endif
        ! build explicit fromp/top for each partition
        self%parts = split_nobjs_even(nmics_tot, params%nparts)
        allocate(self%part_params(params%nparts))
        do ipart = 1, params%nparts
            call self%part_params(ipart)%new(2)
            call self%part_params(ipart)%set('fromp', int2str(self%parts(ipart,1)))
            call self%part_params(ipart)%set('top',   int2str(self%parts(ipart,2)))
        end do
        call self%qenv%new(params, params%nparts)
        call cline%gen_job_descr(self%job_descr)
    end subroutine distr_initialize

    subroutine distr_execute(self, params, cline)
        class(refine_motion_model_distr_strategy), intent(inout) :: self
        type(parameters),                intent(inout) :: params
        class(cmdline),                  intent(inout) :: cline
        if( self%skip_run ) return
        ! schedule & run
        call self%qenv%gen_scripts_and_schedule_jobs( self%job_descr, algnfbody=string(ALGN_FBODY), &
            &part_params=self%part_params, array=L_USE_SLURM_ARR, extra_params=params)
        ! ASSEMBLY - TODO
        call self%spproj%write(params%projfile)
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
