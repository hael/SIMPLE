! Strategy pattern for pick
!
! Hook-less, strategy-heavy, no separate common module.
! Integrated make_pickrefs logic (no commander_make_pickrefs dependency).
!
module simple_pick_strategy
use simple_commanders_api
use simple_parameters,     only: parameters
use simple_cmdline,        only: cmdline
use simple_qsys_env,       only: qsys_env
use simple_sp_project,     only: sp_project
use simple_image,          only: image
use simple_stack_io,       only: stack_io
use simple_default_clines, only: set_automask2D_defaults
use simple_binoris_io,     only: binwrite_oritab
implicit none

public :: pick_strategy
public :: pick_inmem_strategy
public :: pick_distr_strategy
public :: create_pick_strategy
private
#include "simple_local_flags.inc"

! --------------------------------------------------------------------
! Strategy interface
! --------------------------------------------------------------------
type, abstract :: pick_strategy
contains
    procedure(init_interface),           deferred :: initialize
    procedure(exec_interface),           deferred :: execute
    procedure(finalize_interface),       deferred :: finalize_run
    procedure(cleanup_interface),        deferred :: cleanup
    procedure(endmsg_interface),         deferred :: end_message
end type pick_strategy

! Shared-memory / distributed-worker strategy
type, extends(pick_strategy) :: pick_inmem_strategy
contains
    procedure :: initialize     => inmem_initialize
    procedure :: execute        => inmem_execute
    procedure :: finalize_run   => inmem_finalize_run
    procedure :: cleanup        => inmem_cleanup
    procedure :: end_message    => inmem_end_message
end type pick_inmem_strategy

! Distributed-master strategy
type, extends(pick_strategy) :: pick_distr_strategy
    type(sp_project) :: spproj
    type(qsys_env)   :: qenv
    type(chash)      :: job_descr
    integer          :: nmics = 0
    logical          :: templates_provided = .false.
contains
    procedure :: initialize     => distr_initialize
    procedure :: execute        => distr_execute
    procedure :: finalize_run   => distr_finalize_run
    procedure :: cleanup        => distr_cleanup
    procedure :: end_message    => distr_end_message
end type pick_distr_strategy

abstract interface
    subroutine init_interface(self, params, cline)
        import :: pick_strategy, parameters, cmdline
        class(pick_strategy), intent(inout) :: self
        type(parameters),     intent(inout) :: params
        class(cmdline),       intent(inout) :: cline
    end subroutine init_interface

    subroutine exec_interface(self, params, cline)
        import :: pick_strategy, parameters, cmdline
        class(pick_strategy), intent(inout) :: self
        type(parameters),     intent(inout) :: params
        class(cmdline),       intent(inout) :: cline
    end subroutine exec_interface

    subroutine finalize_interface(self, params, cline)
        import :: pick_strategy, parameters, cmdline
        class(pick_strategy), intent(inout) :: self
        type(parameters),     intent(in)    :: params
        class(cmdline),       intent(inout) :: cline
    end subroutine finalize_interface

    subroutine cleanup_interface(self, params, cline)
        import :: pick_strategy, parameters, cmdline
        class(pick_strategy), intent(inout) :: self
        type(parameters),     intent(in)    :: params
        class(cmdline),       intent(inout) :: cline
    end subroutine cleanup_interface

    function endmsg_interface(self) result(msg)
        import :: pick_strategy
        class(pick_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
    end function endmsg_interface
end interface

contains

    ! --------------------------------------------------------------------
    ! Factory
    ! --------------------------------------------------------------------
    function create_pick_strategy(cline) result(strategy)
        class(cmdline), intent(in) :: cline
        class(pick_strategy), allocatable :: strategy
        logical :: is_master
        ! Distributed master when nparts defined but no explicit worker range/part.
        is_master = cline%defined('nparts') .and. (.not.cline%defined('part')) &
                   .and. (.not.cline%defined('fromp')) .and. (.not.cline%defined('top'))
        if( is_master )then
            allocate(pick_distr_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DISTRIBUTED PICK (MASTER)'
        else
            allocate(pick_inmem_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> PICK (SHARED-MEMORY / WORKER)'
        endif
    end function create_pick_strategy

    ! --------------------------------------------------------------------
    ! Shared defaults + validation (private)
    ! --------------------------------------------------------------------
    subroutine set_pick_defaults(cline)
        class(cmdline), intent(inout) :: cline
        type(string) :: which_picker
        call cline%set('oritype', 'mic')
        if( .not. cline%defined('mkdir')        ) call cline%set('mkdir',       'yes')
        if( .not. cline%defined('pcontrast')    ) call cline%set('pcontrast', 'black')
        if( .not. cline%defined('oritype')      ) call cline%set('oritype',     'mic')
        if( .not. cline%defined('thres')        ) call cline%set('thres',         24.)
        if( .not. cline%defined('pick_roi')     ) call cline%set('pick_roi',    'yes')
        if( .not. cline%defined('backgr_subtr') ) call cline%set('backgr_subtr', 'no')
        if( .not. cline%defined('picker')       ) call cline%set('picker',      'new')
        if( .not. cline%defined('lp')           ) call cline%set('lp', PICK_LP_DEFAULT)
        ! ndev default depends on picker choice
        which_picker = cline%get_carg('picker')
        if( which_picker .eq. 'seg' .or. which_picker .eq. 'segdiam' )then
            if( .not. cline%defined('ndev') ) call cline%set('ndev', 1.5)
        else
            if( .not. cline%defined('ndev') ) call cline%set('ndev', 2.0)
        endif
    end subroutine set_pick_defaults

    subroutine validate_pick_mode(params, cline, templates_provided)
        type(parameters), intent(in)    :: params
        class(cmdline),   intent(inout) :: cline
        logical,          intent(in)    :: templates_provided
        select case(trim(params%picker))
            case('new')
                if( templates_provided )then
                    ! reference-based picking OK
                else if( cline%defined('moldiam') )then
                    if( cline%defined('moldiam_max') .and. cline%defined('nmoldiams') )then
                        ! multipick diameter scan OK
                    else if( cline%defined('moldiam_max') .or. cline%defined('nmoldiams') )then
                        THROW_HARD('MOLDIAM, MOLDIAM_MAX & NMOLDIAMS must be provided for determination of optimal diameter!')
                    else
                        ! single diameter picking OK
                    endif
                else if( cline%defined('multi_moldiams') )then
                    ! explicit diameter list OK
                else
                    THROW_HARD('Unsupported new picker mode')
                endif

            case('seg', 'segdiam')
                if( templates_provided ) THROW_HARD('SEGDIAM picker does not requires PICKREFS input')

            case default
                THROW_HARD('Unsupported PICKER: '//trim(params%picker))
        end select
    end subroutine validate_pick_mode

    ! --------------------------------------------------------------------
    ! Integrated make_pickrefs implementation (NO simple_end inside)
    ! --------------------------------------------------------------------
    subroutine make_pickrefs_impl(cline)
        use simple_image_msk, only: automask2D
        class(cmdline), intent(inout) :: cline
        type(parameters)         :: params
        type(stack_io)           :: stkio_r
        type(oris)               :: moldiamori
        type(image)              :: ref2D, ref2D_clip
        type(image), allocatable :: projs(:), masks(:)
        real,        allocatable :: diams(:), shifts(:,:)
        real,    parameter :: MSKDIAM2LP = 0.15, LP_LB = 30., LP_UB = 15.
        integer, parameter :: NREFS = 100
        real    :: ang, rot, lp, diam_max, maxdiam, moldiam, mskdiam
        integer :: nrots, iref, irot, ldim_clip(3), ldim(3), ncavgs, icavg
        integer :: cnt, norefs, box_for_pick, box_for_extract
        ! error check
        if( cline%defined('vol1') ) THROW_HARD('vol1 input no longer supported, use prg=reproject to generate 20 2D references')
        if( .not.cline%defined('pickrefs') ) THROW_HARD('PICKREFS must be informed!')
        ! set defaults for automasking etc.
        call set_automask2D_defaults(cline)
        call cline%set('oritype', 'mic')
        if( .not.cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        ! parse parameters
        call params%new(cline)
        if( params%stream.eq.'yes' ) THROW_HARD('not a streaming application')
        ! read selected cavgs
        call find_ldim_nptcls(params%pickrefs, ldim, ncavgs)
        ldim(3) = 1
        params%msk = real(ldim(1)/2) - COSMSKHALFWIDTH ! for automasking
        allocate( projs(ncavgs), masks(ncavgs) )
        call stkio_r%open(params%pickrefs, params%smpd, 'read', bufsz=ncavgs)
        do icavg = 1, ncavgs
            call projs(icavg)%new(ldim, params%smpd)
            call stkio_r%read(icavg, projs(icavg))
            call masks(icavg)%copy(projs(icavg))
        end do
        call stkio_r%close
        ! automasking
        call automask2D(params, masks, params%ngrow, nint(params%winsz), params%edge, diams, shifts)
        do icavg = 1, ncavgs
            call projs(icavg)%div_below(0., 10.)
            call projs(icavg)%mul(masks(icavg))
            call projs(icavg)%shift([shifts(icavg,1), shifts(icavg,2), 0.])
        end do
        ! estimate new box size and clip
        diam_max        = maxval(diams)
        lp              = min(max(LP_LB, MSKDIAM2LP * diam_max), LP_UB)
        box_for_pick    = min(round2even(diam_max / params%smpd + 2. * COSMSKHALFWIDTH), ldim(1))
        moldiam         = params%smpd * box_for_pick
        mskdiam         = moldiam * MSK_EXP_FAC
        maxdiam         = moldiam + moldiam * BOX_EXP_FAC
        box_for_extract = find_larger_magic_box(round2even(maxdiam / params%smpd))
        write(logfhandle,'(A,1X,I4)') 'ESTIMATED BOX SIZE: ', box_for_pick
        ! set mskdiam in cline for downstream consumers (if any)
        call cline%set('mskdiam', mskdiam)
        ! write info to file
        call moldiamori%new(1, .false.)
        call moldiamori%set(1, 'diam_max',        diam_max)
        call moldiamori%set(1, 'lp',              lp)
        call moldiamori%set(1, 'box_for_pick',    box_for_pick)
        call moldiamori%set(1, 'moldiam',         moldiam)
        call moldiamori%set(1, 'mskdiam',         mskdiam)
        call moldiamori%set(1, 'box_for_extract', box_for_extract)
        if( file_exists(STREAM_MOLDIAM) ) call del_file(STREAM_MOLDIAM)
        call moldiamori%write(1, string(STREAM_MOLDIAM))
        call moldiamori%kill
        ! expand in in-plane rotation, clip and write to file
        ldim_clip = [box_for_pick, box_for_pick, 1]
        do icavg = 1, ncavgs
            call projs(icavg)%bp(0., lp)
        end do
        nrots  = nint( real(NREFS) / real(ncavgs) )
        norefs = ncavgs
        call ref2D_clip%new([ldim_clip(1), ldim_clip(2), 1], params%smpd)
        if( nrots > 1 )then
            call ref2D%new([ldim(1), ldim(2), 1], params%smpd)
            ang = 360. / real(nrots)
            cnt = 0
            do iref = 1, norefs
                rot = 0.
                do irot = 1, nrots
                    cnt = cnt + 1
                    call projs(iref)%rtsq(rot, 0., 0., ref2D)
                    call ref2D%clip(ref2D_clip)
                    call ref2D_clip%write(string(PICKREFS_FBODY)//params%ext, cnt)
                    rot = rot + ang
                end do
            end do
        else
            ! should never happen
            do iref = 1, norefs
                call projs(iref)%clip(ref2D_clip)
                call ref2D_clip%write(string(PICKREFS_FBODY)//params%ext, iref)
            end do
        endif
        ! cleanup
        do icavg = 1, ncavgs
            call masks(icavg)%kill
            call projs(icavg)%kill
        end do
        deallocate(masks, projs)
        call ref2D%kill
        call ref2D_clip%kill
        ! marker file (kept from original behaviour)
        call simple_touch('MAKE_PICKREFS_FINISHED')
    end subroutine make_pickrefs_impl

    ! ====================================================================
    ! PICK (SHARED-MEMORY / WORKER)
    ! ====================================================================

    subroutine inmem_initialize(self, params, cline)
        class(pick_inmem_strategy), intent(inout) :: self
        type(parameters),           intent(inout) :: params
        class(cmdline),             intent(inout) :: cline
        call set_pick_defaults(cline)
        call params%new(cline)
    end subroutine inmem_initialize

    subroutine inmem_execute(self, params, cline)
        use simple_picker_iter, only: picker_iter
        class(pick_inmem_strategy), intent(inout) :: self
        type(parameters),           intent(inout) :: params
        class(cmdline),             intent(inout) :: cline
        type(sp_project)  :: spproj
        type(picker_iter) :: piter
        type(ori)         :: o
        type(string)      :: output_dir, intg_name, imgkind
        type(string)      :: boxfile, thumb_den
        integer           :: fromto(2), imic, ntot, nptcls_out, cnt, state
        
        output_dir = PATH_HERE
        ! loop range
        if( params%stream .eq. 'yes' )then
            fromto(:) = 1
        else
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto = [params%fromp, params%top]
            else
                THROW_HARD('fromp & top args need to be defined in parallel execution; exec_pick')
            endif
        endif
        ntot = fromto(2) - fromto(1) + 1
        call spproj%read(params%projfile)
        if( spproj%get_nintgs() == 0 )then
            THROW_HARD('No integrated micrograph to process!')
        endif
        cnt = 0
        do imic = fromto(1), fromto(2)
            cnt = cnt + 1
            call spproj%os_mic%get_ori(imic, o)
            state = 1
            if( o%isthere('state') ) state = o%get_state()
            if( state == 0 ) cycle
            if( o%isthere('imgkind') )then
                call o%getter('imgkind', imgkind)
                if( imgkind.ne.'mic' )cycle
                call o%getter('intg', intg_name)
                call piter%iterate(params, cline, params%smpd, intg_name, output_dir, boxfile, thumb_den, nptcls_out)
                call spproj%set_boxfile(imic, boxfile, nptcls=nptcls_out)
                if( params%nmoldiams == 1 ) call spproj%os_mic%set(imic, 'thumb_den', thumb_den)
            endif
            write(logfhandle,'(f4.0,1x,a)') 100.*(real(cnt)/real(ntot)), 'percent of the micrographs processed'
        end do
        call binwrite_oritab(params%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        call o%kill
        call spproj%kill
        call piter%kill
    end subroutine inmem_execute

    subroutine inmem_finalize_run(self, params, cline)
        class(pick_inmem_strategy), intent(inout) :: self
        type(parameters),           intent(in)    :: params
        class(cmdline),             intent(inout) :: cline
        call qsys_job_finished(params, string('simple_commanders_pick :: exec_pick'))
    end subroutine inmem_finalize_run

    subroutine inmem_cleanup(self, params, cline)
        class(pick_inmem_strategy), intent(inout) :: self
        type(parameters),           intent(in)    :: params
        class(cmdline),             intent(inout) :: cline
        ! No-op
    end subroutine inmem_cleanup

    function inmem_end_message(self) result(msg)
        class(pick_inmem_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
        msg = '**** SIMPLE_PICK NORMAL STOP ****'
    end function inmem_end_message

    ! ====================================================================
    ! DISTRIBUTED PICK (MASTER)
    ! ====================================================================

    subroutine distr_initialize(self, params, cline)
        class(pick_distr_strategy), intent(inout) :: self
        type(parameters),           intent(inout) :: params
        class(cmdline),             intent(inout) :: cline
        integer :: nmics_local
        call set_pick_defaults(cline)
        call params%new(cline)
        ! sanity check
        call self%spproj%read_segment(params%oritype, params%projfile)
        nmics_local = self%spproj%get_nintgs()
        if( nmics_local == 0 ) THROW_HARD('No micrograph to process! exec_pick_distr')
        call self%spproj%kill
        self%nmics = nmics_local
        ! set mkdir to no (avoid nested directory structure)
        call cline%set('mkdir', 'no')
        if( cline%defined('nparts') )then
            params%nparts = min(params%nparts, self%nmics)
            call cline%set('nparts', params%nparts)
        endif
        params%numlen = len(int2str(params%nparts))
        call cline%set('numlen', params%numlen)
        ! if pick_roi yes -> enforce background subtraction yes
        if( trim(params%pick_roi).eq.'yes' )then
            params%backgr_subtr = 'yes'
            call cline%set('backgr_subtr', params%backgr_subtr)
        endif
        self%templates_provided = cline%defined('pickrefs')
        call validate_pick_mode(params, cline, self%templates_provided)
        ! Ensure merge uses correct entry count
        params%nptcls = self%nmics
        call self%qenv%new(params, params%nparts)
    end subroutine distr_initialize

    subroutine distr_execute(self, params, cline)
        class(pick_distr_strategy), intent(inout) :: self
        type(parameters),           intent(inout) :: params
        class(cmdline),             intent(inout) :: cline
        type(cmdline) :: cline_make_pickrefs
        ! Prepare picking templates if requested (integrated make_pickrefs)
        if( self%templates_provided )then
            cline_make_pickrefs = cline
            call cline_make_pickrefs%set('prg',  'make_pickrefs')
            call cline_make_pickrefs%set('smpd', params%smpd)
            call cline_make_pickrefs%set('mkdir','no')
            if( trim(params%picker).eq.'old' ) call cline_make_pickrefs%set('neg','yes')
            call make_pickrefs_impl(cline_make_pickrefs)
            ! Switch pickrefs to the prepared template stack for the actual picking jobs
            call cline%set('pickrefs', PICKREFS_FBODY//params%ext%to_char())
            params%pickrefs = string(PICKREFS_FBODY)//params%ext
            write(logfhandle,'(A)') '>>> PREPARED PICKING TEMPLATES'
        endif
        call cline%gen_job_descr(self%job_descr)
        call self%qenv%gen_scripts_and_schedule_jobs(self%job_descr, algnfbody=string(ALGN_FBODY), &
            &array=L_USE_SLURM_ARR, extra_params=params)
        call self%spproj%read(params%projfile)
        call self%spproj%merge_algndocs(params%nptcls, params%nparts, 'mic', ALGN_FBODY)
        call self%spproj%kill
    end subroutine distr_execute

    subroutine distr_finalize_run(self, params, cline)
        class(pick_distr_strategy), intent(inout) :: self
        type(parameters),           intent(in)    :: params
        class(cmdline),             intent(inout) :: cline
        ! No-op
    end subroutine distr_finalize_run

    subroutine distr_cleanup(self, params, cline)
        class(pick_distr_strategy), intent(inout) :: self
        type(parameters),           intent(in)    :: params
        class(cmdline),             intent(inout) :: cline
        call qsys_cleanup(params)
        call self%qenv%kill
        call self%job_descr%kill
    end subroutine distr_cleanup

    function distr_end_message(self) result(msg)
        class(pick_distr_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
        msg = '**** SIMPLE_DISTR_PICK NORMAL STOP ****'
    end function distr_end_message

end module simple_pick_strategy
