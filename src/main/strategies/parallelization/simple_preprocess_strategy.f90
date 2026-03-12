! object-oriented strategy pattern for preprocess
!
! Goals
!   - One unified workflow selected via runtime polymorphism
!   - Shared-memory + distributed-worker execution lives in the inmem strategy
!   - Distributed-master execution lives in the distr strategy
!   - Hook-less: strategies implement the work directly
!   - No separate common module: shared helpers are private procedures below
!
module simple_preprocess_strategy
use simple_commanders_api
use simple_parameters, only: parameters
use simple_cmdline,    only: cmdline
use simple_qsys_env,   only: qsys_env
use simple_sp_project, only: sp_project
implicit none

public :: preprocess_strategy, preprocess_inmem_strategy, preprocess_distr_strategy, create_preprocess_strategy
private
#include "simple_local_flags.inc"

! --------------------------------------------------------------------
! Strategy interface
! --------------------------------------------------------------------

type, abstract :: preprocess_strategy
contains
    procedure(init_interface),     deferred :: initialize
    procedure(exec_interface),     deferred :: execute
    procedure(finalize_interface), deferred :: finalize_run
    procedure(cleanup_interface),  deferred :: cleanup
    procedure(endmsg_interface),   deferred :: end_message
end type preprocess_strategy

! Shared-memory + distributed-worker implementation
type, extends(preprocess_strategy) :: preprocess_inmem_strategy
contains
    procedure :: initialize     => inmem_initialize
    procedure :: execute        => inmem_execute
    procedure :: finalize_run   => inmem_finalize_run
    procedure :: cleanup        => inmem_cleanup
    procedure :: end_message    => inmem_end_message
end type preprocess_inmem_strategy

! Distributed-master implementation
type, extends(preprocess_strategy) :: preprocess_distr_strategy
    type(qsys_env)   :: qenv
    type(chash)      :: job_descr
    type(sp_project) :: spproj
contains
    procedure :: initialize     => distr_initialize
    procedure :: execute        => distr_execute
    procedure :: finalize_run   => distr_finalize_run
    procedure :: cleanup        => distr_cleanup
    procedure :: end_message    => distr_end_message
end type preprocess_distr_strategy

abstract interface
    subroutine init_interface(self, params, cline)
        import :: preprocess_strategy, parameters, cmdline
        class(preprocess_strategy), intent(inout) :: self
        type(parameters),           intent(inout) :: params
        class(cmdline),             intent(inout) :: cline
    end subroutine init_interface

    subroutine exec_interface(self, params, cline)
        import :: preprocess_strategy, parameters, cmdline
        class(preprocess_strategy), intent(inout) :: self
        type(parameters),           intent(inout) :: params
        class(cmdline),             intent(inout) :: cline
    end subroutine exec_interface

    subroutine finalize_interface(self, params, cline)
        import :: preprocess_strategy, parameters, cmdline
        class(preprocess_strategy), intent(inout) :: self
        type(parameters),           intent(in)    :: params
        class(cmdline),             intent(inout) :: cline
    end subroutine finalize_interface

    subroutine cleanup_interface(self, params, cline)
        import :: preprocess_strategy, parameters, cmdline
        class(preprocess_strategy), intent(inout) :: self
        type(parameters),           intent(in)    :: params
        class(cmdline),             intent(inout) :: cline
    end subroutine cleanup_interface

    function endmsg_interface(self) result(msg)
        import :: preprocess_strategy
        class(preprocess_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
    end function endmsg_interface
end interface

contains

    ! --------------------------------------------------------------------
    ! Factory
    ! --------------------------------------------------------------------

    function create_preprocess_strategy(cline) result(strategy)
        class(cmdline), intent(in) :: cline
        class(preprocess_strategy), allocatable :: strategy
        logical :: is_master
        ! Master heuristic: nparts is defined, but no explicit worker range/part.
        is_master = cline%defined('nparts') .and. (.not.cline%defined('part')) &
                   .and. (.not.cline%defined('fromp')) .and. (.not.cline%defined('top'))
        if( is_master )then
            allocate(preprocess_distr_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> DISTRIBUTED PREPROCESS (MASTER)'
        else
            allocate(preprocess_inmem_strategy :: strategy)
            if( L_VERBOSE_GLOB ) write(logfhandle,'(A)') '>>> PREPROCESS (SHARED-MEMORY / WORKER)'
        endif
    end function create_preprocess_strategy

    ! --------------------------------------------------------------------
    ! Shared defaults (kept private; no separate common module)
    ! --------------------------------------------------------------------

    subroutine set_preprocess_defaults(cline)
        class(cmdline), intent(inout) :: cline
        ! General
        call cline%set('oritype', 'mic')
        if( .not. cline%defined('stream')  ) call cline%set('stream',  'no')
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',  'yes')
        ! Motion correction
        if( .not. cline%defined('trs')            ) call cline%set('trs',               20.)
        if( .not. cline%defined('lpstart')        ) call cline%set('lpstart',            8.)
        if( .not. cline%defined('lpstop')         ) call cline%set('lpstop',             5.)
        if( .not. cline%defined('bfac')           ) call cline%set('bfac',              50.)
        if( .not. cline%defined('mcconvention')   ) call cline%set('mcconvention', 'simple')
        if( .not. cline%defined('eer_upsampling') ) call cline%set('eer_upsampling',      1)
        if( .not. cline%defined('mcpatch')        ) call cline%set('mcpatch',         'yes')
        if( .not. cline%defined('mcpatch_thres')  ) call cline%set('mcpatch_thres',   'yes')
        if( .not. cline%defined('algorithm')      ) call cline%set('algorithm',     'patch')
        ! CTF estimation
        if( .not. cline%defined('pspecsz')         ) call cline%set('pspecsz',          512)
        if( .not. cline%defined('hp_ctf_estimate') ) call cline%set('hp_ctf_estimate',  30.)
        if( .not. cline%defined('lp_ctf_estimate') ) call cline%set('lp_ctf_estimate',   5.)
        if( .not. cline%defined('dfmin')           ) call cline%set('dfmin',  DFMIN_DEFAULT)
        if( .not. cline%defined('dfmax')           ) call cline%set('dfmax',  DFMAX_DEFAULT)
        if( .not. cline%defined('ctfpatch')        ) call cline%set('ctfpatch',       'yes')
        ! Extraction
        if( .not. cline%defined('pcontrast') ) call cline%set('pcontrast', 'black')
    end subroutine set_preprocess_defaults

    ! ====================================================================
    ! PREPROCESS (SHARED-MEMORY / WORKER)
    ! ====================================================================

    subroutine inmem_initialize(self, params, cline)
        class(preprocess_inmem_strategy), intent(inout) :: self
        type(parameters),                 intent(inout) :: params
        class(cmdline),                   intent(inout) :: cline
        call set_preprocess_defaults(cline)
        call params%new(cline)
        if( params%scale_movies > 1.01 )then
            THROW_HARD('scale_movies cannot be > 1; exec_preprocess')
        endif
    end subroutine inmem_initialize

    subroutine inmem_execute(self, params, cline)
        use FoX_dom
        use simple_motion_correct_iter, only: motion_correct_iter
        use simple_ctf_estimate_iter,   only: ctf_estimate_iter
        use simple_binoris_io,          only: binwrite_oritab
        use simple_mini_stream_utils,   only: segdiampick_preprocess
        class(preprocess_inmem_strategy), intent(inout) :: self
        type(parameters),                 intent(inout) :: params
        class(cmdline),                   intent(inout) :: cline
        type(ori)                     :: o_mov
        type(ctf_estimate_iter)       :: ctfiter
        type(motion_correct_iter)     :: mciter
        type(sp_project)              :: spproj
        type(ctfparams)               :: ctfvars
        type(Node), pointer           :: xmldoc, beamshiftnode, beamshiftnodex, beamshiftnodey
        type(string) :: imgkind, moviename, fbody
        type(string) :: moviename_forctf, output_dir_motion_correct
        type(string) :: output_dir_ctf_estimate, output_dir_inipick_preproc, micname_intg
        type(string) :: eputiltgroup, str_meta
        integer      :: nmovies, fromto(2), imovie, ntot, frame_counter
        logical      :: l_del_forctf
        ! Read in movies/micrographs
        call spproj%read( params%projfile )
        if( spproj%get_nmovies()==0 .and. spproj%get_nintgs()==0 ) THROW_HARD('No movie/micrograph to process!')
        ! Output directories & naming
        output_dir_ctf_estimate        = PATH_HERE
        output_dir_motion_correct      = PATH_HERE
        output_dir_inipick_preproc     = PATH_HERE
        if( params%stream.eq.'yes' )then
            output_dir_ctf_estimate    = DIR_CTF_ESTIMATE
            output_dir_motion_correct  = DIR_MOTION_CORRECT
            output_dir_inipick_preproc = DIR_INIPICK_PREPROC
            if( cline%defined('dir') )then
                output_dir_ctf_estimate    = filepath(params%dir,output_dir_ctf_estimate)//'/'
                output_dir_motion_correct  = filepath(params%dir,output_dir_motion_correct)//'/'
                output_dir_inipick_preproc = filepath(params%dir,output_dir_inipick_preproc)//'/'
            endif
            call simple_mkdir(output_dir_ctf_estimate)
            call simple_mkdir(output_dir_motion_correct)
            if(params%ninipick > 0) then
                call simple_mkdir(output_dir_inipick_preproc)
            endif
        endif
        if( cline%defined('fbody') )then
            fbody = params%fbody
        else
            fbody = ''
        endif
        ! Range
        if( trim(params%stream).eq.'yes' )then
            ! STREAMING MODE
            fromto(:) = 1
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto(1) = params%fromp
                fromto(2) = params%top
            endif
        else
            ! DISTRIBUTED WORKER MODE
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto(1) = params%fromp
                fromto(2) = params%top
            else
                THROW_HARD('fromp & top args need to be defined in parallel execution; exec_preprocess')
            endif
        endif
        ntot = fromto(2) - fromto(1) + 1
        ! numlen (avoid using an uninitialised nmovies)
        nmovies = spproj%os_mic%get_noris()
        if( .not. cline%defined('numlen') )then
            params%numlen = len(int2str(nmovies))
        endif
        frame_counter = 0
        ! Loop over exposures (movies/micrographs)
        do imovie = fromto(1),fromto(2)
            ! Reset deletion flag per exposure (fixes latent carry-over bug)
            l_del_forctf = .false.
            ! Fetch movie orientation
            call spproj%os_mic%get_ori(imovie, o_mov)
            ! Sanity check
            if(.not.o_mov%isthere('imgkind') )cycle
            if(.not.o_mov%isthere('movie') .and. .not.o_mov%isthere('intg'))cycle
            call o_mov%getter('imgkind', imgkind)
            select case(imgkind%to_char())
                case('movie')
                    ! Motion correction
                    ctfvars = spproj%get_micparams(imovie)
                    call o_mov%getter('movie', moviename)
                    if( .not.file_exists(moviename)) cycle
                    if( cline%defined('gainref') )then
                        call mciter%iterate(params, cline, ctfvars, o_mov, fbody, frame_counter, moviename,&
                            &output_dir_motion_correct, gainref_fname=params%gainref)
                    else
                        call mciter%iterate(params, cline, ctfvars, o_mov, fbody, frame_counter, moviename,&
                            &output_dir_motion_correct)
                    endif
                    moviename_forctf = mciter%get_moviename('forctf')
                    l_del_forctf     = .true.
                case('mic')
                    ! Integrated micrograph: use directly for CTF
                    ctfvars = spproj%get_micparams(imovie)
                    call o_mov%getter('intg', moviename_forctf)
                case default
                    cycle
            end select
            ! CTF estimate
            params%hp = params%hp_ctf_estimate
            params%lp = max(params%fny, params%lp_ctf_estimate)
            call ctfiter%iterate(params, ctfvars, moviename_forctf, o_mov, output_dir_ctf_estimate, .false.)
            ! Delete temporary file after estimation
            if( l_del_forctf )then
                call o_mov%delete_entry('forctf')
                call del_file(moviename_forctf)
            endif
            ! Read XML for beamshift
            if(o_mov%isthere('meta'))then
                str_meta = o_mov%get_str('meta')
                if( file_exists(str_meta) ) then
                    xmldoc => parseFile(str_meta%to_char())
                    beamshiftnode  => item(getElementsByTagname(xmldoc, "BeamShift"),   0)
                    beamshiftnodex => item(getElementsByTagname(beamshiftnode, "a:_x"), 0)
                    beamshiftnodey => item(getElementsByTagname(beamshiftnode, "a:_y"), 0)
                    call o_mov%set("shiftx", str2real(getTextContent(beamshiftnodex)))
                    call o_mov%set("shifty", str2real(getTextContent(beamshiftnodey)))
                    call destroy(xmldoc)
                endif
            end if
            ! Assign beamtilt if EPU
            if(o_mov%isthere('intg')) then
                call o_mov%getter('intg', micname_intg)
                micname_intg = basename(micname_intg)
                if(micname_intg%to_char([1,8]) == 'FoilHole') then
                    ! EPU filename
                    eputiltgroup = micname_intg%to_char([micname_intg%substr_ind('Data_') + 5, micname_intg%strlen_trim()])
                    eputiltgroup = eputiltgroup%to_char([1,eputiltgroup%substr_ind('_') - 1])
                    call o_mov%set("tiltgrp", eputiltgroup%to_real())
                else
                    call o_mov%set("tiltgrp", 0.0)
                end if
            end if
            ! Update project
            call spproj%os_mic%set_ori(imovie, o_mov)
        end do
        ! Initial picking preprocessing
        if(params%ninipick > 0 &
        &.and. spproj%os_mic%isthere(1, 'importind') &
        &.and. spproj%os_mic%get(1, 'importind') <= params%ninipick) then
            call segdiampick_preprocess( spproj, params%pcontrast, params%moldiam_max, output_dir_inipick_preproc )
        endif
        if( trim(params%stream).eq.'yes' )then
            call spproj%write_segment_inside(params%oritype)
        else
            call binwrite_oritab(params%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        endif
        ! Cleanup
        call o_mov%kill
        call spproj%kill
    end subroutine inmem_execute

    subroutine inmem_finalize_run(self, params, cline)
        class(preprocess_inmem_strategy), intent(inout) :: self
        type(parameters),                 intent(in)    :: params
        class(cmdline),                   intent(inout) :: cline
        call qsys_job_finished(params, string('simple_commanders_preprocess :: exec_preprocess'))
    end subroutine inmem_finalize_run

    subroutine inmem_cleanup(self, params, cline)
        class(preprocess_inmem_strategy), intent(inout) :: self
        type(parameters),                 intent(in)    :: params
        class(cmdline),                   intent(inout) :: cline
        ! No-op
    end subroutine inmem_cleanup

    function inmem_end_message(self) result(msg)
        class(preprocess_inmem_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
        msg = '**** SIMPLE_PREPROCESS NORMAL STOP ****'
    end function inmem_end_message

    ! ====================================================================
    ! DISTRIBUTED PREPROCESS (MASTER)
    ! ====================================================================

    subroutine distr_initialize(self, params, cline)
        use simple_motion_correct_utils, only: flip_gain
        class(preprocess_distr_strategy), intent(inout) :: self
        type(parameters),                 intent(inout) :: params
        class(cmdline),                   intent(inout) :: cline
        integer :: nmovies
        call set_preprocess_defaults(cline)
        ! Parse parameters
        call params%new(cline)
        ! Set mkdir to no (avoid nested directory structure for workers)
        call cline%set('mkdir', 'no')
        ! Read movie segment only (consistent with other strategies; full re-read in distr_execute)
        call self%spproj%read_segment(params%oritype, params%projfile)
        nmovies = self%spproj%get_nmovies()
        if( nmovies == 0 )then
            THROW_HARD('no movie to process! exec_preprocess_distr')
        endif
        params%nptcls = nmovies
        if( params%nparts > params%nptcls ) THROW_HARD('# partitions (nparts) must be < number of entries in filetable')
        call self%spproj%kill
        ! Deal with numlen so that length matches JOB_FINISHED indicator files
        params%numlen = len(int2str(params%nparts))
        call cline%set('numlen', params%numlen)
        ! Gain reference
        call flip_gain(cline, params%gainref, params%flipgain)
        ! Setup the environment for distributed execution
        call self%qenv%new(params, params%nparts)
        ! Prepare job description
        call cline%gen_job_descr(self%job_descr)
    end subroutine distr_initialize

    subroutine distr_execute(self, params, cline)
        class(preprocess_distr_strategy), intent(inout) :: self
        type(parameters),                 intent(inout) :: params
        class(cmdline),                   intent(inout) :: cline
        ! Schedule & run
        call self%qenv%gen_scripts_and_schedule_jobs(self%job_descr, algnfbody=string(ALGN_FBODY), &
            &array=L_USE_SLURM_ARR, extra_params=params)
        ! Merge docs
        call self%spproj%read(params%projfile)
        call self%spproj%merge_algndocs(params%nptcls, params%nparts, 'mic', ALGN_FBODY)
        call self%spproj%kill
    end subroutine distr_execute

    subroutine distr_finalize_run(self, params, cline)
        class(preprocess_distr_strategy), intent(inout) :: self
        type(parameters),                 intent(in)    :: params
        class(cmdline),                   intent(inout) :: cline
        ! No-op
    end subroutine distr_finalize_run

    subroutine distr_cleanup(self, params, cline)
        class(preprocess_distr_strategy), intent(inout) :: self
        type(parameters),                 intent(in)    :: params
        class(cmdline),                   intent(inout) :: cline
        call qsys_cleanup(params)
        call self%qenv%kill
        call self%job_descr%kill
    end subroutine distr_cleanup

    function distr_end_message(self) result(msg)
        class(preprocess_distr_strategy), intent(in) :: self
        character(len=:), allocatable :: msg
        msg = '**** SIMPLE_DISTR_PREPROCESS NORMAL STOP ****'
    end function distr_end_message

end module simple_preprocess_strategy
