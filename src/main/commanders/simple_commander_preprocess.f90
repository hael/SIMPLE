! concrete commander: pre-processing routines
module simple_commander_preprocess
include 'simple_lib.f08'
use simple_binoris_io
use simple_builder,              only: builder
use simple_cmdline,              only: cmdline
use simple_parameters,           only: parameters, params_glob
use simple_commander_base,       only: commander_base
use simple_image,                only: image
use simple_sp_project,           only: sp_project
use simple_qsys_env,             only: qsys_env
use simple_stack_io,             only: stack_io
use simple_motion_correct_utils, only: flip_gain
use simple_qsys_funs
use simple_progress
implicit none

public :: preprocess_commander_distr
public :: preprocess_commander
public :: motion_correct_commander_distr
public :: motion_correct_commander
public :: gen_pspecs_and_thumbs_commander_distr
public :: gen_pspecs_and_thumbs_commander
public :: ctf_estimate_commander_distr
public :: ctf_estimate_commander
public :: map_cavgs_selection_commander
public :: map_cavgs_states_commander
public :: pick_commander_distr
public :: pick_commander
public :: extract_commander_distr
public :: extract_commander
public :: reextract_commander_distr
public :: reextract_commander
public :: pick_extract_commander
public :: make_pickrefs_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: preprocess_commander_distr
  contains
    procedure :: execute      => exec_preprocess_distr
end type preprocess_commander_distr

type, extends(commander_base) :: preprocess_commander
  contains
    procedure :: execute      => exec_preprocess
end type preprocess_commander

type, extends(commander_base) :: motion_correct_commander_distr
  contains
    procedure :: execute      => exec_motion_correct_distr
end type motion_correct_commander_distr

type, extends(commander_base) :: motion_correct_commander
  contains
    procedure :: execute      => exec_motion_correct
end type motion_correct_commander

type, extends(commander_base) :: gen_pspecs_and_thumbs_commander_distr
  contains
    procedure :: execute      => exec_gen_pspecs_and_thumbs_distr
end type gen_pspecs_and_thumbs_commander_distr

type, extends(commander_base) :: gen_pspecs_and_thumbs_commander
  contains
    procedure :: execute      => exec_gen_pspecs_and_thumbs
end type gen_pspecs_and_thumbs_commander

type, extends(commander_base) :: ctf_estimate_commander_distr
  contains
    procedure :: execute      => exec_ctf_estimate_distr
end type ctf_estimate_commander_distr

type, extends(commander_base) :: ctf_estimate_commander
  contains
    procedure :: execute      => exec_ctf_estimate
end type ctf_estimate_commander

type, extends(commander_base) :: map_cavgs_selection_commander
  contains
    procedure :: execute      => exec_map_cavgs_selection
end type map_cavgs_selection_commander

type, extends(commander_base) :: map_cavgs_states_commander
  contains
    procedure :: execute      => exec_map_cavgs_states
end type map_cavgs_states_commander

type, extends(commander_base) :: pick_commander_distr
  contains
    procedure :: execute      => exec_pick_distr
end type pick_commander_distr

type, extends(commander_base) :: pick_commander
  contains
    procedure :: execute      => exec_pick
end type pick_commander

type, extends(commander_base) :: extract_commander_distr
  contains
    procedure :: execute      => exec_extract_distr
end type extract_commander_distr

type, extends(commander_base) :: extract_commander
  contains
    procedure :: execute      => exec_extract
end type extract_commander

type, extends(commander_base) :: reextract_commander_distr
  contains
    procedure :: execute      => exec_reextract_distr
end type reextract_commander_distr

type, extends(commander_base) :: reextract_commander
  contains
    procedure :: execute      => exec_reextract
end type reextract_commander

type, extends(commander_base) :: pick_extract_commander
  contains
    procedure :: execute      => exec_pick_extract
end type pick_extract_commander

type, extends(commander_base) :: make_pickrefs_commander
  contains
    procedure :: execute      => exec_make_pickrefs
end type make_pickrefs_commander

contains

    subroutine exec_preprocess_distr( self, cline )
        class(preprocess_commander_distr), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(parameters)              :: params
        type(qsys_env)                :: qenv
        type(cmdline)                 :: cline_make_pickrefs
        type(chash)                   :: job_descr
        type(sp_project)              :: spproj
        real    :: pickrefs_smpd
        logical :: l_pick
        if( .not. cline%defined('oritype')         ) call cline%set('oritype',        'mic')
        if( .not. cline%defined('stream')          ) call cline%set('stream',          'no')
        if( .not. cline%defined('mkdir')           ) call cline%set('mkdir',          'yes')
        ! motion correction
        if( .not. cline%defined('trs')             ) call cline%set('trs',              20.)
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart',           8.)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',            5.)
        if( .not. cline%defined('bfac')            ) call cline%set('bfac',             50.)
        if( .not. cline%defined('mcconvention')    ) call cline%set('mcconvention','simple')
        if( .not. cline%defined('eer_upsampling')  ) call cline%set('eer_upsampling',    1.)
        if( .not. cline%defined('mcpatch')         ) call cline%set('mcpatch',        'yes')
        if( .not. cline%defined('mcpatch_thres')   ) call cline%set('mcpatch_thres',  'yes')
        if( .not. cline%defined('algorithm')       ) call cline%set('algorithm',    'patch')
        ! ctf estimation
        if( .not. cline%defined('pspecsz')         ) call cline%set('pspecsz',          512)
        if( .not. cline%defined('hp_ctf_estimate') ) call cline%set('hp_ctf_estimate',  30.)
        if( .not. cline%defined('lp_ctf_estimate') ) call cline%set('lp_ctf_estimate',   5.)
        if( .not. cline%defined('dfmin')           ) call cline%set('dfmin',          DFMIN_DEFAULT)
        if( .not. cline%defined('dfmax')           ) call cline%set('dfmax',          DFMAX_DEFAULT)
        if( .not. cline%defined('ctfpatch')        ) call cline%set('ctfpatch',       'yes')
        ! picking
        if( .not. cline%defined('picker')          ) call cline%set('picker',         'new')
        if( .not. cline%defined('lp_pick')         ) call cline%set('lp_pick',         PICK_LP_DEFAULT )
        if( .not. cline%defined('ndev')            ) call cline%set('ndev',              2.)
        if( .not. cline%defined('thres')           ) call cline%set('thres',            24.)
        if( .not. cline%defined('pick_roi')        ) call cline%set('pick_roi',        'no')
        if( .not. cline%defined('backgr_subtr')    ) call cline%set('backgr_subtr',    'no')
        ! extraction
        if( .not. cline%defined('pcontrast')       ) call cline%set('pcontrast',    'black')
        if( .not. cline%defined('extractfrommov')  ) call cline%set('extractfrommov',  'no')
        call params%new(cline)
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! read in movies
        call spproj%read(params%projfile)
        ! DISTRIBUTED EXECUTION
        params%nptcls = spproj%get_nmovies()
        if( params%nptcls == 0 )then
            THROW_HARD('no movie to process! exec_preprocess_distr')
        endif
        if( params%nparts > params%nptcls ) THROW_HARD('# partitions (nparts) must be < number of entries in filetable')
        ! deal with numlen so that length matches JOB_FINISHED indicator files
        params%numlen = len(int2str(params%nparts))
        call cline%set('numlen', params%numlen)
        ! gain reference
        call flip_gain(cline, params%gainref, params%flipgain)
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepares picking references
        l_pick = .false.
        if( cline%defined('pickrefs') )then
            cline_make_pickrefs = cline
            call cline_make_pickrefs%set('prg','make_pickrefs')
            call cline_make_pickrefs%set('mkdir','no')
            pickrefs_smpd = params%smpd / params%scale
            if( cline_make_pickrefs%defined('eer_upsampling') )then
                pickrefs_smpd = pickrefs_smpd / real(params%eer_upsampling)
            endif
            call cline_make_pickrefs%set('smpd', pickrefs_smpd)
            call qenv%exec_simple_prg_in_queue(cline_make_pickrefs, 'MAKE_PICKREFS_FINISHED')
            call cline%set('pickrefs', trim(PICKREFS_FBODY)//params%ext)
            write(logfhandle,'(A)')'>>> PREPARED PICKING TEMPLATES'
            select case(trim(trim(params%picker)))
                case('old')
                    call cline_make_pickrefs%set('neg','yes')
                case('new')
                    ! nothing to do
                case DEFAULT
                    THROW_HARD('Picker not supported!')
            end select
            call qenv%exec_simple_prg_in_queue(cline_make_pickrefs, 'MAKE_PICKREFS_FINISHED')
            call cline%set('pickrefs', trim(PICKREFS_FBODY)//params%ext)
            write(logfhandle,'(A)')'>>> PREPARED PICKING TEMPLATES'
            l_pick = .true.
        else if( cline%defined('moldiam') )then
            l_pick = .true.
        endif
        ! options compatibility
        if( trim(params%pick_roi).eq.'yes' )then
            params%backgr_subtr = 'yes'
            call cline%set('backgr_subtr', params%backgr_subtr)
        endif
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs(job_descr, algnfbody=trim(ALGN_FBODY), array=L_USE_SLURM_ARR)
        ! merge docs
        call spproj%read(params%projfile)
        call spproj%merge_algndocs(params%nptcls, params%nparts, 'mic', ALGN_FBODY)
        call spproj%kill
        ! cleanup
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_PREPROCESS NORMAL STOP ****')
    end subroutine exec_preprocess_distr

    subroutine exec_preprocess( self, cline )
        use FoX_dom
        use simple_sp_project,          only: sp_project
        use simple_motion_correct_iter, only: motion_correct_iter
        use simple_ctf_estimate_iter,   only: ctf_estimate_iter
        use simple_picker_iter,         only: picker_iter
        class(preprocess_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters)              :: params
        type(ori)                     :: o_mov
        type(ctf_estimate_iter)       :: ctfiter
        type(motion_correct_iter)     :: mciter
        type(picker_iter)             :: piter
        type(extract_commander)       :: xextract
        type(cmdline)                 :: cline_extract
        type(sp_project)              :: spproj
        type(ctfparams)               :: ctfvars
        type(Node), pointer           :: xmldoc, beamshiftnode, beamshiftnodex, beamshiftnodey
        character(len=:), allocatable :: imgkind, moviename, output_dir_picker, fbody
        character(len=:), allocatable :: moviename_forctf, moviename_intg, output_dir_motion_correct
        character(len=:), allocatable :: output_dir_ctf_estimate, output_dir_extract, micname_intg
        character(len=LONGSTRLEN)     :: boxfile, eputiltgroup
        real    :: smpd_pick
        integer :: nmovies, fromto(2), imovie, ntot, frame_counter, nptcls_out
        logical :: l_pick, l_del_forctf, l_skip_pick
        call cline%set('oritype', 'mic')
        call params%new(cline)
        if( params%scale > 1.01 )then
            THROW_HARD('scale cannot be > 1; exec_preprocess')
        endif
        l_pick = .false.
        if( cline%defined('picker') )then
            select case(trim(params%picker))
                case('old')
                    if(.not.cline%defined('pickrefs')) THROW_HARD('PICKREFS required for picker=old')
                case('new')
                    if(cline%defined('pickrefs'))then
                    else
                        if( .not.cline%defined('moldiam') )then
                            THROW_HARD('MOLDIAM required for picker=new')
                        endif
                    endif
            end select
            l_pick = .true.
        endif
        l_del_forctf = .false.
        ! read in movies
        call spproj%read( params%projfile )
        if( spproj%get_nmovies()==0 .and. spproj%get_nintgs()==0 ) THROW_HARD('No movie/micrograph to process!')
        ! output directories & naming
        output_dir_ctf_estimate        = PATH_HERE
        output_dir_motion_correct      = PATH_HERE
        if( l_pick ) output_dir_picker = PATH_HERE
        if( params%stream.eq.'yes' )then
            output_dir_ctf_estimate   = trim(DIR_CTF_ESTIMATE)
            output_dir_motion_correct = trim(DIR_MOTION_CORRECT)
            if( cline%defined('dir') )then
                output_dir_ctf_estimate   = filepath(params%dir,output_dir_ctf_estimate)//'/'
                output_dir_motion_correct = filepath(params%dir,output_dir_motion_correct)//'/'
            endif
            call simple_mkdir(output_dir_ctf_estimate,errmsg="commander_preprocess :: preprocess; ")
            call simple_mkdir(output_dir_motion_correct, errmsg="commander_preprocess :: preprocess;")
            if( l_pick )then
                output_dir_picker  = trim(DIR_PICKER)
                output_dir_extract = trim(DIR_EXTRACT)
                if( cline%defined('dir') )then
                    output_dir_picker  = filepath(params%dir,output_dir_picker)//'/'
                    output_dir_extract = filepath(params%dir,output_dir_extract)//'/'
                endif
                call simple_mkdir(output_dir_picker, errmsg="commander_preprocess :: preprocess; ")
                call simple_mkdir(output_dir_extract, errmsg="commander_preprocess :: preprocess;")
            endif
        endif
        if( cline%defined('fbody') )then
            fbody = trim(params%fbody)
        else
            fbody = ''
        endif
        ! range
        if( trim(params%stream).eq.'yes' )then
            ! STREAMING MODE
            fromto(:) = 1
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto(1) = params%fromp
                fromto(2) = params%top
            endif
        else
            ! DISTRIBUTED MODE
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto(1) = params%fromp
                fromto(2) = params%top
            else
                THROW_HARD('fromp & top args need to be defined in parallel execution; exec_preprocess')
            endif
        endif
        ntot = fromto(2) - fromto(1) + 1
        ! numlen
        if( cline%defined('numlen') )then
            ! nothing to do
        else
            params%numlen = len(int2str(nmovies))
        endif
        frame_counter = 0
        ! loop over exposures (movies)
        do imovie = fromto(1),fromto(2)
            ! fetch movie orientation
            call spproj%os_mic%get_ori(imovie, o_mov)
            ! sanity check
            if(.not.o_mov%isthere('imgkind') )cycle
            if(.not.o_mov%isthere('movie') .and. .not.o_mov%isthere('intg'))cycle
            call o_mov%getter('imgkind', imgkind)
            select case(trim(imgkind))
                case('movie')
                    ! motion_correct
                    ctfvars = spproj%get_micparams(imovie)
                    call o_mov%getter('movie', moviename)
                    if( .not.file_exists(moviename)) cycle
                    if( cline%defined('gainref') )then
                        call mciter%iterate(cline, ctfvars, o_mov, fbody, frame_counter, moviename,&
                            &output_dir_motion_correct, gainref_fname=params%gainref)
                    else
                        call mciter%iterate(cline, ctfvars, o_mov, fbody, frame_counter, moviename,&
                            &output_dir_motion_correct)
                    endif
                    moviename_forctf = mciter%get_moviename('forctf')
                    l_del_forctf     = .true.
                case('mic')
                    ctfvars = spproj%get_micparams(imovie)
                    call o_mov%getter('intg', moviename_forctf)
                case DEFAULT
                    cycle
            end select
            ! ctf_estimate
            params_glob%hp = params%hp_ctf_estimate
            params_glob%lp = max(params%fny, params%lp_ctf_estimate)
            call ctfiter%iterate(ctfvars, moviename_forctf, o_mov, output_dir_ctf_estimate, .false.)
            ! delete file after estimation
            if( l_del_forctf )then
                call o_mov%delete_entry('forctf')
                call del_file(moviename_forctf)
            endif
            ! optional rejection
            l_skip_pick = .false.
            if( trim(params%stream).eq.'yes' .and. trim(params%reject_mics).eq.'yes' )then
                ! based on CTFRES
                if( l_pick .and. o_mov%isthere('ctfres') )then
                    l_skip_pick = o_mov%get('ctfres') > (params%ctfresthreshold-0.001)
                    if( l_skip_pick ) call o_mov%set('nptcls',0)
                end if
                ! based on ice fraction
                if( l_pick .and. .not.l_skip_pick .and. o_mov%isthere('icefrac') )then
                    l_skip_pick = o_mov%get('icefrac') > (params%icefracthreshold-0.001)
                    if( l_skip_pick ) call o_mov%set('nptcls',0)
                endif
            endif
            ! read xml
            if(o_mov%isthere('meta') .and. file_exists(trim(o_mov%get_static('meta')))) then
                xmldoc => parseFile(trim(o_mov%get_static("meta")))
                beamshiftnode  => item(getElementsByTagname(xmldoc, "BeamShift"),   0)
                beamshiftnodex => item(getElementsByTagname(beamshiftnode, "a:_x"), 0)
                beamshiftnodey => item(getElementsByTagname(beamshiftnode, "a:_y"), 0)
                call o_mov%set("shiftx", str2real(getTextContent(beamshiftnodex)))
                call o_mov%set("shifty", str2real(getTextContent(beamshiftnodey)))
                call destroy(xmldoc)
            end if
            ! assign beamtilt if EPU
            if(o_mov%isthere('intg')) then
                call o_mov%getter('intg', micname_intg)
                micname_intg = trim(basename(micname_intg))
                if(micname_intg(:8) == 'FoilHole') then
                    ! EPU filename
                    eputiltgroup = micname_intg(index(micname_intg, 'Data_') + 5:)
                    eputiltgroup = eputiltgroup(:index(eputiltgroup, '_') - 1)
                    call o_mov%set("tiltgrp", str2real(eputiltgroup))
                else
                    call o_mov%set("tiltgrp", 0.0)
                end if
            end if
            ! update project
            call spproj%os_mic%set_ori(imovie, o_mov)
            ! pick
            if( l_pick .and. (.not.l_skip_pick) )then
                smpd_pick      = o_mov%get('smpd')
                params_glob%lp = max(2.*smpd_pick, params%lp_pick)
                call o_mov%getter('intg', moviename_intg)
                call piter%iterate(cline, smpd_pick, moviename_intg, output_dir_picker, boxfile, nptcls_out)
                call o_mov%set('nptcls',  nptcls_out)
                if( nptcls_out > 0 )then
                    call o_mov%set('boxfile', trim(boxfile))
                else
                    call o_mov%set('state',0.)
                endif
                ! update project
                call spproj%os_mic%set_ori(imovie, o_mov)
                ! extract particles
                if( trim(params%stream) .eq. 'yes' )then
                    ! needs to write and re-read project at the end as extract overwrites it
                    call spproj%write_segment_inside(params%oritype)
                    if( nptcls_out > 0 )then
                        cline_extract = cline
                        call cline_extract%set('smpd',      o_mov%get('smpd')) ! in case of scaling
                        call cline_extract%set('dir',       trim(output_dir_extract))
                        call cline_extract%set('pcontrast', params%pcontrast)
                        call cline_extract%delete('msk')
                        if( cline%defined('box_extract') )call cline_extract%set('box', params%box_extract)
                        call xextract%execute(cline_extract)
                        call spproj%kill
                    endif
                endif
            endif
        end do
        if( trim(params%stream).eq.'yes' )then
            if( (.not.l_pick) .or. l_skip_pick)then
                ! because extract performs the writing otherwise
                call spproj%write_segment_inside(params%oritype)
            endif
        else
            call binwrite_oritab(params%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        endif
        call piter%kill
        call o_mov%kill
        ! end gracefully
        call qsys_job_finished( 'simple_commander_preprocess :: exec_preprocess' )
        call simple_end('**** SIMPLE_PREPROCESS NORMAL STOP ****')
    end subroutine exec_preprocess

    subroutine exec_motion_correct_distr( self, cline )
        class(motion_correct_commander_distr), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        type(qsys_env)   :: qenv
        type(chash)      :: job_descr
        if( .not. cline%defined('mkdir')         ) call cline%set('mkdir',           'yes')
        if( .not. cline%defined('trs')           ) call cline%set('trs',               20.)
        if( .not. cline%defined('lpstart')       ) call cline%set('lpstart',            8.)
        if( .not. cline%defined('lpstop')        ) call cline%set('lpstop',             5.)
        if( .not. cline%defined('bfac')          ) call cline%set('bfac',              50.)
        if( .not. cline%defined('mcconvention')  ) call cline%set('mcconvention', 'simple')
        if( .not. cline%defined('wcrit')         ) call cline%set('wcrit',       'softmax')
        if( .not. cline%defined('eer_upsampling')) call cline%set('eer_upsampling',      1)
        if( .not. cline%defined('mcpatch')       ) call cline%set('mcpatch',         'yes')
        if( .not. cline%defined('mcpatch_thres') ) call cline%set('mcpatch_thres',   'yes')
        if( .not. cline%defined('algorithm')     ) call cline%set('algorithm',     'patch')
        call cline%set('oritype', 'mic')
        call params%new(cline)
        call cline%set('numlen', params%numlen)
        ! sanity check
        call spproj%read_segment(params%oritype, params%projfile)
        if( spproj%get_nmovies() ==0 ) THROW_HARD('no movies to process! exec_motion_correct_distr')
        call spproj%kill
        ! gain reference
        call flip_gain(cline, params%gainref, params%flipgain)
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs(job_descr, algnfbody=trim(ALGN_FBODY), array=L_USE_SLURM_ARR)
        ! merge docs
        call spproj%read(params%projfile)
        call spproj%merge_algndocs(params%nptcls, params%nparts, 'mic', ALGN_FBODY)
        call spproj%kill
        ! clean
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_MOTION_CORRECT NORMAL STOP ****')
    end subroutine exec_motion_correct_distr

    subroutine exec_motion_correct( self, cline )
        use simple_sp_project,          only: sp_project
        use simple_motion_correct_iter, only: motion_correct_iter
        class(motion_correct_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline !< command line input
        type(parameters)              :: params
        type(motion_correct_iter)     :: mciter
        type(ctfparams)               :: ctfvars
        type(sp_project)              :: spproj
        type(ori)                     :: o
        character(len=:), allocatable :: output_dir, moviename, fbody
        integer :: nmovies, fromto(2), imovie, ntot, frame_counter, cnt
        call cline%set('oritype', 'mic')
        call params%new(cline)
        call spproj%read(params%projfile)
        ! sanity check
        nmovies = spproj%get_nmovies()
        if( nmovies == 0 )then
            THROW_HARD('No movie to process!')
        endif
        if( params%scale > 1.01 )then
            THROW_HARD('scale cannot be > 1; exec_motion_correct')
        endif
        if( cline%defined('gainref') )then
            if(.not.file_exists(params%gainref) )then
                THROW_HARD('gain reference: '//trim(params%gainref)//' not found; motion_correct')
            endif
        endif
        ! output directory & names
        output_dir = PATH_HERE
        if( cline%defined('fbody') )then
            fbody = trim(params%fbody)
        else
            fbody = ''
        endif
        ! determine loop range & fetch movies oris object
        if( cline%defined('fromp') .and. cline%defined('top') )then
            fromto = [params%fromp, params%top]
        else
            THROW_HARD('fromp & top args need to be defined in parallel execution; motion_correct')
        endif
        ntot = fromto(2) - fromto(1) + 1
        ! align
        frame_counter = 0
        cnt = 0
        do imovie=fromto(1),fromto(2)
            call spproj%os_mic%get_ori(imovie, o)
            if( o%isthere('imgkind') )then
                if( o%isthere('movie') .or. o%isthere('mic') )then
                    cnt = cnt + 1
                    call o%getter('movie', moviename)
                    ctfvars = spproj%get_micparams(imovie)
                    if( cline%defined('gainref') )then
                        call mciter%iterate(cline, ctfvars, o, fbody, frame_counter, moviename, trim(output_dir), gainref_fname=params%gainref)
                    else
                        call mciter%iterate(cline, ctfvars, o, fbody, frame_counter, moviename, trim(output_dir))
                    endif
                    call spproj%os_mic%set_ori(imovie, o)
                    write(logfhandle,'(f4.0,1x,a)') 100.*(real(cnt)/real(ntot)), 'percent of the movies processed'
                endif
            endif
        end do
        ! output
        call binwrite_oritab(params%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        call o%kill
        ! end gracefully
        call qsys_job_finished(  'simple_commander_preprocess :: exec_motion_correct' )
        call simple_end('**** SIMPLE_MOTION_CORRECT NORMAL STOP ****')
    end subroutine exec_motion_correct

    subroutine exec_gen_pspecs_and_thumbs_distr( self, cline )
        class(gen_pspecs_and_thumbs_commander_distr), intent(inout) :: self
        class(cmdline),                               intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        type(qsys_env)   :: qenv
        type(chash)      :: job_descr
        integer          :: nintgs
        call cline%set('oritype', 'mic')
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        params%numlen = len(int2str(params%nparts))
        call cline%set('numlen', params%numlen)
        ! sanity check
        call spproj%read_segment(params%oritype, params%projfile)
        nintgs = spproj%get_nintgs()
        if( nintgs ==0 )then
            THROW_HARD('no integrated movies to process! exec_gen_pspecs_and_thumbs_distr')
        endif
        if( params%nparts > nintgs )then
            call cline%set('nparts', nintgs)
            params%nparts = nintgs
        endif
        call spproj%kill
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs(job_descr, algnfbody=trim(ALGN_FBODY), array=L_USE_SLURM_ARR)
        ! merge docs
        call spproj%read(params%projfile)
        call spproj%merge_algndocs(params%nptcls, params%nparts, 'mic', ALGN_FBODY)
        call spproj%kill
        ! clean
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_GEN_PSPECS_AND_THUMBS NORMAL STOP ****')
    end subroutine exec_gen_pspecs_and_thumbs_distr

    subroutine exec_gen_pspecs_and_thumbs( self, cline )
        use simple_sp_project,       only: sp_project
        use simple_pspec_thumb_iter, only: pspec_thumb_iter
        class(gen_pspecs_and_thumbs_commander), intent(inout) :: self
        class(cmdline),                         intent(inout) :: cline !< command line input
        type(parameters)              :: params
        type(pspec_thumb_iter)        :: ptiter
        type(sp_project)              :: spproj
        type(ori)                     :: o
        character(len=:), allocatable :: output_dir, moviename_intg, imgkind
        integer :: nintgs, fromto(2), iintg, ntot, cnt
        call cline%set('oritype', 'mic')
        call params%new(cline)
        call spproj%read(params%projfile)
        ! sanity check
        nintgs = spproj%get_nintgs()
        if( nintgs == 0 )then
            THROW_HARD('No integrated movies to process!')
        endif
        ! output directory
        output_dir = PATH_HERE
        ! determine loop range & fetch movies oris object
        if( params%l_distr_exec )then
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto = [params%fromp, params%top]
            else
                THROW_HARD('fromp & top args need to be defined in parallel execution; gen_pspecs_and_thumbs')
            endif
        else
            fromto = [1,nintgs]
        endif
        ntot = fromto(2) - fromto(1) + 1
        ! align
        cnt = 0
        do iintg=fromto(1),fromto(2)
            call spproj%os_mic%get_ori(iintg, o)
            if( o%isthere('imgkind').and.o%isthere('intg') )then
                cnt = cnt + 1
                call o%getter('imgkind', imgkind)
                if( imgkind.ne.'mic' )cycle
                call o%getter('intg', moviename_intg)
                call ptiter%iterate(o, moviename_intg, trim(output_dir))
                call spproj%os_mic%set_ori(iintg, o)
                write(logfhandle,'(f4.0,1x,a)') 100.*(real(cnt)/real(ntot)), 'percent of the integrated movies processed'
            endif
        end do
        ! output
        call binwrite_oritab(params%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        call o%kill
        ! end gracefully
        call qsys_job_finished('simple_commander_preprocess :: exec_gen_pspecs_and_thumbs')
        call simple_end('**** SIMPLE_GEN_PSPECS_AND_THUMBS NORMAL STOP ****')
    end subroutine exec_gen_pspecs_and_thumbs

    subroutine exec_ctf_estimate_distr( self, cline )
        class(ctf_estimate_commander_distr), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(chash)                   :: job_descr
        type(qsys_env)                :: qenv
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',  'yes')
        if( .not. cline%defined('pspecsz') ) call cline%set('pspecsz',  512)
        if( .not. cline%defined('hp')      ) call cline%set('hp',       30.)
        if( .not. cline%defined('lp')      ) call cline%set('lp',        5.)
        if( .not. cline%defined('dfmin')   ) call cline%set('dfmin',    DFMIN_DEFAULT)
        if( .not. cline%defined('dfmax')   ) call cline%set('dfmax',    DFMAX_DEFAULT)
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'mic')
        if( .not. cline%defined('ctfpatch')) call cline%set('ctfpatch','yes')
        call params%new(cline)
        ! sanity check
        call spproj%read_segment(params%oritype, params%projfile)
        if( spproj%get_nintgs() ==0 )then
            THROW_HARD('no micrograph to process! exec_ctf_estimate_distr')
        endif
        call spproj%kill
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        params%numlen = len(int2str(params%nparts))
        call cline%set('numlen', params%numlen)
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY), array=L_USE_SLURM_ARR)
        ! merge docs
        call spproj%read(params%projfile)
        call spproj%merge_algndocs(params%nptcls, params%nparts, 'mic', ALGN_FBODY)
        ! cleanup
        call spproj%kill
        call qsys_cleanup
        ! graceful ending
        call simple_end('**** SIMPLE_DISTR_CTF_ESTIMATE NORMAL STOP ****')
    end subroutine exec_ctf_estimate_distr

    subroutine exec_ctf_estimate( self, cline )
        use simple_sp_project,          only: sp_project
        use simple_ctf_estimate_iter,   only: ctf_estimate_iter
        class(ctf_estimate_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline  !< command line input
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(ctf_estimate_iter)       :: ctfiter
        type(ctfparams)               :: ctfvars
        type(ori)                     :: o
        character(len=:), allocatable :: intg_forctf, output_dir, imgkind
        integer                       :: fromto(2), imic, ntot, cnt, state
        logical                       :: l_gen_thumb, l_del_forctf
        call cline%set('oritype', 'mic')
        call params%new(cline)
        call spproj%read(params%projfile)
        ! read in integrated movies
        if( spproj%get_nintgs() == 0 ) THROW_HARD('No integrated micrograph to process!')
        ! output directory
        output_dir = PATH_HERE
        ! parameters & loop range
        if( params%stream .eq. 'yes' )then
            ! determine loop range
            fromto(:) = 1
        else
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto = [params%fromp, params%top]
            else
                THROW_HARD('fromp & top args need to be defined in parallel execution; exec_ctf_estimate')
            endif
        endif
        ntot = fromto(2) - fromto(1) + 1
        ! loop over exposures (movies)
        cnt = 0
        do imic = fromto(1),fromto(2)
            cnt   = cnt + 1
            call spproj%os_mic%get_ori(imic, o)
            state = 1
            if( o%isthere('state') ) state = o%get_state()
            if( state == 0 ) cycle
            if( o%isthere('imgkind') )then
                call o%getter('imgkind', imgkind)
                if( imgkind.ne.'mic' )cycle
                l_del_forctf = .false.
                if( o%isthere('forctf') )then
                    call o%getter('forctf', intg_forctf)
                    if( file_exists(intg_forctf) )then
                        l_del_forctf = .true.
                    else
                        if( o%isthere('intg') )then
                            call o%getter('intg', intg_forctf)
                        endif
                    endif
                else if( o%isthere('intg') )then
                    call o%getter('intg', intg_forctf)
                else
                    THROW_HARD('no image available (forctf|intg) for CTF fittings :: exec_ctf_estimate')
                endif
                l_gen_thumb = .not. o%isthere('thumb')
                ctfvars     = o%get_ctfvars()
                call ctfiter%iterate( ctfvars, intg_forctf, o, trim(output_dir), l_gen_thumb)
                ! delete file after estimation
                if( l_del_forctf )then
                    call o%delete_entry('forctf')
                    call del_file(intg_forctf)
                endif
                ! update project
                call spproj%os_mic%set_ori(imic, o)
            endif
            write(logfhandle,'(f4.0,1x,a)') 100.*(real(cnt)/real(ntot)), 'percent of the micrographs processed'
        end do
        ! output
        call binwrite_oritab(params%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        call o%kill
        ! end gracefully
        call qsys_job_finished(  'simple_commander_preprocess :: exec_ctf_estimate' )
        call simple_end('**** SIMPLE_CTF_ESTIMATE NORMAL STOP ****')
    end subroutine exec_ctf_estimate

    subroutine exec_map_cavgs_selection( self, cline )
        use simple_corrmat,             only: calc_cartesian_corrmat
        class(map_cavgs_selection_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(parameters)              :: params
        type(builder)                 :: build
        type(image),      allocatable :: imgs_sel(:), imgs_all(:)
        integer,          allocatable :: states(:)
        real,             allocatable :: correlations(:,:)
        character(len=:), allocatable :: cavgstk
        integer :: iimg, isel, nall, nsel, loc(1), lfoo(3)
        real    :: smpd
        call cline%set('dir_exec', 'selection')
        call cline%set('mkdir',    'yes')
        if( .not.cline%defined('prune')   ) call cline%set('prune',   'no')
        if( .not.cline%defined('imgkind') ) call cline%set('imgkind', 'cavg')
        call build%init_params_and_build_spproj(cline,params)
        ! find number of selected cavgs
        call find_ldim_nptcls(params%stk2, lfoo, nsel)
        if( cline%defined('ares') ) nsel = int(params%ares)
        ! find number of original cavgs
        if( .not. cline%defined('stk' ) )then
            call build%spproj%get_cavgs_stk(cavgstk, nall, smpd, imgkind=params%imgkind)
            params%stk = trim(cavgstk)
        else
            call find_ldim_nptcls(params%stk, lfoo, nall)
        endif
        ! read images
        allocate(imgs_sel(nsel), imgs_all(nall))
        do isel=1,nsel
            call imgs_sel(isel)%new([params%box,params%box,1], params%smpd)
            call imgs_sel(isel)%read(params%stk2, isel)
        end do
        do iimg=1,nall
            call imgs_all(iimg)%new([params%box,params%box,1], params%smpd)
            call imgs_all(iimg)%read(params%stk, iimg)
        end do
        write(logfhandle,'(a)') '>>> CALCULATING CORRELATIONS'
        call calc_cartesian_corrmat(imgs_sel, imgs_all, correlations)
        ! create the states array for mapping the selection
        allocate(states(nall), source=0)
        do isel=1,nsel
            loc = maxloc(correlations(isel,:))
            states(loc(1)) = 1
        end do
        ! communicate selection to project
        call build%spproj%map_cavgs_selection(states)
        ! optional pruning
        if( trim(params%prune).eq.'yes') call build%spproj%prune_particles
        ! this needs to be a full write as many segments are updated
        call build%spproj%write
        ! end gracefully
        call simple_end('**** SIMPLE_MAP_CAVGS_SELECTION NORMAL STOP ****')
    end subroutine exec_map_cavgs_selection

    subroutine exec_map_cavgs_states( self, cline )
        use simple_corrmat, only: calc_cartesian_corrmat
        class(map_cavgs_states_commander), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline !< command line input
        type(parameters)                   :: params
        type(builder)                      :: build
        type(image),           allocatable :: imgs_sel(:), imgs_all(:)
        integer,               allocatable :: states(:)
        real,                  allocatable :: correlations(:,:)
        character(len=:),      allocatable :: cavgstk, fname
        character(LONGSTRLEN), allocatable :: stkfnames(:)
        integer :: iimg, isel, nall, nsel, loc(1), lfoo(3), s
        real    :: smpd
        call cline%set('dir_exec', 'state_mapping')
        call cline%set('mkdir',    'yes')
        call build%init_params_and_build_spproj(cline,params)
        call read_filetable(params%stktab, stkfnames)
        ! find number of original cavgs
        if( .not. cline%defined('stk' ) )then
            call build%spproj%get_cavgs_stk(cavgstk, nall, smpd)
            params%stk = trim(cavgstk)
        else
            call find_ldim_nptcls(params%stk, lfoo, nall)
        endif
        ! read images
        allocate(imgs_all(nall))
        do iimg=1,nall
            call imgs_all(iimg)%new([params%box,params%box,1], params%smpd)
            call imgs_all(iimg)%read(params%stk, iimg)
        end do
        ! create the states array for mapping the selection
        allocate(states(nall), source=0)
        do s = 1,size(stkfnames)
            ! find number of selected cavgs
            fname = '../'//trim(stkfnames(s))
            call find_ldim_nptcls(fname, lfoo, nsel)
            ! read images
            allocate(imgs_sel(nsel))
            do isel=1,nsel
                call imgs_sel(isel)%new([params%box,params%box,1], params%smpd)
                call imgs_sel(isel)%read(fname, isel)
            end do
            call calc_cartesian_corrmat(imgs_sel, imgs_all, correlations)
            do isel=1,nsel
                loc = maxloc(correlations(isel,:))
                states(loc(1)) = s
            end do
            ! destruct
            do isel=1,nsel
                call imgs_sel(isel)%kill
            end do
            deallocate(imgs_sel)
        end do
        ! communicate selection to project
        call build%spproj%map_cavgs_selection(states)
        ! this needs to be a full write as many segments are updated
        call build%spproj%write
        ! end gracefully
        call simple_end('**** SIMPLE_MAP_CAVGS_SELECTION NORMAL STOP ****')
    end subroutine exec_map_cavgs_states

    subroutine exec_pick_distr( self, cline )
        class(pick_commander_distr), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        type(cmdline)    :: cline_make_pickrefs
        type(qsys_env)   :: qenv
        type(chash)      :: job_descr
        character(len=:), allocatable :: which_picker
        real,             allocatable :: moldiams(:), dists(:), moldiams_prob(:), moldiams_ref(:)
        logical,          allocatable :: mask(:)
        real    :: moldiam_est
        integer :: nmics, loc, i, j, n
        logical :: templates_provided
        if( .not. cline%defined('mkdir')       ) call cline%set('mkdir',       'yes')
        if( .not. cline%defined('pcontrast')   ) call cline%set('pcontrast', 'black')
        if( .not. cline%defined('oritype')     ) call cline%set('oritype',     'mic')
        if( .not. cline%defined('thres')       ) call cline%set('thres',         24.)
        if( .not. cline%defined('pick_roi')    ) call cline%set('pick_roi',     'no')
        if( .not. cline%defined('backgr_subtr')) call cline%set('backgr_subtr', 'no') 
        if( .not. cline%defined('picker')      ) call cline%set('picker',      'new')
        if( .not. cline%defined('lp')          ) call cline%set('lp',PICK_LP_DEFAULT)
        which_picker = cline%get_carg('picker')
        if( trim(which_picker) .eq. 'seg' )then
            if( .not. cline%defined('ndev')        ) call cline%set('ndev',      1.5)
        else
            if( .not. cline%defined('ndev')        ) call cline%set('ndev',       2.)
        endif
        call params%new(cline)
        ! sanity check
        call spproj%read_segment(params%oritype, params%projfile)
        nmics = spproj%get_nintgs()
        if( nmics == 0 ) THROW_HARD('No micrograph to process! exec_pick_distr')
        call spproj%kill
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        if( cline%defined('nparts') )then
            params%nparts = min(params%nparts, nmics)
            call cline%set('nparts', params%nparts)
        endif
        params%numlen = len(int2str(params%nparts))
        call cline%set('numlen', params%numlen)
        ! more sanity checks
        if( trim(params%pick_roi).eq.'yes' )then
            params%backgr_subtr = 'yes'
            call cline%set('backgr_subtr', params%backgr_subtr)
        endif
        templates_provided = cline%defined('pickrefs')
        select case(trim(params%picker))
        case('old')
            if( .not.templates_provided ) THROW_HARD('Old picker requires pickrefs (2D picking references) input')
        case('new')
            if( templates_provided )then
                ! reference-based picking
            else if( cline%defined('moldiam') )then
                if( cline%defined('moldiam_max') .and. cline%defined('nmoldiams') )then
                    ! multipick with nmoldiams diameters ranging from moldiam to moldiam_max
                    ! output is determined diameter moldiam_opt
                else if( cline%defined('moldiam_max') .or. cline%defined('nmoldiams') )then
                    THROW_HARD('MOLDIAM, MOLDIAM_MAX & NMOLDIAMS must be provided for determination of optimal diameter!')
                else
                    ! reference-free single diameter picking
                endif
            elseif( cline%defined('multi_moldiams') )then
                ! reference-free picking with inputted diameters
            else
                THROW_HARD('Unsupported new picker mode')
            endif
        end select
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepares picking references
        if( templates_provided )then
            cline_make_pickrefs = cline
            call cline_make_pickrefs%set('prg',  'make_pickrefs')
            call cline_make_pickrefs%set('smpd', params%smpd)
            call cline_make_pickrefs%set('mkdir','no')
            if( trim(params%picker).eq.'old' ) call cline_make_pickrefs%set('neg','yes')
            call qenv%exec_simple_prg_in_queue(cline_make_pickrefs, 'MAKE_PICKREFS_FINISHED')
            call cline%set('pickrefs', trim(PICKREFS_FBODY)//params%ext)
            write(logfhandle,'(A)')'>>> PREPARED PICKING TEMPLATES'
        endif
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs(job_descr, algnfbody=trim(ALGN_FBODY), array=L_USE_SLURM_ARR)
        ! merge docs
        call spproj%read(params%projfile)
        call spproj%merge_algndocs(params%nptcls, params%nparts, 'mic', ALGN_FBODY)
        if( cline%defined('nmoldiams') )then
            n            = spproj%os_mic%get_noris()
            mask         = spproj%os_mic%get_all('state') > 0.5
            moldiams     = spproj%os_mic%get_all('moldiam')
            moldiams_ref = equispaced_vals(params%moldiam, params%moldiam_max, params%nmoldiams)
            allocate(dists(params%nmoldiams), moldiams_prob(params%nmoldiams), source=0.)
            do i = 1, n
                if( .not.mask(i) ) cycle
                dists = abs(moldiams(i) - moldiams_ref)
                loc   = minloc(dists, dim=1)
                moldiams_prob(loc) = moldiams_prob(loc) + 1.
            end do
            moldiams_prob = moldiams_prob / real(count(mask))
            moldiam_est   = sum(moldiams_ref * moldiams_prob)
            do j = 1, params%nmoldiams
                write(logfhandle,'(a,1x,f6.1,1x,a,1x,f7.2)') 'Molecular diameter:', moldiams_ref(j), 'probability:', moldiams_prob(j)
            end do
            write(logfhandle,'(a,1x,f7.2)') 'Suggested single molecular diameter:', moldiam_est
            call cline%set('moldiam', moldiam_est)
        endif
        ! cleanup
        call qsys_cleanup
        ! graceful exit
        call simple_end('**** SIMPLE_DISTR_PICK NORMAL STOP ****')
    end subroutine exec_pick_distr

    subroutine exec_pick( self, cline )
        use simple_picker_iter, only: picker_iter
        class(pick_commander), intent(inout) :: self
        class(cmdline),        intent(inout) :: cline !< command line input
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(picker_iter)             :: piter
        type(ori)                     :: o
        character(len=:), allocatable :: output_dir, intg_name, imgkind
        character(len=LONGSTRLEN)     :: boxfile
        integer                       :: fromto(2), imic, ntot, nptcls_out, cnt, state
        real                          :: moldiam_opt
        call cline%set('oritype', 'mic')
        call params%new(cline)
        ! output directory
        output_dir = PATH_HERE
        ! parameters & loop range
        if( params%stream .eq. 'yes' )then
            ! determine loop range
            fromto(:) = 1
        else
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto = [params%fromp, params%top]
            else
                THROW_HARD('fromp & top args need to be defined in parallel execution; exec_pick')
            endif
        endif
        ntot = fromto(2) - fromto(1) + 1
        ! read project file
        call spproj%read(params%projfile)
        ! look for movies
        if( spproj%get_nintgs() == 0 )then
            THROW_HARD('No integrated micrograph to process!')
        endif
        ! perform picking
        cnt = 0
        do imic=fromto(1),fromto(2)
            cnt   = cnt + 1
            call spproj%os_mic%get_ori(imic, o)
            state = 1
            if( o%isthere('state') ) state = o%get_state()
            if( state == 0 ) cycle
            if( o%isthere('imgkind') )then
                call o%getter('imgkind', imgkind)
                if( imgkind.ne.'mic' )cycle
                call o%getter('intg', intg_name)
                call piter%iterate(cline, params%smpd, intg_name, output_dir, boxfile, nptcls_out, moldiam_opt=moldiam_opt)
                if( params_glob%nmoldiams > 1 )then
                    call spproj%os_mic%set(imic, 'moldiam', moldiam_opt)
                else
                    call spproj%set_boxfile(imic, boxfile, nptcls=nptcls_out)
                endif
            endif
            write(logfhandle,'(f4.0,1x,a)') 100.*(real(cnt)/real(ntot)), 'percent of the micrographs processed'
        end do
        ! output
        call binwrite_oritab(params%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        ! cleanup
        call o%kill
        call spproj%kill
        call piter%kill
        ! end gracefully
        call qsys_job_finished( 'simple_commander_preprocess :: exec_pick' )
        call simple_end('**** SIMPLE_PICK NORMAL STOP ****')
    end subroutine exec_pick

    subroutine exec_extract_distr( self, cline )
        class(extract_commander_distr), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline !< command line input
        type(parameters)                        :: params
        type(sp_project)                        :: spproj, spproj_part
        type(qsys_env)                          :: qenv
        type(chash)                             :: job_descr
        type(ori)                               :: o_mic, o_tmp
        type(oris)                              :: os_stk
        character(len=LONGSTRLEN),  allocatable :: boxfiles(:), stktab(:), parts_fname(:)
        character(len=:),           allocatable :: mic_name, imgkind, boxfile_name
        real    :: dfx,dfy,ogid,gid
        integer :: boxcoords(2), lfoo(3)
        integer :: nframes,imic,i,nmics_tot,numlen,nmics,cnt,state,istk,nstks,ipart
        if( .not. cline%defined('mkdir')         ) call cline%set('mkdir',          'yes')
        if( .not. cline%defined('outside')       ) call cline%set('outside',         'no')
        if( .not. cline%defined('pcontrast')     ) call cline%set('pcontrast',    'black')
        if( .not. cline%defined('stream')        ) call cline%set('stream',          'no')
        if( .not. cline%defined('extractfrommov')) call cline%set('extractfrommov',  'no')
        if( .not. cline%defined('backgr_subtr')  ) call cline%set('backgr_subtr',    'no')
        if( cline%defined('ctf') )then
            if( cline%get_carg('ctf').ne.'flip' .and. cline%get_carg('ctf').ne.'no' )then
                THROW_HARD('Only CTF=NO/FLIP are allowed')
            endif
        endif
        call cline%set('oritype', 'mic')
        call params%new(cline)
        call cline%set('mkdir', 'no')
        ! read in integrated movies
        call spproj%read(params%projfile)
        if( spproj%get_nintgs() == 0 ) THROW_HARD('No integrated micrograph to process!')
        nmics_tot = spproj%os_mic%get_noris()
        if( nmics_tot < params%nparts ) params%nparts = nmics_tot
        ! wipes previous stacks & particles
        call spproj%os_stk%kill
        call spproj%os_ptcl2D%kill
        call spproj%os_ptcl3D%kill
        call spproj%os_cls2D%kill
        call spproj%os_cls3D%kill
        call spproj%os_out%kill
        ! input directory
        if( cline%defined('dir_box') )then
            if( params%mkdir.eq.'yes' .and. params%dir_box(1:1).ne.'/')then
                params%dir_box = trim(filepath(PATH_PARENT,params%dir_box))
            endif
            params%dir_box = simple_abspath(params%dir_box)
            if( file_exists(params%dir_box) )then
                call simple_list_files_regexp(params%dir_box,'\.box$', boxfiles)
                if(.not.allocated(boxfiles))then
                    write(logfhandle,*)'No box file found in ', trim(params%dir_box), '; simple_commander_preprocess::exec_extract 1'
                    THROW_HARD('No box file found; exec_extract, 1')
                endif
                if(size(boxfiles)==0)then
                    write(logfhandle,*)'No box file found in ', trim(params%dir_box), '; simple_commander_preprocess::exec_extract 2'
                    THROW_HARD('No box file found; exec_extract 2')
                endif
            else
                write(logfhandle,*)'Directory does not exist: ', trim(params%dir_box), 'simple_commander_preprocess::exec_extract'
                THROW_HARD('box directory does not exist; exec_extract')
            endif
            call cline%set('dir_box', params%dir_box)
        endif
        call spproj%write(params%projfile)
        ! sanity checks
        nmics  = 0
        do imic = 1, nmics_tot
            call spproj%os_mic%get_ori(imic, o_mic)
            state = 1
            if( o_mic%isthere('state') ) state = o_mic%get_state()
            if( state == 0 ) cycle
            if( .not. o_mic%isthere('imgkind') )cycle
            if( .not. o_mic%isthere('intg')    )cycle
            call o_mic%getter('imgkind', imgkind)
            if( trim(imgkind).ne.'mic') cycle
            call o_mic%getter('intg', mic_name)
            if( .not.file_exists(mic_name) )cycle
            ! box input
            if( cline%defined('dir_box') )then
                boxfile_name = boxfile_from_mic(mic_name)
                if(trim(boxfile_name).eq.NIL)cycle
            else
                call o_mic%getter('boxfile', boxfile_name)
                if( .not.file_exists(boxfile_name) )cycle
            endif
            ! get number of frames from stack
            call find_ldim_nptcls(mic_name, lfoo, nframes )
            if( nframes > 1 ) THROW_HARD('multi-frame extraction not supported; exec_extract')
            ! update counter
            nmics = nmics + 1
        enddo
        if( nmics == 0 ) THROW_HARD('No particles to extract! exec_extract')
        ! progress
        call progressfile_init_parts(params%nparts) 
        ! DISTRIBUTED EXTRACTION
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY), array=L_USE_SLURM_ARR)
        ! ASSEMBLY
        allocate(parts_fname(params%nparts))
        numlen = len(int2str(params%nparts))
        do ipart = 1,params%nparts
            parts_fname(ipart) = trim(ALGN_FBODY)//int2str_pad(ipart,numlen)//trim(METADATA_EXT)
        enddo
        ! copy updated micrographs
        cnt   = 0
        nstks = 0
        do ipart = 1,params%nparts
            call spproj_part%read_segment('mic',parts_fname(ipart))
            do imic = 1,spproj_part%os_mic%get_noris()
                cnt = cnt + 1
                call spproj_part%os_mic%get_ori(imic, o_mic)
                call spproj%os_mic%set_ori(cnt,o_mic)
                if( o_mic%isthere('nptcls') )then
                    if( o_mic%get_int('nptcls') > 0 ) nstks = nstks + 1
                endif
            enddo
            call spproj_part%kill
        enddo
        if( cnt /= nmics_tot ) THROW_HARD('Inconstistent number of micrographs in individual projects')
        ! fetch stacks table
        if( nstks > 0 )then
            call os_stk%new(nstks, is_ptcl=.false.)
            allocate(stktab(nstks))
            cnt = 0
            do ipart = 1,params%nparts
                call spproj_part%read_segment('stk',parts_fname(ipart))
                do istk = 1,spproj_part%os_stk%get_noris()
                    cnt = cnt + 1
                    call spproj_part%os_stk%get_ori(istk, o_tmp)
                    call os_stk%set_ori(cnt,o_tmp)
                    stktab(cnt) = os_stk%get_static(cnt,'stk')
                enddo
                call spproj_part%kill
            enddo
            ! import stacks into project
            call spproj%add_stktab(stktab,os_stk)
            ! transfer particles locations to ptcl2D & defocus to 2D/3D
            cnt = 0
            do ipart = 1,params%nparts
                call spproj_part%read_segment('ptcl2D',parts_fname(ipart))
                do i = 1,spproj_part%os_ptcl2D%get_noris()
                    cnt = cnt + 1
                    ! picking coordinates
                    call spproj_part%get_boxcoords(i, boxcoords)
                    call spproj%set_boxcoords(cnt, boxcoords)
                    ! defocus from patch-based ctf estimation
                    if( spproj_part%os_ptcl2D%isthere(i,'dfx') )then
                        dfx = spproj_part%os_ptcl2D%get_dfx(i)
                        dfy = spproj_part%os_ptcl2D%get_dfy(i)
                        call spproj%os_ptcl2D%set_dfx(cnt,dfx)
                        call spproj%os_ptcl2D%set_dfy(cnt,dfy)
                        call spproj%os_ptcl3D%set_dfx(cnt,dfx)
                        call spproj%os_ptcl3D%set_dfy(cnt,dfy)
                    endif
                    !optics group id
                    if( spproj_part%os_ptcl2D%isthere(i,'ogid') )then
                        ogid = spproj_part%os_ptcl2D%get(i, 'ogid')
                        call spproj%os_ptcl2D%set(cnt,'ogid',ogid)
                        call spproj%os_ptcl3D%set(cnt,'ogid',ogid)
                    endif
                    !group id
                    if( spproj_part%os_ptcl2D%isthere(i,'gid') )then
                        gid = spproj_part%os_ptcl2D%get(i, 'gid')
                        call spproj%os_ptcl2D%set(cnt,'gid',gid)
                        call spproj%os_ptcl3D%set(cnt,'gid',gid)
                    endif
                enddo
                call spproj_part%kill
            enddo
            call os_stk%kill
        endif
        ! final write
        call spproj%write(params%projfile)
        ! progress
        call progressfile_complete_parts(params%nparts) 
        ! clean
        call spproj%kill
        call o_mic%kill
        call o_tmp%kill
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_EXTRACT_DISTR NORMAL STOP ****')

        contains

            character(len=LONGSTRLEN) function boxfile_from_mic(mic)
                character(len=*), intent(in) :: mic
                character(len=LONGSTRLEN)    :: box_from_mic
                integer :: ibox
                box_from_mic     = fname_new_ext(basename(mic),'box')
                boxfile_from_mic = NIL
                do ibox=1,size(boxfiles)
                    if(trim(basename(boxfiles(ibox))).eq.trim(box_from_mic))then
                        boxfile_from_mic = trim(boxfiles(ibox))
                        return
                    endif
                enddo
            end function boxfile_from_mic

    end subroutine exec_extract_distr

    subroutine exec_extract( self, cline )
        use simple_ctf,                 only: ctf
        use simple_ctf_estimate_fit,    only: ctf_estimate_fit
        use simple_particle_extractor,  only: ptcl_extractor
        class(extract_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline !< command line input
        type(image),                allocatable :: imgs(:)
        type(parameters)                        :: params
        type(ptcl_extractor)                    :: extractor
        type(sp_project)                        :: spproj_in, spproj
        type(nrtxtfile)                         :: boxfile
        type(image)                             :: micrograph, micrograph_pad
        type(oris)                              :: os_mic
        type(ori)                               :: o_mic, o_tmp
        type(ctf)                               :: tfun
        type(ctfparams)                         :: ctfparms
        type(ctf_estimate_fit)                  :: ctffit
        type(stack_io)                          :: stkio_w
        character(len=:),           allocatable :: output_dir, mic_name, imgkind
        real,                       allocatable :: boxdata(:,:)
        integer,                    allocatable :: ptcl_inds(:)
        logical,                    allocatable :: oris_mask(:), mics_mask(:)
        character(len=LONGSTRLEN) :: stack, boxfile_name, box_fname, ctfdoc
        character(len=STDLEN)     :: ext
        real                      :: ptcl_pos(2), stk_mean,stk_sdev,stk_max,stk_min,dfx,dfy,prog
        integer                   :: ldim(3), lfoo(3), fromto(2)
        integer                   :: nframes, imic, iptcl, nptcls,nmics,nmics_here,box, i, iptcl_g
        integer                   :: cnt, nmics_tot, ifoo, state, iptcl_glob, nptcls2extract
        logical                   :: l_ctfpatch, l_gid_present, l_ogid_present,prog_write,prog_part
        call cline%set('oritype', 'mic')
        call cline%set('mkdir',   'no')
        call params%new(cline)
        ! init
        output_dir = PATH_HERE
        fromto(:)  = [params%fromp, params%top]
        nmics_here = fromto(2)-fromto(1)+1
        prog_write = .false.
        prog_part  = .false.
        if( params%stream.eq.'yes' )then
            output_dir = DIR_EXTRACT
            if( cline%defined('dir') ) output_dir = trim(params%dir)//'/'
            ! read in integrated movies
            call spproj%read(params%projfile)
            nmics_tot = spproj%os_mic%get_noris()
            fromto(:)  = [1,1]
            nmics_here = 1
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto(:)  = [params%fromp, params%top]
                nmics_here = params%top-params%fromp+1
            endif
            if( nmics_tot /= nmics_here ) THROW_HARD('Incompatible # of integrated micrograph to process!')
        else
            ! read in integrated movies
            call spproj_in%read_segment(params%oritype, params%projfile)
            nmics_tot = spproj_in%os_mic%get_noris()
            if( spproj_in%get_nintgs() == 0 ) THROW_HARD('No integrated micrograph to process!')
            ! init output project
            call spproj%read_non_data_segments(params%projfile)
            call spproj%projinfo%set(1,'projname', get_fbody(params%outfile,METADATA_EXT,separator=.false.))
            call spproj%projinfo%set(1,'projfile', params%outfile)
            params%projfile = trim(params%outfile)
            call spproj%os_mic%new(nmics_here, is_ptcl=.false.)
            cnt = 0
            do imic = fromto(1),fromto(2)
                cnt = cnt + 1
                call spproj_in%os_mic%get_ori(imic, o_tmp)
                call spproj%os_mic%set_ori(cnt, o_tmp)
            enddo
            prog_write = .true.
            if( cline%defined('part') ) then 
                prog_part = .true.
                call progressfile_init_part(cline%get_iarg('part'))
            else
                call progressfile_init()
            endif
            call spproj_in%kill
        endif
        ! input boxes
        if( cline%defined('dir_box') )then
            if( .not.file_exists(params%dir_box) )then
                write(logfhandle,*)'Directory does not exist: ', trim(params%dir_box), 'simple_commander_preprocess::exec_extract'
                THROW_HARD('box directory does not exist; exec_extract')
            endif
        endif
        ! sanity checks
        allocate(mics_mask(1:nmics_here), source=.false.)
        nmics  = 0
        nptcls = 0
        do imic = 1,nmics_here
            call spproj%os_mic%get_ori(imic, o_mic)
            state = 1
            if( o_mic%isthere('state') ) state = o_mic%get_state()
            if( state == 0 ) cycle
            if( .not. o_mic%isthere('imgkind') )cycle
            if( .not. o_mic%isthere('intg')    )cycle
            call o_mic%getter('imgkind', imgkind)
            if( trim(imgkind).ne.'mic') cycle
            call o_mic%getter('intg', mic_name)
            if( .not.file_exists(mic_name) )cycle
            ! box input
            if( cline%defined('dir_box') )then
                box_fname = trim(params%dir_box)//'/'//fname_new_ext(basename(mic_name),'box')
                if( .not.file_exists(box_fname) )cycle
                boxfile_name = simple_abspath(box_fname, check_exists=.false.)
                call spproj%set_boxfile(imic, boxfile_name)
            else
                boxfile_name = trim(o_mic%get_static('boxfile'))
                if( .not.file_exists(boxfile_name) )cycle
            endif
            ! get number of frames from stack
            call find_ldim_nptcls(mic_name, lfoo, nframes )
            if( nframes > 1 ) THROW_HARD('multi-frame extraction not supported; exec_extract')
            ! update mask
            mics_mask(imic) = .true.
            nmics = nmics + 1
            ! image & box dimensions
            if( nmics == 1 )call find_ldim_nptcls(mic_name, ldim, ifoo)
            if( nptcls == 0 .and. .not.cline%defined('box') )then
                if( nlines(boxfile_name) > 0 )then
                    call boxfile%new(boxfile_name, 1)
                    nptcls = boxfile%get_ndatalines()
                endif
                if( nptcls == 0 )then
                    call spproj%os_mic%set(imic, 'nptcls', 0)
                    cycle
                endif
                allocate( boxdata(boxfile%get_nrecs_per_line(),nptcls) )
                call boxfile%readNextDataLine(boxdata(:,1))
                call boxfile%kill
                params%box = nint(boxdata(3,1))
            endif
        enddo
        ! actual extraction
        if( nmics == 0 )then
            ! done
        else
            if( params%box == 0 )THROW_HARD('box cannot be zero; exec_extract')
            ! init
            call micrograph%new([ldim(1),ldim(2),1], params%smpd)
            if( trim(params%extractfrommov).ne.'yes' ) call extractor%init_mic(params%box, (params%pcontrast .eq. 'black'))
            ! main loop
            iptcl_glob = 0 ! extracted particle index among ALL stacks
            prog = 0.0
            do imic = 1,nmics_here
                if( .not.mics_mask(imic) )then
                    call spproj%os_mic%set(imic, 'nptcls', 0)
                    call spproj%os_mic%set_state(imic, 0)
                    cycle
                endif
                ! fetch micrograph
                call spproj%os_mic%get_ori(imic, o_mic)
                boxfile_name = trim(o_mic%get_static('boxfile'))
                ! box file
                nptcls = 0
                if( nlines(boxfile_name) > 0 )then
                    call boxfile%new(boxfile_name, 1)
                    nptcls = boxfile%get_ndatalines()
                endif
                if( nptcls == 0 ) cycle
                call progress(imic,nmics_tot)
                ! box checks
                if(allocated(oris_mask))deallocate(oris_mask)
                allocate(oris_mask(nptcls), source=.false.)
                ! read box data & update mask
                if(allocated(boxdata))deallocate(boxdata)
                allocate( boxdata(boxfile%get_nrecs_per_line(),nptcls))
                do iptcl=1,nptcls
                    call boxfile%readNextDataLine(boxdata(:,iptcl))
                    box = nint(boxdata(3,iptcl))
                    if( nint(boxdata(3,iptcl)) /= nint(boxdata(4,iptcl)) )then
                        THROW_HARD('only square windows allowed; exec_extract')
                    endif
                    ! modify coordinates if change in box (shift by half the difference)
                    if( box /= params%box ) boxdata(1:2,iptcl) = boxdata(1:2,iptcl) - real(params%box-box)/2.
                    if( .not.cline%defined('box') .and. nint(boxdata(3,iptcl)) /= params%box )then
                        write(logfhandle,*) 'box_current: ', nint(boxdata(3,iptcl)), 'box in params: ', params%box
                        THROW_HARD('inconsistent box sizes in box files; exec_extract')
                    endif
                    ! update particle mask & movie index
                    oris_mask(iptcl)  = (trim(params%outside).eq.'yes') .or. box_inside(ldim, nint(boxdata(1:2,iptcl)), params%box)
                end do
                ! update micrograph field
                nptcls2extract = count(oris_mask)
                call spproj%os_mic%set(imic, 'nptcls', nptcls2extract)
                if( nptcls2extract == 0 )then
                    ! no particles to extract
                    mics_mask(imic) = .false.
                    cycle
                endif
                if(allocated(ptcl_inds))deallocate(ptcl_inds)
                allocate(ptcl_inds(nptcls2extract), source=0)
                cnt = 0
                do iptcl=1,nptcls
                    if( oris_mask(iptcl) )then
                        cnt = cnt + 1
                        ptcl_inds(cnt) = iptcl
                    endif
                enddo
                ! fetch ctf info
                ctfparms      = o_mic%get_ctfvars()
                ctfparms%smpd = params%smpd
                if( o_mic%isthere('dfx') )then
                    if( .not.o_mic%isthere('cs') .or. .not.o_mic%isthere('kv') .or. .not.o_mic%isthere('fraca') )then
                        THROW_HARD('input lacks at least cs, kv or fraca; exec_extract')
                    endif
                endif
                ! output stack name
                call o_mic%getter('intg', mic_name)
                ext   = fname2ext(trim(basename(mic_name)))
                stack = trim(output_dir)//trim(EXTRACT_STK_FBODY)//trim(get_fbody(trim(basename(mic_name)), trim(ext)))//trim(STK_EXT)
                ! init extraction
                call prepimgbatch(nptcls2extract)
                if( trim(params%extractfrommov).eq.'yes' )then
                    ! extraction from movie
                    if( trim(params%ctf).eq.'flip' .and. o_mic%isthere('dfx') )then
                        THROW_HARD('extractfrommov=yes does not support ctf=flip yet')
                    endif
                    call extractor%init_mov(o_mic, params%box, (params%pcontrast .eq. 'black'))
                    call extractor%extract_particles(ptcl_inds, nint(boxdata), imgs, stk_min,stk_max,stk_mean,stk_sdev)
                else
                    ! extraction from micrograph
                    call micrograph%read(mic_name, 1)
                    if( trim(params%backgr_subtr).eq.'yes') call micrograph%subtract_background(HP_BACKGR_SUBTR)
                    ! phase-flip micrograph
                    if( cline%defined('ctf') )then
                        if( trim(params%ctf).eq.'flip' .and. o_mic%isthere('dfx') )then
                            tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                            call micrograph%zero_edgeavg
                            call micrograph%fft
                            call tfun%apply_serial(micrograph, 'flip', ctfparms)
                            call micrograph%ifft
                            ! update stack ctf flag, mic flag unchanged
                            ctfparms%ctfflag = CTFFLAG_FLIP
                        endif
                    endif
                    ! extraction
                    call extractor%extract_particles_from_mic(micrograph, ptcl_inds, nint(boxdata), imgs,&
                        &stk_min,stk_max,stk_mean,stk_sdev)
                endif
                ! write stack
                call stkio_w%open(trim(adjustl(stack)), params%smpd, 'write', box=params%box)
                do i = 1,nptcls2extract
                    call stkio_w%write(i, imgs(i))
                enddo
                call stkio_w%close
                ! update stack stats
                call imgs(1)%update_header_stats(trim(adjustl(stack)), [stk_min, stk_max, stk_mean, stk_sdev])
                ! IMPORT INTO PROJECT (with absolute path)
                call spproj%add_stk(trim(adjustl(stack)), ctfparms)
                ! add box coordinates to ptcl2D field only & updates patch-based defocus
                l_ctfpatch = .false.
                if( o_mic%isthere('ctfdoc') )then
                    ctfdoc = o_mic%get_static('ctfdoc')
                    if( file_exists(ctfdoc) )then
                        call ctffit%read_doc(ctfdoc)
                        l_ctfpatch = .true.
                    endif
                endif
                l_ogid_present = o_mic%isthere('ogid')
                l_gid_present  = o_mic%isthere('gid')
                !$omp parallel do schedule(static) default(shared) proc_bind(close)&
                !$omp private(i,iptcl,iptcl_g,ptcl_pos,dfx,dfy)
                do i = 1,nptcls2extract
                    iptcl    = ptcl_inds(i)
                    iptcl_g  = iptcl_glob + i
                    ptcl_pos = boxdata(1:2,iptcl)
                    ! updates particle position
                    call spproj%set_boxcoords(iptcl_g, nint(ptcl_pos))
                    ! updates particle defocus
                    if( l_ctfpatch )then
                        ptcl_pos = ptcl_pos+1.+real(params%box/2) !  center
                        call ctffit%pix2polyvals(ptcl_pos(1),ptcl_pos(2), dfx,dfy)
                        call spproj%os_ptcl2D%set_dfx(iptcl_g,dfx)
                        call spproj%os_ptcl2D%set_dfy(iptcl_g,dfy)
                        call spproj%os_ptcl3D%set_dfx(iptcl_g,dfx)
                        call spproj%os_ptcl3D%set_dfy(iptcl_g,dfy)
                    endif
                    ! update particle optics group id
                    if( l_ogid_present )then
                        call spproj%os_ptcl2D%set(iptcl_g,'ogid',o_mic%get('ogid'))
                        call spproj%os_ptcl3D%set(iptcl_g,'ogid',o_mic%get('ogid'))
                    endif
                    ! update particle group id
                    if( l_gid_present )then
                        call spproj%os_ptcl2D%set(iptcl_g,'gid',o_mic%get('gid'))
                        call spproj%os_ptcl3D%set(iptcl_g,'gid',o_mic%get('gid'))
                    endif
                end do
                !$omp end parallel do
                ! global particle count
                iptcl_glob = iptcl_glob + nptcls2extract
                ! clean
                call boxfile%kill
                call ctffit%kill
                ! progress
                if(prog_write) then
                    if( (real(imic) / real(nmics_here)) > prog + 0.05 ) then
                        prog = real(imic) / real(nmics_here)
                        if(prog_part) then 
                            call progressfile_update_part(cline%get_iarg('part'), prog)
                        else
                            call progressfile_update(prog)
                        endif
                    endif
                endif
            enddo
            call killimgbatch
        endif
        ! write
        if( trim(params%stream).eq.'yes' )then
            ! purging state=0 and nptcls=0 mics such that all mics (nmics>1)
            ! can be assumed associated with particles in streaming
            nmics = count(mics_mask)
            if( nmics == 0 )then
                call spproj%os_mic%kill
                call spproj%os_stk%kill
                call spproj%os_ptcl2D%kill
                call spproj%os_ptcl3D%kill
            else
                if( nmics < nmics_here )then
                    call os_mic%new(nmics, is_ptcl=.false.)
                    cnt = 0
                    do imic = 1, nmics_here
                        if( mics_mask(imic) )then
                            cnt = cnt+1
                            call os_mic%transfer_ori(cnt, spproj%os_mic, imic)
                        endif
                    enddo
                    spproj%os_mic = os_mic
                endif
            endif
        endif
        call spproj%write(params%projfile)
        ! end gracefully
        call extractor%kill
        call micrograph%kill
        call micrograph_pad%kill
        call o_mic%kill
        call o_tmp%kill
        call os_mic%kill
        if( prog_write ) call progressfile_update(1.0)
        call qsys_job_finished('simple_commander_preprocess :: exec_extract')
        call simple_end('**** SIMPLE_EXTRACT NORMAL STOP ****')
        contains

            subroutine prepimgbatch( batchsz )
                integer,           intent(in) :: batchsz
                integer :: i
                logical :: doprep
                doprep = .false.
                if( .not. allocated(imgs) )then
                    doprep = .true.
                else
                    if( batchsz > size(imgs) ) doprep = .true.
                    if( doprep ) call killimgbatch
                endif
                if( doprep )then
                    allocate(imgs(batchsz))
                    !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
                    do i = 1,batchsz
                        call imgs(i)%new([params%box, params%box, 1], params%smpd, wthreads=.false.)
                    end do
                    !$omp end parallel do
                endif
            end subroutine prepimgbatch

            subroutine killimgbatch
                integer :: i
                if( allocated(imgs) )then
                    do i = 1,size(imgs)
                        call imgs(i)%kill
                    end do
                    deallocate(imgs)
                endif
            end subroutine killimgbatch

    end subroutine exec_extract

    subroutine exec_reextract_distr( self, cline )
        class(reextract_commander_distr), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline !< command line input
        type(parameters)                        :: params
        type(sp_project)                        :: spproj
        type(sp_project),           allocatable :: spproj_parts(:)
        type(qsys_env)                          :: qenv
        type(chash)                             :: job_descr
        type(ori)                               :: o_mic, o
        type(oris)                              :: os_stk
        type(chash),                allocatable :: part_params(:)
        character(len=LONGSTRLEN),  allocatable :: boxfiles(:), stktab(:), parts_fname(:)
        character(len=:),           allocatable :: mic_name, imgkind
        integer,                    allocatable :: parts(:,:)
        integer :: imic,i,nmics_tot,numlen,nmics,cnt,state,istk,nstks,ipart,stkind,nptcls
        if( cline%defined('ctf') )then
            if( cline%get_carg('ctf').ne.'flip' .and. cline%get_carg('ctf').ne.'no' )then
                THROW_HARD('Only CTF=NO/FLIP are allowed')
            endif
        endif
        if( cline%defined('osmpd') )then
            if( .not.cline%defined('box') ) THROW_HARD('BOX must be defined with OSMPD!')
        endif
        if( .not. cline%defined('mkdir')     )     call cline%set('mkdir',          'yes')
        if( .not. cline%defined('pcontrast') )     call cline%set('pcontrast',    'black')
        if( .not. cline%defined('oritype')   )     call cline%set('oritype',     'ptcl3D')
        if( .not. cline%defined('extractfrommov')) call cline%set('extractfrommov',  'no')
        if( .not. cline%defined('backgr_subtr'))   call cline%set('backgr_subtr',    'no')
        call params%new(cline)
        call cline%set('mkdir', 'no')
        ! read in integrated movies
        call spproj%read( params%projfile )
        if( spproj%get_nintgs() == 0 ) THROW_HARD('No integrated micrograph to process!')
        if( spproj%get_nstks() == 0 ) THROW_HARD('This project file does not contain stacks!')
        nmics_tot = spproj%os_mic%get_noris()
        if( nmics_tot < params%nparts )then
            params%nparts = nmics_tot
        endif
        ! sanity checks
        nmics  = 0
        do imic = 1, nmics_tot
            call spproj%os_mic%get_ori(imic, o_mic)
            state = 1
            if( o_mic%isthere('state') ) state = o_mic%get_state()
            if( state == 0 ) cycle
            if( .not. o_mic%isthere('imgkind') )cycle
            if( .not. o_mic%isthere('intg')    )cycle
            call o_mic%getter('imgkind', imgkind)
            if( trim(imgkind).ne.'mic') cycle
            call o_mic%getter('intg', mic_name)
            if( .not.file_exists(mic_name) )cycle
            ! update counter
            nmics = nmics + 1
            ! removes boxfile from micrographs
            call spproj%os_mic%delete_entry(imic,'boxfile')
        enddo
        if( nmics == 0 )then
            THROW_WARN('No particles to re-extract! exec_reextract')
            return
        endif
        call spproj%os_mic%kill
        call spproj%os_stk%kill
        call spproj%os_ptcl2D%kill
        call spproj%os_ptcl3D%kill
        ! DISTRIBUTED EXTRACTION
        ! setup the environment for distributed execution
        parts = split_nobjs_even(nmics_tot, params%nparts)
        allocate(part_params(params%nparts))
        do ipart=1,params%nparts
            call part_params(ipart)%new(2)
            call part_params(ipart)%set('fromp',int2str(parts(ipart,1)))
            call part_params(ipart)%set('top',  int2str(parts(ipart,2)))
        end do
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY), part_params=part_params, array=L_USE_SLURM_ARR)
        ! ASSEMBLY
        allocate(spproj_parts(params%nparts),parts_fname(params%nparts))
        numlen = len(int2str(params%nparts))
        do ipart = 1,params%nparts
            parts_fname(ipart) = trim(ALGN_FBODY)//int2str_pad(ipart,numlen)//trim(METADATA_EXT)
        enddo
        ! copy updated micrographs
        cnt   = 0
        nmics = 0
        do ipart = 1,params%nparts
            call spproj_parts(ipart)%read_segment('mic',parts_fname(ipart))
            nmics = nmics + spproj_parts(ipart)%os_mic%get_noris()
        enddo
        if( nmics > 0 )then
            call spproj%os_mic%new(nmics, is_ptcl=.false.)
            ! transfer stacks
            cnt   = 0
            nstks = 0
            do ipart = 1,params%nparts
                do imic = 1,spproj_parts(ipart)%os_mic%get_noris()
                    cnt = cnt + 1
                    call spproj%os_mic%transfer_ori(cnt, spproj_parts(ipart)%os_mic, imic)
                enddo
                call spproj_parts(ipart)%kill
                call spproj_parts(ipart)%read_segment('stk',parts_fname(ipart))
                nstks = nstks + spproj_parts(ipart)%os_stk%get_noris()
            enddo
            if( nstks /= nmics ) THROW_HARD('Inconstistent number of stacks in individual projects')
            ! generates stacks table
            call os_stk%new(nstks, is_ptcl=.false.)
            allocate(stktab(nstks))
            cnt = 0
            do ipart = 1,params%nparts
                do istk = 1,spproj_parts(ipart)%os_stk%get_noris()
                    cnt = cnt + 1
                    call os_stk%transfer_ori(cnt, spproj_parts(ipart)%os_stk, istk)
                    stktab(cnt) = os_stk%get_static(cnt,'stk')
                enddo
                call spproj_parts(ipart)%kill
            enddo
            ! import stacks into project
            call spproj%add_stktab(stktab,os_stk)
            call os_stk%kill
            ! 2D/3D parameters, transfer everything but stack index
            cnt = 0
            do ipart = 1,params%nparts
                call spproj_parts(ipart)%read_segment('ptcl2D',parts_fname(ipart))
                call spproj_parts(ipart)%read_segment('ptcl3D',parts_fname(ipart))
                nptcls = spproj_parts(ipart)%os_ptcl2D%get_noris()
                if( nptcls /= spproj_parts(ipart)%os_ptcl3D%get_noris())then
                    THROW_HARD('Inconsistent number of particles')
                endif
                do i = 1,nptcls
                    cnt    = cnt + 1
                    stkind = spproj%os_ptcl2D%get_int(cnt,'stkind')
                    call spproj%os_ptcl2D%transfer_ori(cnt, spproj_parts(ipart)%os_ptcl2D, i)
                    call spproj%os_ptcl3D%transfer_ori(cnt, spproj_parts(ipart)%os_ptcl3D, i)
                    call spproj%os_ptcl2D%set(cnt,'stkind',stkind)
                    call spproj%os_ptcl3D%set(cnt,'stkind',stkind)
                enddo
                call spproj_parts(ipart)%kill
            enddo
        endif
        ! final write
        call spproj%write( params%projfile )
        ! clean-up
        call qsys_cleanup
        call spproj%kill
        deallocate(spproj_parts,part_params)
        call o_mic%kill
        call o%kill
        ! end gracefully
        call simple_end('**** SIMPLE_REEXTRACT_DISTR NORMAL STOP ****')
        contains

            character(len=LONGSTRLEN) function boxfile_from_mic(mic)
                character(len=*), intent(in) :: mic
                character(len=LONGSTRLEN)    :: box_from_mic
                integer :: ibox
                box_from_mic     = fname_new_ext(basename(mic),'box')
                boxfile_from_mic = NIL
                do ibox=1,size(boxfiles)
                    if(trim(basename(boxfiles(ibox))).eq.trim(box_from_mic))then
                        boxfile_from_mic = trim(boxfiles(ibox))
                        return
                    endif
                enddo
            end function boxfile_from_mic

    end subroutine exec_reextract_distr

    subroutine exec_reextract( self, cline )
        use simple_ctf,                 only: ctf
        use simple_strategy2D3D_common, only: prepimgbatch, killimgbatch
        use simple_particle_extractor,  only: ptcl_extractor
        class(reextract_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline !< command line input
        type(parameters)              :: params
        type(sp_project)              :: spproj, spproj_in
        type(builder)                 :: build
        type(image)                   :: micrograph, micrograph_sc
        type(ori)                     :: o_mic, o_stk
        type(ctf)                     :: tfun
        type(ctfparams)               :: ctfparms
        type(stack_io)                :: stkio_w
        type(ptcl_extractor)          :: extractor
        character(len=:), allocatable :: mic_name, imgkind, ext
        logical,          allocatable :: mic_mask(:), ptcl_mask(:)
        integer,          allocatable :: mic2stk_inds(:), boxcoords(:,:), ptcl_inds(:), mic_dims(:,:)
        character(len=LONGSTRLEN)     :: stack
        real    :: prev_shift(2),shift2d(2),shift3d(2),prev_shift_sc(2), translation(2), prev_center_sc(2)
        real    :: stk_min, stk_max, stk_mean, stk_sdev, scale
        integer :: prev_pos(2), new_pos(2), ishift(2), ldim(3), prev_center(2), new_center(2), ldim_sc(3)
        integer :: i, nframes, imic, iptcl, nmics, prev_box, box_foo, cnt, nmics_tot, stk_ind
        integer :: fromp, top, istk, nptcls2extract, nptcls
        logical :: l_3d, l_scale_particles, l_movie_frames
        call cline%set('mkdir','no')
        call params%new(cline)
        l_movie_frames    = trim(params%extractfrommov).eq.'yes'
        l_scale_particles = cline%defined('osmpd')
        if( l_scale_particles )then
            if( .not.cline%defined('box') ) THROW_HARD('BOX must be defined with OSMPD!')
            if( l_movie_frames ) THROW_HARD('Particle scaling and extraction of movie frames is not supported!')
        endif
        if( l_movie_frames .and. (trim(params%ctf).eq.'flip') )then
            THROW_HARD('extractfrommov=yes does not support ctf=flip!')
        endif
        ! set normalization radius
        params%msk = RADFRAC_NORM_EXTRACT * real(params%box/2)
        ! whether to use shifts from 2D or 3D
        l_3d = .true.
        if(cline%defined('oritype')) l_3d = trim(params%oritype)=='ptcl3D'
        ! read in integrated movies
        call spproj_in%read_segment('mic', params%projfile)
        nmics_tot = spproj_in%os_mic%get_noris()
        if( spproj_in%get_nintgs() == 0 ) THROW_HARD('No integrated micrograph to process!')
        call spproj_in%read_segment('stk', params%projfile)
        ! sanity checks, dimensions & indexing
        box_foo  = 0
        prev_box = 0
        ldim     = 0
        allocate(mic2stk_inds(nmics_tot), source=0)
        allocate(mic_mask(nmics_tot),     source=.false.)
        stk_ind = 0
        do imic = 1,nmics_tot
            if( imic > params%top ) exit
            call spproj_in%os_mic%get_ori(imic, o_mic)
            if( o_mic%isthere('state') )then
                if( o_mic%get_state() == 0 )cycle
            endif
            if( .not. o_mic%isthere('imgkind') )cycle
            if( .not. o_mic%isthere('intg')    )cycle
            call o_mic%getter('imgkind', imgkind)
            if( trim(imgkind).ne.'mic') cycle
            ! find next selected stack
            do istk=stk_ind,spproj_in%os_stk%get_noris()
                stk_ind = stk_ind+1
                if( spproj_in%os_stk%isthere(stk_ind,'state') )then
                    if( spproj_in%os_stk%get_state(stk_ind) == 1 ) exit
                else
                    exit
                endif
            enddo
            ! update index & mask
            if( imic>=params%fromp .and. imic<=params%top )then
                mic_mask(imic) = .true.
                mic2stk_inds(imic) = stk_ind ! index to os_stk
            endif
        enddo
        nmics = count(mic_mask)
        if( nmics > 0 )then
            call build%build_general_tbox(params, cline, do3d=.false.)
            allocate(mic_dims(3,nmics_tot),source=0)
            ! sanity checks
            do imic = 1,nmics_tot
                if( .not.mic_mask(imic) )cycle
                ! sanity checks
                call spproj_in%os_mic%get_ori(imic, o_mic)
                call o_mic%getter('intg', mic_name)
                if( .not.file_exists(mic_name) )cycle
                 ! micrograph dimensions
                call find_ldim_nptcls(mic_name, ldim, nframes )
                if( nframes > 1 ) THROW_HARD('multi-frame extraction not supported; exec_reextract')
                mic_dims(:,imic) = [ldim(1),ldim(2),1]
                if( l_scale_particles )then
                    ! the following checks are not performed
                else
                    stk_ind = mic2stk_inds(imic)
                    call spproj_in%os_stk%get_ori(stk_ind, o_stk)
                    box_foo = o_stk%get_int('box')
                    if( prev_box == 0 ) prev_box = box_foo
                    if( prev_box /= box_foo ) THROW_HARD('Inconsistent box size; exec_reextract')
                endif
            enddo
            if( .not.cline%defined('box') ) params%box = prev_box
            if( is_odd(params%box) ) THROW_HARD('Box size must be of even dimension! exec_extract')
            ! extraction
            write(logfhandle,'(A)')'>>> EXTRACTING... '
            call spproj_in%read_segment('ptcl2D', params%projfile)
            call spproj_in%read_segment('ptcl3D', params%projfile)
            allocate(ptcl_mask(spproj_in%os_ptcl2D%get_noris()),source=.false.)
            ldim = 0
            do imic = params%fromp,params%top
                if( .not.mic_mask(imic) ) cycle
                ! init micrograph object and extractor
                call spproj_in%os_mic%get_ori(imic, o_mic)
                call o_mic%getter('intg', mic_name)
                ctfparms = o_mic%get_ctfvars()
                if( any(ldim /= mic_dims(:,imic)) )then
                    ! first iteration or different micrograph size
                    ldim = mic_dims(:,imic)
                    call micrograph%new(ldim, ctfparms%smpd)
                    if( .not.l_movie_frames ) call extractor%init_mic(params%box, (params%pcontrast .eq. 'black'))
                    if( l_scale_particles )then
                        scale        = ctfparms%smpd / params%osmpd
                        ldim_sc(1:2) = round2even(scale*real(ldim(1:2)))
                        ldim_sc(3)   = 1
                        call micrograph_sc%new(ldim_sc, params%osmpd)
                    endif
                endif
                ! stack
                stk_ind = mic2stk_inds(imic)
                call spproj_in%os_stk%get_ori(stk_ind, o_stk)
                prev_box = o_stk%get_int('box')
                fromp    = o_stk%get_fromp()
                top      = o_stk%get_top()
                ext      = fname2ext(trim(basename(mic_name)))
                stack    = trim(EXTRACT_STK_FBODY)//trim(get_fbody(trim(basename(mic_name)), trim(ext)))//trim(STK_EXT)
                ! updating shifts, positions, states and doc
                if( allocated(boxcoords) ) deallocate(boxcoords)
                allocate(boxcoords(2,fromp:top),source=0)
                !$omp parallel do default(shared) proc_bind(close) schedule(static)&
                !$omp private(iptcl,prev_pos,prev_shift,prev_center,prev_center_sc,prev_shift_sc,new_center)&
                !$omp private(new_pos,translation,shift2d,shift3d,ishift)
                do iptcl = fromp,top
                    if( spproj_in%os_ptcl2D%get_state(iptcl) == 0 ) cycle
                    if( spproj_in%os_ptcl3D%get_state(iptcl) == 0 ) cycle
                    ! previous position & shift
                    call spproj_in%get_boxcoords(iptcl, prev_pos)
                    if( l_3d )then
                        prev_shift = spproj_in%os_ptcl3D%get_2Dshift(iptcl)
                    else
                        prev_shift = spproj_in%os_ptcl2D%get_2Dshift(iptcl)
                    endif
                    if( l_scale_particles )then
                        ! scale center, shift & positions
                        prev_center      = prev_pos + prev_box/2
                        prev_center_sc   = scale * real(prev_center)
                        prev_shift_sc    = scale * real(prev_shift)
                        new_center       = nint(prev_center_sc - prev_shift_sc)
                        new_pos          = new_center - params%box/2
                        translation      = -(prev_center_sc - real(new_center))
                        ptcl_mask(iptcl) = box_inside(ldim, new_pos, params%box)
                        if( ptcl_mask(iptcl) )then
                            ! updates shifts
                            if( l_3d )then
                                shift2d = scale * spproj_in%os_ptcl2D%get_2Dshift(iptcl) + translation
                                shift3d = scale * prev_shift                             + translation
                            else
                                shift2d = scale * prev_shift                             + translation
                                shift3d = scale * spproj_in%os_ptcl3D%get_2Dshift(iptcl) + translation
                            endif
                        endif
                    else
                        ! calc new position & shift
                        ishift      = nint(prev_shift)
                        new_pos     = prev_pos - ishift
                        translation = -real(ishift)
                        if( prev_box /= params%box ) new_pos = new_pos + (prev_box-params%box)/2
                        ptcl_mask(iptcl) = box_inside(ldim, new_pos, params%box)
                        if( ptcl_mask(iptcl) )then
                            ! updates shifts
                            if( l_3d )then
                                shift2d = spproj_in%os_ptcl2D%get_2Dshift(iptcl) + translation
                                shift3d = prev_shift                             + translation
                            else
                                shift2d = prev_shift                             + translation
                                shift3d = spproj_in%os_ptcl3D%get_2Dshift(iptcl) + translation
                            endif
                        endif
                    endif
                    ! updates document
                    if( ptcl_mask(iptcl) )then
                        ! updates picking position
                        call spproj_in%set_boxcoords(iptcl, new_pos)
                        ! updates shifts
                        call spproj_in%os_ptcl2D%set_shift(iptcl, shift2d)
                        call spproj_in%os_ptcl3D%set_shift(iptcl, shift3d)
                    else
                        ! excluded
                        call spproj_in%os_ptcl2D%set_state(iptcl, 0)
                        call spproj_in%os_ptcl3D%set_state(iptcl, 0)
                    endif
                    ! for actual extraction
                    boxcoords(:,iptcl) = new_pos
                enddo
                !$omp end parallel do
                nptcls2extract = count(ptcl_mask(fromp:top))
                if( nptcls2extract > 0 )then
                    if( allocated(ptcl_inds) ) deallocate(ptcl_inds)
                    allocate(ptcl_inds(nptcls2extract),source=0)
                    cnt = 0
                    do iptcl = fromp,top
                        if( .not.ptcl_mask(iptcl) ) cycle
                        cnt = cnt + 1
                        ptcl_inds(cnt) = iptcl
                        ! updating index of particle in stack
                        call spproj_in%os_ptcl2D%set(iptcl, 'indstk', cnt)
                        call spproj_in%os_ptcl3D%set(iptcl, 'indstk', cnt)
                    enddo
                    ptcl_inds = ptcl_inds -fromp+1 ! because indexing range lost when passed to extractor
                    call prepimgbatch(nptcls2extract)
                    if( l_movie_frames )then
                        ! extraction from movie
                        call extractor%init_mov(o_mic, params%box, (params%pcontrast .eq. 'black'))
                        call extractor%extract_particles(ptcl_inds, boxcoords, build%imgbatch, stk_min,stk_max,stk_mean,stk_sdev)
                    else
                        ! read micrograph
                        call micrograph%read(mic_name)
                        ! preprocess micrograph
                        if( trim(params%backgr_subtr).eq.'yes') call micrograph%subtract_background(HP_BACKGR_SUBTR)
                        if( ctfparms%ctfflag == CTFFLAG_FLIP )then
                            if( o_mic%isthere('dfx') )then
                                ! phase flip micrograph
                                tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                                call micrograph%zero_edgeavg
                                call micrograph%fft
                                call tfun%apply_serial(micrograph, 'flip', ctfparms)
                                call micrograph%ifft
                            endif
                        endif
                        ! Actual extraction
                        if( l_scale_particles )then
                            ! scale
                            if( all(ldim_sc == ldim) )then
                                call micrograph_sc%copy_fast(micrograph)
                            else
                                call micrograph%fft
                                if( any(ldim_sc < ldim) )then
                                    call micrograph%clip(micrograph_sc)
                                else
                                    call micrograph%pad(micrograph_sc, antialiasing=.false.)
                                endif
                                call micrograph_sc%ifft
                                call micrograph%set_ft(.false.)
                            endif
                            call micrograph_sc%set_smpd(params%osmpd) ! safety
                            ! extract
                            call extractor%extract_particles_from_mic(micrograph_sc, ptcl_inds, boxcoords, build%imgbatch,&
                            &stk_min,stk_max,stk_mean,stk_sdev)
                        else
                            ! extract
                            call extractor%extract_particles_from_mic(micrograph, ptcl_inds, boxcoords, build%imgbatch,&
                                &stk_min,stk_max,stk_mean,stk_sdev)
                        endif
                    endif
                    ! write stack
                    if( l_scale_particles )then
                        call stkio_w%open(trim(adjustl(stack)), params%osmpd, 'write', box=params%box)
                    else
                        call stkio_w%open(trim(adjustl(stack)), params%smpd, 'write', box=params%box)
                    endif
                    do i = 1,nptcls2extract
                        call stkio_w%write(i, build%imgbatch(i))
                    enddo
                    call stkio_w%close
                    if( l_scale_particles )then
                        call spproj_in%os_stk%set(stk_ind,'smpd',params%osmpd) !!
                        call micrograph_sc%update_header_stats(trim(adjustl(stack)), [stk_min, stk_max, stk_mean, stk_sdev])
                    else
                        call micrograph%update_header_stats(trim(adjustl(stack)), [stk_min, stk_max, stk_mean, stk_sdev])
                    endif
                    call spproj_in%os_stk%set(stk_ind,'stk',   simple_abspath(stack,check_exists=.false.))
                    call spproj_in%os_stk%set(stk_ind,'box',   params%box)
                    call spproj_in%os_stk%set(stk_ind,'nptcls',nptcls2extract)
                    call spproj_in%os_mic%set(imic,   'nptcls',nptcls2extract)
                    call spproj_in%os_mic%delete_entry(imic,'boxfile')
                else
                    ! all particles in this micrograph excluded
                    call spproj_in%os_stk%set(stk_ind,'state',0)
                    call spproj_in%os_mic%set(imic,'state',0)
                    mic_mask(imic) = .false.
                    mic2stk_inds(imic) = 0
                endif
            enddo
        endif
        call micrograph%kill
        call micrograph_sc%kill
        call extractor%kill
        call killimgbatch
        ! OUTPUT
        call spproj%read_non_data_segments(params%projfile)
        call spproj%projinfo%set(1,'projname', get_fbody(params%outfile,METADATA_EXT,separator=.false.))
        call spproj%projinfo%set(1,'projfile', params%outfile)
        nmics = count(mic_mask)
        ! transfer mics & stk
        call spproj%os_mic%new(nmics, is_ptcl=.false.)
        call spproj%os_stk%new(nmics, is_ptcl=.false.)
        nptcls = count(ptcl_mask)
        cnt = 0
        do imic = params%fromp,params%top
            if( .not.mic_mask(imic) )cycle
            cnt = cnt+1
            call spproj%os_mic%transfer_ori(cnt, spproj_in%os_mic, imic)
            stk_ind = mic2stk_inds(imic)
            call spproj%os_stk%transfer_ori(cnt, spproj_in%os_stk, stk_ind)
        enddo
        ! transfer particles
        nptcls = count(ptcl_mask)
        call spproj%os_ptcl2D%new(nptcls, is_ptcl=.true.)
        call spproj%os_ptcl3D%new(nptcls, is_ptcl=.true.)
        cnt = 0
        do iptcl = 1,size(ptcl_mask)
            if( .not.ptcl_mask(iptcl) )cycle
            cnt = cnt+1
            call spproj%os_ptcl2D%transfer_ori(cnt, spproj_in%os_ptcl2D, iptcl)
            call spproj%os_ptcl3D%transfer_ori(cnt, spproj_in%os_ptcl3D, iptcl)
        enddo
        call spproj_in%kill
        ! final write
        call spproj%write(params%outfile)
        write(logfhandle,'(A,I8)')'>>> RE-EXTRACTED  PARTICLES: ', nptcls
        ! end gracefully
        call qsys_job_finished('simple_commander_preprocess :: exec_reextract')
        call build%kill_general_tbox
        call o_mic%kill
        call o_stk%kill
        call simple_end('**** SIMPLE_REEXTRACT NORMAL STOP ****')
    end subroutine exec_reextract

    ! Stream only application
    subroutine exec_pick_extract( self, cline )
        use simple_sp_project,  only: sp_project
        use simple_picker_iter, only: picker_iter
        class(pick_extract_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters)              :: params
        type(oris)                    :: os_mic
        type(ori)                     :: o_mic
        type(picker_iter)             :: piter
        type(extract_commander)       :: xextract
        type(cmdline)                 :: cline_extract
        type(sp_project)              :: spproj
        character(len=:), allocatable :: micname, output_dir_picker, fbody, output_dir_extract
        character(len=LONGSTRLEN)     :: boxfile
        real    :: moldiam_opt
        integer :: fromto(2), imic, ntot, state, nvalid, i, nptcls
        logical :: l_extract, l_multipick
        ! set oritype
        call cline%set('oritype', 'mic')
        ! parse parameters
        call params%new(cline)
        if( params%stream.ne.'yes' ) THROW_HARD('new streaming only application')
        l_extract   = trim(params%extract).eq.'yes'
        l_multipick = params%nmoldiams > 1
        ! read in movies
        call spproj%read( params%projfile )
        if( spproj%get_nintgs() == 0 ) THROW_HARD('no micrograph to process!')
        params%smpd = spproj%os_mic%get(1,'smpd')
        call cline%set('smpd',params%smpd)
        ! output directories
        output_dir_picker  = trim(DIR_PICKER)
        if( l_extract ) output_dir_extract = trim(DIR_EXTRACT)
        if( cline%defined('dir') )then
            output_dir_picker  = filepath(params%dir,output_dir_picker)//'/'
            if( l_extract) output_dir_extract = filepath(params%dir,output_dir_extract)//'/'
        endif
        call simple_mkdir(output_dir_picker, errmsg="commander_pick_extract; ")
        if( l_extract ) call simple_mkdir(output_dir_extract,errmsg="commander_pick_extract; ")
        ! picker specs
        if( l_multipick )then
            ! nothing to do
        else
            select case(trim(params%picker))
                case('old')
                    if(.not.cline%defined('pickrefs')) THROW_HARD('PICKREFS required for picker=old')
                case('new')
                    if(cline%defined('pickrefs'))then
                    else
                        if( .not.cline%defined('moldiam') )then
                            THROW_HARD('MOLDIAM required for picker=new')
                        endif
                    endif
                case DEFAULT
                    THROW_HARD('Unsupported PICKER: '//trim(params%picker))
            end select
        endif
        ! command lines
        if( l_extract )then
            cline_extract = cline
            call cline_extract%set('dir', trim(output_dir_extract))
            call cline_extract%set('pcontrast', params%pcontrast)
            if( cline%defined('box_extract') ) call cline_extract%set('box', params%box_extract)
            call cline%delete('box')
            call cline_extract%delete('box_extract')
        endif
        ! file name
        if( cline%defined('fbody') )then
            fbody = trim(params%fbody)
        else
            fbody = ''
        endif
        ! range
        fromto(:) = 1
        if( cline%defined('fromp') .and. cline%defined('top') )then
            fromto(:) = [params%fromp, params%top]
        endif
        ntot   = fromto(2) - fromto(1) + 1
        nvalid = 0
        ! main loop
        do imic = fromto(1),fromto(2)
            ! fetch movie orientation
            call spproj%os_mic%get_ori(imic, o_mic)
            ! sanity check
            state = 1
            if( o_mic%isthere('state') ) state = o_mic%get_state()
            if( state == 0 ) cycle
            if( .not.o_mic%isthere('intg')   )cycle
            call o_mic%getter('intg', micname)
            if( .not.file_exists(micname)) cycle
            ! picker
            params_glob%lp = max(params%fny, params%lp_pick)
            call piter%iterate(cline, params%smpd, micname, output_dir_picker, boxfile, nptcls, moldiam_opt=moldiam_opt)
            call o_mic%set('nptcls', nptcls)
            if( nptcls > 0 )then
                call o_mic%set('boxfile', trim(boxfile))
                if( l_multipick ) call o_mic%set('moldiam', moldiam_opt)
            else
                call o_mic%set_state(0)
            endif
            ! update project
            call spproj%os_mic%set_ori(imic, o_mic)
            nvalid = nvalid+1
        enddo
        ! extract particles
        if( l_extract )then
            call spproj%write_segment_inside(params%oritype, params%projfile)
            call xextract%execute(cline_extract)
            ! nothing to write, done by extract
        else
            if( ntot > 1 )then
                ! purging state=0 and nptcls=0 mics such that all mics (nmics>1)
                ! can be assumed valid
                call os_mic%new(nvalid, is_ptcl=.false.)
                i = 0
                do imic = fromto(1),fromto(2)
                    state  = spproj%os_mic%get_state(imic)
                    nptcls = spproj%os_mic%get_int(imic,'nptcls')
                    if( (state == 1) .and. (nptcls > 0) )then
                        i = i+1
                        call os_mic%transfer_ori(i, spproj%os_mic, imic)
                    endif
                enddo
                spproj%os_mic = os_mic
                call os_mic%kill
            endif
            call spproj%write_segment_inside(params%oritype, params%projfile)
        endif
        ! end gracefully
        call qsys_job_finished(  'simple_commander_preprocess :: exec_pick_extract' )
        call o_mic%kill
        call piter%kill
        call simple_end('**** SIMPLE_PICK_EXTRACT NORMAL STOP ****')
    end subroutine exec_pick_extract

    subroutine exec_make_pickrefs( self, cline )
        use simple_masker,   only: automask2D
        use simple_default_clines
        class(make_pickrefs_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters)         :: params
        type(stack_io)           :: stkio_r
        type(oris)               :: moldiamori
        type(image)              :: ref2D, ref2D_clip
        type(image), allocatable :: projs(:), masks(:)
        real,        allocatable :: diams(:), shifts(:,:)
        real,    parameter :: MSKDIAM2LP = 0.15, lP_LB = 30., LP_UB = 15.
        integer, parameter :: NREFS=100
        real    :: ang, rot, lp, diam_max, smpd_here, sc
        integer :: nrots, iref, irot, ldim_clip(3), ldim(3), ldim_sc(3), ncavgs, icavg
        integer :: cnt, norefs, b, new_box
        logical :: l_scale, l_moldiam = .false.
        ! error check
        if( cline%defined('vol1') ) THROW_HARD('vol1 input no longer supported, use prg=reproject to generate 20 2D references')
        if( .not.cline%defined('pickrefs') ) THROW_HARD('PICKREFS must be informed!')
        if( cline%defined('moldiam') ) l_moldiam = .true.
        ! set defaults
        call set_automask2D_defaults(cline)
        call cline%set('oritype', 'mic')
        if( .not.cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( .not.cline%defined('neg')   ) call cline%set('neg',    'no') ! default for new picker
        ! parse parameters
        call params%new(cline)
        if( params%stream.eq.'yes' ) THROW_HARD('not a streaming application')
        ! read selected cavgs
        call find_ldim_nptcls(params%pickrefs, ldim, ncavgs, smpd=smpd_here)
        if( smpd_here < 0.01 ) THROW_HARD('Invalid sampling distance for the cavgs (should be in MRC format)')
        ldim(3) = 1
        ! Pixel size & scaling
        ! smpd_here is the pixel size read in
        ! params%smpd is the target pixel size of the micrograph
        if( .not.cline%defined('smpd') ) params%smpd = smpd_here
        l_scale    = .false.
        params%msk = real(ldim(1)/2) - COSMSKHALFWIDTH ! for automasking
        if( abs(smpd_here-params%smpd) > 0.001 )then
            sc = smpd_here/params%smpd
            b  = round2even(real(ldim(1))*sc)
            if( b /= ldim(1) )then
                ! adjust image size and mask
                l_scale    = .true.
                ldim_sc    = [b,b,1]
                params%msk = real(b/2) - COSMSKHALFWIDTH
            endif
        else
            ldim_sc = ldim
        endif
        ! read
        allocate( projs(ncavgs), masks(ncavgs) )
        call stkio_r%open(params%pickrefs, smpd_here, 'read', bufsz=ncavgs)
        do icavg=1,ncavgs
            call projs(icavg)%new(ldim, smpd_here)
            call stkio_r%read(icavg, projs(icavg))
            if( l_scale )then
                ! scaling to final dimension
                call projs(icavg)%fft
                if( ldim_sc(1) > ldim(1) )then
                    call projs(icavg)%pad_inplace(ldim_sc, antialiasing=.false.)
                else
                    call projs(icavg)%clip_inplace(ldim_sc)
                endif
                call projs(icavg)%ifft
            endif
            call masks(icavg)%copy(projs(icavg))
        end do
        call stkio_r%close
        ! Automasking
        call automask2D(masks, params%ngrow, nint(params%winsz), params%edge, diams, shifts)
        do icavg=1,ncavgs
            call projs(icavg)%div_below(0.,10.)
            call projs(icavg)%mul(masks(icavg))
            call projs(icavg)%shift([shifts(icavg,1),shifts(icavg,2),0.])
        end do
        ! estimate new box size and clip
        if( l_moldiam ) then
            diam_max = params%moldiam
        else
            diam_max = maxval(diams)
        end if
        lp      = min(max(LP_LB,MSKDIAM2LP * diam_max),LP_UB)
        new_box = round2even(diam_max / params%smpd + 2. * COSMSKHALFWIDTH)
        new_box = min(new_box, ldim_sc(1)) ! fail safe: new dimensions cannot be larger than required
        write(logfhandle,'(A,1X,I4)') 'ESTIMATED BOX SIZE: ', new_box
        ldim_clip = [new_box, new_box, 1]
        do icavg=1,ncavgs
            call projs(icavg)%bp(0.,lp)
        end do
        ! write diam_max to file
        call moldiamori%new(1, .false.)
        call moldiamori%set(1, "moldiam", real(round2even(diam_max)))
        call moldiamori%write(1, trim(STREAM_MOLDIAM))
        call moldiamori%kill
        ! expand in in-plane rotation, clip and write to file
        nrots  = nint( real(NREFS)/real(ncavgs) )
        norefs = ncavgs
        call ref2D_clip%new([ldim_clip(1),ldim_clip(2),1], params%smpd)
        if( nrots > 1 )then
            call ref2D%new([ldim(1),ldim(2),1], params%smpd)
            ang = 360./real(nrots)
            cnt = 0
            do iref=1,norefs
                rot = 0.
                do irot=1,nrots
                    cnt = cnt + 1
                    call projs(iref)%rtsq(rot, 0., 0., ref2D)
                    call ref2D%clip(ref2D_clip)
                    if( params%neg .eq. 'yes' ) call ref2D_clip%neg
                    call ref2D_clip%write(trim(PICKREFS_FBODY)//params%ext, cnt)
                    rot = rot + ang
                end do
            end do
        else
            ! should never happen
            do iref=1,norefs
                call projs(iref)%clip(ref2D_clip)
                if( params%neg .eq. 'yes' ) call ref2D_clip%neg
                call ref2D_clip%write(trim(PICKREFS_FBODY)//params%ext, iref)
            end do
        endif
        ! cleanup
        do icavg = 1,ncavgs
            call masks(icavg)%kill
            call projs(icavg)%kill
        enddo
        deallocate(masks,projs)
        call ref2D%kill
        call ref2D_clip%kill
        ! end gracefully
        call simple_touch('MAKE_PICKREFS_FINISHED', errmsg='In: commander_preprocess::exec_make_pickrefs')
        call simple_end('**** SIMPLE_MAKE_PICKREFS NORMAL STOP ****')
    end subroutine exec_make_pickrefs

    ! UTILITIES

    logical function box_inside( ildim, coord, box )
        integer, intent(in) :: ildim(3), coord(2), box
        integer             :: fromc(2), toc(2)
        fromc  = coord+1       ! compensate for the c-range that starts at 0
        toc    = fromc+(box-1) ! the lower left corner is 1,1
        box_inside = .true.    ! box is inside
        if( any(fromc < 1) .or. toc(1) > ildim(1) .or. toc(2) > ildim(2) ) box_inside = .false.
    end function box_inside

end module simple_commander_preprocess
