!@descr: for pre-processing (motion correction, CTF estimation etc.)
module simple_commanders_preprocess
use simple_commander_module_api
use simple_motion_correct_utils, only: flip_gain
use simple_mini_stream_utils,    only: segdiampick_preprocess
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_preprocess_distr
  contains
    procedure :: execute      => exec_preprocess_distr
end type commander_preprocess_distr

type, extends(commander_base) :: commander_preprocess
  contains
    procedure :: execute      => exec_preprocess
end type commander_preprocess

type, extends(commander_base) :: commander_motion_correct_distr
  contains
    procedure :: execute      => exec_motion_correct_distr
end type commander_motion_correct_distr

type, extends(commander_base) :: commander_motion_correct
  contains
    procedure :: execute      => exec_motion_correct
end type commander_motion_correct

type, extends(commander_base) :: commander_gen_pspecs_and_thumbs_distr
  contains
    procedure :: execute      => exec_gen_pspecs_and_thumbs_distr
end type commander_gen_pspecs_and_thumbs_distr

type, extends(commander_base) :: commander_gen_pspecs_and_thumbs
  contains
    procedure :: execute      => exec_gen_pspecs_and_thumbs
end type commander_gen_pspecs_and_thumbs

type, extends(commander_base) :: commander_ctf_estimate_distr
  contains
    procedure :: execute      => exec_ctf_estimate_distr
end type commander_ctf_estimate_distr

type, extends(commander_base) :: commander_ctf_estimate
  contains
    procedure :: execute      => exec_ctf_estimate
end type commander_ctf_estimate

contains

    subroutine exec_preprocess_distr( self, cline )
        class(commander_preprocess_distr), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(parameters) :: params
        type(qsys_env)   :: qenv
        type(chash)      :: job_descr
        type(sp_project) :: spproj
        if( .not. cline%defined('oritype')         ) call cline%set('oritype',        'mic')
        if( .not. cline%defined('stream')          ) call cline%set('stream',          'no')
        if( .not. cline%defined('mkdir')           ) call cline%set('mkdir',          'yes')
        ! motion correction
        if( .not. cline%defined('trs')             ) call cline%set('trs',              20.)
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart',           8.)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',            5.)
        if( .not. cline%defined('bfac')            ) call cline%set('bfac',             50.)
        if( .not. cline%defined('mcconvention')    ) call cline%set('mcconvention','simple')
        if( .not. cline%defined('eer_upsampling')  ) call cline%set('eer_upsampling',     1)
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
        ! extraction
        if( .not. cline%defined('pcontrast')       ) call cline%set('pcontrast',    'black')
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
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs(job_descr, algnfbody=string(ALGN_FBODY), array=L_USE_SLURM_ARR)
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
        use simple_motion_correct_iter, only: motion_correct_iter
        use simple_ctf_estimate_iter,   only: ctf_estimate_iter
        class(commander_preprocess), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters)              :: params
        type(ori)                     :: o_mov
        type(ctf_estimate_iter)       :: ctfiter
        type(motion_correct_iter)     :: mciter
        type(sp_project)              :: spproj
        type(ctfparams)               :: ctfvars
        type(Node), pointer           :: xmldoc, beamshiftnode, beamshiftnodex, beamshiftnodey
        type(string)  :: imgkind, moviename, fbody
        type(string)  :: moviename_forctf, output_dir_motion_correct
        type(string)  :: output_dir_ctf_estimate, output_dir_inipick_preproc, micname_intg
        type(string)  :: eputiltgroup, str_meta
        integer :: nmovies, fromto(2), imovie, ntot, frame_counter
        logical :: l_del_forctf
        call cline%set('oritype', 'mic')
        call params%new(cline)
        if( params%scale_movies > 1.01 )then
            THROW_HARD('scale_movies cannot be > 1; exec_preprocess')
        endif
        l_del_forctf = .false.
        ! read in movies
        call spproj%read( params%projfile )
        if( spproj%get_nmovies()==0 .and. spproj%get_nintgs()==0 ) THROW_HARD('No movie/micrograph to process!')
        ! output directories & naming
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
            select case(imgkind%to_char())
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
            ! read xml
            if(o_mov%isthere('meta'))then
                str_meta   = o_mov%get_str('meta') 
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
            ! assign beamtilt if EPU
            if(o_mov%isthere('intg')) then
                call o_mov%getter('intg', micname_intg)
                micname_intg = basename(micname_intg)
                if(micname_intg%to_char([1,8]) == 'FoilHole') then
                    ! EPU filename
                    eputiltgroup = micname_intg%to_char([micname_intg%substr_ind('Data_') + 5,micname_intg%strlen_trim()])
                    eputiltgroup = eputiltgroup%to_char([1,eputiltgroup%substr_ind('_') - 1])
                    call o_mov%set("tiltgrp", eputiltgroup%to_real())
                else
                    call o_mov%set("tiltgrp", 0.0)
                end if
            end if
            ! update project
            call spproj%os_mic%set_ori(imovie, o_mov)
        end do
        ! do initial picking preprocessing
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
        ! cleanup
        call o_mov%kill
        call spproj%kill
        ! end gracefully
        call qsys_job_finished(string('simple_commanders_preprocess :: exec_preprocess'))
        call simple_end('**** SIMPLE_PREPROCESS NORMAL STOP ****')
    end subroutine exec_preprocess

    subroutine exec_motion_correct_distr( self, cline )
        class(commander_motion_correct_distr), intent(inout) :: self
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
        call qenv%gen_scripts_and_schedule_jobs(job_descr, algnfbody=string(ALGN_FBODY), array=L_USE_SLURM_ARR)
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
        use simple_motion_correct_iter, only: motion_correct_iter
        class(commander_motion_correct), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline !< command line input
        type(parameters)              :: params
        type(motion_correct_iter)     :: mciter
        type(ctfparams)               :: ctfvars
        type(sp_project)              :: spproj
        type(ori)                     :: o
        type(string) :: output_dir, moviename, fbody
        integer :: nmovies, fromto(2), imovie, ntot, frame_counter, cnt
        call cline%set('oritype', 'mic')
        call params%new(cline)
        call spproj%read(params%projfile)
        ! sanity check
        nmovies = spproj%get_nmovies()
        if( nmovies == 0 )then
            THROW_HARD('No movie to process!')
        endif
        if( params%scale_movies > 1.01 )then
            THROW_HARD('scale_movies cannot be > 1; exec_motion_correct')
        endif
        if( cline%defined('gainref') )then
            if(.not.file_exists(params%gainref) )then
                THROW_HARD('gain reference: '//params%gainref%to_char()//' not found; motion_correct')
            endif
        endif
        ! output directory & names
        output_dir = PATH_HERE
        if( cline%defined('fbody') )then
            fbody = params%fbody
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
                        call mciter%iterate(cline, ctfvars, o, fbody, frame_counter, moviename, output_dir, gainref_fname=params%gainref)
                    else
                        call mciter%iterate(cline, ctfvars, o, fbody, frame_counter, moviename, output_dir)
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
        call qsys_job_finished(string('simple_commanders_preprocess :: exec_motion_correct'))
        call simple_end('**** SIMPLE_MOTION_CORRECT NORMAL STOP ****')
    end subroutine exec_motion_correct

    subroutine exec_gen_pspecs_and_thumbs_distr( self, cline )
        class(commander_gen_pspecs_and_thumbs_distr), intent(inout) :: self
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
        call qenv%gen_scripts_and_schedule_jobs(job_descr, algnfbody=string(ALGN_FBODY), array=L_USE_SLURM_ARR)
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
        use simple_pspec_thumb_iter, only: pspec_thumb_iter
        class(commander_gen_pspecs_and_thumbs), intent(inout) :: self
        class(cmdline),                         intent(inout) :: cline !< command line input
        type(parameters)              :: params
        type(pspec_thumb_iter)        :: ptiter
        type(sp_project)              :: spproj
        type(ori)                     :: o
        type(string) :: output_dir, moviename_intg, imgkind
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
                call ptiter%iterate(o, moviename_intg, output_dir)
                call spproj%os_mic%set_ori(iintg, o)
                write(logfhandle,'(f4.0,1x,a)') 100.*(real(cnt)/real(ntot)), 'percent of the integrated movies processed'
            endif
        end do
        ! output
        call binwrite_oritab(params%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        call o%kill
        ! end gracefully
        call qsys_job_finished(string('simple_commanders_preprocess :: exec_gen_pspecs_and_thumbs'))
        call simple_end('**** SIMPLE_GEN_PSPECS_AND_THUMBS NORMAL STOP ****')
    end subroutine exec_gen_pspecs_and_thumbs

    subroutine exec_ctf_estimate_distr( self, cline )
        class(commander_ctf_estimate_distr), intent(inout) :: self
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
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=string(ALGN_FBODY), array=L_USE_SLURM_ARR)
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
        use simple_ctf_estimate_iter,   only: ctf_estimate_iter
        class(commander_ctf_estimate), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline  !< command line input
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(ctf_estimate_iter)       :: ctfiter
        type(ctfparams)               :: ctfvars
        type(ori)                     :: o
        type(string) :: intg_forctf, output_dir, imgkind
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
                call ctfiter%iterate( ctfvars, intg_forctf, o, output_dir, l_gen_thumb)
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
        call qsys_job_finished(string('simple_commanders_preprocess :: exec_ctf_estimate'))
        call simple_end('**** SIMPLE_CTF_ESTIMATE NORMAL STOP ****')
    end subroutine exec_ctf_estimate

end module simple_commanders_preprocess
