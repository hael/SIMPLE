! concrete commander: pre-processing routines
module simple_commander_preprocess
include 'simple_lib.f08'
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_parameters,     only: parameters, params_glob
use simple_commander_base, only: commander_base
use simple_image,          only: image
use simple_ori,            only: ori
use simple_oris,           only: oris
use simple_sp_project,     only: sp_project
use simple_qsys_env,       only: qsys_env
use simple_qsys_funs
implicit none

public :: preprocess_commander_stream
public :: preprocess_commander_distr
public :: preprocess_commander
public :: motion_correct_commander_distr
public :: motion_correct_tomo_commander_distr
public :: motion_correct_commander
public :: gen_pspecs_and_thumbs_commander_distr
public :: gen_pspecs_and_thumbs_commander
public :: ctf_estimate_commander_distr
public :: ctf_estimate_commander
public :: map_cavgs_selection_commander
public :: pick_commander_distr
public :: pick_commander
public :: extract_commander_distr
public :: extract_commander
public :: reextract_commander_distr
public :: reextract_commander
public :: pick_extract_commander_stream
public :: pick_extract_commander
public :: make_pickrefs_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: preprocess_commander_stream
  contains
    procedure :: execute      => exec_preprocess_stream
end type preprocess_commander_stream
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
type, extends(commander_base) :: motion_correct_tomo_commander_distr
  contains
    procedure :: execute      => exec_motion_correct_tomo_distr
end type motion_correct_tomo_commander_distr
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
type, extends(commander_base) :: pick_extract_commander_stream
  contains
    procedure :: execute      => exec_pick_extract_stream
end type pick_extract_commander_stream
type, extends(commander_base) :: pick_extract_commander
  contains
    procedure :: execute      => exec_pick_extract
end type pick_extract_commander
type, extends(commander_base) :: make_pickrefs_commander
  contains
    procedure :: execute      => exec_make_pickrefs
end type make_pickrefs_commander

contains

    subroutine exec_preprocess_stream( self, cline )
        use simple_moviewatcher, only: moviewatcher
        class(preprocess_commander_stream), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(parameters)                       :: params
        integer,                   parameter   :: SHORTTIME = 60   ! folder watched every minute
        integer,                   parameter   :: LONGTIME  = 600  ! time lag after which a movie is processed
        class(cmdline),            allocatable :: completed_jobs_clines(:)
        type(qsys_env)                         :: qenv
        type(cmdline)                          :: cline_make_pickrefs
        type(moviewatcher)                     :: movie_buff
        type(sp_project)                       :: spproj, stream_spproj
        character(len=LONGSTRLEN), allocatable :: movies(:)
        character(len=:),          allocatable :: output_dir, output_dir_ctf_estimate, output_dir_picker
        character(len=:),          allocatable :: output_dir_motion_correct, output_dir_extract, stream_spprojfile
        character(len=LONGSTRLEN)              :: movie
        integer                                :: nmovies, imovie, stacksz, prev_stacksz, iter, icline
        integer                                :: nptcls, nptcls_prev, nmovs, nmovs_prev
        logical                                :: l_pick
        if( .not. cline%defined('oritype')         ) call cline%set('oritype',        'mic')
        if( .not. cline%defined('mkdir')           ) call cline%set('mkdir',          'yes')
        ! mnotion correction
        if( .not. cline%defined('trs')             ) call cline%set('trs',              10.)
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart',           8.)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',            5.)
        if( .not. cline%defined('bfac')            ) call cline%set('bfac',             50.)
        if( .not. cline%defined('groupframes')     ) call cline%set('groupframes',     'no')
        ! ctf estimation
        if( .not. cline%defined('pspecsz')         ) call cline%set('pspecsz',         512.)
        if( .not. cline%defined('hp_ctf_estimate') ) call cline%set('hp_ctf_estimate',  30.)
        if( .not. cline%defined('lp_ctf_estimate') ) call cline%set('lp_ctf_estimate',   5.)
        if( .not. cline%defined('dfmin')           ) call cline%set('dfmin',            0.3)
        if( .not. cline%defined('dfmax')           ) call cline%set('dfmax',            5.0)
        ! picking
        if( .not. cline%defined('lp_pick')         ) call cline%set('lp_pick',          20.)
        ! extraction
        if( .not. cline%defined('pcontrast')       ) call cline%set('pcontrast',    'black')
        if( cline%defined('refs') .and. cline%defined('vol1') )then
            THROW_HARD('REFS and VOL1 cannot be both provided!')
        endif
        call cline%set('numlen', 5.)
        call cline%set('stream','yes')
        ! master parameters
        call params%new(cline)
        params_glob%split_mode = 'stream'
        params_glob%ncunits    = params%nparts
        call cline%set('mkdir', 'no')
        call cline%set('prg',   'preprocess')
        if( cline%defined('dir_prev') .and. .not.file_exists(params%dir_prev) )then
            THROW_HARD('Directory '//trim(params%dir_prev)//' does not exist!')
        endif
        ! read in movies
        call spproj%read( params%projfile )
        ! sanity check
        if( spproj%os_mic%get_noris() /= 0 )then
            THROW_HARD('PREPROCESS_STREAM must always start from an empty project (eg from root proejct folder)')
        endif
        ! picking
        l_pick = .false.
        if( cline%defined('refs') .or. cline%defined('vol1') ) l_pick = .true.
        ! output directories
        output_dir = PATH_HERE
        output_dir_ctf_estimate   = filepath(trim(output_dir), trim(DIR_CTF_ESTIMATE))
        output_dir_motion_correct = filepath(trim(output_dir), trim(DIR_MOTION_CORRECT))
        call simple_mkdir(output_dir_ctf_estimate,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
        call simple_mkdir(output_dir_motion_correct,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
        if( l_pick )then
            output_dir_picker  = filepath(trim(output_dir), trim(DIR_PICKER))
            output_dir_extract = filepath(trim(output_dir), trim(DIR_EXTRACT))
            call simple_mkdir(output_dir_picker,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
            call simple_mkdir(output_dir_extract,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
        endif
        ! setup the environment for distributed execution
        call qenv%new(1,stream=.true.)
        ! prepares picking references
        if( l_pick )then
            cline_make_pickrefs = cline
            call cline_make_pickrefs%set('prg','make_pickrefs')
            call cline_make_pickrefs%set('stream','no')
            call qenv%exec_simple_prg_in_queue(cline_make_pickrefs, 'MAKE_PICKREFS_FINISHED')
            call cline%set('refs', trim(PICKREFS)//params%ext)
            call cline%delete('vol1')
            write(logfhandle,'(A)')'>>> PREPARED PICKING TEMPLATES'
        endif
        ! movie watcher init
        movie_buff = moviewatcher(LONGTIME)
        ! import previous runs
        call import_prev_streams
        ! start watching
        prev_stacksz = 0
        nmovies      = 0
        iter         = 0
        do
            if( file_exists(trim(TERM_STREAM)) )then
                write(logfhandle,'(A)')'>>> TERMINATING PREPROCESS STREAM'
                exit
            endif
            do while( file_exists(trim(PAUSE_STREAM)) )
                if( file_exists(trim(TERM_STREAM)) ) exit
                call write_singlelineoftext(PAUSE_STREAM, 'PAUSED')
                write(logfhandle,'(A,A)')'>>> PREPROCES STREAM PAUSED ',cast_time_char(simple_gettime())
                call simple_sleep(SHORTTIME)
            enddo
            iter = iter + 1
            call movie_buff%watch( nmovies, movies )
            ! append movies to processing stack
            if( nmovies > 0 )then
                do imovie = 1, nmovies
                    movie = trim(adjustl(movies(imovie)))
                    if( .not.file_exists(movie) )cycle ! petty triple checking
                    call create_individual_project
                    call qenv%qscripts%add_to_streaming( cline )
                enddo
            endif
            ! stream scheduling
            call qenv%qscripts%schedule_streaming( qenv%qdescr )
            stacksz = qenv%qscripts%get_stacksz()
            if( stacksz .ne. prev_stacksz )then
                prev_stacksz = stacksz
                write(logfhandle,'(A,I5)')'>>> MOVIES/MICROGRAPHS TO PROCESS: ', stacksz
            endif
            ! completed jobs update the current project
            if( qenv%qscripts%get_done_stacksz() > 0 )then
                ! append new processed movies to project
                call qenv%qscripts%get_stream_done_stack( completed_jobs_clines )
                nptcls_prev = spproj%get_nptcls()
                nmovs_prev  = spproj%os_mic%get_noris()
                do icline=1,size(completed_jobs_clines)
                    stream_spprojfile = completed_jobs_clines(icline)%get_carg('projfile')
                    call stream_spproj%read( stream_spprojfile )
                    call spproj%append_project(stream_spproj, 'mic')
                    if( l_pick )then
                        call spproj%append_project(stream_spproj, 'stk')
                    endif
                    call stream_spproj%kill()
                    deallocate(stream_spprojfile)
                enddo
                nptcls = spproj%get_nptcls()
                nmovs  = spproj%os_mic%get_noris()
                write(logfhandle,'(A,I5)')'>>> # MOVIES PROCESSED:    ',nmovs
                if( l_pick )then
                    write(logfhandle,'(A,I8)')'>>> # PARTICLES EXTRACTED: ',nptcls
                endif
                ! update for 2d streaming
                call update_projects_list
                deallocate(completed_jobs_clines)
                ! exit for trial runs
                if( cline%defined('nmovies_trial') )then
                    if( nmovs  >= params%nmovies_trial )then
                        write(logfhandle,'(A)')'>>> TRIAL # OF MOVIES REACHED'
                        exit
                    endif
                endif
                if( cline%defined('nptcls_trial') )then
                    if( l_pick .and. (nptcls >= params%nptcls_trial) )then
                        write(logfhandle,'(A)')'>>> TRIAL # OF PARTICLES REACHED'
                        exit
                    endif
                endif
                ! write
                call spproj%write
            else
                ! wait
                call simple_sleep(SHORTTIME)
            endif
        end do
        ! termination
        call spproj%write
        call spproj%kill
        ! cleanup
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_PREPROCESS_STREAM NORMAL STOP ****')
        contains

            subroutine update_projects_list
                type(cmdline) :: cline_mov
                character(len=:),          allocatable :: fname, abs_fname
                character(len=LONGSTRLEN), allocatable :: old_fnames(:), fnames(:)
                integer :: i, n_spprojs, n_old
                n_spprojs = size(completed_jobs_clines)
                if( n_spprojs == 0 )return
                if( file_exists(STREAM_SPPROJFILES) )then
                    ! append
                    call read_filetable(STREAM_SPPROJFILES, old_fnames)
                    n_old = size(old_fnames)
                    allocate(fnames(n_spprojs+n_old))
                    fnames(1:n_old) = old_fnames(:)
                    do i=1,n_spprojs
                        cline_mov = completed_jobs_clines(i)
                        fname     = trim(cline_mov%get_carg('projfile'))
                        abs_fname = simple_abspath(fname, errmsg='preprocess_stream :: update_projects_list 1')
                        fnames(n_old+i) = trim(abs_fname)
                        deallocate(abs_fname)
                    enddo
                else
                    ! first write
                    allocate(fnames(n_spprojs))
                    do i=1,n_spprojs
                        cline_mov = completed_jobs_clines(i)
                        fname     = trim(cline_mov%get_carg('projfile'))
                        abs_fname = simple_abspath(fname, errmsg='preprocess_stream :: update_projects_list 2')
                        fnames(i) = trim(abs_fname)
                        deallocate(abs_fname)
                    enddo
                endif
                call write_filetable(STREAM_SPPROJFILES, fnames)
            end subroutine update_projects_list

            subroutine create_individual_project
                type(sp_project)              :: spproj_here
                type(cmdline)                 :: cline_here
                type(ctfparams)               :: ctfvars
                character(len=STDLEN)         :: ext, movie_here
                character(len=LONGSTRLEN)     :: projname, projfile
                movie_here = basename(trim(movie))
                ext        = fname2ext(trim(movie_here))
                projname   = 'preprocess_'//trim(get_fbody(trim(movie_here), trim(ext)))
                projfile   = trim(projname)//trim(METADATA_EXT)
                call cline_here%set('projname', trim(projname))
                call cline_here%set('projfile', trim(projfile))
                call spproj_here%update_projinfo(cline_here)
                spproj_here%compenv  = spproj%compenv
                spproj_here%jobproc  = spproj%jobproc
                ctfvars%ctfflag      = CTFFLAG_YES
                ctfvars%smpd         = params%smpd
                ctfvars%cs           = params%cs
                ctfvars%kv           = params%kv
                ctfvars%fraca        = params%fraca
                ctfvars%l_phaseplate = params%phaseplate.eq.'yes'
                call spproj_here%add_single_movie(trim(movie), ctfvars)
                call spproj_here%write
                call spproj_here%kill
                call cline%set('projname', trim(projname))
                call cline%set('projfile', trim(projfile))
            end subroutine create_individual_project

            !>  import previous run to the current project based on past single project files
            subroutine import_prev_streams
                use simple_ori, only: ori
                type(ori) :: o, o_stk
                character(len=LONGSTRLEN), allocatable :: sp_files(:)
                character(len=:), allocatable :: mic, mov
                integer :: iproj,nprojs,cnt
                logical :: err
                if( .not.cline%defined('dir_prev') ) return
                call simple_list_files(trim(params%dir_prev)//'/preprocess_*.simple', sp_files)
                nprojs = size(sp_files)
                cnt    = 0
                do iproj = 1,nprojs
                    err = .false.
                    call stream_spproj%read( sp_files(iproj) )
                    if( stream_spproj%os_mic%get_noris() /= 1 )then
                        THROW_WARN('Ignoring previous project'//trim(sp_files(iproj)))
                        cycle
                    endif
                    call stream_spproj%os_mic%get_ori(1, o)
                    ! import mic segment
                    call movefile2folder('intg', output_dir_motion_correct, o, err)
                    if( err ) cycle
                    call movefile2folder('forctf', output_dir_motion_correct, o, err)
                    call movefile2folder('thumb', output_dir_motion_correct, o, err)
                    call movefile2folder('mc_starfile', output_dir_motion_correct, o, err)
                    call movefile2folder('mceps', output_dir_motion_correct, o, err)
                    call movefile2folder('ctfjpg', output_dir_ctf_estimate, o, err)
                    call movefile2folder('ctfdoc', output_dir_ctf_estimate, o, err)
                    if( .not.l_pick )then
                        ! import mic segment
                        call stream_spproj%os_mic%set_ori(1, o)
                        call spproj%append_project(stream_spproj, 'mic')
                    else
                        ! import mic & stk segment
                        call movefile2folder('boxfile', output_dir_picker, o, err)
                        call stream_spproj%os_mic%set_ori(1, o)
                        call spproj%append_project(stream_spproj, 'mic')
                        if( .not.err )then
                            if( stream_spproj%os_stk%get_noris() == 1 )then
                                call stream_spproj%os_stk%get_ori(1, o_stk)
                                call movefile2folder('stk', output_dir_extract, o_stk, err)
                                if( .not.err )then
                                    call stream_spproj%os_stk%set_ori(1, o_stk)
                                    call spproj%append_project(stream_spproj, 'stk')
                                endif
                            endif
                        endif
                    endif
                    ! add to history
                    call o%getter('movie', mov)
                    call o%getter('intg', mic)
                    call movie_buff%add2history(mov)
                    call movie_buff%add2history(mic)
                    ! write updated individual project file
                    call stream_spproj%write(basename(sp_files(iproj)))
                    ! count
                    cnt = cnt + 1
                    ! cleanup
                    call stream_spproj%kill
                enddo
                call o%kill
                call o_stk%kill
                write(*,'(A,I3)')'>>> IMPORTED PREVIOUS PROCESSED MOVIES: ', cnt
            end subroutine import_prev_streams

            subroutine movefile2folder(key, folder, o, err)
                use simple_ori, only: ori
                character(len=*), intent(in)    :: key, folder
                class(ori),       intent(inout) :: o
                logical,          intent(out)   :: err
                character(len=:), allocatable :: src
                character(len=LONGSTRLEN) :: dest,reldest
                integer :: iostat
                err = .false.
                if( .not.o%isthere(key) )then
                    err = .true.
                    return
                endif
                call o%getter(key,src)
                if( .not.file_exists(src) )then
                    err = .true.
                    return
                endif
                dest   = trim(folder)//'/'//basename(src)
                iostat = rename(src,dest)
                if( iostat /= 0 )then
                    THROW_WARN('Ignoring '//trim(src))
                    return
                endif
                call make_relativepath(CWD_GLOB,dest,reldest)
                call o%set(key,reldest)
            end subroutine movefile2folder

    end subroutine exec_preprocess_stream

    subroutine exec_preprocess_distr( self, cline )
        class(preprocess_commander_distr), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(parameters)              :: params
        type(qsys_env)                :: qenv
        type(cmdline)                 :: cline_make_pickrefs
        type(chash)                   :: job_descr
        type(sp_project)              :: spproj
        logical                       :: l_pick
        if( .not. cline%defined('oritype')         ) call cline%set('oritype',        'mic')
        if( .not. cline%defined('stream')          ) call cline%set('stream',          'no')
        if( .not. cline%defined('mkdir')           ) call cline%set('mkdir',          'yes')
        ! mnotion correction
        if( .not. cline%defined('trs')             ) call cline%set('trs',              10.)
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart',           8.)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',            5.)
        if( .not. cline%defined('bfac')            ) call cline%set('bfac',             50.)
        if( .not. cline%defined('groupframes')     ) call cline%set('groupframes',     'no')
        ! ctf estimation
        if( .not. cline%defined('pspecsz')         ) call cline%set('pspecsz',         512.)
        if( .not. cline%defined('hp_ctf_estimate') ) call cline%set('hp_ctf_estimate',  30.)
        if( .not. cline%defined('lp_ctf_estimate') ) call cline%set('lp_ctf_estimate',   5.)
        if( .not. cline%defined('dfmin')           ) call cline%set('dfmin',            0.3)
        if( .not. cline%defined('dfmax')           ) call cline%set('dfmax',            5.0)
        ! picking
        if( .not. cline%defined('lp_pick')         ) call cline%set('lp_pick',          20.)
        ! extraction
        if( .not. cline%defined('pcontrast')       ) call cline%set('pcontrast',    'black')
        if( cline%defined('refs') .and. cline%defined('vol1') )then
            THROW_HARD('REFS and VOL1 cannot be both provided!')
        endif
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
        call cline%set('numlen', real(params%numlen))
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepares picking references
        l_pick = .false.
        if( cline%defined('refs') .or. cline%defined('vol1') )then
            l_pick = .true.
            cline_make_pickrefs = cline
            call cline_make_pickrefs%set('prg','make_pickrefs')
            call qenv%exec_simple_prg_in_queue(cline_make_pickrefs, 'MAKE_PICKREFS_FINISHED')
            call cline%set('refs', trim(PICKREFS)//params%ext)
            call cline%delete('vol1')
            write(logfhandle,'(A)')'>>> PREPARED PICKING TEMPLATES'
        endif
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs(job_descr, algnfbody=trim(ALGN_FBODY))
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
        use simple_sp_project,          only: sp_project
        use simple_motion_correct_iter, only: motion_correct_iter
        use simple_ctf_estimate_iter,   only: ctf_estimate_iter
        use simple_picker_iter,         only: picker_iter
        use simple_binoris_io,          only: binwrite_oritab
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
        character(len=:), allocatable :: imgkind, moviename, output_dir_picker, fbody
        character(len=:), allocatable :: moviename_forctf, moviename_intg, output_dir_motion_correct
        character(len=:), allocatable :: output_dir_ctf_estimate, output_dir_extract
        character(len=LONGSTRLEN)     :: boxfile
        integer :: nmovies, fromto(2), imovie, ntot, frame_counter, nptcls_out
        logical :: l_pick
        call cline%set('oritype', 'mic')
        call params%new(cline)
        if( params%scale > 1.05 )then
            THROW_HARD('scale cannot be > 1; exec_preprocess')
        endif
        if( params%tomo .eq. 'yes' )then
            THROW_HARD('tomography mode (tomo=yes) not yet supported!')
        endif
        l_pick = cline%defined('refs')
        ! read in movies
        call spproj%read( params%projfile )
        if( spproj%get_nmovies()==0 .and. spproj%get_nintgs()==0 ) THROW_HARD('No movie/micrograph to process!')
        ! output directories & naming
        output_dir_ctf_estimate   = PATH_HERE
        output_dir_motion_correct = PATH_HERE
        if( l_pick )then
            output_dir_picker  = PATH_HERE
        endif
        if( params%stream.eq.'yes' )then
            output_dir_ctf_estimate   = trim(DIR_CTF_ESTIMATE)
            output_dir_motion_correct = trim(DIR_MOTION_CORRECT)
            call simple_mkdir(output_dir_ctf_estimate,errmsg="commander_preprocess :: preprocess; ")
            call simple_mkdir(output_dir_motion_correct, errmsg="commander_preprocess :: preprocess;")
            if( l_pick )then
                output_dir_picker  = trim(DIR_PICKER)
                output_dir_extract = trim(DIR_EXTRACT)
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
        if( params%stream.eq.'yes' )then
            ! STREAMING MODE
            fromto(:) = 1
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
            call o_mov%delete_entry('forctf')
            call del_file(moviename_forctf)
            ! update project
            call spproj%os_mic%set_ori(imovie, o_mov)
            ! picker
            if( l_pick )then
                params_glob%lp = max(params%fny, params%lp_pick)
                moviename_intg = mciter%get_moviename('intg')
                call piter%iterate(cline,  moviename_intg, boxfile, nptcls_out, output_dir_picker)
                call o_mov%set('boxfile', trim(boxfile)   )
                call o_mov%set('nptcls',  real(nptcls_out))
                ! update project
                call spproj%os_mic%set_ori(imovie, o_mov)
                ! extract particles
                if( trim(params%stream) .eq. 'yes' )then
                    ! needs to write and re-read project at the end as extract overwrites it
                    call spproj%write_segment_inside(params%oritype)
                    cline_extract = cline
                    call cline_extract%set('dir', trim(output_dir_extract))
                    call cline_extract%set('pcontrast', params%pcontrast)
                    call cline_extract%delete('msk')
                    if( cline%defined('box_extract') )call cline_extract%set('box', real(params%box_extract))
                    call xextract%execute(cline_extract)
                    call spproj%kill
                endif
            endif
        end do
        if( params%stream .eq. 'yes' )then
            if( .not.l_pick )then
                ! because extract performs the writing otherwise
                call spproj%write_segment_inside(params%oritype)
            endif
        else
            call binwrite_oritab(params%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        endif
        call o_mov%kill
        ! end gracefully
        call qsys_job_finished(  'simple_commander_preprocess :: exec_preprocess' )
        call simple_end('**** SIMPLE_PREPROCESS NORMAL STOP ****')
    end subroutine exec_preprocess

    subroutine exec_motion_correct_distr( self, cline )
        class(motion_correct_commander_distr), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        type(qsys_env)   :: qenv
        type(chash)      :: job_descr
        if( .not. cline%defined('mkdir')      ) call cline%set('mkdir',       'yes')
        if( .not. cline%defined('trs')        ) call cline%set('trs',           10.)
        if( .not. cline%defined('lpstart')    ) call cline%set('lpstart',        8.)
        if( .not. cline%defined('lpstop')     ) call cline%set('lpstop',         5.)
        if( .not. cline%defined('bfac')       ) call cline%set('bfac',          50.)
        if( .not. cline%defined('groupframes')) call cline%set('groupframes',  'no')
        if( .not. cline%defined('wcrit')      ) call cline%set('wcrit',   'softmax')
        call cline%set('oritype', 'mic')
        call params%new(cline)
        call cline%set('numlen', real(params%numlen))
        ! sanity check
        call spproj%read_segment(params%oritype, params%projfile)
        if( spproj%get_nmovies() ==0 ) THROW_HARD('no movies to process! exec_motion_correct_distr')
        call spproj%kill
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY))
        ! merge docs
        call spproj%read(params%projfile)
        call spproj%merge_algndocs(params%nptcls, params%nparts, 'mic', ALGN_FBODY)
        call spproj%kill
        ! clean
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_MOTION_CORRECT NORMAL STOP ****')
    end subroutine exec_motion_correct_distr

    subroutine exec_motion_correct_tomo_distr( self, cline )
        use simple_oris, only: oris
        class(motion_correct_tomo_commander_distr), intent(inout) :: self
        class(cmdline),                             intent(inout) :: cline
        character(len=LONGSTRLEN), allocatable :: tomonames(:)
        type(parameters)         :: params
        type(oris)               :: exp_doc
        integer                  :: nseries, ipart
        type(qsys_env)           :: qenv
        character(len=KEYLEN)    :: str
        type(chash)              :: job_descr
        type(chash), allocatable :: part_params(:)
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',     'yes')
        if( .not. cline%defined('trs')     ) call cline%set('trs',         10.)
        if( .not. cline%defined('lpstart') ) call cline%set('lpstart',     20.)
        if( .not. cline%defined('lpstop')  ) call cline%set('lpstop',       6.)
        if( .not. cline%defined('tomo')    ) call cline%set('tomo',      'yes')
        if( .not. cline%defined('oritype') ) call cline%set('oritype',   'stk')
        if( .not. cline%defined('wcrit')   ) call cline%set('wcrit', 'softmax')
        call cline%set('prg', 'motion_correct')
        call params%new(cline)
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        if( cline%defined('tomoseries') )then
            call read_filetable(params%tomoseries, tomonames)
        else
            THROW_HARD('need tomoseries (filetable of filetables) to be part of the command line when tomo=yes')
        endif
        nseries = size(tomonames)
        call exp_doc%new(nseries)
        if( cline%defined('exp_doc') )then
            if( file_exists(params%exp_doc) )then
                call exp_doc%read(params%exp_doc)
            else
                THROW_HARD('the required parameter file (flag exp_doc): '//trim(params%exp_doc)//' not in cwd')
            endif
        else
            THROW_HARD('need exp_doc (line: exp_time=X dose_rate=Y) to be part of the command line when tomo=yes')
        endif
        params%nparts = nseries
        params%nptcls = nseries
        ! prepare part-dependent parameters
        allocate(part_params(params%nparts), stat=alloc_stat) ! -1. is default excluded value
        if(alloc_stat.ne.0)call allocchk("exec_motion_correct_tomo_distr ", alloc_stat)
        do ipart=1,params%nparts
            call part_params(ipart)%new(4)
            call part_params(ipart)%set('filetab', trim(tomonames(ipart)))
            call part_params(ipart)%set('fbody', 'tomo'//int2str_pad(ipart,params%numlen_tomo))
            str = real2str(exp_doc%get(ipart,'exp_time'))
            call part_params(ipart)%set('exp_time', trim(str))
            str = real2str(exp_doc%get(ipart,'dose_rate'))
            call part_params(ipart)%set('dose_rate', trim(str))
        end do
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs( job_descr, part_params=part_params)
        call qsys_cleanup
        call simple_end('**** SIMPLE_DISTR_MOTION_CORRECT_TOMO NORMAL STOP ****')
    end subroutine exec_motion_correct_tomo_distr

    subroutine exec_motion_correct( self, cline )
        use simple_sp_project,          only: sp_project
        use simple_binoris_io,          only: binwrite_oritab
        use simple_motion_correct_iter, only: motion_correct_iter
        class(motion_correct_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline !< command line input
        type(parameters)              :: params
        type(motion_correct_iter)     :: mciter
        type(ctfparams)               :: ctfvars
        type(sp_project)              :: spproj
        type(ori)                     :: o
        character(len=:), allocatable :: output_dir, moviename, imgkind, fbody
        integer :: nmovies, fromto(2), imovie, ntot, frame_counter, lfoo(3), nframes, cnt
        call cline%set('oritype', 'mic')
        call params%new(cline)
        call spproj%read(params%projfile)
        ! sanity check
        nmovies = spproj%get_nmovies()
        if( nmovies == 0 )then
            THROW_HARD('No movie to process!')
        endif
        if( params%scale > 1.05 )then
            THROW_HARD('scale cannot be > 1; exec_motion_correct')
        endif
        if( params%tomo .eq. 'yes' )then
            if( .not. params%l_dose_weight )then
                write(logfhandle,*) 'tomo=yes only supported with dose weighting!'
                THROW_HARD('give total exposure time: exp_time (in seconds) and dose_rate (in e/A2/s)')
            endif
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
        if( params%tomo .eq. 'no' )then
            if( cline%defined('fromp') .and. cline%defined('top') )then
                fromto = [params%fromp, params%top]
            else
                THROW_HARD('fromp & top args need to be defined in parallel execution; motion_correct')
            endif
        else
            ! all movies
            fromto(1) = 1
            fromto(2) = spproj%os_mic%get_noris()
        endif
        ntot = fromto(2) - fromto(1) + 1
        ! for series of tomographic movies we need to calculate the time_per_frame
        if( params%tomo .eq. 'yes' )then
            ! get number of frames & dim from stack
            call spproj%os_mic%getter(1, 'movie', moviename)
            call find_ldim_nptcls(moviename, lfoo, nframes)
            ! calculate time_per_frame
            params%time_per_frame = params%exp_time/real(nframes*nmovies)
        endif
        ! align
        frame_counter = 0
        cnt = 0
        do imovie=fromto(1),fromto(2)
            call spproj%os_mic%get_ori(imovie, o)
            if( o%isthere('imgkind').and.o%isthere('movie') )then
                cnt = cnt + 1
                call o%getter('imgkind', imgkind)
                if( imgkind.ne.'movie' )cycle
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
        call cline%set('numlen', real(params%numlen))
        ! sanity check
        call spproj%read_segment(params%oritype, params%projfile)
        nintgs = spproj%get_nintgs()
        if( nintgs ==0 )then
            THROW_HARD('no integrated movies to process! exec_gen_pspecs_and_thumbs_distr')
        endif
        if( params%nparts > nintgs )then
            call cline%set('nparts', real(nintgs))
            params%nparts = nintgs
        endif
        call spproj%kill
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY))
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
        use simple_binoris_io,       only: binwrite_oritab
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
        if( .not. cline%defined('pspecsz') ) call cline%set('pspecsz', 512.)
        if( .not. cline%defined('hp')      ) call cline%set('hp',       30.)
        if( .not. cline%defined('lp')      ) call cline%set('lp',        5.)
        if( .not. cline%defined('dfmin')   ) call cline%set('dfmin',    0.3)
        if( .not. cline%defined('dfmax')   ) call cline%set('dfmax',    5.0)
        if( .not. cline%defined('oritype') ) call cline%set('oritype','mic')
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
        call cline%set('numlen', real(params%numlen))
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY))
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
        use simple_binoris_io,          only: binwrite_oritab
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
        logical                       :: l_gen_thumb
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
                fromto(1) = params%fromp
                fromto(2) = params%top
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
            if( o%isthere('state') ) state = nint(o%get('state'))
            if( state == 0 ) cycle
            if( o%isthere('imgkind') )then
                call o%getter('imgkind', imgkind)
                if( imgkind.ne.'mic' )cycle
                if( o%isthere('forctf') )then
                    call o%getter('forctf', intg_forctf)
                else if( o%isthere('intg') )then
                    call o%getter('intg', intg_forctf)
                else
                    THROW_HARD('no image available (forctf|intg) for CTF fittings :: exec_ctf_estimate')
                endif
                l_gen_thumb = .not. o%isthere('thumb')
                ctfvars     = o%get_ctfvars()
                call ctfiter%iterate( ctfvars, intg_forctf, o, trim(output_dir), l_gen_thumb)
                ! delete file after estimation
                call o%delete_entry('forctf')
                call del_file(intg_forctf)
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
        class(cmdline),                       intent(inout) :: cline !< command line input
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
        call build%init_params_and_build_spproj(cline,params)
        ! find number of selected cavgs
        call find_ldim_nptcls(params%stk2, lfoo, nsel)
        ! find number of original cavgs
        if( .not. cline%defined('stk' ) )then
            call build%spproj%get_cavgs_stk(cavgstk, nall, smpd)
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
        allocate(states(nall))
        states = 0
        do isel=1,nsel
            loc = maxloc(correlations(isel,:))
            states(loc(1)) = 1
        end do
        ! communicate selection to project
        call build%spproj%map_cavgs_selection(states)
        ! this needs to be a full write as many segments are updated
        call build%spproj%write
        ! end gracefully
        call simple_end('**** SIMPLE_MAP_SELECTION NORMAL STOP ****')
    end subroutine exec_map_cavgs_selection

    subroutine exec_pick_distr( self, cline )
        class(pick_commander_distr), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        type(cmdline)    :: cline_make_pickrefs
        type(qsys_env)   :: qenv
        type(chash)      :: job_descr
        logical :: use_segmentation
        if( cline%defined('refs') .and. cline%defined('vol1') )then
            THROW_HARD('REFS and VOL1 cannot be both provided!')
        endif
        use_segmentation = .true.
        if( cline%defined('refs') .or. cline%defined('vol1') ) use_segmentation = .false.
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir',       'yes')
        if( .not. cline%defined('pcontrast') ) call cline%set('pcontrast', 'black')
        if( .not. cline%defined('oritype')   ) call cline%set('oritype',     'mic')
        call params%new(cline)
        ! sanity check
        call spproj%read_segment(params%oritype, params%projfile)
        if( spproj%get_nintgs() ==0 )then
            THROW_HARD('No micrograph to process! exec_pick_distr')
        endif
        call spproj%kill
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        params%numlen = len(int2str(params%nparts))
        call cline%set('numlen', real(params%numlen))
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepares picking references
        if( .not.use_segmentation )then
            cline_make_pickrefs = cline
            call cline_make_pickrefs%set('prg','make_pickrefs')
            call qenv%exec_simple_prg_in_queue(cline_make_pickrefs, 'MAKE_PICKREFS_FINISHED')
            call cline%set('refs', trim(PICKREFS)//params%ext)
            call cline%delete('vol1')
            write(logfhandle,'(A)')'>>> PREPARED PICKING TEMPLATES'
        endif
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY))
        ! merge docs
        call spproj%read(params%projfile)
        call spproj%merge_algndocs(params%nptcls, params%nparts, 'mic', ALGN_FBODY)
        ! cleanup
        call qsys_cleanup
        ! graceful exit
        call simple_end('**** SIMPLE_DISTR_PICK NORMAL STOP ****')
    end subroutine exec_pick_distr

    subroutine exec_pick( self, cline )
        use simple_binoris_io,     only: binwrite_oritab
        use simple_picker_iter,    only: picker_iter
        class(pick_commander), intent(inout) :: self
        class(cmdline),        intent(inout) :: cline !< command line input
        type(parameters)              :: params
        type(sp_project)              :: spproj
        type(picker_iter)             :: piter
        type(ori)                     :: o
        character(len=:), allocatable :: output_dir, intg_name, imgkind
        character(len=LONGSTRLEN)     :: boxfile
        integer :: fromto(2), imic, ntot, nptcls_out, cnt, state
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
                fromto(1) = params%fromp
                fromto(2) = params%top
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
            if( o%isthere('state') ) state = nint(o%get('state'))
            if( state == 0 ) cycle
            if( o%isthere('imgkind') )then
                call o%getter('imgkind', imgkind)
                if( imgkind.ne.'mic' )cycle
                call o%getter('intg', intg_name)
                call piter%iterate(cline, intg_name, boxfile, nptcls_out, output_dir)
                call spproj%os_mic%set_boxfile(imic, boxfile, nptcls=nptcls_out)
            endif
            write(logfhandle,'(f4.0,1x,a)') 100.*(real(cnt)/real(ntot)), 'percent of the micrographs processed'
        end do
        ! output
        call binwrite_oritab(params%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        call o%kill
        ! end gracefully
        call qsys_job_finished(  'simple_commander_preprocess :: exec_pick' )
        call simple_end('**** SIMPLE_PICK NORMAL STOP ****')
    end subroutine exec_pick

    subroutine exec_extract_distr( self, cline )
        use simple_oris,  only: oris
        use simple_ori,   only: ori
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
        real    :: dfx,dfy
        integer :: boxcoords(2), lfoo(3)
        integer :: nframes,imic,i,nmics_tot,numlen,nmics,cnt,state,istk,nstks,ipart
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir',       'yes')
        if( .not. cline%defined('outside')   ) call cline%set('outside',      'no')
        if( .not. cline%defined('pcontrast') ) call cline%set('pcontrast', 'black')
        if( .not. cline%defined('stream')    ) call cline%set('stream',       'no')
        if( cline%defined('ctf') )then
            if( cline%get_carg('ctf').ne.'flip' .and. cline%get_carg('ctf').ne.'no' )then
                THROW_HARD('Only CTF=NO/FLIP are allowed')
            endif
        endif
        call cline%set('nthr', 1.)
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
        call spproj%write
        ! input directory
        if( cline%defined('dir_box') )then
            if( params%mkdir.eq.'yes' .and. params%dir_box(1:1).ne.'/')then
                params%dir_box = trim(filepath(PATH_PARENT,params%dir_box))
            endif
            params%dir_box = simple_abspath(params%dir_box)
            if( file_exists(params%dir_box) )then
                call simple_list_files(trim(params%dir_box)//'/*.box', boxfiles)
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
        ! sanity checks
        nmics  = 0
        do imic = 1, nmics_tot
            call spproj%os_mic%get_ori(imic, o_mic)
            state = 1
            if( o_mic%isthere('state') ) state = nint(o_mic%get('state'))
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
        ! DISTRIBUTED EXTRACTION
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        ! prepare job description
        call cline%gen_job_descr(job_descr)
        ! schedule & clean
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY))
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
                if( nint(o_mic%get('nptcls')) > 0 ) nstks = nstks + 1
            enddo
            call spproj_part%kill
        enddo
        if( cnt /= nmics_tot ) THROW_HARD('Inconstistent number of micrographs in individual projects')
        ! fetch stacks table
        if( nstks > 0 )then
            call os_stk%new(nstks)
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
                        dfx = spproj_part%os_ptcl2D%get(i, 'dfx')
                        dfy = spproj_part%os_ptcl2D%get(i, 'dfy')
                        call spproj%os_ptcl2D%set(cnt,'dfx',dfx)
                        call spproj%os_ptcl2D%set(cnt,'dfy',dfy)
                        call spproj%os_ptcl3D%set(cnt,'dfx',dfx)
                        call spproj%os_ptcl3D%set(cnt,'dfy',dfy)
                    endif
                enddo
                call spproj_part%kill
            enddo
            call os_stk%kill
        endif
        ! final write
        call spproj%write
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
        use simple_ctf,              only: ctf
        use simple_ctf_estimate_fit, only: ctf_estimate_fit
        class(extract_commander), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline !< command line input
        type(builder)                           :: build
        type(parameters)                        :: params
        type(sp_project)                        :: spproj_in, spproj
        type(nrtxtfile)                         :: boxfile
        type(image)                             :: micrograph
        type(ori)                               :: o_mic, o_tmp
        type(ctf)                               :: tfun
        type(ctfparams)                         :: ctfparms
        type(ctf_estimate_fit)                  :: ctffit
        character(len=:),           allocatable :: output_dir, mic_name, imgkind
        real,                       allocatable :: boxdata(:,:)
        logical,                    allocatable :: oris_mask(:), mics_mask(:)
        character(len=LONGSTRLEN) :: stack, boxfile_name, box_fname, ctfdoc
        integer                   :: nframes, imic, iptcl, ldim(3), nptcls,nmics,nmics_here,box, box_first, fromto(2)
        integer                   :: cnt, nmics_tot, lfoo(3), ifoo, noutside, state, iptcl_glob, cnt_stats
        real                      :: ptcl_pos(2), meanv,sddevv,minv,maxv,stk_stats(4),dfx,dfy,sdev_noise
        logical                   :: l_err, l_ctfpatch
        call cline%set('oritype', 'mic')
        call cline%set('mkdir',   'no')
        call params%new(cline)
        ! init
        output_dir = PATH_HERE
        fromto(:)  = [params%fromp, params%top]
        nmics_here = fromto(2)-fromto(1)+1
        if( params%stream.eq.'yes' )then
            output_dir = DIR_EXTRACT
            fromto(:)  = [1,1]
            nmics_here = 1
            ! read in integrated movies, output project = input project
            call spproj%read(params%projfile)
            nmics_tot = spproj_in%os_mic%get_noris()
            if( spproj%get_nintgs() /= 1 ) THROW_HARD('Incompatible # of integrated micrograph to process!')
        else
            ! read in integrated movies
            call spproj_in%read_segment(params%oritype, params%projfile)
            nmics_tot = spproj_in%os_mic%get_noris()
            if( spproj_in%get_nintgs() == 0 ) THROW_HARD('No integrated micrograph to process!')
            ! init output project
            call spproj%read_non_data_segments(params%projfile)
            call spproj%projinfo%set(1,'projname', get_fbody(params%outfile,METADATA_EXT,separator=.false.))
            call spproj%projinfo%set(1,'projfile', params%outfile)
            params%projfile = trim(params%outfile) ! for builder later
            call spproj%os_mic%new(nmics_here)
            cnt = 0
            do imic = fromto(1),fromto(2)
                cnt = cnt + 1
                call spproj_in%os_mic%get_ori(imic, o_tmp)
                call spproj%os_mic%set_ori(cnt, o_tmp)
            enddo
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
            if( o_mic%isthere('state') ) state = nint(o_mic%get('state'))
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
                call make_relativepath(CWD_GLOB,trim(box_fname),boxfile_name)
                call spproj%os_mic%set_boxfile(imic, boxfile_name)
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
                    call spproj%os_mic%set(imic, 'nptcls', 0.)
                    cycle
                endif
                allocate( boxdata(nptcls,boxfile%get_nrecs_per_line()), stat=alloc_stat)
                call boxfile%readNextDataLine(boxdata(1,:))
                call boxfile%kill
                params%box = nint(boxdata(1,3))
                params%boxmatch = params%box
            endif
        enddo
        call spproj%write
        call spproj%kill
        ! actual extraction
        if( nmics == 0 )then
            ! done
        else
            if( params%box == 0 )THROW_HARD('box cannot be zero; exec_extract')
            ! set normalization radius
            params%msk = RADFRAC_NORM_EXTRACT * real(params%box/2)
            ! init
            call build%build_general_tbox(params, cline, do3d=.false.)
            call micrograph%new([ldim(1),ldim(2),1], params%smpd)
            noutside  = 0
            box_first = 0
            ! main loop
            iptcl_glob = 0 ! extracted particle index among ALL stacks
            do imic = 1,nmics_here
                if( .not.mics_mask(imic) )cycle
                ! fetch micrograph
                call build%spproj_field%get_ori(imic, o_mic)
                call o_mic%getter('imgkind', imgkind)
                call o_mic%getter('intg', mic_name)
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
                allocate(oris_mask(nptcls), source=.false., stat=alloc_stat)
                if(alloc_stat.ne.0)call allocchk("In exec_extract oris_mask")
                ! read box data & update mask
                if(allocated(boxdata))deallocate(boxdata)
                allocate( boxdata(nptcls,boxfile%get_nrecs_per_line()), stat=alloc_stat)
                if(alloc_stat.ne.0)call allocchk('In: exec_extract; boxdata etc., 2',alloc_stat)
                do iptcl=1,nptcls
                    call boxfile%readNextDataLine(boxdata(iptcl,:))
                    box = nint(boxdata(iptcl,3))
                    if( nint(boxdata(iptcl,3)) /= nint(boxdata(iptcl,4)) )then
                        THROW_HARD('only square windows allowed; exec_extract')
                    endif
                    ! modify coordinates if change in box (shift by half the difference)
                    if( box /= params%box ) boxdata(iptcl,1:2) = boxdata(iptcl,1:2) - real(params%box-box)/2.
                    if( .not.cline%defined('box') .and. nint(boxdata(iptcl,3)) /= params%box )then
                        write(logfhandle,*) 'box_current: ', nint(boxdata(iptcl,3)), 'box in params: ', params%box
                        THROW_HARD('inconsistent box sizes in box files; exec_extract')
                    endif
                    ! update particle mask & movie index
                    if( box_inside(ldim, nint(boxdata(iptcl,1:2)), params%box) )oris_mask(iptcl) = .true.
                end do
                if( count(oris_mask) == 0 )then
                    ! no particles to extract
                    mics_mask(imic) = .false.
                    cycle
                endif
                ! fetch ctf info
                ctfparms      = o_mic%get_ctfvars()
                ctfparms%smpd = params%smpd
                if( o_mic%isthere('dfx') )then
                    if( .not.o_mic%isthere('cs') .or. .not.o_mic%isthere('kv') .or. .not.o_mic%isthere('fraca') )then
                        THROW_HARD('input lacks at least cs, kv or fraca; exec_extract')
                    endif
                endif
                ! output stack
                stack = trim(output_dir)//trim(EXTRACT_STK_FBODY)//trim(basename(mic_name))
                ! extract windows from integrated movie
                call micrograph%read(mic_name, 1)
                ! phase-flip micrograph
                if( cline%defined('ctf') )then
                    if( trim(params%ctf).eq.'flip' .and. o_mic%isthere('dfx') )then
                        ! phase flip micrograph
                        tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                        call micrograph%zero_edgeavg
                        call micrograph%fft
                        call tfun%apply_serial(micrograph, 'flip', ctfparms)
                        ! update stack ctf flag, mic flag unchanged
                        ctfparms%ctfflag = CTFFLAG_FLIP
                    endif
                endif
                ! filter out frequencies lower than the box can express to avoid aliasing
                call micrograph%ifft ! need to be here in case it was flipped
                ! write new stack
                stk_stats(1)   = huge(stk_stats(1))
                stk_stats(2)   = -stk_stats(1)
                stk_stats(3:4) = 0.
                cnt_stats = 0
                cnt       = 0
                do iptcl=1,nptcls ! loop over boxes
                    if( oris_mask(iptcl) )then
                        cnt = cnt + 1
                        ! extract the window
                        ptcl_pos = boxdata(iptcl,1:2)
                        call micrograph%window(nint(ptcl_pos), params%box, build%img, noutside)
                        if( params%pcontrast .eq. 'black' ) call build%img%neg()
                        call build%img%noise_norm(build%lmsk, sdev_noise)
                        ! keep track of stats
                        call build%img%stats(meanv, sddevv, maxv, minv, errout=l_err)
                        if( .not.l_err )then
                            cnt_stats = cnt_stats + 1
                            stk_stats(1) = min(stk_stats(1),minv)
                            stk_stats(2) = max(stk_stats(2),maxv)
                            stk_stats(3) = stk_stats(3) + meanv
                            stk_stats(4) = stk_stats(4) + sddevv**2.
                        endif
                        ! write
                        call build%img%write(trim(adjustl(stack)), cnt)
                    endif
                end do
                stk_stats(3) = stk_stats(3) / real(cnt_stats)
                stk_stats(4) = sqrt(stk_stats(4) / real(cnt_stats))
                call build%img%update_header_stats(trim(adjustl(stack)), stk_stats)
                ! IMPORT INTO PROJECT
                call build%spproj%add_stk(trim(adjustl(stack)), ctfparms)
                ! add box coordinates to ptcl2D field only & updates patch-based defocus
                l_ctfpatch = .false.
                if( o_mic%isthere('ctfdoc') )then
                    ctfdoc = o_mic%get_static('ctfdoc')
                    if( file_exists(ctfdoc) )then
                        call ctffit%read_doc(ctfdoc)
                        l_ctfpatch = .true.
                    endif
                endif
                do iptcl=1,nptcls
                    if( .not.oris_mask(iptcl) )cycle
                    ! update global counter
                    iptcl_glob = iptcl_glob + 1
                    ptcl_pos   = boxdata(iptcl,1:2)
                    ! updates particle position
                    call build%spproj%set_boxcoords(iptcl_glob, nint(ptcl_pos))
                    ! updates particle defocus
                    if( l_ctfpatch )then
                        ptcl_pos = ptcl_pos+1.+real(params%box/2) !  center
                        call ctffit%pix2polyvals(ptcl_pos(1),ptcl_pos(2), dfx,dfy)
                        call build%spproj%os_ptcl2D%set(iptcl_glob,'dfx',dfx)
                        call build%spproj%os_ptcl3D%set(iptcl_glob,'dfx',dfx)
                        call build%spproj%os_ptcl2D%set(iptcl_glob,'dfy',dfy)
                        call build%spproj%os_ptcl3D%set(iptcl_glob,'dfy',dfy)
                    endif
                end do
                ! clean
                call boxfile%kill
                call ctffit%kill
            enddo
            ! write
            call build%spproj%write
        endif
        ! end gracefully
        call o_mic%kill
        call o_tmp%kill
        call qsys_job_finished('simple_commander_preprocess :: exec_extract')
        call simple_end('**** SIMPLE_EXTRACT NORMAL STOP ****')

        contains

            function box_inside( ldim, coord, box ) result( inside )
                integer, intent(in) :: ldim(3), coord(2), box
                integer             :: fromc(2), toc(2)
                logical             :: inside
                if( params%outside .eq. 'yes' )then
                    inside = .true.
                    return
                endif
                fromc  = coord+1       ! compensate for the c-range that starts at 0
                toc    = fromc+(box-1) ! the lower left corner is 1,1
                inside = .true.        ! box is inside
                if( any(fromc < 1) .or. toc(1) > ldim(1) .or. toc(2) > ldim(2) ) inside = .false.
            end function box_inside
    end subroutine exec_extract

    subroutine exec_reextract_distr( self, cline )
        use simple_oris,  only: oris
        use simple_ori,   only: ori
        class(reextract_commander_distr), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline !< command line input
        type(parameters)                        :: params
        type(sp_project)                        :: spproj
        type(sp_project),           allocatable :: spproj_parts(:)
        type(qsys_env)                          :: qenv
        type(chash)                             :: job_descr
        type(ori)                               :: o_mic, o, o_tmp
        type(oris)                              :: os_stk
        type(chash),                allocatable :: part_params(:)
        character(len=LONGSTRLEN),  allocatable :: boxfiles(:), stktab(:), parts_fname(:)
        character(len=:),           allocatable :: mic_name, imgkind
        integer,                    allocatable :: parts(:,:)
        integer :: boxcoords(2)
        integer :: imic,i,nmics_tot,numlen,nmics,cnt,state,istk,nstks,ipart
        if( cline%defined('ctf') )then
            if( cline%get_carg('ctf').ne.'flip' .and. cline%get_carg('ctf').ne.'no' )then
                THROW_HARD('Only CTF=NO/FLIP are allowed')
            endif
        endif
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir',       'yes')
        if( .not. cline%defined('pcontrast') ) call cline%set('pcontrast', 'black')
        if( .not. cline%defined('oritype')   ) call cline%set('oritype',  'ptcl3D')
        call cline%set('nthr',1.)
        call params%new(cline)
        call cline%set('mkdir', 'no')
        ! read in integrated movies
        call spproj%read(params%projfile)
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
            if( o_mic%isthere('state') ) state = nint(o_mic%get('state'))
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
        call qenv%gen_scripts_and_schedule_jobs( job_descr, algnfbody=trim(ALGN_FBODY),&
            &part_params=part_params)
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
            call spproj%os_mic%new(nmics)
            ! transfer stacks
            cnt   = 0
            nstks = 0
            do ipart = 1,params%nparts
                do imic = 1,spproj_parts(ipart)%os_mic%get_noris()
                    cnt = cnt + 1
                    call spproj_parts(ipart)%os_mic%get_ori(imic,o_tmp)
                    call spproj%os_mic%set_ori(cnt,o_tmp)
                enddo
                call spproj_parts(ipart)%kill
                call spproj_parts(ipart)%read_segment('stk',parts_fname(ipart))
                nstks = nstks + spproj_parts(ipart)%os_stk%get_noris()
            enddo
            if( nstks /= nmics ) THROW_HARD('Inconstistent number of stacks in individual projects')
            ! generates stacks table
            call os_stk%new(nstks)
            allocate(stktab(nstks))
            cnt = 0
            do ipart = 1,params%nparts
                do istk = 1,spproj_parts(ipart)%os_stk%get_noris()
                    cnt = cnt + 1
                    call spproj_parts(ipart)%os_stk%get_ori(istk, o_tmp)
                    call os_stk%set_ori(cnt,o_tmp)
                    stktab(cnt) = os_stk%get_static(cnt,'stk')
                enddo
                call spproj_parts(ipart)%kill
            enddo
            ! import stacks into project
            call spproj%add_stktab(stktab,os_stk)
            ! transfer 2D parameters
            cnt = 0
            do ipart = 1,params%nparts
                call spproj_parts(ipart)%read_segment('ptcl2D',parts_fname(ipart))
                do i = 1,spproj_parts(ipart)%os_ptcl2D%get_noris()
                    cnt = cnt + 1
                    ! particles coordinates
                    call spproj_parts(ipart)%get_boxcoords(i, boxcoords)
                    call spproj%set_boxcoords(cnt, boxcoords)
                    ! search history & parameters
                    call spproj_parts(ipart)%os_ptcl2D%get_ori(i, o)
                    call spproj%os_ptcl2D%transfer_2Dparams(cnt, o)
                enddo
                call spproj_parts(ipart)%kill
            enddo
            ! transfer 3D parameters
            cnt = 0
            do ipart = 1,params%nparts
                call spproj_parts(ipart)%read_segment('ptcl3D',parts_fname(ipart))
                do i = 1,spproj_parts(ipart)%os_ptcl3D%get_noris()
                    cnt = cnt + 1
                    call spproj_parts(ipart)%os_ptcl3D%get_ori(i, o)
                    call spproj%os_ptcl3D%transfer_3Dparams(cnt, o)
                enddo
                call spproj_parts(ipart)%kill
            enddo
            ! clean-up
            call os_stk%kill
        endif
        ! final write
        call spproj%write
        ! clean-up
        call qsys_cleanup
        call spproj%kill
        deallocate(spproj_parts,part_params)
        call o_mic%kill
        call o%kill
        call o_tmp%kill
        call os_stk%kill
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
        use simple_ctf,   only: ctf
        class(reextract_commander), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline !< command line input
        type(parameters)              :: params
        type(sp_project)              :: spproj, spproj_in
        type(image)                   :: micrograph, img, mskimg
        type(ori)                     :: o_mic, o_stk, o_tmp
        type(ctf)                     :: tfun
        type(ctfparams)               :: ctfparms
        character(len=:), allocatable :: mic_name, imgkind
        logical,          allocatable :: pmsk(:,:,:), mic_mask(:), ptcl_mask(:)
        integer,          allocatable :: mic2stk_inds(:)
        character(len=LONGSTRLEN)     :: stack, rel_stack
        integer  :: nframes,imic,iptcl,nmics,prev_box,box_foo,cnt,nmics_tot,nptcls,stk_ind,cnt_stats
        integer :: prev_pos(2),new_pos(2),ishift(2),ldim(3),ldim_foo(3),noutside,fromp,top,istk
        real    :: prev_shift(2), shift2d(2), shift3d(2), minv,maxv,meanv,sddevv,stk_stats(4),sdev_noise
        logical :: l_3d, l_err
        call cline%set('mkdir','no')
        call params%new(cline)
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
        ldim_foo = 0
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
            call spproj_in%read_segment('ptcl2D', params%projfile)
            ! sanity checks
            do imic = 1,nmics_tot
                if( .not.mic_mask(imic) )cycle
                ! sanity checks
                call spproj_in%os_mic%get_ori(imic, o_mic)
                call o_mic%getter('intg', mic_name)
                if( .not.file_exists(mic_name) )cycle
                call find_ldim_nptcls(mic_name, ldim_foo, nframes )
                if( nframes > 1 ) THROW_HARD('multi-frame extraction not supported; exec_reextract')
                if( any(ldim == 0) ) ldim = ldim_foo
                stk_ind = mic2stk_inds(imic)
                call spproj_in%os_stk%get_ori(stk_ind, o_stk)
                fromp   = nint(o_stk%get('fromp'))
                top     = nint(o_stk%get('top'))
                do iptcl=fromp,top
                    if(.not.spproj_in%has_boxcoords(iptcl)) THROW_HARD('missing particle coordinates 2; exec_reextract')
                enddo
                box_foo = nint(o_stk%get('box'))
                if( prev_box == 0 ) prev_box = box_foo
                if( prev_box /= box_foo ) THROW_HARD('Inconsistent box size; exec_reextract')
            enddo
            if( .not.cline%defined('box') ) params%box = prev_box
            if( is_odd(params%box) ) THROW_HARD('Box size must be of even dimension! exec_extract')
            ! extraction
            write(logfhandle,'(A)')'>>> EXTRACTING... '
            call spproj_in%read_segment('ptcl3D', params%projfile)
            allocate(ptcl_mask(spproj_in%os_ptcl2D%get_noris()),source=.false.)
            ! set normalization radius & mask
            params%msk = RADFRAC_NORM_EXTRACT * real(params%box/2)
            call mskimg%disc([params%box,params%box,1], params%smpd, params%msk, pmsk)
            call mskimg%kill
            call micrograph%new([ldim(1),ldim(2),1], params%smpd)
            do imic = params%fromp,params%top
                if( .not.mic_mask(imic) ) cycle
                stk_ind = mic2stk_inds(imic)
                call spproj_in%os_mic%get_ori(imic, o_mic)
                call spproj_in%os_stk%get_ori(stk_ind, o_stk)
                call o_mic%getter('intg', mic_name)
                ctfparms = o_mic%get_ctfvars()
                fromp    = nint(o_stk%get('fromp'))
                top      = nint(o_stk%get('top'))
                call micrograph%read(mic_name)
                ! phase-flip micrograph. To-add as an option to not flip?
                if( ctfparms%ctfflag == CTFFLAG_FLIP )then
                    if( o_mic%isthere('dfx') )then
                        ! phase flip micrograph
                        tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
                        call micrograph%zero_edgeavg
                        call micrograph%fft
                        call tfun%apply_serial(micrograph, 'flip', ctfparms)
                    endif
                endif
                call micrograph%ifft ! need to be here in case it was flipped
                stack = trim(EXTRACT_STK_FBODY)//trim(basename(mic_name))
                ! particles extraction loop
                nptcls         = 0
                cnt_stats      = 0
                stk_stats(1)   = huge(stk_stats(1))
                stk_stats(2)   = -stk_stats(1)
                stk_stats(3:4) = 0.
                do iptcl=fromp,top
                    if( spproj_in%os_ptcl2D%get_state(iptcl) == 0 ) cycle
                    if( spproj_in%os_ptcl3D%get_state(iptcl) == 0 ) cycle
                    ! previous position & shift
                    call spproj_in%get_boxcoords(iptcl, prev_pos)
                    if( l_3d )then
                        prev_shift = spproj_in%os_ptcl3D%get_2Dshift(iptcl)
                    else
                        prev_shift = spproj_in%os_ptcl2D%get_2Dshift(iptcl)
                    endif
                    ! calc new position & shift
                    ishift  = nint(prev_shift)
                    new_pos = prev_pos - ishift
                    if( prev_box /= params%box ) new_pos = new_pos + (prev_box-params%box)/2
                    if( box_inside(ldim, new_pos, params%box) )then
                        ! included
                        nptcls = nptcls+1
                        ptcl_mask(iptcl) = .true.
                        ! updates picking position
                        call spproj_in%set_boxcoords(iptcl, new_pos)
                        ! updates shifts
                        if( l_3d )then
                            shift2d = spproj_in%os_ptcl2D%get_2Dshift(iptcl) - real(ishift)
                            shift3d = prev_shift - real(ishift)
                        else
                            shift2d = prev_shift - real(ishift)
                            shift3d = spproj_in%os_ptcl3D%get_2Dshift(iptcl) - real(ishift)
                        endif
                        call spproj_in%os_ptcl2D%set_shift(iptcl, shift2d)
                        call spproj_in%os_ptcl3D%set_shift(iptcl, shift3d)
                        call spproj_in%os_ptcl2D%set(iptcl,'state',1.)
                        call spproj_in%os_ptcl3D%set(iptcl,'state',1.)
                        ! extracts
                        call micrograph%window(new_pos, params%box, img, noutside)
                        if( params%pcontrast .eq. 'black' ) call img%neg()
                        call img%noise_norm(pmsk, sdev_noise)
                        call img%write(trim(adjustl(stack)), nptcls)
                        ! keep track of stats
                        call img%stats(meanv, sddevv, maxv, minv, errout=l_err)
                        if( .not.l_err )then
                            cnt_stats = cnt_stats + 1
                            stk_stats(1) = min(stk_stats(1),minv)
                            stk_stats(2) = max(stk_stats(2),maxv)
                            stk_stats(3) = stk_stats(3) + meanv
                            stk_stats(4) = stk_stats(4) + sddevv**2.
                        endif
                    else
                        ! excluded
                        call spproj_in%os_ptcl2D%set(iptcl,'state',0.)
                        call spproj_in%os_ptcl3D%set(iptcl,'state',0.)
                    endif
                enddo
                if( nptcls == 0 )then
                    ! all particles in this micrograph excluded
                    call spproj_in%os_stk%set(stk_ind,'state',0.)
                    call spproj_in%os_mic%set(imic,'state',0.)
                    mic_mask(imic) = .false.
                    mic2stk_inds(imic) = 0
                else
                    ! updates header, size, stack & removes box file
                    stk_stats(3) = stk_stats(3) / real(cnt_stats)
                    stk_stats(4) = sqrt(stk_stats(4) / real(cnt_stats))
                    call img%update_header_stats(trim(adjustl(stack)), stk_stats)
                    call make_relativepath(CWD_GLOB, stack, rel_stack)
                    call spproj_in%os_stk%set(stk_ind,'stk',rel_stack)
                    call spproj_in%os_stk%set(stk_ind,'box', real(params%box))
                    call spproj_in%os_mic%set(imic,'nptcls',real(nptcls))
                    call spproj_in%os_mic%delete_entry(imic,'boxfile')
                endif
            enddo
        endif
        ! OUTPUT
        call spproj%read_non_data_segments(params%projfile)
        call spproj%projinfo%set(1,'projname', get_fbody(params%outfile,METADATA_EXT,separator=.false.))
        call spproj%projinfo%set(1,'projfile', params%outfile)
        nmics = count(mic_mask)
        ! transfer mics & stk
        call spproj%os_mic%new(nmics)
        call spproj%os_stk%new(nmics)
        nptcls = count(ptcl_mask)
        cnt = 0
        do imic = params%fromp,params%top
            if( .not.mic_mask(imic) )cycle
            cnt = cnt+1
            call spproj_in%os_mic%get_ori(imic, o_tmp)
            call spproj%os_mic%set_ori(cnt, o_tmp)
            stk_ind = mic2stk_inds(imic)
            call spproj_in%os_stk%get_ori(stk_ind, o_tmp)
            call spproj%os_stk%set_ori(cnt, o_tmp)
        enddo
        ! transfer particles
        nptcls = count(ptcl_mask)
        call spproj%os_ptcl2D%new(nptcls)
        call spproj%os_ptcl3D%new(nptcls)
        cnt = 0
        do iptcl = 1,size(ptcl_mask)
            if( .not.ptcl_mask(iptcl) )cycle
            cnt = cnt+1
            call spproj_in%os_ptcl2D%get_ori(iptcl, o_tmp)
            call spproj%os_ptcl2D%set_ori(cnt, o_tmp)
            call spproj_in%os_ptcl3D%get_ori(iptcl, o_tmp)
            call spproj%os_ptcl3D%set_ori(cnt, o_tmp)
        enddo
        call spproj_in%kill
        ! final write
        call spproj%write(params%outfile)
        write(logfhandle,'(A,I8)')'>>> RE-EXTRACTED  PARTICLES: ', nptcls
        ! end gracefully
        call qsys_job_finished('simple_commander_preprocess :: exec_reextract')
        call o_mic%kill
        call o_stk%kill
        call o_tmp%kill
        call simple_end('**** SIMPLE_REEXTRACT NORMAL STOP ****')

        contains

            function box_inside( ldim, coord, box ) result( inside )
                integer, intent(in) :: ldim(3), coord(2), box
                integer             :: fromc(2), toc(2)
                logical             :: inside
                fromc  = coord+1       ! compensate for the c-range that starts at 0
                toc    = fromc+(box-1) ! the lower left corner is 1,1
                inside = .true.        ! box is inside
                if( any(fromc < 1) .or. toc(1) > ldim(1) .or. toc(2) > ldim(2) ) inside = .false.
            end function box_inside
    end subroutine exec_reextract

    subroutine exec_pick_extract_stream( self, cline )
        class(pick_extract_commander_stream), intent(inout) :: self
        class(cmdline),                             intent(inout) :: cline
        integer,          parameter   :: WAIT_WATCHER        = 10    ! seconds prior to new stack detection
        !integer,          parameter   :: ORIGPROJ_WRITEFREQ  = 600   ! 10mins, Frequency at which the original project file should be updated
        type(parameters)                    :: params
        type(cmdline)                       :: cline_pick_extract, cline_make_pickrefs
        type(sp_project)                    :: orig_proj, stream_proj
        type(qsys_env)                      :: qenv
        class(cmdline),         allocatable :: completed_jobs_clines(:)
        character(LONGSTRLEN),  allocatable :: spproj_list(:)
        character(len=:),       allocatable :: spproj_list_fname, output_dir, output_dir_picker, output_dir_extract, projfname
        integer :: iter, origproj_time, tnow, iproj, icline, nptcls, prev_stacksz, stacksz
        integer :: last_injection, n_spprojs, n_spprojs_prev, n_newspprojs, nmics
        if( cline%defined('refs') .and. cline%defined('vol1') )then
            THROW_HARD('REFS and VOL1 cannot be both provided!')
        endif
        if( .not.cline%defined('refs') .and. .not.cline%defined('vol1') )then
            THROW_HARD('one of REFS and VOL1 must be provided!')
        endif
        if( .not. cline%defined('mkdir')     ) call cline%set('mkdir',       'yes')
        if( .not. cline%defined('pcontrast') ) call cline%set('pcontrast', 'black')
        ! output command line executed
        write(logfhandle,'(a)') '>>> COMMAND LINE EXECUTED'
        write(logfhandle,*) trim(cmdline_glob)
        ! set oritype & defaults
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'mic')
        call cline%set('stream','yes')
        call cline%set('numlen', 5.)
        call params%new(cline)
        params_glob%split_mode = 'stream'
        params_glob%ncunits    = params%nparts
        call cline%set('mkdir', 'no')
        if( .not.file_exists(params%projfile) )THROW_HARD('project file: '//trim(params%projfile)//' does not exist!')
        spproj_list_fname = filepath(trim(params%dir_target),trim(STREAM_SPPROJFILES))
        ! read project info
        call orig_proj%read(params%projfile)
        ! setup the environment for distributed execution
        call qenv%new(1,stream=.true.)
        ! output directories
        output_dir = PATH_HERE
        output_dir_picker  = filepath(trim(output_dir), trim(DIR_PICKER))
        output_dir_extract = filepath(trim(output_dir), trim(DIR_EXTRACT))
        call simple_mkdir(trim(output_dir),errmsg="commander_stream_wflows :: exec_pick_extract_stream;  ")
        call simple_mkdir(trim(output_dir_picker),errmsg="commander_stream_wflows :: exec_pick_extract_stream;  ")
        call simple_mkdir(trim(output_dir_extract),errmsg="commander_stream_wflows :: exec_pick_extract_stream;  ")
        ! init command-lines
        cline_pick_extract  = cline
        cline_make_pickrefs = cline
        call cline_pick_extract%set('prg', 'pick_extract')
        call cline_pick_extract%delete('projname')
        ! prepares picking references
        call cline_make_pickrefs%set('prg','make_pickrefs')
        call cline_make_pickrefs%set('stream','no')
        call qenv%exec_simple_prg_in_queue(cline_make_pickrefs, 'MAKE_PICKREFS_FINISHED')
        call cline_pick_extract%set('refs', trim(PICKREFS)//params%ext)
        call cline_pick_extract%delete('vol1')
        ! wait for the first stacks
        last_injection  = simple_gettime()
        origproj_time   = last_injection
        prev_stacksz    = 0
        n_spprojs       = 0
        n_spprojs_prev  = 0
        do iter = 1,999999
            tnow = simple_gettime()
            if(tnow-last_injection > params%time_inactive .and. stacksz==0)then
                write(logfhandle,*)'>>> TIME LIMIT WITHOUT NEW MICROGRAPHS REACHED'
                exit
            endif
            if( file_exists(spproj_list_fname) )then
                if( .not.is_file_open(spproj_list_fname) )then
                    call read_filetable(spproj_list_fname, spproj_list)
                    if( allocated(spproj_list) )n_spprojs = size(spproj_list)
                    n_newspprojs = n_spprojs - n_spprojs_prev
                    if( n_newspprojs > 0 )then
                        ! copy projects and add to processing stack
                        do iproj = n_spprojs_prev+1, n_spprojs
                            call stream_proj%read(spproj_list(iproj))
                            projfname  = filepath(PATH_HERE, basename(spproj_list(iproj)))
                            call stream_proj%write(projfname)
                            call cline_pick_extract%set('projfile', trim(projfname))
                            call qenv%qscripts%add_to_streaming( cline_pick_extract )
                            call stream_proj%kill
                            deallocate(projfname)
                        enddo
                        n_spprojs_prev = n_spprojs
                    endif
                    ! streaming scheduling
                    call qenv%qscripts%schedule_streaming( qenv%qdescr )
                    stacksz = qenv%qscripts%get_stacksz()
                    if( stacksz .ne. prev_stacksz )then
                        prev_stacksz = stacksz
                        write(logfhandle,'(A,I5)')'>>> MICROGRAPHS TO PROCESS: ', stacksz
                    endif
                    ! completed jobs update the current project
                    if( qenv%qscripts%get_done_stacksz() > 0 )then
                        call qenv%qscripts%get_stream_done_stack( completed_jobs_clines )
                        do icline=1,size(completed_jobs_clines)
                            projfname = completed_jobs_clines(icline)%get_carg('projfile')
                            call stream_proj%read( projfname )
                            call orig_proj%append_project(stream_proj, 'mic')
                            call orig_proj%append_project(stream_proj, 'stk')
                            call stream_proj%kill()
                            deallocate(projfname)
                        enddo
                        nptcls = orig_proj%get_nptcls()
                        nmics  = orig_proj%get_nintgs()
                        write(logfhandle,'(A,I8)')'>>> NEW MICROGRAPHS COUNT: ', nmics
                        write(logfhandle,'(A,I8)')'>>> NEW PARTICLES   COUNT: ', nptcls
                        call orig_proj%write
                        deallocate(completed_jobs_clines)
                    endif
                endif
            endif
            call simple_sleep(WAIT_WATCHER)
        enddo
        ! cleanup
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_PICK_EXTRACT_STREAM NORMAL STOP ****')
    end subroutine exec_pick_extract_stream

    subroutine exec_pick_extract( self, cline )
        use simple_sp_project,          only: sp_project
        use simple_picker_iter,         only: picker_iter
        use simple_binoris_io,          only: binwrite_oritab
        class(pick_extract_commander), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters)              :: params
        type(ori)                     :: o_mic
        type(picker_iter)             :: piter
        type(extract_commander)       :: xextract
        type(cmdline)                 :: cline_extract
        type(sp_project)              :: spproj
        character(len=:), allocatable :: micname, output_dir_picker, fbody, output_dir_extract
        character(len=LONGSTRLEN)     :: boxfile
        integer :: fromto(2), imic, ntot, nptcls_out, state
        ! set oritype
        call cline%set('oritype', 'mic')
        ! parse parameters
        call params%new(cline)
        if( params%scale > 1.01 )then
            THROW_HARD('scale cannot be > 1; exec_preprocess')
        endif
        if( params%stream.ne.'yes' ) THROW_HARD('streaming only application')
        ! read in movies
        call spproj%read( params%projfile )
        if( spproj%get_nintgs() == 0 ) THROW_HARD('no micrograph to process!')
        ! output directories
        if( params%stream.eq.'yes' )then
            output_dir_picker  = trim(DIR_PICKER)
            output_dir_extract = trim(DIR_EXTRACT)
            call simple_mkdir(output_dir_picker, errmsg="commander_preprocess :: preprocess; ")
            call simple_mkdir(output_dir_extract,errmsg="commander_preprocess :: preprocess; ")
        else
            output_dir_picker  = PATH_HERE
            output_dir_extract = PATH_HERE
        endif
        ! command lines
        cline_extract = cline
        call cline_extract%set('dir', trim(output_dir_extract))
        call cline_extract%set('pcontrast', params%pcontrast)
        if( cline%defined('box_extract') )call cline_extract%set('box', real(params%box_extract))
        call cline%delete('box')
        call cline_extract%delete('box_extract')
         ! file name
        if( cline%defined('fbody') )then
            fbody = trim(params%fbody)
        else
            fbody = ''
        endif
        ! range
        if( params%stream.eq.'yes' )then
            fromto(:) = 1
        else
            fromto(:) = [params%fromp, params%top]
        endif
        ntot = fromto(2) - fromto(1) + 1
        ! loop over exposures (movies)
        do imic = fromto(1),fromto(2)
            ! fetch movie orientation
            call spproj%os_mic%get_ori(imic, o_mic)
            ! sanity check
            state = 1
            if( o_mic%isthere('state') ) state = nint(o_mic%get('state'))
            if( state == 0 ) cycle
            if( .not.o_mic%isthere('intg')   )cycle
            call o_mic%getter('intg', micname)
            if( .not.file_exists(micname)) cycle
            ! picker
            params_glob%lp = max(params%fny, params%lp_pick)
            call piter%iterate(cline, micname, boxfile, nptcls_out, output_dir_picker)
            call o_mic%set_boxfile(boxfile, nptcls=nptcls_out)
            ! update project
            call spproj%os_mic%set_ori(imic, o_mic)
            call spproj%write_segment_inside(params%oritype)
            ! extract particles
            call xextract%execute(cline_extract)
            call spproj%kill
        end do
        if( params%stream .eq. 'yes' )then
            ! nothing to do, extract did it
        else
            call binwrite_oritab(params%outfile, spproj, spproj%os_mic, fromto, isegment=MIC_SEG)
        endif
        ! end gracefully
        call qsys_job_finished(  'simple_commander_preprocess :: exec_pick_extract' )
        call o_mic%kill
        call simple_end('**** SIMPLE_PICK_EXTRACT NORMAL STOP ****')
    end subroutine exec_pick_extract

    subroutine exec_make_pickrefs( self, cline )
        use simple_sym,                 only: sym
        use simple_projector_hlev,      only: reproject
        class(make_pickrefs_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters)              :: params
        type(oris)                    :: os
        type(sym)                     :: pgrpsyms
        type(image)                   :: ref3D, ref2D
        type(image),      allocatable :: projs(:)
        integer, parameter :: NREFS=100, NPROJS=20
        real    :: ang, rot, smpd_here
        integer :: nrots, iref, irot, ldim(3), ldim_here(3), ifoo, ncavgs, icavg
        integer :: cnt, norefs
        ! error check
        if( cline%defined('refs') .and. cline%defined('vol1') )then
            THROW_HARD('REFS and VOL1 cannot be both provided!')
        endif
        if( .not.cline%defined('refs') .and. .not.cline%defined('vol1') )then
            THROW_HARD('One of REFS, VOL1 & PROJFILE must be informed!')
        endif
        ! set defaults
        call cline%set('oritype', 'mic')
        if( .not. cline%defined('pcontrast') ) call cline%set('pcontrast','black')
        ! parse parameters
        call params%new(cline)
        if( params%stream.eq.'yes' ) THROW_HARD('not a streaming application')
        if( .not. cline%defined('pgrp') ) params%pgrp = 'd1' ! only northern hemisphere
        ! point-group object
        call pgrpsyms%new(trim(params%pgrp))
        if( cline%defined('vol1') )then
            ! find logical dimension & read reference volume
            call find_ldim_nptcls(params%vols(1), ldim_here, ifoo, smpd=smpd_here)
            if( smpd_here < 0.01 ) THROW_HARD('Invalid sampling distance for the volume (should be in MRC format)')
            call ref3D%new(ldim_here, smpd_here)
            call ref3D%read(params%vols(1))
            call scale_ref(ref3D, params%smpd)
            ! make projection directions
            call os%new(NPROJS)
            call pgrpsyms%build_refspiral(os)
            ! generate reprojections
            projs  = reproject(ref3D, os)
            nrots  = nint( real(NREFS)/real(NPROJS) )
            norefs = NPROJS
        else if( cline%defined('refs') )then
            ! read selected cavgs
            call find_ldim_nptcls(params%refs, ldim_here, ncavgs, smpd=smpd_here)
            if( smpd_here < 0.01 ) THROW_HARD('Invalid sampling distance for the cavgs (should be in MRC format)')
            ldim_here(3) = 1
            allocate( projs(ncavgs) )
            do icavg=1,ncavgs
                call projs(icavg)%new(ldim_here, smpd_here)
                call projs(icavg)%read(params%refs, icavg)
                call scale_ref(projs(icavg), params%smpd)
            end do
            nrots  = nint( real(NREFS)/real(ncavgs) )
            norefs = ncavgs
        else
            THROW_HARD('Missing volume / cavgs for creating picking references!')
        endif
        ! expand in in-plane rotation and write to file
        if( nrots > 1 )then
            call ref2D%new([ldim(1),ldim(2),1], params%smpd)
            ang = 360./real(nrots)
            rot = 0.
            cnt = 0
            do iref=1,norefs
                do irot=1,nrots
                    cnt = cnt + 1
                    call projs(iref)%rtsq(rot, 0., 0., ref2D)
                    if(params%pcontrast .eq. 'black') call ref2D%neg
                    call ref2D%write(trim(PICKREFS)//params%ext, cnt)
                    rot = rot + ang
                end do
            end do
        else
            ! should never happen
            do iref=1,norefs
                if(params%pcontrast .eq. 'black') call projs(iref)%neg
                call projs(iref)%write(trim(PICKREFS)//params%ext, iref)
            end do
        endif
        ! end gracefully
        call simple_touch('MAKE_PICKREFS_FINISHED', errmsg='In: commander_preprocess::exec_make_pickrefs')
        call simple_end('**** SIMPLE_MAKE_PICKREFS NORMAL STOP ****')
        contains

            subroutine scale_ref(refimg, smpd_target)
                class(image), intent(inout) :: refimg
                real,         intent(in)    :: smpd_target
                type(image) :: targetimg
                integer     :: ldim_ref(3), ldim_target(3)
                real        :: smpd_ref, scale
                smpd_ref = refimg%get_smpd()
                scale    = smpd_target / smpd_ref
                if( is_equal(scale,1.) )then
                    ldim = ldim_here
                    return
                endif
                ldim_ref       = refimg%get_ldim()
                ldim_target(1) = round2even(real(ldim_ref(1))/scale)
                ldim_target(2) = ldim_target(1)
                ldim_target(3) = 1
                if( refimg%is_3d() )ldim_target(3) = ldim_target(1)
                call refimg%norm
                if( scale > 1. )then
                    ! downscaling
                    call refimg%fft
                    call refimg%clip_inplace(ldim_target)
                    call refimg%ifft
                else
                    call targetimg%new(ldim_target, smpd_target)
                    call refimg%fft
                    call refimg%pad(targetimg, backgr=0.)
                    call targetimg%ifft
                    refimg = targetimg
                    call targetimg%kill
                endif
                ! updates dimensions
                ldim = ldim_target
            end subroutine

    end subroutine exec_make_pickrefs

end module simple_commander_preprocess
