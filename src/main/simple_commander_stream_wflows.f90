! concrete commander: stream processing routines
module simple_commander_stream_wflows
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_parameters,     only: parameters, params_glob
use simple_sp_project,     only: sp_project
use simple_qsys_env,       only: qsys_env
use simple_qsys_funs
implicit none

! public :: preprocess_commander_stream
! public :: cluster2D_commander_stream
! public :: pick_extract_commander_stream
private
#include "simple_local_flags.inc"

! type, extends(commander_base) :: preprocess_commander_stream
!   contains
!     procedure :: execute      => exec_preprocess_stream
! end type preprocess_commander_stream
! type, extends(commander_base) :: cluster2D_commander_stream
!   contains
!     procedure :: execute      => exec_cluster2D_stream
! end type cluster2D_commander_stream
! type, extends(commander_base) :: pick_extract_commander_stream
!   contains
!     procedure :: execute      => exec_pick_extract_stream
! end type pick_extract_commander_stream

contains

    ! subroutine exec_preprocess_stream( self, cline )
    !     use simple_moviewatcher,         only: moviewatcher
    !     use simple_qsys_funs,            only: qsys_cleanup
    !     use simple_commander_preprocess, only: preprocess_commander
    !     class(preprocess_commander_stream), intent(inout) :: self
    !     class(cmdline),                     intent(inout) :: cline
    !     type(parameters)                       :: params
    !     integer,                   parameter   :: SHORTTIME = 60   ! folder watched every minute
    !     integer,                   parameter   :: LONGTIME  = 600  ! time lag after which a movie is processed
    !     class(cmdline),            allocatable :: completed_jobs_clines(:)
    !     type(qsys_env)                         :: qenv
    !     type(cmdline)                          :: cline_make_pickrefs
    !     type(moviewatcher)                     :: movie_buff
    !     type(sp_project)                       :: spproj, stream_spproj
    !     character(len=LONGSTRLEN), allocatable :: movies(:)
    !     character(len=:),          allocatable :: output_dir, output_dir_ctf_estimate, output_dir_picker
    !     character(len=:),          allocatable :: output_dir_motion_correct, output_dir_extract, stream_spprojfile
    !     character(len=LONGSTRLEN)              :: movie
    !     integer                                :: nmovies, imovie, stacksz, prev_stacksz, iter, icline
    !     integer                                :: nptcls, nptcls_prev, nmovs, nmovs_prev
    !     logical                                :: l_pick
    !     if( .not. cline%defined('oritype')         ) call cline%set('oritype',        'mic')
    !     ! mnotion correction
    !     if( .not. cline%defined('trs')             ) call cline%set('trs',              30.)
    !     if( .not. cline%defined('lpstart')         ) call cline%set('lpstart',           8.)
    !     if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',            5.)
    !     if( .not. cline%defined('bfac')            ) call cline%set('bfac',            150.)
    !     if( .not. cline%defined('groupframes')     ) call cline%set('groupframes',  'patch')
    !     ! ctf estimation
    !     if( .not. cline%defined('pspecsz')         ) call cline%set('pspecsz',         512.)
    !     if( .not. cline%defined('hp_ctf_estimate') ) call cline%set('hp_ctf_estimate',  30.)
    !     if( .not. cline%defined('lp_ctf_estimate') ) call cline%set('lp_ctf_estimate',   5.)
    !     ! picking
    !     if( .not. cline%defined('lp_pick')         ) call cline%set('lp_pick',          20.)
    !     ! extraction
    !     if( .not. cline%defined('pcontrast')       ) call cline%set('pcontrast',    'black')
    !     if( cline%defined('refs') .and. cline%defined('vol1') )then
    !         THROW_HARD('REFS and VOL1 cannot be both provided!')
    !     endif
    !     call cline%set('numlen', 5.)
    !     call cline%set('stream','yes')
    !     ! master parameters
    !     call params%new(cline)
    !     params_glob%split_mode = 'stream'
    !     params_glob%ncunits    = params%nparts
    !     call cline%set('mkdir', 'no')
    !     call cline%set('prg',   'preprocess')
    !     if( cline%defined('dir_prev') .and. .not.file_exists(params%dir_prev) )then
    !         THROW_HARD('Directory '//trim(params%dir_prev)//' does not exist!')
    !     endif
    !     ! read in movies
    !     call spproj%read( params%projfile )
    !     ! sanity check
    !     if( spproj%os_mic%get_noris() /= 0 )then
    !         THROW_HARD('PREPROCESS_STREAM must always start from an empty project (eg from root proejct folder)')
    !     endif
    !     ! picking
    !     l_pick = .false.
    !     if( cline%defined('refs') .or. cline%defined('vol1') ) l_pick = .true.
    !     ! output directories
    !     output_dir = PATH_HERE
    !     output_dir_ctf_estimate   = filepath(trim(output_dir), trim(DIR_CTF_ESTIMATE))
    !     output_dir_motion_correct = filepath(trim(output_dir), trim(DIR_MOTION_CORRECT))
    !     call simple_mkdir(output_dir_ctf_estimate,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
    !     call simple_mkdir(output_dir_motion_correct,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
    !     if( l_pick )then
    !         output_dir_picker  = filepath(trim(output_dir), trim(DIR_PICKER))
    !         output_dir_extract = filepath(trim(output_dir), trim(DIR_EXTRACT))
    !         call simple_mkdir(output_dir_picker,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
    !         call simple_mkdir(output_dir_extract,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
    !     endif
    !     ! setup the environment for distributed execution
    !     call qenv%new(1,stream=.true.)
    !     ! prepares picking references
    !     if( l_pick )then
    !         cline_make_pickrefs = cline
    !         call cline_make_pickrefs%set('prg','make_pickrefs')
    !         call cline_make_pickrefs%set('stream','no')
    !         call qenv%exec_simple_prg_in_queue(cline_make_pickrefs, 'MAKE_PICKREFS_FINISHED')
    !         call cline%set('refs', trim(PICKREFS)//params%ext)
    !         call cline%delete('vol1')
    !         write(logfhandle,'(A)')'>>> PREPARED PICKING TEMPLATES'
    !     endif
    !     ! movie watcher init
    !     movie_buff = moviewatcher(LONGTIME)
    !     ! import previous runs
    !     call import_prev_streams
    !     ! start watching
    !     prev_stacksz = 0
    !     nmovies      = 0
    !     iter         = 0
    !     do
    !         if( file_exists(trim(TERM_STREAM)) )then
    !             write(logfhandle,'(A)')'>>> TERMINATING PREPROCESS STREAM'
    !             exit
    !         endif
    !         do while( file_exists(trim(PAUSE_STREAM)) )
    !             if( file_exists(trim(TERM_STREAM)) ) exit
    !             call write_singlelineoftext(PAUSE_STREAM, 'PAUSED')
    !             write(logfhandle,'(A,A)')'>>> PREPROCES STREAM PAUSED ',cast_time_char(simple_gettime())
    !             call simple_sleep(SHORTTIME)
    !         enddo
    !         iter = iter + 1
    !         call movie_buff%watch( nmovies, movies )
    !         ! append movies to processing stack
    !         if( nmovies > 0 )then
    !             do imovie = 1, nmovies
    !                 movie = trim(adjustl(movies(imovie)))
    !                 if( .not.file_exists(movie) )cycle ! petty triple checking
    !                 call create_individual_project
    !                 call qenv%qscripts%add_to_streaming( cline )
    !             enddo
    !         endif
    !         ! stream scheduling
    !         call qenv%qscripts%schedule_streaming( qenv%qdescr )
    !         stacksz = qenv%qscripts%get_stacksz()
    !         if( stacksz .ne. prev_stacksz )then
    !             prev_stacksz = stacksz
    !             write(logfhandle,'(A,I5)')'>>> MOVIES/MICROGRAPHS TO PROCESS: ', stacksz
    !         endif
    !         ! completed jobs update the current project
    !         if( qenv%qscripts%get_done_stacksz() > 0 )then
    !             ! append new processed movies to project
    !             call qenv%qscripts%get_stream_done_stack( completed_jobs_clines )
    !             nptcls_prev = spproj%get_nptcls()
    !             nmovs_prev  = spproj%os_mic%get_noris()
    !             do icline=1,size(completed_jobs_clines)
    !                 stream_spprojfile = completed_jobs_clines(icline)%get_carg('projfile')
    !                 call stream_spproj%read( stream_spprojfile )
    !                 call spproj%append_project(stream_spproj, 'mic')
    !                 if( l_pick )then
    !                     call spproj%append_project(stream_spproj, 'stk')
    !                 endif
    !                 call stream_spproj%kill()
    !                 deallocate(stream_spprojfile)
    !             enddo
    !             nptcls = spproj%get_nptcls()
    !             nmovs  = spproj%os_mic%get_noris()
    !             write(logfhandle,'(A,I5)')'>>> # MOVIES PROCESSED:    ',nmovs
    !             if( l_pick )then
    !                 write(logfhandle,'(A,I8)')'>>> # PARTICLES EXTRACTED: ',nptcls
    !             endif
    !             ! update for 2d streaming
    !             call update_projects_list
    !             deallocate(completed_jobs_clines)
    !             ! exit for trial runs
    !             if( cline%defined('nmovies_trial') )then
    !                 if( nmovs  >= params%nmovies_trial )then
    !                     write(logfhandle,'(A)')'>>> TRIAL # OF MOVIES REACHED'
    !                     exit
    !                 endif
    !             endif
    !             if( cline%defined('nptcls_trial') )then
    !                 if( l_pick .and. (nptcls >= params%nptcls_trial) )then
    !                     write(logfhandle,'(A)')'>>> TRIAL # OF PARTICLES REACHED'
    !                     exit
    !                 endif
    !             endif
    !             ! write
    !             call spproj%write
    !         else
    !             ! wait
    !             call simple_sleep(SHORTTIME)
    !         endif
    !     end do
    !     ! termination
    !     call spproj%write
    !     call spproj%kill
    !     ! cleanup
    !     call qsys_cleanup
    !     ! end gracefully
    !     call simple_end('**** SIMPLE_PREPROCESS_STREAM NORMAL STOP ****')
    !     contains
    !
    !         subroutine update_projects_list
    !             type(cmdline) :: cline_mov
    !             character(len=:),          allocatable :: fname, abs_fname
    !             character(len=LONGSTRLEN), allocatable :: old_fnames(:), fnames(:)
    !             integer :: i, n_spprojs, n_old
    !             n_spprojs = size(completed_jobs_clines)
    !             if( n_spprojs == 0 )return
    !             if( file_exists(STREAM_SPPROJFILES) )then
    !                 ! append
    !                 call read_filetable(STREAM_SPPROJFILES, old_fnames)
    !                 n_old = size(old_fnames)
    !                 allocate(fnames(n_spprojs+n_old))
    !                 fnames(1:n_old) = old_fnames(:)
    !                 do i=1,n_spprojs
    !                     cline_mov = completed_jobs_clines(i)
    !                     fname     = trim(cline_mov%get_carg('projfile'))
    !                     abs_fname = simple_abspath(fname, errmsg='preprocess_stream :: update_projects_list 1')
    !                     fnames(n_old+i) = trim(abs_fname)
    !                     deallocate(abs_fname)
    !                 enddo
    !             else
    !                 ! first write
    !                 allocate(fnames(n_spprojs))
    !                 do i=1,n_spprojs
    !                     cline_mov = completed_jobs_clines(i)
    !                     fname     = trim(cline_mov%get_carg('projfile'))
    !                     abs_fname = simple_abspath(fname, errmsg='preprocess_stream :: update_projects_list 2')
    !                     fnames(i) = trim(abs_fname)
    !                     deallocate(abs_fname)
    !                 enddo
    !             endif
    !             call write_filetable(STREAM_SPPROJFILES, fnames)
    !         end subroutine update_projects_list
    !
    !         subroutine create_individual_project
    !             type(sp_project)              :: spproj_here
    !             type(cmdline)                 :: cline_here
    !             type(ctfparams)               :: ctfvars
    !             character(len=STDLEN)         :: ext, movie_here
    !             character(len=LONGSTRLEN)     :: projname, projfile
    !             movie_here = basename(trim(movie))
    !             ext        = fname2ext(trim(movie_here))
    !             projname   = 'preprocess_'//trim(get_fbody(trim(movie_here), trim(ext)))
    !             projfile   = trim(projname)//trim(METADATA_EXT)
    !             call cline_here%set('projname', trim(projname))
    !             call cline_here%set('projfile', trim(projfile))
    !             call spproj_here%update_projinfo(cline_here)
    !             spproj_here%compenv  = spproj%compenv
    !             spproj_here%jobproc  = spproj%jobproc
    !             ctfvars%ctfflag      = CTFFLAG_YES
    !             ctfvars%smpd         = params%smpd
    !             ctfvars%cs           = params%cs
    !             ctfvars%kv           = params%kv
    !             ctfvars%fraca        = params%fraca
    !             ctfvars%l_phaseplate = params%phaseplate.eq.'yes'
    !             call spproj_here%add_single_movie(trim(movie), ctfvars)
    !             call spproj_here%write
    !             call spproj_here%kill
    !             call cline%set('projname', trim(projname))
    !             call cline%set('projfile', trim(projfile))
    !         end subroutine create_individual_project
    !
    !         !>  import previous run to the current project based on past single project files
    !         subroutine import_prev_streams
    !             use simple_ori, only: ori
    !             type(ori) :: o, o_stk
    !             character(len=LONGSTRLEN), allocatable :: sp_files(:)
    !             character(len=:), allocatable :: mic, mov
    !             integer :: iproj,nprojs,cnt
    !             logical :: err
    !             if( .not.cline%defined('dir_prev') ) return
    !             call simple_list_files(trim(params%dir_prev)//'/preprocess_*.simple', sp_files)
    !             nprojs = size(sp_files)
    !             cnt    = 0
    !             do iproj = 1,nprojs
    !                 err = .false.
    !                 call stream_spproj%read( sp_files(iproj) )
    !                 if( stream_spproj%os_mic%get_noris() /= 1 )then
    !                     THROW_WARN('Ignoring previous project'//trim(sp_files(iproj)))
    !                     cycle
    !                 endif
    !                 call stream_spproj%os_mic%get_ori(1, o)
    !                 ! import mic segment
    !                 call movefile2folder('intg', output_dir_motion_correct, o, err)
    !                 if( err ) cycle
    !                 call movefile2folder('forctf', output_dir_motion_correct, o, err)
    !                 call movefile2folder('thumb', output_dir_motion_correct, o, err)
    !                 call movefile2folder('mc_starfile', output_dir_motion_correct, o, err)
    !                 call movefile2folder('mceps', output_dir_motion_correct, o, err)
    !                 call movefile2folder('ctfjpg', output_dir_ctf_estimate, o, err)
    !                 call movefile2folder('ctfdoc', output_dir_ctf_estimate, o, err)
    !                 if( .not.l_pick )then
    !                     ! import mic segment
    !                     call stream_spproj%os_mic%set_ori(1, o)
    !                     call spproj%append_project(stream_spproj, 'mic')
    !                 else
    !                     ! import mic & stk segment
    !                     call movefile2folder('boxfile', output_dir_picker, o, err)
    !                     call stream_spproj%os_mic%set_ori(1, o)
    !                     call spproj%append_project(stream_spproj, 'mic')
    !                     if( .not.err )then
    !                         if( stream_spproj%os_stk%get_noris() == 1 )then
    !                             call stream_spproj%os_stk%get_ori(1, o_stk)
    !                             call movefile2folder('stk', output_dir_extract, o_stk, err)
    !                             if( .not.err )then
    !                                 call stream_spproj%os_stk%set_ori(1, o_stk)
    !                                 call spproj%append_project(stream_spproj, 'stk')
    !                             endif
    !                         endif
    !                     endif
    !                 endif
    !                 ! add to history
    !                 call o%getter('movie', mov)
    !                 call o%getter('intg', mic)
    !                 call movie_buff%add2history(mov)
    !                 call movie_buff%add2history(mic)
    !                 ! write updated individual project file
    !                 call stream_spproj%write(basename(sp_files(iproj)))
    !                 ! count
    !                 cnt = cnt + 1
    !                 ! cleanup
    !                 call stream_spproj%kill
    !             enddo
    !             call o%kill
    !             call o_stk%kill
    !             write(*,'(A,I3)')'>>> IMPORTED PREVIOUS PROCESSED MOVIES: ', cnt
    !         end subroutine import_prev_streams
    !
    !         subroutine movefile2folder(key, folder, o, err)
    !             use simple_ori, only: ori
    !             character(len=*), intent(in)    :: key, folder
    !             class(ori),       intent(inout) :: o
    !             logical,          intent(out)   :: err
    !             character(len=:), allocatable :: src
    !             character(len=LONGSTRLEN) :: dest,reldest
    !             integer :: iostat
    !             err = .false.
    !             if( .not.o%isthere(key) )then
    !                 err = .true.
    !                 return
    !             endif
    !             call o%getter(key,src)
    !             if( .not.file_exists(src) )then
    !                 err = .true.
    !                 return
    !             endif
    !             dest   = trim(folder)//'/'//basename(src)
    !             iostat = rename(src,dest)
    !             if( iostat /= 0 )then
    !                 THROW_WARN('Ignoring '//trim(src))
    !                 return
    !             endif
    !             call make_relativepath(CWD_GLOB,dest,reldest)
    !             call o%set(key,reldest)
    !         end subroutine movefile2folder
    !
    ! end subroutine exec_preprocess_stream


    ! subroutine exec_pick_extract_stream( self, cline )
    !     class(pick_extract_commander_stream), intent(inout) :: self
    !     class(cmdline),                             intent(inout) :: cline
    !     integer,          parameter   :: WAIT_WATCHER        = 10    ! seconds prior to new stack detection
    !     !integer,          parameter   :: ORIGPROJ_WRITEFREQ  = 600   ! 10mins, Frequency at which the original project file should be updated
    !     type(parameters)                    :: params
    !     type(cmdline)                       :: cline_pick_extract, cline_make_pickrefs
    !     type(sp_project)                    :: orig_proj, stream_proj
    !     type(qsys_env)                      :: qenv
    !     class(cmdline),         allocatable :: completed_jobs_clines(:)
    !     character(LONGSTRLEN),  allocatable :: spproj_list(:)
    !     character(len=:),       allocatable :: spproj_list_fname, output_dir, output_dir_picker, output_dir_extract, projfname
    !     integer :: iter, origproj_time, tnow, iproj, icline, nptcls, prev_stacksz, stacksz
    !     integer :: last_injection, n_spprojs, n_spprojs_prev, n_newspprojs, nmics
    !     if( cline%defined('refs') .and. cline%defined('vol1') )then
    !         THROW_HARD('REFS and VOL1 cannot be both provided!')
    !     endif
    !     if( .not.cline%defined('refs') .and. .not.cline%defined('vol1') )then
    !         THROW_HARD('one of REFS and VOL1 must be provided!')
    !     endif
    !     if( .not. cline%defined('pcontrast') ) call cline%set('pcontrast', 'black')
    !     ! output command line executed
    !     write(logfhandle,'(a)') '>>> COMMAND LINE EXECUTED'
    !     write(logfhandle,*) trim(cmdline_glob)
    !     ! set oritype & defaults
    !     if( .not. cline%defined('oritype') ) call cline%set('oritype', 'mic')
    !     call cline%set('stream','yes')
    !     call cline%set('numlen', 5.)
    !     call params%new(cline)
    !     params_glob%split_mode = 'stream'
    !     params_glob%ncunits    = params%nparts
    !     call cline%set('mkdir', 'no')
    !     if( .not.file_exists(params%projfile) )THROW_HARD('project file: '//trim(params%projfile)//' does not exist!')
    !     spproj_list_fname = filepath(trim(params%dir_target),trim(STREAM_SPPROJFILES))
    !     ! read project info
    !     call orig_proj%read(params%projfile)
    !     ! setup the environment for distributed execution
    !     call qenv%new(1,stream=.true.)
    !     ! output directories
    !     output_dir = PATH_HERE
    !     output_dir_picker  = filepath(trim(output_dir), trim(DIR_PICKER))
    !     output_dir_extract = filepath(trim(output_dir), trim(DIR_EXTRACT))
    !     call simple_mkdir(trim(output_dir),errmsg="commander_stream_wflows :: exec_pick_extract_stream;  ")
    !     call simple_mkdir(trim(output_dir_picker),errmsg="commander_stream_wflows :: exec_pick_extract_stream;  ")
    !     call simple_mkdir(trim(output_dir_extract),errmsg="commander_stream_wflows :: exec_pick_extract_stream;  ")
    !     ! init command-lines
    !     cline_pick_extract  = cline
    !     cline_make_pickrefs = cline
    !     call cline_pick_extract%set('prg', 'pick_extract')
    !     call cline_pick_extract%delete('projname')
    !     ! prepares picking references
    !     call cline_make_pickrefs%set('prg','make_pickrefs')
    !     call cline_make_pickrefs%set('stream','no')
    !     call qenv%exec_simple_prg_in_queue(cline_make_pickrefs, 'MAKE_PICKREFS_FINISHED')
    !     call cline_pick_extract%set('refs', trim(PICKREFS)//params%ext)
    !     call cline_pick_extract%delete('vol1')
    !     ! wait for the first stacks
    !     last_injection  = simple_gettime()
    !     origproj_time   = last_injection
    !     prev_stacksz    = 0
    !     n_spprojs       = 0
    !     n_spprojs_prev  = 0
    !     do iter = 1,999999
    !         tnow = simple_gettime()
    !         if(tnow-last_injection > params%time_inactive .and. stacksz==0)then
    !             write(logfhandle,*)'>>> TIME LIMIT WITHOUT NEW MICROGRAPHS REACHED'
    !             exit
    !         endif
    !         if( file_exists(spproj_list_fname) )then
    !             if( .not.is_file_open(spproj_list_fname) )then
    !                 call read_filetable(spproj_list_fname, spproj_list)
    !                 if( allocated(spproj_list) )n_spprojs = size(spproj_list)
    !                 n_newspprojs = n_spprojs - n_spprojs_prev
    !                 if( n_newspprojs > 0 )then
    !                     ! copy projects and add to processing stack
    !                     do iproj = n_spprojs_prev+1, n_spprojs
    !                         call stream_proj%read(spproj_list(iproj))
    !                         projfname  = filepath(PATH_HERE, basename(spproj_list(iproj)))
    !                         call stream_proj%write(projfname)
    !                         call cline_pick_extract%set('projfile', trim(projfname))
    !                         call qenv%qscripts%add_to_streaming( cline_pick_extract )
    !                         call stream_proj%kill
    !                         deallocate(projfname)
    !                     enddo
    !                     n_spprojs_prev = n_spprojs
    !                 endif
    !                 ! streaming scheduling
    !                 call qenv%qscripts%schedule_streaming( qenv%qdescr )
    !                 stacksz = qenv%qscripts%get_stacksz()
    !                 if( stacksz .ne. prev_stacksz )then
    !                     prev_stacksz = stacksz
    !                     write(logfhandle,'(A,I5)')'>>> MICROGRAPHS TO PROCESS: ', stacksz
    !                 endif
    !                 ! completed jobs update the current project
    !                 if( qenv%qscripts%get_done_stacksz() > 0 )then
    !                     call qenv%qscripts%get_stream_done_stack( completed_jobs_clines )
    !                     do icline=1,size(completed_jobs_clines)
    !                         projfname = completed_jobs_clines(icline)%get_carg('projfile')
    !                         call stream_proj%read( projfname )
    !                         call orig_proj%append_project(stream_proj, 'mic')
    !                         call orig_proj%append_project(stream_proj, 'stk')
    !                         call stream_proj%kill()
    !                         deallocate(projfname)
    !                     enddo
    !                     nptcls = orig_proj%get_nptcls()
    !                     nmics  = orig_proj%get_nintgs()
    !                     write(logfhandle,'(A,I8)')'>>> NEW MICROGRAPHS COUNT: ', nmics
    !                     write(logfhandle,'(A,I8)')'>>> NEW PARTICLES   COUNT: ', nptcls
    !                     call orig_proj%write
    !                     deallocate(completed_jobs_clines)
    !                 endif
    !             endif
    !         endif
    !         call simple_sleep(WAIT_WATCHER)
    !     enddo
    !     ! cleanup
    !     call qsys_cleanup
    !     ! end gracefully
    !     call simple_end('**** SIMPLE_PICK_EXTRACT_STREAM NORMAL STOP ****')
    ! end subroutine exec_pick_extract_stream


end module simple_commander_stream_wflows
