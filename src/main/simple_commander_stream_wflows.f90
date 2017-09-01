! concrete commander: stream processing routines
module simple_commander_stream_wflows
use simple_defs
use simple_syslib,            only: alloc_errchk,exec_cmdline, simple_sleep 
use simple_cmdline,           only: cmdline
use simple_chash,             only: chash
use simple_params,            only: params
use simple_commander_base,    only: commander_base
use simple_strings,           only: real2str
use simple_commander_preproc, only: preproc_commander
use simple_qsys_env,          only: qsys_env
use simple_qsys_funs          ! use all in there
implicit none

public :: preproc_stream_commander
private

type, extends(commander_base) :: preproc_stream_commander
  contains
    procedure :: execute      => exec_preproc_stream
end type preproc_stream_commander
#include "simple_local_flags.inc"
contains

    ! PRE-PROCESS SINGLE-PARTICLE DDDs IN STREAMING MODE

    subroutine exec_preproc_stream( self, cline )
        use simple_commander_preproc, only: preproc_commander
        use simple_moviewatcher,      only: moviewatcher
        class(preproc_stream_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        integer,               parameter   :: SHORTTIME = 30   ! folder watched every minute
        integer,               parameter   :: LONGTIME  = 15  ! 15 mins before processing a new movie
        character(len=STDLEN), allocatable :: movies(:)
        character(len=STDLEN)    :: movie
        type(qsys_env)           :: qenv
        type(params)             :: p_master
        type(moviewatcher)       :: movie_buff
        integer                  :: nmovies, imovie, stacksz, prev_stacksz
        integer, parameter       :: TRAILING=5
        debug=.true.  ! from local flags
        ! make master parameters
        p_master = params(cline, checkdistr=.false.)
        ! set defaults
        p_master%split_mode = 'stream'
        p_master%numlen     = 5
        call cline%set('numlen', real(p_master%numlen))
        call cline%set('stream', 'yes')
        if( cline%defined('fbody') )then
            call cline%set('fbody', trim(p_master%dir_target)//'/'//trim(p_master%fbody))
        else
            call cline%set('fbody', trim(p_master%dir_target)//'/')
        endif
        ! make target directory
        call exec_cmdline('mkdir -p '//trim(adjustl(p_master%dir_target))//'|| true')
        ! setup the environment for distributed execution
        call qenv%new(p_master, stream=.true.)
        ! movie watcher init
        movie_buff = moviewatcher(p_master, LONGTIME, print=.true.)
        ! start watching
        prev_stacksz = 0
        nmovies      = 0
        do
            ! watch
            call movie_buff%watch( nmovies, movies )
            ! append movies to processing stack
            if( nmovies > 0 )then
                do imovie = 1, nmovies
                    movie = trim(adjustl(movies(imovie)))
                    call create_individual_filetab
                    call qenv%qscripts%add_to_streaming( cline )
                enddo
            endif
            ! manage stream scheduling
            call qenv%qscripts%schedule_streaming( qenv%qdescr )
            stacksz = qenv%qscripts%get_stacksz()
            if( stacksz .ne. prev_stacksz )then
                prev_stacksz = stacksz
                write(*,'(A,I5)')'>>> MOVIES TO PROCESS: ', stacksz
            endif
            ! wait
            call simple_sleep(SHORTTIME)
        end do

        contains

            subroutine create_individual_filetab
                use simple_fileio, only: fopen, fclose, fileio_errmsg
                integer               :: fnr, file_stat
                character(len=STDLEN) :: fname, ext, movie_here
                movie_here = remove_abspath(trim(movie))
                ext   = fname2ext(trim(movie_here))
                fname = 'preproc_'//trim(get_fbody(trim(movie_here), trim(ext)))//'.txt'
                call fopen(fnr, status='replace', action='write', file=trim(fname), iostat=file_stat)
                call fileio_errmsg('exec_preproc_stream :: create_individual_filetab '//trim(fname), file_stat)
                write(fnr,'(a)') trim(movie)
                call fclose(fnr, errmsg='exec_preproc_stream closing filetab '//trim(fname))
                call cline%set('filetab', fname)
            end subroutine create_individual_filetab

    end subroutine exec_preproc_stream

end module simple_commander_stream_wflows
