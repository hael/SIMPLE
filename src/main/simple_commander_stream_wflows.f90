! concrete commander: stream processing routines
module simple_commander_stream_wflows
use simple_defs
use simple_filehandling       ! use all in there
use simple_syscalls           ! use all in there
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

contains

    ! PRE-PROCESS SINGLE-PARTICLE DDDs IN STREAMING MODE

    subroutine exec_preproc_stream( self, cline )
        use simple_commander_preproc, only: preproc_commander
        use simple_moviewatcher,      only: moviewatcher
        class(preproc_stream_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        logical,               parameter   :: DEBUG = .true.
        integer,               parameter   :: SHORTTIME = 30    ! 60 secs for watching folder
        integer,               parameter   :: LONGTIME  = 10   ! 20 mins before processing a new movie
        character(len=STDLEN), allocatable :: movies(:)
        character(len=STDLEN)    :: movie
        type(qsys_env)           :: qenv
        type(params)             :: p_master
        type(moviewatcher)       :: movie_buff
        integer                  :: nmovies, imovie, movie_ind
        integer, parameter       :: TRAILING=5
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
        call exec_cmdline('mkdir -p '//trim(adjustl(p_master%dir_target)))
        ! setup the environment for distributed execution
        call qenv%new(p_master, stream=.true.)
        ! movie watcher init
        movie_buff = moviewatcher(p_master, LONGTIME, print=.true.)
        ! start watching
        nmovies    = 0
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
            ! wait
            call simple_sleep(SHORTTIME)
        end do

        contains

            subroutine create_individual_filetab
                integer               :: fnr, file_stat
                character(len=STDLEN) :: fname, ext, movie_here
                movie_here = remove_abspath(trim(movie))
                ext   = fname2ext(trim(movie_here))
                fname = 'preproc_'//trim(get_fbody(trim(movie_here), trim(ext)))//'.txt'
                fnr   = get_fileunit()
                open(unit=fnr, status='replace', action='write', file=trim(fname), iostat=file_stat)
                call fopen_err('exec_preproc_stream :: create_individual_filetab', file_stat)
                write(fnr,'(a)') trim(movie)
                close(unit=fnr)
                call cline%set('filetab', fname)
            end subroutine create_individual_filetab

    end subroutine exec_preproc_stream

end module simple_commander_stream_wflows
