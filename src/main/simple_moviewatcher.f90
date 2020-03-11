! movie watcher for stream processing
module simple_moviewatcher
include 'simple_lib.f08'
use simple_parameters, only: params_glob
implicit none

public :: moviewatcher
private

type moviewatcher
    private
    character(len=LONGSTRLEN), allocatable :: history(:)         !< history of movies detected
    character(len=LONGSTRLEN)          :: cwd            = ''    !< CWD
    character(len=LONGSTRLEN)          :: watch_dir      = ''    !< movies directory to watch
    character(len=STDLEN)              :: ext            = ''    !< target directory
    character(len=STDLEN)              :: fbody          = ''    !< template name
    integer                            :: n_history      = 0     !< history of movies detected
    integer                            :: report_time    = 600   !< time ellapsed prior to processing
    integer                            :: starttime      = 0     !< time of first watch
    integer                            :: ellapsedtime   = 0     !< time ellapsed between last and first watch
    integer                            :: lastreporttime = 0     !< time ellapsed between last and first watch
    integer                            :: n_watch        = 0     !< number of times the folder has been watched
contains
    ! doers
    procedure          :: watch
    procedure, private :: add2history_1
    procedure, private :: add2history_2
    generic            :: add2history => add2history_1, add2history_2
    procedure, private :: is_past
    ! destructor
    procedure          :: kill
end type

interface moviewatcher
    module procedure constructor
end interface moviewatcher

integer, parameter :: FAIL_THRESH = 50
integer, parameter :: FAIL_TIME   = 7200 ! 2 hours

contains

    !>  \brief  is a constructor
    function constructor( report_time )result( self )
        integer, intent(in) :: report_time  ! in seconds
        type(moviewatcher)  :: self
        call self%kill
        if( .not. file_exists(trim(adjustl(params_glob%dir_movies))) )then
            write(logfhandle,*) 'Directory does not exist: ', trim(adjustl(params_glob%dir_movies))
            stop
        endif
        self%cwd         = trim(params_glob%cwd)
        self%watch_dir   = trim(adjustl(params_glob%dir_movies))
        self%report_time = report_time
        self%ext         = trim(adjustl(params_glob%ext))
        self%fbody       = trim(adjustl(params_glob%fbody))
    end function constructor

    !>  \brief  is the watching procedure
    subroutine watch( self, n_movies, movies )
        class(moviewatcher),           intent(inout) :: self
        integer,                       intent(out)   :: n_movies
        character(len=*), allocatable, intent(out)   :: movies(:)
        character(len=:),          allocatable :: list_glob
        character(len=LONGSTRLEN), allocatable :: farray(:)
        integer,                   allocatable :: fileinfo(:)
        logical,                   allocatable :: is_new_movie(:)
        integer                   :: tnow, last_accessed, last_modified, last_status_change ! in seconds
        integer                   :: i, io_stat, n_lsfiles, cnt, fail_cnt
        character(len=LONGSTRLEN) :: fname
        logical                   :: is_closed
        ! init
        self%n_watch = self%n_watch + 1
        tnow = simple_gettime()
        if( self%n_watch .eq. 1 )then
            self%starttime  = tnow ! first call
        endif
        self%ellapsedtime = tnow - self%starttime
        if(allocated(movies))deallocate(movies)
        n_movies = 0
        fail_cnt = 0
        ! builds files array
        list_glob = trim(self%watch_dir)//PATH_SEPARATOR//'*.mrc'
        list_glob = trim(list_glob)//' '//trim(self%watch_dir)//PATH_SEPARATOR//'*.mrcs'
#ifdef USING_TIFF
        list_glob = trim(list_glob)//' '//trim(self%watch_dir)//PATH_SEPARATOR//'*.tif'
        list_glob = trim(list_glob)//' '//trim(self%watch_dir)//PATH_SEPARATOR//'*.tiff'
#endif
        call simple_list_files(list_glob, farray)
        if( .not.allocated(farray) )return ! nothing to report
        n_lsfiles = size(farray)
        ! identifies closed & untouched files
        allocate(is_new_movie(n_lsfiles), source=.false.)
        do i = 1, n_lsfiles
            fname           = trim(adjustl(farray(i)))
            is_new_movie(i) = .not. self%is_past(fname)
            if( .not.is_new_movie(i) )cycle
            call simple_file_stat(fname, io_stat, fileinfo, doprint=.false.)
            is_closed = .not. is_file_open(fname)
            if( io_stat.eq.0 )then
                ! new movie
                last_accessed      = tnow - fileinfo( 9)
                last_modified      = tnow - fileinfo(10)
                last_status_change = tnow - fileinfo(11)
                if(        (last_accessed      > self%report_time)&
                    &.and. (last_modified      > self%report_time)&
                    &.and. (last_status_change > self%report_time)&
                    &.and. is_closed ) is_new_movie(i) = .true.
            else
                ! some error occured
                fail_cnt = fail_cnt + 1
                write(logfhandle,*)'Error watching file: ', trim(fname), ' with code: ',io_stat
            endif
            if(allocated(fileinfo))deallocate(fileinfo)
        enddo
        ! report
        n_movies = count(is_new_movie)
        if( n_movies > 0 )then
            allocate(movies(n_movies))
            cnt = 0
            do i = 1, n_lsfiles
                if( is_new_movie(i) )then
                    cnt   = cnt + 1
                    fname = trim(adjustl(farray(i)))
                    ! call self%add2history( fname )
                    movies(cnt) = trim(fname)
                endif
            enddo
        endif
    end subroutine watch

    !>  \brief  append to history of previously processed movies/micrographs
    subroutine add2history_1(self, list)
        class(moviewatcher),                    intent(inout) :: self
        character(len=LONGSTRLEN), allocatable, intent(in)    :: list(:)
        integer :: i
        if( allocated(list) )then
            do i=1,size(list)
                call self%add2history_2(trim(list(i)))
            enddo
        endif
    end subroutine add2history_1

    !>  \brief  is for adding to the history of already reported files
    !>          absolute path is implied
    subroutine add2history_2( self, fname )
        class(moviewatcher), intent(inout) :: self
        character(len=*),    intent(in)    :: fname
        character(len=LONGSTRLEN), allocatable :: tmp_farr(:)
        character(len=LONGSTRLEN)              :: abs_fname
        integer :: n, tnow
        if( .not.file_exists(fname) )return ! petty triple checking
        if( .not.allocated(self%history) )then
            n = 0
            allocate(self%history(1), stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk("In: simple_moviewatcher%add2history 1", alloc_stat)
        else
            n = size(self%history)
            allocate(tmp_farr(n), stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk("In: simple_moviewatcher%add2history 2", alloc_stat)
            tmp_farr(:) = self%history
            deallocate(self%history)
            allocate(self%history(n+1), stat=alloc_stat)
            if(alloc_stat.ne.0)call allocchk("In: simple_moviewatcher%add2history 3",alloc_stat)
            self%history(:n) = tmp_farr
            deallocate(tmp_farr)
        endif
        abs_fname = simple_abspath(fname)
        self%history(n+1) = trim(adjustl(abs_fname))
        self%n_history    = self%n_history + 1
        write(logfhandle,'(A,A,A,A)')'>>> NEW MOVIE ADDED: ',trim(adjustl(abs_fname)), '; ', cast_time_char(simple_gettime())
    end subroutine add2history_2

    !>  \brief  is for checking a file has already been reported
    !>          absolute path is implied
    logical function is_past( self, fname )
        class(moviewatcher), intent(inout) :: self
        character(len=*),    intent(in)    :: fname
        integer :: i
        is_past = .false.
        if( allocated(self%history) )then
            do i = 1, size(self%history)
                if( trim(adjustl(fname)) .eq. trim(adjustl(self%history(i))) )then
                    is_past = .true.
                    exit
                endif
            enddo
        endif
    end function is_past

    !>  \brief  is a destructor
    subroutine kill( self )
        class(moviewatcher), intent(inout) :: self
        self%cwd        = ''
        self%watch_dir  = ''
        self%ext        = ''
        if(allocated(self%history))deallocate(self%history)
        self%report_time    = 0
        self%starttime      = 0
        self%ellapsedtime   = 0
        self%lastreporttime = 0
        self%n_watch        = 0
        self%n_history      = 0
    end subroutine kill

end module simple_moviewatcher
