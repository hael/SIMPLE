! movie watcher for stream processing
module simple_moviewatcher
include 'simple_lib.f08'
use simple_parameters, only: params_glob
implicit none

public :: moviewatcher
private

character(len=STDLEN), parameter :: stream_dirs = 'SIMPLE_STREAM_DIRS'

type moviewatcher
    private
    character(len=LONGSTRLEN), allocatable :: history(:)         !< history of movies detected
    character(len=LONGSTRLEN), allocatable :: watch_dirs(:)      !< directories to watch
    character(len=LONGSTRLEN)          :: cwd            = ''    !< CWD
    character(len=LONGSTRLEN)          :: watch_dir      = ''    !< movies directory to watch
    character(len=STDLEN)              :: ext            = ''    !< target directory
    character(len=STDLEN)              :: fbody          = ''    !< template name
    character(len=STDLEN)              :: regexp         = ''    !< movies extensions
    integer                            :: n_history      = 0     !< history of movies detected
    integer                            :: report_time    = 600   !< time ellapsed prior to processing
    integer                            :: starttime      = 0     !< time of first watch
    integer                            :: ellapsedtime   = 0     !< time ellapsed between last and first watch
    integer                            :: lastreporttime = 0     !< time ellapsed between last and first watch
    integer                            :: n_watch        = 0     !< number of times the folder has been watched
contains
    ! doers
    procedure          :: watch
    procedure, private :: watchdirs
    procedure, private :: add2history_1
    procedure, private :: add2history_2
    generic            :: add2history => add2history_1, add2history_2
    procedure          :: is_past
    procedure, private :: add2watchdirs
    procedure, private :: check4dirs
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
        self%regexp = '\.mrc$|\.mrcs$'
#ifdef USING_TIFF
        self%regexp = '\.mrc$|\.mrcs$|\.tif$|\.tiff$|\.eer$'
#endif
    end function constructor

    !>  \brief  is the watching procedure
    subroutine watch( self, n_movies, movies )
        class(moviewatcher),           intent(inout) :: self
        integer,                       intent(out)   :: n_movies
        character(len=*), allocatable, intent(out)   :: movies(:)
        character(len=LONGSTRLEN), allocatable :: farray(:)
        integer,                   allocatable :: fileinfo(:)
        logical,                   allocatable :: is_new_movie(:)
        integer                   :: tnow, last_accessed, last_modified, last_status_change ! in seconds
        integer                   :: i, io_stat, n_lsfiles, cnt, fail_cnt
        character(len=LONGSTRLEN) :: fname
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
        ! get file list
        call self%check4dirs
        call self%watchdirs(farray)
        if( .not.allocated(farray) )return ! nothing to report
        n_lsfiles = size(farray)
        ! identifies closed & untouched files
        allocate(is_new_movie(n_lsfiles), source=.false.)
        do i = 1, n_lsfiles
            fname           = trim(adjustl(farray(i)))
            is_new_movie(i) = .not. self%is_past(fname)
            if( .not.is_new_movie(i) )cycle
            call simple_file_stat(fname, io_stat, fileinfo, doprint=.false.)
            if( io_stat.eq.0 )then
                ! new movie
                last_accessed      = tnow - fileinfo( 9)
                last_modified      = tnow - fileinfo(10)
                last_status_change = tnow - fileinfo(11)
                if(        (last_accessed      > self%report_time)&
                    &.and. (last_modified      > self%report_time)&
                    &.and. (last_status_change > self%report_time) ) is_new_movie(i) = .true.
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
        integer :: n
        if( .not.file_exists(fname) )return ! petty triple checking
        if( .not.allocated(self%history) )then
            n = 0
            allocate(self%history(1))
        else
            n = size(self%history)
            call move_alloc(self%history, tmp_farr)
            allocate(self%history(n+1))
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

    subroutine check4dirs( self )
        class(moviewatcher), intent(inout) :: self
        character(len=LONGSTRLEN), allocatable :: farr(:)
        integer :: i,n
        if( .not.file_exists(stream_dirs))return
        n = nlines(stream_dirs)
        if( n == 0 ) return
        call read_filetable(stream_dirs, farr)
        n = size(farr)
        do i=1,n
            if( trim(self%watch_dir) .ne. trim(farr(i)) ) call self%add2watchdirs(farr(i))
        enddo
    end subroutine check4dirs

    !>  \brief  is for adding a directory to watch
    subroutine add2watchdirs( self, fname )
        class(moviewatcher), intent(inout) :: self
        character(len=*),    intent(in)    :: fname
        character(len=LONGSTRLEN), allocatable :: tmp_farr(:)
        character(len=LONGSTRLEN)              :: abs_fname
        integer :: i,n
        logical :: new
        if( .not.file_exists(fname) )then
            write(logfhandle,'(A)')'>>> Directory does not exist: '//trim(fname)
            return
        endif
        new = .true.
        abs_fname = simple_abspath(fname)
        if( .not.allocated(self%watch_dirs) )then
            allocate(self%watch_dirs(1))
            self%watch_dirs(1) = trim(adjustl(abs_fname))
        else
            n = size(self%watch_dirs)
            do i = 1,n
                if( trim(abs_fname).eq.trim(self%watch_dirs(i)) )then
                    new = .false.
                    exit
                endif
            enddo
            if( new )then
                allocate(tmp_farr(n), source=self%watch_dirs)
                deallocate(self%watch_dirs)
                allocate(self%watch_dirs(n+1))
                self%watch_dirs(:n) = tmp_farr(:)
                self%watch_dirs(n+1) = trim(adjustl(abs_fname))
            endif
        endif
        if( new )then
            write(logfhandle,'(A,A,A,A)')'>>> NEW DIRECTORY ADDED: ',trim(adjustl(abs_fname)), '; ', cast_time_char(simple_gettime())
        endif
    end subroutine add2watchdirs

    !>  \brief  is for watching directories
    subroutine watchdirs( self, farray )
        class(moviewatcher),                    intent(inout) :: self
        character(len=LONGSTRLEN), allocatable, intent(inout) :: farray(:)
        character(len=LONGSTRLEN), allocatable :: tmp_farr(:), tmp_farr2(:)
        character(len=LONGSTRLEN)              :: dir
        integer :: idir,ndirs,n_newfiles,nfiles
        if( allocated(farray) ) deallocate(farray)
        ndirs = 0
        if( allocated(self%watch_dirs) ) ndirs = size(self%watch_dirs)
        do idir = 0,ndirs
            if( idir == 0 )then
                dir = trim(self%watch_dir)
            else
                dir = trim(self%watch_dirs(idir))
            endif
            if(allocated(tmp_farr)) deallocate(tmp_farr)
            call simple_list_files_regexp(dir, self%regexp, tmp_farr)
            if( .not.allocated(tmp_farr) ) cycle
            if( allocated(farray) )then
                n_newfiles = size(tmp_farr)
                nfiles     = size(farray)
                tmp_farr2  = farray(:)
                deallocate(farray)
                allocate(farray(nfiles+n_newfiles))
                farray(1:nfiles) = tmp_farr2(:)
                farray(nfiles+1:nfiles+n_newfiles) = tmp_farr(:)
            else
                farray = tmp_farr(:)
            endif
        enddo
    end subroutine watchdirs

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
