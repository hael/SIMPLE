! movie watcher for stream processing
module simple_moviewatcher
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_progress
implicit none

public :: moviewatcher
private
#include "simple_local_flags.inc"

character(len=STDLEN), parameter :: STREAM_DIRS     = 'SIMPLE_STREAM_DIRS'
character(len=STDLEN), parameter :: WATCHER_HISTORY = 'watcher_history.txt'
character(len=STDLEN), parameter :: WATCHER_DIRS    = 'watcher_dirs.txt'

type moviewatcher
    private
    character(len=LONGSTRLEN), allocatable :: history(:)         !< history of movies detected
    character(len=LONGSTRLEN), allocatable :: watch_dirs(:)      !< directories to watch
    character(len=LONGSTRLEN)          :: cwd            = ''    !< CWD
    character(len=LONGSTRLEN)          :: watch_dir      = ''    !< movies directory to watch
    character(len=STDLEN)              :: ext            = ''    !< target directory
    character(len=STDLEN)              :: regexp         = ''    !< movies extensions
    integer, public                    :: n_history      = 0     !< history of movies detected
    integer                            :: report_time    = 600   !< time ellapsed prior to processing
    integer                            :: starttime      = 0     !< time of first watch
    integer                            :: ellapsedtime   = 0     !< time ellapsed between last and first watch
    integer                            :: lastreporttime = 0     !< time ellapsed between last and first watch
    integer                            :: n_watch        = 0     !< number of times the folder has been watched
    integer                            :: len_ext        = 0
    logical                            :: exists         = .false.
contains
    ! getters
    procedure          :: does_exist
    ! I/O
    procedure          :: write_checkpoint
    ! doers
    procedure          :: watch
    procedure, private :: watchdirs
    procedure, private :: add2history_1
    procedure, private :: add2history_2
    generic            :: add2history => add2history_1, add2history_2
    procedure          :: clear_history
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
    function constructor( report_time, dir, spproj, checkpoint_dir, nretries )result( self )
        integer,                    intent(in) :: report_time  ! in seconds
        character(len=*),           intent(in) :: dir
        logical,          optional, intent(in) :: spproj
        character(len=*), optional, intent(in) :: checkpoint_dir
        integer,          optional, intent(in) :: nretries
        type(moviewatcher)                     :: self
        character(len=LONGSTRLEN), allocatable :: strings(:)
        integer :: i
        logical :: l_movies, l_checkpoint
        call self%kill
        l_movies     = .true.
        l_checkpoint = present(checkpoint_dir)
        if( present(spproj) ) l_movies = .not.spproj
        self%watch_dir   = trim(adjustl(dir))
        self%cwd         = trim(params_glob%cwd)
        self%report_time = report_time
        if( l_checkpoint )then
            if( .not.dir_exists(checkpoint_dir) )then
                THROW_HARD('Could not find checkpoint directory: '//trim(checkpoint_dir))
            endif
            if(         (.not.file_exists(trim(checkpoint_dir)//'/'//trim(WATCHER_HISTORY)))&
                &.and.  (.not.file_exists(trim(checkpoint_dir)//'/'//trim(WATCHER_DIRS))))then
                THROW_HARD('Insufficient information for checkpoint inititalization')
            endif
            ! directories to watch
            call read_filetable(trim(checkpoint_dir)//trim(WATCHER_DIRS), strings)
            if( .not.allocated(strings) ) THROW_HARD('No directory defined in moviewatcher checkpoint')
            do i = 1,size(strings)
                call self%add2watchdirs(strings(i))
            enddo
            deallocate(strings)
            ! history of files watched
            call read_filetable(trim(checkpoint_dir)//trim(WATCHER_HISTORY), strings)
            if( allocated(strings) )then
                call self%add2history(strings)
                deallocate(strings)
                write(logfhandle,'(A,I6,A)')'>>> ADDED ',self%n_history,' FILES TO WATCHING HISTORY'
            endif
        endif
        if( l_movies )then
            ! watching movies
            if( .not.file_exists(self%watch_dir) )then
                THROW_HARD('Directory does not exist: '//trim(self%watch_dir))
            else
                write(logfhandle,'(A,A)')'>>> MOVIES WILL BE DETECTED FROM DIRECTORY: ',trim(self%watch_dir)
            endif
            self%ext    = trim(adjustl(params_glob%ext))
            self%regexp = '\.mrc$|\.mrcs$'
#ifdef USING_TIFF
            self%regexp = '\.mrc$|\.mrcs$|\.tif$|\.tiff$|\.eer$'
#endif
        else
            ! watching simple projects
            if(present(nretries)) then
                do i=1, nretries
                    if(file_exists(self%watch_dir)) then
                        exit
                    endif
                    call sleep(10)
                end do
            endif
            if( .not.file_exists(self%watch_dir) )then
                THROW_HARD('Directory does not exist: '//trim(self%watch_dir))
            else
                write(logfhandle,'(A,A)')'>>> PROJECTS WILL BE DETECTED FROM DIRECTORY: ',trim(self%watch_dir)
            endif
            self%ext         = trim(adjustl(METADATA_EXT))
            self%regexp = '\.simple$'
        endif
        self%len_ext = len_trim(self%ext)
        self%exists  = .true.
    end function constructor

    logical function does_exist( self )
        class(moviewatcher), intent(in) :: self
        does_exist = self%exists
    end function does_exist

    ! I/O

    subroutine write_checkpoint( self )
        class(moviewatcher), intent(in) :: self
        if( self%exists )then
            if( self%n_history == 0 ) return
            call write_filetable(WATCHER_HISTORY, self%history)
            if( allocated(self%watch_dirs))then
                call write_filetable(WATCHER_DIRS, [self%watch_dir, self%watch_dirs(:)])
            else
                call write_singlelineoftext(WATCHER_DIRS, self%watch_dir)
            endif
        endif
    end subroutine write_checkpoint

    ! DOERS

    !>  \brief  is the watching procedure
    subroutine watch( self, n_movies, movies, max_nmovies, chrono )
        class(moviewatcher),           intent(inout) :: self
        integer,                       intent(out)   :: n_movies
        character(len=*), allocatable, intent(out)   :: movies(:)
        integer,          optional,    intent(in)    :: max_nmovies
        logical,          optional,    intent(in)    :: chrono
        character(len=LONGSTRLEN), allocatable :: farray(:)
        integer,                   allocatable :: fileinfo(:)
        logical,                   allocatable :: is_new_movie(:)
        integer                   :: tnow, last_accessed, last_modified, last_status_change ! in seconds
        integer                   :: i, io_stat, n_lsfiles, cnt, fail_cnt
        character(len=LONGSTRLEN) :: fname
        if( allocated(movies) ) deallocate(movies)
        n_movies = 0
        if( .not.self%exists )return
        ! init
        self%n_watch = self%n_watch + 1
        tnow = simple_gettime()
        if( self%n_watch .eq. 1 )then
            self%starttime  = tnow ! first call
        endif
        self%ellapsedtime = tnow - self%starttime
        fail_cnt = 0
        ! get file list
        ! call self%check4dirs ! Only necessary for multiple directories!
        call self%watchdirs(farray, chrono)
        if( .not.allocated(farray) )return ! nothing to report
        n_lsfiles = size(farray)
        ! identifies closed & untouched files
        allocate(is_new_movie(n_lsfiles), source=.false.)
        cnt = 0
        do i = 1, n_lsfiles
            if( present(max_nmovies) )then
                ! maximum required of new movies reached
                if( cnt >= max_nmovies ) exit
            endif
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
                    &.and. (last_status_change > self%report_time) )then
                    is_new_movie(i) = .true.
                    cnt = cnt + 1
                endif
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
                    movies(cnt) = trim(adjustl(farray(i)))
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
                call self%add2history_2(list(i))
            enddo
        endif
    end subroutine add2history_1

    !>  \brief  is for adding to the history of already reported files
    !>          absolute path is implied
    subroutine add2history_2( self, fname )
        class(moviewatcher), intent(inout) :: self
        character(len=*),    intent(in)    :: fname
        character(len=LONGSTRLEN), allocatable :: tmp_farr(:)
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
        self%history(n+1) = trim(adjustl(basename_safe(fname)))
        self%n_history    = self%n_history + 1
        ! write(logfhandle,'(A,A,A,A)')'>>> NEW MOVIE ADDED: ',trim(adjustl(abs_fname)), '; ', cast_time_char(simple_gettime())
        !call lastfoundfile_update()
    end subroutine add2history_2

    !>  \brief  is for clearing the history of imported files
    subroutine clear_history( self )
        class(moviewatcher), intent(inout) :: self
        if( allocated(self%history) ) deallocate(self%history)
        self%n_history    = 0
    end subroutine clear_history

    !>  \brief  is for checking a file has already been reported
    !>          absolute path is implied
    logical function is_past( self, fname )
        class(moviewatcher), intent(inout) :: self
        character(len=*),    intent(in)    :: fname
        character(len=LONGSTRLEN) :: fname1
        integer :: i
        is_past = .false.
        if( allocated(self%history) )then
            ! need to use basename here since if movies are symbolic links ls -1f dereferences the links
            ! which would cause all movies to be declared as new because of the path mismatch
            fname1 = adjustl(basename_safe(fname))
            !$omp parallel do private(i) default(shared) proc_bind(close)
            do i = 1, size(self%history)
                if( .not.is_past )then
                    if( trim(fname1) .eq. trim(self%history(i)) )then
                        !$omp critical
                        is_past = .true.
                        !$omp end critical
                    endif
                endif
            enddo
            !$omp end parallel do
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
        abs_fname = simple_abspath(fname)
        if( trim(abs_fname).eq.trim(self%watch_dir) ) return ! is already base directory
        new = .true.
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
                call move_alloc(self%watch_dirs, tmp_farr)
                allocate(self%watch_dirs(n+1))
                self%watch_dirs(:n) = tmp_farr(:)
                self%watch_dirs(n+1) = trim(adjustl(abs_fname))
            endif
        endif
        if( new )then
            write(logfhandle,'(A,A)')'>>> MOVIES DETECTED FROM NEW DIRECTORY: ',trim(adjustl(abs_fname))
        endif
    end subroutine add2watchdirs

    !>  \brief  is for watching directories
    subroutine watchdirs( self, farray, chrono )
        class(moviewatcher),                    intent(inout) :: self
        character(len=LONGSTRLEN), allocatable, intent(inout) :: farray(:)
        logical,                   optional,    intent(in)    :: chrono
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
            call simple_list_files_regexp(dir, self%regexp, tmp_farr, chronological=chrono)
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
        if( allocated(self%history)    ) deallocate(self%history)
        if( allocated(self%watch_dirs) ) deallocate(self%watch_dirs)
        self%report_time    = 0
        self%starttime      = 0
        self%ellapsedtime   = 0
        self%lastreporttime = 0
        self%n_watch        = 0
        self%n_history      = 0
        self%len_ext        = 0
        self%exists = .false.
    end subroutine kill

end module simple_moviewatcher
