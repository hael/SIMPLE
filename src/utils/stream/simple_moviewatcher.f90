! movie watcher for stream processing
module simple_moviewatcher
include 'simple_lib.f08'
use simple_progress
implicit none

public :: moviewatcher
public :: workout_directory_structure, sniff_folders_SJ
private
#include "simple_local_flags.inc"

character(len=STDLEN), parameter :: WATCHER_HISTORY = 'watcher_history.txt'
character(len=STDLEN), parameter :: WATCHER_DIRS    = 'watcher_dirs.txt'
integer,               parameter :: RATE_INTERVAL   = 3600 ! 1 hour

type moviewatcher
    private
    character(len=LONGSTRLEN), allocatable :: history(:)         !< history of movies detected
    character(len=LONGSTRLEN), allocatable :: watch_dirs(:)      !< directories to watch
    character(len=LONGSTRLEN)          :: watch_dir      = ''    !< movies directory to watch
    character(len=STDLEN)              :: regexp         = ''    !< movies extensions
    integer, public,       allocatable :: ratehistory(:)
    integer, public                    :: n_history      = 0     !< history of movies detected
    integer, public                    :: rate           = 0     !< current rate of movie detection
    integer                            :: report_time    = 600   !< time ellapsed prior to processing
    integer                            :: starttime      = 0     !< time of first watch
    integer                            :: ellapsedtime   = 0     !< time ellapsed between last and first watch
    integer                            :: lastreporttime = 0     !< time ellapsed between last and first watch
    integer                            :: ratetime       = 0     !< time of last rate checkpoint
    integer                            :: raten          = 0     !< number imported at last rate checkpoint
    integer                            :: n_watch        = 0     !< number of times the folder has been watched
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
    procedure          :: detect_and_add_dirs
    procedure          :: add2watchdirs
    ! destructor
    procedure          :: kill
end type

interface moviewatcher
    module procedure constructor
end interface moviewatcher

contains

    !>  \brief  is a constructor
    function constructor( report_time, dir, spproj, nretries, suffix_filter )result( self )
        integer,                    intent(in) :: report_time  ! in seconds
        character(len=*),           intent(in) :: dir
        logical,          optional, intent(in) :: spproj
        integer,          optional, intent(in) :: nretries
        character(len=*), optional, intent(in) :: suffix_filter
        type(moviewatcher)                     :: self
        integer :: i
        logical :: l_movies
        call self%kill
        l_movies = .true.
        if( present(spproj) ) l_movies = .not.spproj
        self%watch_dir   = trim(adjustl(dir))
        self%report_time = report_time
        if( l_movies )then
            ! watching movies
            if( .not.file_exists(self%watch_dir) )then
                THROW_HARD('Directory does not exist: '//trim(self%watch_dir))
            else
                write(logfhandle,'(A,A)')'>>> MOVIES DETECTED FROM: ',trim(self%watch_dir)
            endif
            if( present(suffix_filter) )then
                self%regexp = '\.mrc$|\.mrcs$'
            endif
#ifdef USING_TIFF
            if( present(suffix_filter) )then
                self%regexp = '\'//trim(suffix_filter)//'.mrc$|\'//trim(suffix_filter)//'.mrcs$|\'&
                &//trim(suffix_filter)//'.tif$|\'//trim(suffix_filter)//'.tiff$|\.eer$'
            else
                self%regexp = '\.mrc$|\.mrcs$|\.tif$|\.tiff$|\.eer$'
            endif
#else
            if( present(suffix_filter) )then
                self%regexp = '\'//trim(suffix_filter)//'.mrc$|\'//trim(suffix_filter)//'.mrcs$'
            else
                self%regexp = '\.mrc$|\.mrcs$'
            endif
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
                write(logfhandle,'(A,A)')'>>> PROJECTS DETECTED FROM: ',trim(self%watch_dir)
            endif
            self%regexp = '\.simple$'
        endif
        allocate(self%ratehistory(1))
        self%ratehistory(1) = 0
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
            self%starttime = tnow ! first call
            self%ratetime  = tnow ! first call
        endif
        self%ellapsedtime = tnow - self%starttime
        fail_cnt = 0
        ! get file list
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
        ! update rates
        if(tnow > self%ratetime + RATE_INTERVAL) then
            self%ratehistory = [self%ratehistory, self%rate]
            self%ratetime = self%ratetime + RATE_INTERVAL
            self%raten    = self%n_history
        else
            self%rate = nint(3600.0 * (self%n_history - self%raten) / (tnow - self%ratetime))
            self%ratehistory(size(self%ratehistory)) = self%rate
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
        self%n_history = 0
    end subroutine clear_history

    !>  \brief  is for checking a file has already been reported
    !>          absolute path is implied
    logical function is_past( self, fname )
        class(moviewatcher), intent(in) :: self
        character(len=*),    intent(in) :: fname
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

    subroutine detect_and_add_dirs( self, rootdir, SJdirstruct )
        class(moviewatcher),     intent(inout) :: self
        character(len=*),        intent(in)    :: rootdir
        logical,                 intent(in)    :: SJdirstruct
        character(len=LONGSTRLEN), allocatable :: dir_movies(:)
        integer :: i
        logical :: l_new_dir
        if( SJdirstruct )then
            ! Runtime detection of new grid square directories
            call sniff_folders_SJ( rootdir, l_new_dir, dir_movies )
            if( l_new_dir )then
                do i = 1,size(dir_movies)
                    call self%add2watchdirs(dir_movies(i))
                enddo
                deallocate(dir_movies)
            endif
        else
            ! not relevant yet
        endif
    end subroutine detect_and_add_dirs

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
                self%watch_dirs(:n)  = tmp_farr(:)
                self%watch_dirs(n+1) = trim(adjustl(abs_fname))
            endif
        endif
        if( new )then
            write(logfhandle,'(A,A)')'>>> MOVIES DETECTED FROM: ',trim(adjustl(abs_fname))
        endif
    end subroutine add2watchdirs

    !>  \brief  is for watching directories
    subroutine watchdirs( self, farray, chrono )
        class(moviewatcher),                    intent(in)    :: self
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
        self%watch_dir = ''
        self%regexp    = ''
        if( allocated(self%history)    ) deallocate(self%history)
        if( allocated(self%ratehistory)) deallocate(self%ratehistory)
        if( allocated(self%watch_dirs) ) deallocate(self%watch_dirs)
        self%rate           = 0
        self%report_time    = 0
        self%starttime      = 0
        self%ellapsedtime   = 0
        self%lastreporttime = 0
        self%ratetime       = 0
        self%raten          = 0
        self%n_watch        = 0
        self%n_history      = 0
        self%exists = .false.
    end subroutine kill

    ! PUBLIC UTILITIES

    ! List all directories following the so-called SJ format: directory/xxx/Data
    subroutine sniff_folders_SJ( directory, SJdirstruct, found_directories )
        character(len=*),                       intent(in)    :: directory
        logical,                                intent(inout) :: SJdirstruct
        character(len=LONGSTRLEN), allocatable, intent(inout) :: found_directories(:)
        character(len=:),          allocatable :: dir, absdirectory, subdir
        character(len=LONGSTRLEN), allocatable :: dirs(:), subdirs(:), tmp(:)
        integer :: i, j, nfound
        SJdirstruct  = .false.
        nfound       = 0
        absdirectory = simple_abspath(directory)
        if( allocated(found_directories) ) deallocate(found_directories)
        ! subdirectories, depth=1
        dirs = simple_list_dirs(absdirectory)
        if( .not.allocated(dirs) ) return
        ! subdirectories, depth=2
        allocate(found_directories(0))
        do i = 1,size(dirs)
            dir     = absdirectory//'/'//trim(dirs(i))
            subdirs = simple_list_dirs(dir)
            if( .not.allocated(subdirs) ) cycle
            do j = 1,size(subdirs)
                subdir = dir//'/'//trim(subdirs(j))
                if( trim(subdirs(j)) == 'Data' )then
                    nfound = nfound + 1
                    call move_alloc(found_directories, tmp)
                    allocate(found_directories(nfound))
                    found_directories(1:nfound-1) = tmp(:)
                    found_directories(nfound)     = subdir
                    deallocate(tmp)
                endif
            enddo
            deallocate(subdirs)
        enddo
        SJdirstruct = nfound > 0
    end subroutine sniff_folders_SJ

    ! Determines which directory structure is present:
    ! 1. movies appear in the folder  'directory'
    ! 2. movies appear in the folders 'directory/xx/Data/'
    subroutine workout_directory_structure( directory, found, SJdirstruct )
        character(len=*),          intent(in)  :: directory
        logical,                   intent(out) :: found, SJdirstruct
        character(len=LONGSTRLEN), allocatable :: dirs(:), subdirs(:), tmp(:)
        character(len=:),          allocatable :: regexp, absdirectory, dir
        integer :: i, j, nfound
        found        = .false.
        SJdirstruct  = .false.
        absdirectory = simple_abspath(directory)
        ! Test root folder first
        regexp = '\.mrc$|\.mrcs$|\.tif$|\.tiff$|\.eer$'
        call simple_list_files_regexp(absdirectory, regexp, tmp)
        if( allocated(tmp) )then
            ! Movies are present in the folder: this single folder will be watched
            found = .true.
            return
        endif
        ! Test for subdirectories, depth=2
        nfound = 0
        dirs   = simple_list_dirs(absdirectory)
        if( .not.allocated(dirs) )then
            ! no subfolder, better luck next time
            return
        endif
        do i = 1,size(dirs)
            dir     = absdirectory//'/'//trim(dirs(i))
            subdirs = simple_list_dirs(dir)
            if( .not.allocated(subdirs) ) cycle
            do j = 1,size(subdirs)
                if( trim(subdirs(j)) == 'Data' ) nfound = nfound + 1
            enddo
            deallocate(subdirs)
        enddo
        if( nfound > 0 )then
            ! folder(s) two levels down with suffix Data have been found
            found       = .true.
            SJdirstruct = .true.
        endif
    end subroutine workout_directory_structure

end module simple_moviewatcher
