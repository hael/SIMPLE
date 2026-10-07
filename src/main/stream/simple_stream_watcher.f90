!@descr: movie watcher for stream processing
module simple_stream_watcher
use simple_core_module_api
use simple_progress
implicit none

public :: stream_watcher
public :: workout_directory_structure, sniff_folders_SJ
private
#include "simple_local_flags.inc"

integer,          parameter :: RATE_INTERVAL     = 3600 ! 1 hour
integer,          parameter :: HISTORY_CAPACITY0 = 1024 ! first capacity of the history; doubled when full

type stream_watcher
    private
    type(string),    allocatable :: history(:)             !< basenames of the files reported, in the order added (the first n_history)
    integer,         allocatable :: history_order(:)       !< history indices in lexical order, for is_past's binary search
    type(string),    allocatable :: watch_dirs(:)          !< directories to watch
    type(string)                 :: watch_dir              !< movies directory to watch
    type(string)                 :: regexp                 !< movies extensions
    integer, public, allocatable :: ratehistory(:)
    integer, public              :: n_history      = 0     !< history of movies detected
    integer, public              :: rate           = 0     !< current rate of movie detection
    integer                      :: report_time    = 600   !< a file is reported once untouched this long; negative: as soon as listed
    integer                      :: ratetime       = 0     !< time of last rate checkpoint
    integer                      :: raten          = 0     !< number imported at last rate checkpoint
    integer                      :: n_watch        = 0     !< number of times the folder has been watched
    logical                      :: exists         = .false.
contains
    ! doers
    procedure          :: watch
    procedure, private :: watchdirs
    procedure, private :: add2history_1
    procedure, private :: add2history_2
    generic            :: add2history => add2history_1, add2history_2
    procedure          :: clear_history
    procedure          :: is_past
    procedure, private :: history_pos
    procedure          :: detect_and_add_dirs
    procedure          :: add2watchdirs
    ! destructor
    procedure          :: kill
end type

interface stream_watcher
    module procedure constructor
end interface stream_watcher

contains

    !>  \brief  is a constructor
    function constructor( report_time, dir, spproj, nretries, suffix_filter )result( self )
        integer,                 intent(in) :: report_time  ! in seconds
        class(string),           intent(in) :: dir
        logical,       optional, intent(in) :: spproj
        integer,       optional, intent(in) :: nretries
        class(string), optional, intent(in) :: suffix_filter
        type(stream_watcher) :: self
        integer :: i
        logical :: l_movies
        call self%kill
        l_movies = .true.
        if( present(spproj) ) l_movies = .not.spproj
        self%watch_dir   = dir
        self%report_time = report_time
        if( l_movies )then
            ! watching movies
            if( .not.file_exists(self%watch_dir) )then
                THROW_HARD('Directory does not exist: '//self%watch_dir%to_char())
            else
                write(logfhandle,'(A,A)')'>>> MOVIES DETECTED FROM: ',self%watch_dir%to_char()
            endif
            if( present(suffix_filter) )then
                self%regexp = '\.mrc$|\.mrcs$'
            endif
            if( present(suffix_filter) )then
                self%regexp = '\'//suffix_filter%to_char()//'.mrc$|\'//suffix_filter%to_char()//'.mrcs$|\'&
                &//suffix_filter%to_char()//'.tif$|\'//suffix_filter%to_char()//'.tiff$|\.eer$'
            else
                self%regexp = '\.mrc$|\.mrcs$|\.tif$|\.tiff$|\.eer$'
            endif
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
                THROW_HARD('Directory does not exist: '//self%watch_dir%to_char())
            else
                write(logfhandle,'(A,A)')'>>> PROJECTS DETECTED FROM: ',self%watch_dir%to_char()
            endif
            self%regexp = '\.simple$'
        endif
        allocate(self%ratehistory(1))
        self%ratehistory(1) = 0
        self%exists  = .true.
    end function constructor

    ! DOERS

    !>  \brief  is the watching procedure
    subroutine watch( self, n_movies, movies, max_nmovies )
        class(stream_watcher),       intent(inout) :: self
        integer,                   intent(out)   :: n_movies
        type(string), allocatable, intent(out)   :: movies(:)
        integer, optional,         intent(in)    :: max_nmovies
        type(string), allocatable :: farray(:)
        integer,      allocatable :: fileinfo(:)
        logical,      allocatable :: is_new_movie(:)
        integer                   :: tnow, last_accessed, last_modified, last_status_change ! in seconds
        integer                   :: i, io_stat, n_lsfiles, cnt
        type(string) :: fname
        if( allocated(movies) ) deallocate(movies)
        n_movies = 0
        if( .not.self%exists )return
        ! init
        self%n_watch = self%n_watch + 1
        tnow = simple_gettime()
        if( self%n_watch .eq. 1 )then
            ! first call: the rate counts from here, so the history added before it (a restart's
            ! restored movies) is not taken for movies detected now
            self%ratetime = tnow
            self%raten    = self%n_history
        endif
        ! get file list
        call self%watchdirs(farray)
        if( .not.allocated(farray) )return ! nothing to report
        n_lsfiles = size(farray)
        ! the listing forks a shell; the settle check reads the clock after it, so a file touched
        ! (by any process) while the listing ran cannot have a timestamp past tnow's second
        tnow = simple_gettime()
        ! identifies closed & untouched files
        allocate(is_new_movie(n_lsfiles), source=.false.)
        cnt = 0
        do i = 1, n_lsfiles
            if( present(max_nmovies) )then
                ! maximum required of new movies reached
                if( cnt >= max_nmovies ) exit
            endif
            fname = farray(i)
            if( self%is_past(fname) )cycle
            call simple_file_stat(fname, io_stat, fileinfo)
            if( io_stat.eq.0 )then
                if( self%report_time < 0 )then
                    ! no settle time: taken as soon as listed, whatever its timestamps
                    is_new_movie(i) = .true.
                    cnt = cnt + 1
                else
                    ! not seen before: reported only once untouched for report_time seconds,
                    ! otherwise left for a later watch
                    last_accessed      = tnow - fileinfo( 9)
                    last_modified      = tnow - fileinfo(10)
                    last_status_change = tnow - fileinfo(11)
                    if(        (last_accessed      > self%report_time)&
                        &.and. (last_modified      > self%report_time)&
                        &.and. (last_status_change > self%report_time) )then
                        is_new_movie(i) = .true.
                        cnt = cnt + 1
                    endif
                endif
            else
                ! some error occured
                write(logfhandle,*)'Error watching file: ', fname%to_char(), ' with code: ',io_stat
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
                    movies(cnt) = farray(i)
                endif
            enddo
        endif
        ! update rates
        if(tnow > self%ratetime + RATE_INTERVAL) then
            self%ratehistory = [self%ratehistory, self%rate]
            self%ratetime = self%ratetime + RATE_INTERVAL
            self%raten    = self%n_history
        ! on the first watch, or a second one in the same second, no time has passed since
        ! ratetime: the rate is left as it is rather than divided by zero
        else if( tnow > self%ratetime )then
            self%rate = nint(3600.0 * (self%n_history - self%raten) / (tnow - self%ratetime))
            self%ratehistory(size(self%ratehistory)) = self%rate
        endif
    end subroutine watch

    !>  \brief  append to history of previously processed movies/micrographs
    subroutine add2history_1(self, list)
        class(stream_watcher),       intent(inout) :: self
        type(string), allocatable, intent(in)    :: list(:)
        integer :: i
        if( allocated(list) )then
            do i=1,size(list)
                call self%add2history_2(list(i))
            enddo
        endif
    end subroutine add2history_1

    !>  \brief  is for adding to the history of already reported files
    !>          absolute path is implied; a file already in the history is not added again
    !>  The history grows by doubling and keeps a lexically sorted index (history_order), so
    !>  adding and looking up cost O(log n) comparisons instead of a pass over the history.
    subroutine add2history_2( self, fname )
        class(stream_watcher), intent(inout) :: self
        class(string),         intent(in)    :: fname
        type(string), allocatable :: tmp_farr(:)
        integer,      allocatable :: tmp_order(:)
        type(string) :: name
        integer      :: n, pos
        logical      :: l_found
        if( .not.file_exists(fname) )return ! petty triple checking
        name = basename(fname)
        pos  = self%history_pos(name, l_found)
        if( l_found ) return
        n = self%n_history
        if( .not.allocated(self%history) )then
            allocate(self%history(HISTORY_CAPACITY0), self%history_order(HISTORY_CAPACITY0))
        else if( n == size(self%history) )then
            call move_alloc(self%history, tmp_farr)
            allocate(self%history(2 * n))
            self%history(:n) = tmp_farr
            deallocate(tmp_farr)
            call move_alloc(self%history_order, tmp_order)
            allocate(self%history_order(2 * n))
            self%history_order(:n) = tmp_order
            deallocate(tmp_order)
        endif
        self%history(n+1) = name
        ! insert into the sorted index
        self%history_order(pos+1:n+1) = self%history_order(pos:n)
        self%history_order(pos)       = n + 1
        self%n_history                = n + 1
    end subroutine add2history_2

    !>  \brief  the position of @p name in the sorted history index: where it is (@p l_found) or
    !>          where it would be inserted
    integer function history_pos( self, name, l_found )
        class(stream_watcher), intent(in)  :: self
        class(string),         intent(in)  :: name
        logical,               intent(out) :: l_found
        character(len=:), allocatable :: key, entry
        integer :: lo, hi, mid
        l_found = .false.
        key = name%to_char()
        lo  = 1
        hi  = self%n_history
        do while( lo <= hi )
            mid   = (lo + hi) / 2
            entry = self%history(self%history_order(mid))%to_char()
            if( entry == key )then
                l_found     = .true.
                history_pos = mid
                return
            else if( llt(entry, key) )then
                lo = mid + 1
            else
                hi = mid - 1
            endif
        enddo
        history_pos = lo
    end function history_pos

    !>  \brief  is for clearing the history of imported files
    subroutine clear_history( self )
        class(stream_watcher), intent(inout) :: self
        if( allocated(self%history)       ) deallocate(self%history)
        if( allocated(self%history_order) ) deallocate(self%history_order)
        self%n_history = 0
    end subroutine clear_history

    !>  \brief  is for checking a file has already been reported
    !>          absolute path is implied
    logical function is_past( self, fname )
        class(stream_watcher), intent(in) :: self
        class(string),       intent(in) :: fname
        integer :: pos
        ! need to use basename here since if movies are symbolic links ls -1f dereferences the links
        ! which would cause all movies to be declared as new because of the path mismatch
        pos = self%history_pos(basename(fname), is_past)
    end function is_past

    subroutine detect_and_add_dirs( self, rootdir, SJdirstruct )
        class(stream_watcher), intent(inout) :: self
        type(string),        intent(in)    :: rootdir
        logical,             intent(in)    :: SJdirstruct
        type(string), allocatable :: dir_movies(:)
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
        class(stream_watcher), intent(inout) :: self
        type(string),        intent(in)    :: fname
        type(string), allocatable :: tmp_farr(:)
        type(string)              :: abs_fname
        integer :: i,n
        logical :: new
        if( .not.file_exists(fname) )then
            write(logfhandle,'(A)')'>>> Directory does not exist: '//fname%to_char()
            return
        endif
        abs_fname = simple_abspath(fname)
        if( abs_fname.eq.self%watch_dir ) return ! is already base directory
        new = .true.
        if( .not.allocated(self%watch_dirs) )then
            allocate(self%watch_dirs(1))
            self%watch_dirs(1) = abs_fname
        else
            n = size(self%watch_dirs)
            do i = 1,n
                if( abs_fname.eq.self%watch_dirs(i) )then
                    new = .false.
                    exit
                endif
            enddo
            if( new )then
                call move_alloc(self%watch_dirs, tmp_farr)
                allocate(self%watch_dirs(n+1))
                self%watch_dirs(:n)  = tmp_farr(:)
                self%watch_dirs(n+1) = abs_fname
            endif
        endif
        if( new )then
            write(logfhandle,'(A,A)')'>>> MOVIES DETECTED FROM: ', abs_fname%to_char()
        endif
    end subroutine add2watchdirs

    !>  \brief  is for watching directories
    subroutine watchdirs( self, farray )
        class(stream_watcher),       intent(in)    :: self
        type(string), allocatable, intent(inout) :: farray(:)
        type(string), allocatable :: tmp_farr(:), tmp_farr2(:)
        type(string)              :: dir
        integer :: idir,ndirs,n_newfiles,nfiles,cnt,i
        if( allocated(farray) ) deallocate(farray)
        ndirs = 0
        if( allocated(self%watch_dirs) ) ndirs = size(self%watch_dirs)
        do idir = 0,ndirs
            if( idir == 0 )then
                dir = self%watch_dir
            else
                dir = self%watch_dirs(idir)
            endif
            if(allocated(tmp_farr)) deallocate(tmp_farr)
            call simple_list_files_regexp(dir, self%regexp%to_char(), tmp_farr)
            if( .not.allocated(tmp_farr) ) cycle
            if( size(tmp_farr) == 0 )then
                deallocate(tmp_farr)
                cycle
            endif
            if( allocated(farray) )then
                n_newfiles = size(tmp_farr)
                nfiles     = size(farray)
                allocate(tmp_farr2(nfiles), source=farray)
                call farray(:)%kill
                deallocate(farray)
                allocate(farray(nfiles+n_newfiles))
                do i = 1, nfiles
                    farray(i) = tmp_farr2(i) 
                enddo
                cnt = 0
                do i = nfiles+1, nfiles+n_newfiles
                    cnt = cnt + 1
                    farray(i) = tmp_farr(cnt)
                enddo
                call tmp_farr2(:)%kill
                deallocate(tmp_farr2)
            else
                allocate(farray(size(tmp_farr)), source=tmp_farr)
                do i = 1, size(tmp_farr)
                    farray(i) = tmp_farr(i)
                enddo
            endif
            call tmp_farr(:)%kill
        enddo
    end subroutine watchdirs

    !>  \brief  is a destructor
    subroutine kill( self )
        class(stream_watcher), intent(inout) :: self
        self%watch_dir = ''
        self%regexp    = ''
        if( allocated(self%history)    ) deallocate(self%history)
        if( allocated(self%history_order)) deallocate(self%history_order)
        if( allocated(self%ratehistory)) deallocate(self%ratehistory)
        if( allocated(self%watch_dirs) ) deallocate(self%watch_dirs)
        self%rate           = 0
        self%report_time    = 0
        self%ratetime       = 0
        self%raten          = 0
        self%n_watch        = 0
        self%n_history      = 0
        self%exists = .false.
    end subroutine kill

    ! PUBLIC UTILITIES

    ! List all directories following the so-called SJ format: directory/xxx/Data
    subroutine sniff_folders_SJ( directory, SJdirstruct, found_directories )
        type(string),              intent(in)    :: directory
        logical,                   intent(inout) :: SJdirstruct
        type(string), allocatable, intent(inout) :: found_directories(:)
        type(string)              :: dir, absdirectory, subdir
        type(string), allocatable :: dirs(:), subdirs(:)
        integer :: i, j, nfound
        SJdirstruct  = .false.
        nfound       = 0
        absdirectory = simple_abspath(directory)
        if( allocated(found_directories) ) deallocate(found_directories)
        ! subdirectories, depth=1
        dirs = simple_list_dirs(absdirectory)
        if( .not.allocated(dirs) ) return
        ! subdirectories, depth=2
        do i = 1,size(dirs)
            dir     = absdirectory//'/'//dirs(i)
            subdirs = simple_list_dirs(dir)
            if( .not.allocated(subdirs) ) cycle
            do j = 1,size(subdirs)
                subdir = dir//'/'//subdirs(j)
                if( subdirs(j) == 'Data' )then
                    nfound = nfound + 1
                    if( .not.allocated(found_directories) )then
                        allocate(found_directories(1), source=[subdir])
                    else
                        found_directories = [found_directories(1:nfound-1), subdir]
                    endif
                endif
            enddo
            deallocate(subdirs)
        enddo
        deallocate(dirs)
        SJdirstruct = nfound > 0
    end subroutine sniff_folders_SJ

    ! Determines which directory structure is present:
    ! 1. movies appear in the folder  'directory'
    ! 2. movies appear in the folders 'directory/xx/Data/'
    subroutine workout_directory_structure( directory, found, SJdirstruct )
        class(string), intent(in)  :: directory
        logical,       intent(out) :: found, SJdirstruct
        type(string),  allocatable :: dirs(:), subdirs(:), tmp(:)
        type(string) :: regexp, absdirectory, dir
        integer :: i, j, nfound
        found        = .false.
        SJdirstruct  = .false.
        absdirectory = simple_abspath(directory)
        ! Test root folder first
        regexp = '\.mrc$|\.mrcs$|\.tif$|\.tiff$|\.eer$'
        call simple_list_files_regexp(absdirectory, regexp%to_char(), tmp)
        if( allocated(tmp) .and. size(tmp) > 0 ) then
            ! Movies are present in the folder: this single folder will be watched
            deallocate(tmp)
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
            dir     = absdirectory//'/'//dirs(i)
            subdirs = simple_list_dirs(dir)
            if( .not.allocated(subdirs) ) cycle
            do j = 1,size(subdirs)
                if( subdirs(j) == 'Data' ) nfound = nfound + 1
            enddo
            deallocate(subdirs)
        enddo
        if( nfound > 0 )then
            ! folder(s) two levels down with suffix Data have been found
            found       = .true.
            SJdirstruct = .true.
        endif
    end subroutine workout_directory_structure

end module simple_stream_watcher
