! concrete commander: stream processing routines
module simple_commander_stream_wflows
#include "simple_lib.f08"

use simple_cmdline,           only: cmdline
use simple_chash,             only: chash
use simple_params,            only: params
use simple_commander_base,    only: commander_base
use simple_commander_preproc, only: preproc_commander
use simple_qsys_env,          only: qsys_env
use simple_jiffys,            only: simple_end
use simple_qsys_funs          ! use all in there
implicit none

public :: preproc_stream_commander
public :: prime2D_stream_distr_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: preproc_stream_commander
  contains
    procedure :: execute      => exec_preproc_stream
end type preproc_stream_commander
type, extends(commander_base) :: prime2D_stream_distr_commander
  contains
    procedure :: execute      => exec_prime2D_stream_distr
end type prime2D_stream_distr_commander

contains

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
        debug=.true.  ! from local flags, for debugging
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

    subroutine exec_prime2D_stream_distr( self, cline )
        use simple_commander_distr_wflows, only: prime2D_distr_commander
        use simple_oris,                   only: oris
        use simple_micwatcher        ! use all in there
        use simple_commander_distr   ! use all in there
        use simple_fileio            ! use all in there
        class(prime2D_stream_distr_commander), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        type(prime2D_distr_commander)      :: xprime2D_distr
        type(micwatcher)                   :: mic_watcher
        character(len=STDLEN), allocatable :: newstacks(:), stacktab(:)
        character(len=STDLEN)              :: oritab_glob, str_iter, refs_glob
        type(cmdline)                      :: cline_prime2D
        type(qsys_env)                     :: qenv
        type(params)                       :: p_master
        type(chash)                        :: job_descr
        type(oris)                         :: os
        integer :: ipart, nl, ishift, nparts, iter, ncls, n_newmics
        integer :: nptcls_glob, nstacks_glob, ncls_glob, box_glob, smpd_glob, msk_glob
        character(len=STDLEN), parameter :: DEFTAB         = 'deftab.txt'
        integer,               parameter :: NCLS_INIT_LIM  = 10
        integer,               parameter :: NPTCLS_PER_CLS = 400
        integer,               parameter :: STARTIT        = 1
        character(len=32),     parameter :: ITERFBODY      = 'prime2Ddoc_'
        character(len=STDLEN), parameter :: CAVGNAMES       = 'cavgs_final'//METADATEXT
        character(len=32),     parameter :: CAVGS_ITERFBODY = 'cavgs_iter'
        ! seed the random number generator
        call seed_rnd
        ! output command line executed
        write(*,'(a)') '>>> COMMAND LINE EXECUTED'
        write(*,*) trim(cmdline_glob)
        ! make master parameters
        p_master = params(cline, checkdistr=.false.)
        ! init command-lines
        cline_prime2D = cline
        !call cline_prime2D%set('stream', 'yes')
        call cline_prime2D%set('autoscale', 'no')
        call cline_prime2D%set('startit',   real(STARTIT))
        call cline_prime2D%set('maxits',    real(STARTIT))
        call cline_prime2D%set('deftab',    DEFTAB)
        call cline_prime2D%set('ctf',       'yes')
        call cline_prime2D%set('extr_ietr', 100.) ! guarantees no extremal randomization
        ! init
        p_master%numlen = 5
        ncls_glob    = 0
        nptcls_glob  = 0
        nstacks_glob = 0
        oritab_glob  = 'prime2Ddoc_curr.txt'
        ! Instantiate watcher
        mic_watcher = micwatcher(p_master, 30, print=.true.)
        ! Wait for NCLS_INIT_LIM classes
        do 
            call mic_watcher%watch(n_newmics)
            if(n_newmics > 0)then
                call add_newstacks
                ncls = nint(real(nptcls_glob) / real(NPTCLS_PER_CLS))
                if(ncls > NCLS_INIT_LIM) exit
            endif
            call simple_sleep(60)
        enddo
        ncls_glob = ncls
        call cline_prime2D%set('ncls', real(ncls_glob))
        ! setup the environment for distributed execution
        call qenv%new(p_master)
        ! prepare job description
        call cline_prime2D%gen_job_descr(job_descr)
        ! Main loop
        do iter = 1, 999
            str_iter = int2str_pad(iter,3)
            ! prime2D
            call xprime2D_distr%execute(cline_prime2D)
            oritab_glob = trim(ITERFBODY)//trim(str_iter)//trim(METADATEXT)
            refs_glob   = trim(CAVGS_ITERFBODY)// trim(str_iter) //trim(p_master%ext)
            call job_descr%set('oritab',  trim(oritab_glob))
            call job_descr%set('refs', trim(refs_glob))
            call job_descr%set('startit', int2str(iter))
            ! detects new images
            call mic_watcher%watch(n_newmics)
            if(n_newmics > 0)then
                call add_newstacks
                ncls_glob = nint(real(nptcls_glob) / real(NPTCLS_PER_CLS))
                call remap_classes
                call job_descr%set('ncls', int2str(ncls_glob))
            endif
        enddo

        call simple_end('**** SIMPLE_DISTR_PRIME2D_STREAM NORMAL STOP ****')

        contains

            subroutine add_newstacks
                use simple_commander_imgproc, only: scale_commander
                use simple_binoris_io,        only: binwrite_oritab, binread_oritab
                type(scale_commander)              :: xscale
                type(oris)                         :: deftab_glob
                character(len=STDLEN), allocatable :: new_deftabs(:), new_stacks(:), tmp(:)
                integer                            :: i, nl, nptcls, cnt
                call mic_watcher%get_new_stacks(new_stacks)
                call mic_watcher%get_new_deftabs(new_deftabs)
                ! number of particles
                nptcls = 0
                do i = 1, n_newmics
                    nptcls = nptcls + nlines(new_deftabs(i))
                enddo
                ! consolidate deftabs
                call deftab_glob%new(nptcls_glob+nptcls)
                if( nptcls_glob > 0 )then
                    call binread_oritab(trim(deftab), deftab_glob, [1, nptcls_glob])
                endif
                cnt = nptcls_glob
                do i = 1, n_newmics
                    nl = nlines(new_deftabs(i))
                    call binread_oritab(trim(new_deftabs(i)), deftab_glob, [cnt+1, cnt+nl])
                    cnt = cnt + nl
                enddo
                nptcls_glob = nptcls_glob + nptcls
                call binwrite_oritab(trim(deftab), deftab_glob, [1, nptcls_glob])
                ! update & write stack list
                if( nstacks_glob == 0 )then
                    allocate(stacktab(n_newmics))
                    ! BOX SIZE AND SMPD GLOB NEED TO BE DETERMINED HERE
                else
                    ! updates stacks array
                    tmp = stacktab
                    deallocate(stacktab)
                    allocate(stacktab(nstacks_glob+n_newmics))
                    stacktab(1:nstacks_glob) = tmp(1:nstacks_glob)
                endif
                cnt = 0
                do i = nstacks_glob, nstacks_glob + n_newmics
                    cnt = cnt + 1
                    ! down scale stacks here and update names !!!!!!!!1
                    !!!!!1
                    stacktab(i) = new_stacks(cnt)
                enddo
                nstacks_glob = nstacks_glob + n_newmics
                ! writes downscale stack list here
            end subroutine add_newstacks

            subroutine remap_classes
                ! remap cavgs
                ! new particles randomization & oritab_glob update
            end subroutine remap_classes

    end subroutine exec_prime2D_stream_distr

end module simple_commander_stream_wflows
