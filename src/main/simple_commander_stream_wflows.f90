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
        use simple_image,                  only: image
        use simple_binoris_io,             only: binwrite_oritab, binread_oritab
        use simple_micwatcher        ! use all in there
        use simple_commander_distr   ! use all in there
        use simple_fileio            ! use all in there
        class(prime2D_stream_distr_commander), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        character(len=32), parameter :: ITERFBODY       = 'prime2Ddoc_'
        character(len=32), parameter :: CAVGS_ITERFBODY = 'cavgs_iter'
        character(len=32), parameter :: STK_FILETAB     = 'stktab.txt'
        character(len=32), parameter :: DEFTAB          = 'deftab.txt'
        character(len=32), parameter :: STK_DIR         = './stacks/'
        character(len=32), parameter :: PTCLS_DIR       = 'particles/'
        integer,           parameter :: NCLS_INIT_LIM   = 5    ! FOR TESTING ONLY
        integer,           parameter :: NPTCLS_PER_CLS  = 100
        type(prime2D_distr_commander)      :: xprime2D_distr
        type(micwatcher)                   :: mic_watcher
        character(len=STDLEN), allocatable :: newstacks(:), stktab(:)
        character(len=STDLEN)              :: oritab_glob, str_iter, refs_glob
        type(cmdline)                      :: cline_prime2D, cline_scale
        type(params)                       :: p_master
        real    :: scale, smpd_glob, msk_glob
        integer :: nparts, iter, ncls, n_newmics, ldim(3), nptcls, box_original
        integer :: nptcls_glob, nstacks_glob, ncls_glob, box_glob
        integer :: nptcls_glob_prev, ncls_glob_prev
        ! seed the random number generator
        call seed_rnd
        ! output command line executed
        write(*,'(a)') '>>> COMMAND LINE EXECUTED'
        write(*,*) trim(cmdline_glob)
        ! make target directory
        call exec_cmdline('mkdir -p '//trim(adjustl(STK_DIR))//'|| true')
        ! make master parameters
        p_master = params(cline, checkdistr=.false.)
        ! init command-lines
        cline_scale   = cline
        cline_prime2D = cline
        call cline_prime2D%set('prg',       'prime2D')
        call cline_prime2D%set('autoscale', 'no')
        call cline_prime2D%set('stktab',    STK_FILETAB)
        call cline_prime2D%set('extr_iter', 100.) ! no extremal randomization
        if( p_master%ctf.ne.'no' )then
            call cline_prime2D%set('deftab', DEFTAB)
        endif
        ! init
        ncls_glob        = 0
        ncls_glob_prev   = 0
        nptcls_glob      = 0
        nptcls_glob_prev = 0
        nstacks_glob     = 0
        ! prep scaling parameters
        smpd_glob = LP2SMPDFAC * p_master%lp
        scale     = p_master%smpd / smpd_glob
        ! Instantiate watcher
        mic_watcher = micwatcher(p_master, 30, print=.true.)
        ! Wait for NCLS_INIT_LIM classes
        do 
            call mic_watcher%watch(n_newmics)
            if(n_newmics > 0)then
                if( nptcls_glob .eq. 0 )then
                    ! determines scaling parms
                    call mic_watcher%get_new_stacks(stktab)
                    call find_ldim_nptcls(stktab(1), ldim, nptcls)
                    box_original = ldim(1)
                    box_glob     = find_magic_box(nint(scale*real(ldim(2))))
                    scale        = real(box_glob) / real(box_original)
                    smpd_glob    = p_master%smpd / scale
                    msk_glob     = p_master%msk * scale
                    deallocate(stktab)
                    call cline_scale%set('newbox', real(box_glob))
                    call cline_prime2D%set('box',  real(box_glob))
                    call cline_prime2D%set('smpd', real(smpd_glob))
                    call cline_prime2D%set('msk',  real(msk_glob))
                endif
                nptcls_glob_prev = nptcls_glob
                call add_newstacks
                ncls = nint(real(nptcls_glob) / real(NPTCLS_PER_CLS))
                if(ncls >= NCLS_INIT_LIM) exit
            endif
            call simple_sleep(60) ! parameter instead
        enddo
        ncls_glob = ncls
        call cline_prime2D%set('ncls', real(ncls_glob))
        ! Main loop
        do iter = 1, 999
            str_iter = int2str_pad(iter,3)
            ! prime2D
            call cline_prime2D%set('startit', real(iter))
            call cline_prime2D%set('maxits',  real(iter))
            call xprime2D_distr%execute(cline_prime2D)
            oritab_glob = trim(ITERFBODY)//trim(str_iter)//trim(METADATEXT)
            call rename('prime2Ddoc_final.txt', oritab_glob)
            refs_glob   = trim(CAVGS_ITERFBODY)// trim(str_iter)//trim(p_master%ext)
            call cline_prime2D%set('oritab', trim(oritab_glob))
            call cline_prime2D%set('refs',   trim(refs_glob))
            call remap_empty_classes
            ! detects new images
            call mic_watcher%watch(n_newmics)
            if(n_newmics > 0)then
                nptcls_glob_prev = nptcls_glob
                call add_newstacks
                write(*,'(A,I6)')'>>> NEW PARTICLES COUNT: ', nptcls_glob 
                ! updates particles classes & references
                ncls_glob_prev = ncls_glob
                ncls_glob      = nint(real(nptcls_glob) / real(NPTCLS_PER_CLS))
                call remap_new_classes
                call cline_prime2D%set('ncls', real(ncls_glob))
                write(*,'(A,I6)')'>>> NEW CLASSES COUNT: ', ncls_glob 
            endif
        enddo
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_PRIME2D_STREAM NORMAL STOP ****')

        contains

            subroutine add_newstacks
                use simple_commander_imgproc, only: scale_commander
                type(scale_commander)              :: xscale
                type(oris)                         :: deftab_prev, deftab_glob, deftab_here
                character(len=STDLEN), allocatable :: new_deftabs(:), new_stacks(:), tmp(:)
                character(len=STDLEN) :: stk_scaled, fbody, ext
                integer               :: i, j, nl, nptcls, cnt, cnt2, n, ldim(3), iptcl
                call mic_watcher%get_new_stacks(new_stacks)
                ! number of new particles
                nptcls = 0
                do i = 1, n_newmics
                    call find_ldim_nptcls(new_stacks(i), ldim, n)
                    nptcls = nptcls + n
                enddo
                nptcls_glob = nptcls_glob + nptcls
                if( p_master%ctf.ne.'no' )then
                    ! consolidate deftabs
                    call mic_watcher%get_new_deftabs(new_deftabs)
                    call deftab_glob%new(nptcls_glob)
                    if( nptcls_glob_prev > 0 )then
                        ! transfer previous ctf params to new object
                        call deftab_prev%new(nptcls_glob_prev)
                        call binread_oritab(trim(deftab), deftab_prev, [1, nptcls_glob_prev])
                        do iptcl = 1, nptcls_glob_prev
                            call deftab_glob%set_ori(iptcl, deftab_prev%get_ori(iptcl))
                        enddo 
                        call deftab_prev%kill
                    endif
                    ! transfer new ctf params to new object
                    cnt = nptcls_glob_prev
                    do i = 1, n_newmics
                        nl = nlines(new_deftabs(i))
                        call deftab_here%new(nl)
                        call binread_oritab(trim(new_deftabs(i)), deftab_here, [1, nl])
                        cnt2 = 0
                        do j = cnt+1, cnt+nl
                            cnt2 = cnt2 + 1
                            call deftab_glob%set_ori(j, deftab_here%get_ori(cnt2))
                        enddo
                        cnt = cnt + nl
                    enddo
                    call deftab_here%kill
                    call binwrite_oritab(trim(deftab), deftab_glob, [1, nptcls_glob])
                endif
                ! stacks
                if( nstacks_glob == 0 )then
                    allocate(stktab(n_newmics))
                else
                    ! updates stacks array
                    tmp = stktab
                    deallocate(stktab)
                    allocate(stktab(nstacks_glob+n_newmics))
                    stktab(1:nstacks_glob) = tmp(1:nstacks_glob)
                endif
                ! down-scaling and name updates
                cnt = 0
                do i = nstacks_glob+1, nstacks_glob+n_newmics
                    cnt   = cnt + 1
                    ext   = fname2ext(trim(remove_abspath(trim(new_stacks(cnt)))))
                    fbody = get_fbody(trim(remove_abspath(trim(new_stacks(cnt)))), trim(ext))
                    stk_scaled = trim(STK_DIR)//trim(fbody)// trim('_sc') // trim(p_master%ext)
                    stktab(i)  = trim(stk_scaled)
                    call cline_scale%set('stk', trim(new_stacks(cnt)))
                    call cline_scale%set('outstk', trim(stk_scaled))
                    call xscale%execute(cline_scale)
                enddo
                nstacks_glob = nstacks_glob + n_newmics
                ! writes down-scaled stacks list
                call write_filetable(STK_FILETAB, stktab)
            end subroutine add_newstacks

            subroutine remap_empty_classes
                integer, allocatable :: fromtocls(:,:)
                integer              :: icls
                type(oris)           :: os
                type(image)          :: img_cavg
                call os%new(nptcls_glob)
                call binread_oritab(oritab_glob, os, [1, nptcls_glob])
                call os%fill_empty_classes(ncls_glob, fromtocls)
                if( allocated(fromtocls) )then
                    ! updates document & classes
                    call binwrite_oritab(oritab_glob, os, [1, nptcls_glob])
                    call img_cavg%new([box_glob, box_glob,1], smpd_glob)
                    do icls = 1, size(fromtocls, dim=1)
                        call img_cavg%read(trim(refs_glob), fromtocls(icls, 1))
                        call img_cavg%write(trim(refs_glob), fromtocls(icls, 2))
                    enddo
                    call img_cavg%read(trim(refs_glob), ncls_glob)
                    call img_cavg%write(trim(refs_glob), ncls_glob) ! to preserve size
                endif
            end subroutine remap_empty_classes

            subroutine remap_new_classes
                use simple_ran_tabu, only: ran_tabu
                type(ran_tabu)       :: rt
                type(oris)           :: os, os_prev
                type(image)          :: img_cavg
                integer, allocatable :: fromtocls(:,:), cls(:), pops(:)
                integer              :: icls, ncls_prev, n, iptcl, i
                call os%new(nptcls_glob)
                call os_prev%new(nptcls_glob_prev)
                call binread_oritab(oritab_glob, os_prev, [1, nptcls_glob_prev])
                do iptcl = 1, nptcls_glob_prev
                    call os%set_ori(iptcl, os_prev%get_ori(iptcl))
                enddo
                n = nptcls_glob - nptcls_glob_prev
                ! randomize new ptcls to previous references
                allocate(cls(n))
                rt = ran_tabu(n)
                call rt%balanced(ncls_glob_prev, cls)
                do i = 1, n
                    iptcl = nptcls_glob_prev + i
                    call os%set(iptcl, 'class', real(cls(i)))
                    call os%set(iptcl, 'w',     1.)
                    call os%set(iptcl, 'corr',  0.)
                enddo
                deallocate(cls)
                ! updates references
                if( ncls_glob.eq.ncls_glob_prev )then
                    ! nothing to do
                else
                    call os%fill_empty_classes(ncls_glob, fromtocls)
                    if( allocated(fromtocls) )then
                        call img_cavg%new([box_glob, box_glob,1], smpd_glob)
                        do icls = 1, size(fromtocls, dim=1)
                            call img_cavg%read(trim(refs_glob), fromtocls(icls, 1))
                            call img_cavg%write(trim(refs_glob), fromtocls(icls, 2))
                        enddo
                    endif
                endif
                ! updates document
                call binwrite_oritab(oritab_glob, os, [1, nptcls_glob])
            end subroutine remap_new_classes

    end subroutine exec_prime2D_stream_distr

end module simple_commander_stream_wflows
