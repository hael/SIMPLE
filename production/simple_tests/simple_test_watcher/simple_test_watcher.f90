program simple_test_watcher
use simple_syscalls
use simple_defs
use simple_qsys_funs
use simple_filehandling
implicit none

character(len=STDLEN) :: fbody='./Movies/FoilHole', ext='mrc', filetabname='movies.txt'
integer               :: nmovies, nmovies_prev, i, nmovies_in_queue
logical, allocatable  :: lmask_streaming(:), ltmp(:)
integer, parameter    :: SHORTTIME = 5
logical, parameter    :: DEBUG = .true.

nmovies = 0
do 
    call sys_gen_filetab(fbody, ext, filetabname)
    nmovies_prev = nmovies
    nmovies      = nlines(filetabname)
    if( nmovies > nmovies_prev )then
        ! there are new movies to process
        if( DEBUG ) print *, 'nmovies updated to: ', nmovies
        if( allocated(lmask_streaming) )then
            nmovies_in_queue = size(lmask_streaming)
        else
            nmovies_in_queue = 0
        endif
        ! update the queue
        if( nmovies > nmovies_in_queue )then
            allocate(ltmp(nmovies))
            ltmp = .true.
            ltmp(1:nmovies_in_queue) = lmask_streaming(1:nmovies_in_queue)
            if( allocated(lmask_streaming) ) deallocate(lmask_streaming)
            allocate(lmask_streaming(1:nmovies), source=ltmp)
            deallocate(ltmp)
            if( DEBUG )then
                do i=1,nmovies
                    print *, i, 'lmask: ', lmask_streaming(i)
                end do
            endif
        else
            write(*,'(a)') 'WEIRD! nmovies <= nmovies_in_queue, despite that&
            &there are new movies in filetab'
        endif
    endif
    call sleep(SHORTTIME)
end do




end program simple_test_watcher
