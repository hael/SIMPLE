program simple_test_watcher
use simple_syscalls
use simple_defs
use simple_qsys_funs
use simple_filehandling
implicit none
character(len=STDLEN) :: fbody='./Movies/FoilHole', ext='mrc', filetabname='movies.txt'
integer :: nmovies, nmovies_prev

nmovies = 0
do 
    call sys_gen_filetab( fbody, ext, filetabname )
    nmovies_prev = nmovies
    nmovies      = nlines(filetabname)
    if( nmovies > nmovies_prev )then
        ! there are new movies to process
        print *, 'nmovies updated to: ', nmovies

    endif
end do




end program simple_test_watcher
