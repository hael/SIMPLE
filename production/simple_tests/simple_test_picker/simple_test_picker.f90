program simple_test_picker
!$ use omp_lib
!$ use omp_lib_kinds
use simple_picker ! singleton
implicit none
integer,           parameter   :: NTHR=8, OFFSET=3, MAXKMIT=20
real,              parameter   :: SHRINK=4.0, LP=20.0, SMPD=1.77, MSK=100.0, DISTTHR=160.0
character(len=32), parameter   :: micname='ribo_first_intg011.mrc', refsname='boxrefs.mrc'

! keep the box size to ~60
! no more than 100 references, generated with diverse option in makeoris
!$ call omp_set_num_threads(NTHR)

call init_picker(micname, refsname, SMPD, MSK, OFFSET, SHRINK, LP, DISTTHR)
call pick_particles
call kill_picker

end program simple_test_picker