program simple_plot_textremal
use simple_defs
implicit none
integer, parameter :: NITS  = 20
real,    parameter :: RRATE = 0.8
real    :: extr_thresh
integer :: i
extr_thresh = EXTRINITHRESH/RRATE
do i=1,NITS
    extr_thresh = extr_thresh * RRATE
    print *, real(i), extr_thresh
end do
end program simple_plot_textremal
