program simple_test_corrweights
use simple_fileio
use gnufor2
use simple_math
implicit none
integer, parameter :: NCORRS=38, NDISTR=5
real               :: corrs(NCORRS), x(NCORRS), ws(NCORRS)
character(len=32)  :: fname='corrs4wtst.txt'
integer            :: funit, i, order(NCORRS)

do i=1,NCORRS
    x(i) = real(i)
end do

call fopen(funit, fname)
do i=1,NDISTR

    read(funit,*) corrs
    call corrs2w_softmax
    ! print *, corrs
    call plot(x, ws)

end do
call fclose(funit)

contains

    subroutine corrs2w_factorial
        real    :: logws(NCORRS)
        integer :: order(NCORRS), ipeak
        ws    = exp(-(1.-corrs))
        logws = log(ws)
        order = (/(ipeak,ipeak=1,NCORRS)/)
        call hpsort(logws, order)
        call reverse(order)
        call reverse(logws)
        forall(ipeak=1:NCORRS) ws(order(ipeak)) = exp(sum(logws(:ipeak)))
        ws = ws / sum(ws)
    end subroutine corrs2w_factorial

    subroutine corrs2w_softmax
        real, parameter :: TAU=0.02
        real :: corrmin, corrmax, corrdiff 

        corrmax  = maxval(corrs)
        corrmin  = minval(corrs)
        corrdiff = corrmax - corrmin

        print *, corrdiff

        ws = exp(-(1.-corrs/TAU))
        ws = ws / sum(ws)
    end subroutine corrs2w_softmax


end program simple_test_corrweights
