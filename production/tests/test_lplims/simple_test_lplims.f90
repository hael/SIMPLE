program simple_test_lplims
use simple_estimate_ssnr
implicit none
real :: mskdiam, lpstart,lpstop, lpcen
mskdiam = 300.
do while( mskdiam >= 90.)
    call mskdiam2lplimits( mskdiam, lpstart,lpstop, lpcen )
    print *, 'mskdiam/lpstart/lpstop/lpcen: ', mskdiam, lpstart, lpstop, lpcen
    mskdiam = mskdiam - 10
end do
end program simple_test_lplims
