program simple_test_phaseplate_correct_fsc
use simple_math
implicit none
real :: fsc(72)
integer :: find_plate, k
fsc = 0.
fsc(1) = 0.470
fsc(2) = 0.470
fsc(3) = 0.470
fsc(4) = 0.847
fsc(5) = 0.780
fsc(6) = 0.774
fsc(7) = 0.746
fsc(8) = 0.599
fsc(9) = 0.766
fsc(10) = 0.715
fsc(11) = 0.749
fsc(12) = 0.765
fsc(13) = 0.777
fsc(14) = 0.701
fsc(15) = 0.572
fsc(16) = 0.473
fsc(17) = 0.447
fsc(18) = 0.472
fsc(19) = 0.425
fsc(20) = 0.287
fsc(21) = 0.086
fsc(22) = 0.005
fsc(23) = 0.001
fsc(24) = 0.035
fsc(25) = 0.023
fsc(26) = 0.008
fsc(27) = -0.025
fsc(28) = -0.002
fsc(29) = 0.019
fsc(30) = 0.025
fsc(31) = 0.003
fsc(32) = -0.022
fsc(33) = 0.030
fsc(34) = 0.004
fsc(35) = -0.001
fsc(36) = -0.001
fsc(37) = 0.002
fsc(38) = 0.019
fsc(39) = 0.005
fsc(40) = 0.024
fsc(41) = 0.009
fsc(42) = -0.001
fsc(43) = 0.020
fsc(44) = 0.035
fsc(45) = 0.007
fsc(46) = -0.001
fsc(47) = 0.001
fsc(48) = 0.013
fsc(49) = -0.009
fsc(50) = -0.000
fsc(51) = 0.004
fsc(52) = -0.003
fsc(53) = -0.014
fsc(54) = 0.008
fsc(55) = -0.019
fsc(56) = 0.003
fsc(57) = -0.019
fsc(58) =  -0.015
fsc(59) = -0.016
fsc(60) = -0.062
fsc(61) = -0.023
fsc(62) = 0.008
fsc(63) = -0.011
fsc(64) = 0.073
fsc(65) = 0.005
fsc(66) = 0.041
fsc(67) = 0.034
fsc(68) = 0.132
fsc(69) = -0.133
fsc(70) = -0.034
fsc(71) = -0.070
fsc(72) = -0.046
call phaseplate_correct_fsc( fsc, find_plate )
do k=1,72
    print *, k, fsc(k)
end do
print *, 'find_plate: ', find_plate

end program simple_test_phaseplate_correct_fsc