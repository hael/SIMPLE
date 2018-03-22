module  simple_test_pearsons_as192
implicit none

public ::  test_pearsons

contains

    subroutine test_pearsons
        integer :: itype, ifault
        real    :: xbar, sd, rb1, b2, bndry, sigpt(11)

        do
            write(*, '(a)', advance='no') ' Enter mean & st.devn.: '
            read(*, *) xbar, sd
            write(*, *) 'ITYPE = 1  4 moments used'
            write(*, *) 'ITYPE = 2  3 moments + left boundary used'
            write(*, *) 'ITYPE = 3  3 moments + right boundary used'
            write(*, '(a)', advance='no') ' Enter itype: '
            read(*, *) itype
            write(*, '(a)', advance='no') ' Enter root(beta1) with sign of beta1: '
            read(*, *) rb1
            if (itype == 1) then
                write(*, '(a)', advance='no') ' Enter beta2: '
                read(*, *) b2
            end if
            if (itype > 1) then
                write(*, '(a)', advance='no') ' Enter boundary: '
                read(*, *) bndry
            end if

            call pearsn_192(xbar, sd, itype, rb1, b2, bndry, sigpt, ifault)
            if (ifault /= 0) then
                write(*, *) 'ifault = ', ifault
            else
                write(*, '(" Percentage points"/ " ", 11f7.2)') sigpt
            end if

        end do

        stop
    end subroutine test_pearsons


    subroutine pearsn_192(xbar, sd, itype, rb1, b2, bndry, sigpt, ifault)

        ! Code converted using TO_F90 by Alan Miller
        ! Date: 1999-09-25  Time: 16:22:51
        ! N.B. Array a has been removed from the list of arguments,
        !      and replaced with a PARAMETER statement in this routine.

        !     ALGORITHM AS 192  APPL. STATIST. (1983) VOL.32, NO.3

        !     Computes approximate significance points of a Pearson curve with given
        !     first four moments, or first three moments and left or right boundary.

        ! Error codes
        ! IFAULT = 0  successful completion
        !        = 1 SD < 0.0
        !        = 2 ITYPE < 1 or > 3
        !        = 3 | RB1 | > 2.0
        !        = 4 XBAR value impossible with the value entered for BNDRY
        !        = 5 B2 cannot be computed for the first 3 moments + left boundary
        !        = 6 B2 out of range, that is:
        !            If | RB1 | <= 1.0 then
        !               B2 < 1.5 * | RB1 | + 1.5  or  B2 > 0.2 * | RB1 | + 10.8
        !            If | RB1 | > 1.0 then
        !               B2 < 3.9 * | RB1 | - 0.9  or  B2 > 4.8 * | RB1 | + 6.2
        implicit none
        real, intent(in)      :: xbar !
        real, intent(in)      :: sd
        integer, intent(in)   :: itype
        real, intent(in)      :: rb1
        real, intent(inout)  :: b2
        real, intent(in)      :: bndry
        real, intent(out)     :: sigpt(:)
        integer, intent(out)  :: ifault

        real    :: b(9), rbeta1, ca, ca2, carb1, denom, absca, bnd1, bnd2, sgn, &
            pi1, pi2
        integer :: i, j, k, l
        real, parameter :: zero = 0.0, point2 = 0.2, rthalf = 0.70710678,  &
            point9 = 0.9, one = 1.0, onept5 = 1.5, two = 2.0, &
            three = 3.0, thrpt9 = 3.9, four = 4.0, forpt8 = 4.8, &
            sixpt2 = 6.2, tenpt8 = 10.8
        real, parameter :: a(418) = (/    &
            -1.9336, -1.6061, 2.6955, -1.7036, 1.7236, -1.3209, -2.5464, 1.6812, -0.11812, 0.050875, 0.16042, -1.0616, 0.43929, &
            -0.35799, 0.51040, 1.0958, -0.63453, 0.053175, -0.019945, 0.25597, 0.32520, -0.094571, 2.0418, -2.5786, 0.61148, -0.64242, &
            1.1968, -0.33145, 0.015508, -4.1444, 0.13884, 4.6625, 0.37661, -0.20447, -1.9508, -0.20214, 0.11921, -0.0085268, &

            -1.6453, -1.8494, 2.3915, -1.5844, 1.8290, -1.1167, -3.6091, 2.7150, -0.45857, 0.048102, 0.89117, -1.2700, 0.37653, &
            -0.81542, 0.56815, 1.6562, -1.1908, 0.24574, -0.024375, -1.5063, 4.4876, -0.60765, -6.5584, 2.8944, -0.42381, 2.2664, &
            -1.3425, 0.25806, -0.010421, -1.9526, 0.21332, 2.1317, -1.1996, 0.21033, -0.52154, 0.53597, -0.12355, 0.0060658, &

            -1.1044, -1.1300, 1.7681, 0.10947, 0.30566, -0.90598, -0.98814, 0.49355, 0.30625, -0.018707, 0.74941, -1.3174, 0.23974, &
            -0.12433, 0.59105, 0.85714, -0.47292, -0.16028, 0.012684, -1.4645, 4.5349, -0.71053, -7.9213, 4.0018, -0.61932, 2.7320, &
            -1.7752, 0.35919, -0.014944, -2.0999, 0.30754, 3.6865, -2.3140, 0.39961, -1.1530, 1.0257, -0.22981, 0.011032, &

            -0.85842, 0.78929, 1.2196, -0.55088, -0.96385, -0.57312, -0.076974, 0.34882, 0.47311, -0.056120, -0.63153, -1.2697, 0.64649, &
            0.92712, 0.49088, 0.48118, -0.46407, -0.39753, 0.052413, -0.92873, 2.2039, 0.25828, 0.74312, -2.4307, 0.49934, 0.029833, &
            0.68249, -0.19687, 0.0068136, -2.9758, -0.050738, 0.067069, 2.0835, -0.45211, -0.28149, -0.58012, 0.17621, -0.0067195, &

            -2.4031, -1.4891, 0.75590, -2.3309, 9.2314, -3.8930, -7.0192, 2.3648, 0.44308, -0.029437, 16.616, -4.9826, 1.0626, &
            -18.888, 7.2830, 14.642, -2.6394, -1.2739, 0.060862, -1.6092, 0.81385, -0.059847, -0.18180, 1.2679, -0.36205, 0.29683, &
            -0.83310, 0.26630, -0.015081, 0.43536, 0.27965, 1.2516, -3.0536, 0.65936, -1.0579, 1.7886, -0.47962, 0.025311, &

            -2.0711E-4, 0.0068277, -7.5482E-6, 0.47086, -0.21213, 1.1003E-4, -0.22462, 0.11733, -0.079682, -1.0565E-5, 2.5087, -2.6169, &
            0.94103, -2.3261, 1.6482, 2.3024, -1.4821, -0.020294, 9.9168E-4, -0.097078, 0.58192, -0.019936, -0.7025, 0.11526, 0.0013959, &
            0.1595, -0.062382, 0.01338, 1.8557E-4, -1.8383, -0.54124, 0.68041, 1.7654, -0.37606, -1.113, 0.044447, 0.042456, -0.0017448, &

            0.39672, -0.17140, -0.74246, 0.48553, 0.19650, 0.58942, -1.5768, 0.37082, -0.32481, 0.025296, -0.60818, -1.9005, 0.45451, &
            0.79612, 1.0544, -1.9010, 0.29849, -0.51728, 0.050926, 0.45223, -2.4004, 0.38372, 3.6817, -1.0796, 0.016619, -1.6829, &
            0.81147, -0.10376, 0.0052116, -4.3837, 0.67456, 6.3029, -1.8042, 0.023516, -2.9320, 1.4537, -0.20680, 0.010634, &

            0.77212, 0.028272, -1.3932, 0.95628, -0.37918, 0.85487, 0.38225, -1.3009, 0.42400, -0.051933, 0.17133, -1.5330, 1.2721, &
            -0.58446, 0.78993, 0.33412, -1.2988, 0.39155, -0.046683, 2.2548, -1.7005, -0.65377, -2.9256, 2.1397, 0.093244, 3.2094, &
            -2.2896, 0.24580, -0.0020536, 0.30478, -0.58544, -4.0524, 1.8239, 0.083952, 3.4072, -2.0163, 0.20345, -7.5748E-4, &

            1.1504, 0.22674, -1.7407, 1.1356, -0.58596, 0.88839, 0.34743, -1.2531, 0.42319, -0.014362, 0.17074, -1.2892, 1.1653, &
            -0.46699, 0.58479, 0.33975, -1.0005, 0.25271, -0.0091547, 2.1265, -8.7599, 2.3560, 8.0491, -3.3844, -0.028416, -4.3522, &
            3.9741, -0.92840, 0.055908, -5.5206, 1.6024, 6.1077, -2.4508, -0.012514, -3.7719, 2.8877, -0.59134, 0.034299, &

            1.5989, 1.5353, -2.2124, 0.83145, -1.2543, 0.92247, 2.2076, -1.4147, 0.053879, 0.058513, 0.42581, -1.1636, 0.67165, &
            -0.47849, 0.46159, 0.97601, -0.61165, -0.054400, 0.029952, 10.803, -5.6739, -8.3559, -12.404, 29.846, -7.7439, -6.2673, &
            -1.3612, 1.3966, -0.22122, 2.8678, -3.0580, -14.864, 16.489, -4.0334, -0.49139, -1.7927, 0.98545, -0.10118, &

            2.4201, -1.9281, -3.3357, 2.3720, 1.7318, 1.4149, -0.64616, -2.7558, 0.28474, -0.034471, -0.078628, -1.2092, 0.43924, &
            0.43093, 0.53223, 0.34235, -1.0741, 0.080494, -0.013004, -15.787, -3.9798, 23.933, 24.332, -46.762, 6.0862, 15.874, &
            5.2360, -2.4644, 0.28404, -14.830, 8.5161, 23.701, -19.419, 2.4239, -1.8457, 4.8007, -1.2525, 0.099997 /)

        !     Check for invalid input arguments

        ifault = 1
        if (sd <= zero) return
        ifault = 2
        if (itype < 1 .or. itype > 3) return
        rbeta1 = ABS(rb1)
        ifault = 3
        if (rbeta1 > two) return
        if (itype == 1) go to 10

        !     Compute B2 using XBAR, SD, RB1, and the known boundary

        ifault = 4
        if ((itype == 2 .and. xbar <= bndry) .or. (itype == 3 .and. xbar >= bndry))  &
            return
        ifault = 5
        ca = sd / (xbar - bndry)
        ca2 = ca * ca
        carb1 = ca * rb1
        denom = four * ca2 - carb1 + two
        if (denom <= zero) return
        absca = ABS(ca)
        bnd1 = (ca2 - one) / ca
        if (absca <= rthalf) then
            bnd2 = four * ca / (one - ca2)
        else
            bnd2 = four * ca + two / ca
        end if
        if ((itype == 2 .and. (rb1 < bnd1 .or. rb1 > bnd2)) .or.  &
            (itype == 3 .and. (rb1 < bnd2 .or. rb1 > bnd1))) return
        b2 = three * (carb1 * (carb1 + one) + rb1 * rb1 + two) / denom

10      ifault = 6
        if ((rbeta1 <= one .and. (b2 < onept5 * rbeta1 + onept5 .or.  &
            b2 > point2 * rbeta1 + tenpt8)) .or. (rbeta1 > one .and.  &
            (b2 < thrpt9 * rbeta1 - point9 .or. b2 > forpt8 * rbeta1 + sixpt2))) return
        ifault = 0

        !     Initialization

        b(1) = rbeta1
        b(2) = b2
        b(3) = rbeta1 * rbeta1
        b(4) = rbeta1 * b2
        b(5) = b2 * b2
        b(6) = b(3) * rbeta1
        b(7) = b(3) * b2
        b(8) = b(5) * rbeta1
        b(9) = b(5) * b2

        !     Significance point computation

        k = -19
        if (rbeta1 > one) k = 0
        sgn = one
        if (rb1 < zero) sgn = -one
        do  i = 1, 11
            l = i
            if (rb1 < zero) l = 12 - i
            k = k + 20
            pi1 = a(k)
            do  j = 1, 9
                k = k + 1
                pi1 = pi1 + a(k) * b(j)
            end do
            pi2 = one
            do  j = 1, 9
                k = k + 1
                pi2 = pi2 + a(k) * b(j)
            end do

            !     Compute significance point

            sigpt(l) = xbar + sd * sgn * pi1 / pi2
        end do
        return
    end subroutine pearsn_192

end module simple_test_pearsons_as192
