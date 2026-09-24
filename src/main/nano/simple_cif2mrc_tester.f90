!@descr: analytical all-element tests for CIF-to-MRC atomic density conversion
module simple_cif2mrc_tester
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_atoms,      only: atoms
use simple_defs,       only: dp
use simple_image,      only: image
use simple_imghead,    only: find_ldim_nptcls, find_img_smpd
use simple_string,     only: string
use simple_syslib,     only: del_file, file_exists
use simple_test_utils, only: assert_int, assert_real, assert_true
implicit none
private
public :: run_all_cif2mrc_tests

integer, parameter :: NELEMENTS = 96, EXPECTED_BOX = 52
real,    parameter :: SMPD = 1.0
real(dp), parameter :: MAX_ABS_TOL = 2.0e-3_dp, RMS_TOL = 2.0e-4_dp, REL_SUM_TOL = 5.0e-5_dp
character(len=*), parameter :: CIF_FILE = 'cif2mrc_test.cif'
character(len=*), parameter :: PDB_FILE = 'cif2mrc_test.pdb'
character(len=*), parameter :: MRC_FILE = 'cif2mrc_test.mrc'
character(len=*), parameter :: CENTERED_PDB_FILE = 'cif2mrc_test_centered.pdb'
character(len=2), parameter :: ELEMENTS(NELEMENTS) = [character(len=2) :: &
    &'H ', 'HE', 'LI', 'BE', 'B ', 'C ', 'N ', 'O ', 'F ', 'NE', 'NA', 'MG', 'AL', 'SI', 'P ', 'S ', &
    &'CL', 'AR', 'K ', 'CA', 'SC', 'TI', 'V ', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'GA', 'GE', &
    &'AS', 'SE', 'BR', 'KR', 'RB', 'SR', 'Y ', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'CD', &
    &'IN', 'SN', 'SB', 'TE', 'I ', 'XE', 'CS', 'BA', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', &
    &'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'HF', 'TA', 'W ', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', &
    &'TL', 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', 'AC', 'TH', 'PA', 'U ', 'NP', 'PU', 'AM', 'CM']

contains

    subroutine run_all_cif2mrc_tests
        type(atoms) :: molecule
        type(image) :: volume
        real, pointer :: density(:,:,:)
        real(dp), allocatable :: reference(:,:,:), error(:,:,:)
        real(dp) :: max_abs_error, rms_error, rel_sum_error, local_error
        real :: smpd_out
        integer :: ldim(3), nptcls, z, peak_actual(3), peak_reference(3), center(3), lo(3), hi(3)
        character(len=160) :: message

        call write_all_element_fixture
        call molecule%cif2mrc(string(CIF_FILE), SMPD)
        call assert_int(NELEMENTS, molecule%get_n(), 'cif2mrc parses all 96 supported elements')
        do z = 1, NELEMENTS
            write(message,'(A,A,A,I0)') 'cif2mrc atomic number for ', trim(ELEMENTS(z)), ' is ', z
            call assert_int(z, molecule%get_atomicnumber(z), trim(message))
        enddo
        call assert_true(file_exists(MRC_FILE), 'cif2mrc creates the all-element MRC volume')
        if( file_exists(MRC_FILE) )then
            call find_ldim_nptcls(string(MRC_FILE), ldim, nptcls)
            smpd_out = find_img_smpd(string(MRC_FILE))
            call assert_int(EXPECTED_BOX, ldim(1), 'cif2mrc all-element x dimension')
            call assert_int(EXPECTED_BOX, ldim(2), 'cif2mrc all-element y dimension')
            call assert_int(EXPECTED_BOX, ldim(3), 'cif2mrc all-element z dimension')
            call assert_int(ldim(3), nptcls, 'cif2mrc MRC header reports one section per z plane')
            call assert_real(SMPD, smpd_out, 1.0e-6, 'cif2mrc preserves the requested sampling distance')
            call volume%new(ldim, smpd_out)
            call volume%read(string(MRC_FILE))
            call volume%get_rmat_ptr(density)
            call build_analytical_reference(ldim, reference)
            allocate(error(ldim(1), ldim(2), ldim(3)))
            error = abs(real(density, dp) - reference)
            max_abs_error = maxval(error)
            rms_error = sqrt(sum(error**2) / real(size(error), dp))
            rel_sum_error = abs(sum(real(density, dp)) - sum(reference)) / abs(sum(reference))
            peak_actual = maxloc(density)
            peak_reference = maxloc(reference)
            write(message,'(A,ES12.4,A,ES12.4)') 'all-element maximum voxel error=', max_abs_error, &
                &', tolerance=', MAX_ABS_TOL
            call assert_true(max_abs_error <= MAX_ABS_TOL, trim(message))
            write(message,'(A,ES12.4,A,ES12.4)') 'all-element RMS voxel error=', rms_error, ', tolerance=', RMS_TOL
            call assert_true(rms_error <= RMS_TOL, trim(message))
            write(message,'(A,ES12.4,A,ES12.4)') 'all-element relative density-sum error=', rel_sum_error, &
                &', tolerance=', REL_SUM_TOL
            call assert_true(rel_sum_error <= REL_SUM_TOL, trim(message))
            call assert_true(all(peak_actual == peak_reference), 'cif2mrc density peak matches the analytical reference')
            call assert_true(all(ieee_is_finite(density)), 'cif2mrc density contains only finite values')
            do z = 1, NELEMENTS
                center = nint(expected_position(z)) + 1
                lo = max(1, center - 8)
                hi = min(ldim, center + 8)
                local_error = maxval(error(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)))
                write(message,'(A,A,ES12.4,A,ES12.4)') trim(ELEMENTS(z)), &
                    &' local maximum voxel error=', local_error, ', tolerance=', MAX_ABS_TOL
                call assert_true(local_error <= MAX_ABS_TOL, trim(message))
            enddo
        endif
        if( allocated(reference) ) deallocate(reference)
        if( allocated(error) ) deallocate(error)
        call volume%kill
        call molecule%kill
        call del_file(CIF_FILE)
        call del_file(PDB_FILE)
        call del_file(MRC_FILE)
        call del_file(CENTERED_PDB_FILE)
    end subroutine run_all_cif2mrc_tests

    subroutine write_all_element_fixture
        integer :: file_unit, z
        real(dp) :: xyz(3)
        open(newunit=file_unit, file=CIF_FILE, status='replace', action='write')
        write(file_unit,'(A)') 'data_cif2mrc_all_elements'
        write(file_unit,'(A)') 'loop_'
        write(file_unit,'(A)') '_atom_site.group_PDB'
        write(file_unit,'(A)') '_atom_site.id'
        write(file_unit,'(A)') '_atom_site.type_symbol'
        write(file_unit,'(A)') '_atom_site.label_atom_id'
        write(file_unit,'(A)') '_atom_site.label_alt_id'
        write(file_unit,'(A)') '_atom_site.label_comp_id'
        write(file_unit,'(A)') '_atom_site.label_asym_id'
        write(file_unit,'(A)') '_atom_site.label_seq_id'
        write(file_unit,'(A)') '_atom_site.Cartn_x'
        write(file_unit,'(A)') '_atom_site.Cartn_y'
        write(file_unit,'(A)') '_atom_site.Cartn_z'
        write(file_unit,'(A)') '_atom_site.occupancy'
        write(file_unit,'(A)') '_atom_site.B_iso_or_equiv'
        do z = 1, NELEMENTS
            xyz = fixture_position(z)
            write(file_unit,'(A,1X,I0,1X,A,1X,A,1X,A,1X,A,1X,A,1X,I0,3(1X,F8.3),2(1X,F5.1))') &
                &'ATOM', z, trim(ELEMENTS(z)), trim(ELEMENTS(z)), '.', 'MOL', 'A', z, xyz, 1.0, 0.0
        enddo
        write(file_unit,'(A)') '#'
        close(file_unit)
    end subroutine write_all_element_fixture

    pure function fixture_position( z ) result( xyz )
        integer, intent(in) :: z
        real(dp) :: xyz(3)
        xyz(1) = 4._dp * (real(mod(z - 1, 6), dp) - 2.5_dp)
        xyz(2) = 4._dp * (real(mod((z - 1) / 6, 4), dp) - 1.5_dp)
        xyz(3) = 4._dp * (real((z - 1) / 24, dp) - 1.5_dp)
    end function fixture_position

    pure function expected_position( z ) result( xyz )
        integer, intent(in) :: z
        real(dp) :: xyz(3)
        xyz = fixture_position(z) + 0.5_dp * real(EXPECTED_BOX - 1, dp)
    end function expected_position

    subroutine build_analytical_reference( ldim, reference )
        integer, intent(in) :: ldim(3)
        real(dp), allocatable, intent(out) :: reference(:,:,:)
        real(dp), parameter :: C = 2132.79_dp
        real(dp), parameter :: FOURPISQ = 4._dp * acos(-1._dp)**2
        real(dp), parameter :: BFAC = 64._dp, CUTOFFSQ = 64._dp
        real(dp) :: a(5), b(5), aterm(5), xyz(3), r2
        integer :: z, i, j, k
        allocate(reference(ldim(1), ldim(2), ldim(3)), source=0._dp)
        do z = 1, NELEMENTS
            call reference_coefficients(z, a, b)
            b = b + BFAC
            aterm = a / b**1.5_dp
            xyz = expected_position(z)
            do k = 1, ldim(3)
                do j = 1, ldim(2)
                    do i = 1, ldim(1)
                        r2 = sum((xyz - real([i - 1, j - 1, k - 1], dp))**2)
                        if( r2 > CUTOFFSQ ) cycle
                        reference(i,j,k) = reference(i,j,k) + C * sum(aterm * exp(-FOURPISQ * r2 / b))
                    enddo
                enddo
            enddo
        enddo
    end subroutine build_analytical_reference

    subroutine reference_coefficients( z, a, b )
        integer,  intent(in)  :: z
        real(dp), intent(out) :: a(5), b(5)
            select case(z)
            case(1) ! H
                a = [0.0349, 0.1201, 0.1970, 0.0573, 0.1195]
                b = [0.5347, 3.5867, 12.3471, 18.9525, 38.6269]
            case(2) ! HE
                a = [0.0317, 0.0838, 0.1526, 0.1334, 0.0164]
                b = [0.2507, 1.4751, 4.4938, 12.6646, 31.1653]
            case(3) ! LI
                a = [0.0750, 0.2249, 0.5548, 1.4954, 0.9354]
                b = [0.3864, 2.9383, 15.3829, 53.5545, 138.7337]
            case(4) ! BE
                a = [0.0780, 0.2210, 0.6740, 1.3867, 0.6925]
                b = [0.3131, 2.2381, 10.1517, 30.9061, 78.3273]
            case(5) ! B
                a = [0.0909, 0.2551, 0.7738, 1.2136, 0.4606]
                b = [0.2995, 2.1155, 8.3816, 24.1292, 63.1314]
            case(6) ! C
                a = [0.0893, 0.2563, 0.7570, 1.0487, 0.3575]
                b = [0.2465, 1.7100, 6.4094,18.6113,50.2523]
            case(7) ! N
                a = [0.1022, 0.3219, 0.7982, 0.8197, 0.1715]
                b = [0.2451, 1.7481, 6.1925,17.3894,48.1431]
            case(8) ! O
                a = [0.0974, 0.2921, 0.6910,  0.6990, 0.2039]
                b = [0.2067, 1.3815, 4.6943, 12.7105,32.4726]
            case(9) ! F
                a = [0.1083, 0.3175, 0.6487, 0.5846, 0.1421]
                b = [0.2057, 1.3439, 4.2788, 11.3932, 28.7881]
            case(10) ! NE
                a = [0.1269, 0.3535, 0.5582, 0.4674, 0.1460]
                b = [0.2200, 1.3779, 4.0203, 9.4934, 23.1278]
            case(11) ! NA
                a = [0.2142, 0.6853, 0.7692, 1.6589, 1.4482]
                b = [0.3334, 2.3446, 10.0830, 48.3037, 138.2700]
            case(12) ! MG
                a = [0.2314, 0.6866, 0.9677, 2.1882, 1.1339]
                b = [0.3278, 2.2720, 10.9241, 39.2898, 101.9748]
            case(13) ! AL
                a = [0.2390, 0.6573, 1.2011, 2.5586, 1.2312]
                b = [0.3138, 2.1063, 10.4163, 34.4552, 98.5344]
            case(14) ! SI
                a = [0.2519, 0.6372, 1.3795, 2.5082, 1.0500]
                b = [0.3075, 2.0174, 9.6746, 29.3744, 80.4732]
            case(15) ! P
                a = [0.2548, 0.6106, 1.4541, 2.3204, 0.8477]
                b = [0.2908, 1.8740, 8.5176, 24.3434, 63.2996]
            case(16) ! S
                a = [0.2497, 0.5628, 1.3899, 2.1865, 0.7715]
                b = [0.2681, 1.6711, 7.0267, 19.5377, 50.3888]
            case(17) ! CL
                a = [0.2443, 0.5397, 1.3919, 2.0197, 0.6621]
                b = [0.2468, 1.5242, 6.1537, 16.6687, 42.3086]
            case(18) ! AR
                a = [0.2385, 0.5017, 1.3428, 1.8899, 0.6079]
                b = [0.2289, 1.3694, 5.2561, 14.0928, 35.5361]
            case(19) ! K
                a = [0.4115, 1.4031, 2.2784, 2.6742, 2.2162]
                b = [0.3703, 3.3874, 13.1029, 68.9592, 194.4329]
            case(20) ! CA
                a = [0.4054, 1.3880, 2.1602, 3.7532, 2.2063]
                b = [0.3499, 3.0991, 11.9608, 53.9353, 142.3892]
            case(21) ! SC
                a = [0.3787, 1.2181, 2.0594, 3.2618, 2.3870]
                b = [0.3133, 2.5856, 9.5813, 41.7688, 116.7282]
            case(22) ! TI
                a = [0.3825, 1.2598, 2.0008, 3.0617, 2.0694]
                b = [0.3040, 2.4863, 9.2783, 39.0751, 109.4583]
            case(23) ! V
                a = [0.3876, 1.2750, 1.9109, 2.8314, 1.8979]
                b = [0.2967, 2.3780, 8.7981, 35.9528, 101.7201]
            case(24) ! CR
                a = [0.4046, 1.3696, 1.8941, 2.0800, 1.2196]
                b = [0.2986, 2.3958, 9.1406, 37.4701, 113.7121]
            case(25) ! MN
                a = [0.3796, 1.2094, 1.7815, 2.5420, 1.5937]
                b = [0.2699, 2.0455, 7.4726, 31.0604, 91.5622]
            case(26) ! Fe
                a = [0.3946, 1.2725, 1.7031,  2.3140,  1.4795]
                b = [0.2717, 2.0443, 7.6007, 29.9714, 86.2265]
            case(27) ! CO
                a = [0.4118, 1.3161, 1.6493, 2.1930, 1.2830]
                b = [0.2742, 2.0372, 7.7205, 29.9680, 84.9383]
            case(28) ! NI
                a = [0.3860, 1.1765, 1.5451, 2.0730, 1.3814]
                b = [0.2478, 1.7660, 6.3107, 25.2204, 74.3146]
            case(29) ! CU
                a = [0.4314, 1.3208, 1.5236, 1.4671, 0.8562]
                b = [0.2694, 1.9223, 7.3474, 28.9892, 90.6246]
            case(30) ! CU
                a = [0.4288, 1.2646, 1.4472, 1.8294, 1.0934]
                b = [0.2593, 1.7998, 6.7500, 25.5860, 73.5284]
            case(31) ! GA
                a = [0.4818, 1.4032, 1.6561, 2.4605, 1.1054]
                b = [0.2825, 1.9785, 8.7546, 32.5238, 98.5523]
            case(32) ! GE
                a = [0.4655, 1.3014, 1.6088, 2.6998, 1.3003]
                b = [0.2647, 1.7926, 7.6071, 26.5541, 77.5238]
            case(33) ! AS
                a = [0.4517, 1.2229, 1.5852, 2.7958, 1.2638]
                b = [0.2493, 1.6436, 6.8154, 22.3681, 62.0390]
            case(34) ! SE
                a = [0.4477, 1.1678, 1.5843, 2.8087, 1.1956]
                b = [0.2405, 1.5442, 6.3231, 19.4610, 52.0233]
            case(35) ! BR
                a = [0.4798, 1.1948, 1.8695, 2.6953, 0.8203]
                b = [0.2504, 1.5963, 6.9653, 19.8492, 50.3233]
            case(36) ! KR
                a = [0.4546, 1.0993, 1.7696, 2.7068, 0.8672]
                b = [0.2309, 1.4279, 5.9449, 16.6752, 42.2243]
            case(37) ! RB
                a = [1.0160, 2.8528, 3.5466, -7.7804, 12.1148]
                b = [0.4853, 5.0925, 25.7851, 130.4515, 138.6775]
            case(38) ! SR
                a = [0.6703, 1.4926, 3.3368, 4.4600, 3.1501]
                b = [0.3190, 2.2287, 10.3504, 52.3291, 151.2216]
            case(39) ! Y
                a = [0.6894, 1.5474, 3.2450, 4.2126, 2.9764]
                b = [0.3189, 2.2904, 10.0062, 44.0771, 125.0120]
            case(40) ! ZR
                a = [0.6719, 1.4684, 3.1668, 3.9557, 2.8920]
                b = [0.3036, 2.1249, 8.9236, 36.8458, 108.2049]
            case(41) ! NB
                a = [0.6123, 1.2677, 3.0348, 3.3841, 2.3683]
                b = [0.2709, 1.7683, 7.2489, 27.9465, 98.5624]
            case(42) ! MO
                a = [0.6773, 1.4798, 3.1788, 3.0824, 1.8384]
                b = [0.2920, 2.0606, 8.1129, 30.5336, 100.0658]
            case(43) ! TC
                a = [0.7082, 1.6392, 3.1993, 3.4327, 1.8711]
                b = [0.2976, 2.2106, 8.5246, 33.1456, 96.6377]
            case(44) ! RU
                a = [0.6735, 1.4934, 3.0966, 2.7254, 1.5597]
                b = [0.2773, 1.9716, 7.3249, 26.6891, 90.5581]
            case(45) ! RH
                a = [0.6413, 1.3690, 2.9854, 2.6952, 1.5433]
                b = [0.2580, 1.7721, 6.3854, 23.2549, 85.1517]
            case(46) ! PD
                a = [0.5904, 1.1775, 2.6519, 2.2875, 0.8689]
                b = [0.2324, 1.5019, 5.1591, 15.5428, 46.8213]
            case(47) ! AG
                a = [0.6377, 1.3790, 2.8294, 2.3631, 1.4553]
                b = [0.2466, 1.674, 5.7656, 20.0943, 76.7372]
            case(48) ! CD
                a = [0.6364, 1.4247, 2.7802, 2.5973, 1.7886]
                b = [0.2407, 1.6823, 5.6588, 20.7219, 69.1409]
            case(49) ! IN
                a = [0.6768, 1.6589, 2.7740, 3.1835, 2.1326]
                b = [0.2522, 1.8545, 6.2936, 25.1457, 84.5448]
            case(50) ! SN
                a = [0.7224, 1.9610, 2.7161, 3.5603, 1.8972]
                b = [0.2651, 2.0604, 7.3011, 27.5493, 81.3349]
            case(51) ! SB
                a = [0.7106, 1.9247, 2.6149, 3.8322, 1.8899]
                b = [0.2562, 1.9646, 6.8852, 24.7648, 68.9168]
            case(52) ! TE
                a = [0.6947, 1.8690, 2.5356, 4.0013, 1.8955]
                b = [0.2459, 1.8542, 6.4411, 22.1730, 59.2206]
            case(53) ! I
                a = [0.7047, 1.9484, 2.5940, 4.1526, 1.5057]
                b = [0.2455, 1.8638, 6.7639, 21.8007, 56.4395]
            case(54) ! XE
                a = [0.6737, 1.7908, 2.4129, 4.2100, 1.7058]
                b = [0.2305, 1.6890, 5.8218, 18.3928, 47.2496]
            case(55) ! CS
                a = [1.2704, 3.8018, 5.6618, 0.9205, 4.8105]
                b = [0.4356, 4.2058, 23.4342, 136.778393, 171.7561]
            case(56) ! BA
                a = [0.9049, 2.6076, 4.8498, 5.1603, 4.7388]
                b = [0.3066, 2.4363, 12.1821, 54.6135, 161.9978]
            case(57) ! LA
                a = [0.8405, 2.3863, 4.6139, 5.1514, 4.7949]
                b = [0.2791, 2.1410, 10.3400, 41.9148, 132.0204]
            case(58) ! CE
                a = [0.8551, 2.3915, 4.5772, 5.0278, 4.5118]
                b = [0.2805, 2.1200, 10.1808, 42.0633, 130.9893]
            case(59) ! PR
                a = [0.9096, 2.5313, 4.5266, 4.6376, 4.3690]
                b = [0.2939, 2.2471, 10.8266, 48.8842, 147.6020]
            case(60) ! ND
                a = [0.8807, 2.4183, 4.4448, 4.6858, 4.1725]
                b = [0.2802, 2.0836, 10.0357, 47.4506, 146.9976]
            case(61) ! PM
                a = [0.9471, 2.5463, 4.3523, 4.4789, 3.9080]
                b = [0.2977, 2.2276, 10.5762, 49.3619, 145.3580]
            case(62) ! SM
                a = [0.9699, 2.5837, 4.2778, 4.4575, 3.5985]
                b = [0.3003, 2.2447, 10.6487, 50.7994, 146.4179]
            case(63) ! EU
                a = [0.8694, 2.2413, 3.9196, 3.9694, 4.5498]
                b = [0.2653, 1.8590, 8.3998, 36.7397, 125.7089]
            case(64) ! GD
                a = [0.9673, 2.4702, 4.1148, 4.4972, 3.2099]
                b = [0.2909, 2.1014, 9.7067, 43.4270, 125.9474]
            case(65) ! TB
                a = [0.9325, 2.3673, 3.8791, 3.9674, 3.7996]
                b = [0.2761, 1.9511, 8.9296, 41.5937, 131.0122]
            case(66) ! DY
                a = [0.9505, 2.3705, 3.8218, 4.0471, 3.4451]
                b = [0.2773, 1.9469, 8.8862, 43.0938, 133.1396]
            case(67) ! HO
                a = [0.9248, 2.2428, 3.6182, 3.7910, 3.7912]
                b = [0.2660, 1.8183, 7.9655, 33.1129, 101.8139]
            case(68) ! ER
                a = [1.0373, 2.4824, 3.6558, 3.8925, 3.0056]
                b = [0.2944, 2.0797, 9.4156, 45.8056, 132.7720]
            case(69) ! TM
                a = [1.0075, 2.3787, 3.5440, 3.6932, 3.1759]
                b = [0.2816, 1.9486, 8.7162, 41.8420, 125.0320]
            case(70) ! YB
                a = [1.0347, 2.3911, 3.4619, 3.6556, 3.0052]
                b = [0.2855, 1.9679, 8.7619, 42.3304, 125.6499]
            case(71) ! LU
                a = [0.9927, 2.2436, 3.3554, 3.7813, 3.0994]
                b = [0.2701, 1.8073, 7.8112, 34.4849, 103.3526]
            case(72) ! HF
                a = [1.0295, 2.2911, 3.4110, 3.9497, 2.4925]
                b = [0.2761, 1.8625, 8.0961, 34.2712, 98.5295]
            case(73) ! TA
                a = [1.0190, 2.2291, 3.4097, 3.9252, 2.2679]
                b = [0.2694, 1.7962, 7.6944, 31.0942, 91.1089]
            case(74) ! W
                a = [0.9853, 2.1167, 3.3570, 3.7981, 2.2798]
                b = [0.2569, 1.6745, 7.0098, 26.9234, 81.3910]
            case(75) ! RE
                a = [0.99141, 2.0858, 3.4531, 3.38812, 1.8526]
                b = [0.2548, 1.6518, 6.8845, 26.7234, 81.7215]
            case(76) ! OS
                a = [0.9813, 2.0322, 3.3665, 3.6235, 1.9741]
                b = [0.2487, 1.5973, 6.4737, 23.2817, 70.9254]
            case(77) ! IR
                a = [1.0194, 2.0645, 3.4425, 3.4914, 1.6976]
                b = [0.2554, 1.6475, 6.5966, 23.2269, 70.0272]
            case(78) ! Pt
                a = [0.9148, 1.8096, 3.2134,  3.2953,  1.5754]
                b = [0.2263, 1.3813, 5.3243, 17.5987, 60.0171]
            case(79) ! Au
                a = [0.9674, 1.8916, 3.3993,  3.0524,  1.2607]
                b = [0.2358, 1.4712, 5.6758, 18.7119, 61.5286]
            case(80) ! HG
                a = [1.0033, 1.9469, 3.4396,  3.1548,  1.4180]
                b = [0.2413, 1.5298, 5.8009, 19.4520, 60.5753]
            case(81) ! TL
                a = [1.0689, 2.1038, 3.6039,  3.4927,  1.8283]
                b = [0.2540, 1.6715, 6.3509, 23.1531, 78.7099]
            case(82) ! PB
                a = [1.0891, 2.1867, 3.6160, 3.8031, 1.8994]
                b = [0.2552, 1.7174, 6.5131, 23.9170, 74.7039]
            case(83) ! BI
                a = [1.1007, 2.2306, 3.5689, 4.1549, 2.0382]
                b = [0.2546, 1.7351, 6.4948, 23.6464, 70.3780]
            case(84) ! PO
                a = [1.1568, 2.4353, 3.6459, 4.4064, 1.7179]
                b = [0.2648, 1.8786, 7.1749, 25.1766, 69.2821]
            case(85) ! AT
                a = [1.0909, 2.1976, 3.3831, 4.6700, 2.1277]
                b = [0.2466, 1.6707, 6.0197, 20.7657, 57.2663]
            case(86) ! RN
                a = [1.0756, 2.1630, 3.3178, 4.8852, 2.0489]
                b = [0.2402, 1.6169, 5.7644, 19.4568, 52.5009]
            case(87) ! FR
                a = [1.4282, 3.5081, 5.6767, 4.1964, 3.8946]
                b = [0.3183, 2.6889, 13.4816, 54.3866, 200.8321]
            case(88) ! RA
                a = [1.3127, 3.1243, 5.2988, 5.3891, 5.4133]
                b = [0.2887, 2.2897, 10.8276, 43.5389, 145.6109]
            case(89) ! AC
                a = [1.3128, 3.1021, 5.3385, 5.9611, 4.7562]
                b = [0.2861, 2.2509, 10.5287, 41.7796, 128.2973]
            case(90) ! TH
                a = [1.2553, 2.9178, 5.0862, 6.1206, 4.7122]
                b = [0.2701, 2.0636, 9.3051, 34.5977, 107.9200]
            case(91) ! PA
                a = [1.3218, 3.1444, 5.4371, 5.6444, 4.0107]
                b = [0.2827, 2.2250, 10.2454, 41.1162, 124.4449]
            case(92) ! U
                a = [1.3382, 3.2043, 5.4558, 5.4839, 3.6342]
                b = [0.2838, 2.2452, 10.2519, 41.7251, 124.9023]
            case(93) ! NP
                a = [1.5193, 4.0053, 6.5327, -0.1402, 6.7489]
                b = [0.3213, 2.8206, 14.8878, 68.9103, 81.7257]
            case(94) ! PU
                a = [1.3517, 3.2937, 5.3213, 4.6466, 3.5714]
                b = [0.2813, 2.2418, 9.9952, 42.7939, 132.1739]
            case(95) ! AM
                a = [1.2135, 2.7962, 4.7545, 4.5731, 4.4786]
                b = [0.2483, 1.8437, 7.5421, 29.3841, 112.4579]
            case(96) ! CM
                a = [1.2937, 3.1100, 5.0393, 4.7546, 3.5031]
                b = [0.2638, 2.0341, 8.7101, 35.2992, 109.4972]
            case DEFAULT
                a = 0._dp; b = 1._dp
            end select
    end subroutine reference_coefficients

end module simple_cif2mrc_tester
