module simple_defs_atoms
implicit none

contains

    ! single covalent radii from Cordero, et al., 2008, "Covalent radii revisited"
    ! Dalton Trans. (21): 2832–2838. doi:10.1039/b801115j
    subroutine get_element_Z_and_radius( element_ucase, Z, r )
        character(len=2), intent(in)  :: element_ucase
        integer,          intent(out) :: Z
        real,             intent(out) :: r
        Z = 0
        r = 1.
        select case(element_ucase)
            ! organic
            case('H')
                Z = 1 ; r = 0.31
            case('HE')
                Z = 2 ; r = 0.28
            case('LI')
                Z = 3 ; r = 1.28
            case('BE')
                Z = 4 ; r = 0.96
            case('B')
                Z = 5 ; r = 0.84
            case('C')
                Z = 6 ; r = 0.76 ! sp3
            case('N')
                Z = 7 ; r = 0.71
            case('O')
                Z = 8 ; r = 0.66
            case('F')
                Z = 9 ; r = 0.57
            case('NE')
                Z = 10; r = 0.58
            case('NA')
                Z = 11; r = 1.66
            case('MG')
                Z = 12; r = 1.41
            case('AL')
                Z = 13; r = 1.21
            case('SI')
                Z = 14; r = 1.11
            case('P')
                Z = 15; r = 1.07
            case('S')
                Z = 16; r = 1.05
            case('CL')
                Z = 17; r = 1.02
            case('AR')
                Z = 18; r = 1.06
            case('K')
                Z = 19; r = 2.03
            case('CA')
                Z = 20; r = 1.76
            case('SC')
                Z = 21; r = 1.70
            case('TI')
                Z = 22; r = 1.60
            case('V')
                Z = 23; r = 1.53
            case('CR')
                Z = 24; r = 1.39
            case('MN')
                Z = 25; r = 1.39
            case('FE')
                Z = 26; r = 1.02
                !Z = 26; r = 1.32
            case('CO')
                Z = 27; r = 1.26
            case('NI')
                Z = 28; r = 1.24
            case('CU')
                Z = 29; r = 1.32
            case('ZN')
                Z = 30; r = 1.22
            case('GA')
                Z = 31; r = 1.22
            case('GE')
                Z = 32; r = 1.20
            case('AS')
                Z = 33; r = 1.19
            case('SE')
                Z = 34; r = 1.20
            case('BR')
                Z = 35; r = 1.20
            case('KR')
                Z = 36; r = 1.16
            case('RB')
                Z = 37; r = 2.20
            case('SR')
                Z = 38; r = 1.95
            case('Y')
                Z = 39; r = 1.90
            case('ZR')
                Z = 40; r = 1.75
            case('NB')
                Z = 41; r = 1.64
            case('MO')
                Z = 42; r = 1.54
            case('TC')
                Z = 43; r = 1.47
            case('RU')
                Z = 44; r = 1.46
            case('RH')
                Z = 45; r = 1.42
            case('PD')
                Z = 46; r = 1.12
                !Z = 46; r = 1.39
            case('AG')
                Z = 47; r = 1.45
            case('CD')
                Z = 48; r = 1.44
            case('IN')
                Z = 49; r = 1.42
            case('SN')
                Z = 50; r = 1.39
            case('SB')
                Z = 51; r = 1.39
            case('TE')
                Z = 52; r = 1.38
            case('I')
                Z = 53; r = 1.39
            case('XE')
                Z = 54; r = 1.40
            case('CS')
                Z = 55; r = 2.44
            case('BA')
                Z = 56; r = 2.15
            case('LA')
                Z = 57; r = 2.07
            case('CE')
                Z = 58; r = 2.04
            case('PR')
                Z = 59; r = 2.03
            case('ND')
                Z = 60; r = 2.01
            case('PM')
                Z = 61; r = 1.99
            case('SM')
                Z = 62; r = 1.98
            case('EU')
                Z = 63; r = 1.98
            case('GD')
                Z = 64; r = 1.96
            case('TB')
                Z = 65; r = 1.94
            case('DY')
                Z = 66; r = 1.92
            case('HO')
                Z = 67; r = 1.92
            case('ER')
                Z = 68; r = 1.89
            case('TM')
                Z = 69; r = 1.90
            case('YB')
                Z = 70; r = 1.87
            case('LU')
                Z = 71; r = 1.87
            case('HF')
                Z = 72; r = 1.75
            case('TA')
                Z = 73; r = 1.70
            case('W')
                Z = 74; r = 1.62
            case('RE')
                Z = 75; r = 1.51
            case('OS')
                Z = 76; r = 1.44
            case('IR')
                Z = 77; r = 1.41
            case('PT')
                Z = 78; r = 1.10
                ! Z = 78; r = 1.36
                ! metallic radius below https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
                ! Z = 78; r = 1.385
            case('AU')
                Z = 79; r = 1.23
                ! Z = 79; r = 1.36
                ! metallic radius below https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
                ! Z = 79; r = 1.44
            case('HG')
                Z = 80; r = 1.32
            case('TL')
                Z = 81; r = 1.45
            case('PB')
                Z = 82; r = 1.46
            case('BI')
                Z = 83; r = 1.48
            case('PO')
                Z = 84; r = 1.40
            case('AT')
                Z = 85; r = 1.50
            case('RN')
                Z = 86; r = 1.50
            case('FR')
                Z = 87; r = 2.60
            case('RA')
                Z = 88; r = 2.21
            case('AC')
                Z = 89; r = 2.15
            case('TH')
                Z = 90; r = 2.06
            case('PA')
                Z = 91; r = 2.00
            case('U')
                Z = 92; r = 1.96
            case('NP')
                Z = 93; r = 1.90
            case('PU')
                Z = 94; r = 1.87
            case('AM')
                Z = 95; r = 1.80
            case('CM')
                Z = 96; r = 1.69
        end select
    end subroutine get_element_Z_and_radius

    ! From Wheeler, D, 1925, Physical Review. 25 (6): 753–761, FCC & BCC only.
    subroutine get_lattice_params( element_ucase, crystal_system, a )
        character(len=2), intent(in)    :: element_ucase
        character(len=8), intent(inout) :: crystal_system
        real,             intent(inout) :: a
        crystal_system = 'fcc     ' ! default
        select case( element_ucase )
            case('C')
                a = 3.567 ! diamond
            case('SI')
                a = 5.431020511
            case('GE')
                a = 5.658
            case('AL')
                a = 4.046
            case('NI')
                a = 3.499
            case('CU')
                a = 3.597
            case('PT')
                a = 3.912
            case('AU')
                a = 4.065
            case('AG')
                a = 4.079
            case('PD')
                a = 3.859
            case('PB')
                a = 4.920
            case('FE')
                a = 2.856; crystal_system = 'bcc     '
            case('MO')
                a = 3.142; crystal_system = 'bcc     '
            case('W')
                a = 3.155; crystal_system = 'bcc     '
            case('PBSE')
                a = 6.12;  crystal_system = 'rocksalt'
            case DEFAULT
                a = 3.76
        end select
    end subroutine get_lattice_params

end module simple_defs_atoms
