!@descr: library test of reconstruct3D_wave: a multislice-imaged phantom on a curvature-limited geometry
! A phantom of Gaussian atoms 100 A across, box 128 at 1 A (Nyquist 2 A), imaged by simulate_scattering's
! multislice at 100 kV with Ewald curvature negative (the same lattice and lens as the reconstruction, but
! another model: Rytov reconstructs data multislice made), with Gaussian
! detector noise of 0.002 of the vacuum intensity, converted to SIMPLE's particles (contrast inverted, unit
! noise: x = -(I - 1) / 0.002, the gain 500 of the scale link, given as the dose (500 / smpd)^2):
! sqrt(2 lambda D) = 2.7 A, so the band from 2.7 A to Nyquist is where the flat model loses the curvature.
! reconstruct3D_wave with wave_model=rytov at ewald=neg must beat, in that band (the FSC of the half-map mean
! against the truth, in a soft sphere), both its own flat start (the PCG map under the same prior and support)
! and the solve with the wrong curvature sign (plan section 5.1, library tier).
module simple_rec3D_wave_strategy_tester
use simple_defs,           only: PI
use simple_type_defs,      only: ctfparams, CTFFLAG_YES
use simple_string,         only: string
use simple_syslib,         only: simple_mkdir, simple_chdir
use simple_fileio,         only: simple_copy_file
use simple_image,          only: image
use simple_ori,            only: ori
use simple_oris,           only: oris
use simple_cmdline,        only: cmdline
use simple_sp_project,     only: sp_project
use simple_commanders_sim, only: commander_simulate_scattering
use simple_commanders_rec, only: commander_rec3D_wave
use simple_test_utils,     only: assert_true, fix_test_seed, enter_fixture, leave_fixture, tests_failed
implicit none
private
public :: run_all_rec3D_wave_lib_tests

integer, parameter :: BOX = 128, NPTCLS = 300, NATOMS = 150, NTHR = 4
real,    parameter :: SMPD = 1.0, KV = 100., CS = 2.7, FRACA = 0.07, MSKDIAM = 110., NOISE = 0.002
real,    parameter :: ATOM_PEAK = 15., ATOM_SDEV = 0.7, RADIUS = 50.   !< volts, A, A
real,    parameter :: BAND_LO = 2.72                                   !< sqrt(2 lambda D), A

contains

    subroutine run_all_rec3D_wave_lib_tests()
        type(string) :: cwd_saved, root
        real :: fsc_flat, fsc_neg, fsc_pos
        integer :: nfail0
        write(*,'(A)') '**** running reconstruct3D_wave library tests ****'
        nfail0 = tests_failed
        call fix_test_seed()
        call enter_fixture('rec3D_wave', cwd_saved, root)
        call make_data()
        call reconstruct('neg', fsc_neg, fsc_flat)
        call reconstruct('pos', fsc_pos)
        write(*,'(A,3F8.4)') '    truth FSC beyond the curvature limit: flat start, Rytov ewald=neg, ewald=pos: ', &
            &fsc_flat, fsc_neg, fsc_pos
        call assert_true(fsc_neg > fsc_flat, 'Rytov with the curvature of the data beats the flat map beyond sqrt(2 lambda D)')
        call assert_true(fsc_neg > fsc_pos, 'the curvature sign of the data beats the opposite sign')
        call simple_chdir(root)
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine run_all_rec3D_wave_lib_tests

    !> the phantom, its multislice images and the project of the particles with their truth
    subroutine make_data()
        type(commander_simulate_scattering) :: xsim
        type(cmdline)     :: cline
        type(image)       :: vol
        class(oris),       allocatable :: os
        class(sp_project), allocatable :: project
        type(ori)         :: o
        type(ctfparams)   :: ctf
        type(cmdline)     :: cline_compenv
        real, allocatable :: v(:,:,:)
        real    :: c(3,NATOMS), d2, r(5)
        integer :: i, j, k, a
        ! Gaussian atoms within the sphere
        a = 0
        do while( a < NATOMS )
            call random_number(c(:,a+1))
            c(:,a+1) = (c(:,a+1) - 0.5) * 2. * RADIUS
            if( sum(c(:,a+1)**2) <= RADIUS**2 ) a = a + 1
        end do
        allocate(v(BOX,BOX,BOX), source=0.)
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    do a = 1, NATOMS
                        d2 = (SMPD * real(i - BOX/2 - 1) - c(1,a))**2 + (SMPD * real(j - BOX/2 - 1) - c(2,a))**2 + &
                            &(SMPD * real(k - BOX/2 - 1) - c(3,a))**2
                        v(i,j,k) = v(i,j,k) + ATOM_PEAK * exp(-d2 / (2. * ATOM_SDEV**2))
                    end do
                end do
            end do
        end do
        call vol%new([BOX,BOX,BOX], SMPD)
        call vol%set_rmat(v, .false.)
        call vol%write(string('truth.mrc'))
        call vol%kill
        ! random poses, no shifts, astigmatic underfocus 0.6 to 1.2 um
        allocate(os)
        call os%new(NPTCLS, is_ptcl=.true.)
        call o%new(.true.)
        do i = 1, NPTCLS
            call random_number(r)
            call o%set_euler([360.*r(1), acos(2.*r(2) - 1.) * 180. / PI, 360.*r(3)])
            call o%set_shift([0.,0.])
            call o%set('dfx',    0.6 + 0.6 * r(4))
            call o%set('dfy',    0.95 * (0.6 + 0.6 * r(4)))
            call o%set('angast', 180. * r(5))
            call os%set_ori(i, o)
        end do
        call o%kill
        call os%write(string('truth_oris.txt'), [1,NPTCLS])
        ! noise-free multislice images
        call cline%set('prg',        'simulate_scattering')
        call cline%set('vol1',       'truth.mrc')
        call cline%set('smpd',       SMPD)
        call cline%set('nptcls',     NPTCLS)
        call cline%set('oritab',     'truth_oris.txt')
        call cline%set('wave_model', 'multislice')
        call cline%set('ewald',      'neg')
        call cline%set('kv',         KV)
        call cline%set('cs',         CS)
        call cline%set('fraca',      FRACA)
        call cline%set('bfac',       0.)
        call cline%set('total_dose', 0.)
        call cline%set('noise_sdev', NOISE)
        call cline%set('outstk',     'ptcls.mrcs')
        call cline%set('outfile',    'sim_oris.txt')
        call cline%set('mkdir',      'no')
        call cline%set('nthr',       NTHR)
        call xsim%execute(cline)
        call cline%kill
        ! SIMPLE's particles: contrast inverted, unit noise
        call vol%new([BOX,BOX,1], SMPD)
        do i = 1, NPTCLS
            call vol%read(string('ptcls.mrcs'), i)
            v = vol%get_rmat()
            v = -(v - 1.) / NOISE
            call vol%set_rmat(v, .false.)
            call vol%write(string('particles.mrcs'), i)
        end do
        call vol%kill
        ! the particles' project: truth poses and CTF, alternating halves
        allocate(project)
        ctf%ctfflag = CTFFLAG_YES
        ctf%smpd    = SMPD
        ctf%kv      = KV
        ctf%cs      = CS
        ctf%fraca   = FRACA
        ctf%dfx     = 1.0
        ctf%dfy     = 1.0
        ctf%angast  = 0.0
        ctf%phshift = 0.0
        call project%add_stk(string('particles.mrcs'), ctf)
        do i = 1, NPTCLS
            call project%os_ptcl3D%set_euler(i, os%get_euler(i))
            call project%os_ptcl3D%set_shift(i, [0.,0.])
            call project%os_ptcl3D%set(i, 'dfx',    os%get(i, 'dfx'))
            call project%os_ptcl3D%set(i, 'dfy',    os%get(i, 'dfy'))
            call project%os_ptcl3D%set(i, 'angast', os%get(i, 'angast'))
            call project%os_ptcl3D%set_state(i, 1)
            call project%os_ptcl3D%set(i, 'eo', real(mod(i,2)))
        enddo
        call project%update_projinfo(string('wave.simple'))
        call project%update_compenv(cline_compenv)
        call project%write(string('wave.simple'))
        call project%kill
        call os%kill
        deallocate(project, os)
    end subroutine make_data

    !> reconstruct3D_wave (Rytov, ewald) in its own directory; fsc the band mean of the truth FSC of the half-map
    !! mean, fsc_flat that of the flat start
    subroutine reconstruct( ewald, fsc, fsc_flat )
        character(len=*), intent(in)  :: ewald
        real,             intent(out) :: fsc
        real, optional,   intent(out) :: fsc_flat
        type(commander_rec3D_wave) :: xwave
        type(cmdline) :: cline
        type(string)  :: dir
        integer       :: status
        dir = 'rytov_'//trim(ewald)
        call simple_mkdir(dir)
        call simple_copy_file(string('wave.simple'), dir//'/wave.simple', status)
        call simple_chdir(dir)
        call cline%set('prg',        'reconstruct3D_wave')
        call cline%set('projfile',   'wave.simple')
        call cline%set('mkdir',      'no')
        call cline%set('wave_model', 'rytov')
        call cline%set('ewald',      trim(ewald))
        call cline%set('total_dose', (1. / (NOISE * SMPD))**2)
        call cline%set('objfun',     'cc')
        call cline%set('mskdiam',    MSKDIAM)
        call cline%set('pgrp',       'c1')
        call cline%set('maxits',     4)
        call cline%set('nthr',       NTHR)
        call xwave%execute(cline)
        call cline%kill
        fsc = band_fsc('recvol_state01_even.mrc', 'recvol_state01_odd.mrc')
        if( present(fsc_flat) ) fsc_flat = band_fsc('recvol_state01_even_flat.mrc', 'recvol_state01_odd_flat.mrc')
        call simple_chdir(string('..'))
        call dir%kill
    end subroutine reconstruct

    !> the mean truth FSC of the half-map mean over the shells from BAND_LO to Nyquist, in a soft sphere
    real function band_fsc( even, odd )
        character(len=*), intent(in) :: even, odd
        type(image) :: e, o, t
        real, allocatable :: corrs(:)
        integer :: k, n, klo
        call e%new([BOX,BOX,BOX], SMPD)
        call o%new([BOX,BOX,BOX], SMPD)
        call t%new([BOX,BOX,BOX], SMPD)
        call e%read(string(even))
        call o%read(string(odd))
        call t%read(string('../truth.mrc'))
        call e%add(o)
        call e%mask3D_soft((MSKDIAM / 2. + 4.) / SMPD)
        call t%mask3D_soft((MSKDIAM / 2. + 4.) / SMPD)
        call e%fft
        call t%fft
        n = e%get_filtsz()
        allocate(corrs(n))
        call e%fsc(t, corrs)
        klo = ceiling(real(BOX) * SMPD / BAND_LO)
        band_fsc = 0.
        do k = klo, n
            band_fsc = band_fsc + corrs(k)
        end do
        band_fsc = band_fsc / real(n - klo + 1)
        call e%kill
        call o%kill
        call t%kill
    end function band_fsc

end module simple_rec3D_wave_strategy_tester
