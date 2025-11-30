! Tests anisotropic and isotropic bfactor calculations in atoms_stats
! If all tests pass, then the fitting routines accurately fit ideal isotropic
! and anisotropic trivariate Gaussians.  There could still be errors
! in other aspects of bfactor calculations. 
! Requires an input reconstructed volume, connected component binary volume
! and a PDB file.
! Note that if a test fails, that doesn't mean there's necessarily a bug:
! it could mean that there's a poor fit for one atom.  Better is to analyze
! the errors in the two results output files.
! For assembling the tests in github, using NP7 and NP1 from the Homogeneous
! NP Paper should provide a good range of possible atomic densities.
program simple_test_bfactor
!     include 'simple_lib.f08'
!     use simple_cmdline,          only: cmdline
!     use simple_commanders_atoms, only: commander_atoms_stats
!     use simple_atoms,            only: atoms
!     implicit none

!     ! Output atomic statistics in atoms_stats.csv file (must be in order).
!     ! Used to easily read atoms_stats.csv file line by line.
!     ! Must be kept up-to-date for this test program to work.
!     type :: stats
!         real :: index, nvox, cnstd, nnbondl, cngen, diam, adj_cn13, avgint, maxint
!         real :: cendist, vcorr 
!         real :: uiso, uaniso(3), azimuth, polar ! Relevant
!         real :: doi, doi_min
!         real :: iso_corr, aniso_corr ! Relevant
!         real :: rstrain
!     end type

!     type (cmdline)              :: cline_atoms_stats
!     type(atoms)                 :: inatoms
!     type(commander_atoms_stats) :: xatoms_stats
!     type(stats), allocatable    :: ideal(:), fitiso(:), fitaniso(:)
!     real, parameter             :: SMPD=0.358, MSKDIAM=40., NTHR=16., CORR_THRES=0.95
!     real                        :: mean_abs_error(5), mean_pct_error(5)
!     integer                     :: i, natoms, naniso, fu_ideal, fu_iso, fu_aniso, cc, io, ntests_passed
!     integer                     :: fu_run1, fu_run2
!     logical                     :: test_passed
!     character(len=256)          :: invol, incc, inpdb
!     character(len=2)            :: element
!     character(*), parameter     :: dir='simple_test_bfactor_out'
!     character(*), parameter     :: isodir='test_isotropic_bfactor'
!     character(*), parameter     :: anisodir='test_anisotropic_bfactor'
!     character(*), parameter     :: fn_stats='atoms_stats.csv'
!     character(*), parameter     :: fn_isofit='fit_isotropic.mrc'
!     character(*), parameter     :: fn_anisofit='fit_anisotropic.mrc'
!     character(*), parameter     :: fn_run1='results_simple_test_bfactor_iso_input.csv'
!     character(*), parameter     :: fn_run2='results_simple_test_bfactor_aniso_input.csv'
!     character(*), parameter :: HEADRUN1 = 'INDEX'//CSV_DELIM//'U_ISO_IN'//CSV_DELIM&
!         &//'U_ISO_OUT'//CSV_DELIM//'U_MAJ_OUT'//CSV_DELIM//'U_MED_OUT'//CSV_DELIM//&
!         &'U_MIN_OUT'//CSV_DELIM//'AZIMUTH_OUT'//CSV_DELIM//'POLAR_OUT'//CSV_DELIM//&
!         &'ISO_CORR'//CSV_DELIM//'ANISO_CORR'
!     character(*), parameter :: HEADRUN2 = 'INDEX'//CSV_DELIM//'U_MAJ_IN'//&
!         &CSV_DELIM//'U_MED_IN'//CSV_DELIM//'U_MIN_IN'//CSV_DELIM//'AZIMUTH_IN'//&
!         &CSV_DELIM//'POLAR_IN'//CSV_DELIM//'U_ISO_OUT'//CSV_DELIM//'U_MAJ_OUT'//CSV_DELIM//&
!         &'U_MED_OUT'//CSV_DELIM//'U_MIN_OUT'//CSV_DELIM//'AZIMUTH_OUT'//CSV_DELIM//&
!         &'POLAR_OUT'//CSV_DELIM//'ISO_CORR'//CSV_DELIM//'ANISO_CORR'


!     if( command_argument_count() /= 3 )then
!         write(logfhandle,'(a)') 'Usage: simple_test_bfactor vol.mrc cc.mrc &
!                                 &file.pdb'
!         write(logfhandle,'(a)') 'Note: Paths must be relative'
!         stop
!     else
!         call get_command_argument(1, invol)
!         call get_command_argument(2, incc)
!         call get_command_argument(3, inpdb)
!     endif

!     call inatoms%new(inpdb)
!     element = inatoms%get_element(1)
!     natoms = inatoms%get_n()

!     ! Generate ideal isotropic and anisotropic Gaussian fits from atoms stats
!     write(logfhandle,'(a)') '>>> PASSING INPUT FILES TO ATOMS_STATS'
!     call simple_mkdir(dir)
!     call simple_chdir(dir)
!     call cline_atoms_stats%set('vol1', '../'//invol)
!     call cline_atoms_stats%set('vol2', '../'//incc)
!     call cline_atoms_stats%set('pdbfile', '../'//inpdb)
!     call cline_atoms_stats%set('smpd', SMPD)
!     call cline_atoms_stats%set('element', element)
!     call cline_atoms_stats%set('mskdiam', MSKDIAM)
!     call cline_atoms_stats%set('nthr', NTHR)
!     call xatoms_stats%execute(cline_atoms_stats)

!     ! Run 1: Use a density with ideal isotropic 3D Gaussians
!     write(logfhandle,'(a)') '>>> Run 1: IDEAL ISOTROPIC 3D GAUSSIANS'
!     call simple_mkdir(isodir)
!     call simple_chdir(isodir)
!     call cline_atoms_stats%set('vol1', '../'//fn_isofit)
!     call cline_atoms_stats%set('vol2', '../../'//incc)
!     call cline_atoms_stats%set('pdbfile', '../../'//inpdb)
!     call cline_atoms_stats%set('smpd', SMPD)
!     call cline_atoms_stats%set('element', element)
!     call cline_atoms_stats%set('mskdiam', MSKDIAM)
!     call cline_atoms_stats%set('nthr', NTHR)
!     call xatoms_stats%execute(cline_atoms_stats)
!     call simple_chdir('../')

!     ! Run 2: Use a density with ideal anisotropic 3D Gaussians
!     write(logfhandle,'(a)') '>>> Run 2: IDEAL ANISOTROPIC 3D GAUSSIANS'
!     call simple_mkdir(anisodir)
!     call simple_chdir(anisodir)
!     call cline_atoms_stats%set('vol1', '../'//fn_anisofit)
!     call cline_atoms_stats%set('vol2', '../../'//incc)
!     call cline_atoms_stats%set('pdbfile', '../../'//inpdb)
!     call cline_atoms_stats%set('smpd', SMPD)
!     call cline_atoms_stats%set('element', element)
!     call cline_atoms_stats%set('mskdiam', MSKDIAM)
!     call cline_atoms_stats%set('nthr', NTHR)
!     call xatoms_stats%execute(cline_atoms_stats)
!     call simple_chdir('../')

!     ! Write input/output bfactor stats in two .csv files for in-depth error analysis
!     allocate(ideal(natoms), fitiso(natoms), fitaniso(natoms))
!     call fopen(fu_ideal, FILE=trim(fn_stats), STATUS='OLD', action='READ')
!     call fopen(fu_iso, FILE=trim(isodir//'/'//fn_stats), STATUS='OLD', action='READ')
!     call fopen(fu_aniso, FILE=trim(anisodir//'/'//fn_stats), STATUS='OLD', action='READ')
!     call fopen(fu_run1, FILE=trim(fn_run1), STATUS='REPLACE', action='WRITE', iostat=io)
!     call fopen(fu_run2, FILE=trim(fn_run2), STATUS='REPLACE', action='WRITE')
!     read (fu_ideal, *)      ! Skip header
!     read (fu_iso, *)        ! Skip header
!     read (fu_aniso, *)      ! Skip header
!     write (fu_run1, '(A)') HEADRUN1
!     write (fu_run2, '(A)') HEADRUN2
!     do cc=1,natoms
!         read (fu_ideal, *) ideal(cc)
!         read (fu_iso, *) fitiso(cc)
!         read(fu_aniso, *) fitaniso(cc)
!         ! Write results for both tests
!         call write_test(cc, fu_run1, ideal, fitiso, .true.)
!         call write_test(cc, fu_run2, ideal, fitaniso, .false.)
!     end do
!     call fclose(fu_ideal)
!     call fclose(fu_iso)
!     call fclose(fu_aniso)
!     call fclose(fu_run1)
!     call fclose(fu_run2)

!     ! AUTOMATIC TESTS THAT ALL CORRELATIONS OF INPUT AND FITS >= CORR_THRES
!     ! REPORT MEAN ABSOLUTE AND PERCENT ERROR OF FITS
!     ntests_passed = 0
!     ! TEST 1
!     write(logfhandle,'(a)') '>>> TEST 1: CHECKING ISOTROPIC FITS OF ISOTROPIC INPUT'
!     test_passed = .true.
!     mean_abs_error = 0.
!     mean_pct_error = 0.
!     do cc=1, natoms
!         if (fitiso(cc)%iso_corr < CORR_THRES) then
!             write(logfhandle,'(a,i5,a,f8.4)') 'ERROR AT CC = ', cc, ', ISO_CORR = ',&
!                     &fitiso(cc)%iso_corr
!             test_passed = .false.
!         end if
!         ! Indices of mean_abs_error used later for 3 anisotropic eigenvals and 2 anisotropic angles
!         mean_abs_error(1) = mean_abs_error(1) + abs(ideal(cc)%uiso-fitiso(cc)%uiso)
!         mean_pct_error(1) = mean_pct_error(1) + abs(ideal(cc)%uiso-fitiso(cc)%uiso)/ideal(cc)%uiso
!     end do
!     mean_abs_error = mean_abs_error / natoms
!     mean_pct_error = mean_pct_error / natoms
!     write(logfhandle,'(a, f12.8)') 'mean U_ISO error (Å): ', mean_abs_error(1)
!     write(logfhandle,'(a, f12.8)') 'mean U_ISO (%): ', mean_pct_error(1)*100
!     if (test_passed) then
!         write(logfhandle,'(a)') 'TEST 1 PASSED'
!         ntests_passed = ntests_passed + 1
!     else
!         write(logfhandle,'(a)') 'TEST 1 FAILED'
!     end if

!     ! TEST 2
!     write(logfhandle,'(a)') '>>> TEST 2: CHECKING ANISOTROPIC FITS OF ISOTROPIC INPUT'
!     test_passed = .true.
!     mean_abs_error = 0.
!     mean_pct_error = 0.
!     do cc=1, natoms
!         if (fitiso(cc)%aniso_corr < CORR_THRES) then
!             write(logfhandle,'(a,i5,a,f8.4)') 'ERROR AT CC = ', cc, ', ANISO_CORR = ',&
!                     &fitiso(cc)%aniso_corr
!             test_passed = .false.
!         end if
!         ! Error of 3 eigenvalues
!         do i=1,3
!             mean_abs_error(i) = mean_abs_error(i) + abs(ideal(cc)%uiso-fitiso(cc)%uaniso(i))
!             mean_pct_error(i) = mean_pct_error(i) + abs(ideal(cc)%uiso-fitiso(cc)%uaniso(i))/ideal(cc)%uiso
!         end do
!     end do
!     mean_abs_error = mean_abs_error / natoms
!     mean_pct_error = mean_pct_error / natoms
!     write(logfhandle,'(a, f12.8)') 'U_MAJ mean error (Å): ', mean_abs_error(1)
!     write(logfhandle,'(a, f12.8)') 'U_MAJ mean error (%): ', mean_pct_error(1)*100
!     write(logfhandle,'(a, f12.8)') 'U_MED mean error (Å): ', mean_abs_error(2)
!     write(logfhandle,'(a, f12.8)') 'U_MED mean error (%): ', mean_pct_error(2)*100
!     write(logfhandle,'(a, f12.8)') 'U_MIN mean error (Å): ', mean_abs_error(3)
!     write(logfhandle,'(a, f12.8)') 'U_MIN mean error (%): ', mean_pct_error(3)*100
!     if (test_passed) then
!         write(logfhandle,'(a)') 'TEST 2 PASSED'
!         ntests_passed = ntests_passed + 1
!     else
!         write(logfhandle,'(a)') 'TEST 2 FAILED'
!     end if

!     ! TEST 3
!     write(logfhandle,'(a)') '>>> TEST 3: CHECKING ANISOTROPIC FIT EIGENVALUES OF ANISOTROPIC INPUT'
!     test_passed = .true.
!     mean_abs_error = 0.
!     mean_pct_error = 0.
!     naniso = 0
!     do cc=1, natoms
!         ! Only consider atoms that didn't have aniso fits tossed.
!         if (fitaniso(cc)%uaniso(1) > TINY) then
!             naniso = naniso + 1
!             if (fitaniso(cc)%aniso_corr < CORR_THRES) then
!                 write(logfhandle,'(a,i5,a,f8.4)') 'ERROR AT CC = ', cc, ', ANISO_CORR = ',&
!                         &fitiso(cc)%aniso_corr
!                 test_passed = .false.
!             end if
!             ! Error of 3 eigenvalues
!             do i=1,3
!                 mean_abs_error(i) = mean_abs_error(i) + abs(ideal(cc)%uaniso(i)-fitaniso(cc)%uaniso(i))
!                 mean_pct_error(i) = mean_pct_error(i) + abs(ideal(cc)%uaniso(i)-fitaniso(cc)%uaniso(i))/ideal(cc)%uaniso(i)
!             end do
!             ! Error of 2 orientation angles
!             mean_abs_error(4) = mean_abs_error(4) + abs(ideal(cc)%azimuth-fitaniso(cc)%azimuth)
!             mean_abs_error(5) = mean_abs_error(5) + abs(ideal(cc)%polar-fitaniso(cc)%polar)
!         end if
!     end do
!     mean_abs_error = mean_abs_error / naniso
!     mean_pct_error = mean_pct_error / naniso
!     write(logfhandle,'(a, f12.8)') 'U_MAJ mean error (Å): ', mean_abs_error(1)
!     write(logfhandle,'(a, f12.8)') 'U_MAJ mean error (%): ', mean_pct_error(1)*100
!     write(logfhandle,'(a, f12.8)') 'U_MED mean error (Å): ', mean_abs_error(2)
!     write(logfhandle,'(a, f12.8)') 'U_MED mean error (%): ', mean_pct_error(2)*100
!     write(logfhandle,'(a, f12.8)') 'U_MIN mean error (Å): ', mean_abs_error(3)
!     write(logfhandle,'(a, f12.8)') 'U_MIN mean error (%): ', mean_pct_error(3)*100
!     write(logfhandle,'(a, f12.8)') 'AZIMUTH mean error (rad): ', mean_abs_error(4)
!     write(logfhandle,'(a, f12.8)') 'POLAR mean error (rad): ', mean_abs_error(5)
!     if (test_passed) then
!         write(logfhandle,'(a)') 'TEST 3 PASSED'
!         ntests_passed = ntests_passed + 1
!     else
!         write(logfhandle,'(a)') 'TEST 3 FAILED'
!     end if

!     ! Sanity check for isotropic fits of anisotropic input: ensure U_MIN <= U_ISO <= U_MAJ
!     write(logfhandle,'(a)') '>>> TEST 4: CHECKING REASONABLE ISOTROPIC FITS OF ANISOTROPIC INPUT'
!     do cc=1, natoms
!         if (fitaniso(cc)%uaniso(1) > TINY) then
!             if (fitaniso(cc)%uiso < fitaniso(cc)%uaniso(3) .or. fitaniso(cc)%uiso > fitaniso(cc)%uaniso(1)) then
!                 write(logfhandle,'(a,i5,a,f8.4,a,f8.4,a,f8.4)') 'ERROR AT CC = ', cc, &
!                         &', U_ISO = ', fitaniso(cc)%uiso, ', U_MAJ = ', fitaniso(cc)%uaniso(1), &
!                         &'U_MIN = ', fitaniso(cc)%uaniso(3)
!                 test_passed = .false.
!             end if
!         end if
!     end do
!     if (test_passed) then
!         write(logfhandle,'(a)') 'TEST 4 PASSED'
!         ntests_passed = ntests_passed + 1
!     else
!         write(logfhandle,'(a)') 'TEST 4 FAILED'
!     end if

!     deallocate(ideal, fitiso, fitaniso)
!     if (ntests_passed == 4) then 
!         write(logfhandle,'(a)') 'SUCCESS: 4/4 TESTS PASSED'
!     else
!         write(logfhandle,'(a,i1,a)') 'FAILURE: ', ntests_passed, '/4 TESTS PASSED'
!     endif
!     write(logfhandle, '(a)') '***END SIMPLE_TEST_BFACTOR***'

! contains

!     subroutine write_test(cc, fu, stats_in, stats_out, iso_in)
!         integer, intent(in) :: cc, fu
!         type(stats), allocatable, intent(in) :: stats_in(:), stats_out(:)
!         logical, intent(in) :: iso_in

!         601 format(F8.4,A2)
!         602 format(F8.4)
!         ! various per-atom parameters
!         write(fu,601,advance='no') real(cc),                    CSV_DELIM ! INDEX
!         if (iso_in) then
!             write(fu,601,advance='no') stats_in(cc)%uiso,       CSV_DELIM ! BFAC_IN
!         else
!             write(fu,601,advance='no') stats_in(cc)%uaniso(1),  CSV_DELIM ! MAJSAX_IN
!             write(fu,601,advance='no') stats_in(cc)%uaniso(2),  CSV_DELIM ! MEDSAX_IN
!             write(fu,601,advance='no') stats_in(cc)%uaniso(3),  CSV_DELIM ! MINSAX_IN
!             write(fu,601,advance='no') stats_in(cc)%azimuth,    CSV_DELIM ! AZMTH_IN
!             write(fu,601,advance='no') stats_in(cc)%polar,      CSV_DELIM ! POLAR_IN
!         end if
!         write(fu,601,advance='no') stats_out(cc)%uiso,          CSV_DELIM ! BFAC_OUT
!         write(fu,601,advance='no') stats_out(cc)%uaniso(1),     CSV_DELIM ! MAJSAX_OUT
!         write(fu,601,advance='no') stats_out(cc)%uaniso(2),     CSV_DELIM ! MEDSAX_OUT
!         write(fu,601,advance='no') stats_out(cc)%uaniso(3),     CSV_DELIM ! MINSAX_OUT
!         write(fu,601,advance='no') stats_out(cc)%azimuth,       CSV_DELIM ! AZIMUTH_OUT
!         write(fu,601,advance='no') stats_out(cc)%polar,         CSV_DELIM ! POLAR_OUT
!         write(fu,601,advance='no') stats_out(cc)%iso_corr,      CSV_DELIM ! ISO_CORR
!         write(fu,602)              stats_out(cc)%aniso_corr               ! ANISO_CORR
!     end subroutine write_test


end program simple_test_bfactor
