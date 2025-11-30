! This program tests the fitting of NP atomic positions with a crystal lattice
! Takes in a PDB file to generate an ideal lattice with potentially missing atoms.
! For a range of lattice parameters, fits an ideal lattice and a lattice with simulated
! displacements.  The test succeeds if the error in the mean fitted parameters
! is within the error thresholds (ideal_err_thres and displ_err_thres) for all input
! lattice parameters. Reports results in a .csv file.
program simple_test_fit_lattice

!     include 'simple_lib.f08'
!     use simple_atoms,              only: atoms
!     use simple_nanoparticle_utils, only: fit_lattice, find_rMax
!     implicit none

!     type(atoms)             :: atomsin
!     real, allocatable       :: displ(:,:), centers(:,:), ideal(:,:)
!     real, parameter         :: displ_mean=0., aStep=0.2, ideal_err_thres=0.01, displ_err_thres=0.1
!     real                    :: a0, a(3), a_fit(3), aMin, aMax, distTemp, dist, cMin(3)
!     real                    :: cMax(3), cMid(3), displ_std, limits(2)
!     integer, allocatable    :: lat_pos(:,:)
!     integer                 :: funit, natoms, cc, centerAtom
!     logical                 :: test_passed = .true.
!     character(len=2)        :: el_ucase
!     character(len=8)        :: crystal_system
!     character(len=256)      :: pdbin
!     character(*), parameter :: fn_results='results.csv'
!     character(*), parameter :: header='INPUT_LAT_PARAM'//CSV_DELIM//'DISPL_SDEV'//CSV_DELIM//&
!                                 &'FIT_LAT_PARAM'//CSV_DELIM//'TEST_SUCCESS'
    
!     if( command_argument_count() /= 1 ) then
!         write(logfhandle,'(a)') '>>> Usage: simple_test_fit_lattice file.pdb'
!         write(logfhandle,'(2a)') 'file.pdb contains the input atomic elements and coordinates.',&
!                 & NEW_LINE('a')
!         stop
!     else
!         call get_command_argument(1, pdbin)
!     endif

!     ! Read info from PDB (This allows us to test lattice fitting on lattices with gaps)
!     write(logfhandle,'(a)') '>>> FINDING IDEAL LATTICE'
!     call atomsin%new(string(trim(pdbin)))
!     natoms = atomsin%get_n()
!     write(logfhandle,'(a, i5, a)') '>>> Lattice contains ', natoms, ' atoms.'
!     allocate(centers(3, natoms), source=0.)
!     el_ucase = uppercase(trim(adjustl(atomsin%get_element(1))))
!     call get_lattice_params(el_ucase, crystal_system, a0)
!     do cc=1, natoms
!         centers(1:3, cc) = atomsin%get_coord(cc)
!     end do

!     ! Get center atom and generate integer lattice
!     centerAtom = 0
!     distTemp = HUGE(distTemp)
!     cMin       = minval(centers,dim=2)
!     cMax       = maxval(centers,dim=2)
!     cMid       = (cMax-cMin)/2.+cMin
!     do cc=1, natoms
!         dist = euclid(cMid, centers(1:3,cc))
!         if (dist < distTemp) then
!             centerAtom = cc
!             distTemp = dist
!         end if
!     end do
!     allocate(lat_pos(3, natoms), source=0)
!     do cc=1, natoms
!         lat_pos(1:3,cc) = nint((centers(1:3,cc) - centers(1:3,centerAtom)) / (a0 / 2.))
!     end do

!     allocate(ideal(3, natoms), displ(3, natoms), source=0.)
!     call fopen(funit, FILE=string(fn_results), STATUS='REPLACE', action='WRITE')
!     write(funit, '(a)') header
!     ! Lattice fitting assumes that the experimental lattice parameter a falls in the range
!     ! rMax < a < rMax * sqrt(2).
!     aMin = find_rMax(el_ucase)
!     aMax = sqrt(2.) * aMin
!     a = aMin + aStep
!     do while(a(1) < aMax)

!         ! Test 1: Ideal lattice
!         do cc=1, natoms
!             ideal(1:3,cc) = centers(centerAtom,1:3) + a/2.*lat_pos(1:3,cc)
!         end do
!         call fit_lattice(el_ucase, ideal, a_fit)
!         if (abs(sum(a_fit)/3. - sum(a)/3.) < ideal_err_thres) then
!             call write_result(funit, a, 0., a_fit, .true.)
!         else
!             call write_result(funit, a, 0., a_fit, .false.)
!             test_passed = .false.
!         end if

!         ! Test 2: Lattice with simulated random disiplacements
!         displ_std = a(1) / 20.
!         limits = [-2.*displ_std, 2*displ_std] ! So there aren't extreme outliers
!         do cc=1, natoms
!             displ(1:3,cc) = ideal(1:3,cc) + gasdev(displ_mean, displ_std, limits)
!         end do
!         call fit_lattice(el_ucase, displ, a_fit)
!         if (abs(sum(a_fit)/3. - sum(a)/3.) < displ_err_thres) then
!             call write_result(funit, a, displ_std, a_fit, .true.)
!         else
!             call write_result(funit, a, displ_std, a_fit, .false.)
!             test_passed = .false.
!         end if

!         a = a + aStep
!     end do

!     if (test_passed) then
!         write(logfhandle,'(a)') '***TEST SUCCESS***'
!     else
!         write(logfhandle,'(a)') '***TEST FAILURE (SEE RESULTS.CSV)***'
!     end if
!     write(logfhandle,'(a)') '***END SIMPLE_TEST_FIT_LATTICE***'

! contains

!     ! Writes the input lattice parameter a, standard deviation of input lattice 
!     ! displacements displ_std, the output lattice parameters a_fit, and whether 
!     ! the particular test was successful to the .csv file associated with funit.
!     subroutine write_result(funit, a, displ_std, a_fit, success)
!         integer, intent(in) :: funit
!         real, intent(in)    :: a(3), displ_std, a_fit(3)
!         logical, intent(in) :: success
!             write(funit, '(f8.4,a)', advance='no')  sum(a)/3., CSV_DELIM     ! INPUT_LAT_PARAM
!             write(funit, '(f8.4,a)', advance='no')  displ_std, CSV_DELIM     ! DISPL_SDEV
!             write(funit, '(f8.4,a)', advance='no')  sum(a_fit)/3., CSV_DELIM ! FIT_LAT_PARAM
!             if (success) then
!                 write(funit, '(f8.4,a)')            1., CSV_DELIM            ! FIT_LAT_PARAM
!             else
!                 write(funit, '(f8.4,a)')            0., CSV_DELIM            ! FIT_LAT_PARAM
!             end if
!     end subroutine write_result

end program simple_test_fit_lattice