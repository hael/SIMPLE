!*******************************************************************************
!                                                                              !
! Name:                                                                        !
! simple_defs - basic definitions used in all modules.                         !
!                                                                              !
! Description:                                                                 !
! simple_defs provides basic definitions for the types and declarations        !
! and all the modules used through out the code. It is based on numerical      !
! nrtypes.f90 module.                                                          !
!*******************************************************************************
!
module simple_defs
use, intrinsic :: iso_c_binding
implicit none
character(len=1), parameter :: default_file_format = 'M' ! 'I', 'M' or 'S' for imagic, mrc, spider
integer, parameter  :: IMPORTANT=10 ! number of solutions considered important
integer, parameter  :: MAXS=99      ! maximum number of states
integer, parameter  :: STDLEN=256   ! standard string length
integer, parameter  :: short = selected_int_kind(4)
integer, parameter  :: long  = selected_int_kind(9)
integer, parameter  :: longer  = selected_int_kind(16)
integer, parameter  :: I4B = SELECTED_INT_KIND(9)
integer, parameter  :: I2B = SELECTED_INT_KIND(4)
integer, parameter  :: I1B = SELECTED_INT_KIND(2)
integer, parameter  :: SP = KIND(1.0)
integer, parameter  :: DP = KIND(1.0D0)
integer, parameter  :: DOUBLE = KIND(1.0D0)
integer, parameter  :: SPC = KIND((1.0,1.0))
integer, parameter  :: DPC = KIND((1.0D0,1.0D0))
integer, parameter  :: LGT = KIND(.true.)
integer, parameter  :: line_max_len = 8192 !< Max number of characters on line
real(sp), parameter :: PI=acos(-1.)
real(sp), parameter :: PIO2=acos(-1.)/2.
real(sp), parameter :: TWOPI=2.*acos(-1.)
real(sp), parameter :: FOURPI=4.*acos(-1.)
real(sp), parameter :: SQRT2=sqrt(2.)
real(sp), parameter :: EUL=0.5772156649015328606065120900824024310422_sp
real(sp), parameter :: um2a = 10000.
real(sp), parameter :: TINY=1e-10
real(sp), parameter :: SMALL=1e-6
real(sp), parameter :: MINEULSDEV=3.
real(sp), parameter :: MINTRSSDEV=0.5
real(sp), parameter :: FTOL=1e-4
real(dp), parameter :: DTINY=1e-10
real(dp), parameter :: DSMALL=1e-6
real(dp), parameter :: pisqr = PI*PI   ! PI^2.

! plan for the CTF
type :: ctfplan
    character(len=STDLEN) :: mode='' !< astig/noastig
    character(len=STDLEN) :: flag='' !< flag: <mul|flip|no> 
end type ctfplan

! constants for picker
real,    parameter :: PICKER_SHRINK        = 4.
real,    parameter :: PICKER_SHRINK_REFINE = 2.
integer, parameter :: PICKER_OFFSET        = 3

! constants for interpolation
real, parameter :: KBWINSZ = 1.5
real, parameter :: KBALPHA = 2.0

! CONSTANTS THAT CONTROL SEARCH AND CONVERGENCE
real,    parameter :: FRAC_SH_LIM     = 80.0 ! at what frac to turn on the shift search
real,    parameter :: SHW_FRAC_LIM    = 50.  ! at what frac to turn on shell-weights
real,    parameter :: EXTRINITHRESH   = 0.5
integer, parameter :: LPLIM1ITERBOUND = 5
integer, parameter :: LPLIM3ITERBOUND = 7

! endianness conversion
character(len=:), allocatable :: endconv

! number of threads global variable
integer(kind=c_int):: nthr_glob

! global distributed execution flag
logical :: l_distr_exec_glob

! global executable absolute path
character(len=STDLEN) :: exec_abspath_glob

end module simple_defs
