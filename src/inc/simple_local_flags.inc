!------------------------------------------------------------------------------!
! SIMPLE v3.0 Elmlund & Elmlund Lab simplecryoem.com !
!------------------------------------------------------------------------------!
!> Simple debug and vebose flags
!!
!! simple_local_flags.inc is for insertion into modules or subroutines. The
!! flags debug, warn, verbose can dictate how much is printed to the terminal
!
!! @author Michael Eager
!
! The code is distributed with the hope that it will be useful, but WITHOUT ANY
! WARRANTY. Redistribution and modification is regulated by the GNU General
! Public License.
!-----------------------------------------------------------------------------!
#if defined(_DEBUG)
    logical :: debug=.true.
    logical :: warn=.true.
    logical :: verbose=.true.
#else
    logical :: debug=.false.
    logical :: warn=.false.
    logical :: verbose=.false.
#endif

#define _WPRINT_FILE_LINE_(FFFF,LLLL) if(debug.or.global_debug)write(*,'(A,A,I0,A)',advance='no') FFFF,':',LLLL,': ';
#define _PRINT_FILE_LINE_(myFILE,myLine) if(debug.or.global_debug) print *, myFILE//':myLine',

#ifdef __PRETTY_FUNCTION__
# ifndef DebugPrint
#  define DebugPrint _PRINT_FILE_LINE_( __FILENAME__  , __LINE__ ); if(debug.or.global_debug) print *,__PRETTY_FUNCTION__,':',
# endif
#else
# ifndef DebugPrint
#  if defined(_DEBUG) || ! defined(PGI)
#   define DebugPrint _WPRINT_FILE_LINE_(__FILENAME__,__LINE__)if(debug.or.global_debug) print *,
#  else
#   define DebugPrint ! /* Disabled debug printing*/
#  endif
# endif
#endif
#ifndef WarnPrint
# define WarnPrint if(warn.or.debug.or.global_debug) print *, __FILENAME__,':',__LINE__,':Warning:',
#endif
#ifndef VerbosePrint
# define VerbosePrint if(verbose.or.global_verbose.or.debug.or.global_debug) print *, 
#endif
#ifndef ErrorPrint
# define ErrorPrint print *, __FILENAME__,':',__LINE__,':Error:',
#endif
