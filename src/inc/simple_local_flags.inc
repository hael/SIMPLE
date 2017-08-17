!------------------------------------------------------------------------------!
! SIMPLE v3.0 Elmlund & Elmlund Lab simplecryoem.com !
!------------------------------------------------------------------------------!
!> Simple debug and vebose flags
!! simple_local_flags.inc is for insertion into modules or subroutines. The
!! flags debug, warn, verbose can dictate how much is printed to the terminal
!! \author Michael Eager
! The code is distributed with the hope that it will be useful, but WITHOUT ANY
! WARRANTY. Redistribution and modification is regulated by the GNU General
! Public License.
!-----------------------------------------------------------------------------!
    logical :: debug=.false.
    logical :: warn=.false.
    logical :: verbose=.false.
#define _WPRINT_FILE_LINE_(FFFF,LLLL) if(debug.or.global_debug)write(*,'(A,A,I0,A)',advance='no') FFFF,':',LLLL,': '
#define __WPRINT_FILE_LINE_(FFFF,LLLL) write(*,'("DEBUG:",/,A,/,":",I0,":",A)',advance='no') FFFF,LLLL
#define _PRINT_FILE_LINE_(myFILE,myLine) if(debug.or.global_debug) print *, myFILE//':myLine',
#ifdef __PRETTY_FUNCTION__
# ifndef DebugPrint
#  define DebugPrint _PRINT_FILE_LINE_( __FILENAME__  , __LINE__ ); if(debug.or.global_debug) print *,__PRETTY_FUNCTION__,':',
# endif
#else
# ifndef DebugPrint
#  if defined(PGI)
#    define DebugPrint if(debug.or.global_debug) print *,
#  elif defined(Intel)
#    define LinePrint(N) #N
#    define DebugPrint if(debug.or.global_debug) print *,__FILENAME__//':',__LINE__,':',
#  else
#   define DebugPrint  _WPRINT_FILE_LINE_(__FILENAME__,__LINE__);if(debug.or.global_debug) print *,
#  endif
# endif
#endif
#ifndef WarnPrint
# define WarnPrint if(warn.or.debug.or.global_debug) print *, __FILENAME__,':',__LINE__,':Warning:',
#endif
#ifndef VerbosePrint
# define VerbosePrint if(verbose.or.global_verbose.or.debug.or.global_debug) print *, 
#endif

#ifndef ErrorMsg
#define ErrorMsg(msg) Error_Msg ( __FILENAME__, __LINE__, msg )
#endif