!! -*- mode: f90 -*-
!!
!! simple_local_flags.inc is for insertion
!! into modules or subroutines.  The flags:
!! debug, warn, verbose
#if defined(_DEBUG)
    logical :: debug=.true.
    logical :: warn=.true.
    logical :: verbose=.true.
#else
    logical :: debug=.false.
    logical :: warn=.false.
    logical :: verbose=.false.
#endif  ! /* _DEBUG */

#ifdef __PRETTY_FUNCTION__
#  ifndef DebugPrint
#    define DebugPrint(...) if( debug.or.global_debug ) print *, __FILE__,':',__LINE__,':',__PRETTY_FUNCTION__,':', __VA_ARGS__
#  endif
#else
#  ifndef DebugPrint
#    define DebugPrint(...) if( debug.or.global_debug ) print *, __FILE__,':',__LINE__,':', __VA_ARGS__
#  endif
#endif
#ifndef WarnPrint
#  define WarnPrint(...) if( warn.or.debug.or.global_warn.or.global_debug ) print *, __FILE__,':',__LINE__,':Warning:', __VA_ARGS__
#endif
#ifndef VerbosePrint
#  define VerbosePrint(...) if( verbose.or.warn.or.debug.or.global_verbose.or.global_warn.or.global_debug ) print *,  __VA_ARGS__
#endif
