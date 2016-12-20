/* mpe_logging_conf.h.  Generated from mpe_logging_conf.h.in by configure.  */
/* mpe_logging_conf.h.in.  Generated from configure.ac by autoheader.  */

/* Define if building universal (internal helper macro) */
/* #undef AC_APPLE_UNIVERSAL_BUILD */

/* The configuration header enables CLOG2 implementation related API */
#define CLOG_IMPL 1

/* Define if stdint.h should be before system headers */
/* #undef FIX_STDINT_ORDER */

/* Define if CRAY's FCD string is found */
/* #undef HAVE_CRAY_FCD_STRING */

/* Define to 1 if you have the <fcntl.h> header file. */
#define HAVE_FCNTL_H 1

/* Define to 1 if you have the `gethostname' function. */
#define HAVE_GETHOSTNAME 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define if libpthread.a(so) is defined */
#define HAVE_LIBPTHREAD 1

/* Define to 1 if you have the `lrand48' function. */
#define HAVE_LRAND48 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the `mkstemp' function. */
#define HAVE_MKSTEMP 1

/* Define to 1 if O_BINARY flag for open() exists. */
/* #undef HAVE_O_BINARY */

/* Define to 1 if O_LARGEFILE flag for open() exists. */
#define HAVE_O_LARGEFILE 1

/* Define if PMPI_Comm_Create_keyval() available */
#define HAVE_PMPI_COMM_CREATE_KEYVAL 1

/* Define if PMPI_Comm_free_keyval() available */
#define HAVE_PMPI_COMM_FREE_KEYVAL 1

/* Define if PMPI_Comm_get_attr() available */
#define HAVE_PMPI_COMM_GET_ATTR 1

/* Define if PMPI_Comm_set_attr() available */
#define HAVE_PMPI_COMM_SET_ATTR 1

/* Define to 1 if you have the <pthread.h> header file. */
#define HAVE_PTHREAD_H 1

/* Define to 1 if you have the `snprintf' function. */
#define HAVE_SNPRINTF 1

/* Define to 1 if you have the `srand48' function. */
#define HAVE_SRAND48 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdio.h> header file. */
#define HAVE_STDIO_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/bitypes.h> header file. */
#define HAVE_SYS_BITYPES_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define if gethostname needs a declaration */
/* #undef NEEDS_GETHOSTNAME_DECL */

/* Define if lrand48 needs a declaration */
/* #undef NEEDS_LRAND48_DECL */

/* Define if mkstemp needs a declaration */
/* #undef NEEDS_MKSTEMP_DECL */

/* Define if snprintf needs a declaration */
/* #undef NEEDS_SNPRINTF_DECL */

/* Define if srand48 needs a declaration */
/* #undef NEEDS_SRAND48_DECL */

/* Define to 1 if your C compiler doesn't accept -c and -o together. */
/* #undef NO_MINUS_C_MINUS_O */

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT ""

/* Define to the full name of this package. */
#define PACKAGE_NAME ""

/* Define to the full name and version of this package. */
#define PACKAGE_STRING ""

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME ""

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION ""

/* The size of `char', as computed by sizeof. */
/* #undef SIZEOF_CHAR */

/* The size of `int', as computed by sizeof. */
/* #undef SIZEOF_INT */

/* define if sizeof(int) = sizeof(void*) */
/* #undef SIZEOF_INT_IS_AINT */

/* The size of `long', as computed by sizeof. */
/* #undef SIZEOF_LONG */

/* The size of `long long', as computed by sizeof. */
/* #undef SIZEOF_LONG_LONG */

/* The size of `short', as computed by sizeof. */
/* #undef SIZEOF_SHORT */

/* The size of `void *', as computed by sizeof. */
#define SIZEOF_VOID_P 8

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define WORDS_BIGENDIAN to 1 if your processor stores words with the most
   significant byte first (like Motorola and SPARC, unlike Intel). */
#if defined AC_APPLE_UNIVERSAL_BUILD
# if defined __BIG_ENDIAN__
#  define WORDS_BIGENDIAN 1
# endif
#else
# ifndef WORDS_BIGENDIAN
/* #  undef WORDS_BIGENDIAN */
# endif
#endif

/* Enable large inode numbers on Mac OS X 10.5.  */
#ifndef _DARWIN_USE_64_BIT_INODE
# define _DARWIN_USE_64_BIT_INODE 1
#endif

/* Number of bits in a file offset, on hosts where this is settable. */
/* #undef _FILE_OFFSET_BITS */

/* Define if 64-bit file address support in 32-bit OS. */
#define _LARGEFILE64_SOURCE 1

/* Define for large files, on AIX-style hosts. */
/* #undef _LARGE_FILES */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */


/* Temporary fix of undefined int64_t on SLES8 with "icc -no-gcc" */
#if defined(FIX_STDINT_ORDER)
#include "clog_inttypes.h"
#endif

/* Define WINDOWS specific features */
/*
   Windows' open() opens an ASCII file by default, add Windows specific
   flag O_BINARY to open()'s argument
*/
#if !defined( OPEN )

#if defined(HAVE_WINDOWS_H)

#define OPEN( a, b, c )    open( a, b | O_BINARY, c )

#else

#if defined(HAVE_O_LARGEFILE)

#if defined(HAVE_O_BINARY)
#define OPEN( a, b, c )    open( a, b | O_LARGEFILE | O_BINARY, c )
#else
#define OPEN( a, b, c )    open( a, b | O_LARGEFILE, c )
#endif    /* HAVE_O_BINARY */

#else

#if defined(HAVE_O_BINARY)
#define OPEN( a, b, c )    open( a, b | O_BINARY, c )
#else
#define OPEN( a, b, c )    open( a, b, c )
#endif    /* HAVE_O_BINARY */

#endif    /* HAVE_O_LARGEFILE */

#endif    /* HAVE_WINDOWS_H */

#endif    /* OPEN */

