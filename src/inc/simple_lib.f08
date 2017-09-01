    !! # define ABORT() call simple_abort (__FILENAME__, __LINE__)
#ifndef HALT
# define HALT(X) call simple_stop (X,__FILENAME__, __LINE__)
#endif

#ifndef ABANDON
# define ABANDON(usermsg) print *,usermsg,__FILENAME__, __LINE__;call abort
#endif

#ifndef allocchk
#define allocchk( X ) call alloc_errchk(X ,alloc_stat,__FILENAME__,__LINE__)
#endif

