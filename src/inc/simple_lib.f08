    !! # define ABORT() call simple_abort (__FILENAME__, __LINE__)
#ifndef HALT
#define HALT(X) call simple_stop (X,__FILENAME__, __LINE__)
#endif

#ifndef SysError
# define SysError(usermsg,ioerr,iomsg)  call simple_stop(usermsg, __FILENAME__, __LINE__)
#endif
    
