
#define CUDASIZE_OF( element ) SIZE( TRANSFER( element, ['a'] ) )

#define HANDLE_CUDAERROR( return_err, bool , msg)  \
call my_cudaErrorCheck( return_err, bool);\
if(bool) call simple_exception(msg, __FILENAME__ , __LINE__)
