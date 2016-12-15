 
/* For now, the GPU only supports a 32-bit address space, so device pointers
   can be represented as INTEGER*4 in Fortran. In the future, device pointers
   may become 64-bit pointers, and will have to be represented as INTEGER*8 in
   Fortran, at which point devptr_t needs to be typedef'ed as long long.
*/
typedef size_t devptr_t;

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */
  int SIMPLE_CUDNN_CREATE(cudnnHandle_t *handle);
  int SIMPLE_CUDNN_DESTROY(cudnnHandle_t *handle);
#if defined(__cplusplus)
}
#endif /* __cplusplus */

/*
 * Fortran callable thin wrappers. Fortran application must allocate and
 * deallocate GPU memory, and copy data up and down.
 */
#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

  //Insert future wrappers here.

#if defined(__cplusplus)
}
#endif /* __cplusplus */
