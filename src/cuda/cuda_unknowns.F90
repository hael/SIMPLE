module cuda_unknowns
  use, intrinsic :: ISO_C_BINDING
implicit none
  type, bind(C) :: cudaStream_t
     type (c_ptr) :: stream
   end type cudaStream_t

   type, bind(C) :: cudaEvent_t
      type (c_ptr) :: event
   end type cudaEvent_t

   type, bind(C) :: cudaArray_t
      type (c_ptr) :: array
   end type cudaArray_t

   type, bind(C) :: cudaContext_t
      type (c_ptr) :: context
   end type cudaContext_t

   type, bind(C) :: cudaModule_t
      type (c_ptr) :: mod
   end type cudaModule_t

   type, bind(C) :: cudaFunction_t
      type (c_ptr) :: func
   end type cudaFunction_t

   type, bind(C) :: cudaSurfRef_t
      type (c_ptr) :: surf_ref
   end type cudaSurfRef_t

   type, bind(C) :: cudaTexRef_t
      type (c_ptr) :: tex_ref
   end type cudaTexRef_t

   type, bind(C) :: cudaGraphicsResource_t
      type (c_ptr) :: graphics_resource
   end type cudaGraphicsResource_t

end module cuda_unknowns
