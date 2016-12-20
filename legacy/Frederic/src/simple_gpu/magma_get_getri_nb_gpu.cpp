/*
  -- MAGMA (version 1.0) --
  Returns the number of blocks
*/

#if defined (MAGMA) /*preprossing for the OPENCL environment */

// ==== Definition of blocking sizes for Tesla ===============================
#if (GPUSHMEM < 200)

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for zgetri based on n
*/
extern "C" int magma_get_getri_nb_gpu(int n)
{
  return 128;
  //return 64;
}
// ==== End Definition of blocking sizes for Tesla ===========================

#else

// ====     Definition of blocking sizes for Fermi ===========================

/* ////////////////////////////////////////////////////////////////////////////
   -- Return nb for zgetri based on n
*/
extern "C" int magma_get_getri_nb_gpu(int n)
{
  if (n<=3072)
    return 64;
  else
    return 128;
}
// ==== End Definition of blocking sizes for Fermi ===========================
#endif

#endif /* MAGMA */
