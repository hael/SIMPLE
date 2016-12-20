// max GPU's to manage for multi-GPU parallel compute
#if defined (OPENCL) /*preprossing for the OPENCL environment */
const unsigned int MAX_GPU_COUNT = 8;
#endif /* OPENCL */
