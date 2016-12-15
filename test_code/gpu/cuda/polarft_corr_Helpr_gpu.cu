/*******************************************************************************
 * SIMPLE routine for CUDA BLAS, including a CUDA kernel and wappers           *
 *******************************************************************************
 *
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 5th Feb 2016
 *      Monash University
 *      February 2016
 *
 *      Helper functions to assist the other kernels and debugging.
 *
 * @precisions normal z -> s d c
 */
#include "common_magma.h"
#include "commonblas_zz2d.h"
#include "polarft_gpu.h"
#include "simple.h"

#if defined (OPENMP) /*preprossing for the OpenMP environment */
#include <omp.h>
#endif

#define debug false
#define debug_high false
#define debug_write false
#define debug_write_AandB false
#define debug_write_C true
#define nthr 12

#if defined (CUDA) /*preprossing for the CUDA environment */

/* getting the error id from the cuda */
extern "C" int
get_error_id_corr_Hadmr_gpu(cudaError_t error_id) {
  int rc = RC_SUCCESS;
  printf("cudaDriverGetVersion returned %d\n-> %s\n", 
         (int)error_id, cudaGetErrorString(error_id));
  printf("Result = FAIL\n");
  printf("Line %i:%s ---> %s\n",__LINE__,__FILE__,__FUNCTION__);
  rc = (int)error_id;
  exit(EXIT_FAILURE);
  return rc;
}
/* getting the warning message from the cuda */
extern "C"
int get_warning_message_corr_Hadmr_gpu() {
  int rc = RC_SUCCESS;

  printf("***************************WARNING*****************************\n");
  printf("You need to compile with -DCUDA to acces the CUDA environment  \n");
  printf("computation using GPU                                          \n");
  printf("***************************************************************\n");
  printf("\n");
  printf("Line %i:%s ---> %s\n",__LINE__,__FILE__,__FUNCTION__);
  printf("\n");
  rc = RC_FAIL;

  return rc;
}
/* warning message for the missmatch in the number of threadsPerBlock */
extern "C"
int get_warning_threadsPerBlock(polar_corr_calc_t *s_polar,
                                int threadsPerBlock) {
  int rc = RC_SUCCESS;

  printf(ANSI_COLOR_BRIGHT_CYAN);
  printf("***************************WARNING*****************************\n");
  printf("Number of threads from the s_polar structure is not the same   \n");
  printf("as the number of threadsPerBlock. The size of the share memory \n");
  printf("array shared_A in kernel: "ANSI_COLOR_BRIGHT_RED"%i "
         ANSI_COLOR_BRIGHT_CYAN"is not the same as entered     \n",
         s_polar->ikrnl);
  printf("threadsPerBlock check that computation is correct from the GPU \n");
  printf(ANSI_COLOR_BRIGHT_GREEN"        threadsPerBlock: %i, "
         ANSI_COLOR_BRIGHT_BLUE"s_polar->threadsPerBlock: %i"
         ANSI_COLOR_BRIGHT_CYAN"\n",threadsPerBlock,s_polar->threadsPerBlock);
  printf("Line %i:%s ---> %s\n",__LINE__,__FILE__,__FUNCTION__);
  printf("***************************************************************\n");
  printf(ANSI_COLOR_RESET);
  rc = RC_FAIL;

  return rc;
}
/* getting the warning message from the cuda */
extern "C"
int print_summing_method(int isum) {
  int rc = RC_SUCCESS;

  printf(ANSI_COLOR_BRIGHT_CYAN"Summing method  : "
         ANSI_COLOR_RESET);
  switch (isum) {
  case 0: printf(ANSI_COLOR_BRIGHT_YELLOW"Summing with standard vectorisation\n"
                 ANSI_COLOR_RESET); break;
  case 1: printf(ANSI_COLOR_BRIGHT_GREEN"Summing with pointer arithmetic\n"
                 ANSI_COLOR_RESET); break;
  default:
    printf(ANSI_COLOR_BRIGHT_CYAN"Unknown summation technique to be decided\n"
           ANSI_COLOR_RESET);
    rc = RC_FAIL;
  }
   
  return rc;
}
/* printing out the true dimension of the problem */ 
extern "C" int
print_npart_nrot_nk(int npart, int nrot, int nk) {
  int rc = RC_SUCCESS;

  int npts = npart * nrot * nk;

  printf(ANSI_COLOR_BRIGHT_YELLOW"npart=%i, "
         ANSI_COLOR_BRIGHT_GREEN"nrot=%i, "
         ANSI_COLOR_BRIGHT_BLUE  "nk=%i, "
         ANSI_COLOR_BRIGHT_CYAN  "npts=%i, npart*nrot*nk=%i\n",
         npart, nrot, nk, npts, npart*nrot*nk);

  printf(ANSI_COLOR_BRIGHT_RED
         "npart*nrot*nk*sizeof(cuFloatComplex)=%lu, Memory=%.2f MB\n"
         ANSI_COLOR_RESET,
         npart*nrot*nk*sizeof(cuFloatComplex),
         npart*nrot*nk*sizeof(cuFloatComplex)/1.e6);

  return rc;
}
/* printing out the true dimension of the problem */ 
extern "C" int
print_vx_vy_vz(int vx, int vy, int vz) {
  int rc = RC_SUCCESS;

  int npts = vx * vy * vz;

  printf(ANSI_COLOR_BRIGHT_YELLOW"vx=%i, "
         ANSI_COLOR_BRIGHT_GREEN"vy=%i, "
         ANSI_COLOR_BRIGHT_BLUE  "vz=%i, "
         ANSI_COLOR_BRIGHT_CYAN  "npts=%i, vx*vy*vz=%i\n",
         vx, vy, vz, npts, vx*vy*vz);

  printf(ANSI_COLOR_BRIGHT_RED
         "vx*vy*vz*sizeof(cuFloatComplex)=%lu, Memory=%.2f MB\n"
         ANSI_COLOR_RESET,
         vx*vy*vz*sizeof(cuFloatComplex),
         vx*vy*vz*sizeof(cuFloatComplex)/1.e6);

  return rc;
}
/* getting info for the 1D mesh */
extern "C"
int print_1D_mesh(int s_cnstr_mesh_1D, mesh_1D *p_mesh_1D, int N, int threadsPerBlock, int blocksPerGrid) {
  int rc = RC_SUCCESS;

  printf(ANSI_COLOR_BRIGHT_CYAN"Dim 1D          : "
         ANSI_COLOR_BRIGHT_YELLOW"N=%i, "
         ANSI_COLOR_BRIGHT_GREEN "threadsPerBlock=%i, "
         ANSI_COLOR_BRIGHT_BLUE  "blocksPerGrid=%i\n"
         ANSI_COLOR_RESET,
         N, threadsPerBlock, blocksPerGrid);
  if ( s_cnstr_mesh_1D == 1) {
  printf(ANSI_COLOR_BRIGHT_CYAN"mesh1D Objct    : "
         ANSI_COLOR_BRIGHT_YELLOW"p_mesh_1D->get_mesh1D_chunk=%i, "
         ANSI_COLOR_BRIGHT_GREEN "p_mesh_1D->get_mesh1D_threadsPerBlock=%i, "
         ANSI_COLOR_BRIGHT_BLUE  "p_mesh_1D->get_mesh1D_blocksPerGrid=%i\n"
         ANSI_COLOR_RESET,
         p_mesh_1D->get_mesh1D_chunk(), p_mesh_1D->get_mesh1D_threadsPerBlock(),
         p_mesh_1D->get_mesh1D_blocksPerGrid());
    
  }  
  return rc;
}
/* getting the warning message from the cuda */
extern "C"
int print_3D_mesh(int s_cnstr_mesh_3D, mesh_3D *p_mesh_3D,
                  polar_corr_calc_t *s_polar, pfts_Sizes *p_pfts_Sizes,
                  int nx, int ny, int nz) {
  int rc = RC_SUCCESS;

  int gridx,gridy,gridz;
  
  printf(ANSI_COLOR_BRIGHT_CYAN"mesh3D Objct: "
         ANSI_COLOR_BRIGHT_YELLOW"p_mesh_3D->get_mesh3D_nx()=%i, "
         ANSI_COLOR_BRIGHT_GREEN "p_mesh_3D->get_mesh3D_ny()=%i, "
         ANSI_COLOR_BRIGHT_BLUE  "p_mesh_3D->get_mesh3D_nz()=%i, "
                                 "total(nx*ny*nz)=%i\n",
         p_mesh_3D->get_mesh3D_nx(), p_mesh_3D->get_mesh3D_ny(),
         p_mesh_3D->get_mesh3D_nz(), p_mesh_3D->get_mesh3D_nx()*
         p_mesh_3D->get_mesh3D_ny()* p_mesh_3D->get_mesh3D_nz() );
  if ( s_cnstr_mesh_3D == 1) {
  printf(ANSI_COLOR_BRIGHT_CYAN"              "
         ANSI_COLOR_BRIGHT_YELLOW"p_mesh_3D->get_mesh3D_gridx()=%i, "
         ANSI_COLOR_BRIGHT_GREEN "p_mesh_3D->get_mesh3D_gridy()=%i, "
         ANSI_COLOR_BRIGHT_BLUE  "p_mesh_3D->get_mesh3D_gridz()=%i, "
                                 "total(*)=%i\n",
         p_mesh_3D->get_mesh3D_gridx(), p_mesh_3D->get_mesh3D_gridy(),
         p_mesh_3D->get_mesh3D_gridz(),
         p_mesh_3D->get_mesh3D_gridx()*p_mesh_3D->get_mesh3D_gridy()*
         p_mesh_3D->get_mesh3D_gridz() );
  }
  printf(ANSI_COLOR_BRIGHT_CYAN"Dim_3D      : "
         ANSI_COLOR_BRIGHT_YELLOW"(float)nx=%f, "
         ANSI_COLOR_BRIGHT_GREEN "(float)ny=%f, "
         ANSI_COLOR_BRIGHT_BLUE  "(float)nz=%f\n",
         (float)nx, (float)ny, (float)nz);
  printf(ANSI_COLOR_BRIGHT_CYAN"Threads     : "
         ANSI_COLOR_BRIGHT_YELLOW"nx=%i, "
         ANSI_COLOR_BRIGHT_GREEN "ny=%i, "
         ANSI_COLOR_BRIGHT_BLUE  "nz=%i, total(nx*ny*nz)=%i\n",
         nx,ny,nz,nx*ny*nz);

  gridx =    p_pfts_Sizes->get_npart()/(float)s_polar->nx +
         (   p_pfts_Sizes->get_npart()%s_polar->nx!=0     );
  gridy =  p_pfts_Sizes->get_nhlfrot()/(float)s_polar->ny +
         ( p_pfts_Sizes->get_nhlfrot()%s_polar->ny!=0     );
  gridz =       p_pfts_Sizes->get_nk()/(float)s_polar->nz +
         (      p_pfts_Sizes->get_nk()%s_polar->nz!=0     );

  printf(ANSI_COLOR_BRIGHT_CYAN"Grid3D      : "
         ANSI_COLOR_BRIGHT_YELLOW"npart/(float)nx+(npart mod(nx)!=0)= %i, "
         ANSI_COLOR_BRIGHT_GREEN"nhlfrot/(float)ny+(nhlfrot mod(ny)!=0)=%i, "
         ANSI_COLOR_BRIGHT_BLUE"nk/(float)nz+(nk mod(nz)!=0)=%i, total(*)=%i\n"
         ANSI_COLOR_RESET,
         gridx,gridy,gridz, gridx*gridy*gridz);
  
  return rc;
}
/* getting the warning message from the cuda */
extern "C"
int print_3D_V_mesh(int s_cnstr_mesh_3D, mesh_3D *p_mesh_3D,
                    polar_corr_calc_t *s_carte,
                    img_2D_cart_Sizes *p_img_2D_cart_Sizes,
                    int nx, int ny, int nz) {
  int rc = RC_SUCCESS;

  //pfts_Sizes *p_pfts_Sizes
  
  int gridx,gridy,gridz;
  
  printf(ANSI_COLOR_BRIGHT_CYAN"mesh3D Objct: "
         ANSI_COLOR_BRIGHT_YELLOW"p_mesh_3D->get_mesh3D_nx()=%i, "
         ANSI_COLOR_BRIGHT_GREEN "p_mesh_3D->get_mesh3D_ny()=%i, "
         ANSI_COLOR_BRIGHT_BLUE  "p_mesh_3D->get_mesh3D_nz()=%i, "
                                 "total(nx*ny*nz)=%i\n",
         p_mesh_3D->get_mesh3D_nx(), p_mesh_3D->get_mesh3D_ny(),
         p_mesh_3D->get_mesh3D_nz(), p_mesh_3D->get_mesh3D_nx()*
         p_mesh_3D->get_mesh3D_ny()* p_mesh_3D->get_mesh3D_nz() );
  if ( s_cnstr_mesh_3D == 1) {
  printf(ANSI_COLOR_BRIGHT_CYAN"              "
         ANSI_COLOR_BRIGHT_YELLOW"p_mesh_3D->get_mesh3D_gridx()=%i, "
         ANSI_COLOR_BRIGHT_GREEN "p_mesh_3D->get_mesh3D_gridy()=%i, "
         ANSI_COLOR_BRIGHT_BLUE  "p_mesh_3D->get_mesh3D_gridz()=%i, "
                                 "total(*)=%i\n",
         p_mesh_3D->get_mesh3D_gridx(), p_mesh_3D->get_mesh3D_gridy(),
         p_mesh_3D->get_mesh3D_gridz(),
         p_mesh_3D->get_mesh3D_gridx()*p_mesh_3D->get_mesh3D_gridy()*
         p_mesh_3D->get_mesh3D_gridz() );
  }
  printf(ANSI_COLOR_BRIGHT_CYAN"Dim_3D      : "
         ANSI_COLOR_BRIGHT_YELLOW"(float)nx=%f, "
         ANSI_COLOR_BRIGHT_GREEN "(float)ny=%f, "
         ANSI_COLOR_BRIGHT_BLUE  "(float)nz=%f\n",
         (float)nx, (float)ny, (float)nz);
  printf(ANSI_COLOR_BRIGHT_CYAN"Threads     : "
         ANSI_COLOR_BRIGHT_YELLOW"nx=%i, "
         ANSI_COLOR_BRIGHT_GREEN "ny=%i, "
         ANSI_COLOR_BRIGHT_BLUE  "nz=%i, total(nx*ny*nz)=%i\n",
         nx,ny,nz,nx*ny*nz);

  gridx =    p_img_2D_cart_Sizes->get_2D_vx()/(float)s_carte->nx +
         (   p_img_2D_cart_Sizes->get_2D_vx()%s_carte->nx!=0     );
  gridy =  p_img_2D_cart_Sizes->get_2D_vy()/(float)s_carte->ny +
         ( p_img_2D_cart_Sizes->get_2D_vy()%s_carte->ny!=0     );
  gridz =       p_img_2D_cart_Sizes->get_2D_vz()/(float)s_carte->nz +
         (      p_img_2D_cart_Sizes->get_2D_vz()%s_carte->nz!=0     );

  printf(ANSI_COLOR_BRIGHT_CYAN"Grid3D      : "
         ANSI_COLOR_BRIGHT_YELLOW"vx/(float)nx+(vx mod(nx)!=0)= %i, "
         ANSI_COLOR_BRIGHT_GREEN"vy/(float)ny+(vy mod(ny)!=0)=%i, "
         ANSI_COLOR_BRIGHT_BLUE"vz/(float)nz+(vz mod(nz)!=0)=%i, total(*)=%i\n"
         ANSI_COLOR_RESET,
         gridx,gridy,gridz, gridx*gridy*gridz);
  
  return rc;
}


/* printing info for s_polar struct */
extern "C"
int print_s_polar_struct(polar_corr_calc_t *s_polar) {
  int rc = RC_SUCCESS;

  printf(ANSI_COLOR_BRIGHT_CYAN"s_polar     : "
         ANSI_COLOR_BRIGHT_YELLOW"s_polar->r_polar=%f, "
         ANSI_COLOR_BRIGHT_YELLOW"s_polar->sumasq_polar=%f, "
         ANSI_COLOR_BRIGHT_YELLOW"s_polar->sumbsq_polar=%f, \n"
         ANSI_COLOR_BRIGHT_GREEN "              s_polar->ikrnl=%i, "
         ANSI_COLOR_BRIGHT_BLUE  "s_polar->threadsPerBlock=%i, \n"
         ANSI_COLOR_BRIGHT_YELLOW"              s_polar->nx=%i, "
         ANSI_COLOR_BRIGHT_GREEN "s_polar->ny=%i, "
         ANSI_COLOR_BRIGHT_BLUE  "s_polar->nz=%i, "
         ANSI_COLOR_BRIGHT_BLUE  "at Line %i %s\n"
         ANSI_COLOR_RESET,
         s_polar->r_polar, s_polar->sumasq_polar, s_polar->sumbsq_polar,
         s_polar->ikrnl, s_polar->threadsPerBlock,
         s_polar->nx,s_polar->ny,s_polar->nz,
         __LINE__,__FUNCTION__);
  
  return rc;
}
/* printing info for s_debug struct */
extern "C"
int print_s_debug_struct(debug_gpu_t *s_debug_gpu) {
  int rc = RC_SUCCESS;

  printf(ANSI_COLOR_BRIGHT_CYAN"s_debug_gpu : "

         ANSI_COLOR_BRIGHT_YELLOW "get_bool(debug_i): %s, "
         ANSI_COLOR_BRIGHT_GREEN  "get_bool(debug_cpu_i): %s, "
         ANSI_COLOR_BRIGHT_BLUE   "get_bool(debug_high_i): %s, "
         ANSI_COLOR_BRIGHT_MAGENTA"get_bool(debug_write_i): %s, "
         ANSI_COLOR_BRIGHT_RED    "get_bool(debug_write_C_i): %s, \n"

         ANSI_COLOR_BRIGHT_YELLOW "                       debug_i : %i, "
         ANSI_COLOR_BRIGHT_GREEN  "            debug_cpu_i : %i, "
         ANSI_COLOR_BRIGHT_BLUE   "             debug_high_i : %i, "
         ANSI_COLOR_BRIGHT_MAGENTA"            debug_write_i : %i, "
         ANSI_COLOR_BRIGHT_RED    "             debug_write_C_i : %i, \n"

         ANSI_COLOR_BRIGHT_YELLOW "                       debug_i : %s, "
         ANSI_COLOR_BRIGHT_GREEN  "         debug_cpu_i : %s, "
         ANSI_COLOR_BRIGHT_BLUE   "         debug_high_i : %s, "
         ANSI_COLOR_BRIGHT_MAGENTA"         debug_write_i : %s, "
         ANSI_COLOR_BRIGHT_RED    "         debug_write_C_i : %s, \n"

         ANSI_COLOR_RESET,

         get_bool(s_debug_gpu->debug_i),
         get_bool(s_debug_gpu->debug_cpu_i),
         get_bool(s_debug_gpu->debug_high_i),
         get_bool(s_debug_gpu->debug_write_i),
         get_bool(s_debug_gpu->debug_write_C_i),

         s_debug_gpu->debug_i,
         s_debug_gpu->debug_cpu_i,
         s_debug_gpu->debug_high_i,
         s_debug_gpu->debug_write_i,
         s_debug_gpu->debug_write_C_i,
         
         s_debug_gpu->debug_i?"true":"false",
         s_debug_gpu->debug_cpu_i?"true":"false",
         s_debug_gpu->debug_high_i?"true":"false",
         s_debug_gpu->debug_write_i?"true":"false",
         s_debug_gpu->debug_write_C_i?"true":"false"

         );

  return rc;
}
/* printing info for s_debug struct */
extern "C"
int print_s_bench_struct(bench_t *s_bench) {
  int rc = RC_SUCCESS;

  printf(ANSI_COLOR_BRIGHT_CYAN"s_bench     : "
         ANSI_COLOR_BRIGHT_YELLOW "bench_i: %i, "
         ANSI_COLOR_BRIGHT_GREEN  "bench_write_i: %i\n"
         ANSI_COLOR_RESET,
         s_bench->bench_i,
         s_bench->bench_write_i);

  return rc;
}
/* printing info for s_dev[] struct */
//extern "C"
int print_s_devD_struct(deviceDetails_t * s_devD) {
  int rc = RC_SUCCESS;
  /*indexers*/
  int i;
  cudaDeviceProp deviceProp;
  cudaError_t err;

  printf(ANSI_COLOR_BRIGHT_CYAN"s_devD[]    : \n");
  //         ANSI_COLOR_BRIGHT_BLUE"Line %i %s\n"
  //         ANSI_COLOR_BRIGHT_CYAN,__LINE__,__FUNCTION__);
  printf("+-----------------------------------------------------------+");
  printf("\n| NVIDIA-GPU SYSTEM CUDA: "
         ANSI_COLOR_BRIGHT_YELLOW"%4.2f"
         ANSI_COLOR_BRIGHT_CYAN"      Drivers Version: "
         ANSI_COLOR_BRIGHT_YELLOW"%4.2f"
         ANSI_COLOR_BRIGHT_CYAN"   |",s_devD[0].d_runver,s_devD[0].d_ver);
  printf("\n| Number of devices on Sys: "
         ANSI_COLOR_BRIGHT_YELLOW"%i"
         ANSI_COLOR_BRIGHT_CYAN"                               |",s_devD[0].ndev);
  printf("\n+--------------------------+--------------------------------+");
  printf("\n| Device Names             |  ");
  for(i=0;i<s_devD[0].ndev;i++){
    err = cudaGetDeviceProperties(&deviceProp, i);
    if (err != cudaSuccess) { rc = get_error_id_corr_Hadmr_gpu(err); }
    s_devD[i].dev_name = deviceProp.name;
    printf(ANSI_COLOR_BRIGHT_GREEN"%s  ",s_devD[i].dev_name);}
  printf(ANSI_COLOR_BRIGHT_CYAN"         |");
  printf("\n| Total global Mem MB      |  ");
  for(i=0;i<s_devD[0].ndev;i++){printf(ANSI_COLOR_BRIGHT_GREEN
                                       "%6.2f   "
                                       ,s_devD[i].tot_global_mem_MB);}
  printf(ANSI_COLOR_BRIGHT_CYAN"         |");
  printf("\n| N MultiProc              |  ");
  for(i=0;i<s_devD[0].ndev;i++){printf(ANSI_COLOR_BRIGHT_GREEN
                                       "%5d      ",s_devD[i].nMultiProc);}
  printf(ANSI_COLOR_BRIGHT_CYAN"        |");
  printf("\n| N cudacores/MultiProc    |  ");
  for(i=0;i<s_devD[0].ndev;i++){printf(ANSI_COLOR_BRIGHT_GREEN
                                       "%5d      ",s_devD[i].ncudacores_per_MultiProc);}
  printf(ANSI_COLOR_BRIGHT_CYAN"        |");
  printf("\n| N cudacores              |  ");
  for(i=0;i<s_devD[0].ndev;i++){printf(ANSI_COLOR_BRIGHT_GREEN
                                       "%5d      ",s_devD[i].ncudacores);}
  printf(ANSI_COLOR_BRIGHT_CYAN"        |");
  printf("\n| Is SMsuitable 0:no 1:yes |  ");
  for(i=0;i<s_devD[0].ndev;i++){printf(ANSI_COLOR_BRIGHT_GREEN
                                       "%5d      ",s_devD[i].is_SMsuitable);}
  printf(ANSI_COLOR_BRIGHT_CYAN"        |");
  printf("\n| N registers/blk          |  ");
  for(i=0;i<s_devD[0].ndev;i++){printf(ANSI_COLOR_BRIGHT_GREEN
                                       "%5d      ",s_devD[i].nregisters_per_blk);}
  printf(ANSI_COLOR_BRIGHT_CYAN"        |");
  printf("\n| WarpSze                  |  ");
  for(i=0;i<s_devD[0].ndev;i++){printf(ANSI_COLOR_BRIGHT_GREEN
                                       "%5d      ",s_devD[i].warpSze);}
  printf(ANSI_COLOR_BRIGHT_CYAN"        |");
  printf("\n| Maxthreads/MultiProc     |  ");
  for(i=0;i<s_devD[0].ndev;i++){printf(ANSI_COLOR_BRIGHT_GREEN
                                       "%5d      ",s_devD[i].maxthreads_per_mp);}
  printf(ANSI_COLOR_BRIGHT_CYAN"        |");
  printf("\n| Maxthreads/blk           |  ");
  for(i=0;i<s_devD[0].ndev;i++){printf(ANSI_COLOR_BRIGHT_GREEN
                                       "%5d      ",s_devD[i].maxthreads_per_blk);}
  printf(ANSI_COLOR_BRIGHT_CYAN"        |");
  printf("\n| Is ecc: err corr code Mem|  ");
  for(i=0;i<s_devD[0].ndev;i++){printf(ANSI_COLOR_BRIGHT_GREEN
                                       "%5d      ",s_devD[i].is_ecc);}
  printf(ANSI_COLOR_BRIGHT_CYAN"        |");
  printf("\n| Is p2p                   |");
  printf("\n+--------------------------+--------------------------------+\n");

  printf(ANSI_COLOR_RESET);
  
  return rc;
}
/* Checker for the mesh3D via mesh3D getting the kernel call info */
extern "C"
int check_grid3D(int gridx, int gridy, int gridz, polar_corr_calc_t *s_polar,
                 pfts_Sizes *p_pfts_Sizes) {
  int rc = RC_SUCCESS;

  if ( gridx != ( p_pfts_Sizes->get_npart()/(float)s_polar->nx   +
                  ( p_pfts_Sizes->get_npart()%s_polar->nx!=0 )   ) )
    {rc = RC_FAIL;}
  if ( gridy != ( p_pfts_Sizes->get_nhlfrot()/(float)s_polar->ny +
                  ( p_pfts_Sizes->get_nhlfrot()%s_polar->ny!=0 ) ) )
    {rc = RC_FAIL;}
  if ( gridz != ( p_pfts_Sizes->get_nk()/(float)s_polar->nz      +
                  ( p_pfts_Sizes->get_nk()%s_polar->nz!=0      ) ) )
    {rc = RC_FAIL;}
  
  return rc;
}
/* Checker for the mesh3D via mesh3D getting the kernel call info */
extern "C"
int check_grid3D_V(int gridx, int gridy, int gridz, polar_corr_calc_t *s_carte,
                 img_2D_cart_Sizes *p_img_2D_cart_Sizes) {
  int rc = RC_SUCCESS;

  if ( gridx != ( p_img_2D_cart_Sizes->get_2D_vx()/(float)s_carte->nx   +
                ( p_img_2D_cart_Sizes->get_2D_vx()%s_carte->nx!=0 )   ) )
    {rc = RC_FAIL;}
  if ( gridy != ( p_img_2D_cart_Sizes->get_2D_vy()/(float)s_carte->ny +
                ( p_img_2D_cart_Sizes->get_2D_vy()%s_carte->ny!=0 ) ) )
    {rc = RC_FAIL;}
  if ( gridz != ( p_img_2D_cart_Sizes->get_2D_vz()/(float)s_carte->nz      +
                ( p_img_2D_cart_Sizes->get_2D_vz()%s_carte->nz!=0      ) ) )
    {rc = RC_FAIL;}
  
  return rc;
}
/* getting the kernel call info */
extern "C"
int print_iKernel_threadsPerBlock(polar_corr_calc_t *s_polar,
                                 int N, int threadsPerBlock,
                                 int blocksPerGrid) {
  int rc = RC_SUCCESS;

  printf(ANSI_COLOR_BRIGHT_CYAN"Kernel details  : "
         ANSI_COLOR_BRIGHT_YELLOW"N=%i, "
         ANSI_COLOR_BRIGHT_GREEN "iKernel=%i, "
         ANSI_COLOR_BRIGHT_BLUE  "s_polar->threadsPerBlock=%i\n"
         ANSI_COLOR_RESET,
         N, s_polar->ikrnl, s_polar->threadsPerBlock);
  if (s_polar->threadsPerBlock != threadsPerBlock){
    rc = get_warning_threadsPerBlock(s_polar, threadsPerBlock);
  }
  return rc;
}
/* //////////////////////////////////////////////////////////////////////////// 
 -- Polarft_gpu kernel N_N helpers
*/
/* printing out the input values */ 
extern "C" int
print_function_header_N_N_info(float *r,
                               cuFloatComplex *C,
                               const cuFloatComplex *A,
                               const cuFloatComplex *B,
                               int npart, int nrot, int nk,
                               float alpha)
{
  int rc = RC_SUCCESS;
  int i,j,k;
  float h_r = *r;

  int npts = npart * nrot * nk;

  printf(ANSI_COLOR_BRIGHT_YELLOW"h_r=%f, ", h_r);

  printf(ANSI_COLOR_BRIGHT_YELLOW"npart=%i, "
         ANSI_COLOR_BRIGHT_GREEN"nrot=%i, "
         ANSI_COLOR_BRIGHT_BLUE  "nk=%i, "
         ANSI_COLOR_BRIGHT_CYAN  "npts=%i, npart*nrot*nk=%i\n",
         npart, nrot, nk, npts, npart*nrot*nk);

  printf(ANSI_COLOR_BRIGHT_RED
         "npart*nrot*nk*sizeof(cuFloatComplex)=%lu, Memory=%.2f MB\n"
         ANSI_COLOR_RESET,
         npart*nrot*nk*sizeof(cuFloatComplex),
         npart*nrot*nk*sizeof(cuFloatComplex)/1.e6);

  FILE * pftFile;
  pftFile = fopen("pFTs_gpu_CUDA.log","w");
  
  printf("alpha: %f\n",alpha);
  for (i=0; i<npart/*-(npart-2)*/ ; i++) {
    for (j=0; j<nrot/*-(nrot-2)*/ ; j++) {
      for (k=0; k<nk/*-(nk-2)*/ ; k++) {

        if ( debug_write == true ) {
          fprintf(pftFile,"%i %i %i %15.8f %15.8f %i %i %i %15.8f %15.8f %i %i %i %15.8f %15.8f\n",
                  i,j,k, cuCrealf(A[(j+nrot*k)*npart+i]),
                  cuCimagf(A[(j+nrot*k)*npart+i]),
                  i,j,k, cuCrealf(B[(j+nrot*k)*npart+i]),
                  cuCimagf(B[(j+nrot*k)*npart+i]),
                  i,j,k, cuCrealf(cuCstarf(B[(j+nrot*k)*npart+i])),
                  cuCimagf(cuCstarf(B[(j+nrot*k)*npart+i])) );
        }

        printf(ANSI_COLOR_GREEN"A[%i][%i][%i]=(%15.8f,%15.8f)"
               ANSI_COLOR_BLUE" B[%i][%i][%i]=(%15.8f,%15.8f)"
               ANSI_COLOR_YELLOW" Bstar[%i][%i][%i]=(%15.8f,%15.8f)\n"
               ANSI_COLOR_RESET,
               i,j,k, cuCrealf(A[(j+nrot*k)*npart+i]),
               cuCimagf(A[(j+nrot*k)*npart+i]),
               i,j,k, cuCrealf(B[(j+nrot*k)*npart+i]),
               cuCimagf(B[(j+nrot*k)*npart+i]),
               i,j,k, cuCrealf(cuCstarf(B[(j+nrot*k)*npart+i])),
               cuCimagf(cuCstarf(B[(j+nrot*k)*npart+i])) );
        C[(j+nrot*k)*npart+i] = make_cuFloatComplex ( 1.0, 0.0 );
      }
    }
  }
  
  fclose(pftFile);

  float sumr = 0.0;
  float sumA = 0.0;
  float sumB = 0.0;
  for (int ipart=0; ipart<npart ; ipart++) {
    for (int irot=0; irot<nrot ; irot++) {
      for (int ik=0; ik<nk ; ik++) {
        sumr += cuReCCstarmulf( A[(irot+nrot*ik)*npart+ipart], 
                                B[(irot+nrot*ik)*npart+ipart] ); 
        sumA += cuReCCstarf( A[(irot+nrot*ik)*npart+ipart]);
        sumB += cuReCCstarf( B[(irot+nrot*ik)*npart+ipart]);
      }
    }
  }
  printf("nvcc CPU sums:\n"
         ANSI_COLOR_GREEN"sum(Re(A*cong(B))[1:%i][1:%i][1:%i]))=%15.8f\n"
         ANSI_COLOR_GREEN"sum(Re(A*cong(A))[1:%i][1:%i][1:%i]))=%15.8f\n"
         ANSI_COLOR_GREEN"sum(Re(B*cong(B))[1:%i][1:%i][1:%i]))=%15.8f\n"
         ANSI_COLOR_RESET,
         npart, nrot, nk, sumr,
         npart, nrot, nk, sumA,
         npart, nrot, nk, sumB );

  return rc;

} /* end of print_function_header_N_N_info */

/* //////////////////////////////////////////////////////////////////////////// 
 -- Polarft_gpu kernel P_N helpers
*/

/* printing out the input values */ 
extern "C" int
print_function_header_P_N_info(polar_corr_calc_t *s_polar,
                              float *r,
                              cuFloatComplex *C,
                              const cuFloatComplex *A,
                              const cuFloatComplex *B,
                              int npart, int nrot, int nk,
                              float alpha)
{
  int rc = RC_SUCCESS;
  int i,j,k;
  float h_r = *r;

  int npts = npart * nrot * nk;

  printf(ANSI_COLOR_BRIGHT_YELLOW"h_r=%f, at Line %i %s\n"
         ANSI_COLOR_RESET,
         h_r,__LINE__,__FUNCTION__);
  
  printf(ANSI_COLOR_BRIGHT_YELLOW"s_polar->r_polar=%f, "
         ANSI_COLOR_BRIGHT_GREEN "s_polar->sumasq_polar=%f, "
         ANSI_COLOR_BRIGHT_BLUE  "s_polar->sumbsq_polar=%f, at Line %i %s\n"
         ANSI_COLOR_RESET,
         s_polar->r_polar, s_polar->sumasq_polar, s_polar->sumbsq_polar,
         __LINE__,__FUNCTION__);

  printf(ANSI_COLOR_BRIGHT_YELLOW"npart=%i, "
         ANSI_COLOR_BRIGHT_GREEN"nrot=%i, "
         ANSI_COLOR_BRIGHT_BLUE  "nk=%i, "
         ANSI_COLOR_BRIGHT_CYAN  "npts=%i, npart*nrot*nk=%i\n",
         npart, nrot, nk, npts, npart*nrot*nk);

  printf(ANSI_COLOR_BRIGHT_RED
         "npart*nrot*nk*sizeof(cuFloatComplex)=%lu, Memory=%.2f MB\n"
         ANSI_COLOR_RESET,
         npart*nrot*nk*sizeof(cuFloatComplex),
         npart*nrot*nk*sizeof(cuFloatComplex)/1.e6);

  FILE * pftFile;
  pftFile = fopen("pFTs_gpu_CUDA.log","a");
  
  printf("alpha: %f\n",alpha);
  for (i=0; i<npart/*-(npart-2)*/ ; i++) {
    for (j=0; j<nrot/*-(nrot-2)*/ ; j++) {
      for (k=0; k<nk/*-(nk-2)*/ ; k++) {

        if ( debug_write == true ) {
          fprintf(pftFile,"%i %i %i %15.8f %15.8f %i %i %i %15.8f %15.8f %i %i %i %15.8f %15.8f\n",
                  i,j,k, cuCrealf(A[(j+nrot*k)*npart+i]),
                  cuCimagf(A[(j+nrot*k)*npart+i]),
                  i,j,k, cuCrealf(B[(j+nrot*k)*npart+i]),
                  cuCimagf(B[(j+nrot*k)*npart+i]),
                  i,j,k, cuCrealf(cuCstarf(B[(j+nrot*k)*npart+i])),
                  cuCimagf(cuCstarf(B[(j+nrot*k)*npart+i])) );
        }

        /*
        printf(ANSI_COLOR_GREEN"A[%i][%i][%i]=(%15.8f,%15.8f)"
               ANSI_COLOR_BLUE" B[%i][%i][%i]=(%15.8f,%15.8f)"
               ANSI_COLOR_YELLOW" Bstar[%i][%i][%i]=(%15.8f,%15.8f)\n" 
               ANSI_COLOR_RESET,
               i,j,k, cuCrealf(A[(j+nrot*k)*npart+i]), 
               cuCimagf(A[(j+nrot*k)*npart+i]),
               i,j,k, cuCrealf(B[(j+nrot*k)*npart+i]), 
               cuCimagf(B[(j+nrot*k)*npart+i]),
               i,j,k, cuCrealf(cuCstarf(B[(j+nrot*k)*npart+i])), 
               cuCimagf(cuCstarf(B[(j+nrot*k)*npart+i])) );
        */
        //C[(j+nrot*k)*npart+i] = make_cuFloatComplex ( 1.0, 0.0 );
      }
    }
  }
  
  fclose(pftFile);

  float sumr = 0.0;
  float sumA = 0.0;
  float sumB = 0.0;
  for (int ipart=0; ipart<npart ; ipart++) {
    for (int irot=0; irot<nrot ; irot++) {
      for (int ik=0; ik<nk ; ik++) {
        sumr += cuReCCstarmulf( A[(irot+nrot*ik)*npart+ipart], 
                                B[(irot+nrot*ik)*npart+ipart] ); 
        sumA += cuReCCstarf( A[(irot+nrot*ik)*npart+ipart]);
        sumB += cuReCCstarf( B[(irot+nrot*ik)*npart+ipart]);
      }
    }
  }
  printf("nvcc CPU sums:\n"
         ANSI_COLOR_GREEN"sum(Re(A*cong(B))[1:%i][1:%i][1:%i]))=%15.8f\n"
         ANSI_COLOR_GREEN"sum(Re(A*cong(A))[1:%i][1:%i][1:%i]))=%15.8f\n"
         ANSI_COLOR_GREEN"sum(Re(B*cong(B))[1:%i][1:%i][1:%i]))=%15.8f\n"
         ANSI_COLOR_RESET,
         npart, nrot, nk, sumr,
         npart, nrot, nk, sumA,
         npart, nrot, nk, sumB );

  return rc;

} /* end of print_function_header_P_N_info */

/* //////////////////////////////////////////////////////////////////////////// 
 -- Polarft_gpu kernel X_N helpers
*/
/* printing out the input values */ 
extern "C" int
print_function_header_X_N_info(polar_corr_calc_t *s_polar,
                               float *r,
                               cuFloatComplex *C,
                               const cuFloatComplex *A,
                               const cuFloatComplex *B,
                               int npart, int nrot, int nk,
                               float alpha)
{
  int rc = RC_SUCCESS;
  int i,j,k;

  int npts = npart * nrot * nk;

  printf(ANSI_COLOR_BRIGHT_YELLOW"s_polar->r_polar=%f, "
         ANSI_COLOR_BRIGHT_GREEN "s_polar->sumasq_polar=%f, "
         ANSI_COLOR_BRIGHT_BLUE  "s_polar->sumbsq_polar=%f, at Line %i %s\n"
         ANSI_COLOR_RESET,
         s_polar->r_polar, s_polar->sumasq_polar, s_polar->sumbsq_polar,
         __LINE__,__FUNCTION__);

  printf(ANSI_COLOR_BRIGHT_YELLOW"npart=%i, "
         ANSI_COLOR_BRIGHT_GREEN"nrot=%i, "
         ANSI_COLOR_BRIGHT_BLUE  "nk=%i, "
         ANSI_COLOR_BRIGHT_CYAN  "npts=%i, npart*nrot*nk=%i\n",
         npart, nrot, nk, npts, npart*nrot*nk);

  printf(ANSI_COLOR_BRIGHT_RED
         "npart*nrot*nk*sizeof(cuFloatComplex)=%lu, Memory=%.2f MB\n"
         ANSI_COLOR_RESET,
         npart*nrot*nk*sizeof(cuFloatComplex),
         npart*nrot*nk*sizeof(cuFloatComplex)/1.e6);

  FILE * pftFile;
  pftFile = fopen("pFTs_gpu_CUDA.log","a");
  
  printf("alpha: %f\n",alpha);
  if ( debug_write == true ) {
    for (i=0; i<npart ; i++) {
      for (j=0; j<nrot ; j++) {
        for (k=0; k<nk ; k++) {

          fprintf(pftFile,"%i %i %i %15.8f %15.8f %i %i %i %15.8f %15.8f %i %i %i %15.8f %15.8f\n",
                  i,j,k, cuCrealf(A[(j+nrot*k)*npart+i]),
                  cuCimagf(A[(j+nrot*k)*npart+i]),
                  i,j,k, cuCrealf(B[(j+nrot*k)*npart+i]),
                  cuCimagf(B[(j+nrot*k)*npart+i]),
                  i,j,k, cuCrealf(cuCstarf(B[(j+nrot*k)*npart+i])),
                  cuCimagf(cuCstarf(B[(j+nrot*k)*npart+i])) );
        }
      }
    }
  }

  for (i=0; i<npart-(npart-2) ; i++) {
    for (j=0; j<nrot-(nrot-2) ; j++) {
      for (k=0; k<nk-(nk-2) ; k++) {
        printf(ANSI_COLOR_GREEN "A[%i][%i][%i]=(%15.8f,%15.8f) "
               ANSI_COLOR_BLUE  "B[%i][%i][%i]=(%15.8f,%15.8f) "
               ANSI_COLOR_YELLOW"Bstar[%i][%i][%i]=(%15.8f,%15.8f)\n"
               ANSI_COLOR_RESET,
               i,j,k, cuCrealf(A[(j+nrot*k)*npart+i]),
               cuCimagf(A[(j+nrot*k)*npart+i]),
               i,j,k, cuCrealf(B[(j+nrot*k)*npart+i]),
               cuCimagf(B[(j+nrot*k)*npart+i]),
               i,j,k, cuCrealf(cuCstarf(B[(j+nrot*k)*npart+i])),
               cuCimagf(cuCstarf(B[(j+nrot*k)*npart+i])) );
        //C[(j+nrot*k)*npart+i] = make_cuFloatComplex ( 1.0, 0.0 );
      }
    }
  }
  
  fclose(pftFile);

  float sumr = 0.0;
  float sumA = 0.0;
  float sumB = 0.0;
  for (int ipart=0; ipart<npart ; ipart++) {
    for (int irot=0; irot<nrot ; irot++) {
      for (int ik=0; ik<nk ; ik++) {
        sumr += cuReCCstarmulf( A[(irot+nrot*ik)*npart+ipart], 
                                B[(irot+nrot*ik)*npart+ipart] ); 
        sumA += cuReCCstarf( A[(irot+nrot*ik)*npart+ipart]);
        sumB += cuReCCstarf( B[(irot+nrot*ik)*npart+ipart]);
      }
    }
  }
  printf("nvcc CPU sums:\n"
         ANSI_COLOR_GREEN"sum(Re(A*cong(B))[1:%i][1:%i][1:%i]))=%15.8f\n"
         ANSI_COLOR_GREEN"sum(Re(A*cong(A))[1:%i][1:%i][1:%i]))=%15.8f\n"
         ANSI_COLOR_GREEN"sum(Re(B*cong(B))[1:%i][1:%i][1:%i]))=%15.8f\n"
         ANSI_COLOR_RESET,
         npart, nrot, nk, sumr,
         npart, nrot, nk, sumA,
         npart, nrot, nk, sumB );

  return rc;

} /* end of print_function_header_X_N_info */
/* printing out the input values */ 
extern "C" int
print_sumVecs_X_N(int npart, float *sum_vec, float *h_ptr_vec,
                  float *r, double *sum_dble_vec, int depth) {
  int rc = RC_SUCCESS;
  int i;
  printf(ANSI_COLOR_BRIGHT_RED"%22s%30s%22s%33s\n","CPU","CPU","GPU","CPU");
  for (i=0; i<npart-(npart-depth) ; i++) {
    printf(ANSI_COLOR_BRIGHT_CYAN    "cc vector: "
           ANSI_COLOR_BRIGHT_GREEN   "sum_vec[%i]=%15.8f, "
           ANSI_COLOR_BRIGHT_YELLOW  "h_ptr_vec[%i]=%15.8f, "
           ANSI_COLOR_BRIGHT_BLUE    "r[%i]=%15.8f, "
           ANSI_COLOR_BRIGHT_MAGENTA "sum_dble_vec[%i]=%15.8f\n"
           ANSI_COLOR_RESET,
           i,sum_vec[i],
           i,h_ptr_vec[i],
           i,r[i],
           i,sum_dble_vec[i]);
  }

  return rc;
}
/* printing out sanity check on the d{1,2}lim arrays */ 
extern "C" int
print_d12lim(int d1lim[], int d2lim[])
{
  int rc = RC_SUCCESS;

  printf(ANSI_COLOR_BRIGHT_CYAN"d{1,2}lim arrays: "
         ANSI_COLOR_BRIGHT_YELLOW"d1lim[0] = %i, d1lim[1] = %i, "
         ANSI_COLOR_BRIGHT_GREEN"d2lim[0] = %i, d2lim[1] = %i\n"
         ANSI_COLOR_RESET,
         d1lim[0], d1lim[1], d2lim[0], d2lim[1]);

  return rc;

}

/* //////////////////////////////////////////////////////////////////////////// 
 -- Polarft_gpu kernel Z_N helpers
*/
/* printing out the input values */ 
extern "C" int
print_function_header_Z_N_info(polar_corr_calc_t *s_polar,
                               float *r, float *cormat3d,
                               const cuFloatComplex *A,
                               const cuFloatComplex *B,
                               const float *sqsums_A,
                               const float *sqsums_B,
                               pfts_Sizes *p_pfts_Sizes,
                               int npart, int nrot, int nk,
                               float alpha, int depth) {
  int rc = RC_SUCCESS;
  
  rc = print_s_polar_struct(s_polar);
  rc = print_npart_nrot_nk(npart, nrot, nk);
  rc = print_A_and_B_depth(A, B, p_pfts_Sizes, npart, nrot, nk, depth);
  rc = print_sizeA_B_Z_N_info(A, B, p_pfts_Sizes);
  
  return rc;
}
/* //////////////////////////////////////////////////////////////////////////// 
 -- Polarft_gpu kernel C_N helpers
*/
/* printing out the input values */ 
extern "C" int
print_function_header_C_N_info(deviceDetails_t * s_devD,
                               polar_corr_calc_t *s_carte,
                               float *r,
                               cuFloatComplex *shmat,
                               const cuFloatComplex *A,
                               const cuFloatComplex *B,
                               img_2D_cart_Sizes *p_img_2D_cart_Sizes,
                               int vx, int vy, int vz,
                               float alpha, int depth) {
  int rc = RC_SUCCESS;
  
  rc = print_s_polar_struct(s_carte);
  rc = print_vx_vy_vz(vx,vy,vz);
  rc = print_A_and_B_V_depth(shmat,A,B,p_img_2D_cart_Sizes, vx, vy, vz, depth);
  rc = print_sizeA_B_C_N_info(A, B, p_img_2D_cart_Sizes);
  
  return rc;
}
/* printing out the first few elements of the matrices A and B */
extern "C" int
print_A_and_B_V_depth(cuFloatComplex *shmat,
                      const cuFloatComplex *A,
                      const cuFloatComplex *B,
                      img_2D_cart_Sizes *p_img_2D_cart_Sizes,
                      int vx, int vy, int vz,
                      int depth) {
  //
  //pfts_Sizes *p_pfts_Sizes

  int rc = RC_SUCCESS;
  int i,j,k;

  /* the loops go from */

  for (i=0; i<vx-(vx-depth) ; i++) {
    for (j=0; j<vy-(vy-depth) ; j++) {
      for (k=0; k<vz-(vz-1) ; k++) {
        printf(ANSI_COLOR_GREEN "A[%i][%i][%i]=(%15.8f,%15.8f) "
               ANSI_COLOR_BLUE  "B[%i][%i][%i]=(%15.8f,%15.8f) "
               ANSI_COLOR_YELLOW"Bstar[%i][%i][%i]=(%15.8f,%15.8f) "
               ANSI_COLOR_MAGENTA"shmat[%i][%i][%i]=(%15.8f,%15.8f)\n"
               ANSI_COLOR_RESET,
               i,j,k, cuCrealf(A[(j+vy*k)*vx+i]),
                      cuCimagf(A[(j+vy*k)*vx+i]),
               
               i,j,k, cuCrealf(B[(j+vy*k)*vx+i]),
                      cuCimagf(B[(j+vy*k)*vx+i]),

               i,j,k, cuCrealf(cuCstarf(B[(j+vy*k)*vx+i])),
                      cuCimagf(cuCstarf(B[(j+vy*k)*vx+i])),

               i,j,k, cuCrealf(shmat[(j+vy*k)*vx+i]),
                      cuCimagf(shmat[(j+vy*k)*vx+i])  );
      }
    }
  }

  /* the loops go from */
  if ( debug_write_AandB == true ) {

    FILE * AFile;
    AFile = fopen("carte2D_A_gpu_CUDA.log","a");

    for (i=0; i<  p_img_2D_cart_Sizes->get_2D_vx() ; i++) {
      for (j=0; j<  p_img_2D_cart_Sizes->get_2D_vy() ; j++) {
        for (k=0; k<  p_img_2D_cart_Sizes->get_2D_vz() ; k++) {
          fprintf(AFile,"A[%i][%i][%i]=(%15.8f,%15.8f) \n",
                  i,j,k, cuCrealf(A[(j+ p_img_2D_cart_Sizes->get_2D_vy()*k)* p_img_2D_cart_Sizes->get_2D_vx()+i]),
                         cuCimagf(A[(j+ p_img_2D_cart_Sizes->get_2D_vy()*k)* p_img_2D_cart_Sizes->get_2D_vx()+i]) );
          
        }
      }
    }
    fclose(AFile);

    FILE * BFile;
    BFile = fopen("carte2D_B_gpu_CUDA.log","a");

    for (i=0; i<  p_img_2D_cart_Sizes->get_2D_vx() ; i++) {
      for (j=0; j<  p_img_2D_cart_Sizes->get_2D_vy() ; j++) {
        for (k=0; k<  p_img_2D_cart_Sizes->get_2D_vz() ; k++) {
          fprintf(BFile,"B[%i][%i][%i]=(%15.8f,%15.8f) "
                 "Bstar[%i][%i][%i]=(%15.8f,%15.8f)\n",

                 i,j,k, cuCrealf(B[(j+p_img_2D_cart_Sizes->get_2D_vy()*k)*p_img_2D_cart_Sizes->get_2D_vx()+i]),
                        cuCimagf(B[(j+p_img_2D_cart_Sizes->get_2D_vy()*k)*p_img_2D_cart_Sizes->get_2D_vx()+i]),

                 i,j,k, cuCrealf(cuCstarf(B[(j+p_img_2D_cart_Sizes->get_2D_vy()*k)*p_img_2D_cart_Sizes->get_2D_vx()+i])),
                        cuCimagf(cuCstarf(B[(j+p_img_2D_cart_Sizes->get_2D_vy()*k)*p_img_2D_cart_Sizes->get_2D_vx()+i])) );
        }
      }
    }
    fclose(BFile);
    
  }
  
  return rc;
}
/* printing out the first few elements of the matrices A and B */
extern "C" int
print_cmat2sh_V_depth(const cuFloatComplex *d_cmat2sh,
                      img_2D_cart_Sizes *p_img_2D_cart_Sizes,
                      int vx, int vy, int vz,
                      int depth) {
  /*return code*/
  int rc = RC_SUCCESS;
  /* the error handlers from the cuda library */
  cudaError_t err;
  /* size of the element in consideration */
  int size_cmat2sh = p_img_2D_cart_Sizes->get_npts()*sizeof(cuFloatComplex);
  /*indexers*/
  int i,j,k;

  /*retrieving the matrix from device*/
  cuFloatComplex *h_cmat2sh = (cuFloatComplex*)malloc(size_cmat2sh);
  err = cudaMemcpy(h_cmat2sh, d_cmat2sh, size_cmat2sh, cudaMemcpyDeviceToHost);
  if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_corr_Hadmr_gpu(err); }
  for (i=0/*vx-depth*/; i<vx-(vx-depth) ; i++) {
    for (j=0/*vy-depth*/; j<vy-(vy-depth) ; j++) {
      for (k=0/*vz-1*/; k<vz-(vz-1) ; k++) {
        printf(ANSI_COLOR_GREEN "h_cmat2sh[%i][%i][%i]=(%15.8f,%15.8f)\n"
               ANSI_COLOR_RESET,
               i,j,k, cuCrealf(h_cmat2sh[(j+vy*k)*vx+i]),
               cuCimagf(h_cmat2sh[(j+vy*k)*vx+i]) );
      }
    }
  }
  /* Wrtitting the matrix to diesk for debug*/
  if ( debug_write_AandB == true ) {
    FILE * h_cmat2shFile;
    h_cmat2shFile = fopen("carte2D_h_cmat2sh_gpu_CUDA.log","a");
    for (i=0; i<  p_img_2D_cart_Sizes->get_2D_vx() ; i++) {
      for (j=0; j<  p_img_2D_cart_Sizes->get_2D_vy() ; j++) {
        for (k=0; k<  p_img_2D_cart_Sizes->get_2D_vz() ; k++) {
          fprintf(h_cmat2shFile,"h_cmat2sh[%i][%i][%i]=(%15.8f,%15.8f)\n",
                  i,j,k,
                  cuCrealf(h_cmat2sh
                           [(j+ p_img_2D_cart_Sizes->get_2D_vy()*k)*
                            p_img_2D_cart_Sizes->get_2D_vx()+i]),
                  cuCimagf(h_cmat2sh
                           [(j+ p_img_2D_cart_Sizes->get_2D_vy()*k)*
                            p_img_2D_cart_Sizes->get_2D_vx()+i]) );
        }
      }
    }
    /*closingh file */
    fclose(h_cmat2shFile);
  }
  /*freeing the memory*/
  free(h_cmat2sh);
  
  return rc;
}
/* printing out the first few elements of the matrices A and B */
extern "C" int
print_A_and_B_depth(const cuFloatComplex *A,
                    const cuFloatComplex *B,
                    pfts_Sizes *p_pfts_Sizes,
                    int npart, int nrot, int nk,
                    int depth) {
  int rc = RC_SUCCESS;
  int i,j,k;
  /* Here nhlfrot-->refsz, n2rot-->ptclsz, n2part-->2*ipart in SIMPLE */ 
  int nhlfrot = nrot / 2; // Here the actual size of the arrays for both
  int n2rot = nrot * 2;   // A and B are not [npart][nrot][nk] but
  int n2part = npart * 2; // A[npart][nhlfrot][nk] and B[n2part][n2rot][nk] 

  /* the loops go from */
  for (i=0; i<npart-(npart-depth) ; i++) {
    for (j=0; j<nrot-(nrot-depth) ; j++) {
      for (k=0; k<nk-(nk-depth) ; k++) {
        printf(ANSI_COLOR_GREEN "A[%i][%i][%i]=(%15.8f,%15.8f) "
               ANSI_COLOR_BLUE  "B[%i][%i][%i]=(%15.8f,%15.8f) "
               ANSI_COLOR_YELLOW"Bstar[%i][%i][%i]=(%15.8f,%15.8f)\n"
               ANSI_COLOR_RESET,
               i,j,k, cuCrealf(A[(j+nhlfrot*k)*npart+i]),
               cuCimagf(A[(j+nhlfrot*k)*npart+i]),

               i,j,k, cuCrealf(B[(j+n2rot*k)*n2part+i]),
               cuCimagf(B[(j+n2rot*k)*n2part+i]),

               i,j,k, cuCrealf(cuCstarf(B[(j+n2rot*k)*n2part+i])),
               cuCimagf(cuCstarf(B[(j+n2rot*k)*n2part+i])) );
      }
    }
  }

  /* the loops go from */

  /* the loops go from */
  if ( debug_write_AandB == true ) {

    FILE * AFile;
    AFile = fopen("pFTs_A_gpu_CUDA.log","a");

    for (i=0; i<  p_pfts_Sizes->get_npart() ; i++) {
      for (j=0; j<  p_pfts_Sizes->get_nhlfrot() ; j++) {
        for (k=0; k<  p_pfts_Sizes->get_nk() ; k++) {
          fprintf(AFile,"A[%i][%i][%i]=(%15.8f,%15.8f) \n",
                  i,j,k, cuCrealf(A[(j+ p_pfts_Sizes->get_nhlfrot()*k)* p_pfts_Sizes->get_npart()+i]),
                  cuCimagf(A[(j+ p_pfts_Sizes->get_nhlfrot()*k)* p_pfts_Sizes->get_npart()+i]) );
          
        }
      }
    }
    fclose(AFile);

    FILE * BFile;
    BFile = fopen("pFTs_B_gpu_CUDA.log","a");

    for (i=0; i<  p_pfts_Sizes->get_n2part() ; i++) {
      for (j=0; j<  p_pfts_Sizes->get_n2rot() ; j++) {
        for (k=0; k<  p_pfts_Sizes->get_nk() ; k++) {
          fprintf(BFile,"B[%i][%i][%i]=(%15.8f,%15.8f) "
                 "Bstar[%i][%i][%i]=(%15.8f,%15.8f)\n",

                 i,j,k, cuCrealf(B[(j+p_pfts_Sizes->get_n2rot()*k)*n2part+i]),
                 cuCimagf(B[(j+p_pfts_Sizes->get_n2rot()*k)*n2part+i]),

                 i,j,k, cuCrealf(cuCstarf(B[(j+p_pfts_Sizes->get_n2rot()*k)*p_pfts_Sizes->get_n2part()+i])),
                 cuCimagf(cuCstarf(B[(j+p_pfts_Sizes->get_n2rot()*k)*p_pfts_Sizes->get_n2part()+i])) );
        }
      }
    }
    fclose(BFile);
    
  }
  
  return rc;
}
/* printing out to file the hadamard product C */
extern "C" int
print_C_HadrProd(const float *C,
                 pfts_Sizes *p_pfts_Sizes) {

  int rc = RC_SUCCESS;
  int i,j,k;
  int jpart,jrot,jk;

  /* the loops go from */
  FILE * hdrFile;
  hdrFile = fopen("pFTs_Hdr_gpu_CUDA.log","w");
  i = 0;
  for (jpart = 0 ; jpart < p_pfts_Sizes->get_npart() ; jpart++ ) {
    j = 0;
    for (jrot = 0 ; jrot < p_pfts_Sizes->get_nhlfrot() ; jrot++ ) {
      k = 0;
      for ( jk = 0 ; jk < p_pfts_Sizes->get_nk() ; jk++ ) {
        fprintf(hdrFile,"C[%i][%i][%i]=%15.8f \n",
               i,j,k, C[ (j + p_pfts_Sizes->get_nhlfrot()*k) * p_pfts_Sizes->get_npart() + i] );
        k++;
      }
      j++;
    }
    i++;
  }
  fclose(hdrFile);

  return rc;
}
/* printing out to file the corrmat3D from the cpu calculation */
extern "C" int
print_C_cormat3d(const float *cormat3d,
                 pfts_Sizes *p_pfts_Sizes, const char *filename) {

  int rc = RC_SUCCESS;
  int i,j,k;
  int jpart,jrot,jk;

  /* the loops go from */
  FILE * c3dFilecpu;
  c3dFilecpu = fopen(filename,"w");
  rewind(c3dFilecpu);
  char buff[50];
  i=0;
  for (jpart = 0 ; jpart < p_pfts_Sizes->get_npart() ; jpart++ ) {
    j = 0;
    for (jrot = 0 ; jrot < p_pfts_Sizes->get_npart() ; jrot++ ) {
      k = 0;
      for ( jk = 0 ; jk < p_pfts_Sizes->get_nk() ; jk++ ) {
        sprintf(buff,"cormat3d[%i][%i][%i]=%15.8f \n",i,j,k,
                cormat3d[ (j + p_pfts_Sizes->get_npart()*k) *
                          p_pfts_Sizes->get_npart() + i]);
        fwrite(buff,strlen(buff),1,c3dFilecpu);
        k++;
      }
      j++;
    }
    i++;
  }
  fclose(c3dFilecpu);

  return rc;
}
/* printing out to file the hadamard product C */
extern "C" int
print_dC_HadrProd(float *d_C,
                  pfts_Sizes *p_pfts_Sizes) {

  int rc = RC_SUCCESS;
  int i,j,k;
  int jpart,jrot,jk;
  /* indexers */
  int c_indx;
  /* debugging variables */
  float *C = NULL;
  int size_C = p_pfts_Sizes->get_npts_A() * sizeof(float);
  /* the error handlers from the cuda library */
  cudaError_t err;

  /* allocating the resources for Real(A*Bstar) */
  C = (float*)malloc(size_C);
  i = 0;
  for (jpart = 0 ; jpart < p_pfts_Sizes->get_npart() ; jpart++ ) {
    j = 0;
    for (jrot = 0 ; jrot < p_pfts_Sizes->get_nhlfrot() ; jrot++ ) {
      k = 0;
      for (jk = 0 ; jk < p_pfts_Sizes->get_nk() ; jk++ ) {
        c_indx = (j + p_pfts_Sizes->get_nhlfrot()*k) * p_pfts_Sizes->get_npart() + i;
        C[c_indx] = 0.0;
        k++;
      }
      j++;
    }
    i++;
  }

  err = cudaMemcpy(C, d_C, size_C, cudaMemcpyDeviceToHost);
  if ( (int)err != CUDA_SUCCESS ) {rc = get_error_id_corr_Hadmr_gpu(err); }
  rc = print_C_HadrProd(C, p_pfts_Sizes);
  free(C);
 
  return rc;
}
/* printing out the true dimension of the problem */ 
extern "C" int
print_sizeA_B_Z_N_info(const cuFloatComplex *A,
                       const cuFloatComplex *B,
                       pfts_Sizes *p_pfts_Sizes) {
  int rc = RC_SUCCESS;

  int npart = p_pfts_Sizes->get_npart();
  int  nrot = p_pfts_Sizes->get_nrot();
  int    nk = p_pfts_Sizes->get_nk();
  
  /* Here nhlfrot-->refsz, n2rot-->ptclsz, n2part-->2*ipart in SIMPLE
   * nhlfrot = nrot / 2; Here the actual size of the arrays for both
   *   n2rot = nrot * 2; A and B are not [npart][nrot][nk] but
   *  n2part = npart* 2; A[npart][nhlfrot][nk] and B[n2part][n2rot][nk]
  */ 
  int nhlfrot   = p_pfts_Sizes->get_nhlfrot();
  int   n2rot   = p_pfts_Sizes->get_n2rot();
  int  n2part   = p_pfts_Sizes->get_n2part();

  int npts      = p_pfts_Sizes->get_npts();  //npart * nrot    * nk;
  int npts_A    = p_pfts_Sizes->get_npts_A();//npart * nhlfrot * nk;
  int npts_B    = p_pfts_Sizes->get_npts_B();//n2part* n2rot   * nk;
  int npts_cr3D = p_pfts_Sizes->get_npts_cr3D();//npart**2 * nrot;

  printf(ANSI_COLOR_BRIGHT_CYAN"True            : "
         ANSI_COLOR_BRIGHT_YELLOW" npart=%i, "
         ANSI_COLOR_BRIGHT_GREEN"   nrot=%i, "
         ANSI_COLOR_BRIGHT_BLUE  "nk=%i, "
         ANSI_COLOR_BRIGHT_CYAN  "  npts=%i,   npart*nrot*nk=%i\n"
         ANSI_COLOR_RESET,
         npart, nrot, nk, npts, npart*nrot*nk);

  printf(ANSI_COLOR_BRIGHT_CYAN"PFT_ref   pft1 A: "
         ANSI_COLOR_BRIGHT_YELLOW" npart=%i, "
         ANSI_COLOR_BRIGHT_GREEN"nhlfrot=%i, "
         ANSI_COLOR_BRIGHT_BLUE  "nk=%i, "
         ANSI_COLOR_BRIGHT_CYAN  "npts_A=%i, npart*nhlfrot*nk=%i\n"
         ANSI_COLOR_RESET,
         npart, nhlfrot, nk, npts_A, npart*nhlfrot*nk);

  printf(ANSI_COLOR_BRIGHT_CYAN"PFT_ptcls pft2 B: "
         ANSI_COLOR_BRIGHT_YELLOW"n2part=%i, "
         ANSI_COLOR_BRIGHT_GREEN" n2rot=%i, "
         ANSI_COLOR_BRIGHT_BLUE  "nk=%i, "
         ANSI_COLOR_BRIGHT_CYAN  "npts_B=%i,  npart*n2rot*nk=%i\n"
         ANSI_COLOR_RESET,
         n2part, n2rot, nk, npts_B, n2part*n2rot*nk);
  
  printf(ANSI_COLOR_BRIGHT_RED
         "Total true size :    npart*nrot*nk*sizeof(cuFloatComplex) = %lu, "
         " Memory = %.2f MB\n"
         ANSI_COLOR_RESET,
         npart*nrot*nk*sizeof(cuFloatComplex),
         npart*nrot*nk*sizeof(cuFloatComplex)/1.e6);

  printf(ANSI_COLOR_BRIGHT_RED
         "Size of A       : npart*nhlfrot*nk*sizeof(cuFloatComplex) = %lu, "
         " Memory = %.2f MB\n"
         ANSI_COLOR_RESET,
         npts_A*sizeof(cuFloatComplex),
         npts_A*sizeof(cuFloatComplex)/1.e6);

  printf(ANSI_COLOR_BRIGHT_RED
         "Size of B       :  n2part*n2rot*nk*sizeof(cuFloatComplex) = %lu, "
         "Memory = %.2f MB\n"
         ANSI_COLOR_RESET,
         npts_B*sizeof(cuFloatComplex),
         npts_B*sizeof(cuFloatComplex)/1.e6);
  
  printf(ANSI_COLOR_BRIGHT_RED
         "--------------------------------> Ratio of Sizes of B / A : %i\n"
         ANSI_COLOR_RESET,
         (npts_B)/(npts_A));

  printf(ANSI_COLOR_BRIGHT_RED
         "Size of cormat3D: npart*npart*nrot*sizeof(float)          = %lu, "
         "Memory = %.2f MB\n"
         ANSI_COLOR_RESET,
         npts_cr3D*sizeof(cuFloatComplex),
         npts_cr3D*sizeof(cuFloatComplex)/1.e6);

  return rc;
}
/* printing out the true dimension of the problem */ 
extern "C" int
print_sizeA_B_C_N_info(const cuFloatComplex *A,
                       const cuFloatComplex *B,
                       img_2D_cart_Sizes *p_img_2D_cart_Sizes) {
  //  pfts_Sizes *p_pfts_Sizes
  int rc = RC_SUCCESS;

  int vx = p_img_2D_cart_Sizes->get_2D_vx();
  int vy = p_img_2D_cart_Sizes->get_2D_vy();
  int vz = p_img_2D_cart_Sizes->get_2D_vz();

  int npts      = p_img_2D_cart_Sizes->get_npts();  //vx * vy * vz;
  int npts_A    = p_img_2D_cart_Sizes->get_npts_A();//vx * vy * vz;
  int npts_B    = p_img_2D_cart_Sizes->get_npts_B();//vx*  vy * vz;
  //  int npts_cr3D = p_img_2D_cart_Sizes->get_npts_cr3D();//vx**2 * vy;

  printf(ANSI_COLOR_BRIGHT_CYAN"True            : "
         ANSI_COLOR_BRIGHT_YELLOW"vx=%i, "
         ANSI_COLOR_BRIGHT_GREEN "vy=%i, "
         ANSI_COLOR_BRIGHT_BLUE  "vz=%i, "
         ANSI_COLOR_BRIGHT_CYAN  "  npts=%i, vx*vy*vz=%i\n"
         ANSI_COLOR_RESET,
         vx, vy, vz, npts, vx*vy*vz);

  printf(ANSI_COLOR_BRIGHT_CYAN"PFT_ref   pft1 A: "
         ANSI_COLOR_BRIGHT_YELLOW"vx=%i, "
         ANSI_COLOR_BRIGHT_GREEN "vy=%i, "
         ANSI_COLOR_BRIGHT_BLUE  "vz=%i, "
         ANSI_COLOR_BRIGHT_CYAN  "npts_A=%i, vx*vy*vz=%i\n"
         ANSI_COLOR_RESET,
         vx, vy, vz, npts_A, vx*vy*vz);

  printf(ANSI_COLOR_BRIGHT_CYAN"PFT_ptcls pft2 B: "
         ANSI_COLOR_BRIGHT_YELLOW"vx=%i, "
         ANSI_COLOR_BRIGHT_GREEN "vy=%i, "
         ANSI_COLOR_BRIGHT_BLUE  "vz=%i, "
         ANSI_COLOR_BRIGHT_CYAN  "npts_B=%i, vx*vy*vz=%i\n"
         ANSI_COLOR_RESET,
         vx, vy, vz, npts_B, vx*vy*vz);
  
  printf(ANSI_COLOR_BRIGHT_RED
         "Total true size : vx*vy*vz*sizeof(cuFloatComplex) = %lu, "
         "Memory = %.2f MB\n"
         ANSI_COLOR_RESET,
         vx*vy*vz*sizeof(cuFloatComplex),
         vx*vy*vz*sizeof(cuFloatComplex)/1.e6);

  printf(ANSI_COLOR_BRIGHT_RED
         "Size of A       : vx*vy*vz*sizeof(cuFloatComplex) = %lu, "
         "Memory = %.2f MB\n"
         ANSI_COLOR_RESET,
         npts_A*sizeof(cuFloatComplex),
         npts_A*sizeof(cuFloatComplex)/1.e6);

  printf(ANSI_COLOR_BRIGHT_RED
         "Size of B       : vx*vy*vz*sizeof(cuFloatComplex) = %lu, "
         "Memory = %.2f MB\n"
         ANSI_COLOR_RESET,
         npts_B*sizeof(cuFloatComplex),
         npts_B*sizeof(cuFloatComplex)/1.e6);
  
  printf(ANSI_COLOR_BRIGHT_RED
         "------------------------> Ratio of Sizes of B / A : %i\n"
         ANSI_COLOR_RESET,
         (npts_B)/(npts_A));

  return rc;
}
/* printing out the true dimension of the problem */ 
extern "C" int
polarft_gencorrAll_Z_N_cpu(polar_corr_calc_t *s_polar,
                           float *r, float *cormat3d,
                           const cuFloatComplex *A,
                           const cuFloatComplex *B,
                           const float *sqsums_A,
                           const float *sqsums_B,
                           pfts_Sizes *p_pfts_Sizes,
                           int npart, int nrot, int nk,
                           float alpha)
{
  int rc = RC_SUCCESS;
  /* array dimension */
  int d1lim[2], d2lim[2];
  /* creating the work space for the CPU */
  float *C = NULL;
  /* sizes of the element in consideration */
  int size_C = p_pfts_Sizes->get_npts_A() * sizeof(float);
  //now summing the over the chunk and mapping it to vector r  
  int chunk = p_pfts_Sizes->get_nhlfrot() *
              p_pfts_Sizes->get_nk()      ;  
  /* indexer */
  int ipart, irot;               // main loop
  int jpart,jrot,jk;             //for internal sums
  int i,j,k;                     //movers
  int a_indx, b_indx, c_indx;    //locators
  int cr3d_indx;
  int sqB_indx;
  /* printer index depth */
  float depth = 2;

#if defined (OPENMP) /*preprossing for the OpenMP environment */
  omp_set_num_threads(nthr);
#endif
  
  /* allocating the resources for Real(A*Bstar) */
  C = (float*)malloc(size_C);
  
  if (debug == true ) {
    rc = print_function_header_Z_N_info(s_polar,
                                        r, cormat3d,
                                        A, B, 
                                        sqsums_A, sqsums_B,
                                        p_pfts_Sizes,
                                        npart, nrot, nk,
                                        alpha, depth);
  }

  /* the pft cross correlator calculation */
  /* initilizing the dim arrays */
  for (i = 0 ; i < 2 ; i++) {d1lim[i] = 0;} 
  for (i = 0 ; i < 2 ; i++) {d2lim[i] = 0;}
  /* start of the loop ipart irot */
  for (ipart = 0 ; ipart < p_pfts_Sizes->get_npart() ; ipart++) {
    d1lim[0] = ipart;
    d1lim[1] = ipart + p_pfts_Sizes->get_npart();//-1;
    for (irot = 0 ; irot < p_pfts_Sizes->get_nrot() ; irot++) {
      d2lim[0] = irot;
      d2lim[1] = irot + p_pfts_Sizes->get_nhlfrot();//-1;

      i = 0;
      //#pragma omp parallel default(shared) private (j,k,jpart,jrot,jk)
      //{
      //#pragma omp for schedule(auto)

      for (jpart = d1lim[0] ; jpart < d1lim[1] ; jpart++ ) {
        j = 0;
        for (jrot = d2lim[0] ; jrot < d2lim[1] ; jrot++ ) {
          k = 0;
          for ( jk = 0 ; jk < p_pfts_Sizes->get_nk() ; jk++ ) {
            
            a_indx = (j    + p_pfts_Sizes->get_nhlfrot() * k  ) * p_pfts_Sizes->get_npart() + i;
            c_indx = (j    + p_pfts_Sizes->get_nhlfrot() * k  ) * p_pfts_Sizes->get_npart() + i;
            b_indx = (jrot +   p_pfts_Sizes->get_n2rot() * jk ) * p_pfts_Sizes->get_n2part() + jpart;

            C[c_indx] = cuReCCstarmulf(A[a_indx],
                                       B[b_indx] );
            k++;
          }
          j++;
        }
        i++;
      }

      //Start loop here
      //#pragma omp for schedule(auto)
      for (jpart = 0 ; jpart < p_pfts_Sizes->get_npart() ; jpart++ ) {
        float *startptr = C + jpart;
        r[jpart] = 0.0;
        //looping over the partial chunks and summing.
        for (int ichunk = 0 ; ichunk < chunk ; ichunk++ ) {
          float *endptr = startptr + ichunk * p_pfts_Sizes->get_npart();
          r[jpart] += *endptr;
        }
        cr3d_indx = (ipart+p_pfts_Sizes->get_npart()*irot) *
          p_pfts_Sizes->get_npart()+jpart;
        sqB_indx = d1lim[0]+jpart;
        cormat3d[cr3d_indx] = r[jpart] /
          sqrt(sqsums_A[jpart] * sqsums_B[sqB_indx]);
      }//finish the loop here.

      
      //} /* end of omp parallel loop */
      
    }
  } /* end of the loop ipart irot */

  if (debug_write_C == true){
    rc = print_C_HadrProd(C, p_pfts_Sizes);
    printf("printing before: print_C_cormat3d, line: %i, in %s\n", __LINE__,__FUNCTION__);
    const char *filename_cpu = "cormat3d_cpu_CUDA.log";
    rc = print_C_cormat3d(cormat3d, p_pfts_Sizes, filename_cpu);
  }

  /* sanity check on the d{1,2}lim arrays */
  rc = print_d12lim(d1lim, d2lim);
  
  /* freeing the resources on host CPU */
  free(C);
  
  return rc;
}

/* the aliases for external access */
extern "C" int print_s_devd_struct_() __attribute__((weak,alias("print_s_devD_struct")));

#endif /* CUDA */
