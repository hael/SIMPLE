/*
 *   -- Polarft kernel addon for polarft calculation on GPU 
 *      Author: Frederic Bonnet, Date: 27th April 2015
 *      Monash University
 *      April 2015
 *
 * @precisions normal z -> s d c
 */
#include "simple.h"
#include "get_fft123D_gpu.h"
#include "simple_cuDoubleComplex.h"

#ifndef _POLARFT_GPU_H_
#define _POLARFT_GPU_H_

#ifdef __cplusplus
extern "C" {
#endif


/* //////////////////////////////////////////////////////////////////////////// 
 -- Polarft_gpu kernel function definitions / Data on GPU
*/

  void zz2dgemm_ElmtWs_tesla_gpu_( char TRANSA, char TRANSB,
				   int m, int n, int k, 
				   double alpha, 
				   const cuDoubleComplex *A, int lda, 
				   const cuDoubleComplex *B, int ldb,
				   double beta, 
				   double *C, int ldc);
  
  void zz2dgemm_ElmtWs_tesla_sumsq_gpu_( char TRANSA, char TRANSB,
					 int m, int n, int k,
					 double alpha, 
					 const cuDoubleComplex *A, int lda, 
					 const cuDoubleComplex *B, int ldb,
					 double beta,
					 double *C, int ldc);
  /* GETTERS */  
  int get_polarft_corr_gpu(deviceDetails_t * s_devD,
                           polar_corr_calc_t *s_polar,
                           char *TRANSA, char *TRANSB,
                           float *r,
                           const float complex *A_in,
                           const float complex *B_in,
                           int *npart, int *nrot, int *nk,
                           int *lda, int *ldb, int *ldc,
                           float *alpha,
                           bench_t *s_bench, debug_gpu_t *s_debug_gpu);
  
 int get_polarft_gencorrall_gpu(deviceDetails_t * s_devD,
                                polar_corr_calc_t *s_polar,
                                char *TRANSA, char *TRANSB,
                                float *r, float *cormat3d_in,
                                const float complex *A_in, 
                                const float complex *B_in,
                                const float *sqsums_A_in,
                                const float *sqsums_B_in,
                                int *npart, int *nrot, int *nk,
                                int *lda, int *ldb, int *ldc,
                                float *alpha,
                                bench_t *s_bench, debug_gpu_t *s_debug_gpu);

  int get_polarft_multi_GPUs_gpu(deviceDetails_t * s_devD,
                                 polar_corr_calc_t *s_polar,
                                 char *TRANSA, char *TRANSB,
                                 float *r, float *cormat3d_in,
                                 const float complex *A_in, 
                                 const float complex *B_in,
                                 const float *sqsums_A_in,
                                 const float *sqsums_B_in,
                                 int *npart, int *nrot, int *nk,
                                 int *lda, int *ldb, int *ldc,
                                 float *alpha,
                                 bench_t *s_bench, debug_gpu_t *s_debug_gpu);

 int get_carte2d_ftExt_corr_gpu(deviceDetails_t * s_devD,
                                 polar_corr_calc_t *s_carte,
                                 char *TRANSA, char *TRANSB,
                                 float *r,
                                 float complex *shmat_in,
                                 const float complex *A_in, 
                                 const float complex *B_in,
                                 int *vx, int *vy, int *vz,
                                 int *lda, int *ldb, int *ldc,
                                 float *alpha,
                                 bench_t *s_bench, debug_gpu_t *s_debug_gpu);
  /* ACTIVATORS */  
  int polarft_corr_Hadmr_gpu_(deviceDetails_t * s_devD,
                              polar_corr_calc_t *s_polar,
                              char TRANSA, char TRANSB,
                              float *r,
                              const cuFloatComplex *A,
                              const cuFloatComplex *B,
                              int npart, int nrot, int nk,
                              int lda, int ldb, int ldc,
                              float alpha,
                              bench_t *s_bench, debug_gpu_t *s_debug_gpu);

  int polarft_gencorrAll_gpu_(deviceDetails_t * s_devD,
                              polar_corr_calc_t *s_polar,
                              char TRANSA, char TRANSB,
                              float *r, float *cormat3d,
                              const cuFloatComplex *A,
                              const cuFloatComplex *B,
                              const float *sqsums_A,
                              const float *sqsums_B,
                              int npart, int nrot, int nk,
                              int lda, int ldb, int ldc,
                              float alpha,
                              bench_t *s_bench, debug_gpu_t *s_debug_gpu);

  int polarft_multi_GPUs_gpu_(deviceDetails_t * s_devD,
                              polar_corr_calc_t *s_polar,
                              char TRANSA, char TRANSB,
                              float *r, float *cormat3d,
                              const cuFloatComplex *A,
                              const cuFloatComplex *B,
                              const float *sqsums_A,
                              const float *sqsums_B,
                              int npart, int nrot, int nk,
                              int lda, int ldb, int ldc,
                              float alpha,
                              bench_t *s_bench, debug_gpu_t *s_debug_gpu);

  int carte2d_ftExt_corr_gpu_(deviceDetails_t * s_devD,
                              polar_corr_calc_t *s_carte,
                              char TRANSA, char TRANSB,
                              float *r,
                              cuFloatComplex *shmat,
                              const cuFloatComplex *A,
                              const cuFloatComplex *B,
                              int vx, int vy, int vz,
                              int lda, int ldb, int ldc,
                              float alpha,
                              bench_t *s_bench, debug_gpu_t *s_debug_gpu);

  /*OPTIMSATION*/
  int get_polarft_krnl_Opti_gpu(deviceDetails_t * s_devD,
                                polar_corr_calc_t *s_polar,
                                char *TRANSA, char *TRANSB,
                                float *r, float *cormat3d_in,
                                const float complex *A_in, 
                                const float complex *B_in,
                                 const float *sqsums_A_in,
                                const float *sqsums_B_in,
                                int *npart, int *nrot, int *nk,
                                int *lda, int *ldb, int *ldc,
                                float *alpha,
                                bench_t *s_bench, debug_gpu_t *s_debug_gpu);
  int polarft_krnl_Opti_gpu_(deviceDetails_t * s_devD,
                             polar_corr_calc_t *s_polar,
                             char TRANSA, char TRANSB,
                             float *r, float *cormat3d,
                             const cuFloatComplex *A,
                             const cuFloatComplex *B,
                             const float *sqsums_A,
                             const float *sqsums_B,
                             int npart, int nrot, int nk,
                             int lda, int ldb, int ldc,
                             float alpha,
                             bench_t *s_bench, debug_gpu_t *s_debug_gpu);
  
  /* helper methods from CUDA */
  /* class for the 2D_cart_Sizes */
  class img_2D_cart_Sizes {
  private:
  public:
    /*global variables */
    int vx_in;
    int vy_in;
    int vz_in;
    /*constructor*/
    img_2D_cart_Sizes(int vx_in, int vy_in, int vz_in);
    //setters
    void set_2D_vx(int vx);
    void set_2D_vy(int vy);
    void set_2D_vz(int vz);
    //getters
    int get_2D_vx();
    int get_2D_vy();
    int get_2D_vz();
    int get_npts();
    int get_npts_A();
    int get_npts_B();
    /*destructors*/
    ~img_2D_cart_Sizes();
  };
  /* class for the pfts_Sizes */
  class pfts_Sizes {
  private:
  public:
    /*global variables*/
    int npart_in;
    int nrot_in;
    int nk_in;
    /*constructor*/
    pfts_Sizes(int npart_in, int nrot_in, int nk_in);
    /*setters*/
    void set_npart(int npart);
    void  set_nrot(int nrot);
    void    set_nk(int nk);
    /*getters*/
    int   get_npart();   int get_nrot();     int get_nk();
    int get_nhlfrot();  int get_n2rot(); int get_n2part();
    int    get_npts(); int get_npts_A(); int get_npts_B();
    int get_npts_cr3D();
    /*destructors*/
    ~pfts_Sizes();
  };
  /* class for the 3D mesh used on the 3D data structure on device */
  class mesh_3D {
  private:
  public:
    /*global variables*/
    polar_corr_calc *s_polar_in;
    int n_particles_in;
    int n_rotations_in;
    int n_rezstions_in;
    /*constructor*/
    mesh_3D(polar_corr_calc *s_polar_in); /*overloading the constructor*/
    mesh_3D(polar_corr_calc *s_polar_in,
            int n_particles_in, int n_rotations_in, int n_rezstions_in);
    /*setters*/
    void set_nx(polar_corr_calc *s_polar);
    void set_ny(polar_corr_calc *s_polar);
    void set_nz(polar_corr_calc *s_polar);
    /*getters*/
    int get_mesh3D_nx(); int get_mesh3D_ny(); int get_mesh3D_nz();
    int get_mesh3D_gridx();
    int get_mesh3D_gridy();
    int get_mesh3D_gridz();
    /*destructor*/
    ~mesh_3D();
  };
  /* class for the 3D mesh used on the 3D data structure on device */
  class mesh_1D {
  private:
  public:
    /*global variables */
    polar_corr_calc *s_polar_in;
    int n_particles_in;
    int n_rotations_in;
    int n_rezstions_in;
    /*constructor*/
    mesh_1D(polar_corr_calc *s_polar_in);
    mesh_1D(polar_corr_calc *s_polar_in,/* overloading the constructors*/
            int n_particles_in, int n_rotations_in, int n_rezstions_in);
    /*setters*/
    void set_nx(polar_corr_calc *s_polar);
    void set_ny(polar_corr_calc *s_polar);
    void set_nz(polar_corr_calc *s_polar);
    /*getters*/
    int get_mesh1D_nx();
    int get_mesh1D_ny();
    int get_mesh1D_nz();
    int get_mesh1D_chunk();
    int get_mesh1D_threadsPerBlock();
    int get_mesh1D_blocksPerGrid();
    int get_mesh1D_size_pB();
    /*destructor*/
    ~mesh_1D();
  };
  /* class for the 1DV mesh used on the 3D data structure on device */
  class mesh_1DV {
  private:
  public:
    /*global variables */
    polar_corr_calc *s_carte_in;
    int vx_in;
    int vy_in;
    int vz_in;
    /*constructor*/
    mesh_1DV(polar_corr_calc *s_carte_in);
    mesh_1DV(polar_corr_calc *s_carte_in,/* overloading the constructors*/
                  int vx_in, int vy_in, int vz_in);
    /*setters*/
    void set_nx(polar_corr_calc *s_carte);
    void set_ny(polar_corr_calc *s_carte);
    void set_nz(polar_corr_calc *s_carte);
    /*getters*/
    int get_mesh1DV_nx();
    int get_mesh1DV_ny();
    int get_mesh1DV_nz();
    int get_mesh1DV_chunk();
    int get_mesh1DV_threadsPerBlock();
    int get_mesh1DV_blocksPerGrid();
    int get_mesh1DV_size_pB();
    /*the destructor*/
    ~mesh_1DV();
  };
  /*error and warning handlers methods */
  int get_warning_message_corr_Hadmr_gpu();
  int get_error_id_corr_Hadmr_gpu(cudaError_t error_id);
  int get_warning_threadsPerBlock(polar_corr_calc_t *s_polar,
                                  int threadsPerBlock);
  /*checkers*/
  int check_grid3D(int gridx, int gridy, int gridz, polar_corr_calc_t *s_polar,
                   pfts_Sizes *p_pfts_Sizes);
  int check_grid3D_V(int gridx, int gridy, int gridz,
                     polar_corr_calc_t *s_carte,
                     img_2D_cart_Sizes *p_img_2D_cart_Sizes);
/*printers*/
  int print_summing_method(int isum);
  int print_1D_mesh(int s_cnstr_mesh_1D, mesh_1D *p_mesh_1D,
                    int N, int threadsPerBlock, int blocksPerGrid);

  int print_3D_mesh(int s_cnstr_mesh_3D, mesh_3D *p_mesh_3D,
                    polar_corr_calc_t *s_polar, pfts_Sizes *p_pfts_Sizes,
                    int nx, int ny, int nz);
  int print_3D_V_mesh(int s_cnstr_mesh_3D, mesh_3D *p_mesh_3D,
                      polar_corr_calc_t *s_carte,
                      img_2D_cart_Sizes *p_img_2D_cart_Sizes,
                      int nx, int ny, int nz);
  int print_s_polar_struct(polar_corr_calc_t *s_polar);
  int print_s_debug_struct(debug_gpu_t *s_debug_gpu);
  int print_s_bench_struct(bench_t *s_bench);
  int print_s_devD_struct(deviceDetails_t * s_devD);
  int print_npart_nrot_nk(int npart, int nrot, int nk);
  int print_vx_vy_vz(int vx, int vy, int vz);
  int print_d12lim(int d1lim[], int d2lim[]);
  int print_iKernel_threadsPerBlock(polar_corr_calc_t *s_polar,
                                    int N, int threadsPerBlock,
                                    int blocksPerGrid);
  int print_function_header_N_N_info(float *r,
                                     cuFloatComplex *C,
                                     const cuFloatComplex *A,
                                     const cuFloatComplex *B,
                                     int npart, int nrot, int nk,
                                     float alpha);
  int print_function_header_P_N_info(polar_corr_calc_t *s_polar,
                                     float *r,
                                     cuFloatComplex *C,
                                     const cuFloatComplex *A,
                                     const cuFloatComplex *B,
                                     int npart, int nrot, int nk,
                                     float alpha);
  int print_function_header_X_N_info(polar_corr_calc_t *s_polar,
                                     float *r,
                                     cuFloatComplex *C,
                                     const cuFloatComplex *A,
                                     const cuFloatComplex *B,
                                     int npart, int nrot, int nk,
                                     float alpha);
  int print_function_header_Z_N_info(polar_corr_calc_t *s_polar,
                                     float *r, float *cormat3d,
                                     const cuFloatComplex *A,
                                     const cuFloatComplex *B,
                                     const float *sqsums_A,
                                     const float *sqsums_B,
                                     pfts_Sizes *p_pfts_Sizes,
                                     int npart, int nrot, int nk,
                                     float alpha, int depth);
  int print_function_header_C_N_info(deviceDetails_t * s_devD,
                                     polar_corr_calc_t *s_carte,
                                     float *r, cuFloatComplex *shmat,
                                     const cuFloatComplex *A,
                                     const cuFloatComplex *B,
                                     img_2D_cart_Sizes *p_img_2D_cart_Sizes,
                                     int vx, int vy, int vz,
                                     float alpha, int depth);
  int print_sumVecs_X_N(int npart, float *sum_vec, float *h_ptr_vec,
                        float *r, double *sum_dble_vec, int depth);
  int print_sizeA_B_Z_N_info(const cuFloatComplex *A,
                             const cuFloatComplex *B,
                             pfts_Sizes *p_pfts_Sizes);
  int print_sizeA_B_C_N_info(const cuFloatComplex *A,
                             const cuFloatComplex *B,
                             img_2D_cart_Sizes *p_img_2D_cart_Sizes);
  int print_A_and_B_depth(const cuFloatComplex *A,
                          const cuFloatComplex *B,
                          pfts_Sizes *p_pfts_Sizes,
                          int npart, int nrot, int nk,
                          int depth);
  int print_A_and_B_V_depth(cuFloatComplex *shmat,
                            const cuFloatComplex *A,
                            const cuFloatComplex *B,
                            img_2D_cart_Sizes *p_img_2D_cart_Sizes,
                            int vx, int vy, int vz,
                            int depth);
  int print_cmat2sh_V_depth(const cuFloatComplex *d_cmat2sh,
                            img_2D_cart_Sizes *p_img_2D_cart_Sizes,
                            int vx, int vy, int vz,
                            int depth);
  int print_C_HadrProd(const float *C,
                       pfts_Sizes *p_pfts_Sizes);
  int print_C_cormat3d(const float *cormat3d,
                       pfts_Sizes *p_pfts_Sizes, const char *filename);
  int print_dC_HadrProd(float *d_C,
                        pfts_Sizes *p_pfts_Sizes);
  /*calculators*/
  int polarft_gencorrAll_Z_N_cpu(polar_corr_calc_t *s_polar,
                                 float *r, float *cormat3d,
                                 const cuFloatComplex *A,
                                 const cuFloatComplex *B,
                                 const float *sqsums_A,
                                 const float *sqsums_B,
                                 pfts_Sizes *p_pfts_Sizes,
                                 int npart, int nrot, int nk,
                                 float alpha);

#ifdef __cplusplus
}
#endif

#endif /* _POLARFT_GPU_H_ */
