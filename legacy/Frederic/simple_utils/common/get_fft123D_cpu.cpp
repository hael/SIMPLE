/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 17th June 2015
 *
 *      June 2015
 *
 *   code to obtain fftw on cpu
 *
 * @precisions normal z -> s d c
 */

#include "get_fft123D_cpu.h"
/* CPU header files */
#include <fftw3.h>

#ifdef __cplusplus
extern "C" {
#endif
/***************************************************************************/
/**
 *  FORTRAN API - math functions (simple interface)
 **/

  /* 3D fourier transform */
  /* Z2Z */
  int gather_fft_3D_Z2Z_cpu_cpp(int *nx, int *ny, int *nz,
				double complex *in, double complex *out,
				int *sign) {
    int rc = 0;
    int nptx = *nx;
    int npty = *ny;
    int nptz = *nz;

    int npts = nptx * npty * nptz;

    int direction = *sign;
    /* the plan for the fftw */
    fftw_plan plan_fwd;
    fftw_plan plan_bwd;

    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_BLUE  "nptz=%i, "ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty*nptz=%i "
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_out=%i, nptx*npty*nptz=%i\n",
	   nptx, npty, nptz, npts, nptx*npty*nptz, npts, nptx*npty*nptz);
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*npty*nptz*sizeof(double complex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*npty*nptz*sizeof(double complex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*npty*nptz*sizeof(double complex), nptx*npty*nptz*sizeof(double complex)/1.e6,
	   nptx*npty*nptz*sizeof(double complex), nptx*npty*nptz*sizeof(double complex)/1.e6);

    if ( direction == FFTW_FORWARD ) {
      plan_fwd = fftw_plan_dft_3d(nptx, npty, nptz, in, out, direction, FFTW_ESTIMATE);
      fftw_execute(plan_fwd);
      fftw_destroy_plan ( plan_fwd );
    } else if ( direction == FFTW_BACKWARD) {
      plan_bwd = fftw_plan_dft_3d(nptx, npty, nptz, in, out, direction, FFTW_ESTIMATE);
      fftw_execute(plan_bwd);
      fftw_destroy_plan ( plan_bwd );
    } else { rc = get_warning_message_fft123D_cpu(); }

    return rc;
  }
  /* C2C */
  int gather_fft_3D_C2C_cpu_cpp(int *nx, int *ny, int *nz,
				fftwf_complex *in, fftwf_complex *out,
				int *sign) {
    int rc = 0;
    int nptx = *nx;
    int npty = *ny;
    int nptz = *nz;

    int npts = nptx * npty * nptz;

    int direction = *sign;
    /* the plan for the fftw */
    fftwf_plan plan_fwd;
    fftwf_plan plan_bwd;

    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_BLUE  "nptz=%i, "ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty*nptz=%i "
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_out=%i, nptx*npty*nptz=%i\n",
	   nptx, npty, nptz, npts, nptx*npty*nptz, npts, nptx*npty*nptz);
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*npty*nptz*sizeof(fftwf_complex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*npty*nptz*sizeof(fftwf_complex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*npty*nptz*sizeof(fftwf_complex), nptx*npty*nptz*sizeof(fftwf_complex)/1.e6,
	   nptx*npty*nptz*sizeof(fftwf_complex), nptx*npty*nptz*sizeof(fftwf_complex)/1.e6);

    if ( direction == FFTW_FORWARD ) {
      plan_fwd = fftwf_plan_dft_3d(nptx, npty, nptz, in, out, direction, FFTW_ESTIMATE);
      fftwf_execute(plan_fwd);
      fftwf_destroy_plan ( plan_fwd );
    } else if ( direction == FFTW_BACKWARD) {
      plan_bwd = fftwf_plan_dft_3d(nptx, npty, nptz, in, out, direction, FFTW_ESTIMATE);
      fftwf_execute(plan_bwd);
      fftwf_destroy_plan ( plan_bwd );
    } else { rc = get_warning_message_fft123D_cpu(); }

    return rc;
  }
  /* D2Z */
  int gather_fft_3D_D2Z_cpu_cpp(int *nx, int *ny, int *nz,
				double *in, double complex *out) {
    int rc = 0;
    int nptx = *nx;
    int npty = *ny;
    int nptz = *nz;

    int npts = nptx * npty * nptz;
    int npts_out = nptx * npty * (nptz / 2 + 1);

    /* the plan for the fftw */
    fftw_plan plan_fwd;
    
    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_BLUE  "nptz=%i, "ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty*nptz=%i "
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_out=%i, nptx*npty*(nptz/2+1)=%i\n",
	   nptx, npty, nptz, npts, nptx*npty*nptz, npts_out, nptx*npty*(nptz/2+1));
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*npty*nptz*sizeof(double complex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*npty*(nptz/2+1)*sizeof(double complex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*npty*nptz*sizeof(double complex), nptx*npty*nptz*sizeof(double complex)/1.e6,
	   nptx*npty*(nptz/2+1)*sizeof(double complex), nptx*npty*(nptz/2+1)*sizeof(double complex)/1.e6);

    plan_fwd = fftw_plan_dft_r2c_3d(nptx, npty, nptz, in, out, FFTW_ESTIMATE);
    fftw_execute(plan_fwd);
    fftw_destroy_plan ( plan_fwd );

    return rc;
  }
  /* S2C */
  int gather_fft_3D_S2C_cpu_cpp(int *nx, int *ny, int *nz,
				float *in, fftwf_complex *out) {
    int rc = 0;
    int nptx = *nx;
    int npty = *ny;
    int nptz = *nz;

    int npts = nptx * npty * nptz;
    int npts_out = nptx * npty * (nptz / 2 + 1);

    /* the plan for the fftwf single precision */
    fftwf_plan plan_fwd;

    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_BLUE  "nptz=%i, "ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty*nptz=%i "
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_out=%i, nptx*npty*(nptz/2+1)=%i\n",
	   nptx, npty, nptz, npts, nptx*npty*nptz, npts_out, nptx*npty*(nptz/2+1));
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*npty*nptz*sizeof(fftwf_complex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*npty*(nptz/2+1)*sizeof(fftwf_complex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*npty*nptz*sizeof(fftwf_complex), nptx*npty*nptz*sizeof(fftwf_complex)/1.e6,
	   nptx*npty*(nptz/2+1)*sizeof(fftwf_complex), nptx*npty*(nptz/2+1)*sizeof(fftwf_complex)/1.e6);

    plan_fwd = fftwf_plan_dft_r2c_3d(nptx, npty, nptz, in, out, FFTW_ESTIMATE);
    fftwf_execute(plan_fwd);
    fftwf_destroy_plan ( plan_fwd );

    return rc;
  }
  /* Z2D */
  int gather_fft_3D_Z2D_cpu_cpp(int *nx, int *ny, int *nz,
				double complex *in, double *out) {
    int rc = 0;
    int nptx = *nx;
    int npty = *ny;
    int nptz = *nz;

    int npts = nptx * npty * nptz;
    int npts_in = nptx * npty * (nptz / 2 + 1);

    /* the plan for the fftw */
    fftw_plan plan_fwd;

    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_BLUE  "nptz=%i, "ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty*nptz=%i "
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_in=%i, nptx*npty*(nptz/2+1)=%i\n",
	   nptx, npty, nptz, npts, nptx*npty*nptz, npts_in, nptx*npty*(nptz/2+1));
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*npty*(nptz/2+1)*sizeof(double complex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*npty*nptz*sizeof(double complex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*npty*(nptz/2+1)*sizeof(double complex), nptx*npty*(nptz/2+1)*sizeof(double complex)/1.e6,
	   nptx*npty*nptz*sizeof(double complex), nptx*npty*nptz*sizeof(double complex)/1.e6);

    plan_fwd = fftw_plan_dft_c2r_3d(nptx, npty, nptz, in, out, FFTW_ESTIMATE);
    fftw_execute(plan_fwd);
    fftw_destroy_plan ( plan_fwd );

    return rc;
  }
  /* C2S */
  int gather_fft_3D_C2S_cpu_cpp(int *nx, int *ny, int *nz,
				fftwf_complex *in, float *out) {
    int rc = 0;
    int nptx = *nx;
    int npty = *ny;
    int nptz = *nz;

    int npts = nptx * npty * nptz;
    int npts_in = nptx * npty * (nptz / 2 + 1);

    /* the plan for the fftwf for single precision*/
    fftwf_plan plan_fwd;

    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_BLUE  "nptz=%i, "ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty*nptz=%i "
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_in=%i, nptx*npty*(nptz/2+1)=%i\n",
	   nptx, npty, nptz, npts, nptx*npty*nptz, npts_in, nptx*npty*(nptz/2+1));
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*npty*(nptz/2+1)*sizeof(fftwf_complex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*npty*nptz*sizeof(fftwf_complex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*npty*(nptz/2+1)*sizeof(fftwf_complex), nptx*npty*(nptz/2+1)*sizeof(fftwf_complex)/1.e6,
	   nptx*npty*nptz*sizeof(fftwf_complex), nptx*npty*nptz*sizeof(fftwf_complex)/1.e6);

    plan_fwd = fftwf_plan_dft_c2r_3d(nptx, npty, nptz, in, out, FFTW_ESTIMATE);
    fftwf_execute(plan_fwd);
    fftwf_destroy_plan ( plan_fwd );

    return rc;
  }
  /* 2D fourier transform */
  /* Z2Z */
  int gather_fft_2D_Z2Z_cpu_cpp(int *nx, int *ny,
				double complex *in, double complex *out,
				int *sign) {
    int rc = 0;
    int nptx = *nx;
    int npty = *ny;

    int npts = nptx * npty;

    int direction = *sign;
    /* the plan for the fftw */
    fftw_plan plan_fwd;
    fftw_plan plan_bwd;

    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty=%i "
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_out=%i, nptx*npty=%i\n",
	   nptx, npty, npts, nptx*npty, npts, nptx*npty);
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*npty*sizeof(double complex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*npty*sizeof(double complex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*npty*sizeof(double complex), nptx*npty*sizeof(double complex)/1.e6,
	   nptx*npty*sizeof(double complex), nptx*npty*sizeof(double complex)/1.e6);

    if ( direction == FFTW_FORWARD ) {
      plan_fwd = fftw_plan_dft_2d(nptx, npty, in, out, direction, FFTW_ESTIMATE);
      fftw_execute(plan_fwd);
      fftw_destroy_plan ( plan_fwd );
    } else if ( direction == FFTW_BACKWARD) {
      plan_bwd = fftw_plan_dft_2d(nptx, npty, in, out, direction, FFTW_ESTIMATE);
      fftw_execute(plan_bwd);
      fftw_destroy_plan ( plan_bwd );
    } else { rc = get_warning_message_fft123D_cpu(); }

    return rc;
  }
  /* C2C */
  int gather_fft_2D_C2C_cpu_cpp(int *nx, int *ny,
                                fftwf_complex *in, fftwf_complex *out,
                                int *sign) {
    int rc = 0;
    int nptx = *nx;
    int npty = *ny;

    int npts = nptx * npty;

    int direction = *sign;
    /* the plan for the fftw */
    fftwf_plan plan_fwd;
    fftwf_plan plan_bwd;

    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty=%i "
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_out=%i, nptx*npty=%i\n",
	   nptx, npty, npts, nptx*npty, npts, nptx*npty);
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*npty*sizeof(fftwf_complex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*npty*sizeof(fftwf_complex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*npty*sizeof(fftwf_complex), nptx*npty*sizeof(fftwf_complex)/1.e6,
	   nptx*npty*sizeof(fftwf_complex), nptx*npty*sizeof(fftwf_complex)/1.e6);

    if ( direction == FFTW_FORWARD ) {
      plan_fwd = fftwf_plan_dft_2d(nptx, npty, in, out, direction, FFTW_ESTIMATE);
      fftwf_execute(plan_fwd);
      fftwf_destroy_plan ( plan_fwd );
    } else if ( direction == FFTW_BACKWARD) {
      plan_bwd = fftwf_plan_dft_2d(nptx, npty, in, out, direction, FFTW_ESTIMATE);
      fftwf_execute(plan_bwd);
      fftwf_destroy_plan ( plan_bwd );
    } else { rc = get_warning_message_fft123D_cpu(); }

    return rc;
  }
  /* S2C */
  int gather_fft_2D_S2C_cpu_cpp(int *nx, int *ny,
				float *in, fftwf_complex *out) {
    int rc = 0;
    int nptx = *nx;
    int npty = *ny;

    int npts = nptx * npty;
    int npts_out = nptx * (npty / 2 + 1);

    /* the plan for the fftw */
    fftwf_plan plan_fwd;

    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty=%i "
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_out=%i, nptx*(npty/2+1)=%i\n",
	   nptx, npty, npts, nptx*npty, npts_out, nptx*(npty/2+1));
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*npty*sizeof(fftwf_complex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*(npty/2+1)*sizeof(fftwf_complex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*npty*sizeof(fftwf_complex), nptx*npty*sizeof(fftwf_complex)/1.e6,
	   nptx*(npty/2+1)*sizeof(fftwf_complex), nptx*(npty/2+1)*sizeof(fftwf_complex)/1.e6);

    plan_fwd = fftwf_plan_dft_r2c_2d(nptx, npty, in, out, FFTW_ESTIMATE);
    fftwf_execute(plan_fwd);
    fftwf_destroy_plan ( plan_fwd );

    return rc;
  }
  /* D2Z */
  int gather_fft_2D_D2Z_cpu_cpp(int *nx, int *ny,
				double *in, double complex *out) {
    int rc = 0;
    int nptx = *nx;
    int npty = *ny;

    int npts = nptx * npty;
    int npts_out = nptx * (npty / 2 + 1);

    /* the plan for the fftw */
    fftw_plan plan_fwd;

    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty=%i "
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_out=%i, nptx*(npty/2+1)=%i\n",
	   nptx, npty, npts, nptx*npty, npts_out, nptx*(npty/2+1));
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*npty*sizeof(double complex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*(npty/2+1)*sizeof(double complex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*npty*sizeof(double complex), nptx*npty*sizeof(double complex)/1.e6,
	   nptx*(npty/2+1)*sizeof(double complex), nptx*(npty/2+1)*sizeof(double complex)/1.e6);

    plan_fwd = fftw_plan_dft_r2c_2d(nptx, npty, in, out, FFTW_ESTIMATE);
    fftw_execute(plan_fwd);
    fftw_destroy_plan ( plan_fwd );

    return rc;
  }
  /* Z2D */
  int gather_fft_2D_Z2D_cpu_cpp(int *nx, int *ny,
				double complex *in, double *out) {
    int rc = 0;
    int nptx = *nx;
    int npty = *ny;

    int npts = nptx * npty;
    int npts_in = nptx * (npty / 2 + 1);

    /* the plan for the fftw */
    fftw_plan plan_fwd;

    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty=%i "
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_in=%i, nptx*npty*(nptz/2+1)=%i\n",
	   nptx, npty, npts, nptx*npty, npts_in, nptx*(npty/2+1));
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*(npty/2+1)*sizeof(double complex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*npty*sizeof(double complex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*(npty/2+1)*sizeof(double complex), nptx*(npty/2+1)*sizeof(double complex)/1.e6,
	   nptx*npty*sizeof(double complex), nptx*npty*sizeof(double complex)/1.e6);

    plan_fwd = fftw_plan_dft_c2r_2d(nptx, npty, in, out, FFTW_ESTIMATE);
    fftw_execute(plan_fwd);
    fftw_destroy_plan ( plan_fwd );

    return rc;
  }
  /* C2S */
  int gather_fft_2D_C2S_cpu_cpp(int *nx, int *ny,
				fftwf_complex *in, float *out) {
    int rc = 0;
    int nptx = *nx;
    int npty = *ny;

    int npts = nptx * npty;
    int npts_in = nptx * (npty / 2 + 1);

    /* the plan for the fftw */
    fftwf_plan plan_fwd;

    printf(ANSI_COLOR_BRIGHT_YELLOW"nptx=%i, "ANSI_COLOR_BRIGHT_GREEN"npty=%i, "
	   ANSI_COLOR_BRIGHT_CYAN  "npts=%i, nptx*npty=%i "
	   ANSI_COLOR_BOLD_BRIGHT_CYAN  "npts_in=%i, nptx*npty*(nptz/2+1)=%i\n",
	   nptx, npty, npts, nptx*npty, npts_in, nptx*(npty/2+1));
    printf(ANSI_COLOR_BRIGHT_RED"size of in=nptx*(npty/2+1)*sizeof(fftwf_complex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_BOLD_BRIGHT_RED"size of out=nptx*npty*sizeof(fftwf_complex)=%lu, Memory=%.2f MB\n"
	   ANSI_COLOR_RESET,
	   nptx*(npty/2+1)*sizeof(fftwf_complex), nptx*(npty/2+1)*sizeof(fftwf_complex)/1.e6,
	   nptx*npty*sizeof(fftwf_complex), nptx*npty*sizeof(fftwf_complex)/1.e6);

    plan_fwd = fftwf_plan_dft_c2r_2d(nptx, npty, in, out, FFTW_ESTIMATE);
    fftwf_execute(plan_fwd);
    fftwf_destroy_plan ( plan_fwd );

    return rc;
  }
  /* 1D fourier transform */
  /* Z2Z */
  int gather_fft_1D_Z2Z_cpu_cpp(int *nx, 
				double complex *in, double complex *out,
				int *sign) {
    int rc = 0;

    int nptx = *nx;
    int direction = *sign;

    fftw_plan plan_fwd;
    fftw_plan plan_bwd;

    if ( direction == FFTW_FORWARD ) {
      plan_fwd = fftw_plan_dft_1d(nptx, in, out, direction, FFTW_ESTIMATE);
      fftw_execute(plan_fwd);
      fftw_destroy_plan ( plan_fwd );
    } else if ( direction == FFTW_BACKWARD) {
      plan_bwd = fftw_plan_dft_1d(nptx, in, out, direction, FFTW_ESTIMATE);
      fftw_execute(plan_bwd);
      fftw_destroy_plan ( plan_bwd );
    } else { rc = get_warning_message_fft123D_cpu(); }

    return rc;
  }

  /* test code to test the FFTW3 libtrary on CPU must have 
   * #include <fftw3.h> header and not the CUDA ones
   */
  int test_fft_1D_cpu_cpp(double *xi, double *xf, int *nstep,
			  double *x, double *fx, double complex *fh) {
    int rc = 0;

    int istep;
    int nstepi = *nstep;

    fftw_complex *in_c;
    fftw_complex *in_back_c;

    fftw_plan plan_fwd;
    fftw_plan plan_bwd;

    in_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nstepi);
    in_back_c = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*nstepi);

    for (istep = 0 ; istep < nstepi ; istep++ ) { 
      in_c[istep] = fx[istep]; 
    }

    plan_fwd = fftw_plan_dft_1d(nstepi, in_c, fh, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan_fwd);

    plan_bwd = fftw_plan_dft_1d(nstepi, fh, in_back_c, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan_bwd);

    FILE * testfile;
    testfile = fopen("Testfile.log","a");

    for (istep = 0 ; istep < nstepi ; istep++ ) {
      //fprintf(testfile,"%15.8f %15.8f\n",x[istep], fx[istep]);
      //fprintf(testfile,"%15.5f %15.5f\n", x[istep], creal(in_c[istep]));
      fprintf(testfile,"%15.5f %15.5f\n", x[istep], creal(in_back_c[istep])/nstepi);
    }

    fclose(testfile);

    fftw_destroy_plan ( plan_fwd );
    fftw_destroy_plan ( plan_bwd );
    fftw_free(in_c);
    fftw_free(in_back_c);

    return rc;
  }

  /* getting the warning message from the cuda */
  int get_warning_message_fft123D_cpu() {
    int rc = 0;

    printf("***************************WARNING*****************************\n");
    printf("You need to set the direction of the transform correctly       \n");
    printf("***************************************************************\n");
    printf("\n");
    printf("Exit at Line %i in file %s %s\n",__LINE__,__FILE__,__FUNCTION__);
    printf("\n");
    rc = -1;

    return rc;
  }

  /* the aliases for external access */

  /* 1D */
  //double precision
  extern "C" int gather_fft_1d_z2z_cpu_cpp_() __attribute__((weak,alias("gather_fft_1D_Z2Z_cpu_cpp")));
  /* 2D */
  extern "C" int gather_fft_2d_c2c_cpu_cpp_() __attribute__((weak,alias("gather_fft_2D_C2C_cpu_cpp")));
  extern "C" int gather_fft_2d_s2c_cpu_cpp_() __attribute__((weak,alias("gather_fft_2D_S2C_cpu_cpp")));
  extern "C" int gather_fft_2d_c2s_cpu_cpp_() __attribute__((weak,alias("gather_fft_2D_C2S_cpu_cpp")));
  //double precision
  extern "C" int gather_fft_2d_z2z_cpu_cpp_() __attribute__((weak,alias("gather_fft_2D_Z2Z_cpu_cpp")));
  extern "C" int gather_fft_2d_d2z_cpu_cpp_() __attribute__((weak,alias("gather_fft_2D_D2Z_cpu_cpp")));
  extern "C" int gather_fft_2d_z2d_cpu_cpp_() __attribute__((weak,alias("gather_fft_2D_Z2D_cpu_cpp")));
  /* 3D */
  extern "C" int gather_fft_3d_c2c_cpu_cpp_() __attribute__((weak,alias("gather_fft_3D_C2C_cpu_cpp")));
  extern "C" int gather_fft_3d_s2c_cpu_cpp_() __attribute__((weak,alias("gather_fft_3D_S2C_cpu_cpp")));
  extern "C" int gather_fft_3d_c2s_cpu_cpp_() __attribute__((weak,alias("gather_fft_3D_C2S_cpu_cpp")));
  //double precision
  extern "C" int gather_fft_3d_z2z_cpu_cpp_() __attribute__((weak,alias("gather_fft_3D_Z2Z_cpu_cpp")));
  extern "C" int gather_fft_3d_d2z_cpu_cpp_() __attribute__((weak,alias("gather_fft_3D_D2Z_cpu_cpp")));
  extern "C" int gather_fft_3d_z2d_cpu_cpp_() __attribute__((weak,alias("gather_fft_3D_Z2D_cpu_cpp")));

#ifdef __cplusplus
}
#endif
