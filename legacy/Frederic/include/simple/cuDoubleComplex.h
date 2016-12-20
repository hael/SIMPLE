/*
 */

#if !defined(CU_DOUBLECOMPLEX_H_)
#define CU_DOUBLECOMPLEX_H_

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

#include <math.h>       /* import fabsf, sqrt */

#if defined (CUDA) /*preprossing for the CUDA environment */

#include "vector_types.h"

  typedef float2 cuFloatComplex;

  __host__ __device__ static __inline__ cuFloatComplex cuFCaddf (float x, cuFloatComplex y)
  {
    return make_cuFloatComplex ( x + cuCrealf(y), cuCimagf(y) );
  }

  __host__ __device__ static __inline__ cuFloatComplex cuCFaddf (cuFloatComplex x, float y)
  {
    return make_cuFloatComplex (cuCrealf(x) + y, cuCimagf(x) );
  }

  __host__ __device__ static __inline__ cuFloatComplex cuFCsubf (float x, cuFloatComplex y)
  {
    return make_cuFloatComplex ( x - cuCrealf(y), - cuCimagf(y));
  }

  __host__ __device__ static __inline__ cuFloatComplex cuCFsubf (cuFloatComplex x, float y)
  {
    return make_cuFloatComplex ( cuCrealf(x) - y, cuCimagf(x) );
  }

  /* This implementation could suffer from intermediate overflow even though
   * the final result would be in range. However, various implementations do
   * not guard against this (presumably to avoid losing performance), so we 
   * don't do it either to stay competitive.
   */

  __host__ __device__ static __inline__ cuFloatComplex cuFCmulf (float x, cuFloatComplex y)
  {
    cuFloatComplex prod;
    prod = make_cuFloatComplex  (x * cuCrealf(y) ,
				 x * cuCimagf(y) );
    return prod;
  }


  __host__ __device__ static __inline__ cuFloatComplex cuCFmulf (cuFloatComplex x, float y)
  {
    cuFloatComplex prod;
    prod = make_cuFloatComplex  (cuCrealf(x) * y ,
                                 cuCimagf(x) * y );
    return prod;
  }

  /* This implementation guards against intermediate underflow and overflow
   * by scaling. Such guarded implementations are usually the default for
   * complex library implementations, with some also offering an unguarded,
   * faster version.
   */

  __host__ __device__ static __inline__ cuFloatComplex cuFCdivf (float x, cuFloatComplex y)
  {
    cuFloatComplex quot;
    float s = fabsf(cuCrealf(y)) + fabsf(cuCimagf(y));
    float oos = 1.0f / s;
    float ars = x * oos;
    float brs = cuCrealf(y) * oos;
    float bis = cuCimagf(y) * oos;
    s = (brs * brs) + (bis * bis);
    oos = 1.0f / s;
    quot = make_cuFloatComplex ( ( (ars * brs) ) * oos ,
				 (-(ars * bis) ) * oos );
    return quot;
  }
  
  __host__ __device__ static __inline__ cuFloatComplex cuCFdivf (cuFloatComplex x, float y)
  {
    cuFloatComplex quot;
    float s = fabsf( y );
    float oos = 1.0f / s;
    float ars = cuCrealf(x) * oos;
    float ais = cuCimagf(x) * oos;
    float brs = y * oos;

    quot = make_cuFloatComplex ( (ars * brs) , 
				 (ais * brs) );
    return quot;
  }

  /* Double precision */
  typedef double2 cuDoubleComplex;

  __host__ __device__ static __inline__ cuDoubleComplex cuDCadd(double x, cuDoubleComplex y)
  {
    return make_cuDoubleComplex ( x + cuCreal(y), cuCimag(y) );
  }

  __host__ __device__ static __inline__ cuDoubleComplex cuCDadd(cuDoubleComplex x, double y)
  {
    return make_cuDoubleComplex (cuCreal(x) + y, cuCimag(x) );
  }

  __host__ __device__ static __inline__ cuDoubleComplex cuDCsub(double x, cuDoubleComplex y)
  {
    return make_cuDoubleComplex ( x - cuCreal(y), - cuCimag(y) );
  }

  __host__ __device__ static __inline__ cuDoubleComplex cuCDsub(cuDoubleComplex x, double y)
  {
    return make_cuDoubleComplex (cuCreal(x) - y, cuCimag(x) );
  }

  /* This implementation could suffer from intermediate overflow even though
   * the final result would be in range. However, various implementations do
   * not guard against this (presumably to avoid losing performance), so we 
   * don't do it either to stay competitive.
   */

  __host__ __device__ static __inline__ cuDoubleComplex cuDCmul(double x, cuDoubleComplex y)
  {
    cuDoubleComplex prod;
    prod = make_cuDoubleComplex ( x * cuCreal(y) ,
				  x * cuCimag(y) );
    return prod;
  }

  __host__ __device__ static __inline__ cuDoubleComplex cuCDmul(cuDoubleComplex x, double y)
  {
    cuDoubleComplex prod;
    prod = make_cuDoubleComplex ( cuCreal(x) * y ,
				  cuCimag(x) * y );
    return prod;
  }

  /* This implementation guards against intermediate underflow and overflow
   * by scaling. Such guarded implementations are usually the default for
   * complex library implementations, with some also offering an unguarded,
   * faster version.
   */

  __host__ __device__ static __inline__ cuDoubleComplex cuDCdiv(double x, cuDoubleComplex y)
  {
    cuDoubleComplex quot;
    double s = (fabs(cuCreal(y))) + (fabs(cuCimag(y)));
    double oos = 1.0 / s;
    double ars = x * oos;
    double brs = cuCreal(y) * oos;
    double bis = cuCimag(y) * oos;
    s = (brs * brs) + (bis * bis);
    oos = 1.0 / s;
    quot = make_cuDoubleComplex ( ( (ars * brs) ) * oos ,
				  (-(ars * bis) ) * oos );
    return quot;
  }

  __host__ __device__ static __inline__ cuDoubleComplex cuCDdiv(cuDoubleComplex x, double y)
  {
    cuDoubleComplex quot;
    double s = fabs( y );
    double oos = 1.0 / s;
    double ars = cuCreal(x) * oos;
    double ais = cuCimag(x) * oos;
    double brs = y * oos;
    //double bis = cuCimag(y) * oos;
    //s = (brs * brs) + (bis * bis);
    //oos = 1.0 / s;

    quot = make_cuDoubleComplex ( (ars * brs) ,
				  (ais * brs) );
    return quot;
  }


#if defined(__cplusplus)
}
#endif /* __cplusplus */

/* aliases */
typedef cuFloatComplex cuComplex;

/* float-to-double promotion */

__host__ __device__ static __inline__  cuComplex cuDCfmaf( float x, cuComplex y, cuComplex d)
{
  float real_res;
  float imag_res;
    
  real_res = (x *  cuCrealf(y)) + cuCrealf(d);
  imag_res = (x *  cuCimagf(y)) + cuCimagf(d);
     
  return make_cuComplex(real_res, imag_res);
}

__host__ __device__ static __inline__  cuComplex cuCDfmaf( cuComplex x, float y, cuComplex d)
{
  float real_res;
  float imag_res;
    
  real_res = (cuCrealf(x) * y) + cuCrealf(d);
  imag_res = (cuCimagf(x) * y) + cuCimagf(d);          
     
  return make_cuComplex(real_res, imag_res);
}

__host__ __device__ static __inline__  cuDoubleComplex cuDCfma( double x, cuDoubleComplex y, cuDoubleComplex d)
{
  double real_res;
  double imag_res;

  real_res = (x *  cuCreal(y)) + cuCreal(d);
  imag_res = (x *  cuCimag(y)) + cuCimag(d);
     
  return make_cuDoubleComplex(real_res, imag_res);
}

__host__ __device__ static __inline__  cuDoubleComplex cuCDfma( cuDoubleComplex x, double y, cuDoubleComplex d)
{
  double real_res;
  double imag_res;

  real_res = (cuCreal(x) * y) + cuCreal(d);
  imag_res = (cuCimag(x) * y) + cuCimag(d);

  return make_cuDoubleComplex(real_res, imag_res);
}

/* product of double time complex scaled by a constant alpha 
 *
 * alpha: complex
 *     x: double
 *     y: cuDoubleComplex
 * returns: alpha * x * y
 */
__host__ __device__ static __inline__  cuDoubleComplex cuDCmula( double x, cuDoubleComplex y, cuDoubleComplex alpha)
{
  double real_res;
  double imag_res;

  real_res = (x * cuCreal(y));
  imag_res = (x * cuCimag(y));
 
  real_res = (cuCreal(alpha) * real_res) - (cuCimag(alpha) * imag_res);
  imag_res = (cuCreal(alpha) * imag_res) + (cuCimag(alpha) * real_res);
    
  return make_cuDoubleComplex(real_res, imag_res);
}
/* product of double time complex scaled by a constant alpha 
 *
 * alpha: complex
 *     x: cuDoubleComplex
 *     y: double
 * returns: alpha * x * y
 */
__host__ __device__ static __inline__  cuDoubleComplex cuCDmula( cuDoubleComplex x, double y, cuDoubleComplex alpha)
{
  double real_res;
  double imag_res;

  real_res = (cuCreal(x) * y);
  imag_res = (cuCimag(x) * y);
 
  real_res = (cuCreal(alpha) * real_res) - (cuCimag(alpha) * imag_res);
  imag_res = (cuCreal(alpha) * imag_res) + (cuCimag(alpha) * real_res);
    
  return make_cuDoubleComplex(real_res, imag_res);
}

#endif /* CUDA */

#endif /* !defined(CU_DOUBLECOMPLEX_H_) */
