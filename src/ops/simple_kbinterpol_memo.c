#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static double alpha = 0., beta = 0., betasq = 0., oneoW = 0., piW = 0., twooW = 0., W = 0., Whalf = 0.;
static double dx = 0.;
static int Nx = -1;

static double* apod_mem         = NULL;
static double* apod_deriv_mem   = NULL;
static double* dapod_mem        = NULL;
static double* dapod_deriv_mem  = NULL;

#ifndef M_PI
const double M_PI = 3.14159265358979323846;
#endif

void kbinterp_memo_set(double Whalf_in, double alpha_in, int Nx_in);
void kbinterp_memo_memoize();
void kbinterp_memo_kill();
double  apod_nointerp(double x);
double  apod_lininterp(double x);
double dapod_nointerp(double x);
double dapod_lininterp(double x);
float  apod_nointerp_sp(float x);
float  apod_lininterp_sp(float x);
float  apod_nointerp_sp(float x);
float  apod_lininterp_sp(float x);
static double apod( double x );
static double dapod( double x ) ;
static double bessi0( double x );
static double bessi1( double x );
static double deallocate();

void
kbinterp_memo_set( double Whalf_in, double alpha_in, int Nx_in )
{
    Whalf  = Whalf_in;
    alpha  = alpha_in;
    W      = 2.0 * Whalf;
    piW    = M_PI * W;
    if( Whalf <= 1.5 )
        beta = 7.4;
    else
        beta = M_PI * sqrt((W*W / alpha / alpha) *
                           (alpha - 0.5)*(alpha - 0.5) - 0.8);
    betasq = beta * beta;
    twooW  = 2.0 / W;
    oneoW  = 1.0 / W;
    if (Nx_in <= 1)
    {
	printf("error in simple_kbinterpol_memo :: new: Nx <= 1");
	exit(1);
    }
    Nx = Nx_in;
    dx = Whalf / (Nx - 1);
}

void
kbinterp_memo_memoize()
{
    int i;
    double x;
    deallocate();
    apod_mem        = (double*)malloc(sizeof(double)*Nx);
    apod_deriv_mem  = (double*)malloc(sizeof(double)*Nx);
    dapod_mem       = (double*)malloc(sizeof(double)*Nx);
    dapod_deriv_mem = (double*)malloc(sizeof(double)*Nx);
    for (i = 0; i < Nx; ++i)
    {
	x = dx * i;
	apod_mem[i]  = apod(x);
	dapod_mem[i] = dapod(x);
    }
     apod_deriv_mem[0] = ( apod_mem[1]    -  apod_mem[0]   ) / dx; // forward deriv
    dapod_deriv_mem[0] = (dapod_mem[1]    - dapod_mem[0]   ) / dx; // forward deriv
     apod_deriv_mem[1] = 2.0 * ( apod_mem[3]    -  apod_mem[0]   ) / (2. * dx); // 3p stencil
    dapod_deriv_mem[1] = 2.0 * (dapod_mem[3]    - dapod_mem[0]   ) / (2. * dx); // 3p stencil

    for (i = 2; i < Nx-3; ++i)
    {
	/* 5p stencil */
	 apod_deriv_mem[i] = (- apod_mem[i+2]+ 8. *  apod_mem[i+1] - 8. *  apod_mem[i-1] +  apod_mem[i-2]) / (12. * dx);
	dapod_deriv_mem[i] = (-dapod_mem[i+2]+ 8. * dapod_mem[i+1] - 8. * dapod_mem[i-1] + dapod_mem[i-2]) / (12. * dx);
    }

     apod_deriv_mem[Nx-3] = ( apod_mem[Nx-3] -  apod_mem[Nx-4]) / dx; // forward deriv
    dapod_deriv_mem[Nx-3] = (dapod_mem[Nx-3] - dapod_mem[Nx-4]) / dx; // forward deriv
     apod_deriv_mem[Nx-2] = ( apod_mem[Nx-2] -  apod_mem[Nx-3]) / dx; // forward deriv
    dapod_deriv_mem[Nx-2] = (dapod_mem[Nx-2] - dapod_mem[Nx-3]) / dx; // forward deriv
     apod_deriv_mem[Nx-1] = 0.;
    dapod_deriv_mem[Nx-1] = 0.;
}

void
kbinterp_memo_kill()
{
    deallocate();
    alpha = 0.;
    beta = 0.;
    betasq = 0.;
    oneoW = 0.;
    piW = 0.;
    twooW = 0.;
    W = 0.;
    Whalf = 0.;
    dx = 0.;
    Nx = -1;
}

double
apod_nointerp(double x)
{
    int i;
    i = (int)(round(x / dx));
    if (( i < 0 ) || ( i > Nx-1 ))
	return 0.;
    return apod_mem[i];
}

double
apod_lininterp(double x)
{
    int i;
    double x0;
    double delta;
    i = (int)(round(x / dx));
    if (( i < 0 ) || ( i > Nx-1 ))
	return 0.;
    x0 = i * dx;
    delta = x - x0;
    return apod_mem[i] + delta * apod_deriv_mem[i];
}

double
dapod_nointerp(double x)
{
    int i;
    i = (int)(round(x / dx));
    if (( i < 0 ) || ( i > Nx-1 ))
	return 0.;
    return dapod_mem[i];
}

double
dapod_lininterp(double x)
{
    int i;
    double x0;
    double delta;
    i = (int)(round(x / dx));
    if (( i < 0 ) || ( i > Nx-1 ))
	return 0.;
    x0 = i * dx;
    delta = x - x0;
    return dapod_mem[i] + delta * dapod_deriv_mem[i];
}

float
apod_nointerp_sp(float x)
{
    int i;
    i = (int)(round(x / dx));
    if (( i < 0 ) || ( i > Nx-1 ))
	return 0.;
    return (float)(apod_mem[i]);
}

float
apod_lininterp_sp(float x)
{
    int i;
    double x0;
    double delta;
    i = (int)(round(x / dx));
    if (( i < 0 ) || ( i > Nx-1 ))
	return 0.;
    x0 = i * dx;
    delta = x - x0;
    return (float)(apod_mem[i] + delta * apod_deriv_mem[i]);
}

float
dapod_nointerp_sp(float x)
{
    int i;
    i = (int)(round(x / dx));
    if (( i < 0 ) || ( i > Nx-1 ))
	return 0.;
    return (float)(dapod_mem[i]);
}

float
dapod_lininterp_sp(float x)
{
    int i;
    double x0;
    double delta;
    i = (int)(round(x / dx));
    if (( i < 0 ) || ( i > Nx-1 ))
	return 0.;
    x0 = i * dx;
    delta = x - x0;
    return (float)(dapod_mem[i] + delta * dapod_deriv_mem[i]);
}

static double
apod( double x )
{
    double r, arg;
    arg = twooW * x;
    arg = 1. - arg * arg;
    if (arg < 0.) return 0.;
    return oneoW * bessi0(beta * sqrt(arg));
}

static double
dapod( double x )
{
    double r, arg, sqrtarg;
    arg  = twooW * x;
    arg  = 1. - arg * arg;
    if (arg <= 0.) return 0.;
    sqrtarg = sqrt(arg);
    return - 4. * beta * x * bessi1(beta * sqrtarg) /
        sqrtarg / (W * W * W);
}

static double
bessi0( double x )
{
    double y, ax;
    ax = x;
    if ( ax < 3.75 )
    {
        y = x / 3.75;
        y = y * y;
        return 1.0 +
                y * (3.5156229 + y * (3.0899424 + y * (1.2067492 +
                y * (0.2659732 + y * (0.0360768 + y *  0.0045813)))));
    }
    else
    {
        y = 3.75 / ax;
        return ( 0.39894228 + y * ( 0.01328592 +
                y * ( 0.00225319 + y * ( -0.00157565 + y * ( 0.00916281 +
                y * (-0.02057706 + y * (  0.02635537 + y * (-0.01647633 +
                y *  0.00392377)))))))) * exp( ax ) / sqrt( ax );
    }
}

static double
bessi1( double x )
{
    double y, ax, bx;
    ax = x;
    if ( ax < 3.75)
    {
        y = x / 3.75;
        y = y * y;
        return x * (
            0.5 + y * (0.87890594 + y * (0.51498869 +
                  y * (0.15084934 + y * (0.2658733E-1 + y * (0.301532E-2 +
                  y * 0.32411E-3))))));
    }
    else
    {
        y   = 3.75 / ax;
        bx  = exp(ax) / sqrt(ax);
        ax  = 0.39894228   + y * (-0.3988024E-1 + y * (-0.362018E-2 +
                             y * ( 0.163801E-2  + y * (-0.1031555E-1 + y * (0.2282967E-1 +
                             y * (-0.2895312E-1 + y * ( 0.1787654E-1 + y * (-0.420059E-2))))))));
        return ax * bx;
    }
}

static double
deallocate()
{
    if (apod_mem != NULL)
	free(apod_mem);
    if (apod_deriv_mem != NULL)
	free(apod_deriv_mem);
    if (dapod_mem != NULL)
	free(dapod_mem);
    if (dapod_deriv_mem != NULL)
	free(dapod_deriv_mem);
}
