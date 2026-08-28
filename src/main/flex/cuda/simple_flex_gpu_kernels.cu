// CUDA-C kernels for the flex_pca insertion family (P1 of the polar/GPU plan).
//
// Mirrors insert_planes_oversamp_multi_scaled_batch (simple_flex_reconstructor_latent_ops.f90)
// exactly: same Friedel half-plane reads, same nyq-disk and window bounds tests (incl. the
// NONREDUNDANT_HMIN gate and the INSIDE per-sample bounds test), same KB apodization
// (kbinterpol apod_fast 14-term polynomial at W=1.5, alpha=2) with the same per-axis sum
// normalization, same dscale/rscale accumulation into per-component cmat/rho. The only
// permitted difference is fp32 atomic summation ORDER, which the statistical gates own.
//
// Lifecycle (host Fortran side: simple_flex_gpu.f90): insert_begin allocates and zeroes the
// device accumulators; insert_batch uploads one packed particle batch and launches; insert_fetch
// copies one component's accumulators back for the host to add into the reconstructor pair;
// insert_free releases. One active accumulation at a time (the flex stages are sequential).
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cufft.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#define WDIM   3       // LATENT_WDIM at KBWINSZ=1.5
#define OSMPL_PAD_FAC_DEV 2
#define IWINSZ 1

namespace {

float2 *d_cmat = nullptr;   // [ncomp][nvox]
float  *d_rho  = nullptr;   // [ncomp][nvox]
size_t  g_nvox = 0;
int     g_ncomp = 0;

// coupled-stats accumulation (probe M-step shape): ncomp cmplx accumulators + nrho
// density accumulators in the HOST layout (rho fastest: flat = ir + nrho*vox)
float2 *d_cmat2 = nullptr;  // [ncomp][nvox]
float  *d_rho2  = nullptr;  // [nvox][nrho], rho index fastest
size_t  g2_nvox = 0;
int     g2_ncomp = 0, g2_nrho = 0;

// kbinterpol apod_fast at KBWINSZ: Whalf = 1.5, so W = 2*Whalf = 3.0 (kbinterpol%new).
// r = (1/W) * P(u), u = 1 - (2x/W)^2; zero outside |x| > Whalf.
__device__ __forceinline__ float kb_apod_fast(float x)
{
    const float coeffs[15] = {
        1.0000000f, 1.0517297e1f, 2.7653385e1f, 3.2315430e1f, 2.1241936e1f,
        8.9363103f, 2.6107175f, 5.6036106e-1f, 9.2085685e-2f, 1.1956698e-2f,
        1.2575214e-3f, 1.0930353e-4f, 7.9831782e-6f, 4.9681336e-7f, 2.6658846e-8f };
    if (fabsf(x) > 1.5f) return 0.f;
    const float q = (2.0f/3.0f)*x;
    const float u = fmaxf(0.f, 1.f - q*q);
    float p = coeffs[14];
    #pragma unroll
    for (int i = 13; i >= 0; --i) p = p*u + coeffs[i];
    return (1.0f/3.0f)*p;
}

__global__ void insert_multi_kernel(
    const float2* __restrict__ pl_c,    // packed half-planes [nrec][k<=0 rows][Hpd] (h fastest)
    const float*  __restrict__ pl_ct,
    const float*  __restrict__ rotm,    // [nrec][nsym][9] column-major 3x3
    const float*  __restrict__ dscale,  // [nrec][ncomp]
    const float*  __restrict__ rscale,  // [nrec][ncomp]
    const int*    __restrict__ act,     // [nrec][ncomp] live component indices (1-based)
    const int*    __restrict__ nact,    // [nrec]
    const int*    __restrict__ nyqdisk, // [nrec]
    int nrec, int ncomp, int nsym,
    int hpd_lo, int hpd_hi, int kpd_lo,          // packed plane bounds (padded units, k <= 0 half)
    int hlo, int hhi, int klo, int khi, int pf,  // unpadded sample grid
    int lb1, int lb2, int lb3, int ub1, int ub2, int ub3,
    float tiny, size_t nvox,
    float2* __restrict__ cmat, float* __restrict__ rho)
{
    const int nh  = hhi - hlo + 1;
    const int nk  = khi - klo + 1;
    const size_t tid = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (tid >= (size_t)nrec*nh*nk) return;
    const int i  = (int)(tid / ((size_t)nh*nk));
    const size_t r2 = tid - (size_t)i*nh*nk;
    const int h  = hlo + (int)(r2 % nh);
    const int k  = klo + (int)(r2 / nh);
    if (h*h + k*k > nyqdisk[i]) return;
    const int na = nact[i];
    if (na == 0) return;
    // Friedel half-plane read. Planes arrive PACKED on the unpadded sample lattice (the kernels
    // only ever read pf-multiples of the padded storage, so the wrappers ship 1/pf^2 of it):
    // storage bounds hpd_lo/hpd_hi/kpd_lo are now in SAMPLE units, k <= 0 half.
    const int hpdn = hpd_hi - hpd_lo + 1;
    const size_t plsz = (size_t)(1 - kpd_lo)*hpdn;
    float2 raw; float ct;
    if (k <= 0) {
        const size_t pidx = (size_t)i*plsz + (size_t)(k - kpd_lo)*hpdn + (h - hpd_lo);
        raw = pl_c[pidx]; ct = pl_ct[pidx];
    } else {
        const size_t pidx = (size_t)i*plsz + (size_t)(-k - kpd_lo)*hpdn + (-h - hpd_lo);
        const float2 c = pl_c[pidx];
        raw = make_float2(c.x, -c.y); ct = pl_ct[pidx];
    }
    if (fabsf(raw.x) + fabsf(raw.y) <= tiny && ct <= tiny) return;
    const float pf2 = (float)(pf*pf);
    const float cbre = pf2*raw.x, cbim = pf2*raw.y;
    const int d1 = ub1-lb1+1, d2 = ub2-lb2+1;
    for (int isym = 0; isym < nsym; ++isym) {
        const float* m = rotm + ((size_t)i*nsym + isym)*9;
        // CPU convention: loc = h*(r11,r12,r13) + k*(r21,r22,r23); column-major m[(col)*3+(row)]
        const float lx = h*m[0] + k*m[1];
        const float ly = h*m[3] + k*m[4];
        const float lz = h*m[6] + k*m[7];
        const int wx0 = __float2int_rn(lx) - IWINSZ;
        const int wy0 = __float2int_rn(ly) - IWINSZ;
        const int wz0 = __float2int_rn(lz) - IWINSZ;
        if (wx0 + 2*IWINSZ < 0) continue;                          // NONREDUNDANT_HMIN gate
        if (wx0 < lb1 || wy0 < lb2 || wz0 < lb3) continue;         // bounds INSIDE, per sample
        if (wx0+2*IWINSZ > ub1 || wy0+2*IWINSZ > ub2 || wz0+2*IWINSZ > ub3) continue;
        float wx[WDIM], wy[WDIM], wz[WDIM];
        float sx = 0.f, sy = 0.f, sz = 0.f;
        #pragma unroll
        for (int j = 0; j < WDIM; ++j) {
            wx[j] = kb_apod_fast((float)(wx0+j) - lx); sx += wx[j];
            wy[j] = kb_apod_fast((float)(wy0+j) - ly); sy += wy[j];
            wz[j] = kb_apod_fast((float)(wz0+j) - lz); sz += wz[j];
        }
        const float eps = 1.1920929e-07f;   // epsilon(1.0) sp, as the CPU normalization uses
        if (fabsf(sx) > eps) { const float s = 1.f/sx; wx[0]*=s; wx[1]*=s; wx[2]*=s; }
        else                 { wx[0]=wx[1]=wx[2]=1.f/3.f; }
        if (fabsf(sy) > eps) { const float s = 1.f/sy; wy[0]*=s; wy[1]*=s; wy[2]*=s; }
        else                 { wy[0]=wy[1]=wy[2]=1.f/3.f; }
        if (fabsf(sz) > eps) { const float s = 1.f/sz; wz[0]*=s; wz[1]*=s; wz[2]*=s; }
        else                 { wz[0]=wz[1]=wz[2]=1.f/3.f; }
        for (int iq = 0; iq < na; ++iq) {
            const int q  = act[(size_t)i*ncomp + iq] - 1;
            const float ds = dscale[(size_t)i*ncomp + q];
            const float rs = rscale[(size_t)i*ncomp + q];
            float2* cq = cmat + (size_t)q*nvox;
            float*  rq = rho  + (size_t)q*nvox;
            #pragma unroll
            for (int iz = 0; iz < WDIM; ++iz) {
                const size_t zoff = (size_t)(wz0+iz-lb3)*d2*d1;
                #pragma unroll
                for (int iy = 0; iy < WDIM; ++iy) {
                    const size_t yoff = zoff + (size_t)(wy0+iy-lb2)*d1;
                    #pragma unroll
                    for (int ix = 0; ix < WDIM; ++ix) {
                        const float ww = wx[ix]*(wy[iy]*wz[iz]);
                        const size_t vidx = yoff + (size_t)(wx0+ix-lb1);
                        atomicAdd(&cq[vidx].x, ds*cbre*ww);
                        atomicAdd(&cq[vidx].y, ds*cbim*ww);
                        atomicAdd(&rq[vidx],   rs*ct*ww);
                    }
                }
            }
        }
    }
}

// Coupled-stats kernel (insert_planes_oversamp_coupled_batch_scaled): the plane value is the
// HOST-premultiplied conj(transfer)*cmplx product (Friedel of the product = its conjugate, so
// the read logic is unchanged); every component takes dscale(q)*value, every density slot takes
// rscale(ir)*ctfsq — nrho = 1 (shared), ncomp (diagonal) or npairs (packed) with the scales
// built host-side, so one kernel serves all three CPU modes.
__global__ void insert_coupled_kernel(
    const float2* __restrict__ pl_c, const float* __restrict__ pl_ct,
    const float*  __restrict__ rotm, const float* __restrict__ dscale,
    const float*  __restrict__ rscale, const int* __restrict__ nyqdisk,
    const unsigned char* __restrict__ vmask,
    int nrec, int ncomp, int nrho, int nsym,
    int hpd_lo, int hpd_hi, int kpd_lo,
    int hlo, int hhi, int klo, int khi, int pf,
    int lb1, int lb2, int lb3, int ub1, int ub2, int ub3,
    float tiny, size_t nvox,
    float2* __restrict__ cmat, float* __restrict__ rho)
{
    const int nh  = hhi - hlo + 1;
    const int nk  = khi - klo + 1;
    const size_t tid = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (tid >= (size_t)nrec*nh*nk) return;
    const int i  = (int)(tid / ((size_t)nh*nk));
    if (!vmask[i]) return;
    const size_t r2 = tid - (size_t)i*nh*nk;
    const int h  = hlo + (int)(r2 % nh);
    const int k  = klo + (int)(r2 / nh);
    if (h*h + k*k > nyqdisk[i]) return;
    // packed unpadded read (see insert_multi_kernel)
    const int hpdn = hpd_hi - hpd_lo + 1;
    const size_t plsz = (size_t)(1 - kpd_lo)*hpdn;
    float2 raw; float ct;
    if (k <= 0) {
        const size_t pidx = (size_t)i*plsz + (size_t)(k - kpd_lo)*hpdn + (h - hpd_lo);
        raw = pl_c[pidx]; ct = pl_ct[pidx];
    } else {
        const size_t pidx = (size_t)i*plsz + (size_t)(-k - kpd_lo)*hpdn + (-h - hpd_lo);
        const float2 c = pl_c[pidx];
        raw = make_float2(c.x, -c.y); ct = pl_ct[pidx];
    }
    if (fabsf(raw.x) + fabsf(raw.y) <= tiny && ct <= tiny) return;
    const float pf2 = (float)(pf*pf);
    const float cbre = pf2*raw.x, cbim = pf2*raw.y;
    const int d1 = ub1-lb1+1, d2 = ub2-lb2+1;
    for (int isym = 0; isym < nsym; ++isym) {
        const float* m = rotm + ((size_t)i*nsym + isym)*9;
        const float lx = h*m[0] + k*m[1];
        const float ly = h*m[3] + k*m[4];
        const float lz = h*m[6] + k*m[7];
        const int wx0 = __float2int_rn(lx) - IWINSZ;
        const int wy0 = __float2int_rn(ly) - IWINSZ;
        const int wz0 = __float2int_rn(lz) - IWINSZ;
        if (wx0 + 2*IWINSZ < 0) continue;
        if (wx0 < lb1 || wy0 < lb2 || wz0 < lb3) continue;
        if (wx0+2*IWINSZ > ub1 || wy0+2*IWINSZ > ub2 || wz0+2*IWINSZ > ub3) continue;
        float wx[WDIM], wy[WDIM], wz[WDIM];
        float sx = 0.f, sy = 0.f, sz = 0.f;
        #pragma unroll
        for (int j = 0; j < WDIM; ++j) {
            wx[j] = kb_apod_fast((float)(wx0+j) - lx); sx += wx[j];
            wy[j] = kb_apod_fast((float)(wy0+j) - ly); sy += wy[j];
            wz[j] = kb_apod_fast((float)(wz0+j) - lz); sz += wz[j];
        }
        const float eps = 1.1920929e-07f;
        if (fabsf(sx) > eps) { const float s = 1.f/sx; wx[0]*=s; wx[1]*=s; wx[2]*=s; }
        else                 { wx[0]=wx[1]=wx[2]=1.f/3.f; }
        if (fabsf(sy) > eps) { const float s = 1.f/sy; wy[0]*=s; wy[1]*=s; wy[2]*=s; }
        else                 { wy[0]=wy[1]=wy[2]=1.f/3.f; }
        if (fabsf(sz) > eps) { const float s = 1.f/sz; wz[0]*=s; wz[1]*=s; wz[2]*=s; }
        else                 { wz[0]=wz[1]=wz[2]=1.f/3.f; }
        #pragma unroll
        for (int iz = 0; iz < WDIM; ++iz) {
            const size_t zoff = (size_t)(wz0+iz-lb3)*d2*d1;
            #pragma unroll
            for (int iy = 0; iy < WDIM; ++iy) {
                const size_t yoff = zoff + (size_t)(wy0+iy-lb2)*d1;
                #pragma unroll
                for (int ix = 0; ix < WDIM; ++ix) {
                    const float ww  = wx[ix]*(wy[iy]*wz[iz]);
                    const size_t vidx = yoff + (size_t)(wx0+ix-lb1);
                    for (int q = 0; q < ncomp; ++q) {
                        const float ds = dscale[(size_t)i*ncomp + q];
                        atomicAdd(&cmat[(size_t)q*nvox + vidx].x, ds*cbre*ww);
                        atomicAdd(&cmat[(size_t)q*nvox + vidx].y, ds*cbim*ww);
                    }
                    float* rv = rho + vidx*(size_t)nrho;
                    for (int ir = 0; ir < nrho; ++ir) {
                        atomicAdd(&rv[ir], rscale[(size_t)i*nrho + ir]*ct*ww);
                    }
                }
            }
        }
    }
}

// Column-accumulation phase-B kernel (accumulate_covariance_columns / scatter_one_particle):
// consumes the CPU-built per-particle sample caches VERBATIM (windows, per-axis weights, plane
// values, ctf^2) so the interpolation weights are bit-identical to the CPU's exact-apod path.
// Per (particle, cached sample) thread: for every column s, replicate the dead-column skip, the
// noise-debias term at the column's q-voxel, and the 27-tap scatter into that halfset's
// Bcol/Hcol lattice. Layout matches Fortran (h,k,l,s): flat = vox + nvox*s.
__global__ void cols_scatter_kernel(
    const int*    __restrict__ cwin,   // [nb][maxc][3]
    const float*  __restrict__ cwx,    // [nb][maxc][3]
    const float*  __restrict__ cwy, const float* __restrict__ cwz,
    const float2* __restrict__ cpl,    // [nb][maxc]
    const float*  __restrict__ cct,    // [nb][maxc]
    const int*    __restrict__ ncache, // [nb]
    const unsigned char* __restrict__ vmask,  // [nb] 0 invalid, 1 even, 2 odd
    const float2* __restrict__ gcol,   // [nb][ncol]
    const float*  __restrict__ hcolv,  // [nb][ncol]
    const int*    __restrict__ qhkl,   // [ncol][3]
    int nb, int maxc, int ncol,
    int lb1, int lb2, int lb3, int ub1, int ub2, int ub3,
    int deadskip, int debias, float pf2, size_t nvox,
    float2* __restrict__ Bcol_e, float* __restrict__ Hcol_e,
    float2* __restrict__ Bcol_o, float* __restrict__ Hcol_o)
{
    const size_t tid = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (tid >= (size_t)nb*maxc) return;
    const int ib = (int)(tid / maxc);
    const int j  = (int)(tid % maxc);
    const unsigned char vm = vmask[ib];
    if (vm == 0 || j >= ncache[ib]) return;
    float2* Bcol = (vm == 1) ? Bcol_e : Bcol_o;
    float*  Hcol = (vm == 1) ? Hcol_e : Hcol_o;
    const int w1 = cwin[((size_t)ib*maxc + j)*3 + 0];
    const int w2 = cwin[((size_t)ib*maxc + j)*3 + 1];
    const int w3 = cwin[((size_t)ib*maxc + j)*3 + 2];
    const float* wx = cwx + ((size_t)ib*maxc + j)*3;
    const float* wy = cwy + ((size_t)ib*maxc + j)*3;
    const float* wz = cwz + ((size_t)ib*maxc + j)*3;
    const float2 pl = cpl[(size_t)ib*maxc + j];
    const float  ct = cct[(size_t)ib*maxc + j];
    const int d1 = ub1-lb1+1, d2 = ub2-lb2+1;
    for (int s = 0; s < ncol; ++s) {
        const float2 g = gcol[(size_t)ib*ncol + s];
        const float  h = hcolv[(size_t)ib*ncol + s];
        if (deadskip && g.x == 0.f && g.y == 0.f && h == 0.f) continue;
        // cg = pf2 * conj(g)
        const float cgre = pf2*g.x, cgim = -pf2*g.y;
        // noise-debias: stencil weight AT the column's q-voxel, when inside this sample's window
        float wqs = 0.f;
        if (debias) {
            const int qh = qhkl[3*s+0], qk = qhkl[3*s+1], ql = qhkl[3*s+2];
            if (qh >= w1 && qh <= w1+2 && qk >= w2 && qk <= w2+2 && ql >= w3 && ql <= w3+2)
                wqs = pf2 * wx[qh-w1]*wy[qk-w2]*wz[ql-w3] * ct;
        }
        // cnum = cg * pl
        const float cnre = cgre*pl.x - cgim*pl.y;
        const float cnim = cgre*pl.y + cgim*pl.x;
        for (int iz = 0; iz < 3; ++iz) {
            const size_t zoff = (size_t)(w3+iz-lb3)*d2*d1;
            for (int iy = 0; iy < 3; ++iy) {
                const size_t yoff = zoff + (size_t)(w2+iy-lb2)*d1;
                const float wyz = wy[iy]*wz[iz];
                for (int ix = 0; ix < 3; ++ix) {
                    const float w = wx[ix]*wyz;
                    const size_t idx = (size_t)s*nvox + yoff + (size_t)(w1+ix-lb1);
                    atomicAdd(&Bcol[idx].x, cnre*w - wqs*w);
                    atomicAdd(&Bcol[idx].y, cnim*w);
                    atomicAdd(&Hcol[idx],   h*w*ct);
                }
            }
        }
    }
}

// ---- fused columns: cache build + column gather on device ----
// The adjoint residual is the RESIDENT packed pair (cb_c/cb_ct after the E-step residual),
// so the cache entries are born in the same co_* pools the scatter kernel reads — the whole
// CPU cache build (mean projection, residual formation, KB stencils) and its 4x-contended
// uploads disappear. Entry ORDER within a particle's cache is nondeterministic (atomic
// append); every consumer is an order-free sum, so the gate is statistical, as for insert.
__global__ void cols_cache_kernel(
    const float2* __restrict__ plc, const float* __restrict__ plct,
    const float*  __restrict__ rotm, const unsigned char* __restrict__ vmask,
    int nb, int hs_lo, int hs_hi, int ks_lo, int khi,
    int nyq_disk, int iwinsz,
    int lb1, int lb2, int lb3, int ub1, int ub2, int ub3,
    float tiny, int maxc,
    int* __restrict__ cwin, float* __restrict__ cwx, float* __restrict__ cwy,
    float* __restrict__ cwz, float2* __restrict__ cpl, float* __restrict__ cct,
    int* __restrict__ ncache, int* __restrict__ ovfl)
{
    const int nh = hs_hi - hs_lo + 1;
    const int nk = khi  - ks_lo + 1;
    const size_t npts = (size_t)nh*nk;
    const size_t tid  = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (tid >= (size_t)nb*npts) return;
    const int ib = (int)(tid / npts);
    if (!vmask[ib]) return;
    const size_t pt = tid - (size_t)ib*npts;
    const int h = hs_lo + (int)(pt % nh);
    const int k = ks_lo + (int)(pt / nh);
    if (h*h + k*k > nyq_disk) return;
    // adjoint value: stored k<=0 half, Friedel for k>0 (packed unpadded layout)
    const size_t plsz = (size_t)(1 - ks_lo)*nh;
    float2 pv; float cv;
    if (k <= 0) {
        const size_t pidx = (size_t)ib*plsz + (size_t)(k - ks_lo)*nh + (h - hs_lo);
        pv = plc[pidx]; cv = plct[pidx];
    } else {
        const size_t pidx = (size_t)ib*plsz + (size_t)(-k - ks_lo)*nh + (-h - hs_lo);
        pv = plc[pidx]; pv.y = -pv.y; cv = plct[pidx];
    }
    if (fabsf(pv.x) + fabsf(pv.y) <= tiny && cv <= tiny) return;
    // loc(j) = h*R(1,j) + k*R(2,j), R Fortran column-major
    const float* R = rotm + (size_t)ib*9;
    const float lx = (float)h*R[0] + (float)k*R[1];
    const float ly = (float)h*R[3] + (float)k*R[4];
    const float lz = (float)h*R[6] + (float)k*R[7];
    const int wl1 = (int)nearbyintf(lx) - iwinsz;
    const int wl2 = (int)nearbyintf(ly) - iwinsz;
    const int wl3 = (int)nearbyintf(lz) - iwinsz;
    const int tw  = 2*iwinsz;
    if (wl1 < lb1 || wl2 < lb2 || wl3 < lb3 ||
        wl1 + tw > ub1 || wl2 + tw > ub2 || wl3 + tw > ub3) return;
    // per-axis normalized KB apod, the cov_kb_weights scheme
    float wxx[WDIM], wyy[WDIM], wzz[WDIM];
    float sx = 0.f, sy = 0.f, sz = 0.f;
    for (int t = 0; t < WDIM; ++t) {
        wxx[t] = kb_apod_fast((float)(wl1+t) - lx); sx += wxx[t];
        wyy[t] = kb_apod_fast((float)(wl2+t) - ly); sy += wyy[t];
        wzz[t] = kb_apod_fast((float)(wl3+t) - lz); sz += wzz[t];
    }
    const float eps = 1.1920929e-7f;
    if (fabsf(sx) > eps) { for (int t=0;t<WDIM;++t) wxx[t] /= sx; }
    else                 { for (int t=0;t<WDIM;++t) wxx[t] = 1.f/3.f; }
    if (fabsf(sy) > eps) { for (int t=0;t<WDIM;++t) wyy[t] /= sy; }
    else                 { for (int t=0;t<WDIM;++t) wyy[t] = 1.f/3.f; }
    if (fabsf(sz) > eps) { for (int t=0;t<WDIM;++t) wzz[t] /= sz; }
    else                 { for (int t=0;t<WDIM;++t) wzz[t] = 1.f/3.f; }
    const int j = atomicAdd(&ncache[ib], 1);
    if (j >= maxc) { atomicMax(ovfl, 1); return; }
    const size_t e = (size_t)ib*maxc + j;
    cwin[e*3+0] = wl1; cwin[e*3+1] = wl2; cwin[e*3+2] = wl3;
    for (int t = 0; t < WDIM; ++t) {
        cwx[e*3+t] = wxx[t]; cwy[e*3+t] = wyy[t]; cwz[e*3+t] = wzz[t];
    }
    cpl[e] = pv; cct[e] = cv;
}

__global__ void cols_gather_kernel(
    const int* __restrict__ cwin, const float* __restrict__ cwx,
    const float* __restrict__ cwy, const float* __restrict__ cwz,
    const float2* __restrict__ cpl, const float* __restrict__ cct,
    const int* __restrict__ ncache, const unsigned char* __restrict__ vmask,
    const int* __restrict__ lookup,
    int nb, int maxc, int ncol,
    int lb1, int lb2, int lb3, int ub1, int ub2, int ub3,
    float2* __restrict__ gcol, float* __restrict__ hcolv)
{
    const size_t tid = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (tid >= (size_t)nb*maxc) return;
    const int ib = (int)(tid / maxc);
    const int j  = (int)(tid % maxc);
    if (!vmask[ib] || j >= ncache[ib]) return;
    const size_t e = (size_t)ib*maxc + j;
    const int w1 = cwin[e*3+0], w2 = cwin[e*3+1], w3 = cwin[e*3+2];
    const float* wx = cwx + e*3;
    const float* wy = cwy + e*3;
    const float* wz = cwz + e*3;
    const float2 pl = cpl[e];
    const float  ct = cct[e];
    const int d1 = ub1-lb1+1, d2 = ub2-lb2+1;
    for (int iz = 0; iz < 3; ++iz) {
        const size_t zoff = (size_t)(w3+iz-lb3)*d2*d1;
        for (int iy = 0; iy < 3; ++iy) {
            const size_t yoff = zoff + (size_t)(w2+iy-lb2)*d1;
            const float wyz = wy[iy]*wz[iz];
            for (int ix = 0; ix < 3; ++ix) {
                const int s = lookup[yoff + (size_t)(w1+ix-lb1)];
                if (s == 0) continue;
                const float w = wx[ix]*wyz;
                atomicAdd(&gcol[(size_t)ib*ncol + (s-1)].x, w*pl.x);
                atomicAdd(&gcol[(size_t)ib*ncol + (s-1)].y, w*pl.y);
                atomicAdd(&hcolv[(size_t)ib*ncol + (s-1)],  w*ct);
            }
        }
    }
}

} // namespace

#define CKI(x) do { cudaError_t e = (x); if (e != cudaSuccess) { \
    fprintf(stderr, "flex_gpu CUDA error %s at %s:%d\n", cudaGetErrorString(e), __FILE__, __LINE__); \
    return 1; } } while (0)

// Grow-only per-stage device buffer pools. The batch entry points used to cudaMalloc/cudaFree
// their staging buffers on EVERY batch; with the heavy kernels gone (banked adjoint) that
// allocator churn dominated: cudaFree synchronizes device-wide, and with 4 worker processes
// contending, ~13 alloc/free pairs x ~170 batches cost ~20 s of a 29 s stage. Buffers now grow
// to the high-water mark and are released once in the stage's *_free().
namespace {
struct DBuf { void* p = nullptr; size_t cap = 0; };
inline int dbuf_ensure(DBuf& b, size_t bytes)
{
    if (bytes <= b.cap) return 0;
    if (b.p) cudaFree(b.p);
    b.p = nullptr; b.cap = 0;
    cudaError_t e = cudaMalloc(&b.p, bytes);
    if (e != cudaSuccess) {
        fprintf(stderr, "flex_gpu CUDA error %s (dbuf_ensure %zu bytes)\n", cudaGetErrorString(e), bytes);
        b.p = nullptr; return 1;
    }
    b.cap = bytes; return 0;
}
inline void dbuf_release(DBuf& b) { if (b.p) cudaFree(b.p); b.p = nullptr; b.cap = 0; }
// pinned host mirrors for the big plane uploads: pageable H2D staged by the driver runs at
// ~2-4 GB/s and serializes the four workers on the PCIe link; memcpy into a pinned mirror
// (host-bandwidth) + async DMA restores full link speed on 1/pf^2 of the old payload
struct PBuf { void* p = nullptr; size_t cap = 0; };
inline int pbuf_ensure(PBuf& b, size_t bytes)
{
    if (bytes <= b.cap) return 0;
    if (b.p) cudaFreeHost(b.p);
    b.p = nullptr; b.cap = 0;
    cudaError_t e = cudaMallocHost(&b.p, bytes);
    if (e != cudaSuccess) {
        fprintf(stderr, "flex_gpu CUDA error %s (pbuf_ensure %zu bytes)\n", cudaGetErrorString(e), bytes);
        b.p = nullptr; return 1;
    }
    b.cap = bytes; return 0;
}
inline void pbuf_release(PBuf& b) { if (b.p) cudaFreeHost(b.p); b.p = nullptr; b.cap = 0; }
// Segment-chunk width for the banked M-step. This is a pure BLOCKING factor: the inner loop runs
// one GEMM pair PER SEGMENT for any value, so the arithmetic is identical and only the number of
// splat launches (and the parallelism inside each) changes. It is also the DOMINANT VRAM term --
// cb_cmbr = NSEGC*nrho*npts*4 -- and both grow with the box and quadratically with neigs, so a
// hardcoded 64 is a footprint that scales with the science settings and OOMs without warning
// (measured: 9.4 GB for that one buffer at npts=65792/nrho=600, killing nparts=4 on a 32 GB card).
// Size it from a byte BUDGET instead, so per-worker VRAM is bounded and predictable when several
// workers share one card. Even one segment per chunk leaves npts work items per splat launch
// (~65 k here = 257 blocks), which is ample occupancy, so the budget costs launches, not throughput.
constexpr int NSEGC_MAX = 64;
inline int nsegc_for(size_t npts, int ncomp, int nrho)
{
    static long budget_mb = -1;
    if (budget_mb < 0) {
        budget_mb = 768;                       // per worker, for cb_cmbc + cb_cmbr together
        const char* e = getenv("SIMPLE_FLEX_GPU_BANK_MB");
        if (e && *e) { long v = atol(e); if (v > 0) budget_mb = v; }
    }
    const size_t per = ((size_t)ncomp*sizeof(float2) + (size_t)nrho*sizeof(float)) * npts;
    if (per == 0) return NSEGC_MAX;
    size_t n = ((size_t)budget_mb << 20) / per;
    if (n < 1) n = 1;
    if (n > (size_t)NSEGC_MAX) n = (size_t)NSEGC_MAX;
    static long last = -1;
    if ((long)n != last) {
        last = (long)n;
        fprintf(stderr, "flex_gpu bank chunk: NSEGC=%d (npts=%zu ncomp=%d nrho=%d) -> "
            "cb_cmbc+cb_cmbr %.2f GB [budget %ld MB, SIMPLE_FLEX_GPU_BANK_MB]\n",
            (int)n, npts, ncomp, nrho, (double)(n*per)/1073741824.0, budget_mb);
        fflush(stderr);
    }
    return (int)n;
}
inline int hstage(PBuf& pin, void* dev, const void* src, size_t bytes)
{
    if (pbuf_ensure(pin, bytes)) return 1;
    memcpy(pin.p, src, bytes);
    cudaError_t e = cudaMemcpyAsync(dev, pin.p, bytes, cudaMemcpyHostToDevice);
    if (e != cudaSuccess) {
        fprintf(stderr, "flex_gpu CUDA error %s (hstage %zu bytes)\n", cudaGetErrorString(e), bytes);
        return 1;
    }
    return 0;
}
// insert (states) pool
DBuf in_c, in_ct, in_rm, in_ds, in_rs, in_act, in_na, in_nd;
// per-particle coupled pool
DBuf cp_c, cp_ct, cp_rm, cp_ds, cp_rs, cp_vm, cp_nd;
// banked coupled pool
DBuf cb_c, cb_ct, cb_ca, cb_sa, cb_ds, cb_rs, cb_vm, cb_al, cb_alct, cb_cmbc, cb_cmbr, cb_sd;
DBuf cb_cy, cb_t;   // device prep: separated whitened pair (image plane, CTF plane)
// packed-plane geometry registered by the device prep / E-step uploads; shared by every
// resident consumer (fused E-step, banked M-step, resident insert)
int ge_hs_lo = 0, ge_hs_hi = 0, ge_ks_lo = 0, ge_khi = 0, ge_nrec = 0, ge_nyq2 = 0;
// columns phase-B pool
DBuf co_cw, co_wx, co_wy, co_wz, co_pl, co_ct, co_nc, co_vm, co_gc, co_hv, co_qh;
DBuf co_rm, co_ovfl;          // fused columns: rotation matrices + cache-overflow flag
int *dc_lookup = nullptr;     // fused columns: column lookup volume (uploaded once per stage)
// shared pinned mirrors for the plane payloads (stages are sequential within a process, and the
// high-water size is ~tens of MB, so these live for the process lifetime)
PBuf pin_c, pin_ct;
// pinned mirrors for the columns cache payload (the fattest upload in the pipeline)
PBuf pin_cw, pin_wx, pin_wy, pin_wz, pin_pl, pin_pct;
}
#define DBUF(b, n) do { if (dbuf_ensure(b, (n))) return 4; } while (0)

extern "C" {

int flex_gpu_device_count(void)
{
    int n = 0;
    if (cudaGetDeviceCount(&n) != cudaSuccess) return 0;
    return n;
}

int flex_gpu_insert_begin(int ncomp, long nvox)
{
    if (d_cmat) { cudaFree(d_cmat); d_cmat = nullptr; }
    if (d_rho)  { cudaFree(d_rho);  d_rho  = nullptr; }
    g_nvox = (size_t)nvox; g_ncomp = ncomp;
    // insertion accumulators: ncomp x nvox x 12 B. In the state stage ncomp is 2*nstates
    // (even+odd), so this scales with the STATE COUNT and box_rec^3 and is independent of neigs
    // -- it is the peak-VRAM term of the whole program at large box_rec (measured 10076
    // box_rec=256, 24 states: ~4.9 GB of a 7.4 GB per-worker peak). Unlike the M-step bank it
    // can only be reduced by processing states in batches, which costs extra particle passes.
    fprintf(stderr, "flex_gpu insert accumulators: d_cmat %.2f GB + d_rho %.2f GB "
        "(ncomp=%d nvox=%zu)\n",
        (double)((size_t)ncomp*g_nvox*sizeof(float2))/1073741824.0,
        (double)((size_t)ncomp*g_nvox*sizeof(float))/1073741824.0, ncomp, g_nvox);
    fflush(stderr);
    CKI(cudaMalloc(&d_cmat, (size_t)ncomp*g_nvox*sizeof(float2)));
    CKI(cudaMalloc(&d_rho,  (size_t)ncomp*g_nvox*sizeof(float)));
    CKI(cudaMemset(d_cmat, 0, (size_t)ncomp*g_nvox*sizeof(float2)));
    CKI(cudaMemset(d_rho,  0, (size_t)ncomp*g_nvox*sizeof(float)));
    return 0;
}

int flex_gpu_insert_batch(
    const void* pl_c, const void* pl_ct, const void* rotm,
    const void* dscale, const void* rscale, const void* act, const void* nact,
    const void* nyqdisk, int nrec, int ncomp, int nsym,
    int hpd_lo, int hpd_hi, int kpd_lo,
    int hlo, int hhi, int klo, int khi, int pf,
    const int* lb, const int* ub, float tiny)
{
    if (!d_cmat || ncomp != g_ncomp) return 2;
    const int hpdn = hpd_hi - hpd_lo + 1;
    const size_t plsz = (size_t)(1 - kpd_lo)*hpdn;
    DBUF(in_c,  (size_t)nrec*plsz*sizeof(float2));
    DBUF(in_ct, (size_t)nrec*plsz*sizeof(float));
    DBUF(in_rm, (size_t)nrec*nsym*9*sizeof(float));
    DBUF(in_ds, (size_t)nrec*ncomp*sizeof(float));
    DBUF(in_rs, (size_t)nrec*ncomp*sizeof(float));
    DBUF(in_act,(size_t)nrec*ncomp*sizeof(int));
    DBUF(in_na, (size_t)nrec*sizeof(int));
    DBUF(in_nd, (size_t)nrec*sizeof(int));
    float2 *dp_c = (float2*)in_c.p; float *dp_ct = (float*)in_ct.p, *dp_rm = (float*)in_rm.p,
        *dp_ds = (float*)in_ds.p, *dp_rs = (float*)in_rs.p;
    int *dp_act = (int*)in_act.p, *dp_na = (int*)in_na.p, *dp_nd = (int*)in_nd.p;
    if (hstage(pin_c,  dp_c,  pl_c,  (size_t)nrec*plsz*sizeof(float2))) return 5;
    if (hstage(pin_ct, dp_ct, pl_ct, (size_t)nrec*plsz*sizeof(float)))  return 5;
    CKI(cudaMemcpy(dp_rm, rotm,  (size_t)nrec*nsym*9*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dp_ds, dscale,(size_t)nrec*ncomp*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dp_rs, rscale,(size_t)nrec*ncomp*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dp_act,act,   (size_t)nrec*ncomp*sizeof(int),   cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dp_na, nact,  (size_t)nrec*sizeof(int),         cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dp_nd, nyqdisk,(size_t)nrec*sizeof(int),        cudaMemcpyHostToDevice));
    const int nh = hhi - hlo + 1, nk = khi - klo + 1;
    const size_t nwork = (size_t)nrec*nh*nk;
    const int nblk = (int)((nwork + 255)/256);
    insert_multi_kernel<<<nblk, 256>>>(dp_c, dp_ct, dp_rm, dp_ds, dp_rs, dp_act, dp_na,
        dp_nd, nrec, ncomp, nsym, hpd_lo, hpd_hi, kpd_lo, hlo, hhi, klo, khi, pf,
        lb[0], lb[1], lb[2], ub[0], ub[1], ub[2], tiny, g_nvox, d_cmat, d_rho);
    CKI(cudaGetLastError());
    CKI(cudaDeviceSynchronize());
    return 0;
}

// resident variant: consume the device-prepped packed planes (cb_c/cb_ct, geometry as
// registered by flex_gpu_prep_batch) directly -- zero plane traffic. Records are NOT
// compacted: dead records carry nact=0 and vmask'd zero planes.
int flex_gpu_insert_batch_resident(
    const void* rotm, const void* dscale, const void* rscale, const void* act,
    const void* nact, const void* nyqdisk, int nrec, int ncomp, int nsym,
    const int* lb, const int* ub, float tiny)
{
    if (!d_cmat || ncomp != g_ncomp) return 2;
    if (!cb_c.p || !cb_ct.p || ge_nrec < nrec) return 3;
    const int hs_lo = ge_hs_lo, hs_hi = ge_hs_hi, ks_lo = ge_ks_lo, khi = ge_khi;
    DBUF(in_rm, (size_t)nrec*nsym*9*sizeof(float));
    DBUF(in_ds, (size_t)nrec*ncomp*sizeof(float));
    DBUF(in_rs, (size_t)nrec*ncomp*sizeof(float));
    DBUF(in_act,(size_t)nrec*ncomp*sizeof(int));
    DBUF(in_na, (size_t)nrec*sizeof(int));
    DBUF(in_nd, (size_t)nrec*sizeof(int));
    float *dp_rm = (float*)in_rm.p, *dp_ds = (float*)in_ds.p, *dp_rs = (float*)in_rs.p;
    int *dp_act = (int*)in_act.p, *dp_na = (int*)in_na.p, *dp_nd = (int*)in_nd.p;
    CKI(cudaMemcpy(dp_rm, rotm,  (size_t)nrec*nsym*9*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dp_ds, dscale,(size_t)nrec*ncomp*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dp_rs, rscale,(size_t)nrec*ncomp*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dp_act,act,   (size_t)nrec*ncomp*sizeof(int),   cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dp_na, nact,  (size_t)nrec*sizeof(int),         cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dp_nd, nyqdisk,(size_t)nrec*sizeof(int),        cudaMemcpyHostToDevice));
    const int nh = hs_hi - hs_lo + 1, nk = khi - ks_lo + 1;
    const size_t nwork = (size_t)nrec*nh*nk;
    const int nblk = (int)((nwork + 255)/256);
    insert_multi_kernel<<<nblk, 256>>>((const float2*)cb_c.p, (const float*)cb_ct.p,
        dp_rm, dp_ds, dp_rs, dp_act, dp_na, dp_nd, nrec, ncomp, nsym,
        hs_lo, hs_hi, ks_lo, hs_lo, hs_hi, ks_lo, khi, OSMPL_PAD_FAC_DEV,
        lb[0], lb[1], lb[2], ub[0], ub[1], ub[2], tiny, g_nvox, d_cmat, d_rho);
    CKI(cudaGetLastError());
    CKI(cudaDeviceSynchronize());
    return 0;
}

int flex_gpu_insert_fetch(int q1, void* cmat_host, void* rho_host)
{
    if (!d_cmat || q1 < 1 || q1 > g_ncomp) return 2;
    const size_t q = (size_t)(q1 - 1);
    CKI(cudaMemcpy(cmat_host, d_cmat + q*g_nvox, g_nvox*sizeof(float2), cudaMemcpyDeviceToHost));
    CKI(cudaMemcpy(rho_host,  d_rho  + q*g_nvox, g_nvox*sizeof(float),  cudaMemcpyDeviceToHost));
    return 0;
}

int flex_gpu_insert_free(void)
{
    if (d_cmat) { cudaFree(d_cmat); d_cmat = nullptr; }
    if (d_rho)  { cudaFree(d_rho);  d_rho  = nullptr; }
    g_nvox = 0; g_ncomp = 0;
    dbuf_release(in_c); dbuf_release(in_ct); dbuf_release(in_rm); dbuf_release(in_ds);
    dbuf_release(in_rs); dbuf_release(in_act); dbuf_release(in_na); dbuf_release(in_nd);
    return 0;
}

int flex_gpu_coupled_begin(int ncomp, int nrho, long nvox)
{
    if (d_cmat2) { cudaFree(d_cmat2); d_cmat2 = nullptr; }
    if (d_rho2)  { cudaFree(d_rho2);  d_rho2  = nullptr; }
    g2_nvox = (size_t)nvox; g2_ncomp = ncomp; g2_nrho = nrho;
    // resident accumulators: the OTHER half of the per-worker footprint, and unlike the bank
    // they cannot be chunked (they ARE the output). Report so a VRAM budget can be planned.
    {
        const double GB = 1.0/1073741824.0;
        fprintf(stderr, "flex_gpu resident accumulators: d_cmat2 %.2f GB + d_rho2 %.2f GB "
            "(ncomp=%d nrho=%d nvox=%zu)\n",
            (double)((size_t)ncomp*g2_nvox*sizeof(float2))*GB,
            (double)((size_t)nrho*g2_nvox*sizeof(float))*GB, ncomp, nrho, g2_nvox);
        fflush(stderr);
    }
    CKI(cudaMalloc(&d_cmat2, (size_t)ncomp*g2_nvox*sizeof(float2)));
    CKI(cudaMalloc(&d_rho2,  (size_t)nrho*g2_nvox*sizeof(float)));
    CKI(cudaMemset(d_cmat2, 0, (size_t)ncomp*g2_nvox*sizeof(float2)));
    CKI(cudaMemset(d_rho2,  0, (size_t)nrho*g2_nvox*sizeof(float)));
    return 0;
}

int flex_gpu_coupled_batch(
    const void* pl_c, const void* pl_ct, const void* rotm,
    const void* dscale, const void* rscale, const void* vmask, const void* nyqdisk,
    int nrec, int ncomp, int nrho, int nsym,
    int hpd_lo, int hpd_hi, int kpd_lo,
    int hlo, int hhi, int klo, int khi, int pf,
    const int* lb, const int* ub, float tiny)
{
    if (!d_cmat2 || ncomp != g2_ncomp || nrho != g2_nrho) return 2;
    const int hpdn = hpd_hi - hpd_lo + 1;
    const size_t plsz = (size_t)(1 - kpd_lo)*hpdn;
    DBUF(cp_c,  (size_t)nrec*plsz*sizeof(float2));
    DBUF(cp_ct, (size_t)nrec*plsz*sizeof(float));
    DBUF(cp_rm, (size_t)nrec*nsym*9*sizeof(float));
    DBUF(cp_ds, (size_t)nrec*ncomp*sizeof(float));
    DBUF(cp_rs, (size_t)nrec*nrho*sizeof(float));
    DBUF(cp_vm, (size_t)nrec*sizeof(unsigned char));
    DBUF(cp_nd, (size_t)nrec*sizeof(int));
    float2 *dp_c = (float2*)cp_c.p; float *dp_ct = (float*)cp_ct.p, *dp_rm = (float*)cp_rm.p,
        *dp_ds = (float*)cp_ds.p, *dp_rs = (float*)cp_rs.p;
    int *dp_nd = (int*)cp_nd.p; unsigned char *dp_vm = (unsigned char*)cp_vm.p;
    if (hstage(pin_c,  dp_c,  pl_c,  (size_t)nrec*plsz*sizeof(float2))) return 5;
    if (hstage(pin_ct, dp_ct, pl_ct, (size_t)nrec*plsz*sizeof(float)))  return 5;
    CKI(cudaMemcpy(dp_rm, rotm,  (size_t)nrec*nsym*9*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dp_ds, dscale,(size_t)nrec*ncomp*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dp_rs, rscale,(size_t)nrec*nrho*sizeof(float),  cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dp_vm, vmask, (size_t)nrec*sizeof(unsigned char), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dp_nd, nyqdisk,(size_t)nrec*sizeof(int),        cudaMemcpyHostToDevice));
    const int nh = hhi - hlo + 1, nk = khi - klo + 1;
    const size_t nwork = (size_t)nrec*nh*nk;
    const int nblk = (int)((nwork + 255)/256);
    insert_coupled_kernel<<<nblk, 256>>>(dp_c, dp_ct, dp_rm, dp_ds, dp_rs, dp_nd, dp_vm,
        nrec, ncomp, nrho, nsym, hpd_lo, hpd_hi, kpd_lo, hlo, hhi, klo, khi, pf,
        lb[0], lb[1], lb[2], ub[0], ub[1], ub[2], tiny, g2_nvox, d_cmat2, d_rho2);
    CKI(cudaGetLastError());
    CKI(cudaDeviceSynchronize());
    return 0;
}

int flex_gpu_coupled_fetch_cmat(int q1, void* cmat_host)
{
    if (!d_cmat2 || q1 < 1 || q1 > g2_ncomp) return 2;
    CKI(cudaMemcpy(cmat_host, d_cmat2 + (size_t)(q1-1)*g2_nvox, g2_nvox*sizeof(float2),
        cudaMemcpyDeviceToHost));
    return 0;
}

int flex_gpu_coupled_fetch_rho(void* rho_host)
{
    if (!d_rho2) return 2;
    CKI(cudaMemcpy(rho_host, d_rho2, (size_t)g2_nrho*g2_nvox*sizeof(float),
        cudaMemcpyDeviceToHost));
    return 0;
}

int flex_gpu_coupled_free(void)
{
    if (d_cmat2) { cudaFree(d_cmat2); d_cmat2 = nullptr; }
    if (d_rho2)  { cudaFree(d_rho2);  d_rho2  = nullptr; }
    g2_nvox = 0; g2_ncomp = 0; g2_nrho = 0;
    dbuf_release(cp_c); dbuf_release(cp_ct); dbuf_release(cp_rm); dbuf_release(cp_ds);
    dbuf_release(cp_rs); dbuf_release(cp_vm); dbuf_release(cp_nd);
    return 0;
}

}   // extern "C" (reopened after the banked-adjoint section)

// ---- Banked coupled adjoint (probe M-step) ----
// The per-particle coupled insert is atomic-bound: every record splats into ALL ncomp cmplx +
// nrho density rows (dense E[zz'] coupling), ~1.4e12 atomicAdds per probe iteration at 84k
// records x 304 rows. Reformulation: factor each record's rotation R_i = R_bank(dir) * Rz(psi),
// resample the record's plane by psi IN-PLANE (2D KB gather on the same unit-tap lattice the
// polar former uses -- polar_interp_plane's exact scheme), GEMM-combine the aligned planes of
// each bank direction's records (scales become the GEMM's A-matrix), and splat each direction's
// combined planes ONCE at R_bank. Atomic volume drops by the mean records-per-direction; the
// heavy work becomes sgemm. The direction quantization is the SAME approximation family as the
// gated banked polar forward, and is science-gated (Ribosembly ARI band), not byte-compared.
// Callers must pre-sort records by direction id so segments are contiguous within a batch;
// a direction spanning two batches just splats twice (additive, exact).
namespace {
float *d_cb_rmat = nullptr;   // [ndir][9] bank rotations, same column-major layout as rotm
int    g_cb_ndir = 0;
cublasHandle_t g_hb_cb = nullptr;
// stage profile across a bank's lifetime (upload/align/gemm/splat in ms of stream time,
// segment stats); printed and reset by flex_gpu_coupled_bank_free
double g_cbt_up = 0., g_cbt_al = 0., g_cbt_gm = 0., g_cbt_sp = 0.;
long   g_cb_nbatch = 0, g_cb_nseg_tot = 0; int g_cb_nseg_max = 0;

// 2D in-plane KB resample of the packed half-planes, mirroring polar_interp_plane: taps at UNIT
// spacing in (h,k) around the rotated position, values read at pf-multiples on the stored padded
// half-plane with per-tap Friedel, per-axis normalized apod weights, OOB taps dropped WITHOUT
// renormalization. pf^2 is folded into the cmplx output (the insert convention); ctfsq is not.
__global__ void align_rot2d_kernel(
    const float2* __restrict__ pl_c, const float* __restrict__ pl_ct,
    const float*  __restrict__ cav, const float* __restrict__ sav,
    const unsigned char* __restrict__ vmask,
    int nrec, int hpd_lo, int hpd_hi, int kpd_lo,
    int hlo, int hhi, int klo, int khi, int pf, int nyq2,
    float2* __restrict__ al_c, float* __restrict__ al_ct)
{
    const int nh = hhi - hlo + 1;
    const int nk = khi - klo + 1;
    const size_t npts = (size_t)nh*nk;
    const size_t tid  = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (tid >= (size_t)nrec*npts) return;
    const int i  = (int)(tid / npts);
    const size_t pt = tid - (size_t)i*npts;
    const int h  = hlo + (int)(pt % nh);
    const int k  = klo + (int)(pt / nh);
    float2 acc = make_float2(0.f, 0.f); float acct = 0.f;
    al_c[tid] = acc; al_ct[tid] = acct;
    if (!vmask[i]) return;
    if (h*h + k*k > nyq2) return;
    // sample the particle plane at the bank-frame position rotated by the record's in-plane
    // angle: hu = h*ca + k*sa, ku = -h*sa + k*ca (the polar former's phi+alpha convention)
    const float ca = cav[i], sa = sav[i];
    const float hu =  h*ca + k*sa;
    const float ku = -h*sa + k*ca;
    const int wx0 = __float2int_rn(hu) - IWINSZ;
    const int wy0 = __float2int_rn(ku) - IWINSZ;
    float wx[WDIM], wy[WDIM];
    float sx = 0.f, sy = 0.f;
    #pragma unroll
    for (int j = 0; j < WDIM; ++j) {
        wx[j] = kb_apod_fast((float)(wx0+j) - hu); sx += wx[j];
        wy[j] = kb_apod_fast((float)(wy0+j) - ku); sy += wy[j];
    }
    const float eps = 1.1920929e-07f;
    if (fabsf(sx) > eps) { const float s = 1.f/sx; wx[0]*=s; wx[1]*=s; wx[2]*=s; }
    else                 { wx[0]=wx[1]=wx[2]=1.f/3.f; }
    if (fabsf(sy) > eps) { const float s = 1.f/sy; wy[0]*=s; wy[1]*=s; wy[2]*=s; }
    else                 { wy[0]=wy[1]=wy[2]=1.f/3.f; }
    // packed unpadded taps: unit spacing in (h,k), Friedel + stored-range bounds per tap,
    // exactly polar_interp_plane's scheme on the packed lattice
    const int hpdn = hpd_hi - hpd_lo + 1;
    const size_t plsz = (size_t)(1 - kpd_lo)*hpdn;
    const float2* pli = pl_c  + (size_t)i*plsz;
    const float*  pci = pl_ct + (size_t)i*plsz;
    #pragma unroll
    for (int iy = 0; iy < WDIM; ++iy) {
        const int ky = wy0 + iy;
        #pragma unroll
        for (int ix = 0; ix < WDIM; ++ix) {
            const int hx = wx0 + ix;
            const float w = wx[ix]*wy[iy];
            float2 cv; float cc;
            if (ky <= 0) {
                if (hx < hpd_lo || hx > hpd_hi || ky < kpd_lo) continue;
                const size_t pidx = (size_t)(ky - kpd_lo)*hpdn + (hx - hpd_lo);
                cv = pli[pidx]; cc = pci[pidx];
            } else {
                if (-hx < hpd_lo || -hx > hpd_hi || -ky < kpd_lo) continue;
                const size_t pidx = (size_t)(-ky - kpd_lo)*hpdn + (-hx - hpd_lo);
                const float2 c = pli[pidx];
                cv = make_float2(c.x, -c.y); cc = pci[pidx];
            }
            acc.x += w*cv.x; acc.y += w*cv.y; acct += w*cc;
        }
    }
    const float pf2 = (float)(pf*pf);
    al_c[tid]  = make_float2(pf2*acc.x, pf2*acc.y);
    al_ct[tid] = acct;
}

// Splat one chunk of direction segments: thread per (segment, plane point), tap geometry
// computed ONCE per point at the bank rotation, then all ncomp+nrho combined values accumulated
// through the same 27 taps. Identical bounds/window/apod arithmetic to insert_coupled_kernel.
__global__ void splat_bank_kernel(
    const float2* __restrict__ cmb_c,  // [nsegc][ncomp][npts]
    const float*  __restrict__ cmb_r,  // [nsegc][nrho][npts]
    const float*  __restrict__ bank_rm,
    const int*    __restrict__ seg_dir, // [nsegc] 1-based bank ids
    int nsegc, int ncomp, int nrho,
    int hlo, int hhi, int klo, int khi, int nyq2,
    int lb1, int lb2, int lb3, int ub1, int ub2, int ub3,
    size_t nvox,
    float2* __restrict__ cmat, float* __restrict__ rho)
{
    const int nh = hhi - hlo + 1;
    const int nk = khi - klo + 1;
    const size_t npts = (size_t)nh*nk;
    const size_t tid  = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (tid >= (size_t)nsegc*npts) return;
    const int js = (int)(tid / npts);
    const size_t pt = tid - (size_t)js*npts;
    const int h  = hlo + (int)(pt % nh);
    const int k  = klo + (int)(pt / nh);
    if (h*h + k*k > nyq2) return;
    const float* m = bank_rm + (size_t)(seg_dir[js] - 1)*9;
    const float lx = h*m[0] + k*m[1];
    const float ly = h*m[3] + k*m[4];
    const float lz = h*m[6] + k*m[7];
    const int wx0 = __float2int_rn(lx) - IWINSZ;
    const int wy0 = __float2int_rn(ly) - IWINSZ;
    const int wz0 = __float2int_rn(lz) - IWINSZ;
    if (wx0 + 2*IWINSZ < 0) return;
    if (wx0 < lb1 || wy0 < lb2 || wz0 < lb3) return;
    if (wx0+2*IWINSZ > ub1 || wy0+2*IWINSZ > ub2 || wz0+2*IWINSZ > ub3) return;
    float wx[WDIM], wy[WDIM], wz[WDIM];
    float sx = 0.f, sy = 0.f, sz = 0.f;
    #pragma unroll
    for (int j = 0; j < WDIM; ++j) {
        wx[j] = kb_apod_fast((float)(wx0+j) - lx); sx += wx[j];
        wy[j] = kb_apod_fast((float)(wy0+j) - ly); sy += wy[j];
        wz[j] = kb_apod_fast((float)(wz0+j) - lz); sz += wz[j];
    }
    const float eps = 1.1920929e-07f;
    if (fabsf(sx) > eps) { const float s = 1.f/sx; wx[0]*=s; wx[1]*=s; wx[2]*=s; }
    else                 { wx[0]=wx[1]=wx[2]=1.f/3.f; }
    if (fabsf(sy) > eps) { const float s = 1.f/sy; wy[0]*=s; wy[1]*=s; wy[2]*=s; }
    else                 { wy[0]=wy[1]=wy[2]=1.f/3.f; }
    if (fabsf(sz) > eps) { const float s = 1.f/sz; wz[0]*=s; wz[1]*=s; wz[2]*=s; }
    else                 { wz[0]=wz[1]=wz[2]=1.f/3.f; }
    const int d1 = ub1-lb1+1, d2 = ub2-lb2+1;
    const float2* vc = cmb_c + (size_t)js*ncomp*npts;
    const float*  vr = cmb_r + (size_t)js*nrho*npts;
    #pragma unroll
    for (int iz = 0; iz < WDIM; ++iz) {
        const size_t zoff = (size_t)(wz0+iz-lb3)*d2*d1;
        #pragma unroll
        for (int iy = 0; iy < WDIM; ++iy) {
            const size_t yoff = zoff + (size_t)(wy0+iy-lb2)*d1;
            #pragma unroll
            for (int ix = 0; ix < WDIM; ++ix) {
                const float ww  = wx[ix]*(wy[iy]*wz[iz]);
                const size_t vidx = yoff + (size_t)(wx0+ix-lb1);
                for (int q = 0; q < ncomp; ++q) {
                    const float2 v = vc[(size_t)q*npts + pt];
                    if (v.x != 0.f || v.y != 0.f) {
                        atomicAdd(&cmat[(size_t)q*nvox + vidx].x, v.x*ww);
                        atomicAdd(&cmat[(size_t)q*nvox + vidx].y, v.y*ww);
                    }
                }
                float* rv = rho + vidx*(size_t)nrho;
                for (int ir = 0; ir < nrho; ++ir) {
                    const float v = vr[(size_t)ir*npts + pt];
                    if (v != 0.f) atomicAdd(&rv[ir], v*ww);
                }
            }
        }
    }
}
}   // anon namespace

extern "C" {

int flex_gpu_coupled_bank(const void* rmat, int ndir)
{
    if (d_cb_rmat) { cudaFree(d_cb_rmat); d_cb_rmat = nullptr; }
    g_cb_ndir = ndir;
    CKI(cudaMalloc(&d_cb_rmat, (size_t)ndir*9*sizeof(float)));
    CKI(cudaMemcpy(d_cb_rmat, rmat, (size_t)ndir*9*sizeof(float), cudaMemcpyHostToDevice));
    if (!g_hb_cb) {
        if (cublasCreate(&g_hb_cb) != CUBLAS_STATUS_SUCCESS) { g_hb_cb = nullptr; return 3; }
    }
    return 0;
}

int flex_gpu_coupled_batch_banked(
    const void* pl_c, const void* pl_ct, const void* cav, const void* sav,
    const void* dscale, const void* rscale, const void* vmask, const int* dirid,
    int nrec, int ncomp, int nrho,
    int hpd_lo, int hpd_hi, int kpd_lo,
    int hlo, int hhi, int klo, int khi, int pf, int nyq2,
    const int* lb, const int* ub)
{
    if (!d_cmat2 || ncomp != g2_ncomp || nrho != g2_nrho) return 2;
    if (!d_cb_rmat || !g_hb_cb) return 2;
    const int hpdn = hpd_hi - hpd_lo + 1;
    const size_t plsz = (size_t)(1 - kpd_lo)*hpdn;
    const int nh = hhi - hlo + 1, nk = khi - klo + 1;
    const size_t npts = (size_t)nh*nk;
    // segments of equal dirid (records pre-sorted by the caller)
    const int NSEGC = nsegc_for(npts, ncomp, nrho);
    int nseg = 0;
    int *seg_dir_h = (int*)malloc((size_t)nrec*sizeof(int));
    int *seg_beg_h = (int*)malloc((size_t)nrec*sizeof(int));
    int *seg_cnt_h = (int*)malloc((size_t)nrec*sizeof(int));
    for (int i = 0; i < nrec; ++i) {
        if (nseg > 0 && dirid[i] == seg_dir_h[nseg-1]) { seg_cnt_h[nseg-1]++; continue; }
        seg_dir_h[nseg] = dirid[i]; seg_beg_h[nseg] = i; seg_cnt_h[nseg] = 1; nseg++;
    }
    g_cb_nbatch++; g_cb_nseg_tot += nseg;
    if (nseg > g_cb_nseg_max) g_cb_nseg_max = nseg;
    cudaEvent_t ev0, ev1, ev2, ev3, ev4;
    cudaEventCreate(&ev0); cudaEventCreate(&ev1); cudaEventCreate(&ev2);
    cudaEventCreate(&ev3); cudaEventCreate(&ev4);
    cudaEventRecord(ev0);
    DBUF(cb_c,   (size_t)nrec*plsz*sizeof(float2));
    DBUF(cb_ct,  (size_t)nrec*plsz*sizeof(float));
    DBUF(cb_ca,  (size_t)nrec*sizeof(float));
    DBUF(cb_sa,  (size_t)nrec*sizeof(float));
    DBUF(cb_ds,  (size_t)nrec*ncomp*sizeof(float));
    DBUF(cb_rs,  (size_t)nrec*nrho*sizeof(float));
    DBUF(cb_vm,  (size_t)nrec*sizeof(unsigned char));
    DBUF(cb_al,  (size_t)nrec*npts*sizeof(float2));
    DBUF(cb_alct,(size_t)nrec*npts*sizeof(float));
    DBUF(cb_cmbc,(size_t)NSEGC*ncomp*npts*sizeof(float2));
    DBUF(cb_cmbr,(size_t)NSEGC*nrho*npts*sizeof(float));
    DBUF(cb_sd,  (size_t)NSEGC*sizeof(int));
    float2 *dp_c = (float2*)cb_c.p, *dp_al = (float2*)cb_al.p, *dp_cmbc = (float2*)cb_cmbc.p;
    float *dp_ct = (float*)cb_ct.p, *dp_alct = (float*)cb_alct.p, *dp_cmbr = (float*)cb_cmbr.p,
        *dp_ca = (float*)cb_ca.p, *dp_sa = (float*)cb_sa.p, *dp_ds = (float*)cb_ds.p,
        *dp_rs = (float*)cb_rs.p;
    unsigned char *dp_vm = (unsigned char*)cb_vm.p; int *dp_sd = (int*)cb_sd.p;
    if (hstage(pin_c,  dp_c,  pl_c,  (size_t)nrec*plsz*sizeof(float2))) return 5;
    if (hstage(pin_ct, dp_ct, pl_ct, (size_t)nrec*plsz*sizeof(float)))  return 5;
    CKI(cudaMemcpy(dp_ca, cav,   (size_t)nrec*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dp_sa, sav,   (size_t)nrec*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dp_ds, dscale,(size_t)nrec*ncomp*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dp_rs, rscale,(size_t)nrec*nrho*sizeof(float),  cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dp_vm, vmask, (size_t)nrec*sizeof(unsigned char), cudaMemcpyHostToDevice));
    cudaEventRecord(ev1);
    {
        const size_t nwork = (size_t)nrec*npts;
        const int nblk = (int)((nwork + 255)/256);
        align_rot2d_kernel<<<nblk, 256>>>(dp_c, dp_ct, dp_ca, dp_sa, dp_vm, nrec,
            hpd_lo, hpd_hi, kpd_lo, hlo, hhi, klo, khi, pf, nyq2, dp_al, dp_alct);
        CKI(cudaGetLastError());
    }
    cudaEventRecord(ev2);
    const float one = 1.f, zero = 0.f;
    for (int s0 = 0; s0 < nseg; s0 += NSEGC) {
        const int nsc = (nseg - s0 < NSEGC) ? (nseg - s0) : NSEGC;
        for (int js = 0; js < nsc; ++js) {
            const int beg = seg_beg_h[s0+js], cnt = seg_cnt_h[s0+js];
            if (cublasSgemm(g_hb_cb, CUBLAS_OP_N, CUBLAS_OP_T, 2*(int)npts, ncomp, cnt,
                    &one, (const float*)dp_al + (size_t)beg*2*npts, 2*(int)npts,
                    dp_ds + (size_t)beg*ncomp, ncomp,
                    &zero, (float*)dp_cmbc + (size_t)js*ncomp*2*npts, 2*(int)npts)
                != CUBLAS_STATUS_SUCCESS) return 3;
            if (cublasSgemm(g_hb_cb, CUBLAS_OP_N, CUBLAS_OP_T, (int)npts, nrho, cnt,
                    &one, dp_alct + (size_t)beg*npts, (int)npts,
                    dp_rs + (size_t)beg*nrho, nrho,
                    &zero, dp_cmbr + (size_t)js*nrho*npts, (int)npts)
                != CUBLAS_STATUS_SUCCESS) return 3;
        }
        CKI(cudaMemcpy(dp_sd, seg_dir_h + s0, (size_t)nsc*sizeof(int), cudaMemcpyHostToDevice));
        cudaEventRecord(ev3);
        const size_t nwork = (size_t)nsc*npts;
        const int nblk = (int)((nwork + 255)/256);
        splat_bank_kernel<<<nblk, 256>>>(dp_cmbc, dp_cmbr, d_cb_rmat, dp_sd,
            nsc, ncomp, nrho, hlo, hhi, klo, khi, nyq2,
            lb[0], lb[1], lb[2], ub[0], ub[1], ub[2], g2_nvox, d_cmat2, d_rho2);
        CKI(cudaGetLastError());
    }
    cudaEventRecord(ev4);
    CKI(cudaDeviceSynchronize());
    {
        float ms;
        cudaEventElapsedTime(&ms, ev0, ev1); g_cbt_up += ms;
        cudaEventElapsedTime(&ms, ev1, ev2); g_cbt_al += ms;
        cudaEventElapsedTime(&ms, ev2, ev3); g_cbt_gm += ms;   // GEMMs of the LAST chunk only when multi-chunk
        cudaEventElapsedTime(&ms, ev3, ev4); g_cbt_sp += ms;
        cudaEventDestroy(ev0); cudaEventDestroy(ev1); cudaEventDestroy(ev2);
        cudaEventDestroy(ev3); cudaEventDestroy(ev4);
    }
    free(seg_dir_h); free(seg_beg_h); free(seg_cnt_h);
    return 0;
}

int flex_gpu_coupled_bank_free(void);   // forward (defined below)

}   // extern "C" (psample section follows)

// ---- Device particle prep (norm_noise_mask_pad_fft + gen_fplane4rec, observation model) ----
// Raw box images go up (16-64 KB each, ~4-16x less than the packed planes they replace); the
// noise-normalization, soft mask, pad+fftshift, batched cuFFT R2C and the packed plane build
// (shift phase + CTF transfer) all run on device, writing conj(T)y and ctfsq DIRECTLY into the
// resident cb_c/cb_ct buffers the fused E-step and banked M-step consume. CPU prep remains the
// reference and the fallback (sigma2 whitening, taper path, nsym>1, non-fused stages).
namespace {
float *dpr_raw = nullptr;      // [nrec][n1*n1] raw images
float *dpr_pad = nullptr;      // [nrec][n1o*n1o] padded/shifted real
cufftComplex *dpr_spec = nullptr; // [nrec][(n1o/2+1)*n1o]
float *dpr_cs2 = nullptr;      // [n1] (i-1-n1/2)^2
unsigned char *dpr_lmsk = nullptr; // [n1*n1] 1 = inside resolution mask
float *dpr_ctf = nullptr;      // [nrec][8]: sum_df, diff_df, angast, phshift, ampc, wl, hwl2cs, flip
float *dpr_psh = nullptr;      // [nrec][2]
float *dpr_stat = nullptr;     // [nrec][4]: mean, invstd, ave, edge_mean
cufftHandle gpr_plan = 0;
int gpr_n1 = 0, gpr_n1o = 0, gpr_maxrec = 0, gpr_lctf = 1;
float gpr_mskrad = 0.f, gpr_minlen = 0.f;
int gpr_wbg = 1;

__device__ __forceinline__ float pr_cosedge(float r2, float minlen, float mskrad)
{
    const float maxrad  = 0.5f*minlen;
    const float width   = 2.f*(maxrad - mskrad);
    const float r0      = maxrad - width;
    const float maxrad2 = maxrad*maxrad;
    const float r02     = r0*r0;
    if (r2 >= maxrad2) return 0.f;
    if (r2 <= r02)     return 1.f;
    const float t = (sqrtf(r2) - r0) / width;
    return 0.5f*(cosf(t*3.14159265358979323846f) + 1.f);
}

// TAPER variant (the live flex path: COV_MASK_IMAGES=.false.): per image, subtract ramped
// smoothed edge profiles in x then y (the y pass reads the x-tapered values, exactly as the
// CPU does), then border mean + noise stats over !lmsk on the TAPERED image. One block per
// image; profiles live in shared memory; the raw buffer is mutated in place.
__global__ void prep_taper_kernel(
    float* __restrict__ raw, const unsigned char* __restrict__ lmsk,
    int n1, int wtap, int wavg, int wbg,
    float* __restrict__ stat)
{
    const int i = blockIdx.x;
    float* im = raw + (size_t)i*n1*n1;
    extern __shared__ float shp[];       // 4*n1 floats: prof_a, prof_b, smooth_a, smooth_b
    float* pa = shp; float* pb = shp + n1; float* sa = shp + 2*n1; float* sb = shp + 3*n1;
    __shared__ double sh_a[256], sh_b[256];
    __shared__ int    sh_c[256];
    const int t = threadIdx.x, nt = blockDim.x;
    const float inv_wavg = 1.f/(float)wavg;
    const float inv_wtap = 1.f/(float)wtap;
    if (wtap > 0) {
        // ---- x edges: profiles over j (column index) ----
        for (int j = t; j < n1; j += nt) {
            float s0 = 0.f, s1 = 0.f;
            for (int a = 0; a < wavg; ++a) {
                s0 += im[(size_t)j*n1 + a];
                s1 += im[(size_t)j*n1 + (n1-1-a)];
            }
            s0 *= inv_wavg; s1 *= inv_wavg;
            const float av = 0.5f*(s0 + s1);
            pa[j] = s0 - av; pb[j] = s1 - av;
        }
        __syncthreads();
        for (int j = t; j < n1; j += nt) {
            const int j1 = (j > 0) ? j-1 : 0;
            const int j2 = (j < n1-1) ? j+1 : n1-1;
            const float inv = 1.f/(float)(j2 - j1 + 1);
            float s0 = 0.f, s1 = 0.f;
            for (int a = j1; a <= j2; ++a) { s0 += pa[a]; s1 += pb[a]; }
            sa[j] = s0*inv; sb[j] = s1*inv;
        }
        __syncthreads();
        for (int p = t; p < n1*wtap; p += nt) {
            const int ii = p % wtap;             // 0-based row into the ramp
            const int j  = p / wtap;
            const float al0 = (float)(wtap - ii)*inv_wtap;      // rows 1..wtap
            const float al1 = (float)(ii + 1)*inv_wtap;         // rows n-wtap+1..n
            im[(size_t)j*n1 + ii]          -= sa[j]*al0;
            im[(size_t)j*n1 + (n1-wtap+ii)] -= sb[j]*al1;
        }
        __syncthreads();
        // ---- y edges on the x-tapered image: profiles over i (row index) ----
        for (int ii = t; ii < n1; ii += nt) {
            float s0 = 0.f, s1 = 0.f;
            for (int a = 0; a < wavg; ++a) {
                s0 += im[(size_t)a*n1 + ii];
                s1 += im[(size_t)(n1-1-a)*n1 + ii];
            }
            s0 *= inv_wavg; s1 *= inv_wavg;
            const float av = 0.5f*(s0 + s1);
            pa[ii] = s0 - av; pb[ii] = s1 - av;
        }
        __syncthreads();
        for (int ii = t; ii < n1; ii += nt) {
            const int j1 = (ii > 0) ? ii-1 : 0;
            const int j2 = (ii < n1-1) ? ii+1 : n1-1;
            const float inv = 1.f/(float)(j2 - j1 + 1);
            float s0 = 0.f, s1 = 0.f;
            for (int a = j1; a <= j2; ++a) { s0 += pa[a]; s1 += pb[a]; }
            sa[ii] = s0*inv; sb[ii] = s1*inv;
        }
        __syncthreads();
        for (int p = t; p < n1*wtap; p += nt) {
            const int jj = p % wtap;
            const int ii = p / wtap;
            const float al0 = (float)(wtap - jj)*inv_wtap;
            const float al1 = (float)(jj + 1)*inv_wtap;
            im[(size_t)jj*n1 + ii]           -= sa[ii]*al0;
            im[(size_t)(n1-wtap+jj)*n1 + ii] -= sb[ii]*al1;
        }
        __syncthreads();
    }
    // ---- border (edge) mean over the tapered image ----
    double a = 0.;
    for (int p = t; p < n1*n1; p += nt) {
        const int px = p % n1, py = p / n1;
        const bool border = (px < wbg) || (px >= n1-wbg) ||
                            ((py < wbg || py >= n1-wbg) && px >= wbg && px < n1-wbg);
        if (border) a += (double)im[p];
    }
    sh_a[t] = a;
    __syncthreads();
    for (int sdim = nt/2; sdim > 0; sdim >>= 1) {
        if (t < sdim) sh_a[t] += sh_a[t+sdim];
        __syncthreads();
    }
    __shared__ float s_edgem;
    if (t == 0) {
        const int bc = 4*wbg*n1 - 4*wbg*wbg;
        s_edgem = (float)(sh_a[0]/(double)bc);
    }
    __syncthreads();
    // ---- noise stats over !lmsk on the tapered image ----
    a = 0.; double b = 0.; int cnt = 0;
    for (int p = t; p < n1*n1; p += nt) {
        if (!lmsk[p]) { const double v = im[p]; a += v; b += v*v; cnt++; }
    }
    sh_a[t] = a; sh_b[t] = b; sh_c[t] = cnt;
    __syncthreads();
    for (int sdim = nt/2; sdim > 0; sdim >>= 1) {
        if (t < sdim) { sh_a[t]+=sh_a[t+sdim]; sh_b[t]+=sh_b[t+sdim]; sh_c[t]+=sh_c[t+sdim]; }
        __syncthreads();
    }
    if (t == 0) {
        const int np = sh_c[0];
        float mean = 0.f, invstd = 1.f;
        if (np > 1) {
            const double m  = sh_a[0]/np;
            const double vr = (sh_b[0] - sh_a[0]*sh_a[0]/np) / (np - 1);
            if (vr > 0. && isfinite(vr)) { mean = (float)m; invstd = (float)(1./sqrt(vr)); }
        }
        stat[(size_t)i*4+0] = mean;
        stat[(size_t)i*4+1] = invstd;
        stat[(size_t)i*4+2] = 0.f;
        stat[(size_t)i*4+3] = s_edgem;
    }
}

__global__ void prep_fill_taper_kernel(
    const float* __restrict__ raw, const float* __restrict__ stat,
    int nrec, int n1, int n1o,
    float* __restrict__ pad)
{
    const size_t tid = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    const size_t npo = (size_t)n1o*n1o;
    if (tid >= (size_t)nrec*npo) return;
    const int i = (int)(tid / npo);
    const size_t po = tid - (size_t)i*npo;
    const int xo = (int)(po % n1o), yo = (int)(po / n1o);
    const float mean = stat[(size_t)i*4+0], invstd = stat[(size_t)i*4+1];
    const float edgem = stat[(size_t)i*4+3];
    const int h1o = n1o/2, start = (n1o - n1)/2;
    const int x0 = (xo + h1o) % n1o;
    const int y0 = (yo + h1o) % n1o;
    const int ii = x0 - start, jj = y0 - start;
    float v;
    const bool donorm = (mean != 0.f || invstd != 1.f);
    if (ii < 0 || ii >= n1 || jj < 0 || jj >= n1) {
        v = donorm ? (edgem - mean)*invstd : edgem;
    } else {
        v = raw[(size_t)i*n1*n1 + (size_t)jj*n1 + ii];
        if (donorm) v = (v - mean)*invstd;
    }
    pad[(size_t)i*npo + po] = v;
}

// one block per image: three shared reductions (noise stats over !lmsk, outside-mask mean of
// normalized values, border mean of normalized+masked values)
__global__ void prep_stats_kernel(
    const float* __restrict__ raw, const unsigned char* __restrict__ lmsk,
    const float* __restrict__ cs2, int n1, float mskrad, float minlen, int wbg,
    float* __restrict__ stat)
{
    const int i = blockIdx.x;
    const float* im = raw + (size_t)i*n1*n1;
    __shared__ double sh_a[256], sh_b[256];
    __shared__ float  s_mean, s_invstd, s_ave;
    const int t = threadIdx.x, nt = blockDim.x;
    const int npx = n1*n1;
    double a = 0., b = 0.;
    int cnt = 0;
    for (int p = t; p < npx; p += nt) {
        if (!lmsk[p]) { const double v = im[p]; a += v; b += v*v; cnt++; }
    }
    sh_a[t] = a; sh_b[t] = b;
    __shared__ int sh_c[256];
    sh_c[t] = cnt;
    __syncthreads();
    for (int sdim = nt/2; sdim > 0; sdim >>= 1) {
        if (t < sdim) { sh_a[t]+=sh_a[t+sdim]; sh_b[t]+=sh_b[t+sdim]; sh_c[t]+=sh_c[t+sdim]; }
        __syncthreads();
    }
    if (t == 0) {
        const int np = sh_c[0];
        float mean = 0.f, invstd = 1.f;
        if (np > 1) {
            const double m  = sh_a[0]/np;
            const double vr = (sh_b[0] - sh_a[0]*sh_a[0]/np) / (np - 1);
            if (vr > 0. && isfinite(vr)) { mean = (float)m; invstd = (float)(1./sqrt(vr)); }
        }
        s_mean = mean; s_invstd = invstd;
    }
    __syncthreads();
    const float mean = s_mean, invstd = s_invstd;
    const float rad2 = mskrad*mskrad;
    a = 0.; cnt = 0;
    for (int p = t; p < npx; p += nt) {
        const int px = p % n1, py = p / n1;
        if (cs2[px] + cs2[py] > rad2) { a += (double)((im[p]-mean)*invstd); cnt++; }
    }
    sh_a[t] = a; sh_c[t] = cnt;
    __syncthreads();
    for (int sdim = nt/2; sdim > 0; sdim >>= 1) {
        if (t < sdim) { sh_a[t]+=sh_a[t+sdim]; sh_c[t]+=sh_c[t+sdim]; }
        __syncthreads();
    }
    if (t == 0) s_ave = (sh_c[0] > 0) ? (float)(sh_a[0]/sh_c[0]) : 0.f;
    __syncthreads();
    const float ave = s_ave;
    a = 0.; cnt = 0;
    for (int p = t; p < npx; p += nt) {
        const int px = p % n1, py = p / n1;
        const bool border = (px < wbg) || (px >= n1-wbg) ||
                            ((py < wbg || py >= n1-wbg) && px >= wbg && px < n1-wbg);
        if (!border) continue;
        float v = (im[p]-mean)*invstd;
        const float e = pr_cosedge(cs2[px]+cs2[py], minlen, mskrad);
        if (e < 1.e-4f) v = ave;
        else if (e < 1.f-1.e-4f) v = e*v + (1.f-e)*ave;
        a += (double)v; cnt++;
    }
    sh_a[t] = a; sh_c[t] = cnt;
    __syncthreads();
    for (int sdim = nt/2; sdim > 0; sdim >>= 1) {
        if (t < sdim) { sh_a[t]+=sh_a[t+sdim]; sh_c[t]+=sh_c[t+sdim]; }
        __syncthreads();
    }
    if (t == 0) {
        stat[(size_t)i*4+0] = mean;
        stat[(size_t)i*4+1] = invstd;
        stat[(size_t)i*4+2] = ave;
        stat[(size_t)i*4+3] = (sh_c[0] > 0) ? (float)(sh_a[0]/sh_c[0]) : ave;
    }
}

__global__ void prep_fill_kernel(
    const float* __restrict__ raw, const float* __restrict__ cs2,
    const float* __restrict__ stat, int nrec, int n1, int n1o,
    float mskrad, float minlen,
    float* __restrict__ pad)
{
    const size_t tid = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    const size_t npo = (size_t)n1o*n1o;
    if (tid >= (size_t)nrec*npo) return;
    const int i = (int)(tid / npo);
    const size_t po = tid - (size_t)i*npo;
    const int xo = (int)(po % n1o), yo = (int)(po / n1o);
    const float mean = stat[(size_t)i*4+0], invstd = stat[(size_t)i*4+1];
    const float ave  = stat[(size_t)i*4+2], edgem  = stat[(size_t)i*4+3];
    // invert pad+fftshift: source (i0,j0) of this padded position, or background
    const int h1o = n1o/2, start = (n1o - n1)/2;   // 0-based start
    const int x0 = (xo + h1o) % n1o;               // undo fftshift (0-based)
    const int y0 = (yo + h1o) % n1o;
    const int ii = x0 - start, jj = y0 - start;    // 0-based input coords
    float v;
    if (ii < 0 || ii >= n1 || jj < 0 || jj >= n1) {
        v = edgem;
    } else {
        v = (raw[(size_t)i*n1*n1 + (size_t)jj*n1 + ii] - mean)*invstd;
        const float e = pr_cosedge(cs2[ii]+cs2[jj], minlen, mskrad);
        if (e < 1.e-4f) v = ave;
        else if (e < 1.f-1.e-4f) v = e*v + (1.f-e)*ave;
    }
    pad[(size_t)i*npo + po] = v;
}

// packed plane build straight into the resident cb_c/cb_ct: c = tval*y*phase (conj(T)y with T
// real), ctfsq = tval^2; spectrum read with the FFTW/cuFFT r2c half-layout + Friedel
__global__ void prep_plane_kernel(
    const cufftComplex* __restrict__ spec, const float* __restrict__ ctfp,
    const float* __restrict__ psh, const unsigned char* __restrict__ vmask,
    const float* __restrict__ sig2, int nsig,
    int nrec, int n1o, int hs_lo, int hs_hi, int ks_lo, int nyq_pd, int pf, int l_ctf,
    float2* __restrict__ plc, float* __restrict__ plct,
    float2* __restrict__ plcy, float* __restrict__ plt)
{
    const int nh = hs_hi - hs_lo + 1;
    const int nk = 1 - ks_lo;
    const size_t npts = (size_t)nh*nk;
    const size_t tid  = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (tid >= (size_t)nrec*npts) return;
    const int i = (int)(tid / npts);
    const size_t pt = tid - (size_t)i*npts;
    const int h = hs_lo + (int)(pt % nh);
    const int k = ks_lo + (int)(pt / nh);
    const size_t pidx = (size_t)i*npts + pt;
    if (!vmask[i]) {
        plc[pidx] = make_float2(0.f,0.f); plct[pidx] = 0.f;
        if (plcy) { plcy[pidx] = make_float2(0.f,0.f); plt[pidx] = 0.f; }
        return;
    }
    const int hp = h*pf, kp = k*pf;
    const int shell2 = hp*hp + kp*kp;
    if (shell2 > nyq_pd*(nyq_pd+1)) {
        plc[pidx] = make_float2(0.f,0.f); plct[pidx] = 0.f;
        if (plcy) { plcy[pidx] = make_float2(0.f,0.f); plt[pidx] = 0.f; }
        return;
    }
    // r2c physical address on the padded spectrum (x-fastest, nxc = n1o/2+1 columns)
    const int nxc = n1o/2 + 1;
    int ah = hp, ak = kp; bool cj = false;
    if (ah < 0) { ah = -ah; ak = -ak; cj = true; }
    if (ak < 0) ak += n1o;
    cufftComplex c = spec[(size_t)i*nxc*n1o + (size_t)ak*nxc + ah];
    if (cj) c.y = -c.y;
    const float scale = 1.f/((float)n1o*(float)n1o);
    c.x *= scale; c.y *= scale;
    // shift phase exp(i(hp*ps1 + kp*ps2)) -- pshift given in PADDED index units
    const float ph = hp*psh[(size_t)i*2+0] + kp*psh[(size_t)i*2+1];
    const float cph = cosf(ph), sph = sinf(ph);
    const float yre = c.x*cph - c.y*sph;
    const float yim = c.x*sph + c.y*cph;
    float tval = 1.f;
    if (l_ctf) {
        const float* cp = ctfp + (size_t)i*8;
        const float s2  = ((float)(hp*hp) + (float)(kp*kp))/((float)n1o*(float)n1o);
        const float cterm = cosf(2.f*(atan2f((float)kp,(float)hp) - cp[2]));
        const float df  = 0.5f*(cp[0] + cterm*cp[1]);
        const float pws = 3.14159265358979323846f*cp[5]*s2;
        tval = sinf(pws*(df - cp[6]*s2) + cp[3] + cp[4]);
        // flag: 0=no ctf, 1=normal (signed), 2=phase-flipped (|ctf|); abs applies to flip ONLY
        if (cp[7] > 1.5f) tval = fabsf(tval);
    }
    float w = 1.f;
    if (nsig > 0) {
        // observation model with sigma2 whitening: y and T each carry 1/sqrt(s), so the
        // premultiplied product and ctfsq both carry 1/s
        int sh = __float2int_rn(sqrtf((float)shell2));
        if (sh >= nsig) sh = nsig - 1;
        const float sv = sig2[(size_t)i*nsig + sh];
        if (sv > 1.e-20f) w = 1.f/sv;
    }
    plc[pidx]  = make_float2(w*tval*yre, w*tval*yim);
    plct[pidx] = w*tval*tval;
    if (plcy) {
        // separated observation-model pair: y and T each carry 1/sqrt(s)
        const float sw = sqrtf(w);
        plcy[pidx] = make_float2(sw*yre, sw*yim);
        plt[pidx]  = sw*tval;
    }
}
}   // anon namespace

extern "C" {

int flex_gpu_prep_begin(
    const void* lmsk, const void* cs2, int n1, int n1o, int maxrec,
    float mskrad, float minlen, int wbg, int l_ctf)
{
    gpr_n1 = n1; gpr_n1o = n1o; gpr_maxrec = maxrec; gpr_lctf = l_ctf;
    gpr_mskrad = mskrad; gpr_minlen = minlen; gpr_wbg = wbg;
    CKI(cudaMalloc(&dpr_raw,  (size_t)maxrec*n1*n1*sizeof(float)));
    CKI(cudaMalloc(&dpr_pad,  (size_t)maxrec*n1o*n1o*sizeof(float)));
    CKI(cudaMalloc(&dpr_spec, (size_t)maxrec*(n1o/2+1)*n1o*sizeof(cufftComplex)));
    CKI(cudaMalloc(&dpr_cs2,  (size_t)n1*sizeof(float)));
    CKI(cudaMalloc(&dpr_lmsk, (size_t)n1*n1*sizeof(unsigned char)));
    CKI(cudaMalloc(&dpr_ctf,  (size_t)maxrec*8*sizeof(float)));
    CKI(cudaMalloc(&dpr_psh,  (size_t)maxrec*2*sizeof(float)));
    CKI(cudaMalloc(&dpr_stat, (size_t)maxrec*4*sizeof(float)));
    CKI(cudaMemcpy(dpr_cs2,  cs2,  (size_t)n1*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dpr_lmsk, lmsk, (size_t)n1*n1*sizeof(unsigned char), cudaMemcpyHostToDevice));
    if (cufftPlan1d(&gpr_plan, 0, CUFFT_R2C, 1) == CUFFT_SUCCESS) cufftDestroy(gpr_plan);
    {
        int n[2] = {n1o, n1o};
        if (cufftPlanMany(&gpr_plan, 2, n, NULL, 1, n1o*n1o, NULL, 1, (n1o/2+1)*n1o,
                CUFFT_R2C, maxrec) != CUFFT_SUCCESS) { gpr_plan = 0; return 7; }
    }
    return 0;
}

// runs the whole chain and leaves the packed planes RESIDENT in cb_c/cb_ct with the estep
// geometry registered, so flex_gpu_estep_batch_resident and the banked M-step follow directly
int flex_gpu_prep_batch(
    const void* raw, const void* ctfp, const void* psh, const void* vmask,
    const void* sig2, int nsig,
    int nrec, int hs_lo, int hs_hi, int ks_lo, int nyq_pd)
{
    if (!dpr_raw || nrec > gpr_maxrec) return 2;
    const int n1 = gpr_n1, n1o = gpr_n1o;
    const int hpdn = hs_hi - hs_lo + 1;
    const size_t plsz = (size_t)(1 - ks_lo)*hpdn;
    DBUF(cb_c,  (size_t)nrec*plsz*sizeof(float2));
    DBUF(cb_ct, (size_t)nrec*plsz*sizeof(float));
    DBUF(cb_cy, (size_t)nrec*plsz*sizeof(float2));
    DBUF(cb_t,  (size_t)nrec*plsz*sizeof(float));
    DBUF(cb_vm, (size_t)nrec*sizeof(unsigned char));
    if (hstage(pin_c, dpr_raw, raw, (size_t)nrec*n1*n1*sizeof(float))) return 5;
    static DBuf dpr_sig2;
    if (nsig > 0) {
        if (dbuf_ensure(dpr_sig2, (size_t)nrec*nsig*sizeof(float))) return 4;
        CKI(cudaMemcpy(dpr_sig2.p, sig2, (size_t)nrec*nsig*sizeof(float), cudaMemcpyHostToDevice));
    }
    CKI(cudaMemcpy(dpr_ctf, ctfp, (size_t)nrec*8*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dpr_psh, psh,  (size_t)nrec*2*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(cb_vm.p, vmask,(size_t)nrec*sizeof(unsigned char), cudaMemcpyHostToDevice));
    if (gpr_mskrad > 0.f) {
        prep_stats_kernel<<<nrec, 256>>>(dpr_raw, dpr_lmsk, dpr_cs2, n1,
            gpr_mskrad, gpr_minlen, gpr_wbg, dpr_stat);
        CKI(cudaGetLastError());
        const size_t nwork = (size_t)nrec*n1o*n1o;
        prep_fill_kernel<<<(int)((nwork+255)/256), 256>>>(dpr_raw, dpr_cs2, dpr_stat,
            nrec, n1, n1o, gpr_mskrad, gpr_minlen, dpr_pad);
        CKI(cudaGetLastError());
    } else {
        // taper path: winsz = nint(COSMSKHALFWIDTH)=6 -> wtap=6 (n>=12), wavg=max(1,wtap/8)=1,
        // wbg=max(1,winsz/2)=3 -- passed in via gpr_wbg/gpr_minlen slots at begin
        const int wtap = (int)gpr_minlen;   // repurposed: taper width
        prep_taper_kernel<<<nrec, 256, 4*n1*sizeof(float)>>>(dpr_raw, dpr_lmsk,
            n1, wtap, 1, gpr_wbg, dpr_stat);
        CKI(cudaGetLastError());
        const size_t nwork = (size_t)nrec*n1o*n1o;
        prep_fill_taper_kernel<<<(int)((nwork+255)/256), 256>>>(dpr_raw, dpr_stat,
            nrec, n1, n1o, dpr_pad);
        CKI(cudaGetLastError());
    }
    if (cufftExecR2C(gpr_plan, dpr_pad, dpr_spec) != CUFFT_SUCCESS) return 7;
    {
        const size_t nwork = (size_t)nrec*plsz;
        prep_plane_kernel<<<(int)((nwork+255)/256), 256>>>(dpr_spec, dpr_ctf, dpr_psh,
            (const unsigned char*)cb_vm.p,
            (nsig > 0) ? (const float*)dpr_sig2.p : (const float*)dpr_ctf, nsig,
            nrec, n1o, hs_lo, hs_hi, ks_lo, nyq_pd,
            OSMPL_PAD_FAC_DEV, gpr_lctf, (float2*)cb_c.p, (float*)cb_ct.p,
            (float2*)cb_cy.p, (float*)cb_t.p);
        CKI(cudaGetLastError());
    }
    // register the packed geometry for resident consumers (fused E-step re-registers with
    // its own band; the resident insert relies on THIS in stages that never run an E-step).
    // Planes are square, so the +k packed limit equals the +h one.
    ge_hs_lo = hs_lo; ge_hs_hi = hs_hi; ge_ks_lo = ks_lo; ge_khi = hs_hi; ge_nrec = nrec;
    return 0;
}

// fetch the resident packed planes (validation cross-check only)
int flex_gpu_prep_fetch(void* plc_host, void* plct_host, int nrec, int hs_lo, int hs_hi, int ks_lo)
{
    const int hpdn = hs_hi - hs_lo + 1;
    const size_t plsz = (size_t)(1 - ks_lo)*hpdn;
    if (!cb_c.p || !cb_ct.p) return 2;
    CKI(cudaDeviceSynchronize());
    CKI(cudaMemcpy(plc_host,  cb_c.p,  (size_t)nrec*plsz*sizeof(float2), cudaMemcpyDeviceToHost));
    CKI(cudaMemcpy(plct_host, cb_ct.p, (size_t)nrec*plsz*sizeof(float),  cudaMemcpyDeviceToHost));
    return 0;
}

// fetch the separated observation-model planes (image, CTF, ctfsq) for host-side consumers:
// the general-stage path unpacks these into fplane_type
int flex_gpu_prep_fetch_sep(void* plcy_host, void* plt_host, void* plct_host,
    int nrec, int hs_lo, int hs_hi, int ks_lo)
{
    const int hpdn = hs_hi - hs_lo + 1;
    const size_t plsz = (size_t)(1 - ks_lo)*hpdn;
    if (!cb_cy.p || !cb_t.p || !cb_ct.p) return 2;
    CKI(cudaDeviceSynchronize());
    CKI(cudaMemcpy(plcy_host, cb_cy.p, (size_t)nrec*plsz*sizeof(float2), cudaMemcpyDeviceToHost));
    CKI(cudaMemcpy(plt_host,  cb_t.p,  (size_t)nrec*plsz*sizeof(float),  cudaMemcpyDeviceToHost));
    CKI(cudaMemcpy(plct_host, cb_ct.p, (size_t)nrec*plsz*sizeof(float),  cudaMemcpyDeviceToHost));
    return 0;
}

int flex_gpu_prep_free(void)
{
    if (dpr_raw)  { cudaFree(dpr_raw);  dpr_raw  = nullptr; }
    if (dpr_pad)  { cudaFree(dpr_pad);  dpr_pad  = nullptr; }
    if (dpr_spec) { cudaFree(dpr_spec); dpr_spec = nullptr; }
    if (dpr_cs2)  { cudaFree(dpr_cs2);  dpr_cs2  = nullptr; }
    if (dpr_lmsk) { cudaFree(dpr_lmsk); dpr_lmsk = nullptr; }
    if (dpr_ctf)  { cudaFree(dpr_ctf);  dpr_ctf  = nullptr; }
    if (dpr_psh)  { cudaFree(dpr_psh);  dpr_psh  = nullptr; }
    if (dpr_stat) { cudaFree(dpr_stat); dpr_stat = nullptr; }
    if (gpr_plan) { cufftDestroy(gpr_plan); gpr_plan = 0; }
    gpr_n1 = 0; gpr_n1o = 0; gpr_maxrec = 0;
    return 0;
}

}   // extern "C"

// ---- Fused probe E-step (exact per-particle directions, volumes resident) ----
// The banked-forward experiment established that the E-step LATENTS need exact per-particle
// projection directions. This keeps them: the mean + basis expanded lattices live on device
// (~43 MB), and per (record, plane point) the kernel gathers all ncomp+1 volumes at the
// record's own rotation (the Cartesian projector's exact convention: NONREDUNDANT half-volume
// with conjugation at loc_x < 0) and accumulates the five stat groups the posterior solve
// needs -- e_mm, myv, b_q, c_q, G_qr -- as DOUBLE atomics. The payload is the SAME packed
// (conj(T)y, ctfsq) batch the M-step uploads, so the marginal transfer is the 170 doubles per
// record coming back. After the host solves the 16x16 systems, a second kernel forms the
// premultiplied residual conj(T)y - a*ctfsq*proj0 on device and it is fetched packed.
namespace {
float2 *de_vols = nullptr;     // [ncomp+1][nvox_exp] expanded volume lattices
int     ge_ncomp1 = 0;
size_t  ge_nvox = 0;
int     ge_lb[3] = {0,0,0}, ge_ub[3] = {0,0,0};
DBuf    de_rm, de_proj0, de_stats, de_avec;

__global__ void estep_stats_kernel(
    const float2* __restrict__ pl_c, const float* __restrict__ pl_ct,
    const float*  __restrict__ rotm, const unsigned char* __restrict__ vmask,
    const float2* __restrict__ vols, int ncomp1,
    int nrec, int hs_lo, int hs_hi, int ks_lo, int khi, int nyq2,
    int lb1, int lb2, int lb3, int ub1, int ub2, int ub3, size_t nvox,
    float2* __restrict__ proj0, double* __restrict__ stats, int nst)
{
    const int nh = hs_hi - hs_lo + 1;
    const int nk = khi  - ks_lo + 1;
    const size_t npts = (size_t)nh*nk;
    const size_t tid  = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (tid >= (size_t)nrec*npts) return;
    const int i  = (int)(tid / npts);
    const size_t pt = tid - (size_t)i*npts;
    if (!vmask[i]) return;
    const int h = hs_lo + (int)(pt % nh);
    const int k = ks_lo + (int)(pt / nh);
    if (k > 0) return;                     // stats live on the stored k<=0 half
    if (h*h + k*k > nyq2) return;
    if (k == 0 && h > 0) return;           // cov_herm_inner: k=0 keeps only h<=0
    const int hpdn = nh;
    const size_t plsz = (size_t)(1 - ks_lo)*hpdn;
    const size_t pidx = (size_t)i*plsz + (size_t)(k - ks_lo)*hpdn + (h - hs_lo);
    const float2 plc = pl_c[pidx];
    const float  ct  = pl_ct[pidx];
    // exact projection location, the Cartesian projector's convention
    const float* m = rotm + (size_t)i*9;
    float lx = h*m[0] + k*m[1];
    float ly = h*m[3] + k*m[4];
    float lz = h*m[6] + k*m[7];
    bool cj = false;
    if (lx < 0.f) { lx = -lx; ly = -ly; lz = -lz; cj = true; }
    const int wx0 = __float2int_rn(lx) - IWINSZ;
    const int wy0 = __float2int_rn(ly) - IWINSZ;
    const int wz0 = __float2int_rn(lz) - IWINSZ;
    const bool inb = !(wx0 < lb1 || wy0 < lb2 || wz0 < lb3 ||
                       wx0+2*IWINSZ > ub1 || wy0+2*IWINSZ > ub2 || wz0+2*IWINSZ > ub3);
    // ncomp1 <= 64 guarded host-side. Indexed by a runtime loop variable, so these live in
    // LOCAL memory (not registers) regardless of size -- raising the bound costs small-rank
    // users nothing; it only extends the per-thread stack frame actually touched (512 B at 64).
    float u_re[64], u_im[64];
    if (inb) {
        float wx[WDIM], wy[WDIM], wz[WDIM];
        float sx = 0.f, sy = 0.f, sz = 0.f;
        #pragma unroll
        for (int j = 0; j < WDIM; ++j) {
            wx[j] = kb_apod_fast((float)(wx0+j) - lx); sx += wx[j];
            wy[j] = kb_apod_fast((float)(wy0+j) - ly); sy += wy[j];
            wz[j] = kb_apod_fast((float)(wz0+j) - lz); sz += wz[j];
        }
        const float eps = 1.1920929e-07f;
        if (fabsf(sx) > eps) { const float sc = 1.f/sx; wx[0]*=sc; wx[1]*=sc; wx[2]*=sc; }
        else                 { wx[0]=wx[1]=wx[2]=1.f/3.f; }
        if (fabsf(sy) > eps) { const float sc = 1.f/sy; wy[0]*=sc; wy[1]*=sc; wy[2]*=sc; }
        else                 { wy[0]=wy[1]=wy[2]=1.f/3.f; }
        if (fabsf(sz) > eps) { const float sc = 1.f/sz; wz[0]*=sc; wz[1]*=sc; wz[2]*=sc; }
        else                 { wz[0]=wz[1]=wz[2]=1.f/3.f; }
        const int d1 = ub1-lb1+1, d2 = ub2-lb2+1;
        for (int q = 0; q < ncomp1; ++q) { u_re[q] = 0.f; u_im[q] = 0.f; }
        #pragma unroll
        for (int iz = 0; iz < WDIM; ++iz) {
            const size_t zoff = (size_t)(wz0+iz-lb3)*d2*d1;
            #pragma unroll
            for (int iy = 0; iy < WDIM; ++iy) {
                const size_t yoff = zoff + (size_t)(wy0+iy-lb2)*d1;
                #pragma unroll
                for (int ix = 0; ix < WDIM; ++ix) {
                    const float ww = wx[ix]*(wy[iy]*wz[iz]);
                    const size_t vidx = yoff + (size_t)(wx0+ix-lb1);
                    for (int q = 0; q < ncomp1; ++q) {
                        const float2 v = vols[(size_t)q*nvox + vidx];
                        u_re[q] += ww*v.x; u_im[q] += ww*v.y;
                    }
                }
            }
        }
        if (cj) for (int q = 0; q < ncomp1; ++q) u_im[q] = -u_im[q];
    } else {
        for (int q = 0; q < ncomp1; ++q) { u_re[q] = 0.f; u_im[q] = 0.f; }
    }
    proj0[(size_t)i*npts + pt] = make_float2(u_re[0], u_im[0]);
    // five stat groups (REAL parts, as the CPU takes them); layout per record:
    // [0]=e_mm  [1]=myv  [2..nc+1]=b_q  [nc+2..2nc+1]=c_q  [2nc+2..]=G upper (q<=r)
    double* st = stats + (size_t)i*nst;
    const int nc = ncomp1 - 1;
    atomicAdd(&st[0], (double)ct*((double)u_re[0]*u_re[0] + (double)u_im[0]*u_im[0]));
    atomicAdd(&st[1], (double)u_re[0]*plc.x + (double)u_im[0]*plc.y);
    for (int q = 1; q <= nc; ++q) {
        atomicAdd(&st[1+q],    (double)u_re[q]*plc.x + (double)u_im[q]*plc.y);
        atomicAdd(&st[1+nc+q], (double)ct*((double)u_re[q]*u_re[0] + (double)u_im[q]*u_im[0]));
    }
    int off = 2 + 2*nc;
    for (int q = 1; q <= nc; ++q)
        for (int r = q; r <= nc; ++r, ++off)
            atomicAdd(&st[off], (double)ct*((double)u_re[q]*u_re[r] + (double)u_im[q]*u_im[r]));
}

__global__ void estep_resid_kernel(
    float2* __restrict__ pl_c, const float* __restrict__ pl_ct,
    const float2* __restrict__ proj0, const float* __restrict__ avec,
    const unsigned char* __restrict__ vmask,
    int nrec, int hs_lo, int hs_hi, int ks_lo, int khi)
{
    const int nh = hs_hi - hs_lo + 1;
    const int nk = khi  - ks_lo + 1;
    const size_t npts = (size_t)nh*nk;
    const size_t tid  = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (tid >= (size_t)nrec*npts) return;
    const int i  = (int)(tid / npts);
    const size_t pt = tid - (size_t)i*npts;
    if (!vmask[i]) return;
    const int k = ks_lo + (int)(pt / nh);
    if (k > 0) return;
    const int hpdn = nh;
    const size_t plsz = (size_t)(1 - ks_lo)*hpdn;
    const int h = hs_lo + (int)(pt % nh);
    const size_t pidx = (size_t)i*plsz + (size_t)(k - ks_lo)*hpdn + (h - hs_lo);
    const float a = avec[i];
    const float2 p0 = proj0[(size_t)i*npts + pt];
    pl_c[pidx].x -= a*pl_ct[pidx]*p0.x;
    pl_c[pidx].y -= a*pl_ct[pidx]*p0.y;
}
}   // anon namespace

extern "C" {

int flex_gpu_estep_vols(const void* vol, int q1, int ncomp1, const int* lb, const int* ub)
{
    if (q1 == 1) {
        if (ncomp1 > 64) return 6;   // must match u_re/u_im extents in estep kernels
        for (int j = 0; j < 3; ++j) { ge_lb[j] = lb[j]; ge_ub[j] = ub[j]; }
        ge_nvox = (size_t)(ub[0]-lb[0]+1)*(ub[1]-lb[1]+1)*(ub[2]-lb[2]+1);
        ge_ncomp1 = ncomp1;
        if (de_vols) { cudaFree(de_vols); de_vols = nullptr; }
        CKI(cudaMalloc(&de_vols, (size_t)ncomp1*ge_nvox*sizeof(float2)));
    }
    if (!de_vols || q1 < 1 || q1 > ge_ncomp1) return 2;
    CKI(cudaMemcpy(de_vols + (size_t)(q1-1)*ge_nvox, vol, ge_nvox*sizeof(float2),
        cudaMemcpyHostToDevice));
    return 0;
}

// consumes the SAME packed (conj(T)y, ctfsq) payload layout the coupled wrappers build;
// planes stay resident in the cb_* pools for the residual kernel and the M-step
int flex_gpu_estep_batch(
    const void* pl_c, const void* pl_ct, const void* rotm, const void* vmask,
    int nrec, int hs_lo, int hs_hi, int ks_lo, int khi, int nyq2,
    void* stats_host)
{
    if (!de_vols) return 2;
    const int hpdn = hs_hi - hs_lo + 1;
    const size_t plsz = (size_t)(1 - ks_lo)*hpdn;
    const size_t npts = (size_t)hpdn*(khi - ks_lo + 1);
    const int nc  = ge_ncomp1 - 1;
    const int nst = 2 + 2*nc + (nc*(nc+1))/2;
    DBUF(cb_c,  (size_t)nrec*plsz*sizeof(float2));
    DBUF(cb_ct, (size_t)nrec*plsz*sizeof(float));
    DBUF(cb_vm, (size_t)nrec*sizeof(unsigned char));
    DBUF(de_rm,    (size_t)nrec*9*sizeof(float));
    DBUF(de_proj0, (size_t)nrec*npts*sizeof(float2));
    DBUF(de_stats, (size_t)nrec*nst*sizeof(double));
    ge_hs_lo = hs_lo; ge_hs_hi = hs_hi; ge_ks_lo = ks_lo; ge_khi = khi; ge_nrec = nrec; ge_nyq2 = nyq2;
    if (hstage(pin_c,  cb_c.p,  pl_c,  (size_t)nrec*plsz*sizeof(float2))) return 5;
    if (hstage(pin_ct, cb_ct.p, pl_ct, (size_t)nrec*plsz*sizeof(float)))  return 5;
    CKI(cudaMemcpy(cb_vm.p, vmask, (size_t)nrec*sizeof(unsigned char), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(de_rm.p, rotm,  (size_t)nrec*9*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemset(de_stats.p, 0, (size_t)nrec*nst*sizeof(double)));
    const size_t nwork = (size_t)nrec*npts;
    estep_stats_kernel<<<(int)((nwork+255)/256), 256>>>(
        (const float2*)cb_c.p, (const float*)cb_ct.p, (const float*)de_rm.p,
        (const unsigned char*)cb_vm.p, de_vols, ge_ncomp1,
        nrec, hs_lo, hs_hi, ks_lo, khi, nyq2,
        ge_lb[0], ge_lb[1], ge_lb[2], ge_ub[0], ge_ub[1], ge_ub[2], ge_nvox,
        (float2*)de_proj0.p, (double*)de_stats.p, nst);
    CKI(cudaGetLastError());
    CKI(cudaMemcpy(stats_host, de_stats.p, (size_t)nrec*nst*sizeof(double),
        cudaMemcpyDeviceToHost));
    return 0;
}

// estep over planes ALREADY RESIDENT (device prep filled cb_c/cb_ct and cb_vm): only the
// rotations go up and the stat vectors come back
int flex_gpu_estep_batch_resident(
    const void* rotm, int nrec, int hs_lo, int hs_hi, int ks_lo, int khi, int nyq2,
    void* stats_host)
{
    if (!de_vols) return 2;
    const int hpdn = hs_hi - hs_lo + 1;
    const size_t npts = (size_t)hpdn*(khi - ks_lo + 1);
    const int nc  = ge_ncomp1 - 1;
    const int nst = 2 + 2*nc + (nc*(nc+1))/2;
    ge_hs_lo = hs_lo; ge_hs_hi = hs_hi; ge_ks_lo = ks_lo; ge_khi = khi; ge_nrec = nrec; ge_nyq2 = nyq2;
    DBUF(de_rm,    (size_t)nrec*9*sizeof(float));
    DBUF(de_proj0, (size_t)nrec*npts*sizeof(float2));
    DBUF(de_stats, (size_t)nrec*nst*sizeof(double));
    CKI(cudaMemcpy(de_rm.p, rotm, (size_t)nrec*9*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemset(de_stats.p, 0, (size_t)nrec*nst*sizeof(double)));
    const size_t nwork = (size_t)nrec*npts;
    estep_stats_kernel<<<(int)((nwork+255)/256), 256>>>(
        (const float2*)cb_c.p, (const float*)cb_ct.p, (const float*)de_rm.p,
        (const unsigned char*)cb_vm.p, de_vols, ge_ncomp1,
        nrec, hs_lo, hs_hi, ks_lo, khi, nyq2,
        ge_lb[0], ge_lb[1], ge_lb[2], ge_ub[0], ge_ub[1], ge_ub[2], ge_nvox,
        (float2*)de_proj0.p, (double*)de_stats.p, nst);
    CKI(cudaGetLastError());
    CKI(cudaMemcpy(stats_host, de_stats.p, (size_t)nrec*nst*sizeof(double),
        cudaMemcpyDeviceToHost));
    return 0;
}

// forms the premultiplied residual on device and fetches it PACKED; the host unpacks into
// the fplanes so the (unchanged) M-step path consumes it with a premultiplied flag
int flex_gpu_estep_resid(
    const void* avec, const void* vmask, int nrec,
    int hs_lo, int hs_hi, int ks_lo, int khi, void* resid_host)
{
    if (!de_vols) return 2;
    // sentinel bounds (lo > hi) = resident-only call: use the geometry registered by the
    // device prep / E-step upload. The host may have no built fplane to read bounds from
    // (the empty-frlims read here is what silently DISABLED mean deflation on the fused path).
    // khi stays 0: the residual lives on the stored k<=0 half-plane.
    if (hs_lo > hs_hi) {
        hs_lo = ge_hs_lo; hs_hi = ge_hs_hi; ks_lo = ge_ks_lo; khi = 0;
        if (hs_lo > hs_hi) return 3;
    }
    const int hpdn = hs_hi - hs_lo + 1;
    const size_t plsz = (size_t)(1 - ks_lo)*hpdn;
    const size_t npts = (size_t)hpdn*(khi - ks_lo + 1);
    DBUF(de_avec, (size_t)nrec*sizeof(float));
    CKI(cudaMemcpy(de_avec.p, avec, (size_t)nrec*sizeof(float), cudaMemcpyHostToDevice));
    const size_t nwork = (size_t)nrec*npts;
    estep_resid_kernel<<<(int)((nwork+255)/256), 256>>>(
        (float2*)cb_c.p, (const float*)cb_ct.p, (const float2*)de_proj0.p,
        (const float*)de_avec.p, (const unsigned char*)cb_vm.p,
        nrec, hs_lo, hs_hi, ks_lo, khi);
    CKI(cudaGetLastError());
    if (resid_host) {
        CKI(cudaMemcpy(resid_host, cb_c.p, (size_t)nrec*plsz*sizeof(float2),
            cudaMemcpyDeviceToHost));
    }
    return 0;
}

// resident banked M-step: consumes the residual the E-step left in cb_c/cb_ct -- NO plane
// re-upload (the double-upload made fused v1 a net loss; the four workers serialize on the
// link, so each avoided upload pays its full wall cost back)
int flex_gpu_coupled_batch_banked_res(
    const void* cav, const void* sav, const void* dscale, const void* rscale,
    const void* vmask, const int* dirid, int nrec, int ncomp, int nrho,
    const int* lb, const int* ub)
{
    if (!d_cmat2 || ncomp != g2_ncomp || nrho != g2_nrho) return 2;
    if (!d_cb_rmat || !g_hb_cb || ge_nrec < nrec) return 2;
    const int nyq2 = ge_nyq2;
    const int hs_lo = ge_hs_lo, hs_hi = ge_hs_hi, ks_lo = ge_ks_lo, khi = ge_khi;
    const int hpdn = hs_hi - hs_lo + 1;
    const int nh = hpdn, nk = khi - ks_lo + 1;
    const size_t npts = (size_t)nh*nk;
    const int NSEGC = nsegc_for(npts, ncomp, nrho);
    // TEMPORARY VRAM audit: report the geometry and every device buffer this call demands,
    // once per process, so an OOM can be attributed to a term instead of inferred from the
    // failing byte count. Remove once the footprint is bounded.
    {
        static bool audited = false;
        if (!audited) {
            audited = true;
            const double MB = 1.0/1048576.0;
            fprintf(stderr,
                "flex_gpu VRAM audit: h=[%d,%d] k=[%d,%d] nh=%d nk=%d npts=%zu | "
                "nrec=%d ncomp=%d nrho=%d NSEGC=%d\n",
                hs_lo, hs_hi, ks_lo, khi, nh, nk, npts, nrec, ncomp, nrho, NSEGC);
            fprintf(stderr,
                "flex_gpu VRAM audit: cb_al=%.1f cb_alct=%.1f cb_ds=%.1f cb_rs=%.1f "
                "cb_cmbc=%.1f cb_cmbr=%.1f TOTAL=%.1f MB\n",
                (double)nrec*npts*sizeof(float2)*MB,
                (double)nrec*npts*sizeof(float)*MB,
                (double)nrec*ncomp*sizeof(float)*MB,
                (double)nrec*nrho*sizeof(float)*MB,
                (double)NSEGC*ncomp*npts*sizeof(float2)*MB,
                (double)NSEGC*nrho*npts*sizeof(float)*MB,
                ((double)nrec*npts*(sizeof(float2)+sizeof(float))
                 + (double)nrec*((size_t)ncomp+nrho)*sizeof(float)
                 + (double)NSEGC*npts*((size_t)ncomp*sizeof(float2)+(size_t)nrho*sizeof(float)))*MB);
            fflush(stderr);
        }
    }
    int nseg = 0;
    int *seg_dir_h = (int*)malloc((size_t)nrec*sizeof(int));
    int *seg_beg_h = (int*)malloc((size_t)nrec*sizeof(int));
    int *seg_cnt_h = (int*)malloc((size_t)nrec*sizeof(int));
    for (int i = 0; i < nrec; ++i) {
        if (nseg > 0 && dirid[i] == seg_dir_h[nseg-1]) { seg_cnt_h[nseg-1]++; continue; }
        seg_dir_h[nseg] = dirid[i]; seg_beg_h[nseg] = i; seg_cnt_h[nseg] = 1; nseg++;
    }
    DBUF(cb_ca,  (size_t)nrec*sizeof(float));
    DBUF(cb_sa,  (size_t)nrec*sizeof(float));
    DBUF(cb_ds,  (size_t)nrec*ncomp*sizeof(float));
    DBUF(cb_rs,  (size_t)nrec*nrho*sizeof(float));
    DBUF(cb_al,  (size_t)nrec*npts*sizeof(float2));
    DBUF(cb_alct,(size_t)nrec*npts*sizeof(float));
    DBUF(cb_cmbc,(size_t)NSEGC*ncomp*npts*sizeof(float2));
    DBUF(cb_cmbr,(size_t)NSEGC*nrho*npts*sizeof(float));
    DBUF(cb_sd,  (size_t)NSEGC*sizeof(int));
    CKI(cudaMemcpy(cb_ca.p, cav,   (size_t)nrec*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(cb_sa.p, sav,   (size_t)nrec*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(cb_ds.p, dscale,(size_t)nrec*ncomp*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(cb_rs.p, rscale,(size_t)nrec*nrho*sizeof(float),  cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(cb_vm.p, vmask, (size_t)nrec*sizeof(unsigned char), cudaMemcpyHostToDevice));
    {
        const size_t nwork = (size_t)nrec*npts;
        align_rot2d_kernel<<<(int)((nwork+255)/256), 256>>>(
            (const float2*)cb_c.p, (const float*)cb_ct.p, (const float*)cb_ca.p,
            (const float*)cb_sa.p, (const unsigned char*)cb_vm.p, nrec,
            hs_lo, hs_hi, ks_lo, hs_lo, hs_hi, ks_lo, khi, 2, nyq2,
            (float2*)cb_al.p, (float*)cb_alct.p);
        CKI(cudaGetLastError());
    }
    const float one = 1.f, zero = 0.f;
    for (int s0 = 0; s0 < nseg; s0 += NSEGC) {
        const int nsc = (nseg - s0 < NSEGC) ? (nseg - s0) : NSEGC;
        for (int js = 0; js < nsc; ++js) {
            const int beg = seg_beg_h[s0+js], cnt = seg_cnt_h[s0+js];
            if (cublasSgemm(g_hb_cb, CUBLAS_OP_N, CUBLAS_OP_T, 2*(int)npts, ncomp, cnt,
                    &one, (const float*)cb_al.p + (size_t)beg*2*npts, 2*(int)npts,
                    (const float*)cb_ds.p + (size_t)beg*ncomp, ncomp,
                    &zero, (float*)cb_cmbc.p + (size_t)js*ncomp*2*npts, 2*(int)npts)
                != CUBLAS_STATUS_SUCCESS) return 3;
            if (cublasSgemm(g_hb_cb, CUBLAS_OP_N, CUBLAS_OP_T, (int)npts, nrho, cnt,
                    &one, (const float*)cb_alct.p + (size_t)beg*npts, (int)npts,
                    (const float*)cb_rs.p + (size_t)beg*nrho, nrho,
                    &zero, (float*)cb_cmbr.p + (size_t)js*nrho*npts, (int)npts)
                != CUBLAS_STATUS_SUCCESS) return 3;
        }
        CKI(cudaMemcpy(cb_sd.p, seg_dir_h + s0, (size_t)nsc*sizeof(int), cudaMemcpyHostToDevice));
        const size_t nwork = (size_t)nsc*npts;
        splat_bank_kernel<<<(int)((nwork+255)/256), 256>>>(
            (const float2*)cb_cmbc.p, (const float*)cb_cmbr.p, d_cb_rmat, (const int*)cb_sd.p,
            nsc, ncomp, nrho, hs_lo, hs_hi, ks_lo, khi, nyq2,
            lb[0], lb[1], lb[2], ub[0], ub[1], ub[2], g2_nvox, d_cmat2, d_rho2);
        CKI(cudaGetLastError());
    }
    CKI(cudaDeviceSynchronize());
    free(seg_dir_h); free(seg_beg_h); free(seg_cnt_h);
    return 0;
}

int flex_gpu_estep_free(void)
{
    if (de_vols) { cudaFree(de_vols); de_vols = nullptr; }
    ge_ncomp1 = 0; ge_nvox = 0;
    dbuf_release(de_rm); dbuf_release(de_proj0); dbuf_release(de_stats); dbuf_release(de_avec);
    return 0;
}

}   // extern "C"

// ---- Polar particle sampler (solve/embed PASS-A) ----
// Device port of polar_sample_particle_packed: per (record, band sample) rotate the ring angle
// by the record's in-plane (ca,sa), interpolate the PACKED y and transfer planes with the
// polar former's unit-tap 2D KB scheme (per-axis normalized apod, per-tap Friedel; parity
// sub-kernels renormalized by their weight sum exactly as polar_interp_plane does), emit
// xw = conj(T)y packed as sqwq*(re,im), ring |T|^2 sums for wr/tazim, and the noise-ring
// power for sig2. The apod is kb_apod_fast (the 14-term polynomial); the CPU sampler uses the
// exact apod -- the difference is ~1e-7 per tap and the statistical gates own it.
namespace {
float *dps_rad = nullptr, *dps_cs = nullptr, *dps_sn = nullptr, *dps_sqwq = nullptr;
int   *dps_ring = nullptr, *dps_nang = nullptr;
float *dps_nrad = nullptr, *dps_ncs = nullptr, *dps_nsn = nullptr, *dps_nwq = nullptr;
int    gps_nsamp = 0, gps_nk = 0, gps_nsampn = 0;
DBuf   ps_y, ps_t, ps_ca, ps_sa, ps_vm, ps_xw, ps_xw1, ps_xw2, ps_wrs, ps_wrq, ps_pw, ps_taz;
DBuf   ps_wro, ps_pwf, ps_tazf;   // float outputs (ring stats accumulate in double: the
                                  // Sum(t2^2)/n - mean^2 cancellation is too lossy in fp32)
PBuf   ps_pin_y, ps_pin_t;

__device__ __forceinline__ void ps_interp(
    const float2* __restrict__ pl, int hpd_lo, int hpd_hi, int kpd_lo, int hpdn,
    float hu, float ku,
    float2& full, float2& p1, float2& p2, float& ws1, float& ws2)
{
    const int wx0 = __float2int_rn(hu) - IWINSZ;
    const int wy0 = __float2int_rn(ku) - IWINSZ;
    float wx[WDIM], wy[WDIM];
    float sx = 0.f, sy = 0.f;
    #pragma unroll
    for (int j = 0; j < WDIM; ++j) {
        wx[j] = kb_apod_fast((float)(wx0+j) - hu); sx += wx[j];
        wy[j] = kb_apod_fast((float)(wy0+j) - ku); sy += wy[j];
    }
    const float eps = 1.1920929e-07f;
    if (fabsf(sx) > eps) { const float s = 1.f/sx; wx[0]*=s; wx[1]*=s; wx[2]*=s; }
    else                 { wx[0]=wx[1]=wx[2]=1.f/3.f; }
    if (fabsf(sy) > eps) { const float s = 1.f/sy; wy[0]*=s; wy[1]*=s; wy[2]*=s; }
    else                 { wy[0]=wy[1]=wy[2]=1.f/3.f; }
    full = make_float2(0.f,0.f); p1 = full; p2 = full; ws1 = 0.f; ws2 = 0.f;
    #pragma unroll
    for (int iy = 0; iy < WDIM; ++iy) {
        const int ky = wy0 + iy;
        #pragma unroll
        for (int ix = 0; ix < WDIM; ++ix) {
            const int hx = wx0 + ix;
            const float w = wx[ix]*wy[iy];
            float2 cv;
            if (ky <= 0) {
                if (hx < hpd_lo || hx > hpd_hi || ky < kpd_lo) continue;
                cv = pl[(size_t)(ky - kpd_lo)*hpdn + (hx - hpd_lo)];
            } else {
                if (-hx < hpd_lo || -hx > hpd_hi || -ky < kpd_lo) continue;
                const float2 c = pl[(size_t)(-ky - kpd_lo)*hpdn + (-hx - hpd_lo)];
                cv = make_float2(c.x, -c.y);
            }
            full.x += w*cv.x; full.y += w*cv.y;
            const int par = ((hx + ky) % 2 + 2) % 2 + 1;   // modulo(hx+ky,2)+1
            if (par == 1) { p1.x += w*cv.x; p1.y += w*cv.y; ws1 += w; }
            else          { p2.x += w*cv.x; p2.y += w*cv.y; ws2 += w; }
        }
    }
}

__global__ void psample_kernel(
    const float2* __restrict__ ply, const float2* __restrict__ plt,
    const float* __restrict__ cav, const float* __restrict__ sav,
    const unsigned char* __restrict__ vmask,
    const float* __restrict__ rad, const float* __restrict__ cs, const float* __restrict__ sn,
    const float* __restrict__ sqwq, const int* __restrict__ ring_of,
    int nrec, int nsamp, int want_halves,
    int hpd_lo, int hpd_hi, int kpd_lo, int nk,
    float* __restrict__ xw, float* __restrict__ xw1, float* __restrict__ xw2,
    double* __restrict__ wrs_sum, double* __restrict__ wrs_sq)
{
    const size_t tid = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (tid >= (size_t)nrec*nsamp) return;
    const int i = (int)(tid / nsamp);
    const int j = (int)(tid - (size_t)i*nsamp);
    if (!vmask[i]) return;
    const float ca = cav[i], sa = sav[i];
    const float c1 = cs[j]*ca - sn[j]*sa;
    const float s1 = sn[j]*ca + cs[j]*sa;
    const float hu =  rad[j]*c1;
    const float ku = -rad[j]*s1;
    const int hpdn = hpd_hi - hpd_lo + 1;
    const size_t plsz = (size_t)(1 - kpd_lo)*hpdn;
    float2 yf, y1, y2, tf, t1, t2c; float wsy1, wsy2, wst1, wst2;
    ps_interp(ply + (size_t)i*plsz, hpd_lo, hpd_hi, kpd_lo, hpdn, hu, ku, yf, y1, y2, wsy1, wsy2);
    ps_interp(plt + (size_t)i*plsz, hpd_lo, hpd_hi, kpd_lo, hpdn, hu, ku, tf, t1, t2c, wst1, wst2);
    // xw = conj(T) * y, packed sqwq*(re, im)
    const float sq = sqwq[j];
    const float xr =  tf.x*yf.x + tf.y*yf.y;
    const float xi =  tf.x*yf.y - tf.y*yf.x;
    xw[(size_t)i*2*nsamp + 2*j    ] = sq*xr;
    xw[(size_t)i*2*nsamp + 2*j + 1] = sq*xi;
    if (want_halves) {
        const float e12 = 1.e-12f;
        float2 a = y1; if (fabsf(wsy1) > e12) { a.x /= wsy1; a.y /= wsy1; }
        float2 b = y2; if (fabsf(wsy2) > e12) { b.x /= wsy2; b.y /= wsy2; }
        xw1[(size_t)i*2*nsamp + 2*j    ] = sq*( tf.x*a.x + tf.y*a.y);
        xw1[(size_t)i*2*nsamp + 2*j + 1] = sq*( tf.x*a.y - tf.y*a.x);
        xw2[(size_t)i*2*nsamp + 2*j    ] = sq*( tf.x*b.x + tf.y*b.y);
        xw2[(size_t)i*2*nsamp + 2*j + 1] = sq*( tf.x*b.y - tf.y*b.x);
    }
    const double t2 = (double)tf.x*tf.x + (double)tf.y*tf.y;
    atomicAdd(&wrs_sum[(size_t)i*nk + ring_of[j]], t2);
    atomicAdd(&wrs_sq [(size_t)i*nk + ring_of[j]], t2*t2);
}

__global__ void psample_noise_kernel(
    const float2* __restrict__ ply, const unsigned char* __restrict__ vmask,
    const float* __restrict__ nrad, const float* __restrict__ ncs, const float* __restrict__ nsn,
    const float* __restrict__ nwq,
    const float* __restrict__ cav, const float* __restrict__ sav,
    int nrec, int nsampn, int hpd_lo, int hpd_hi, int kpd_lo,
    double* __restrict__ pw)
{
    const size_t tid = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (tid >= (size_t)nrec*nsampn) return;
    const int i = (int)(tid / nsampn);
    const int j = (int)(tid - (size_t)i*nsampn);
    if (!vmask[i]) return;
    const float ca = cav[i], sa = sav[i];
    const float c1 = ncs[j]*ca - nsn[j]*sa;
    const float s1 = nsn[j]*ca + ncs[j]*sa;
    const float hu =  nrad[j]*c1;
    const float ku = -nrad[j]*s1;
    const int hpdn = hpd_hi - hpd_lo + 1;
    const size_t plsz = (size_t)(1 - kpd_lo)*hpdn;
    float2 yf, p1, p2; float w1, w2;
    ps_interp(ply + (size_t)i*plsz, hpd_lo, hpd_hi, kpd_lo, hpdn, hu, ku, yf, p1, p2, w1, w2);
    atomicAdd(&pw[i], (double)nwq[j]*((double)yf.x*yf.x + (double)yf.y*yf.y));
}

__global__ void psample_finalize_kernel(
    const double* __restrict__ wrs_sum, const double* __restrict__ wrs_sq,
    const int* __restrict__ nang, const unsigned char* __restrict__ vmask,
    int nrec, int nk,
    float* __restrict__ wr, double* __restrict__ taz)
{
    const size_t tid = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (tid >= (size_t)nrec*nk) return;
    const int i  = (int)(tid / nk);
    const int ir = (int)(tid - (size_t)i*nk);
    if (!vmask[i]) return;
    const double na = (double)nang[ir];
    const double tm = wrs_sum[(size_t)i*nk + ir] / na;
    double tv = wrs_sq[(size_t)i*nk + ir] / na - tm*tm;
    if (tv < 0.) tv = 0.;
    wr[(size_t)i*nk + ir] = (float)tm;
    if (tm > 1.e-8) atomicAdd(&taz[i], sqrt(tv)/tm);
}

__global__ void psample_tofloat_kernel(
    const double* __restrict__ pw, const double* __restrict__ taz, int nrec,
    float* __restrict__ pwf, float* __restrict__ tazf)
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= nrec) return;
    pwf[i]  = (float)pw[i];
    tazf[i] = (float)taz[i];
}

// resident psample bridge: the device prep stores the transfer plane REAL (cb_t); the
// sampler kernels read it as float2
__global__ void real_to_c_kernel(const float* __restrict__ src, float2* __restrict__ dst,
    size_t n)
{
    const size_t i = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= n) return;
    dst[i] = make_float2(src[i], 0.f);
}

// Cartesian-lattice noise diagnostics over the resident whitened plane (cb_cy): per-shell
// power/count plus the high-frequency band sums the solve debias consumes. Accumulated on
// device across the whole stage, fetched once. Layout of acc:
// [0..nyq]=pwsh, [nyq+1..2nyq+1]=cntsh, [2nyq+2]=hfpw, [2nyq+3]=hfcnt
__global__ void psample_diag_kernel(
    const float2* __restrict__ ply, const unsigned char* __restrict__ vmask,
    int nrec, int hs_lo, int hs_hi, int ks_lo, int nyq, int sh_lo,
    double* __restrict__ acc)
{
    const int nh = hs_hi - hs_lo + 1;
    const int nk = 1 - ks_lo;   // stored k<=0 half, as the CPU diagnostics iterate
    const size_t npts = (size_t)nh*nk;
    const size_t tid  = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (tid >= (size_t)nrec*npts) return;
    const int ib = (int)(tid / npts);
    if (!vmask[ib]) return;
    const size_t pt = tid - (size_t)ib*npts;
    const int h = hs_lo + (int)(pt % nh);
    const int k = ks_lo + (int)(pt / nh);
    const int sh = (int)nearbyintf(sqrtf((float)(h*h + k*k)));
    if (sh > nyq) return;
    const float2 v = ply[(size_t)ib*npts + pt];
    const double p = (double)v.x*v.x + (double)v.y*v.y;
    atomicAdd(&acc[sh], p);
    atomicAdd(&acc[nyq+1+sh], 1.0);
    if (sh >= sh_lo) {
        atomicAdd(&acc[2*nyq+2], p);
        atomicAdd(&acc[2*nyq+3], 1.0);
    }
}
}   // anon namespace

extern "C" {

int flex_gpu_psample_begin(
    const void* rad, const void* cs, const void* sn, const void* sqwq, const int* ring_of,
    const int* nang, int nsamp, int nk,
    const void* nrad, const void* ncs, const void* nsn, const void* nwq, int nsampn)
{
    gps_nsamp = nsamp; gps_nk = nk; gps_nsampn = nsampn;
    CKI(cudaMalloc(&dps_rad,  (size_t)nsamp*sizeof(float)));
    CKI(cudaMalloc(&dps_cs,   (size_t)nsamp*sizeof(float)));
    CKI(cudaMalloc(&dps_sn,   (size_t)nsamp*sizeof(float)));
    CKI(cudaMalloc(&dps_sqwq, (size_t)nsamp*sizeof(float)));
    CKI(cudaMalloc(&dps_ring, (size_t)nsamp*sizeof(int)));
    CKI(cudaMalloc(&dps_nang, (size_t)nk*sizeof(int)));
    CKI(cudaMemcpy(dps_rad,  rad,  (size_t)nsamp*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dps_cs,   cs,   (size_t)nsamp*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dps_sn,   sn,   (size_t)nsamp*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dps_sqwq, sqwq, (size_t)nsamp*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dps_ring, ring_of, (size_t)nsamp*sizeof(int), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dps_nang, nang, (size_t)nk*sizeof(int), cudaMemcpyHostToDevice));
    if (nsampn > 0) {
        CKI(cudaMalloc(&dps_nrad, (size_t)nsampn*sizeof(float)));
        CKI(cudaMalloc(&dps_ncs,  (size_t)nsampn*sizeof(float)));
        CKI(cudaMalloc(&dps_nsn,  (size_t)nsampn*sizeof(float)));
        CKI(cudaMalloc(&dps_nwq,  (size_t)nsampn*sizeof(float)));
        CKI(cudaMemcpy(dps_nrad, nrad, (size_t)nsampn*sizeof(float), cudaMemcpyHostToDevice));
        CKI(cudaMemcpy(dps_ncs,  ncs,  (size_t)nsampn*sizeof(float), cudaMemcpyHostToDevice));
        CKI(cudaMemcpy(dps_nsn,  nsn,  (size_t)nsampn*sizeof(float), cudaMemcpyHostToDevice));
        CKI(cudaMemcpy(dps_nwq,  nwq,  (size_t)nsampn*sizeof(float), cudaMemcpyHostToDevice));
    }
    return 0;
}

int flex_gpu_psample_batch(
    const void* ply, const void* plt, const void* cav, const void* sav, const void* vmask,
    int nrec, int want_halves, int hpd_lo, int hpd_hi, int kpd_lo,
    void* xw_host, void* xw1_host, void* xw2_host, void* wr_host, void* pw_host, void* taz_host)
{
    if (gps_nsamp <= 0) return 2;
    const int hpdn = hpd_hi - hpd_lo + 1;
    const size_t plsz = (size_t)(1 - kpd_lo)*hpdn;
    const int ns2 = 2*gps_nsamp;
    DBUF(ps_y,  (size_t)nrec*plsz*sizeof(float2));
    DBUF(ps_t,  (size_t)nrec*plsz*sizeof(float2));
    DBUF(ps_ca, (size_t)nrec*sizeof(float));
    DBUF(ps_sa, (size_t)nrec*sizeof(float));
    DBUF(ps_vm, (size_t)nrec*sizeof(unsigned char));
    DBUF(ps_xw, (size_t)nrec*ns2*sizeof(float));
    if (want_halves) {
        DBUF(ps_xw1, (size_t)nrec*ns2*sizeof(float));
        DBUF(ps_xw2, (size_t)nrec*ns2*sizeof(float));
    }
    DBUF(ps_wrs, (size_t)nrec*gps_nk*sizeof(double));
    DBUF(ps_wrq, (size_t)nrec*gps_nk*sizeof(double));
    DBUF(ps_pw,  (size_t)nrec*sizeof(double));
    DBUF(ps_taz, (size_t)nrec*sizeof(double));
    DBUF(ps_wro, (size_t)nrec*gps_nk*sizeof(float));
    DBUF(ps_pwf, (size_t)nrec*sizeof(float));
    DBUF(ps_tazf,(size_t)nrec*sizeof(float));
    if (hstage(ps_pin_y, ps_y.p, ply, (size_t)nrec*plsz*sizeof(float2))) return 5;
    if (hstage(ps_pin_t, ps_t.p, plt, (size_t)nrec*plsz*sizeof(float2))) return 5;
    CKI(cudaMemcpy(ps_ca.p, cav, (size_t)nrec*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(ps_sa.p, sav, (size_t)nrec*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(ps_vm.p, vmask, (size_t)nrec*sizeof(unsigned char), cudaMemcpyHostToDevice));
    CKI(cudaMemset(ps_wrs.p, 0, (size_t)nrec*gps_nk*sizeof(double)));
    CKI(cudaMemset(ps_wrq.p, 0, (size_t)nrec*gps_nk*sizeof(double)));
    CKI(cudaMemset(ps_pw.p,  0, (size_t)nrec*sizeof(double)));
    CKI(cudaMemset(ps_taz.p, 0, (size_t)nrec*sizeof(double)));
    CKI(cudaMemset(ps_wro.p, 0, (size_t)nrec*gps_nk*sizeof(float)));
    CKI(cudaMemset(ps_xw.p,  0, (size_t)nrec*ns2*sizeof(float)));
    if (want_halves) {
        CKI(cudaMemset(ps_xw1.p, 0, (size_t)nrec*ns2*sizeof(float)));
        CKI(cudaMemset(ps_xw2.p, 0, (size_t)nrec*ns2*sizeof(float)));
    }
    {
        const size_t nwork = (size_t)nrec*gps_nsamp;
        psample_kernel<<<(int)((nwork+255)/256), 256>>>(
            (const float2*)ps_y.p, (const float2*)ps_t.p, (const float*)ps_ca.p,
            (const float*)ps_sa.p, (const unsigned char*)ps_vm.p,
            dps_rad, dps_cs, dps_sn, dps_sqwq, dps_ring,
            nrec, gps_nsamp, want_halves, hpd_lo, hpd_hi, kpd_lo, gps_nk,
            (float*)ps_xw.p, want_halves ? (float*)ps_xw1.p : (float*)ps_xw.p,
            want_halves ? (float*)ps_xw2.p : (float*)ps_xw.p,
            (double*)ps_wrs.p, (double*)ps_wrq.p);
        CKI(cudaGetLastError());
    }
    if (gps_nsampn > 0) {
        const size_t nwork = (size_t)nrec*gps_nsampn;
        psample_noise_kernel<<<(int)((nwork+255)/256), 256>>>(
            (const float2*)ps_y.p, (const unsigned char*)ps_vm.p,
            dps_nrad, dps_ncs, dps_nsn, dps_nwq,
            (const float*)ps_ca.p, (const float*)ps_sa.p,
            nrec, gps_nsampn, hpd_lo, hpd_hi, kpd_lo, (double*)ps_pw.p);
        CKI(cudaGetLastError());
    }
    {
        const size_t nwork = (size_t)nrec*gps_nk;
        psample_finalize_kernel<<<(int)((nwork+255)/256), 256>>>(
            (const double*)ps_wrs.p, (const double*)ps_wrq.p, dps_nang,
            (const unsigned char*)ps_vm.p, nrec, gps_nk,
            (float*)ps_wro.p, (double*)ps_taz.p);
        CKI(cudaGetLastError());
        psample_tofloat_kernel<<<(nrec+255)/256, 256>>>(
            (const double*)ps_pw.p, (const double*)ps_taz.p, nrec,
            (float*)ps_pwf.p, (float*)ps_tazf.p);
        CKI(cudaGetLastError());
    }
    CKI(cudaMemcpy(xw_host, ps_xw.p, (size_t)nrec*ns2*sizeof(float), cudaMemcpyDeviceToHost));
    if (want_halves) {
        CKI(cudaMemcpy(xw1_host, ps_xw1.p, (size_t)nrec*ns2*sizeof(float), cudaMemcpyDeviceToHost));
        CKI(cudaMemcpy(xw2_host, ps_xw2.p, (size_t)nrec*ns2*sizeof(float), cudaMemcpyDeviceToHost));
    }
    CKI(cudaMemcpy(wr_host,  ps_wro.p, (size_t)nrec*gps_nk*sizeof(float), cudaMemcpyDeviceToHost));
    CKI(cudaMemcpy(pw_host,  ps_pwf.p, (size_t)nrec*sizeof(float), cudaMemcpyDeviceToHost));
    CKI(cudaMemcpy(taz_host, ps_tazf.p,(size_t)nrec*sizeof(float), cudaMemcpyDeviceToHost));
    return 0;
}

// stage-persistent diagnostics accumulator for the resident sampler path
namespace { DBuf ps_diag; int gps_diag_nyq = -1; }

int flex_gpu_psample_diag_begin(int nyq)
{
    const size_t n = (size_t)(2*(nyq+1) + 2);
    if (dbuf_ensure(ps_diag, n*sizeof(double))) return 4;
    CKI(cudaMemset(ps_diag.p, 0, n*sizeof(double)));
    gps_diag_nyq = nyq;
    return 0;
}

// runs on the CURRENT resident whitened planes with the CURRENT ps_vm mask: call right
// after flex_gpu_psample_batch_resident for the same batch
int flex_gpu_psample_diag_batch(int nrec, int sh_lo)
{
    if (gps_diag_nyq < 0 || !cb_cy.p || !ps_vm.p) return 2;
    const int nh = ge_hs_hi - ge_hs_lo + 1;
    const int nk = 1 - ge_ks_lo;
    const size_t nwork = (size_t)nrec*nh*nk;
    psample_diag_kernel<<<(int)((nwork+255)/256), 256>>>(
        (const float2*)cb_cy.p, (const unsigned char*)ps_vm.p,
        nrec, ge_hs_lo, ge_hs_hi, ge_ks_lo, gps_diag_nyq, sh_lo,
        (double*)ps_diag.p);
    CKI(cudaGetLastError());
    return 0;
}

int flex_gpu_psample_diag_fetch(void* acc_host)
{
    if (gps_diag_nyq < 0) return 2;
    const size_t n = (size_t)(2*(gps_diag_nyq+1) + 2);
    CKI(cudaDeviceSynchronize());
    CKI(cudaMemcpy(acc_host, ps_diag.p, n*sizeof(double), cudaMemcpyDeviceToHost));
    gps_diag_nyq = -1;
    return 0;
}

// resident variant: sample straight off the separated device-prep pair (cb_cy = whitened
// image plane, cb_t = real transfer plane) with the registered packed geometry — the
// standalone port was transport-bound and measured NEGATIVE; with the planes already on
// device the same kernels never touch the bus.
int flex_gpu_psample_batch_resident(
    const void* cav, const void* sav, const void* vmask,
    int nrec, int want_halves,
    void* xw_host, void* xw1_host, void* xw2_host, void* wr_host, void* pw_host, void* taz_host)
{
    if (gps_nsamp <= 0) return 2;
    if (!cb_cy.p || !cb_t.p || ge_nrec < nrec) return 3;
    const int hpd_lo = ge_hs_lo, hpd_hi = ge_hs_hi, kpd_lo = ge_ks_lo;
    const int hpdn = hpd_hi - hpd_lo + 1;
    const size_t plsz = (size_t)(1 - kpd_lo)*hpdn;
    const int ns2 = 2*gps_nsamp;
    DBUF(ps_t,  (size_t)nrec*plsz*sizeof(float2));
    DBUF(ps_ca, (size_t)nrec*sizeof(float));
    DBUF(ps_sa, (size_t)nrec*sizeof(float));
    DBUF(ps_vm, (size_t)nrec*sizeof(unsigned char));
    DBUF(ps_xw, (size_t)nrec*ns2*sizeof(float));
    if (want_halves) {
        DBUF(ps_xw1, (size_t)nrec*ns2*sizeof(float));
        DBUF(ps_xw2, (size_t)nrec*ns2*sizeof(float));
    }
    DBUF(ps_wrs, (size_t)nrec*gps_nk*sizeof(double));
    DBUF(ps_wrq, (size_t)nrec*gps_nk*sizeof(double));
    DBUF(ps_pw,  (size_t)nrec*sizeof(double));
    DBUF(ps_taz, (size_t)nrec*sizeof(double));
    DBUF(ps_wro, (size_t)nrec*gps_nk*sizeof(float));
    DBUF(ps_pwf, (size_t)nrec*sizeof(float));
    DBUF(ps_tazf,(size_t)nrec*sizeof(float));
    {
        const size_t n = (size_t)nrec*plsz;
        real_to_c_kernel<<<(int)((n+255)/256), 256>>>((const float*)cb_t.p, (float2*)ps_t.p, n);
        CKI(cudaGetLastError());
    }
    CKI(cudaMemcpy(ps_ca.p, cav, (size_t)nrec*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(ps_sa.p, sav, (size_t)nrec*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(ps_vm.p, vmask, (size_t)nrec*sizeof(unsigned char), cudaMemcpyHostToDevice));
    CKI(cudaMemset(ps_wrs.p, 0, (size_t)nrec*gps_nk*sizeof(double)));
    CKI(cudaMemset(ps_wrq.p, 0, (size_t)nrec*gps_nk*sizeof(double)));
    CKI(cudaMemset(ps_pw.p,  0, (size_t)nrec*sizeof(double)));
    CKI(cudaMemset(ps_taz.p, 0, (size_t)nrec*sizeof(double)));
    CKI(cudaMemset(ps_wro.p, 0, (size_t)nrec*gps_nk*sizeof(float)));
    CKI(cudaMemset(ps_xw.p,  0, (size_t)nrec*ns2*sizeof(float)));
    if (want_halves) {
        CKI(cudaMemset(ps_xw1.p, 0, (size_t)nrec*ns2*sizeof(float)));
        CKI(cudaMemset(ps_xw2.p, 0, (size_t)nrec*ns2*sizeof(float)));
    }
    {
        const size_t nwork = (size_t)nrec*gps_nsamp;
        psample_kernel<<<(int)((nwork+255)/256), 256>>>(
            (const float2*)cb_cy.p, (const float2*)ps_t.p, (const float*)ps_ca.p,
            (const float*)ps_sa.p, (const unsigned char*)ps_vm.p,
            dps_rad, dps_cs, dps_sn, dps_sqwq, dps_ring,
            nrec, gps_nsamp, want_halves, hpd_lo, hpd_hi, kpd_lo, gps_nk,
            (float*)ps_xw.p, want_halves ? (float*)ps_xw1.p : (float*)ps_xw.p,
            want_halves ? (float*)ps_xw2.p : (float*)ps_xw.p,
            (double*)ps_wrs.p, (double*)ps_wrq.p);
        CKI(cudaGetLastError());
    }
    if (gps_nsampn > 0) {
        const size_t nwork = (size_t)nrec*gps_nsampn;
        psample_noise_kernel<<<(int)((nwork+255)/256), 256>>>(
            (const float2*)cb_cy.p, (const unsigned char*)ps_vm.p,
            dps_nrad, dps_ncs, dps_nsn, dps_nwq,
            (const float*)ps_ca.p, (const float*)ps_sa.p,
            nrec, gps_nsampn, hpd_lo, hpd_hi, kpd_lo, (double*)ps_pw.p);
        CKI(cudaGetLastError());
    }
    {
        const size_t nwork = (size_t)nrec*gps_nk;
        psample_finalize_kernel<<<(int)((nwork+255)/256), 256>>>(
            (const double*)ps_wrs.p, (const double*)ps_wrq.p, dps_nang,
            (const unsigned char*)ps_vm.p, nrec, gps_nk,
            (float*)ps_wro.p, (double*)ps_taz.p);
        CKI(cudaGetLastError());
        psample_tofloat_kernel<<<(nrec+255)/256, 256>>>(
            (const double*)ps_pw.p, (const double*)ps_taz.p, nrec,
            (float*)ps_pwf.p, (float*)ps_tazf.p);
        CKI(cudaGetLastError());
    }
    CKI(cudaMemcpy(xw_host, ps_xw.p, (size_t)nrec*ns2*sizeof(float), cudaMemcpyDeviceToHost));
    if (want_halves) {
        CKI(cudaMemcpy(xw1_host, ps_xw1.p, (size_t)nrec*ns2*sizeof(float), cudaMemcpyDeviceToHost));
        CKI(cudaMemcpy(xw2_host, ps_xw2.p, (size_t)nrec*ns2*sizeof(float), cudaMemcpyDeviceToHost));
    }
    CKI(cudaMemcpy(wr_host,  ps_wro.p, (size_t)nrec*gps_nk*sizeof(float), cudaMemcpyDeviceToHost));
    CKI(cudaMemcpy(pw_host,  ps_pwf.p, (size_t)nrec*sizeof(float), cudaMemcpyDeviceToHost));
    CKI(cudaMemcpy(taz_host, ps_tazf.p,(size_t)nrec*sizeof(float), cudaMemcpyDeviceToHost));
    return 0;
}

int flex_gpu_psample_free(void)
{
    if (dps_rad)  { cudaFree(dps_rad);  dps_rad  = nullptr; }
    if (dps_cs)   { cudaFree(dps_cs);   dps_cs   = nullptr; }
    if (dps_sn)   { cudaFree(dps_sn);   dps_sn   = nullptr; }
    if (dps_sqwq) { cudaFree(dps_sqwq); dps_sqwq = nullptr; }
    if (dps_ring) { cudaFree(dps_ring); dps_ring = nullptr; }
    if (dps_nang) { cudaFree(dps_nang); dps_nang = nullptr; }
    if (dps_nrad) { cudaFree(dps_nrad); dps_nrad = nullptr; }
    if (dps_ncs)  { cudaFree(dps_ncs);  dps_ncs  = nullptr; }
    if (dps_nsn)  { cudaFree(dps_nsn);  dps_nsn  = nullptr; }
    if (dps_nwq)  { cudaFree(dps_nwq);  dps_nwq  = nullptr; }
    dbuf_release(ps_y); dbuf_release(ps_t); dbuf_release(ps_ca); dbuf_release(ps_sa);
    dbuf_release(ps_vm); dbuf_release(ps_xw); dbuf_release(ps_xw1); dbuf_release(ps_xw2);
    dbuf_release(ps_wrs); dbuf_release(ps_wrq); dbuf_release(ps_pw); dbuf_release(ps_taz);
    dbuf_release(ps_wro); dbuf_release(ps_pwf); dbuf_release(ps_tazf);
    gps_nsamp = 0; gps_nk = 0; gps_nsampn = 0;
    return 0;
}

int flex_gpu_coupled_bank_free(void)
{
    if (g_cb_nbatch > 0) {
        printf(">>> FLEX_GPU BANKED PROFILE (stream seconds): upload=%7.2f align=%7.2f gemm=%7.2f splat=%7.2f  batches=%ld  nseg mean=%.1f max=%d\n",
            g_cbt_up*1e-3, g_cbt_al*1e-3, g_cbt_gm*1e-3, g_cbt_sp*1e-3,
            g_cb_nbatch, (double)g_cb_nseg_tot/(double)g_cb_nbatch, g_cb_nseg_max);
        fflush(stdout);
    }
    g_cbt_up = 0.; g_cbt_al = 0.; g_cbt_gm = 0.; g_cbt_sp = 0.;
    g_cb_nbatch = 0; g_cb_nseg_tot = 0; g_cb_nseg_max = 0;
    if (d_cb_rmat) { cudaFree(d_cb_rmat); d_cb_rmat = nullptr; }
    g_cb_ndir = 0;
    dbuf_release(cb_c); dbuf_release(cb_ct); dbuf_release(cb_ca); dbuf_release(cb_sa);
    dbuf_release(cb_ds); dbuf_release(cb_rs); dbuf_release(cb_vm); dbuf_release(cb_al);
    dbuf_release(cb_alct); dbuf_release(cb_cmbc); dbuf_release(cb_cmbr); dbuf_release(cb_sd);
    return 0;
}

}   // extern "C"

extern "C" {

// ---- SNR variance (slab design) ----
// var(x) += (sum of a particle's inserted taps at x)^2 is nonlinear per particle, so each
// particle needs its own completed splat before squaring. Chunk the batch: NSLAB per-particle
// slabs (cmplx+rho lattices) live in device scratch; kernel 1 scatters a chunk's particles into
// their slabs (same Friedel/window/apod-fast arithmetic as the multi kernel, no scales), kernel 2
// squares each slab cell into the var/dens accumulators and zeroes the slab for the next chunk.
namespace {
float2 *ds_slabC = nullptr; float *ds_slabR = nullptr;
float  *ds_var = nullptr, *ds_dens = nullptr;
size_t  gs_nvox = 0; int gs_nslab = 0;

__global__ void snr_insert_kernel(
    const float2* __restrict__ pl_c, const float* __restrict__ pl_ct,
    const float* __restrict__ rotm, const int* __restrict__ nyqdisk,
    const unsigned char* __restrict__ vmask,
    int nrec, int nsym, int hpd_lo, int hpd_hi, int kpd_lo,
    int hlo, int hhi, int klo, int khi, int pf,
    int lb1, int lb2, int lb3, int ub1, int ub2, int ub3,
    float tiny, size_t nvox,
    float2* __restrict__ slabC, float* __restrict__ slabR)
{
    const int nh = hhi - hlo + 1, nk = khi - klo + 1;
    const size_t tid = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (tid >= (size_t)nrec*nh*nk) return;
    const int i = (int)(tid / ((size_t)nh*nk));
    if (!vmask[i]) return;
    const size_t r2 = tid - (size_t)i*nh*nk;
    const int h = hlo + (int)(r2 % nh);
    const int k = klo + (int)(r2 / nh);
    if (h*h + k*k > nyqdisk[i]) return;
    // packed unpadded read (see insert_multi_kernel)
    const int hpdn = hpd_hi - hpd_lo + 1;
    const size_t plsz = (size_t)(1 - kpd_lo)*hpdn;
    float2 raw; float ct;
    if (k <= 0) {
        const size_t pidx = (size_t)i*plsz + (size_t)(k - kpd_lo)*hpdn + (h - hpd_lo);
        raw = pl_c[pidx]; ct = pl_ct[pidx];
    } else {
        const size_t pidx = (size_t)i*plsz + (size_t)(-k - kpd_lo)*hpdn + (-h - hpd_lo);
        const float2 c = pl_c[pidx];
        raw = make_float2(c.x, -c.y); ct = pl_ct[pidx];
    }
    if (fabsf(raw.x) + fabsf(raw.y) <= tiny && ct <= tiny) return;
    const float pf2 = (float)(pf*pf);
    const float cbre = pf2*raw.x, cbim = pf2*raw.y;
    const int d1 = ub1-lb1+1, d2 = ub2-lb2+1;
    for (int isym = 0; isym < nsym; ++isym) {
        const float* m = rotm + ((size_t)i*nsym + isym)*9;
        const float lx = h*m[0] + k*m[1];
        const float ly = h*m[3] + k*m[4];
        const float lz = h*m[6] + k*m[7];
        const int wx0 = __float2int_rn(lx) - IWINSZ;
        const int wy0 = __float2int_rn(ly) - IWINSZ;
        const int wz0 = __float2int_rn(lz) - IWINSZ;
        if (wx0 + 2*IWINSZ < 0) continue;
        if (wx0 < lb1 || wy0 < lb2 || wz0 < lb3) continue;
        if (wx0+2*IWINSZ > ub1 || wy0+2*IWINSZ > ub2 || wz0+2*IWINSZ > ub3) continue;
        float wx[WDIM], wy[WDIM], wz[WDIM];
        float sx = 0.f, sy = 0.f, sz = 0.f;
        #pragma unroll
        for (int j = 0; j < WDIM; ++j) {
            wx[j] = kb_apod_fast((float)(wx0+j) - lx); sx += wx[j];
            wy[j] = kb_apod_fast((float)(wy0+j) - ly); sy += wy[j];
            wz[j] = kb_apod_fast((float)(wz0+j) - lz); sz += wz[j];
        }
        const float eps = 1.1920929e-07f;
        if (fabsf(sx) > eps) { const float s = 1.f/sx; wx[0]*=s; wx[1]*=s; wx[2]*=s; }
        else                 { wx[0]=wx[1]=wx[2]=1.f/3.f; }
        if (fabsf(sy) > eps) { const float s = 1.f/sy; wy[0]*=s; wy[1]*=s; wy[2]*=s; }
        else                 { wy[0]=wy[1]=wy[2]=1.f/3.f; }
        if (fabsf(sz) > eps) { const float s = 1.f/sz; wz[0]*=s; wz[1]*=s; wz[2]*=s; }
        else                 { wz[0]=wz[1]=wz[2]=1.f/3.f; }
        float2* sc = slabC + (size_t)i*nvox;
        float*  sr = slabR + (size_t)i*nvox;
        #pragma unroll
        for (int iz = 0; iz < WDIM; ++iz) {
            const size_t zoff = (size_t)(wz0+iz-lb3)*d2*d1;
            #pragma unroll
            for (int iy = 0; iy < WDIM; ++iy) {
                const size_t yoff = zoff + (size_t)(wy0+iy-lb2)*d1;
                #pragma unroll
                for (int ix = 0; ix < WDIM; ++ix) {
                    const float w = wx[ix]*(wy[iy]*wz[iz]);
                    const size_t vidx = yoff + (size_t)(wx0+ix-lb1);
                    atomicAdd(&sc[vidx].x, cbre*w);
                    atomicAdd(&sc[vidx].y, cbim*w);
                    atomicAdd(&sr[vidx],   ct*w);
                }
            }
        }
    }
}

__global__ void snr_sweep_kernel(
    const unsigned char* __restrict__ vmask, int nrec, size_t nvox, float inv_pf2,
    float2* __restrict__ slabC, float* __restrict__ slabR,
    float* __restrict__ var, float* __restrict__ dens)
{
    const size_t tid = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (tid >= (size_t)nrec*nvox) return;
    const int    i = (int)(tid / nvox);
    const size_t v = tid % nvox;
    if (!vmask[i]) return;
    float2 c = slabC[(size_t)i*nvox + v];
    float  r = slabR[(size_t)i*nvox + v];
    if (c.x != 0.f || c.y != 0.f || r != 0.f) {
        const float re = c.x*inv_pf2, im = c.y*inv_pf2;
        atomicAdd(&var[v],  re*re + im*im);
        atomicAdd(&dens[v], r);
        slabC[(size_t)i*nvox + v] = make_float2(0.f, 0.f);
        slabR[(size_t)i*nvox + v] = 0.f;
    }
}
} // namespace

// ---- columns phase B ----
static float2 *dc_Be = nullptr, *dc_Bo = nullptr;
static float  *dc_He = nullptr, *dc_Ho = nullptr;
static size_t gc_nvox = 0;
static int    gc_ncol = 0;

int flex_gpu_cols_begin(int ncol, long nvox)
{
    if (dc_Be) { cudaFree(dc_Be); dc_Be = nullptr; }
    if (dc_Bo) { cudaFree(dc_Bo); dc_Bo = nullptr; }
    if (dc_He) { cudaFree(dc_He); dc_He = nullptr; }
    if (dc_Ho) { cudaFree(dc_Ho); dc_Ho = nullptr; }
    gc_nvox = (size_t)nvox; gc_ncol = ncol;
    CKI(cudaMalloc(&dc_Be, (size_t)ncol*gc_nvox*sizeof(float2)));
    CKI(cudaMalloc(&dc_Bo, (size_t)ncol*gc_nvox*sizeof(float2)));
    CKI(cudaMalloc(&dc_He, (size_t)ncol*gc_nvox*sizeof(float)));
    CKI(cudaMalloc(&dc_Ho, (size_t)ncol*gc_nvox*sizeof(float)));
    CKI(cudaMemset(dc_Be, 0, (size_t)ncol*gc_nvox*sizeof(float2)));
    CKI(cudaMemset(dc_Bo, 0, (size_t)ncol*gc_nvox*sizeof(float2)));
    CKI(cudaMemset(dc_He, 0, (size_t)ncol*gc_nvox*sizeof(float)));
    CKI(cudaMemset(dc_Ho, 0, (size_t)ncol*gc_nvox*sizeof(float)));
    return 0;
}

int flex_gpu_cols_batch(
    const void* cwin, const void* cwx, const void* cwy, const void* cwz,
    const void* cpl, const void* cct, const void* ncache, const void* vmask,
    const void* gcol, const void* hcolv, const void* qhkl,
    int nb, int maxc, int ncol,
    const int* lb, const int* ub, int deadskip, int debias, float pf2)
{
    if (!dc_Be || ncol != gc_ncol) return 2;
    const size_t nmc = (size_t)nb*maxc;
    DBUF(co_cw, nmc*3*sizeof(int));
    DBUF(co_wx, nmc*3*sizeof(float));
    DBUF(co_wy, nmc*3*sizeof(float));
    DBUF(co_wz, nmc*3*sizeof(float));
    DBUF(co_pl, nmc*sizeof(float2));
    DBUF(co_ct, nmc*sizeof(float));
    DBUF(co_nc, (size_t)nb*sizeof(int));
    DBUF(co_vm, (size_t)nb*sizeof(unsigned char));
    DBUF(co_gc, (size_t)nb*ncol*sizeof(float2));
    DBUF(co_hv, (size_t)nb*ncol*sizeof(float));
    DBUF(co_qh, (size_t)ncol*3*sizeof(int));
    int *d_cw = (int*)co_cw.p, *d_nc = (int*)co_nc.p, *d_qh = (int*)co_qh.p;
    float *d_wx = (float*)co_wx.p, *d_wy = (float*)co_wy.p, *d_wz = (float*)co_wz.p,
        *d_ct = (float*)co_ct.p, *d_hv = (float*)co_hv.p;
    float2 *d_pl = (float2*)co_pl.p, *d_gc = (float2*)co_gc.p;
    unsigned char *d_vm = (unsigned char*)co_vm.p;
    // pageable on purpose: the driver's chunked staging pipelines memcpy+DMA; a full
    // memcpy-into-pinned THEN DMA serializes the two and measured SLOWER (gpu13: 11.0->13.4 s)
    CKI(cudaMemcpy(d_cw, cwin,  nmc*3*sizeof(int),   cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(d_wx, cwx,   nmc*3*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(d_wy, cwy,   nmc*3*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(d_wz, cwz,   nmc*3*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(d_pl, cpl,   nmc*sizeof(float2),  cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(d_ct, cct,   nmc*sizeof(float),   cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(d_nc, ncache, nb*sizeof(int),     cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(d_vm, vmask, nb*sizeof(unsigned char), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(d_gc, gcol,  (size_t)nb*ncol*sizeof(float2), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(d_hv, hcolv, (size_t)nb*ncol*sizeof(float),  cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(d_qh, qhkl,  (size_t)ncol*3*sizeof(int),     cudaMemcpyHostToDevice));
    const int nblk = (int)((nmc + 127)/128);
    cols_scatter_kernel<<<nblk, 128>>>(d_cw, d_wx, d_wy, d_wz, d_pl, d_ct, d_nc, d_vm,
        d_gc, d_hv, d_qh, nb, maxc, ncol, lb[0], lb[1], lb[2], ub[0], ub[1], ub[2],
        deadskip, debias, pf2, gc_nvox, dc_Be, dc_He, dc_Bo, dc_Ho);
    CKI(cudaGetLastError());
    CKI(cudaDeviceSynchronize());
    return 0;
}

int flex_gpu_cols_fetch(int which, void* B_host, void* H_host)
{
    if (!dc_Be) return 2;
    const float2* B = (which == 0) ? dc_Be : dc_Bo;
    const float*  H = (which == 0) ? dc_He : dc_Ho;
    CKI(cudaMemcpy(B_host, B, (size_t)gc_ncol*gc_nvox*sizeof(float2), cudaMemcpyDeviceToHost));
    CKI(cudaMemcpy(H_host, H, (size_t)gc_ncol*gc_nvox*sizeof(float),  cudaMemcpyDeviceToHost));
    return 0;
}

// fused columns: upload the column lookup volume + column coordinates once per stage
int flex_gpu_cols_lookup(const void* lookup, const void* qhkl, int ncol,
    const int* lb, const int* ub)
{
    const size_t nvx = (size_t)(ub[0]-lb[0]+1)*(size_t)(ub[1]-lb[1]+1)*(size_t)(ub[2]-lb[2]+1);
    if (dc_lookup) { cudaFree(dc_lookup); dc_lookup = nullptr; }
    CKI(cudaMalloc(&dc_lookup, nvx*sizeof(int)));
    CKI(cudaMemcpy(dc_lookup, lookup, nvx*sizeof(int), cudaMemcpyHostToDevice));
    DBUF(co_qh, (size_t)ncol*3*sizeof(int));
    CKI(cudaMemcpy(co_qh.p, qhkl, (size_t)ncol*3*sizeof(int), cudaMemcpyHostToDevice));
    return 0;
}

// fused columns batch: cache build + gather + scatter, all on the resident adjoint pair
// (cb_c/cb_ct after estep_resid, geometry as registered). Zero cache traffic over PCIe.
int flex_gpu_cols_fused_batch(
    const void* rotm, const void* vmask, int nb, int maxc,
    int nyq_disk, int iwinsz, const int* lb, const int* ub,
    int deadskip, int debias, float pf2, float tiny)
{
    if (!dc_Be || !dc_lookup) return 2;
    if (!cb_c.p || !cb_ct.p || ge_nrec < nb) return 3;
    const size_t nmc = (size_t)nb*maxc;
    DBUF(co_cw, nmc*3*sizeof(int));
    DBUF(co_wx, nmc*3*sizeof(float));
    DBUF(co_wy, nmc*3*sizeof(float));
    DBUF(co_wz, nmc*3*sizeof(float));
    DBUF(co_pl, nmc*sizeof(float2));
    DBUF(co_ct, nmc*sizeof(float));
    DBUF(co_nc, (size_t)nb*sizeof(int));
    DBUF(co_vm, (size_t)nb*sizeof(unsigned char));
    DBUF(co_gc, (size_t)nb*gc_ncol*sizeof(float2));
    DBUF(co_hv, (size_t)nb*gc_ncol*sizeof(float));
    DBUF(co_rm, (size_t)nb*9*sizeof(float));
    DBUF(co_ovfl, sizeof(int));
    CKI(cudaMemset(co_nc.p, 0, (size_t)nb*sizeof(int)));
    CKI(cudaMemset(co_gc.p, 0, (size_t)nb*gc_ncol*sizeof(float2)));
    CKI(cudaMemset(co_hv.p, 0, (size_t)nb*gc_ncol*sizeof(float)));
    CKI(cudaMemset(co_ovfl.p, 0, sizeof(int)));
    CKI(cudaMemcpy(co_rm.p, rotm, (size_t)nb*9*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(co_vm.p, vmask, (size_t)nb*sizeof(unsigned char), cudaMemcpyHostToDevice));
    const int nh = ge_hs_hi - ge_hs_lo + 1;
    const int nk = ge_khi  - ge_ks_lo + 1;
    {
        const size_t nwork = (size_t)nb*nh*nk;
        cols_cache_kernel<<<(int)((nwork+255)/256), 256>>>(
            (const float2*)cb_c.p, (const float*)cb_ct.p, (const float*)co_rm.p,
            (const unsigned char*)co_vm.p, nb, ge_hs_lo, ge_hs_hi, ge_ks_lo, ge_khi,
            nyq_disk, iwinsz, lb[0], lb[1], lb[2], ub[0], ub[1], ub[2],
            tiny, maxc, (int*)co_cw.p, (float*)co_wx.p, (float*)co_wy.p, (float*)co_wz.p,
            (float2*)co_pl.p, (float*)co_ct.p, (int*)co_nc.p, (int*)co_ovfl.p);
        CKI(cudaGetLastError());
    }
    {
        const size_t nwork = nmc;
        cols_gather_kernel<<<(int)((nwork+127)/128), 128>>>(
            (const int*)co_cw.p, (const float*)co_wx.p, (const float*)co_wy.p,
            (const float*)co_wz.p, (const float2*)co_pl.p, (const float*)co_ct.p,
            (const int*)co_nc.p, (const unsigned char*)co_vm.p, dc_lookup,
            nb, maxc, gc_ncol, lb[0], lb[1], lb[2], ub[0], ub[1], ub[2],
            (float2*)co_gc.p, (float*)co_hv.p);
        CKI(cudaGetLastError());
    }
    {
        const int nblk = (int)((nmc + 127)/128);
        cols_scatter_kernel<<<nblk, 128>>>((const int*)co_cw.p, (const float*)co_wx.p,
            (const float*)co_wy.p, (const float*)co_wz.p, (const float2*)co_pl.p,
            (const float*)co_ct.p, (const int*)co_nc.p, (const unsigned char*)co_vm.p,
            (const float2*)co_gc.p, (const float*)co_hv.p, (const int*)co_qh.p,
            nb, maxc, gc_ncol, lb[0], lb[1], lb[2], ub[0], ub[1], ub[2],
            deadskip, debias, pf2, gc_nvox, dc_Be, dc_He, dc_Bo, dc_Ho);
        CKI(cudaGetLastError());
    }
    CKI(cudaDeviceSynchronize());
    {
        int ovfl = 0;
        CKI(cudaMemcpy(&ovfl, co_ovfl.p, sizeof(int), cudaMemcpyDeviceToHost));
        if (ovfl) return 5;
    }
    return 0;
}

// unit-test helper: inject synthetic packed adjoint planes as the resident pair and
// register their geometry, so the fused columns path can run without the full prep chain
int flex_gpu_cols_test_upload(const void* plc, const void* plct,
    int nrec, int hs_lo, int hs_hi, int ks_lo, int khi)
{
    const int hpdn = hs_hi - hs_lo + 1;
    const size_t plsz = (size_t)(1 - ks_lo)*hpdn;
    DBUF(cb_c,  (size_t)nrec*plsz*sizeof(float2));
    DBUF(cb_ct, (size_t)nrec*plsz*sizeof(float));
    CKI(cudaMemcpy(cb_c.p,  plc,  (size_t)nrec*plsz*sizeof(float2), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(cb_ct.p, plct, (size_t)nrec*plsz*sizeof(float),  cudaMemcpyHostToDevice));
    ge_hs_lo = hs_lo; ge_hs_hi = hs_hi; ge_ks_lo = ks_lo; ge_khi = khi; ge_nrec = nrec;
    return 0;
}

int flex_gpu_cols_free(void)
{
    if (dc_Be) { cudaFree(dc_Be); dc_Be = nullptr; }
    if (dc_Bo) { cudaFree(dc_Bo); dc_Bo = nullptr; }
    if (dc_He) { cudaFree(dc_He); dc_He = nullptr; }
    if (dc_Ho) { cudaFree(dc_Ho); dc_Ho = nullptr; }
    if (dc_lookup) { cudaFree(dc_lookup); dc_lookup = nullptr; }
    gc_nvox = 0; gc_ncol = 0;
    dbuf_release(co_cw); dbuf_release(co_wx); dbuf_release(co_wy); dbuf_release(co_wz);
    dbuf_release(co_pl); dbuf_release(co_ct); dbuf_release(co_nc); dbuf_release(co_vm);
    dbuf_release(co_gc); dbuf_release(co_hv); dbuf_release(co_qh);
    dbuf_release(co_rm); dbuf_release(co_ovfl);
    return 0;
}

int flex_gpu_snr_begin(long nvox, int nslab)
{
    if (ds_slabC) { cudaFree(ds_slabC); ds_slabC = nullptr; }
    if (ds_slabR) { cudaFree(ds_slabR); ds_slabR = nullptr; }
    if (ds_var)   { cudaFree(ds_var);   ds_var   = nullptr; }
    if (ds_dens)  { cudaFree(ds_dens);  ds_dens  = nullptr; }
    gs_nvox = (size_t)nvox; gs_nslab = nslab;
    CKI(cudaMalloc(&ds_slabC, (size_t)nslab*gs_nvox*sizeof(float2)));
    CKI(cudaMalloc(&ds_slabR, (size_t)nslab*gs_nvox*sizeof(float)));
    CKI(cudaMalloc(&ds_var,   gs_nvox*sizeof(float)));
    CKI(cudaMalloc(&ds_dens,  gs_nvox*sizeof(float)));
    CKI(cudaMemset(ds_slabC, 0, (size_t)nslab*gs_nvox*sizeof(float2)));
    CKI(cudaMemset(ds_slabR, 0, (size_t)nslab*gs_nvox*sizeof(float)));
    CKI(cudaMemset(ds_var,   0, gs_nvox*sizeof(float)));
    CKI(cudaMemset(ds_dens,  0, gs_nvox*sizeof(float)));
    return 0;
}

int flex_gpu_snr_batch(
    const void* pl_c, const void* pl_ct, const void* rotm, const void* nyqdisk,
    const void* vmask, int nrec, int nsym,
    int hpd_lo, int hpd_hi, int kpd_lo,
    int hlo, int hhi, int klo, int khi, int pf,
    const int* lb, const int* ub, float tiny, float inv_pf2)
{
    if (!ds_slabC) return 2;
    const int hpdn = hpd_hi - hpd_lo + 1;
    const size_t plsz = (size_t)(1 - kpd_lo)*hpdn;
    float2 *dp_c; float *dp_ct, *dp_rm; int *dp_nd; unsigned char *dp_vm;
    CKI(cudaMalloc(&dp_c,  (size_t)nrec*plsz*sizeof(float2)));
    CKI(cudaMalloc(&dp_ct, (size_t)nrec*plsz*sizeof(float)));
    CKI(cudaMalloc(&dp_rm, (size_t)nrec*nsym*9*sizeof(float)));
    CKI(cudaMalloc(&dp_nd, (size_t)nrec*sizeof(int)));
    CKI(cudaMalloc(&dp_vm, (size_t)nrec*sizeof(unsigned char)));
    if (hstage(pin_c,  dp_c,  pl_c,  (size_t)nrec*plsz*sizeof(float2))) return 5;
    if (hstage(pin_ct, dp_ct, pl_ct, (size_t)nrec*plsz*sizeof(float)))  return 5;
    CKI(cudaMemcpy(dp_rm, rotm,  (size_t)nrec*nsym*9*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dp_nd, nyqdisk, (size_t)nrec*sizeof(int),       cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dp_vm, vmask, (size_t)nrec*sizeof(unsigned char), cudaMemcpyHostToDevice));
    const int nh = hhi - hlo + 1, nk = khi - klo + 1;
    for (int off = 0; off < nrec; off += gs_nslab) {
        const int nc = (nrec - off < gs_nslab) ? (nrec - off) : gs_nslab;
        const size_t nwork1 = (size_t)nc*nh*nk;
        snr_insert_kernel<<<(int)((nwork1+255)/256), 256>>>(
            dp_c + (size_t)off*plsz, dp_ct + (size_t)off*plsz, dp_rm + (size_t)off*nsym*9,
            dp_nd + off, dp_vm + off, nc, nsym, hpd_lo, hpd_hi, kpd_lo,
            hlo, hhi, klo, khi, pf, lb[0], lb[1], lb[2], ub[0], ub[1], ub[2],
            tiny, gs_nvox, ds_slabC, ds_slabR);
        CKI(cudaGetLastError());
        const size_t nwork2 = (size_t)nc*gs_nvox;
        snr_sweep_kernel<<<(int)((nwork2+255)/256), 256>>>(
            dp_vm + off, nc, gs_nvox, inv_pf2, ds_slabC, ds_slabR, ds_var, ds_dens);
        CKI(cudaGetLastError());
    }
    CKI(cudaDeviceSynchronize());
    cudaFree(dp_c); cudaFree(dp_ct); cudaFree(dp_rm); cudaFree(dp_nd); cudaFree(dp_vm);
    return 0;
}

// resident variant: the adjoint residual pair is already on device (cb_c/cb_ct after
// estep_resid) with registered geometry — the slab insertion consumes it in place. The
// upload version above was the day-1 measured negative; residency is what flips it.
int flex_gpu_snr_batch_resident(
    const void* rotm, const void* nyqdisk, const void* vmask, int nrec, int nsym,
    const int* lb, const int* ub, float tiny, float inv_pf2)
{
    if (!ds_slabC) return 2;
    if (!cb_c.p || !cb_ct.p || ge_nrec < nrec) return 3;
    const int hs_lo = ge_hs_lo, hs_hi = ge_hs_hi, ks_lo = ge_ks_lo, khi = ge_khi;
    const int hpdn = hs_hi - hs_lo + 1;
    const size_t plsz = (size_t)(1 - ks_lo)*hpdn;
    DBUF(in_rm, (size_t)nrec*nsym*9*sizeof(float));
    DBUF(in_nd, (size_t)nrec*sizeof(int));
    DBUF(co_vm, (size_t)nrec*sizeof(unsigned char));
    float *dp_rm = (float*)in_rm.p;
    int   *dp_nd = (int*)in_nd.p;
    unsigned char *dp_vm = (unsigned char*)co_vm.p;
    CKI(cudaMemcpy(dp_rm, rotm,   (size_t)nrec*nsym*9*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dp_nd, nyqdisk,(size_t)nrec*sizeof(int),          cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dp_vm, vmask,  (size_t)nrec*sizeof(unsigned char),cudaMemcpyHostToDevice));
    const int nh = hpdn, nk = khi - ks_lo + 1;
    for (int off = 0; off < nrec; off += gs_nslab) {
        const int nc = (nrec - off < gs_nslab) ? (nrec - off) : gs_nslab;
        const size_t nwork1 = (size_t)nc*nh*nk;
        snr_insert_kernel<<<(int)((nwork1+255)/256), 256>>>(
            (const float2*)cb_c.p + (size_t)off*plsz, (const float*)cb_ct.p + (size_t)off*plsz,
            dp_rm + (size_t)off*nsym*9, dp_nd + off, dp_vm + off, nc, nsym,
            hs_lo, hs_hi, ks_lo, hs_lo, hs_hi, ks_lo, khi, OSMPL_PAD_FAC_DEV,
            lb[0], lb[1], lb[2], ub[0], ub[1], ub[2],
            tiny, gs_nvox, ds_slabC, ds_slabR);
        CKI(cudaGetLastError());
        const size_t nwork2 = (size_t)nc*gs_nvox;
        snr_sweep_kernel<<<(int)((nwork2+255)/256), 256>>>(
            dp_vm + off, nc, gs_nvox, inv_pf2, ds_slabC, ds_slabR, ds_var, ds_dens);
        CKI(cudaGetLastError());
    }
    CKI(cudaDeviceSynchronize());
    return 0;
}

} // extern "C"

// ---- polar reduced-solve particle pipeline (P2+P3) ----
// Replaces the blas_particle block of reduced_solve_accumulate_polar: the per-particle sample
// cache lives device-resident; per direction the CPU-built tables (Us, Cpk, Cm0, c00) come up
// with the chunk and the per-slice small GEMMs become chunk-wide cuBLAS calls; Vg never touches
// the host and the d^4 rank-k update runs as ONE resident FP64-EMULATED DGEMM per chunk
// (SYRK is not emulation-eligible; P0 measured the emulated GEMM at ~4x the native syrk).
namespace {
cublasHandle_t g_hb = nullptr, g_hb_emu = nullptr;
float  *dv_xws = nullptr, *dv_wrs = nullptr;          // resident cache: [np][nsamp2], [np][nk]
double *dv_Mspk = nullptr, *dv_Sbb = nullptr;         // resident accumulators
float  *dv_Us = nullptr;                              // chunk tables: [ndc][nsamp2][dt1]
double *dv_Cpk = nullptr, *dv_Cm0 = nullptr, *dv_c00 = nullptr;
int    *dv_dlist = nullptr;                           // chunk particle indices (cache columns)
float  *dv_X = nullptr;                               // gathered [mchunk][nsamp2]
double *dv_W = nullptr;                               // gathered+cast [mchunk][nk]
float  *dv_Reb = nullptr;                             // [mchunk][dt1]
double *dv_Vg = nullptr, *dv_Corr = nullptr, *dv_den = nullptr, *dv_B = nullptr, *dv_wrsum = nullptr;
long long gv_np = 0; int gv_nsamp2 = 0, gv_nk = 0, gv_dt = 0, gv_npk = 0;
int gv_maxdir = 0; long long gv_maxm = 0;
int gv_unitc = 0; float gv_alo = 0.1f, gv_ahi = 5.0f;

__global__ void solve_gather_kernel(const float* __restrict__ xws, const float* __restrict__ wrs,
    const int* __restrict__ dlist, int m, int nsamp2, int nk,
    float* __restrict__ X, double* __restrict__ W)
{
    const size_t tid = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    const size_t tot = (size_t)m*(nsamp2 + nk);
    if (tid >= tot) return;
    const int ii = (int)(tid / (nsamp2 + nk));
    const int j  = (int)(tid % (nsamp2 + nk));
    const int ip = dlist[ii] - 1;                      // 1-based cache column
    if (j < nsamp2) X[(size_t)ii*nsamp2 + j] = xws[(size_t)ip*nsamp2 + j];
    else            W[(size_t)ii*nk + (j - nsamp2)] = (double)wrs[(size_t)ip*nk + (j - nsamp2)];
}

// per-particle contrast + B assembly + per-direction W row-sums
__global__ void solve_contrast_kernel(const float* __restrict__ Reb, const double* __restrict__ Corr,
    const double* __restrict__ den, const double* __restrict__ W, int m, int dt, int nk,
    int unitc, float alo, float ahi, int diroff_of_first, const int* __restrict__ dir_of,
    double* __restrict__ B, double* __restrict__ wrsum)
{
    const int ii = blockIdx.x*blockDim.x + threadIdx.x;
    if (ii >= m) return;
    double cc = 1.0;
    if (!unitc) {
        const double d = den[ii];
        cc = (double)Reb[(size_t)ii*(dt+1)] / (d > 1e-300 ? d : 1e-300);
        if (cc < (double)alo) cc = (double)alo;
        if (cc > (double)ahi) cc = (double)ahi;
    }
    for (int q = 1; q <= dt; ++q)
        B[(size_t)ii*dt + (q-1)] = (double)Reb[(size_t)ii*(dt+1) + q] - cc*Corr[(size_t)ii*dt + (q-1)];
    const int idl = dir_of[ii] - diroff_of_first;      // local direction slot in chunk
    for (int k = 0; k < nk; ++k)
        atomicAdd(&wrsum[(size_t)idl*nk + k], W[(size_t)ii*nk + k]);
}
} // namespace

extern "C" {

int flex_gpu_solve_begin(long long np, int nsamp2, int nk, int dt, int npk,
                         int maxdir, long long maxm, int unitc, float alo, float ahi)
{
    gv_np = np; gv_nsamp2 = nsamp2; gv_nk = nk; gv_dt = dt; gv_npk = npk;
    gv_maxdir = maxdir; gv_maxm = maxm; gv_unitc = unitc; gv_alo = alo; gv_ahi = ahi;
    if (!g_hb)     { if (cublasCreate(&g_hb) != CUBLAS_STATUS_SUCCESS) return 3; }
    if (!g_hb_emu) {
        if (cublasCreate(&g_hb_emu) != CUBLAS_STATUS_SUCCESS) return 3;
#if defined(CUDART_VERSION) && (CUDART_VERSION >= 13000)
        // Blackwell FP64 emulation (CUDA 13+ only). On pre-13 toolkits / pre-Blackwell
        // GPUs the handle stays a plain one and the Dgemm below runs native FP64.
        cublasSetMathMode(g_hb_emu, (cublasMath_t)8);   // CUBLAS_FP64_EMULATED_FIXEDPOINT_MATH
        cublasSetEmulationStrategy(g_hb_emu, CUBLAS_EMULATION_STRATEGY_EAGER);
#endif
    }
    CKI(cudaMalloc(&dv_xws, (size_t)np*nsamp2*sizeof(float)));
    CKI(cudaMalloc(&dv_wrs, (size_t)np*nk*sizeof(float)));
    CKI(cudaMalloc(&dv_Mspk, (size_t)npk*npk*sizeof(double)));
    CKI(cudaMalloc(&dv_Sbb,  (size_t)dt*dt*sizeof(double)));
    CKI(cudaMemset(dv_Mspk, 0, (size_t)npk*npk*sizeof(double)));
    CKI(cudaMemset(dv_Sbb,  0, (size_t)dt*dt*sizeof(double)));
    CKI(cudaMalloc(&dv_Us,  (size_t)maxdir*nsamp2*(dt+1)*sizeof(float)));
    CKI(cudaMalloc(&dv_Cpk, (size_t)maxdir*npk*nk*sizeof(double)));
    CKI(cudaMalloc(&dv_Cm0, (size_t)maxdir*dt*nk*sizeof(double)));
    CKI(cudaMalloc(&dv_c00, (size_t)maxdir*nk*sizeof(double)));
    CKI(cudaMalloc(&dv_dlist, (size_t)maxm*sizeof(int)));
    CKI(cudaMalloc(&dv_X,   (size_t)maxm*nsamp2*sizeof(float)));
    CKI(cudaMalloc(&dv_W,   (size_t)maxm*nk*sizeof(double)));
    CKI(cudaMalloc(&dv_Reb, (size_t)maxm*(dt+1)*sizeof(float)));
    CKI(cudaMalloc(&dv_Vg,  (size_t)maxm*npk*sizeof(double)));
    CKI(cudaMalloc(&dv_Corr,(size_t)maxm*dt*sizeof(double)));
    CKI(cudaMalloc(&dv_den, (size_t)maxm*sizeof(double)));
    CKI(cudaMalloc(&dv_B,   (size_t)maxm*dt*sizeof(double)));
    CKI(cudaMalloc(&dv_wrsum, (size_t)maxdir*nk*sizeof(double)));
    return 0;
}

int flex_gpu_solve_cache(const void* xws, const void* wrs)
{
    if (!dv_xws) return 2;
    CKI(cudaMemcpy(dv_xws, xws, (size_t)gv_np*gv_nsamp2*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dv_wrs, wrs, (size_t)gv_np*gv_nk*sizeof(float),    cudaMemcpyHostToDevice));
    return 0;
}

// One chunk: ndc directions with per-direction particle counts dcnt[], concatenated particle
// index list dlist (1-based cache columns), per-direction tables. Returns wrsum (nk x ndc).
int flex_gpu_solve_chunk(const void* Us, const void* Cpk, const void* Cm0, const void* c00,
                         const void* dlist, const void* dcnt_v, int ndc, long long mtot,
                         void* wrsum_host)
{
    if (!dv_xws || ndc > gv_maxdir || mtot > gv_maxm) return 2;
    const int dt = gv_dt, dt1 = gv_dt + 1, nk = gv_nk, npk = gv_npk, ns2 = gv_nsamp2;
    const int* dcnt = (const int*)dcnt_v;
    CKI(cudaMemcpy(dv_Us,  Us,  (size_t)ndc*ns2*dt1*sizeof(float),  cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dv_Cpk, Cpk, (size_t)ndc*npk*nk*sizeof(double),  cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dv_Cm0, Cm0, (size_t)ndc*dt*nk*sizeof(double),   cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dv_c00, c00, (size_t)ndc*nk*sizeof(double),      cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dv_dlist, dlist, (size_t)mtot*sizeof(int),       cudaMemcpyHostToDevice));
    CKI(cudaMemset(dv_wrsum, 0, (size_t)ndc*nk*sizeof(double)));
    // gather the chunk's particles once
    {
        const size_t tot = (size_t)mtot*(ns2 + nk);
        solve_gather_kernel<<<(int)((tot+255)/256), 256>>>(dv_xws, dv_wrs, dv_dlist,
            (int)mtot, ns2, nk, dv_X, dv_W);
        CKI(cudaGetLastError());
    }
    const float  onef = 1.f, zerof = 0.f;
    const double oned = 1.0, zerod = 0.0;
    long long off = 0;
    for (int id = 0; id < ndc; ++id) {
        const int m = dcnt[id];
        if (m <= 0) continue;
        // Reb = Us^T X   (dt1 x m, SP)
        cublasSgemm(g_hb, CUBLAS_OP_T, CUBLAS_OP_N, dt1, m, ns2, &onef,
            dv_Us + (size_t)id*ns2*dt1, ns2, dv_X + off*ns2, ns2, &zerof,
            dv_Reb + off*dt1, dt1);
        // Vg = Cpk W   (npk x m, DP)
        cublasDgemm(g_hb, CUBLAS_OP_N, CUBLAS_OP_N, npk, m, nk, &oned,
            dv_Cpk + (size_t)id*npk*nk, npk, dv_W + off*nk, nk, &zerod,
            dv_Vg + off*npk, npk);
        // Corr = Cm0 W   (dt x m, DP)
        cublasDgemm(g_hb, CUBLAS_OP_N, CUBLAS_OP_N, dt, m, nk, &oned,
            dv_Cm0 + (size_t)id*dt*nk, dt, dv_W + off*nk, nk, &zerod,
            dv_Corr + off*dt, dt);
        // den = W^T c00   (m, DP)
        cublasDgemv(g_hb, CUBLAS_OP_T, nk, m, &oned, dv_W + off*nk, nk,
            dv_c00 + (size_t)id*nk, 1, &zerod, dv_den + off, 1);
        off += m;
    }
    // contrast + B + wrsum in one pass; dir_of encoded via a device array built from dcnt
    {
        // build dir_of on host (mtot ints) is wasteful; reuse dv_dlist slot trick instead:
        // we overwrite dv_dlist with the local direction id per particle (host-computed once).
        // (caller passes dlist already consumed by the gather above)
        int* h_dirof = (int*)malloc((size_t)mtot*sizeof(int));
        long long o = 0;
        for (int id = 0; id < ndc; ++id)
            for (int k = 0; k < dcnt[id]; ++k) h_dirof[o++] = id + 1;
        CKI(cudaMemcpy(dv_dlist, h_dirof, (size_t)mtot*sizeof(int), cudaMemcpyHostToDevice));
        free(h_dirof);
        solve_contrast_kernel<<<(int)((mtot+255)/256), 256>>>(dv_Reb, dv_Corr, dv_den, dv_W,
            (int)mtot, dt, nk, gv_unitc, gv_alo, gv_ahi, 1, dv_dlist, dv_B, dv_wrsum);
        CKI(cudaGetLastError());
    }
    // Sbb += B B^T (dt x dt), native DP
    cublasDsyrk(g_hb, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, dt, (int)mtot, &oned,
        dv_B, dt, &oned, dv_Sbb, dt);
    // Mspk += Vg Vg^T (npk x npk) -- the d^4 update, FP64-emulated GEMM
    cublasDgemm(g_hb_emu, CUBLAS_OP_N, CUBLAS_OP_T, npk, npk, (int)mtot, &oned,
        dv_Vg, npk, dv_Vg, npk, &oned, dv_Mspk, npk);
    CKI(cudaDeviceSynchronize());
    CKI(cudaMemcpy(wrsum_host, dv_wrsum, (size_t)ndc*nk*sizeof(double), cudaMemcpyDeviceToHost));
    return 0;
}

int flex_gpu_solve_end(void* Sbb_host, void* Mspk_host)
{
    if (!dv_Mspk) return 2;
    CKI(cudaMemcpy(Sbb_host,  dv_Sbb,  (size_t)gv_dt*gv_dt*sizeof(double),   cudaMemcpyDeviceToHost));
    CKI(cudaMemcpy(Mspk_host, dv_Mspk, (size_t)gv_npk*gv_npk*sizeof(double), cudaMemcpyDeviceToHost));
    cudaFree(dv_xws); dv_xws = nullptr; cudaFree(dv_wrs); dv_wrs = nullptr;
    cudaFree(dv_Mspk); dv_Mspk = nullptr; cudaFree(dv_Sbb); dv_Sbb = nullptr;
    cudaFree(dv_Us); dv_Us = nullptr; cudaFree(dv_Cpk); dv_Cpk = nullptr;
    cudaFree(dv_Cm0); dv_Cm0 = nullptr; cudaFree(dv_c00); dv_c00 = nullptr;
    cudaFree(dv_dlist); dv_dlist = nullptr; cudaFree(dv_X); dv_X = nullptr;
    cudaFree(dv_W); dv_W = nullptr; cudaFree(dv_Reb); dv_Reb = nullptr;
    cudaFree(dv_Vg); dv_Vg = nullptr; cudaFree(dv_Corr); dv_Corr = nullptr;
    cudaFree(dv_den); dv_den = nullptr; cudaFree(dv_B); dv_B = nullptr;
    cudaFree(dv_wrsum); dv_wrsum = nullptr;
    gv_np = 0;
    return 0;
}

int flex_gpu_snr_end(void* var_host, void* dens_host)
{
    if (!ds_var) return 2;
    CKI(cudaMemcpy(var_host,  ds_var,  gs_nvox*sizeof(float), cudaMemcpyDeviceToHost));
    CKI(cudaMemcpy(dens_host, ds_dens, gs_nvox*sizeof(float), cudaMemcpyDeviceToHost));
    cudaFree(ds_slabC); ds_slabC = nullptr;
    cudaFree(ds_slabR); ds_slabR = nullptr;
    cudaFree(ds_var);   ds_var   = nullptr;
    cudaFree(ds_dens);  ds_dens  = nullptr;
    gs_nvox = 0; gs_nslab = 0;
    return 0;
}

} // extern "C"

// ---- POLAR SHARED-DIRECTION E-STEP, device stage 2 (SIMPLE_COV_POLAR_GPU=1) ----
// Device port of the l_pol_es per-particle former in probe_subspace_iteration. The bank
// (UsallE + ring Gram tables) is HOST-built once per EM iteration and uploaded here once per
// iteration; per batch the raw (y, T) plane sections come up, the packed (conj(T)y, |T|^2)
// payload is derived on device, and the per-particle statistics are formed by
//   1. the RING part: a psample-clone gathers xw = sqwq*conj(T)y at the (bank direction,
//      relative in-plane) ring samples and the per-ring mean |T|^2; then Reb = Us^T xw
//      (float accumulate, the CPU path's sgemv) and the ring-table contractions
//      G/c/e_mm = tables . wrd (double, the CPU path's dgemv);
//   2. the HYBRID exact part: the UNMODIFIED estep_stats_kernel at the particle's own
//      exact direction with the disc capped at rhyb*(rhyb+1) -- its half-plane rule
//      (k<=0; on k=0 only h<=0) and per-axis-normalized KB taps are exactly the lattice
//      enumeration and interpolation of polar_hybrid_exact_accum. It reads the resident
//      volumes uploaded via flex_gpu_estep_vols (de_vols).
// Both parts accumulate into ONE per-record stat vector in the fused E-step's layout
// [e_mm, myv, b(1..nc), c(1..nc), G upper], which the host consumes exactly as it consumes
// the fused device stats (shared solve, mixture, ONLINE windows, e/o untouched).
// All state here is per-stage (ring geometry), per-iteration (bank) or per-batch (planes),
// so the path is window-safe under SIMPLE_COV_ONLINE by construction.
namespace {
float *dpe_rad = nullptr, *dpe_cs = nullptr, *dpe_sn = nullptr, *dpe_sqwq = nullptr;
int   *dpe_ring = nullptr, *dpe_nang = nullptr;
int    gpe_nsamp = 0, gpe_nk = 0, gpe_nc = 0, gpe_ndir = 0, gpe_nst = 0;
DBuf   pe_Us, pe_Cf, pe_Cm0, pe_c00;                       // the per-iteration bank
DBuf   pe_y, pe_t, pe_plc, pe_ct, pe_rm, pe_ca, pe_sa, pe_dir, pe_vm;
DBuf   pe_xw, pe_wrs, pe_wrd, pe_reb, pe_st, pe_p0;
PBuf   pe_pin_y, pe_pin_t;

// packed payload from the raw pair: plc = conj(T)*y, ct = |T|^2 (the coupled wrappers'
// host-side packing, moved on device so y and T stay available for the ring sampler)
__global__ void poles_payload_kernel(const float2* __restrict__ y, const float2* __restrict__ t,
    size_t n, float2* __restrict__ plc, float* __restrict__ ct)
{
    const size_t i = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= n) return;
    const float2 yv = y[i], tv = t[i];
    plc[i] = make_float2(tv.x*yv.x + tv.y*yv.y, tv.x*yv.y - tv.y*yv.x);
    ct[i]  = tv.x*tv.x + tv.y*tv.y;
}

// ring sampler: polar_sample_particle_fused minus taz/noise/halves. Interpolates y and T
// SEPARATELY (product-of-interpolants, the CPU convention) via ps_interp.
__global__ void poles_sample_kernel(
    const float2* __restrict__ ply, const float2* __restrict__ plt,
    const float* __restrict__ cav, const float* __restrict__ sav,
    const unsigned char* __restrict__ vmask,
    const float* __restrict__ rad, const float* __restrict__ cs, const float* __restrict__ sn,
    const float* __restrict__ sqwq, const int* __restrict__ ring_of,
    int nrec, int nsamp, int hpd_lo, int hpd_hi, int kpd_lo, int nk,
    float* __restrict__ xw, double* __restrict__ wrs_sum)
{
    const size_t tid = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (tid >= (size_t)nrec*nsamp) return;
    const int i = (int)(tid / nsamp);
    const int j = (int)(tid - (size_t)i*nsamp);
    if (!vmask[i]) return;
    const float ca = cav[i], sa = sav[i];
    const float c1 = cs[j]*ca - sn[j]*sa;
    const float s1 = sn[j]*ca + cs[j]*sa;
    const float hu =  rad[j]*c1;
    const float ku = -rad[j]*s1;
    const int hpdn = hpd_hi - hpd_lo + 1;
    const size_t plsz = (size_t)(1 - kpd_lo)*hpdn;
    float2 yf, p1, p2, tf, t1, t2c; float w1, w2, w3, w4;
    ps_interp(ply + (size_t)i*plsz, hpd_lo, hpd_hi, kpd_lo, hpdn, hu, ku, yf, p1, p2, w1, w2);
    ps_interp(plt + (size_t)i*plsz, hpd_lo, hpd_hi, kpd_lo, hpdn, hu, ku, tf, t1, t2c, w3, w4);
    const float sq = sqwq[j];
    xw[(size_t)i*2*nsamp + 2*j    ] = sq*( tf.x*yf.x + tf.y*yf.y);
    xw[(size_t)i*2*nsamp + 2*j + 1] = sq*( tf.x*yf.y - tf.y*yf.x);
    const double t2 = (double)tf.x*tf.x + (double)tf.y*tf.y;
    atomicAdd(&wrs_sum[(size_t)i*nk + ring_of[j]], t2);
}

__global__ void poles_wr_kernel(const double* __restrict__ wrs_sum,
    const int* __restrict__ nang, const unsigned char* __restrict__ vmask,
    int nrec, int nk, double* __restrict__ wrd)
{
    const size_t tid = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (tid >= (size_t)nrec*nk) return;
    const int i  = (int)(tid / nk);
    const int ir = (int)(tid - (size_t)i*nk);
    if (!vmask[i]) return;
    // round through float: the CPU forms wr in single (fused sampler) and casts to double
    wrd[(size_t)i*nk + ir] = (double)(float)(wrs_sum[(size_t)i*nk + ir] / (double)nang[ir]);
}

// Reb(0:nc) = Us(:,:,dir)^T xw -- one thread per (record, component), float accumulate
// (the CPU path's sgemv); Us layout (nsamp2, 0:nc, ndir) column-major from Fortran
__global__ void poles_reb_kernel(const float* __restrict__ Us, const float* __restrict__ xw,
    const int* __restrict__ dirid, const unsigned char* __restrict__ vmask,
    int nrec, int nc1, int ns2, float* __restrict__ Reb)
{
    const size_t tid = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (tid >= (size_t)nrec*nc1) return;
    const int i = (int)(tid / nc1);
    const int q = (int)(tid - (size_t)i*nc1);
    if (!vmask[i]) return;
    const int id = dirid[i] - 1;
    const float* u = Us + ((size_t)id*nc1 + q)*ns2;
    const float* x = xw + (size_t)i*ns2;
    float acc = 0.f;
    for (int j = 0; j < ns2; ++j) acc += u[j]*x[j];
    Reb[(size_t)i*nc1 + q] = acc;
}

// ring contributions += into the fused-stats layout; one thread per (record, stat index),
// each index written by exactly one thread AFTER the exact kernel finished (stream order),
// so plain read-modify-write. Table layouts are the Fortran column-major ones:
// Cf (nc*nc, nk, ndir) with (q,r) at (q-1)+(r-1)*nc; Cm0 (nc, nk, ndir); c00 (nk, ndir).
__global__ void poles_stats_kernel(const double* __restrict__ Cf, const double* __restrict__ Cm0,
    const double* __restrict__ c00, const double* __restrict__ wrd, const float* __restrict__ Reb,
    const int* __restrict__ dirid, const unsigned char* __restrict__ vmask,
    int nrec, int nc, int nk, int nst, double* __restrict__ stats)
{
    const size_t tid = (size_t)blockIdx.x*blockDim.x + threadIdx.x;
    if (tid >= (size_t)nrec*nst) return;
    const int i = (int)(tid / nst);
    const int s = (int)(tid - (size_t)i*nst);
    if (!vmask[i]) return;
    const int id = dirid[i] - 1;
    const int nc1 = nc + 1;
    const double* w = wrd + (size_t)i*nk;
    double acc = 0.0;
    if (s == 0) {                       // e_mm = c00 . wrd
        const double* c = c00 + (size_t)id*nk;
        for (int k = 0; k < nk; ++k) acc += c[k]*w[k];
    } else if (s == 1) {                // myv = Reb(0)
        acc = (double)Reb[(size_t)i*nc1];
    } else if (s <= 1 + nc) {           // b_q = Reb(q)
        acc = (double)Reb[(size_t)i*nc1 + (s - 1)];
    } else if (s <= 1 + 2*nc) {         // c_q = Cm0 . wrd
        const int q = s - 1 - nc;
        const double* m = Cm0 + ((size_t)id*nk)*nc + (q - 1);
        for (int k = 0; k < nk; ++k) acc += m[(size_t)k*nc]*w[k];
    } else {                            // G upper (q<=r), row-major pair order of the layout
        int g0 = s - 2 - 2*nc;
        int q = 1;
        for (; q <= nc; ++q) {
            const int len = nc - q + 1;
            if (g0 < len) break;
            g0 -= len;
        }
        const int r = q + g0;
        const double* cf = Cf + ((size_t)id*nk)*nc*nc + (size_t)(r-1)*nc + (q-1);
        for (int k = 0; k < nk; ++k) acc += cf[(size_t)k*nc*nc]*w[k];
    }
    stats[(size_t)i*nst + s] += acc;
}
}   // anon namespace

extern "C" {

// per-stage ring geometry (the E-step polar grid); the bank dims arrive with each bank upload
int flex_gpu_poles_begin(const void* rad, const void* cs, const void* sn, const void* sqwq,
    const int* ring_of, const int* nang, int nsamp, int nk)
{
    if (dpe_rad)  { cudaFree(dpe_rad);  dpe_rad  = nullptr; }
    if (dpe_cs)   { cudaFree(dpe_cs);   dpe_cs   = nullptr; }
    if (dpe_sn)   { cudaFree(dpe_sn);   dpe_sn   = nullptr; }
    if (dpe_sqwq) { cudaFree(dpe_sqwq); dpe_sqwq = nullptr; }
    if (dpe_ring) { cudaFree(dpe_ring); dpe_ring = nullptr; }
    if (dpe_nang) { cudaFree(dpe_nang); dpe_nang = nullptr; }
    gpe_nsamp = nsamp; gpe_nk = nk;
    CKI(cudaMalloc(&dpe_rad,  (size_t)nsamp*sizeof(float)));
    CKI(cudaMalloc(&dpe_cs,   (size_t)nsamp*sizeof(float)));
    CKI(cudaMalloc(&dpe_sn,   (size_t)nsamp*sizeof(float)));
    CKI(cudaMalloc(&dpe_sqwq, (size_t)nsamp*sizeof(float)));
    CKI(cudaMalloc(&dpe_ring, (size_t)nsamp*sizeof(int)));
    CKI(cudaMalloc(&dpe_nang, (size_t)nk*sizeof(int)));
    CKI(cudaMemcpy(dpe_rad,  rad,  (size_t)nsamp*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dpe_cs,   cs,   (size_t)nsamp*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dpe_sn,   sn,   (size_t)nsamp*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dpe_sqwq, sqwq, (size_t)nsamp*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dpe_ring, ring_of, (size_t)nsamp*sizeof(int), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(dpe_nang, nang, (size_t)nk*sizeof(int), cudaMemcpyHostToDevice));
    return 0;
}

// once per EM iteration, after the host bank build (rank may change between iterations)
int flex_gpu_poles_bank(const void* Us, const void* Cf, const void* Cm0, const void* c00,
    int nc, int ndir)
{
    if (gpe_nsamp <= 0) return 2;
    const int ns2 = 2*gpe_nsamp;
    gpe_nc = nc; gpe_ndir = ndir;
    gpe_nst = 2 + 2*nc + (nc*(nc+1))/2;
    DBUF(pe_Us,  (size_t)ndir*(nc+1)*ns2*sizeof(float));
    DBUF(pe_Cf,  (size_t)ndir*gpe_nk*nc*nc*sizeof(double));
    DBUF(pe_Cm0, (size_t)ndir*gpe_nk*nc*sizeof(double));
    DBUF(pe_c00, (size_t)ndir*gpe_nk*sizeof(double));
    CKI(cudaMemcpy(pe_Us.p,  Us,  (size_t)ndir*(nc+1)*ns2*sizeof(float),  cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(pe_Cf.p,  Cf,  (size_t)ndir*gpe_nk*nc*nc*sizeof(double), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(pe_Cm0.p, Cm0, (size_t)ndir*gpe_nk*nc*sizeof(double),  cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(pe_c00.p, c00, (size_t)ndir*gpe_nk*sizeof(double),     cudaMemcpyHostToDevice));
    return 0;
}

// one E-step batch: raw (y, T) plane sections + per-record geometry up, per-record stat
// vectors (ring + optional hybrid-exact contributions) back. nyq2_ex = rhyb*(rhyb+1) caps
// the exact-lattice disc; 0 = pure rings (no volume access, de_vols not required).
int flex_gpu_poles_batch(const void* ply, const void* plt, const void* rotm,
    const void* cav, const void* sav, const int* dirid, const void* vmask,
    int nrec, int hs_lo, int hs_hi, int ks_lo, int nyq2_ex, void* stats_host)
{
    if (gpe_nsamp <= 0 || gpe_nc <= 0) return 2;
    const int hpdn = hs_hi - hs_lo + 1;
    const size_t plsz = (size_t)(1 - ks_lo)*hpdn;
    const int nc = gpe_nc, nc1 = nc + 1, ns2 = 2*gpe_nsamp, nst = gpe_nst;
    if (nyq2_ex > 0) {
        if (!de_vols) return 2;
        if (ge_ncomp1 != nc1) return 6;   // resident volumes at a different rank
    }
    DBUF(pe_y,   (size_t)nrec*plsz*sizeof(float2));
    DBUF(pe_t,   (size_t)nrec*plsz*sizeof(float2));
    DBUF(pe_plc, (size_t)nrec*plsz*sizeof(float2));
    DBUF(pe_ct,  (size_t)nrec*plsz*sizeof(float));
    DBUF(pe_rm,  (size_t)nrec*9*sizeof(float));
    DBUF(pe_ca,  (size_t)nrec*sizeof(float));
    DBUF(pe_sa,  (size_t)nrec*sizeof(float));
    DBUF(pe_dir, (size_t)nrec*sizeof(int));
    DBUF(pe_vm,  (size_t)nrec*sizeof(unsigned char));
    DBUF(pe_xw,  (size_t)nrec*ns2*sizeof(float));
    DBUF(pe_wrs, (size_t)nrec*gpe_nk*sizeof(double));
    DBUF(pe_wrd, (size_t)nrec*gpe_nk*sizeof(double));
    DBUF(pe_reb, (size_t)nrec*nc1*sizeof(float));
    DBUF(pe_st,  (size_t)nrec*nst*sizeof(double));
    if (hstage(pe_pin_y, pe_y.p, ply, (size_t)nrec*plsz*sizeof(float2))) return 5;
    if (hstage(pe_pin_t, pe_t.p, plt, (size_t)nrec*plsz*sizeof(float2))) return 5;
    CKI(cudaMemcpy(pe_rm.p,  rotm,  (size_t)nrec*9*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(pe_ca.p,  cav,   (size_t)nrec*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(pe_sa.p,  sav,   (size_t)nrec*sizeof(float), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(pe_dir.p, dirid, (size_t)nrec*sizeof(int), cudaMemcpyHostToDevice));
    CKI(cudaMemcpy(pe_vm.p,  vmask, (size_t)nrec*sizeof(unsigned char), cudaMemcpyHostToDevice));
    CKI(cudaMemset(pe_st.p,  0, (size_t)nrec*nst*sizeof(double)));
    CKI(cudaMemset(pe_wrs.p, 0, (size_t)nrec*gpe_nk*sizeof(double)));
    CKI(cudaMemset(pe_xw.p,  0, (size_t)nrec*ns2*sizeof(float)));
    {
        const size_t n = (size_t)nrec*plsz;
        poles_payload_kernel<<<(int)((n+255)/256), 256>>>(
            (const float2*)pe_y.p, (const float2*)pe_t.p, n,
            (float2*)pe_plc.p, (float*)pe_ct.p);
        CKI(cudaGetLastError());
    }
    if (nyq2_ex > 0) {
        // hybrid-exact low-k statistics: the unmodified fused kernel with the disc capped
        const size_t npts = (size_t)hpdn*(0 - ks_lo + 1);
        DBUF(pe_p0, (size_t)nrec*npts*sizeof(float2));
        const size_t nwork = (size_t)nrec*npts;
        estep_stats_kernel<<<(int)((nwork+255)/256), 256>>>(
            (const float2*)pe_plc.p, (const float*)pe_ct.p, (const float*)pe_rm.p,
            (const unsigned char*)pe_vm.p, de_vols, nc1,
            nrec, hs_lo, hs_hi, ks_lo, 0, nyq2_ex,
            ge_lb[0], ge_lb[1], ge_lb[2], ge_ub[0], ge_ub[1], ge_ub[2], ge_nvox,
            (float2*)pe_p0.p, (double*)pe_st.p, nst);
        CKI(cudaGetLastError());
    }
    {
        const size_t nwork = (size_t)nrec*gpe_nsamp;
        poles_sample_kernel<<<(int)((nwork+255)/256), 256>>>(
            (const float2*)pe_y.p, (const float2*)pe_t.p, (const float*)pe_ca.p,
            (const float*)pe_sa.p, (const unsigned char*)pe_vm.p,
            dpe_rad, dpe_cs, dpe_sn, dpe_sqwq, dpe_ring,
            nrec, gpe_nsamp, hs_lo, hs_hi, ks_lo, gpe_nk,
            (float*)pe_xw.p, (double*)pe_wrs.p);
        CKI(cudaGetLastError());
    }
    {
        const size_t nwork = (size_t)nrec*gpe_nk;
        poles_wr_kernel<<<(int)((nwork+255)/256), 256>>>(
            (const double*)pe_wrs.p, dpe_nang, (const unsigned char*)pe_vm.p,
            nrec, gpe_nk, (double*)pe_wrd.p);
        CKI(cudaGetLastError());
    }
    {
        const size_t nwork = (size_t)nrec*nc1;
        poles_reb_kernel<<<(int)((nwork+255)/256), 256>>>(
            (const float*)pe_Us.p, (const float*)pe_xw.p, (const int*)pe_dir.p,
            (const unsigned char*)pe_vm.p, nrec, nc1, ns2, (float*)pe_reb.p);
        CKI(cudaGetLastError());
    }
    {
        const size_t nwork = (size_t)nrec*nst;
        poles_stats_kernel<<<(int)((nwork+255)/256), 256>>>(
            (const double*)pe_Cf.p, (const double*)pe_Cm0.p, (const double*)pe_c00.p,
            (const double*)pe_wrd.p, (const float*)pe_reb.p, (const int*)pe_dir.p,
            (const unsigned char*)pe_vm.p, nrec, nc, gpe_nk, nst, (double*)pe_st.p);
        CKI(cudaGetLastError());
    }
    CKI(cudaMemcpy(stats_host, pe_st.p, (size_t)nrec*nst*sizeof(double),
        cudaMemcpyDeviceToHost));
    return 0;
}

int flex_gpu_poles_free(void)
{
    if (dpe_rad)  { cudaFree(dpe_rad);  dpe_rad  = nullptr; }
    if (dpe_cs)   { cudaFree(dpe_cs);   dpe_cs   = nullptr; }
    if (dpe_sn)   { cudaFree(dpe_sn);   dpe_sn   = nullptr; }
    if (dpe_sqwq) { cudaFree(dpe_sqwq); dpe_sqwq = nullptr; }
    if (dpe_ring) { cudaFree(dpe_ring); dpe_ring = nullptr; }
    if (dpe_nang) { cudaFree(dpe_nang); dpe_nang = nullptr; }
    dbuf_release(pe_Us); dbuf_release(pe_Cf); dbuf_release(pe_Cm0); dbuf_release(pe_c00);
    dbuf_release(pe_y); dbuf_release(pe_t); dbuf_release(pe_plc); dbuf_release(pe_ct);
    dbuf_release(pe_rm); dbuf_release(pe_ca); dbuf_release(pe_sa); dbuf_release(pe_dir);
    dbuf_release(pe_vm); dbuf_release(pe_xw); dbuf_release(pe_wrs); dbuf_release(pe_wrd);
    dbuf_release(pe_reb); dbuf_release(pe_st); dbuf_release(pe_p0);
    gpe_nsamp = 0; gpe_nk = 0; gpe_nc = 0; gpe_ndir = 0; gpe_nst = 0;
    return 0;
}

} // extern "C"
