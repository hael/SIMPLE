!@descr: measure-only benchmark of the dense probability table's euclid scoring layouts (polar_matching_bench)
! Synthetic references and CTF-weighted particles at production polar sizes. Layout 1 is the production per-pair
! routine as fill_tab_range calls it; the replica and layouts 2-4 score the same data from a memo rebuilt here
! with FFTW plans of production's geometry. See doc/implementation_notes/planned/polar_matching_layout_bench.md.
module simple_polar_matching_bench
use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated, c_float, c_float_complex
use simple_defs,         only: sp, dp
use simple_error,        only: simple_exception
use simple_parameters,   only: parameters
use simple_polarft_calc, only: polarft_calc
!$ use omp_lib,          only: omp_get_wtime, omp_get_thread_num, omp_get_max_threads, omp_set_num_threads
implicit none
#include "simple_local_flags.inc"

public :: polar_matching_bench
public :: PMB_PER_PAIR, PMB_REPLICA, PMB_BLOCKED, PMB_STATES_INNER, PMB_MATRIX, PMB_MATRIX_T, PMB_PER_PAIR_CC
public :: PMB_NKINDS
private

integer, parameter :: PMB_PER_PAIR     = 1   !< layout 1: production gen_likelihood_val, particle against every reference
integer, parameter :: PMB_REPLICA      = 2   !< layout 1's sums and FFT on the benchmark memo (the identity baseline)
integer, parameter :: PMB_BLOCKED      = 3   !< layout 2: B particles x R references per tile, same per-pair sums
integer, parameter :: PMB_STATES_INNER = 4   !< layout 3: every state of a direction, one batched inverse FFT
integer, parameter :: PMB_MATRIX       = 5   !< layout 4: cgemm per angular bin over a tile, batched inverse FFT
integer, parameter :: PMB_MATRIX_T     = 6   !< layout 4 with the products transposed to contiguous FFT input
integer, parameter :: PMB_PER_PAIR_CC  = 7   !< layout 1 with objfun=cc (gen_corrs)
integer, parameter :: PMB_NKINDS       = 4   !< production kernels of prod_kernels: unshifted, shifted, bank draw, bank values
real(sp), parameter :: PMB_SHIFT(2)    = [1.3, -0.7]  !< the shift of the shifted and bank kernels (pixels)
integer, parameter :: PMB_BOX      = 180           !< native box (pixels)
integer, parameter :: PMB_BOX_CROP = 160
real,    parameter :: PMB_SMPD     = 2.055         !< native sampling (A)
real,    parameter :: PMB_MSKDIAM  = 270.          !< A
real,    parameter :: PMB_LP       = 6.            !< matching limit (A)
real,    parameter :: PMB_FALLOFF  = 40.           !< e-fold of the synthetic spectra's radial amplitude (shells)
real,    parameter :: PMB_KV = 300., PMB_CS = 2.7, PMB_FRACA = 0.1
integer, parameter :: PMB_SEED     = 20261021

! the BLAS of the build, called inside the OpenMP team: run with a single-threaded BLAS (OPENBLAS_NUM_THREADS=1)
interface
    subroutine cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        import :: sp
        character(len=1), intent(in)    :: transa, transb
        integer,          intent(in)    :: m, n, k, lda, ldb, ldc
        complex(sp),      intent(in)    :: alpha, beta
        complex(sp),      intent(in)    :: a(lda,*), b(ldb,*)
        complex(sp),      intent(inout) :: c(ldc,*)
    end subroutine cgemm
end interface

!> one thread's FFT buffers and draw work arrays
type :: pmb_buf
    type(c_ptr) :: p_vec = c_null_ptr, p_bat = c_null_ptr, p_bins = c_null_ptr, p_out = c_null_ptr, p_ct = c_null_ptr
    complex(c_float_complex), pointer :: cvec(:)     => null()   !< (pftsz+1), in place with rvec
    real(c_float),            pointer :: rvec(:)     => null()   !< (nrots+2)
    complex(c_float_complex), pointer :: cbat(:,:)   => null()   !< (pftsz+1, nstates), in place with rbat
    real(c_float),            pointer :: rbat(:,:)   => null()   !< (nrots+2, nstates)
    complex(c_float_complex), pointer :: cbins(:,:,:) => null()  !< (R, B, pftsz+1) cgemm output
    real(c_float),            pointer :: rout(:,:)   => null()   !< (nrots, R*B)
    complex(c_float_complex), pointer :: ct(:,:)     => null()   !< (pftsz+1, R*B) transposed products
    complex(sp), allocatable :: q(:,:,:)                         !< (2nk, B, pftsz) weighted particle block
    real(sp),    allocatable :: dists(:), pvec(:), cdf(:)
    integer,     allocatable :: inds(:)
    real(dp) :: t_parts(3) = 0._dp                               !< thread seconds: sums, inverse FFTs, distances and draw
end type pmb_buf

type :: polar_matching_bench
    private
    class(parameters),   allocatable :: params
    class(polarft_calc), allocatable :: pftc
    real,                allocatable :: sigma2(:,:)          !< noise power (shell, half), pointed at by pftc
    real(sp),            allocatable :: dists(:,:)           !< table rows (nrefs, nptcls): the drawn distance
    integer,             allocatable :: irots(:,:)           !< and its in-plane rotation
    ! the benchmark memo: production's arrays rebuilt (first index the angular bin, shells kfromto(1):kfromto(2))
    complex(sp),         allocatable :: ft_ref(:,:,:,:), ft_ref2(:,:,:,:)   !< (pftsz+1, k, nrefs, half)
    complex(sp),         allocatable :: ft_ptcl_ctf(:,:,:), ft_ctf2(:,:,:)  !< (pftsz+1, k, nptcls)
    complex(sp),         allocatable :: gmat(:,:,:,:)        !< layout 4 bin-major (2nk, nrefs, pftsz, half)
    real(dp),            allocatable :: wk8(:,:)             !< k/sigma2 (shell, half)
    real(dp),            allocatable :: norms(:)             !< wsqsum * 2 nrots per particle
    logical,             allocatable :: iseven(:)
    ! check records of particles 1..nkeep: full distance vectors and the nsample lowest with their rotations
    real(sp),            allocatable :: vecs(:,:,:)          !< (nrots, nrefs, nkeep)
    real(sp),            allocatable :: tops_d(:,:,:)        !< (nsample, nrefs, nkeep)
    integer,             allocatable :: tops_i(:,:,:)
    type(c_ptr) :: plan_c2r = c_null_ptr, plan_bat = c_null_ptr, plan_bins = c_null_ptr
    integer :: plan_bins_shape(3) = 0                      !< R, B and layout the plan was made for
    integer :: nspace = 0, nstates = 0, nrefs = 0, nptcls = 0, nthr_max = 0, nkeep = 0
    integer :: pftsz = 0, nrots = 0, kfromto(2) = 0, nk = 0, interpklim = 0, nsample = 0
    real(dp) :: t_gmat = 0._dp                                !< seconds spent building gmat
    real(dp) :: t_parts(3) = 0._dp                            !< last score's thread seconds per part (replica, matrix)
    logical :: with_prod_refs = .false., with_memo = .false., with_gmat = .false., l_cc = .false.
    logical :: exists = .false.
contains
    procedure :: new
    procedure :: kill
    procedure :: score
    procedure :: get_sizes
    procedure :: get_vecs
    procedure :: get_tops
    procedure :: get_gmat_build_time
    procedure :: get_parts
    procedure :: ref_bytes_per_pair
    procedure :: memo_gbytes
    procedure :: prod_kernels
end type polar_matching_bench

contains

    !> Synthetic references (nspace directions x nstates states, both halves) and nptcls CTF-weighted particles
    !! for teams of up to nthr_max threads: in the production calculator (prod_refs), in the benchmark memo
    !! (memo) and in layout 4's bin-major copy (gmat); the data are the same in every combination.
    !! objfun (euclid by default) is the objective of the production calculator.
    subroutine new( self, nspace, nstates, nptcls, nthr_max, prod_refs, memo, gmat, nkeep, objfun )
        use simple_cmdline,    only: cmdline
        use simple_sp_project, only: sp_project
        use simple_math_ft,    only: calc_fourier_index
        use simple_rnd,        only: ran3
        use simple_test_utils, only: set_fixed_seed
        use simple_fftw3,      only: fftwf_alloc_complex, fftwf_free, fftwf_plan_many_dft, fftwf_plan_many_dft_r2c,&
            &fftwf_execute_dft, fftwf_execute_dft_r2c, fftwf_destroy_plan, FFTW_FORWARD, FFTW_ESTIMATE
        use, intrinsic :: iso_c_binding, only: c_f_pointer, c_int, c_size_t
        class(polar_matching_bench), target, intent(inout) :: self
        integer,                             intent(in)    :: nspace, nstates, nptcls, nthr_max, nkeep
        logical,                             intent(in)    :: prod_refs, memo, gmat
        character(len=*), optional,          intent(in)    :: objfun
        class(cmdline),    allocatable :: cline
        class(sp_project), allocatable :: spproj
        complex(sp),       allocatable :: pft(:,:)
        real,              allocatable :: dfxs(:), dfys(:), angasts(:)
        complex(c_float_complex), pointer :: cmat(:,:), cmat_out(:,:), crc(:,:)
        real(c_float),            pointer :: crr(:,:)
        type(c_ptr) :: p_cmat, p_crmat, plan_fwd, plan_r2c
        real(dp)    :: t0
        integer     :: pdim(3), iref, ieo, i, k
        call self%kill
        if( nspace < 1 .or. nstates < 1 .or. nptcls < 1 .or. nthr_max < 1 ) THROW_HARD('invalid benchmark case')
        if( nkeep < 0 .or. nkeep > nptcls ) THROW_HARD('invalid number of kept particles')
        self%nspace   = nspace
        self%nstates  = nstates
        self%nrefs    = nspace * nstates
        self%nptcls   = nptcls
        self%nthr_max = nthr_max
        self%nkeep    = nkeep
        self%with_prod_refs = prod_refs
        self%with_memo      = memo
        self%with_gmat      = gmat
        self%l_cc           = .false.
        if( present(objfun) ) self%l_cc = objfun == 'cc'
        allocate(cline, self%params, self%pftc, spproj)
        call cline%set('box',      real(PMB_BOX))
        call cline%set('box_crop', real(PMB_BOX_CROP))
        call cline%set('smpd',     PMB_SMPD)
        call cline%set('mskdiam',  PMB_MSKDIAM)
        call cline%set('lp',       PMB_LP)
        call cline%set('nspace',   real(nspace))
        call cline%set('nstates',  real(nstates))
        call cline%set('nptcls',   real(nptcls))
        call cline%set('nthr',     real(nthr_max))
        call cline%set('ctf',      'no')
        if( self%l_cc )then
            call cline%set('objfun', 'cc')
        else
            call cline%set('objfun', 'euclid')
        endif
        call self%params%new(cline, silent=.true.)
        ! the shell range of set_bp_range at an explicit lp
        self%kfromto(1) = max(2, calc_fourier_index(self%params%hp, self%params%box, self%params%smpd))
        self%kfromto(2) = calc_fourier_index(self%params%lp, self%params%box, self%params%smpd)
        ! a fixed seed after parameters%new: reproducible data and FFTW_ESTIMATE plans (no wisdom file)
        call set_fixed_seed(PMB_SEED)
        ! without production references the calculator holds one (it serves the particle side)
        if( prod_refs )then
            call self%pftc%new(self%params, self%nrefs, [1,nptcls], self%kfromto)
        else
            call self%pftc%new(self%params, 1, [1,nptcls], self%kfromto)
        endif
        call self%pftc%set_with_ctf(.true.)
        self%pftsz      = self%pftc%get_pftsz()
        self%nrots      = self%pftc%get_nrots()
        self%nk         = self%kfromto(2) - self%kfromto(1) + 1
        pdim            = self%pftc%get_pdim_interp()
        self%interpklim = pdim(3)
        ! fill_tab_range's in-plane sample size on a fresh project (inpl_athres = prob_athres)
        self%nsample    = min(self%nrots, max(1, int(self%params%prob_athres * real(self%nrots) / 180.)))
        ! production's memo transforms, on buffers of its alignment
        p_cmat  = fftwf_alloc_complex(int(self%nrots * self%nk, c_size_t))
        call c_f_pointer(p_cmat, cmat, [self%nrots, self%nk])
        call c_f_pointer(p_cmat, cmat_out, [self%nrots, self%nk])
        p_crmat = fftwf_alloc_complex(int((self%pftsz+1) * self%nk, c_size_t))
        call c_f_pointer(p_crmat, crc, [self%pftsz+1, self%nk])
        call c_f_pointer(p_crmat, crr, [self%nrots+2, self%nk])
        plan_fwd = fftwf_plan_many_dft(1_c_int, [int(self%nrots,c_int)], int(self%nk,c_int), cmat, [int(self%nrots,c_int)],&
            &1_c_int, int(self%nrots,c_int), cmat_out, [int(self%nrots,c_int)], 1_c_int, int(self%nrots,c_int), FFTW_FORWARD,&
            &FFTW_ESTIMATE)
        plan_r2c = fftwf_plan_many_dft_r2c(1_c_int, [int(self%nrots,c_int)], int(self%nk,c_int), crr,&
            &[int(self%nrots,c_int)], 1_c_int, int(self%nrots+2,c_int), crc, [int(self%pftsz+1,c_int)], 1_c_int,&
            &int(self%pftsz+1,c_int), FFTW_ESTIMATE)
        if( memo ) allocate(self%ft_ref( self%pftsz+1,self%kfromto(1):self%kfromto(2),self%nrefs,2),&
                           &self%ft_ref2(self%pftsz+1,self%kfromto(1):self%kfromto(2),self%nrefs,2))
        if( gmat ) allocate(self%gmat(2*self%nk,self%nrefs,self%pftsz,2))
        self%t_gmat = 0._dp
        ! references: random complex values of unit scale with a smooth radial falloff in the matched band
        allocate(pft(self%pftsz,self%kfromto(1):self%interpklim))
        pft = cmplx(0.,0.,kind=sp)
        do iref = 1,self%nrefs
            do ieo = 1,2
                call fill_band(pft)
                if( prod_refs ) call self%pftc%set_ref_pft(iref, pft, iseven=(ieo == 1))
                if( memo .or. gmat ) call memoize_ref(iref, ieo)
            enddo
        enddo
        allocate(self%iseven(nptcls))
        do i = 1,nptcls
            call fill_band(pft)
            call self%pftc%set_ptcl_pft(i, pft)
            self%iseven(i) = mod(i,2) == 1
            call self%pftc%set_eo(i, self%iseven(i))
        enddo
        ! per-particle astigmatic CTFs through the production path: a one-stack project
        allocate(dfxs(nptcls), dfys(nptcls), angasts(nptcls))
        call spproj%os_stk%new(1, is_ptcl=.false.)
        call spproj%os_stk%set(1, 'fromp', 1)
        call spproj%os_stk%set(1, 'top',   nptcls)
        call spproj%os_stk%set(1, 'nptcls_stk', nptcls)
        call spproj%os_stk%set(1, 'box',   PMB_BOX)
        call spproj%os_stk%set(1, 'smpd',  PMB_SMPD)
        call spproj%os_stk%set(1, 'ctf',   'yes')
        call spproj%os_stk%set(1, 'kv',    PMB_KV)
        call spproj%os_stk%set(1, 'cs',    PMB_CS)
        call spproj%os_stk%set(1, 'fraca', PMB_FRACA)
        call spproj%os_ptcl3D%new(nptcls, is_ptcl=.true.)
        do i = 1,nptcls
            dfxs(i)    = 0.8 + 1.7 * ran3()
            dfys(i)    = dfxs(i) + 0.05 + 0.1 * ran3()
            angasts(i) = 180. * ran3()
            call spproj%os_ptcl3D%set(i, 'stkind', 1)
            call spproj%os_ptcl3D%set(i, 'indstk', i)
            call spproj%os_ptcl3D%set(i, 'dfx',    dfxs(i))
            call spproj%os_ptcl3D%set(i, 'dfy',    dfys(i))
            call spproj%os_ptcl3D%set(i, 'angast', angasts(i))
            call spproj%os_ptcl3D%set_state(i, 1)
        enddo
        call self%pftc%create_polar_absctfmats(spproj, 'ptcl3D')
        ! noise power of realistic shape: falling with resolution to a floor, the halves slightly apart
        allocate(self%sigma2(self%kfromto(1):self%kfromto(2),2), self%wk8(self%kfromto(1):self%kfromto(2),2))
        do k = self%kfromto(1),self%kfromto(2)
            self%sigma2(k,1) = 0.05 + exp(-real(k) / 20.)
            self%sigma2(k,2) = 1.02 * self%sigma2(k,1)
            self%wk8(k,:)    = real(k,dp) / real(self%sigma2(k,:),dp)
        enddo
        call self%pftc%assign_sigma2_noise(self%sigma2)
        call self%pftc%memoize_refs
        call self%pftc%memoize_ptcls
        allocate(self%norms(nptcls))
        do i = 1,nptcls
            self%norms(i) = self%pftc%wsqsum_ptcl_now(i) * real(2*self%nrots, dp)
        enddo
        if( memo .or. gmat ) call memoize_ptcls_here
        allocate(self%dists(self%nrefs,nptcls), source=0._sp)
        allocate(self%irots(self%nrefs,nptcls), source=0)
        if( nkeep > 0 )then
            allocate(self%vecs(self%nrots,self%nrefs,nkeep),     source=0._sp)
            allocate(self%tops_d(self%nsample,self%nrefs,nkeep), source=0._sp)
            allocate(self%tops_i(self%nsample,self%nrefs,nkeep), source=0)
        endif
        call fftwf_destroy_plan(plan_fwd)
        call fftwf_destroy_plan(plan_r2c)
        call fftwf_free(p_cmat)
        call fftwf_free(p_crmat)
        call spproj%kill
        call cline%kill
        deallocate(spproj, cline, pft, dfxs, dfys, angasts)
        self%exists = .true.

      contains

        subroutine fill_band( pft_band )
            complex(sp), intent(inout) :: pft_band(:,self%kfromto(1):)
            real    :: amp
            integer :: kk, p
            do kk = self%kfromto(1),self%kfromto(2)
                amp = exp(-real(kk) / PMB_FALLOFF)
                do p = 1,self%pftsz
                    pft_band(p,kk) = amp * cmplx(2. * ran3() - 1., 2. * ran3() - 1., kind=sp)
                enddo
            enddo
        end subroutine fill_band

        ! memoize_refs' transforms of pft (FT(REF), FT(REF2)) into the memo and the bin-major copy
        subroutine memoize_ref( iref_here, ieo_here )
            integer, intent(in) :: iref_here, ieo_here
            integer :: kk, k0, kh, j
            k0 = self%kfromto(1)
            do kh = self%kfromto(1), self%kfromto(2)
                kk = kh - k0 + 1
                cmat(1:self%pftsz,kk) = pft(:,kh)
                cmat(self%pftsz+1:self%nrots,kk) = conjg(cmat(1:self%pftsz,kk))
            enddo
            call fftwf_execute_dft(plan_fwd, cmat, cmat)
            do kh = self%kfromto(1), self%kfromto(2)
                kk = kh - k0 + 1
                crr(1:self%pftsz,kk) = real(pft(:,kh) * conjg(pft(:,kh)))
                crr(self%pftsz+1:self%nrots,kk) = crr(1:self%pftsz,kk)
            enddo
            call fftwf_execute_dft_r2c(plan_r2c, crr, crc)
            if( memo )then
                do kh = self%kfromto(1), self%kfromto(2)
                    kk = kh - k0 + 1
                    self%ft_ref(:,kh,iref_here,ieo_here)  = cmat(1:self%pftsz+1,kk)
                    self%ft_ref2(:,kh,iref_here,ieo_here) = conjg(crc(:,kk))
                enddo
            endif
            if( gmat )then
                t0 = 0._dp
                !$ t0 = omp_get_wtime()
                ! conjugated FT(REF2) so that cgemm's 'C' restores it; FT(REF) is conjugated by 'C' as the cross term needs
                do j = 1,self%pftsz
                    do kk = 1,self%nk
                        self%gmat(kk,iref_here,j,ieo_here)         = crc(j,kk)
                        self%gmat(self%nk+kk,iref_here,j,ieo_here) = cmat(j,kk)
                    enddo
                enddo
                !$ self%t_gmat = self%t_gmat + (omp_get_wtime() - t0)
            endif
        end subroutine memoize_ref

        ! memoize_ptcls' FT(X.CTF) (memoize_shifted_ptcl at zero shift: the same products and plan) and FT(CTF2)
        ! from create_polar_absctfmats' CTF matrices, recomputed here
        subroutine memoize_ptcls_here
            use simple_ctf, only: ctf
            type(ctf) :: tfun
            real(sp), allocatable :: px(:,:), py(:,:), spafreqsq(:,:), angs(:,:), ctfmat(:,:)
            complex(sp), allocatable :: bank(:,:)
            real    :: hinv, kinv
            integer :: irot, kh, kk, k0, ip
            k0 = self%kfromto(1)
            allocate(self%ft_ptcl_ctf(self%pftsz+1,self%kfromto(1):self%kfromto(2),nptcls),&
                    &self%ft_ctf2(self%pftsz+1,self%kfromto(1):self%kfromto(2),nptcls))
            allocate(px(self%kfromto(1):self%kfromto(2),self%pftsz), py(self%kfromto(1):self%kfromto(2),self%pftsz),&
                &spafreqsq(self%pftsz,self%kfromto(1):self%kfromto(2)), angs(self%pftsz,self%kfromto(1):self%kfromto(2)),&
                &ctfmat(self%pftsz,self%kfromto(1):self%kfromto(2)), bank(self%pftsz+1,self%nk))
            call self%pftc%get_polar_coords(1, px)
            call self%pftc%get_polar_coords(2, py)
            do irot = 1,self%pftsz
                do kh = self%kfromto(1),self%kfromto(2)
                    hinv = px(kh,irot) / PMB_BOX
                    kinv = py(kh,irot) / PMB_BOX
                    spafreqsq(irot,kh) = hinv*hinv+kinv*kinv
                    angs(irot,kh)      = atan2(py(kh,irot),px(kh,irot))
                end do
            end do
            do ip = 1,nptcls
                call self%pftc%memoize_shifted_ptcl(ip, [0._sp,0._sp], bank)
                self%ft_ptcl_ctf(:,:,ip) = bank
                tfun = ctf(PMB_SMPD, PMB_KV, PMB_CS, PMB_FRACA)
                call tfun%init(dfxs(ip), dfys(ip), angasts(ip))
                ctfmat = abs(tfun%eval_canonical(spafreqsq, angs, 0.))
                do kh = self%kfromto(1), self%kfromto(2)
                    kk = kh - k0 + 1
                    crr(1:self%pftsz,kk)            = ctfmat(:,kh) * ctfmat(:,kh)
                    crr(self%pftsz+1:self%nrots,kk) = crr(1:self%pftsz,kk)
                end do
                call fftwf_execute_dft_r2c(plan_r2c, crr, crc)
                do kh = self%kfromto(1), self%kfromto(2)
                    kk = kh - k0 + 1
                    self%ft_ctf2(:,kh,ip) = crc(:,kk)
                end do
            enddo
            deallocate(px, py, spafreqsq, angs, ctfmat, bank)
        end subroutine memoize_ptcls_here

    end subroutine new

    subroutine kill( self )
        use simple_fftw3, only: fftwf_destroy_plan
        class(polar_matching_bench), intent(inout) :: self
        if( allocated(self%pftc) )then
            call self%pftc%kill
            deallocate(self%pftc)
        endif
        if( allocated(self%params)      ) deallocate(self%params)
        if( allocated(self%sigma2)      ) deallocate(self%sigma2)
        if( allocated(self%dists)       ) deallocate(self%dists)
        if( allocated(self%irots)       ) deallocate(self%irots)
        if( allocated(self%ft_ref)      ) deallocate(self%ft_ref, self%ft_ref2)
        if( allocated(self%ft_ptcl_ctf) ) deallocate(self%ft_ptcl_ctf, self%ft_ctf2)
        if( allocated(self%gmat)        ) deallocate(self%gmat)
        if( allocated(self%wk8)         ) deallocate(self%wk8)
        if( allocated(self%norms)       ) deallocate(self%norms)
        if( allocated(self%iseven)      ) deallocate(self%iseven)
        if( allocated(self%vecs)        ) deallocate(self%vecs, self%tops_d, self%tops_i)
        if( c_associated(self%plan_c2r)  ) call fftwf_destroy_plan(self%plan_c2r)
        if( c_associated(self%plan_bat)  ) call fftwf_destroy_plan(self%plan_bat)
        if( c_associated(self%plan_bins) ) call fftwf_destroy_plan(self%plan_bins)
        self%plan_c2r  = c_null_ptr
        self%plan_bat  = c_null_ptr
        self%plan_bins = c_null_ptr
        self%plan_bins_shape = 0
        self%exists = .false.
    end subroutine kill

    !> Fills the table rows of particles 1..nptcls_use with the given layout on a team of nthr threads
    !! (B particles x R references per tile for layouts 2 and 4); returns the wall time of the scoring region
    function score( self, layout, nthr, nptcls_use, bsz, rsz ) result( wall )
        class(polar_matching_bench), intent(inout) :: self
        integer,                     intent(in)    :: layout, nthr, nptcls_use, bsz, rsz
        real(dp) :: wall
        type(pmb_buf), allocatable :: bufs(:)
        integer :: nteam, it
        if( .not. self%exists ) THROW_HARD('polar_matching_bench: score before new')
        if( nthr < 1 .or. nthr > self%nthr_max ) THROW_HARD('polar_matching_bench: team larger than the calculator')
        if( nptcls_use < 1 .or. nptcls_use > self%nptcls ) THROW_HARD('polar_matching_bench: invalid particle count')
        if( bsz < 1 .or. rsz < 1 ) THROW_HARD('polar_matching_bench: invalid tile')
        select case(layout)
            case(PMB_PER_PAIR, PMB_PER_PAIR_CC)
                if( .not. self%with_prod_refs ) THROW_HARD('polar_matching_bench: layout 1 needs production references')
                if( self%l_cc .neqv. (layout == PMB_PER_PAIR_CC) ) THROW_HARD('polar_matching_bench: objective mismatch')
            case(PMB_REPLICA, PMB_BLOCKED, PMB_STATES_INNER)
                if( .not. self%with_memo ) THROW_HARD('polar_matching_bench: layout needs the benchmark memo')
            case(PMB_MATRIX, PMB_MATRIX_T)
                if( .not. self%with_gmat ) THROW_HARD('polar_matching_bench: layout 4 needs the bin-major copy')
            case DEFAULT
                THROW_HARD('polar_matching_bench: unknown layout')
        end select
        !$ call omp_set_num_threads(nthr)
        nteam = 1
        !$ nteam = omp_get_max_threads()
        if( layout /= PMB_PER_PAIR .and. layout /= PMB_PER_PAIR_CC ) call make_bufs(self, layout, nteam, bsz, rsz, bufs)
        select case(layout)
            case(PMB_PER_PAIR, PMB_PER_PAIR_CC)
                wall = score_per_pair(self, nptcls_use)
            case(PMB_REPLICA)
                wall = score_replica(self, nptcls_use, bufs)
            case(PMB_BLOCKED)
                wall = score_blocked(self, nptcls_use, bsz, rsz, bufs)
            case(PMB_STATES_INNER)
                wall = score_states_inner(self, nptcls_use, bufs)
            case(PMB_MATRIX)
                wall = score_matrix(self, nptcls_use, bsz, rsz, bufs, .false.)
            case(PMB_MATRIX_T)
                wall = score_matrix(self, nptcls_use, bsz, rsz, bufs, .true.)
        end select
        self%t_parts = 0._dp
        if( allocated(bufs) )then
            do it = 1,size(bufs)
                self%t_parts = self%t_parts + bufs(it)%t_parts
            enddo
            call kill_bufs(bufs)
        endif
        !$ call omp_set_num_threads(self%nthr_max)
    end function score

    !> per-thread buffers for a team of nteam and the plans of the layout (made outside the parallel region)
    subroutine make_bufs( self, layout, nteam, bsz, rsz, bufs )
        use simple_fftw3, only: fftwf_alloc_complex, fftwf_plan_dft_c2r_1d, fftwf_plan_many_dft_c2r,&
            &fftwf_destroy_plan, FFTW_ESTIMATE
        use, intrinsic :: iso_c_binding, only: c_f_pointer, c_int, c_size_t
        class(polar_matching_bench), intent(inout) :: self
        integer,                     intent(in)    :: layout, nteam, bsz, rsz
        type(pmb_buf), allocatable,  intent(inout) :: bufs(:)
        integer :: ithr, n1, nst
        allocate(bufs(nteam))
        n1  = self%pftsz + 1
        nst = self%nstates
        do ithr = 1,nteam
            bufs(ithr)%p_vec = fftwf_alloc_complex(int(n1, c_size_t))
            call c_f_pointer(bufs(ithr)%p_vec, bufs(ithr)%cvec, [n1])
            call c_f_pointer(bufs(ithr)%p_vec, bufs(ithr)%rvec, [self%nrots+2])
            allocate(bufs(ithr)%dists(self%nrots), bufs(ithr)%pvec(self%nrots), bufs(ithr)%cdf(self%nrots),&
                &bufs(ithr)%inds(self%nrots))
            bufs(ithr)%cvec = cmplx(0.,0.,kind=c_float_complex)
            if( layout == PMB_STATES_INNER )then
                bufs(ithr)%p_bat = fftwf_alloc_complex(int(n1 * nst, c_size_t))
                call c_f_pointer(bufs(ithr)%p_bat, bufs(ithr)%cbat, [n1, nst])
                call c_f_pointer(bufs(ithr)%p_bat, bufs(ithr)%rbat, [self%nrots+2, nst])
                bufs(ithr)%cbat = cmplx(0.,0.,kind=c_float_complex)
            endif
            if( layout == PMB_MATRIX .or. layout == PMB_MATRIX_T )then
                bufs(ithr)%p_bins = fftwf_alloc_complex(int(rsz * bsz * n1, c_size_t))
                call c_f_pointer(bufs(ithr)%p_bins, bufs(ithr)%cbins, [rsz, bsz, n1])
                bufs(ithr)%p_out  = fftwf_alloc_complex(int((self%nrots * rsz * bsz + 1) / 2, c_size_t))
                call c_f_pointer(bufs(ithr)%p_out, bufs(ithr)%rout, [self%nrots, rsz * bsz])
                bufs(ithr)%cbins = cmplx(0.,0.,kind=c_float_complex)
                allocate(bufs(ithr)%q(2*self%nk, bsz, self%pftsz), source=cmplx(0.,0.,kind=sp))
            endif
            if( layout == PMB_MATRIX_T )then
                bufs(ithr)%p_ct = fftwf_alloc_complex(int(n1 * rsz * bsz, c_size_t))
                call c_f_pointer(bufs(ithr)%p_ct, bufs(ithr)%ct, [n1, rsz * bsz])
                bufs(ithr)%ct = cmplx(0.,0.,kind=c_float_complex)
            endif
        enddo
        ! production's plan_bwd1_single
        if( .not. c_associated(self%plan_c2r) ) self%plan_c2r = fftwf_plan_dft_c2r_1d(int(self%nrots, c_int),&
            &bufs(1)%cvec, bufs(1)%rvec, FFTW_ESTIMATE)
        if( layout == PMB_STATES_INNER .and. .not. c_associated(self%plan_bat) )then
            ! one in-place inverse transform per state, as production's plan_bwd_many_refs lays them out
            self%plan_bat = fftwf_plan_many_dft_c2r(1_c_int, [int(self%nrots,c_int)], int(nst,c_int), bufs(1)%cbat,&
                &[int(n1,c_int)], 1_c_int, int(n1,c_int), bufs(1)%rbat, [int(self%nrots,c_int)], 1_c_int,&
                &int(self%nrots+2,c_int), FFTW_ESTIMATE)
        endif
        if( (layout == PMB_MATRIX .or. layout == PMB_MATRIX_T) .and. any(self%plan_bins_shape /= [rsz, bsz, layout]) )then
            if( c_associated(self%plan_bins) ) call fftwf_destroy_plan(self%plan_bins)
            if( layout == PMB_MATRIX )then
                ! R*B transforms whose bins are R*B apart (cgemm's output layout), out of place
                self%plan_bins = fftwf_plan_many_dft_c2r(1_c_int, [int(self%nrots,c_int)], int(rsz*bsz,c_int),&
                    &bufs(1)%cbins, [int(n1,c_int)], int(rsz*bsz,c_int), 1_c_int, bufs(1)%rout,&
                    &[int(self%nrots,c_int)], 1_c_int, int(self%nrots,c_int), FFTW_ESTIMATE)
            else
                ! R*B contiguous transforms, out of place
                self%plan_bins = fftwf_plan_many_dft_c2r(1_c_int, [int(self%nrots,c_int)], int(rsz*bsz,c_int),&
                    &bufs(1)%ct, [int(n1,c_int)], 1_c_int, int(n1,c_int), bufs(1)%rout,&
                    &[int(self%nrots,c_int)], 1_c_int, int(self%nrots,c_int), FFTW_ESTIMATE)
            endif
            self%plan_bins_shape = [rsz, bsz, layout]
        endif
    end subroutine make_bufs

    subroutine kill_bufs( bufs )
        use simple_fftw3, only: fftwf_free
        type(pmb_buf), allocatable, intent(inout) :: bufs(:)
        integer :: ithr
        do ithr = 1,size(bufs)
            if( c_associated(bufs(ithr)%p_vec)  ) call fftwf_free(bufs(ithr)%p_vec)
            if( c_associated(bufs(ithr)%p_bat)  ) call fftwf_free(bufs(ithr)%p_bat)
            if( c_associated(bufs(ithr)%p_bins) ) call fftwf_free(bufs(ithr)%p_bins)
            if( c_associated(bufs(ithr)%p_out)  ) call fftwf_free(bufs(ithr)%p_out)
            if( c_associated(bufs(ithr)%p_ct)   ) call fftwf_free(bufs(ithr)%p_ct)
        enddo
        deallocate(bufs)
    end subroutine kill_bufs

    !> Layout 1: fill_tab_range's loop, a thread takes a particle and scores it against every reference in turn
    function score_per_pair( self, nptcls_use ) result( wall )
        class(polar_matching_bench), intent(inout) :: self
        integer,                     intent(in)    :: nptcls_use
        real(dp) :: wall
        real(sp), allocatable :: pvec_sorted(:,:)
        integer,  allocatable :: sorted_inds(:,:)
        real(sp) :: dist, corr
        integer  :: nteam, i, ithr, iref, irot
        nteam = 1
        !$ nteam = omp_get_max_threads()
        allocate(pvec_sorted(self%nrots,nteam), source=0._sp)
        allocate(sorted_inds(self%nrots,nteam), source=0)
        wall = 0._dp
        !$ wall = omp_get_wtime()
        !$omp parallel do default(shared) private(i,ithr,iref,dist,corr,irot) proc_bind(close) schedule(dynamic)
        do i = 1,nptcls_use
            ithr = 1
            !$ ithr = omp_get_thread_num() + 1
            do iref = 1,self%nrefs
                call self%pftc%gen_likelihood_val(iref, i, [0._sp,0._sp], self%nsample, dist, corr, irot,&
                    &pvec_sorted(:,ithr), sorted_inds(:,ithr))
                self%dists(iref,i) = dist
                self%irots(iref,i) = irot
                if( i <= self%nkeep )then
                    self%tops_d(:,iref,i) = pvec_sorted(1:self%nsample,ithr)
                    self%tops_i(:,iref,i) = sorted_inds(1:self%nsample,ithr)
                endif
            enddo
        enddo
        !$omp end parallel do
        !$ wall = omp_get_wtime() - wall
        deallocate(pvec_sorted, sorted_inds)
    end function score_per_pair

    !> The replica of layout 1 on the benchmark memo: particle-major, gen_euclid_crvec's sums and FFT per pair
    function score_replica( self, nptcls_use, bufs ) result( wall )
        use simple_fftw3, only: fftwf_execute_dft_c2r
        class(polar_matching_bench), intent(inout) :: self
        integer,                     intent(in)    :: nptcls_use
        type(pmb_buf),               intent(inout) :: bufs(:)
        real(dp) :: wall, t0, t1, t2, t3
        integer  :: i, ithr, iref
        wall = 0._dp
        t0 = 0._dp
        t1 = 0._dp
        t2 = 0._dp
        t3 = 0._dp
        !$ wall = omp_get_wtime()
        !$omp parallel do default(shared) private(i,ithr,iref,t0,t1,t2,t3) proc_bind(close) schedule(dynamic)
        do i = 1,nptcls_use
            ithr = 1
            !$ ithr = omp_get_thread_num() + 1
            do iref = 1,self%nrefs
                !$ t0 = omp_get_wtime()
                call pair_sums(self, i, iref, bufs(ithr)%cvec)
                !$ t1 = omp_get_wtime()
                call fftwf_execute_dft_c2r(self%plan_c2r, bufs(ithr)%cvec, bufs(ithr)%rvec)
                !$ t2 = omp_get_wtime()
                call pair_finish(self, i, iref, bufs(ithr)%rvec, bufs(ithr))
                !$ t3 = omp_get_wtime()
                bufs(ithr)%t_parts = bufs(ithr)%t_parts + [t1-t0, t2-t1, t3-t2]
            enddo
        enddo
        !$omp end parallel do
        !$ wall = omp_get_wtime() - wall
    end function score_replica

    !> Layout 2: a thread takes a block of B same-half particles; per tile of R references, each reference against
    !! every particle of the block, the per-pair computation unchanged
    function score_blocked( self, nptcls_use, bsz, rsz, bufs ) result( wall )
        class(polar_matching_bench), intent(inout) :: self
        integer,                     intent(in)    :: nptcls_use, bsz, rsz
        type(pmb_buf),               intent(inout) :: bufs(:)
        real(dp) :: wall
        integer, allocatable :: porder(:), blocks(:,:)
        integer  :: ib, ithr, r0, iref, jj, i
        call make_blocks(self, nptcls_use, bsz, porder, blocks)
        wall = 0._dp
        !$ wall = omp_get_wtime()
        !$omp parallel do default(shared) private(ib,ithr,r0,iref,jj,i) proc_bind(close) schedule(dynamic)
        do ib = 1,size(blocks,2)
            ithr = 1
            !$ ithr = omp_get_thread_num() + 1
            do r0 = 1,self%nrefs,rsz
                do iref = r0,min(r0+rsz-1,self%nrefs)
                    do jj = blocks(1,ib),blocks(2,ib)
                        i = porder(jj)
                        call pair_crvec(self, i, iref, bufs(ithr))
                        call pair_finish(self, i, iref, bufs(ithr)%rvec, bufs(ithr))
                    enddo
                enddo
            enddo
        enddo
        !$omp end parallel do
        !$ wall = omp_get_wtime() - wall
        deallocate(porder, blocks)
    end function score_blocked

    !> Layout 3: per particle and direction, the sums of every state's reference, then one batched inverse FFT
    function score_states_inner( self, nptcls_use, bufs ) result( wall )
        use simple_fftw3, only: fftwf_execute_dft_c2r
        class(polar_matching_bench), intent(inout) :: self
        integer,                     intent(in)    :: nptcls_use
        type(pmb_buf),               intent(inout) :: bufs(:)
        real(dp) :: wall
        integer  :: i, ithr, iproj, istate, iref
        wall = 0._dp
        !$ wall = omp_get_wtime()
        !$omp parallel do default(shared) private(i,ithr,iproj,istate,iref) proc_bind(close) schedule(dynamic)
        do i = 1,nptcls_use
            ithr = 1
            !$ ithr = omp_get_thread_num() + 1
            do iproj = 1,self%nspace
                do istate = 1,self%nstates
                    iref = (istate-1)*self%nspace + iproj
                    call pair_sums(self, i, iref, bufs(ithr)%cbat(:,istate))
                enddo
                call fftwf_execute_dft_c2r(self%plan_bat, bufs(ithr)%cbat, bufs(ithr)%rbat)
                do istate = 1,self%nstates
                    iref = (istate-1)*self%nspace + iproj
                    call pair_finish(self, i, iref, bufs(ithr)%rbat(:,istate), bufs(ithr))
                enddo
            enddo
        enddo
        !$omp end parallel do
        !$ wall = omp_get_wtime() - wall
    end function score_states_inner

    !> Layout 4: per block of B same-half particles and tile of R references, both sums of every angular bin as one
    !! cgemm over the shells (2nk deep), then R*B inverse FFTs in one plan
    function score_matrix( self, nptcls_use, bsz, rsz, bufs, l_transpose ) result( wall )
        use simple_fftw3, only: fftwf_execute_dft_c2r
        class(polar_matching_bench), intent(inout) :: self
        integer,                     intent(in)    :: nptcls_use, bsz, rsz
        type(pmb_buf),               intent(inout) :: bufs(:)
        logical,                     intent(in)    :: l_transpose
        real(dp) :: wall
        complex(sp), parameter :: ONE = cmplx(1.,0.,kind=sp), ZERO = cmplx(0.,0.,kind=sp)
        integer, allocatable :: porder(:), blocks(:,:)
        integer  :: ib, ithr, r0, nr, nb, j, ih, b, ir, i, kk, kh, n2k
        real(sp) :: wk
        real(dp) :: t0, t1, t2, t3
        call make_blocks(self, nptcls_use, bsz, porder, blocks)
        n2k  = 2 * self%nk
        wall = 0._dp
        !$ wall = omp_get_wtime()
        t0 = 0._dp
        t1 = 0._dp
        t2 = 0._dp
        t3 = 0._dp
        !$omp parallel do default(shared) private(ib,ithr,r0,nr,nb,j,ih,b,ir,i,kk,kh,wk,t0,t1,t2,t3) proc_bind(close)&
        !$omp& schedule(dynamic)
        do ib = 1,size(blocks,2)
            ithr = 1
            !$ ithr = omp_get_thread_num() + 1
            !$ t0 = omp_get_wtime()
            nb = blocks(2,ib) - blocks(1,ib) + 1
            ih = merge(1, 2, self%iseven(porder(blocks(1,ib))))
            ! the block's weighted particle terms: shells 1..nk the CTF2 term, nk+1..2nk the cross term
            do b = 1,nb
                i = porder(blocks(1,ib) + b - 1)
                do j = 1,self%pftsz
                    do kh = self%kfromto(1),self%kfromto(2)
                        kk = kh - self%kfromto(1) + 1
                        wk = real(self%wk8(kh,ih), sp)
                        bufs(ithr)%q(kk,b,j)         = wk * self%ft_ctf2(j,kh,i)
                        bufs(ithr)%q(self%nk+kk,b,j) = -2. * wk * self%ft_ptcl_ctf(j,kh,i)
                    enddo
                enddo
            enddo
            !$ t1 = omp_get_wtime()
            bufs(ithr)%t_parts(1) = bufs(ithr)%t_parts(1) + (t1 - t0)
            do r0 = 1,self%nrefs,rsz
                !$ t0 = omp_get_wtime()
                nr = min(rsz, self%nrefs - r0 + 1)
                do j = 1,self%pftsz
                    call cgemm('C', 'N', nr, nb, n2k, ONE, self%gmat(:,r0:r0+nr-1,j,ih), n2k, bufs(ithr)%q(:,:,j), n2k,&
                        &ZERO, bufs(ithr)%cbins(:,:,j), rsz)
                enddo
                !$ t1 = omp_get_wtime()
                if( l_transpose )then
                    ! each pair's bins contiguous; the angular Nyquist bin stays zero as in the per-pair sums
                    do b = 1,nb
                        do ir = 1,nr
                            bufs(ithr)%ct(1:self%pftsz,(b-1)*rsz+ir) = bufs(ithr)%cbins(ir,b,1:self%pftsz)
                            bufs(ithr)%ct(self%pftsz+1,(b-1)*rsz+ir) = ZERO
                        enddo
                    enddo
                    call fftwf_execute_dft_c2r(self%plan_bins, bufs(ithr)%ct, bufs(ithr)%rout)
                else
                    ! the angular Nyquist bin is dropped as in the per-pair sums; the c2r may have overwritten it
                    bufs(ithr)%cbins(:,:,self%pftsz+1) = ZERO
                    call fftwf_execute_dft_c2r(self%plan_bins, bufs(ithr)%cbins, bufs(ithr)%rout)
                endif
                !$ t2 = omp_get_wtime()
                do b = 1,nb
                    i = porder(blocks(1,ib) + b - 1)
                    do ir = 1,nr
                        call pair_finish(self, i, r0+ir-1, bufs(ithr)%rout(:,(b-1)*rsz+ir), bufs(ithr))
                    enddo
                enddo
                !$ t3 = omp_get_wtime()
                bufs(ithr)%t_parts = bufs(ithr)%t_parts + [t1-t0, t2-t1, t3-t2]
            enddo
        enddo
        !$omp end parallel do
        !$ wall = omp_get_wtime() - wall
        deallocate(porder, blocks)
    end function score_matrix

    !> particles 1..nptcls_use ordered even half first, cut into blocks of at most bsz that never mix halves
    subroutine make_blocks( self, nptcls_use, bsz, porder, blocks )
        class(polar_matching_bench), intent(in)    :: self
        integer,                     intent(in)    :: nptcls_use, bsz
        integer, allocatable,        intent(inout) :: porder(:), blocks(:,:)
        integer :: i, n, neven, nblk, j0
        allocate(porder(nptcls_use))
        n = 0
        do i = 1,nptcls_use
            if( self%iseven(i) )then
                n = n + 1
                porder(n) = i
            endif
        enddo
        neven = n
        do i = 1,nptcls_use
            if( .not. self%iseven(i) )then
                n = n + 1
                porder(n) = i
            endif
        enddo
        nblk = (neven + bsz - 1) / bsz + (nptcls_use - neven + bsz - 1) / bsz
        allocate(blocks(2,nblk))
        nblk = 0
        do j0 = 1,neven,bsz
            nblk = nblk + 1
            blocks(:,nblk) = [j0, min(j0+bsz-1, neven)]
        enddo
        do j0 = neven+1,nptcls_use,bsz
            nblk = nblk + 1
            blocks(:,nblk) = [j0, min(j0+bsz-1, nptcls_use)]
        enddo
    end subroutine make_blocks

    !> gen_euclid_crvec's sums (unshifted branch, angular bins 1..pftsz) and inverse FFT into cvec/rvec
    subroutine pair_crvec( self, i, iref, buf )
        use simple_fftw3, only: fftwf_execute_dft_c2r
        class(polar_matching_bench), intent(in)    :: self
        integer,                     intent(in)    :: i, iref
        type(pmb_buf),               intent(inout) :: buf
        call pair_sums(self, i, iref, buf%cvec)
        call fftwf_execute_dft_c2r(self%plan_c2r, buf%cvec, buf%rvec)
    end subroutine pair_crvec

    subroutine pair_sums( self, i, iref, cvec )
        class(polar_matching_bench), intent(in)    :: self
        integer,                     intent(in)    :: i, iref
        complex(c_float_complex),    intent(inout) :: cvec(:)
        complex(sp) :: c
        real(dp)    :: wk
        integer     :: k, p, ieo, ih
        ih  = merge(1, 2, self%iseven(i))
        ieo = ih
        cvec = cmplx(0.,0.,kind=c_float_complex)
        do k = self%kfromto(1), self%kfromto(2)
            wk = self%wk8(k,ih)
            do p = 1,self%pftsz
                c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%ft_ref(p,k,iref,ieo))
                cvec(p) = cmplx(cvec(p) + wk * c, kind=kind(cvec))
            enddo
        enddo
    end subroutine pair_sums

    !> euclid_dists_from_crvec over all rotations, then the production likelihood draw; stores the row entry
    !! and, for kept particles, the check records
    subroutine pair_finish( self, i, iref, rv, buf )
        use simple_defs,               only: TINY
        use simple_eul_prob_tab_utils, only: sample_likelihood_dists
        class(polar_matching_bench), intent(inout) :: self
        integer,                     intent(in)    :: i, iref
        real(c_float),               intent(in)    :: rv(:)
        type(pmb_buf),               intent(inout) :: buf
        real(dp) :: v
        real(sp) :: dist, corr
        integer  :: p, irot
        do p = 1,self%nrots
            v = 1.d0 + real(rv(p), dp) / self%norms(i)
            if( v > -log(real(TINY,dp)) )then
                buf%dists(p) = huge(dist)
            else
                buf%dists(p) = real(v, sp)
            endif
        enddo
        call sample_likelihood_dists(self%nrots, buf%dists, self%nsample, dist, corr, irot, buf%pvec, buf%inds, buf%cdf)
        self%dists(iref,i) = dist
        self%irots(iref,i) = irot
        if( i <= self%nkeep )then
            self%vecs(:,iref,i)   = buf%dists
            self%tops_d(:,iref,i) = buf%pvec(1:self%nsample)
            self%tops_i(:,iref,i) = buf%inds(1:self%nsample)
        endif
    end subroutine pair_finish

    !> The production kernels on particles 1..nptcls_use against every reference, team of nthr: the full distance
    !! vectors (nrots, nrefs, nptcls_use, PMB_NKINDS) of the unshifted kernel, the kernel at PMB_SHIFT, the bank draw
    !! and the bank values at PMB_SHIFT, the nsample lowest of gen_likelihood_val and, for euclid, gen_best_euclid_val.
    !! Euclid: gen_likelihood_val and gen_likelihood_val_bank with nsample = nrots sort every rotation's distance,
    !! gen_objfun_vals_bank gives exp(-loss). cc: gen_objfun_vals and gen_objfun_vals_bank give the correlations,
    !! gen_likelihood_val_bank the distances 1 - cc.
    subroutine prod_kernels( self, nthr, nptcls_use, vecs, tops_d, tops_i, best_d, best_i )
        class(polar_matching_bench), intent(inout) :: self
        integer,                     intent(in)    :: nthr, nptcls_use
        real(sp), allocatable,       intent(inout) :: vecs(:,:,:,:), tops_d(:,:,:), best_d(:,:)
        integer,  allocatable,       intent(inout) :: tops_i(:,:,:), best_i(:,:)
        real(sp),    allocatable :: pvec(:,:)
        integer,     allocatable :: inds(:,:)
        complex(sp), allocatable :: bank(:,:,:)
        real(sp) :: dist, corr
        integer  :: nteam, i, ithr, iref, irot, j
        if( .not. self%with_prod_refs ) THROW_HARD('polar_matching_bench: prod_kernels needs production references')
        if( nthr < 1 .or. nthr > self%nthr_max ) THROW_HARD('polar_matching_bench: team larger than the calculator')
        if( nptcls_use < 1 .or. nptcls_use > self%nptcls ) THROW_HARD('polar_matching_bench: invalid particle count')
        if( allocated(vecs)   ) deallocate(vecs)
        if( allocated(tops_d) ) deallocate(tops_d)
        if( allocated(tops_i) ) deallocate(tops_i)
        if( allocated(best_d) ) deallocate(best_d)
        if( allocated(best_i) ) deallocate(best_i)
        allocate(vecs(self%nrots,self%nrefs,nptcls_use,PMB_NKINDS), source=0._sp)
        allocate(tops_d(self%nsample,self%nrefs,nptcls_use), source=0._sp)
        allocate(tops_i(self%nsample,self%nrefs,nptcls_use), source=0)
        allocate(best_d(self%nrefs,nptcls_use), source=0._sp)
        allocate(best_i(self%nrefs,nptcls_use), source=0)
        !$ call omp_set_num_threads(nthr)
        nteam = 1
        !$ nteam = omp_get_max_threads()
        allocate(pvec(self%nrots,nteam), source=0._sp)
        allocate(inds(self%nrots,nteam), source=0)
        allocate(bank(self%pftsz+1,self%nk,nteam), source=cmplx(0.,0.,kind=sp))
        !$omp parallel do default(shared) private(i,ithr,iref,irot,j,dist,corr) proc_bind(close) schedule(dynamic)
        do i = 1,nptcls_use
            ithr = 1
            !$ ithr = omp_get_thread_num() + 1
            call self%pftc%memoize_shifted_ptcl(i, PMB_SHIFT, bank(:,:,ithr))
            do iref = 1,self%nrefs
                if( self%l_cc )then
                    call self%pftc%gen_objfun_vals(iref, i, [0._sp,0._sp], vecs(:,iref,i,1))
                    call self%pftc%gen_objfun_vals(iref, i, PMB_SHIFT, vecs(:,iref,i,2))
                else
                    call self%pftc%gen_likelihood_val(iref, i, [0._sp,0._sp], self%nrots, dist, corr, irot,&
                        &pvec(:,ithr), inds(:,ithr))
                    do j = 1,self%nrots
                        vecs(inds(j,ithr),iref,i,1) = pvec(j,ithr)
                    enddo
                    call self%pftc%gen_likelihood_val(iref, i, PMB_SHIFT, self%nrots, dist, corr, irot,&
                        &pvec(:,ithr), inds(:,ithr))
                    do j = 1,self%nrots
                        vecs(inds(j,ithr),iref,i,2) = pvec(j,ithr)
                    enddo
                    call self%pftc%gen_best_objfun_val(iref, i, [0._sp,0._sp], best_d(iref,i), best_i(iref,i))
                endif
                call self%pftc%gen_likelihood_val_bank(iref, i, bank(:,:,ithr), self%nrots, dist, corr, irot,&
                    &pvec(:,ithr), inds(:,ithr))
                do j = 1,self%nrots
                    vecs(inds(j,ithr),iref,i,3) = pvec(j,ithr)
                enddo
                call self%pftc%gen_objfun_vals_bank(iref, i, bank(:,:,ithr), vecs(:,iref,i,4))
                call self%pftc%gen_likelihood_val(iref, i, [0._sp,0._sp], self%nsample, dist, corr, irot,&
                    &pvec(:,ithr), inds(:,ithr))
                tops_d(:,iref,i) = pvec(1:self%nsample,ithr)
                tops_i(:,iref,i) = inds(1:self%nsample,ithr)
            enddo
        enddo
        !$omp end parallel do
        !$ call omp_set_num_threads(self%nthr_max)
        deallocate(pvec, inds, bank)
    end subroutine prod_kernels

    !> pftsz, nrots, kfromto, nsample, nrefs
    subroutine get_sizes( self, pftsz, nrots, kfromto, nsample, nrefs )
        class(polar_matching_bench), intent(in)  :: self
        integer,                     intent(out) :: pftsz, nrots, kfromto(2), nsample, nrefs
        pftsz   = self%pftsz
        nrots   = self%nrots
        kfromto = self%kfromto
        nsample = self%nsample
        nrefs   = self%nrefs
    end subroutine get_sizes

    !> the distance vectors of the kept particles from the last score
    subroutine get_vecs( self, vecs )
        class(polar_matching_bench), intent(in)    :: self
        real(sp), allocatable,       intent(inout) :: vecs(:,:,:)
        if( .not. allocated(self%vecs) ) THROW_HARD('polar_matching_bench: no kept particles')
        vecs = self%vecs
    end subroutine get_vecs

    !> the nsample lowest distances and their rotations of the kept particles from the last score
    subroutine get_tops( self, tops_d, tops_i )
        class(polar_matching_bench), intent(in)    :: self
        real(sp), allocatable,       intent(inout) :: tops_d(:,:,:)
        integer,  allocatable,       intent(inout) :: tops_i(:,:,:)
        if( .not. allocated(self%tops_d) ) THROW_HARD('polar_matching_bench: no kept particles')
        tops_d = self%tops_d
        tops_i = self%tops_i
    end subroutine get_tops

    !> thread seconds of the last score (replica, matrix): sums or cgemm, inverse FFTs, distances and draw
    function get_parts( self ) result( parts )
        class(polar_matching_bench), intent(in) :: self
        real(dp) :: parts(3)
        parts = self%t_parts
    end function get_parts

    real(dp) function get_gmat_build_time( self )
        class(polar_matching_bench), intent(in) :: self
        get_gmat_build_time = self%t_gmat
    end function get_gmat_build_time

    !> bytes of memoized reference spectra one pair reads: ft_ref and ft_ref2 of one half
    real(dp) function ref_bytes_per_pair( self )
        class(polar_matching_bench), intent(in) :: self
        ref_bytes_per_pair = 2._dp * real(self%pftsz+1,dp) * real(self%nk,dp) * 8._dp
    end function ref_bytes_per_pair

    !> GB of the reference and particle spectra a layout scores from
    real(dp) function memo_gbytes( self, layout )
        class(polar_matching_bench), intent(in) :: self
        integer,                     intent(in) :: layout
        real(dp) :: refs, ptcls
        refs  = 2._dp * self%ref_bytes_per_pair() * real(self%nrefs,dp)
        ptcls = self%ref_bytes_per_pair() * real(self%nptcls,dp)
        if( layout == PMB_MATRIX .or. layout == PMB_MATRIX_T )then
            ! the bin-major copy keeps bins 1..pftsz
            refs = refs * real(self%pftsz,dp) / real(self%pftsz+1,dp)
        endif
        memo_gbytes = (refs + ptcls) / 1.e9_dp
    end function memo_gbytes

end module simple_polar_matching_bench
