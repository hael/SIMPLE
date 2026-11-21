!@descr: measure-only benchmark of the dense probability table's euclid scoring layouts (polar_matching_bench)
! Synthetic references and CTF-weighted particles at production polar sizes. Layout 1 is the production per-pair
! routine as fill_tab_range calls it; the replica and layouts 2-4 score the same data from a memo rebuilt here
! with FFTW plans of production's geometry. See doc/implementation_notes/planned/polar_matching_layout_bench.md.
! The neighbourhood route (polar_neigh_bench) scores the neighbourhood table's candidates by the per-pair routine
! and by the block kernel's paths (doc/implementation_notes/planned/prob_table_matrix_scoring.md, decision 19).
module simple_polar_matching_bench
use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated, c_float, c_float_complex
use, intrinsic :: iso_fortran_env, only: int64
use simple_defs,         only: sp, dp
use simple_error,        only: simple_exception
use simple_parameters,   only: parameters
use simple_polarft_calc, only: polarft_calc, PFTC_BLOCK_NB, PFTC_BLOCK_NR
use simple_linalg,       only: cgemm
use simple_builder,      only: builder
use simple_oris,         only: oris
use simple_sym,          only: sym
!$ use omp_lib,          only: omp_get_wtime, omp_get_thread_num, omp_get_max_threads, omp_set_num_threads
implicit none
#include "simple_local_flags.inc"

public :: polar_matching_bench
public :: PMB_PER_PAIR, PMB_REPLICA, PMB_BLOCKED, PMB_STATES_INNER, PMB_MATRIX, PMB_MATRIX_T, PMB_PER_PAIR_CC, PMB_BLOCK_PROD
public :: PMB_NKINDS
public :: pmb_fill_result, PMB_ROUTE_3D_DENSE, PMB_ROUTE_3D_NEIGH, PMB_ROUTE_2D_DENSE, PMB_ROUTE_2D_SPARSE, PMB_NROUTES
public :: PMB_ROUTE_3D_SEEDED, PMB_ROUTE_3D_NOSHIFT, PMB_ROUTE_3D_STATE, PMB_ROUTE_3D_GEOM
public :: PMB_NDIFF, PMB_DIFF_NAMES, pmb_route_name, kill_fill_result, compare_fill_results
public :: polar_neigh_bench, pnb_pass, PNB_PATH_REPS, PNB_PATH_MEMBERS, PNB_PATH_CONE, PNB_NPATHS
private

integer, parameter :: PMB_PER_PAIR     = 1   !< layout 1: production gen_posterior_val, particle against every reference
integer, parameter :: PMB_REPLICA      = 2   !< layout 1's sums and FFT on the benchmark memo (the identity baseline)
integer, parameter :: PMB_BLOCKED      = 3   !< layout 2: B particles x R references per tile, same per-pair sums
integer, parameter :: PMB_STATES_INNER = 4   !< layout 3: every state of a direction, one batched inverse FFT
integer, parameter :: PMB_MATRIX       = 5   !< layout 4: cgemm per angular bin over a tile, batched inverse FFT
integer, parameter :: PMB_MATRIX_T     = 6   !< layout 4 with the products transposed to contiguous FFT input
integer, parameter :: PMB_PER_PAIR_CC  = 7   !< layout 1 with objfun=cc (gen_corrs)
integer, parameter :: PMB_BLOCK_PROD   = 8   !< production's block kernel (score_block_posterior) on blocks of B particles
integer, parameter :: PMB_NKINDS       = 4   !< production kernels of prod_kernels: unshifted, shifted, bank row, bank values
real(sp), parameter :: PMB_SHIFT(2)    = [1.3, -0.7]  !< the shift of the shifted and bank kernels (pixels)
integer, parameter :: PMB_BOX      = 180           !< native box (pixels)
integer, parameter :: PMB_BOX_CROP = 160
real,    parameter :: PMB_SMPD     = 2.055         !< native sampling (A)
real,    parameter :: PMB_MSKDIAM  = 270.          !< A
real,    parameter :: PMB_LP       = 6.            !< matching limit (A)
real,    parameter :: PMB_FALLOFF  = 40.           !< e-fold of the synthetic spectra's radial amplitude (shells)
real,    parameter :: PMB_KV = 300., PMB_CS = 2.7, PMB_FRACA = 0.1
integer, parameter :: PMB_SEED     = 20261021
! the fill harness (suite=fills of polar_matching_bench; prob_align_batch_pipeline.md, level 2 of section 8): each
! probabilistic route's batch pipeline (the fill, the selection and the extras of the matchers) on the synthetic
! data, in parts as the distributed parts run it, compared by particle and reference identity between thread
! counts, batch sizes and part counts
integer, parameter :: PMB_ROUTE_3D_DENSE  = 1   !< refine=prob on ptcl3D: the dense table (eul_prob_tab)
integer, parameter :: PMB_ROUTE_3D_NEIGH  = 2   !< refine=prob_neigh prob_neigh_mode=snhc on ptcl3D (eul_prob_tab_neigh)
integer, parameter :: PMB_ROUTE_2D_DENSE  = 3   !< refine=prob on ptcl2D: the dense rows (eul_prob_tab2D)
integer, parameter :: PMB_ROUTE_2D_SPARSE = 4   !< refine=prob_neigh on ptcl2D: the candidate stream (eul_prob_tab2D)
integer, parameter :: PMB_ROUTE_3D_SEEDED  = 5   !< the dense table without the polar shift map (shift_sigma=0): the seed and its bank
integer, parameter :: PMB_ROUTE_3D_NOSHIFT = 6   !< the dense table without shifts (trs=0)
integer, parameter :: PMB_ROUTE_3D_STATE   = 7   !< refine=prob_neigh prob_neigh_mode=state: the representatives, the peaks, the pooled cells
integer, parameter :: PMB_ROUTE_3D_GEOM    = 8   !< refine=prob_neigh prob_neigh_mode=geom: the cones about the previous direction
integer, parameter :: PMB_NROUTES         = 8
real,    parameter :: PMB_TRS   = 3.            !< the fills' shift window (pixels); the polar shift map is on, as in production
real,    parameter :: PMB_RES05 = 300.          !< the insertion cone's FSC=0.5 resolution of every state (A): above the mask diameter, so that every refined entry is an eligible extra on random rows (the extras path is exercised)
integer, parameter :: PMB_NDIFF = 7             !< the categories compare_fill_results counts
! the neighbourhood route (decision 19): bgal's late stages (stage 8 of solve3D: box 256 at 1.275 A cropped to 162,
! mask 180 A, lp 4.03 A, shells 2-81)
integer, parameter :: PNB_PATH_REPS    = 1     !< (a) state's coarse pass: a block against every state's representatives
integer, parameter :: PNB_PATH_MEMBERS = 2     !< (b) state's member pass: a pooled cell's members against its sub-block
integer, parameter :: PNB_PATH_CONE    = 3     !< (c) geom's cones: a reachable cell's members against its sub-block
integer, parameter :: PNB_NPATHS       = 3
integer, parameter :: PNB_BOX      = 256
integer, parameter :: PNB_BOX_CROP = 162
real,    parameter :: PNB_SMPD     = 1.275
real,    parameter :: PNB_MSKDIAM  = 180.
real,    parameter :: PNB_LP       = 4.03
real,    parameter :: PNB_SMOOTH_CUT  = 2.5    !< the references' smoothing reaches this many cell spacings
real,    parameter :: PNB_STATE_CORR  = 0.9    !< a second state's anchors: this share of the first's
real,    parameter :: PNB_NOISE       = 1.     !< the particles' noise amplitude relative to the references' (per shell)
character(len=*), parameter :: PMB_DIFF_NAMES(PMB_NDIFF) = [character(len=10) :: 'dists', 'inpls', 'shifts', 'temps',&
    &'seeds', 'assignment', 'extras']

!> the outputs of one fill-and-select run (run_fill), by particle and full reference identity ((state-1)*nspace+proj
!! in 3D, the class in 2D); huge, zero or false where a row holds no entry (the sparse routes' candidate lists)
type :: pmb_fill_result
    real(sp), allocatable :: dists(:,:)       !< (nrefs, nptcls) the rows' distances
    integer,  allocatable :: inpls(:,:)       !< (nrefs, nptcls) the drawn in-plane indices
    real(sp), allocatable :: shifts(:,:,:)    !< (2, nrefs, nptcls) the refined entries' shifts
    real(sp), allocatable :: temps(:)         !< (nptcls) the particle temperatures, which no draw enters
    real(sp), allocatable :: seed_shifts(:,:) !< (2, nptcls) the seeded shifts
    logical,  allocatable :: seed_has_sh(:)   !< (nptcls) their flags
    integer,  allocatable :: assign_ref(:)    !< (nptcls) the assigned reference
    integer,  allocatable :: assign_inpl(:)   !< (nptcls) its in-plane index
    real(sp), allocatable :: assign_vals(:,:) !< (4, nptcls) its distance, shift and posterior shift width
    real(sp), allocatable :: extras(:,:,:)    !< (5, topk, nptcls) the weighted insertion's extras: reference, in-plane, shift, weight (3D)
    integer :: nptcls = 0, nrefs = 0, nentries = 0 !< entries of the rows (dense: references x particles; sparse: the candidates)
    integer :: projs_ns = 0                         !< the dense 3D fill's refined-reference budget (the iteration's model)
end type pmb_fill_result

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
    complex(sp), allocatable :: q(:,:,:)                         !< (2nk, B, pftsz+1) weighted particle block
    real(sp),    allocatable :: dists(:), w(:)                   !< one row over the rotations, the draw's weights
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
    complex(sp),         allocatable :: gmat(:,:,:,:)        !< layout 4 bin-major (2nk, nrefs, pftsz+1, half)
    real(dp),            allocatable :: wk8(:,:)             !< k/sigma2 (shell, half)
    real(dp),            allocatable :: norms(:)             !< wsqsum * 2 nrots per particle
    real(sp),            allocatable :: temps(:)             !< the temperature of every draw (ptcl_posterior_temp)
    logical,             allocatable :: iseven(:)
    ! check records of particles 1..nkeep: the full rows over the rotations (the distances; the correlations with cc)
    real(sp),            allocatable :: vecs(:,:,:)          !< (nrots, nrefs, nkeep)
    type(c_ptr) :: plan_c2r = c_null_ptr, plan_bat = c_null_ptr, plan_bins = c_null_ptr
    integer :: plan_bins_shape(3) = 0                      !< R, B and layout the plan was made for
    integer :: nspace = 0, nstates = 0, nrefs = 0, nptcls = 0, nthr_max = 0, nkeep = 0
    integer :: pftsz = 0, nrots = 0, kfromto(2) = 0, nk = 0, interpklim = 0
    real(dp) :: t_gmat = 0._dp                                !< seconds spent building gmat
    real(dp) :: t_parts(3) = 0._dp                            !< last score's thread seconds per part (replica, matrix)
    logical :: with_prod_refs = .false., with_memo = .false., with_gmat = .false., l_cc = .false.
    logical :: with_ref_bins = .false.   !< the calculator's bin-major memo is made (layout PMB_BLOCK_PROD)
    logical :: exists = .false.
    class(builder), allocatable :: build  !< the fill harness's toolbox (fill_route > 0): the same references and particles in build%pftc, the particle field of previous poses, the direction space
    integer :: fill_route = 0            !< the harness's route (PMB_ROUTE_*; 0: no harness)
contains
    procedure :: new
    procedure :: run_fill
    procedure :: kill
    procedure :: score
    procedure :: get_sizes
    procedure :: get_vecs
    procedure :: get_entries
    procedure :: get_temps
    procedure :: scale_temps
    procedure :: get_gmat_build_time
    procedure :: get_parts
    procedure :: ref_bytes_per_pair
    procedure :: memo_gbytes
    procedure :: prod_kernels
end type polar_matching_bench

!> one pass of the neighbourhood route: the entries of the path's pairs (npairs, in the order of production's fill)
!! and their cost
type :: pnb_pass
    real(sp), allocatable :: dists(:)
    integer,  allocatable :: irots(:)
    real(dp) :: wall = 0._dp, tsec = 0._dp   !< the region's wall time; the threads' time in the scoring, summed
    integer(int64) :: nscored = 0_int64      !< pairs scored (the path's needed pairs and its wasted ones)
    integer(int64) :: nwasted = 0_int64      !< pairs scored outside the particle's cone (path c)
    integer(int64) :: ncalls  = 0_int64      !< block kernel calls
    integer(int64) :: nrows   = 0_int64      !< their particles, summed
    integer(int64) :: ncols   = 0_int64      !< their references, summed
    integer(int64) :: nslots  = 0_int64      !< the pairs their products computed (16 x 32 per tile)
    integer :: bsz(2) = 0                    !< the blocks of the halves (decision 3)
end type pnb_pass

!> The neighbourhood route (decision 19 of prob_table_matrix_scoring.md): the candidates of fill_tab_subspace_range's
!! state and geom modes on a direction space with its subspace map, scored by the per-pair routine and by the block
!! kernel over the same pairs: (a) a block's particles against every state's representatives; (b) each pooled cell's
!! members (every state's) against the particles of the block that pooled it; (c) each reachable cell's members
!! against the particles of the block whose cone reaches it, the entries inside a particle's cone recorded and the
!! others wasted. A cell's members are the list of one call, so that each reference keeps its column (decision 2's
!! revision) whatever the block. The references are smooth over the directions and every particle is its previous
!! direction's reference, in-plane rotated, under its CTF and noise, so that the peaks lie about that direction
type :: polar_neigh_bench
    private
    class(parameters),   allocatable :: params
    class(polarft_calc), allocatable :: pftc
    type(oris) :: eulspace                              !< the direction space (the point group's spiral)
    type(sym)  :: pgrpsyms
    real,     allocatable :: sigma2(:,:)
    real(sp), allocatable :: temps(:)                   !< the particles' posterior temperatures
    logical,  allocatable :: iseven(:)
    integer,  allocatable :: prev_proj(:)               !< the particles' previous directions
    integer,  allocatable :: sub_inds(:), full2sub(:)   !< the representatives' directions; every direction's cell
    integer,  allocatable :: cell_off(:), cell_proj(:)  !< the cells' member directions, ascending (CSR over cells)
    real,     allocatable :: cell_radius(:)             !< geom's cell radii (degrees)
    integer,  allocatable :: rep_refs(:)                !< (a) every state's representatives, state-major in subspace order
    integer,  allocatable :: npooled(:), pooled(:,:)    !< state mode: the particles' pooled cells (ncap, nptcls)
    integer,  allocatable :: mem_off(:), mem_ref(:)     !< (b) the members each particle scores (CSR, production's order)
    integer,  allocatable :: mem_cell(:,:,:)            !< (b) where a pooled cell's members of a state start (ncap, nstates, nptcls)
    logical,  allocatable :: reach(:,:)                 !< geom: the cells each particle's cone reaches (nsubs, nptcls)
    integer,  allocatable :: cone_off(:), cone_ref(:)   !< (c) the references inside each particle's cone (CSR, ascending)
    integer :: nspace = 0, nstates = 0, nrefs = 0, nptcls = 0, nsubs = 0, npeak = 0, ncap = 0, nthr_max = 0
    integer :: pftsz = 0, nrots = 0, kfromto(2) = 0, maxcell = 0
    real    :: cone = 0.
    logical :: exists = .false.
contains
    procedure :: new  => new_neigh
    procedure :: kill => kill_neigh
    procedure :: score_path
    procedure :: npairs
    procedure :: get_neigh_sizes
    procedure :: pooled_changes
end type polar_neigh_bench

contains

    !> Synthetic references (nspace directions x nstates states, both halves) and nptcls CTF-weighted particles
    !! for teams of up to nthr_max threads: in the production calculator (prod_refs), in the benchmark memo
    !! (memo) and in layout 4's bin-major copy (gmat); the data are the same in every combination.
    !! objfun (euclid by default) is the objective of the production calculator. fills_route (PMB_ROUTE_*) builds
    !! the fill harness's toolbox beside: the same data in a builder's calculator and a particle field (run_fill).
    subroutine new( self, nspace, nstates, nptcls, nthr_max, prod_refs, memo, gmat, nkeep, objfun, fills_route )
        use simple_cmdline,    only: cmdline
        use simple_sp_project, only: sp_project
        use simple_math_ft,    only: calc_fourier_index
        use simple_rnd,        only: ran3, rnd_scope_key, RND_PROG_REFINE3D, RND_PASS_ALIGN
        use simple_test_utils, only: set_fixed_seed
        use simple_fftw3,      only: fftwf_alloc_complex, fftwf_free, fftwf_plan_many_dft, fftwf_plan_many_dft_r2c,&
            &fftwf_execute_dft, fftwf_execute_dft_r2c, fftwf_destroy_plan, FFTW_FORWARD, FFTW_ESTIMATE
        use, intrinsic :: iso_c_binding, only: c_f_pointer, c_int, c_size_t
        class(polar_matching_bench), target, intent(inout) :: self
        integer,                             intent(in)    :: nspace, nstates, nptcls, nthr_max, nkeep
        logical,                             intent(in)    :: prod_refs, memo, gmat
        character(len=*), optional,          intent(in)    :: objfun
        integer,          optional,          intent(in)    :: fills_route
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
        self%fill_route = 0
        if( present(fills_route) ) self%fill_route = fills_route
        if( self%fill_route < 0 .or. self%fill_route > PMB_NROUTES ) THROW_HARD('invalid fill route')
        if( self%fill_route > 0 .and. (self%l_cc .or. .not. prod_refs) )&
            &THROW_HARD('the fill harness needs production references and objfun=euclid')
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
        call cline%set('ctf',      'yes')
        if( self%l_cc )then
            call cline%set('objfun', 'cc')
        else
            call cline%set('objfun', 'euclid')
        endif
        if( self%fill_route > 0 )then
            ! the fills' keys: the shift window, the point group of the direction space and the route (the polar
            ! shift map is on by default, shift_sigma, as in production)
            call cline%set('trs',  PMB_TRS)
            call cline%set('pgrp', 'c1')
            select case(self%fill_route)
                case(PMB_ROUTE_3D_DENSE)
                    call cline%set('oritype', 'ptcl3D')
                    call cline%set('refine',  'prob')
                case(PMB_ROUTE_3D_SEEDED)
                    call cline%set('oritype',     'ptcl3D')
                    call cline%set('refine',      'prob')
                    call cline%set('shift_sigma', 0.)
                case(PMB_ROUTE_3D_NOSHIFT)
                    call cline%set('oritype', 'ptcl3D')
                    call cline%set('refine',  'prob')
                    call cline%set('trs',     0.)
                case(PMB_ROUTE_3D_NEIGH)
                    call cline%set('oritype', 'ptcl3D')
                    call cline%set('refine',  'prob_neigh')
                    call cline%set('prob_neigh_mode', 'snhc')
                case(PMB_ROUTE_3D_STATE, PMB_ROUTE_3D_GEOM)
                    ! a subspace of a tenth of the directions (cells of about ten)
                    call cline%set('oritype', 'ptcl3D')
                    call cline%set('refine',  'prob_neigh')
                    call cline%set('nspace_sub', real(max(2, nspace / 10)))
                    if( self%fill_route == PMB_ROUTE_3D_STATE )then
                        call cline%set('prob_neigh_mode', 'state')
                    else
                        call cline%set('prob_neigh_mode', 'geom')
                    endif
                case(PMB_ROUTE_2D_DENSE)
                    call cline%set('oritype', 'ptcl2D')
                    call cline%set('refine',  'prob')
                    call cline%set('ncls',    real(self%nrefs))
                case(PMB_ROUTE_2D_SPARSE)
                    call cline%set('oritype', 'ptcl2D')
                    call cline%set('refine',  'prob_neigh')
                    call cline%set('ncls',    real(self%nrefs))
            end select
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
        ! the fill harness's calculator, given the same references and particles
        if( self%fill_route > 0 )then
            allocate(self%build)
            call self%build%pftc%new(self%params, self%nrefs, [1,nptcls], self%kfromto)
        endif
        self%pftsz      = self%pftc%get_pftsz()
        self%nrots      = self%pftc%get_nrots()
        self%nk         = self%kfromto(2) - self%kfromto(1) + 1
        pdim            = self%pftc%get_pdim_interp()
        self%interpklim = pdim(3)
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
        if( gmat ) allocate(self%gmat(2*self%nk,self%nrefs,self%pftsz+1,2))
        self%t_gmat = 0._dp
        ! references: random complex values of unit scale with a smooth radial falloff in the matched band
        allocate(pft(self%pftsz,self%kfromto(1):self%interpklim))
        pft = cmplx(0.,0.,kind=sp)
        do iref = 1,self%nrefs
            do ieo = 1,2
                call fill_band(pft)
                if( prod_refs ) call self%pftc%set_ref_pft(iref, pft, iseven=(ieo == 1))
                if( self%fill_route > 0 ) call self%build%pftc%set_ref_pft(iref, pft, iseven=(ieo == 1))
                if( memo .or. gmat ) call memoize_ref(iref, ieo)
            enddo
        enddo
        allocate(self%iseven(nptcls))
        do i = 1,nptcls
            call fill_band(pft)
            call self%pftc%set_ptcl_pft(i, pft)
            self%iseven(i) = mod(i,2) == 1
            call self%pftc%set_eo(i, self%iseven(i))
            if( self%fill_route > 0 )then
                call self%build%pftc%set_ptcl_pft(i, pft)
                call self%build%pftc%set_eo(i, self%iseven(i))
            endif
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
        call self%pftc%set_rnd_scope(rnd_scope_key(int(PMB_SEED, int64), RND_PROG_REFINE3D, RND_PASS_ALIGN, 1))
        allocate(self%norms(nptcls), self%temps(nptcls))
        do i = 1,nptcls
            self%norms(i) = self%pftc%wsqsum_ptcl_now(i) * real(2*self%nrots, dp)
            self%temps(i) = self%pftc%ptcl_posterior_temp(i)
        enddo
        if( memo .or. gmat ) call memoize_ptcls_here
        if( self%fill_route > 0 ) call build_fill_toolbox
        allocate(self%dists(self%nrefs,nptcls), source=0._sp)
        allocate(self%irots(self%nrefs,nptcls), source=0)
        if( nkeep > 0 ) allocate(self%vecs(self%nrots,self%nrefs,nkeep), source=0._sp)
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
                do j = 1,self%pftsz+1
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

        ! the harness's toolbox: the calculator the tables reach through the builder holds the same references,
        ! particles, CTFs and noise; the particle field carries the CTF records with the state (round robin, so
        ! that every state is populated), a previous pose of random direction, a previous class, the sampling
        ! bookkeeping calc_num2sample reads and the centring status of the polar shift map (odd particles not yet
        ! centred take the disk, even ones the window); the direction space is the point group's spiral
        subroutine build_fill_toolbox
            integer :: istate
            call self%build%pftc%create_polar_absctfmats(spproj, 'ptcl3D')
            call self%build%pftc%assign_sigma2_noise(self%sigma2)
            call self%build%pftc%memoize_refs
            call self%build%pftc%memoize_ptcls
            call self%build%spproj%os_ptcl3D%copy(spproj%os_ptcl3D, is_ptcl=.true.)
            do i = 1,nptcls
                istate = mod(i-1, nstates) + 1
                call self%build%spproj%os_ptcl3D%set_state(i, istate)
                call self%build%spproj%os_ptcl3D%set_euler(i, [360. * ran3(), 180. * ran3(), 360. * ran3()])
                call self%build%spproj%os_ptcl3D%set_shift(i, [0., 0.])
                call self%build%spproj%os_ptcl3D%set(i, 'corr',     0.5)
                call self%build%spproj%os_ptcl3D%set(i, 'sampled',  1.)
                call self%build%spproj%os_ptcl3D%set(i, 'dist',     20. * ran3())
                call self%build%spproj%os_ptcl3D%set(i, 'shincarg', ran3())
                call self%build%spproj%os_ptcl3D%set(i, 'shsig',    merge(-1., 0.5, mod(i,2) == 1))
                call self%build%spproj%os_ptcl3D%set(i, 'class',    real(min(self%nrefs, 1 + int(real(self%nrefs) * ran3()))))
                call self%build%spproj%os_ptcl3D%set(i, 'eo',       merge(0., 1., self%iseven(i)))
            enddo
            call self%build%spproj%os_ptcl2D%copy(self%build%spproj%os_ptcl3D, is_ptcl=.true.)
            call self%build%pgrpsyms%new('c1')
            call self%build%build_eulspace(self%params)
            if( any(self%fill_route == [PMB_ROUTE_3D_STATE, PMB_ROUTE_3D_GEOM]) ) call build_fill_subspaces
        end subroutine build_fill_toolbox

        ! the subspace map of the state and geom routes, as the builder makes it for them
        subroutine build_fill_subspaces
            use simple_srchspace_map, only: srchspace_map
            type(oris)          :: eulspace_sub
            type(srchspace_map) :: sub_mapper
            real, allocatable   :: sub_distmat(:,:)
            call eulspace_sub%new(self%params%nspace_sub, is_ptcl=.false.)
            call self%build%pgrpsyms%build_refspiral(eulspace_sub)
            allocate(sub_distmat(self%params%nspace_sub, self%params%nspace))
            call self%build%pgrpsyms%sym_dist_table(eulspace_sub, self%build%eulspace, sub_distmat)
            call sub_mapper%new(self%params%nspace, self%params%nspace_sub, sub_distmat)
            self%build%subspace_inds         = sub_mapper%get_sub2full_map()
            self%build%subspace_full2sub_map = sub_mapper%get_full2sub_map()
            call sub_mapper%kill
            call eulspace_sub%kill
            deallocate(sub_distmat)
        end subroutine build_fill_subspaces

    end subroutine new

    subroutine kill( self )
        use simple_fftw3, only: fftwf_destroy_plan
        class(polar_matching_bench), intent(inout) :: self
        if( allocated(self%build) )then
            nullify(self%build%spproj_field)
            call self%build%pftc%kill
            call self%build%spproj%kill
            call self%build%eulspace%kill
            call self%build%pgrpsyms%kill
            deallocate(self%build)
        endif
        self%fill_route = 0
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
        if( allocated(self%temps)       ) deallocate(self%temps)
        if( allocated(self%iseven)      ) deallocate(self%iseven)
        if( allocated(self%vecs)        ) deallocate(self%vecs)
        if( c_associated(self%plan_c2r)  ) call fftwf_destroy_plan(self%plan_c2r)
        if( c_associated(self%plan_bat)  ) call fftwf_destroy_plan(self%plan_bat)
        if( c_associated(self%plan_bins) ) call fftwf_destroy_plan(self%plan_bins)
        self%plan_c2r  = c_null_ptr
        self%plan_bat  = c_null_ptr
        self%plan_bins = c_null_ptr
        self%plan_bins_shape = 0
        self%with_ref_bins   = .false.
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
            case(PMB_BLOCK_PROD)
                if( .not. self%with_prod_refs .or. self%l_cc ) THROW_HARD('polar_matching_bench: the block kernel needs production references and euclid')
                if( bsz > PFTC_BLOCK_NB ) THROW_HARD('polar_matching_bench: the block kernel takes at most PFTC_BLOCK_NB particles')
                ! the kernel's memo, made once outside the timed region as a table's new_pipeline makes it
                if( .not. self%with_ref_bins ) call self%pftc%memoize_ref_bins
                self%with_ref_bins = .true.
            case DEFAULT
                THROW_HARD('polar_matching_bench: unknown layout')
        end select
        !$ call omp_set_num_threads(nthr)
        nteam = 1
        !$ nteam = omp_get_max_threads()
        if( layout /= PMB_PER_PAIR .and. layout /= PMB_PER_PAIR_CC .and. layout /= PMB_BLOCK_PROD )then
            call make_bufs(self, layout, nteam, bsz, rsz, bufs)
        endif
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
            case(PMB_BLOCK_PROD)
                wall = score_block_prod(self, nptcls_use, bsz, nteam)
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
            allocate(bufs(ithr)%dists(self%nrots), bufs(ithr)%w(self%nrots))
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
                allocate(bufs(ithr)%q(2*self%nk, bsz, self%pftsz+1), source=cmplx(0.,0.,kind=sp))
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
    !! (the posterior in-plane draw at the particle's temperature; with cc the correlations and their best rotation).
    !! A kept particle's check record is its full row, read again through the same kernel after the draw.
    function score_per_pair( self, nptcls_use ) result( wall )
        class(polar_matching_bench), intent(inout) :: self
        integer,                     intent(in)    :: nptcls_use
        real(dp) :: wall
        real(sp), allocatable :: row(:,:)
        real(sp) :: dist
        integer  :: nteam, i, ithr, iref, irot
        nteam = 1
        !$ nteam = omp_get_max_threads()
        allocate(row(self%nrots,nteam), source=0._sp)
        wall = 0._dp
        !$ wall = omp_get_wtime()
        !$omp parallel do default(shared) private(i,ithr,iref,dist,irot) proc_bind(close) schedule(dynamic)
        do i = 1,nptcls_use
            ithr = 1
            !$ ithr = omp_get_thread_num() + 1
            do iref = 1,self%nrefs
                if( self%l_cc )then
                    call self%pftc%gen_objfun_vals(iref, i, [0._sp,0._sp], row(:,ithr))
                    irot = maxloc(row(:,ithr), dim=1)
                    dist = 1. - row(irot,ithr)
                    if( i <= self%nkeep ) self%vecs(:,iref,i) = row(:,ithr)
                else
                    call self%pftc%gen_posterior_val(iref, i, [0._sp,0._sp], self%temps(i), dist, irot)
                    if( i <= self%nkeep ) call self%pftc%gen_euclid_dists(iref, i, [0._sp,0._sp], self%vecs(:,iref,i))
                endif
                self%dists(iref,i) = dist
                self%irots(iref,i) = irot
            enddo
        enddo
        !$omp end parallel do
        !$ wall = omp_get_wtime() - wall
        deallocate(row)
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
                do j = 1,self%pftsz+1
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
                do j = 1,self%pftsz+1
                    call cgemm('C', 'N', nr, nb, n2k, ONE, self%gmat(:,r0:r0+nr-1,j,ih), n2k, bufs(ithr)%q(:,:,j), n2k,&
                        &ZERO, bufs(ithr)%cbins(:,:,j), rsz)
                enddo
                !$ t1 = omp_get_wtime()
                if( l_transpose )then
                    ! each pair's bins contiguous, the angular Nyquist bin included as in the per-pair sums
                    do b = 1,nb
                        do ir = 1,nr
                            bufs(ithr)%ct(1:self%pftsz+1,(b-1)*rsz+ir) = bufs(ithr)%cbins(ir,b,1:self%pftsz+1)
                        enddo
                    enddo
                    call fftwf_execute_dft_c2r(self%plan_bins, bufs(ithr)%ct, bufs(ithr)%rout)
                else
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

    !> Production's block kernel as the dense fill calls it: a thread takes a block of B same-half particles and scores
    !! it against every reference (score_block_posterior: one cgemm per angular bin per tile of PFTC_BLOCK_NR
    !! references, the per-pair kernels' distances and keyed draws)
    function score_block_prod( self, nptcls_use, bsz, nteam ) result( wall )
        class(polar_matching_bench), intent(inout) :: self
        integer,                     intent(in)    :: nptcls_use, bsz, nteam
        real(dp) :: wall
        integer,  allocatable :: porder(:), blocks(:,:), irefs(:), irots_t(:,:)
        real(sp), allocatable :: temps_po(:), dists_t(:,:)
        integer  :: ib, ithr, nb, b, ir, i, j0, j
        call make_blocks(self, nptcls_use, bsz, porder, blocks)
        irefs    = [(j, j=1,self%nrefs)]
        temps_po = self%temps(porder)
        allocate(dists_t(self%nrefs*PFTC_BLOCK_NB,nteam), irots_t(self%nrefs*PFTC_BLOCK_NB,nteam))
        wall = 0._dp
        !$ wall = omp_get_wtime()
        !$omp parallel do default(shared) private(ib,ithr,nb,b,ir,i,j0) proc_bind(close) schedule(dynamic)
        do ib = 1,size(blocks,2)
            ithr = 1
            !$ ithr = omp_get_thread_num() + 1
            j0 = blocks(1,ib)
            nb = blocks(2,ib) - j0 + 1
            call self%pftc%score_block_posterior(nb, porder(j0:blocks(2,ib)), self%nrefs, irefs, temps_po(j0:blocks(2,ib)),&
                &dists_t(:,ithr), irots_t(:,ithr))
            do b = 1,nb
                i = porder(j0+b-1)
                do ir = 1,self%nrefs
                    self%dists(ir,i) = dists_t((b-1)*self%nrefs+ir,ithr)
                    self%irots(ir,i) = irots_t((b-1)*self%nrefs+ir,ithr)
                enddo
            enddo
        enddo
        !$omp end parallel do
        !$ wall = omp_get_wtime() - wall
        deallocate(porder, blocks, irefs, temps_po, dists_t, irots_t)
    end function score_block_prod

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

    !> gen_euclid_crvec's sums (unshifted branch, angular bins 1..pftsz+1) and inverse FFT into cvec/rvec
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
        real(sp)    :: wk
        integer     :: k, p, ieo, ih
        ih  = merge(1, 2, self%iseven(i))
        ieo = ih
        cvec = cmplx(0.,0.,kind=c_float_complex)
        do k = self%kfromto(1), self%kfromto(2)
            ! single precision, as production's sigma2_wk
            wk = real(k,sp) / self%sigma2(k,ih)
            do p = 1,self%pftsz+1
                c = self%ft_ctf2(p,k,i) * self%ft_ref2(p,k,iref,ieo)
                c = c - 2.0 * self%ft_ptcl_ctf(p,k,i) * conjg(self%ft_ref(p,k,iref,ieo))
                cvec(p) = cvec(p) + wk * c
            enddo
        enddo
    end subroutine pair_sums

    !> euclid_dists_from_crvec over all rotations, then the production posterior draw at the particle's temperature
    !! (draw_posterior_row, as draw_posterior_crvec); stores the row entry and, for kept particles, the check record
    subroutine pair_finish( self, i, iref, rv, buf )
        use simple_defs,               only: TINY
        use simple_eul_prob_tab_utils, only: draw_posterior_row
        use simple_rnd,                only: ran3
        class(polar_matching_bench), intent(inout) :: self
        integer,                     intent(in)    :: i, iref
        real(c_float),               intent(in)    :: rv(:)
        type(pmb_buf),               intent(inout) :: buf
        real(dp) :: v
        real(sp) :: dist
        integer  :: p, irot
        do p = 1,self%nrots
            v = 1.d0 + real(rv(p), dp) / self%norms(i)
            if( v > -log(real(TINY,dp)) )then
                buf%dists(p) = huge(dist)
            else
                buf%dists(p) = real(v, sp)
            endif
        enddo
        call draw_posterior_row(self%nrots, buf%dists, self%temps(i), buf%w, real(ran3(),dp), irot, dist)
        if( irot < 1 )then
            irot = 1
            dist = huge(dist)
        endif
        self%dists(iref,i) = dist
        self%irots(iref,i) = irot
        if( i <= self%nkeep ) self%vecs(:,iref,i) = buf%dists
    end subroutine pair_finish

    !> The production kernels on particles 1..nptcls_use against every reference, team of nthr: the full rows
    !! (nrots, nrefs, nptcls_use, PMB_NKINDS) of the unshifted kernel, the kernel at PMB_SHIFT, the bank entry at
    !! PMB_SHIFT and the bank values, and, for euclid, gen_best_euclid_val. Euclid: gen_euclid_dists and
    !! gen_euclid_dists_bank give the losses the posterior draws run over, gen_objfun_vals_bank exp(-loss). cc:
    !! gen_objfun_vals and gen_objfun_vals_bank give the correlations (kinds 1-3; kind 4 stays zero).
    subroutine prod_kernels( self, nthr, nptcls_use, vecs, best_d, best_i )
        class(polar_matching_bench), intent(inout) :: self
        integer,                     intent(in)    :: nthr, nptcls_use
        real(sp), allocatable,       intent(inout) :: vecs(:,:,:,:), best_d(:,:)
        integer,  allocatable,       intent(inout) :: best_i(:,:)
        complex(sp), allocatable :: bank(:,:,:)
        integer  :: nteam, i, ithr, iref
        if( .not. self%with_prod_refs ) THROW_HARD('polar_matching_bench: prod_kernels needs production references')
        if( nthr < 1 .or. nthr > self%nthr_max ) THROW_HARD('polar_matching_bench: team larger than the calculator')
        if( nptcls_use < 1 .or. nptcls_use > self%nptcls ) THROW_HARD('polar_matching_bench: invalid particle count')
        if( allocated(vecs)   ) deallocate(vecs)
        if( allocated(best_d) ) deallocate(best_d)
        if( allocated(best_i) ) deallocate(best_i)
        allocate(vecs(self%nrots,self%nrefs,nptcls_use,PMB_NKINDS), source=0._sp)
        allocate(best_d(self%nrefs,nptcls_use), source=0._sp)
        allocate(best_i(self%nrefs,nptcls_use), source=0)
        !$ call omp_set_num_threads(nthr)
        nteam = 1
        !$ nteam = omp_get_max_threads()
        allocate(bank(self%pftsz+1,self%nk,nteam), source=cmplx(0.,0.,kind=sp))
        !$omp parallel do default(shared) private(i,ithr,iref) proc_bind(close) schedule(dynamic)
        do i = 1,nptcls_use
            ithr = 1
            !$ ithr = omp_get_thread_num() + 1
            call self%pftc%memoize_shifted_ptcl(i, PMB_SHIFT, bank(:,:,ithr))
            do iref = 1,self%nrefs
                if( self%l_cc )then
                    call self%pftc%gen_objfun_vals(iref, i, [0._sp,0._sp], vecs(:,iref,i,1))
                    call self%pftc%gen_objfun_vals(iref, i, PMB_SHIFT, vecs(:,iref,i,2))
                    call self%pftc%gen_objfun_vals_bank(iref, i, bank(:,:,ithr), vecs(:,iref,i,3))
                else
                    call self%pftc%gen_euclid_dists(iref, i, [0._sp,0._sp], vecs(:,iref,i,1))
                    call self%pftc%gen_euclid_dists(iref, i, PMB_SHIFT, vecs(:,iref,i,2))
                    call self%pftc%gen_euclid_dists_bank(iref, i, bank(:,:,ithr), vecs(:,iref,i,3))
                    call self%pftc%gen_objfun_vals_bank(iref, i, bank(:,:,ithr), vecs(:,iref,i,4))
                    call self%pftc%gen_best_objfun_val(iref, i, [0._sp,0._sp], best_d(iref,i), best_i(iref,i))
                endif
            enddo
        enddo
        !$omp end parallel do
        !$ call omp_set_num_threads(self%nthr_max)
        deallocate(bank)
    end subroutine prod_kernels

    !> pftsz, nrots, kfromto, nrefs
    subroutine get_sizes( self, pftsz, nrots, kfromto, nrefs )
        class(polar_matching_bench), intent(in)  :: self
        integer,                     intent(out) :: pftsz, nrots, kfromto(2), nrefs
        pftsz   = self%pftsz
        nrots   = self%nrots
        kfromto = self%kfromto
        nrefs   = self%nrefs
    end subroutine get_sizes

    !> the rows of the kept particles from the last score
    subroutine get_vecs( self, vecs )
        class(polar_matching_bench), intent(in)    :: self
        real(sp), allocatable,       intent(inout) :: vecs(:,:,:)
        if( .not. allocated(self%vecs) ) THROW_HARD('polar_matching_bench: no kept particles')
        vecs = self%vecs
    end subroutine get_vecs

    !> the table entries of the last score: the drawn distance and in-plane index of every reference and particle
    subroutine get_entries( self, dists, irots )
        class(polar_matching_bench), intent(in)    :: self
        real(sp), allocatable,       intent(inout) :: dists(:,:)
        integer,  allocatable,       intent(inout) :: irots(:,:)
        dists = self%dists
        irots = self%irots
    end subroutine get_entries

    !> the temperatures of the posterior draws (ptcl_posterior_temp, times the last scale_temps factor)
    function get_temps( self ) result( temps )
        class(polar_matching_bench), intent(in) :: self
        real(sp), allocatable :: temps(:)
        temps = self%temps
    end function get_temps

    !> the draws' temperatures at factor times the particles' own (a huge factor draws the argmin rotation)
    subroutine scale_temps( self, factor )
        class(polar_matching_bench), intent(inout) :: self
        real,                        intent(in)    :: factor
        integer :: i
        do i = 1,self%nptcls
            self%temps(i) = factor * self%pftc%ptcl_posterior_temp(i)
        enddo
    end subroutine scale_temps

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
    real(dp) function memo_gbytes( self )
        class(polar_matching_bench), intent(in) :: self
        real(dp) :: refs, ptcls
        refs  = 2._dp * self%ref_bytes_per_pair() * real(self%nrefs,dp)
        ptcls = self%ref_bytes_per_pair() * real(self%nptcls,dp)
        memo_gbytes = (refs + ptcls) / 1.e9_dp
    end function memo_gbytes

    ! ---- the fill harness ----

    !> One fill-and-select of the harness's route on the pipeline of the matchers: the particles in nparts parts, each
    !! part a pipeline table over its index range (as the distributed parts are), the iteration's model computed once
    !! from the complete incoming field and set on every part (as each part computes it from the field the master
    !! wrote), the cone of the weighted insertion (3D); then batch by batch the fill and the selection on a team of
    !! nthr threads, the rows gathered from the batch's store before it is released; the temperatures, seeds,
    !! assignment and extras at the close. The seed is reset before the parts. The outputs are gathered by particle
    !! and full reference identity. With l_perturb, after each part's fill the field of the particles already filled
    !! is changed as a commit would change it (their distances, states, classes and shift increments), so that the
    !! later parts' rows show whether a fill reads other particles' outputs (the input-boundary regression of the
    !! iteration's model); the field is restored afterwards.
    subroutine run_fill( self, nthr, batchsz, nparts, res, l_perturb )
        use simple_eul_prob_tab,       only: eul_prob_tab
        use simple_eul_prob_tab_neigh, only: eul_prob_tab_neigh
        use simple_eul_prob_tab2D,     only: eul_prob_tab2D
        use simple_eul_prob_tab_utils, only: prob_insert_cone, prob_insert_pose, prob_insert_poses, calc_num2sample, &
            &shift_prior_sigma2, prob_candidate_buffer
        use simple_map_reduce,         only: split_nobjs_even
        use simple_test_utils,         only: set_fixed_seed
        use simple_oris,               only: oris
        class(polar_matching_bench), target, intent(inout) :: self
        integer,                             intent(in)    :: nthr, batchsz, nparts
        type(pmb_fill_result),               intent(inout) :: res
        logical, optional,                   intent(in)    :: l_perturb
        type(eul_prob_tab),       allocatable :: tab3D
        type(eul_prob_tab_neigh), allocatable :: tabnb
        type(eul_prob_tab2D),     allocatable :: tab2D
        type(prob_insert_cone) :: cone
        type(oris)             :: field_saved
        integer,  allocatable  :: parts(:,:), pinds(:)
        integer :: ipart, nloc, first, last, i, ioff, projs_ns_model
        real    :: sigma2_model
        logical :: l_pert
        l_pert = .false.
        if( present(l_perturb) ) l_pert = l_perturb
        if( self%fill_route < 1 ) THROW_HARD('polar_matching_bench: no fill harness (new with fills_route)')
        if( nthr < 1 .or. nthr > self%nthr_max ) THROW_HARD('polar_matching_bench: team larger than the calculator')
        if( batchsz < 1 .or. nparts < 1 .or. nparts > self%nptcls ) THROW_HARD('polar_matching_bench: invalid batch size or part count')
        call kill_fill_result(res)
        select case(self%fill_route)
            case(PMB_ROUTE_3D_DENSE, PMB_ROUTE_3D_NEIGH, PMB_ROUTE_3D_SEEDED, PMB_ROUTE_3D_NOSHIFT, PMB_ROUTE_3D_STATE,&
                &PMB_ROUTE_3D_GEOM)
                self%build%spproj_field => self%build%spproj%os_ptcl3D
            case DEFAULT
                self%build%spproj_field => self%build%spproj%os_ptcl2D
        end select
        !$ call omp_set_num_threads(nthr)
        call set_fixed_seed(PMB_SEED, propagate=.true.)
        if( l_pert ) call field_saved%copy(self%build%spproj_field, .true.)
        parts = split_nobjs_even(self%nptcls, nparts)
        ! the iteration's model once, from the complete incoming field, as every part computes it from the field
        ! the master wrote (the perturbation of the earlier parts' rows must not enter it)
        sigma2_model   = iteration_shift_prior()
        projs_ns_model = 0
        if( any(self%fill_route == [PMB_ROUTE_3D_DENSE, PMB_ROUTE_3D_SEEDED, PMB_ROUTE_3D_NOSHIFT]) )then
            allocate(tab3D)
            call tab3D%new_pipeline(self%params, self%build, [(i, i=1,self%nptcls)])
            projs_ns_model = dense_projs_ns(tab3D)
            call tab3D%kill
            deallocate(tab3D)
        endif
        call cone%new([(PMB_RES05, i=1,self%nstates)], PMB_MSKDIAM, PMB_SMPD)
        call alloc_result
        res%projs_ns = projs_ns_model
        do ipart = 1,nparts
            pinds = [(i, i=parts(ipart,1),parts(ipart,2))]
            nloc  = size(pinds)
            ioff  = parts(ipart,1) - 1
            select case(self%fill_route)
                case(PMB_ROUTE_3D_DENSE, PMB_ROUTE_3D_SEEDED, PMB_ROUTE_3D_NOSHIFT)
                    allocate(tab3D)
                    call tab3D%new_pipeline(self%params, self%build, pinds, block_memo=.true.)
                    call tab3D%set_iteration_model(projs_ns_model, sigma2_model)
                    call tab3D%set_insert_cone(cone)
                    call tab3D%open_assignment(nglob=self%nptcls, rank0=ioff)
                    first = 1
                    do while( first <= nloc )
                        last = min(nloc, first + batchsz - 1)
                        call tab3D%fill_tab_range(first, last)
                        call gather_dense_batch3D(tab3D, first, last)
                        call tab3D%assign_batch(first, last)
                        first = last + 1
                    enddo
                    call tab3D%close_assignment
                    call gather_tab3D(tab3D)
                    call tab3D%kill
                    deallocate(tab3D)
                case(PMB_ROUTE_3D_NEIGH, PMB_ROUTE_3D_STATE, PMB_ROUTE_3D_GEOM)
                    allocate(tabnb)
                    call tabnb%new_neigh_pipeline(self%params, self%build, pinds, block_memo=.true.)
                    call tabnb%set_iteration_model(0, sigma2_model)
                    call tabnb%set_insert_cone(cone)
                    call tabnb%open_assignment(nglob=self%nptcls, rank0=ioff)
                    first = 1
                    do while( first <= nloc )
                        last = min(nloc, first + batchsz - 1)
                        call tabnb%fill_tab_range(first, last)
                        call gather_buffers(tabnb%candidate_buffers, .true.)
                        call tabnb%assign_batch(first, last)
                        first = last + 1
                    enddo
                    call tabnb%close_assignment
                    call gather_tab3D(tabnb)
                    call tabnb%kill
                    deallocate(tabnb)
                case DEFAULT
                    allocate(tab2D)
                    call tab2D%new_pipeline(self%params, self%build, pinds)
                    call tab2D%set_iteration_model(sigma2_model)
                    call tab2D%open_assignment(nglob=self%nptcls, rank0=ioff)
                    first = 1
                    do while( first <= nloc )
                        last = min(nloc, first + batchsz - 1)
                        call tab2D%fill_tab_range(first, last)
                        call gather_buffers(tab2D%candidate_buffers, self%fill_route == PMB_ROUTE_2D_SPARSE)
                        call tab2D%assign_batch(first, last)
                        first = last + 1
                    enddo
                    call tab2D%close_assignment
                    call gather_tab2D(tab2D)
                    call tab2D%kill
                    deallocate(tab2D)
            end select
            if( l_pert ) call perturb_committed(parts(ipart,1), parts(ipart,2))
        enddo
        if( l_pert ) call self%build%spproj_field%copy(field_saved, .true.)
        call cone%kill
        !$ call omp_set_num_threads(self%nthr_max)

    contains

        !> the batch's dense 3D rows and refined shifts, from the store the fill kept
        subroutine gather_dense_batch3D( tab, ifirst, ilast )
            class(eul_prob_tab), intent(in) :: tab
            integer,             intent(in) :: ifirst, ilast
            integer :: ri, iref, k
            do i = ifirst,ilast
                do ri = 1,tab%nrefs
                    iref = tab%ref_full(ri)
                    res%dists(iref,ioff+i) = tab%batch_rows(ri,i)%dist
                    res%inpls(iref,ioff+i) = tab%batch_rows(ri,i)%inpl
                enddo
                do k = 1,tab%batch_nrefined(i)
                    iref = tab%ref_full(tab%batch_refined(k,i)%ri)
                    res%shifts(:,iref,ioff+i) = [tab%batch_refined(k,i)%x, tab%batch_refined(k,i)%y]
                enddo
                res%nentries = res%nentries + tab%nrefs
            enddo
        end subroutine gather_dense_batch3D

        !> the candidates of a batch as the row, from the thread buffers the fill kept: the sparse routes take every
        !! candidate's shift, as the batch's store holds it (l_all_shifts); the 2D dense route the refined entries'
        subroutine gather_buffers( buffers, l_all_shifts )
            type(prob_candidate_buffer), intent(in) :: buffers(:)
            logical,                     intent(in) :: l_all_shifts
            integer :: ithr, j, iref
            do ithr = 1,size(buffers)
                do j = 1,buffers(ithr)%nused
                    i    = buffers(ithr)%particle_indices(j)
                    iref = buffers(ithr)%candidates(j)%iref
                    if( iref < 1 .or. iref > self%nrefs ) cycle
                    res%dists(iref,ioff+i) = buffers(ithr)%candidates(j)%dist
                    res%inpls(iref,ioff+i) = buffers(ithr)%candidates(j)%inpl
                    if( l_all_shifts .or. buffers(ithr)%candidates(j)%has_sh )&
                        &res%shifts(:,iref,ioff+i) = [buffers(ithr)%candidates(j)%x, buffers(ithr)%candidates(j)%y]
                    res%nentries = res%nentries + 1
                enddo
            enddo
        end subroutine gather_buffers

        !> the temperatures, seeds, assignment and extras of a 3D pipeline table at its close
        subroutine gather_tab3D( tab )
            class(eul_prob_tab), intent(in) :: tab
            type(prob_insert_pose)  :: pose
            type(prob_insert_poses) :: poses
            integer :: k
            do i = 1,nloc
                res%temps(ioff+i)         = tab%temps(i)
                res%seed_shifts(:,ioff+i) = tab%seed_shifts(:,i)
                res%seed_has_sh(ioff+i)   = tab%seed_has_sh(i)
                res%assign_ref(ioff+i)    = (tab%assgn_map(i)%istate - 1) * self%nspace + tab%assgn_map(i)%iproj
                res%assign_inpl(ioff+i)   = tab%assgn_map(i)%inpl
                res%assign_vals(:,ioff+i) = [tab%assgn_map(i)%dist, tab%assgn_map(i)%x, tab%assgn_map(i)%y, tab%assgn_map(i)%sdev]
            enddo
            if( tab%insert_tab%get_nptcls() == nloc )then
                ! the extras as the search reads them, through the part's layout
                call poses%setup_range(tab%insert_tab, pinds, 1, nloc)
                do i = 1,nloc
                    do k = 1,min(poses%nextra(i), size(res%extras,2))
                        pose = tab%insert_tab%get_extra(i, k)
                        res%extras(:,k,ioff+i) = [real(pose%iref), real(pose%inpl), pose%x, pose%y, pose%w]
                    enddo
                enddo
                call poses%kill
            endif
        end subroutine gather_tab3D

        !> the temperatures, seeds and assignment of a 2D pipeline table at its close
        subroutine gather_tab2D( tab )
            type(eul_prob_tab2D), intent(in) :: tab
            do i = 1,nloc
                res%temps(ioff+i)         = tab%temps(i)
                res%seed_shifts(:,ioff+i) = tab%seed_shifts(:,i)
                res%seed_has_sh(ioff+i)   = tab%seed_has_sh(i)
                res%assign_ref(ioff+i)    = tab%assgn_map(i)%icls
                res%assign_inpl(ioff+i)   = tab%assgn_map(i)%inpl
                res%assign_vals(:,ioff+i) = [tab%assgn_map(i)%dist, tab%assgn_map(i)%x, tab%assgn_map(i)%y, tab%assgn_map(i)%sdev]
            enddo
        end subroutine gather_tab2D

        !> what a commit changes in the field of the particles of a part already filled: their distances,
        !! states, classes and shift increments move, and the later parts must not see it
        subroutine perturb_committed( ifirst, ilast )
            integer, intent(in) :: ifirst, ilast
            integer :: j, istate, icls
            do j = ifirst, ilast
                call self%build%spproj_field%set(j, 'dist',     self%build%spproj_field%get(j, 'dist') + 7.5)
                call self%build%spproj_field%set(j, 'shincarg', 3. * self%build%spproj_field%get(j, 'shincarg') + 0.5)
                istate = self%build%spproj_field%get_state(j)
                call self%build%spproj_field%set_state(j, mod(istate, self%nstates) + 1)
                icls = nint(self%build%spproj_field%get(j, 'class'))
                call self%build%spproj_field%set(j, 'class', real(mod(icls, self%nrefs) + 1))
            enddo
        end subroutine perturb_committed

        !> the iteration's model as a part computes it from the complete incoming field: the dense table's
        !! refined-reference budget over the active states
        integer function dense_projs_ns( tab )
            class(eul_prob_tab), intent(in) :: tab
            integer :: si, n
            dense_projs_ns = 0
            do si = 1, tab%nstates
                call calc_num2sample(self%build%spproj_field, self%params%nspace, 'dist', n, self%params%prob_athres, &
                    &state=tab%ssinds(si))
                dense_projs_ns = max(dense_projs_ns, n)
            enddo
        end function dense_projs_ns

        !> and the shift prior's variance per axis
        real function iteration_shift_prior()
            iteration_shift_prior = 0.
            if( self%params%l_shift_prior )&
                &iteration_shift_prior = shift_prior_sigma2(self%build%spproj_field, self%params%lp / self%params%smpd)
        end function iteration_shift_prior

        subroutine alloc_result
            res%nptcls   = self%nptcls
            res%nrefs    = self%nrefs
            res%nentries = 0
            allocate(res%dists(self%nrefs,self%nptcls),    source=huge(1._sp))
            allocate(res%inpls(self%nrefs,self%nptcls),    source=0)
            allocate(res%shifts(2,self%nrefs,self%nptcls), source=huge(1._sp))
            allocate(res%temps(self%nptcls),               source=0._sp)
            allocate(res%seed_shifts(2,self%nptcls),       source=0._sp)
            allocate(res%seed_has_sh(self%nptcls),         source=.false.)
            allocate(res%assign_ref(self%nptcls),          source=0)
            allocate(res%assign_inpl(self%nptcls),         source=0)
            allocate(res%assign_vals(4,self%nptcls),       source=0._sp)
            allocate(res%extras(5,max(1,self%params%topk),self%nptcls), source=0._sp)
        end subroutine alloc_result

    end subroutine run_fill

    subroutine kill_fill_result( res )
        type(pmb_fill_result), intent(inout) :: res
        if( allocated(res%dists) )then
            deallocate(res%dists, res%inpls, res%shifts, res%temps, res%seed_shifts, res%seed_has_sh, res%assign_ref,&
                &res%assign_inpl, res%assign_vals, res%extras)
        endif
        res%nptcls   = 0
        res%nrefs    = 0
        res%nentries = 0
        res%projs_ns = 0
    end subroutine kill_fill_result

    !> the counts of differing entries between two fill results, by particle and reference identity: the rows'
    !! distances (value or presence), their in-plane indices, the refined shifts, the temperatures, the seeds (shift
    !! or flag), the assignments (reference, in-plane index, distance, shift or posterior shift width) and the extras
    subroutine compare_fill_results( a, b, ndiff )
        type(pmb_fill_result), intent(in)  :: a, b
        integer,               intent(out) :: ndiff(PMB_NDIFF)
        integer :: i, iref, k
        ndiff = 0
        if( a%nptcls /= b%nptcls .or. a%nrefs /= b%nrefs ) THROW_HARD('fill results of different sizes; compare_fill_results')
        do i = 1,a%nptcls
            do iref = 1,a%nrefs
                if( a%dists(iref,i) /= b%dists(iref,i) )         ndiff(1) = ndiff(1) + 1
                if( a%inpls(iref,i) /= b%inpls(iref,i) )         ndiff(2) = ndiff(2) + 1
                if( any(a%shifts(:,iref,i) /= b%shifts(:,iref,i)) ) ndiff(3) = ndiff(3) + 1
            enddo
            if( a%temps(i) /= b%temps(i) ) ndiff(4) = ndiff(4) + 1
            if( any(a%seed_shifts(:,i) /= b%seed_shifts(:,i)) .or. (a%seed_has_sh(i) .neqv. b%seed_has_sh(i)) )&
                &ndiff(5) = ndiff(5) + 1
            if( a%assign_ref(i) /= b%assign_ref(i) .or. a%assign_inpl(i) /= b%assign_inpl(i) .or.&
                &any(a%assign_vals(:,i) /= b%assign_vals(:,i)) ) ndiff(6) = ndiff(6) + 1
            do k = 1,min(size(a%extras,2), size(b%extras,2))
                if( any(a%extras(:,k,i) /= b%extras(:,k,i)) ) ndiff(7) = ndiff(7) + 1
            enddo
        enddo
    end subroutine compare_fill_results

    function pmb_route_name( route ) result( name )
        integer, intent(in) :: route
        character(len=:), allocatable :: name
        select case(route)
            case(PMB_ROUTE_3D_DENSE)
                name = '3D_dense'
            case(PMB_ROUTE_3D_SEEDED)
                name = '3D_dense_seeded'
            case(PMB_ROUTE_3D_NOSHIFT)
                name = '3D_dense_noshift'
            case(PMB_ROUTE_3D_NEIGH)
                name = '3D_neigh'
            case(PMB_ROUTE_3D_STATE)
                name = '3D_neigh_state'
            case(PMB_ROUTE_3D_GEOM)
                name = '3D_neigh_geom'
            case(PMB_ROUTE_2D_DENSE)
                name = '2D_dense'
            case(PMB_ROUTE_2D_SPARSE)
                name = '2D_sparse'
            case DEFAULT
                name = 'unknown'
        end select
    end function pmb_route_name

    ! ---- the neighbourhood route (decision 19) ----

    !> The route's case: nspace directions of point group pgrp with nsubs representatives (solve3D's subspace map), nstates
    !! states, nptcls particles for teams of up to nthr_max threads; the candidates of the state mode (the per-pair
    !! routine's entries of the representatives at the particles' temperatures, the peaks and the pool of
    !! fill_tab_subspace_range) and of the geom mode (the cone of athres_cont about the previous direction)
    subroutine new_neigh( self, nspace, nstates, nsubs, pgrp, nptcls, nthr_max )
        use simple_cmdline,            only: cmdline
        use simple_sp_project,         only: sp_project
        use simple_math_ft,            only: calc_fourier_index
        use simple_defs,               only: PI
        use simple_rnd,                only: ran3, rnd_scope_key, RND_PROG_REFINE3D, RND_PASS_ALIGN
        use simple_test_utils,         only: set_fixed_seed
        use simple_srchspace_map,      only: srchspace_map
        use simple_eul_prob_tab_utils, only: cone_cell_radii, cone_reachable_cells, cone_view_dist
        class(polar_neigh_bench), intent(inout) :: self
        integer,                  intent(in)    :: nspace, nstates, nsubs, nptcls, nthr_max
        character(len=*),         intent(in)    :: pgrp
        class(cmdline),    allocatable :: cline
        class(sp_project), allocatable :: spproj
        type(oris)          :: eulspace_sub
        type(srchspace_map) :: sub_mapper
        type(pnb_pass)      :: reps
        complex(sp), allocatable :: anchors(:,:,:,:), band(:,:,:), refb(:,:), ptclb(:,:), ring(:)
        real,        allocatable :: sub_distmat(:,:), ctfmat(:,:)
        logical,     allocatable :: cone_in(:,:)
        integer,     allocatable :: nfill(:)
        real    :: sigma_deg, euls(3), dfx
        integer :: pdim(3), interpklim, i, k, p, s, c, a, ithr, nteam, istate, irot0, q, n, isub
        call self%kill
        if( nspace < 1 .or. nstates < 1 .or. nsubs < 1 .or. nsubs > nspace .or. nptcls < 1 .or. nthr_max < 1 )&
            &THROW_HARD('invalid neighbourhood case')
        self%nspace   = nspace
        self%nstates  = nstates
        self%nrefs    = nspace * nstates
        self%nsubs    = nsubs
        self%nptcls   = nptcls
        self%nthr_max = nthr_max
        allocate(cline, self%params, self%pftc, spproj)
        call cline%set('box',      real(PNB_BOX))
        call cline%set('box_crop', real(PNB_BOX_CROP))
        call cline%set('smpd',     PNB_SMPD)
        call cline%set('mskdiam',  PNB_MSKDIAM)
        call cline%set('lp',       PNB_LP)
        call cline%set('nspace',   real(nspace))
        call cline%set('nstates',  real(nstates))
        call cline%set('nptcls',   real(nptcls))
        call cline%set('nthr',     real(nthr_max))
        call cline%set('ctf',      'yes')
        call cline%set('objfun',   'euclid')
        call cline%set('pgrp',     trim(pgrp))
        call self%params%new(cline, silent=.true.)
        self%kfromto(1) = max(2, calc_fourier_index(self%params%hp, self%params%box, self%params%smpd))
        self%kfromto(2) = calc_fourier_index(self%params%lp, self%params%box, self%params%smpd)
        ! the production defaults of the candidate rules: npeaks peaks per state, the cone of athres_cont
        self%npeak = min(max(1, self%params%npeaks), nsubs)
        self%ncap  = min(nsubs, max(1, self%npeak * max(1, nstates)))
        self%cone  = self%params%athres_cont
        call set_fixed_seed(PMB_SEED)
        call self%pftc%new(self%params, self%nrefs, [1,nptcls], self%kfromto)
        self%pftsz = self%pftc%get_pftsz()
        self%nrots = self%pftc%get_nrots()
        pdim       = self%pftc%get_pdim_interp()
        interpklim = pdim(3)
        ! the direction space and its subspace map, as the builder makes them for prob_neigh's state and geom modes
        call self%pgrpsyms%new(trim(pgrp))
        call self%eulspace%new(nspace, is_ptcl=.false.)
        call self%pgrpsyms%build_refspiral(self%eulspace)
        call eulspace_sub%new(nsubs, is_ptcl=.false.)
        call self%pgrpsyms%build_refspiral(eulspace_sub)
        allocate(sub_distmat(nsubs,nspace))
        call self%pgrpsyms%sym_dist_table(eulspace_sub, self%eulspace, sub_distmat)
        call sub_mapper%new(nspace, nsubs, sub_distmat)
        self%sub_inds = sub_mapper%get_sub2full_map()
        self%full2sub = sub_mapper%get_full2sub_map()
        call sub_mapper%kill
        call eulspace_sub%kill
        ! the cells' members, ascending (sub_ref_list of the neighbourhood table's workspace)
        allocate(self%cell_off(nsubs+1), source=0)
        do p = 1,nspace
            self%cell_off(self%full2sub(p)+1) = self%cell_off(self%full2sub(p)+1) + 1
        enddo
        self%maxcell = maxval(self%cell_off(2:nsubs+1))
        self%cell_off(1) = 1
        do c = 1,nsubs
            self%cell_off(c+1) = self%cell_off(c) + self%cell_off(c+1)
        enddo
        allocate(self%cell_proj(nspace), nfill(nsubs))
        nfill = 0
        do p = 1,nspace
            c = self%full2sub(p)
            self%cell_proj(self%cell_off(c) + nfill(c)) = p
            nfill(c) = nfill(c) + 1
        enddo
        call cone_cell_radii(self%eulspace, self%pgrpsyms, self%sub_inds, self%full2sub, self%cell_radius)
        ! references smooth over the directions: per state, random anchors at the subspace spiral's directions blended
        ! with a Gaussian of the point group's view distance at the cells' spacing (both halves alike)
        sigma_deg = (180. / PI) * sqrt(4. * PI / real(self%pgrpsyms%get_nsym() * nsubs))
        allocate(anchors(self%pftsz,self%kfromto(1):self%kfromto(2),nsubs,nstates))
        allocate(band(self%pftsz,self%kfromto(1):interpklim,nthr_max), source=cmplx(0.,0.,kind=sp))
        do s = 1,nstates
            do a = 1,nsubs
                call fill_anchor(anchors(:,:,a,s))
                if( s > 1 ) anchors(:,:,a,s) = PNB_STATE_CORR * anchors(:,:,a,1) +&
                    &sqrt(1. - PNB_STATE_CORR**2) * anchors(:,:,a,s)
            enddo
        enddo
        nteam = 1
        !$ nteam = omp_get_max_threads()
        if( nteam > nthr_max ) THROW_HARD('the team exceeds the case''s threads; polar_neigh_bench%new')
        !$omp parallel do default(shared) private(p,s,ithr) proc_bind(close) schedule(static)
        do p = 1,nspace
            ithr = 1
            !$ ithr = omp_get_thread_num() + 1
            do s = 1,nstates
                call smooth_band(p, s, band(:,:,ithr))
                call self%pftc%set_ref_pft((s-1)*nspace + p, band(:,:,ithr), iseven=.true.)
                call self%pftc%set_ref_pft((s-1)*nspace + p, band(:,:,ithr), iseven=.false.)
            enddo
        enddo
        !$omp end parallel do
        ! the CTFs through the production path (a one-stack project), then the particles: the previous direction's
        ! reference of the particle's state (round robin), rotated in plane, under the |CTF|, with noise
        call spproj%os_stk%new(1, is_ptcl=.false.)
        call spproj%os_stk%set(1, 'fromp', 1)
        call spproj%os_stk%set(1, 'top',   nptcls)
        call spproj%os_stk%set(1, 'nptcls_stk', nptcls)
        call spproj%os_stk%set(1, 'box',   PNB_BOX)
        call spproj%os_stk%set(1, 'smpd',  PNB_SMPD)
        call spproj%os_stk%set(1, 'ctf',   'yes')
        call spproj%os_stk%set(1, 'kv',    PMB_KV)
        call spproj%os_stk%set(1, 'cs',    PMB_CS)
        call spproj%os_stk%set(1, 'fraca', PMB_FRACA)
        call spproj%os_ptcl3D%new(nptcls, is_ptcl=.true.)
        do i = 1,nptcls
            call spproj%os_ptcl3D%set(i, 'stkind', 1)
            call spproj%os_ptcl3D%set(i, 'indstk', i)
            dfx = 0.8 + 1.7 * ran3()
            call spproj%os_ptcl3D%set(i, 'dfx',    dfx)
            call spproj%os_ptcl3D%set(i, 'dfy',    dfx + 0.05 + 0.1 * ran3())
            call spproj%os_ptcl3D%set(i, 'angast', 180. * ran3())
            call spproj%os_ptcl3D%set_state(i, 1)
        enddo
        call self%pftc%create_polar_absctfmats(spproj, 'ptcl3D')
        allocate(self%prev_proj(nptcls), self%iseven(nptcls), ctfmat(self%pftsz,self%kfromto(1):self%kfromto(2)),&
            &ring(self%nrots))
        allocate(refb(self%pftsz,self%kfromto(1):interpklim), ptclb(self%pftsz,self%kfromto(1):interpklim),&
            &source=cmplx(0.,0.,kind=sp))
        do i = 1,nptcls
            istate = mod(i-1, nstates) + 1
            self%prev_proj(i) = min(nspace, 1 + int(real(nspace) * ran3()))
            self%iseven(i)    = mod(i,2) == 1
            irot0 = min(self%nrots-1, int(real(self%nrots) * ran3()))
            call smooth_band(self%prev_proj(i), istate, refb)
            call self%pftc%get_ptcl_ctfmat(i, ctfmat)
            do k = self%kfromto(1),self%kfromto(2)
                ring(1:self%pftsz)            = refb(:,k)
                ring(self%pftsz+1:self%nrots) = conjg(refb(:,k))
                do p = 1,self%pftsz
                    q = mod(p - 1 + irot0, self%nrots) + 1
                    ptclb(p,k) = ctfmat(p,k) * ring(q) + PNB_NOISE * exp(-real(k) / PMB_FALLOFF) *&
                        &cmplx(2. * ran3() - 1., 2. * ran3() - 1., kind=sp)
                enddo
            enddo
            call self%pftc%set_ptcl_pft(i, ptclb)
            call self%pftc%set_eo(i, self%iseven(i))
        enddo
        ! the noise power of the dense benchmark
        allocate(self%sigma2(self%kfromto(1):self%kfromto(2),2))
        do k = self%kfromto(1),self%kfromto(2)
            self%sigma2(k,1) = 0.05 + exp(-real(k) / 20.)
            self%sigma2(k,2) = 1.02 * self%sigma2(k,1)
        enddo
        call self%pftc%assign_sigma2_noise(self%sigma2)
        call self%pftc%memoize_refs
        call self%pftc%memoize_ptcls
        call self%pftc%memoize_ref_bins
        call self%pftc%set_rnd_scope(rnd_scope_key(int(PMB_SEED, int64), RND_PROG_REFINE3D, RND_PASS_ALIGN, 1))
        allocate(self%temps(nptcls))
        do i = 1,nptcls
            self%temps(i) = self%pftc%ptcl_posterior_temp(i)
        enddo
        ! (a) every state's representatives, state-major in subspace order (find_peak_subspaces)
        allocate(self%rep_refs(nsubs*nstates))
        do s = 1,nstates
            do a = 1,nsubs
                self%rep_refs((s-1)*nsubs + a) = (s-1)*nspace + self%sub_inds(a)
            enddo
        enddo
        ! state mode: the pool of every particle from the per-pair routine's entries of the representatives
        self%exists = .true.
        call self%score_path(PNB_PATH_REPS, .false., nthr_max, 1., reps)
        allocate(self%npooled(nptcls), self%pooled(self%ncap,nptcls))
        call pool_particles(self, reps%dists, self%npooled, self%pooled)
        ! (b) the members of the pooled cells, per state then per pooled cell (evaluate_neighborhood's order)
        allocate(self%mem_off(nptcls+1), self%mem_cell(self%ncap,nstates,nptcls), source=0)
        n = 0
        do i = 1,nptcls
            do s = 1,nstates
                do k = 1,self%npooled(i)
                    c = self%pooled(k,i)
                    n = n + self%cell_off(c+1) - self%cell_off(c)
                enddo
            enddo
        enddo
        allocate(self%mem_ref(n))
        n = 0
        do i = 1,nptcls
            self%mem_off(i) = n + 1
            do s = 1,nstates
                do k = 1,self%npooled(i)
                    c = self%pooled(k,i)
                    self%mem_cell(k,s,i) = n + 1
                    do a = self%cell_off(c),self%cell_off(c+1)-1
                        n = n + 1
                        self%mem_ref(n) = (s-1)*nspace + self%cell_proj(a)
                    enddo
                enddo
            enddo
        enddo
        self%mem_off(nptcls+1) = n + 1
        ! geom mode: the cells the cone reaches and the directions inside it (evaluate_cone: the previous direction,
        ! no mirror view; the reach filter, then the view distance)
        allocate(self%reach(nsubs,nptcls), cone_in(nspace,nptcls))
        !$omp parallel do default(shared) private(i,euls,p,isub) proc_bind(close) schedule(static)
        do i = 1,nptcls
            euls = self%eulspace%get_euler(self%prev_proj(i))
            call cone_reachable_cells(self%eulspace, self%pgrpsyms, self%sub_inds, self%cell_radius, euls, euls,&
                &.false., self%cone, self%reach(:,i))
            do p = 1,nspace
                cone_in(p,i) = .false.
                isub = self%full2sub(p)
                if( isub >= 1 .and. isub <= nsubs )then
                    if( .not. self%reach(isub,i) ) cycle
                endif
                cone_in(p,i) = cone_view_dist(self%pgrpsyms, euls, euls, .false., self%eulspace%get_euler(p)) <= self%cone
            enddo
        enddo
        !$omp end parallel do
        n = count(cone_in) * nstates
        allocate(self%cone_off(nptcls+1), self%cone_ref(n))
        n = 0
        do i = 1,nptcls
            self%cone_off(i) = n + 1
            do s = 1,nstates
                do p = 1,nspace
                    if( .not. cone_in(p,i) ) cycle
                    n = n + 1
                    self%cone_ref(n) = (s-1)*nspace + p
                enddo
            enddo
        enddo
        self%cone_off(nptcls+1) = n + 1
        call kill_pass(reps)
        call spproj%kill
        call cline%kill
        deallocate(spproj, cline, anchors, band, refb, ptclb, ring, sub_distmat, ctfmat, cone_in, nfill)

      contains

        subroutine fill_anchor( pft_band )
            complex(sp), intent(inout) :: pft_band(:,self%kfromto(1):)
            integer :: kk, pp
            do kk = self%kfromto(1),self%kfromto(2)
                do pp = 1,self%pftsz
                    pft_band(pp,kk) = exp(-real(kk) / PMB_FALLOFF) * cmplx(2. * ran3() - 1., 2. * ran3() - 1., kind=sp)
                enddo
            enddo
        end subroutine fill_anchor

        ! direction p's band of state s: the anchors within PNB_SMOOTH_CUT spacings, Gaussian weights normalized in
        ! power (the shells past kfromto(2) zero)
        subroutine smooth_band( p_here, s_here, pft_band )
            integer,     intent(in)    :: p_here, s_here
            complex(sp), intent(inout) :: pft_band(:,self%kfromto(1):)
            real    :: w, wsq
            integer :: a_here
            pft_band = cmplx(0.,0.,kind=sp)
            wsq = 0.
            do a_here = 1,self%nsubs
                if( sub_distmat(a_here,p_here) > PNB_SMOOTH_CUT * sigma_deg ) cycle
                w   = exp(-0.5 * (sub_distmat(a_here,p_here) / sigma_deg)**2)
                wsq = wsq + w * w
                pft_band(:,self%kfromto(1):self%kfromto(2)) = pft_band(:,self%kfromto(1):self%kfromto(2)) +&
                    &w * anchors(:,:,a_here,s_here)
            enddo
            if( wsq > 0. )then
                pft_band = pft_band / sqrt(wsq)
            else
                pft_band(:,self%kfromto(1):self%kfromto(2)) = anchors(:,:,self%full2sub(p_here),s_here)
            endif
        end subroutine smooth_band

    end subroutine new_neigh

    subroutine kill_neigh( self )
        class(polar_neigh_bench), intent(inout) :: self
        if( allocated(self%pftc) )then
            call self%pftc%kill
            deallocate(self%pftc)
        endif
        if( allocated(self%params) ) deallocate(self%params)
        call self%eulspace%kill
        call self%pgrpsyms%kill
        if( allocated(self%sigma2)      ) deallocate(self%sigma2)
        if( allocated(self%temps)       ) deallocate(self%temps)
        if( allocated(self%iseven)      ) deallocate(self%iseven)
        if( allocated(self%prev_proj)   ) deallocate(self%prev_proj)
        if( allocated(self%sub_inds)    ) deallocate(self%sub_inds)
        if( allocated(self%full2sub)    ) deallocate(self%full2sub)
        if( allocated(self%cell_off)    ) deallocate(self%cell_off)
        if( allocated(self%cell_proj)   ) deallocate(self%cell_proj)
        if( allocated(self%cell_radius) ) deallocate(self%cell_radius)
        if( allocated(self%rep_refs)    ) deallocate(self%rep_refs)
        if( allocated(self%npooled)     ) deallocate(self%npooled)
        if( allocated(self%pooled)      ) deallocate(self%pooled)
        if( allocated(self%mem_off)     ) deallocate(self%mem_off)
        if( allocated(self%mem_ref)     ) deallocate(self%mem_ref)
        if( allocated(self%mem_cell)    ) deallocate(self%mem_cell)
        if( allocated(self%reach)       ) deallocate(self%reach)
        if( allocated(self%cone_off)    ) deallocate(self%cone_off)
        if( allocated(self%cone_ref)    ) deallocate(self%cone_ref)
        self%nspace = 0
        self%nptcls = 0
        self%exists = .false.
    end subroutine kill_neigh

    subroutine kill_pass( res )
        type(pnb_pass), intent(inout) :: res
        if( allocated(res%dists) ) deallocate(res%dists)
        if( allocated(res%irots) ) deallocate(res%irots)
        res%wall    = 0._dp
        res%tsec    = 0._dp
        res%nscored = 0_int64
        res%nwasted = 0_int64
        res%ncalls  = 0_int64
        res%nrows   = 0_int64
        res%ncols   = 0_int64
        res%nslots  = 0_int64
        res%bsz     = 0
    end subroutine kill_pass

    !> the pairs a path's particles need: (a) every state's representatives, (b) their pooled cells' members, (c) the
    !! references inside their cones
    integer function npairs( self, path )
        class(polar_neigh_bench), intent(in) :: self
        integer,                  intent(in) :: path
        select case(path)
            case(PNB_PATH_REPS)
                npairs = self%nptcls * size(self%rep_refs)
            case(PNB_PATH_MEMBERS)
                npairs = self%mem_off(self%nptcls+1) - 1
            case(PNB_PATH_CONE)
                npairs = self%cone_off(self%nptcls+1) - 1
            case DEFAULT
                npairs = 0
                THROW_HARD('unknown path; polar_neigh_bench%npairs')
        end select
    end function npairs

    subroutine get_neigh_sizes( self, pftsz, nrots, kfromto, nsubs, ncap, maxcell, cone )
        class(polar_neigh_bench), intent(in)  :: self
        integer,                  intent(out) :: pftsz, nrots, kfromto(2), nsubs, ncap, maxcell
        real,                     intent(out) :: cone
        pftsz   = self%pftsz
        nrots   = self%nrots
        kfromto = self%kfromto
        nsubs   = self%nsubs
        ncap    = self%ncap
        maxcell = self%maxcell
        cone    = self%cone
    end subroutine get_neigh_sizes

    !> the particles whose pool differs when it is formed from the given entries of the representatives (path a's
    !! layout) instead of the per-pair routine's
    integer function pooled_changes( self, rep_dists )
        class(polar_neigh_bench), intent(in) :: self
        real(sp),                 intent(in) :: rep_dists(:)
        integer, allocatable :: npooled(:), pooled(:,:)
        integer :: i
        allocate(npooled(self%nptcls), pooled(self%ncap,self%nptcls))
        call pool_particles(self, rep_dists, npooled, pooled)
        pooled_changes = 0
        do i = 1,self%nptcls
            if( npooled(i) /= self%npooled(i) )then
                pooled_changes = pooled_changes + 1
            else if( any(pooled(1:npooled(i),i) /= self%pooled(1:npooled(i),i)) )then
                pooled_changes = pooled_changes + 1
            endif
        enddo
        deallocate(npooled, pooled)
    end function pooled_changes

    !> state mode's pool per particle from the entries of the representatives (rep_dists in path a's layout):
    !! consider_peak_subspace's sorted top npeak per state, then build_pooled_neighborhood's merge with the cell of
    !! the previous direction (fill_tab_subspace_range, replicated)
    subroutine pool_particles( self, rep_dists, npooled, pooled )
        class(polar_neigh_bench), intent(in)    :: self
        real(sp),                 intent(in)    :: rep_dists(:)
        integer,                  intent(inout) :: npooled(:), pooled(:,:)
        real    :: pk_dists(self%npeak,self%nstates)
        integer :: pk_inds(self%npeak,self%nstates), pk_count(self%nstates)
        integer :: i, s, isub, nrep, nfound, insert, k, prev_isub, ncopy, cnt
        logical :: has_prev
        real    :: dist
        nrep = size(self%rep_refs)
        if( size(rep_dists) /= nrep * self%nptcls ) THROW_HARD('wrong entries of the representatives; pool_particles')
        do i = 1,self%nptcls
            pk_dists = huge(1.0)
            pk_inds  = 0
            pk_count = 0
            do s = 1,self%nstates
                do isub = 1,self%nsubs
                    dist   = rep_dists((i-1)*nrep + (s-1)*self%nsubs + isub)
                    nfound = pk_count(s)
                    if( nfound >= self%npeak )then
                        if( dist >= pk_dists(self%npeak,s) ) cycle
                        insert = self%npeak
                    else
                        nfound      = nfound + 1
                        pk_count(s) = nfound
                        insert      = nfound
                    endif
                    do k = 1,nfound
                        if( dist < pk_dists(k,s) )then
                            insert = k
                            exit
                        endif
                    enddo
                    do k = nfound,insert+1,-1
                        pk_dists(k,s) = pk_dists(k-1,s)
                        pk_inds(k,s)  = pk_inds(k-1,s)
                    enddo
                    pk_dists(insert,s) = dist
                    pk_inds(insert,s)  = isub
                enddo
            enddo
            cnt = 0
            prev_isub = self%full2sub(max(1, min(self%nspace, self%prev_proj(i))))
            if( prev_isub < 1 .or. prev_isub > self%nsubs ) prev_isub = 1
            do s = 1,self%nstates
                has_prev = any(pk_inds(1:pk_count(s),s) == prev_isub)
                ncopy    = pk_count(s)
                if( (.not. has_prev) .and. pk_count(s) >= self%npeak ) ncopy = max(0, self%npeak - 1)
                do k = 1,ncopy
                    call add_pooled(pk_inds(k,s))
                enddo
                if( .not. has_prev ) call add_pooled(prev_isub)
            enddo
            npooled(i) = cnt
        enddo

      contains

        subroutine add_pooled( isub_here )
            integer, intent(in) :: isub_here
            if( isub_here < 1 .or. isub_here > self%nsubs ) return
            if( any(pooled(1:cnt,i) == isub_here) ) return
            if( cnt >= self%ncap ) return
            cnt = cnt + 1
            pooled(cnt,i) = isub_here
        end subroutine add_pooled

    end subroutine pool_particles

    !> decision 18's order (by half, even first, then by the cell of the previous direction, then by index) cut into
    !! decision 3's blocks per half, B = min(PFTC_BLOCK_NB, max(4, n_half / (2 nteam)))
    subroutine neigh_blocks( self, nteam, porder, blocks, bsz )
        class(polar_neigh_bench), intent(in)    :: self
        integer,                  intent(in)    :: nteam
        integer, allocatable,     intent(inout) :: porder(:), blocks(:,:)
        integer,                  intent(out)   :: bsz(2)
        integer :: nhalf(2), ih, c, i, n, nblk, j0, jstart
        allocate(porder(self%nptcls))
        n = 0
        do ih = 1,2
            jstart = n
            do c = 1,self%nsubs
                do i = 1,self%nptcls
                    if( self%iseven(i) .neqv. (ih == 1) ) cycle
                    if( self%full2sub(self%prev_proj(i)) /= c ) cycle
                    n = n + 1
                    porder(n) = i
                enddo
            enddo
            nhalf(ih) = n - jstart
            bsz(ih)   = min(PFTC_BLOCK_NB, max(4, nhalf(ih) / (2 * nteam)))
        enddo
        if( n /= self%nptcls ) THROW_HARD('a particle without a cell; neigh_blocks')
        nblk = (nhalf(1) + bsz(1) - 1) / bsz(1) + (nhalf(2) + bsz(2) - 1) / bsz(2)
        allocate(blocks(2,nblk))
        nblk = 0
        do j0 = 1,nhalf(1),bsz(1)
            nblk = nblk + 1
            blocks(:,nblk) = [j0, min(j0+bsz(1)-1, nhalf(1))]
        enddo
        do j0 = nhalf(1)+1,self%nptcls,bsz(2)
            nblk = nblk + 1
            blocks(:,nblk) = [j0, min(j0+bsz(2)-1, self%nptcls)]
        enddo
    end subroutine neigh_blocks

    !> One pass of a path on a team of nthr threads at the particles' temperatures times tscale: the per-pair routine
    !! (a thread takes a particle and scores its pairs in production's order, gen_posterior_val) or the block kernel
    !! (a thread takes a block of decision 18's order and scores it as the path does); the entries in res (npairs).
    !! With l_cellorder the per-pair routine takes the blocks of decision 18's order instead of single particles, as a
    !! fill in blocks would score a path left on it (reported, not judged)
    subroutine score_path( self, path, l_block, nthr, tscale, res, l_cellorder )
        class(polar_neigh_bench), intent(inout) :: self
        integer,                  intent(in)    :: path, nthr
        logical,                  intent(in)    :: l_block
        real,                     intent(in)    :: tscale
        type(pnb_pass),           intent(out)   :: res
        logical, optional,        intent(in)    :: l_cellorder
        integer,        allocatable :: porder(:), blocks(:,:), wlist(:,:), wsub(:,:), wslot(:,:), wirot(:,:)
        real(sp),       allocatable :: temps(:), wdist(:,:), wtemp(:,:)
        real(dp),       allocatable :: tthr(:)
        integer(int64), allocatable :: cnt(:,:)
        real(dp) :: t0
        integer  :: nteam, ib, i, ithr, nwork, j
        logical  :: l_order
        if( .not. self%exists ) THROW_HARD('polar_neigh_bench: score_path before new')
        if( nthr < 1 .or. nthr > self%nthr_max ) THROW_HARD('polar_neigh_bench: team larger than the calculator')
        if( path < 1 .or. path > PNB_NPATHS ) THROW_HARD('polar_neigh_bench: unknown path')
        allocate(res%dists(self%npairs(path)), source=huge(1._sp))
        allocate(res%irots(self%npairs(path)), source=0)
        l_order = .false.
        if( present(l_cellorder) ) l_order = l_cellorder .and. .not. l_block
        temps = self%temps * tscale
        !$ call omp_set_num_threads(nthr)
        nteam = 1
        !$ nteam = omp_get_max_threads()
        allocate(tthr(nteam), source=0._dp)
        allocate(cnt(6,nteam), source=0_int64)
        t0 = 0._dp
        if( l_block )then
            call neigh_blocks(self, nteam, porder, blocks, res%bsz)
            nwork = max(size(self%rep_refs), self%maxcell * self%nstates)
            allocate(wlist(nwork,nteam), wsub(PFTC_BLOCK_NB,nteam), wslot(PFTC_BLOCK_NB,nteam),&
                &wtemp(PFTC_BLOCK_NB,nteam), wdist(nwork*PFTC_BLOCK_NB,nteam), wirot(nwork*PFTC_BLOCK_NB,nteam))
            !$ res%wall = omp_get_wtime()
            ! dynamic over the blocks; a block is processed whole by one thread in its work columns
            !$omp parallel do default(shared) private(ib,ithr,t0) proc_bind(close) schedule(dynamic)
            do ib = 1,size(blocks,2)
                ithr = 1
                !$ ithr = omp_get_thread_num() + 1
                !$ t0 = omp_get_wtime()
                select case(path)
                    case(PNB_PATH_REPS)
                        call block_reps(blocks(1,ib), blocks(2,ib), ithr)
                    case(PNB_PATH_MEMBERS)
                        call block_cells(blocks(1,ib), blocks(2,ib), ithr, .false.)
                    case(PNB_PATH_CONE)
                        call block_cells(blocks(1,ib), blocks(2,ib), ithr, .true.)
                end select
                !$ tthr(ithr) = tthr(ithr) + (omp_get_wtime() - t0)
            enddo
            !$omp end parallel do
            !$ res%wall = omp_get_wtime() - res%wall
            deallocate(porder, blocks, wlist, wsub, wslot, wtemp, wdist, wirot)
        else if( l_order )then
            call neigh_blocks(self, nteam, porder, blocks, res%bsz)
            !$ res%wall = omp_get_wtime()
            !$omp parallel do default(shared) private(ib,j,ithr,t0) proc_bind(close) schedule(dynamic)
            do ib = 1,size(blocks,2)
                ithr = 1
                !$ ithr = omp_get_thread_num() + 1
                !$ t0 = omp_get_wtime()
                do j = blocks(1,ib),blocks(2,ib)
                    call particle_pairs(porder(j), ithr)
                enddo
                !$ tthr(ithr) = tthr(ithr) + (omp_get_wtime() - t0)
            enddo
            !$omp end parallel do
            !$ res%wall = omp_get_wtime() - res%wall
            deallocate(porder, blocks)
        else
            !$ res%wall = omp_get_wtime()
            ! dynamic over the particles, as the neighbourhood fill takes them
            !$omp parallel do default(shared) private(i,ithr,t0) proc_bind(close) schedule(dynamic)
            do i = 1,self%nptcls
                ithr = 1
                !$ ithr = omp_get_thread_num() + 1
                !$ t0 = omp_get_wtime()
                call particle_pairs(i, ithr)
                !$ tthr(ithr) = tthr(ithr) + (omp_get_wtime() - t0)
            enddo
            !$omp end parallel do
            !$ res%wall = omp_get_wtime() - res%wall
        endif
        res%tsec    = sum(tthr)
        res%nscored = sum(cnt(1,:))
        res%nwasted = sum(cnt(2,:))
        res%ncalls  = sum(cnt(3,:))
        res%nrows   = sum(cnt(4,:))
        res%ncols   = sum(cnt(5,:))
        res%nslots  = sum(cnt(6,:))
        deallocate(temps, tthr, cnt)
        !$ call omp_set_num_threads(self%nthr_max)

      contains

        ! the per-pair routine over particle i's pairs of the path, in production's order
        subroutine particle_pairs( i_here, ithr_here )
            integer, intent(in) :: i_here, ithr_here
            integer :: j, nrep, jfrom, jto
            select case(path)
                case(PNB_PATH_REPS)
                    nrep = size(self%rep_refs)
                    do j = 1,nrep
                        call self%pftc%gen_posterior_val(self%rep_refs(j), i_here, [0._sp,0._sp], temps(i_here),&
                            &res%dists((i_here-1)*nrep + j), res%irots((i_here-1)*nrep + j))
                    enddo
                    cnt(1,ithr_here) = cnt(1,ithr_here) + int(nrep, int64)
                case(PNB_PATH_MEMBERS)
                    jfrom = self%mem_off(i_here)
                    jto   = self%mem_off(i_here+1) - 1
                    do j = jfrom,jto
                        call self%pftc%gen_posterior_val(self%mem_ref(j), i_here, [0._sp,0._sp], temps(i_here),&
                            &res%dists(j), res%irots(j))
                    enddo
                    cnt(1,ithr_here) = cnt(1,ithr_here) + int(jto - jfrom + 1, int64)
                case(PNB_PATH_CONE)
                    jfrom = self%cone_off(i_here)
                    jto   = self%cone_off(i_here+1) - 1
                    do j = jfrom,jto
                        call self%pftc%gen_posterior_val(self%cone_ref(j), i_here, [0._sp,0._sp], temps(i_here),&
                            &res%dists(j), res%irots(j))
                    enddo
                    cnt(1,ithr_here) = cnt(1,ithr_here) + int(jto - jfrom + 1, int64)
            end select
        end subroutine particle_pairs

        ! (a) the block's particles against every state's representatives, one listed tile
        subroutine block_reps( j0, j1, ithr_here )
            integer, intent(in) :: j0, j1, ithr_here
            integer :: nb, nrep, b, r, i_here
            nb   = j1 - j0 + 1
            nrep = size(self%rep_refs)
            do b = 1,nb
                wtemp(b,ithr_here) = temps(porder(j0+b-1))
            enddo
            call self%pftc%score_block_posterior(nb, porder(j0:j1), nrep, self%rep_refs, wtemp(1:nb,ithr_here),&
                &wdist(:,ithr_here), wirot(:,ithr_here))
            do b = 1,nb
                i_here = porder(j0+b-1)
                do r = 1,nrep
                    res%dists((i_here-1)*nrep + r) = wdist((b-1)*nrep + r,ithr_here)
                    res%irots((i_here-1)*nrep + r) = wirot((b-1)*nrep + r,ithr_here)
                enddo
            enddo
            call count_call(ithr_here, nb, nrep, 0)
        end subroutine block_reps

        ! (b) per cell, every state's members against the block's particles that pooled it; (c) against those whose
        ! cone reaches it, the entries inside a particle's cone recorded, the others wasted
        subroutine block_cells( j0, j1, ithr_here, l_cone )
            integer, intent(in) :: j0, j1, ithr_here
            logical, intent(in) :: l_cone
            integer :: c, b, k, nsb, ncm, nr, s, m, i_here, slot, iref, nwaste
            do c = 1,self%nsubs
                ncm = self%cell_off(c+1) - self%cell_off(c)
                if( ncm < 1 ) cycle
                nsb = 0
                do b = j0,j1
                    i_here = porder(b)
                    if( l_cone )then
                        if( .not. self%reach(c,i_here) ) cycle
                        nsb = nsb + 1
                        wsub(nsb,ithr_here) = i_here
                    else
                        do k = 1,self%npooled(i_here)
                            if( self%pooled(k,i_here) == c )then
                                nsb = nsb + 1
                                wsub(nsb,ithr_here)  = i_here
                                wslot(nsb,ithr_here) = k
                                exit
                            endif
                        enddo
                    endif
                enddo
                if( nsb == 0 ) cycle
                nr = ncm * self%nstates
                do s = 1,self%nstates
                    do m = 1,ncm
                        wlist((s-1)*ncm + m,ithr_here) = (s-1)*self%nspace + self%cell_proj(self%cell_off(c) + m - 1)
                    enddo
                enddo
                do b = 1,nsb
                    wtemp(b,ithr_here) = temps(wsub(b,ithr_here))
                enddo
                call self%pftc%score_block_posterior(nsb, wsub(1:nsb,ithr_here), nr, wlist(1:nr,ithr_here),&
                    &wtemp(1:nsb,ithr_here), wdist(:,ithr_here), wirot(:,ithr_here))
                nwaste = 0
                do b = 1,nsb
                    i_here = wsub(b,ithr_here)
                    do s = 1,self%nstates
                        do m = 1,ncm
                            if( l_cone )then
                                iref = wlist((s-1)*ncm + m,ithr_here)
                                slot = cone_slot(i_here, iref)
                                if( slot < 1 )then
                                    nwaste = nwaste + 1
                                    cycle
                                endif
                            else
                                slot = self%mem_cell(wslot(b,ithr_here),s,i_here) + m - 1
                            endif
                            res%dists(slot) = wdist((b-1)*nr + (s-1)*ncm + m,ithr_here)
                            res%irots(slot) = wirot((b-1)*nr + (s-1)*ncm + m,ithr_here)
                        enddo
                    enddo
                enddo
                call count_call(ithr_here, nsb, nr, nwaste)
            enddo
        end subroutine block_cells

        ! where reference iref lies in particle i's cone entries (ascending), 0 outside the cone
        integer function cone_slot( i_here, iref )
            integer, intent(in) :: i_here, iref
            integer :: lo, hi, mid
            cone_slot = 0
            lo = self%cone_off(i_here)
            hi = self%cone_off(i_here+1) - 1
            do while( lo <= hi )
                mid = (lo + hi) / 2
                if( self%cone_ref(mid) == iref )then
                    cone_slot = mid
                    return
                else if( self%cone_ref(mid) < iref )then
                    lo = mid + 1
                else
                    hi = mid - 1
                endif
            enddo
        end function cone_slot

        subroutine count_call( ithr_here, nb, nr, nwaste )
            integer, intent(in) :: ithr_here, nb, nr, nwaste
            cnt(1,ithr_here) = cnt(1,ithr_here) + int(nb, int64) * int(nr, int64)
            cnt(2,ithr_here) = cnt(2,ithr_here) + int(nwaste, int64)
            cnt(3,ithr_here) = cnt(3,ithr_here) + 1_int64
            cnt(4,ithr_here) = cnt(4,ithr_here) + int(nb, int64)
            cnt(5,ithr_here) = cnt(5,ithr_here) + int(nr, int64)
            cnt(6,ithr_here) = cnt(6,ithr_here) + int(PFTC_BLOCK_NB * PFTC_BLOCK_NR, int64) *&
                &int((nr + PFTC_BLOCK_NR - 1) / PFTC_BLOCK_NR, int64)
        end subroutine count_call

    end subroutine score_path

end module simple_polar_matching_bench
