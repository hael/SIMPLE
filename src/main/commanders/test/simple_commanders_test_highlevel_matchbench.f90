!@descr: polar matching layout benchmark of simple_commanders_test_highlevel (measure only, no CTest entry)
submodule (simple_commanders_test_highlevel) simple_commanders_test_highlevel_matchbench
use simple_polar_matching_bench, only: polar_matching_bench, PMB_PER_PAIR, PMB_REPLICA, PMB_BLOCKED, PMB_STATES_INNER,&
    &PMB_MATRIX, PMB_MATRIX_T, PMB_PER_PAIR_CC, PMB_BLOCK_PROD, PMB_NKINDS, pmb_fill_result, PMB_NROUTES, PMB_NDIFF,&
    &PMB_DIFF_NAMES,&
    &pmb_route_name, kill_fill_result, compare_fill_results, PMB_ROUTE_3D_DENSE, PMB_ROUTE_3D_SEEDED,&
    &polar_neigh_bench, pnb_pass, PNB_PATH_REPS, PNB_PATH_MEMBERS, PNB_PATH_CONE, PNB_NPATHS
use simple_polarft_calc,         only: SHIFTMAP_SCREEN, SHIFTMAP_WINDOW
implicit none
#include "simple_local_flags.inc"

character(len=*), parameter :: MATCHBENCH_TSV   = 'polar_matching_bench.tsv'
character(len=*), parameter :: MATCHCHECK_TSV   = 'polar_matching_check.tsv'
character(len=*), parameter :: BLOCKCHECK_TSV   = 'polar_block_check.tsv'
character(len=*), parameter :: IDENTITY_TSV     = 'polar_kernel_identity.tsv'
character(len=*), parameter :: FILLS_TSV        = 'polar_fill_identity.tsv'
character(len=*), parameter :: NEIGH_TSV        = 'polar_neigh_bench.tsv'
character(len=*), parameter :: NEIGHCHECK_TSV   = 'polar_neigh_check.tsv'
character(len=*), parameter :: SHIFTMAP_TSV     = 'polar_shiftmap_bench.tsv'
integer,          parameter :: NSPACE_FILLS     = 200    !< directions of suite fills unless nspace is given
integer,          parameter :: NSTATES_FILLS(2) = [1, 2] !< its states unless nstates is given
integer,          parameter :: NPTCLS_FILLS     = 64     !< its particles unless nptcls is given
integer,          parameter :: NKEEP_IDENTITY   = 32     !< particles of suite identity unless nptcls is given
integer,          parameter :: NREPS            = 3      !< timing: best of three after one warm-up
integer,          parameter :: NREPS_CHECK      = 5      !< identity: repeats at the full team
integer,          parameter :: NPTCLS_ONE_THREAD = 200
integer,          parameter :: BLOCKED_BSZ(4)   = [4, 8, 16, 32]
integer,          parameter :: BLOCKED_RSZ      = 8      !< 8 references of one half: 1.4 MB, within a performance core's share of L2
integer,          parameter :: MATRIX_RSZ(3)    = [32, 64, 128]
integer,          parameter :: MATRIX_BSZ(4)    = [4, 8, 16, 32]
integer,          parameter :: BLOCK_BSZ(3)     = [4, 8, 16]  !< the block kernel's blocks (at most PFTC_BLOCK_NB)
integer,          parameter :: BLOCK_RSZ        = 32          !< its tile (PFTC_BLOCK_NR, fixed in the kernel)
real,             parameter :: BLOCK_FLIP_RATE  = 2.e-5       !< its drawn rotations may differ in max(1, this share) of pairs
real,             parameter :: BLOCK_RTOL       = 1.e-6       !< its drawn losses within this relative difference
real,             parameter :: BLOCK_TSCALES(2) = [1.e25, 0.1] !< the draws' temperature scales of the sensitivity check
integer,          parameter :: BLOCK_BSZ_SENS = 16       !< the block of the sensitivity check
integer,          parameter :: NSPACE_NEIGH     = 5000        !< suite neigh: the late stages' directions unless nspace is given
integer,          parameter :: NSTATES_NEIGH(2) = [1, 2]      !< its states unless nstates is given
character(len=2), parameter :: NEIGH_PGRPS(2)   = ['c1', 'd2']
integer,          parameter :: NEIGH_NREPS      = 2           !< timed repeats per team
integer,          parameter :: NEIGH_TEAM_ID    = 4           !< the team of the identity pass
integer,          parameter :: NSPACE_MAPS      = 2500        !< suite shiftmap: directions unless nspace is given (block_prod's largest case)
integer,          parameter :: NSTATES_MAPS(1)  = [4]         !< its states unless nstates is given
integer,          parameter :: NCAND_MAPS       = 16          !< its candidate references per particle (at most nrefs)
integer,          parameter :: NCHECK_MAPS      = 64          !< its team identity: particles compared at one thread and the full team

contains

module subroutine exec_test_polar_matching_bench( self, cline )
    !$ use omp_lib, only: omp_get_num_procs
    class(commander_test_polar_matching_bench), intent(inout) :: self
    class(cmdline),                             intent(inout) :: cline
    integer, parameter :: NSPACE_ALL(2) = [1000, 2500], NSTATES_ALL(3) = [1, 2, 4]
    integer, allocatable :: nspaces(:), nstates_list(:), nthrs(:), layouts(:)
    type(string) :: suite, fftw_plans
    integer      :: nptcls, nthr_max, ispace, istate
    logical      :: l_check, l_identity, l_fills, l_neigh, l_maps, l_measure, l_planner
    integer      :: ipg
    suite    = ''
    nptcls   = 2000
    nthr_max = 1
    !$ nthr_max = omp_get_num_procs()
    if( cline%defined('suite')  ) suite    = cline%get_carg('suite')
    if( cline%defined('nptcls') ) nptcls   = cline%get_iarg('nptcls')
    if( cline%defined('nthr')   ) nthr_max = cline%get_iarg('nthr')
    nspaces      = NSPACE_ALL
    nstates_list = NSTATES_ALL
    if( cline%defined('nspace')  ) nspaces      = [cline%get_iarg('nspace')]
    if( cline%defined('nstates') ) nstates_list = [cline%get_iarg('nstates')]
    l_check    = .false.
    l_identity = .false.
    l_fills    = .false.
    l_neigh    = .false.
    l_maps     = .false.
    select case(suite%to_char())
        case('list')
            write(logfhandle,'(a)') 'Available suites for polar_matching_bench: per_pair, replica, blocked, states_inner,&
                & matrix, matrix_t, block_prod, per_pair_cc, check, identity, fills, neigh, shiftmap'
            return
        case('')
            layouts = [PMB_PER_PAIR, PMB_REPLICA, PMB_BLOCKED, PMB_STATES_INNER, PMB_MATRIX, PMB_MATRIX_T, PMB_BLOCK_PROD]
        case('per_pair')
            layouts = [PMB_PER_PAIR]
        case('replica')
            layouts = [PMB_REPLICA]
        case('blocked')
            layouts = [PMB_BLOCKED]
        case('states_inner')
            layouts = [PMB_STATES_INNER]
        case('matrix')
            layouts = [PMB_MATRIX]
        case('matrix_t')
            layouts = [PMB_MATRIX_T]
        case('block_prod')
            layouts = [PMB_BLOCK_PROD]
        case('per_pair_cc')
            layouts = [PMB_PER_PAIR_CC]
        case('check')
            l_check = .true.
        case('identity')
            l_identity = .true.
        case('fills')
            l_fills = .true.
            if( .not. cline%defined('nspace')  ) nspaces      = [NSPACE_FILLS]
            if( .not. cline%defined('nstates') ) nstates_list = NSTATES_FILLS
            if( .not. cline%defined('nptcls')  ) nptcls       = NPTCLS_FILLS
        case('neigh')
            ! the matcher's batch at the full team (BATCHTHRSZ particles per thread)
            l_neigh = .true.
            if( .not. cline%defined('nspace')  ) nspaces      = [NSPACE_NEIGH]
            if( .not. cline%defined('nstates') ) nstates_list = NSTATES_NEIGH
            if( .not. cline%defined('nptcls')  ) nptcls       = nthr_max * BATCHTHRSZ
        case('shiftmap')
            ! the dense fill's screen and window maps, by default on block_prod's largest case, so that both
            ! suites run on the same data and reference memo
            l_maps = .true.
            if( .not. cline%defined('nspace')  ) nspaces      = [NSPACE_MAPS]
            if( .not. cline%defined('nstates') ) nstates_list = NSTATES_MAPS
        case default
            THROW_HARD('no suite '//suite%to_char()//' in polar_matching_bench; use suite=list')
    end select
    ! block_prod and shiftmap plan the calculator's FFTs as production does, measured, unless fftw_plans=estimate
    ! asks for the estimated plans, whose output is reproducible for comparing builds; the other suites estimate
    l_planner  = suite%to_char() == 'block_prod' .or. suite%to_char() == 'shiftmap'
    fftw_plans = 'estimate'
    if( l_planner ) fftw_plans = 'measure'
    if( cline%defined('fftw_plans') )then
        if( .not. l_planner ) THROW_HARD('fftw_plans is supported only for block_prod and shiftmap')
        fftw_plans = cline%get_carg('fftw_plans')
    endif
    select case(fftw_plans%to_char())
        case('estimate')
            l_measure = .false.
        case('measure')
            l_measure = .true.
        case default
            THROW_HARD('fftw_plans must be estimate or measure')
    end select
    if( nptcls < 1 .or. nthr_max < 1 .or. any(nspaces < 1) .or. any(nstates_list < 1) )then
        THROW_HARD('polar_matching_bench needs nptcls, nthr, nspace and nstates >= 1')
    endif
    ! the teams, the full one first (the tile sweeps run there): all threads, half of them, one thread
    ! (half of 11 is the Mac's five performance cores, half of 24 the Dell's 12)
    nthrs = [nthr_max]
    if( nthr_max > 3 ) nthrs = [nthrs, nthr_max / 2]
    if( nthr_max > 1 ) nthrs = [nthrs, 1]
    do ispace = 1,size(nspaces)
        do istate = 1,size(nstates_list)
            if( l_check )then
                call check_case(nspaces(ispace), nstates_list(istate), nthr_max)
            else if( l_neigh )then
                do ipg = 1,size(NEIGH_PGRPS)
                    call neigh_case(nspaces(ispace), nstates_list(istate), NEIGH_PGRPS(ipg), nptcls, nthr_max)
                enddo
            else if( l_fills )then
                call fills_case(nspaces(ispace), nstates_list(istate), nptcls, nthr_max)
            else if( l_maps )then
                call shiftmap_case(nspaces(ispace), nstates_list(istate), nptcls, nthrs, l_measure)
            else if( l_identity )then
                if( cline%defined('nptcls') )then
                    call identity_case(nspaces(ispace), nstates_list(istate), nptcls, nthr_max)
                else
                    call identity_case(nspaces(ispace), nstates_list(istate), NKEEP_IDENTITY, nthr_max)
                endif
            else
                call time_case(nspaces(ispace), nstates_list(istate), nptcls, nthrs, layouts, l_measure)
            endif
        enddo
    enddo
    call simple_end('**** SIMPLE_TEST_POLAR_MATCHING_BENCH NORMAL STOP ****')
end subroutine exec_test_polar_matching_bench

!> The timings of one case: per layout and team, best of NREPS after a warm-up; the tile sizes of layouts 2 and 4
!! are swept at the full team and the fastest is timed at the others
subroutine time_case( nspace, nstates, nptcls, nthrs, layouts, l_measure )
    use simple_rnd, only: has_global_seed
    integer, intent(in) :: nspace, nstates, nptcls, nthrs(:), layouts(:)
    logical, intent(in) :: l_measure
    class(polar_matching_bench), allocatable, target :: bench
    character(len=:), allocatable :: plans
    integer, allocatable :: bszs(:), rszs(:)
    real(dp) :: wall, wall_best
    integer  :: ilay, it, nt, best_b, best_r, ib, ir
    logical  :: l_prod, l_memo, l_gmat, l_cc
    ! the planner of the calculator's FFTs, for the TSV
    plans = 'measure'
    if( has_global_seed() .and. .not. l_measure ) plans = 'estimate'
    l_prod = any(layouts == PMB_PER_PAIR) .or. any(layouts == PMB_PER_PAIR_CC) .or. any(layouts == PMB_BLOCK_PROD)
    l_cc   = any(layouts == PMB_PER_PAIR_CC)
    if( l_cc .and. size(layouts) > 1 ) THROW_HARD('per_pair_cc is timed on its own')
    l_memo = any(layouts == PMB_REPLICA) .or. any(layouts == PMB_BLOCKED) .or. any(layouts == PMB_STATES_INNER)
    l_gmat = any(layouts == PMB_MATRIX) .or. any(layouts == PMB_MATRIX_T)
    allocate(bench)
    if( l_cc )then
        call bench%new(nspace, nstates, nptcls, maxval(nthrs), l_prod, l_memo, l_gmat, 0, objfun='cc', fftw_measure=l_measure)
    else
        call bench%new(nspace, nstates, nptcls, maxval(nthrs), l_prod, l_memo, l_gmat, 0, fftw_measure=l_measure)
    endif
    call report_case(bench, nspace, nstates)
    do ilay = 1,size(layouts)
        select case(layouts(ilay))
            case(PMB_BLOCKED)
                bszs = BLOCKED_BSZ
                rszs = [BLOCKED_RSZ]
            case(PMB_MATRIX, PMB_MATRIX_T)
                bszs = MATRIX_BSZ
                rszs = MATRIX_RSZ
            case(PMB_BLOCK_PROD)
                bszs = BLOCK_BSZ
                rszs = [BLOCK_RSZ]
            case DEFAULT
                bszs = [1]
                rszs = [1]
        end select
        best_b = bszs(1)
        best_r = rszs(1)
        do it = 1,size(nthrs)
            nt = nthrs(it)
            if( it == 1 )then
                wall_best = huge(wall_best)
                do ir = 1,size(rszs)
                    do ib = 1,size(bszs)
                        wall = time_one(bench, layouts(ilay), nspace, nstates, nptcls, nt, bszs(ib), rszs(ir), plans)
                        if( wall < wall_best )then
                            wall_best = wall
                            best_b    = bszs(ib)
                            best_r    = rszs(ir)
                        endif
                    enddo
                enddo
            else
                wall = time_one(bench, layouts(ilay), nspace, nstates, nptcls, nt, best_b, best_r, plans)
            endif
        enddo
    enddo
    call bench%kill
    deallocate(bench)
end subroutine time_case

!> one timing row: warm-up, best of NREPS, appended to the TSV with the calculator's FFTW planner; returns the best
!! wall time
function time_one( bench, layout, nspace, nstates, nptcls, nthr, bsz, rsz, plans ) result( wall_best )
    class(polar_matching_bench), intent(inout) :: bench
    integer,                     intent(in)    :: layout, nspace, nstates, nptcls, nthr, bsz, rsz
    character(len=*),            intent(in)    :: plans
    real(dp) :: wall_best, walls(NREPS), pairs, us_pair, pairs_per_s, ref_gbps, parts(3)
    real     :: load1
    integer  :: nptcls_use, nwarm, irep, funit, pftsz, nrots, kfromto(2), nrefs
    logical  :: l_new
    call bench%get_sizes(pftsz, nrots, kfromto, nrefs)
    nptcls_use = nptcls
    if( nthr == 1 ) nptcls_use = min(nptcls, NPTCLS_ONE_THREAD)
    nwarm = min(nptcls_use, max(2*nthr, nptcls_use/10))
    load1 = load_average()
    walls(1) = bench%score(layout, nthr, nwarm, bsz, rsz)
    do irep = 1,NREPS
        walls(irep) = bench%score(layout, nthr, nptcls_use, bsz, rsz)
    enddo
    parts       = bench%get_parts() * 1.e6_dp / (real(nptcls_use,dp) * real(nrefs,dp))
    wall_best   = minval(walls)
    pairs       = real(nptcls_use,dp) * real(nrefs,dp)
    us_pair     = wall_best * real(nthr,dp) * 1.e6_dp / pairs
    pairs_per_s = pairs / wall_best
    ref_gbps    = pairs_per_s * bench%ref_bytes_per_pair() / 1.e9_dp
    inquire(file=MATCHBENCH_TSV, exist=l_new)
    l_new = .not. l_new
    open(newunit=funit, file=MATCHBENCH_TSV, position='append', action='write', status='unknown')
    if( l_new ) write(funit,'(a)') 'layout'//achar(9)//'nspace'//achar(9)//'nstates'//achar(9)//'nrefs'//achar(9)//&
        &'nptcls'//achar(9)//'nthr'//achar(9)//'B'//achar(9)//'R'//achar(9)//'pftsz'//achar(9)//'nrots'//achar(9)//&
        &'kfrom'//achar(9)//'kto'//achar(9)//'load1'//achar(9)//'wall_s'//achar(9)//&
        &'wall_rep1'//achar(9)//'wall_rep2'//achar(9)//'wall_rep3'//achar(9)//'us_pair_thread'//achar(9)//&
        &'pairs_per_s'//achar(9)//'ref_GBps'//achar(9)//'memo_GB'//achar(9)//'check'//achar(9)//&
        &'us_sums'//achar(9)//'us_fft'//achar(9)//'us_draw'//achar(9)//'fftw_plans'
    write(funit,'(a,11(a,i0),a,f6.2,4(a,f10.4),a,f9.3,a,es12.5,a,f8.3,a,f8.3,a,a,3(a,f9.3),a,a)') layout_name(layout),&
        &achar(9), nspace, achar(9), nstates, achar(9), nrefs, achar(9), nptcls_use, achar(9), nthr, achar(9), bsz,&
        &achar(9), rsz, achar(9), pftsz, achar(9), nrots, achar(9), kfromto(1), achar(9), kfromto(2),&
        &achar(9), load1, achar(9), wall_best, achar(9), walls(1), achar(9), walls(2), achar(9), walls(3), achar(9),&
        &us_pair, achar(9), pairs_per_s, achar(9), ref_gbps, achar(9), bench%memo_gbytes(), achar(9), '-',&
        &achar(9), parts(1), achar(9), parts(2), achar(9), parts(3), achar(9), plans
    close(funit)
    write(logfhandle,'(4x,a,a,i0,a,i0,a,i0,a,i0,a,f6.2,a,f9.3,a,f8.2,a,f8.2,a)') layout_name(layout), ' nthr=', nthr,&
        &' B=', bsz, ' R=', rsz, ' nptcls=', nptcls_use, ' load=', load1, ' wall=', wall_best, ' s, ', us_pair,&
        &' us/pair/thread, ', ref_gbps, ' GB/s'
end function time_one

!> The shift maps of one case: the window and the screen map, per team, best of NREPS after a warm-up; then each
!! map's results on NCHECK_MAPS particles at one thread, held to the full team's and, under estimated plans, written
!! for comparing builds
subroutine shiftmap_case( nspace, nstates, nptcls, nthrs, l_measure )
    integer, intent(in) :: nspace, nstates, nptcls, nthrs(:)
    logical, intent(in) :: l_measure
    class(polar_matching_bench), allocatable, target :: bench
    integer,          parameter :: MAPS(2) = [SHIFTMAP_WINDOW, SHIFTMAP_SCREEN]
    character(len=6), parameter :: MAPNAMES(2) = ['window', 'screen']
    real(dp), allocatable :: res_one(:,:,:), res_team(:,:,:)
    real(dp) :: wall
    integer  :: imap, it, pftsz, nrots, kfromto(2), nrefs, ncand, nchk, ndiff, funit
    allocate(bench)
    call bench%new(nspace, nstates, nptcls, maxval(nthrs), .true., .false., .false., 0, fftw_measure=l_measure)
    call report_case(bench, nspace, nstates)
    call bench%get_sizes(pftsz, nrots, kfromto, nrefs)
    ncand = min(NCAND_MAPS, nrefs)
    nchk  = min(NCHECK_MAPS, nptcls)
    do imap = 1,size(MAPS)
        do it = 1,size(nthrs)
            call time_maps(bench, MAPS(imap), nspace, nstates, nptcls, nthrs(it), l_measure)
        enddo
        ! the one-thread results (loss, width, rotation, radius, direction per map) for comparing builds byte for
        ! byte; measured plans may differ between runs in the last bits, so only estimated plans write the file
        wall = bench%score_maps(MAPS(imap), 1, nchk, ncand)
        call bench%get_map_results(res_one)
        if( l_measure )then
            write(logfhandle,'(4x,a,a,a)') 'shiftmap ', MAPNAMES(imap),&
                &' identity file not written under measured plans; fftw_plans=estimate writes it'
        else
            open(newunit=funit, file='polar_shiftmap_identity_'//MAPNAMES(imap)//'.bin', access='stream',&
                &form='unformatted', status='replace', action='write')
            write(funit) res_one
            close(funit)
        endif
        if( maxval(nthrs) > 1 )then
            wall = bench%score_maps(MAPS(imap), maxval(nthrs), nchk, ncand)
            call bench%get_map_results(res_team)
            ndiff = count(any(abs(res_one - res_team) > 0._dp, dim=1))
            write(logfhandle,'(4x,a,a,a,i0,a,i0,a,i0,a,i0,a)') 'shiftmap ', MAPNAMES(imap), ' team identity: ', nchk,&
                &' particles x ', ncand, ' candidates, 1 against ', maxval(nthrs), ' threads: ', ndiff, ' maps differ'
        endif
    enddo
    call bench%kill
    deallocate(bench)
end subroutine shiftmap_case

!> one shift-map row: warm-up, best of NREPS, appended to the TSV with the map's geometry and the last repeat's
!! thread microseconds per map of construction, cell selection and posterior reduction
subroutine time_maps( bench, imap, nspace, nstates, nptcls, nthr, l_measure )
    use simple_rnd, only: has_global_seed
    class(polar_matching_bench), intent(inout) :: bench
    integer,                     intent(in)    :: imap, nspace, nstates, nptcls, nthr
    logical,                     intent(in)    :: l_measure
    character(len=:), allocatable :: mapname, plans
    real(dp) :: walls(NREPS), wall_best, nmaps, us_map, maps_per_s, parts(3), rho_max, prior_w
    real     :: load1
    integer  :: nptcls_use, nwarm, ncand, irep, funit, pftsz, nrots, kfromto(2), nrefs, geom(4)
    logical  :: l_new
    call bench%get_sizes(pftsz, nrots, kfromto, nrefs)
    mapname = 'window'
    if( imap == SHIFTMAP_SCREEN ) mapname = 'screen'
    ! report the planner used for this run
    plans = 'measure'
    if( has_global_seed() .and. .not. l_measure ) plans = 'estimate'
    ncand      = min(NCAND_MAPS, nrefs)
    nptcls_use = nptcls
    if( nthr == 1 ) nptcls_use = min(nptcls, NPTCLS_ONE_THREAD)
    nwarm = min(nptcls_use, max(2*nthr, nptcls_use/10))
    load1 = load_average()
    walls(1) = bench%score_maps(imap, nthr, nwarm, ncand)
    do irep = 1,NREPS
        walls(irep) = bench%score_maps(imap, nthr, nptcls_use, ncand)
    enddo
    nmaps      = real(nptcls_use,dp) * real(ncand,dp)
    parts      = bench%get_map_parts() * 1.e6_dp / nmaps
    wall_best  = minval(walls)
    us_map     = wall_best * real(nthr,dp) * 1.e6_dp / nmaps
    maps_per_s = nmaps / wall_best
    call bench%get_map_geometry(imap, geom, rho_max, prior_w)
    inquire(file=SHIFTMAP_TSV, exist=l_new)
    l_new = .not. l_new
    open(newunit=funit, file=SHIFTMAP_TSV, position='append', action='write', status='unknown')
    if( l_new ) write(funit,'(a)') 'map'//achar(9)//'nspace'//achar(9)//'nstates'//achar(9)//'nrefs'//achar(9)//&
        &'nptcls'//achar(9)//'ncand'//achar(9)//'nthr'//achar(9)//'pftsz'//achar(9)//'nrots'//achar(9)//&
        &'kfrom'//achar(9)//'kscreen'//achar(9)//'nrho'//achar(9)//'nalpha'//achar(9)//'norder'//achar(9)//&
        &'rho_max_px'//achar(9)//'prior_w'//achar(9)//'load1'//achar(9)//'wall_s'//achar(9)//&
        &'wall_rep1'//achar(9)//'wall_rep2'//achar(9)//'wall_rep3'//achar(9)//'us_map_thread'//achar(9)//&
        &'maps_per_s'//achar(9)//'us_construct'//achar(9)//'us_select'//achar(9)//'us_reduce'//achar(9)//'fftw_plans'
    write(funit,'(a,13(a,i0),a,f10.4,a,es12.5,a,f6.2,4(a,f10.4),a,f10.3,a,es12.5,3(a,f10.3),a,a)') mapname,&
        &achar(9), nspace, achar(9), nstates, achar(9), nrefs, achar(9), nptcls_use,&
        &achar(9), ncand, achar(9), nthr, achar(9), pftsz, achar(9), nrots, achar(9), kfromto(1), achar(9), geom(4),&
        &achar(9), geom(1), achar(9), geom(2), achar(9), geom(3), achar(9), rho_max, achar(9), prior_w,&
        &achar(9), load1, achar(9), wall_best, achar(9), walls(1), achar(9), walls(2), achar(9), walls(3),&
        &achar(9), us_map, achar(9), maps_per_s, achar(9), parts(1), achar(9), parts(2), achar(9), parts(3),&
        &achar(9), plans
    close(funit)
    write(logfhandle,'(4x,a,a,a,i0,a,i0,a,i0,a,f6.2,a,f9.3,a,f10.2,a,es12.5,a)') 'shiftmap ', mapname, ' nthr=', nthr,&
        &' nptcls=', nptcls_use, ' ncand=', ncand, ' load=', load1, ' wall=', wall_best, ' s, ', us_map,&
        &' us/map/thread, ', maps_per_s, ' maps/s'
end subroutine time_maps

!> The identity checks of one case on nkeep particles, on the full rows over the rotations (the losses the posterior
!! draws run over): production (layout 1) at one thread is the reference every layout is held to byte for byte, the
!! replica at one thread a second reference for layouts 2-4, and each layout's own one-thread run its thread identity.
!! Production's block kernel is held to production's per-pair draws (rows of polar_block_check.tsv).
subroutine check_case( nspace, nstates, nthr_max )
    integer, intent(in) :: nspace, nstates, nthr_max
    class(polar_matching_bench), allocatable, target :: bench
    real(sp), allocatable :: vecs_prod(:,:,:), vecs_ref(:,:,:), vecs_one(:,:,:), dists_pp(:,:), dists_first(:,:), temps(:)
    integer,  allocatable :: irots_pp(:,:), irots_first(:,:)
    real(dp) :: wall_pp
    integer  :: nkeep, irep, ib, itemp
    nkeep = 64
    if( nspace * nstates > 5000 ) nkeep = 32
    allocate(bench)
    ! production
    call bench%new(nspace, nstates, nkeep, nthr_max, .true., .false., .false., nkeep)
    call report_case(bench, nspace, nstates)
    call check_run(PMB_PER_PAIR, 1, 0, 1, 1, .true.)
    do irep = 1,NREPS_CHECK
        call check_run(PMB_PER_PAIR, nthr_max, irep, 1, 1, .false.)
    enddo
    ! production's block kernel against production's per-pair draws (the same keyed variates), its thread identity
    ! and its independence of the block size, on the same calculator
    wall_pp = bench%score(PMB_PER_PAIR, 1, nkeep, 1, 1)
    call bench%get_entries(dists_pp, irots_pp)
    temps = bench%get_temps()
    write(logfhandle,'(a,3es11.3)') '    posterior temperatures min/median/max: ', minval(temps), median_nocopy(temps),&
        &maxval(temps)
    do ib = 1,size(BLOCK_BSZ)
        call block_run(1, 0, BLOCK_BSZ(ib), 1.)
        do irep = 1,NREPS_CHECK
            call block_run(nthr_max, irep, BLOCK_BSZ(ib), 1.)
        enddo
    enddo
    ! the draws' sensitivity to the last bits of the losses: the argmin (a huge temperature) and a ten times cooler draw
    do itemp = 1,size(BLOCK_TSCALES)
        call bench%scale_temps(BLOCK_TSCALES(itemp))
        wall_pp = bench%score(PMB_PER_PAIR, 1, nkeep, 1, 1)
        call bench%get_entries(dists_pp, irots_pp)
        deallocate(dists_first, irots_first)
        call block_run(1, 0, BLOCK_BSZ_SENS, BLOCK_TSCALES(itemp))
        call block_run(nthr_max, 1, BLOCK_BSZ_SENS, BLOCK_TSCALES(itemp))
    enddo
    call bench%scale_temps(1.)
    call bench%kill
    ! the same data in the benchmark memo
    call bench%new(nspace, nstates, nkeep, nthr_max, .false., .true., .true., nkeep)
    write(logfhandle,'(a,f8.3,a)') '    bin-major copy built in ', bench%get_gmat_build_time(), ' s'
    call check_run(PMB_REPLICA, 1, 0, 1, 1, .true.)
    do irep = 1,NREPS_CHECK
        call check_run(PMB_REPLICA, nthr_max, irep, 1, 1, .false.)
    enddo
    do ib = 1,size(BLOCKED_BSZ)
        call check_run(PMB_BLOCKED, 1, 0, BLOCKED_BSZ(ib), BLOCKED_RSZ, .false.)
        do irep = 1,NREPS_CHECK
            call check_run(PMB_BLOCKED, nthr_max, irep, BLOCKED_BSZ(ib), BLOCKED_RSZ, .false.)
        enddo
    enddo
    call check_run(PMB_STATES_INNER, 1, 0, 1, 1, .false.)
    do irep = 1,NREPS_CHECK
        call check_run(PMB_STATES_INNER, nthr_max, irep, 1, 1, .false.)
    enddo
    do ib = 1,size(MATRIX_BSZ)
        call check_run(PMB_MATRIX, 1, 0, MATRIX_BSZ(ib), 64, .false.)
        do irep = 1,NREPS_CHECK
            call check_run(PMB_MATRIX, nthr_max, irep, MATRIX_BSZ(ib), 64, .false.)
        enddo
        call check_run(PMB_MATRIX_T, 1, 0, MATRIX_BSZ(ib), 64, .false.)
        do irep = 1,NREPS_CHECK
            call check_run(PMB_MATRIX_T, nthr_max, irep, MATRIX_BSZ(ib), 64, .false.)
        enddo
    enddo
    call bench%kill
    deallocate(bench)

  contains

    ! one scoring of the kept particles, compared with the references (or made the reference)
    subroutine check_run( layout, nthr, irep, bsz, rsz, l_make_ref )
        integer, intent(in) :: layout, nthr, irep, bsz, rsz
        logical, intent(in) :: l_make_ref
        real(sp), allocatable :: vecs(:,:,:)
        real(dp) :: wall, max_abs, max_rel, d, max_prod_abs
        integer  :: nrefs, p, iref, i, n_vec_diff, n_argmin_diff, n_prod_diff, n_prod_ptcls, n_vec_diff_one, funit
        logical  :: l_ptcl_diff
        logical  :: l_new
        wall = bench%score(layout, nthr, nkeep, bsz, rsz)
        call bench%get_vecs(vecs)
        nrefs = size(vecs,2)
        if( l_make_ref )then
            if( layout == PMB_PER_PAIR )then
                vecs_prod = vecs
            else
                vecs_ref = vecs
            endif
        endif
        ! the rows against production's, byte for byte
        n_prod_diff  = 0
        n_prod_ptcls = 0
        max_prod_abs = 0._dp
        do i = 1,nkeep
            l_ptcl_diff = .false.
            do iref = 1,nrefs
                if( any(vecs(:,iref,i) /= vecs_prod(:,iref,i)) )then
                    n_prod_diff = n_prod_diff + 1
                    l_ptcl_diff = .true.
                endif
                max_prod_abs = max(max_prod_abs, maxval(abs(real(vecs(:,iref,i),dp) - real(vecs_prod(:,iref,i),dp))))
            enddo
            if( l_ptcl_diff ) n_prod_ptcls = n_prod_ptcls + 1
        enddo
        ! the rows against this layout's own one-thread run (irep 0), thread identity
        n_vec_diff_one = -1
        if( irep == 0 )then
            vecs_one = vecs
        else
            n_vec_diff_one = 0
            do i = 1,nkeep
                do iref = 1,nrefs
                    if( any(vecs(:,iref,i) /= vecs_one(:,iref,i)) ) n_vec_diff_one = n_vec_diff_one + 1
                enddo
            enddo
        endif
        ! the rows against the replica's (layouts 2-4)
        n_vec_diff    = -1
        n_argmin_diff = -1
        max_abs       = -1._dp
        max_rel       = -1._dp
        if( layout /= PMB_PER_PAIR .and. allocated(vecs_ref) )then
            n_vec_diff    = 0
            n_argmin_diff = 0
            max_abs       = 0._dp
            max_rel       = 0._dp
            do i = 1,nkeep
                do iref = 1,nrefs
                    if( any(vecs(:,iref,i) /= vecs_ref(:,iref,i)) ) n_vec_diff = n_vec_diff + 1
                    if( minloc(vecs(:,iref,i),dim=1) /= minloc(vecs_ref(:,iref,i),dim=1) ) n_argmin_diff = n_argmin_diff + 1
                    do p = 1,size(vecs,1)
                        d       = abs(real(vecs(p,iref,i),dp) - real(vecs_ref(p,iref,i),dp))
                        max_abs = max(max_abs, d)
                        max_rel = max(max_rel, d / max(abs(real(vecs_ref(p,iref,i),dp)), 1.e-30_dp))
                    enddo
                enddo
            enddo
        endif
        inquire(file=MATCHCHECK_TSV, exist=l_new)
        l_new = .not. l_new
        open(newunit=funit, file=MATCHCHECK_TSV, position='append', action='write', status='unknown')
        if( l_new ) write(funit,'(a)') 'layout'//achar(9)//'nspace'//achar(9)//'nstates'//achar(9)//'nrefs'//achar(9)//&
            &'nptcls'//achar(9)//'nthr'//achar(9)//'rep'//achar(9)//'B'//achar(9)//'R'//achar(9)//'pairs'//achar(9)//&
            &'vec_diff_vs_prod'//achar(9)//'vec_diff_vs_replica'//achar(9)//'argmin_diff_vs_replica'//achar(9)//&
            &'max_abs_diff'//achar(9)//'max_rel_diff'//achar(9)//'wall_s'//achar(9)//'vec_diff_ptcls_vs_prod'//achar(9)//&
            &'max_abs_diff_vs_prod'//achar(9)//'vec_diff_vs_own_1thr'
        write(funit,'(a,12(a,i0),2(a,es12.5),a,f10.4,a,i0,a,es12.5,a,i0)') layout_name(layout), achar(9), nspace, achar(9), nstates,&
            &achar(9), nrefs, achar(9), nkeep, achar(9), nthr, achar(9), irep, achar(9), bsz, achar(9), rsz, achar(9),&
            &nrefs*nkeep, achar(9), n_prod_diff, achar(9), n_vec_diff, achar(9), n_argmin_diff, achar(9), max_abs,&
            &achar(9), max_rel, achar(9), wall, achar(9), n_prod_ptcls, achar(9), max_prod_abs, achar(9), n_vec_diff_one
        close(funit)
        write(logfhandle,'(4x,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,es10.3)') layout_name(layout), ' nthr=', nthr,&
            &' rep=', irep, ' B=', bsz, ' R=', rsz, ' vec_diff_vs_prod=', n_prod_diff, ' vec_diff=', n_vec_diff,&
            &' argmin_diff=', n_argmin_diff, ' max_rel=', max_rel
    end subroutine check_run

    ! one scoring of the kept particles by production's block kernel against the per-pair draws: rotations differing,
    ! the largest relative difference of the drawn loss where the rotation agrees, and the entries against the first
    ! one-thread block run (thread identity; independence of the block size)
    subroutine block_run( nthr, irep, bsz, tscale )
        integer, intent(in) :: nthr, irep, bsz
        real,    intent(in) :: tscale
        real(sp), allocatable :: dists(:,:)
        integer,  allocatable :: irots(:,:)
        real(dp) :: wall, max_rel, d
        integer  :: nrefs, iref, i, nrot_diff, ninvalid, ndiff_first, nallowed, funit
        logical  :: l_new, l_ok
        wall = bench%score(PMB_BLOCK_PROD, nthr, nkeep, bsz, BLOCK_RSZ)
        call bench%get_entries(dists, irots)
        nrefs     = size(dists,1)
        nrot_diff = 0
        ninvalid  = 0
        max_rel   = 0._dp
        do i = 1,nkeep
            do iref = 1,nrefs
                if( irots(iref,i) /= irots_pp(iref,i) )then
                    nrot_diff = nrot_diff + 1
                else if( dists_pp(iref,i) > 0.5 * huge(1.) )then
                    if( .not. (dists(iref,i) > 0.5 * huge(1.)) ) ninvalid = ninvalid + 1
                else
                    d       = abs(real(dists(iref,i),dp) - real(dists_pp(iref,i),dp))
                    max_rel = max(max_rel, d / max(abs(real(dists_pp(iref,i),dp)), 1.e-30_dp))
                endif
            enddo
        enddo
        ndiff_first = 0
        if( .not. allocated(dists_first) )then
            dists_first = dists
            irots_first = irots
        else
            ndiff_first = count(dists /= dists_first) + count(irots /= irots_first)
        endif
        nallowed = max(1, int(BLOCK_FLIP_RATE * real(nrefs*nkeep)))
        l_ok     = nrot_diff <= nallowed .and. max_rel <= real(BLOCK_RTOL,dp) .and. ninvalid == 0 .and. ndiff_first == 0
        inquire(file=BLOCKCHECK_TSV, exist=l_new)
        l_new = .not. l_new
        open(newunit=funit, file=BLOCKCHECK_TSV, position='append', action='write', status='unknown')
        if( l_new ) write(funit,'(a)') 'layout'//achar(9)//'nspace'//achar(9)//'nstates'//achar(9)//'nrefs'//achar(9)//&
            &'nptcls'//achar(9)//'nthr'//achar(9)//'rep'//achar(9)//'B'//achar(9)//'pairs'//achar(9)//&
            &'rot_diff_vs_per_pair'//achar(9)//'rot_diff_allowed'//achar(9)//'max_rel_diff_vs_per_pair'//achar(9)//&
            &'invalid_mismatch'//achar(9)//'entries_diff_vs_first_1thr'//achar(9)//'wall_s'//achar(9)//'within_tolerance'//&
            &achar(9)//'temp_scale'
        write(funit,'(a,10(a,i0),a,es12.5,2(a,i0),a,f10.4,a,a,a,es10.3)') layout_name(PMB_BLOCK_PROD), achar(9), nspace,&
            &achar(9), nstates, achar(9), nrefs, achar(9), nkeep, achar(9), nthr, achar(9), irep, achar(9), bsz, achar(9),&
            &nrefs*nkeep, achar(9), nrot_diff, achar(9), nallowed, achar(9), max_rel, achar(9), ninvalid, achar(9),&
            &ndiff_first, achar(9), wall, achar(9), trim(merge('yes','no ',l_ok)), achar(9), tscale
        close(funit)
        write(logfhandle,'(4x,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,es10.3,a,i0,a,es10.3)') layout_name(PMB_BLOCK_PROD), ' nthr=', nthr,&
            &' rep=', irep, ' B=', bsz, ' rot_diff_vs_per_pair=', nrot_diff, ' of ', nrefs*nkeep, ' max_rel=', max_rel,&
            &' entries_diff_vs_first_1thr=', ndiff_first, ' temp_scale=', tscale
    end subroutine block_run

end subroutine check_case

!> The production kernels of one case on nkeep particles, euclid and cc (prod_kernels): the one-thread result is
!! written to polar_kernel_identity_<objfun>_<nspace>x<nstates>.bin for a comparison between builds, and NREPS_CHECK
!! runs at the full team are compared with it byte for byte (rows of polar_kernel_identity.tsv)
subroutine identity_case( nspace, nstates, nkeep, nthr_max )
    !$ use omp_lib, only: omp_get_wtime
    integer, intent(in) :: nspace, nstates, nkeep, nthr_max
    character(len=*), parameter :: OBJFUNS(2) = ['euclid', 'cc    ']
    class(polar_matching_bench), allocatable, target :: bench
    real(sp), allocatable :: v1(:,:,:,:), bd1(:,:), v(:,:,:,:), bd(:,:)
    integer,  allocatable :: bi1(:,:), bi(:,:)
    character(len=:), allocatable :: fname
    character(len=32) :: casestr
    real(dp) :: wall
    integer  :: iobj, irep, funit, ikind, iref, i, ndiff(PMB_NKINDS), ndiff_best
    integer  :: pftsz, nrots, kfromto(2), nrefs
    logical  :: l_new
    allocate(bench)
    write(casestr,'(i0,a,i0)') nspace, 'x', nstates
    do iobj = 1,size(OBJFUNS)
        call bench%new(nspace, nstates, nkeep, nthr_max, .true., .false., .false., 0, objfun=trim(OBJFUNS(iobj)))
        call report_case(bench, nspace, nstates)
        call bench%get_sizes(pftsz, nrots, kfromto, nrefs)
        call bench%prod_kernels(1, nkeep, v1, bd1, bi1)
        fname = 'polar_kernel_identity_'//trim(OBJFUNS(iobj))//'_'//trim(casestr)//'.bin'
        open(newunit=funit, file=fname, access='stream', form='unformatted', action='write', status='replace')
        write(funit) nrots, nrefs, nkeep, PMB_NKINDS
        write(funit) v1, bd1, bi1
        close(funit)
        write(logfhandle,'(4x,a,a)') 'one-thread kernels written to ', fname
        do irep = 1,NREPS_CHECK
            wall = 0._dp
            !$ wall = omp_get_wtime()
            call bench%prod_kernels(nthr_max, nkeep, v, bd, bi)
            !$ wall = omp_get_wtime() - wall
            ndiff      = 0
            ndiff_best = 0
            do i = 1,nkeep
                do iref = 1,nrefs
                    do ikind = 1,PMB_NKINDS
                        if( any(v(:,iref,i,ikind) /= v1(:,iref,i,ikind)) ) ndiff(ikind) = ndiff(ikind) + 1
                    enddo
                    if( bd(iref,i) /= bd1(iref,i) .or. bi(iref,i) /= bi1(iref,i) ) ndiff_best = ndiff_best + 1
                enddo
            enddo
            inquire(file=IDENTITY_TSV, exist=l_new)
            l_new = .not. l_new
            open(newunit=funit, file=IDENTITY_TSV, position='append', action='write', status='unknown')
            if( l_new ) write(funit,'(a)') 'objfun'//achar(9)//'nspace'//achar(9)//'nstates'//achar(9)//'nrefs'//&
                &achar(9)//'nptcls'//achar(9)//'nthr'//achar(9)//'rep'//achar(9)//'pairs'//achar(9)//&
                &'diff_unshifted'//achar(9)//'diff_shifted'//achar(9)//'diff_bank_row'//achar(9)//&
                &'diff_bank_vals'//achar(9)//'diff_best'//achar(9)//'wall_s'
            write(funit,'(a,12(a,i0),a,f10.4)') trim(OBJFUNS(iobj)), achar(9), nspace, achar(9), nstates, achar(9),&
                &nrefs, achar(9), nkeep, achar(9), nthr_max, achar(9), irep, achar(9), nrefs*nkeep, achar(9), ndiff(1),&
                &achar(9), ndiff(2), achar(9), ndiff(3), achar(9), ndiff(4), achar(9), ndiff_best, achar(9), wall
            close(funit)
            write(logfhandle,'(4x,a,a,i0,a,i0,a,4(1x,i0),a,i0)') trim(OBJFUNS(iobj)), ' nthr=', nthr_max,&
                &' rep=', irep, ' rows differing from one thread:', ndiff, ' best: ', ndiff_best
        enddo
        call bench%kill
    enddo
    deallocate(bench)
end subroutine identity_case

!> The frozen-input harness of the probabilistic routes (prob_align_batch_pipeline.md, level 2 of section 8): each
!! route's batch pipeline (the fill, the selection and the extras the matchers run, in parts as the distributed
!! parts run it) on the bench's synthetic references and particles at one thread, one batch and one part (the
!! reference), then repeated as it is, at the full team, in batches of a third, in three parts, all three together,
!! and in three parts with the earlier parts' committed outputs perturbed before the later parts' fills (the
!! input-boundary regression of the iteration's model); the entries differing from the reference are counted by
!! particle and reference identity per category (rows of polar_fill_identity.tsv: the rows' distances, in-plane
!! indices, refined shifts, temperatures, seeds, assignments, extras). Every draw is keyed on the run, the pass,
!! the iteration, the particle and the draw, so every configuration must give the reference's entries: the run
!! fails on any differing entry.
subroutine fills_case( nspace, nstates, nptcls, nthr_max )
    !$ use omp_lib, only: omp_get_wtime
    integer, intent(in) :: nspace, nstates, nptcls, nthr_max
    integer, parameter :: NCONF = 6
    class(polar_matching_bench), allocatable, target :: bench
    type(pmb_fill_result) :: ref, res
    character(len=:), allocatable :: route
    integer  :: conf_nthr(NCONF), conf_bsz(NCONF), conf_nparts(NCONF)
    logical  :: conf_pert(NCONF)
    integer  :: iroute, iconf, ndiff(PMB_NDIFF), nfail, funit, bsz3, nparts3, idiff
    real(dp) :: wall
    logical  :: l_new
    bsz3        = max(1, (nptcls + 2) / 3)
    nparts3     = min(3, nptcls)
    ! the pipeline at one and nthr_max threads, two batch sizes, one and three parts, the earlier parts perturbed
    conf_nthr   = [1,       nthr_max, 1,       1,       nthr_max, 1      ]
    conf_bsz    = [nptcls,  nptcls,   bsz3,    nptcls,  bsz3,     nptcls ]
    conf_nparts = [1,       1,        1,       nparts3, nparts3,  nparts3]
    conf_pert   = [.false., .false.,  .false., .false., .false.,  .true. ]
    nfail = 0
    allocate(bench)
    do iroute = 1,PMB_NROUTES
        route = pmb_route_name(iroute)
        call bench%new(nspace, nstates, nptcls, nthr_max, .true., .false., .false., 0, fills_route=iroute)
        call report_case(bench, nspace, nstates)
        write(logfhandle,'(4x,a)') 'route '//route//': the reference fill at one thread, one batch, one part'
        call bench%run_fill(1, nptcls, 1, ref)
        if( iroute == PMB_ROUTE_3D_DENSE .or. iroute == PMB_ROUTE_3D_SEEDED ) call check_dense_refinement(ref, nspace, nstates)
        do iconf = 1,NCONF
            wall = 0._dp
            !$ wall = omp_get_wtime()
            call bench%run_fill(conf_nthr(iconf), conf_bsz(iconf), conf_nparts(iconf), res, l_perturb=conf_pert(iconf))
            !$ wall = omp_get_wtime() - wall
            call compare_fill_results(ref, res, ndiff)
            inquire(file=FILLS_TSV, exist=l_new)
            l_new = .not. l_new
            open(newunit=funit, file=FILLS_TSV, position='append', action='write', status='unknown')
            if( l_new )then
                write(funit,'(a)',advance='no') 'route'//achar(9)//'nspace'//achar(9)//'nstates'//achar(9)//'nptcls'//&
                    &achar(9)//'nthr'//achar(9)//'batchsz'//achar(9)//'nparts'//achar(9)//'perturbed'//achar(9)//'entries'
                do idiff = 1,PMB_NDIFF
                    write(funit,'(a)',advance='no') achar(9)//'diff_'//trim(PMB_DIFF_NAMES(idiff))
                enddo
                write(funit,'(a)') achar(9)//'wall_s'
            endif
            write(funit,'(a,8(a,i0))',advance='no') route, achar(9), nspace, achar(9), nstates, achar(9), nptcls,&
                &achar(9), conf_nthr(iconf), achar(9), conf_bsz(iconf), achar(9), conf_nparts(iconf), achar(9), &
                &merge(1, 0, conf_pert(iconf)), achar(9), res%nentries
            do idiff = 1,PMB_NDIFF
                write(funit,'(a,i0)',advance='no') achar(9), ndiff(idiff)
            enddo
            write(funit,'(a,f10.4)') achar(9), wall
            close(funit)
            write(logfhandle,'(4x,a,a,i0,a,i0,a,i0,a,l1,a,i0,a,7(1x,i0),a,f8.3,a)') route, ' nthr=', conf_nthr(iconf),&
                &' batchsz=', conf_bsz(iconf), ' nparts=', conf_nparts(iconf), ' perturbed=', conf_pert(iconf),&
                &' entries=', res%nentries,&
                &' differing from the reference (dists inpls shifts temps seeds assignment extras):', ndiff, ' in ', wall, ' s'
            if( any(ndiff /= 0) )then
                nfail = nfail + 1
                write(logfhandle,'(4x,a)') 'FAIL: a configuration of route '//route//' differs from its reference'
            endif
        enddo
        call bench%kill
    enddo
    call kill_fill_result(ref)
    call kill_fill_result(res)
    deallocate(bench)
    if( nfail > 0 ) THROW_HARD('the fill harness differs between configurations; see '//FILLS_TSV)
    write(logfhandle,'(a)') '>>> the fill harness gives identical rows, selections and extras in every configuration'
end subroutine fills_case

!> the dense 3D fill's refinement on the reference fill: the budget projs_ns split over the states' blocks, so no block
!! holds more refined entries than its share (ceil(projs_ns / nstates)) and no particle more than projs_ns up to
!! rounding; the shift map records a refinement only where it beats the stored shift, so the particles left without a
!! refined entry in some block are reported, not held to zero
subroutine check_dense_refinement( res, nspace, nstates )
    type(pmb_fill_result), intent(in) :: res
    integer,               intent(in) :: nspace, nstates
    integer :: i, s, nlack, nover, nref_s, ntot, nmax, nper, nmax_blk
    nper = res%projs_ns
    if( nstates > 1 ) nper = min(nspace, (res%projs_ns + nstates - 1) / nstates)
    nlack    = 0
    nover    = 0
    nmax     = 0
    nmax_blk = 0
    do i = 1,res%nptcls
        ntot = 0
        do s = 1,nstates
            nref_s = count(res%shifts(1,(s-1)*nspace+1:s*nspace,i) < 0.5 * huge(1._sp))
            if( nref_s == 0 ) nlack = nlack + 1
            if( nref_s > nper ) nover = nover + 1
            nmax_blk = max(nmax_blk, nref_s)
            ntot = ntot + nref_s
        enddo
        nmax = max(nmax, ntot)
        if( ntot > res%projs_ns + nstates - 1 ) nover = nover + 1
    enddo
    write(logfhandle,'(4x,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') '3D_dense refinement: budget ', res%projs_ns, &
        &', states ', nstates, ', share per state ', nper, ', most refined per particle ', nmax, ' and per state ', &
        &nmax_blk, ', over the share or budget ', nover, ', particles without a refined entry in a state ', nlack
    if( nover > 0 ) THROW_HARD('the dense fill refines more than a state''s share or the budget')
end subroutine check_dense_refinement

!> The neighbourhood route of one case (decision 19 of prob_table_matrix_scoring.md): per path and team (the full one,
!! then one thread), NEIGH_NREPS repeats of the per-pair routine and of the block path at the particles' temperatures
!! (timed; their entries against the first repeat at the full team), and the per-pair routine over the blocks of
!! decision 18's order (reported: how a fill in blocks would score a path left per pair); the block path once at a
!! team of NEIGH_TEAM_ID (identity); then the check (rulings R7): at the temperature scale that makes the draw the best
!! rotation, the best rotations differ in at most max(1, BLOCK_FLIP_RATE of the pairs) and the distances lie
!! within BLOCK_RTOL where they agree; the drawn rotations' differences at the particles' temperatures are reported
subroutine neigh_case( nspace, nstates, pgrp, nptcls, nthr_max )
    use, intrinsic :: iso_fortran_env, only: int64
    integer,          intent(in) :: nspace, nstates, nptcls, nthr_max
    character(len=*), intent(in) :: pgrp
    class(polar_neigh_bench), allocatable :: bench
    type(pnb_pass), allocatable :: pp, bk, po, pp_first, bk_first, bk_team, pp_best, bk_best
    integer, allocatable :: teams(:)
    real(dp) :: best_rel, drawn_rel
    real     :: cone, load
    integer  :: nsubs, path, it, irep, pftsz, nrots, kfromto(2), ncap, maxcell, ndiff_pp, ndiff_bk, ndiff_team
    integer  :: best_rot, drawn_rot, nallowed, npooled_diff, funit, npairs
    logical  :: l_new, l_ok, l_ident
    ! solve3D's subspace: 126 representatives per 2,500 directions, even
    nsubs = max(2, round2even(real(nint(126. * real(nspace) / 2500.))))
    allocate(bench, pp, bk, po, pp_first, bk_first, bk_team, pp_best, bk_best)
    call bench%new(nspace, nstates, nsubs, pgrp, nptcls, nthr_max)
    call bench%get_neigh_sizes(pftsz, nrots, kfromto, nsubs, ncap, maxcell, cone)
    write(logfhandle,'(a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,f5.1,a,3(1x,f7.1))') '>>> POLAR_NEIGH_BENCH case pgrp=',&
        &trim(pgrp), ' nspace=', nspace, ' nstates=', nstates, ' nsubs=', nsubs, ' nptcls=', nptcls, ' pftsz=', pftsz,&
        &' k=', kfromto(1), '..', kfromto(2), ' pool<=', ncap, ' cell<=', maxcell, ' cone=', cone,&
        &' pairs/particle (reps, members, cone):', real(bench%npairs(PNB_PATH_REPS)) / real(nptcls),&
        &real(bench%npairs(PNB_PATH_MEMBERS)) / real(nptcls), real(bench%npairs(PNB_PATH_CONE)) / real(nptcls)
    teams = [nthr_max]
    if( nthr_max > 1 ) teams = [teams, 1]
    do path = 1,PNB_NPATHS
        npairs  = bench%npairs(path)
        l_ident = .true.
        do it = 1,size(teams)
            do irep = 1,NEIGH_NREPS
                load = load_average()
                call bench%score_path(path, .false., teams(it), 1., pp)
                call bench%score_path(path, .true.,  teams(it), 1., bk)
                call bench%score_path(path, .false., teams(it), 1., po, l_cellorder=.true.)
                ndiff_pp = 0
                ndiff_bk = 0
                if( it == 1 .and. irep == 1 )then
                    pp_first = pp
                    bk_first = bk
                else
                    ndiff_pp = ndiff_pass(pp, pp_first)
                    ndiff_bk = ndiff_pass(bk, bk_first)
                endif
                ndiff_pp = ndiff_pp + ndiff_pass(po, pp_first)
                l_ident = l_ident .and. ndiff_pp == 0 .and. ndiff_bk == 0
                inquire(file=NEIGH_TSV, exist=l_new)
                l_new = .not. l_new
                open(newunit=funit, file=NEIGH_TSV, position='append', action='write', status='unknown')
                if( l_new ) write(funit,'(a)') 'path'//achar(9)//'pgrp'//achar(9)//'nspace'//achar(9)//'nstates'//&
                    &achar(9)//'nsubs'//achar(9)//'nptcls'//achar(9)//'nthr'//achar(9)//'rep'//achar(9)//'B_even'//&
                    &achar(9)//'B_odd'//achar(9)//'needed_pairs'//achar(9)//'scored_pairs'//achar(9)//'wasted_pairs'//&
                    &achar(9)//'calls'//achar(9)//'rows_per_call'//achar(9)//'refs_per_call'//achar(9)//'slot_fill'//&
                    &achar(9)//'per_pair_tsec'//achar(9)//'per_pair_wall'//achar(9)//'block_tsec'//achar(9)//&
                    &'block_wall'//achar(9)//'ratio_tsec'//achar(9)//'ratio_wall'//achar(9)//'per_pair_us'//achar(9)//&
                    &'block_us'//achar(9)//'load'//achar(9)//'per_pair_diff_vs_first'//achar(9)//'block_diff_vs_first'//&
                    &achar(9)//'per_pair_cellorder_tsec'//achar(9)//'ratio_cellorder'
                write(funit,'(a,a,a,9(a,i0),3(a,i0),3(a,f10.4),4(a,es12.5),4(a,f10.4),a,f6.2,2(a,i0),a,es12.5,a,f10.4)')&
                    &path_name(path),&
                    &achar(9), trim(pgrp), achar(9), nspace, achar(9), nstates, achar(9), nsubs, achar(9), nptcls,&
                    &achar(9), teams(it), achar(9), irep, achar(9), bk%bsz(1), achar(9), bk%bsz(2), achar(9), npairs,&
                    &achar(9), bk%nscored, achar(9), bk%nwasted, achar(9), bk%ncalls,&
                    &achar(9), real(bk%nrows) / real(max(1_int64,bk%ncalls)),&
                    &achar(9), real(bk%ncols) / real(max(1_int64,bk%ncalls)),&
                    &achar(9), real(bk%nscored) / real(max(1_int64,bk%nslots)),&
                    &achar(9), pp%tsec, achar(9), pp%wall, achar(9), bk%tsec, achar(9), bk%wall,&
                    &achar(9), bk%tsec / max(pp%tsec, 1.e-12_dp), achar(9), bk%wall / max(pp%wall, 1.e-12_dp),&
                    &achar(9), 1.e6_dp * pp%tsec / real(npairs,dp), achar(9), 1.e6_dp * bk%tsec / real(npairs,dp),&
                    &achar(9), load, achar(9), ndiff_pp, achar(9), ndiff_bk, achar(9), po%tsec,&
                    &achar(9), po%tsec / max(pp%tsec, 1.e-12_dp)
                close(funit)
                write(logfhandle,'(4x,a,a,i0,a,i0,a,f8.3,a,f8.3,a,f6.3,a,i0,a,i0)') path_name(path), ' nthr=', teams(it),&
                    &' rep=', irep, ' per_pair_tsec=', pp%tsec, ' block_tsec=', bk%tsec, ' ratio=',&
                    &bk%tsec / max(pp%tsec, 1.e-12_dp), ' diff_vs_first=', ndiff_pp, '/', ndiff_bk
            enddo
        enddo
        ndiff_team = -1
        if( nthr_max > NEIGH_TEAM_ID )then
            call bench%score_path(path, .true., NEIGH_TEAM_ID, 1., bk_team)
            ndiff_team = ndiff_pass(bk_team, bk_first)
            l_ident    = l_ident .and. ndiff_team == 0
        endif
        ! the check: the best rotations and their distances, then the drawn rotations at the particles' temperatures
        call bench%score_path(path, .false., nthr_max, BLOCK_TSCALES(1), pp_best)
        call bench%score_path(path, .true.,  nthr_max, BLOCK_TSCALES(1), bk_best)
        call compare_pass(bk_best,  pp_best,  best_rot,  best_rel)
        call compare_pass(bk_first, pp_first, drawn_rot, drawn_rel)
        npooled_diff = -1
        if( path == PNB_PATH_REPS ) npooled_diff = bench%pooled_changes(bk_first%dists)
        nallowed = max(1, int(BLOCK_FLIP_RATE * real(npairs)))
        l_ok     = best_rot <= nallowed .and. best_rel <= real(BLOCK_RTOL,dp)
        inquire(file=NEIGHCHECK_TSV, exist=l_new)
        l_new = .not. l_new
        open(newunit=funit, file=NEIGHCHECK_TSV, position='append', action='write', status='unknown')
        if( l_new ) write(funit,'(a)') 'path'//achar(9)//'pgrp'//achar(9)//'nspace'//achar(9)//'nstates'//achar(9)//&
            &'nptcls'//achar(9)//'nthr'//achar(9)//'pairs'//achar(9)//'best_rot_diff'//achar(9)//'allowed'//achar(9)//&
            &'best_max_rel'//achar(9)//'drawn_rot_diff'//achar(9)//'drawn_max_rel'//achar(9)//'team4_diff_vs_first'//&
            &achar(9)//'pooled_changes'//achar(9)//'within_tolerance'//achar(9)//'identical'
        write(funit,'(a,a,a,7(a,i0),a,es12.5,a,i0,a,es12.5,2(a,i0),2(a,a))') path_name(path), achar(9), trim(pgrp),&
            &achar(9), nspace, achar(9), nstates, achar(9), nptcls, achar(9), nthr_max, achar(9), npairs, achar(9),&
            &best_rot, achar(9), nallowed, achar(9), best_rel, achar(9), drawn_rot, achar(9), drawn_rel, achar(9),&
            &ndiff_team, achar(9), npooled_diff, achar(9), trim(merge('yes','no ',l_ok)), achar(9),&
            &trim(merge('yes','no ',l_ident))
        close(funit)
        write(logfhandle,'(4x,a,a,i0,a,i0,a,es10.3,a,i0,a,es10.3,a,i0,a,i0,a,a,a,a)') path_name(path), ' best_rot_diff=',&
            &best_rot, ' of ', npairs, ' best_max_rel=', best_rel, ' drawn_rot_diff=', drawn_rot, ' drawn_max_rel=',&
            &drawn_rel, ' team4_diff=', ndiff_team, ' pooled_changes=', npooled_diff, ' within_tolerance=',&
            &trim(merge('yes','no ',l_ok)), ' identical=', trim(merge('yes','no ',l_ident))
    enddo
    call bench%kill
    deallocate(bench, pp, bk, po, pp_first, bk_first, bk_team, pp_best, bk_best)

contains

    ! entries differing bit for bit
    integer function ndiff_pass( a, b )
        type(pnb_pass), intent(in) :: a, b
        ndiff_pass = count(a%dists /= b%dists) + count(a%irots /= b%irots)
    end function ndiff_pass

    ! the block path against the per-pair routine: rotations differing; the largest relative difference of the
    ! distance where they agree
    subroutine compare_pass( a, ref, nrot, max_rel )
        type(pnb_pass), intent(in)  :: a, ref
        integer,        intent(out) :: nrot
        real(dp),       intent(out) :: max_rel
        integer :: j
        nrot    = 0
        max_rel = 0._dp
        do j = 1,size(ref%dists)
            if( a%irots(j) /= ref%irots(j) )then
                nrot = nrot + 1
            else
                max_rel = max(max_rel, abs(real(a%dists(j),dp) - real(ref%dists(j),dp)) /&
                    &max(abs(real(ref%dists(j),dp)), 1.e-30_dp))
            endif
        enddo
    end subroutine compare_pass

    function path_name( path_here ) result( name )
        integer, intent(in) :: path_here
        character(len=:), allocatable :: name
        select case(path_here)
            case(PNB_PATH_REPS)
                name = 'a_reps'
            case(PNB_PATH_MEMBERS)
                name = 'b_members'
            case(PNB_PATH_CONE)
                name = 'c_cone'
            case DEFAULT
                name = 'unknown'
        end select
    end function path_name

end subroutine neigh_case

subroutine report_case( bench, nspace, nstates )
    class(polar_matching_bench), intent(in) :: bench
    integer,                     intent(in) :: nspace, nstates
    integer :: pftsz, nrots, kfromto(2), nrefs
    call bench%get_sizes(pftsz, nrots, kfromto, nrefs)
    write(logfhandle,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') '>>> POLAR_MATCHING_BENCH case nspace=', nspace,&
        &' nstates=', nstates, ' nrefs=', nrefs, ' pftsz=', pftsz, ' nrots=', nrots, ' k=', kfromto(1), '..',&
        &kfromto(2)
end subroutine report_case

function layout_name( layout ) result( name )
    integer, intent(in) :: layout
    character(len=:), allocatable :: name
    select case(layout)
        case(PMB_PER_PAIR)
            name = 'per_pair'
        case(PMB_PER_PAIR_CC)
            name = 'per_pair_cc'
        case(PMB_REPLICA)
            name = 'replica'
        case(PMB_BLOCKED)
            name = 'blocked'
        case(PMB_STATES_INNER)
            name = 'states_inner'
        case(PMB_MATRIX)
            name = 'matrix'
        case(PMB_BLOCK_PROD)
            name = 'block_prod'
        case(PMB_MATRIX_T)
            name = 'matrix_t'
        case DEFAULT
            name = 'unknown'
    end select
end function layout_name

!> the one-minute load average (-1 when unavailable)
real function load_average()
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    interface
        integer(c_int) function getloadavg( loadavg, nelem ) bind(C, name='getloadavg')
            import :: c_int, c_double
            real(c_double)        :: loadavg(*)
            integer(c_int), value :: nelem
        end function getloadavg
    end interface
    real(c_double) :: loads(3)
    load_average = -1.
    if( getloadavg(loads, 3_c_int) >= 1 ) load_average = real(loads(1))
end function load_average

end submodule simple_commanders_test_highlevel_matchbench
