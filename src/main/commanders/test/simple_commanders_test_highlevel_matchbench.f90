!@descr: polar matching layout benchmark of simple_commanders_test_highlevel (measure only, no CTest entry)
submodule (simple_commanders_test_highlevel) simple_commanders_test_highlevel_matchbench
use simple_polar_matching_bench, only: polar_matching_bench, PMB_PER_PAIR, PMB_REPLICA, PMB_BLOCKED, PMB_STATES_INNER,&
    &PMB_MATRIX, PMB_MATRIX_T, PMB_PER_PAIR_CC, PMB_NKINDS
implicit none
#include "simple_local_flags.inc"

character(len=*), parameter :: MATCHBENCH_TSV   = 'polar_matching_bench.tsv'
character(len=*), parameter :: MATCHCHECK_TSV   = 'polar_matching_check.tsv'
character(len=*), parameter :: IDENTITY_TSV     = 'polar_kernel_identity.tsv'
integer,          parameter :: NKEEP_IDENTITY   = 32     !< particles of suite identity unless nptcls is given
integer,          parameter :: NREPS            = 3      !< timing: best of three after one warm-up
integer,          parameter :: NREPS_CHECK      = 5      !< identity: repeats at the full team
integer,          parameter :: NPTCLS_ONE_THREAD = 200
integer,          parameter :: BLOCKED_BSZ(4)   = [4, 8, 16, 32]
integer,          parameter :: BLOCKED_RSZ      = 8      !< 8 references of one half: 1.4 MB, within a performance core's share of L2
integer,          parameter :: MATRIX_RSZ(3)    = [32, 64, 128]
integer,          parameter :: MATRIX_BSZ(4)    = [4, 8, 16, 32]

contains

module subroutine exec_test_polar_matching_bench( self, cline )
    !$ use omp_lib, only: omp_get_num_procs
    class(commander_test_polar_matching_bench), intent(inout) :: self
    class(cmdline),                             intent(inout) :: cline
    integer, parameter :: NSPACE_ALL(2) = [1000, 2500], NSTATES_ALL(3) = [1, 2, 4]
    integer, allocatable :: nspaces(:), nstates_list(:), nthrs(:), layouts(:)
    type(string) :: suite
    integer      :: nptcls, nthr_max, ispace, istate
    logical      :: l_check, l_identity
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
    select case(suite%to_char())
        case('list')
            write(logfhandle,'(a)') 'Available suites for polar_matching_bench: per_pair, replica, blocked, states_inner,&
                & matrix, matrix_t, per_pair_cc, check, identity'
            return
        case('')
            layouts = [PMB_PER_PAIR, PMB_REPLICA, PMB_BLOCKED, PMB_STATES_INNER, PMB_MATRIX, PMB_MATRIX_T]
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
        case('per_pair_cc')
            layouts = [PMB_PER_PAIR_CC]
        case('check')
            l_check = .true.
        case('identity')
            l_identity = .true.
        case default
            THROW_HARD('no suite '//suite%to_char()//' in polar_matching_bench; use suite=list')
    end select
    if( nptcls < 1 .or. nthr_max < 1 .or. any(nspaces < 1) .or. any(nstates_list < 1) )then
        THROW_HARD('polar_matching_bench needs nptcls, nthr, nspace and nstates >= 1')
    endif
    ! the teams, the full one first (the tile sweeps run there): all threads, the five performance cores, one thread
    nthrs = [nthr_max]
    if( nthr_max > 5 ) nthrs = [nthrs, 5]
    if( nthr_max > 1 ) nthrs = [nthrs, 1]
    do ispace = 1,size(nspaces)
        do istate = 1,size(nstates_list)
            if( l_check )then
                call check_case(nspaces(ispace), nstates_list(istate), nthr_max)
            else if( l_identity )then
                if( cline%defined('nptcls') )then
                    call identity_case(nspaces(ispace), nstates_list(istate), nptcls, nthr_max)
                else
                    call identity_case(nspaces(ispace), nstates_list(istate), NKEEP_IDENTITY, nthr_max)
                endif
            else
                call time_case(nspaces(ispace), nstates_list(istate), nptcls, nthrs, layouts)
            endif
        enddo
    enddo
    call simple_end('**** SIMPLE_TEST_POLAR_MATCHING_BENCH NORMAL STOP ****')
end subroutine exec_test_polar_matching_bench

!> The timings of one case: per layout and team, best of NREPS after a warm-up; the tile sizes of layouts 2 and 4
!! are swept at the full team and the fastest is timed at the others
subroutine time_case( nspace, nstates, nptcls, nthrs, layouts )
    integer, intent(in) :: nspace, nstates, nptcls, nthrs(:), layouts(:)
    class(polar_matching_bench), allocatable, target :: bench
    integer, allocatable :: bszs(:), rszs(:)
    real(dp) :: wall, wall_best
    integer  :: ilay, it, nt, best_b, best_r, ib, ir
    logical  :: l_prod, l_memo, l_gmat, l_cc
    l_prod = any(layouts == PMB_PER_PAIR) .or. any(layouts == PMB_PER_PAIR_CC)
    l_cc   = any(layouts == PMB_PER_PAIR_CC)
    if( l_cc .and. size(layouts) > 1 ) THROW_HARD('per_pair_cc is timed on its own')
    l_memo = any(layouts == PMB_REPLICA) .or. any(layouts == PMB_BLOCKED) .or. any(layouts == PMB_STATES_INNER)
    l_gmat = any(layouts == PMB_MATRIX) .or. any(layouts == PMB_MATRIX_T)
    allocate(bench)
    if( l_cc )then
        call bench%new(nspace, nstates, nptcls, maxval(nthrs), l_prod, l_memo, l_gmat, 0, objfun='cc')
    else
        call bench%new(nspace, nstates, nptcls, maxval(nthrs), l_prod, l_memo, l_gmat, 0)
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
                        wall = time_one(bench, layouts(ilay), nspace, nstates, nptcls, nt, bszs(ib), rszs(ir))
                        if( wall < wall_best )then
                            wall_best = wall
                            best_b    = bszs(ib)
                            best_r    = rszs(ir)
                        endif
                    enddo
                enddo
            else
                wall = time_one(bench, layouts(ilay), nspace, nstates, nptcls, nt, best_b, best_r)
            endif
        enddo
    enddo
    call bench%kill
    deallocate(bench)
end subroutine time_case

!> one timing row: warm-up, best of NREPS, appended to the TSV; returns the best wall time
function time_one( bench, layout, nspace, nstates, nptcls, nthr, bsz, rsz ) result( wall_best )
    class(polar_matching_bench), intent(inout) :: bench
    integer,                     intent(in)    :: layout, nspace, nstates, nptcls, nthr, bsz, rsz
    real(dp) :: wall_best, walls(NREPS), pairs, us_pair, pairs_per_s, ref_gbps, parts(3)
    real     :: load1
    integer  :: nptcls_use, nwarm, irep, funit, pftsz, nrots, kfromto(2), nsample, nrefs
    logical  :: l_new
    call bench%get_sizes(pftsz, nrots, kfromto, nsample, nrefs)
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
        &'kfrom'//achar(9)//'kto'//achar(9)//'nsample'//achar(9)//'load1'//achar(9)//'wall_s'//achar(9)//&
        &'wall_rep1'//achar(9)//'wall_rep2'//achar(9)//'wall_rep3'//achar(9)//'us_pair_thread'//achar(9)//&
        &'pairs_per_s'//achar(9)//'ref_GBps'//achar(9)//'memo_GB'//achar(9)//'check'//achar(9)//&
        &'us_sums'//achar(9)//'us_fft'//achar(9)//'us_draw'
    write(funit,'(a,12(a,i0),a,f6.2,4(a,f10.4),a,f9.3,a,es12.5,a,f8.3,a,f8.3,a,a,3(a,f9.3))') layout_name(layout),&
        &achar(9), nspace, achar(9), nstates, achar(9), nrefs, achar(9), nptcls_use, achar(9), nthr, achar(9), bsz,&
        &achar(9), rsz, achar(9), pftsz, achar(9), nrots, achar(9), kfromto(1), achar(9), kfromto(2), achar(9), nsample,&
        &achar(9), load1, achar(9), wall_best, achar(9), walls(1), achar(9), walls(2), achar(9), walls(3), achar(9),&
        &us_pair, achar(9), pairs_per_s, achar(9), ref_gbps, achar(9), bench%memo_gbytes(layout), achar(9), '-',&
        &achar(9), parts(1), achar(9), parts(2), achar(9), parts(3)
    close(funit)
    write(logfhandle,'(4x,a,a,i0,a,i0,a,i0,a,i0,a,f6.2,a,f9.3,a,f8.2,a,f8.2,a)') layout_name(layout), ' nthr=', nthr,&
        &' B=', bsz, ' R=', rsz, ' nptcls=', nptcls_use, ' load=', load1, ' wall=', wall_best, ' s, ', us_pair,&
        &' us/pair/thread, ', ref_gbps, ' GB/s'
end function time_one

!> The identity checks of one case on nkeep particles. Production (layout 1) at one thread is the reference of
!! the nsample lowest distances and their rotations; the replica at one thread that of the full distance vectors.
subroutine check_case( nspace, nstates, nthr_max )
    integer, intent(in) :: nspace, nstates, nthr_max
    class(polar_matching_bench), allocatable, target :: bench
    real(sp), allocatable :: prod_d(:,:,:), vecs_ref(:,:,:), vecs_one(:,:,:)
    integer,  allocatable :: prod_i(:,:,:)
    integer :: nkeep, irep, ib
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
        real(sp), allocatable :: vecs(:,:,:), tops_d(:,:,:)
        integer,  allocatable :: tops_i(:,:,:)
        real(dp) :: wall, max_abs, max_rel, d, max_top_abs
        integer  :: nrefs, p, iref, i, n_vec_diff, n_argmin_diff, n_top_diff, n_top_ptcls, n_vec_diff_one, funit
        logical  :: l_ptcl_diff
        logical  :: l_new, l_vecs
        wall = bench%score(layout, nthr, nkeep, bsz, rsz)
        call bench%get_tops(tops_d, tops_i)
        nrefs = size(tops_d,2)
        l_vecs = layout /= PMB_PER_PAIR
        if( l_vecs ) call bench%get_vecs(vecs)
        if( l_make_ref )then
            if( layout == PMB_PER_PAIR )then
                prod_d = tops_d
                prod_i = tops_i
            else
                vecs_ref = vecs
            endif
        endif
        ! the nsample lowest distances and their rotations against production's, byte for byte
        n_top_diff  = 0
        n_top_ptcls = 0
        max_top_abs = 0._dp
        do i = 1,nkeep
            l_ptcl_diff = .false.
            do iref = 1,nrefs
                if( any(tops_d(:,iref,i) /= prod_d(:,iref,i)) .or. any(tops_i(:,iref,i) /= prod_i(:,iref,i)) )then
                    n_top_diff  = n_top_diff + 1
                    l_ptcl_diff = .true.
                endif
                max_top_abs = max(max_top_abs, maxval(abs(real(tops_d(:,iref,i),dp) - real(prod_d(:,iref,i),dp))))
            enddo
            if( l_ptcl_diff ) n_top_ptcls = n_top_ptcls + 1
        enddo
        ! the full vectors against this layout's own one-thread run (irep 0), thread identity
        n_vec_diff_one = -1
        if( l_vecs )then
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
        endif
        ! the full vectors against the replica's
        n_vec_diff    = -1
        n_argmin_diff = -1
        max_abs       = -1._dp
        max_rel       = -1._dp
        if( l_vecs .and. allocated(vecs_ref) )then
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
            &'top_diff_vs_prod'//achar(9)//'vec_diff_vs_replica'//achar(9)//'argmin_diff_vs_replica'//achar(9)//&
            &'max_abs_diff'//achar(9)//'max_rel_diff'//achar(9)//'wall_s'//achar(9)//'top_diff_ptcls'//achar(9)//&
            &'max_top_abs_diff'//achar(9)//'vec_diff_vs_own_1thr'
        write(funit,'(a,12(a,i0),2(a,es12.5),a,f10.4,a,i0,a,es12.5,a,i0)') layout_name(layout), achar(9), nspace, achar(9), nstates,&
            &achar(9), nrefs, achar(9), nkeep, achar(9), nthr, achar(9), irep, achar(9), bsz, achar(9), rsz, achar(9),&
            &nrefs*nkeep, achar(9), n_top_diff, achar(9), n_vec_diff, achar(9), n_argmin_diff, achar(9), max_abs,&
            &achar(9), max_rel, achar(9), wall, achar(9), n_top_ptcls, achar(9), max_top_abs, achar(9), n_vec_diff_one
        close(funit)
        write(logfhandle,'(4x,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,es10.3)') layout_name(layout), ' nthr=', nthr,&
            &' rep=', irep, ' B=', bsz, ' R=', rsz, ' top_diff_vs_prod=', n_top_diff, ' vec_diff=', n_vec_diff,&
            &' argmin_diff=', n_argmin_diff, ' max_rel=', max_rel
    end subroutine check_run

end subroutine check_case

!> The production kernels of one case on nkeep particles, euclid and cc (prod_kernels): the one-thread result is
!! written to polar_kernel_identity_<objfun>_<nspace>x<nstates>.bin for a comparison between builds, and NREPS_CHECK
!! runs at the full team are compared with it byte for byte (rows of polar_kernel_identity.tsv)
subroutine identity_case( nspace, nstates, nkeep, nthr_max )
    !$ use omp_lib, only: omp_get_wtime
    integer, intent(in) :: nspace, nstates, nkeep, nthr_max
    character(len=*), parameter :: OBJFUNS(2) = ['euclid', 'cc    ']
    class(polar_matching_bench), allocatable, target :: bench
    real(sp), allocatable :: v1(:,:,:,:), td1(:,:,:), bd1(:,:), v(:,:,:,:), td(:,:,:), bd(:,:)
    integer,  allocatable :: ti1(:,:,:), bi1(:,:), ti(:,:,:), bi(:,:)
    character(len=:), allocatable :: fname
    character(len=32) :: casestr
    real(dp) :: wall
    integer  :: iobj, irep, funit, ikind, iref, i, ndiff(PMB_NKINDS), ndiff_tops, ndiff_best
    integer  :: pftsz, nrots, kfromto(2), nsample, nrefs
    logical  :: l_new
    allocate(bench)
    write(casestr,'(i0,a,i0)') nspace, 'x', nstates
    do iobj = 1,size(OBJFUNS)
        call bench%new(nspace, nstates, nkeep, nthr_max, .true., .false., .false., 0, objfun=trim(OBJFUNS(iobj)))
        call report_case(bench, nspace, nstates)
        call bench%get_sizes(pftsz, nrots, kfromto, nsample, nrefs)
        call bench%prod_kernels(1, nkeep, v1, td1, ti1, bd1, bi1)
        fname = 'polar_kernel_identity_'//trim(OBJFUNS(iobj))//'_'//trim(casestr)//'.bin'
        open(newunit=funit, file=fname, access='stream', form='unformatted', action='write', status='replace')
        write(funit) nrots, nrefs, nkeep, PMB_NKINDS, nsample
        write(funit) v1, td1, ti1, bd1, bi1
        close(funit)
        write(logfhandle,'(4x,a,a)') 'one-thread kernels written to ', fname
        do irep = 1,NREPS_CHECK
            wall = 0._dp
            !$ wall = omp_get_wtime()
            call bench%prod_kernels(nthr_max, nkeep, v, td, ti, bd, bi)
            !$ wall = omp_get_wtime() - wall
            ndiff      = 0
            ndiff_tops = 0
            ndiff_best = 0
            do i = 1,nkeep
                do iref = 1,nrefs
                    do ikind = 1,PMB_NKINDS
                        if( any(v(:,iref,i,ikind) /= v1(:,iref,i,ikind)) ) ndiff(ikind) = ndiff(ikind) + 1
                    enddo
                    if( any(td(:,iref,i) /= td1(:,iref,i)) .or. any(ti(:,iref,i) /= ti1(:,iref,i)) )then
                        ndiff_tops = ndiff_tops + 1
                    endif
                    if( bd(iref,i) /= bd1(iref,i) .or. bi(iref,i) /= bi1(iref,i) ) ndiff_best = ndiff_best + 1
                enddo
            enddo
            inquire(file=IDENTITY_TSV, exist=l_new)
            l_new = .not. l_new
            open(newunit=funit, file=IDENTITY_TSV, position='append', action='write', status='unknown')
            if( l_new ) write(funit,'(a)') 'objfun'//achar(9)//'nspace'//achar(9)//'nstates'//achar(9)//'nrefs'//&
                &achar(9)//'nptcls'//achar(9)//'nthr'//achar(9)//'rep'//achar(9)//'pairs'//achar(9)//&
                &'diff_unshifted'//achar(9)//'diff_shifted'//achar(9)//'diff_bank_draw'//achar(9)//&
                &'diff_bank_vals'//achar(9)//'diff_tops'//achar(9)//'diff_best'//achar(9)//'wall_s'
            write(funit,'(a,13(a,i0),a,f10.4)') trim(OBJFUNS(iobj)), achar(9), nspace, achar(9), nstates, achar(9),&
                &nrefs, achar(9), nkeep, achar(9), nthr_max, achar(9), irep, achar(9), nrefs*nkeep, achar(9), ndiff(1),&
                &achar(9), ndiff(2), achar(9), ndiff(3), achar(9), ndiff(4), achar(9), ndiff_tops, achar(9),&
                &ndiff_best, achar(9), wall
            close(funit)
            write(logfhandle,'(4x,a,a,i0,a,i0,a,4(1x,i0),a,i0,a,i0)') trim(OBJFUNS(iobj)), ' nthr=', nthr_max,&
                &' rep=', irep, ' vectors differing from one thread:', ndiff, ' tops: ', ndiff_tops, ' best: ',&
                &ndiff_best
        enddo
        call bench%kill
    enddo
    deallocate(bench)
end subroutine identity_case

subroutine report_case( bench, nspace, nstates )
    class(polar_matching_bench), intent(in) :: bench
    integer,                     intent(in) :: nspace, nstates
    integer :: pftsz, nrots, kfromto(2), nsample, nrefs
    call bench%get_sizes(pftsz, nrots, kfromto, nsample, nrefs)
    write(logfhandle,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') '>>> POLAR_MATCHING_BENCH case nspace=', nspace,&
        &' nstates=', nstates, ' nrefs=', nrefs, ' pftsz=', pftsz, ' nrots=', nrots, ' k=', kfromto(1), '..',&
        &kfromto(2), ' nsample=', nsample
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
