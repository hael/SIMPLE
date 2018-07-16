! projection-matching based on Hadamard products, high-level search routines for REFINE3D
module simple_strategy3D_matcher
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_strategy3D_alloc ! singleton s3D
use simple_timer
use simple_qsys_funs,                only: qsys_job_finished
use simple_binoris_io,               only: binwrite_oritab
use simple_kbinterpol,               only: kbinterpol
use simple_ori,                      only: ori
use simple_sym,                      only: sym
use simple_image,                    only: image
use simple_cmdline,                  only: cmdline
use simple_parameters,               only: params_glob
use simple_builder,                  only: build_glob
use simple_polarft_corrcalc,         only: polarft_corrcalc
use simple_strategy2D3D_common,      only: killrecvols, set_bp_range, preprecvols,&
    prepimgbatch, grid_ptcl, read_imgbatch, eonorm_struct_facts,norm_struct_facts
use simple_strategy3D_cluster,       only: strategy3D_cluster
use simple_strategy3D_cluster_snhc,  only: strategy3D_cluster_snhc
use simple_strategy3D_single,        only: strategy3D_single
use simple_strategy3D_multi,         only: strategy3D_multi
use simple_strategy3D_snhc_single,   only: strategy3D_snhc_single
use simple_strategy3D_greedy_single, only: strategy3D_greedy_single
use simple_strategy3D_hard_single,   only: strategy3D_hard_single
use simple_strategy3D_greedy_multi,  only: strategy3D_greedy_multi
use simple_strategy3D_hard_multi,    only: strategy3D_hard_multi
use simple_strategy3D_cont_single,   only: strategy3D_cont_single
use simple_strategy3D,               only: strategy3D
use simple_strategy3D_srch,          only: strategy3D_spec
use simple_convergence,              only: convergence
implicit none

public :: refine3D_exec, preppftcc4align, pftcc
private

logical, parameter             :: L_BENCH = .false., DEBUG = .false.
type(polarft_corrcalc), target :: pftcc
integer, allocatable           :: pinds(:)
logical, allocatable           :: ptcl_mask(:)
integer                        :: nptcls2update
integer(timer_int_kind)        :: t_init, t_prep_pftcc, t_align, t_rec, t_tot, t_prep_primesrch3D
real(timer_int_kind)           :: rt_init, rt_prep_pftcc, rt_align, rt_rec, rt_prep_primesrch3D
real(timer_int_kind)           :: rt_tot
character(len=STDLEN)          :: benchfname

contains

    subroutine refine3D_exec( cline, which_iter, converged )
        class(cmdline),        intent(inout) :: cline
        integer,               intent(in)    :: which_iter
        logical,               intent(inout) :: converged
        type(image),     allocatable :: rec_imgs(:)
        integer, target, allocatable :: symmat(:,:)
        logical,         allocatable :: het_mask(:)
        !---> The below is to allow particle-dependent decision about which 3D strategy to use
        type :: strategy3D_per_ptcl
            class(strategy3D), pointer :: ptr  => null()
        end type strategy3D_per_ptcl
        type(strategy3D_per_ptcl), allocatable :: strategy3Dsrch(:)
        !<---- hybrid or combined search strategies can then be implemented as extensions of the
        !      relevant strategy3D base class
        type(strategy3D_spec) :: strategy3Dspec
        type(ori)             :: orientation
        type(kbinterpol)      :: kbwin
        type(sym)             :: c1_symop
        type(ctfparams)       :: ctfvars
        type(convergence)     :: conv
        type(image)           :: mskimg
        real    :: frac_srch_space, extr_thresh, extr_score_thresh
        integer :: iptcl, iextr_lim, i, zero_pop, fnr, cnt, i_batch
        integer :: ibatch, npeaks, batchlims(2), updatecnt
        logical :: doprint, do_extr

        if( L_BENCH )then
            t_init = tic()
            t_tot  = t_init
        endif

        if( params_glob%dryrun.eq.'yes' )then
            ! TBD
        endif

        ! CHECK THAT WE HAVE AN EVEN/ODD PARTITIONING
        if( params_glob%eo .ne. 'no' )then
            if( build_glob%spproj_field%get_nevenodd() == 0 ) &
                call simple_stop('ERROR! no eo partitioning available; strategy3D_matcher :: refine3D_exec')
        else
            call build_glob%spproj_field%set_all2single('eo', -1.)
        endif

        ! SET FOURIER INDEX RANGE
        call set_bp_range(cline)

        ! DETERMINE THE NUMBER OF PEAKS
        select case(params_glob%refine)
            case('cluster', 'snhc', 'clustersym', 'cluster_snhc', 'cont_single')
                npeaks = 1
            case('hard_single','hard_multi')
                npeaks = MAXNPEAKS
                ! command line overrides
                if( cline%defined('npeaks') ) npeaks = params_glob%npeaks
            case DEFAULT
                npeaks = min(build_glob%eulspace%find_npeaks_from_athres(NPEAKSATHRES), MAXNPEAKS)
                ! command line overrides
                if( cline%defined('npeaks') ) npeaks = params_glob%npeaks
        end select
        if( DEBUG ) print *, '*** strategy3D_matcher ***: determined the number of peaks'

        ! SET FRACTION OF SEARCH SPACE
        frac_srch_space = build_glob%spproj_field%get_avg('frac')

        ! SETUP WEIGHTS
        if( params_glob%weights3D.eq.'yes' )then
            if( params_glob%nptcls <= SPECWMINPOP )then
                call build_glob%spproj_field%calc_hard_weights(params_glob%frac)
            else
                call build_glob%spproj_field%calc_spectral_weights(params_glob%frac)
            endif
        else
            call build_glob%spproj_field%calc_hard_weights(params_glob%frac)
        endif

        ! READ FOURIER RING CORRELATIONS
        select case(params_glob%refine)
            case('cluster', 'clustersym', 'cluster_snhc')
                if( file_exists(CLUSTER3D_FRCS) ) call build_glob%projfrcs%read(CLUSTER3D_FRCS)
            case DEFAULT
                if( file_exists(params_glob%frcs) ) call build_glob%projfrcs%read(params_glob%frcs)
        end select

        ! PARTICLE INDEX SAMPLING FOR FRACTIONAL UPDATE (OR NOT)
        if( allocated(pinds) )     deallocate(pinds)
        if( allocated(ptcl_mask) ) deallocate(ptcl_mask)
        if( params_glob%l_frac_update )then
            allocate(ptcl_mask(params_glob%fromp:params_glob%top))
            call build_glob%spproj_field%sample4update_and_incrcnt([params_glob%fromp,params_glob%top],&
            &params_glob%update_frac, nptcls2update, pinds, ptcl_mask)
            ! correct convergence stats
            do iptcl=params_glob%fromp,params_glob%top
                if( .not. ptcl_mask(iptcl) )then
                    ! these are not updated
                    call build_glob%spproj_field%set(iptcl, 'mi_proj',     1.0)
                    call build_glob%spproj_field%set(iptcl, 'mi_inpl',     1.0)
                    call build_glob%spproj_field%set(iptcl, 'mi_state',    1.0)
                    call build_glob%spproj_field%set(iptcl, 'mi_joint',    1.0)
                    call build_glob%spproj_field%set(iptcl, 'dist',        0.0)
                    call build_glob%spproj_field%set(iptcl, 'dist_inpl',   0.0)
                    call build_glob%spproj_field%set(iptcl, 'frac',      100.0)
                endif
            end do
        else
            nptcls2update = params_glob%top - params_glob%fromp + 1
            allocate(pinds(nptcls2update), ptcl_mask(params_glob%fromp:params_glob%top))
            pinds = (/(i,i=params_glob%fromp,params_glob%top)/)
            ptcl_mask = .true.
            call build_glob%spproj_field%incr_updatecnt([params_glob%fromp,params_glob%top])
        endif

        ! B-factor weighted reconstruction
        if( params_glob%shellw.eq.'yes' ) call build_glob%spproj_field%calc_bfac_rec

        ! EXTREMAL LOGICS
        do_extr = .false.
        select case(trim(params_glob%refine))
            case('cluster_snhc')
                if(allocated(het_mask))deallocate(het_mask)
                allocate(het_mask(params_glob%fromp:params_glob%top), source=ptcl_mask)
                zero_pop = count(.not.build_glob%spproj_field%included(consider_w=.false.))
                extr_score_thresh = -huge(extr_score_thresh)
                if( params_glob%l_frac_update )then
                    iextr_lim = ceiling(2.*log(real(params_glob%nptcls-zero_pop)) * (2.-params_glob%update_frac))
                    if( which_iter==1 .or.(frac_srch_space <= 99. .and. params_glob%extr_iter <= iextr_lim) )&
                        &do_extr = .true.
                else
                    iextr_lim = ceiling(2.*log(real(params_glob%nptcls-zero_pop)))
                    if( which_iter==1 .or.(frac_srch_space <= 98. .and. params_glob%extr_iter <= iextr_lim) )&
                        &do_extr = .true.
                endif
                if( do_extr )then
                    !extr_score_thresh = EXTRINITHRESH * cos(PI/2. * real(params_glob%extr_iter-1)/real(iextr_lim))
                    extr_thresh = EXTRINITHRESH * cos(PI/2. * real(params_glob%extr_iter-1)/real(iextr_lim))
                    extr_score_thresh = build_glob%spproj_field%extremal_bound(extr_thresh, 'corr')
                endif
            case('cluster','clustersym')
                if(allocated(het_mask))deallocate(het_mask)
                allocate(het_mask(params_glob%fromp:params_glob%top), source=ptcl_mask)
                zero_pop          = count(.not.build_glob%spproj_field%included(consider_w=.false.))
                extr_score_thresh = -huge(extr_score_thresh)
                if( params_glob%l_frac_update )then
                    ptcl_mask = .true.
                    iextr_lim = ceiling(2.*log(real(params_glob%nptcls-zero_pop)) * (2.-params_glob%update_frac))
                    if( which_iter==1 .or.(frac_srch_space <= 99. .and. params_glob%extr_iter <= iextr_lim) )&
                        &do_extr = .true.
                else
                    iextr_lim = ceiling(2.*log(real(params_glob%nptcls-zero_pop)))
                    if( which_iter==1 .or.(frac_srch_space <= 98. .and. params_glob%extr_iter <= iextr_lim) )&
                        &do_extr = .true.
                endif
                if( do_extr )then
                    extr_thresh = EXTRINITHRESH * cos(PI/2. * real(params_glob%extr_iter-1)/real(iextr_lim))
                    extr_score_thresh = build_glob%spproj_field%extremal_bound(extr_thresh, 'corr')
                endif
                if(trim(params_glob%refine).eq.'clustersym')then
                   ! symmetry pairing matrix
                    c1_symop = sym('c1')
                    params_glob%nspace = min( params_glob%nspace*build_glob%pgrpsyms%get_nsym(), 3000 )
                    call build_glob%eulspace%new( params_glob%nspace )
                    call build_glob%eulspace%spiral
                    call build_glob%pgrpsyms%nearest_sym_neighbors(build_glob%eulspace, symmat)
                endif
            case DEFAULT
                ! nothing to do
        end select
        if( L_BENCH ) rt_init = toc(t_init)

        ! PREPARE THE POLARFT_CORRCALC DATA STRUCTURE
        if( L_BENCH ) t_prep_pftcc = tic()
        call preppftcc4align(cline)
        if( L_BENCH ) rt_prep_pftcc = toc(t_prep_pftcc)

        write(*,'(A,1X,I3)') '>>> REFINE3D SEARCH, ITERATION:', which_iter

        ! STOCHASTIC IMAGE ALIGNMENT
        if( L_BENCH ) t_prep_primesrch3D = tic()
        ! clean big objects before starting to allocate new big memory chunks
        ! cannot kill build_glob%vol since used in continuous search
        call build_glob%vol2%kill

        ! array allocation for strategy3D
        if( DEBUG ) print *, '*** strategy3D_matcher ***: array allocation for strategy3D'
        call prep_strategy3D( ptcl_mask, npeaks )  ! allocate s3D singleton
        if( DEBUG ) print *, '*** strategy3D_matcher ***: array allocation for strategy3D, DONE'
        if( L_BENCH ) rt_prep_primesrch3D = toc(t_prep_primesrch3D)
        ! switch for per-particle polymorphic strategy3D construction
        allocate(strategy3Dsrch(params_glob%fromp:params_glob%top), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("In simple_strategy3D_matcher::refine3D_exec strategy3Dsrch",alloc_stat)
        select case(trim(params_glob%refine))
            case('snhc')
                do iptcl=params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) ) then
                        allocate(strategy3D_snhc_single :: strategy3Dsrch(iptcl)%ptr, stat=alloc_stat)
                        if(alloc_stat.ne.0)&
                            &call allocchk("In simple_strategy3D_matcher::refine3D_exec strategy3Dsrch snhc",alloc_stat)
                    endif
                end do
            case('single')
                do iptcl=params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) )then
                        updatecnt = nint(build_glob%spproj_field%get(iptcl,'updatecnt'))
                        if( .not.build_glob%spproj_field%has_been_searched(iptcl) .or. updatecnt == 1 )then
                            allocate(strategy3D_greedy_single :: strategy3Dsrch(iptcl)%ptr, stat=alloc_stat)
                        else
                            allocate(strategy3D_single        :: strategy3Dsrch(iptcl)%ptr, stat=alloc_stat)
                        endif
                        if(alloc_stat.ne.0)&
                           &call allocchk("In simple_strategy3D_matcher::refine3D_exec strategy3Dsrch single",alloc_stat)
                    endif
                end do
            case('hard_single')
                do iptcl=params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) )then
                        allocate(strategy3D_hard_single :: strategy3Dsrch(iptcl)%ptr, stat=alloc_stat)
                        if(alloc_stat.ne.0)&
                            &call allocchk("In simple_strategy3D_matcher::refine3D_exec strategy3Dsrch hard_single",alloc_stat)
                    endif
                end do
            case('greedy_single')
                do iptcl=params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) ) allocate(strategy3D_greedy_single :: strategy3Dsrch(iptcl)%ptr)
                end do
            case('cont_single')
                do iptcl=params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) )then
                        allocate(strategy3D_cont_single   :: strategy3Dsrch(iptcl)%ptr, stat=alloc_stat)
                        if(alloc_stat.ne.0)&
                            &call allocchk("In simple_strategy3D_matcher::refine3D_exec strategy3Dsrch snhc",alloc_stat)
                    endif
                end do
            case('multi')
                do iptcl=params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) )then
                        updatecnt = nint(build_glob%spproj_field%get(iptcl,'updatecnt'))
                        if( updatecnt == 1 )then
                            allocate(strategy3D_greedy_multi :: strategy3Dsrch(iptcl)%ptr, stat=alloc_stat)
                        else
                            allocate(strategy3D_multi        :: strategy3Dsrch(iptcl)%ptr, stat=alloc_stat)
                        endif
                        if(alloc_stat.ne.0)&
                            &call allocchk("In simple_strategy3D_matcher::refine3D_exec strategy3Dsrch multi",alloc_stat)
                    endif
                end do
            case('hard_multi')
                do iptcl=params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) )then
                        allocate(strategy3D_hard_multi :: strategy3Dsrch(iptcl)%ptr, stat=alloc_stat)
                        if(alloc_stat.ne.0)&
                            &call allocchk("In simple_strategy3D_matcher::refine3D_exec strategy3Dsrch hard_multi",alloc_stat)
                    endif
                end do
            case('greedy_multi')
                do iptcl=params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) )then
                        allocate(strategy3D_greedy_multi  :: strategy3Dsrch(iptcl)%ptr, stat=alloc_stat)
                        if(alloc_stat.ne.0)&
                            &call allocchk("In simple_strategy3D_matcher::refine3D_exec strategy3Dsrch snhc",alloc_stat)
                    endif
                end do
            case('cluster','clustersym')
                do iptcl=params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) )then
                        allocate(strategy3D_cluster       :: strategy3Dsrch(iptcl)%ptr, stat=alloc_stat)
                        if(alloc_stat.ne.0)&
                             call allocchk("In simple_strategy3D_matcher::refine3D_exec strategy3Dsrch cluster",alloc_stat)
                    endif
                end do
            case('cluster_snhc')
                do iptcl=params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) )then
                        allocate(strategy3D_cluster_snhc :: strategy3Dsrch(iptcl)%ptr, stat=alloc_stat)
                        if(alloc_stat.ne.0)call allocchk("In simple_strategy3D_matcher::refine3D_exec strategy3Dsrch cluster_neigh",alloc_stat)
                    endif
                end do
            case DEFAULT
                write(*,*) 'refine flag: ', trim(params_glob%refine)
                stop 'Refinement mode unsupported'
        end select
        ! construct search objects
        cnt = 0
        do iptcl=params_glob%fromp,params_glob%top
            if( ptcl_mask(iptcl) )then
                cnt = cnt + 1
                ! search spec
                strategy3Dspec%iptcl     =  iptcl
                strategy3Dspec%iptcl_map =  cnt
                strategy3Dspec%szsn      =  params_glob%szsn
                strategy3Dspec%extr_score_thresh =  extr_score_thresh
                if( allocated(het_mask) )  strategy3Dspec%do_extr =  het_mask(iptcl)
                if( allocated(symmat) )    strategy3Dspec%symmat  => symmat
                ! search object
                call strategy3Dsrch(iptcl)%ptr%new(strategy3Dspec, npeaks)
            endif
        end do
        if( DEBUG ) print *, '*** strategy3D_matcher ***: search object construction, DONE'
        ! memoize CTF matrices
        if( trim(params_glob%oritype) .eq. 'ptcl3D' )then
            if( build_glob%spproj%get_ctfflag('ptcl3D').ne.'no' )&
                &call pftcc%create_polar_absctfmats(build_glob%spproj, 'ptcl3D')
        else
            ! class averages have no CTF
        endif
        ! memoize FFTs
        call pftcc%memoize_ffts

        ! SEARCH
        if( params_glob%dryrun.eq.'yes' )then
            ! TBD
        else
            if( L_BENCH ) t_align = tic()
            !$omp parallel do default(shared) private(i,iptcl) schedule(static) proc_bind(close)
            do i=1,nptcls2update
                iptcl = pinds(i)
                call strategy3Dsrch(iptcl)%ptr%srch
            end do
            !$omp end parallel do
        endif
        ! clean
        call clean_strategy3D  ! deallocate s3D singleton
        call pftcc%kill
        call build_glob%vol%kill
        call build_glob%vol_odd%kill
        do iptcl = params_glob%fromp,params_glob%top
            if( ptcl_mask(iptcl) )then
                call strategy3Dsrch(iptcl)%ptr%kill
                strategy3Dsrch(iptcl)%ptr => null()
            endif
        end do
        deallocate(strategy3Dsrch)
        if( L_BENCH ) rt_align = toc(t_align)
        if( allocated(symmat)   ) deallocate(symmat)
        if( allocated(het_mask) ) deallocate(het_mask)

        ! OUTPUT ORIENTATIONS
        select case(trim(params_glob%oritype))
            case('ptcl3D')
                call binwrite_oritab(params_glob%outfile, build_glob%spproj, &
                    &build_glob%spproj_field, [params_glob%fromp,params_glob%top], isegment=PTCL3D_SEG)
            case('cls3D')
                call binwrite_oritab(params_glob%outfile, build_glob%spproj, &
                    &build_glob%spproj_field, [params_glob%fromp,params_glob%top], isegment=CLS3D_SEG)
            case DEFAULT
                write(*,*) 'oritype: ', trim(params_glob%oritype)
                stop 'Unsupported oritype; strategy3D_matcher :: refine3D_exec'
        end select
        params_glob%oritab = params_glob%outfile

        ! VOLUMETRIC 3D RECONSTRUCTION
        if( L_BENCH ) t_rec = tic()
        ! make the gridding prepper
        if( params_glob%eo .ne. 'no' )then
            kbwin = build_glob%eorecvols(1)%get_kbwin()
        else
            kbwin = build_glob%recvols(1)%get_kbwin()
        endif
        ! init volumes
        call preprecvols
        ! prep rec imgs
        allocate(rec_imgs(MAXIMGBATCHSZ))
        do i=1,MAXIMGBATCHSZ
            call rec_imgs(i)%new([params_glob%boxpd, params_glob%boxpd, 1], params_glob%smpd)
        end do
        ! prep batch imgs
        call prepimgbatch( MAXIMGBATCHSZ)
        ! make mask image
        call mskimg%disc([params_glob%box,params_glob%box,1], params_glob%smpd, params_glob%msk)
        ! gridding batch loop
        do i_batch=1,nptcls2update,MAXIMGBATCHSZ
            batchlims = [i_batch,min(nptcls2update,i_batch + MAXIMGBATCHSZ - 1)]
            call read_imgbatch( nptcls2update, pinds, batchlims)
            !$omp parallel do default(shared) private(i,ibatch) schedule(static) proc_bind(close)
            do i=batchlims(1),batchlims(2)
                ibatch = i - batchlims(1) + 1
                call build_glob%imgbatch(ibatch)%norm_subtr_backgr_pad_serial(mskimg, rec_imgs(ibatch))
            end do
            !$omp end parallel do
            ! gridding
            do i=batchlims(1),batchlims(2)
                iptcl       = pinds(i)
                ibatch      = i - batchlims(1) + 1
                orientation = build_glob%spproj_field%get_ori(iptcl)
                ctfvars     = build_glob%spproj%get_ctfparams(params_glob%oritype, iptcl)
                if( orientation%isstatezero() ) cycle
                select case(trim(params_glob%refine))
                    case('clustersym')
                        ! always C1 reconstruction
                        call grid_ptcl( rec_imgs(ibatch), c1_symop, orientation, s3D%o_peaks(iptcl), ctfvars)
                    case('hard_multi', 'hard_single')
                        ! npeaks > 1 but only one is selected for reconstruction by probabilistic logic
                        call grid_ptcl( rec_imgs(ibatch), build_glob%pgrpsyms, orientation, ctfvars)
                    case DEFAULT
                        call grid_ptcl( rec_imgs(ibatch), build_glob%pgrpsyms, orientation, s3D%o_peaks(iptcl), ctfvars)
                end select
            end do
        end do
        ! normalise structure factors
        if( params_glob%eo .ne. 'no' )then
            call eonorm_struct_facts( cline, which_iter)
        else
            call norm_struct_facts( which_iter)
        endif
        ! destruct
        call killrecvols()
        do ibatch=1,MAXIMGBATCHSZ
            call rec_imgs(ibatch)%kill
            call build_glob%imgbatch(ibatch)%kill
        end do
        deallocate(rec_imgs, build_glob%imgbatch)
        call mskimg%kill
        if( L_BENCH ) rt_rec = toc(t_rec)

        ! REPORT CONVERGENCE
        call qsys_job_finished(  'simple_strategy3D_matcher :: refine3D_exec')
        if( .not. params_glob%l_distr_exec ) converged = conv%check_conv3D(cline, params_glob%msk)
        if( L_BENCH )then
            rt_tot  = toc(t_tot)
            doprint = .true.
            if( params_glob%part /= 1 ) doprint = .false.
            if( doprint )then
                benchfname = 'HADAMARD3D_BENCH_ITER'//int2str_pad(which_iter,3)//'.txt'
                call fopen(fnr, FILE=trim(benchfname), STATUS='REPLACE', action='WRITE')
                write(fnr,'(a)') '*** TIMINGS (s) ***'
                write(fnr,'(a,1x,f9.2)') 'initialisation          : ', rt_init
                write(fnr,'(a,1x,f9.2)') 'pftcc preparation       : ', rt_prep_pftcc
                write(fnr,'(a,1x,f9.2)') 'primesrch3D preparation : ', rt_prep_primesrch3D
                write(fnr,'(a,1x,f9.2)') 'stochastic alignment    : ', rt_align
                write(fnr,'(a,1x,f9.2)') 'reconstruction          : ', rt_rec
                write(fnr,'(a,1x,f9.2)') 'total time              : ', rt_tot
                write(fnr,'(a)') ''
                write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
                write(fnr,'(a,1x,f9.2)') 'initialisation          : ', (rt_init/rt_tot)             * 100.
                write(fnr,'(a,1x,f9.2)') 'pftcc preparation       : ', (rt_prep_pftcc/rt_tot)       * 100.
                write(fnr,'(a,1x,f9.2)') 'primesrch3D preparation : ', (rt_prep_primesrch3D/rt_tot) * 100.
                write(fnr,'(a,1x,f9.2)') 'stochastic alignment    : ', (rt_align/rt_tot)            * 100.
                write(fnr,'(a,1x,f9.2)') 'reconstruction          : ', (rt_rec/rt_tot)              * 100.
                write(fnr,'(a,1x,f9.2)') '% accounted for         : ',&
                    &((rt_init+rt_prep_pftcc+rt_prep_primesrch3D+rt_align+rt_rec)/rt_tot) * 100.
                call fclose(fnr)
            endif
        endif
    end subroutine refine3D_exec

    !> Prepare alignment search using polar projection Fourier cross correlation
    subroutine preppftcc4align( cline )
        use simple_polarizer,             only: polarizer
        use simple_cmdline,               only: cmdline
        use simple_strategy2D3D_common,   only: calcrefvolshift_and_mapshifts2ptcls, preprefvol, build_pftcc_particles
        class(cmdline), intent(inout) :: cline !< command line
        type(polarizer),  allocatable :: match_imgs(:)
        real      :: xyz(3)
        integer   :: cnt, s, ind, iref, nrefs, imatch
        logical   :: do_center, has_been_searched
        nrefs             = params_glob%nspace * params_glob%nstates
        has_been_searched = .not.build_glob%spproj%is_virgin_field(params_glob%oritype)
        ! must be done here since params_glob%kfromto is dynamically set
        if( params_glob%eo .ne. 'no' )then
            call pftcc%new(nrefs, ptcl_mask, nint(build_glob%spproj_field%get_all('eo', [params_glob%fromp,params_glob%top])))
        else
            call pftcc%new(nrefs, ptcl_mask)
        endif

        ! PREPARATION OF REFERENCES IN PFTCC
        ! read reference volumes and create polar projections
        cnt = 0
        do s=1,params_glob%nstates
            if( has_been_searched )then
                if( build_glob%spproj_field%get_pop(s, 'state') == 0 )then
                    ! empty state
                    cnt = cnt + params_glob%nspace
                    call progress(cnt, nrefs)
                    cycle
                endif
            endif
            call calcrefvolshift_and_mapshifts2ptcls( cline, s, params_glob%vols(s), do_center, xyz)
            if( params_glob%eo .ne. 'no' )then
                if( params_glob%nstates.eq.1 )then
                    ! PREPARE ODD REFERENCES
                    call preprefvol(cline, s, params_glob%vols_odd(s), do_center, xyz)
                    !$omp parallel do default(shared) private(iref) schedule(static) proc_bind(close)
                    do iref=1,params_glob%nspace
                        call build_glob%vol%fproject_polar((s - 1) * params_glob%nspace + iref, &
                            &build_glob%eulspace%get_ori(iref), pftcc, iseven=.false.)
                    end do
                    !$omp end parallel do
                    ! copy odd volume
                    build_glob%vol_odd = build_glob%vol
                    ! expand for fast interpolation
                    call build_glob%vol_odd%expand_cmat(params_glob%alpha)
                    ! PREPARE EVEN REFERENCES
                    call preprefvol( cline, s, params_glob%vols_even(s), do_center, xyz)
                    !$omp parallel do default(shared) private(iref) schedule(static) proc_bind(close)
                    do iref=1,params_glob%nspace
                        call build_glob%vol%fproject_polar((s - 1) * params_glob%nspace + iref, &
                            &build_glob%eulspace%get_ori(iref), pftcc, iseven=.true.)
                    end do
                    !$omp end parallel do
                else
                    call preprefvol( cline, s, params_glob%vols(s), do_center, xyz)
                    !$omp parallel do default(shared) private(iref, ind) schedule(static) proc_bind(close)
                    do iref=1,params_glob%nspace
                        ind = (s - 1) * params_glob%nspace + iref
                        call build_glob%vol%fproject_polar(ind, build_glob%eulspace%get_ori(iref), &
                            &pftcc, iseven=.true.)
                        call pftcc%cp_even2odd_ref(ind)
                    end do
                    !$omp end parallel do
                endif
            else
                ! low-pass set or multiple states
                call preprefvol( cline, s, params_glob%vols(s), do_center, xyz)
                !$omp parallel do default(shared) private(iref) schedule(static) proc_bind(close)
                do iref=1,params_glob%nspace
                    call build_glob%vol%fproject_polar((s - 1) * params_glob%nspace + iref, &
                        &build_glob%eulspace%get_ori(iref), pftcc, iseven=.true.)
                end do
                !$omp end parallel do
            endif
        end do

        ! PREPARATION OF PARTICLES IN PFTCC
        ! prepare the polarizer images
        call build_glob%img_match%init_polarizer(pftcc, params_glob%alpha)
        allocate(match_imgs(MAXIMGBATCHSZ))
        do imatch=1,MAXIMGBATCHSZ
            call match_imgs(imatch)%new([params_glob%boxmatch, params_glob%boxmatch, 1], params_glob%smpd)
            call match_imgs(imatch)%copy_polarizer(build_glob%img_match)
        end do
        call build_pftcc_particles(pftcc, MAXIMGBATCHSZ, match_imgs, .true., ptcl_mask)

        ! DESTRUCT
        do imatch=1,MAXIMGBATCHSZ
            call match_imgs(imatch)%kill_polarizer
            call match_imgs(imatch)%kill
            call build_glob%imgbatch(imatch)%kill
        end do
        deallocate(match_imgs, build_glob%imgbatch)
        if( DEBUG ) print *, '*** strategy3D_matcher ***: finished preppftcc4align'
    end subroutine preppftcc4align

end module simple_strategy3D_matcher
