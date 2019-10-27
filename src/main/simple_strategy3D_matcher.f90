! projection-matching based on Hadamard products, high-level search routines for REFINE3D
module simple_strategy3D_matcher
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_strategy3D_alloc ! singleton s3D
use simple_timer
use simple_o_peaks_io
use simple_oris,                     only: oris
use simple_qsys_funs,                only: qsys_job_finished
use simple_binoris_io,               only: binwrite_oritab
use simple_kbinterpol,               only: kbinterpol
use simple_ori,                      only: ori
use simple_sym,                      only: sym
use simple_image,                    only: image
use simple_cmdline,                  only: cmdline
use simple_parameters,               only: params_glob
use simple_builder,                  only: build_glob
use simple_polarizer,                only: polarizer
use simple_polarft_corrcalc,         only: polarft_corrcalc
use simple_strategy2D3D_common,      only: killrecvols, set_bp_range, preprecvols,&
    prepimgbatch, grid_ptcl, read_imgbatch, norm_struct_facts
use simple_strategy3D_cluster,       only: strategy3D_cluster
use simple_strategy3D_clustersoft,   only: strategy3D_clustersoft
use simple_strategy3D_single,        only: strategy3D_single
use simple_strategy3D_multi,         only: strategy3D_multi
use simple_strategy3D_snhc_single,   only: strategy3D_snhc_single
use simple_strategy3D_greedy_single, only: strategy3D_greedy_single
use simple_strategy3D_greedy_multi,  only: strategy3D_greedy_multi
use simple_strategy3D_neigh_single,  only: strategy3D_neigh_single
use simple_strategy3D_neigh_multi,   only: strategy3D_neigh_multi
use simple_strategy3D_cont_single,   only: strategy3D_cont_single
use simple_strategy3D,               only: strategy3D
use simple_strategy3D_srch,          only: strategy3D_spec, set_ptcl_stats, eval_ptcl
use simple_convergence,              only: convergence
use simple_euclid_sigma2,            only: euclid_sigma2
implicit none

public :: refine3D_exec, preppftcc4align, pftcc, setup_weights_read_o_peaks
public :: calc_global_ori_weights, calc_3Drec, calc_proj_weights
private
#include "simple_local_flags.inc"

logical, parameter              :: L_BENCH = .false., DEBUG_HERE = .false.
logical                         :: has_been_searched
type(polarft_corrcalc),  target :: pftcc
type(polarizer),    allocatable :: match_imgs(:)
integer,            allocatable :: pinds(:), prev_states(:)
logical,            allocatable :: ptcl_mask(:)
type(sym)                       :: c1_symop
integer                         :: nptcls2update
integer                         :: npeaks
integer(timer_int_kind)         :: t_init, t_prep_pftcc, t_align, t_rec, t_tot, t_prep_primesrch3D
real(timer_int_kind)            :: rt_init, rt_prep_pftcc, rt_align, rt_rec, rt_prep_primesrch3D
real(timer_int_kind)            :: rt_tot
character(len=STDLEN)           :: benchfname
type(euclid_sigma2)             :: eucl_sigma

contains

    subroutine refine3D_exec( cline, which_iter, converged )
        class(cmdline),        intent(inout) :: cline
        integer,               intent(in)    :: which_iter
        logical,               intent(inout) :: converged
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
        type(convergence)     :: conv
        type(oris)            :: o_peak_prev
        real, allocatable     :: resarr(:)
        real    :: frac_srch_space, extr_thresh, extr_score_thresh, mi_proj, anneal_ratio
        integer :: iptcl, i, fnr, ithr, updatecnt, state, n_nozero
        integer :: ibatch, iextr_lim, lpind_anneal, lpind_start
        logical :: doprint, do_extr

        if( L_BENCH )then
            t_init = tic()
            t_tot  = t_init
        endif

        ! CHECK THAT WE HAVE AN EVEN/ODD PARTITIONING
        if( build_glob%spproj_field%get_nevenodd() == 0 ) THROW_HARD('no eo partitioning available; refine3D_exec')

        ! CHECK WHETHER WE HAVE PREVIOUS 3D ORIENTATIONS
        has_been_searched = .not.build_glob%spproj%is_virgin_field(params_glob%oritype)

        ! SET FOURIER INDEX RANGE
        call set_bp_range(cline)

        ! DETERMINE THE NUMBER OF PEAKS
        select case(params_glob%refine)
            case('cluster', 'snhc', 'clustersym', 'cont_single', 'eval')
                npeaks = 1
            case DEFAULT
                npeaks = NPEAKS2REFINE
                ! command line overrides
                if( cline%defined('npeaks') ) npeaks = params_glob%npeaks
        end select
        if( DEBUG_HERE ) write(logfhandle,*) '*** strategy3D_matcher ***: determined the number of peaks'

        ! SET FRACTION OF SEARCH SPACE
        frac_srch_space = build_glob%spproj_field%get_avg('frac')

        ! SETUP WEIGHTS
        if( params_glob%l_ptclw )then
            call build_glob%spproj_field%calc_soft_weights(params_glob%frac)
        else
            call build_glob%spproj_field%calc_hard_weights(params_glob%frac)
        endif

        ! READ FOURIER RING CORRELATIONS
        if( params_glob%nstates.eq.1 )then
            if( file_exists(params_glob%frcs) ) call build_glob%projfrcs%read(params_glob%frcs)
        else
            if( file_exists(CLUSTER3D_FRCS) )then
                call build_glob%projfrcs%read(CLUSTER3D_FRCS)
            else
                if( file_exists(params_glob%frcs) )then
                    call build_glob%projfrcs%read(params_glob%frcs)
                endif
            endif
        endif

        ! PARTICLE INDEX SAMPLING FOR FRACTIONAL UPDATE (OR NOT)
        if( allocated(pinds) )     deallocate(pinds)
        if( allocated(ptcl_mask) ) deallocate(ptcl_mask)
        if( params_glob%l_frac_update )then
            allocate(ptcl_mask(params_glob%fromp:params_glob%top))
            call build_glob%spproj_field%sample4update_and_incrcnt([params_glob%fromp,params_glob%top],&
            &params_glob%update_frac, nptcls2update, pinds, ptcl_mask)
        else
            allocate(ptcl_mask(params_glob%fromp:params_glob%top))
            call build_glob%spproj_field%sample4update_and_incrcnt_nofrac([params_glob%fromp,params_glob%top],&
            nptcls2update, pinds, ptcl_mask)
        endif

        ! EXTREMAL LOGICS
        do_extr           = .false.
        extr_score_thresh = -huge(extr_score_thresh)
        select case(trim(params_glob%refine))
            case('cluster','clustersym','clustersoft')
                ! general logics
                if(allocated(het_mask))deallocate(het_mask)
                allocate(het_mask(params_glob%fromp:params_glob%top), source=ptcl_mask)
                call build_glob%spproj_field%set_extremal_vars(params_glob%extr_init, params_glob%extr_iter,&
                    &which_iter, frac_srch_space, do_extr, iextr_lim, update_frac=params_glob%update_frac)
                if( do_extr )then
                    anneal_ratio      = max(0., cos(PI/2.*real(params_glob%extr_iter-1)/real(iextr_lim)))
                    extr_thresh       = params_glob%extr_init * anneal_ratio
                    extr_score_thresh = build_glob%spproj_field%extremal_bound(extr_thresh, 'corr')
                    if( cline%defined('lpstart') )then
                        ! resolution limit update
                        lpind_start       = calc_fourier_index(params_glob%lpstart,params_glob%boxmatch,params_glob%smpd)
                        lpind_anneal      = nint(real(lpind_start) + (1.-anneal_ratio)*real(params_glob%kstop-lpind_start))
                        params_glob%kstop = min(lpind_anneal, params_glob%kstop)
                        resarr            = build_glob%img%get_res()
                        params_glob%lp    = resarr(params_glob%kstop)
                        if( params_glob%cc_objfun .ne. OBJFUN_EUCLID ) params_glob%kfromto(2) = params_glob%kstop
                        call build_glob%spproj_field%set_all2single('lp',params_glob%lp)
                        deallocate(resarr)
                    endif
                else
                    het_mask = .false.
                endif
                ! refinement mode specifics
                select case(trim(params_glob%refine))
                case('clustersym')
                   ! symmetry pairing matrix
                    c1_symop = sym('c1')
                    params_glob%nspace = min(params_glob%nspace*build_glob%pgrpsyms%get_nsym(), 2500)
                    call build_glob%eulspace%new( params_glob%nspace )
                    call build_glob%eulspace%spiral
                    call build_glob%pgrpsyms%nearest_sym_neighbors(build_glob%eulspace, symmat)
                case('clustersoft')
                    prev_states = nint(build_glob%spproj_field%get_all('state',[params_glob%fromp,params_glob%top]))
                end select
        end select
        if( L_BENCH ) rt_init = toc(t_init)

        ! PREPARE THE POLARFT_CORRCALC DATA STRUCTURE
        if( L_BENCH ) t_prep_pftcc = tic()
        call preppftcc4align(cline)
        if( L_BENCH ) rt_prep_pftcc = toc(t_prep_pftcc)

        ! STOCHASTIC IMAGE ALIGNMENT
        if( L_BENCH ) t_prep_primesrch3D = tic()
        write(logfhandle,'(A,1X,I3)') '>>> REFINE3D SEARCH, ITERATION:', which_iter
        ! clean big objects before starting to allocate new big memory chunks
        ! cannot kill build_glob%vol since used in continuous search
        call build_glob%vol2%kill
        ! array allocation for strategy3D
        if( DEBUG_HERE ) write(logfhandle,*) '*** strategy3D_matcher ***: array allocation for strategy3D'
        call prep_strategy3D( ptcl_mask, npeaks ) ! allocate s3D singleton
        if( DEBUG_HERE ) write(logfhandle,*) '*** strategy3D_matcher ***: array allocation for strategy3D, DONE'
        if( L_BENCH ) rt_prep_primesrch3D = toc(t_prep_primesrch3D)
        ! switch for per-particle polymorphic strategy3D construction
        allocate(strategy3Dsrch(params_glob%fromp:params_glob%top), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("In simple_strategy3D_matcher::refine3D_exec strategy3Dsrch",alloc_stat)
        select case(trim(params_glob%refine))
            case('snhc')
                do iptcl=params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) ) allocate(strategy3D_snhc_single :: strategy3Dsrch(iptcl)%ptr)
                end do
            case('single')
                do iptcl=params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) )then
                        if( .not.build_glob%spproj_field%has_been_searched(iptcl) .or. ran3() < GREEDY_FREQ )then
                            allocate(strategy3D_greedy_single :: strategy3Dsrch(iptcl)%ptr)
                        else
                            allocate(strategy3D_single :: strategy3Dsrch(iptcl)%ptr)
                        endif
                    endif
                end do
            case('greedy_single')
                do iptcl=params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) ) allocate(strategy3D_greedy_single :: strategy3Dsrch(iptcl)%ptr)
                end do
            case('cont_single')
                do iptcl=params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) ) allocate(strategy3D_cont_single :: strategy3Dsrch(iptcl)%ptr)
                end do
            case('multi')
                do iptcl=params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) )then
                        updatecnt = nint(build_glob%spproj_field%get(iptcl,'updatecnt'))
                        if( .not.build_glob%spproj_field%has_been_searched(iptcl) .or. updatecnt == 1 )then
                            allocate(strategy3D_greedy_multi :: strategy3Dsrch(iptcl)%ptr)
                        else
                            allocate(strategy3D_multi :: strategy3Dsrch(iptcl)%ptr)
                        endif
                    endif
                end do
            case('greedy_multi')
                do iptcl=params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) ) allocate(strategy3D_greedy_multi :: strategy3Dsrch(iptcl)%ptr)
                end do
            case('cluster','clustersym')
                do iptcl=params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) ) allocate(strategy3D_cluster :: strategy3Dsrch(iptcl)%ptr)
                end do
            case('clustersoft')
                call open_o_peaks_io(trim(params_glob%o_peaks_file))
                do iptcl=params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) )then
                        allocate(strategy3D_clustersoft :: strategy3Dsrch(iptcl)%ptr)
                        call read_o_peak(s3D%o_peaks(iptcl), [params_glob%fromp,params_glob%top], iptcl, n_nozero)
                    endif
                end do
                call close_o_peaks_io
            case('neigh_single')
                do iptcl=params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) ) allocate(strategy3D_neigh_single :: strategy3Dsrch(iptcl)%ptr)
                end do
            case('neigh_multi')
                do iptcl=params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) ) allocate(strategy3D_neigh_multi :: strategy3Dsrch(iptcl)%ptr)
                end do
            case('eval')
                ! nothing to do
            case DEFAULT
                THROW_HARD('refinement mode: '//trim(params_glob%refine)//' unsupported')
        end select
        ! construct search objects
        select case(trim(params_glob%refine))
            case('eval')
                ! nothing to do
            case DEFAULT
                do iptcl=params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) )then
                        ! search spec
                        strategy3Dspec%iptcl =  iptcl
                        strategy3Dspec%szsn  =  params_glob%szsn
                        strategy3Dspec%extr_score_thresh = extr_score_thresh
                        if( allocated(het_mask) ) strategy3Dspec%do_extr =  het_mask(iptcl)
                        if( allocated(symmat) )   strategy3Dspec%symmat  => symmat
                        ! search object
                        call strategy3Dsrch(iptcl)%ptr%new(strategy3Dspec, npeaks)
                    endif
                end do
        end select
        if( DEBUG_HERE ) write(logfhandle,*) '*** strategy3D_matcher ***: search object construction, DONE'
        ! memoize CTF matrices
        if( trim(params_glob%oritype) .eq. 'ptcl3D' )then
            if( build_glob%spproj%get_ctfflag('ptcl3D').ne.'no' )&
                &call pftcc%create_polar_absctfmats(build_glob%spproj, 'ptcl3D')
        else
            ! class averages have no CTF
        endif
        ! memoize FFTs
        call pftcc%memoize_ffts
        ! read o_peaks for neigh refinement modes
        if( str_has_substr(params_glob%refine, 'neigh') ) call read_o_peaks

        ! SEARCH
        select case(trim(params_glob%refine))
            case('eval')
                !$omp parallel do default(shared) private(i) schedule(static) proc_bind(close)
                do i=1,nptcls2update
                    call eval_ptcl(pftcc, pinds(i), pinds(i))
                end do
                !$omp end parallel do
            case DEFAULT
                if( L_BENCH ) t_align = tic()
                !$omp parallel do default(shared) private(i,iptcl,ithr) schedule(static) proc_bind(close)
                do i=1,nptcls2update
                    iptcl = pinds(i)
                    ithr  = omp_get_thread_num() + 1
                    call strategy3Dsrch(iptcl)%ptr%srch(ithr)
                end do
                !$omp end parallel do
        end select

        ! CALCULATE AND WRITE SIGMAS FOR ML-BASED REFINEMENT
        if ( params_glob%l_needs_sigma ) then
            call eucl_sigma%calc_and_write_sigmas(build_glob%spproj_field, s3D%o_peaks, ptcl_mask)
            call eucl_sigma%kill_ptclsigma2
        end if

        ! UPDATE PARTICLE STATS
        call calc_ptcl_stats

        ! O_PEAKS I/O & CONVERGENCE STATS
        ! here we read all peaks to allow deriving statistics based on the complete set
        ! this is needed for deriving projection direction weights
        select case(trim(params_glob%refine))
        case('eval','cluster','clustersym') ! 'neigh_multi','neigh_single'
                ! nothing to do
            case DEFAULT
                if( .not. file_exists(trim(params_glob%o_peaks_file)) )then
                    ! write an empty one to be filled in
                    call write_empty_o_peaks_file(params_glob%o_peaks_file, [params_glob%fromp,params_glob%top])
                endif
                call open_o_peaks_io(trim(params_glob%o_peaks_file))
                do iptcl=params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) )then
                        state = build_glob%spproj_field%get_state(iptcl)
                        call read_o_peak(o_peak_prev, [params_glob%fromp,params_glob%top], iptcl, n_nozero)
                        if( n_nozero == 0 )then
                            ! there's nothing to compare, set overlap to zero
                            call build_glob%spproj_field%set(iptcl, 'mi_proj', 0.)
                        else
                            mi_proj = s3D%o_peaks(iptcl)%overlap(o_peak_prev, 'proj', state)
                            call build_glob%spproj_field%set(iptcl, 'mi_proj', mi_proj)
                        endif
                        ! replace the peak on disc
                        call write_o_peak(s3D%o_peaks(iptcl), [params_glob%fromp,params_glob%top], iptcl)
                    else
                        call read_o_peak(s3D%o_peaks(iptcl), [params_glob%fromp,params_glob%top], iptcl, n_nozero)
                    endif
                end do
                call close_o_peaks_io
        end select

        ! CALCULATE GLOBAL ORIENTATION WEIGHTS
        if( params_glob%l_wglob ) call calc_global_ori_weights

        ! CALCULATE PROJECTION DIRECTION WEIGHTS
        ! call calc_proj_weights !!!!!!!!!! turned off 4 now, needs integration and testing

        ! CLEAN
        call clean_strategy3D ! deallocate s3D singleton
        call pftcc%kill
        call build_glob%vol%kill
        call build_glob%vol_odd%kill
        ! if( params_glob%cc_objfun /= OBJFUN_EUCLID ) call eucl_sigma%kill
        select case(trim(params_glob%refine))
            case('eval')
                ! nothing to do
            case DEFAULT
                do iptcl = params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) )then
                        call strategy3Dsrch(iptcl)%ptr%kill
                        strategy3Dsrch(iptcl)%ptr => null()
                    endif
                end do
        end select
        do ibatch=1,MAXIMGBATCHSZ
            call match_imgs(ibatch)%kill_polarizer
            call match_imgs(ibatch)%kill
        end do
        deallocate(strategy3Dsrch, match_imgs)
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
                THROW_HARD('unsupported oritype: '//trim(params_glob%oritype)//'; refine3D_exec')
        end select
        params_glob%oritab = params_glob%outfile

        ! VOLUMETRIC 3D RECONSTRUCTION
        call calc_3Drec( cline, which_iter )
        call eucl_sigma%kill

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
        use simple_polarizer,           only: polarizer
        use simple_cmdline,             only: cmdline
        use simple_strategy2D3D_common, only: calcrefvolshift_and_mapshifts2ptcls, preprefvol, build_pftcc_particles
        class(cmdline), intent(inout) :: cline !< command line
        type(ori) :: o_tmp
        real      :: xyz(3)
        integer   :: cnt, s, ind, iref, nrefs, imatch
        logical   :: do_center
        nrefs = params_glob%nspace * params_glob%nstates
        ! must be done here since params_glob%kfromto is dynamically set
        if( params_glob%l_lpset )then
            call pftcc%new(nrefs, [params_glob%fromp,params_glob%top], ptcl_mask)
        else
            call pftcc%new(nrefs, [params_glob%fromp,params_glob%top], ptcl_mask,&
                &nint(build_glob%spproj_field%get_all('eo', [params_glob%fromp,params_glob%top])))
        endif
        if ( params_glob%l_needs_sigma ) then
            call eucl_sigma%new(SIGMA2_FBODY//int2str_pad(params_glob%part,params_glob%numlen)//'.dat')
            call eucl_sigma%read(build_glob%spproj_field, ptcl_mask)
        end if
        ! PREPARATION OF REFERENCES IN PFTCC
        ! read reference volumes and create polar projections
        cnt = 0
        do s=1,params_glob%nstates
            if( str_has_substr(params_glob%refine,'greedy') )then
                if( .not.file_exists(params_glob%vols(s)) )then
                    cnt = cnt + params_glob%nspace
                    call progress(cnt, nrefs)
                    cycle
                endif
            else
                if( has_been_searched )then
                    if( build_glob%spproj_field%get_pop(s, 'state') == 0 )then
                        ! empty state
                        cnt = cnt + params_glob%nspace
                        call progress(cnt, nrefs)
                        cycle
                    endif
                endif
            endif
            call calcrefvolshift_and_mapshifts2ptcls( cline, s, params_glob%vols(s), do_center, xyz)
            if( params_glob%l_lpset )then
                ! low-pass set or multiple states
                call preprefvol(pftcc, cline, s, params_glob%vols(s), do_center, xyz, .true.)
                !$omp parallel do default(shared) private(iref, o_tmp) schedule(static) proc_bind(close)
                do iref=1,params_glob%nspace
                    call build_glob%eulspace%get_ori(iref, o_tmp)
                    call build_glob%vol%fproject_polar((s - 1) * params_glob%nspace + iref, &
                        &o_tmp, pftcc, iseven=.true., mask=build_glob%l_resmsk)
                    call o_tmp%kill
                end do
                !$omp end parallel do
            else
                if( params_glob%nstates.eq.1 )then
                    ! PREPARE ODD REFERENCES
                    call preprefvol(pftcc, cline, s, params_glob%vols_odd(s), do_center, xyz, .false.)
                    !$omp parallel do default(shared) private(iref, o_tmp) schedule(static) proc_bind(close)
                    do iref=1,params_glob%nspace
                        call build_glob%eulspace%get_ori(iref, o_tmp)
                        call build_glob%vol%fproject_polar((s - 1) * params_glob%nspace + iref, &
                            &o_tmp, pftcc, iseven=.false., mask=build_glob%l_resmsk)
                        call o_tmp%kill
                    end do
                    !$omp end parallel do
                    ! copy odd volume
                    build_glob%vol_odd = build_glob%vol
                    ! expand for fast interpolation
                    call build_glob%vol_odd%expand_cmat(params_glob%alpha,norm4proj=.true.)
                    ! PREPARE EVEN REFERENCES
                    call preprefvol(pftcc,  cline, s, params_glob%vols_even(s), do_center, xyz, .true.)
                    !$omp parallel do default(shared) private(iref, o_tmp) schedule(static) proc_bind(close)
                    do iref=1,params_glob%nspace
                        call build_glob%eulspace%get_ori(iref, o_tmp)
                        call build_glob%vol%fproject_polar((s - 1) * params_glob%nspace + iref, &
                            &o_tmp, pftcc, iseven=.true., mask=build_glob%l_resmsk)
                        call o_tmp%kill
                    end do
                    !$omp end parallel do
                else
                    call preprefvol(pftcc, cline, s, params_glob%vols(s), do_center, xyz, .true.)
                    !$omp parallel do default(shared) private(iref, ind, o_tmp) schedule(static) proc_bind(close)
                    do iref=1,params_glob%nspace
                        ind = (s - 1) * params_glob%nspace + iref
                        call build_glob%eulspace%get_ori(iref, o_tmp)
                        call build_glob%vol%fproject_polar(ind, o_tmp, pftcc, iseven=.true., mask=build_glob%l_resmsk)
                        call pftcc%cp_even2odd_ref(ind)
                        call o_tmp%kill
                    end do
                    !$omp end parallel do
                endif
            endif
        end do
        if( params_glob%l_needs_sigma .and. params_glob%cc_objfun /= OBJFUN_EUCLID ) then
            ! When calculating sigma2 prior to OBJFUN_EUCLID the refernces are zeroed out
            !  beyond the resolution limit such that sigma2 is the weighted sum of particle power spectrum
            call pftcc%zero_refs_beyond_kstop
        endif
        ! PREPARATION OF PARTICLES IN PFTCC
        ! prepare the polarizer images
        call build_glob%img_match%init_polarizer(pftcc, params_glob%alpha)
        allocate(match_imgs(MAXIMGBATCHSZ))
        do imatch=1,MAXIMGBATCHSZ
            call match_imgs(imatch)%new([params_glob%boxmatch, params_glob%boxmatch, 1], params_glob%smpd)
            call match_imgs(imatch)%copy_polarizer(build_glob%img_match)
        end do
        call build_pftcc_particles(pftcc, MAXIMGBATCHSZ, match_imgs, ptcl_mask)
        if( DEBUG_HERE ) write(logfhandle,*) '*** strategy3D_matcher ***: finished preppftcc4align'
    end subroutine preppftcc4align

    !> Prepare alignment search using polar projection Fourier cross correlation
    subroutine calc_ptcl_stats
        use simple_strategy2D3D_common, only: prepimg4align
        type(polarft_corrcalc) :: pftcc_here
        logical, allocatable   :: ptcl_mask_not(:)
        integer :: nptcls, nrefs, iptcl_batch, batchlims(2), iptcl, imatch, eoarr(MAXIMGBATCHSZ)
        if( .not.params_glob%l_frac_update ) return
        select case(params_glob%refine)
        case('cluster','clustersym','clustersoft','eval')
                return
            case DEFAULT
                ! all good
        end select
        ! init local mask with states
        allocate(ptcl_mask_not(params_glob%fromp:params_glob%top))
        ptcl_mask_not = .not.ptcl_mask
        do iptcl=params_glob%fromp,params_glob%top
            if( ptcl_mask_not(iptcl) )ptcl_mask_not(iptcl) = build_glob%spproj_field%get_state(iptcl)>0
        enddo
        ! create a local pftcc object with identical references
        nrefs = pftcc%get_nrefs()
        call pftcc_here%new(nrefs, [1,MAXIMGBATCHSZ])
        call pftcc_here%cp_refs(pftcc)
        ! this is the image batch for reading / polarization
        call prepimgbatch(MAXIMGBATCHSZ)
        ! parallelise over batches of size MAXIMGBATCHSZ
        do iptcl_batch=params_glob%fromp,params_glob%top,MAXIMGBATCHSZ
            ! read images
            batchlims = [iptcl_batch,min(params_glob%top,iptcl_batch + MAXIMGBATCHSZ - 1)]
            nptcls    = batchlims(2)-batchlims(1)+1
            call read_imgbatch(batchlims, ptcl_mask_not)
            ! get eo-arrays
            eoarr = -1
            eoarr(1:nptcls) = nint(build_glob%spproj_field%get_all('eo', batchlims))
            call pftcc_here%set_eos(eoarr)
            ! memoize CTF matrices
            if( trim(params_glob%oritype) .eq. 'ptcl3D' )then
                if( build_glob%spproj%get_ctfflag('ptcl3D').ne.'no' )&
                    &call pftcc_here%create_polar_absctfmats(build_glob%spproj, 'ptcl3D', [1,nptcls])
            endif
            !$omp parallel do default(shared) private(iptcl,imatch) schedule(static) proc_bind(close)
            do iptcl=batchlims(1),batchlims(2)
                if( .not.ptcl_mask_not(iptcl) )then
                    ! stuff already calculated or particle excluded
                else
                    imatch = iptcl - batchlims(1) + 1
                    call prepimg4align(iptcl, build_glob%imgbatch(imatch), match_imgs(imatch))
                    ! transfer to polar coordinates
                    call match_imgs(imatch)%polarize(pftcc_here, imatch, .true., .true., mask=build_glob%l_resmsk)
                    ! calc stats
                    call set_ptcl_stats(pftcc_here, iptcl, imatch)
                endif
            enddo
            !$omp end parallel do
        end do
        ! cleanup
        call pftcc_here%kill
        deallocate(ptcl_mask_not)
    end subroutine calc_ptcl_stats

    !> volumetric 3d reconstruction
    subroutine calc_3Drec( cline, which_iter )
        use simple_fplane, only: fplane
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: which_iter
        type(fplane),    allocatable :: fpls(:)
        type(ctfparams), allocatable :: ctfparms(:)
        type(ori)                :: orientation
        type(kbinterpol)         :: kbwin
        real    :: sdev_noise
        integer :: batchlims(2), iptcl, i, i_batch, ibatch
        select case(trim(params_glob%refine))
            case('eval')
                ! no reconstruction
            case DEFAULT
                c1_symop = sym('c1')
                ! make the gridding prepper
                kbwin = build_glob%eorecvols(1)%get_kbwin()
                ! init volumes
                call preprecvols
                ! prep batch imgs
                call prepimgbatch(MAXIMGBATCHSZ)
                ! allocate array
                allocate(fpls(MAXIMGBATCHSZ),ctfparms(MAXIMGBATCHSZ))
                ! prep batch imgs
                call prepimgbatch(MAXIMGBATCHSZ)
                ! gridding batch loop
                do i_batch=1,nptcls2update,MAXIMGBATCHSZ
                    batchlims = [i_batch,min(nptcls2update,i_batch + MAXIMGBATCHSZ - 1)]
                    call read_imgbatch( nptcls2update, pinds, batchlims)
                    !$omp parallel do default(shared) private(i,iptcl,ibatch) schedule(static) proc_bind(close)
                    do i=batchlims(1),batchlims(2)
                        iptcl  = pinds(i)
                        ibatch = i - batchlims(1) + 1
                        if( .not.fpls(ibatch)%does_exist() ) call fpls(ibatch)%new(build_glob%imgbatch(1), build_glob%spproj)
                        call build_glob%imgbatch(ibatch)%noise_norm(build_glob%lmsk, sdev_noise)
                        call build_glob%imgbatch(ibatch)%fft
                        ctfparms(ibatch) = build_glob%spproj%get_ctfparams(params_glob%oritype, iptcl)
                        call fpls(ibatch)%gen_planes(build_glob%imgbatch(ibatch), ctfparms(ibatch))
                    end do
                    !$omp end parallel do
                    ! gridding
                    do i=batchlims(1),batchlims(2)
                        iptcl       = pinds(i)
                        ibatch      = i - batchlims(1) + 1
                        call build_glob%spproj_field%get_ori(iptcl, orientation)
                        if( orientation%isstatezero() ) cycle
                        call eucl_sigma%set_sigma2(iptcl)
                        select case(trim(params_glob%refine))
                            case('clustersym')
                                ! always C1 reconstruction
                                call grid_ptcl(fpls(ibatch), c1_symop, orientation, s3D%o_peaks(iptcl))
                            case DEFAULT
                                call grid_ptcl(fpls(ibatch), build_glob%pgrpsyms, orientation, s3D%o_peaks(iptcl))
                        end select
                    end do
                end do
                ! normalise structure factors
                call norm_struct_facts( cline, which_iter)
                ! destruct
                call killrecvols()
                do ibatch=1,MAXIMGBATCHSZ
                    call fpls(ibatch)%kill
                end do
                deallocate(fpls,ctfparms)
       end select
       call orientation%kill
    end subroutine calc_3Drec

    subroutine setup_weights_read_o_peaks
        ! set npeaks
        npeaks = NPEAKS2REFINE
        ! particle weights
        if( params_glob%l_ptclw )then
            call build_glob%spproj_field%calc_soft_weights(params_glob%frac)
        else
            call build_glob%spproj_field%calc_hard_weights(params_glob%frac)
        endif
        ! prepare particle mask
        allocate(ptcl_mask(params_glob%fromp:params_glob%top))
        call build_glob%spproj_field%sample4update_and_incrcnt_nofrac([params_glob%fromp,params_glob%top],&
        nptcls2update, pinds, ptcl_mask)
        ! allocate s3D singleton
        call prep_strategy3D(ptcl_mask, npeaks)
        ! read peaks
        call read_o_peaks
    end subroutine setup_weights_read_o_peaks

    subroutine read_o_peaks
        use simple_strategy3D_utils, only: update_softmax_weights
        integer :: iptcl, n_nozero
        if( .not. file_exists(trim(params_glob%o_peaks_file)) )then
            THROW_HARD(trim(params_glob%o_peaks_file)//' file does not exist')
        endif
        call open_o_peaks_io(trim(params_glob%o_peaks_file))
        do iptcl=params_glob%fromp,params_glob%top
            if( ptcl_mask(iptcl) )then
                call read_o_peak(s3D%o_peaks(iptcl), [params_glob%fromp,params_glob%top], iptcl, n_nozero)
                call update_softmax_weights(iptcl, npeaks)
            endif
        end do
        call close_o_peaks_io
    end subroutine read_o_peaks

    subroutine calc_global_ori_weights
        real,    allocatable :: weights_glob(:), weights(:), ws_nonzero(:)
        real    :: weight_thres, wsum
        integer :: nweights, cnt, iptcl, i, nw
        select case(params_glob%refine)
        case('cluster', 'snhc', 'clustersym', 'cont_single', 'eval')
                ! nothing to do
            case DEFAULT
                ! extract weights
                nweights = npeaks * nptcls2update
                allocate(weights_glob(nweights), source=0.)
                cnt = 0
                do iptcl=params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) )then
                        weights = s3D%o_peaks(iptcl)%get_all('ow')
                        do i=1,size(weights)
                            cnt = cnt + 1
                            weights_glob(cnt) = weights(i)
                        end do
                    endif
                end do
                ! define a threshold using Otsu's algorithm
                ws_nonzero = pack(weights_glob(:cnt), mask=weights_glob(:cnt) > TINY)
                call otsu(ws_nonzero, weight_thres)
                deallocate(ws_nonzero)
                if( weight_thres <= TINY ) weight_thres = SMALL
                ! zero and normalize weights, apply rank-based weighting scheme if so specified
                do iptcl=params_glob%fromp,params_glob%top
                    if( ptcl_mask(iptcl) )then
                        weights = s3D%o_peaks(iptcl)%get_all('ow')
                        nw = size(weights)
                        where( weights <= weight_thres ) weights = 0.
                        wsum = sum(weights)
                        if( wsum > TINY )then
                            weights = weights / wsum
                            select case(params_glob%wcrit_enum)
                                case(RANK_SUM_CRIT,RANK_CEN_CRIT,RANK_EXP_CRIT,RANK_INV_CRIT)
                                    if( params_glob%wcrit_enum == RANK_EXP_CRIT )then
                                        call conv2rank_weights(nw, weights, params_glob%wcrit_enum, RANKW_EXP)
                                    else
                                        call conv2rank_weights(nw, weights, params_glob%wcrit_enum)
                                    endif
                                    call s3D%o_peaks(iptcl)%set_all('ow', weights)
                                    call build_glob%spproj_field%set(iptcl, 'npeaks', real(count(weights > TINY)))
                                    call build_glob%spproj_field%set(iptcl, 'ow', maxval(weights))
                                    deallocate(weights)
                                case DEFAULT
                                    call s3D%o_peaks(iptcl)%set_all('ow', weights)
                                    call build_glob%spproj_field%set(iptcl, 'npeaks', real(count(weights > TINY)))
                                    call build_glob%spproj_field%set(iptcl, 'ow',     maxval(weights))
                                    deallocate(weights)
                            end select
                        else
                            call s3D%o_peaks(iptcl)%set_all2single('ow', 0.)
                            call build_glob%spproj_field%set(iptcl, 'npeaks', 0.)
                        endif
                    endif
                end do
                deallocate(weights_glob)
        end select
    end subroutine calc_global_ori_weights

    subroutine calc_proj_weights
        real, allocatable :: weights(:), projs(:)
        integer :: i, ind, iptcl
        real    :: pw, proj_weights(params_glob%nspace), minw
        select case(params_glob%refine)
            case('cluster', 'snhc', 'clustersym', 'cont_single', 'eval')
                ! nothing to do
            case DEFAULT
                if( build_glob%spproj_field%get_avg('updatecnt') < 1.0 )then
                    ! nothing to do
                else
                    ! calculate the weight strenght per projection direction
                    proj_weights = 0.
                    do iptcl=params_glob%fromp,params_glob%top
                        pw = 1.0
                        if( build_glob%spproj_field%isthere(iptcl, 'w') ) pw = build_glob%spproj_field%get(iptcl, 'w')
                        if( s3D%o_peaks(iptcl)%isthere('ow') )then
                            weights = s3D%o_peaks(iptcl)%get_all('ow')
                            projs   = s3D%o_peaks(iptcl)%get_all('proj')
                            do i=1,size(projs)
                                ind = nint(projs(i))
                                if( weights(i) > TINY ) proj_weights(ind) = proj_weights(ind) + weights(i) * pw
                            end do
                            deallocate(weights, projs)
                        endif
                    end do
                    minw = minval(proj_weights, mask=proj_weights > TINY)
                    where( proj_weights < TINY ) proj_weights = minw ! to prevent division with zero
                endif
                call arr2file(proj_weights, params_glob%proj_weights_file)
        end select
    end subroutine calc_proj_weights

end module simple_strategy3D_matcher
