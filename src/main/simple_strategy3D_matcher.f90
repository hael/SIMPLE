! projection-matching based on Hadamard products, high-level search routines for REFINE3D
module simple_strategy3D_matcher
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_strategy3D_alloc  ! singleton s3D
use simple_singletons
use simple_timer
implicit none

public :: refine3D_exec
public :: preppftcc4align, pftcc
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
        use simple_qsys_funs,  only: qsys_job_finished
        use simple_binoris_io, only: binwrite_oritab
        use simple_kbinterpol, only: kbinterpol
        use simple_ori,        only: ori
        use simple_sym,        only: sym
        use simple_prep4cgrid, only: prep4cgrid
        use simple_image,      only: image
        use simple_cmdline,    only: cmdline
        use simple_strategy2D3D_common,   only: killrecvols, set_bp_range, preprecvols,&
            prepimgbatch, grid_ptcl, read_imgbatch, eonorm_struct_facts,norm_struct_facts
        use simple_strategy3D_cluster,       only: strategy3D_cluster
        use simple_strategy3D_single,        only: strategy3D_single
        use simple_strategy3D_multi,         only: strategy3D_multi
        use simple_strategy3D_snhc_single,   only: strategy3D_snhc_single
        use simple_strategy3D_greedy_single, only: strategy3D_greedy_single
        use simple_strategy3D_greedy_multi,  only: strategy3D_greedy_multi
        use simple_strategy3D_cont_single,   only: strategy3D_cont_single
        use simple_strategy3D,               only: strategy3D
        use simple_strategy3D_srch,          only: strategy3D_spec
        use simple_convergence,              only: convergence
        class(cmdline),        intent(inout) :: cline
        integer,               intent(in)    :: which_iter
        logical,               intent(inout) :: converged
        type(image),     allocatable :: rec_imgs(:)
        integer, target, allocatable :: symmat(:,:)
        logical,         allocatable :: het_mask(:)
        !---> The below is to allow particle-dependent decision about which 3D strategy to use
        type :: strategy3D_per_ptcl
            class(strategy3D), pointer :: ptr => null()
        end type strategy3D_per_ptcl
        type(strategy3D_per_ptcl), allocatable :: strategy3Dsrch(:)
        !<---- hybrid or combined search strategies can then be implemented as extensions of the
        !      relevant strategy3D base class
        type(strategy3D_spec) :: strategy3Dspec
        type(ori)             :: orientation
        type(kbinterpol)      :: kbwin
        type(sym)             :: c1_symop
        type(prep4cgrid)      :: gridprep
        type(ctfparams)       :: ctfvars
        type(convergence)     :: conv
        real    :: frac_srch_space, extr_thresh, corr_thresh, bfac_rec, specscore_avg
        integer :: iptcl, iextr_lim, i, zero_pop, fnr, cnt, i_batch
        integer :: ibatch, npeaks, batchlims(2), updatecnt
        logical :: doprint, do_extr

        if( L_BENCH )then
            t_init = tic()
            t_tot  = t_init
        endif

        if( p%dryrun.eq.'yes' )then
            ! TBD
        endif

        ! CHECK THAT WE HAVE AN EVEN/ODD PARTITIONING
        if( p%eo .ne. 'no' )then
            if( b%a%get_nevenodd() == 0 ) stop 'ERROR! no eo partitioning available; strategy3D_matcher :: refine3D_exec'
        else
            call b%a%set_all2single('eo', -1.)
        endif

        ! SET FOURIER INDEX RANGE
                call set_bp_range(cline)

        ! DETERMINE THE NUMBER OF PEAKS
        select case(p%refine)
            case('cluster', 'snhc', 'clustersym', 'clusterdev', 'cont_single')
                npeaks = 1
            case DEFAULT
                if( p%eo .ne. 'no' )then
                    npeaks = min(b%e%find_npeaks_from_athres(NPEAKSATHRES), MAXNPEAKS)
                else
                    npeaks = min(10,b%e%find_npeaks(p%lp, p%moldiam))
                endif
        end select
        if( DEBUG ) print *, '*** strategy3D_matcher ***: determined the number of peaks'

        ! SET FRACTION OF SEARCH SPACE
        frac_srch_space = b%a%get_avg('frac')

        ! SETUP WEIGHTS
        if( p%weights3D.eq.'yes' )then
            if( p%nptcls <= SPECWMINPOP )then
                call b%a%calc_hard_weights(p%frac)
            else
                call b%a%calc_spectral_weights(p%frac)
            endif
        else
            call b%a%calc_hard_weights(p%frac)
        endif

        ! READ FOURIER RING CORRELATIONS
        if( file_exists(p%frcs) ) call b%projfrcs%read(p%frcs)

        ! PARTICLE INDEX SAMPLING FOR FRACTIONAL UPDATE (OR NOT)
        if( allocated(pinds) )     deallocate(pinds)
        if( allocated(ptcl_mask) ) deallocate(ptcl_mask)
        if( p%l_frac_update )then
            allocate(ptcl_mask(p%fromp:p%top))
            call b%a%sample4update_and_incrcnt([p%fromp,p%top], p%update_frac, nptcls2update, pinds, ptcl_mask)
            ! correct convergence stats
            do iptcl=p%fromp,p%top
                if( .not. ptcl_mask(iptcl) )then
                    ! these are not updated
                    call b%a%set(iptcl, 'mi_proj',     1.0)
                    call b%a%set(iptcl, 'mi_inpl',     1.0)
                    call b%a%set(iptcl, 'mi_state',    1.0)
                    call b%a%set(iptcl, 'mi_joint',    1.0)
                    call b%a%set(iptcl, 'dist',        0.0)
                    call b%a%set(iptcl, 'dist_inpl',   0.0)
                    call b%a%set(iptcl, 'frac',      100.0)
                endif
            end do
        else
            nptcls2update = p%top - p%fromp + 1
            allocate(pinds(nptcls2update), ptcl_mask(p%fromp:p%top))
            pinds = (/(i,i=p%fromp,p%top)/)
            ptcl_mask = .true.
        endif

        ! B-factor weighted reconstruction
        if( p%shellw.eq.'yes' ) call b%a%calc_bfac_rec

        ! EXTREMAL LOGICS
        do_extr = .false.
        select case(trim(p%refine))
            case('cluster','clusterdev','clustersym')
                if(allocated(het_mask))deallocate(het_mask)
                allocate(het_mask(p%fromp:p%top), source=ptcl_mask)
                zero_pop    = count(.not.b%a%included(consider_w=.false.))
                corr_thresh = -huge(corr_thresh)
                if( p%l_frac_update )then
                    ptcl_mask = .true.
                    iextr_lim = ceiling(2.*log(real(p%nptcls-zero_pop)) * (2.-p%update_frac))
                    if( which_iter==1 .or.(frac_srch_space <= 99. .and. p%extr_iter <= iextr_lim) )&
                        &do_extr = .true.
                else
                    iextr_lim = ceiling(2.*log(real(p%nptcls-zero_pop)))
                    if( which_iter==1 .or.(frac_srch_space <= 98. .and. p%extr_iter <= iextr_lim) )&
                        &do_extr = .true.
                endif
                if( do_extr )then
                    extr_thresh = EXTRINITHRESH * cos(PI/2. * real(p%extr_iter-1)/real(iextr_lim))
                    corr_thresh = b%a%extremal_bound(extr_thresh)
                endif
                if(trim(p%refine).eq.'clustersym')then
                   ! symmetry pairing matrix
                    c1_symop = sym('c1')
                    p%nspace = min( p%nspace*b%se%get_nsym(), 3000 )
                    call b%e%new( p%nspace )
                    call b%e%spiral
                    call b%se%nearest_sym_neighbors( b%e, symmat )
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
        ! cannot kill b%vol since used in continuous search
        call b%vol2%kill

        ! array allocation for strategy3D
        if( DEBUG ) print *, '*** strategy3D_matcher ***: array allocation for strategy3D'
        call prep_strategy3D( ptcl_mask, npeaks )  ! allocate s3D singleton
        if( DEBUG ) print *, '*** strategy3D_matcher ***: array allocation for strategy3D, DONE'
        if( L_BENCH ) rt_prep_primesrch3D = toc(t_prep_primesrch3D)
        ! switch for per-particle polymorphic strategy3D construction
        allocate(strategy3Dsrch(p%fromp:p%top), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("In simple_strategy3D_matcher::refine3D_exec strategy3Dsrch",alloc_stat)
        select case(trim(p%refine))
            case('snhc')
                do iptcl=p%fromp,p%top
                    if( ptcl_mask(iptcl) ) then
                        allocate(strategy3D_snhc_single :: strategy3Dsrch(iptcl)%ptr, stat=alloc_stat)
                        if(alloc_stat.ne.0)call allocchk("In simple_strategy3D_matcher::refine3D_exec strategy3Dsrch snhc",alloc_stat)
                    endif
                end do
            case('single')
                ! polymorphic assignment based on is_virgin and updatecnt
                do iptcl=p%fromp,p%top
                    if( ptcl_mask(iptcl) )then
                        updatecnt = nint(b%a%get(iptcl,'updatecnt'))
                        if( .not.b%a%has_been_searched(iptcl) .or. updatecnt == 1 )then
                            allocate(strategy3D_greedy_single :: strategy3Dsrch(iptcl)%ptr, stat=alloc_stat)
                        else
                            allocate(strategy3D_single        :: strategy3Dsrch(iptcl)%ptr, stat=alloc_stat)
                        endif
                        if(alloc_stat.ne.0)call allocchk("In simple_strategy3D_matcher::refine3D_exec strategy3Dsrch single",alloc_stat)
                    endif
                end do
            case('multi')
                do iptcl=p%fromp,p%top
                    if( ptcl_mask(iptcl) )then
                        updatecnt = nint(b%a%get(iptcl,'updatecnt'))
                        if( updatecnt == 1 )then
                            allocate(strategy3D_greedy_multi :: strategy3Dsrch(iptcl)%ptr, stat=alloc_stat)
                        else
                            allocate(strategy3D_multi        :: strategy3Dsrch(iptcl)%ptr, stat=alloc_stat)
                        endif
                        if(alloc_stat.ne.0)call allocchk("In simple_strategy3D_matcher::refine3D_exec strategy3Dsrch multi",alloc_stat)
                    endif
                end do
            case('cont_single')
                do iptcl=p%fromp,p%top
                    if( ptcl_mask(iptcl) )then
                        allocate(strategy3D_cont_single   :: strategy3Dsrch(iptcl)%ptr, stat=alloc_stat)
                        if(alloc_stat.ne.0)call allocchk("In simple_strategy3D_matcher::refine3D_exec strategy3Dsrch snhc",alloc_stat)
                    endif
                end do
            case('greedy_single')
                do iptcl=p%fromp,p%top
                    if( ptcl_mask(iptcl) ) allocate(strategy3D_greedy_single :: strategy3Dsrch(iptcl)%ptr)
                end do
            case('greedy_multi')
                do iptcl=p%fromp,p%top
                    if( ptcl_mask(iptcl) )then
                        allocate(strategy3D_greedy_multi  :: strategy3Dsrch(iptcl)%ptr, stat=alloc_stat)
                        if(alloc_stat.ne.0)call allocchk("In simple_strategy3D_matcher::refine3D_exec strategy3Dsrch snhc",alloc_stat)
                    endif
                end do
            case('cluster','clustersym','clusterdev')
                do iptcl=p%fromp,p%top
                    if( ptcl_mask(iptcl) )then
                        allocate(strategy3D_cluster       :: strategy3Dsrch(iptcl)%ptr, stat=alloc_stat)
                        if(alloc_stat.ne.0)call allocchk("In simple_strategy3D_matcher::refine3D_exec strategy3Dsrch snhc",alloc_stat)
                    endif
                end do
            case DEFAULT
                write(*,*) 'refine flag: ', trim(p%refine)
                stop 'Refinement mode unsupported'
        end select
        ! construct search objects
        cnt = 0
        do iptcl=p%fromp,p%top
            if( ptcl_mask(iptcl) )then
                cnt = cnt + 1
                ! search spec
                strategy3Dspec%iptcl       =  iptcl
                strategy3Dspec%iptcl_map   =  cnt
                strategy3Dspec%szsn        =  p%szsn
                strategy3Dspec%corr_thresh =  corr_thresh
                strategy3Dspec%ppftcc      => pftcc
                if( allocated(het_mask) )     strategy3Dspec%do_extr    =  het_mask(iptcl)
                if( allocated(symmat) )       strategy3Dspec%symmat     => symmat
                ! search object
                call strategy3Dsrch(iptcl)%ptr%new(strategy3Dspec, npeaks)
            endif
        end do
        if( DEBUG ) print *, '*** strategy3D_matcher ***: search object construction, DONE'
        ! memoize CTF matrices
        if( trim(p%oritype) .eq. 'ptcl3D' )then
            if( b%spproj%get_ctfflag('ptcl3D').ne.'no' ) call pftcc%create_polar_absctfmats(b%spproj, 'ptcl3D')
        else
            ! class averages have no CTF
        endif
        ! memoize FFTs
        call pftcc%memoize_ffts

        ! SEARCH
        if( p%dryrun.eq.'yes' )then
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
        call b%vol%kill
        call b%vol_odd%kill
        do iptcl = p%fromp,p%top
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
        select case(trim(p%oritype))
            case('ptcl3D')
                call binwrite_oritab(p%outfile, b%spproj, b%a, [p%fromp,p%top], isegment=PTCL3D_SEG)
            case('cls3D')
                call binwrite_oritab(p%outfile, b%spproj, b%a, [p%fromp,p%top], isegment=CLS3D_SEG)
            case DEFAULT
                write(*,*) 'oritype: ', trim(p%oritype)
                stop 'Unsupported oritype; strategy3D_matcher :: refine3D_exec'
        end select
        p%oritab = p%outfile

        ! VOLUMETRIC 3D RECONSTRUCTION
        if( L_BENCH ) t_rec = tic()
        ! make the gridding prepper
        if( p%eo .ne. 'no' )then
            kbwin = b%eorecvols(1)%get_kbwin()
        else
            kbwin = b%recvols(1)%get_kbwin()
        endif
        call gridprep%new(b%img, kbwin, [p%boxpd,p%boxpd,1])
        ! init volumes
        call preprecvols
        ! prep rec imgs
        allocate(rec_imgs(MAXIMGBATCHSZ))
        do i=1,MAXIMGBATCHSZ
            call rec_imgs(i)%new([p%boxpd, p%boxpd, 1], p%smpd)
        end do
        ! prep batch imgs
        call prepimgbatch( MAXIMGBATCHSZ)
        !call prepimgbatch(b, p, MAXIMGBATCHSZ)
        ! gridding batch loop
        do i_batch=1,nptcls2update,MAXIMGBATCHSZ
            batchlims = [i_batch,min(nptcls2update,i_batch + MAXIMGBATCHSZ - 1)]
            call read_imgbatch( nptcls2update, pinds, batchlims)
            ! parallel gridprep
            !$omp parallel do default(shared) private(i,ibatch) schedule(static) proc_bind(close)
            do i=batchlims(1),batchlims(2)
                ibatch = i - batchlims(1) + 1
                ! normalise (read_imgbatch does not normalise)
                call b%imgbatch(ibatch)%norm()
                ! in dev=yes code, we filter before inserting into 3D vol
                call gridprep%prep_serial_no_fft(b%imgbatch(ibatch), rec_imgs(ibatch))
            end do
            !$omp end parallel do
            ! gridding
            do i=batchlims(1),batchlims(2)
                iptcl       = pinds(i)
                ibatch      = i - batchlims(1) + 1
                orientation = b%a%get_ori(iptcl)
                ctfvars     = b%spproj%get_ctfparams(p%oritype, iptcl)
                if( orientation%isstatezero() ) cycle
                if( trim(p%refine).eq.'clustersym' )then
                    ! always C1 reconstruction
                    call grid_ptcl( rec_imgs(ibatch), c1_symop, orientation, s3D%o_peaks(iptcl), ctfvars)
                else
                    call grid_ptcl( rec_imgs(ibatch), b%se, orientation, s3D%o_peaks(iptcl), ctfvars)
                endif
            end do
        end do
        ! normalise structure factors
        if( p%eo .ne. 'no' )then
            call eonorm_struct_facts( cline, which_iter)
        else
            call norm_struct_facts( which_iter)
        endif
        ! destruct
        call killrecvols()
        call gridprep%kill
        do ibatch=1,MAXIMGBATCHSZ
            call rec_imgs(ibatch)%kill
            call b%imgbatch(ibatch)%kill
        end do
        deallocate(rec_imgs, b%imgbatch)
        call gridprep%kill
        if( L_BENCH ) rt_rec = toc(t_rec)

        ! REPORT CONVERGENCE
        call qsys_job_finished( 'simple_strategy3D_matcher :: refine3D_exec')
        if( .not. p%l_distr_exec ) converged = conv%check_conv3D(cline)
        if( L_BENCH )then
            rt_tot  = toc(t_tot)
            doprint = .true.
            if( p%part /= 1 ) doprint = .false.
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
        use simple_polarizer, only: polarizer
        use simple_cmdline,   only: cmdline
        use simple_strategy2D3D_common,   only: cenrefvol_and_mapshifts2ptcls, preprefvol, build_pftcc_particles
        class(cmdline),             intent(inout) :: cline !< command line
        type(polarizer), allocatable :: match_imgs(:)
        integer   :: cnt, s, ind, iref, nrefs
        integer   :: imatch
        logical   :: do_center
        real      :: xyz(3)
        nrefs  = p%nspace * p%nstates
        ! must be done here since p%kfromto is dynamically set based on FSC from previous round
        ! or based on dynamic resolution limit update
        if( p%eo .ne. 'no' )then
            call pftcc%new(nrefs,  ptcl_mask, nint(b%a%get_all('eo', [p%fromp,p%top])))
        else
            call pftcc%new(nrefs,  ptcl_mask)
        endif

        ! PREPARATION OF REFERENCES IN PFTCC
        ! read reference volumes and create polar projections
        cnt = 0
        do s=1,p%nstates
            if( p%oritab .ne. '' )then
                if( b%a%get_pop(s, 'state') == 0 )then
                    ! empty state
                    cnt = cnt + p%nspace
                    call progress(cnt, nrefs)
                    cycle
                endif
            endif
            call cenrefvol_and_mapshifts2ptcls( cline, s, p%vols(s), do_center, xyz)
            if( p%eo .ne. 'no' )then
                if( p%nstates.eq.1 )then
                    ! PREPARE ODD REFERENCES
                    call preprefvol(cline, s, p%vols_odd(s), do_center, xyz)
                    !$omp parallel do default(shared) private(iref) schedule(static) proc_bind(close)
                    do iref=1,p%nspace
                        call b%vol%fproject_polar((s - 1) * p%nspace + iref, b%e%get_ori(iref), pftcc, iseven=.false.)
                    end do
                    !$omp end parallel do
                    ! copy odd volume
                    b%vol_odd = b%vol
                    ! expand for fast interpolation
                    call b%vol_odd%expand_cmat(p%alpha)
                    ! PREPARE EVEN REFERENCES
                    call preprefvol( cline, s, p%vols_even(s), do_center, xyz)
                    !$omp parallel do default(shared) private(iref) schedule(static) proc_bind(close)
                    do iref=1,p%nspace
                        call b%vol%fproject_polar((s - 1) * p%nspace + iref, b%e%get_ori(iref), pftcc, iseven=.true.)
                    end do
                    !$omp end parallel do
                else
                    call preprefvol( cline, s, p%vols(s), do_center, xyz)
                    !$omp parallel do default(shared) private(iref, ind) schedule(static) proc_bind(close)
                    do iref=1,p%nspace
                        ind = (s - 1) * p%nspace + iref
                        call b%vol%fproject_polar(ind, b%e%get_ori(iref), pftcc, iseven=.true.)
                        call pftcc%cp_even2odd_ref(ind)
                    end do
                    !$omp end parallel do
                endif
            else
                ! low-pass set or multiple states
                call preprefvol( cline, s, p%vols(s), do_center, xyz)
                !$omp parallel do default(shared) private(iref) schedule(static) proc_bind(close)
                do iref=1,p%nspace
                    call b%vol%fproject_polar((s - 1) * p%nspace + iref, b%e%get_ori(iref), pftcc, iseven=.true.)
                end do
                !$omp end parallel do
            endif
        end do

        ! PREPARATION OF PARTICLES IN PFTCC
        ! prepare the polarizer images
        call b%img_match%init_polarizer(pftcc, p%alpha)
        allocate(match_imgs(MAXIMGBATCHSZ))
        do imatch=1,MAXIMGBATCHSZ
            call match_imgs(imatch)%new([p%boxmatch, p%boxmatch, 1], p%smpd)
            call match_imgs(imatch)%copy_polarizer(b%img_match)
        end do
        call build_pftcc_particles(pftcc, MAXIMGBATCHSZ, match_imgs, .true., ptcl_mask)

        ! DESTRUCT
        do imatch=1,MAXIMGBATCHSZ
            call match_imgs(imatch)%kill_polarizer
            call match_imgs(imatch)%kill
            call b%imgbatch(imatch)%kill
        end do
        deallocate(match_imgs, b%imgbatch)
        if( DEBUG ) print *, '*** strategy3D_matcher ***: finished preppftcc4align'
    end subroutine preppftcc4align

end module simple_strategy3D_matcher
