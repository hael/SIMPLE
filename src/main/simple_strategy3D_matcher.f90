! projection-matching based on Hadamard products, high-level search routines for PRIME3D
module simple_strategy3D_matcher
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'

use simple_ori,                      only: ori
use simple_oris,                     only: oris
use simple_build,                    only: build
use simple_params,                   only: params
use simple_cmdline,                  only: cmdline
use simple_binoris_io,               only: binwrite_oritab
use simple_qsys_funs,                only: qsys_job_finished
use simple_kbinterpol,               only: kbinterpol
use simple_prep4cgrid,               only: prep4cgrid
use simple_polarft_corrcalc,         only: polarft_corrcalc
use simple_strategy3D,               only: strategy3D
use simple_strategy3D_srch,          only: strategy3D_spec
use simple_strategy3D_alloc,         only: o_peaks, clean_strategy3D, prep_strategy3D
use simple_strategy3D_cluster,       only: strategy3D_cluster
use simple_strategy3D_single,        only: strategy3D_single
use simple_strategy3D_multi,         only: strategy3D_multi
use simple_strategy3D_snhc_single,   only: strategy3D_snhc_single
use simple_strategy3D_greedy_single, only: strategy3D_greedy_single
use simple_strategy3D_greedy_multi,  only: strategy3D_greedy_multi
use simple_strategy3D_cont_single,   only: strategy3D_cont_single
use simple_strategy2D3D_common       ! use all in there
implicit none

public :: refine3D_exec
public :: preppftcc4align, pftcc
private
! #include "simple_local_flags.inc"

logical, parameter             :: L_BENCH = .false., DEBUG = .true.
type(polarft_corrcalc), target :: pftcc
integer, allocatable           :: pinds(:)
logical, allocatable           :: ptcl_mask(:)
integer                        :: nptcls2update
integer(timer_int_kind)        :: t_init, t_prep_pftcc, t_align, t_rec, t_tot, t_prep_primesrch3D
real(timer_int_kind)           :: rt_init, rt_prep_pftcc, rt_align, rt_rec, rt_prep_primesrch3D
real(timer_int_kind)           :: rt_tot
character(len=STDLEN)          :: benchfname

contains

    subroutine refine3D_exec( b, p, cline, which_iter, converged )
        use simple_qsys_funs, only: qsys_job_finished
        use simple_sym,       only: sym
        use simple_image,     only: image
        class(build),  target, intent(inout) :: b
        class(params), target, intent(inout) :: p
        class(cmdline),        intent(inout) :: cline
        integer,               intent(in)    :: which_iter
        logical,               intent(inout) :: converged
        type(image),     allocatable :: rec_imgs(:)
        integer, target, allocatable :: symmat(:,:)
        logical,         allocatable :: het_mask(:)
        class(strategy3D), pointer   :: strategy3Dsrch(:)
        type(strategy3D_spec) :: strategy3Dspec
        type(ori)             :: orientation
        type(kbinterpol)      :: kbwin
        type(sym)             :: c1_symop
        type(prep4cgrid)      :: gridprep
        type(ctfparams)       :: ctfvars
        real    :: frac_srch_space, reslim, extr_thresh, corr_thresh
        real    :: bfac_rec, specscore_avg
        integer :: iptcl, iextr_lim, i, zero_pop, fnr, cnt, i_batch, ibatch, npeaks
        integer :: batchlims(2)
        logical :: doprint, do_extr, is_virgin

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
        call set_bp_range( b, p, cline )

        ! CALCULATE ANGULAR THRESHOLD (USED BY THE SPARSE WEIGHTING SCHEME)
        reslim   = p%lp
        if( DEBUG ) print *, '*** strategy3D_matcher ***: calculated angular threshold (used by the sparse weighting scheme)'

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

        ! SPECSCORE & B-FACTOR RECONSTRUCTION
        if( b%a%isthere('specscore') )then
            specscore_avg = b%a%get_avg('specscore')
            do iptcl = p%fromp,p%top
                if( .not.ptcl_mask(iptcl) ) cycle
                bfac_rec  = BSC * (b%a%get(iptcl,'specscore') - specscore_avg)
                call b%a%set(iptcl, 'bfac_rec', bfac_rec)
            enddo
        else
            specscore_avg = 0.
            call b%a%set_all2single('bfac_rec', 0.)
        endif

        ! EXTREMAL LOGICS
        do_extr  = .false.
        select case(trim(p%refine))
            case('cluster','clusterdev','clustersym')
                if(allocated(het_mask))deallocate(het_mask)
                allocate(het_mask(p%fromp:p%top), source=ptcl_mask)
                zero_pop    = count(.not.b%a%included(consider_w=.false.))
                corr_thresh = -huge(corr_thresh)
                if(p%l_frac_update) then
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
        call preppftcc4align( b, p, cline )
        if( L_BENCH ) rt_prep_pftcc = toc(t_prep_pftcc)

        write(*,'(A,1X,I3)') '>>> REFINE3D SEARCH, ITERATION:', which_iter

        ! STOCHASTIC IMAGE ALIGNMENT
        if( L_BENCH ) t_prep_primesrch3D = tic()
        ! clean big objects before starting to allocate new big memory chunks
        ! cannot kill b%vol since used in continuous search
        call b%vol2%kill

        ! array allocation for strategy3D
        if( DEBUG ) print *, '*** strategy3D_matcher ***: array allocation for strategy3D'
        call prep_strategy3D( b, p, ptcl_mask, npeaks )
        if( DEBUG ) print *, '*** strategy3D_matcher ***: array allocation for strategy3D, DONE'
        if( L_BENCH ) rt_prep_primesrch3D = toc(t_prep_primesrch3D)
        ! switch for polymorphic strategy3D construction
        select case(trim(p%refine))
            case('snhc')
                allocate(strategy3D_snhc_single :: strategy3Dsrch(p%fromp:p%top), stat=alloc_stat)
            case('single')
                ! check if virgin
                select case(trim(p%oritype))
                    case('ptcl3D')
                        is_virgin = b%spproj%is_virgin_field('ptcl3D')
                    case('cls3D')
                        is_virgin = b%spproj%is_virgin_field('cls3D')
                    case DEFAULT
                        write(*,*) 'oritype: ', trim(p%oritype)
                        stop 'Unsupported oritype; strategy3D_matcher :: refine3D_exec'
                end select
                ! allocate the appropriate polymorphic type
                if( is_virgin )then
                    allocate(strategy3D_greedy_single :: strategy3Dsrch(p%fromp:p%top))
                else
                    allocate(strategy3D_single        :: strategy3Dsrch(p%fromp:p%top))
                endif
            case('multi')
                allocate(strategy3D_multi         :: strategy3Dsrch(p%fromp:p%top), stat=alloc_stat)
            case('cont_single')
                allocate(strategy3D_cont_single   :: strategy3Dsrch(p%fromp:p%top), stat=alloc_stat)
            case('greedy_single')
                allocate(strategy3D_greedy_single :: strategy3Dsrch(p%fromp:p%top), stat=alloc_stat)
            case('greedy_multi')
                allocate(strategy3D_greedy_multi  :: strategy3Dsrch(p%fromp:p%top), stat=alloc_stat)
            case('cluster','clustersym','clusterdev')
                allocate(strategy3D_cluster       :: strategy3Dsrch(p%fromp:p%top), stat=alloc_stat)
            case DEFAULT
                write(*,*) 'refine flag: ', trim(p%refine)
                stop 'Refinement mode unsupported'
        end select
        if(alloc_stat.ne.0)call allocchk("In simple_strategy3D_matcher::prime3D_exec strategy3D objects ",alloc_stat)
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
                strategy3Dspec%pb          => b
                strategy3Dspec%pp          => p
                strategy3Dspec%ppftcc      => pftcc
                strategy3Dspec%pa          => b%a
                strategy3Dspec%pse         => b%se
                if( allocated(b%nnmat) )      strategy3Dspec%nnmat      => b%nnmat
                if( allocated(b%grid_projs) ) strategy3Dspec%grid_projs => b%grid_projs
                if( allocated(het_mask) )     strategy3Dspec%do_extr    =  het_mask(iptcl)
                if( allocated(symmat) )       strategy3Dspec%symmat     => symmat
                ! search object
                call strategy3Dsrch(iptcl)%new(strategy3Dspec, npeaks)
            endif
        end do
        if( DEBUG ) print *, '*** strategy3D_matcher ***: search object construction, DONE'
        ! memoize CTF matrices
        if( trim(p%oritype) .eq. 'ptcl3D' )then
            if( b%spproj%get_ctfflag('ptcl3D').ne.'no' ) call pftcc%create_polar_ctfmats(b%spproj, 'ptcl3D')
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
                call strategy3Dsrch(iptcl)%srch
            end do
            !$omp end parallel do
        endif
        ! clean
        call clean_strategy3D()
        call pftcc%kill
        call b%vol%kill
        call b%vol_even%kill
        do iptcl = p%fromp,p%top
            if( ptcl_mask(iptcl) ) call strategy3Dsrch(iptcl)%kill
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
        call preprecvols(b, p)
        ! prep rec imgs
        allocate(rec_imgs(MAXIMGBATCHSZ))
        do i=1,MAXIMGBATCHSZ
            call rec_imgs(i)%new([p%boxpd, p%boxpd, 1], p%smpd)
        end do
        ! prep batch imgs
        call prepimgbatch(b, p, MAXIMGBATCHSZ)
        ! gridding batch loop
        do i_batch=1,nptcls2update,MAXIMGBATCHSZ
            batchlims = [i_batch,min(nptcls2update,i_batch + MAXIMGBATCHSZ - 1)]
            call read_imgbatch(b, p, nptcls2update, pinds, batchlims)
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
                    call grid_ptcl(b, p, rec_imgs(ibatch), c1_symop, orientation, o_peaks(iptcl), ctfvars)
                else
                    call grid_ptcl(b, p, rec_imgs(ibatch), b%se, orientation, o_peaks(iptcl), ctfvars)
                endif
            end do
        end do
        ! normalise structure factors
        if( p%eo .ne. 'no' )then
            call eonorm_struct_facts(b, p, cline, reslim, which_iter)
        else
            call norm_struct_facts(b, p, which_iter)
        endif
        ! destruct
        call killrecvols(b, p)
        call gridprep%kill
        do ibatch=1,MAXIMGBATCHSZ
            call rec_imgs(ibatch)%kill
            call b%imgbatch(ibatch)%kill
        end do
        deallocate(rec_imgs, b%imgbatch)
        call gridprep%kill
        if( L_BENCH ) rt_rec = toc(t_rec)

        ! REPORT CONVERGENCE
        call qsys_job_finished( p, 'simple_strategy3D_matcher :: refine3D_exec')
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
    subroutine preppftcc4align( b, p, cline )
        use simple_polarizer, only: polarizer
        class(build),               intent(inout) :: b     !< build object
        class(params),              intent(inout) :: p     !< param object
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
            call pftcc%new(nrefs, p, ptcl_mask, nint(b%a%get_all('eo', [p%fromp,p%top])))
        else
            call pftcc%new(nrefs, p, ptcl_mask)
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
            call cenrefvol_and_mapshifts2ptcls(b, p, cline, s, p%vols(s)%str, do_center, xyz)
            if( p%eo .ne. 'no' )then
                if( p%nstates.eq.1 )then
                    call preprefvol(b, p, cline, s, p%vols_even(s)%str, do_center, xyz)
                    !$omp parallel do default(shared) private(iref) schedule(static) proc_bind(close)
                    do iref=1,p%nspace
                        call b%vol%fproject_polar((s - 1) * p%nspace + iref, b%e%get_ori(iref), pftcc, iseven=.true.)
                    end do
                    !$omp end parallel do
                    ! copy even volume
                    b%vol_even = b%vol
                    ! expand for fast interpolation
                    call b%vol_even%expand_cmat(p%alpha)
                    call preprefvol(b, p, cline, s, p%vols_odd(s)%str, do_center, xyz)
                    !$omp parallel do default(shared) private(iref) schedule(static) proc_bind(close)
                    do iref=1,p%nspace
                        call b%vol%fproject_polar((s - 1) * p%nspace + iref, b%e%get_ori(iref), pftcc, iseven=.false.)
                    end do
                    !$omp end parallel do
                    ! odd volume in b%vol
                else
                    call preprefvol(b, p, cline, s, p%vols(s)%str, do_center, xyz)
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
                call preprefvol(b, p, cline, s, p%vols(s)%str, do_center, xyz)
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
        call build_pftcc_particles( b, p, pftcc, MAXIMGBATCHSZ, match_imgs, .true., ptcl_mask)

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
