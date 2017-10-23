! projection-matching based on Hadamard products, high-level search routines for PRIME3D
module simple_hadamard3D_matcher
!$ use omp_lib
!$ use omp_lib_kinds
#include "simple_lib.f08"
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_prime3D_srch,     only: prime3D_srch
use simple_ori,              only: ori
use simple_build,            only: build
use simple_params,           only: params
use simple_cmdline,          only: cmdline
use simple_gridding,         only: prep4cgrid
use simple_binoris_io,       only: binwrite_oritab
use simple_hadamard_common   ! use all in there
use simple_timer             ! use all in there
implicit none

public :: prime3D_find_resrange, prime3D_exec, gen_random_model
public :: preppftcc4align, prep_refs_pftcc4align, pftcc

private
#include "simple_local_flags.inc"

logical, parameter              :: L_BENCH = .false.
type(polarft_corrcalc)          :: pftcc
type(prime3D_srch), allocatable :: primesrch3D(:)
integer(timer_int_kind)         :: t_init, t_prep_pftcc, t_align, t_rec, t_tot
real(timer_int_kind)            :: rt_init, rt_prep_pftcc, rt_align, rt_rec
real(timer_int_kind)            :: rt_tot
character(len=STDLEN)           :: benchfname

contains

    subroutine prime3D_find_resrange( b, p, lp_start, lp_finish )
        use simple_oris, only: oris
        class(build),  intent(inout) :: b
        class(params), intent(inout) :: p
        real,          intent(out)   :: lp_start, lp_finish
        real, allocatable :: peaks(:)
        type(oris)        :: o
        integer :: lfny, k, pos10, pos6
        call o%new(p%nspace)
        call o%spiral
        lfny = b%img_match%get_lfny(1)
        allocate( peaks(lfny), stat=alloc_stat )
        allocchk("In: prime3D_find_resrange, simple_hadamard3D_matcher")
        do k=2,b%img_match%get_lfny(1)
            peaks(k) = real(o%find_npeaks(b%img_match%get_lp(k), p%moldiam))
        end do
        peaks(1)  = peaks(2)
        pos10     = locate(peaks, lfny, 10.)
        pos6      = locate(peaks, lfny,  6.)
        lp_start  = b%img_match%get_lp(pos10)
        lp_finish = b%img_match%get_lp(pos6)
        deallocate(peaks)
        call o%kill
    end subroutine prime3D_find_resrange

    subroutine prime3D_exec( b, p, cline, which_iter, update_res, converged )
        use simple_qsys_funs, only: qsys_job_finished
        use simple_oris,      only: oris
        use simple_fileio,    only: del_file
        use simple_rnd,       only: irnd_uni
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: which_iter
        logical,        intent(inout) :: update_res, converged
        integer,          allocatable :: proj_space_inds(:,:), proj_srch_order(:,:)
        logical,          allocatable :: state_exists(:), to_update(:)
        type(oris),       allocatable :: reforis(:)
        type(oris) :: prime3D_oris
        type(ori)  :: orientation
        real       :: norm, corr_thresh, skewness, frac_srch_space
        real       :: extr_thresh, update_frac, reslim
        integer    :: iptcl, inptcls, istate, iextr_lim
        integer    :: update_ind, nupdates_target, nupdates, fnr
        integer    :: statecnt(p%nstates)
        logical    :: doprint

        if( L_BENCH )then
            t_init = tic()
            t_tot  = t_init
        endif

        ! CHECK THAT WE HAVE AN EVEN/ODD PARTITIONING
        if( p%eo .ne. 'no' )then
            if( p%l_distr_exec )then
                if( b%a%get_nevenodd() == 0 ) stop 'ERROR! no eo partitioning available; hadamard3D_matcher :: prime2D_exec'
            else
                if( b%a%get_nevenodd() == 0 ) call b%a%partition_eo
            endif
        endif

        inptcls = p%top - p%fromp + 1

        ! SET FOURIER INDEX RANGE
        call set_bp_range( b, p, cline )

        ! CALCULATE ANGULAR THRESHOLD (USED BY THE SPARSE WEIGHTING SCHEME)
        p%athres = rad2deg( atan(max(p%fny,p%lp)/(p%moldiam/2.) ))
        reslim   = p%lp
        DebugPrint '*** hadamard3D_matcher ***: calculated angular threshold (used by the sparse weighting scheme)'

        ! DETERMINE THE NUMBER OF PEAKS
        if( .not. cline%defined('npeaks') )then
            select case(p%refine)
                case('no', 'neigh', 'greedy', 'greedyneigh', 'states', 'tseries')
                    if( p%eo .ne. 'no' )then
                        p%npeaks = min(b%e%find_npeaks_from_athres(NPEAKSATHRES), MAXNPEAKS)
                    else
                        p%npeaks = min(10,b%e%find_npeaks(p%lp, p%moldiam))
                    endif
                case DEFAULT
                    p%npeaks = 1
            end select
            DebugPrint '*** hadamard3D_matcher ***: determined the number of peaks'
        endif

        ! RANDOM MODEL GENERATION
        if( p%vols(1) .eq. '' .and. p%nstates == 1 )then
            if( p%nptcls > 1000 )then
                call gen_random_model(b, p, 1000)
            else
                call gen_random_model(b, p)
            endif
            DebugPrint '*** hadamard3D_matcher ***: generated random model'
        endif

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

        ! POPULATION BALANCING LOGICS
        ! this needs to be done prior to search such that each part
        ! sees the same information in distributed execution
        if( p%balance > 0 )then
            call b%a%balance( p%balance, NSPACE_BALANCE, p%nsym, p%eullims, skewness )
            write(*,'(A,F8.2)') '>>> PROJECTION DISTRIBUTION SKEWNESS(%):', 100. * skewness
        else
            call b%a%set_all2single('state_balance', 1.0)
        endif

        ! EXTREMAL LOGICS
        if( p%refine.eq.'het' )then
            iextr_lim = ceiling(2.*log(real(p%nptcls)))
            if( frac_srch_space < 98. .or. p%extr_iter <= iextr_lim )then
                extr_thresh = EXTRINITHRESH * cos(PI/2. * real(p%extr_iter-1)/real(iextr_lim)) ! cosine decay
                extr_thresh = min(EXTRINITHRESH, max(0., extr_thresh))
                corr_thresh = b%a%extremal_bound(extr_thresh)
                statecnt(:) = 0
            else
                corr_thresh = -huge(corr_thresh)
            endif
        endif

        ! FRACTIONAL UPDATE
        allocate( to_update(p%fromp:p%top), source=.true. )
        if( p%l_frac_update )then
        ! Soon to come
        !     nupdates = inptcls
        !     nupdates_target = nint(p%update_frac * real(inptcls))
        !     do while( nupdates > nupdates_target )
        !         update_ind = irnd_uni(inptcls)
        !         if( to_update(update_ind) )then
        !             to_update(update_ind) = .false.
        !             nupdates = nupdates - 1
        !         else
        !             ! better luck next time
        !         endif
        !     enddo
        else
            ! all done
        endif
        if( L_BENCH ) rt_init = toc(t_init)

        ! PREPARE THE POLARFT_CORRCALC DATA STRUCTURE
        if( L_BENCH ) t_prep_pftcc = tic()
        call preppftcc4align( b, p, cline )
        if( L_BENCH ) rt_prep_pftcc = toc(t_prep_pftcc)

        ! INITIALIZE
        write(*,'(A,1X,I3)') '>>> PRIME3D DISCRETE STOCHASTIC SEARCH, ITERATION:', which_iter
        if( .not. p%l_distr_exec )then
            if( p%refine .eq. 'snhc')then
                p%outfile = SNHCDOC
            else
                p%outfile = 'prime3Ddoc_'//int2str_pad(which_iter,3)//METADATEXT
            endif
        endif

        ! STOCHASTIC IMAGE ALIGNMENT
        ! reference projections indices, here to avoid online allocation in prime3d_srch
        allocate(reforis(p%fromp:p%top))
        do iptcl = p%fromp, p%top
            call reforis(iptcl)%new(p%nstates*p%nspace)
        enddo
        ! reference projections indices, here to avoid online allocation in prime3d_srch
        allocate(proj_space_inds(p%fromp:p%top,1:p%nspace*p%nstates), source=0)
        ! reference projection search order, here to avoid online allocation in prime3d_srch
        allocate(proj_srch_order(p%fromp:p%top,1:p%nspace*p%nstates), source=0)
        ! states existence
        if( p%oritab.ne.'' )then
            state_exists = b%a%states_exist(p%nstates)
        else
            allocate(state_exists(p%nstates), source=.true.)
        endif
        ! create the search objects, need to re-create every round because parameters are changing
        allocate( primesrch3D(p%fromp:p%top) , stat=alloc_stat)
        allocchk("In hadamard3D_matcher::prime3D_exec primesrch3D objects ")
        do iptcl=p%fromp,p%top
            call primesrch3D(iptcl)%new(iptcl, pftcc, b%a, b%e, p, b%se,&
                &proj_space_inds(iptcl,1:p%nspace*p%nstates),&
                &proj_srch_order(iptcl,1:p%nspace*p%nstates),&
                &reforis(iptcl), state_exists)
        end do
        ! apply CTF to particles
        if( p%ctf .ne. 'no' ) call pftcc%apply_ctf_to_ptcls(b%a)
        ! memoize FFTs for improved performance
        call pftcc%memoize_ffts
        ! execute the search
        call del_file(p%outfile)
        if( L_BENCH ) t_align = tic()
        select case(p%refine)
            case( 'snhc' )
                !$omp parallel do default(shared) schedule(guided) private(iptcl) proc_bind(close)
                do iptcl=p%fromp,p%top
                    if( to_update(iptcl) )then
                        call primesrch3D(iptcl)%exec_prime3D_srch(p%lp, szsn=p%szsn)
                    endif
                end do
                !$omp end parallel do
            case( 'no','shc' )
                if( p%oritab .eq. '' )then
                    !$omp parallel do default(shared) schedule(guided) private(iptcl) proc_bind(close)
                    do iptcl=p%fromp,p%top
                        call primesrch3D(iptcl)%exec_prime3D_srch(p%lp, greedy=.true.)
                    end do
                    !$omp end parallel do
                else
                    !$omp parallel do default(shared) schedule(guided) private(iptcl) proc_bind(close)
                    do iptcl=p%fromp,p%top
                        if( to_update(iptcl) )then
                            call primesrch3D(iptcl)%exec_prime3D_srch(p%lp)
                        endif
                    end do
                    !$omp end parallel do
                endif
            case('neigh','shcneigh')
                if( p%oritab .eq. '' ) stop 'cannot run the refine=neigh mode without input oridoc (oritab)'
                !$omp parallel do default(shared) schedule(guided) private(iptcl) proc_bind(close)
                do iptcl=p%fromp,p%top
                    if( to_update(iptcl) )then
                        call primesrch3D(iptcl)%exec_prime3D_srch(p%lp, nnmat=b%nnmat, grid_projs=b%grid_projs)
                    endif
                end do
                !$omp end parallel do
            case('greedy')
                if( p%oritab .eq. '' )then
                    !$omp parallel do default(shared) schedule(guided) private(iptcl) proc_bind(close)
                    do iptcl=p%fromp,p%top
                        call primesrch3D(iptcl)%exec_prime3D_srch(p%lp, greedy=.true.)
                    end do
                    !$omp end parallel do
                else
                    !$omp parallel do default(shared) schedule(guided) private(iptcl) proc_bind(close)
                    do iptcl=p%fromp,p%top
                        if( to_update(iptcl) )then
                            call primesrch3D(iptcl)%exec_prime3D_srch(p%lp, greedy=.true.)
                        endif
                    end do
                    !$omp end parallel do
                endif
            case('greedyneigh')
                if( p%oritab .eq. '' )then                
                    !$omp parallel do default(shared) schedule(guided) private(iptcl) proc_bind(close)
                    do iptcl=p%fromp,p%top
                        call primesrch3D(iptcl)%exec_prime3D_srch(p%lp,&
                            &greedy=.true., nnmat=b%nnmat, grid_projs=b%grid_projs)
                    end do
                    !$omp end parallel do
                else
                    !$omp parallel do default(shared) schedule(guided) private(iptcl) proc_bind(close)
                    do iptcl=p%fromp,p%top
                        if( to_update(iptcl) )then
                            call primesrch3D(iptcl)%exec_prime3D_srch(p%lp,&
                                &greedy=.true., nnmat=b%nnmat, grid_projs=b%grid_projs)
                        endif
                    end do
                    !$omp end parallel do
                endif
            case('het')
                if(p%oritab .eq. '') stop 'cannot run the refine=het mode without input oridoc (oritab)'
                if( corr_thresh > TINY )then
                    write(*,'(A,F8.2)') '>>> PARTICLE RANDOMIZATION(%):', 100.*extr_thresh
                    write(*,'(A,F8.2)') '>>> CORRELATION THRESHOLD:    ', corr_thresh
                endif
                !$omp parallel do default(shared) schedule(guided) private(iptcl) reduction(+:statecnt) proc_bind(close)
                do iptcl=p%fromp,p%top
                    call primesrch3D(iptcl)%exec_prime3D_srch_het(corr_thresh, statecnt)
                end do
                !$omp end parallel do
                if( corr_thresh > TINY )then
                    norm = real(sum(statecnt))
                    do istate=1,p%nstates
                        print *,'% randomized ptcls for state ',istate,' is ',100.*(real(statecnt(istate))/norm),&
                            &'; pop=',statecnt(istate)
                    end do
                endif
            case ('states')
                if(p%oritab .eq. '') stop 'cannot run the refine=states mode without input oridoc (oritab)'
                !$omp parallel do default(shared) schedule(guided) private(iptcl) proc_bind(close)
                do iptcl=p%fromp,p%top
                    if( to_update(iptcl) )then
                        call primesrch3D(iptcl)%exec_prime3D_srch(p%lp, greedy=.true., nnmat=b%nnmat)
                    endif
                end do
                !$omp end parallel do
            case ('tseries')
                if(p%oritab .eq. '') stop 'cannot run the refine=tseries mode without input oridoc (oritab)'
                !$omp parallel do default(shared) schedule(guided) private(iptcl) proc_bind(close)
                do iptcl=p%fromp,p%top
                    if( to_update(iptcl) )then
                        call primesrch3D(iptcl)%exec_prime3D_srch(p%lp, greedy=.false.)
                    endif
                end do
                !$omp end parallel do
            case DEFAULT
                write(*,*) 'The refinement mode: ', trim(p%refine), ' is unsupported'
                stop
        end select
        if( L_BENCH ) rt_align = toc(t_align)
        

        ! PARTICLE REJECTION BASED ON ALIGNABILITY (SDEV OF ANGULAR ORIS)
        if( cline%defined('sdev_thres') )then
            call b%a%reject_above('sdev', p%sdev_thres)
        endif

        ! OUTPUT ORIENTATIONS
        call binwrite_oritab(p%outfile, b%a, [p%fromp,p%top])
        p%oritab = p%outfile

        ! VOLUMETRIC 3D RECONSTRUCTION
        if( L_BENCH ) t_rec = tic()
        if( p%norec .ne. 'yes' )then
            ! init volumes
            call preprecvols(b, p)
            if( p%l_frac_update )then
                ! need to read in part volumes & rho
                ! todo in simple_hadamard_common
                ! call readrecvols_for_update(p)
                ! discarding contribution to volume should be done within grid_ptcl
            endif
            ! reconstruction
            do iptcl=p%fromp,p%top
                if( to_update(iptcl) )then
                    orientation = b%a%get_ori(iptcl)
                    if( nint(orientation%get('state')) == 0 .or.&
                       &nint(orientation%get('state_balance')) == 0 ) cycle
                    call read_img_and_norm( b, p, iptcl )
                    if( p%npeaks > 1 )then
                        call primesrch3D(iptcl)%get_oris(prime3D_oris, orientation)
                        call grid_ptcl(b, p, orientation, os=prime3D_oris)
                    else
                        call grid_ptcl(b, p, orientation)
                    endif
                else
                    ! reconstruction contribution already part of recvols
                endif
            end do
            ! normalise structure factors
            if( p%eo .ne. 'no' )then
                call eonorm_struct_facts(b, p, cline, reslim, which_iter)
            else
                call norm_struct_facts(b, p, which_iter)
            endif
        endif
        if( L_BENCH ) rt_rec = toc(t_rec)

        ! DESTRUCT
        if( .not. p%l_distr_exec )then
            do iptcl=p%fromp,p%top
                call primesrch3D(iptcl)%kill
                call reforis(iptcl)%kill
            end do
            call pftcc%kill
            deallocate( primesrch3D, state_exists, proj_space_inds, reforis )
            call prime3D_oris%kill
            call killrecvols(b, p)
        endif

        ! REPORT CONVERGENCE
        if( p%l_distr_exec )then
            call qsys_job_finished( p, 'simple_hadamard3D_matcher :: prime3D_exec')
        else
            if( p%refine .eq. 'het' )then
                converged = b%conv%check_conv_het()
            else
                converged = b%conv%check_conv3D(update_res)
            endif
        endif
        if( L_BENCH )then
            rt_tot  = toc(t_tot)
            doprint = .true.
            if( p%l_distr_exec .and. p%part /= 1 ) doprint = .false.
            if( doprint )then
                benchfname = 'HADAMARD3D_BENCH_ITER'//int2str_pad(which_iter,3)//'.txt'
                call fopen(fnr, FILE=trim(benchfname), STATUS='REPLACE', action='WRITE')
                write(fnr,'(a)') '*** TIMINGS (s) ***'
                write(fnr,'(a,1x,f9.2)') 'initialisation       : ', rt_init
                write(fnr,'(a,1x,f9.2)') 'pftcc preparation    : ', rt_prep_pftcc
                write(fnr,'(a,1x,f9.2)') 'stochastic alignment : ', rt_align
                write(fnr,'(a,1x,f9.2)') 'reconstruction       : ', rt_rec
                write(fnr,'(a,1x,f9.2)') 'total time           : ', rt_tot
                write(fnr,'(a)') ''
                write(fnr,'(a)') '*** REATIVE TIMINGS (%) ***'
                write(fnr,'(a,1x,f9.2)') 'initialisation       : ', (rt_init/rt_tot)        * 100.
                write(fnr,'(a,1x,f9.2)') 'pftcc preparation    : ', (rt_prep_pftcc/rt_tot)  * 100.
                write(fnr,'(a,1x,f9.2)') 'stochastic alignment : ', (rt_align/rt_tot)       * 100.
                write(fnr,'(a,1x,f9.2)') 'reconstruction       : ', (rt_rec/rt_tot)        * 100.
                write(fnr,'(a,1x,f9.2)') '% accounted for      : ',&
                    &((rt_init+rt_prep_pftcc+rt_align+rt_rec)/rt_tot) * 100.
                call fclose(fnr)
            endif
        endif
    end subroutine prime3D_exec

    subroutine gen_random_model( b, p, nsamp_in )
        use simple_ran_tabu,   only: ran_tabu
        use simple_kbinterpol, only: kbinterpol
        class(build),      intent(inout) :: b         !< build object
        class(params),     intent(inout) :: p         !< param object
        integer, optional, intent(in)    :: nsamp_in  !< num input samples
        type(ran_tabu)       :: rt
        type(ori)            :: orientation, o_sym
        integer, allocatable :: sample(:)
        integer              :: i, k, nsamp, alloc_stat
        type(kbinterpol)     :: kbwin
        if( p%vols(1) == '' )then
            ! init volumes
            call preprecvols(b, p)
            p%oritab = 'prime3D_startdoc'//METADATEXT
            if( trim(p%refine).eq.'tseries' )then
                call b%a%spiral
            else
                call b%a%rnd_oris
                call b%a%zero_shifts
            endif
            if( p%l_distr_exec .and. p%part.ne.1 )then
                ! so random oris only written once in distributed mode
            else
                call binwrite_oritab(p%oritab, b%a, [1,p%nptcls])
            endif
            p%vols(1) = 'startvol'//p%ext
            if( p%noise .eq. 'yes' )then
                call b%vol%ran
                call b%vol%write(p%vols(1), del_if_exists=.true.)
                return
            endif
            nsamp = p%top - p%fromp + 1
            if( present(nsamp_in) ) nsamp = nsamp_in
            allocate( sample(nsamp), stat=alloc_stat )
            call alloc_errchk("In: gen_random_model; simple_hadamard3D_matcher", alloc_stat)
            if( present(nsamp_in) )then
                rt = ran_tabu(p%top - p%fromp + 1)
                call rt%ne_ran_iarr(sample)
                call rt%kill
            else
                forall(i=1:nsamp) sample(i) = i
            endif
            write(*,'(A)') '>>> RECONSTRUCTING RANDOM MODEL'
            kbwin = b%recvols(1)%get_kbwin()
            do i=1,nsamp
                call progress(i, nsamp)
                orientation = b%a%get_ori(sample(i) + p%fromp - 1)
                call read_img_and_norm(b, p, sample(i) + p%fromp - 1)
                call prep4cgrid(b%img, b%img_pad, p%msk, kbwin)
                if( p%pgrp == 'c1' )then
                    call b%recvols(1)%inout_fplane(orientation, .true., b%img_pad, pwght=1.0)
                else
                    do k=1,b%se%get_nsym()
                        o_sym = b%se%apply(orientation, k)
                        call b%recvols(1)%inout_fplane(o_sym, .true., b%img_pad, pwght=1.0)
                    end do
                endif
            end do
            deallocate(sample)
            call norm_struct_facts(b, p)
            call killrecvols(b, p)
        endif
    end subroutine gen_random_model
    
    !> Prepare alignment search using polar projection Fourier cross correlation
    subroutine preppftcc4align( b, p, cline, ppfts_fname )
        class(build),               intent(inout) :: b       !< build object
        class(params),              intent(inout) :: p       !< param object
        class(cmdline),             intent(inout) :: cline   !< command line
        character(len=*), optional, intent(in)    :: ppfts_fname
        integer :: nrefs
        if( .not. p%l_distr_exec ) write(*,'(A)') '>>> BUILDING PRIME3D SEARCH ENGINE'
        ! must be done here since p%kfromto is dynamically set based on FSC from previous round
        ! or based on dynamic resolution limit update
        nrefs = p%nspace * p%nstates
        if( p%eo .ne. 'no' )then
            call pftcc%new(nrefs, p, nint(b%a%get_all('eo', [p%fromp,p%top])))
        else
            call pftcc%new(nrefs, p)
        endif
        call prep_refs_pftcc4align( b, p, cline )
        call prep_ptcls_pftcc4align( b, p, ppfts_fname )
        DebugPrint '*** hadamard3D_matcher ***: finished preppftcc4align'
    end subroutine preppftcc4align

    !> Prepare reference images and create polar projections
    subroutine prep_refs_pftcc4align( b, p, cline )
        class(build),   intent(inout) :: b          !< build object
        class(params),  intent(inout) :: p          !< param object
        class(cmdline), intent(inout) :: cline      !< command line
        type(ori) :: o
        integer   :: cnt, s, iref, nrefs, ldim(3)
        logical   :: do_center
        real      :: xyz(3)
        ! PREPARATION OF REFERENCES IN PFTCC
        ! read reference volumes and create polar projections
        nrefs = p%nspace * p%nstates
        cnt   = 0
        if( .not. p%l_distr_exec ) write(*,'(A)') '>>> BUILDING REFERENCES'
        do s=1,p%nstates
            if( p%oritab .ne. '' )then
                if( b%a%get_pop(s, 'state') == 0 )then
                    ! empty state
                    cnt = cnt + p%nspace
                    call progress(cnt, nrefs)
                    cycle
                endif
            endif
            call cenrefvol_and_mapshifts2ptcls(b, p, cline, s, p%vols(s), do_center, xyz)
            if( p%eo .ne. 'no' )then
                if( .not. p%l_distr_exec ) write(*,'(A)') '>>> BUILDING EVEN REFERENCES'
                call preprefvol(b, p, cline, s, p%vols_even(s), do_center, xyz)
                cnt = 0
                do iref=1,p%nspace
                    cnt = cnt + 1
                    call progress(cnt, p%nspace)
                    o = b%e%get_ori(iref)
                    call b%vol%fproject_polar(cnt, o, pftcc, iseven=.true.) ! polar central sections
                end do
                if( .not. p%l_distr_exec ) write(*,'(A)') '>>> BUILDING ODD REFERENCES'
                call preprefvol(b, p, cline, s, p%vols_odd(s), do_center, xyz)
                cnt = 0
                do iref=1,p%nspace
                    cnt = cnt + 1
                    call progress(cnt, p%nspace)
                    o = b%e%get_ori(iref)
                    call b%vol%fproject_polar(cnt, o, pftcc, iseven=.false.) ! polar central sections
                end do
            else
                call preprefvol(b, p, cline, s, p%vols(s), do_center, xyz )
                do iref=1,p%nspace
                    cnt = cnt + 1
                    call progress(cnt, nrefs)
                    o = b%e%get_ori(iref)
                    call b%vol%fproject_polar(cnt, o, pftcc, iseven=.true.) ! polar central sections
                end do
            endif
        end do
        ! cleanup
        call b%vol%kill_expanded
        ! bring back the original b%vol size for clean exit
        call b%vol%new([p%box,p%box,p%box], p%smpd)
    end subroutine prep_refs_pftcc4align

    !> Prepare particle images and create polar projections
    subroutine prep_ptcls_pftcc4align( b, p, ppfts_fname )
        use simple_fileio, only: file_exists
        class(build),               intent(inout) :: b          !< build object
        class(params),              intent(inout) :: p          !< param object
        character(len=*), optional, intent(in)    :: ppfts_fname
        type(ori) :: o
        integer   :: cnt, s, iptcl, istate, ntot
        if( .not. p%l_distr_exec ) write(*,'(A)') '>>> BUILDING PARTICLES'
        ! initialize
        call b%img_match%init_polarizer(pftcc, p%alpha)
        ntot = (p%top-p%fromp+1) * p%nstates
        cnt  = 0
        do s=1,p%nstates
            if( b%a%get_pop(s, 'state') == 0 )then
                ! empty state
                cycle
            endif
            do iptcl=p%fromp,p%top
                o      = b%a%get_ori(iptcl)
                istate = nint(o%get('state'))
                cnt = cnt + 1
                if( istate /= s ) cycle
                call progress(cnt, ntot)
                call read_img_and_norm( b, p, iptcl )
                call prepimg4align(b, p, o, is3D=.true.)
                call b%img_match%polarize(pftcc, iptcl)
            end do
        end do
    end subroutine prep_ptcls_pftcc4align

end module simple_hadamard3D_matcher
