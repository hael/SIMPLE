! projection-matching based on Hadamard products, high-level search routines for PRIME3D
module simple_hadamard3D_matcher
!$ use omp_lib
!$ use omp_lib_kinds
use simple_defs
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_prime3D_srch,     only: prime3D_srch
use simple_ori,              only: ori
use simple_build,            only: build
use simple_params,           only: params
use simple_cmdline,          only: cmdline
use simple_gridding,         only: prep4cgrid
use simple_strings,          only: str_has_substr, int2str_pad
use simple_cont3D_matcher    ! use all in there
use simple_hadamard_common   ! use all in there
use simple_math              ! use all in there
implicit none

public :: prime3D_find_resrange, prime3D_exec, gen_random_model
public :: preppftcc4align, prep_refs_pftcc4align, pftcc

private
#include "simple_local_flags.inc"

type(polarft_corrcalc)          :: pftcc
type(prime3D_srch), allocatable :: primesrch3D(:)
real                            :: reslim

type(ori)                       :: orientation, o_sym
character(len=:), allocatable   :: ppfts_fname

contains

    !> Find resolution range in Prime3D search
    subroutine prime3D_find_resrange( b, p, lp_start, lp_finish )
        use simple_oris, only: oris
        class(build),  intent(inout) :: b
        class(params), intent(inout) :: p
        real,          intent(out)   :: lp_start, lp_finish
        real, allocatable :: peaks(:)
        type(oris)        :: o
        integer :: lfny, alloc_stat, k, pos10, pos6
        call o%new(p%nspace)
        call o%spiral
        lfny = b%img_match%get_lfny(1)
        allocate( peaks(lfny), stat=alloc_stat )
        call alloc_err("In: prime3D_find_resrange, simple_hadamard3D_matcher", alloc_stat)
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

    !> Execute prime3D (Hadamard method)
    subroutine prime3D_exec( b, p, cline, which_iter, update_res, converged )
        use simple_qsys_funs, only: qsys_job_finished
        use simple_oris,      only: oris
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: which_iter
        logical,        intent(inout) :: update_res, converged
        type(oris) :: prime3D_oris
        real       :: norm, corr_thresh, skewness, frac_srch_space, extr_thresh
        integer    :: iptcl, inptcls, istate
        integer    :: statecnt(p%nstates)
        inptcls = p%top - p%fromp + 1

        ! SET BAND-PASS LIMIT RANGE
        call set_bp_range( b, p, cline )

        ! CALCULATE ANGULAR THRESHOLD (USED BY THE SPARSE WEIGHTING SCHEME)
        p%athres = rad2deg( atan(max(p%fny,p%lp)/(p%moldiam/2.) ))
        reslim   = p%lp
        DebugPrint '*** hadamard3D_matcher ***: calculated angular threshold (used by the sparse weighting scheme)'

        ! DETERMINE THE NUMBER OF PEAKS
        if( .not. cline%defined('npeaks') )then
            select case(p%refine)
                case('no', 'neigh', 'greedy', 'greedyneigh', 'exp')
                    if( p%eo .eq. 'yes' )then
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
        if( p%vols(1) .eq. '' .and. p%nstates==1 )then
            if( p%nptcls > 1000 )then
                call gen_random_model(b, p, 1000)
            else
                call gen_random_model(b, p)
            endif
            DebugPrint '*** hadamard3D_matcher ***: generated random model'
        endif

        ! SET FRACTION OF SEARCH SPACE
        frac_srch_space = b%a%get_avg('frac')

        ! EXTREMAL LOGICS
        if( p%refine.eq.'het' )then
            if( frac_srch_space < 98. .or. which_iter <= 15 )then
            ! if( frac_srch_space < 98. .or. extr_thresh > 0.025 )then
                !extr_thresh = p%extr_thresh                                          ! factorial decay: the old way
                !extr_thresh = EXTRINITHRESH * exp(-(real(which_iter-1)/6.)**2. / 2.) ! gaussian decay: untested
                extr_thresh = EXTRINITHRESH * cos(PI/2. * real(which_iter-1)/15.)   ! cosine decay
                extr_thresh = max(0., extr_thresh)
                extr_thresh = min(EXTRINITHRESH, extr_thresh)
                corr_thresh = b%a%extremal_bound(extr_thresh)
                statecnt(:) = 0
            else
                corr_thresh = -huge(corr_thresh)
            endif
        endif

        ! PREPARE THE POLARFT_CORRCALC DATA STRUCTURE
        if( p%refine.eq.'het' )then
            ! generate filename for memoization of particle pfts
            if( allocated(ppfts_fname) ) deallocate(ppfts_fname)
            if( p%l_distr_exec )then
                allocate( ppfts_fname, source='ppfts_memoized_part'//int2str_pad(p%part,p%numlen)//'.bin' )
            else
                allocate( ppfts_fname, source='ppfts_memoized.bin' )
            endif
            ! generate projections (polar FTs)
            call preppftcc4align( b, p, cline, ppfts_fname )
        else
            ! generate projections (polar FTs)
            call preppftcc4align( b, p, cline )
        endif

        ! INITIALIZE
        write(*,'(A,1X,I3)') '>>> PRIME3D DISCRETE STOCHASTIC SEARCH, ITERATION:', which_iter
        if( .not. p%l_distr_exec )then
            if( p%refine .eq. 'snhc')then
                p%outfile = SNHCDOC
            else
                p%outfile = 'prime3Ddoc_'//int2str_pad(which_iter,3)//'.txt'
            endif
        endif

        ! STOCHASTIC IMAGE ALIGNMENT
        ! create the search objects, need to re-create every round because parameters are changing
        allocate( primesrch3D(p%fromp:p%top) )
        do iptcl=p%fromp,p%top
            call primesrch3D(iptcl)%new(b%a, p, pftcc)
        end do
        ! prep ctf
        if(p%ctf .ne. 'no') call pftcc%create_polar_ctfmats(b%a)
        ! execute the search
        call del_file(p%outfile)
        select case(p%refine)
            case( 'snhc' )
                !$omp parallel do default(shared) schedule(guided) private(iptcl) proc_bind(close)
                do iptcl=p%fromp,p%top
                    call primesrch3D(iptcl)%exec_prime3D_srch(pftcc, iptcl, b%a, b%e, p%lp, szsn=p%szsn)
                end do
                !$omp end parallel do
            case( 'no','shc' )
                if( p%oritab .eq. '' )then
                    !$omp parallel do default(shared) schedule(guided) private(iptcl) proc_bind(close)
                    do iptcl=p%fromp,p%top
                        call primesrch3D(iptcl)%exec_prime3D_srch(pftcc, iptcl, b%a, b%e, p%lp, greedy=.true.)
                    end do
                    !$omp end parallel do
                else
                    !$omp parallel do default(shared) schedule(guided) private(iptcl) proc_bind(close)
                    do iptcl=p%fromp,p%top
                        call primesrch3D(iptcl)%exec_prime3D_srch(pftcc, iptcl, b%a, b%e, p%lp)
                    end do
                    !$omp end parallel do
                endif
            case('neigh','shcneigh')
                if( p%oritab .eq. '' ) stop 'cannot run the refine=neigh mode without input oridoc (oritab)'
                !$omp parallel do default(shared) schedule(guided) private(iptcl) proc_bind(close)
                do iptcl=p%fromp,p%top
                    call primesrch3D(iptcl)%exec_prime3D_srch(pftcc, iptcl, b%a,&
                        b%e, p%lp, nnmat=b%nnmat, grid_projs=b%grid_projs)
                end do
                !$omp end parallel do
            case('greedy')
                !$omp parallel do default(shared) schedule(guided) private(iptcl) proc_bind(close)
                do iptcl=p%fromp,p%top
                    call primesrch3D(iptcl)%exec_prime3D_srch(pftcc, iptcl, b%a, b%e, p%lp, greedy=.true.)
                end do
                !$omp end parallel do
            case('greedyneigh')
                !$omp parallel do default(shared) schedule(guided) private(iptcl) proc_bind(close)
                do iptcl=p%fromp,p%top
                    call primesrch3D(iptcl)%exec_prime3D_srch(pftcc, iptcl, b%a, b%e, p%lp,&
                        greedy=.true., nnmat=b%nnmat, grid_projs=b%grid_projs)
                end do
                !$omp end parallel do
            case('het')
                if(p%oritab .eq. '') stop 'cannot run the refine=het mode without input oridoc (oritab)'
                if( corr_thresh > TINY )then
                    ! write(*,'(A,F8.2)') '>>> PARTICLE RANDOMIZATION(%):', 100.*p%extr_thresh
                    write(*,'(A,F8.2)') '>>> PARTICLE RANDOMIZATION(%):', 100.*extr_thresh
                    write(*,'(A,F8.2)') '>>> CORRELATION THRESHOLD:    ', corr_thresh
                endif
                !$omp parallel do default(shared) schedule(guided) private(iptcl) reduction(+:statecnt) proc_bind(close)
                do iptcl=p%fromp,p%top
                    call primesrch3D(iptcl)%exec_prime3D_srch_het(pftcc, iptcl, b%a, b%e, corr_thresh, statecnt)
                end do
                !$omp end parallel do
                if( corr_thresh > TINY )then
                    norm = real(sum(statecnt))
                    do istate=1,p%nstates
                        print *,'% randomized ptcls for state ',istate,' is ',100.*(real(statecnt(istate))/norm),&
                            &'; pop=',statecnt(istate)
                    end do
                endif
            case ('exp')
                !$omp parallel do default(shared) schedule(guided) private(iptcl) proc_bind(close)
                do iptcl=p%fromp,p%top
                    call primesrch3D(iptcl)%exec_prime3D_srch(pftcc, iptcl, b%a, b%e, p%lp,&
                        greedy=.true., nnmat=b%nnmat, grid_projs=b%grid_projs)
                end do
                !$omp end parallel do
            case DEFAULT
                write(*,*) 'The refinement mode: ', trim(p%refine), ' is unsupported'
                stop
        end select
        call pftcc%kill

        ! SETUP WEIGHTS
        if( p%nptcls <= SPECWMINPOP )then
            call b%a%calc_hard_ptcl_weights(p%frac)
        else
            call b%a%calc_spectral_weights(p%frac)
        endif

        ! POPULATION BALANCING LOGICS
        if( p%balance > 0 )then
            call b%a%balance('proj', p%balance, skewness)
            write(*,'(A,F8.2)') '>>> PROJECTION DISTRIBUTION SKEWNESS(%):', 100. * skewness
        else
            call b%a%set_all2single('state_balance', 1.0)
        endif

        ! OUTPUT ORIENTATIONS
        call b%a%write(p%outfile, [p%fromp,p%top])
        p%oritab = p%outfile

        ! VOLUMETRIC 3D RECONSTRUCTION
        if( p%norec .eq. 'no' )then
            ! init volumes
            call preprecvols(b, p)
            ! reconstruction
            do iptcl=p%fromp,p%top
                orientation = b%a%get_ori(iptcl)
                if( nint(orientation%get('state')) == 0 .or.&
                   &nint(orientation%get('state_balance')) == 0 ) cycle
                call read_img_from_stk( b, p, iptcl )
                if( p%npeaks > 1 )then
                    call primesrch3D(iptcl)%get_oris(prime3D_oris, orientation)
                    call grid_ptcl(b, p, orientation, prime3D_oris)
                else
                    call grid_ptcl(b, p, orientation)
                endif
            end do
            ! normalise structure factors
            if( p%eo .eq. 'yes' )then
                call eonorm_struct_facts(b, p, reslim, which_iter)
            else
                call norm_struct_facts(b, p, which_iter)
            endif
            ! destruct volumes
            call killrecvols(b, p)
        endif

        ! DESTRUCT
        do iptcl=p%fromp,p%top
            call primesrch3D(iptcl)%kill
        end do
        deallocate( primesrch3D )
        call prime3D_oris%kill

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
    end subroutine prime3D_exec

    subroutine gen_random_model( b, p, nsamp_in )
        use simple_ran_tabu,   only: ran_tabu
        use simple_kbinterpol, only: kbinterpol
        class(build),      intent(inout) :: b         !< build object
        class(params),     intent(inout) :: p         !< param object
        integer, optional, intent(in)    :: nsamp_in  !< num input samples
        type(ran_tabu)       :: rt
        integer, allocatable :: sample(:)
        integer              :: i, k, nsamp, alloc_stat
        type(kbinterpol)     :: kbwin
        if( p%vols(1) == '' )then
            ! init volumes
            call preprecvols(b, p)
            p%oritab = 'prime3D_startdoc.txt'
            call b%a%rnd_oris
            call b%a%zero_shifts
            if( p%l_distr_exec .and. p%part.ne.1 )then
                ! so random oris only written once in distributed mode
            else
                call b%a%write( p%oritab )
            endif
            p%vols(1) = 'startvol'//p%ext
            if( p%noise .eq. 'yes' )then
                call b%vol%ran
                call b%vol%write(p%vols(1), del_if_exists=.true.)
                return
            endif
            nsamp = p%nptcls
            if( present(nsamp_in) ) nsamp = nsamp_in
            allocate( sample(nsamp), stat=alloc_stat )
            call alloc_err("In: gen_random_model; simple_hadamard3D_matcher", alloc_stat)
            if( present(nsamp_in) )then
                rt = ran_tabu(p%nptcls)
                call rt%ne_ran_iarr(sample)
                call rt%kill
            else
                forall(i=1:nsamp) sample(i) = i
            endif
            write(*,'(A)') '>>> RECONSTRUCTING RANDOM MODEL'
            kbwin = b%recvols(1)%get_kbwin()
            do i=1,nsamp
                call progress(i, nsamp)
                orientation = b%a%get_ori(sample(i))
                call b%img%read(p%stk, sample(i))
                call prep4cgrid(b%img, b%img_pad, p%msk, kbwin)
                if( p%pgrp == 'c1' )then
                    call b%recvols(1)%inout_fplane(orientation, .true., b%img_pad)
                else
                    do k=1,b%se%get_nsym()
                        o_sym = b%se%apply(orientation, k)
                        call b%recvols(1)%inout_fplane(o_sym, .true., b%img_pad)
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
        nrefs = p%nspace*p%nstates
        call pftcc%new(nrefs, [p%fromp,p%top], [p%boxmatch,p%boxmatch,1],&
        p%smpd, p%kfromto, p%ring2, p%ctf)
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
        integer   :: cnt, s, iref, nrefs
        ! PREPARATION OF REFERENCES IN PFTCC
        ! read reference volumes and create polar projections
        nrefs = p%nspace*p%nstates
        cnt   = 0
        if( .not. p%l_distr_exec ) write(*,'(A)') '>>> BUILDING REFERENCES'
        do s=1,p%nstates
            if( p%oritab .ne. '' )then
                ! greedy start
                if( b%a%get_state_pop(s) == 0 )then
                    ! empty state
                    cnt = cnt + p%nspace
                    call progress(cnt, nrefs)
                    cycle
                endif
            endif
            call preprefvol( b, p, cline, s )
            ! generate discrete projections
            do iref=1,p%nspace
                cnt = cnt+1
                call progress(cnt, nrefs)
                o = b%e%get_ori(iref)
                call b%vol%fproject_polar(cnt, o, pftcc)
            end do
        end do
        ! cleanup
        call b%vol%kill_expanded
        ! bring back the original b%vol size for clean exit
        if( p%boxmatch < p%box )call b%vol%new([p%box,p%box,p%box], p%smpd)
    end subroutine prep_refs_pftcc4align

    !> Prepare particle images and create polar projections
    subroutine prep_ptcls_pftcc4align( b, p, ppfts_fname )
        class(build),               intent(inout) :: b          !< build object
        class(params),              intent(inout) :: p          !< param object
        character(len=*), optional, intent(in)    :: ppfts_fname
        ! read particle images and create polar projections
        if( present(ppfts_fname) )then
            if( file_exists(ppfts_fname) )then
                call pftcc%read_pfts_ptcls(ppfts_fname)
            else
                call prep_pftcc_local
                call pftcc%write_pfts_ptcls(ppfts_fname)
            endif
        else
            call prep_pftcc_local
        endif

        contains

            subroutine prep_pftcc_local
                type(ori) :: o
                integer   :: cnt, s, iptcl, istate, ntot
                if( .not. p%l_distr_exec ) write(*,'(A)') '>>> BUILDING PARTICLES'
                ! initialize
                call b%img_match%init_polarizer(pftcc)
                ntot = (p%top-p%fromp+1) * p%nstates
                cnt  = 0
                do s=1,p%nstates
                    if( b%a%get_state_pop(s) == 0 )then
                        ! empty state
                        cycle
                    endif
                    do iptcl=p%fromp,p%top
                        o      = b%a%get_ori(iptcl)
                        istate = nint(o%get('state'))
                        cnt = cnt + 1
                        if( istate /= s ) cycle
                        call progress(cnt, ntot)
                        call read_img_from_stk( b, p, iptcl )
                        call prepimg4align(b, p, o)
                        call b%img_match%polarize(pftcc, iptcl)
                    end do
                end do
            end subroutine prep_pftcc_local

    end subroutine prep_ptcls_pftcc4align

end module simple_hadamard3D_matcher
