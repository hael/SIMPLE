module simple_hadamard3D_matcher
use simple_defs
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_prime3D_srch,     only: prime3D_srch
use simple_ori,              only: ori
use simple_build,            only: build
use simple_params,           only: params
use simple_cmdline,          only: cmdline
use simple_gridding,         only: prep4cgrid
use simple_masker,           only: automask
use simple_rnd,              only: ran3
use simple_strings,          only: str_has_substr
use simple_cont3D_matcher    ! use all in there
use simple_hadamard_common   ! use all in there
use simple_math              ! use all in there
implicit none

public :: prime3D_exec, gen_random_model, prime3D_find_resrange, pftcc, primesrch3D
public :: preppftcc4align, prep_refs_pftcc4align
private

integer, parameter            :: MAXNPEAKS=10
logical, parameter            :: debug=.false.

type(polarft_corrcalc)        :: pftcc
type(prime3D_srch)            :: primesrch3D
real                          :: reslim
real                          :: frac_srch_space
type(ori)                     :: orientation, o_sym, o_tmp
integer                       :: cnt_glob=0
character(len=:), allocatable :: ppfts_fname

contains

    !>  \brief  to find the resolution range for initial prime3D alignment
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
        lfny = b%img%get_lfny(1)
        allocate( peaks(lfny), stat=alloc_stat )
        call alloc_err("In: prime3D_find_resrange, simple_hadamard3D_matcher", alloc_stat)
        do k=2,b%img%get_lfny(1)
            peaks(k) = real(o%find_npeaks(b%img%get_lp(k), p%moldiam))
        end do
        peaks(1)  = peaks(2)
        pos10     = locate(peaks, lfny, 10.)
        pos6      = locate(peaks, lfny,  6.)
        lp_start  = b%img%get_lp(pos10)
        lp_finish = b%img%get_lp(pos6)
        deallocate(peaks)
        call o%kill
    end subroutine prime3D_find_resrange

    !>  \brief  is the prime3D algorithm
    subroutine prime3D_exec( b, p, cline, which_iter, update_res, converged )
        use simple_filterer,  only: resample_filter
        use simple_qsys_funs, only: qsys_job_finished
        use simple_oris,      only: oris
        use simple_strings,   only: int2str_pad
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: which_iter
        logical,        intent(inout) :: update_res, converged
        type(oris)                    :: prime3D_oris, exist_oris
        real, allocatable             :: wmat(:,:), wresamp(:), res(:), res_pad(:), corrs(:), corrs_incl(:)
        integer, allocatable          :: smpl_inds(:), inds(:), inds_incl(:)
        logical, allocatable          :: is_rnd(:), incl(:)
        real                          :: w, lptmp, norm, het_corr_thresh
        integer                       :: iptcl, fnr, file_stat, s, inptcls, prev_state, istate
        integer                       :: statecnt(p%nstates), ind, thresh_ind, n_samples
        logical                       :: doshellweight

        inptcls = p%top - p%fromp + 1

        ! SET FRACTION OF SEARCH SPACE
        frac_srch_space = b%a%get_avg('frac')

        ! SET BAND-PASS LIMIT RANGE
        call set_bp_range( b, p, cline )

        ! CALCULATE ANGULAR THRESHOLD (USED BY THE SPARSE WEIGHTING SCHEME)
        p%athres = rad2deg(atan(max(p%fny,p%lp)/(p%moldiam/2.)))
        reslim   = p%lp
        if( debug ) write(*,*) '*** hadamard3D_matcher ***: calculated angular threshold (used by the sparse weighting scheme)'

        ! DETERMINE THE NUMBER OF PEAKS
        if( .not. cline%defined('npeaks') )then
            select case(p%refine)
                case('no', 'neigh', 'adasym')
                    p%npeaks = min(MAXNPEAKS,b%e%find_npeaks(p%lp, p%moldiam))
                    if( str_has_substr(p%refine,'adasym') ) p%npeaks = p%npeaks * p%nsym
                case DEFAULT
                    p%npeaks = 1
            end select
            if( debug ) write(*,*) '*** hadamard3D_matcher ***: determined the number of peaks'
        endif

        ! RANDOM MODEL GENERATION
        if( p%vols(1) .eq. '' .and. p%nstates==1 )then
            if( p%nptcls > 1000 )then
                call gen_random_model(b, p, 1000)
            else
                call gen_random_model(b, p)
            endif
            if( debug ) write(*,*) '*** hadamard3D_matcher ***: generated random model'
        endif

        ! SETUP WEIGHTS FOR THE 3D RECONSTRUCTION
        if( p%oritab .ne. '' .and. p%frac < 0.99 ) call b%a%calc_hard_ptcl_weights(p%frac)
        if( p%l_distr_exec )then
            ! nothing to do
        else
            if( p%l_shellw .and. frac_srch_space >= 50. ) call cont3D_shellweight(b, p, cline)
        endif
        call setup_shellweights(b, p, doshellweight, wmat, res, res_pad)

        ! HETEROGEINITY
        if( p%refine.eq.'het' )then
            het_corr_thresh = -1.
            if( frac_srch_space<0.98 .or. p%het_thresh>0.02 )then
                write(*,'(A,F8.2)')'>>> STATE RANDOMIZATION %:', 100.*p%het_thresh
                corrs      = b%a%get_all('corr')
                incl       = b%a%included()
                inds       = b%a%order()
                inds_incl  = pack(inds, mask=incl)
                corrs_incl = pack(corrs, mask=incl)
                thresh_ind = nint( size(inds_incl)*p%het_thresh )
                het_corr_thresh = corrs_incl( inds_incl(thresh_ind) )
                write(*,'(A,F8.2)')'>>> CORRELATION THRESHOLD:', het_corr_thresh
                allocate(is_rnd(p%nptcls))
                is_rnd = .false.
                where( corrs <= het_corr_thresh )is_rnd = .true.
                deallocate(corrs, inds, incl, inds_incl, corrs_incl)
            endif

            ! generate filename for memoization of particle pfts
            if( allocated(ppfts_fname) ) deallocate(ppfts_fname)
            if( p%l_distr_exec )then
                allocate( ppfts_fname, source='ppfts_memoized_part'//int2str_pad(p%part,p%numlen)//'.bin' )
            else
                allocate( ppfts_fname, source='ppfts_memoized'//'.bin' )
            endif
            ! generate projections (polar FTs)
            call preppftcc4align( b, p, cline, ppfts_fname )
        else
            ! generate projections (polar FTs)
            call preppftcc4align( b, p, cline )
        endif

        ! INITIALIZE
        if( which_iter <= 0 )then
            write(*,'(A)') '>>> PRIME3D DISCRETE STOCHASTIC SEARCH'
        else
            write(*,'(A,1X,I3)') '>>> PRIME3D DISCRETE STOCHASTIC SEARCH, ITERATION:', which_iter
        endif
        if( which_iter > 0 ) p%outfile = 'prime3Ddoc_'//int2str_pad(which_iter,3)//'.txt'

        ! FOR BENCHMARKING (GPU LOGICS ON CPU) OR USE OF GPU WE CALCULATE CORRS BEFORHAND
        if( p%use_gpu .eq. 'yes' .or. p%bench_gpu .eq. 'yes' )then
           if( p%bench_gpu .eq. 'yes' .and. p%use_gpu .eq. 'no' )then
              call primesrch3D%calc_corrs(pftcc, mode='bench')
           else
              call primesrch3D%calc_corrs(pftcc)
           endif
        endif

        ! RESET RECVOLS
        do s=1,p%nstates
            if( p%eo .eq. 'yes' )then
                call b%eorecvols(s)%reset_all
            else
                call b%recvols(s)%reset
            endif
        end do
        if( debug ) write(*,*) '*** hadamard3D_matcher ***: did reset recvols'

        ! ALIGN & GRID
        call del_file(p%outfile)
        cnt_glob = 0
        statecnt = 0
        if( debug ) write(*,*) '*** hadamard3D_matcher ***: loop fromp/top:', p%fromp, p%top
        do iptcl=p%fromp,p%top
            cnt_glob = iptcl - p%fromp + 1
            call progress( cnt_glob, inptcls )
            orientation = b%a%get_ori(iptcl)
            prev_state  = nint( orientation%get('state') )
            if( prev_state > 0 )then
                call preprefs4align(b, p, iptcl, pftcc)
                ! execute the high-level routines in prime3D_srch
                if( p%use_gpu .eq. 'yes' .or. p%bench_gpu .eq. 'yes' )then
                    select case(p%refine)
                        case('no')
                            if( p%oritab .eq. '' )then
                                call primesrch3D%exec_prime3D_srch(pftcc, iptcl, p%lp, cnt_glob=cnt_glob)
                            else
                                call primesrch3D%exec_prime3D_srch(pftcc, iptcl, p%lp, orientation, cnt_glob=cnt_glob)
                            endif
                        case('neigh')
                            stop 'refine=neigh mode not currently implemented on GPU'
                         case('shc')
                            stop 'refine=shc mode not currently implemented on GPU'
                            if( p%oritab .eq. '' )then
                                call primesrch3D%exec_prime3D_shc_srch(pftcc, iptcl,  p%lp, cnt_glob=cnt_glob)
                            else
                                call primesrch3D%exec_prime3D_shc_srch(pftcc, iptcl,  p%lp, orientation, cnt_glob=cnt_glob)
                            endif
                        case('shcneigh')
                            stop 'refine=shcneigh mode not currently implemented on GPU'
                        case DEFAULT
                            write(*,*) 'The refinement mode: ', trim(p%refine), ' is unsupported on GPU'
                            stop 
                    end select
                else
                   select case(p%refine)
                        case('no')
                            if( p%oritab .eq. '' )then
                                call primesrch3D%exec_prime3D_srch(pftcc, iptcl, p%lp)
                            else
                                call primesrch3D%exec_prime3D_srch(pftcc, iptcl, p%lp, orientation)
                            endif
                        case('neigh')
                            if( p%oritab .eq. '' ) stop 'cannot run the refine=neigh mode without input oridoc (oritab)'
                            call primesrch3D%exec_prime3D_srch(pftcc, iptcl, p%lp, orientation, nnmat=b%nnmat)
                        case('shc')
                            if( p%oritab .eq. '' )then
                                call primesrch3D%exec_prime3D_shc_srch(pftcc, iptcl, p%lp)
                            else
                                call primesrch3D%exec_prime3D_shc_srch(pftcc, iptcl, p%lp, orientation)
                            endif
                        case('shcneigh')
                            if( p%oritab .eq. '' ) stop 'cannot run the refine=shcneigh mode without input oridoc (oritab)'
                            call primesrch3D%exec_prime3D_shc_srch(pftcc, iptcl, p%lp, orientation, nnmat=b%nnmat)
                        case('shift')
                            if( p%oritab .eq. '' ) stop 'cannot run the refine=shift mode without input oridoc (oritab)'
                            call primesrch3D%exec_prime3D_inpl_srch(pftcc, iptcl, p%lp, orientation, greedy=.false.)
                        case('het')
                            if( p%oritab .eq. '' ) stop 'cannot run the refine=het mode without input oridoc (oritab)'
                            ! if(orientation%get('corr') < het_corr_thresh)then
                            !     call primesrch3D%exec_prime3D_het_srch(pftcc, iptcl, orientation, statecnt, do_rnd=.true.)
                            ! else
                            !     call primesrch3D%exec_prime3D_het_srch(pftcc, iptcl, orientation, statecnt, do_rnd=.false.)
                            ! endif
                            call primesrch3D%exec_prime3D_het_srch(pftcc, iptcl, orientation, statecnt, do_rnd=is_rnd(iptcl))
                        case('adasym')
                            if( p%oritab .eq. '' )then
                                call primesrch3D%exec_prime3D_srch(pftcc, iptcl, p%lp)
                            else
                                call primesrch3D%exec_prime3D_srch(pftcc, iptcl, p%lp, orientation)
                            endif
                        case DEFAULT
                            write(*,*) 'The refinement mode: ', trim(p%refine), ' is unsupported on CPU'
                            stop 
                    end select
                endif
                call primesrch3D%get_ori_best(orientation)
                call b%a%set_ori(iptcl,orientation)
                if( p%npeaks>1 )then
                    call primesrch3D%get_oris(prime3D_oris, orientation)
                    if( doshellweight )then
                        wresamp = resample_filter(wmat(iptcl,:), res, res_pad)
                        call grid_ptcl(b, p, iptcl, cnt_glob, orientation, prime3D_oris, shellweights=wresamp)
                    else
                        call grid_ptcl(b, p, iptcl, cnt_glob, orientation, prime3D_oris)
                    endif
                else
                    if( doshellweight )then
                        wresamp = resample_filter(wmat(iptcl,:), res, res_pad)
                        call grid_ptcl(b, p, iptcl, cnt_glob, orientation, shellweights=wresamp)
                    else
                        call grid_ptcl(b, p, iptcl, cnt_glob, orientation)
                    endif
                endif
            else
                call orientation%reject
                call b%a%set_ori(iptcl,orientation)
            endif
        end do
        ! DEV
        if( p%refine.eq.'het')then
            norm = real(sum(statecnt))
            do istate=1,p%nstates
               print *, '% state ', istate, ' is ', 100.*(real(statecnt(istate))/norm)
            end do
        endif
        ! END DEV
        ! orientations output
        call b%a%write(p%outfile, [p%fromp,p%top])
        p%oritab = p%outfile
        call pftcc%kill
        ! NORMALIZE STRUCTURE FACTORS
        if( p%eo .eq. 'yes' )then
            call eonorm_struct_facts(b, p, reslim, which_iter)
        else
            call norm_struct_facts(b, p, which_iter)
        endif
        ! REPORT CONVERGENCE
        if( p%l_distr_exec )then
            call qsys_job_finished( p, 'simple_hadamard3D_matcher :: prime3D_exec')
        else
            ! CONVERGENCE TEST
            if( p%refine .eq. 'het' )then
                converged = b%conv%check_conv_het()
            else
                converged = b%conv%check_conv3D(update_res)
            endif
        endif
        ! DEALLOCATE
        if( allocated(wmat)    ) deallocate(wmat)
        if( allocated(wresamp) ) deallocate(wresamp)
        if( allocated(res)     ) deallocate(res)
        if( allocated(res_pad) ) deallocate(res_pad)
    end subroutine prime3D_exec

    subroutine gen_random_model( b, p, nsamp_in )
        use simple_ran_tabu, only: ran_tabu
        class(build),      intent(inout) :: b
        class(params),     intent(inout) :: p
        integer, optional, intent(in)    :: nsamp_in
        type(ran_tabu)       :: rt
        integer              :: i, j, k, nsamp, alloc_stat
        integer, allocatable :: sample(:)
        if( p%vols(1) == '' )then
            p%oritab = 'prime3D_startdoc.txt'
            if( p%refine .eq. 'no' .or. p%refine.eq.'adasym')then
                call b%a%rnd_oris
                call b%a%zero_shifts
                if( p%l_distr_exec .and. p%part.ne.1 )then
                    ! so random oris only written once in distributed mode
                else
                    call b%a%write( p%oritab )
                endif
            else
                stop 'Unsupported refine option; simple_hadamard3D_matcher::gen_random_model'
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
            do i=1,nsamp
                call progress(i, nsamp)
                orientation = b%a%get_ori(sample(i))
                call b%img%read(p%stk, sample(i), p%l_xfel)
                if( p%l_xfel )then
                    call b%img%pad(b%img_pad)
                else
                    call prep4cgrid(b%img, b%img_pad, p%msk, b%recvols(1)%get_wfuns())
                endif
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
            ! NORMALIZE STRUCTURE FACTORS
            call norm_struct_facts(b, p)
        endif
    end subroutine gen_random_model

    subroutine preppftcc4align( b, p, cline, ppfts_fname )
        class(build),               intent(inout) :: b
        class(params),              intent(inout) :: p
        class(cmdline),             intent(inout) :: cline
        character(len=*), optional, intent(in)    :: ppfts_fname
        integer :: nrefs
        if( .not. p%l_distr_exec ) write(*,'(A)') '>>> BUILDING PRIME3D SEARCH ENGINE'
        ! must be done here since constants in p are dynamically set
        call primesrch3D%new( b%a, b%e, p )
        ! must be done here since p%kfromto is dynamically set based on FSC from previous round
        ! or based on dynamic resolution limit update
        nrefs = p%nspace*p%nstates
        if( p%l_xfel )then
            call pftcc%new(nrefs, [p%fromp,p%top], [p%boxmatch,p%boxmatch,1],&
            p%kfromto, p%ring2, p%ctf, isxfel='yes')
        else
            call pftcc%new(nrefs, [p%fromp,p%top], [p%boxmatch,p%boxmatch,1],&
            p%kfromto, p%ring2, p%ctf)
        endif
        ! PREPARATION OF REFERENCES IN PFTCC
        call prep_refs_pftcc4align( b, p, cline )
        ! PREPARATION OF PARTICLES IN PFTCC
        call prep_ptcls_pftcc4align( b, p, cline, ppfts_fname )
        ! subtract the mean shell values for xfel correlations
        if( p%l_xfel ) call pftcc%xfel_subtract_shell_mean()
        if( debug ) write(*,*) '*** hadamard3D_matcher ***: finished preppftcc4align'
    end subroutine preppftcc4align

    subroutine prep_refs_pftcc4align( b, p, cline )
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        type(ori)                     :: o
        integer :: cnt, s, iref, nrefs
        ! PREPARATION OF REFERENCES IN PFTCC
        ! read reference volumes and create polar projections
        nrefs = p%nspace*p%nstates
        cnt   = 0
        if( .not. p%l_distr_exec ) write(*,'(A)') '>>> BUILDING REFERENCES'
        do s=1,p%nstates
            if( p%oritab .ne. '' )then
                ! greedy start
                if( b%a%get_statepop(s) == 0 )then
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
                call b%proj%fproject_polar(cnt, b%vol, o, p, pftcc)
            end do
        end do
        ! bring back the original b%vol size
        if( p%boxmatch < p%box ) call b%vol%new([p%box,p%box,p%box], p%smpd)
    end subroutine prep_refs_pftcc4align

    subroutine prep_ptcls_pftcc4align( b, p, cline, ppfts_fname )
        class(build),               intent(inout) :: b
        class(params),              intent(inout) :: p
        class(cmdline),             intent(inout) :: cline
        character(len=*), optional, intent(in)    :: ppfts_fname
        type(ori) :: o
        integer   :: cnt, s, iptcl, istate, ntot, progress_cnt
        ! PREPARATION OF PARTICLES IN PFTCC
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
                ! read particle images and create polar projections
                if( .not. p%l_distr_exec ) write(*,'(A)') '>>> BUILDING PARTICLES'
                progress_cnt = 0
                cnt_glob     = 0
                ntot         = p%top-p%fromp+1
                do s=1,p%nstates
                    if( b%a%get_statepop(s) ==0 )then
                        ! empty state
                        cnt_glob = cnt_glob + p%top-p%fromp + 1
                        cycle
                    endif
                    if( p%doautomsk )then
                        ! read & pre-process mask volume
                        call b%mskvol%read(p%masks(s))
                        if( p%eo .eq. 'yes' )then
                            call prep4cgrid(b%mskvol, b%vol_pad, p%msk, b%eorecvols(1)%get_wfuns())
                        else
                            call prep4cgrid(b%mskvol, b%vol_pad, p%msk, b%recvols(1)%get_wfuns())
                        endif
                    endif
                    cnt = 0
                    do iptcl=p%fromp,p%top
                        cnt      = cnt + 1
                        o        = b%a%get_ori(iptcl)
                        istate   = nint(o%get('state'))
                        cnt_glob = cnt_glob+1
                        if( istate /= s ) cycle
                        progress_cnt = progress_cnt + 1
                        call progress( progress_cnt, ntot )
                        if( p%boxmatch < p%box )then
                            ! back to the original size as b%img has been
                            ! and will be clipped/padded in prepimg4align
                            call b%img%new([p%box,p%box,1],p%smpd)
                        endif
                        if( p%l_distr_exec )then
                            call b%img%read(p%stk_part, cnt, p%l_xfel)
                        else
                            call b%img%read(p%stk, iptcl, p%l_xfel)
                        endif
                        call prepimg4align(b, p, o)
                        call b%proj%img2polarft(iptcl, b%img, pftcc)
                    end do
                end do
                ! restores b%img dimensions for clean exit
                if( p%boxmatch < p%box )call b%img%new([p%box,p%box,1],p%smpd)
            end subroutine prep_pftcc_local
        
    end subroutine prep_ptcls_pftcc4align

end module simple_hadamard3D_matcher
