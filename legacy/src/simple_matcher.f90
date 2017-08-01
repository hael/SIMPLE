module simple_matcher
use simple_defs         ! singleton
use simple_polarft_srch ! singleton
use simple_onflyft_srch ! singleton
use simple_cmdline      ! singleton
use simple_ori,         only: ori
use simple_pori,        only: pori
use simple_build,       only: build
use simple_params,      only: params
use simple_jiffys,      only: fopen_err, alloc_err, progress, get_fileunit, del_txtfile
use simple_gridding,    only: prep4cgrid
use simple_masker,      only: automask
use simple_rnd,         only: ran3
use simple_math,        only: deg2rad, rad2deg, file2rarr,&
get_lplim,  ssnr2wiener, get_resolution, locate ! fsc2ssnr, fsc2wiener
implicit none

public :: prime_exec, gen_random_model, prime_find_resrange
private

real                        :: res
logical, parameter          :: report=.false.
type(soft_ori), allocatable :: weights(:)
type(pori)                  :: po
type(ori)                   :: orientation, orisoft, o_sym, o_trial

contains

    !>  \brief  to find the resolution range for initial prime alignment
    subroutine prime_find_resrange( b, p, lp_start, lp_finish )
        use simple_oris, only: oris
        class(build), intent(inout)  :: b
        class(params), intent(inout) :: p
        real, intent(out)            :: lp_start, lp_finish
        real, allocatable            :: peaks(:)
        type(oris)                   :: o
        integer :: lfny, alloc_stat, k, pos10, pos6
        call o%new(p%nspace)
        call o%spiral
        lfny = b%img%get_lfny(1)
        allocate( peaks(lfny), stat=alloc_stat )
        call alloc_err("In: prime_find_resrange, simple_matcher", alloc_stat)
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
    end subroutine
    
    !>  \brief  is the prime algorithm
    subroutine prime_exec( b, p, which_iter, update_res, converged )
        class(build), intent(inout)  :: b
        class(params), intent(inout) :: p        
        integer, intent(in)          :: which_iter
        logical, intent(inout)       :: update_res, converged
        real                         :: corr, wcorr, dist
        real                         :: sdev, frac, w, nmin
        real                         :: nbetter, frac_better, lptmp, mi_hard
        real                         :: fsc0143, fsc05, mapres(p%nstates)
        character(len=STDLEN)        :: dig, dig2
        integer                      :: cnt, i, j, s, k, fnr
        integer                      :: file_stat, loc(1)
        logical                      :: didsrch
        real, allocatable            :: resarr(:)

        if( p%eo .eq. 'yes' )then
            ! we need the worst resolved fsc
            resarr = b%img%get_res()
            do s=1,p%nstates
                write(dig,*) s
                b%fsc(s,:) = file2rarr('fsc_state'//trim(adjustl(dig))//'.bin')
                call get_resolution(b%fsc(s,:), resarr, fsc05, fsc0143)
                mapres(s) = fsc0143
            end do
            if( report ) write(*,*) '*** matcher ***: extracted FSC info'
            loc = maxloc(mapres)
            ! set Fourier index range
            p%kfromto(1) = max(2,b%img%get_find(p%hp))
            if( defined_cmd_arg('lpstop'))then
                p%kfromto(2) = min(get_lplim(b%fsc(loc(1),:)),b%img%get_find(p%lpstop))
            else
                p%kfromto(2) = get_lplim(b%fsc(loc(1),:))
            endif
            p%lp = b%img%get_lp(p%kfromto(2))
            p%lp_dyn = p%lp
            call b%a%set_all('lp',p%lp)
            if( report ) write(*,*) '*** matcher ***: did set Fourier index range'
        else
            ! SET FOURIER INDEX RANGE
            p%kfromto(1) = max(2,b%img%get_find(p%hp))
            if( defined_cmd_arg('lpstop'))then
                p%kfromto(2) = min(b%img%get_find(p%lp),b%img%get_find(p%lpstop))
            else
                p%kfromto(2) = b%img%get_find(p%lp)
            endif
            call b%a%set_all('lp',p%lp)
            if( report ) write(*,*) '*** matcher ***: did set Fourier index range'
        endif
        
        ! CALCULATE ANGULAR TRESHOLD (USED BY THE SPARSE WEIGHTING SCHEME)        
        p%tres = rad2deg(atan(max(p%fny,p%lp)/(p%moldiam/2.)))
        res    = p%lp
        if( report ) write(*,*) '*** matcher ***: calculated angular treshold (used by the sparse weighting scheme)'
        
        ! DETERMINE THE NUMBER OF PEAKS
        if( .not. defined_cmd_arg('npeaks') )then
            p%npeaks = 1
            if( p%refine .eq. 'soft' )then
                p%npeaks = IMPORTANT
            else if( p%refine.eq.'no' .or. p%refine.eq.'het' )then
                if( p%lp >= 10. ) p%npeaks = min(10,b%e%find_npeaks(p%lp, p%moldiam))
            endif
            if( report ) write(*,*) '*** matcher ***: determined the number of peaks'
        endif
        
        ! RANDOM MODEL GENERATION
        if( p%vols(1) .eq. '' )then
            if( p%nptcls > 1000 )then
                call gen_random_model(b, p, 1000)
            else
                call gen_random_model(b, p)
            endif
            if( report ) write(*,*) '*** matcher ***: generated random model'
        endif
        
        ! GENERATE PROJECTIONS (POLAR FTs)
        if( p%oritab .ne. '' .and. p%frac < 0.99 ) call b%a%calc_hard_ptcl_weights(p%frac)
        if( p%refine .eq. 'no' )then
            ! read reference volumes and create polar projections
            do s=1,p%nstates
                call b%refvols(s)%read(p%vols(s))
                if( p%boxmatch < p%box )then
                    call b%refvols(s)%clip_inplace([p%boxmatch,p%boxmatch,p%boxmatch]) ! SQUARE DIMS ASSUMED
                endif
                call b%proj%fprojvol_polar(b%refvols(s), b%e, p, b%refimgs, s)
                call b%refvols(s)%fwd_ft
                b%sigpow(s,:) = b%refvols(s)%spectrum('power')
            end do
            if( report ) write(*,*) '*** matcher ***: read reference volumes and created polar projections'
        else
            ! read reference volumes
            do s=1,p%nstates
                call b%vol%read(p%vols(s))
                b%refvols(s) = b%vol
                if( p%boxmatch < p%box )then
                    call b%refvols(s)%clip_inplace([p%boxmatch,p%boxmatch,p%boxmatch]) ! SQUARE DIMS ASSUMED
                endif
                call b%refvols(s)%fwd_ft
                b%sigpow(s,:) = b%refvols(s)%spectrum('power')
            end do
            if( report )then
                if( p%boxmatch < p%box )then
                    write(*,*) '*** matcher ***: read & clipped reference volumes'
                else
                    write(*,*) '*** matcher ***: read reference volumes'
                endif
            endif
        endif
        
        ! INITIALIZE
        if( p%refine .eq. 'no' )then    
            call polarft_srch_init(p, b) ! initialize polarft_srch object
            if( which_iter <= 0 )then
                write(*,'(A)') '>>> PRIME DISCRETE STOCHASTIC SEARCH'
            else
                write(*,'(A,1X,I3)') '>>> PRIME DISCRETE STOCHASTIC SEARCH, ITERATION:', which_iter
            endif
            write(dig2,*) which_iter
            if( which_iter > 0 ) p%outfile = 'prime_doc'//trim(adjustl(dig2))//'.txt'
        else
            call onflyft_srch_init(b, p) ! initialize onflyft_srch object
            if( report ) write(*,*) '*** matcher ***: initialized onflyft srch'
            if( which_iter <= 0 )then
                if( p%refine .eq. 'het' )then
                    write(*,'(A)') '>>> OASIS HETEROGENEITY ANALYSIS'
                else
                    write(*,'(A)') '>>> OASIS REFINEMENT'
                endif
            else
                if( p%refine .eq. 'het' )then
                    write(*,'(A,1X,I3)') '>>> OASIS HETEROGENEITY ANALYSIS, ITERATION:', which_iter
                else
                    write(*,'(A,1X,I3)') '>>> OASIS REFINEMENT, ITERATION:', which_iter
                endif
            endif
            write(dig2,*) which_iter
            if( which_iter > 0 )then
                if( p%refine .eq. 'het' )then
                    p%outfile = 'oasis_het_doc'//trim(adjustl(dig2))//'.txt'
                else
                    p%outfile = 'oasis_refine_doc'//trim(adjustl(dig2))//'.txt'
                endif
            endif
        endif
        
        ! RESET RECVOLS
        if( p%norec .eq. 'yes' )then
            ! do nothing
        else
            do s=1,p%nstates
                if( p%eo .eq. 'yes' )then
                    call b%eorecvols(s)%reset_all
                else
                    call b%recvols(s)%reset
                endif
            end do
            if( report ) write(*,*) '*** matcher ***: reseted recvols'
        endif

        ! ALIGN & GRID
        call del_txtfile(p%outfile)
        cnt = 0
        if( report ) write(*,*) '*** matcher ***: loop fromp/top:', p%fromp, p%top
        do i=p%fromp,p%top
            cnt = cnt+1
            call progress(cnt, p%top-p%fromp+1)
            if( report ) write(*,*) '*** matcher ***: trying to get orientation:', i
            if( report ) write(*,*) '*** matcher ***: size of oris:', b%a%get_noris()
            orientation = b%a%get_ori(i)
            if( report ) write(*,*) '*** matcher ***: got orientation'     
            o_trial = orientation ! transfer of parameters
            if( report ) write(*,*) '*** matcher ***: transferred ori parameters'
            call prep4align(b, p, i, orientation)
            if( p%refine .eq. 'no' )then
                if( p%oritab .eq. '' )then
                    call polarft_srch_align_ptcl(b%img)
                else
                    call polarft_srch_align_ptcl(b%img, orientation)
                endif
                call polarft_srch_calc_weights(weights, wcorr)
                call polarft_srch_get_ori_best(orientation)
                call orientation%set('corr',wcorr)
                call orientation%set('sdev',polarft_srch_ang_sdev(weights))
                call b%a%set_ori(i,orientation)
            else if( p%refine .eq. 'het' )then
                ! the b%a object is automatically updated
                call onflyft_srch_align(i, po)
            else
                call onflyft_srch_align(i, po)
                if( report ) write(*,*) '*** matcher ***: finished onflyft_srch_align'
            endif
            if( b%pvalid(i) )then
                if( p%refine .ne. 'no' )then                    
                    orientation = b%a%get_ori(i)
                    lptmp = p%lp
                    p%lp  = p%lpvalid
                    call b%a%set(i, 'corr', onflyft_srch_corr_volimg(orientation, nint(orientation%get('state'))))
                    p%lp  = lptmp
                endif
            else
                if( p%norec .eq. 'yes' )then
                    ! no reconstruction
                else
                    call grid_ptcl( b, p, i )
                 endif
            endif
            call b%a%write(i, p%outfile)
            call orientation%kill
        end do
        p%oritab  = p%outfile
        
        ! NORMALIZE STRUCTURE FACTORS
        if( p%norec .eq. 'yes' )then
            ! no normalization needed
        else
            if( p%eo .eq. 'yes' )then
                call eonorm_struct_facts(b, p, which_iter)
            else
                call norm_struct_facts(b, p, which_iter)
            endif
        endif
        
        ! CONVERGENCE TEST
        if( .not. defined_cmd_arg('part') )then
            if( p%refine .eq. 'no' )then
                corr = b%a%get_avg('corr')
                dist = b%a%get_avg('dist')
                frac = b%a%get_avg('frac')
                sdev = b%a%get_avg('sdev')
                write(*,'(A,11X,F5.1)') '>>> ANGLE OF FEASIBLE REGION:', p%tres
                write(*,'(A,2X,F5.1)')  '>>> AVERAGE ANGULAR DISTANCE BTW ORIS:', dist
                write(*,'(A,1X,F5.1)')  '>>> PERCENTAGE OF SEARCH SPACE SCANNED:', frac
                write(*,'(A,25X,F7.4)') '>>> CORRELATION:', corr
                write(*,'(A,11X,F9.2)') '>>> ANGULAR SDEV OF MODEL:', sdev
                ! automatic resolution stepping
                if( update_res )then
                    ! the previous round updated the resolution limit, so
                    ! don't update this round
                    update_res = .false.
                else
                    update_res = .false.
                    if( p%dynlp .eq. 'yes' .and. .not. defined_cmd_arg('lp') .and. dist <= p%tres/2. ) update_res = .true.
                endif    
                converged = .false.
                if( dist < p%tres/5. .and. frac > 99. ) converged = .true.
            else if( p%refine .eq. 'het' )then
                corr = b%a%get_avg('corr')
                mi_hard = b%a%get_avg('mi_hard')
                write(*,'(A,16X,F7.4)') '>>> DISTRIBUTION OVERLAP:', mi_hard
                write(*,'(A,25X,F7.4)') '>>> CORRELATION:', corr
                if( p%fracvalid > 0. )then
                    corr = b%a%get_avg('corr',mask=b%pvalid)
                    write(*,'(A,20X,F7.4)') '>>> FREE CORRELATION:', corr
                endif
                if( mi_hard > 0.99 ) converged = .true.
            else
                corr = b%a%get_avg('corr')
                dist = b%a%get_nonzero_avg('dist')
                nbetter = b%a%get_sum('nbetter')
                frac_better = 100.*(nbetter/real(cnt))
                sdev = b%a%get_avg('sdev')
                nmin = b%a%get_avg('nmin')
                write(*,'(A,11X,F5.1)') '>>> ANGLE OF FEASIBLE REGION:', p%tres
                write(*,'(A,2X,F5.1)')  '>>> AVERAGE ANGULAR DISTANCE BTW ORIS:', dist
                write(*,'(A,4X,F5.1)')  '>>> FOUND BETTER ORIS FOR THIS FRAC:', frac_better
                write(*,'(A,25X,F7.4)') '>>> CORRELATION:', corr
                write(*,'(A,11X,F9.2)') '>>> ANGULAR SDEV OF MODEL:', sdev
                write(*,'(A,11X,F9.2)') '>>> AVERAGE NR OF MINIMAS:', nmin
                if( p%fracvalid > 0. )then
                    corr = b%a%get_avg('corr',mask=b%pvalid)
                    write(*,'(A,20X,F7.4)') '>>> EARLY STOPPING CORRELATION:', corr
                endif
                ! check 4 convergence
                converged = .false.
                if( frac_better <= 15. )then
                    converged = .true.
                else if( dist < p%tres/5. .or. res <= p%fny )then
                    converged = .true.
                endif
            endif
        else
            if( p%refine .eq. 'no' )then
                corr = b%a%get_avg('corr',fromto=[p%fromp,p%top])
                dist = b%a%get_avg('dist',fromto=[p%fromp,p%top])
                frac = b%a%get_avg('frac',fromto=[p%fromp,p%top])
                sdev = b%a%get_avg('sdev',fromto=[p%fromp,p%top])
                fnr  = get_fileunit()
                write(dig,*) p%part
                open(unit=fnr, FILE='JOB_FINISHED_'//trim(adjustl(dig)), STATUS='REPLACE', action='WRITE', iostat=file_stat)
                call fopen_err( 'In: prime_ini; simple_matcher.f90', file_stat )
                write(fnr,'(A,11X,F5.1)') '>>> ANGLE OF FEASIBLE REGION:', p%tres
                write(fnr,'(A,2X,F5.1)')  '>>> AVERAGE ANGULAR DISTANCE BTW ORIS:', dist
                write(fnr,'(A,1X,F5.1)')  '>>> PERCENTAGE OF SEARCH SPACE SCANNED:', frac
                write(fnr,'(A,25X,F7.4)') '>>> CORRELATION:', corr
                write(fnr,'(A,11X,F9.2)') '>>> ANGULAR SDEV OF MODEL:', sdev
                ! automatic resolution stepping
                update_res = .false.
                if( p%dynlp .eq. 'yes' .and. .not. defined_cmd_arg('lp') .and. dist <= p%tres/2. )then 
                    update_res = .true.
                    write(fnr,'(A)') '>>> UPDATE LOW-PASS LIMIT: .YES.'
                else
                    write(fnr,'(A)') '>>> UPDATE LOW-PASS LIMIT: .NO.'
                endif
                converged = .false.
                if( dist < p%tres/5. .and. frac > 99. )then
                    converged = .true.
                    write(fnr,'(A)') '>>> CONVERGED: .YES.'
                else
                    write(fnr,'(A)') '>>> CONVERGED: .NO.'
                endif
                close(fnr)
            else if( p%refine .eq. 'het' )then
                corr    = b%a%get_avg('corr',fromto=[p%fromp,p%top])
                mi_hard = b%a%get_avg('mi_hard',fromto=[p%fromp,p%top])
                fnr     = get_fileunit()
                write(dig,*) p%part
                open(unit=fnr, FILE='JOB_FINISHED_'//trim(adjustl(dig)), STATUS='REPLACE', action='WRITE', iostat=file_stat)
                call fopen_err( 'In: statematch; simple_matcher.f90', file_stat )
                write(fnr,'(A,16X,F7.4)') '>>> DISTRIBUTION OVERLAP:', mi_hard
                write(fnr,'(A,25X,F7.4)') '>>> CORRELATION:', corr
                if( p%fracvalid > 0. )then
                    corr = b%a%get_avg('corr',fromto=[p%fromp,p%top],mask=b%pvalid)
                    write(fnr,'(A,20X,F7.4)') '>>> FREE CORRELATION:', corr
                endif
            else
                dist = b%a%get_nonzero_avg('dist',fromto=[p%fromp,p%top])
                nbetter = b%a%get_sum('nbetter',fromto=[p%fromp,p%top])
                frac_better = 100.*(nbetter/real(cnt))
                sdev = b%a%get_avg('sdev',fromto=[p%fromp,p%top])
                nmin = b%a%get_avg('nmin',fromto=[p%fromp,p%top])
                fnr  = get_fileunit()
                write(dig,*) p%part
                open(unit=fnr, FILE='JOB_FINISHED_'//trim(adjustl(dig)), STATUS='REPLACE', action='WRITE', iostat=file_stat)
                call fopen_err( 'In: prime_refine; simple_matcher.f90', file_stat )
                write(fnr,'(A,11X,F5.1)') '>>> ANGLE OF FEASIBLE REGION:', p%tres
                write(fnr,'(A,2X,F5.1)')  '>>> AVERAGE ANGULAR DISTANCE BTW ORIS:', dist
                write(fnr,'(A,4X,F5.1)')  '>>> FOUND BETTER ORIS FOR THIS FRAC:', frac_better
                write(fnr,'(A,25X,F7.4)') '>>> CORRELATION:', corr
                write(fnr,'(A,11X,F9.2)') '>>> ANGULAR SDEV OF MODEL:', sdev
                write(fnr,'(A,11X,F9.2)') '>>> AVERAGE NR OF MINIMAS:', nmin
                if( p%fracvalid > 0. )then
                    corr = b%a%get_avg('corr',fromto=[p%fromp,p%top],mask=b%pvalid)
                    write(fnr,'(A,20X,F7.4)') '>>> EARLY STOPPING CORRELATION:', corr
                endif
                converged = .false.
                if( frac_better <= 15. )then
                    converged = .true.
                    write(fnr,'(A)') '>>> CONVERGED: .YES.'
                else if( dist < p%tres/5. .or. res <= p%fny )then
                    converged = .true.
                    write(fnr,'(A)') '>>> CONVERGED: .YES.'
                else
                    write(fnr,'(A)') '>>> CONVERGED: .NO.'
                endif
                close(fnr)
            endif
        endif
    end subroutine
    
    !>  \brief  grids one particle image to the volume
    subroutine grid_ptcl( b, p, i )
        class(build), intent(inout)  :: b
        class(params), intent(inout) :: p        
        integer, intent(in)          :: i
        real      :: pw, ran, w
        integer   :: j, s, k, endit
        pw = orientation%get('w')
        if( pw > 0. )then
            ! prepare image for gridding
            ! using the uncorrected/unmodified image as input
            if( p%eo .eq. 'yes' )then
                call prep4cgrid(b%img_copy, b%img_pad, p%msk, b%eorecvols(1)%get_wfuns())
            else
                call prep4cgrid(b%img_copy, b%img_pad, p%msk, b%recvols(1)%get_wfuns())
            endif
            if( report ) write(*,*) '*** matcher ***: prepared image for gridding'
            ran = ran3()
            if( p%refine .eq. 'soft' .or. p%refine .eq. 'het' )then
                endit = po%get_nbetter()
            else
                endit = p%npeaks
            endif
            orisoft = orientation
            do j=1,endit
                if( report ) write(*,*) '*** matcher ***: gridding, iteration:', j
                ! get ori info
                if( p%refine .eq. 'no' )then
                    call polarft_srch_get_ori(weights(j), orisoft)
                    w = weights(j)%w
                    s = weights(j)%state
                else if( p%refine .eq. 'soft' .or. p%refine .eq. 'het' )then
                    call po%minimum2ori(j,orisoft)
                    w = orisoft%get('ow')
                    s = nint(orisoft%get('state'))
                else
                    orisoft = b%a%get_ori(i)
                    w = 1.
                    s = nint(orisoft%get('state'))
                endif
                if( report ) write(*,*) '*** matcher ***: got orientation'
                if( p%frac < 0.99 ) w = w*pw
                if( w > 0. )then
                    ! grid
                    if( p%pgrp == 'c1' )then
                        if( p%eo .eq. 'yes' )then
                            call b%eorecvols(s)%grid_fplane(orisoft, p%doshift, b%img_pad, pwght=w, ran=ran)
                        else
                            call b%recvols(s)%inout_fplane(orisoft, .true., p%doshift, b%img_pad, pwght=w)
                        endif
                    else
                        do k=1,b%se%get_nsym()
                            o_sym = b%se%apply(orisoft, k)
                            if( p%eo .eq. 'yes' )then
                                call b%eorecvols(s)%grid_fplane(o_sym, p%doshift, b%img_pad, pwght=w, ran=ran)
                            else
                                call b%recvols(s)%inout_fplane(o_sym, .true., p%doshift, b%img_pad, pwght=w)
                            endif
                        end do
                    endif
                endif
                if( report ) write(*,*) '*** matcher ***: gridded ptcl'
            end do
            call orisoft%kill
        endif
    end subroutine
     
    subroutine gen_random_model( b, p, nsamp_in )
        use simple_ran_tabu, only: ran_tabu
        class(build), intent(inout)   :: b
        class(params), intent(inout)  :: p
        integer, intent(in), optional :: nsamp_in
        type(ran_tabu)                :: rt
        integer                       :: i, j, s, k, nsamp, alloc_stat, en
        character(len=STDLEN)         :: dig
        integer, allocatable          :: sample(:)
        real                          :: euls(3), ran
        if( p%vols(1) == '' )then
            p%oritab = 'rndoris.txt'
            if( p%refine .eq. 'no')then
                call b%a%rnd_oris(p%eullims, 0.)
                call b%a%zero_shifts
                if( p%npeaks > 1 ) call b%a%rnd_weights(p%nspace, p%npeaks, p%nstates)
                call b%a%rnd_states(p%nstates)
                call b%a%write(p%oritab)
            else if( p%refine .eq. 'het' )then
                if( .not. defined_cmd_arg('lp'))then
                    stop 'Need low-pass limit 4 state randomization; gen_random_model; simple_matcher'
                endif
                if( .not. defined_cmd_arg('oritab'))then
                    stop 'Need input oris 4 state randomization; gen_random_model; simple_matcher'
                endif
                if( p%nstates < 2 )then
                    stop 'Too few states 4 state randomization (<2); gen_random_model; simple_matcher'
                endif
                p%tres = rad2deg(atan(max(p%fny,p%lp)/(p%moldiam/2.)))
                p%npeaks = 1
                call b%a%rnd_states(p%nstates)
                call b%a%write(p%oritab)
                call b%a%introd_alig_err(2.*p%tres, p%trs)
            else
                stop 'Unknown refine option; gen_random_model; simple_matcher'
            endif
            if( p%noise .eq. 'yes' )then
                do s=1,p%nstates
                    write(dig,*) s
                    p%vols(s) = 'startvol'//'_state'//trim(adjustl(dig))//p%ext
                    call b%vol%ran
                    if( defined_cmd_arg('inner') )then
                        call b%vol%mask(p%msk,'soft',inner=p%inner,width=p%width)
                    else
                        call b%vol%mask(p%msk,'soft')
                    endif
                    call b%vol%write(p%vols(s), del_if_exists=.true.)
                end do
            else
                if( present(nsamp_in) )then
                    nsamp = nsamp_in
                else
                    nsamp = p%nptcls
                endif
                allocate( sample(nsamp), stat=alloc_stat )
                call alloc_err("In: gen_random_model; simple_matcher", alloc_stat)
                if( present(nsamp_in) )then
                    rt = ran_tabu(p%nptcls)
                    call rt%ne_ran_iarr(sample)
                    call rt%kill
                else
                    do i=1,nsamp
                        sample(i) = i
                    end do
                endif
                if( p%nstates > 1 )then
                    write(*,'(A)') '>>> RECONSTRUCTING RANDOM MODELS'
                else
                    write(*,'(A)') '>>> RECONSTRUCTING RANDOM MODEL'
                endif
                do i=1,nsamp
                    call progress(i, nsamp)
                    orisoft = b%a%get_ori(sample(i))
                    if( p%refine .eq. 'no' .and. p%npeaks > 1 )then
                        weights = orisoft%get_weights()
                        en = size(weights)
                    else
                        if( allocated(weights) ) deallocate(weights)
                        allocate( weights(1) )
                        weights(1)%w = 1.
                        en = 1
                        weights(1)%state = nint(orisoft%get('state'))
                    endif
                    call b%img%read(p%stk, sample(i))
                    if( p%eo .eq. 'yes' )then
                        call prep4cgrid(b%img, b%img_pad, p%msk, b%eorecvols(1)%get_wfuns())
                    else
                        call prep4cgrid(b%img, b%img_pad, p%msk, b%recvols(1)%get_wfuns())
                    endif
                    ran = ran3()
                    do j=1,en
                        if( p%refine .eq. 'no')then
                            euls = b%e%get_euler(weights(j)%proj)
                            call orisoft%set_euler(euls)
                        endif
                        if( p%pgrp == 'c1' )then
                            if( p%eo .eq. 'yes' )then      
                                call b%eorecvols(weights(j)%state)%grid_fplane(orisoft,&
                                p%doshift, b%img_pad, pwght=weights(j)%w, ran=ran)
                            else
                                call b%recvols(weights(j)%state)%inout_fplane(orisoft,&
                                .true., p%doshift, b%img_pad, pwght=weights(j)%w)
                            endif
                        else
                            do k=1,b%se%get_nsym()
                                o_sym = b%se%apply(orisoft, k)
                                if( p%eo .eq. 'yes' )then
                                    call b%eorecvols(s)%grid_fplane(o_sym, p%doshift, b%img_pad,&
                                    pwght=weights(j)%w, ran=ran)
                                else
                                    call b%recvols(s)%inout_fplane(o_sym, .true., p%doshift,&
                                    b%img_pad, pwght=weights(j)%w)
                                endif
                            end do
                        endif
                    end do
                    call orisoft%kill
                    deallocate(weights)
                end do
                deallocate( sample )
            endif            
            ! NORMALIZE STRUCTURE FACTORS
            if( p%eo .eq. 'yes' )then
                call eonorm_struct_facts(b, p)
            else
                call norm_struct_facts(b, p)
            endif
        endif
    end subroutine
    
    subroutine prep4align( b, p, i, o )
        class(build), intent(inout)  :: b
        class(params), intent(inout) :: p
        integer, intent(in)          :: i
        class(ori), intent(inout)    :: o
        real                         :: x, y
        integer                      :: s, j
        ! create new image if needed (if the image 4 alignment is auto-clipped later; see below)
        if( p%boxmatch < p%box ) call b%img%new([p%box,p%box,1],p%smpd)
        call b%img%read(p%stk, i)
        if( report ) write(*,*) '*** matcher ***: read image'
        ! save a copy for later
        b%img_copy = b%img
        ! parse ori
        x = o%get('x')                            
        y = o%get('y')
        s = nint(o%get('state'))
        if( p%ctf .eq. 'yes' )then ! multiply ptcl image with CTF and make CTF**2 image
            call b%img%fwd_ft
            if( p%ctfmode .eq. 'astig' )then ! astigmatic CTF
                if( p%ctf .eq. 'flip' )then
                    call b%tfun%apply(b%img, b%a%get(i,'dfx'), 'abs', b%a%get(i,'dfy'), b%a%get(i,'angast'))
                else if( p%ctf .eq. 'yes' )then
                    call b%tfun%apply(b%img, b%a%get(i,'dfx'), 'ctf', b%a%get(i,'dfy'), b%a%get(i,'angast'))
                endif
            else if( p%ctfmode .eq. 'noastig' )then ! non-astigmatic CTF
                if( p%ctf .eq. 'flip' )then
                    call b%tfun%apply(b%img, b%a%get(i,'dfx'), 'abs')
                else if( p%ctf .eq. 'yes' )then
                    call b%tfun%apply(b%img, b%a%get(i,'dfx'), 'ctf')
                endif
            endif
            ! shift image to rotational origin
            call b%img%shift(-x, -y)
            ! apply a soft-edged mask
            call b%img%bwd_ft
            ! clip image if needed
            if( p%boxmatch < p%box ) call b%img%clip_inplace([p%boxmatch,p%boxmatch,1]) ! SQUARE DIMS ASSUMED
            call b%img%mask(p%msk, 'soft')
        else if( p%ctf .eq. 'refine' )then
            ! do nothing
        else ! no CTF
            call b%img%fwd_ft
            call b%img%shift(-x, -y)
            if( report ) write(*,*) '*** matcher ***: shifted image to previously found origin'
            call b%img%bwd_ft
            ! clip image if needed
            if( p%boxmatch < p%box )then
                call b%img%clip_inplace([p%boxmatch,p%boxmatch,1]) ! SQUARE DIMS ASSUMED
                if( report ) write(*,*) '*** matcher ***: clipped inplace'
            endif
            call b%img%mask(p%msk, 'soft') ! apply a soft-edged mask
        endif
        if( report ) write(*,*) '*** matcher ***: finished prep4align'
    end subroutine
    
    subroutine eonorm_struct_facts( b, p, which_iter )
        class(build), intent(inout)   :: b
        class(params), intent(inout)  :: p
        integer, intent(in), optional :: which_iter
        integer                       :: s
        character(len=STDLEN)         :: dig, dig2, dig3
        real                          :: res05s(p%nstates), res0143s(p%nstates)
        do s=1,p%nstates
            ! create new vol if needed (if the refvol 4 alignment was auto-clipped later; see above)
            if( p%boxmatch < p%box ) call b%refvols(s)%new([p%box,p%box,p%box],p%smpd)
            write(dig,*) s
            if( defined_cmd_arg('part') )then
                write(dig3,*) p%part
                call b%eorecvols(s)%write_eos('recvol'//'_state'//trim(adjustl(dig))//'_part'//trim(adjustl(dig3)))
            else
                if( present(which_iter) )then
                    if( which_iter <= 0 )then
                        p%vols(s)     = 'recvol'//'_state'//trim(adjustl(dig))//p%ext
                        p%vols_msk(s) = 'recvol'//'_state'//trim(adjustl(dig))//'msk'//p%ext
                        p%masks(s)    = 'automask'//'_state'//trim(adjustl(dig))//p%ext
                    else
                        write(dig2,*) which_iter
                        p%vols(s)     = 'recvol'//'_state'//trim(adjustl(dig))//'_iter'//trim(adjustl(dig2))//p%ext
                        p%vols_msk(s) = 'recvol'//'_state'//trim(adjustl(dig))//'_iter'//trim(adjustl(dig2))//'msk'//p%ext
                        p%masks(s)    = 'automask'//'_state'//trim(adjustl(dig))//'_iter'//trim(adjustl(dig2))//p%ext
                    endif
                else
                     p%vols(s)     = 'startvol'//'_state'//trim(adjustl(dig))//p%ext
                     p%vols_msk(s) = 'startvol'//'_state'//trim(adjustl(dig))//'msk'//p%ext
                     p%masks(s)    = 'startmask'//'_state'//trim(adjustl(dig))//p%ext
                endif
                call b%eorecvols(s)%sum_eos
                call b%eorecvols(s)%sampl_dens_correct_eos(p%msk,s)
                call b%eorecvols(s)%get_res(res05s(s), res0143s(s))
                call b%eorecvols(s)%sampl_dens_correct_sum(b%refvols(s))
                if( defined_cmd_arg('inner') )then
                    call b%refvols(s)%mask(p%msk, 'soft', inner=p%inner, width=p%width)
                else
                    call b%refvols(s)%mask(p%msk, 'soft')
                endif
                call b%refvols(s)%write(p%vols(s), del_if_exists=.true.)
            endif
        end do
        if( .not. defined_cmd_arg('part') )then
            ! set the resolution limit according to the worst resolved model
            res      = maxval(res0143s)
            p%lp     = min(p%lp,max(p%lpstop,res))
            if( .not. defined_cmd_arg('amsklp') )then
                p%amsklp = max(b%img%get_lp(b%img%get_find(res)-2),maxval(res05s))
            endif
            do s=1,p%nstates                
                if( (defined_cmd_arg('mw') .or. defined_cmd_arg('nvox'))&
                .or. defined_cmd_arg('mskfile') )then
                    call automask(b, p, b%refvols(s), b%mskvol, p%vols_msk(s), p%masks(s))
                    p%vols(s) = p%vols_msk(s)
                endif
            end do
        endif
    end subroutine

    subroutine norm_struct_facts( b, p, which_iter )
        class(build), intent(inout)   :: b
        class(params), intent(inout)  :: p
        integer, intent(in), optional :: which_iter
        integer                       :: s
        character(len=STDLEN)         :: dig, dig2, dig3, fbody
        do s=1,p%nstates
            ! create new vol if needed (if the refvol 4 alignment was auto-clipped later; see above)
            if( p%boxmatch < p%box ) call b%refvols(s)%new([p%box,p%box,p%box],p%smpd)
            write(dig,*) s
            if( defined_cmd_arg('part') )then
                write(dig3,*) p%part
                fbody = 'recvol_state'//trim(adjustl(dig))//'_part'//trim(adjustl(dig3))
                p%vols(s)  = trim(adjustl(fbody))//p%ext
                p%masks(s) = 'rho_'//trim(adjustl(fbody))//p%ext
                call b%recvols(s)%write(p%vols(s), del_if_exists=.true.)
                call b%recvols(s)%write_rho(p%masks(s))
            else
                if( present(which_iter) )then
                    if( which_iter <= 0 )then
                        p%vols(s)     = 'recvol'//'_state'//trim(adjustl(dig))//p%ext
                        p%vols_msk(s) = 'recvol'//'_state'//trim(adjustl(dig))//'msk'//p%ext
                        p%masks(s)    = 'automask'//'_state'//trim(adjustl(dig))//p%ext
                    else
                        write(dig2,*) which_iter
                        p%vols(s)     = 'recvol'//'_state'//trim(adjustl(dig))//'_iter'//trim(adjustl(dig2))//p%ext
                        p%vols_msk(s) = 'recvol'//'_state'//trim(adjustl(dig))//'_iter'//trim(adjustl(dig2))//'msk'//p%ext
                        p%masks(s)    = 'automask'//'_state'//trim(adjustl(dig))//'_iter'//trim(adjustl(dig2))//p%ext
                    endif
                else
                     p%vols(s)     = 'startvol'//'_state'//trim(adjustl(dig))//p%ext
                     p%vols_msk(s) = 'startvol'//'_state'//trim(adjustl(dig))//'msk'//p%ext
                     p%masks(s)    = 'startmask'//'_state'//trim(adjustl(dig))//p%ext
                endif
                call b%recvols(s)%sampl_dens_correct(self_out=b%vol_pad) ! this preserves the recvol for online update
                call b%vol_pad%bwd_ft
                call b%vol_pad%clip(b%refvols(s))
                call b%refvols(s)%bp(0.,p%lp)
                if( defined_cmd_arg('inner') )then
                    call b%refvols(s)%mask(p%msk, 'soft', inner=p%inner, width=p%width)
                else
                    call b%refvols(s)%mask(p%msk, 'soft')
                endif
                call b%refvols(s)%write(p%vols(s), del_if_exists=.true.)
                if( (defined_cmd_arg('mw') .or. defined_cmd_arg('nvox'))&
                .or. defined_cmd_arg('mskfile') )then
                    call automask(b, p, b%refvols(s), b%mskvol, p%vols_msk(s), p%masks(s))
                    p%vols(s) = p%vols_msk(s)
                endif
            endif
        end do
    end subroutine
    
end module simple_matcher
