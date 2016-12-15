module simple_hadamard3D_matcher
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_prime3D_srch,     only: prime3D_srch
use simple_ori,              only: ori
use simple_build,            only: build
use simple_params,           only: params
use simple_gridding,         only: prep4cgrid
use simple_masker,           only: automask
use simple_rnd,              only: ran3
use simple_hadamard_common   ! singleton
use simple_defs              ! singleton
use simple_cuda_defs         ! singleton
use simple_jiffys,           ! singleton
use simple_cmdline           ! singleton
use simple_math              ! singleton
use simple_timing            ! singleton
implicit none

public :: prime3D_exec, gen_random_model, prime3D_find_resrange, preppftcc4align
private

type(polarft_corrcalc) :: pftcc
type(prime3D_srch)     :: primesrch3D
real                   :: res
type(ori)              :: orientation, orisoft, o_sym, o_tmp
integer                :: cnt_glob=0
integer, parameter     :: MAXNPEAKS=10
logical, parameter     :: debug=.false.
!timer string
character(len=10) :: char_out
character(len=80) :: tmr_name
!CUDA err variable for the return function calls
integer      :: err
!integer,external      :: strlen
!function calls
!integer,external      :: get_length_of_string_c
!integer,external      :: convert_int2char_pos_c
integer      :: length

interface external_c_function_timer
   function get_length_of_string_c(int_in)
     integer :: int_in
     integer :: get_length_of_string_c
   end function get_length_of_string_c
   function convert_int2char_pos_c(char_out, int_in)
     integer :: int_in
     character(len=3) :: char_out
     integer :: convert_int2char_pos_c
   end function convert_int2char_pos_c
   function convert_int2char_indexed_c(char_out, iter, iint, max_iter)
     integer :: iter,iint,max_iter
     character(len=3) :: char_out
     integer :: convert_int2char_indexed_c
   end function convert_int2char_indexed_c
   function strlen(string)
     character(len=*) string
     integer :: strlen
   end function strlen
end interface external_c_function_timer

contains

    !>  \brief  to find the resolution range for initial prime3D alignment
    subroutine prime3D_find_resrange( b, p, lp_start, lp_finish )
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
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        integer, intent(in)           :: which_iter
        logical, intent(inout)        :: update_res, converged
        real                          :: wcorr, w, lptmp
        integer                       :: iptcl, fnr, file_stat, s
        double precision              :: elps_r
        double precision,dimension(2) :: st_r, et_r
        !local variables
        integer                       :: rc
        !fucntion call
        integer                       :: setdevice_c

        ! SET BAND-PASS LIMIT RANGE 
        call set_bp_range( b, p, cline )

        ! CALCULATE ANGULAR THRESHOLD (USED BY THE SPARSE WEIGHTING SCHEME)
        p%athres = rad2deg(atan(max(p%fny,p%lp)/(p%moldiam/2.)))
        res      = p%lp
        if( debug ) write(*,*) '*** hadamard3D_matcher ***: calculated angular threshold (used by the sparse weighting scheme)'

        ! DETERMINE THE NUMBER OF PEAKS
        if( .not. cline%defined('npeaks') )then
            select case(p%refine)
                case('no', 'qcont', 'qcontneigh')
                    p%npeaks = min(MAXNPEAKS,b%e%find_npeaks(p%lp, p%moldiam))
                case DEFAULT
                    p%npeaks = 1
            end select
            if( debug ) write(*,*) '*** hadamard3D_matcher ***: determined the number of peaks'
        endif

        ! RANDOM MODEL GENERATION
        if( p%vols(1) .eq. '' )then
            if( p%nptcls > 1000 )then
                call gen_random_model(b, p, 1000)
            else
                call gen_random_model(b, p)
            endif
            if( debug ) write(*,*) '*** hadamard3D_matcher ***: generated random model'
        endif

        ! GENERATE PROJECTIONS (POLAR FTs)
        if( p%oritab .ne. '' .and. p%frac < 0.99 ) call b%a%calc_hard_ptcl_weights(p%frac)

        call preppftcc4align( b, p, cline )

        ! INITIALIZE
        if( which_iter <= 0 )then
            if( str_has_substr(p%refine,'qcont') )then
                write(*,'(A)') '>>> PRIME3D QUASI-CONTINUOUS STOCHASTIC SEARCH'
            else
                write(*,'(A)') '>>> PRIME3D DISCRETE STOCHASTIC SEARCH'
            endif
        else
            if( str_has_substr(p%refine,'qcont') )then
                write(*,'(A,1X,I3)') '>>> PRIME3D QUASI-CONTINUOUS STOCHASTIC SEARCH, ITERATION:', which_iter
            else
                write(*,'(A,1X,I3)') '>>> PRIME3D DISCRETE STOCHASTIC SEARCH, ITERATION:', which_iter
            endif
        endif
        if( which_iter > 0 ) p%outfile = 'prime3Ddoc_'//int2str_pad(which_iter,3)//'.txt'

        length = get_length_of_string_c(which_iter)
        !err = convert_int2char_pos_c(char_out,which_iter)
        err = convert_int2char_indexed_c(char_out,which_iter,1,p%maxits)
        tmr_name = 'y_primesrch3D_iter'
        tmr_name = tmr_name(1:strlen(tmr_name))//char_out
        tmr_name = tmr_name(1:strlen(tmr_name))
!         write(*,*) "tmr_name: ",tmr_name
!         write(*,*) "tmr_name: ",tmr_name
!         write(*,*) "p%use_gpu: ",p%use_gpu
!         write(*,*) "p%bench_gpu: ",p%bench_gpu
!         write(*,*) "p%fix_gpu: ",p%fix_gpu
!         write(*,*) "p%set_gpu: ",p%set_gpu
!         write(*,*) "has_gpu: ",has_gpu
!         write(*,*) "use_gpu: ",use_gpu
!         write(*,*) "has_multi_gpu: ",has_multi_gpu
#if defined (CUDA)
        rc = setdevice_c(p%set_gpu)
#endif
        if( p%use_gpu .eq. 'yes' .or. p%bench_gpu .eq. 'yes' )then
           if (b%s_bench%bench_i==1 .or. ibench .eqv. .true. ) call start_timer_cpu(tmr_name)
           if (b%s_bench%bench_i==1 .or. ibench .eqv. .true. ) call gettimeofday_c(st_r)
           ! CALCULATE ALL CORRELATIONS BEFORE THE SEARCH BEGINS
!!!!           
           !TODO: need to fix the logix between the use_gpu and bench_gpu
!!!!
           if( p%bench_gpu .eq. 'yes' .and. p%use_gpu .eq. 'no')then
              call primesrch3D%calc_corrs_on_gpu(pftcc, b%a, .true.)
           else
              call primesrch3D%calc_corrs_on_gpu(pftcc, b%a)
           endif
           if (b%s_bench%bench_i==1 .or. ibench .eqv. .true. ) then
              call stop_timer_cpu(tmr_name)
              call gettimeofday_c(et_r)
              call elapsed_time_c(st_r,et_r,elps_r)
              write(*,'(a,x,f15.6,x,a)') "Elapsed time for primesrch3D: ",elps_r," secs"
              write(100+which_iter,'(f7.4,x,f8.3,x,i5,x,i5,x,i5,x,i3)') &
                   p%smpd,p%lp,p%ring2,p%nspace,round2even(twopi*real(p%ring2)),p%kfromto(2)-p%kfromto(1)+1
              write(100+which_iter,'(a,x,f12.4,x,a)') "Elapsed time for primesrch3D: ",elps_r," secs"
           end if
           
           if (b%s_bench%bench_i==1 .or. ibench .eqv. .true. ) call elapsed_time_c(st_r,et_r,elps_r)
           !write(1500,*) elps_r
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
        call del_txtfile(p%outfile)
        cnt_glob = 0
        if( debug ) write(*,*) '*** hadamard3D_matcher ***: loop fromp/top:', p%fromp, p%top
        if( p%use_gpu .eq. 'no' .and. p%bench_gpu .eq. 'no' )then
           if (b%s_bench%bench_i==1 .or. ibench .eqv. .true. ) call start_timer_cpu(tmr_name)
           if (b%s_bench%bench_i==1 .or. ibench .eqv. .true. ) call gettimeofday_c(st_r)
        end if
        do iptcl=p%fromp,p%top
            cnt_glob = cnt_glob+1
            call progress(cnt_glob, p%top-p%fromp+1)
            orientation = b%a%get_ori(iptcl)
            if( nint(orientation%get('state')) > 0 )then
                if( .not. str_has_substr(p%refine,'qcont') ) call preprefs4align(b, p, iptcl, pftcc)
                ! execute the high-level routines in prime3D_srch
                if(  p%use_gpu .eq. 'yes' .or. p%bench_gpu .eq. 'yes' )then
                    select case(p%refine)
                        case('no')
                            if( p%oritab .eq. '' )then
                                call primesrch3D%exec_prime3D_srch(pftcc, iptcl, p%lp, wcorr, cnt_glob=cnt_glob)
                            else
                                call primesrch3D%exec_prime3D_srch(pftcc, iptcl, p%lp, wcorr, orientation, cnt_glob=cnt_glob)
                            endif
                        case('neigh')
                            stop 'refine=neigh mode not currently implemented on GPU'
                         case('shc')
                            stop 'refine=shc mode not currently implemented on GPU'
                            if( p%oritab .eq. '' )then
                                call primesrch3D%exec_prime3D_shc_srch(pftcc, iptcl,  p%lp, wcorr, cnt_glob=cnt_glob)
                            else
                                call primesrch3D%exec_prime3D_shc_srch(pftcc, iptcl,  p%lp, wcorr, orientation, cnt_glob=cnt_glob)
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
                                call primesrch3D%exec_prime3D_srch(pftcc, iptcl, p%lp, wcorr)
                            else
                                call primesrch3D%exec_prime3D_srch(pftcc, iptcl, p%lp, wcorr, orientation)
                            endif
                        case('neigh')
                            if( p%oritab .eq. '' ) stop 'cannot run the refine=neigh mode without input oridoc (oritab)'
                            call primesrch3D%exec_prime3D_srch(pftcc, iptcl, p%lp, wcorr, orientation, nnmat=b%nnmat)
                        case('shc')
                            if( p%oritab .eq. '' )then
                                call primesrch3D%exec_prime3D_shc_srch(pftcc, iptcl,  p%lp, wcorr)
                            else
                                call primesrch3D%exec_prime3D_shc_srch(pftcc, iptcl,  p%lp, wcorr, orientation)
                            endif
                        case('shcneigh')
                            if( p%oritab .eq. '' ) stop 'cannot run the refine=shcneigh mode without input oridoc (oritab)'
                            call primesrch3D%exec_prime3D_shc_srch(pftcc,&
                            iptcl, p%lp, wcorr, orientation, nnmat=b%nnmat)
                        case('qcont')
                            if( p%oritab .eq. '' )then
                                if( p%oritab .eq. '' ) stop 'cannot run the refine=qcont mode without input oridoc (oritab)'
                            else
                                call primesrch3D%exec_prime3D_qcont_srch(b%refvols,&
                                b%proj, b%tfun, pftcc, iptcl, p%lp, wcorr, orientation)
                            endif
                        case('qcontneigh')
                            if( p%oritab .eq. '' )then
                                if( p%oritab .eq. '' ) stop 'cannot run the refine=qcontneigh mode without input oridoc (oritab)'
                            else
                                call primesrch3D%exec_prime3D_qcont_srch(b%refvols, b%proj, b%tfun,&
                                pftcc, iptcl, p%lp, wcorr, orientation, athres=max(10.,2.*p%athres))
                            endif
                        case DEFAULT
                            write(*,*) 'The refinement mode: ', trim(p%refine), ' is unsupported on CPU'
                            stop 
                    end select
                endif
                call primesrch3D%get_ori_best(orientation)
                call orientation%set('corr', wcorr)
                call b%a%set_ori(iptcl,orientation)
                call grid_ptcl( b, p, iptcl, cnt_glob, orientation, primesrch3D )
            else
                call orientation%reject
                call b%a%set_ori(iptcl,orientation)
            endif
            call b%a%write(iptcl, p%outfile)
            call orientation%kill
        end do
        if( p%use_gpu .eq. 'no' .and. p%bench_gpu .eq. 'no' )then
           if (b%s_bench%bench_i==1 .or. ibench .eqv. .true. ) then
              call stop_timer_cpu(tmr_name)
              call gettimeofday_c(et_r)
              call elapsed_time_c(st_r,et_r,elps_r)
              write(*,'(a,x,f15.6,x,a)') "Elapsed time for primesrch3D: ",elps_r," secs"
              write(100+which_iter,'(f7.4,x,f8.3,x,i5,x,i5,x,i5,x,i3)') &
                   p%smpd,p%lp,p%ring2,p%nspace,round2even(twopi*real(p%ring2)),p%kfromto(2)-p%kfromto(1)+1
              write(100+which_iter,'(a,x,f12.4,x,a)') "Elapsed time for primesrch3D: ",elps_r," secs"
           end if
        end if
        p%oritab = p%outfile
        call pftcc%kill

        ! NORMALIZE STRUCTURE FACTORS
        if( p%eo .eq. 'yes' )then
            call eonorm_struct_facts(b, p, res, which_iter)
        else
            call norm_struct_facts(b, p, which_iter)
        endif

        if( p%l_distr_exec )then
            ! generation of this file marks completion of the partition
            ! this file is empty 4 now but may contain run stats etc.
            fnr  = get_fileunit()
            open(unit=fnr, FILE='JOB_FINISHED_'//int2str_pad(p%part,p%numlen), STATUS='REPLACE', action='WRITE', iostat=file_stat)
            call fopen_err( 'In: prime3D_exec; simple_hadamard3D_matcher.f90', file_stat )
            close(fnr)
        else
            ! CONVERGENCE TEST
            converged = b%conv%check_conv3D(update_res)
        endif
        
    end subroutine prime3D_exec

    subroutine gen_random_model( b, p, nsamp_in )
        use simple_ran_tabu, only: ran_tabu
        class(build), intent(inout)   :: b
        class(params), intent(inout)  :: p
        integer, intent(in), optional :: nsamp_in
        type(ran_tabu)                :: rt
        integer                       :: i, j, k, nsamp, alloc_stat, en
        integer, allocatable          :: sample(:)
        real                          :: euls(3), ran
        type(soft_ori), allocatable   :: soris(:)
        if( p%vols(1) == '' )then
            p%oritab = 'prime3D_startdoc.txt'
            if( p%refine .eq. 'no')then
                call b%a%rnd_oris
                call b%a%zero_shifts
                if( p%npeaks > 1 ) call b%a%rnd_weights(p%nspace, p%npeaks, p%nstates)
                call b%a%write(p%oritab)
            else
                stop 'Unsupported refine option; simple_hadamard3D_matcher::gen_random_model'
            endif
            p%vols(1) = 'startvol'//p%ext
            if( p%noise .eq. 'yes' )then
                call b%vol%ran
                call b%vol%write(p%vols(1), del_if_exists=.true.)
            else
                if( present(nsamp_in) )then
                    nsamp = nsamp_in
                else
                    nsamp = p%nptcls
                endif
                allocate( sample(nsamp), stat=alloc_stat )
                call alloc_err("In: gen_random_model; simple_hadamard3D_matcher", alloc_stat)
                if( present(nsamp_in) )then
                    rt = ran_tabu(p%nptcls)
                    call rt%ne_ran_iarr(sample)
                    call rt%kill
                else
                    do i=1,nsamp
                        sample(i) = i
                    end do
                endif
                write(*,'(A)') '>>> RECONSTRUCTING RANDOM MODEL'
                do i=1,nsamp
                    call progress(i, nsamp)
                    orisoft = b%a%get_ori(sample(i))
                    if( p%refine .eq. 'no' .and. p%npeaks > 1 )then
                        soris = orisoft%get_weights()
                        en = size(soris)
                    else
                        if( allocated(soris) ) deallocate(soris)
                        allocate( soris(1) )
                        soris(1)%w = 1.
                        en = 1
                        soris(1)%state = 1
                    endif
                    call b%img%read(p%stk, sample(i), p%l_xfel)
                    if( p%l_xfel )then
                        call b%img%pad(b%img_pad)
                    else
                        call prep4cgrid(b%img, b%img_pad, p%msk, b%recvols(1)%get_wfuns())
                    endif
                    ran = ran3()
                    do j=1,en
                        if( p%refine .eq. 'no' .and. p%npeaks > 1 )then
                            euls = b%e%get_euler(soris(j)%proj)
                            call orisoft%set_euler(euls)
                        endif
                        if( p%pgrp == 'c1' )then
                            call b%recvols(1)%inout_fplane(orisoft, .true., b%img_pad, pwght=soris(j)%w)
                        else
                            do k=1,b%se%get_nsym()
                                o_sym = b%se%apply(orisoft, k)
                                call b%recvols(1)%inout_fplane(o_sym, .true., b%img_pad, pwght=soris(j)%w)
                            end do
                        endif
                    end do
                    call orisoft%kill
                    deallocate(soris)
                end do
                deallocate(sample)
            endif
            ! NORMALIZE STRUCTURE FACTORS
            if( p%noise .ne. 'yes' ) call norm_struct_facts(b, p)
        endif
    end subroutine gen_random_model

    subroutine preppftcc4align( b, p, cline )
        use simple_image, only: image
        class(build), intent(inout)   :: b
        class(params), intent(inout)  :: p
        class(cmdline), intent(inout) :: cline
        type(image)                   :: vol4grid
        type(ori)                     :: o
        real, allocatable             :: sqrtssnr(:)
        integer :: cnt, s, iref, iptcl, istate, nrefs, ntot
        write(*,'(A)') '>>> BUILDING PRIME3D SEARCH ENGINE'
        ! must be done here since constants in p are dynamically set
        call primesrch3D%new( b%a, b%e, p )
        write(*,'(A)') '>>> BUILDING DATA STRUCTURE FOR POLARFT CORRELATION CALCULATION'
        ! must be done here since p%kfromto is dynamically set based on FSC from previous round
        ! or based on dynamic resolution limit update
        nrefs = p%nspace*p%nstates
        if( str_has_substr(p%refine,'qcont') ) nrefs = 1 ! on-line update of the reference section
        if( p%l_xfel )then
            call pftcc%new(nrefs, [p%fromp,p%top], [p%boxmatch,p%boxmatch,1],&
            p%kfromto, p%ring2, p%ctf, isxfel='yes')
        else
            call pftcc%new(nrefs, [p%fromp,p%top], [p%boxmatch,p%boxmatch,1],&
            p%kfromto, p%ring2, p%ctf)
        endif
        ! PREPARATION OF REFERENCES IN PFTCC
        ! read reference volumes and create polar projections
        cnt = 0
        write(*,'(A)') '>>> BUILDING REFERENCES'
        do s=1,p%nstates
            ! read & clip
            if( p%boxmatch < p%box ) call b%vol%new([p%box,p%box,p%box],p%smpd) ! ensure correct dim
            call b%vol%read(p%vols(s), isxfel=p%l_xfel)
            if( p%boxmatch < p%box )then
                call b%vol%clip_inplace([p%boxmatch,p%boxmatch,p%boxmatch]) ! SQUARE DIMS ASSUMED
            endif
            ! take care of masking
            if( p%l_xfel )then
                ! no masking
            else
                ! mask volume using a spherical soft-edged mask
                p%vols_msk(s) = add2fbody(p%vols(s), p%ext, 'msk')
                if( p%l_innermsk )then
                    call b%vol%mask(p%msk, 'soft', inner=p%inner, width=p%width)
                else
                    call b%vol%mask(p%msk, 'soft')
                endif
                ! mask using a molecular envelope
                if( p%doautomsk )then
                    p%masks(s) = 'automask_state'//int2str_pad(s,2)//p%ext
                    if( p%l_distr_exec )then
                        if( p%part == 1 )then
                            ! automask & write files
                            call automask(b, p, cline, b%vol, b%mskvol, p%vols_msk(s), p%masks(s))
                        else
                            ! automask & DO NOT write files
                            call automask(b, p, cline, b%vol, b%mskvol)
                        endif
                    else
                        ! automask & write files
                        call automask(b, p, cline, b%vol, b%mskvol, p%vols_msk(s), p%masks(s))
                    endif
                endif
            endif
            ! FT volume
            call b%vol%fwd_ft
!             if( p%eo .eq. 'yes' )then
!                 ! shell normalise
!                 call b%vol%shellnorm
!                 ! scale amplitudes according to PSSNR
!                 allocate(sqrtssnr(size(b%pssnr3D(s,:))))
!                 where(b%pssnr3D(s,:) > 0.)
!                     sqrtssnr = sqrt(b%pssnr3D(s,:))
!                 else where
!                     sqrtssnr = 0.
!                 end where
!                 call b%vol%apply_filter(sqrtssnr)
!                 deallocate(sqrtssnr)
!             else
!                 ! scale amplitudes with a low-pass filter
!                 call b%vol%bp(0.,p%lp)
!             endif
            if( str_has_substr(p%refine,'qcont') )then
                ! store reference volume
                b%refvols(s) = b%vol
            else 
                ! generate discrete projections
                do iref=1,p%nspace
                    cnt = cnt+1
                    call progress(cnt, p%nstates*p%nspace)
                    o = b%e%get_ori(iref)
                    call b%proj%fproject_polar(cnt, b%vol, o, p, pftcc)
                end do
            endif
        end do
        ! bring back the original b%vol size
        if( p%boxmatch < p%box ) call b%vol%new([p%box,p%box,p%box], p%smpd)
        ! PREPARATION OF PARTICLES IN PFTCC
        ! read particle images and create polar projections
        write(*,'(A)') '>>> BUILDING PARTICLES'
        cnt_glob = 0
        ntot     = (p%top-p%fromp+1)*p%nstates
        do s=1,p%nstates
            if( p%doautomsk )then
                ! read & pre-process mask volume
                if( p%boxmatch < p%box ) call b%mskvol%new([p%box,p%box,p%box], p%smpd)
                call b%mskvol%read(p%masks(s))
                if( p%eo .eq. 'yes' )then
                    call prep4cgrid(b%mskvol, b%vol_pad, p%msk, b%eorecvols(1)%get_wfuns())
                else
                    call prep4cgrid(b%mskvol, b%vol_pad, p%msk, b%recvols(1)%get_wfuns())
                endif
            endif
            cnt = 0
            do iptcl=p%fromp,p%top
                cnt    = cnt+1
                o      = b%a%get_ori(iptcl)
                istate = nint(o%get('state'))
                cnt_glob = cnt_glob+1
                call progress(cnt_glob, p%nstates*(p%top-p%fromp+1))
                if( istate /= s )cycle
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
                if( p%eo .eq. 'yes' )then
                    call prepimg4align(b, p, iptcl, pssnr=b%pssnr2D(istate,:))
                else
                    call prepimg4align(b, p, iptcl)
                endif
                call b%proj%img2polarft(iptcl, b%img, pftcc)
            end do
        end do
        ! restores b%img dimensions for clean exit
        if( p%boxmatch < p%box )call b%img%new([p%box,p%box,1],p%smpd)
        if( debug ) write(*,*) '*** hadamard3D_matcher ***: finished preppftcc4align'
    end subroutine preppftcc4align

end module simple_hadamard3D_matcher
