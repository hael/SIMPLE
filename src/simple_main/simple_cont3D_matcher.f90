module simple_cont3D_matcher
use simple_defs
use simple_cartft_corrcalc, only: cartft_corrcalc
use simple_ori,             only: ori
use simple_build,           only: build
use simple_params,          only: params
use simple_masker,          only: automask
use simple_cmdline,         only: cmdline
use simple_qsys_funs,       only: qsys_job_finished
use simple_strings,         only: int2str_pad
use simple_hadamard_common  ! us all in there
use simple_math             ! us all in there
use simple_cftcc_srch       ! us all in there
implicit none

public :: cont3D_exec, cont3D_shellweight
private

type(cartft_corrcalc) :: cftcc
real                  :: res
integer               :: cnt_glob=0
integer, parameter    :: MAXNPEAKS=10, NRESTARTS=3
logical, parameter    :: debug=.false.
character(len=STDLEN) :: OPT_STR='simplex'

contains

    subroutine cont3D_shellweight( b, p, cline )
        use simple_ori,  only: ori
        use simple_stat, only: corrs2weights, normalize_sigm
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        real, allocatable             :: res(:), corrs(:), weights_tmp(:), wmat(:,:), wsums(:)
        character(len=:), allocatable :: fname
        type(ori)                     :: orientation
        integer                       :: iptcl, state, filtsz, alloc_stat
        integer                       :: ishell, io_stat, filnum, fnr, file_stat
        filtsz = b%img%get_filtsz()
        allocate( wmat(p%top-p%fromp+1,filtsz), corrs(filtsz), stat=alloc_stat)
        call alloc_err('In: simple_cont3D_matcher :: cont3D_shellweight', alloc_stat)
        wmat  = 0.
        corrs = 0.
        call cftcc%new( b, p, cline )
        cnt_glob = 0
        write(*,'(A)') '>>> CALCULATING FOURIER RING CORRELATIONS'
        do iptcl=p%fromp,p%top
            cnt_glob = cnt_glob + 1
            call progress(cnt_glob, p%top-p%fromp+1)
            orientation = b%a%get_ori(iptcl)
            if( nint(orientation%get('state')) > 0 )then
                if( p%boxmatch < p%box )call b%img%new([p%box,p%box,1],p%smpd) ! ensures correct dimensions
                if( p%l_distr_exec )then
                    call b%img%read(p%stk_part, cnt_glob, p%l_xfel)
                else
                    call b%img%read(p%stk, iptcl, p%l_xfel)
                endif
                call prepimg4align(b, p, orientation)
                call cftcc%frc(orientation, 1, b%img, res, corrs)
                wmat(cnt_glob,:) = corrs(:)
            else
                ! weights should be zero
                if( .not. allocated(corrs) ) allocate(corrs(filtsz))
                corrs = 0.
                wmat(cnt_glob,:) = corrs
            endif
        end do
        ! write the weight matrix
        filnum = get_fileunit()
        if( p%l_distr_exec )then
            allocate(fname, source='shellweights_part'//int2str_pad(p%part,p%numlen)//'.bin')
        else
            allocate(fname, source='shellweights.bin')
        endif
        open(unit=filnum, status='REPLACE', action='WRITE', file=fname, access='STREAM')
        write(unit=filnum,pos=1,iostat=io_stat) wmat
        ! check if the write was successful
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(cont3D_shellweight): I/O error ',&
            io_stat, ' when writing to cont3D_shellweights_partX.bin'
            stop 'I/O error; cont3D_shellweight_part*; simple_cont3D_matcher'
        endif
        close(filnum)
        deallocate(wmat, fname)
        call qsys_job_finished( p, 'simple_cont3D_matcher :: cont3D_shellweight' )
    end subroutine cont3D_shellweight
    
    !>  \brief  is the continuous refinement algorithm
    subroutine cont3D_exec( b, p, cline, which_iter, converged )
        use simple_ori,       only: ori
        use simple_filterer,  only: resample_filter
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        integer, intent(in)           :: which_iter
        logical, intent(inout)        :: converged
        type(ori)         :: orientation
        real, allocatable :: wmat(:,:), wresamp(:), res(:), res_pad(:)
        real              :: corr, dist, frac_srch_space, reslim
        integer           :: state, file_stat, fnr, iptcl, s
        logical           :: doshellweight
        
        ! CREATE THE CARTESIAN CORRELATOR
        call cftcc%new( b, p, cline )
        if( p%boxmatch < p%box )call b%vol%new([p%box,p%box,p%box],p%smpd) ! ensures correct dimensions
        
        ! SET BAND-PASS LIMIT RANGE 
        call set_bp_range( b, p, cline )
        
        ! CALCULATE ANGULAR THRESHOLD
        p%athres = rad2deg(atan(max(p%fny,p%lp)/(p%moldiam/2.)))
        reslim   = p%lp

        ! SET FRACTION OF SEARCH SPACE
        frac_srch_space = b%a%get_avg('frac')
        
        ! GENERATE PARTICLE WEIGHTS
        if( p%oritab .ne. '' .and. p%frac < 0.99 ) call b%a%calc_hard_ptcl_weights(p%frac)
        if( p%l_distr_exec )then
            ! nothing to do
        else
            if( p%l_shellw .and. frac_srch_space >= 50. )then
                call cont3D_shellweight(b, p, cline)
            endif
        endif
        call setup_shellweights(b, p, doshellweight, wmat, res, res_pad)

        ! INITIALIZE
        call cftcc_srch_init(cftcc, b%img, OPT_STR, p%optlims(:5,:), NRESTARTS)
        
        if( which_iter <= 0 )then
            write(*,'(A)') '>>> CONTINUOUS ORIENTATION SEARCH'
        else
            write(*,'(A,1X,I3)') '>>> CONTINUOUS ORIENTATION SEARCH, ITERATION:', which_iter
        endif
        if( which_iter > 0 ) p%outfile = 'cont3Ddoc_'//int2str_pad(which_iter,3)//'.txt'
        
        ! RESET RECVOLS
        do s=1,p%nstates
            if( p%eo .eq. 'yes' )then
                call b%eorecvols(s)%reset_all
            else
                call b%recvols(s)%reset
            endif
        end do
        if( debug ) write(*,*) '*** cont3D_matcher ***: did reset recvols'
        
        ! ALIGN & GRID
        call del_file(p%outfile)
        cnt_glob = 0
        if( debug ) write(*,*) '*** cont3D_matcher ***: loop fromp/top:', p%fromp, p%top
        do iptcl=p%fromp,p%top
            cnt_glob = cnt_glob+1
            call progress(cnt_glob, p%top-p%fromp+1)
            orientation = b%a%get_ori(iptcl)
            state = nint(orientation%get('state'))
            if( s > 0 )then
                if( p%boxmatch < p%box )call b%img%new([p%box,p%box,1],p%smpd) ! ensures correct dimensions
                if( p%l_distr_exec )then
                    call b%img%read(p%stk_part, cnt_glob, p%l_xfel)
                else
                    call b%img%read(p%stk, iptcl, p%l_xfel)
                endif
                call prepimg4align(b, p, orientation)
                call cftcc_srch_set_state(state)
                call cftcc_srch_minimize(orientation)
                call b%a%set_ori(iptcl,orientation)
                if( doshellweight )then
                    wresamp = resample_filter(wmat(iptcl,:), res, res_pad)
                    call grid_ptcl(b, p, iptcl, cnt_glob, orientation, shellweights=wresamp)
                else
                    call grid_ptcl(b, p, iptcl, cnt_glob, orientation)
                endif
            else
                call orientation%reject
                call b%a%set_ori(iptcl,orientation)
            endif
            call b%a%write(iptcl, p%outfile)
        end do
        p%oritab = p%outfile
        
        ! NORMALIZE STRUCTURE FACTORS
        if( p%eo .eq. 'yes' )then
            call eonorm_struct_facts(b, p, reslim, which_iter)
        else
            call norm_struct_facts(b, p, which_iter)
        endif

        if( p%l_distr_exec )then
            call qsys_job_finished( p, 'simple_hadamard3D_matcher.f90 :: prime3D_exec' )
        else
            ! CONVERGENCE TEST
            converged = b%conv%check_conv3D()
        endif

        ! DEALLOCATE
        if( allocated(wmat)    ) deallocate(wmat)
        if( allocated(wresamp) ) deallocate(wresamp)
        if( allocated(res)     ) deallocate(res)
        if( allocated(res_pad) ) deallocate(res_pad)
    end subroutine cont3D_exec

end module simple_cont3D_matcher