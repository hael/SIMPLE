module simple_cont3D_matcher
use simple_cartft_corrcalc, only: cartft_corrcalc
use simple_ori,             only: ori
use simple_build,           only: build
use simple_params,          only: params
use simple_masker,          only: automask
use simple_cmdline,         only: cmdline
use simple_hadamard_common  ! singleton
use simple_defs             ! singleton
use simple_jiffys           ! singleton
use simple_math             ! singleton
use simple_cftcc_srch       ! singleton
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
        use simple_stat, only: corrs2weights, moment, normalize_sigm
        use gnufor2,     only: plot
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        real, allocatable             :: res(:), corrs(:), weights_tmp(:), wmat(:,:), wsums(:)
        character(len=:), allocatable :: fname
        type(ori)                     :: orientation
        integer                       :: iptcl, state, filtsz, recsz, alloc_stat
        integer                       :: ishell, io_stat, filnum, i, fnr, file_stat
        real                          :: shellstats(4), var
        logical                       :: err, doprint=.false.
        if( p%l_distr_exec )then
            filtsz = b%img%get_filtsz()
            allocate( wmat(p%top-p%fromp+1,filtsz), corrs(filtsz), stat=alloc_stat)
            call alloc_err('In: simple_cont3D_matcher :: cont3D_shellweight, 1', alloc_stat)
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
            open(unit=filnum, status='REPLACE', action='WRITE',&
            file='cont3D_shellweights_part'//int2str_pad(p%part,p%numlen)//'.bin', access='STREAM')
            write(unit=filnum,pos=1,iostat=io_stat) wmat
            ! check if the write was successful
            if( io_stat .ne. 0 )then
                write(*,'(a,i0,2a)') '**ERROR(cont3D_shellweight): I/O error ',&
                io_stat, ' when writing to cont3D_shellweights_partX.bin'
                stop 'I/O error; cont3D_shellweight_part*; simple_cont3D_matcher'
            endif
            close(filnum)
            ! generation of this file marks completion of the partition
            ! this file is empty 4 now but may contain run stats etc.
            fnr  = get_fileunit()
            open(unit=fnr, FILE='JOB_FINISHED_'//int2str_pad(p%part,p%numlen),&
            STATUS='REPLACE', action='WRITE', iostat=file_stat)
            call fopen_err( 'In: cont3D_shellweight; simple_cont3D_matcher.f90', file_stat )
            close(fnr)
            deallocate(wmat)
        else
            write(*,*) 'simple_cont3D_matcher :: cont3D_shellweight'
            stop 'not implemented for non-distributed mode'
        endif
    end subroutine cont3D_shellweight
    
    !>  \brief  is the prime3D algorithm
    subroutine cont3D_exec( b, p, cline, which_iter, converged )
        use simple_ori, only: ori
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        integer, intent(in)           :: which_iter
        logical, intent(inout)        :: converged
        type(ori) :: orientation
        integer   :: state, file_stat, fnr, iptcl, s
        real      :: corr, dist
        
        ! CREATE THE CARTESIAN CORRELATOR
        call cftcc%new( b, p, cline )
        
        ! SET BAND-PASS LIMIT RANGE 
        call set_bp_range( b, p, cline )
        
        ! CALCULATE ANGULAR THRESHOLD
        p%athres = rad2deg(atan(max(p%fny,p%lp)/(p%moldiam/2.)))
        
        ! GENERATE PARTICLE WEIGHTS
        if( p%oritab .ne. '' .and. p%frac < 0.99 ) call b%a%calc_hard_ptcl_weights(p%frac)
        
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
        call del_txtfile(p%outfile)
        cnt_glob = 0
        if( debug ) write(*,*) '*** cont3D_matcher ***: loop fromp/top:', p%fromp, p%top
        do iptcl=p%fromp,p%top
            cnt_glob = cnt_glob+1
            call progress(cnt_glob, p%top-p%fromp+1)
            orientation = b%a%get_ori(iptcl)
            if( nint(orientation%get('state')) > 0 )then
                if( p%l_distr_exec )then
                    call b%img%read(p%stk_part, cnt_glob, p%l_xfel)
                else
                    call b%img%read(p%stk, iptcl, p%l_xfel)
                endif
                call prepimg4align(b, p, orientation)
                state = nint(orientation%get('state'))
                call cftcc_srch_set_state(state)
                call cftcc_srch_minimize(orientation)
                call b%a%set_ori(iptcl,orientation)
                call grid_ptcl( b, p, iptcl, cnt_glob, orientation)
            else
                call orientation%reject
                call b%a%set_ori(iptcl,orientation)
            endif
            call b%a%write(iptcl, p%outfile)
        end do
        p%oritab = p%outfile
        
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
            converged = b%conv%check_conv3D()
        endif

    end subroutine cont3D_exec

end module simple_cont3D_matcher