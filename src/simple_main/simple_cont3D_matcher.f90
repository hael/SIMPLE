module simple_cont3D_matcher
use simple_defs
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_ori,              only: ori
use simple_oris,             only: oris
use simple_build,            only: build
use simple_params,           only: params
!use simple_masker,           only: automask
use simple_cmdline,          only: cmdline
use simple_qsys_funs,        only: qsys_job_finished
use simple_strings,          only: int2str_pad
use simple_hadamard_common  ! use all in there
use simple_math             ! use all in there
use simple_pcont3D_srch,     only: pcont3D_srch

implicit none

public :: cont3D_exec
private

type(polarft_corrcalc)  :: pftcc
integer, parameter      :: NRESTARTS = 5
character(len=STDLEN)   :: OPT_STR   = 'simplex'

type(oris)             :: orefs                   !< per particle projection direction search space
logical, allocatable   :: state_exists(:)
integer                :: nptcls          = 0
integer                :: nrefs_per_ptcl  = 0
integer                :: neff_states     = 0

integer, parameter     :: BATCHSZ_MUL = 10   ! particles per thread
integer, parameter     :: MAXNPEAKS   = 10
integer, parameter     :: NREFS       = 50
logical, parameter     :: debug = .false.

contains
    
    subroutine cont3D_exec(b, p, cline, which_iter, converged)
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        integer,        intent(inout) :: which_iter
        logical,        intent(inout) :: converged
        real    :: reslim
        integer :: state
        ! SET BAND-PASS LIMIT RANGE 
        call set_bp_range( b, p, cline )
        reslim = p%lp

        ! RESET RECVOLS
        if(p%norec .eq. 'no')then
            do state=1,p%nstates
                if( p%eo .eq. 'yes' )then
                    call b%eorecvols(state)%reset_all
                else
                    call b%recvols(state)%reset
                endif
            end do
            if( debug ) write(*,*) '*** cont3D_matcher ***: did reset recvols'
        endif

        ! GENERATE PARTICLE WEIGHTS
        if( p%nptcls <= SPECWMINPOP )then
            call b%a%calc_hard_ptcl_weights(p%frac)
        else
            call b%a%calc_spectral_weights(p%frac)
        endif

        ! POLAR/CARTESIAN FORK
        select case(p%refine)
            case('greedy')
                call cont3D_polar_exec(b, p, cline, which_iter, converged)
            case DEFAULT
                stop 'Unknown refinment mode; simple_cont3D_matcher%cont3D_exec'
        end select

        ! NORMALIZE STRUCTURE FACTORS
        if(p%norec .eq. 'no')then
            if( p%eo .eq. 'yes' )then
                call eonorm_struct_facts(b, p, reslim, which_iter)
            else
                call norm_struct_facts(b, p, which_iter)
            endif
        endif

        ! FINISHING
        if( p%l_distr_exec )then
            call qsys_job_finished( p, 'simple_cont3D_matcher.f90 :: cont3D_exec' )
        else
            ! CONVERGENCE TEST
            converged = b%conv%check_conv3D()
        endif
    end subroutine cont3D_exec

    !>  \brief  is the continuous refinement algorithm
    subroutine cont3D_polar_exec( b, p, cline, which_iter, converged )
        use simple_pftcc_contsrch    ! use all in there
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: which_iter
        logical,        intent(inout) :: converged
        type(ori)              :: orientation
        integer                :: state, iptcl, cnt_glob

        ! CREATE CALCULATOR
        if( p%l_xfel )then
            call pftcc%new(1, [1,1], [p%boxmatch,p%boxmatch,1], p%kfromto, p%ring2, p%ctf, isxfel='yes')
        else
            call pftcc%new(1, [1,1], [p%boxmatch,p%boxmatch,1], p%kfromto, p%ring2, p%ctf)
        endif

        ! INITIALIZE
        call pftcc_contsrch_init(b, p, cline, pftcc, b%img, OPT_STR, NRESTARTS)
        if(p%boxmatch < p%box)call b%vol%new([p%box,p%box,p%box],p%smpd) ! ensures correct dimensions
        write(*,'(A)',advance='no') '>>> POLAR-FT '
        if( which_iter <= 0 )then
            write(*,'(A)')'CONTINUOUS ORIENTATION SEARCH'
        else
            write(*,'(A,1X,I3)')'CONTINUOUS ORIENTATION SEARCH, ITERATION:', which_iter
            if( which_iter > 0 ) p%outfile = 'cont3Ddoc_'//int2str_pad(which_iter,3)//'.txt'
        endif
        
        ! ALIGN & GRID
        call del_file(p%outfile)
        cnt_glob = 0
        if( debug ) write(*,*) '*** cont3D_matcher ***: loop fromp/top:', p%fromp, p%top
        do iptcl=p%fromp,p%top
            cnt_glob = cnt_glob+1
            call progress(cnt_glob, p%top-p%fromp+1)
            orientation = b%a%get_ori(iptcl)
            state = nint(orientation%get('state'))
            if(state > 0)then
                call read_img_from_stk(b, p, iptcl)
                b%img_copy = b%img ! stash raw image for rec
                ! ALIGN
                call prepimg4align(b, p, orientation)
                call pftcc_contsrch_set_state(state)
                call pftcc_contsrch_minimize(orientation)
                ! GRID
                if(p%norec.eq.'no')then
                    b%img = b%img_copy
                    call grid_ptcl(b, p, orientation)
                endif
            else
                call orientation%reject
            endif
            call b%a%set_ori(iptcl, orientation)
            call b%a%write(iptcl, p%outfile)
        end do
        p%oritab = p%outfile
        
        ! cleanup
        if( p%boxmatch < p%box )call b%img%new([p%boxmatch,p%boxmatch,1],p%smpd) ! for next iteration in local execution
        call b%refvols(1)%kill_expanded
        call pftcc_contsrch_reset
        call pftcc%kill
    end subroutine cont3D_polar_exec

    !>  \brief  is the prime3D algorithm
    subroutine pcont3D_exec( b, p, cline, which_iter, converged )
        use simple_map_reduce, only: split_nobjs_even
        use simple_qsys_funs,  only: qsys_job_finished
        use simple_strings,    only: int2str_pad
        use simple_projector,  only: projector
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: which_iter
        logical,        intent(inout) :: converged
        ! batches-related variables
        type(projector),        allocatable :: batch_imgs(:)
        type(polarft_corrcalc), allocatable :: pftccs(:)
        !type(pcont3D_srch),     allocatable :: pcont3Dsrchs(:)
        integer,                allocatable :: batches(:,:)
        ! other variables
        type(ori)  :: orientation
        real       :: reslim, frac_srch_space
        integer    :: nbatches, batch, fromp, top, iptcl, state, alloc_stat
        logical    :: update_res
        ! AUTOMASKING DEACTIVATED FOR NOW
        ! INIT
        nptcls          = p%top - p%fromp + 1       ! number of particles processed
        allocate(state_exists(p%nstates))
        state_exists    = b%a%get_state_exist(p%nstates)    ! state existence
        neff_states     = count(state_exists)       ! number of non-empty states
        nrefs_per_ptcl  = NREFS*neff_states         ! number of references per particle
        frac_srch_space = b%a%get_avg('frac')       ! Fraction of the search space
        ! batches
        nbatches = ceiling(real(nptcls)/real(p%nthr*BATCHSZ_MUL)) ! congruence with number of threads
        batches  = split_nobjs_even(nptcls, nbatches)

        ! SET BAND-PASS LIMIT RANGE
        call set_bp_range( b, p, cline )
        reslim = p%lp

        ! CALCULATE ANGULAR THRESHOLD (USED BY THE SPARSE WEIGHTING SCHEME)
        p%athres = rad2deg(atan(max(p%fny,p%lp)/(p%moldiam/2.)))
        write(*,'(A,F6.2)')'>>> ANGULAR THRESHOLD: ', p%athres

        ! DETERMINE THE NUMBER OF PEAKS
        if( .not. cline%defined('npeaks') )then
            p%npeaks = min(MAXNPEAKS,b%e%find_npeaks(p%lp, p%moldiam))
        endif
        write(*,'(A,I2)')'>>> NPEAKS: ', p%npeaks

        ! SETUP WEIGHTS FOR THE 3D RECONSTRUCTION
        if( p%nptcls <= SPECWMINPOP )then
            call b%a%calc_hard_ptcl_weights(p%frac)
        else
            call b%a%calc_spectral_weights(p%frac)
        endif

        ! PREPARE REFVOLS
        call prep_vols(b, p, cline)

        ! RESET RECVOLS
        if(p%norec .eq. 'no')then
            do state=1,p%nstates
                if( state_exists(state) )then
                    if( p%eo .eq. 'yes' )then
                        call b%eorecvols(state)%reset_all
                    else
                        call b%recvols(state)%reset
                    endif
                endif
            end do
            if(debug)write(*,*)'*** pcont3D_matcher ***: did reset recvols'
        endif

        ! INIT IMGPOLARIZER
        ! dummy pftcc is only init here so the img polarizer can be initialized
        call pftcc%new(nrefs_per_ptcl, [1,1], [p%boxmatch,p%boxmatch,1],p%kfromto, p%ring2, p%ctf)
        call b%img%init_imgpolarizer(pftcc)
        call pftcc%kill

        ! INITIALIZE
        if( which_iter <= 0 )then
            write(*,'(A)')'>>> CONTINUOUS POLAR-FT ORIENTATION SEARCH'
        else
            write(*,'(A,1X,I3)')'>>> CONTINUOUS POLAR-FT ORIENTATION SEARCH, ITERATION:', which_iter
            p%outfile = 'cont3Ddoc_'//int2str_pad(which_iter,3)//'.txt'
        endif

        ! BATCH PROCESSING
        call del_file(p%outfile)      
        do batch = 1, nbatches
            ! BATCH INDICES
            fromp = p%fromp-1 + batches(batch,1)
            top   = p%fromp-1 + batches(batch,2)
            ! PREP BATCH
            allocate(pftccs(fromp:top),& !pcont3Dsrchs(fromp:top),&
                &batch_imgs(fromp:top), stat=alloc_stat)
            call alloc_err('In pcont3D_matcher::pcont3D_exec_single',alloc_stat)
            do iptcl = fromp, top
                state = nint(b%a%get(iptcl, 'state'))
                if(state == 0)cycle
                ! stash raw image for rec
                call read_img_from_stk(b, p, iptcl)
                batch_imgs(iptcl) = b%img
                ! prep pftccs & ctf
                call prep_pftcc(b, p, iptcl, pftccs(iptcl))
                if( p%ctf.ne.'no' )call pftccs(iptcl)%create_polar_ctfmats(p%smpd, b%a)
                !call pcont3Dsrchs(iptcl)%new(p, b%a, orefs, pftccs(iptcl), iptcl)
            enddo
            ! SERIAL SEARCHES
            !$omp parallel do default(shared) schedule(guided) private(iptcl)
            do iptcl = fromp, top
                !call pcont3Dsrchs(iptcl)%do_srch(b%a)
            enddo
            !$omp end parallel do
            ! GRID & 3D REC
            if(p%norec .eq. 'no')then
                do iptcl = fromp, top
                    orientation = b%a%get_ori(iptcl)
                    state       = nint(orientation%get('state'))
                    if(state == 0)cycle
                    b%img = batch_imgs(iptcl)
                    call grid_ptcl(b, p, orientation)
                enddo
            endif
            ! ORIENTATIONS OUPTUT: only here for debug
            do iptcl = fromp, top
                call b%a%write(iptcl, p%outfile)
            enddo
            ! CLEANUP BATCH
            do iptcl = fromp, top
                !call pcont3Dsrchs(iptcl)%kill
                call pftccs(iptcl)%kill
                call batch_imgs(iptcl)%kill
            enddo
            deallocate(pftccs, batch_imgs)
            !deallocate(pftccs, pcont3Dsrchs, batch_imgs)
        enddo
        ! orientations output
        !call b%a%write(p%outfile, [p%fromp,p%top])
        p%oritab = p%outfile

        ! NORMALIZE STRUCTURE FACTORS
        if(p%norec .eq. 'no')then
            if( p%eo .eq. 'yes' )then
                call eonorm_struct_facts(b, p, reslim, which_iter)
            else
                call norm_struct_facts(b, p, which_iter)
            endif
        endif

        ! REPORT CONVERGENCE
        if( p%l_distr_exec )then
            call qsys_job_finished( p, 'simple_pcont3D_matcher :: cont3D_exec')
        else
            ! CONVERGENCE TEST
            converged = b%conv%check_conv3D(update_res)
        endif
        
        ! DEALLOCATE
        deallocate(batches)
        do state=1,p%nstates
            if(state_exists(state))then
                call b%refvols(state)%kill_expanded
                call b%refvols(state)%kill
            endif
        enddo
        deallocate(state_exists)
        call b%img%kill_expanded
    end subroutine pcont3D_exec

    subroutine prep_vols( b, p, cline )
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        integer :: state
        ! PREPARATION OF VOLUMES FOR PROJECTION
        do state=1,p%nstates
            if( state_exists(state) )then
                call preprefvol( b, p, cline, state, doexpand=.false. )
                b%refvols(state) = b%vol
                call b%refvols(state)%expand_cmat
            endif
        enddo
        if( debug )write(*,*)'prep volumes done'
        ! bring back the original b%vol size
        if( p%boxmatch < p%box )call b%vol%new([p%box,p%box,p%box], p%smpd) ! to double check
    end subroutine prep_vols

    subroutine prep_pftcc(b, p, iptcl, pftcc)
        class(build),               intent(inout) :: b
        class(params),              intent(inout) :: p
        integer,                    intent(in)    :: iptcl
        class(polarft_corrcalc),    intent(inout) :: pftcc
        type(oris) :: cone
        type(ori)  :: optcl, oref
        real       :: eullims(3,2)
        integer    :: state, iref, cnt
        optcl = b%a%get_ori(iptcl)
        ! RE-INIT PFTCC
        if( p%l_xfel )then
            call pftcc%new(nrefs_per_ptcl, [iptcl,iptcl], [p%boxmatch,p%boxmatch,1],p%kfromto, p%ring2, p%ctf, isxfel='yes')
        else
            call pftcc%new(nrefs_per_ptcl, [iptcl,iptcl], [p%boxmatch,p%boxmatch,1],p%kfromto, p%ring2, p%ctf)
        endif
        ! PREP PARTICLE
        call read_img_from_stk(b, p, iptcl)
        call prepimg4align(b, p, optcl)
        call b%img%imgpolarizer(pftcc, iptcl, isptcl=.true.)
        ! restores b%img dimensions for clean exit
        if(p%boxmatch < p%box)call b%img%new([p%box,p%box,1],p%smpd)
    end subroutine prep_pftcc

end module simple_cont3D_matcher