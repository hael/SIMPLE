module simple_cont3D_matcher
use simple_defs
use simple_cartft_corrcalc,  only: cartft_corrcalc
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
implicit none

public :: cont3D_exec
private

type(cartft_corrcalc)   :: cftcc
type(polarft_corrcalc)  :: pftcc
integer, parameter      :: NRESTARTS = 5
integer, parameter      :: MAXNPEAKS = 10
character(len=STDLEN)   :: OPT_STR   = 'simplex'
logical, parameter      :: debug = .false.

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
            case('cart')
                call cont3D_cart_exec(b, p, cline, which_iter, converged)
            case('polar')
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
    subroutine cont3D_cart_exec( b, p, cline, which_iter, converged )
        use simple_cftcc_srch        ! use all in there
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: which_iter
        logical,        intent(inout) :: converged
        type(ori)              :: orientation
        integer                :: state, iptcl, cnt_glob

        ! CREATE CALCULATOR
        call cftcc%new( b, p, cline )
        if(p%boxmatch < p%box)call b%vol%new([p%box,p%box,p%box],p%smpd) ! ensures correct dimensions

        ! INITIALIZE
        call cftcc_srch_init(cftcc, b%img, OPT_STR, p%optlims(:5,:), NRESTARTS)
        write(*,'(A)',advance='no') '>>> CARTESIAN-FT '
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
                call cftcc_srch_set_state(state)
                call cftcc_srch_minimize(orientation)
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
        call cftcc_srch_reset
        call cftcc%kill
    end subroutine cont3D_cart_exec

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

    ! !>  \brief  is the continuous refinement algorithm
    ! subroutine cont3D_polar_exec_dev( b, p, cline, which_iter, converged )
    !     use simple_map_reduce, only: split_nobjs_even
    !     use simple_projector,  only: projector
    !     use simple_pftcc_contsrch    ! use all in there
    !     class(build),   intent(inout) :: b
    !     class(params),  intent(inout) :: p
    !     class(cmdline), intent(inout) :: cline
    !     integer,        intent(in)    :: which_iter
    !     logical,        intent(inout) :: converged
    !     ! batch related variables
    !     type(projector),        allocatable :: batch_imgs(:)
    !     type(polarft_corrcalc), allocatable :: pftccs(:)
    !     !type(pcont3D_srch),     allocatable :: pcont3Dsrchs(:)
    !     integer,                allocatable :: batches(:,:)
    !     ! other variables
    !     type(polarft_corrcalc) :: pftcc
    !     type(ori)              :: orientation
    !     integer                :: state, iptcl, cnt_glob, nptcls
    !     ! INIT
    !     nptcls = p%top - p%fromp + 1       ! number of particles processed
    !     ! batches
    !     nbatches = ceiling(real(nptcls)/real(p%nthr*BATCHSZ_MUL)) ! congruence with number of threads
    !     batches  = split_nobjs_even(nptcls, nbatches)


    !     ! INIT IMAGE PROJECTION
    !     ! dummy pftcc is only init here so the img polarizer can be initialized
    !     if( p%l_xfel )then
    !         call pftcc%new(1, [1,1], [p%boxmatch,p%boxmatch,1], p%kfromto, p%ring2, p%ctf, isxfel='yes')
    !     else
    !         call pftcc%new(1, [1,1], [p%boxmatch,p%boxmatch,1], p%kfromto, p%ring2, p%ctf)
    !     endif
    !     call b%img%init_imgpolarizer(pftcc)
    !     call pftcc%kill

    !     ! VOLS NEED TO BE PREPPED HERE TOGETHER WITH ORIENTATIONS CENTERING

    !     ! INITIALIZE
    !     !call pftcc_contsrch_init(b, p, cline, pftcc, b%img, OPT_STR, NRESTARTS)
    !     !if(p%boxmatch < p%box)call b%vol%new([p%box,p%box,p%box],p%smpd) ! ensures correct dimensions
    !     write(*,'(A)',advance='no') '>>> POLAR-FT '
    !     if( which_iter <= 0 )then
    !         write(*,'(A)')'CONTINUOUS ORIENTATION SEARCH'
    !     else
    !         write(*,'(A,1X,I3)')'CONTINUOUS ORIENTATION SEARCH, ITERATION:', which_iter
    !         if( which_iter > 0 ) p%outfile = 'cont3Ddoc_'//int2str_pad(which_iter,3)//'.txt'
    !     endif

    !     ! BATCH PROCESSING
    !     call del_file(p%outfile)      
    !     do batch = 1, nbatches
    !         ! BATCH INDICES
    !         fromp = p%fromp-1 + batches(batch,1)
    !         top   = p%fromp-1 + batches(batch,2)
    !         ! PREP BATCH
    !         allocate(pftccs(fromp:top), &!pcont3Dsrchs(fromp:top),&
    !             &batch_imgs(fromp:top), stat=alloc_stat)
    !         call alloc_err('In cont3D_matcher::cont3D_polar_exec_single',alloc_stat)
    !         do iptcl = fromp, top
    !             orientation = b%a%get_ori(iptcl)              
    !             state = nint(orientation%get('state'))
    !             if(state == 0)cycle
    !             ! stash raw image for rec
    !             call read_img_from_stk(b, p, iptcl)
    !             batch_imgs(iptcl) = b%img
    !             ! prep pftccs & search
    !     if( p%l_xfel )then
    !         call pftcc%new(nrefs_per_ptcl, [1,1], [p%boxmatch,p%boxmatch,1],p%kfromto, p%ring2, p%ctf, isxfel='yes')
    !     else
    !         call pftcc%new(nrefs_per_ptcl, [1,1], [p%boxmatch,p%boxmatch,1],p%kfromto, p%ring2, p%ctf)
    !     endif



    !             call prep_pftcc(b, p, iptcl, pftccs(iptcl), ptcl_img=b%img)
    !             !call pcont3Dsrchs(iptcl)%new(p, orientation, orefs, pftccs(iptcl))
    !         enddo
    !         ! SEARCH
    !         !$omp parallel do default(shared) schedule(guided) private(iptcl)
    !         do iptcl = fromp, top
    !             !call pcont3Dsrchs(iptcl)%do_srch
    !             call pftccs(iptcl)%kill ! cleanup
    !         enddo
    !         !$omp end parallel do
    !         ! GRID & 3D REC
    !         do iptcl = fromp, top
    !             orientation = b%a%get_ori(iptcl)              
    !             state       = nint(orientation%get('state'))
    !             if(state == 0)then
    !                 call orientation%reject
    !             else
    !                 orientation = pcont3Dsrchs(iptcl)%get_best_ori()
    !                 if(p%norec .eq. 'no')then
    !                     b%img = batch_imgs(iptcl)
    !                     call grid_ptcl(b, p, orientation)
    !                 endif
    !             endif
    !             ! set output
    !             call b%a%set_ori(iptcl,orientation)
    !             call b%a%write(iptcl, p%outfile)
    !             ! cleanup
    !             !call pcont3Dsrchs(iptcl)%kill
    !             call batch_imgs(iptcl)%kill
    !         enddo
    !         ! cleanup         
    !         deallocate(pftccs, pcont3Dsrchs, batch_imgs)
    !     enddo
    !     p%oritab = p%outfile
        
    !     ! cleanup
    !     if( p%boxmatch < p%box )call b%img%new([p%boxmatch,p%boxmatch,1],p%smpd) ! for next iteration in local execution
    !     call b%refvols(1)%kill_expanded
    !     call pftcc_contsrch_reset
    !     call pftcc%kill
    ! end subroutine cont3D_polar_exec_dev

    ! subroutine prep_pftcc(b, p, iptcl, pftcc)
    !     use simple_projector, only: projector
    !     class(build),               intent(inout) :: b
    !     class(params),              intent(inout) :: p
    !     integer,                    intent(in)    :: iptcl
    !     class(polarft_corrcalc),    intent(inout) :: pftcc
    !     type(ori) :: optcl
    !     optcl = b%a%get_ori(iptcl)
    !     ! INIT PFTCC
    !     if( p%l_xfel )then
    !         call pftcc%new(nrefs_per_ptcl, [1,1], [p%boxmatch,p%boxmatch,1],p%kfromto, p%ring2, p%ctf, isxfel='yes')
    !     else
    !         call pftcc%new(nrefs_per_ptcl, [1,1], [p%boxmatch,p%boxmatch,1],p%kfromto, p%ring2, p%ctf)
    !     endif
    !     ! PREP PARTICLE
    !     call read_img_from_stk(b, p, iptcl)
    !     call prepimg4align(b, p, optcl)
    !     call b%img%imgpolarizer(pftcc, 1, isptcl=.true.)
    !     ! restores b%img dimensions for clean exit
    !     if(p%boxmatch < p%box)call b%img%new([p%box,p%box,1],p%smpd)
    ! end subroutine prep_pftcc

end module simple_cont3D_matcher