module simple_pcont3D_matcher
use simple_defs
use simple_build,            only: build
use simple_params,           only: params
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_oris,             only: oris
use simple_ori,              only: ori
use simple_cmdline,          only: cmdline
!use simple_masker,           only: automask
use simple_pcont3D_srch,     only: pcont3D_srch
use simple_hadamard_common,  ! use all in there
use simple_math              ! use all in there
implicit none

public :: pcont3D_exec, pcont3D_exec_single
private

type(polarft_corrcalc) :: pftcc
type(oris)             :: orefs                   !< per particle projection direction search space
logical, allocatable   :: state_exists(:)
integer                :: nptcls          = 0
integer                :: nrefs_per_ptcl  = 0
integer                :: neff_states     = 0

integer, parameter     :: BATCHSZ_MUL = 5   ! particles per thread
integer, parameter     :: MAXNPEAKS   = 10
integer, parameter     :: NREFS       = 50
logical, parameter     :: debug = .false.

contains

    !>  \brief  is the prime3D algorithm
    subroutine pcont3D_exec( b, p, cline, which_iter, converged )
        use simple_qsys_funs, only: qsys_job_finished
        use simple_strings,   only: int2str_pad
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: which_iter
        logical,        intent(inout) :: converged
        type(pcont3D_srch)            :: pcont3Dsrch
        type(oris)                    :: softoris
        type(ori)                     :: orientation
        real                          :: reslim, frac_srch_space
        integer                       :: iptcl, state, cnt_glob
        logical                       :: update_res
        ! AUTOMASKING DEACTIVATED FOR NOW
        ! INIT
        nptcls          = p%top - p%fromp + 1       ! number of particles processed
        state_exists    = b%a%get_state_exist(p%nstates)    ! state existence
        neff_states     = count(state_exists)       ! number of non-empty states
        nrefs_per_ptcl  = NREFS*neff_states         ! number of references per particle
        frac_srch_space = b%a%get_avg('frac')       ! Fraction of the search space

        ! SET BAND-PASS LIMIT RANGE
        call set_bp_range( b, p, cline )
        reslim = p%lp

        ! CALCULATE ANGULAR THRESHOLD (USED BY THE SPARSE WEIGHTING SCHEME)
        p%athres = rad2deg(atan(max(p%fny,p%lp)/(p%moldiam/2.)))
        write(*,'(A,F8.2)')'>>> ANGULAR THRESHOLD:', p%athres

        ! DETERMINE THE NUMBER OF PEAKS
        if( .not. cline%defined('npeaks') )then
            p%npeaks = min(MAXNPEAKS,b%e%find_npeaks(p%lp, p%moldiam))
        endif
        write(*,'(A,I3)')'>>> NPEAKS:', p%npeaks

        ! SETUP WEIGHTS FOR THE 3D RECONSTRUCTION
        if( p%frac < 0.99 ) call b%a%calc_hard_ptcl_weights(p%frac)

        ! PREPARE REFVOLS
        call prep_vols(b, p, cline)

        ! RESET RECVOLS
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

        ! INIT PFTCC & IMGPOLARIZER
        if( p%l_xfel )then
            call pftcc%new(nrefs_per_ptcl, [1,1], [p%boxmatch,p%boxmatch,1],p%kfromto, p%ring2, p%ctf, isxfel='yes')
        else
            call pftcc%new(nrefs_per_ptcl, [1,1], [p%boxmatch,p%boxmatch,1],p%kfromto, p%ring2, p%ctf)
        endif
        ! the pftcc is only initialized here so the img polarizer can be
        call b%img%init_imgpolarizer(pftcc)

        ! INITIALIZE
        if( which_iter <= 0 )then
            write(*,'(A)')'>>> CONTINUOUS POLAR-FT ORIENTATION SEARCH'
        else
            write(*,'(A,1X,I3)')'>>> CONTINUOUS POLAR-FT ORIENTATION SEARCH, ITERATION:', which_iter
            p%outfile = 'cont3Ddoc_'//int2str_pad(which_iter,3)//'.txt'
        endif

        ! ALIGN & GRID
        call del_file(p%outfile)
        cnt_glob = 0
        do iptcl=p%fromp,p%top
            cnt_glob = cnt_glob + 1
            call progress(cnt_glob, nptcls)
            ! orientation to align
            orientation = b%a%get_ori(iptcl)
            state = nint(orientation%get('state'))
            if(state == 0) then
                call orientation%reject
                call b%a%set_ori(iptcl,orientation)
            else
                ! re-fills pftcc
                call prep_pftcc(b, p, iptcl, pftcc)
                ! CTF matrices
                call apply_ctf(p, orientation, pftcc)
                ! align
                call pcont3Dsrch%new(p, orientation, orefs, pftcc)
                call pcont3Dsrch%do_srch
                orientation = pcont3Dsrch%get_best_ori()
                call b%a%set_ori(iptcl, orientation)
                ! grid
                call read_img_from_stk( b, p, iptcl )
                if(p%npeaks == 1)then
                    call grid_ptcl(b, p, iptcl, orientation)
                else
                    softoris = pcont3Dsrch%get_softoris()
                    call grid_ptcl(b, p, iptcl, orientation, os=softoris)
                endif
            endif
            ! output orientation
            call b%a%write(iptcl, p%outfile)
        enddo

        ! orientations output
        !call b%a%write(p%outfile, [p%fromp,p%top])
        p%oritab = p%outfile

        ! NORMALIZE STRUCTURE FACTORS
        if( p%eo .eq. 'yes' )then
            call eonorm_struct_facts(b, p, reslim, which_iter)
        else
            call norm_struct_facts(b, p, which_iter)
        endif
        ! REPORT CONVERGENCE
        if( p%l_distr_exec )then
            call qsys_job_finished( p, 'simple_pcont3D_matcher :: cont3D_exec')
        else
            ! CONVERGENCE TEST
            converged = b%conv%check_conv3D(update_res)
        endif
        ! DEALLOCATE
        call pftcc%kill
        do state=1,p%nstates
            if(state_exists(state))call b%refvols(state)%kill_expanded
        enddo
        !call b%img%kill_imgpolarizer is private a the moment
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

    subroutine prep_pftcc(b, p, iptcl, pftcc, ptcl_img)
        use simple_jiffys,    only: alloc_err
        use simple_projector, only: projector
        class(build),               intent(inout) :: b
        class(params),              intent(inout) :: p
        integer,                    intent(in)    :: iptcl
        class(polarft_corrcalc),    intent(inout) :: pftcc
        class(projector), optional, intent(in)    :: ptcl_img 
        type(oris) :: cone
        type(ori)  :: optcl, oref
        real       :: eullims(3,2)
        integer    :: state, iref, cnt
        optcl = b%a%get_ori(iptcl)
        ! RE-INIT PFTCC
        if( p%l_xfel )then
            call pftcc%new(nrefs_per_ptcl, [1,1], [p%boxmatch,p%boxmatch,1],p%kfromto, p%ring2, p%ctf, isxfel='yes')
        else
            call pftcc%new(nrefs_per_ptcl, [1,1], [p%boxmatch,p%boxmatch,1],p%kfromto, p%ring2, p%ctf)
        endif
        ! SEARCH SPACE PREP
        eullims = b%se%srchrange()
        call cone%rnd_proj_space(NREFS, optcl, p%athres, eullims)
        call cone%set_euler(1, optcl%get_euler()) ! previous best is the first
        do iref = 1, NREFS
            call cone%e3set(iref, 0.)
        enddo
        ! replicates to states
        if( p%nstates==1 )then
            call cone%set_all2single('state', 1.)
            orefs = cone
        else
            call orefs%new(NREFS*neff_states)
            cnt = 0
            do state = 1,p%nstates
                if(state_exists(state))then
                    call cone%set_all2single('state',real(state))
                    do iref=1,NREFS
                        cnt = cnt+1
                        call orefs%set_ori(cnt, cone%get_ori(iref))
                    enddo
                endif
            enddo
        endif
        ! REFERENCES PROJECTION
        do iref=1,nrefs_per_ptcl
            oref  = orefs%get_ori(iref)
            state = nint(oref%get('state'))
            call b%refvols(state)%fproject_polar(iref, oref, pftcc, expanded=.true.)
        enddo
        ! PREP PARTICLE
        if(present(ptcl_img))then
            if(any( ptcl_img%get_ldim()-b%img%get_ldim() .ne. 0))& ! dimension check
                &stop 'Incompatible image diensions in pcont3D_matcher::prep_pftcc'
            b%img = ptcl_img
        else
            call read_img_from_stk(b, p, iptcl)
        endif
        call prepimg4align(b, p, optcl)
        call b%img%imgpolarizer(pftcc, 1, isptcl=.true.)
        ! restores b%img dimensions for clean exit
        if(p%boxmatch < p%box)call b%img%new([p%box,p%box,1],p%smpd)
    end subroutine prep_pftcc

    subroutine apply_ctf(p, o, pftcc)
        use simple_ctf, only: ctf
        class(params),           intent(inout) :: p
        class(ori),              intent(inout) :: o
        class(polarft_corrcalc), intent(inout) :: pftcc
        type(oris) :: weirdos
        type(ctf)  :: tfun
        real       :: kV, cs, fraca, dfx, dfy, angast
        if( p%ctf .ne. 'no' )then
            ! grab CTF parameters
            dfx = o%get('dfx')
            if( p%tfplan%mode .eq. 'astig' )then ! astigmatic CTF
                dfy    = o%get('dfy')
                angast = o%get('angast')
            else if( p%tfplan%mode .eq. 'noastig' )then
                dfy    = dfx
                angast = 0.
            else
                stop 'Unsupported ctf mode; apply_ctf; simple_pcont3D_matcher'
            endif
            kV    = o%get('kv')
            cs    = o%get('cs')
            fraca = o%get('fraca')
            ! transfer function
            tfun = ctf(p%smpd, kV, cs, fraca)
            ! create CTF polar matrices for online application
            call weirdos%new(1)
            call weirdos%set_ori(1, o)
            call pftcc%create_polar_ctfmats(p%smpd, weirdos)
            call weirdos%kill
        endif
    end subroutine apply_ctf

    !>  \brief  is the prime3D algorithm
    subroutine pcont3D_exec_single( b, p, cline, which_iter, converged )
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

        type(projector),        allocatable :: batch_imgs(:)
        type(polarft_corrcalc), allocatable :: pftccs(:)
        type(pcont3D_srch),     allocatable :: pcont3Dsrchs(:)
        integer,                allocatable :: batches(:,:)
        type(oris) :: softoris
        type(ori)  :: orientation
        real       :: reslim, frac_srch_space
        integer    :: nbatches, batch, fromp, top, iptcl, state, alloc_stat
        logical    :: update_res
        ! AUTOMASKING DEACTIVATED FOR NOW
        ! INIT
        nptcls          = p%top - p%fromp + 1       ! number of particles processed
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
        write(*,'(A,F8.2)')'>>> ANGULAR THRESHOLD:', p%athres

        ! DETERMINE THE NUMBER OF PEAKS
        if( .not. cline%defined('npeaks') )then
            p%npeaks = min(MAXNPEAKS,b%e%find_npeaks(p%lp, p%moldiam))
        endif
        write(*,'(A,I3)')'>>> NPEAKS:', p%npeaks

        ! SETUP WEIGHTS FOR THE 3D RECONSTRUCTION
        if( p%nptcls <= SPECWMINPOP )then
            call b%a%calc_hard_ptcl_weights(p%frac)
        else
            call b%a%calc_spectral_weights(p%frac)
        endif

        ! PREPARE REFVOLS
        call prep_vols(b, p, cline)

        ! RESET RECVOLS
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

        ! INIT IMGPOLARIZER
        ! dummy pftcc is only init here so the img polarizer can be initialized
        if( p%l_xfel )then
            call pftcc%new(nrefs_per_ptcl, [1,1], [p%boxmatch,p%boxmatch,1],p%kfromto, p%ring2, p%ctf, isxfel='yes')
        else
            call pftcc%new(nrefs_per_ptcl, [1,1], [p%boxmatch,p%boxmatch,1],p%kfromto, p%ring2, p%ctf)
        endif
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
            allocate(pftccs(fromp:top), pcont3Dsrchs(fromp:top),&
                &batch_imgs(fromp:top), stat=alloc_stat)
            call alloc_err('In pcont3D_matcher::pcont3D_exec_single',alloc_stat)
            do iptcl = fromp, top
                orientation = b%a%get_ori(iptcl)              
                state = nint(orientation%get('state'))
                if(state == 0)cycle
                ! stash raw image for rec
                call read_img_from_stk(b, p, iptcl)
                batch_imgs(iptcl) = b%img
                ! prep pftccs & search
                call prep_pftcc(b, p, iptcl, pftccs(iptcl), ptcl_img=b%img)
                call apply_ctf(p, orientation, pftccs(iptcl))
                call pcont3Dsrchs(iptcl)%new(p, orientation, orefs, pftccs(iptcl))
            enddo
            ! SEARCH
            !$omp parallel do default(shared) schedule(auto) private(iptcl)
            do iptcl = fromp, top
                call pcont3Dsrchs(iptcl)%do_srch
                call pftccs(iptcl)%kill ! cleanup
            enddo
            !$omp end parallel do
            ! GRID & 3D REC
            do iptcl = fromp, top
                orientation = b%a%get_ori(iptcl)              
                state       = nint(orientation%get('state'))
                if(state == 0)then
                    call orientation%reject
                else
                    orientation = pcont3Dsrchs(iptcl)%get_best_ori()
                    if(p%norec .eq. 'no')then
                        b%img = batch_imgs(iptcl)
                        if(p%npeaks == 1)then
                            call grid_ptcl(b, p, iptcl, orientation)
                        else
                            softoris = pcont3Dsrchs(iptcl)%get_softoris()
                            call grid_ptcl(b, p, iptcl, orientation, os=softoris)
                        endif
                    endif
                endif
                ! set output
                call b%a%set_ori(iptcl,orientation)
                call b%a%write(iptcl, p%outfile)
                ! cleanup
                call pcont3Dsrchs(iptcl)%kill
                call batch_imgs(iptcl)%kill
            enddo
            ! cleanup         
            deallocate(pftccs, pcont3Dsrchs, batch_imgs)
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
            if(state_exists(state))call b%refvols(state)%kill_expanded
        enddo
        !call b%img%kill_imgpolarizer is private a the moment
    end subroutine pcont3D_exec_single

end module simple_pcont3D_matcher

