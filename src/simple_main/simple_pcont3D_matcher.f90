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

public :: pcont3D_exec
private

type(polarft_corrcalc) :: pftcc
type(oris)             :: orefs                   !< per particle projection direction search space
logical, allocatable   :: state_exists(:)
integer                :: nptcls          = 0
integer                :: nrefs_per_ptcl  = 0
integer                :: neff_states     = 0

integer, parameter     :: MAXNPEAKS = 10
integer, parameter     :: NREFS     = 50
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
        integer                       :: iptcl, state, alloc_stat, cnt_glob
        logical                       :: update_res
        ! AUTOMASKING DEACTIVATED FOR NOW
        ! INIT
        nptcls = p%top - p%fromp + 1                ! number of particles processed
        ! state existence
        allocate(state_exists(p%nstates), stat=alloc_stat)
        do state = 1, p%nstates
            state_exists(state) = (b%a%get_statepop(state) > 0)
        enddo
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
                call prep_pftcc(b, p, iptcl)
                ! CTF matrices
                call apply_ctf(p, orientation, pftcc)
                ! align
                call pcont3Dsrch%new(p, orientation, orefs, pftcc)
                call pcont3Dsrch%do_srch
                orientation = pcont3Dsrch%get_best_ori()
                call b%a%set_ori(iptcl, orientation)
                ! grid
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

    subroutine prep_pftcc(b, p, iptcl)
        use simple_rnd, only: ran3
        class(build),  intent(inout) :: b
        class(params), intent(inout) :: p
        integer,       intent(in)    :: iptcl
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
        call read_img_from_stk(b, p, iptcl)
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

end module simple_pcont3D_matcher

