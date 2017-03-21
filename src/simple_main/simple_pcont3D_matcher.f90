module simple_pcont3D_matcher
use simple_defs
use simple_build,            only: build
use simple_params,           only: params
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_oris,             only: oris
use simple_ori,              only: ori
use simple_cmdline,          only: cmdline
!use simple_masker,           only: automask
use simple_cont3D_matcher,   only: cont3D_shellweight
use simple_pcont3D_srch,     only: pcont3D_srch
use simple_hadamard_common,  only: set_bp_range, norm_struct_facts, eonorm_struct_facts,&
                                &prepimg4align, preprefvol, setup_shellweights, grid_ptcl
use simple_math              ! use all in there
implicit none

public :: pcont3D_exec
private

type(polarft_corrcalc) :: pftcc
type(oris)             :: orefs                   !< per particle projection direction search space
type(oris)             :: spiral                  !< projection direction search space
logical, allocatable   :: state_exists(:)
integer                :: nptcls          = 0
integer                :: nrefs_per_ptcl  = 0
integer                :: neff_states     = 0

integer, parameter     :: MAXNPEAKS = 10
integer, parameter     :: NREFS     = 200
logical, parameter     :: debug = .false.

contains

    !>  \brief  is the prime3D algorithm
    subroutine pcont3D_exec( b, p, cline, which_iter, converged )
        use simple_filterer,  only: resample_filter
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
        real, allocatable             :: wmat(:,:), wresamp(:), res(:), res_pad(:)
        real                          :: reslim, frac_srch_space
        integer                       :: iptcl, state, alloc_stat, cnt_glob
        logical                       :: doshellweight, update_res
        ! AUTOMASKING DEACTIVATED FOR NOW
        ! INIT
        nptcls = p%top - p%fromp + 1                ! number of particles processed
        ! state existence
        allocate(state_exists(p%nstates), stat=alloc_stat)
        state_exists = .true.
        do state = 1, p%nstates
            state_exists(state) = (b%a%get_statepop(state) > 0)
        enddo
        neff_states = count(state_exists)           ! number of non-empty states
        nrefs_per_ptcl = NREFS*neff_states          ! number of references per particle
        frac_srch_space = b%a%get_avg('frac')       ! Fraction of the search space

        ! SET BAND-PASS LIMIT RANGE
        call set_bp_range( b, p, cline )
        reslim = p%lp

        ! CALCULATE ANGULAR THRESHOLD (USED BY THE SPARSE WEIGHTING SCHEME)
        ! Unused for now
        ! p%athres = rad2deg(atan(max(p%fny,p%lp)/(p%moldiam/2.)))
        ! write(*,'(A,F8.2)')'>>> ANGULAR THRESHOLD:', p%athres

        ! DETERMINE THE NUMBER OF PEAKS
        if( .not. cline%defined('npeaks') )then
            p%npeaks = min(MAXNPEAKS,b%e%find_npeaks(p%lp, p%moldiam)) ! to update dependence on b%e?
        endif
        write(*,'(A,I3)')'>>> NPEAKS:', p%npeaks

        ! SETUP WEIGHTS FOR THE 3D RECONSTRUCTION
        if( p%frac < 0.99 ) call b%a%calc_hard_ptcl_weights(p%frac)
        if( p%l_distr_exec )then
            ! nothing to do
        else
            if( p%l_shellw .and. frac_srch_space>=50. )call cont3D_shellweight(b, p, cline)
        endif
        call setup_shellweights(b, p, doshellweight, wmat, res, res_pad)

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

        ! SEARCH SPACE PREP
        call spiral%new(p%nspace * b%se%get_nsym())
        call spiral%spiral

        ! INIT PFTCC & IMGPOLARIZER
        if( p%l_xfel )then
            call pftcc%new(nrefs_per_ptcl, [1,1], [p%boxmatch,p%boxmatch,1],p%kfromto, p%ring2, p%ctf, isxfel='yes')
        else
            call pftcc%new(nrefs_per_ptcl, [1,1], [p%boxmatch,p%boxmatch,1],p%kfromto, p%ring2, p%ctf)
        endif
        call b%img%init_imgpolarizer(pftcc, p%smpd)

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
            orientation = b%a%get_ori(iptcl)
            state = nint(orientation%get('state'))
            if(state==0) then
                call orientation%reject
                call b%a%set_ori(iptcl,orientation)
                cycle
            endif
            ! re-fills pftcc
            call prep_pftcc(b, p, iptcl, cnt_glob )
            ! multiply refs with CTF
            call apply_ctf(p, orientation, pftcc)
            ! align
            call pcont3Dsrch%new(p, orientation, orefs, pftcc)
            call pcont3Dsrch%do_srch
            orientation = pcont3Dsrch%get_best_ori()
            call b%a%set_ori(iptcl, orientation)
            softoris = pcont3Dsrch%get_softoris()
            ! grid
            if( doshellweight )then
                wresamp = resample_filter(wmat(iptcl,:), res, res_pad)
                call grid_ptcl(b, p, iptcl, cnt_glob, orientation, softoris, shellweights=wresamp)
            else
                call grid_ptcl(b, p, iptcl, cnt_glob, orientation, softoris )
            endif
            call b%a%write(iptcl, p%outfile)
        enddo
        ! cleanup (mostly for debug purposes)
        call spiral%kill
        call pftcc%kill
        do state=1,p%nstates
            if(state_exists(state))call b%refvols(state)%kill_expanded
        enddo
        !call b%img%kill_imgpolarizer is private a the moment

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
        if( allocated(wmat)    ) deallocate(wmat)
        if( allocated(wresamp) ) deallocate(wresamp)
        if( allocated(res)     ) deallocate(res)
        if( allocated(res_pad) ) deallocate(res_pad)
    end subroutine pcont3D_exec

    subroutine prep_vols( b, p, cline )
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline
        integer :: state
        ! PREPARATION OF VOLUMES FOR PROJECTION
        do state=1,p%nstates
            if( state_exists(state) )then
                call preprefvol( b, p, cline, state )
                b%refvols(state) = b%vol
                call b%refvols(state)%expand_cmat
            endif
        enddo
        if( debug )write(*,*)'prep volumes done'
        ! bring back the original b%vol size
        if( p%boxmatch < p%box ) call b%vol%new([p%box,p%box,p%box], p%smpd) ! to double check
    end subroutine prep_vols

    subroutine prep_pftcc(b, p, iptcl, cnt_glob)
        use simple_rnd, only: ran3
        class(build),            intent(inout) :: b
        class(params),           intent(inout) :: p
        integer,                 intent(in)    :: iptcl
        integer,                 intent(in)    :: cnt_glob
        type(oris) :: cone
        type(ori)  :: optcl, oref, ostoch
        integer    :: inds(p%nspace*b%se%get_nsym())
        integer    :: state, iref, cnt, iref_start, half_nrefs
        optcl = b%a%get_ori(iptcl)
        ! RE-INIT PFTCC
        if( p%l_xfel )then
            call pftcc%new(nrefs_per_ptcl, [1,1], [p%boxmatch,p%boxmatch,1],p%kfromto, p%ring2, p%ctf, isxfel='yes')
        else
            call pftcc%new(nrefs_per_ptcl, [1,1], [p%boxmatch,p%boxmatch,1],p%kfromto, p%ring2, p%ctf)
        endif
        ! SEARCH SPACE PREP
        ! random rotation of spiral
        call ostoch%rnd_ori
        call ostoch%e3set(0.)
        call spiral%rot(ostoch)
        half_nrefs = ceiling(real(NREFS) / 2.)
        call spiral%find_closest_projs(optcl, inds)
        call cone%new(NREFS)
        ! previous orientation goes first
        call oref%new
        call oref%set_euler(optcl%get_euler())
        call oref%e3set(0.)
        call cone%set_ori(1, oref)
        ! fine grained first-half
        cnt = 1
        do iref = 1,p%nspace*b%se%get_nsym()
            oref = spiral%get_ori(inds(iref))
            if( b%se%within_asymunit(oref) )then
                cnt = cnt + 1
                if(cnt > half_nrefs)exit
                call oref%e3set(0.)
                call cone%set_ori(cnt, oref)
            endif
        enddo
        ! coarse grained second-half
        iref_start = iref
        do iref = iref_start,p%nspace*b%se%get_nsym()
            oref = spiral%get_ori(inds(iref))
            if( b%se%within_asymunit(oref) )then
                if(ran3() >= 0.5)cycle ! skips every second on average
                cnt = cnt + 1
                if(cnt > NREFS)exit
                call oref%e3set(0.)
                call cone%set_ori(cnt, oref)
            endif
        enddo
        ! replicates to states
        if( p%nstates==1 )then
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
            call b%refvols(state)%fproject_polar(iref, oref, p, pftcc, expanded=.true.)
        enddo
        ! PREP PARTICLE
        if( p%boxmatch < p%box )then
            ! back to the original size as b%img has been
            ! and will be clipped/padded in prepimg4align
            call b%img%new([p%box,p%box,1],p%smpd)
        endif
        if( p%l_distr_exec )then
            call b%img%read(p%stk_part, cnt_glob, p%l_xfel)
        else
            call b%img%read(p%stk, iptcl, p%l_xfel)
        endif
        call prepimg4align(b, p, optcl)
        call b%img%imgpolarizer(pftcc, 1, isptcl=.true.)
        ! restores b%img dimensions for clean exit
        if( p%boxmatch < p%box )call b%img%new([p%box,p%box,1],p%smpd)
    end subroutine prep_pftcc

    subroutine apply_ctf(p, o, pftcc)
        use simple_ctf, only: ctf
        class(params),           intent(inout) :: p
        class(ori),              intent(inout) :: o
        class(polarft_corrcalc), intent(inout) :: pftcc
        type(ctf) :: tfun
        real      :: kV, cs, fraca, dfx, dfy, angast
        if( p%ctf .ne. 'no' )then
            ! grabf CTF parameters
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
            ! CTF multiplication of references
            call pftcc%apply_ctf(tfun, dfx, dfy=dfy, angast=angast)
        endif
    end subroutine apply_ctf

end module simple_pcont3D_matcher

