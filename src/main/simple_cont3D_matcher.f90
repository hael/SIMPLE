! high-level continuous 3D projection matching
module simple_cont3D_matcher
use simple_build,             only: build
use simple_params,            only: params
use simple_cmdline,           only: cmdline
use simple_polarft_corrcalc,  only: polarft_corrcalc
use simple_oris,              only: oris
use simple_ori,               only: ori
use simple_cont3D_ada_srch,   only: cont3D_ada_srch
use simple_hadamard_common   ! use all in there
use simple_math              ! use all in there
use simple_binoris_io        ! use all in there
use simple_defs              ! use all in there
implicit none

public :: cont3D_exec
private
#include "simple_local_flags.inc"
integer,                   parameter :: BATCHSZ_MUL = 10   !< particles per thread
integer,                   parameter :: NREFS       = 100  !< Max num of references
type(cont3D_ada_srch),   allocatable :: cont3Dadasrch(:)   !< pftcc array for refine=ada
logical, allocatable                 :: state_exists(:)
integer                              :: nptcls          = 0 !< number of particle images per part
integer                              :: nrefs_per_ptcl  = 0 !< number of references per particle 
integer                              :: neff_states     = 0

contains

    !>  \brief  is the 3D continous algorithm
    subroutine cont3D_exec( b, p, cline, which_iter, converged )
        use simple_map_reduce, only: split_nobjs_even
        use simple_qsys_funs,  only: qsys_job_finished
        use simple_strings,    only: int2str_pad
        use simple_projector,  only: projector
        use simple_fileio,     only: del_file
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(build),   intent(inout) :: b               !< build object
        class(params),  intent(inout) :: p               !< params object
        class(cmdline), intent(inout) :: cline           !< command line input
        integer,        intent(in)    :: which_iter      !< iteration
        logical,        intent(inout) :: converged       !< have converged
        ! batches-related variables
        type(projector),        allocatable :: batch_imgs(:)
        type(polarft_corrcalc), allocatable :: pftccs(:)
        integer,                allocatable :: batches(:,:)
        integer                             :: nbatches, batch
        ! other variables
        type(polarft_corrcalc) :: pftcc            !< dummy convenience pftcc
        real, allocatable      :: eopart(:)
        type(oris)             :: softoris
        type(ori)              :: orientation
        real                   :: reslim
        integer                :: fromp, top, iptcl, state, alloc_stat,iptcl_tmp
        ! MULTIPLE STATES DEACTIVATED FOR NOW
        if(p%nstates>1)stop 'MULTIPLE STATES DEACTIVATED FOR NOW; cont3D_matcher::cont3Dexec'
        ! INIT
        nptcls = p%top - p%fromp + 1               ! number of particles processed
        ! states
        allocate(state_exists(p%nstates))
        state_exists = b%a%states_exist(p%nstates) ! state existence
        neff_states  = count(state_exists)         ! number of non-empty states
        ! number of references per particle
        select case(p%refine)
            case('yes')
                nrefs_per_ptcl = NREFS * neff_states
            case DEFAULT
                stop 'Unknown refinement mode; pcont3D_matcher::cont3D_exec'
        end select
        ! batches
        nbatches = ceiling(real(nptcls)/real(p%nthr*BATCHSZ_MUL))
        batches  = split_nobjs_even(nptcls, nbatches)

        ! SET BAND-PASS LIMIT RANGE
        call set_bp_range( b, p, cline )
        reslim = p%lp

        ! CALCULATE ANGULAR THRESHOLD
        if( .not.cline%defined('athres') ) p%athres = max(p%lp, ATHRES_LIM)
        write(*,'(A,F6.2)')'>>> ANGULAR THRESHOLD: ', p%athres

        ! DETERMINE THE NUMBER OF PEAKS
        select case(p%refine)
            case('yes')
                if( .not. cline%defined('npeaks') )then
                    if( p%eo .eq. 'yes' )then
                        p%npeaks = min(b%e%find_npeaks_from_athres(NPEAKSATHRES), MAXNPEAKS)
                    else
                        p%npeaks = min(10,b%e%find_npeaks(p%lp, p%moldiam))
                    endif
                endif
                p%npeaks = min(NREFS, p%npeaks)
            case DEFAULT
                stop 'Unknown refinement mode; pcont3D_matcher::cont3D_exec'
        end select
        write(*,'(A,I2)')'>>> NPEAKS: ', p%npeaks

        ! SETUP WEIGHTS FOR THE 3D RECONSTRUCTION
        if( p%nptcls <= SPECWMINPOP )then
            call b%a%calc_hard_weights(p%frac)
        else
            call b%a%calc_spectral_weights(p%frac)
        endif

        ! PREPARE REFERENCE & RECONSTRUCTION VOLUMES
        call prep_vols(b, p, cline)
        if(p%norec .eq. 'no')then
            call preprecvols(b, p)
            if( p%eo.eq.'yes')call prep_eopairs(b, p, eopart)
        endif

        ! INIT IMGPOLARIZER
        call pftcc%new(nrefs_per_ptcl, [1,1], [p%boxmatch,p%boxmatch,1], p%smpd, p%kfromto, p%ring2, p%ctf)
        call b%img_match%init_polarizer(pftcc)
        call pftcc%kill

        ! INITIALIZE
        write(*,'(A,1X,I3)')'>>> CONTINUOUS ORIENTATION SEARCH, ITERATION:', which_iter
        if( .not. p%l_distr_exec )then
            p%outfile = 'cont3Ddoc_'//int2str_pad(which_iter,3)//METADATEXT
        endif

        ! BATCH PROCESSING
        call del_file(p%outfile)
        do batch = 1, nbatches
            ! BATCH INDICES
            fromp = p%fromp-1 + batches(batch,1)
            top   = p%fromp-1 + batches(batch,2)
            ! PREP BATCH
            allocate(pftccs(fromp:top), cont3Dadasrch(fromp:top),&
            &batch_imgs(fromp:top), stat=alloc_stat)
            call alloc_errchk('In pcont3D_matcher::pcont3D_exec',alloc_stat)
            do iptcl = fromp, top
                state = nint(b%a%get(iptcl, 'state'))
                if(state == 0)cycle
                ! stash raw image for rec
                call read_img_from_stk(b, p, iptcl)
                call batch_imgs(iptcl)%copy(b%img)
                ! prep pftccs & ctf
                call init_pftcc(p, iptcl, pftccs(iptcl))
                iptcl_tmp = iptcl
                 call prep_pftcc_ptcl(b, p, iptcl_tmp, pftccs(iptcl_tmp))
                if( p%ctf.ne.'no' )call pftccs(iptcl_tmp)%create_polar_ctfmats(b%a)
                select case(p%refine)
                    case('yes')
                        call cont3Dadasrch(iptcl)%new(p, pftccs(iptcl), b%refvols, b%se)
                    case DEFAULT
                        stop 'Unknown refinement mode; pcont3D_matcher::cont3D_exec'
                end select
            enddo
            ! SERIAL SEARCHES
            !$omp parallel do default(shared) schedule(guided) private(iptcl) proc_bind(close)
            do iptcl = fromp, top
                call cont3Dadasrch(iptcl)%exec_srch(b%a, iptcl)
            enddo
            !$omp end parallel do
            ! ORIENTATIONS OUTPUT: only here for now
            call binwrite_oritab(p%outfile, b%a, [fromp,top])

            ! GRID & 3D REC
            if(p%norec .eq. 'no')then
                do iptcl = fromp, top
                    orientation = b%a%get_ori(iptcl)
                    state       = nint(orientation%get('state'))
                    if(state == 0)cycle
                    call b%img%copy(batch_imgs(iptcl))
                    if(p%npeaks == 1)then
                        if( p%eo.eq.'yes' )then
                            ! call grid_ptcl(b, p, orientation, ran_eo=eopart(iptcl) )
                            call grid_ptcl(b, p, orientation )
                        else
                            call grid_ptcl(b, p, orientation)
                        endif
                    else
                        softoris = cont3Dadasrch(iptcl)%get_peaks()
                        if( p%eo.eq.'yes' )then
                            !call grid_ptcl(b, p, orientation, os=softoris, ran_eo=eopart(iptcl) )
                            call grid_ptcl(b, p, orientation, os=softoris)
                        else
                            call grid_ptcl(b, p, orientation, os=softoris )
                        endif
                    endif
                enddo
            endif
            ! CLEANUP BATCH
            do iptcl = fromp, top
                call cont3Dadasrch(iptcl)%kill
                call pftccs(iptcl)%kill
                call batch_imgs(iptcl)%kill
            enddo
            deallocate(pftccs, cont3Dadasrch, batch_imgs)
        enddo

        ! CLEANUP SEARCH
        call b%img_match%kill_polarizer
        do state=1,p%nstates
            call b%refvols(state)%kill_expanded
            call b%refvols(state)%kill
        enddo
        if(allocated(eopart))deallocate(eopart)
        deallocate(batches, state_exists)

        ! ORIENTATIONS OUTPUT
        ! call binwrite_oritab(p%outfile, b%a, [p%fromp,p%top])
        p%oritab = p%outfile

        ! NORMALIZE STRUCTURE FACTORS
        if(p%norec .eq. 'no')then
            call b%vol%new([p%box, p%box,p%box], p%smpd)
            if( p%eo .eq. 'yes' )then
                call eonorm_struct_facts(b, p, reslim, which_iter)
            else
                call norm_struct_facts(b, p, which_iter)
            endif
            call killrecvols(b, p)
        endif

        ! REPORT CONVERGENCE
        if( p%l_distr_exec )then
            call qsys_job_finished( p, 'simple_pcont3D_matcher :: cont3D_exec')
        else
            converged = b%conv%check_conv_cont3D()
        endif
    end subroutine cont3D_exec

    !>  \brief  preps volumes for projection
    subroutine prep_vols( b, p, cline )
        class(build),   intent(inout) :: b           !< build object
        class(params),  intent(inout) :: p           !< params object
        class(cmdline), intent(inout) :: cline       !< command line input
        integer :: state
        do state=1,p%nstates
            if( state_exists(state) )then
                call b%vol%new([p%box,p%box,p%box], p%smpd) 
                call preprefvol( b, p, cline, state, doexpand=.false. )
                b%refvols(state) = b%vol
                call b%refvols(state)%expand_cmat
            endif
        enddo
        call b%vol%kill
        ! bring back the original b%vol size
        if( p%boxmatch < p%box )call b%vol%new([p%box,p%box,p%box], p%smpd) ! to double check
        DebugPrint 'prep volumes done'
    end subroutine prep_vols

    !>  \brief  initialize pftcc
    subroutine init_pftcc(p, iptcl, pftcc)
        class(params),              intent(inout) :: p        !< params object
        integer,                    intent(in)    :: iptcl    !< index to particle
        class(polarft_corrcalc),    intent(inout) :: pftcc    !< calculator object
        call pftcc%new(nrefs_per_ptcl, [iptcl,iptcl], [p%boxmatch,p%boxmatch,1], p%smpd, p%kfromto, p%ring2, p%ctf)
    end subroutine init_pftcc

    !>  \brief  particle projection into pftcc
    subroutine prep_pftcc_ptcl(b, p, iptcl, pftcc)
        class(build),               intent(inout) :: b        !< build object
        class(params),              intent(inout) :: p        !< params object
        integer,                    intent(in)    :: iptcl       !< index to particle
        class(polarft_corrcalc),    intent(inout) :: pftcc    !< calculator object
        type(ori)  :: optcl
        optcl = b%a%get_ori(iptcl)
        call prepimg4align(b, p, optcl)
        call b%img_match%polarize(pftcc, iptcl, isptcl=.true.)
    end subroutine prep_pftcc_ptcl

    !>  \brief  is for balanced distribution of the even/odd pairs
    subroutine prep_eopairs(b, p, eopart)
        use simple_ran_tabu, only: ran_tabu
        class(build),      intent(inout) :: b
        class(params),     intent(inout) :: p
        real, allocatable, intent(out)   :: eopart(:)
        type(ran_tabu)       :: rt
        real,    allocatable :: rec_weights(:)
        integer, allocatable :: part(:)
        integer              :: iptcl, cnt, n_recptcls
        allocate(eopart(p%fromp:p%top), source=-1.) ! -1. is default excluded value
        rec_weights = b%a%get_all('w')
        n_recptcls  = count(rec_weights(p%fromp:p%top) > TINY)
        rt          = ran_tabu( n_recptcls )
        allocate(part(n_recptcls), source=0)
        call rt%balanced(2, part)
        part = part - 1 ! 0/1 outcome
        cnt = 0
        do iptcl = p%fromp, p%top
            if(rec_weights(iptcl) > TINY)then
                cnt = cnt + 1
                eopart(iptcl) = real(part(cnt))
            endif
        enddo
        if(count(eopart(p%fromp:p%top)>=0.).ne.n_recptcls)stop 'indexing error; simple_cont3Dmatcher%prep_eopairs'
        call rt%kill
        deallocate(rec_weights, part)
    end subroutine prep_eopairs

end module simple_cont3D_matcher
