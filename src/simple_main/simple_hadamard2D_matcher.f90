module simple_hadamard2D_matcher
use simple_defs
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_prime_srch,       only: prime_srch
use simple_prime2D_srch,     only: prime2D_srch
use simple_ori,              only: ori
use simple_build,            only: build
use simple_params,           only: params
use simple_cmdline,          only: cmdline
use simple_strings,          only: int2str_pad
use simple_jiffys,           only: progress
use simple_filehandling      ! use all in there
use simple_hadamard_common   ! use all in there
use simple_filterer          ! use all in there
implicit none

public :: prime2D_exec, prime2D_assemble_sums, prime2D_norm_sums, prime2D_assemble_sums_from_parts,&
prime2D_write_sums, preppftcc4align, pftcc, primesrch2D, prime2D_read_sums, prime2D_write_partial_sums
private

logical, parameter              :: DEBUG = .false.
type(polarft_corrcalc)          :: pftcc
type(prime2D_srch), allocatable :: primesrch2D(:)

contains
    
    !>  \brief  is the prime2D algorithm
    subroutine prime2D_exec( b, p, cline, which_iter, converged )
        use simple_qsys_funs,   only: qsys_job_finished
        use simple_strings,     only: str_has_substr
        use simple_procimgfile, only: random_selection_from_imgfile
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline     
        integer,        intent(in)    :: which_iter
        logical,        intent(inout) :: converged
        integer, allocatable :: ptcls2process(:) ! ptcl inds 2 process
        logical   :: srch_shifts(p%fromp:p%top)  ! per-particle shift flags
        integer   :: iptcl, fnr, icls, io_stat, inorm, cands(3), pop, iproc, nproc
        real      :: corr_thresh, frac_srch_space
        type(ori) :: orientation
        
        ! SET FRACTION OF SEARCH SPACE
        frac_srch_space = b%a%get_avg('frac')

        ! PER-PARTICLE CONVERGENCE
        if( file_exists(p%ppconvfile) )then
            call b%ppconv%read(p%ppconvfile)
        else
            call b%ppconv%zero_joint_distr_olap
        endif
        call b%ppconv%set_shift_larr(srch_shifts)
        ptcls2process = b%ppconv%identify_ptcls2process()
        nproc = size(ptcls2process)
        write(*,'(A,F8.2)') '>>> NON-CONVERGED PARTICLES(%):', 100.*(real(nproc) / real(p%top - p%fromp + 1))

        ! PREP REFERENCES
        if( p%l_distr_exec )then
            if( .not. cline%defined('refs') )then
                stop 'need refs to be part of command line for distributed prime2D execution'
            else if( cline%defined('refs') )then
                if( .not. file_exists(p%refs) ) stop 'input references (refs) does not exist in cwd'
            endif
            call prime2D_read_sums( b, p )
        else
            ! for shared-memory or chunk-based parallellisation we need initial references for iter=1 only
            if( which_iter == p%startit )then
                if( cline%defined('refs') )then
                    if( .not. file_exists(p%refs) ) stop 'input references (refs) does not exist in cwd'
                    call prime2D_read_sums( b, p )
                else
                    ! we need to make references
                    if( cline%defined('oritab') )then
                        ! we make class averages
                        call prime2D_assemble_sums(b, p)
                    else
                        ! we randomly select particle images as initial references
                        p%refs = 'start2Drefs'//p%ext
                        if( p%chunktag .ne. '' ) p%refs = trim(p%chunktag)//trim(p%refs)
                        call random_selection_from_imgfile(p%stk, p%refs, p%ncls, p%smpd)
                        call prime2D_read_sums( b, p )
                    endif
                endif
            endif 
        endif

        ! SETUP WEIGHTS
        if( p%nptcls <= SPECWMINPOP )then
            call b%a%calc_hard_ptcl_weights(p%frac)
        else
            call b%a%calc_spectral_weights(p%frac)
        endif

        ! EXTREMAL LOGICS
        if( frac_srch_space < 0.98 .or. p%extr_thresh > 0.025 )then
            corr_thresh = b%a%extremal_bound(p%extr_thresh)
        else
            corr_thresh = -huge(corr_thresh)
        endif
        
        ! SET FOURIER INDEX RANGE
        call set_bp_range2D( b, p, cline, which_iter, frac_srch_space )
        
        ! GENERATE REFERENCE & PARTICLE POLAR FTs
        call preppftcc4align( b, p )

        ! INITIALIZE
        if( which_iter <= 0 )then
            write(*,'(A)')       '>>> PRIME2D DISCRETE STOCHASTIC SEARCH'
        else
            write(*,'(A,1X,I3)') '>>> PRIME2D DISCRETE STOCHASTIC SEARCH, ITERATION:', which_iter
        endif
        if( .not. p%l_distr_exec .and. which_iter > 0 )then
            p%outfile = 'prime2Ddoc_'//int2str_pad(which_iter,3)//'.txt'
            if( p%chunktag .ne. '' ) p%outfile= trim(p%chunktag)//trim(p%outfile)
        endif
        call prime2D_init_sums( b, p )

        ! STOCHASTIC IMAGE ALIGNMENT
        allocate( primesrch2D(nproc) )
        do iproc=1,nproc
            call primesrch2D(iproc)%new(p, pftcc) 
        end do
        ! calculate CTF matrices
        if( p%ctf .ne. 'no' ) call pftcc%create_polar_ctfmats(p%smpd, b%a)
        ! execute the search
        if( p%refine .eq. 'neigh' )then
            call del_file(p%outfile)
            !$omp parallel do default(shared) schedule(auto) private(iproc,iptcl)
            do iproc=1,nproc
                iptcl = ptcls2process(iproc)
                call primesrch2D(iproc)%nn_srch(pftcc, iptcl, b%a, b%nnmat, srch_shifts(iptcl))
            end do
            !$omp end parallel do
        else
            ! execute the search
            call del_file(p%outfile)
            if( p%oritab .eq. '' )then
                !$omp parallel do default(shared) schedule(auto) private(iproc,iptcl)
                do iproc=1,nproc
                    iptcl = ptcls2process(iproc)
                    call primesrch2D(iproc)%exec_prime2D_srch(pftcc, iptcl, b%a, srch_shifts(iptcl), greedy=.true.)
                end do
                !$omp end parallel do
            else
                if(corr_thresh > 0.)then
                    write(*,'(A,F8.2)') '>>> PARTICLE RANDOMIZATION(%):', 100.*p%extr_thresh
                    write(*,'(A,F8.2)') '>>> CORRELATION THRESHOLD:    ', corr_thresh
                endif
                !$omp parallel do default(shared) schedule(auto) private(iproc,iptcl)
                do iproc=1,nproc
                    iptcl = ptcls2process(iproc)
                    call primesrch2D(iproc)%exec_prime2D_srch(pftcc, iptcl, b%a, srch_shifts(iptcl), extr_bound=corr_thresh)
                end do
                !$omp end parallel do
            endif
        endif
        if( DEBUG ) print *, 'DEBUG, hadamard2D_matcher; completed alignment'
        ! output orientations
        call b%a%write(p%outfile, [p%fromp,p%top])
        p%oritab = p%outfile
        
        ! WIENER RESTORATION OF CLASS AVERAGES
        do iptcl=p%fromp,p%top
            orientation = b%a%get_ori(iptcl)
            if( nint(orientation%get('state')) > 0 )then
                call read_img_from_stk( b, p, iptcl )
                icls = nint(orientation%get('class'))
                call wiener_restore2D_online(b%img, orientation, p%tfplan,&
                &b%cavgs(icls), b%ctfsqsums(icls), p%msk)
            endif
        end do
        if( DEBUG ) print *, 'DEBUG, hadamard2D_matcher; generated class averages'
        
        ! WRITE CLASS AVERAGES
        if( p%l_distr_exec )then
            call prime2D_write_partial_sums( b, p )
        else
            call prime2D_norm_sums( b, p )
            call prime2D_write_sums( b, p, which_iter )
        endif

        ! DESTRUCT
        do iproc=1,nproc
            call primesrch2D(iproc)%kill
        end do
        deallocate( primesrch2D )
        call pftcc%kill

        ! REPORT CONVERGENCE
        call b%ppconv%update_joint_distr_olap
        call b%ppconv%write(p%ppconvfile)
        if( p%l_distr_exec )then
            call qsys_job_finished(p, 'simple_hadamard2D_matcher :: prime2D_exec')
        else
            ! CONVERGENCE TEST
            converged = b%conv%check_conv2D()
        endif
    end subroutine prime2D_exec
    
    subroutine prime2D_read_sums( b, p )
        class(build),  intent(inout) :: b
        class(params), intent(inout) :: p
        integer :: icls
        if( file_exists(p%refs) )then
            do icls=1,p%ncls
                call b%cavgs(icls)%read(p%refs, icls)                
            end do
        else
            write(*,*) 'File does not exists: ', trim(p%refs)
            stop 'In: simple_hadamard2D_matcher :: prime2D_read_sums'
        endif
    end subroutine prime2D_read_sums
    
    subroutine prime2D_init_sums( b, p )
        class(build),  intent(inout) :: b
        class(params), intent(inout) :: p
        integer :: icls
        !$omp parallel do schedule(auto) default(shared) private(icls)
        do icls=1,p%ncls
            b%cavgs(icls) = 0.
            b%ctfsqsums(icls) = cmplx(0.,0.)
        end do
        !$omp end parallel do 
    end subroutine prime2D_init_sums
    
    subroutine prime2D_assemble_sums( b, p )
        use simple_ctf, only: ctf
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        type(ori) :: orientation
        integer   :: icls, iptcl, cnt, istart, iend, filtsz
        if( .not. p%l_distr_exec ) write(*,'(a)') '>>> ASSEMBLING CLASS SUMS'
        filtsz = b%img_pad%get_filtsz()
        call prime2D_init_sums( b, p )
        if( p%l_distr_exec )then
            istart  = p%fromp
            iend    = p%top
        else
            istart  = 1
            iend    = p%nptcls
        endif
        cnt = 0
        do iptcl=istart,iend
            cnt = cnt+1
            call progress( cnt, iend-istart+1 )
            orientation = b%a%get_ori(iptcl)
            if( nint(orientation%get('state')) > 0 )then
                call read_img_from_stk( b, p, iptcl )
                icls = nint(orientation%get('class'))
                call wiener_restore2D_online(b%img, orientation, p%tfplan,&
                &b%cavgs(icls), b%ctfsqsums(icls), p%msk)
            endif
        end do
        if( .not. p%l_distr_exec ) call prime2D_norm_sums( b, p )
    end subroutine prime2D_assemble_sums
    
    subroutine prime2D_assemble_sums_from_parts( b, p )
        class(build),  intent(inout) :: b
        class(params), intent(inout) :: p
        character(len=STDLEN) :: fname_cavgs, fname_ctfsqsums
        integer :: ipart, icls
        call prime2D_init_sums( b, p )
        do ipart=1,p%nparts
            fname_cavgs     = 'cavgs_part'//int2str_pad(ipart,p%numlen)//p%ext
            fname_ctfsqsums = 'ctfsqsums_part'//int2str_pad(ipart,p%numlen)//p%ext
            ! read & sum partial class averages
            if( file_exists(fname_cavgs) )then
                do icls=1,p%ncls
                    call b%img%read(fname_cavgs, icls)
                    ! add subaverage to class
                    call b%cavgs(icls)%add(b%img)
                end do
            else
                write(*,*) 'File does not exists: ', trim(fname_cavgs)
                stop 'In: simple_hadamard2D_matcher :: prime2D_assemble'
            endif
            ! read & sum partial ctfsqsums
            if( file_exists(fname_ctfsqsums) )then
                do icls=1,p%ncls
                    call b%img%read(fname_ctfsqsums, icls)
                    ! add subaverage to class
                    call b%ctfsqsums(icls)%add(b%img)
                end do
            else
                write(*,*) 'File does not exists: ', trim(fname_ctfsqsums)
                stop 'In: simple_hadamard2D_matcher :: prime2D_assemble'
            endif
        end do
        call prime2D_norm_sums( b, p )
    end subroutine prime2D_assemble_sums_from_parts
    
    subroutine prime2D_write_partial_sums( b, p )
        class(build),  intent(inout) :: b
        class(params), intent(inout) :: p
        integer :: icls
        do icls=1,p%ncls
            call b%cavgs(icls)%write('cavgs_part'//int2str_pad(p%part,p%numlen)//p%ext, icls)
            call b%ctfsqsums(icls)%write('ctfsqsums_part'//int2str_pad(p%part,p%numlen)//p%ext, icls)
        end do 
    end subroutine prime2D_write_partial_sums

    subroutine prime2D_norm_sums( b, p )
        class(build),  intent(inout) :: b
        class(params), intent(inout) :: p
        integer :: icls, pop
        do icls=1,p%ncls
            pop = b%a%get_cls_pop(icls)
            if( pop > 1 )then
                call b%cavgs(icls)%fwd_ft
                call b%cavgs(icls)%ctf_dens_correct(b%ctfsqsums(icls))
                call b%cavgs(icls)%bwd_ft
            endif
        end do
    end subroutine prime2D_norm_sums
    
    subroutine prime2D_write_sums( b, p, which_iter, fname )
        class(build),               intent(inout) :: b
        class(params),              intent(inout) :: p
        integer,          optional, intent(in)    :: which_iter
        character(len=*), optional, intent(in)    :: fname
        integer :: icls
        if( present(which_iter) )then
            if( present(fname) ) stop&
            &'fname cannot be present together with which_iter; simple_hadamard2D_matcher :: prime2D_write_sums'
            if( which_iter <= 0 )then
                p%refs = 'cavgs'//p%ext
            else
                p%refs = 'cavgs_iter'//int2str_pad(which_iter,3)//p%ext
            endif
        else
            if( present(fname) )then
                p%refs = fname
            else
                p%refs = 'startcavgs'//p%ext
            endif
        endif
        if( p%chunktag .ne. '' ) p%refs = trim(p%chunktag)//trim(p%refs)
        ! write to disk
        do icls=1,p%ncls
            call b%cavgs(icls)%write(p%refs, icls)
        end do
    end subroutine prime2D_write_sums

    !>  \brief  prepares the polarft corrcalc object for search
    subroutine preppftcc4align( b, p )
        use simple_masker,       only: automask2D
        use simple_jiffys,       only: alloc_err
        class(build),  intent(inout) :: b
        class(params), intent(inout) :: p
        type(ori) :: o
        integer   :: cnt, iptcl, icls, sz, pop, istate
        integer   :: filtsz, alloc_stat, filnum, io_stat
        if( .not. p%l_distr_exec ) write(*,'(A)') '>>> BUILDING PRIME2D SEARCH ENGINE'
        ! must be done here since constants in p are dynamically set
        call pftcc%new(p%ncls, [p%fromp,p%top], [p%box,p%box,1], p%kfromto, p%ring2, p%ctf)
        ! prepare the polarizers
        call b%img%init_imgpolarizer(pftcc)
        ! PREPARATION OF REFERENCES IN PFTCC
        ! read references and transform into polar coordinates
        if( .not. p%l_distr_exec ) write(*,'(A)') '>>> BUILDING REFERENCES'
        sz = b%img%get_lfny(1)
        do icls=1,p%ncls
            call progress(icls, p%ncls)
            pop = 2 
            if( p%oritab /= '' ) pop = b%a%get_cls_pop(icls)
            if( pop > 1 )then
                ! prepare the reference
                b%img = b%cavgs(icls)
                call prep2Dref(p, b%img, b%a, icls)
                b%refs(icls) = b%img
                ! transfer to polar coordinates
                call b%img%imgpolarizer(pftcc, icls, isptcl=.false.)
            endif
        end do
        ! PREPARATION OF PARTICLES IN PFTCC
        ! read particle images and create polar projections
        if( .not. p%l_distr_exec ) write(*,'(A)') '>>> BUILDING PARTICLES'
        cnt = 0
        do iptcl=p%fromp,p%top
            cnt = cnt+1
            call progress(cnt, p%top-p%fromp+1)
            call read_img_from_stk( b, p, iptcl )
            o      = b%a%get_ori(iptcl)
            icls   = nint(o%get('class'))
            istate = nint(o%get('state'))
            if( istate == 0 ) icls = 0
            call prepimg4align(b, p, o)
            ! transfer to polar coordinates
            call b%img%imgpolarizer(pftcc, iptcl)
        end do
        if( debug ) write(*,*) '*** hadamard2D_matcher ***: finished preppftcc4align'
    end subroutine preppftcc4align
    
end module simple_hadamard2D_matcher
