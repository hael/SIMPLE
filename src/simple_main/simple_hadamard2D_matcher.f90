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

type(polarft_corrcalc) :: pftcc       ! need to be revealed to the outside world for testing purposes
type(prime2D_srch)     :: primesrch2D ! need to be revealed to the outside world for testing purposes
type(ori)              :: orientation
integer                :: cnt_glob = 0
real                   :: frac_srch_space = 0.
integer, allocatable   :: defgroups(:)
real,    allocatable   :: wmat(:,:)
real,    allocatable   :: ctfparams(:,:)
logical, parameter     :: DEBUG = .false.
real,    parameter     :: prime2Deps = 0.3

contains
    
    !>  \brief  is the prime2D algorithm
    subroutine prime2D_exec( b, p, cline, which_iter, converged )
        use simple_qsys_funs, only: qsys_job_finished
        use simple_strings,   only: str_has_substr
        use simple_ran_tabu,  only: ran_tabu
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline     
        integer,        intent(in)    :: which_iter
        logical,        intent(inout) :: converged
        real, allocatable :: res(:), res_pad(:), wresamp(:)
        integer :: iptcl, fnr, icls, io_stat, inorm, cands(3)
        logical :: doshellweight
        
        ! SET FRACTION OF SEARCH SPACE
        frac_srch_space = b%a%get_avg('frac')

        ! READ IMAGES
        call read_imgs_from_stk( b, p )

        ! SETUP SHELLWEIGHTS
        call setup_shellweights( b, p, doshellweight, wmat, res, res_pad )
        
        ! SET FOURIER INDEX RANGE
        call set_bp_range( b, p, cline )
         
        ! READ REFERENCES
        call prime2D_read_sums( b, p )
        
        ! GENERATE REFERENCE & PARTICLE POLAR FTs
        call preppftcc4align( b, p )

        ! INITIALIZE
        if( which_iter <= 0 )then
            write(*,'(A)') '>>> PRIME2D DISCRETE STOCHASTIC SEARCH'
        else
            write(*,'(A,1X,I3)') '>>> PRIME2D DISCRETE STOCHASTIC SEARCH, ITERATION:', which_iter
        endif
        if( which_iter > 0 ) p%outfile = 'prime2Ddoc_'//int2str_pad(which_iter,3)//'.txt'

        ! INITIALISE SUMS
        call prime2D_init_sums( b, p )

        ! ALIGN & GRID
        call del_file(p%outfile)
        if( p%ctf .ne. 'no' ) call pftcc%create_polar_ctfmats(p%smpd, b%a)
        select case(p%refine)
            case('no','greedy')
                if( p%oritab .eq. '' )then
                    call primesrch2D%exec_prime2D_srch(pftcc, b%a, [p%fromp,p%top],&
                    &frac_srch_space, greedy=.true.)
                else
                    call primesrch2D%exec_prime2D_srch(pftcc, b%a, [p%fromp,p%top],&
                    &frac_srch_space)
                endif
            case DEFAULT
                write(*,*) 'The refinement mode: ', trim(p%refine), ' is unsupported'
                stop
        end select
        
        ! WIENER RESTORATION OF CLASS AVERAGES
        cnt_glob = 0
        do iptcl=p%fromp,p%top
            cnt_glob = cnt_glob + 1
            orientation = b%a%get_ori(iptcl)            
            if( nint(orientation%get('state')) > 0 )then
                b%img = b%imgs(iptcl) ! put the original image back
                icls = nint(orientation%get('class'))
                if( allocated(wmat) )then
                    wresamp = resample_filter(wmat(iptcl,:), res, res_pad)
                    call wiener_restore2D_online_fast(b%img, orientation, p%tfplan,&
                    &b%cavgs(icls), b%ctfsqsums(icls), p%msk, wresamp)
                else
                    call wiener_restore2D_online_fast(b%img, orientation, p%tfplan,&
                    &b%cavgs(icls), b%ctfsqsums(icls), p%msk)
                endif
            endif
        end do

        ! orientations output
        call b%a%write(p%outfile, [p%fromp,p%top])
        p%oritab = p%outfile
        call pftcc%kill
        
        ! status here: sums have been assembled but not normalised
        ! in parallel setting: write partial files and move on
        ! in serial setting: normalise sums and write to disk
        
        ! WRITE CLASS AVERAGES
        if( p%l_distr_exec )then
            call prime2D_write_partial_sums( b, p )
        else
            call prime2D_norm_sums( b, p )
            call prime2D_write_sums( b, p, which_iter )
        endif

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
        integer   :: icls, iptcl, cnt, istart, iend
        if( .not. p%l_distr_exec ) write(*,'(a)') '>>> ASSEMBLING CLASS SUMS'
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
                b%img = b%imgs(iptcl) ! put the original image back
                icls = nint(orientation%get('class'))
                call wiener_restore2D_online_fast(b%img, orientation, p%tfplan,&
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
            pop = b%a%get_clspop(icls)
            if( pop > 1 )then
                call b%cavgs(icls)%fwd_ft
                call b%cavgs(icls)%ctf_dens_correct(b%ctfsqsums(icls))
                call b%cavgs(icls)%bwd_ft
            endif
        end do
    end subroutine prime2D_norm_sums
    
    subroutine prime2D_write_sums( b, p, which_iter )
        class(build),      intent(inout) :: b
        class(params),     intent(inout) :: p
        integer, optional, intent(in)    :: which_iter
        integer                          :: icls
        if( present(which_iter) )then
            if( which_iter <= 0 )then
                p%refs = 'cavgs'//p%ext
            else
                p%refs = 'cavgs_iter'//int2str_pad(which_iter,3)//p%ext
            endif
        else
            p%refs = 'startcavgs'//p%ext
        endif
        ! write to disk
        do icls=1,p%ncls
            call b%cavgs(icls)%write(p%refs, icls)
        end do
    end subroutine prime2D_write_sums

    !>  \brief  prepares the polarft corrcalc object for search
    subroutine preppftcc4align( b, p )
        use simple_image,        only: image
        use simple_masker,       only: automask2D
        use simple_jiffys,       only: alloc_err
        class(build),  intent(inout) :: b
        class(params), intent(inout) :: p
        type(ori) :: o
        integer   :: cnt, iptcl, icls, sz, pop, istate
        integer   :: filtsz, alloc_stat, filnum, io_stat
        if( .not. p%l_distr_exec ) write(*,'(A)') '>>> BUILDING PRIME2D SEARCH ENGINE'
        if( frac_srch_space >= SHW_FRAC_LIM .and. p%oritab .ne. '' )then
            filtsz = b%img%get_filtsz()
            if( allocated(wmat) ) deallocate(wmat)
            allocate(wmat(p%top-p%fromp+1,filtsz), stat=alloc_stat)
            call alloc_err("In simple_hadamard2D_matcher :: preppftcc4align", alloc_stat)
            wmat = 1.0
        endif
        ! must be done here since constants in p are dynamically set
        call primesrch2D%new(p)
        call pftcc%new(p%ncls, [p%fromp,p%top], [p%box,p%box,1], p%kfromto, p%ring2, p%nthr, p%ctf)
        ! prepare the polarizers
        call b%img%init_imgpolarizer(pftcc)
        ! PREPARATION OF REFERENCES IN PFTCC
        ! read references and transform into polar coordinates
        if( .not. p%l_distr_exec ) write(*,'(A)') '>>> BUILDING REFERENCES'
        sz = b%img%get_lfny(1)
        do icls=1,p%ncls
            call progress(icls, p%ncls)
            pop = 2 
            if( p%oritab /= '' ) pop = b%a%get_clspop(icls)
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
            b%img  = b%imgs(iptcl) ! put the original image back
            o      = b%a%get_ori(iptcl)
            icls   = nint(o%get('class'))
            istate = nint(o%get('state'))
            if( istate == 0 ) icls = 0
            call prepimg4align(b, p, o)
            if( allocated(wmat) ) call calc_frc( b, p, o, icls, cnt, wmat )
            ! transfer to polar coordinates
            call b%img%imgpolarizer(pftcc, iptcl)
        end do
        if( allocated(wmat) )then
            filnum = get_fileunit()
            if( p%l_distr_exec )then  
                open(unit=filnum, status='REPLACE', action='WRITE',&
                file='shellweights_part'//int2str_pad(p%part,p%numlen)//'.bin', access='STREAM')
            else
                open(unit=filnum, status='REPLACE', action='WRITE',&
                file='shellweights.bin', access='STREAM')
            endif
            write(unit=filnum,pos=1,iostat=io_stat) wmat
            ! check if the write was successful
            if( io_stat .ne. 0 )then
                write(*,'(a,i0,2a)') '**ERROR(preppftcc4align): I/O error ',&
                io_stat, ' when writing shellweights*.bin'
                stop 'I/O error; preppftcc4align; simple_hadamard2D_matcher'
            endif
            close(filnum)
            deallocate(wmat)
        endif
        if( debug ) write(*,*) '*** hadamard2D_matcher ***: finished preppftcc4align'
    end subroutine preppftcc4align

    !>  \brief  calculates the FRC between the prepared reference
    !!          image and the prepared particle image
    subroutine calc_frc( b, p, o, icls, cnt_glob, wmat )
        use simple_ctf, only: ctf
        class(build),  intent(inout) :: b
        class(params), intent(in)    :: p 
        integer,       intent(in)    :: icls, cnt_glob
        class(ori),    intent(inout) :: o
        real,          intent(inout) :: wmat(:,:)
        real, allocatable :: res(:), corrs(:)
        type(image)       :: ref_local
        type(ctf)         :: tfun
        real              :: dfx, dfy, angast
        if( icls == 0 )then
            wmat(cnt_glob,:) = -1.
            return
        endif
        ref_local = b%refs(icls)
        if( p%tfplan%flag .ne. 'no' )then
            if( p%tfplan%mode .eq. 'astig' )then ! astigmatic CTF
                dfx    = o%get('dfx')
                dfy    = o%get('dfy')
                angast = o%get('angast')
            else if( p%tfplan%mode .eq. 'noastig' )then
                dfx    = o%get('dfx')
                dfy    = dfx
                angast = 0.
            else
                stop 'Unsupported ctf mode; simple_hadamard2D_matcher :: calc_frc'
            endif
            tfun = ctf(p%smpd, o%get('kv'), o%get('cs'), o%get('fraca'))
            call tfun%apply(ref_local, dfx, 'ctf', dfy, angast)
        endif
        ! calculate FRC    
        call ref_local%fsc(b%img, res, corrs)
        wmat(cnt_glob,:) = corrs
        call ref_local%kill
        deallocate(res, corrs)
    end subroutine calc_frc
    
end module simple_hadamard2D_matcher