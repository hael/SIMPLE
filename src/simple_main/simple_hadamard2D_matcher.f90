module simple_hadamard2D_matcher
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_prime_srch,       only: prime_srch
use simple_prime2D_srch,     only: prime2D_srch
use simple_ori,              only: ori
use simple_build,            only: build
use simple_params,           only: params
use simple_cmdline,          only: cmdline
use simple_hadamard_common   ! singleton
use simple_defs              ! singleton
use simple_jiffys            ! singleton
use simple_filterer          ! singleton
implicit none

public :: prime2D_exec, prime2D_assemble_sums, prime2D_norm_sums, prime2D_assemble_sums_from_parts,&
prime2D_write_sums, preppftcc4align, pftcc, primesrch2D, prime2D_read_sums, prime2D_write_partial_sums
private

type(polarft_corrcalc) :: pftcc       ! need to be revealed to the outside world for testing purposes
type(prime2D_srch)     :: primesrch2D ! need to be revealed to the outside world for testing purposes
type(ori)              :: orientation
integer                :: cnt_glob = 0
logical, parameter     :: debug=.false.

contains
    
    !>  \brief  is the prime2D algorithm
    subroutine prime2D_exec( b, p, cline, which_iter, converged )
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline     
        integer,        intent(in)    :: which_iter
        logical,        intent(inout) :: converged
        integer :: iptcl, fnr, icls, file_stat

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
        
        ! IN GPU MODE WE CALCULATE ALL THE CORRELATIONS BEFORE THE SEARCH BEGINS
        if( p%use_gpu .eq. 'yes' .or. p%bench_gpu .eq. 'yes' )then
            if( p%bench_gpu .eq. 'yes' .and. p%use_gpu .eq. 'no')then
                if( p%oritab .eq. '' )then
                    call primesrch2D%calc_corrs(pftcc, mode='bench')
                else
                     call primesrch2D%calc_corrs(pftcc, b%a, [p%fromp,p%top], 'bench')
                endif
            else
                if( p%oritab .eq. '' )then
                    call primesrch2D%calc_corrs(pftcc, mode='cpu') ! 4 NOW
                    ! call primesrch2D%calc_corrs(pftcc, mode='gpu')
                else
                    call primesrch2D%calc_corrs(pftcc, b%a, [p%fromp,p%top], 'cpu') ! 4 NOW
                    ! call primesrch2D%calc_corrs(pftcc, b%a, [p%fromp,p%top], 'gpu')
                endif
            endif
        endif
        
        ! INITIALISE SUMS
        call prime2D_init_sums( b, p )

        ! ALIGN & GRID
        call del_txtfile(p%outfile)
        cnt_glob = 0
        if( debug ) write(*,*) '*** hadamard2D_matcher ***: loop fromp/top:', p%fromp, p%top
        do iptcl=p%fromp,p%top
            cnt_glob = cnt_glob+1
            call progress(cnt_glob, p%top-p%fromp+1)        
            orientation = b%a%get_ori(iptcl)
            if( nint(orientation%get('state')) > 0 )then
                call preprefs4align(b, p, iptcl, pftcc)
                ! execute the high-level routines in prime2D_srch
                if(  p%use_gpu .eq. 'yes' .or. p%bench_gpu .eq. 'yes' )then
                    if( str_has_substr(p%refine,'neigh') ) stop 'refine=neigh modes not currently implemented on GPU'
                    if( p%oritab .eq. '' )then
                        call primesrch2D%exec_prime2D_srch(pftcc, iptcl, p%lp, cnt_glob=cnt_glob)
                    else
                        call primesrch2D%exec_prime2D_srch(pftcc, iptcl, p%lp, orientation, cnt_glob=cnt_glob)
                    endif
                else             
                    select case(p%refine)
                        case('no')
                            if( p%oritab .eq. '' )then
                                call primesrch2D%exec_prime2D_srch(pftcc, iptcl, p%lp)
                            else
                                call primesrch2D%exec_prime2D_srch(pftcc, iptcl, p%lp, orientation)
                            endif
                        case('neigh')
                            if( p%oritab .eq. '' )then
                                call primesrch2D%exec_prime2D_srch(pftcc, iptcl, p%lp, nnmat=b%nnmat)
                            else
                                call primesrch2D%exec_prime2D_srch(pftcc, iptcl, p%lp, orientation, nnmat=b%nnmat)
                            endif
                        case DEFAULT
                            write(*,*) 'The refinement mode: ', trim(p%refine), ' is unsupported on CPU'
                            stop
                    end select
                endif
                call primesrch2D%get_cls(orientation)
            else
                call orientation%reject
            endif
            call b%a%set_ori(iptcl,orientation)
            call b%a%write(iptcl, p%outfile)
            ! read back the image again 4 shift, rotation and cavg update
            if( p%l_distr_exec )then
                call b%img%read(p%stk_part, cnt_glob)
            else
                call b%img%read(p%stk, iptcl)
            endif
            if( nint(orientation%get('state')) > 0 )then
                icls = nint(orientation%get('class'))
                ! the wiener_restore2D_online modifies b%img
                call wiener_restore2D_online(b%img, orientation,&
                b%tfun, p%tfplan, b%cavgs(icls), p%msk)
                call assemble_ctfsqsum_online(b%img, orientation,&
                b%tfun, p%tfplan, b%ctfsqsums(icls))
            endif
        end do
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
            ! generation of this file marks completion of the partition
            ! this file is empty 4 now but may contain run stats etc.
            fnr = get_fileunit()
            open(unit=fnr, FILE='JOB_FINISHED_'//int2str_pad(p%part,p%numlen),&
            STATUS='REPLACE', action='WRITE', iostat=file_stat)
            call fopen_err( 'In: prime_ini; simple_hadamard2D_matcher.f90', file_stat )
            close(fnr)
        else
            ! CONVERGENCE TEST
            converged = b%conv%check_conv2D()
        endif
        
    end subroutine prime2D_exec
    
    subroutine prime2D_read_sums( b, p )
        use simple_jiffys, only: file_exists
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
            b%cavgs(icls)     = 0.
            b%ctfsqsums(icls) = cmplx(0.,0.)
        end do
        !$omp end parallel do 
    end subroutine prime2D_init_sums
    
    subroutine prime2D_assemble_sums( b, p )
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        type(ori) :: orientation
        integer   :: icls, iptcl, cnt, istart, iend
        write(*,'(a)') '>>> ASSEMBLING CLASS SUMS'
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
                icls = nint(orientation%get('class'))
                call orientation%print
                if( p%l_distr_exec )then
                    call b%img%read( p%stk_part, cnt )
                else
                    call b%img%read(p%stk, iptcl)
                endif
                call wiener_restore2D_online(b%img, orientation,&
                b%tfun, p%tfplan, b%cavgs(icls), p%msk)
                call assemble_ctfsqsum_online(b%img, orientation,&
                b%tfun, p%tfplan, b%ctfsqsums(icls))             
            endif
        end do
        if( .not. p%l_distr_exec ) call prime2D_norm_sums( b, p )
    end subroutine prime2D_assemble_sums
    
    subroutine prime2D_assemble_sums_from_parts( b, p )
        use simple_jiffys, only: file_exists
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
        use simple_image,  only: image
        use simple_masker, only: automask2D
        class(build),  intent(inout) :: b
        class(params), intent(inout) :: p
        type(image)                  :: ref
        type(ori) :: o
        integer   :: cnt, iptcl, icls, sz, pop
        write(*,'(A)') '>>> BUILDING PRIME2D SEARCH ENGINE'
        ! must be done here since constants in p are dynamically set
        call primesrch2D%new(p)
        write(*,'(A)') '>>> BUILDING DATA STRUCTURE FOR POLARFT CORRELATION CALCULATION'
        call pftcc%new(p%ncls, [p%fromp,p%top], [p%box,p%box,1], p%kfromto, p%ring2, p%ctf)
        ! PREPARATION OF REFERENCES IN PFTCC
        ! read references and transform into polar coordinates
        write(*,'(A)') '>>> BUILDING REFERENCES'
        sz = b%img%get_lfny(1)
        do icls=1,p%ncls
            call progress(icls, p%ncls)
            pop = 2 
            if( p%oritab /= '' ) pop = b%a%get_clspop(icls)
            if( pop > 1 )then
                ref = b%cavgs(icls)
                call prep2Dref( p, ref )
                ! transfer to polar coordinates
                call b%proj%img2polarft(icls, ref, pftcc, isptcl=.false.)
                call ref%kill
            endif
        end do
        ! PREPARATION OF PARTICLES IN PFTCC
        ! read particle images and create polar projections
        write(*,'(A)') '>>> BUILDING PARTICLES'
        cnt = 0
        do iptcl=p%fromp,p%top
            cnt = cnt+1
            call progress(cnt, p%top-p%fromp+1)
            if( p%l_distr_exec )then
                call b%img%read(p%stk_part, cnt)
            else
                call b%img%read(p%stk, iptcl)
            endif
            o = b%a%get_ori(iptcl)
            call prepimg4align(b, p, o)
            call b%proj%img2polarft(iptcl, b%img, pftcc)
        end do
        if( debug ) write(*,*) '*** hadamard2D_matcher ***: finished preppftcc4align'
    end subroutine preppftcc4align
    
end module simple_hadamard2D_matcher
