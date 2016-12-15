module simple_hadamard2D_matcher
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_prime2D_srch,     only: prime2D_srch
use simple_ori,              only: ori
use simple_build,            only: build
use simple_params,           only: params
use simple_estimate_ssnr     ! singleton
use simple_hadamard_common   ! singleton
use simple_defs              ! singleton
use simple_cmdline           ! singleton
use simple_jiffys            ! singleton
use simple_filterer          ! singleton
implicit none

public :: prime2D_exec, prime2D_assemble_sums, prime2D_norm_sums, &
&prime2D_assemble_sums_from_parts, prime2D_write_sums, preppftcc4align
private

type(polarft_corrcalc) :: pftcc
type(prime2D_srch)     :: primesrch2D
type(ori)              :: orientation
integer                :: cnt_glob = 0
logical, parameter     :: debug=.false.
logical, parameter     :: test=.true.
real, parameter        :: SHTHRESH=0.0001

contains
    
    !>  \brief  is the prime2D algorithm
    subroutine prime2D_exec( b, p, which_iter, converged )
        class(build), intent(inout)  :: b
        class(params), intent(inout) :: p        
        integer, intent(in)          :: which_iter
        logical, intent(inout)       :: converged
        real                         :: corr
        real                         :: frac
        real                         :: frac_better, mi_hard, x, y
        integer                      :: iptcl, j, k, fnr, icls, file_stat
        logical                      :: didsrch
        character(len=STDLEN)        :: fname

        ! SET FOURIER INDEX RANGE
        call set_bp_range(b, p)
         
        ! READ REFERENCES
        call prime2D_read_sums( b, p )
        
!         ! ESTIMATE SSNR FOR THE CLUSTERS
!         b%ssnr = 1.
!         call estimate_cluster_ssnrs(b, p)
        
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
                if( p%oritab .eq. '' )then
                    call primesrch2D%exec_prime2D_srch(pftcc, iptcl, p%lp)
                else
                    call primesrch2D%exec_prime2D_srch(pftcc, iptcl, p%lp, orientation)
                endif
                call primesrch2D%get_cls(orientation)
            else
                call orientation%reject
            endif
            call b%a%set_ori(iptcl,orientation)
            call b%a%write(iptcl, p%outfile)
            ! read back the image again 4 shift, rotation and cavg update
            if( defined_cmd_arg('part') )then
                call b%img%read(p%stk_part, cnt_glob)
            else
                call b%img%read(p%stk, iptcl)
            endif
            if( nint(orientation%get('state')) > 0 )then
                icls = nint(orientation%get('class'))
                call wiener_restore2D_online(b%img, orientation,&
                b%tfun, p%ctfmode, p%ctf, b%refs(icls))
                ! The wiener_restore2D_online modifies the image so 
                ! that it is ready for noise power estimation
                call assemble_ctfsqsum_online(b%img, orientation,&
                b%tfun, p%ctfmode, p%ctf, b%ctfsqsums(icls))
!                 if( defined_cmd_arg('inner') )then
!                     call estimate_specnoise_online(b%refs(icls),&
!                     b%img, p%msk, b%spec_noise2D(icls,:), inner_width=[p%inner,p%width])
!                 else
!                     call estimate_specnoise_online(b%refs(icls),&
!                     b%img, p%msk, b%spec_noise2D(icls,:))
!                 endif
            endif            
            call orientation%kill
        end do
        p%oritab = p%outfile
        
        ! status here: sums have been assembled but not normalised
        ! in parallel setting: write partial files and move on
        ! in serial setting: normalise sums and write to disk
        
        ! WRITE CLASS AVERAGES
        if( defined_cmd_arg('part') )then
            call prime2D_write_partial_sums( b, p )
        else
            call prime2D_norm_sums( b, p )
            call prime2D_write_sums( b, p, which_iter )
        endif
        
        ! CONVERGENCE TEST
        if( .not. defined_cmd_arg('part') )then
            corr    = b%a%get_avg('corr')
            frac    = b%a%get_avg('frac')
            mi_hard = b%a%get_avg('mi_hard')
            write(*,'(A,19X,F7.4)') '>>> DISTRIBUTION OVERLAP:', mi_hard
            write(*,'(A,4X,F5.1)')  '>>> PERCENTAGE OF SEARCH SPACE SCANNED:', frac
            write(*,'(A,28X,F7.4)') '>>> CORRELATION:', corr
            ! dynamic shift search range update
            if( frac >= 90.0 )then
                if( .not. defined_cmd_arg('trs') .or. p%trs < 2.0 )then
                    p%trs     = 0.025*real(p%box)
                    p%trs     = max(2.0,p%trs) ! minimum 2 pixels shift
                    p%trs     = min(6.0,p%trs) ! maximum 6 pixels shift
                    p%doshift = .true.
                endif
            endif
            converged = .false.
            if( mi_hard > 0.95 .and. frac > 98. ) converged = .true.
        else
            corr = b%a%get_avg('corr',       fromto=[p%fromp,p%top])
            frac = b%a%get_avg('frac',       fromto=[p%fromp,p%top])
            mi_hard = b%a%get_avg('mi_hard', fromto=[p%fromp,p%top])
            fnr  = get_fileunit()
            open(unit=fnr, FILE='JOB_FINISHED_'//int2str_pad(p%part,p%numlen),&
            STATUS='REPLACE', action='WRITE', iostat=file_stat)
            call fopen_err( 'In: prime_ini; simple_hadamard2D_matcher.f90', file_stat )
            write(fnr,'(A,19X,F7.4)') '>>> DISTRIBUTION OVERLAP:', mi_hard
            write(fnr,'(A,4X,F5.1)')  '>>> PERCENTAGE OF SEARCH SPACE SCANNED:', frac
            write(fnr,'(A,28X,F7.4)') '>>> CORRELATION:', corr
            converged = .false.
            if( mi_hard > 0.95 .and. frac > 98. )then
                converged = .true.
                write(fnr,'(A)') '>>> CONVERGED: .YES.'
            else
                write(fnr,'(A)') '>>> CONVERGED: .NO.'
            endif
            close(fnr)
        endif
    end subroutine
    
    subroutine prime2D_read_sums( b, p )
        use simple_jiffys, only: file_exists
        class(build), intent(inout)   :: b
        class(params), intent(inout)  :: p
        character(len=:), allocatable :: fname_refs, fname_ctfsqsums, fname_noisespecs
        real, allocatable             :: noise_specs(:,:)
        integer                       :: icls
        allocate(fname_refs,       source=trim(p%refs))
        allocate(fname_ctfsqsums,  source='ctfsqsums'//p%ext)
        allocate(fname_noisespecs, source='noisespecs'//'.bin')
        ! read class averages
        if( file_exists(fname_refs) )then
            ! loop over clusters
            do icls=1,p%ncls
                call b%refs(icls)%read(fname_refs, icls)
            end do
        else
            write(*,*) 'File does not exists: ', trim(fname_refs)
            stop 'In: simple_hadamard2D_matcher :: prime2D_read_sums'
        endif
        ! read ctfsqsums
        if( file_exists(fname_ctfsqsums) )then
            ! loop over clusters
            do icls=1,p%ncls
                call b%ctfsqsums(icls)%read(fname_ctfsqsums, icls)
            end do
        else
            write(*,*) 'File does not exists: ', trim(fname_ctfsqsums)
            stop 'In: simple_hadamard2D_matcher :: prime2D_read_sums'
        endif        
        ! read noise power spectra
!         if( file_exists(fname_noisespecs) )then
!             if( allocated(b%spec_noise2D) ) deallocate(b%spec_noise2D)
!             b%spec_noise2D = file2arr2D(fname_noisespecs)
!         else
!             write(*,*) 'WARNING! File does not exists: ', trim(fname_noisespecs)
!             write(*,*) 'If this is the first iteration, all is ok!'
!             allocate( b%spec_noise2D(p%ncls,b%img%get_lfny(1)) )
!             b%spec_noise2D = 1.0
!         endif
        deallocate(fname_refs, fname_ctfsqsums, fname_noisespecs)
    end subroutine
    
!     subroutine estimate_cluster_ssnrs( b, p )
!         class(build), intent(inout)  :: b
!         class(params), intent(inout) :: p
!         integer           :: icls, pop
!         real, allocatable :: spec_sig2D(:)
!         write(*,'(A)') '>>> ESTIMATING CLUSTER SSNRS'
!         b%ssnr = 0.
!         do icls=1,p%ncls
!             call progress(icls,p%ncls)
!             pop = b%a%get_clspop(icls)
!             if( pop > 1 )then
!                 spec_sig2D = b%refs(icls)%spectrum('power')
!                 where( b%spec_noise2D(icls,:) > 1e-6 )
!                     b%ssnr(icls,:) = spec_sig2D/b%spec_noise2D(icls,:)
!                 end where
!                 deallocate(spec_sig2D)
!             endif
!         end do
!     end subroutine
    
!     subroutine print_average_frc( b, p )
!         class(build), intent(inout)  :: b
!         class(params), intent(inout) :: p
!         real, allocatable :: ssnr_avg(:), frc_avg(:)
!         integer           :: dim2sz, icls, pop, cnt, k
!         dim2sz = size(b%ssnr, dim=2)
!         allocate( ssnr_avg(dim2sz) )
!         ssnr_avg = 0.
!         cnt      = 0
!         do icls=1,p%ncls
!             pop = b%a%get_clspop(icls)
!             if( pop > 1 )then
!                 cnt = cnt+1
!                 ssnr_avg = ssnr_avg+b%ssnr(icls,:)
!             endif
!         end do
!         ssnr_avg = ssnr_avg/real(cnt)
!         frc_avg  = ssnr2fsc(ssnr_avg)
!         write(*,'(A)') '>>> AVERAGE FRC'
!         do k=1,dim2sz
!             print *, 'RES/E{FRC}: ', b%img%get_lp(k), frc_avg(k)
!         end do
!         deallocate(ssnr_avg,frc_avg)
!     end subroutine
    
    subroutine prime2D_init_sums( b, p )
        class(build), intent(inout)   :: b
        class(params), intent(inout)  :: p
        integer :: icls
        do icls=1,p%ncls
            b%refs(icls)      = 0.
!             b%cavgs(icls)     = 0.
            b%ctfsqsums(icls) = cmplx(0.,0.)
        end do
!         if( allocated(b%spec_noise2D) ) b%spec_noise2D = 0.
    end subroutine
    
    subroutine prime2D_assemble_sums( b, p, mul )
        class(build),  intent(inout)  :: b
        class(params),  intent(inout) :: p
        real, optional, intent(in)    :: mul
        type(ori) :: orientation
        integer   :: icls, iptcl
        write(*,'(a)') '>>> ASSEMBLING CLASS SUMS'
        call prime2D_init_sums( b, p )
        do iptcl=1,p%nptcls
            call progress(iptcl,p%nptcls)
            orientation = b%a%get_ori(iptcl)
            if( nint(orientation%get('state')) > 0 )then
                icls = nint(orientation%get('class'))
                call b%img%read(p%stk, iptcl)
                call wiener_restore2D_online(b%img, orientation,&
                b%tfun, p%ctfmode, p%ctf, b%refs(icls), mul)
                call assemble_ctfsqsum_online(b%img, orientation,&
                b%tfun, p%ctfmode, p%ctf, b%ctfsqsums(icls))             
            endif
        end do
        call prime2D_norm_sums( b, p )
    end subroutine
    
    subroutine prime2D_assemble_sums_from_parts( b, p )
        use simple_jiffys, only: file_exists
        class(build), intent(inout)   :: b
        class(params), intent(inout)  :: p
        character(len=:), allocatable :: fname_refs, fname_ctfsqsums, fname_noisespecs
        real, allocatable             :: noise_specs(:,:)
        integer                       :: ipart, icls, dim1sz, dim2sz
        call prime2D_init_sums( b, p )
        do ipart=1,p%npart
            allocate(fname_refs,       source='refs_part'//int2str_pad(ipart,p%numlen)//p%ext)
            allocate(fname_ctfsqsums,  source='ctfsqsums_part'//int2str_pad(ipart,p%numlen)//p%ext)
            allocate(fname_noisespecs, source='noisespecs_part'//int2str_pad(ipart,p%numlen)//'.bin')
            ! read & sum partial class averages
            if( file_exists(fname_refs) )then
                do icls=1,p%ncls
                    call b%img%read(fname_refs, icls)
                    ! add subaverage to class
                    call b%refs(icls)%add(b%img)
                end do
            else
                write(*,*) 'File does not exists: ', trim(fname_refs)
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
!             ! read & sum partial noise power spectra
!             if( file_exists(fname_noisespecs) )then
!                 noise_specs = file2arr2D(fname_noisespecs)
!                 dim1sz = size(noise_specs, dim=1)
!                 dim2sz = size(noise_specs, dim=2)
!                 if( .not. allocated(b%spec_noise2D) )then
!                     allocate(b%spec_noise2D(dim1sz,dim2sz), source=noise_specs)
!                 else
!                     b%spec_noise2D = b%spec_noise2D+noise_specs
!                 endif
!                 deallocate(noise_specs)
!             else
!                 write(*,*) 'File does not exists: ', trim(fname_noisespecs)
!                 stop 'In: simple_hadamard2D_matcher :: prime2D_assemble'
!             endif
            deallocate(fname_refs, fname_ctfsqsums, fname_noisespecs)
        end do
        call prime2D_norm_sums( b, p )
    end subroutine
    
    subroutine prime2D_write_partial_sums( b, p )
        class(build), intent(inout)   :: b
        class(params), intent(inout)  :: p
        character(len=:), allocatable :: fname_refs, fname_ctfsqsums, fname_noisespecs
        real, allocatable             :: noise_specs(:,:)
        integer                       :: icls
        allocate(fname_refs,       source='refs_part'//int2str_pad(p%part,p%numlen)//p%ext)
        allocate(fname_ctfsqsums,  source='ctfsqsums_part'//int2str_pad(p%part,p%numlen)//p%ext)
        allocate(fname_noisespecs, source='noisespecs_part'//int2str_pad(p%part,p%numlen)//'.bin')
        do icls=1,p%ncls
            call b%refs(icls)%write(fname_refs, icls)
            call b%ctfsqsums(icls)%write(fname_ctfsqsums, icls)
        end do
!         call arr2D2file(b%spec_noise2D, fname_noisespecs)
        deallocate(fname_refs, fname_ctfsqsums, fname_noisespecs) 
    end subroutine

    subroutine prime2D_norm_sums( b, p )
        class(build),  intent(inout) :: b
        class(params), intent(inout) :: p
        integer :: icls, pop
        do icls=1,p%ncls
            pop = b%a%get_clspop(icls)
            if( pop > 1 )then
                
!                 call b%refs(icls)%div(real(pop))
!                 b%cavgs(icls) = b%refs(icls)
!                 call b%cavgs(icls)%fwd_ft
!                 call b%cavgs(icls)%ctf_dens_correct(b%ctfsqsums(icls))
!                 call b%cavgs(icls)%bwd_ft
!                 call b%cavgs(icls)%norm
!                 call b%refs(icls)%norm
                call b%refs(icls)%fwd_ft
                call b%refs(icls)%ctf_dens_correct(b%ctfsqsums(icls))
                call b%refs(icls)%bwd_ft
!                 if( allocated(b%spec_noise2D) )&
!                 &b%spec_noise2D(icls,:) = b%spec_noise2D(icls,:)/real(pop)
            endif
        end do
    end subroutine
    
    subroutine prime2D_write_sums( b, p, which_iter )
        class(build),      intent(inout) :: b
        class(params),     intent(inout) :: p
        integer, optional, intent(in)    :: which_iter
        integer                          :: icls, pop
        character(len=:), allocatable    :: fname_ctfsqsums, fname_cavgs, fname_noisespecs
        if( present(which_iter) )then
            if( which_iter <= 0 )then
                p%refs = 'refs'//p%ext
                allocate(fname_cavgs, source='cavgs'//p%ext)
            else
                p%refs = 'refs_iter'//int2str_pad(which_iter,3)//p%ext
                allocate(fname_cavgs, source='cavgs_iter'//int2str_pad(which_iter,3)//p%ext)
            endif
        else
             p%refs = 'startrefs'//p%ext
             allocate(fname_cavgs, source='startcavgs'//p%ext)
        endif
        allocate(fname_ctfsqsums,  source='ctfsqsums'//p%ext)
        allocate(fname_noisespecs, source='noisespecs'//'.bin')
        ! write to disk
        do icls=1,p%ncls
            call b%refs(icls)%write(p%refs, icls)
!             call b%cavgs(icls)%write(fname_cavgs, icls)
            call b%ctfsqsums(icls)%write(fname_ctfsqsums, icls)
        end do        
!         if( allocated(b%spec_noise2D) )then
!             call arr2D2file(b%spec_noise2D, fname_noisespecs)
!         endif
        deallocate(fname_ctfsqsums,fname_cavgs, fname_noisespecs)
    end subroutine
    
    !>  \brief  prepares the polarft corrcalc object for search
    subroutine preppftcc4align( b, p )
        use simple_image,  only: image
        class(build), intent(inout)  :: b
        class(params), intent(inout) :: p
        real, allocatable            :: sqrtssnr(:)
        type(image)                  :: ref
        integer :: cnt, iptcl, icls, sz
        write(*,'(A)') '>>> BUILDING PRIME2D SEARCH ENGINE'
        ! must be done here since constants in p are dynamically set
        call primesrch2D%new(p)
        write(*,'(A)') '>>> BUILDING DATA STRUCTURE FOR POLARFT CORRELATION CALCULATION'
        
!         call pftcc%new(p%ncls, [p%fromp,p%top], [p%box,p%box,1], p%kfromto, p%ring2, 'no')

        call pftcc%new(p%ncls, [p%fromp,p%top], [p%box,p%box,1], p%kfromto, p%ring2, p%ctf)
        ! PREPARATION OF REFERENCES IN PFTCC
        ! read references and transform into polar coordinates
        write(*,'(A)') '>>> BUILDING REFERENCES'
        sz = b%img%get_lfny(1)
        do icls=1,p%ncls
            call progress(icls, p%ncls)
            ref = b%refs(icls)
            ! apply a soft-edged mask
            if( defined_cmd_arg('inner') )then
                call ref%mask(p%msk, 'soft', inner=p%inner, width=p%width)
            else 
                call ref%mask(p%msk, 'soft')
            endif
            ! move to Fourier space
            call ref%fwd_ft
!             if( test )then
!                 ! scale amplitudes according to SSNR
!                 allocate(sqrtssnr(sz))
!                 sqrtssnr = 0.
!                 where( b%ssnr(icls,:) > 0. )
!                     sqrtssnr = sqrt(b%ssnr(icls,:))
!                 end where
!                 call ref%shellnorm
!                 call ref%apply_filter(sqrtssnr)
!                 deallocate(sqrtssnr)
!             endif
            ! transfer to polar coordinates
            call b%proj%img2polarft(icls,ref, pftcc, isptcl=.false.)
            call ref%kill
        end do
        ! PREPARATION OF PARTICLES IN PFTCC
        ! read particle images and create polar projections
        write(*,'(A)') '>>> BUILDING PARTICLES'
        cnt = 0
        do iptcl=p%fromp,p%top
            cnt = cnt+1
            call progress(cnt, p%top-p%fromp+1)
            if( defined_cmd_arg('part') )then
                call b%img%read(p%stk_part, cnt)
            else
                call b%img%read(p%stk, iptcl)
            endif
            call prepimg4align2D(b, p, iptcl)
            ! transfer to polar coordinates
            call b%proj%img2polarft(iptcl, b%img, pftcc)
        end do
        if( debug ) write(*,*) '*** hadamard2D_matcher ***: finished preppftcc4align'
    end subroutine
    
    subroutine prepimg4align2D( b, p, iptcl )
        class(build), intent(inout)  :: b
        class(params), intent(inout) :: p
        integer, intent(in)          :: iptcl
        real                         :: x, y, dfx, dfy, angast
        integer                      :: icls, sz
        real, allocatable            :: filter(:)
        if( p%lxfel )then
            ! nothing to do 4 now
            return
        else
            ! parse ori
            x    = b%a%get(iptcl, 'x')
            y    = b%a%get(iptcl, 'y')
            icls = nint(b%a%get(iptcl, 'class'))          
            ! move to Fourier space
            call b%img%fwd_ft
!             if( test )then
!                 ! scale amplitudes according to SSNR
!                 sz = b%img%get_lfny(1)
!                 allocate(filter(sz))
!                 filter = sqrt(b%ssnr(icls,:)+1.0)
!                 call b%img%shellnorm
!                 call b%img%apply_filter(filter)
!                 deallocate(filter)
!             endif
            ! deal with CTF
            select case(p%ctf)
                case('mul')  ! images have been multiplied with the CTF, no CTF-dependent weighting of the correlations
                    stop 'CTF multiplied particle images are not supported; simple_hadamard2D_matcher :: prepimg4align2D' 
                case('no')   ! do nothing
                case('yes')  ! do nothing
                case('flip') ! flip back
                    ! set CTF parameters
                    select case(p%ctfmode)
                        case('astig') ! astigmatic CTF
                            dfx    = b%a%get(iptcl,'dfx')
                            dfy    = b%a%get(iptcl,'dfy')
                            angast = b%a%get(iptcl,'angast')
                        case('noastig') ! non-astigmatic CTF
                            dfx    = b%a%get(iptcl,'dfx')
                            dfy    = dfx
                            angast = 0.
                        case DEFAULT
                            write(*,*) 'Unsupported p%ctfmode: ', trim(p%ctfmode)
                            stop 'simple_hadamard2D_matcher :: prepimg4align2D'
                    end select
                    ! modify the image
                    call b%tfun%apply(b%img, dfx, 'flip ', dfy, angast)
                case DEFAULT
                    stop 'Unsupported ctf mode; prepimg4align; simple_hadamard2D_matcher :: prepimg4align2D'
            end select
            ! shift image to rotational origin
            if( abs(x) > SHTHRESH .or. abs(y) > SHTHRESH ) call b%img%shift(-x, -y)
            ! apply a soft-edged mask
            call b%img%bwd_ft
            ! apply a soft-edged mask
            if( defined_cmd_arg('inner') )then
                call b%img%mask(p%msk, 'soft', inner=p%inner, width=p%width)
            else 
                call b%img%mask(p%msk, 'soft')
            endif            
            ! return in Fourier space
            call b%img%fwd_ft
        endif
        if( debug ) write(*,*) '*** simple_hadamard2D_matcher ***: finished prepimg4align2D'
    end subroutine
    
end module simple_hadamard2D_matcher
