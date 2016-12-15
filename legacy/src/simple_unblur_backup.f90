module simple_unblur
use simple_opt_factory, only: opt_factory
use simple_opt_spec,    only: opt_spec
use simple_optimizer,   only: optimizer
use simple_image,       only: image
use simple_imgfile,     only: imgfile
use simple_params,      only: params, find_ldim
use simple_jiffys,      only: alloc_err
use simple_defs         ! singleton
implicit none

public :: unblur_step1, unblur_step2
private

type(opt_factory)         :: ofac                     !< optimizer factory
type(opt_spec)            :: ospec                    !< optimizer specification object
class(optimizer), pointer :: nlopt=>null()            !< pointer to nonlinear optimizer
type(imgfile)             :: imgstk                   !< object to deal with image stacks
type(image), allocatable  :: movie_frames(:)          !< movie frames
type(image), allocatable  :: shmovie_frames(:)        !< shifted movie frames
type(image)               :: movie_sum_global         !< global movie sum for refinement
real, allocatable         :: lims(:,:)                !< limits for optimization (enforced by the barrier method)
real, allocatable         :: corrmat(:,:)             !< matrix of correlations (to solve the exclusion problem)
real, allocatable         :: corrs(:)                 !< per-frame correlations
real, allocatable         :: frameweights(:)          !< array of frameweights
real, allocatable         :: opt_shifts(:,:)          !< optimal shifts identified
integer                   :: ndim         = 0         !< dimensionality of the problem
integer                   :: npairs       = 0         !< number of frame pairs
integer                   :: nframes      = 0         !< number of frames
integer                   :: fixed_frame  = 0         !< fixed frame of reference (0,0)
integer                   :: ldim(3)      = [0,0,0]   !< logical dimension of frame
integer                   :: frame_global = 0         !< global frame index (used in avg-based refinement)
integer                   :: mits         = 0         !< maximum number of iterations
real                      :: maxshift     = 0.        !< maximum halfwidth shift
real                      :: hp           = 0.        !< high-pass limit
real                      :: lp           = 0.        !< low-pass limit
real                      :: smpd         = 0.        !< sampling distance
logical                   :: refine       = .false.   !< refinement flag
logical                   :: doprint      = .true.    !< print out correlations
logical                   :: debug        = .false.   !< debug or not
logical                   :: existence    = .false.   !< to indicate existence

integer, parameter :: MITSREF   = 50 !< max nr iterations of refinement optimisation
integer, parameter :: NRESTARTS = 3  !< nr restarts 4 initial simplex opt with randomised bounds

interface shift_frames
    module procedure shift_frames_1
    module procedure shift_frames_2
end interface

contains

    subroutine init_step1( movie_stack_fname, p, opt, nrestarts )
        character(len=*), intent(in) :: movie_stack_fname
        class(params), intent(inout) :: p
        character(len=*), intent(in) :: opt
        integer, intent(in)          :: nrestarts
        real                         :: moldiam, dimo4
        real, allocatable            :: vec(:)
        integer                      :: i, alloc_stat
        call unblur_kill  
        ! SET SAMPLING DISTANCE
        smpd = p%smpd
        ! GET NUMBER OF FRAMES & DIM FROM STACK
        call imgstk%open(movie_stack_fname)
        nframes = imgstk%getStackSz()
        call imgstk%close
        if( debug ) write(*,*) 'number of frames: ', nframes
        ldim = find_ldim(movie_stack_fname)
        ldim(3) = 1 ! to correct for the stupide 3:d dim of mrc stacks
        if( debug ) write(*,*) 'logical dimension: ', ldim
        if( debug ) write(*,*) 'logical dimension of frame: ', ldim
        ndim   = 2*(nframes-1) ! because the center frame is kept fixed
        npairs = ((nframes)*(nframes-1))/2 
        mits   = 10*ndim ! maximum number of iterations for the refinement optimiser
        if( debug ) print *, 'ndim: ', ndim
        ! ALLOCATE
        allocate( vec(ndim), movie_frames(nframes), shmovie_frames(nframes),&
        lims(ndim,2), corrmat(nframes,nframes), corrs(nframes),&
        opt_shifts(nframes,2), stat=alloc_stat )
        call alloc_err('init; simple_unblur', alloc_stat)
        vec     = 0.
        lims    = 0.
        corrmat = 0.
        corrs   = 0.        
        ! read and FT frames
        do i=1,nframes
            call movie_frames(i)%new(ldim, smpd)
            call movie_frames(i)%read(movie_stack_fname,i)
            call movie_frames(i)%fwd_ft
            shmovie_frames(i) = movie_frames(i)
        end do
        ! INITIALIZE
        ! set optimization limits
        do i=1,ndim
            lims(i,1) = -p%trs
            lims(i,2) = p%trs
        end do
        maxshift = p%trs
        ! set fixed frame (all others are shifted by reference to this at 0,0)
        fixed_frame = nint(real(nframes)/2.)
        ! SET RESOLUTION LIMITS
        dimo4   = real(minval(ldim(1:2)))/4.
        moldiam = 0.7*real(p%box)*p%smpd
        hp      = dimo4
        lp      = moldiam/4.
        ! MAKE THE OPTIMIZER READY
        call ospec%specify(opt, ndim, ftol=1e-4, gtol=1e-4, limits=lims, nrestarts=nrestarts)
        call ospec%set_costfun(costfun)
        if( debug ) write(*,'(a,1x,f7.4)') '>>> START CORRELATION: ', -costfun(vec, ndim)
        ! generate optimizer object with the factory
        call ofac%new(ospec, nlopt)
        deallocate(vec)
        existence = .true.
    end subroutine
    
    function get_shifts() result( shifts )
        real, allocatable :: shifts(:,:)
        shifts = vec2shifts( ospec%x, ndim ) 
    end function
    
    subroutine center_shifts( shifts )
        real, intent(inout) :: shifts(nframes,2)
        real    :: xsh, ysh
        integer :: iframe
        xsh = -shifts(fixed_frame,1)
        ysh = -shifts(fixed_frame,2)
        do iframe=1,nframes
            shifts(iframe,1) = shifts(iframe,1)+xsh
            shifts(iframe,2) = shifts(iframe,2)+ysh
            if( abs(shifts(iframe,1)) < 1e-6 ) shifts(iframe,1) = 0.
            if( abs(shifts(iframe,2)) < 1e-6 ) shifts(iframe,2) = 0.
        end do
    end subroutine
    
    subroutine unblur_step1( movie_stack_fname, p, oind, oout, corr, movie_sum )
        use simple_oris,   only: oris
        use simple_jiffys, only: int2str
        character(len=*), intent(in) :: movie_stack_fname
        class(params), intent(inout) :: p
        integer, intent(in)          :: oind
        class(oris), intent(inout)   :: oout
        real, intent(out)            :: corr
        type(image), intent(out)     :: movie_sum
        integer :: iframe
        call init_step1(movie_stack_fname, p, 'simplex', NRESTARTS)
!         corr = minimize()
        ospec%x = 0.
        opt_shifts = get_shifts()
        do iframe=1,size(opt_shifts,1)
            call oout%set(oind, 'x'//int2str(iframe), opt_shifts(iframe,1))
            call oout%set(oind, 'y'//int2str(iframe), opt_shifts(iframe,2))
        end do
        call shift_frames(opt_shifts)
        call calc_corrmat
        call corrmat2weights
        ! movie sum for refinement
        call wsum_movie_frames
        ! outputted movie sum (always in real space)
        movie_sum = movie_sum_global
        call movie_sum%bwd_ft 
    end subroutine
    
    subroutine unblur_step2( oind, oout, corr, movie_sum )
        use simple_oris,     only: oris
        use simple_jiffys,   only: int2str
        use simple_ft_shsrch ! singleton
        use simple_stat,     only: corrs2weights
        integer, intent(in)          :: oind
        class(oris), intent(inout)   :: oout
        real, intent(out)            :: corr
        type(image), intent(out)     :: movie_sum
        real    :: cxy(3), lims(2,2), corr_prev
        integer :: iframe, iter, i, nimproved
        ! make search object ready
        lims(:,1) = -maxshift
        lims(:,2) = maxshift
        call ft_shsrch_init(movie_sum_global, movie_frames(1), lims, lp, hp)
        ! calc avg corr to weighted avg
        call shift_frames(opt_shifts)
        call calc_corrs
        corr = sum(corrs)/real(nframes)
        if( doprint ) write(*,'(a,1x,f7.4)') '>>> INITIAL CORRELATION:', corr
        iter = 0
        do i=1,MITSREF
            iter = iter+1
            nimproved = 0
            do iframe=1,nframes
                call ft_shsrch_reset_ptrs(movie_sum_global, movie_frames(iframe))
                cxy = ft_shsrch_minimize(corrs(iframe), opt_shifts(iframe,:))
                if( cxy(1) > corrs(iframe) )  nimproved = nimproved+1
                opt_shifts(iframe,:) = cxy(2:3)
                corrs(iframe)        = cxy(1)
            end do
            if( doprint ) print *, 'This % of frames improved their alignment: ',&
            real(nimproved)/real(nframes)*100.
            call center_shifts(opt_shifts)
            frameweights = corrs2weights(corrs)
            call shift_frames(opt_shifts)
            call wsum_movie_frames
            corr_prev = corr
            corr = sum(corrs)/real(nframes)
            if( doprint )  write(*,'(a,1x,f7.4)') '>>> OPTIMAL CORRELATION:', corr
            if( (nimproved == 0 .and. i > 2) .or. corr_prev/corr > 0.9999 ) exit
        end do
        ! outputted shift parameters
        do iframe=1,size(opt_shifts,1)
            call oout%set(oind, 'x'//int2str(iframe), opt_shifts(iframe,1))
            call oout%set(oind, 'y'//int2str(iframe), opt_shifts(iframe,2))
        end do
        ! outputted movie sum (always in realspace)
        movie_sum = movie_sum_global
        call movie_sum%bwd_ft
    end subroutine
    
    function minimize( ) result( corr )
        real :: corr
        ospec%x = 0.
        call nlopt%minimize(ospec, corr)
        corr  = -corr
        if( debug ) write(*,'(a,1x,f7.4)') '>>> OPTIMAL CORRELATION:', corr
    end function

    function costfun( vec, D ) result( cost )
        integer, intent(in) :: D
        real, intent(in)    :: vec(D)
        real                :: corr, corrsum, cost
        integer             :: i, j, ncorrs
        real, allocatable   :: shifts(:,:)
        call shift_frames( vec, D, shifts )
        ! IMPLEMENT THE BARRIER METHOD
        if( any(abs(shifts(:,:)) > maxshift) )then
            deallocate(shifts)
            cost = 1.
            return
        endif
        ncorrs = (nframes*(nframes-1))/2
        corrsum = 0.
        do i=1,nframes-1
            do j=i+1,nframes
                corr = shmovie_frames(i)%corr(shmovie_frames(j), lp, hp)
                corrsum = corrsum+corr
            end do
        end do
        cost = -corrsum/real(ncorrs)
        deallocate(shifts)
    end function
    
    subroutine shift_frames_1( shifts )
        real, intent(in) :: shifts(nframes,2)
        integer :: frame
        do frame=1,nframes
            call movie_frames(frame)%shift(-shifts(frame,1), -shifts(frame,2),&
            lp_dyn=lp, imgout=shmovie_frames(frame))
        end do
    end subroutine
    
    subroutine shift_frames_2( vec, D, shifts_out )
        integer, intent(in)                      :: D
        real, intent(in)                         :: vec(D)
        real, allocatable, intent(out), optional :: shifts_out(:,:)
        real, allocatable                        :: shifts(:,:)
        integer                                  :: alloc_stat        
        shifts = vec2shifts( vec, D )
        call shift_frames( shifts )
        if( present(shifts_out) )then
            if( allocated(shifts_out) ) deallocate(shifts_out)
            allocate( shifts_out(nframes,2), source=shifts, stat=alloc_stat )
            call alloc_err('In: simple_unblur::shift_frames_3', alloc_stat)
        endif
        deallocate(shifts)        
    end subroutine
    
    subroutine calc_corrmat
        integer :: iframe, jframe
        corrmat = 1. ! diagonal elements are 1
        do iframe=1,nframes-1
            do jframe=iframe+1,nframes
                corrmat(iframe,jframe) = shmovie_frames(iframe)%corr(shmovie_frames(jframe), lp, hp)
                corrmat(jframe,iframe) = corrmat(iframe,jframe)
            end do
        end do
    end subroutine
    
    subroutine calc_corrs
        integer :: frame
        do frame=1,nframes
            corrs(frame) = movie_sum_global%corr(shmovie_frames(frame), lp, hp)
        end do
    end subroutine

    subroutine corrmat2weights
        use simple_stat, only: corrs2weights
        integer :: iframe, jframe
        do iframe=1,nframes
            corrs(iframe) = 0.
            do jframe=1,nframes
                if( jframe == iframe ) cycle
                corrs(iframe) = corrs(iframe)+corrmat(iframe,jframe)
            end do
            corrs(iframe) = corrs(iframe)/real(nframes-1)
        end do
        frameweights = corrs2weights(corrs)
    end subroutine
    
    subroutine wsum_movie_frames
        integer     :: frame
        call movie_sum_global%new(ldim, smpd)
        call movie_sum_global%set_ft(.true.)
        do frame=1,nframes
            if( frameweights(frame) > 0. )then
                call movie_sum_global%add(shmovie_frames(frame), w=frameweights(frame))
            endif
        end do
    end subroutine
    
    function vec2shifts( vec, D ) result( shifts )
        integer, intent(in) :: D
        real, intent(in)    :: vec(D)
        real, allocatable   :: shifts(:,:)
        integer :: frame, i, alloc_stat
        logical :: reached_fixed_frame
        allocate( shifts(nframes,2), stat=alloc_stat )
        call alloc_err('In: simple_unblur::vec2shifts', alloc_stat)
        frame = 0
        reached_fixed_frame = .false.
        do i=1,ndim+2,2
            frame = frame+1
            if( frame == fixed_frame )then
                shifts(frame,:) = 0.
                reached_fixed_frame = .true.
                cycle
            endif
            if( reached_fixed_frame )then
                shifts(frame,:) = [vec(i-2),vec(i-1)]
            else
                shifts(frame,:) = [vec(i),vec(i+1)]
            endif
            if( abs(shifts(frame,1)) < 1e-6 ) shifts(frame,1) = 0.
            if( abs(shifts(frame,2)) < 1e-6 ) shifts(frame,2) = 0.
        end do
        shifts = -shifts ! 2 fit shifting convention
    end function

    subroutine unblur_kill
        integer :: i
        if( existence )then
            call ofac%kill
            call ospec%kill
            nlopt => null()
            do i=1,nframes
                call movie_frames(i)%kill
                call shmovie_frames(i)%kill
            end do
            deallocate( movie_frames, shmovie_frames,&
            lims, corrmat, frameweights, corrs, opt_shifts )
            refine    = .false.
            existence = .false.
        endif
    end subroutine

end module simple_unblur