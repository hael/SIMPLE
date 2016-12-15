module simple_unblur
use simple_image,  only: image
use simple_params, only: params
use simple_jiffys   ! singleton
use simple_filterer ! singleton
use simple_defs     ! singleton
implicit none

public :: unblur_movie, unblur_calc_sums, unblur_calc_sums_tomo, test_unblur
private

type(image), allocatable  :: movie_frames(:)            !< movie frames
type(image), allocatable  :: shmovie_frames(:)          !< shifted movie frames
type(image)               :: movie_sum_global           !< global movie sum for refinement
real, allocatable         :: corrmat(:,:)               !< matrix of correlations (to solve the exclusion problem)
real, allocatable         :: corrs(:)                   !< per-frame correlations
real, allocatable         :: frameweights(:)            !< array of frameweights
real, allocatable         :: frameweights_saved(:)      !< array of frameweights
real, allocatable         :: opt_shifts(:,:)            !< optimal shifts identified
real, allocatable         :: opt_shifts_saved(:,:)      !< optimal shifts for local opt saved
real, allocatable         :: acc_doses(:)               !< accumulated doses
integer                   :: nframes        = 0         !< number of frames
integer                   :: fixed_frame    = 0         !< fixed frame of reference (0,0)
integer                   :: ldim(3)        = [0,0,0]   !< logical dimension of frame
integer                   :: ldim_scaled(3) = [0,0,0]   !< shrunken logical dimension of frame
real                      :: maxshift       = 0.        !< maximum halfwidth shift
real                      :: hp             = 0.        !< high-pass limit
real                      :: lp             = 0.        !< low-pass limit
real                      :: resstep        = 0.        !< resolution step size (in Å)
real                      :: smpd           = 0.        !< sampling distance
real                      :: corr_saved     = 0.        !< opt corr for local opt saved
real                      :: kV             = 300.      !< acceleration voltage
real                      :: dose_rate      = 0.        !< dose rate
logical                   :: do_dose_weight = .false.   !< dose weight or not
logical                   :: doscale        = .false.   !< scale or not
logical                   :: doprint        = .true.    !< print out correlations
logical                   :: debug          = .false.   !< debug or not
logical                   :: existence      = .false.   !< to indicate existence

integer, parameter :: MITSREF    = 30 !< max nr iterations of refinement optimisation
integer, parameter :: SMALLSHIFT = 2. !< small initial shift to blur out fixed pattern noise

contains
    
    subroutine unblur_movie( movie_stack_fname, p, corr )
        use simple_oris,     only: oris
        use simple_jiffys,   only: int2str
        use simple_ft_shsrch ! singleton
        use simple_rnd,      only: ran3
        use simple_stat,     only: corrs2weights, moment
        character(len=*), intent(in)    :: movie_stack_fname
        class(params),    intent(inout) :: p
        real,             intent(out)   :: corr
        real    :: ave, sdev, var, minw, maxw
        real    :: cxy(3), lims(2,2), corr_prev, frac_improved, corrfrac
        integer :: iframe, iter, nimproved, ires, updateres, i
        logical :: didsave, didupdateres, err
        ! initialise
        call unblur_init(movie_stack_fname, p)
        ! make search object ready
        lims(:,1) = -maxshift
        lims(:,2) =  maxshift
        call ft_shsrch_init(movie_sum_global, movie_frames(1), lims, lp, hp)
        ! initialise with small random shifts (to average out dead/hot pixels)
        do iframe=1,nframes
            opt_shifts(iframe,1) = ran3()*2.*SMALLSHIFT-SMALLSHIFT
            opt_shifts(iframe,2) = ran3()*2.*SMALLSHIFT-SMALLSHIFT
        end do
        ! generate movie sum for refinement
        call shift_frames(opt_shifts)
        call calc_corrmat
        call corrmat2weights ! this should remove any possible extreme outliers
        call wsum_movie_frames
        ! calc avg corr to weighted avg
        call calc_corrs
        corr = sum(corrs)/real(nframes)
        if( doprint ) write(*,'(a)') '>>> WEIGHTED AVERAGE-BASED REFINEMENT'
        if( doprint ) write(*,'(a,7x,f7.4)') '>>> INITIAL CORRELATION:', corr
        iter       = 0
        corr_saved = -1.
        didsave    = .false.
        updateres  = 0
        do i=1,MITSREF
            iter = iter+1
            nimproved = 0
            do iframe=1,nframes
                ! subtract the movie frame being aligned to reduce bias
                call subtract_movie_frame( iframe )
                call ft_shsrch_reset_ptrs(movie_sum_global, movie_frames(iframe))
                cxy = ft_shsrch_minimize(corrs(iframe), opt_shifts(iframe,:))
                if( cxy(1) > corrs(iframe) ) nimproved = nimproved+1
                opt_shifts(iframe,:) = cxy(2:3)
                corrs(iframe)        = cxy(1)
                ! add the subtracted movie frame back to the weighted sum
                call add_movie_frame( iframe )
            end do            
            frac_improved = real(nimproved)/real(nframes)*100.
            if( doprint ) write(*,'(a,1x,f4.0)') 'This % of frames improved their alignment: ', frac_improved
            call center_shifts(opt_shifts)
            frameweights = corrs2weights(corrs)
            call shift_frames(opt_shifts)
            call wsum_movie_frames
            corr_prev = corr
            corr = sum(corrs)/real(nframes)
            if( corr >= corr_saved )then ! save the local optimum
                corr_saved         = corr
                opt_shifts_saved   = opt_shifts
                frameweights_saved = frameweights
                didsave = .true.
            endif
            if( doprint )  write(*,'(a,7x,f7.4)') '>>> OPTIMAL CORRELATION:', corr
            corrfrac = corr_prev/corr
            didupdateres = .false.
            select case(updateres)
                case(0)
                    call update_res( 0.96, 70., updateres )
                case(1)
                    call update_res( 0.97, 60., updateres )
                case(2)
                    call update_res( 0.98, 50., updateres )
                case DEFAULT
                    ! nothing to do
            end select
            if( updateres > 2 .and. .not. didupdateres )then ! at least one iteration with new lim
                if( nimproved == 0 .and. i > 2 )  exit
                if( i > 10 .and. corrfrac > 0.9999 ) exit
            endif
        end do
        ! put the best local optimum back
        corr         = corr_saved
        opt_shifts   = opt_shifts_saved
        frameweights = frameweights_saved        
        call shift_frames(opt_shifts)
        ! print
        if( doprint ) write(*,'(a,7x,f7.4)') '>>> OPTIMAL CORRELATION:', corr
        call moment(frameweights, ave, sdev, var, err)
        minw = minval(frameweights)
        maxw = maxval(frameweights)
        if( doprint ) write(*,'(a,7x,f7.4)') '>>> AVERAGE WEIGHT     :', ave
        if( doprint ) write(*,'(a,7x,f7.4)') '>>> SDEV OF WEIGHTS    :', sdev
        if( doprint ) write(*,'(a,7x,f7.4)') '>>> MIN WEIGHT         :', minw
        if( doprint ) write(*,'(a,7x,f7.4)') '>>> MAX WEIGHT         :', maxw
        
        contains
            
            subroutine update_res( thres_corrfrac, thres_frac_improved, which_update )
                real,    intent(in) :: thres_corrfrac, thres_frac_improved
                integer, intent(in) :: which_update
                if( corrfrac > thres_corrfrac .and. frac_improved <= thres_frac_improved&
                .and. updateres == which_update )then
                    lp = lp - resstep
                    if( doprint )  write(*,'(a,1x,f7.4)') '>>> LOW-PASS LIMIT UPDATED TO:', lp
                    ! need to update correlation values
                    call calc_corrs
                    ! need to indicate that we updated resolution limit
                    updateres  = updateres + 1
                    ! need to destroy all previous knowledge about correlations
                    corr       = sum(corrs)/real(nframes)
                    corr_prev  = corr
                    corr_saved = corr
                    ! indicate that reslim was updated
                    didupdateres = .true.
                endif
            end subroutine update_res
                    
    end subroutine unblur_movie

    subroutine unblur_calc_sums( movie_sum, movie_sum_corrected, movie_sum_ctf )
        type(image), intent(out) :: movie_sum, movie_sum_corrected, movie_sum_ctf
        integer :: iframe
        ! calculate the sum for CTF estimation
        call sum_movie_frames
        movie_sum_ctf = movie_sum_global
        call movie_sum_ctf%bwd_ft
        ! re-calculate the weighted sum
        if( do_dose_weight ) call dose_weight_movie_frames
        call wsum_movie_frames
        movie_sum_corrected = movie_sum_global
        call movie_sum_corrected%bwd_ft
        ! generate straight integrated movie frame for comparison
        do iframe=1,nframes
            shmovie_frames(iframe) = movie_frames(iframe)
        end do
        call sum_movie_frames
        movie_sum = movie_sum_global
        call movie_sum%bwd_ft
    end subroutine unblur_calc_sums

    subroutine unblur_calc_sums_tomo( frame_counter, time_per_frame, movie_sum, movie_sum_corrected, movie_sum_ctf )
        integer,     intent(inout) :: frame_counter
        real,        intent(in)    :: time_per_frame
        type(image), intent(out)   :: movie_sum, movie_sum_corrected, movie_sum_ctf
        integer :: iframe
        ! calculate the sum for CTF estimation
        call sum_movie_frames
        movie_sum_ctf = movie_sum_global
        call movie_sum_ctf%bwd_ft
        ! re-calculate the weighted sum
        call dose_weight_movie_frames_tomo( frame_counter, time_per_frame )
        call wsum_movie_frames
        movie_sum_corrected = movie_sum_global
        call movie_sum_corrected%bwd_ft
        ! generate straight integrated movie frame for comparison
        do iframe=1,nframes
            shmovie_frames(iframe) = movie_frames(iframe)
        end do
        call sum_movie_frames
        movie_sum = movie_sum_global
        call movie_sum%bwd_ft
    end subroutine unblur_calc_sums_tomo
     
    subroutine unblur_init( movie_stack_fname, p )
        character(len=*), intent(in)    :: movie_stack_fname
        class(params),    intent(inout) :: p
        real        :: moldiam, dimo4
        integer     :: alloc_stat, iframe
        type(image) :: frame_tmp
        real        :: time_per_frame, current_time
        call unblur_kill  
        ! GET NUMBER OF FRAMES & DIM FROM STACK
        call find_ldim_nptcls(movie_stack_fname, ldim, nframes, endconv=endconv)
        if( debug ) write(*,*) 'logical dimension: ', ldim
        ldim(3) = 1 ! to correct for the stupide 3:d dim of mrc stacks
        if( p%scale < 0.99 )then
            ldim_scaled(1) = nint(real(ldim(1))*p%scale)
            ldim_scaled(2) = nint(real(ldim(2))*p%scale)
            ldim_scaled(3) = 1
            doscale        = .true.
        else
            ldim_scaled = ldim
            doscale     = .false.
        endif
        ! SET SAMPLING DISTANCE
        smpd = p%smpd/p%scale
        if( debug ) write(*,*) 'logical dimension of frame: ',        ldim
        if( debug ) write(*,*) 'scaled logical dimension of frame: ', ldim_scaled
        if( debug ) write(*,*) 'number of frames: ',                  nframes
        ! ALLOCATE
        allocate( movie_frames(nframes), shmovie_frames(nframes),&
        corrs(nframes), opt_shifts(nframes,2), opt_shifts_saved(nframes,2),&
        corrmat(nframes,nframes), frameweights(nframes),&
        frameweights_saved(nframes), stat=alloc_stat )
        call alloc_err('unblur_init; simple_unblur', alloc_stat)
        corrmat = 0.
        corrs   = 0.
        ! read and FT frames
        do iframe=1,nframes
            call frame_tmp%new(ldim, smpd)
            call frame_tmp%read(movie_stack_fname,iframe)
            call frame_tmp%fwd_ft
            call movie_frames(iframe)%new(ldim_scaled, smpd)
            call movie_frames(iframe)%set_ft(.true.)
            call frame_tmp%clip(movie_frames(iframe))
            call shmovie_frames(iframe)%copy(movie_frames(iframe))
        end do
        maxshift = p%trs/p%scale
        ! set fixed frame (all others are shifted by reference to this at 0,0)
        fixed_frame = nint(real(nframes)/2.)
        ! set reslims
        dimo4     = (real(minval(ldim_scaled(1:2)))*p%smpd)/4.
        moldiam   = 0.7*real(p%box)*p%smpd
        hp        = min(dimo4,2000.)
        lp        = p%lpstart
        resstep   = (p%lpstart-p%lpstop)/3.
        ! check if we are doing dose weighting
        if( p%l_dose_weight )then
            do_dose_weight = .true.
            allocate( acc_doses(nframes), stat=alloc_stat )
            call alloc_err('unblur_init; simple_unblur, 2', alloc_stat)
            kV = p%kV
            time_per_frame = p%exp_time/real(nframes)           ! unit: s
            dose_rate      = p%dose_rate
            do iframe=1,nframes
                current_time      = real(iframe)*time_per_frame ! unit: s
                acc_doses(iframe) = dose_rate*current_time      ! unit: e/A2/s * s = e/A2
            end do
        endif
        existence = .true.
        if( debug ) write(*,*) 'unblur_init, done'
    end subroutine unblur_init

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
    end subroutine center_shifts
    
    subroutine shift_frames( shifts )
        real, intent(in) :: shifts(nframes,2)
        integer :: iframe
        do iframe=1,nframes
            call movie_frames(iframe)%shift(-shifts(iframe,1), -shifts(iframe,2),&
            lp_dyn=lp, imgout=shmovie_frames(iframe))
        end do
    end subroutine shift_frames
    
    subroutine calc_corrmat
        integer :: iframe, jframe
        corrmat = 1. ! diagonal elements are 1
        do iframe=1,nframes-1
            do jframe=iframe+1,nframes
                corrmat(iframe,jframe) = shmovie_frames(iframe)%corr(shmovie_frames(jframe), lp, hp)
                corrmat(jframe,iframe) = corrmat(iframe,jframe)
            end do
        end do
    end subroutine calc_corrmat

    subroutine calc_corrs
        integer :: iframe
        real    :: old_corr
        do iframe=1,nframes
            ! subtract the movie frame being correlated to reduce bias
            call subtract_movie_frame( iframe )
            corrs(iframe) = movie_sum_global%corr(shmovie_frames(iframe), lp, hp)
            ! add the subtracted movie frame back to the weighted sum
            call add_movie_frame( iframe )
        end do
    end subroutine calc_corrs

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
    end subroutine corrmat2weights

    subroutine dose_weight_movie_frames
        integer           :: iframe
        real, allocatable :: filter(:)
        do iframe=1,nframes
            filter = acc_dose2filter(shmovie_frames(iframe), acc_doses(iframe), kV)
            call shmovie_frames(iframe)%apply_filter(filter)
            deallocate(filter)
        end do
    end subroutine dose_weight_movie_frames

    subroutine dose_weight_movie_frames_tomo( frame_counter, time_per_frame )
        integer, intent(inout) :: frame_counter
        real,    intent(in)    :: time_per_frame
        real, allocatable      :: filter(:)
        integer                :: iframe
        real                   :: current_time, acc_dose
        do iframe=1,nframes
            frame_counter = frame_counter + 1
            current_time  = real(frame_counter)*time_per_frame ! unit: s
            acc_dose      = dose_rate*current_time             ! unit e/A2
            filter        = acc_dose2filter(shmovie_frames(iframe), acc_dose, kV)
            call shmovie_frames(iframe)%apply_filter(filter)
            deallocate(filter)
        end do
    end subroutine dose_weight_movie_frames_tomo

    subroutine sum_movie_frames
        integer :: iframe
        integer :: ldim_here(3)
        real    :: w
        ldim_here = shmovie_frames(1)%get_ldim()
        call movie_sum_global%new(ldim_here, smpd)
        call movie_sum_global%set_ft(.true.)
        w = 1./real(nframes)
        do iframe=1,nframes
            call movie_sum_global%add(shmovie_frames(iframe), w=w)
        end do
    end subroutine sum_movie_frames
    
    subroutine wsum_movie_frames
        integer :: iframe
        integer :: ldim_here(3)
        ldim_here = shmovie_frames(1)%get_ldim()
        call movie_sum_global%new(ldim_here, smpd)
        call movie_sum_global%set_ft(.true.)
        do iframe=1,nframes
            if( frameweights(iframe) > 0. )&
            &call movie_sum_global%add(shmovie_frames(iframe), w=frameweights(iframe))
        end do
    end subroutine wsum_movie_frames

    subroutine add_movie_frame( iframe, w )
        integer,        intent(in) :: iframe
        real, optional, intent(in) :: w
        real :: ww
        ww = frameweights(iframe)
        if( present(w) ) ww = w
        if( frameweights(iframe) > 0. )&
        &call movie_sum_global%add(shmovie_frames(iframe), w=ww)
    end subroutine add_movie_frame

    subroutine subtract_movie_frame( iframe )
        integer, intent(in) :: iframe
        if( frameweights(iframe) > 0. )&
        &call movie_sum_global%subtr(shmovie_frames(iframe), w=frameweights(iframe))
    end subroutine subtract_movie_frame
    
    subroutine test_unblur
        use simple_oris,    only: oris
        use simple_cmdline, only: cmdline
        real          :: shifts(7,2)
        type(image)   :: squares(7), straight_sum, corrected , sum4ctf
        type(params)  :: p_here 
        type(cmdline) :: cline
        real          :: corr
        integer       :: i
        shifts(1,1) = -3.
        shifts(1,2) = -3.
        shifts(2,1) = -2.
        shifts(2,2) = -2.
        shifts(3,1) = -1.
        shifts(3,2) = -1.
        shifts(4,1) =  0.
        shifts(4,2) =  0.
        shifts(5,1) =  1.
        shifts(5,2) =  1.
        shifts(6,1) =  2.
        shifts(6,2) =  2.
        shifts(7,1) =  3.
        shifts(7,2) =  3.
        do i=1,7
            call squares(i)%new([100,100,1], 2.)
            call squares(i)%square(30)
            call squares(i)%shift(shifts(i,1),shifts(i,2))
            call squares(i)%write('sugarcubes.mrc', i)
        end do
        p_here         = params(cline)
        p_here%smpd    = 2.0
        p_here%trs     = 3.5
        p_here%box     = 100
        p_here%lpstart = 10.
        p_here%lpstop  = 4.0
        call unblur_movie('sugarcubes.mrc', p_here, corr)
        call unblur_calc_sums(straight_sum, corrected, sum4ctf)
        if( sum(abs(opt_shifts-shifts)) < 0.5 )then
            write(*,'(a)') 'SIMPLE_UNBLUR_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
        else
            write(*,'(a)') 'SIMPLE_UNBLUR_UNIT_TEST FAILED :-('
        endif
        do i=1,7
            call squares(i)%kill
        end do
        call straight_sum%kill
        call corrected%kill
        call sum4ctf%kill
    end subroutine test_unblur

    subroutine unblur_kill
        integer :: iframe
        if( existence )then
            do iframe=1,nframes
                call movie_frames(iframe)%kill
                call shmovie_frames(iframe)%kill
            end do
            call movie_sum_global%kill
            deallocate( movie_frames, shmovie_frames, frameweights,&
            frameweights_saved, corrs, corrmat, opt_shifts, opt_shifts_saved )
            if( allocated(acc_doses) ) deallocate(acc_doses)
            existence = .false.
        endif
    end subroutine unblur_kill

end module simple_unblur