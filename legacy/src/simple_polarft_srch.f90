module simple_polarft_srch
use simple_defs        ! singleton
use simple_polarft,    only: polarft
use simple_oris,       only: oris
use simple_align_pair, only: align_pair
use simple_ran_tabu,   only: ran_tabu
use simple_jiffys,     only: alloc_err
implicit none

public :: polarft_srch_init, polarft_srch_align_ptcl, polarft_srch_get_ori, polarft_srch_get_ori_best,&
polarft_srch_get_state, polarft_srch_get_proj, polarft_srch_calc_weights, polarft_srch_get_class,&
polarft_srch_get_euler, polarft_srch_ang_sdev
private

type proj_dir
    real, allocatable    :: inpl_corrs(:)
    integer, allocatable :: inpl_inds(:)
    integer              :: proj=1, state=1
    real                 :: x=0., y=0.
end type

type polarft_ref
    class(polarft), pointer :: p=>null()
end type    

! CLASS PARAMETERS/VARIABLES
class(oris), pointer           :: o_refs=>null()        !< even reference orientations
type(polarft_ref), allocatable :: refs(:)               !< polar reference sections
integer, allocatable           :: ref_srch_order(:)     !< stochastic search order            
integer, allocatable           :: proj_space_inds(:)    !< projection space index array
real, allocatable              :: proj_space_corrs(:)   !< projection space correlations
type(proj_dir), allocatable    :: proj_space(:)         !< search space
integer, allocatable           :: refinds(:,:)          !< reference section indices
real, allocatable              :: probs(:)              !< multinomal distribution for the states
real, allocatable              :: probs_prev(:)         !< previous multinomal distribution for the states
type(align_pair)               :: ap                    !< align pair object
type(ran_tabu)                 :: rt                    !< random number generator
character(len=STDLEN)          :: diversify             !< search mode
integer                        :: nstates, nrefs, nobjs !< constants
integer                        :: kfromto(2)            !< resolution range
integer                        :: nsamp                 !< nr of rotations in polar image
integer                        :: npeaks                !< nr of peaks for soft projection matching
integer                        :: nbetter               !< nr of better orientations
integer                        :: nrefs_eval            !< nr of references evaluated
real                           :: corr_best_glob        !< globally best correlation
real                           :: corr_prev             !< previous best correlation
real                           :: rotxy(3)              !< global rotxy
real                           :: eps=0.5               !< learning rate for the state probabilities
integer                        :: best_ref              !< best reference found so far
logical                        :: debug=.false.         !< debugging mode or not
real                           :: xold=0.,yold=0.       !< old origin shift

contains

    ! CONSTRUCTOR

    !>  \brief  is a constructor
    subroutine polarft_srch_init( p, b, sh )
        use simple_build,  only: build  
        use simple_params, only: params
        class(params), intent(in)        :: p  !< parameters
        class(build), intent(in), target :: b  !< buld instance 
        integer, intent(in), optional    :: sh !< shift (2 be able to zero it)
        integer                          :: i, s, alloc_stat, cnt 
        ! set constants
        o_refs    => b%e
        kfromto   = p%kfromto
        nstates   = p%nstates
        nrefs     = p%nspace
        nobjs     = nstates*nrefs
        npeaks    = p%npeaks
        nsamp     = b%refimgs(1,1)%get_nrots()
        diversify = p%diversify
        ! construct objects
        rt = ran_tabu(nobjs)
        ! make polar references
        if( .not. allocated(refinds) )then
            allocate(refinds(nstates,nrefs), ref_srch_order(nobjs),&
            probs(nstates), probs_prev(nstates), refs(nobjs), proj_space(nobjs),&
            proj_space_corrs(nobjs), proj_space_inds(nobjs), stat=alloc_stat)
            call alloc_err( 'polarft_srch_init; simple_polarft_srch, 1', alloc_stat )
        endif
        ! make object for pair-wise alignment
        if( .not. ap%exists() ) call ap%new(p, sh)
        cnt = 0
        do i=1,nrefs
            do s=1,nstates
                cnt = cnt+1
                if( b%refimgs(s,i)%exists() )then
                    refs(cnt)%p => b%refimgs(s,i)
                else
                    write(*,*) 'state:', s, 'proj:', i
                    stop 'reference image does not exists; polarft_srch_init; simple_polarft_srch'
                endif
            end do
        end do
        cnt = 0
        do i=1,nrefs
            do s=1,nstates
                cnt = cnt+1
                refinds(s,i) = cnt
                proj_space(cnt)%state = s
                proj_space(cnt)%proj  = i
                if( allocated(proj_space(cnt)%inpl_corrs) )then
                    deallocate(proj_space(cnt)%inpl_corrs, proj_space(cnt)%inpl_inds)
                endif
                allocate( proj_space(cnt)%inpl_corrs(nsamp), proj_space(cnt)%inpl_inds(nsamp), stat=alloc_stat ) 
                call alloc_err( 'polarft_srch_init; simple_polarft_srch, 2', alloc_stat )
            end do
        end do
        write(*,'(A)') '>>> DONE INITIALIZING POLARFT_SRCH'
    end subroutine
    
    ! HIGH-LEVEL FUNCTIONALITY (IN ORDER OF EXECUTION)

    !>  \brief  rotationally aligns the ptcl to the reference sections, fills the search space
    !!          shift aligns the ptcl if sh /= 0
    subroutine polarft_srch_align_ptcl( img, o_prev )
        use simple_ori,   only: ori 
        use simple_image, only: image
        use simple_math,  only: put_last, hpsort
        class(image), intent(inout)        :: img
        type(ori), intent(inout), optional :: o_prev
        integer                            :: s, i, j, cnt, which_ref
        integer                            :: endit, projinpl(2)
        real                               :: xy(2)
        ! prepare
        if( present(o_prev) )then
            call prep4srch(img,o_prev,projinpl)
            xy(1) = o_prev%get('x')
            xy(2) = o_prev%get('y')
        else
            call prep4srch(img)
            xy = 0.
        endif
        ! make 'random' projection direction order
        call rt%reset
        call rt%ne_ran_iarr(ref_srch_order)
        if( present(o_prev) )then
            ! put projinpl(1) last to avoid cycling
            call put_last(projinpl(1), ref_srch_order)
        endif
        endit = nrefs
        ! do stochastic rotational search
        cnt = 0
        nbetter = 0
        corr_best_glob = -1.
        proj_space_corrs = -1.
        nrefs_eval = 0
        do i=1,endit
            do s=1,nstates
                cnt = cnt+1
                call ap%set_refimg(refs(ref_srch_order(cnt))%p)
                do j=1,nsamp
                    proj_space(ref_srch_order(cnt))%inpl_inds(j) = j
                    proj_space(ref_srch_order(cnt))%inpl_corrs(j) = ap%rocorr(j)
                    if( proj_space(ref_srch_order(cnt))%inpl_corrs(j) > corr_best_glob )then
                        corr_best_glob = proj_space(ref_srch_order(cnt))%inpl_corrs(j)
                        best_ref = ref_srch_order(cnt)
                    endif  
                end do
                call hpsort(nsamp, proj_space(ref_srch_order(cnt))%inpl_corrs,&
                proj_space(ref_srch_order(cnt))%inpl_inds)
                proj_space_inds(ref_srch_order(cnt))  = ref_srch_order(cnt)
                proj_space_corrs(ref_srch_order(cnt)) = proj_space(ref_srch_order(cnt))%inpl_corrs(nsamp)
                if( npeaks == 1 )then
                    if( proj_space_corrs(ref_srch_order(cnt)) > corr_prev ) nbetter = nbetter+1
                else
                    if( proj_space_corrs(ref_srch_order(cnt)) >= corr_prev ) nbetter = nbetter+1
                endif
                if( diversify .eq. 'yes' .and. nbetter == npeaks ) exit
            end do
            nrefs_eval = nrefs_eval+1
            if( diversify .eq. 'yes' .and. nbetter == npeaks ) exit
        end do
        call hpsort(nobjs, proj_space_corrs, proj_space_inds)
        if( ap%get_sh() == 0 ) return
        ! search the shifts for the npeaks best refs
        do i=nobjs,nobjs-npeaks+1,-1
            which_ref = proj_space(proj_space_inds(i))%proj
            print *, 'searching the shifts for reference: ', which_ref
            print *, 'previous correlation: ', proj_space_corrs(i)
            print *, 'previous origin shift: ', xy
            call ap%set_refimg(refs(which_ref)%p)
            rotxy = ap%align()
            proj_space(which_ref)%x=rotxy(2)
            proj_space(which_ref)%y=rotxy(3)
            print *, 'new correlation: ', ap%get_corr()
            print *, 'new origin shift: ', rotxy(2:3)
            print *, '*********************************'
        end do
    end subroutine
    
    !>  \brief  prepares the align_pair object for alignment
    subroutine prep4srch( img, o_prev, projinpl )
        use simple_ori,   only: ori
        use simple_image, only: image
        class(image), intent(inout)        :: img
        type(ori), intent(inout), optional :: o_prev
        integer, intent(out), optional     :: projinpl(2)
        integer                            :: proj, inpl, state, s
        character(len=STDLEN)              :: dig, key
        ! set ptcl image in align pair object (which takes care of masking)
        call ap%set_pimg(img)
        if( present(o_prev) )then   
            ! check for presence of state probabilities and 
            ! extract the multimodal distribution
            if( nstates == 1 )then
                probs_prev = 1.
            else if( o_prev%isthere('sprob1') )then
                do s=1,nstates
                    write(dig,*) s
                    key = 'sprob'//trim(adjustl(dig))
                    probs_prev(s) = o_prev%get(key)
                end do
            else
                probs_prev = 1./real(nstates)
            endif
            ! find previous discrete alignment parameters
            proj = o_refs%find_closest_proj(o_prev)         ! projection
            inpl = refs(1)%p%get_roind(360.-o_prev%e3get()) ! in-plane angle
            state = nint(o_prev%get('state'))
            if( state > nstates )then
                stop 'previous best state outside boundary; prep4srch; simple_polarft_srch'
            endif
            ! set reference image
            call ap%set_refimg(refs(refinds(state,proj))%p) 
            ! calculate previous correlation
            corr_prev = ap%rocorr(inpl)
            if( present(projinpl) )then
                projinpl(1) = proj
                projinpl(2) = inpl
            endif
        else
            if( present(projinpl) ) stop 'need o_prev to be present for setting projinpl; prep4srch; simple_polarft_srch'    
            corr_prev = 1.
        endif
    end subroutine
    
    !>  \brief  calculates orientation weights directly proportional to 
    !!          the correlation in the feasible subset
    subroutine polarft_srch_calc_weights( weights, wcorr )
        type(soft_ori), intent(out), allocatable :: weights(:)
        real, intent(out)                        :: wcorr
        real, allocatable                        :: corrs(:)
        integer, allocatable                     :: ind_inpls(:) 
        real                                     :: wsum
        integer                                  :: i, j, nbetter_inpl, cnt, alloc_stat
        if( allocated(weights) ) deallocate(weights)
        allocate(weights(npeaks), corrs(npeaks),&
        ind_inpls(nobjs-npeaks+1:nobjs), stat=alloc_stat)
        call alloc_err("In: polarft_srch_calc_weights; simple_polarft_srch", alloc_stat)
        ! update parameters
        cnt   = 0
        probs = 0.
        do i=nobjs,nobjs-npeaks+1,-1
            call param_update
        end do
        ! normalize the state probs
        if( nstates == 1 )then
            probs(1) = 1.
        else
            probs = probs/sum(probs)
            probs = (1.-eps)*probs_prev+eps*(probs)
        endif
        ! update weights
        cnt = 0
        wsum = 0.
        do i=nobjs,nobjs-npeaks+1,-1
            call weight_update
        end do
        wcorr = 0.     
        do i=1,npeaks
            weights(i)%w = weights(i)%w/wsum
            wcorr = wcorr+corrs(i)*weights(i)%w
        end do
        deallocate(corrs, ind_inpls)
        
        contains
            
            subroutine param_update
                use simple_rnd, only: irnd_uni
                cnt = cnt+1
                nbetter_inpl = 1
                if( diversify .eq. 'yes' )then
                    ! count the number of better in-plane rotations
                    nbetter_inpl = 0
                    do j=nsamp,1,-1
                        if( proj_space(proj_space_inds(i))%inpl_corrs(j) >= corr_prev )then
                            nbetter_inpl = nbetter_inpl+1
                            cycle
                        else
                            exit
                        endif
                    end do
                endif
                if( nbetter_inpl > 1 )then
                    ! select a random improving inplane
                    ind_inpls(i) = nsamp-irnd_uni(nbetter_inpl)+1
                else
                    ! select the best in-plane
                    ind_inpls(i) = nsamp
                endif
                ! calculate weights
                weights(cnt)%proj         = proj_space(proj_space_inds(i))%proj
                weights(cnt)%state        = proj_space(proj_space_inds(i))%state
                weights(cnt)%inpl         = proj_space(proj_space_inds(i))%inpl_inds(ind_inpls(i))
                weights(cnt)%x            = proj_space(proj_space_inds(i))%x
                weights(cnt)%y            = proj_space(proj_space_inds(i))%y
                corrs(cnt)                = proj_space(proj_space_inds(i))%inpl_corrs(ind_inpls(i))
                weights(cnt)%corr         = corrs(cnt)  
                probs(weights(cnt)%state) = probs(weights(cnt)%state)+1.                
            end subroutine
            
            subroutine weight_update
                cnt = cnt+1
                if( corrs(cnt) > 0. )then
                    weights(cnt)%w = exp(corrs(cnt)*probs(weights(cnt)%state))
                else
                    weights(cnt)%w = 0.
                endif
                wsum = wsum+weights(cnt)%w
            end subroutine
        
    end subroutine
    
    !>  \brief  to get one orientation from the discrete space
    subroutine polarft_srch_get_ori_best( o )
        use simple_ori,  only: ori
        use simple_math, only: rad2deg
        class(ori), intent(inout) :: o
        type(ori)                 :: o_copy
        real                      :: euls(3), mi, frac, dist
        character(len=STDLEN)     :: dig, key
        integer                   :: proj1, proj2, state1, state2, inpl1, inpl2, class, s
        real                      :: x, y 
        o_copy = o
        ! get old indices
        proj1 = o_refs%find_closest_proj(o)
        inpl1 = refs(1)%p%get_roind(360.-o%e3get())
        state1 = nint(o%get('state'))
        ! get new indices
        proj2  = proj_space(proj_space_inds(nobjs))%proj
        inpl2  = proj_space(proj_space_inds(nobjs))%inpl_inds(nsamp)
        state2 = proj_space(proj_space_inds(nobjs))%state
        euls(1) = o_refs%e1get(proj2)
        euls(2) = o_refs%e2get(proj2)
        euls(3) = 360.-refs(1)%p%get_rot(inpl2) ! change sgn to fit convention
        if( euls(3) == 360. ) euls(3) = 0.
        class = proj_space_inds(nobjs)
        ! calculate overlap between distributions
        mi = 0.
        if( proj1  ==  proj2 ) mi = mi+1
        if( inpl1  ==  inpl2 ) mi = mi+1
        if( state1 == state2 ) mi = mi+1
        mi = mi/3.
        ! set parameters
        call o%set_euler(euls)
        ! shifts must be obtained by vector addition
        xold = o%get('x')
        yold = o%get('y')
        x = xold+proj_space(proj_space_inds(nobjs))%x
        y = yold+proj_space(proj_space_inds(nobjs))%y
        call o%set('x',x)
        call o%set('y',y)
        frac = 100.*(real(nrefs_eval)/real(nrefs))
        dist = 0.5*rad2deg(o_copy.euldist.o)+0.5*o%get('dist')
        call o%set_list(['state  ','class  ','corr   ','dist   ','mi_hard','frac   '],&
        [real(state2),real(class),corr_best_glob,dist,mi,frac])
        ! output multinomal state distribution
        if( nstates > 1 )then
            do s=1,nstates
                write(dig,*) s
                key = 'sprob'//trim(adjustl(dig))
                call o%set(trim(adjustl(key)), probs(s))
            end do
        endif
    end subroutine
    
    !>  \brief  get the state corresponding to class
    function polarft_srch_get_state( class ) result( s )
        integer, intent(in) :: class
        integer             :: s
        s = proj_space(class)%state
    end function
    
    !>  \brief  get the projection direction corresponding to class
    function polarft_srch_get_proj( class ) result( p )
        integer, intent(in) :: class
        integer             :: p
        p = proj_space(class)%proj
    end function
    
    !>  \brief  get the projection direction corresponding to class
    function polarft_srch_get_class( state, proj ) result( class )
        integer, intent(in) :: state, proj
        integer             :: class
        class = refinds(state,proj)
    end function

    !>  \brief  to get one orientation from the discrete space
    subroutine polarft_srch_get_ori( so, o )
        use simple_ori, only: ori
        type(soft_ori), intent(in) :: so
        type(ori), intent(inout)   :: o
        real                       :: euls(3)
        euls(1) = o_refs%e1get(so%proj)
        euls(2) = o_refs%e2get(so%proj)
        euls(3) = 360.-refs(1)%p%get_rot(so%inpl) ! change sgn to fit convention  
        call o%set_euler(euls)
        call o%set('state',real(so%state))
        call o%set('corr',so%corr)
        ! shifts must be obtained by vector addition
        call o%set('x',xold+so%x)
        call o%set('y',yold+so%y)
    end subroutine
    
    !>  \brief  to get one orientation from the discrete space
    subroutine polarft_srch_get_euler( so, o )
        use simple_ori, only: ori
        type(soft_ori), intent(in) :: so
        type(ori), intent(inout)   :: o
        real                       :: euls(3)
        euls(1) = o_refs%e1get(so%proj)
        euls(2) = o_refs%e2get(so%proj)
        euls(3) = 360.-refs(1)%p%get_rot(so%inpl) ! change sgn to fit convention  
        call o%set_euler(euls)
        call o%set('state',real(so%state))
    end subroutine
    
    !>  \brief  standard deviation of soft orientations
    function polarft_srch_ang_sdev( sos ) result( sdev )
        use simple_ori, only: ori
        type(soft_ori), intent(in) :: sos(:)
        real                       :: sdev
        integer                    :: nstates, s
        nstates = maxval(sos(:)%state)
        sdev = 0.
        do s=1,nstates
            sdev = sdev+polarft_srch_ang_sdev_state()
        end do
        sdev = sdev/real(nstates)
    
        contains

            function polarft_srch_ang_sdev_state( ) result( sdev )
                use simple_math, only: rad2deg
                use simple_stat, only: moment
                real                        :: ave, sdev, var
                type(ori)                   :: o_best, o
                integer                     :: loc(1), alloc_stat, i, cnt
                real, allocatable           :: dists(:)
                type(soft_ori), allocatable :: sos_state(:)
                logical                     :: err
                call o_best%new
                call o%new
                ! count the number of occurances of state
                cnt = 0
                do i=1,size(sos)
                    if( sos(i)%state == s ) cnt = cnt+1
                end do
                if( cnt <= 3 ) then ! because one is excluded in the next step & moment needs at least 2 objs
                    sdev = 0.
                    return
                endif
                ! allocate state copy
                allocate( sos_state(cnt), dists(cnt-1), stat=alloc_stat )
                call alloc_err( 'polarft_srch_ang_sdev_state', alloc_stat )
                ! copy data
                cnt = 0
                do i=1,size(sos)
                    if( sos(i)%state == s )then
                        cnt = cnt+1
                        sos_state(cnt) = sos(i)
                    endif
                end do
                loc = maxloc(sos_state(:)%w)
                call polarft_srch_get_euler(sos_state(loc(1)), o_best)       
                cnt = 0
                do i=1,size(sos_state)
                    if( i == loc(1) ) cycle
                    cnt = cnt+1
                    call polarft_srch_get_euler( sos_state(i), o )
                    dists(cnt) = rad2deg(o_best.euldist.o) 
                end do
                call moment(dists, ave, sdev, var, err )
                deallocate(sos_state, dists)
                call o_best%kill
                call o%kill
            end function
            
    end function
    
end module simple_polarft_srch
