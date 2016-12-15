module simple_refinecluster
use simple_params,           only: params
use simple_build,            only: build
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_oris,             only: oris
use simple_ori,              only: ori
use simple_image,            only: image
use simple_jiffys,           only: alloc_err, progress
use simple_math,             only: rad2deg, round2even
use simple_pftcc_shsrch      ! singleton
use simple_defs              ! singleton

implicit none

public :: refinecluster
private

type refinecluster

    private
    class(params), pointer, private :: pp =>null()  !< pointer to params
    class(build),  pointer, private :: pb =>null()  !< pointer to builder
    type(oris), pointer             :: pos=>null()  !< pointer to original oris 
    type(oris)               :: os                  !< local copy of oris 
    type(polarft_corrcalc)   :: pftcc               !< PFTCC object
    type(image)              :: avg                 !< local cavg
    type(image), allocatable :: ptcls(:)            !< local copy of ptcl stack
    real,        allocatable :: angtab(:)           !< list of in-plane angles
    real                     :: lims(2,2)           !< shift srch boundaries
    integer                  :: n         = 0       !< number of ptcls
    integer                  :: nrots     = 0       !< number of in-plane rotations in polar representation
    integer                  :: class     = 0       !< class being processed
    integer                  :: maxits    = 50      !< maximum number of iterations
    real                     :: trs       = 0.      !< shift range parameter [-trs,trs]
    real                     :: angthresh = 1.      !< inpl angular threshold for search convergence  
    real                     :: trsthresh = 1.      !< inpl shift threshold for search convergence  
    logical                  :: doshift   = .true.  !< origin shift search indicator
    logical                  :: inplace   = .false.
    logical                  :: exists    = .false.
    logical                  :: doprint   = .true. 

  contains
    ! CONSTRUCTOR & INITIALIZER
    procedure :: new
    ! CALCULATORS
    procedure, private :: make_avg
    procedure, private :: calc_roind
    procedure, private :: prep_ptcls
    ! GETTERS & SETTERS
    procedure :: get_cavg
    ! SEARCH ROUTINE
    procedure :: exec
    procedure, private :: greedy_srch
    ! DESTRUCTOR
    procedure :: kill

end type

contains

    !>  \brief  is a constructor
    subroutine new( self, b, p, cls, os, imgs, doprint )
        class(refinecluster), intent(inout)            :: self
        class(build),  intent(inout), target           :: b
        class(params), intent(in)   , target           :: p
        type(oris), intent(inout),    target, optional :: os
        type(image), intent(in),              optional :: imgs(:)
        integer, intent(in),                  optional :: cls
        logical, intent(in),                  optional :: doprint
        type(ori)                            :: o
        integer, allocatable                 :: inds(:)
        integer                              :: alloc_stat, i, state
        real                                 :: dang
        call self%kill
        ! pointers
        self%pb => b
        self%pp => p
        if( .not.present( cls) )then
            self%inplace = .true.
            if( .not.present(os)   )stop 'Orientations not provided in simple_refinecluster%new'
            if( .not.present(imgs) )stop 'Images not provided in simple_refinecluster%new'
            self%pos   => os
            self%class = nint(os%get(1,'class'))
            self%n     = os%get_noris()
            if( self%n.ne.size(imgs) ) &
                stop 'Inconsistent number of images and oris in refinecluster%new'
        else
            self%inplace = .false.
            self%class   = cls
            self%n       = b%a%get_clspop( cls )
            inds = b%a%get_cls_pinds( cls )
        endif
        ! set constants
        self%nrots   = round2even(twopi*real(p%ring2))
        self%trs     = p%trs
        self%doshift = p%doshift
        self%maxits  = p%maxits 
        if( present(doprint) )self%doprint=doprint
        ! retrieve oris
        if( self%n==0 )then
            write(*,'(A,I6)')'Empty class in new; simple_cluster_refine',self%class
            return
        else if( self%n==1 )then
            write(*,'(A,I6)')'Nothing to refine in new; simple_cluster_refine for class',self%class
            return
        endif
        call self%os%new( self%n )
        do i=1,self%n
            if( .not.self%inplace )then 
                o = b%a%get_ori( inds(i) )
            else
                o = os%get_ori( i )
            endif
            if( p%rnd=='yes' )then
                ! randomised starting cavg
                state = nint( o%get('state') )
                call o%rnd_ori(p%trs)
                call o%set('state',real(state))
            endif
            call self%os%set_ori( i,o )
        enddo
        ! construct composites
        allocate( self%angtab(self%nrots), self%ptcls(self%n), stat=alloc_stat )
        call alloc_err('In: new; simple_prime2d_srch, 1', alloc_stat)
        ! generate the list of in-plane angles and indices
        dang = twopi/real(self%nrots)
        do i=1,self%nrots
            self%angtab(i) = rad2deg((i-1)*dang)
        end do
        self%angthresh = rad2deg(dang)/pi  ! angular threshold for convergence
        self%trsthresh = 1./p%smpd         ! shift threshold for convergence
        ! pftcc
        call self%pftcc%new(self%n, [1,self%n], [p%box,p%box,1], p%kfromto, p%ring2, 'no')
        ! gather images
        do i=1,self%n
            call self%ptcls(i)%new( [p%box,p%box,1],p%smpd,p%imgkind )
            if( .not.self%inplace )then
                call self%ptcls(i)%read( p%stk,inds(i) )
            else
                call self%ptcls(i)%copy( imgs(i) )
            endif                
        enddo
        ! cavg
        call self%avg%new( [p%box,p%box,1],p%smpd,p%imgkind )
        self%avg = 0.
        ! shift search
        if( self%doshift )then
            self%lims(:,1) = -self%trs
            self%lims(:,2) = self%trs
        endif
        if( allocated(inds) )deallocate( inds )
        self%exists = .true.
    end subroutine

    subroutine exec( self )
        class(refinecluster), intent(inout) :: self
        type(ori)                           :: o
        integer, allocatable                :: inds(:)
        real                    :: corr, ang_dist, e3_prev, e3diff
        real                    :: x_prev, y_prev, sh_dist, shdiff
        integer                 :: iptcl, ind, roind, cnt
        logical                 :: converged
        write(*,'(A,I6,A,I6)')'>>> REFINING CLUSTER',self%class,'; POP=',self%n
        !if( .not.self%doprint )call progress( 0,self%maxits )
        ! INIT
        ! init ptcls
        call self%prep_ptcls
        ! init avg
        call self%make_avg
        call self%avg%mask(self%pp%msk, 'soft')
        call self%pb%proj%img2polarft(1, self%avg, self%pftcc, isptcl=.false.)
        ! init correlations
        do iptcl=1,self%n
            roind = self%calc_roind(360.-self%os%e3get(iptcl)) ! in-plane angle index
            corr  = self%pftcc%corr( 1,iptcl,roind )
            call self%os%set( iptcl,'corr',corr )
        enddo
        ! MASTER LOOP
        converged = .false.
        cnt       = 0
        corr      = self%os%get_avg('corr')
        if( self%doprint )write(*,'(A,I6,A,F8.3,A,F8.3)')'Step',cnt,'; corr=',corr
        do while( .not.converged )
            cnt = cnt+1
            ! per ptcl shc srch
            ang_dist  = 0.
            sh_dist   = 0.
            do iptcl=1,self%n
                o       = self%os%get_ori( iptcl )
                e3_prev = o%e3get()
                x_prev  = o%get('x')
                y_prev  = o%get('y')
                call self%greedy_srch( iptcl,o )
                call self%os%set_ori( iptcl,o )
                e3diff = abs(e3_prev-o%e3get())
                if( e3diff>180. )e3diff=360.-e3diff
                ang_dist = ang_dist+e3diff/real(self%n)
                shdiff   = sqrt( (o%get('x')-x_prev)**2+(o%get('y')-y_prev)**2 )
                sh_dist  = sh_dist+shdiff/real(self%n)
            enddo
            corr = self%os%get_avg('corr')
            ! convergence & output
            if( (ang_dist<self%angthresh.and.sh_dist<self%trsthresh) .or. cnt==self%maxits )converged=.true.
            if( self%doprint )then
                write(*,'(A,I6,A,F8.3,A,F8.3,A,F8.3)')'Step',cnt,'; corr=',corr,'; avg ang dist=', ang_dist, &
                    '; avg trs=', sh_dist
            ! else
            !     call progress( cnt,self%maxits )
            endif
            ! updates ptcls & cavg
            call self%prep_ptcls
            call self%make_avg
            call self%avg%mask(self%pp%msk, 'soft')
            call self%pb%proj%img2polarft(1, self%avg, self%pftcc, isptcl=.false.)
        enddo
        ! THE END
        if( self%inplace )then
            do iptcl=1,self%n
                call self%pos%set_ori( iptcl,self%os%get_ori( iptcl ) )
            enddo            
        else
            inds = self%pb%a%get_cls_pinds( self%class )
            do iptcl=1,self%n
                ind = inds( iptcl )
                call self%pb%a%set_ori( ind,self%os%get_ori( iptcl ) )
            enddo
            deallocate(inds)
        endif
        ! if( .not.self%doprint )call progress( 1,1 )
    end subroutine exec

    function get_cavg( self )result( img )
        class(refinecluster), intent(inout) :: self
        type(image) :: img
        call img%copy( self%avg )
    end function get_cavg

    !>  \brief  translates & rotates & mask ptcls prior to each srch iteration.
    subroutine prep_ptcls( self )
        class(refinecluster), intent(inout) :: self
        type(ori)                           :: o
        integer                             :: iptcl
        do iptcl=1,self%n
            call self%pb%img%copy( self%ptcls(iptcl) )
            if( self%doshift )then
                o = self%os%get_ori( iptcl )
                call self%pb%img%fwd_ft
                call self%pb%img%shift(-o%get('x'), -o%get('y'))
                call self%pb%img%bwd_ft
            endif
            call self%pb%img%mask(self%pp%msk, 'soft')
            call self%pb%proj%img2polarft(iptcl, self%pb%img, self%pftcc, .true.)
        enddo
    end subroutine prep_ptcls

    !>  \brief  executes the stochastic rotational search
    subroutine greedy_srch( self, iptcl, o )
        class(refinecluster), intent(inout)    :: self
        type(ori), intent(inout)               :: o
        integer, intent(in)                    :: iptcl
        real                :: corr_prev, inpl_corr, x,y
        real                :: shvec(2), cxy(3), corrs(self%nrots)
        integer             :: loc(1), rot, inpl_ind
        corr_prev    = o%get('corr')
        rot          = self%calc_roind( 360.-o%e3get() )
        corrs = self%pftcc%gencorrs(1, iptcl)
        if( any(corrs > corr_prev) )then
            loc       = maxloc(corrs)          
            inpl_corr = corrs( loc(1) )
            inpl_ind  = loc(1)
        else
            inpl_corr = corr_prev
            inpl_ind  = rot
        endif
        ! search shifts
        x = o%get('x')
        y = o%get('y')
        if( self%doshift )then
            call pftcc_shsrch_init( self%pftcc,self%lims )
            call pftcc_shsrch_set_indices(1, iptcl, inpl_ind)
            cxy       = pftcc_shsrch_minimize()
            shvec     = cxy(2:3)
            inpl_corr = cxy(1)
            shvec(1)  = shvec(1)+x
            shvec(2)  = shvec(2)+y
        else
            shvec = [x,y]
        endif
        ! best in plane parms
        call o%e3set( 360.-self%angtab(inpl_ind) )
        call o%set( 'x',shvec(1) )
        call o%set( 'y',shvec(2) )
        call o%set( 'corr',inpl_corr )
    end subroutine

    subroutine make_avg( self )
        class(refinecluster), intent(inout)    :: self
        type(ori)   :: o
        integer     :: iptcl
        real        :: x, y
        self%avg = 0.
        ! loop over all particles
        do iptcl=1,self%n
            o = self%os%get_ori(iptcl)
            call self%pb%img%copy( self%ptcls(iptcl) )
            ! shift & rotate
            call self%pb%img%fwd_ft
            x = o%get('x')
            y = o%get('y')
            call self%pb%img%shift(-x, -y)
            call self%pb%img%bwd_ft
            call self%pb%img%rtsq(-o%e3get(), 0., 0.)
            ! sum
            call self%avg%add( self%pb%img )
        end do
        ! div and noise norm
        call self%avg%div( real(self%n) )
        call self%avg%noise_norm( self%pp%msk )
    end subroutine

    !>  \brief calculate the in-plane rotational index for the rot in-plane angle
    function calc_roind( self, rot )result( ind )
        class(refinecluster), intent(inout) :: self
        real, intent(in)                    :: rot
        integer           :: ind, i, loc(1)
        real              :: dists_sq(self%nrots)
        do i=1,self%nrots
            dists_sq(i) = (self%angtab(i)-rot)**2.
        end do
        loc = minloc(dists_sq)
        ind = loc(1)
    end function

    ! DESTRUCTOR
    subroutine kill( self )
        class(refinecluster), intent(inout)    :: self
        integer :: i
        self%pb  => null()
        self%pp  => null()
        self%pos => null()
        call self%pftcc%kill
        call self%os%kill
        if( self%avg%exists() )call self%avg%kill
        if( allocated(self%ptcls) )then
            do i=1,self%n
                if( self%ptcls(i)%exists() )call self%ptcls(i)%kill
            enddo
            deallocate( self%ptcls )
        endif
        if( allocated(self%angtab) )deallocate( self%angtab )
        self%inplace = .false.
        self%exists  = .false.
    end subroutine

end module simple_refinecluster
