! cluster validation
module simple_cluster_valid
    use simple_defs
    use simple_jiffys, only: alloc_err
implicit none
#include "simple_local_flags.inc"

!>  \brief is the type definition for the validation class, driven by pointers
!>  since the data most likley will anyway be allocated elsewhere
type cluster_valid
    integer              :: nptcls     = 0   !< number of particles being clustered
    integer              :: npairs     = 0   !< number of particle pairs
    integer              :: ncls       = 0   !< number of clusters
    integer              :: ninds      = 0   !< number of indices available
    integer, allocatable :: labels(:)        !< class labels
    integer, allocatable :: centers(:)       !< class centers
    real                 :: dmin_cen   = 0.  !< minimum distance between cluster centers     (separation metric)
    real                 :: sepavg_cen = 0.  !< average distance between cluster centers     (separation metric)
    real, allocatable    :: sepmat(:,:)      !< matrix of average distances between clusters (separation metric)
    real, allocatable    :: mindists(:,:)    !< minimum distance between clusters            (separation metric)
    real, allocatable    :: maxdists(:)      !< maximum distance within clusters             (cohesion metric)
    real, allocatable    :: avgdists(:)      !< average of distances within cluster          (cohesion metric)
    real, allocatable    :: avgdists_cen(:)  !< average of distances memebers and centers    (cohesion metric)
    real, allocatable    :: maxdists_cen(:)  !< maximum distance between members and centers (cohesion metric)
    real, pointer        :: S(:,:)           !< similarity matrix
    real, pointer        :: D(:,:)           !< distance matrix
    logical              :: exists = .false. !< to indicate existence
  contains
    ! CONSTRUCTOR
    procedure :: new
    ! GETTERS
    procedure, private :: get_clspop
    procedure, private :: get_clsarr
    ! CALCULATORS
    procedure, private :: define_centers
    procedure, private :: gen_separation_metrics
    procedure, private :: gen_cohesion_metrics
    ! CLUSTER VALIDITY INDICES
    procedure, private :: cohesion
    procedure, private :: separation
    procedure :: ratio_index
    ! DESTRUCTOR
    procedure :: kill
end type cluster_valid

contains

    ! CONSTRUCTOR

    !>  \brief  is a constructor
    subroutine new( self, o, ncls, which, S, D )
        use simple_oris, only: oris
        class(cluster_valid), intent(inout) :: self   !< instance
        class(oris), intent(inout)          :: o      !< clustering solution
        integer, intent(in)                 :: ncls   !< number of clusters
        character(len=*), intent(in)        :: which  !< tag which=<state|class|subclass>
        real, intent(in), optional, target  :: S(:,:) !< similarity matrix
        real, intent(in), optional, target  :: D(:,:) !< distance matrix
        integer :: alloc_stat, iptcl
        logical :: spresen, dpresen
        call self%kill
        ! set constants
        self%nptcls = o%get_noris()
        self%npairs = (self%nptcls*(self%nptcls-1))/2
        self%ncls   = ncls
        self%ninds  = 2
        ! error check
        if( self%nptcls < 2 ) stop 'number of objects (nptcls) has to be > 1; simple_cluster_valid::new'
        if( self%ncls  < 2 )  stop 'number of classes (ncls) has to be > 1; simple_cluster_valid::new'
        select case(which)
            case('state')
            case('class')
            case('subclass')
            case DEFAULT
                stop 'unknown which tag, which=<state|class|subclass>; simple_cluster_valid::new'
        end select
        ! allocate
        allocate( self%labels(self%nptcls), self%centers(self%ncls), self%maxdists(self%ncls),&
        self%avgdists(self%ncls), self%sepmat(self%ncls,self%ncls), self%mindists(self%ncls,self%ncls),&
        self%avgdists_cen(self%ncls), self%maxdists_cen(self%ncls), stat=alloc_stat )
        call alloc_err( 'In: simple_cluster_valid::new', alloc_stat )
        ! fill-up labels
        do iptcl=1,self%nptcls
            self%labels(iptcl) = nint(o%get(iptcl, which))
        end do
        ! set similarity matrix
        spresen = .false.
        if( present(S) )then
            self%S => S
            spresen = .true.
        endif
        ! set distance matrix
        dpresen = .false.
        if( present(D) )then
            self%D => D
            dpresen = .true.
        endif
        ! check presence or absence of S and D
        if( spresen .and. dpresen )then
            stop 'Cannot use both similarity matrix (S) and distance marix (D)&
            &for validation; simple_cluster_valid::new'
        else if( spresen .or. dpresen )then
            ! all fine
        else
            stop 'Must have either similarity matrix (S) or distance marix (D)&
            &for validation; simple_cluster_valid::new'
        endif
        self%exists = .true.
    end subroutine new

    ! GETTERS

    !>  \brief  get population of a cluster
    function get_clspop( self, class ) result( pop )
        class(cluster_valid), intent(inout) :: self
        integer, intent(in)                 :: class
        integer :: pop, iptcl
        pop = 0
        do iptcl=1,self%nptcls
            if( self%labels(iptcl) == class ) pop = pop+1
        end do
    end function get_clspop

    !>  \brief  get an array of particle indices within a cluster (class)
    function get_clsarr( self, class, pop ) result( arr )
        class(cluster_valid), intent(inout) :: self
        integer, intent(in)                 :: class, pop
        integer, allocatable :: arr(:)
        integer :: alloc_stat, cnt, iptcl
        allocate( arr(pop), stat=alloc_stat )
        call alloc_err( 'In: simple_cluster_valid::get_clsarr', alloc_stat )
        cnt = 0
        do iptcl=1,self%nptcls
            if( self%labels(iptcl) == class )then
                cnt = cnt+1
                arr(cnt) = iptcl
            endif
        end do
    end function get_clsarr

    ! CALCULATORS

    !>  \brief  calculates cluster separation metrics
    subroutine define_centers( self )
        class(cluster_valid), intent(inout) :: self
        integer, allocatable :: karr(:)
        integer :: k, sz, i, j
        real :: sumd, mind, x, d
        do k=1,self%ncls
            sz = self%get_clspop(k)
            if( sz > 0 )then
                karr = self%get_clsarr(k, sz)
                if( sz == 1 )then
                    self%centers(k) = karr(1)
                else
                    mind = huge(x)
                    do i=1,sz
                        sumd = 0.
                        do j=1,sz
                            if( j == i )then
                                d = 0.
                            else
                                if( associated(self%S) )then
                                    d = 1.-self%S(karr(i),karr(j))
                                else
                                    d = self%D(karr(i),karr(j))
                                endif
                            endif
                            sumd = sumd+d
                        end do
                        if( sumd < mind )then
                            self%centers(k) = i
                            mind = sumd
                        endif
                    end do
                endif
            endif
        end do
    end subroutine define_centers

    !>  \brief  calculates cluster separation metrics
    subroutine gen_separation_metrics( self )
        class(cluster_valid), intent(inout) :: self
        integer :: alloc_stat, icls, jcls, i, j, isz, jsz, ncl_pairs
        integer, allocatable :: iclsarr(:), jclsarr(:)
        real    :: d, dmin, x
        ! initialize
        self%sepmat   = 0.
        self%mindists = 0.
        self%dmin_cen = huge(x)
        self%sepavg_cen = 0.
        ! loop over all cluster pairs
        do icls=1,self%ncls-1
            ! get the particle lables for icls
            isz = self%get_clspop(icls)
            if( isz > 0 )then
                iclsarr = self%get_clsarr(icls, isz)
                do jcls=icls+1,self%ncls
                    ! get the particle lables for jcls
                    jsz = self%get_clspop(jcls)
                    if( jsz > 0 )then
                        jclsarr = self%get_clsarr(jcls, jsz)
                        dmin = huge(x)
                        do i=1,isz
                            do j=1,jsz
                                if( associated(self%S) )then
                                    d = 1.-self%S(iclsarr(i),jclsarr(j))
                                else
                                    d = self%D(iclsarr(i),jclsarr(j))
                                endif
                                self%sepmat(icls,jcls) = self%sepmat(icls,jcls)+d
                                if( d < dmin ) dmin = d
                            end do
                        end do
                        self%mindists(icls,jcls) = dmin
                        self%mindists(jcls,icls) = dmin
                        self%sepmat(icls,jcls)   = self%sepmat(icls,jcls)/real(isz*jsz)
                        self%sepmat(jcls,icls)   = self%sepmat(icls,jcls)
                        deallocate(jclsarr)
                        ! generate the center-based metrics
                        if( associated(self%S) )then
                            d = 1.-self%S(self%centers(icls),self%centers(jcls))
                        else
                            d = self%D(self%centers(icls),self%centers(jcls))
                        endif
                        self%sepavg_cen = self%sepavg_cen+d
                        if( d < self%dmin_cen ) self%dmin_cen = d
                    endif
                end do
                deallocate(iclsarr)
            endif
        end do
        ncl_pairs     = (self%ncls*(self%ncls-1))/2
        self%sepavg_cen = self%sepavg_cen/real(ncl_pairs)
    end subroutine gen_separation_metrics

    !>  \brief  calculates an array of all cluster cohesions
    subroutine gen_cohesion_metrics( self )
        class(cluster_valid), intent(inout) :: self
        integer, allocatable :: karr(:)
        integer :: sz, k, alloc_stat, i, j, npairs
        real    :: d, dmax, dsum, dmax_cen
        self%maxdists     = 0.
        self%avgdists     = 0.
        self%avgdists_cen = 0.
        do k=1,self%ncls
            sz = self%get_clspop(k)
            if( sz > 1 )then
                karr = self%get_clsarr(k,sz)
                ! calculate sum and maximum distance between all points assigned to the cluster
                dmax = 0.
                dsum = 0.
                npairs = (sz*(sz-1))/2
                do i=1,sz-1
                    do j=i+1,sz
                        if( associated(self%S) )then
                            d = 1.-self%S(karr(i),karr(j))
                        else
                            d = self%D(karr(i),karr(j))
                        endif
                        dsum = dsum+d
                        if( d > dmax ) dmax = d
                    end do
                end do
                self%maxdists(k) = dmax
                self%avgdists(k) = dsum/real(npairs)
                ! generate the center-based metrics
                dmax_cen = 0.
                do i=1,sz
                    if( associated(self%S) )then
                        d = 1.-self%S(karr(i),self%centers(k))
                    else
                        d = self%D(karr(i),self%centers(k))
                    endif
                    self%avgdists_cen(k) = self%avgdists_cen(k)+d
                    if( d > dmax_cen ) dmax_cen = d
                end do
                self%avgdists_cen(k) = self%avgdists_cen(k)/real(sz)
                self%maxdists_cen(k) = dmax_cen
                deallocate(karr)
            endif
        end do
        DebugPrint  'finished gen_cohesion_metrics'
    end subroutine gen_cohesion_metrics

    !>  \brief  calculates cohesion of the clusters in terms of distance
    !!          we want to minimize cohesion
    function cohesion( self ) result( coh )
        class(cluster_valid), intent(inout) :: self
        integer :: k
        real :: coh
        call self%gen_cohesion_metrics
        coh = 0.
        do k=1,self%ncls
            coh = coh+self%maxdists(k)+self%avgdists(k)
        end do
        coh = coh/real(2*self%ncls)
        DebugPrint  'calculated cohesion'
    end function cohesion

    !>  \brief  calculates separation of the clusters in terms of distance
    !!          we want to maximise separaration
    function separation( self ) result( sep )
        class(cluster_valid), intent(inout) :: self
        integer :: npairs, icls, jcls
        real :: sep
        call self%gen_separation_metrics
        npairs = self%ncls*(self%ncls-1)/2
        sep = 0.
        do icls=1,self%ncls-1
            do jcls=icls+1,self%ncls
                sep = sep+self%mindists(icls,jcls)+self%sepmat(icls,jcls)
            end do
        end do
        sep = sep/real(2*npairs)
        DebugPrint  'calculated separation'
    end function separation

    !>  \brief  calculates the ratio index for cluster validation
    !!          which tries to balance between cohesion and separation
    !!          we want to minimize this index
    function ratio_index( self ) result( eratio )
        class(cluster_valid), intent(inout) :: self
        real :: eratio
        call self%define_centers
        eratio = self%cohesion()/self%separation()
    end function ratio_index

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill( self )
        class(cluster_valid), intent(inout) :: self
        if( self%exists )then
            self%nptcls = 0
            self%npairs = 0
            self%ncls   = 0
            deallocate(self%labels, self%centers, self%sepmat, self%avgdists_cen,&
            self%mindists, self%maxdists, self%avgdists, self%maxdists_cen)
            self%S      => null()
            self%D      => null()
            self%exists =.false.
        endif
    end subroutine kill

end module simple_cluster_valid
