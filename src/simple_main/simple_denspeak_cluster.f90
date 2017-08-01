!------------------------------------------------------------------------------!
! SIMPLE v3.0         Elmlund & Elmlund Lab          simplecryoem.com          !
!------------------------------------------------------------------------------!
!> Simple cluster class : peak density
module simple_denspeak_cluster
use simple_jiffys, only: alloc_err
implicit none

type denspeak_cluster
    private
    integer              :: N=0
    real,    pointer     :: dmat(:,:)=>null()
    real                 :: d_crit=0.
    real,    allocatable :: rhos(:)
    real,    allocatable :: deltas(:)
    integer, allocatable :: centers(:)
    integer, allocatable :: labels(:)
    logical              :: exists=.false.
  contains
    procedure          :: new
    procedure          :: cluster
    procedure          :: get_ncls
    procedure          :: get_centers
    procedure          :: get_labels
    procedure, private :: calc_rhos_and_deltas
    procedure, private :: identify_centers
    procedure, private :: rho
    procedure, private :: delta
    procedure          :: kill
end type denspeak_cluster

interface denspeak_cluster
    module procedure constructor
end interface denspeak_cluster

contains

    function constructor( N, dmat, d_crit ) result( self )
        integer,      intent(in) :: N
        real, target, intent(in) :: dmat(N,N)
        real,         intent(in) :: d_crit
        type(denspeak_cluster)   :: self
        call self%new(N, dmat, d_crit)
    end function constructor

    subroutine new( self, N, dmat, d_crit )
        class(denspeak_cluster), intent(inout) :: self
        integer,                 intent(in)    :: N
        real, target,            intent(in)    :: dmat(N,N)
        real,                    intent(in)    :: d_crit
        integer :: alloc_stat
        self%N      =  N
        self%dmat   => dmat
        self%d_crit =  d_crit
        allocate( self%rhos(N), self%deltas(N), stat=alloc_stat )
        call alloc_err("In: simple_denspeak_cluster :: new", alloc_stat)
        self%exists = .true.
    end subroutine new

    subroutine cluster( self )
        class(denspeak_cluster), intent(inout) :: self
        call self%calc_rhos_and_deltas
        call self%identify_centers
    end subroutine cluster

    integer function get_ncls( self )
        class(denspeak_cluster), intent(inout) :: self
        if( allocated(self%centers) )then
            get_ncls = size(self%centers)
        else
            stop 'centers not allocated; simple_denspeak_cluster :: get_ncls'
        endif
    end function get_ncls

    function get_labels( self ) result( labels )
        class(denspeak_cluster), intent(inout) :: self
        integer, allocatable :: labels(:)
        if( allocated(self%labels) )then
            allocate(labels(self%N), source=self%labels)
        else
            stop 'labels not allocated; simple_denspeak_cluster :: get_labels'
        endif
    end function get_labels

    function get_centers( self ) result( centers )
        class(denspeak_cluster), intent(inout) :: self
        integer, allocatable :: centers(:)
        integer :: ncls
        if( allocated(self%centers) )then
            ncls = size(self%centers)
            allocate(centers(ncls), source=self%centers)
        else
            stop 'centers not allocated; simple_denspeak_cluster :: get_centers'
        endif
    end function get_centers

    subroutine calc_rhos_and_deltas( self )
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(denspeak_cluster), intent(inout) :: self
        integer :: i, loc(1)
        !$omp parallel default(shared) private(i) proc_bind(close)
        !$omp do schedule(static)
        do i=1,self%N
            self%rhos(i) = self%rho(i)
        end do
        !$omp end do
        !$omp do schedule(static)
        do i=1,self%N
            self%deltas(i) = self%delta(i)
        end do
        !$omp end do nowait
        !$omp end parallel
        ! convention
        loc = maxloc(self%rhos)
        self%deltas(loc(1)) = maxval(self%dmat(loc(1),:))
    end subroutine calc_rhos_and_deltas

    subroutine identify_centers( self )
        use simple_math, only: hpsort
        class(denspeak_cluster), intent(inout) :: self
        real    :: delta_peak, delta_fuzz, delta_peak_new, delta_fuzz_new
        real    :: dist_peak, dist_fuzz
        real    :: deltas_sorted(self%N)
        integer :: kmit, idelta, pop_peak, pop_fuzz, cnt, j, k, loc(1), alloc_stat
        logical :: j_is_center
        integer, parameter :: MAXKMIT = 20
        real, allocatable  :: distances(:)
        ! binary delta quantisation
        deltas_sorted = self%deltas
        call hpsort(self%N, deltas_sorted)
        delta_fuzz = sum(deltas_sorted(:self%N/2))/real(self%N/2)
        delta_peak = sum(deltas_sorted(self%N/2+1:))/real(self%N/2)
        do kmit=1,MAXKMIT
            pop_peak       = 0
            pop_fuzz       = 0
            delta_peak_new = 0.
            delta_fuzz_new = 0.
            do idelta=1,self%N
                dist_peak = sqrt((delta_peak-self%deltas(idelta))**2.0)
                dist_fuzz = sqrt((delta_fuzz-self%deltas(idelta))**2.0)
                if( dist_peak <= dist_fuzz )then
                    delta_peak_new = delta_peak_new+self%deltas(idelta)
                    pop_peak       = pop_peak + 1
                else
                    delta_fuzz_new = delta_fuzz_new+self%deltas(idelta)
                    pop_fuzz       = pop_fuzz + 1
                endif
            end do
            delta_peak = delta_peak_new/real(pop_peak)
            delta_fuzz = delta_fuzz_new/real(pop_fuzz)
        end do
        ! re-calculate pop peak (since the centers have been updated)
        pop_peak = 0
        do idelta=1,self%N
            dist_peak = sqrt((delta_peak-self%deltas(idelta))**2.0)
            dist_fuzz = sqrt((delta_fuzz-self%deltas(idelta))**2.0)
            if( dist_peak <= dist_fuzz )then
                pop_peak = pop_peak + 1
            endif
        end do
        ! identify centers & create labels
        allocate(self%centers(pop_peak), self%labels(self%N), distances(pop_peak), stat=alloc_stat)
        call alloc_err("In: simple_denspeak_cluster :: identify_centers", alloc_stat)
        cnt = 0
        do idelta=1,self%N
            dist_peak = sqrt((delta_peak-self%deltas(idelta))**2.0)
            dist_fuzz = sqrt((delta_fuzz-self%deltas(idelta))**2.0)
            if( dist_peak <= dist_fuzz )then
                cnt = cnt + 1
                self%centers(cnt) = idelta
            endif
        end do
        ! if( cnt /= pop_peak ) stop 'cnt /= pop_peak; simple_denspeak_cluster :: identify_centers'
        ! if( cnt == 0 )            stop 'no clusters identified; simple_denspeak_cluster :: identify_centers'
        ! create labels
        do j=1,self%N
            j_is_center = .false.
            do k=1,pop_peak
                if( self%centers(k) == j )then
                    self%labels(j) = k
                    j_is_center = .true.
                endif
            end do
            if( j_is_center ) cycle
            do k=1,pop_peak
                ! if( self%centers(k) < 1 .or. self%centers(k) > self%N )then
                !     write(*,*) 'center nr: ', k, ' is: ', self%centers(k)
                !     stop 'center out of bound; simple_denspeak_cluster :: identify_centers'
                ! endif
                distances(k) = self%dmat(self%centers(k),j)
            end do
            loc = minloc(distances)
            self%labels(j) = loc(1)
        end do
    end subroutine identify_centers

    real function rho( self, i )
        class(denspeak_cluster), intent(in) :: self
        integer,                 intent(in) :: i
        real :: diff(self%N)
        real :: binarr(self%N)
        diff   = self%dmat(i,:)-self%d_crit
        binarr = 0.
        where(diff < 0.) binarr = 1.
        rho = sum(binarr)
    end function rho

    real function delta( self, i )
        class(denspeak_cluster), intent(in) :: self
        integer,                 intent(in) :: i
        logical :: mask(self%N)
        mask = .false.
        where( self%rhos > self%rhos(i) ) mask = .true.
        if( all(.not. mask))then
            delta = 0.
        else
            delta = minval(self%dmat(i,:), mask=mask)
        endif
    end function delta

    subroutine kill( self )
        class(denspeak_cluster), intent(inout) :: self
        if( self%exists )then
            self%dmat => null()
            deallocate(self%rhos, self%deltas)
            if( allocated(self%centers) ) deallocate(self%centers)
            if( allocated(self%labels)  ) deallocate(self%labels)
            self%exists = .false.
        endif
    end subroutine kill

end module simple_denspeak_cluster
