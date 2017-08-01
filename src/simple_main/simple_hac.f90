!------------------------------------------------------------------------------!
! SIMPLE v3.0         Elmlund & Elmlund Lab          simplecryoem.com          !
!------------------------------------------------------------------------------!
!> Simple hierarchical clustering module
!
!! simple_hac is the SIMPLE class for hierarchical clustering.
!
! The code is
! distributed with the hope that it will be useful, but _WITHOUT_
! _ANY_ _WARRANTY_. Redistribution or modification is regulated by the
! GNU General Public License. *Author:* Hans Elmlund, 2012-02-02.
!
!==Changes are documented below
!
module simple_hac
use simple_sll,    only: sll
use simple_oris,   only: oris
use simple_jiffys, only: alloc_err
implicit none

public :: hac, test_hac
private

type hac
    private
    integer                 :: nobj           !< number of objects in total
    integer                 :: nran           !< number of ojects to cluster
    class(sll), allocatable :: clusters(:)    !< singly linked list data structure
    class(oris), pointer    :: o_ptr=>null()  !< for storing clustering solution
    logical                 :: exists=.false. !< to check existence
  contains
    procedure          :: new
    procedure          :: get_node
    procedure          :: cluster
    procedure, private :: convert_sll_rep
    procedure          :: kill
end type

interface hac
    module procedure constructor
end interface

contains

    !>  \brief  is a constructor
    function constructor( nobj, o, nran ) result( self )
        integer, intent(in)             :: nobj !< nr of objects to cluster
        class(oris), intent(in), target :: o    !< for storing cluster info
        integer, intent(in), optional   :: nran !< size of random sample to cluster
        type(hac)                       :: self
        if( present(nran) )then
            call self%new( nobj, o, nran )
        else
            call self%new( nobj, o )
        endif
    end function

    !>  \brief  is a constructor
    subroutine new( self, nobj, o, nran )
        use simple_ran_tabu,  only: ran_tabu
        class(hac), intent(inout)       :: self !< object
        integer, intent(in)             :: nobj !< nr of objects to cluster
        class(oris), intent(in), target :: o    !< for storing cluster info
        integer, intent(in), optional   :: nran !< size of random sample to cluster
        type(ran_tabu)                  :: rt
        integer                         :: alloc_stat, irnd, i
        self%nobj = nobj
        self%nran = nobj
        if( present(nran) ) self%nran = nran
        self%o_ptr => o
        allocate( self%clusters(self%nran), stat=alloc_stat )
        call alloc_err('new; simple_hac', alloc_stat)
        call self%o_ptr%zero('class') ! zero class vars in oris
        do i=1,self%nran
            call self%clusters(i)%new ! one list per object to cluster
        end do
        if( present(nran) )then       ! nran is a random sample from the larger set of nobj objs
            rt = ran_tabu(nobj)
            do i=1,nran
                irnd = rt%irnd()
                call rt%insert(irnd)
                call self%clusters(i)%add(iarr=[irnd])
            end do
            call rt%kill
        else
            do i=1,self%nobj          ! all included
                call self%clusters(i)%add(iarr=[i])
            end do
        endif
        self%exists = .false.
    end subroutine

    !>  \brief  is a getter
    function get_node( self, a ) result( ax )
        class(hac), intent(inout) :: self
        integer, intent(in)       :: a
        integer, allocatable      :: axa(:)
        integer                   :: ax
        if( a > self%nran .or. a < 0 ) stop 'index out of range; get_node; simple_hac'
        call self%clusters(a)%get(1, iarr=axa)
        ax = axa(1)
    end function

    !>  \brief  Hierachical Agglomerative Clustering of a pair-distance table that is cleared
    !!          after this operation and stored on disk name: pdfile
    subroutine cluster( self, pd_tab, pdfile, ncls, minpop )
        use simple_pair_dtab, only: pair_dtab
        class(hac), intent(inout)       :: self    !< object
        class(pair_dtab), intent(inout) :: pd_tab  !< pairwise distance table
        integer, intent(in)             :: ncls    !< nr of clusters
        character(len=*), intent(in)    :: pdfile  !< pairwise distance table, file
        integer, intent(in), optional   :: minpop  !< minimum class population
        real    :: dist
        integer :: i, j
        integer :: round, nclasses
        ! initialize parameters
        nclasses = self%nran
        round    = 0
        ! make a copy of the pd table before molesting it
        call pd_tab%write(pdfile)
        write(*,'(A)') '>>> HIERARCHICAL CLUSTERING USING THE GROUP AVERAGE METHOD'
        do
            round = round+1
            ! get the best merge out
            call pd_tab%get_dmin_pair(i,j,dist)
            ! merge the cluster pair
            self%clusters(i) = self%clusters(i)%append(self%clusters(j))
            ! merge the entries in the pd-tab
            call pd_tab%merge_pair_d(i,j)
            ! update the number of classes
            nclasses = nclasses-1
            if( round == 1 .or. mod(round,500) == 0 ) then
                write(*,"(1X,A)", advance="no") 'Iter:'
                write(*,"(1X,I7)", advance="no") round
                write(*,"(1X,A)", advance="no") 'Ncls:'
                write(*,"(1X,I7)", advance="no") nclasses
                write(*,"(1X,A)", advance="no") 'Dist:'
                write(*,"(1X,F9.4)" ) dist
            endif
            if( nclasses == ncls ) exit
        end do
        write(*,"(1X,A)", advance="no") 'Iter:'
        write(*,"(1X,I7)", advance="no") round
        write(*,"(1X,A)", advance="no") 'Ncls:'
        write(*,"(1X,I7)") nclasses
        call pd_tab%clear
        call self%convert_sll_rep(minpop)
    end subroutine

    function find_centroids( self, pd_tab, pdfile ) result( centroids )
        use simple_pair_dtab, only: pair_dtab
        class(hac), intent(inout)       :: self   !< object
        class(pair_dtab), intent(inout) :: pd_tab !< pairwise distance table
        character(len=*), intent(in)    :: pdfile !< pairwise distance table, file
        integer, allocatable            :: centroids(:)
        integer :: nonempty, alloc_stat, i, j, k, cnt
        real    :: dsum, dmin, x
        ! read back the distance data
        call pd_tab%read(pdfile)
        ! find number of nonempty classes
        nonempty = 0
        do i=1,self%nobj
            if( self%clusters(i)%size() > 0 ) nonempty = nonempty+1
        end do
        allocate( centroids(nonempty), stat=alloc_stat )
        call alloc_err('In: find_centroids; simple_hac', alloc_stat)
        ! find the centroids
        cnt = 0
        do i=1,self%nobj
            if( self%clusters(i)%size() == 0 )then
                cycle
            else if( self%clusters(i)%size() == 1 )then
                cnt = cnt+1
                centroids(cnt) = self%get_node(1)
            else
                cnt  = cnt+1
                dmin = huge(x)
                dsum = 0.
                do j=1,self%clusters(i)%size()
                    do k=1,self%clusters(i)%size()
                        if( k /= j )then
                            dsum = dsum+pd_tab%get_pair_d(j, k)
                        endif
                    end do
                    if( dsum < dmin )then
                        dmin = dsum
                        centroids(cnt) = j
                    endif
                end do
            endif
        end do
    end function

    !>  \brief  communicate the sll representation to the oris obj
    subroutine convert_sll_rep( self, minpop )
        class(hac), intent(inout)     :: self
        integer, intent(in), optional :: minpop
        integer :: i, j, clsnr, lim
        integer, allocatable :: val(:)
        call self%o_ptr%zero('class')
        lim = 1
        if( present(minpop) ) lim = minpop
        clsnr = 1
        do i=1,self%nran
            if( self%clusters(i)%size() >= lim )then
                do j=1,self%clusters(i)%size()
                    call self%clusters(i)%get(j, iarr=val)
                    call self%o_ptr%set(val(1), 'class', real(clsnr))
                end do
                clsnr = clsnr+1
            endif
        end do
    end subroutine

    ! HAC UNIT TEST

    !>  \brief  is the hac unit test
    subroutine test_hac
        use simple_pair_dtab, only: pair_dtab
        use simple_math
        real :: datavecs(90,5)
        type(pair_dtab) :: pd
        type(oris)      :: o
        type(hac)       :: hacls
        integer :: i, j
        write(*,'(a)') '**info(simple_hac_unit_test): testing all functionality'
        call o%new(90)
        pd = pair_dtab(90)
        call hacls%new(90, o)
        ! make data
        do i=1,30
            datavecs(i,:) = 1.
        end do
        do i=31,60
            datavecs(i,:) = 3.
        end do
        do i=61,90
            datavecs(i,:) = 5.
        end do
        ! build pair_dtab
        do i=1,89
            do j=i+1,90
                call pd%set_pair_d(i, j, euclid(datavecs(i,:),datavecs(j,:)) )
            end do
        end do
        call hacls%cluster(pd, 'pdfile.bin', 3, 1 )
        call o%write( 'test_hac_clsdoc.txt' )
        call hacls%kill
        call o%kill
        call pd%kill
        write(*,'(a)') 'SIMPLE_HAC_UNIT_TEST COMPLETED WITHOUT TERMINAL BUGS ;-)'
        write(*,'(a)') 'PLEASE, INSPECT THE RESULTS'
    end subroutine

    !>  \brief  is a destructor
    subroutine kill( self )
        class(hac), intent(inout) :: self !< object
        integer :: i
        if( self%exists )then
            do i=1,self%nran
                call self%clusters(i)%kill
            end do
            self%o_ptr => null()
        endif
    end subroutine

end module simple_hac
