module simple_fsc_compare
implicit none
    
type :: fsc_compare
    private
    real,    pointer     :: pfscs(:,:) => null()
    integer, allocatable :: order(:)
    integer              :: n = 0, filtsz = 0
    logical              :: ranked = .false.
  contains
    procedure          :: best_match
    procedure          :: shc_mask
    procedure, private :: rank

end type

interface fsc_compare
    module procedure constructor
end interface fsc_compare

contains

    !>  \brief  is a constructor
    function constructor( fscs ) result( self )
        use simple_jiffys, only: alloc_err
        real, target, intent(in) :: fscs(:,:)
        type(fsc_compare) :: self
        integer :: alloc_stat
        self%n      =  size(fscs,1)
        self%filtsz =  size(fscs,2)
        self%pfscs  => fscs
        if( allocated(self%order) ) deallocate(self%order)
        allocate( self%order(self%n), stat=alloc_stat )
        call alloc_err('constructor; fsc_compare', alloc_stat)
        self%ranked = .false.
    end function constructor

    !>  \brief  is for finding the best match
    function best_match( self ) result( best )
        class(fsc_compare), intent(inout) :: self
        integer :: best
        if( .not. self%ranked ) call self%rank
        best = self%order(1)
    end function best_match

    !>  \brief  generates a mask for the stochastic hill-climbing
    function shc_mask( self, prev_state ) result( mask )
        class(fsc_compare), intent(inout) :: self
        integer,            intent(in)    :: prev_state
        logical :: mask(self%n)
        integer :: i
        if( .not. self%ranked ) call self%rank
        mask = .false.
        do i=1,self%n
            if( fsc1_gt_fsc2(i,prev_state) ) mask(i) = .true.
        end do

        contains

            logical function fsc1_gt_fsc2( fsc1, fsc2 )
                integer, intent(in) :: fsc1, fsc2
                integer :: count1, count2
                count1 = count(self%pfscs(fsc1,:) > self%pfscs(fsc2,:))
                count2 = self%filtsz - count1
                if( count1 > count2 )then
                    fsc1_gt_fsc2 = .true.
                else
                    fsc1_gt_fsc2 = .false.
                endif
            end function fsc1_gt_fsc2

    end function shc_mask

    !>  \brief  orders FSC:s
    subroutine rank( self )
        use simple_math, only: hpsort
        class(fsc_compare), intent(inout) :: self
        integer :: i, alloc_stat
        self%order = (/(i,i=1,self%n)/)
        call hpsort( self%n, self%order, fsc1_gt_fsc2 )
        self%ranked = .true.

        contains

            logical function fsc1_gt_fsc2( fsc1, fsc2 )
                integer, intent(in) :: fsc1, fsc2
                integer :: count1, count2
                count1 = count(self%pfscs(fsc1,:) > self%pfscs(fsc2,:))
                count2 = self%filtsz - count1
                if( count1 > count2 )then
                    fsc1_gt_fsc2 = .true.
                else
                    fsc1_gt_fsc2 = .false.
                endif
            end function fsc1_gt_fsc2

    end subroutine rank

end module simple_fsc_compare