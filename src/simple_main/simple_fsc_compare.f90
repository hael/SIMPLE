module simple_fsc_compare
implicit none
    
type :: fsc_compare
    private
    real,    pointer     :: pfscs(:,:) => null()
    integer, allocatable :: order(:)
    real,    allocatable :: state_cnt_diff_wtab(:,:)
    integer              :: n = 0, filtsz = 0
    logical              :: exists = .false.
    logical              :: ranked = .false.
  contains
    procedure          :: best_match
    procedure          :: shc_mask
    procedure, private :: calc_weights_1
    procedure, private :: calc_weights_2
    generic            :: calc_weights => calc_weights_1, calc_weights_2
    procedure, private :: rank
    procedure          :: kill
end type

interface fsc_compare
    module procedure constructor
end interface fsc_compare

integer, parameter :: CNT_THRESH = 6

contains

    !>  \brief  is a constructor
    function constructor( fscs ) result( self )
        use simple_jiffys, only: alloc_err
        real, target, intent(in) :: fscs(:,:)
        type(fsc_compare) :: self
        integer :: alloc_stat
        call self%kill
        self%n      =  size(fscs,1)
        self%filtsz =  size(fscs,2)
        self%pfscs  => fscs
        allocate( self%order(self%n), self%state_cnt_diff_wtab(self%n,3), stat=alloc_stat )
        call alloc_err('constructor; fsc_compare', alloc_stat)
        self%ranked = .false.
        self%exists = .true.
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

    subroutine calc_weights_1( self, weights )
        class(fsc_compare), intent(inout) :: self
        real,               intent(out)   :: weights(self%n)
        integer :: s, i, iplus, ione
        if( .not. self%ranked ) call self%rank
        self%state_cnt_diff_wtab = 0.
        ione = self%order(1)
        do s=1,self%n - 1
            i     = self%order(s)
            iplus = self%order(s + 1)
            self%state_cnt_diff_wtab(i,1) = real(count(self%pfscs(i,:) > self%pfscs(iplus,:)))
            self%state_cnt_diff_wtab(i,2) = 2.*self%state_cnt_diff_wtab(i,1) - real(self%filtsz)
            self%state_cnt_diff_wtab(i,3) = sum(self%pfscs(i,:)) / real(self%filtsz)
        end do
        if( nint(self%state_cnt_diff_wtab(ione,3)) <=  CNT_THRESH )then
            ! we weight
            self%state_cnt_diff_wtab(:,3) = self%state_cnt_diff_wtab(:,3) / sum(self%state_cnt_diff_wtab(:,3))
        else
            self%state_cnt_diff_wtab(:,3)    = 0.
            self%state_cnt_diff_wtab(ione,3) = 1.
        endif
        weights = self%state_cnt_diff_wtab(:,3)
    end subroutine calc_weights_1

    subroutine calc_weights_2( self, mask, weights )
        class(fsc_compare), intent(inout) :: self
        logical,            intent(in)    :: mask(self%n)
        real,               intent(out)   :: weights(self%n)
        integer :: s
        do s=1,self%n
            if( mask(s) ) weights(s) = sum(self%pfscs(s,:)) / real(self%filtsz)
        end do
        weights = weights / sum(weights)
    end subroutine calc_weights_2

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

    subroutine kill( self )
        class(fsc_compare), intent(inout) :: self
        if( self%exists )then
            deallocate( self%order, self%state_cnt_diff_wtab )
            self%pfscs => null()
            self%exists = .false.
        endif
    end subroutine kill

end module simple_fsc_compare