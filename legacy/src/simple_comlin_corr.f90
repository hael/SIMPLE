module simple_comlin_corr
use simple_build, only: build
use simple_ori,   only: ori
implicit none

class(build), pointer, private :: bp=>null()
integer,      pointer, private :: pptcl=>null()
real,         pointer, private :: plp_dyn=>null()
type(ori),             private :: o_old
logical,               private :: assoc=.false.
integer                        :: varwhich

interface pcorr_comlin
    module procedure pcorr_comlin_1
    module procedure pcorr_comlin_2
end interface pcorr_comlin

contains

    !>  \brief  initializes the comlin_corr singleton
    subroutine comlin_corr_init( b, ptcl, lp_dyn )
        use simple_build, only: build
        class(build), target, intent(in) :: b
        integer,      target, intent(in) :: ptcl
        real,         target, intent(in) :: lp_dyn
        bp      => b
        pptcl   => ptcl
        plp_dyn => lp_dyn
        assoc   = .true.
    end subroutine comlin_corr_init
    
    !>  \brief  is for calculating the per-ptcl continuous common line correlation coefficient
    !!          modifies the clins data structure (make sure to store the one_ptcl contribution)
    function pcorr_comlin_1() result( corr )
        real :: corr
        corr = -1.
        if( assoc )then
            corr = bp%clins%pcorr( pptcl, plp_dyn )
        else
            write(*,'(a)') 'Associate comlin_corr pointers before calculating correlation!'
            write(*,'(a)') 'pcorr_comlin; simple_comlin_corr'
            stop
        endif
    end function pcorr_comlin_1
    
    !>  \brief  is for calculating the per-ptcl continuous pairwise common line correlation coefficient
    function pcorr_comlin_2( iptcl, jptcl ) result( corr )
        integer, intent(in) :: iptcl, jptcl
        real :: corr
        corr = -1.
        if( assoc )then
            corr = bp%clins%pcorr(iptcl, jptcl, plp_dyn)
        else
            write(*,'(a)') 'Associate comlin_corr pointers before calculating correlation!'
            write(*,'(a)') 'pcorr_comlin; simple_comlin_corr'
            stop
        endif
    end function pcorr_comlin_2
    
end module simple_comlin_corr
