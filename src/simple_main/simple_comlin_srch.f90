module simple_comlin_srch
use simple_comlin_corr ! use all in there
use simple_build,      only: build
use simple_jiffys,     only: progress
implicit none

public :: comlin_srch_init, comlin_srch_symaxis
private

class(build), pointer :: bp=>null()      !< pointer to builder
integer,      pointer :: ptcl=>null()    !< ptcl index
real,         pointer :: dynlim=>null()  !< dynamic lowpass limit
integer               :: nptcls=0        !< nr of particles

contains

    !>  \brief  is a constructor
    subroutine comlin_srch_init( b, p )
        use simple_build,  only: build
        use simple_params, only: params
        class(build), target,  intent(in) :: b
        class(params), target, intent(in) :: p
        nptcls = p%nptcls
        bp     => b
        ptcl   => p%ptcl
        dynlim => p%lp_dyn
        ! make comlin_corr functionality:
        call comlin_corr_init( b, ptcl, dynlim )
    end subroutine comlin_srch_init
    
    !>  \brief  is for calculating the joint common line correlation
    function comlin_srch_corr() result( corr )
        real    :: corr
        integer :: i
        corr = 0.
        do i=1,nptcls
            ptcl = i
            corr = corr+pcorr_comlin()
        end do
        corr = corr/real(nptcls)
    end function comlin_srch_corr
    
    !>  \brief  is for finding the symmetry axis give an aligned set of images
    subroutine comlin_srch_symaxis( orientation_best, doprint )
        use simple_ori,  only: ori
        use simple_oris, only: oris
        class(ori), intent(inout) :: orientation_best
        logical,    intent(in)    :: doprint
        integer    :: i, j, k
        type(ori)  :: orientation
        type(oris) :: a_copy
        real       :: corr, corr_best
        integer    :: ntot, cnt, lims(3,2)
        if( doprint ) write(*,'(A)') '>>> GLOBAL SYMMETRY AXIS SEARCH'
        a_copy    = bp%a
        corr_best = -1.
        ntot      = 24624
        cnt       = 0
        lims(1,1) = 0
        lims(1,2) = 359
        lims(2,1) = 0
        lims(2,2) = 180
        lims(3,1) = 0
        lims(3,2) = 359
        call orientation%new
        do i=lims(1,1),lims(1,2),10
            call orientation%e1set(real(i))
            do j=lims(2,1),lims(2,2),10
                call orientation%e2set(real(j))
                do k=lims(3,1),lims(3,2),10
                    cnt = cnt+1
                    if( doprint ) call progress(cnt, ntot)
                    call orientation%e3set(real(k))
                    bp%a = a_copy
                    call bp%a%rot(orientation)
                    call bp%se%apply2all(bp%a)
                    corr = comlin_srch_corr()
                    call orientation%set('corr',corr)
                    if( corr > corr_best )then
                        corr_best = corr
                        orientation_best = orientation
                    endif
                end do
            end do
        end do
        if( doprint )then
            write(*,'(A)') '>>> FOUND FIRST SYMMETRY AXIS ORIENTATION'
            call orientation_best%print
            write(*,'(A)') '>>> REFINED SYMMETRY AXIS SEARCH'
        endif
        lims(1,1) = max(nint(orientation_best%e1get()-10.),lims(1,1))
        lims(1,2) = min(nint(orientation_best%e1get()+10.),lims(1,2))
        lims(2,1) = max(nint(orientation_best%e2get()-10.),lims(2,1))
        lims(2,2) = min(nint(orientation_best%e2get()+10.),lims(2,2))
        lims(3,1) = max(nint(orientation_best%e3get()-10.),lims(3,1))
        lims(3,2) = min(nint(orientation_best%e3get()+10.),lims(3,2))
        ntot      = (lims(1,2)-lims(1,1)+1)*(lims(2,2)-lims(2,1)+1)*(lims(3,2)-lims(3,1)+1)
        cnt       = 0 
        do i=lims(1,1),lims(1,2)
            call orientation%e1set(real(i))
            do j=lims(2,1),lims(2,2)
                call orientation%e2set(real(j))
                do k=lims(3,1),lims(3,2)
                    cnt = cnt+1
                    if( doprint ) call progress(cnt, ntot)
                    call orientation%e3set(real(k))
                    bp%a = a_copy
                    call bp%a%rot(orientation)
                    call bp%se%apply2all(bp%a)
                    corr = comlin_srch_corr()
                    call orientation%set('corr',corr)
                    if( corr > corr_best )then
                        corr_best = corr
                        orientation_best = orientation
                    endif
                end do
            end do
        end do
        if( doprint )then
            write(*,'(A)') '>>> FOUND REFINED SYMMETRY AXIS ORIENTATION'
            call orientation_best%print
        endif
        bp%a = a_copy
        call a_copy%kill
    end subroutine comlin_srch_symaxis

end module simple_comlin_srch
