module simple_crossover
use simple_ori,  only: ori
use simple_oris, only: oris
use simple_rnd ! singleton
implicit none

contains

    !>  \brief  most common GA crossover
    subroutine blx_alpha( D, xparent, yparent, alpha, offspring )
        integer, intent(in)  :: D                      !< solution vector dimension
        real,    intent(in)  :: xparent(D), yparent(D) !< parents
        real,    intent(in)  :: alpha                  !< positive real parameter, typically 0.5
        real,    intent(out) :: offspring(2,D)         !< the two children
        real                 :: diffs(D), harvest(D), lb(D), ub(D)
        integer              :: i
        diffs = abs(xparent - yparent)
        do i=1,D
            lb(i) = min(xparent(i),yparent(i)) - alpha * diffs(i)
            ub(i) = max(xparent(i),yparent(i)) + alpha * diffs(i)
        end do
        call ran3arr(harvest)
        offspring(1,:) = harvest * (ub - lb) + lb
        call ran3arr(harvest)
        offspring(2,:) = harvest * (ub - lb) + lb
    end subroutine blx_alpha

    !>  \brief  for generation a subset of orientations using randomised crossover
    function rnd_ori_xover_subset( os_in ) result( os_out )
        class(oris), intent(inout) :: os_in
        type(ori)  :: offspring, oi, oj
        type(oris) :: os_out
        integer    :: iref, jref, nsub, cnt, noris
        noris = os_in%get_noris()
        nsub  = (noris*(noris-1))/2
        call os_out%new(nsub)
        cnt = 0
        do iref=1,noris-1
            oi = os_in%get_ori(iref)
            do jref=iref+1,noris
                cnt = cnt + 1
                oj  = os_in%get_ori(jref)
                call rnd_ori_xover(oi, oj, offspring)
                call os_out%set_ori(cnt, offspring)
            end do
        end do
    end function rnd_ori_xover_subset

    !>  \brief  randomised crossover of orientations
    subroutine rnd_ori_xover( xparent, yparent, offspring )
        class(ori), intent(in)  :: xparent, yparent
        type(ori),  intent(out) :: offspring
        type(oris) :: opair
        real       :: w(2)
        ! generate stochastic weights
        call ran3arr(w)
        w = w/sum(w)
        ! prepare for geodesic optimisation
        call opair%new(2)
        call opair%set_ori(1,xparent)
        call opair%set_ori(2,yparent)
        ! optimise
        offspring = opair%ori_generator('median', weights=w)
        call opair%kill
    end subroutine rnd_ori_xover

end module simple_crossover
