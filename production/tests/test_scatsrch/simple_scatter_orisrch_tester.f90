module simple_scatter_orisrch_tester
use simple_scatter_orisrch, only: scatter_orisrch
use simple_oris,            only: oris
use simple_ori,             only: ori
implicit none

public :: exec_scatter_orisrch_test
private

type(scatter_orisrch) :: oscat
type(ori)             :: facit
type(oris)            :: pool
integer               :: nevals

contains

    subroutine exec_scatter_orisrch_test
        !type(oris)         :: subset_div
        type(ori)          :: o, o_prev
        integer            :: npool, i, t, conv_cnt
        integer, parameter :: MAXITS=15
        real               :: score_best, score, score_prev
        call setup_testenv
        ! generate the first pool
        pool       = oscat%gen_initial_pool()
        ! initialise
        score_best = -1.0*sqrt(2.0)
        nevals     = 0
        conv_cnt = 0
        do t=1,MAXITS  
            ! score pool
            npool = pool%get_noris()
            do i=1,npool
                o     = pool%get_ori(i)
                score = scorefun(o)
                call pool%set(i, 'corr', score)
            end do
            ! generate reference sets from pool
            call oscat%gen_ref_sets_from_pool(pool)
            call pool%kill
            o_prev = o
            o = oscat%get_best_quality_ref_ori()
            score_prev = score_best
            score_best = o%get('corr')
            write(*,'(a,1x,i2,1x,a,1x,f7.3,1x,a,1x,i3,1x)') 'ITERATION:', t,&
            'SCORE:', score_best, 'NOPT:', oscat%get_nopt()
            if( abs(score_best-score_prev) < 1e-4 )then
                conv_cnt = conv_cnt + 1
            else
                conv_cnt = 0
            endif
            if( conv_cnt == 3 ) exit
            pool = oscat%update_pool()
        end do
        write(*,*) 'NEVALS: ', nevals
    end subroutine exec_scatter_orisrch_test

    subroutine setup_testenv
        call facit%rnd_euler
        call oscat%new
    end subroutine setup_testenv

    real function scorefun( otrial )
        class(ori), intent(in) :: otrial
        scorefun = -(otrial.geod.facit)
        nevals = nevals + 1
    end function scorefun

end module simple_scatter_orisrch_tester
