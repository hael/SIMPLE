program simple_test_tab_align
include 'simple_lib.f08'
use simple_cmdline,    only: cmdline
use simple_parameters, only: parameters
implicit none
type(parameters)   :: p
type(cmdline)      :: cline
integer, parameter :: N_P = 100000, N_R = 10
integer :: sorted_inds(N_R), assigned_ptcl, assigned_iref, iref
logical :: ptcl_avail(N_P)
real    :: tab(N_P, N_R), ref_dist(N_R), sorted_dist(N_R)
if( command_argument_count() < 2 )then
    write(logfhandle,'(a)') 'Usage: simple_test_tab_align smpd=xx nthr=yy'
    stop
else
    call cline%parse_oldschool
endif
call cline%checkvar('smpd', 1)
call cline%checkvar('nthr', 1)
call cline%check
call p%new(cline)
call seed_rnd
call random_number(tab)
ptcl_avail = .true.
do while( any(ptcl_avail) )
    ref_dist = 0.
    !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iref)
    do iref = 1, N_R
        ref_dist(iref) = minval(tab(:,iref), dim=1, mask=ptcl_avail)
    enddo
    !$omp end parallel do
    assigned_iref = ref_multinomal(ref_dist)
    assigned_ptcl = minloc(tab(:,assigned_iref), dim=1, mask=ptcl_avail)
    ptcl_avail(assigned_ptcl) = .false.
enddo
contains
    function ref_multinomal( pvec ) result( which )
        real, intent(in) :: pvec(:) !< probabilities
        integer :: i, which
        real    :: rnd, bound, sum_refs_corr
        rnd         = ran3()
        sorted_dist = pvec
        sorted_inds = (/(i,i=1,N_R)/)
        call hpsort(sorted_dist, sorted_inds)
        sum_refs_corr = sum(sorted_dist)
        if( sum_refs_corr < TINY )then
            ! uniform sampling
            which = 1 + floor(real(N_R) * rnd)
        else
            ! normalizing within the hard-limit
            sorted_dist = sorted_dist / sum_refs_corr
            bound = 0.
            do which=1,N_R-1
                bound = bound + sorted_dist(which)
                if( rnd >= bound )exit
            enddo
        endif
        which = sorted_inds(which)
    end function ref_multinomal
end program simple_test_tab_align
