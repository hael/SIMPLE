program simple_test_ref_assign
include 'simple_lib.f08'
use simple_cmdline,    only: cmdline
use simple_parameters, only: parameters
implicit none
type(parameters)   :: p
type(cmdline)      :: cline
integer, parameter :: N_P = 6000, N_R = 100
integer(8) :: t_cur
integer :: sorted_inds(N_R), assigned_ptcl, assigned_iref, iref, iptcl
integer :: stab_inds(N_P, N_R), ref_dist_inds(N_R), cnt, cur_map(N_P), imp_map(N_P), sorted_vec_inds(N_P)
logical :: ptcl_avail(N_P)
real    :: tab(N_P, N_R), ref_dist(N_R), sorted_dist(N_R), sorted_tab(N_P, N_R), rnd_list(N_P), sorted_vec(N_P)
if( command_argument_count() < 2 )then
    write(logfhandle,'(a)') 'ERROR! Usage: simple_test_ref_assign smpd=xx nthr=yy'
    write(logfhandle,'(a)') 'Example: simple_test_ref_assign smpd=1. nthr=4'
    write(logfhandle,'(a)') 'DEFAULT TEST (example above) is running now...'
else
    call cline%parse_oldschool
endif
call cline%checkvar('smpd', 1)
call cline%checkvar('nthr', 2)
call cline%check
call p%new(cline)
call seed_rnd
call random_number(tab)
do iptcl = 1, N_P
    rnd_list(iptcl) = ran3()
enddo
ptcl_avail = .true.
cnt        = 1
t_cur      = tic()
do while( any(ptcl_avail) )
    ref_dist = 0.
    !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iref)
    do iref = 1, N_R
        ref_dist(iref) = minval(tab(:,iref), dim=1, mask=ptcl_avail)
    enddo
    !$omp end parallel do
    assigned_iref = ref_multinomal(ref_dist, cnt)
    assigned_ptcl = minloc(tab(:,assigned_iref), dim=1, mask=ptcl_avail)
    ptcl_avail(assigned_ptcl) = .false.
    cur_map(assigned_ptcl)    = assigned_iref
    cnt = cnt + 1
enddo
print *, 'current timing = ', toc(t_cur)
print *, '-----------'
t_cur = tic()
! sorting each columns
sorted_tab = tab
!$omp parallel do default(shared) proc_bind(close) schedule(static) private(iref,iptcl)
do iref = 1, N_R
    stab_inds(:,iref) = (/(iptcl,iptcl=1,N_P)/)
    call hpsort(sorted_tab(:,iref), stab_inds(:,iref))
enddo
!$omp end parallel do
ptcl_avail    = .true.
ref_dist_inds = 1
ref_dist      = sorted_tab(1,:)
cnt           = 1
do while( any(ptcl_avail) )
    assigned_iref = ref_multinomal(ref_dist, cnt)
    assigned_ptcl = stab_inds(ref_dist_inds(assigned_iref), assigned_iref)
    ptcl_avail(assigned_ptcl) = .false.
    imp_map(assigned_ptcl)    = assigned_iref
    ! update the ref_dist and ref_dist_inds
    !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iref)
    do iref = 1, N_R
        do while( ref_dist_inds(iref) < N_P .and. .not.(ptcl_avail(stab_inds(ref_dist_inds(iref), iref))) )
            ref_dist_inds(iref) = ref_dist_inds(iref) + 1
            ref_dist(iref)      = sorted_tab(ref_dist_inds(iref), iref)
        enddo
    enddo
    !$omp end parallel do
    cnt = cnt + 1
enddo
print *, 'improved timing = ', toc(t_cur)
print *, '-----------'
if( all(cur_map == imp_map) ) print *, 'mappings match! PASSED 1!'
print *, '-----------'
imp_map = 0
t_cur = tic()
! sorting each columns
sorted_tab = tab
!$omp parallel do default(shared) proc_bind(close) schedule(static) private(iref,iptcl)
do iref = 1, N_R
    stab_inds(:,iref) = (/(iptcl,iptcl=1,N_P)/)
    call hpsort(sorted_tab(:,iref), stab_inds(:,iref))
enddo
!$omp end parallel do
ptcl_avail    = .true.
ref_dist_inds = 1
ref_dist      = sorted_tab(1,:)
cnt           = 1
do while( any(ptcl_avail) )
    assigned_iref = ref_multinomal2(ref_dist, cnt)
    assigned_ptcl = stab_inds(ref_dist_inds(assigned_iref), assigned_iref)
    ptcl_avail(assigned_ptcl) = .false.
    imp_map(assigned_ptcl)    = assigned_iref
    ! update the ref_dist and ref_dist_inds
    !!$omp parallel do default(shared) proc_bind(close) schedule(static) private(iref)
    do iref = 1, N_R
        do while( ref_dist_inds(iref) < N_P .and. .not.(ptcl_avail(stab_inds(ref_dist_inds(iref), iref))) )
            ref_dist_inds(iref) = ref_dist_inds(iref) + 1
            ref_dist(iref)      = sorted_tab(ref_dist_inds(iref), iref)
        enddo
    enddo
    !!$omp end parallel do
    cnt = cnt + 1
enddo
print *, 'improved timing = ', toc(t_cur)
print *, '-----------'
if( all(cur_map == imp_map) ) print *, 'mappings match! PASSED 2!'

contains
    function ref_multinomal( pvec, cnt ) result( which )
        real,    intent(in) :: pvec(:) !< probabilities
        integer, intent(in) :: cnt
        integer :: i, which
        real    :: rnd, bound, sum_refs_corr
        rnd         = rnd_list(cnt)
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

    function ref_multinomal2( pvec, cnt ) result( which )
        real,    intent(in) :: pvec(:) !< probabilities
        integer, intent(in) :: cnt
        integer :: i, which
        real    :: rnd, bound, sum_refs_corr
        rnd           = rnd_list(cnt)
        sorted_dist   = pvec
        sum_refs_corr = sum(sorted_dist)
        if( sum_refs_corr < TINY )then
            ! uniform sampling
            which = 1 + floor(real(N_R) * rnd)
        else
            rnd   = rnd * sum_refs_corr
            which = minloc(sorted_dist,dim=1)
            bound = sorted_dist(which)
            if( rnd >= bound )then
                ! done
            else
                sorted_inds = (/(i,i=1,N_R)/)
                call hpsort(sorted_dist, sorted_inds)
                bound = 0.
                do which=1,N_R-1
                    bound = bound + sorted_dist(which)
                    if( rnd >= bound )exit
                enddo
                which = sorted_inds(which)
            endif
        endif
    end function ref_multinomal2

end program simple_test_ref_assign
