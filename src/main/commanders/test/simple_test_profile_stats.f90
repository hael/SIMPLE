!@descr: robustness-profile statistics: the gate tables (metrics.tsv) of many executions summarized per metric, and two profiles compared
! A metric's distribution over the executions of a profile: runs, median, robust spread s (1.4826
! times the median absolute deviation, or the range over 4 when that is zero), the observed range,
! the worst run in the sense of its floor, and the failure rate: the fraction of all executions in
! which it failed or was never written (a failed or timed-out execution counts).
module simple_test_profile_stats
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
use simple_stat, only: median, mad_gau
implicit none

public :: gate_row, metric_dist, read_gate_table, robust_spread, summarize_metric, profile_table
public :: write_profile, read_profile, compare_profiles
private

integer, parameter :: NAMELEN = 160

!> one row of a gate table
type :: gate_row
    character(len=NAMELEN) :: name  = ''
    real                   :: value = 0.
    real                   :: floor = 0.
    logical                :: pass  = .false.
    character(len=8)       :: sense = ''  !< min, max, check or report
end type gate_row

!> one metric over the executions of a profile
type :: metric_dist
    character(len=NAMELEN) :: name   = ''
    character(len=8)       :: sense  = ''
    integer :: runs         = 0   !< executions that wrote it
    integer :: nexec        = 0   !< executions of the profile
    real    :: floor        = 0.  !< median of the floors written
    real    :: med          = 0.
    real    :: spread       = 0.
    real    :: vmin         = 0.
    real    :: vmax         = 0.
    real    :: worst        = 0.
    real    :: failure_rate = 0.
end type metric_dist

contains

    !> the rows of a metrics.tsv; none when the file is absent or unreadable
    subroutine read_gate_table( fname, rows )
        character(len=*),            intent(in)  :: fname
        type(gate_row), allocatable, intent(out) :: rows(:)
        type(gate_row), allocatable :: tmp(:)
        character(len=1024) :: line
        character(len=NAMELEN) :: fields(5)
        integer :: funit, ios, n, nf
        logical :: exists
        allocate(rows(0))
        inquire(file=fname, exist=exists)
        if( .not. exists ) return
        open(newunit=funit, file=fname, action='read', status='old', iostat=ios)
        if( ios /= 0 ) return
        read(funit,'(a)',iostat=ios) line ! header
        allocate(tmp(64))
        n = 0
        do
            read(funit,'(a)',iostat=ios) line
            if( ios /= 0 ) exit
            if( len_trim(line) == 0 ) cycle
            call split_tabs(line, fields, nf)
            if( nf < 4 ) cycle
            if( n == size(tmp) ) tmp = [tmp, tmp]
            n = n + 1
            tmp(n)%name = fields(1)
            read(fields(2),*,iostat=ios) tmp(n)%value
            if( ios /= 0 ) tmp(n)%value = ieee_value(0., ieee_quiet_nan) ! an unreadable value does not pass
            read(fields(3),*,iostat=ios) tmp(n)%floor
            tmp(n)%pass  = trim(adjustl(fields(4))) == 'yes'
            tmp(n)%sense = 'min'
            if( nf >= 5 ) tmp(n)%sense = fields(5)(1:len(tmp(n)%sense))
        enddo
        close(funit)
        rows = tmp(1:n)
    end subroutine read_gate_table

    !> 1.4826 times the median absolute deviation; the range over 4 when that is zero
    real function robust_spread( vals )
        real, intent(in) :: vals(:)
        robust_spread = 0.
        if( size(vals) < 2 ) return
        robust_spread = mad_gau(vals, median(vals))
        if( robust_spread <= 0. ) robust_spread = (maxval(vals) - minval(vals)) / 4.
    end function robust_spread

    !> the distribution of one metric: vals over its rows and passes over the executions that wrote it
    !! (a non-finite value has failed, and stays out of the statistics), nexec executions in all
    function summarize_metric( name, sense, vals, floors, passes, nexec ) result( d )
        character(len=*), intent(in) :: name, sense
        real,             intent(in) :: vals(:), floors(:)
        logical,          intent(in) :: passes(:)
        integer,          intent(in) :: nexec
        type(metric_dist) :: d
        logical, allocatable :: finite(:)
        real,    allocatable :: v(:)
        d%name  = name
        d%sense = sense
        d%runs  = size(passes)
        d%nexec = nexec
        finite  = ieee_is_finite(vals)
        d%failure_rate = 1.
        if( nexec > 0 ) d%failure_rate = real(nexec - count(passes)) / real(nexec)
        if( size(floors) > 0 ) d%floor = median(floors)
        v = pack(vals, finite)
        if( size(v) == 0 ) return
        d%med    = median(v)
        d%spread = robust_spread(v)
        d%vmin   = minval(v)
        d%vmax   = maxval(v)
        if( trim(sense) == 'max' )then
            d%worst = d%vmax
        else
            d%worst = d%vmin
        endif
    end function summarize_metric

    !> the distribution of every metric of the gate tables of nexec executions (one table per file;
    !! a missing file is an execution that wrote nothing), in order of first appearance; a metric
    !! written several times by one execution passes there when every row passes and the execution
    !! itself passed (exec_ok: it ended with exit code 0 inside its deadline), and takes the sense of
    !! its bounded rows
    subroutine profile_table( fnames, exec_ok, dists )
        character(len=*),               intent(in)  :: fnames(:)
        logical,                        intent(in)  :: exec_ok(:)
        type(metric_dist), allocatable, intent(out) :: dists(:)
        type(gate_row), allocatable :: rows(:), allrows(:)
        character(len=NAMELEN), allocatable :: names(:)
        character(len=8),       allocatable :: senses(:)
        integer,                allocatable :: iexec(:)
        logical, allocatable :: sel(:), bounded(:), passes(:)
        integer :: i, j, k, n
        allocate(allrows(0), names(0), senses(0), iexec(0))
        do i = 1, size(fnames)
            call read_gate_table(trim(fnames(i)), rows)
            allrows = [allrows, rows]
            iexec   = [iexec, spread(i, 1, size(rows))]
            do j = 1, size(rows)
                k = findloc(names, rows(j)%name, dim=1)
                if( k == 0 )then
                    names  = [names,  rows(j)%name]
                    senses = [senses, rows(j)%sense]
                ! a value reported and then held to a floor takes the sense of the floor
                else if( trim(senses(k)) == 'report' )then
                    senses(k) = rows(j)%sense
                endif
            enddo
        enddo
        n = size(names)
        allocate(dists(n))
        do j = 1, n
            sel = allrows%name == names(j)
            ! one verdict per execution that wrote the metric: its rows and its own outcome
            passes = [(exec_ok(k) .and. all(allrows%pass .or. .not. (sel .and. iexec == k)), k = 1, size(fnames))]
            passes = pack(passes, [(any(sel .and. iexec == k), k = 1, size(fnames))])
            ! the floor of a reported and gated name is that of its bounded rows
            bounded = sel .and. (allrows%sense == 'min' .or. allrows%sense == 'max')
            if( .not. any(bounded) ) bounded = sel
            dists(j) = summarize_metric(trim(names(j)), trim(senses(j)), pack(allrows%value, sel), &
                &pack(allrows%floor, bounded), passes, size(fnames))
        enddo
    end subroutine profile_table

    subroutine write_profile( fname, dists )
        character(len=*),  intent(in) :: fname
        type(metric_dist), intent(in) :: dists(:)
        character(len=1), parameter :: T = achar(9)
        integer :: funit, i
        open(newunit=funit, file=fname, action='write', status='replace')
        write(funit,'(a)') 'name'//T//'sense'//T//'floor'//T//'runs'//T//'executions'//T//'median'//T//'spread'//T//&
            &'min'//T//'max'//T//'worst'//T//'failure_rate'
        do i = 1, size(dists)
            write(funit,'(a,a,a,a,es13.6,a,i0,a,i0,5(a,es13.6),a,f8.4)') trim(dists(i)%name), T, trim(dists(i)%sense), T, &
                &dists(i)%floor, T, dists(i)%runs, T, dists(i)%nexec, T, dists(i)%med, T, dists(i)%spread, T, &
                &dists(i)%vmin, T, dists(i)%vmax, T, dists(i)%worst, T, dists(i)%failure_rate
        enddo
        close(funit)
    end subroutine write_profile

    !> a table written by write_profile
    subroutine read_profile( fname, dists )
        character(len=*),               intent(in)  :: fname
        type(metric_dist), allocatable, intent(out) :: dists(:)
        type(metric_dist)      :: d
        character(len=2048)    :: line
        character(len=NAMELEN) :: f(11)
        integer :: funit, ios, nf
        allocate(dists(0))
        open(newunit=funit, file=fname, action='read', status='old', iostat=ios)
        if( ios /= 0 ) return
        read(funit,'(a)',iostat=ios) line
        do
            read(funit,'(a)',iostat=ios) line
            if( ios /= 0 ) exit
            call split_tabs(line, f, nf)
            if( nf < 11 ) cycle
            d%name  = f(1)
            d%sense = f(2)(1:len(d%sense))
            read(f(3),*)  d%floor
            read(f(4),*)  d%runs
            read(f(5),*)  d%nexec
            read(f(6),*)  d%med
            read(f(7),*)  d%spread
            read(f(8),*)  d%vmin
            read(f(9),*)  d%vmax
            read(f(10),*) d%worst
            read(f(11),*) d%failure_rate
            dists = [dists, d]
        enddo
        close(funit)
    end subroutine read_profile

    !> the metrics of profile a and profile b side by side: the shift of the median in units of the
    !! larger spread (signed so that positive is better), and both failure rates
    subroutine compare_profiles( a, b, fname )
        type(metric_dist), intent(in) :: a(:), b(:)
        character(len=*),  intent(in) :: fname
        character(len=1), parameter :: T = achar(9)
        real    :: s, shift
        integer :: funit, i, j
        open(newunit=funit, file=fname, action='write', status='replace')
        write(funit,'(a)') 'name'//T//'sense'//T//'median_a'//T//'median_b'//T//'shift_in_spreads'//T//&
            &'failure_rate_a'//T//'failure_rate_b'
        do i = 1, size(a)
            j = findloc(b%name, a(i)%name, dim=1)
            if( j == 0 ) cycle
            s     = max(a(i)%spread, b(j)%spread)
            shift = 0.
            if( s > 0. ) shift = (b(j)%med - a(i)%med) / s
            if( trim(a(i)%sense) == 'max' .and. shift /= 0. ) shift = -shift
            write(funit,'(a,a,a,3(a,es13.6),2(a,f8.4))') trim(a(i)%name), T, trim(a(i)%sense), T, a(i)%med, T, &
                &b(j)%med, T, shift, T, a(i)%failure_rate, T, b(j)%failure_rate
        enddo
        close(funit)
    end subroutine compare_profiles

    subroutine split_tabs( line, fields, nf )
        character(len=*), intent(in)  :: line
        character(len=*), intent(out) :: fields(:)
        integer,          intent(out) :: nf
        integer :: first, i
        fields = ''
        nf     = 0
        first  = 1
        do i = 1, len_trim(line) + 1
            if( i > len_trim(line) .or. line(min(i,len(line)):min(i,len(line))) == achar(9) )then
                if( nf == size(fields) ) return
                nf = nf + 1
                fields(nf) = adjustl(line(first:i-1))
                first = i + 1
            endif
        enddo
    end subroutine split_tabs

end module simple_test_profile_stats
