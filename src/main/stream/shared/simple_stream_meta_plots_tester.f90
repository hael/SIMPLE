!@descr: unit tests for the GUI plots stream p01 fills (simple_stream_meta_plots)
! The windowed time plots and the rate plot, within and beyond a plot's MAX_TIMEPLOT_POINTS: a long
! run widens the window instead of overflowing the plot, a long rate history keeps its newest hours.
! Sub-suite "meta plots" of unit_stream.
module simple_stream_meta_plots_tester
use simple_test_utils
use simple_string,                 only: string
use simple_oris,                   only: oris
use simple_gui_metadata_types,     only: GUI_METADATA_TIMEPLOT_TYPE
use simple_gui_metadata_timeplot,  only: gui_metadata_timeplot, MAX_TIMEPLOT_POINTS
use simple_stream_meta_plots,      only: set_timeplot_from_oris, set_rate_timeplot
implicit none
private
public :: run_all_stream_meta_plots_tests

contains

    subroutine run_all_stream_meta_plots_tests()
        write(*,'(A)') '**** running all stream meta plots tests ****'
        call test_timeplot_windows()
        call test_timeplot_widens()
        call test_rate_timeplot()
        call test_rate_timeplot_newest()
    end subroutine run_all_stream_meta_plots_tests

    !> windows of the given size, the last one partial, labelled 1, 2, ...
    subroutine test_timeplot_windows()
        type(gui_metadata_timeplot) :: meta
        type(oris)                  :: os
        real, allocatable           :: labels(:), data(:), data2(:)
        write(*,'(A)') 'test_timeplot_windows'
        call make_oris(os, 1200)
        call meta%new(GUI_METADATA_TIMEPLOT_TYPE)
        call set_timeplot_from_oris(meta, os, 'ctfres', 500)
        call get_plot(meta, labels, data, data2)
        call assert_int(3, size(labels), 'three windows of 500, 500 and 200 micrographs')
        call assert_real(1.,     labels(1), 1.e-6, 'labelled from 1')
        call assert_real(250.5,  data(1),   1.e-3, 'the first window''s mean')
        call assert_real(1100.5, data(3),   1.e-3, 'the partial last window''s mean')
        call meta%kill()
        call os%kill()
    end subroutine test_timeplot_windows

    !> a run longer than the plot holds widens the window: the whole run, in at most
    !! MAX_TIMEPLOT_POINTS points
    subroutine test_timeplot_widens()
        type(gui_metadata_timeplot) :: meta
        type(oris)                  :: os
        real, allocatable           :: labels(:), data(:), data2(:)
        integer :: n
        write(*,'(A)') 'test_timeplot_widens'
        n = 2 * MAX_TIMEPLOT_POINTS
        call make_oris(os, n)
        call meta%new(GUI_METADATA_TIMEPLOT_TYPE)
        call set_timeplot_from_oris(meta, os, 'ctfres', 1)
        call get_plot(meta, labels, data, data2)
        call assert_int(MAX_TIMEPLOT_POINTS, size(labels), 'the plot is full, not overflowed')
        call assert_real(1.5,              data(1),                   1.e-3, 'each point covers two micrographs')
        call assert_real(real(n) - 0.5,    data(MAX_TIMEPLOT_POINTS), 1.e-3, 'and the last point the run''s end')
        call meta%kill()
        call os%kill()
    end subroutine test_timeplot_widens

    !> one point per rate interval, labelled with its interval
    subroutine test_rate_timeplot()
        type(gui_metadata_timeplot) :: meta
        real, allocatable           :: labels(:), data(:), data2(:)
        write(*,'(A)') 'test_rate_timeplot'
        call meta%new(GUI_METADATA_TIMEPLOT_TYPE)
        call set_rate_timeplot(meta, [10, 20, 30])
        call get_plot(meta, labels, data, data2)
        call assert_int(3, size(labels), 'one point per interval')
        call assert_real(3.,  labels(3), 1.e-6, 'labelled with its interval')
        call assert_real(30., data(3),   1.e-6, 'holding its rate')
        call meta%kill()
    end subroutine test_rate_timeplot

    !> a rate history longer than the plot holds keeps its newest intervals, with their numbers
    subroutine test_rate_timeplot_newest()
        type(gui_metadata_timeplot) :: meta
        real,    allocatable        :: labels(:), data(:), data2(:)
        integer, allocatable        :: rates(:)
        integer :: i, n
        write(*,'(A)') 'test_rate_timeplot_newest'
        n = MAX_TIMEPLOT_POINTS + 10
        rates = [(i, i = 1,n)]
        call meta%new(GUI_METADATA_TIMEPLOT_TYPE)
        call set_rate_timeplot(meta, rates)
        call get_plot(meta, labels, data, data2)
        call assert_int(MAX_TIMEPLOT_POINTS, size(labels), 'the plot is full, not overflowed')
        call assert_real(11.,     labels(1),                   1.e-6, 'the oldest intervals are left out')
        call assert_real(real(n), labels(MAX_TIMEPLOT_POINTS), 1.e-6, 'the newest is kept, with its number')
        call assert_real(real(n), data(MAX_TIMEPLOT_POINTS),   1.e-6, 'and its rate')
        call meta%kill()
    end subroutine test_rate_timeplot_newest

    ! ---- fixtures ------------------------------------------------------------

    ! n micrographs, all selected, with ctfres i for micrograph i
    subroutine make_oris( os, n )
        type(oris), intent(inout) :: os
        integer,    intent(in)    :: n
        integer :: i
        call os%new(n, is_ptcl=.false.)
        do i = 1,n
            call os%set(i, 'state',  1.)
            call os%set(i, 'ctfres', real(i))
        enddo
    end subroutine make_oris

    subroutine get_plot( meta, labels, data, data2 )
        type(gui_metadata_timeplot), intent(in)  :: meta
        real, allocatable,           intent(out) :: labels(:), data(:), data2(:)
        type(string) :: name
        call assert_true(meta%get(name, labels, data, data2), 'the plot is set')
    end subroutine get_plot

end module simple_stream_meta_plots_tester
