!@descr: fills GUI histogram and time-plot metadata from a micrograph or particle segment
!==============================================================================
! MODULE: simple_stream_meta_plots
!
! PURPOSE:
!   The histogram / windowed-statistics / rate plots that stream p01 built
!   inline three times over (ctfres, icefrac, astig, df). Each helper only
!   fills a GUI metadata object; sending it is the caller's business.
!
!   Lives on the stream side rather than in src/utils/gui/metadata because it
!   reads oris: the GUI metadata types import nothing that reaches src/main, so
!   the stages and the master can use them cheaply.
!
!   A plot holds MAX_TIMEPLOT_POINTS points: the windowed plots widen their
!   window to keep the whole run, the rate plot keeps its newest intervals.
!==============================================================================
module simple_stream_meta_plots
use simple_oris,                    only: oris
use simple_string,                  only: string
use simple_histogram,               only: histogram
use simple_gui_metadata_histogram,  only: gui_metadata_histogram
use simple_gui_metadata_timeplot,   only: gui_metadata_timeplot, MAX_TIMEPLOT_POINTS
implicit none

public :: set_histogram_from_oris, set_timeplot_from_oris, set_rate_timeplot, recent_shifts_by_optics_group
private

contains

    !> Counts @p key over every entry of @p os into bins centred on @p bin_centres.
    subroutine set_histogram_from_oris( meta, os, key, bin_centres )
        class(gui_metadata_histogram), intent(inout) :: meta
        class(oris),                   intent(in)    :: os
        character(len=*),              intent(in)    :: key
        real,                          intent(in)    :: bin_centres(:)
        type(histogram)      :: hist
        real,    allocatable :: labels(:)
        integer, allocatable :: counts(:)
        integer :: i
        labels = bin_centres
        call hist%new(labels)
        call hist%zero()
        do i = 1,os%get_noris()
            call hist%update(os%get(i, key))
        enddo
        allocate(counts(size(labels)))
        do i = 1,size(labels)
            counts(i) = int(hist%get(i))
        enddo
        call meta%set(name=string(key), labels=labels, data=counts)
        call hist%kill()
    end subroutine set_histogram_from_oris

    !> Mean and standard deviation of @p key over consecutive windows of @p window entries of @p os,
    !! or of more entries when the windows would not fit in a plot; labelled 1, 2, ...
    subroutine set_timeplot_from_oris( meta, os, key, window )
        class(gui_metadata_timeplot), intent(inout) :: meta
        class(oris),                  intent(inout) :: os
        character(len=*),             intent(in)    :: key
        integer,                      intent(in)    :: window
        real, allocatable :: labels(:), avgs(:), sdevs(:)
        real    :: ave, sdev, var
        integer :: n, win, nwin, iwin, fromto(2)
        logical :: err
        n = os%get_noris()
        if( n == 0 .or. window <= 0 ) return
        win  = max(window, ceiling(real(n) / real(MAX_TIMEPLOT_POINTS)))
        nwin = ceiling(real(n) / real(win))
        allocate(labels(nwin), avgs(nwin), sdevs(nwin))
        do iwin = 1,nwin
            fromto(1) = (iwin - 1) * win + 1
            fromto(2) = min(iwin * win, n)
            call os%stats(key, ave, sdev, var, err, fromto)
            labels(iwin) = real(iwin)
            avgs(iwin)   = ave
            sdevs(iwin)  = sdev
        enddo
        call meta%set(name=string(key), labels=labels, data=avgs, data2=sdevs)
    end subroutine set_timeplot_from_oris

    !> One point per rate interval of a stream_watcher (movies per hour), for the newest
    !! MAX_TIMEPLOT_POINTS intervals, each labelled with its interval's number.
    subroutine set_rate_timeplot( meta, rates )
        class(gui_metadata_timeplot), intent(inout) :: meta
        integer,                      intent(in)    :: rates(:)
        real, allocatable :: labels(:), data(:)
        integer :: i, first, n
        if( size(rates) == 0 ) return
        first = max(1, size(rates) - MAX_TIMEPLOT_POINTS + 1)
        n     = size(rates) - first + 1
        allocate(labels(n), data(n))
        do i = 1,n
            labels(i) = real(first + i - 1)
            data(i)   = real(rates(first + i - 1))
        enddo
        call meta%set(name=string('rate'), labels=labels, data=data)
    end subroutine set_rate_timeplot

    !> Beam-image shifts of the newest micrographs of each optics group of @p os_optics, at most
    !! @p max_points per group, in one pass from the newest micrograph back: column i of
    !! @p xshifts / @p yshifts holds group i, with @p npoints(i) entries filled.
    subroutine recent_shifts_by_optics_group( os_mic, os_optics, max_points, xshifts, yshifts, npoints )
        class(oris),          intent(in)  :: os_mic, os_optics
        integer,              intent(in)  :: max_points
        real,    allocatable, intent(out) :: xshifts(:,:), yshifts(:,:)
        integer, allocatable, intent(out) :: npoints(:)
        integer, allocatable :: ogids(:)
        integer :: ngroups, igroup, imic
        ngroups = os_optics%get_noris()
        allocate(xshifts(max_points,ngroups), yshifts(max_points,ngroups), source=0.)
        allocate(npoints(ngroups), source=0)
        if( ngroups == 0 .or. max_points <= 0 ) return
        allocate(ogids(ngroups))
        do igroup = 1,ngroups
            ogids(igroup) = os_optics%get_int(igroup, 'ogid')
        enddo
        do imic = os_mic%get_noris(),1,-1
            if( all(npoints >= max_points) ) exit
            if( .not. os_mic%isthere(imic, 'ogid') ) cycle
            igroup = findloc(ogids, os_mic%get_int(imic, 'ogid'), 1)
            if( igroup == 0 ) cycle
            if( npoints(igroup) >= max_points ) cycle
            npoints(igroup) = npoints(igroup) + 1
            xshifts(npoints(igroup),igroup) = os_mic%get(imic, 'shiftx')
            yshifts(npoints(igroup),igroup) = os_mic%get(imic, 'shifty')
        enddo
    end subroutine recent_shifts_by_optics_group

end module simple_stream_meta_plots
