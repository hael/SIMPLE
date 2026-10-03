!@descr: GUI messages several stream stages send: class averages or picking references as sprite-sheet tiles, and the latest micrographs
!==============================================================================
! MODULE: simple_stream_gui_senders
!
! PURPOSE:
!   The initial-analysis stage and the reference-picking stage both show the
!   GUI a set of class averages (or picking references) as tiles of one
!   sprite sheet, and the thumbnails of the most recent micrographs with their
!   CTF values and picked positions; the 3D stage shows each state's three
!   orthogonal reprojections as tiles. These are those messages, built from
!   the stage's metadata object and sent through its pipe; none keeps state
!   between calls.
!
!   It could move beside the metadata types in src/utils/gui, which would then
!   depend on simple_stream_pipe.
!==============================================================================
module simple_stream_gui_senders
use simple_defs,                    only: logfhandle
use simple_string,                  only: string
use simple_fileio,                  only: file_exists
use simple_nrtxtfile,               only: nrtxtfile
use simple_oris,                    only: oris
use simple_gui_metadata_cavg2D,     only: gui_metadata_cavg2D, sprite_sheet_pos
use simple_gui_metadata_micrograph, only: gui_metadata_micrograph, MAX_MIC_COORDINATES
use simple_stream_pipe,             only: stream_pipe
implicit none

public :: send_cavgs, send_recent_micrographs, send_reproj_tiles
private

contains

    !> One message per class average @p inds(i) (or picking reference) as tile i of the sprite
    !! sheet @p jpg, row by row in an @p xtiles x @p ytiles grid; @p stk is the stack the tiles
    !! come from. The messages carry each class's resolution and population when given, either
    !! from @p os_cls2D (indexed by class) or as @p res and @p pop (indexed like @p inds).
    subroutine send_cavgs( pipe, meta, jpg, inds, stk, xtiles, ytiles, os_cls2D, res, pop )
        class(stream_pipe),        intent(inout) :: pipe
        type(gui_metadata_cavg2D), intent(inout) :: meta
        class(string),             intent(in)    :: jpg, stk
        integer,                   intent(in)    :: inds(:), xtiles, ytiles
        class(oris), optional,     intent(in)    :: os_cls2D
        real,        optional,     intent(in)    :: res(:)
        integer,     optional,     intent(in)    :: pop(:)
        type(string) :: jpg_here, stk_here
        integer      :: n, i, idx, xtile, ytile
        n        = size(inds)
        jpg_here = jpg
        stk_here = stk
        write(logfhandle,*) '>>> SENDING', n, ' CLASS AVERAGES TO GUI'
        xtile = 0
        ytile = 0
        do i = 1, n
            idx = inds(i)
            if( present(os_cls2D) )then
                call meta%set(path=jpg_here, mrcpath=stk_here, i=i, i_max=n, res=os_cls2D%get(idx, 'res'),&
                    &pop=os_cls2D%get_int(idx, 'pop'), idx=idx, sprite=tile_pos(xtile, ytile))
            else if( present(res) .and. present(pop) )then
                call meta%set(path=jpg_here, mrcpath=stk_here, i=i, i_max=n, res=res(i), pop=pop(i), idx=idx,&
                    &sprite=tile_pos(xtile, ytile))
            else
                call meta%set(path=jpg_here, mrcpath=stk_here, i=i, i_max=n, idx=idx, sprite=tile_pos(xtile, ytile))
            endif
            call pipe%send_meta(meta)
            xtile = xtile + 1
            if( xtile == xtiles )then
                xtile = 0
                ytile = ytile + 1
            endif
        end do

    contains

        ! a single row or column sits at 0%; the step is computed only when there are two or more
        ! tiles, as merge() would evaluate the division by zero (trapped in Debug builds)
        type(sprite_sheet_pos) function tile_pos( xt, yt )
            integer, intent(in) :: xt, yt
            real :: xstep, ystep
            xstep = 0.
            ystep = 0.
            if( xtiles > 1 ) xstep = 100.0 / real(xtiles - 1)
            if( ytiles > 1 ) ystep = 100.0 / real(ytiles - 1)
            tile_pos = sprite_sheet_pos(x = real(xt) * xstep, y = real(yt) * ystep, h = 100 * ytiles, w = 100 * xtiles)
        end function tile_pos

    end subroutine send_cavgs

    !> The three orthogonal reprojections of state @p istate (of @p nstates) as tiles 1-3 of the
    !! one-row sprite sheet @p jpg of volume @p vol, one message each: index @p istate, position
    !! (istate - 1) * 3 + tile of 3 * @p nstates, and the state's population @p pop.
    subroutine send_reproj_tiles( pipe, meta, jpg, vol, istate, nstates, pop )
        class(stream_pipe),        intent(inout) :: pipe
        type(gui_metadata_cavg2D), intent(inout) :: meta
        class(string),             intent(in)    :: jpg, vol
        integer,                   intent(in)    :: istate, nstates, pop
        integer, parameter :: NTILES = 3
        type(string) :: jpg_here, vol_here
        integer      :: itile
        jpg_here = jpg
        vol_here = vol
        do itile = 1,NTILES
            call meta%set(path=jpg_here, mrcpath=vol_here, idx=istate, i=(istate - 1) * NTILES + itile,&
                &i_max=nstates * NTILES, pop=pop, sprite=sprite_sheet_pos(x=real(itile - 1) * (100.0 / real(NTILES - 1)),&
                &y=0.0, h=100, w=100 * NTILES))
            call pipe%send_meta(meta)
        enddo
    end subroutine send_reproj_tiles

    !> Thumbnail, CTF values and picked positions of the @p nmax most recent micrographs of
    !! @p os_mic, one message each; nothing when the segment has no thumbnails or box files yet.
    !! A micrograph whose box file is missing is sent without positions.
    subroutine send_recent_micrographs( pipe, meta, os_mic, nmax )
        class(stream_pipe),            intent(inout) :: pipe
        type(gui_metadata_micrograph), intent(inout) :: meta
        class(oris),                   intent(in)    :: os_mic
        integer,                       intent(in)    :: nmax
        type(nrtxtfile)   :: boxfile
        type(string)      :: boxpath
        real, allocatable :: boxdata(:)
        integer :: nmics, i_max, iori, imic, i, nrecs, nlines, xdim, ydim
        if( .not. (os_mic%isthere('thumb') .and. os_mic%isthere('xdim') .and. os_mic%isthere('ydim')&
            &.and. os_mic%isthere('smpd') .and. os_mic%isthere('boxfile')) ) return
        nmics = os_mic%get_noris()
        i_max = min(nmics, nmax)
        do iori = 1, i_max
            imic = nmics - i_max + iori
            call meta%set(path=os_mic%get_str(imic, 'thumb'), dfx=os_mic%get(imic, 'dfx'), dfy=os_mic%get(imic, 'dfy'),&
                &ctfres=os_mic%get(imic, 'ctfres'), i=iori, i_max=i_max)
            call meta%clear_coordinates()
            boxpath = os_mic%get_str(imic, 'boxfile')
            if( boxpath%strlen() > 0 .and. file_exists(boxpath) )then
                call boxfile%new(boxpath, 1)
                nrecs  = boxfile%get_nrecs_per_line()
                nlines = boxfile%get_ndatalines()
                xdim   = os_mic%get_int(imic, 'xdim')
                ydim   = os_mic%get_int(imic, 'ydim')
                if( nrecs >= 4 )then
                    allocate(boxdata(nrecs))
                    ! the thumbnail shows at most the picks the metadata holds
                    do i = 1, min(nlines, MAX_MIC_COORDINATES)
                        call boxfile%readNextDataLine(boxdata)
                        call meta%set_coordinate(i, nint(boxdata(1) + boxdata(3)/2), nint(boxdata(2) + boxdata(4)/2), xdim, ydim)
                    enddo
                    deallocate(boxdata)
                endif
                call boxfile%kill()
            endif
            call pipe%send_meta(meta)
        end do
    end subroutine send_recent_micrographs

end module simple_stream_gui_senders
