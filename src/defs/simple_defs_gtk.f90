!@descr: settings of SIMPLE's GTK code: colour theme, fonts, symbols, key codes, display limits and widget geometry
! Shared by src/main/gtk and the image-stack display model; not re-exported by simple_core_module_api, so a
! change here recompiles only those. Names carry GTKUI_ (GUI_ is NICE's, GTK_ is gtk-fortran's).
module simple_defs_gtk
use, intrinsic :: iso_c_binding, only: c_int
use simple_defs,                 only: dp
implicit none

! colour theme (RGB in 0..1): text, histogram and limit colours from select2d.py, background and accent chosen
real(dp), parameter :: GTKUI_COLOR_BG(3)       = [246._dp, 246._dp, 246._dp] / 255._dp
real(dp), parameter :: GTKUI_COLOR_FG(3)       = [ 32._dp,  36._dp,  43._dp] / 255._dp
real(dp), parameter :: GTKUI_COLOR_ACCENT(3)   = [ 53._dp, 132._dp, 228._dp] / 255._dp
real(dp), parameter :: GTKUI_COLOR_WHITE(3)    = [1._dp, 1._dp, 1._dp]
real(dp), parameter :: GTKUI_COLOR_BACKDROP(3) = [ 32._dp,  34._dp,  38._dp] / 255._dp
real(dp), parameter :: GTKUI_COLOR_BARS(3)     = [116._dp, 139._dp, 167._dp] / 255._dp
real(dp), parameter :: GTKUI_COLOR_MIN(3)      = [ 51._dp, 120._dp, 189._dp] / 255._dp
real(dp), parameter :: GTKUI_COLOR_MAX(3)      = [192._dp, 101._dp,  50._dp] / 255._dp
real(dp), parameter :: GTKUI_COLOR_RANGE(3)    = [139._dp,  74._dp, 158._dp] / 255._dp
real(dp), parameter :: GTKUI_COLOR_SURFACE(3)  = GTKUI_COLOR_BG + (GTKUI_COLOR_FG - GTKUI_COLOR_BG) * 0.025_dp
real(dp), parameter :: GTKUI_COLOR_BORDER(3)   = GTKUI_COLOR_BG + (GTKUI_COLOR_FG - GTKUI_COLOR_BG) * 0.17_dp
real(dp), parameter :: GTKUI_COLOR_HOVER(3)    = GTKUI_COLOR_BG + (GTKUI_COLOR_FG - GTKUI_COLOR_BG) * 0.38_dp
real(dp), parameter :: GTKUI_COLOR_MUTED(3)    = GTKUI_COLOR_SURFACE + (GTKUI_COLOR_FG - GTKUI_COLOR_SURFACE) * 0.65_dp
real(dp), parameter :: GTKUI_COLOR_SELECTED(3) = GTKUI_COLOR_SURFACE + (GTKUI_COLOR_ACCENT - GTKUI_COLOR_SURFACE) * 0.09_dp

! fonts: the family, Tk point sizes of the script in pixels, small print, Pango's units per pixel
character(len=*), parameter :: GTKUI_FONT        = 'Sans'
real(dp),         parameter :: GTKUI_PT2PX       = 1.25_dp
real(dp),         parameter :: GTKUI_SMALL_PX    = 10._dp
real(dp),         parameter :: GTKUI_MESSAGE_PX  = 15._dp
real(dp),         parameter :: GTKUI_PANGO_UNITS = 1024._dp

! UTF-8 of symbols in labels
character(len=*), parameter :: GTKUI_U_ELLIPSIS = char(226)//char(128)//char(166)
character(len=*), parameter :: GTKUI_U_MIDDOT   = char(194)//char(183)
character(len=*), parameter :: GTKUI_U_TIMES    = char(195)//char(151)
character(len=*), parameter :: GTKUI_U_UP       = char(226)//char(134)//char(145)
character(len=*), parameter :: GTKUI_U_DOWN     = char(226)//char(134)//char(147)
character(len=*), parameter :: GTKUI_U_LEFT     = char(226)//char(134)//char(144)
character(len=*), parameter :: GTKUI_U_RIGHT    = char(226)//char(134)//char(146)
character(len=*), parameter :: GTKUI_U_CMD      = char(226)//char(140)//char(152)

! GDK key codes gtk-fortran does not provide (gdkkeysyms.h)
integer(c_int), parameter :: GTKUI_KEY_F10   = int(z'FFC7', c_int)
integer(c_int), parameter :: GTKUI_KEY_MENU  = int(z'FF67', c_int)
integer(c_int), parameter :: GTKUI_KEY_LEFT  = int(z'FF51', c_int)
integer(c_int), parameter :: GTKUI_KEY_UP    = int(z'FF52', c_int)
integer(c_int), parameter :: GTKUI_KEY_RIGHT = int(z'FF53', c_int)
integer(c_int), parameter :: GTKUI_KEY_DOWN  = int(z'FF54', c_int)

! image-stack display (select2d.py): brightness in percent, contrast factor, automatic-window percentiles, bins
real(dp), parameter :: GTKUI_BRIGHTNESS_LIM(2) = [-100._dp, 100._dp]
real(dp), parameter :: GTKUI_CONTRAST_LIM(2)   = [0.1_dp, 5._dp]
real(dp), parameter :: GTKUI_AUTO_PCT(2)       = [1._dp, 99.9_dp]
integer,  parameter :: GTKUI_HIST_NBINS        = 256

! image grid (select2d.py): margins and spacing, image border, tile sizes, zoom per wheel step
integer,  parameter :: GTKUI_GRID_MARGIN    = 12
integer,  parameter :: GTKUI_GRID_SPACING   = 12
integer,  parameter :: GTKUI_GRID_BORDER    = 2
integer,  parameter :: GTKUI_TILE_START     = 240
integer,  parameter :: GTKUI_PICK_BOX_START = 64
integer,  parameter :: GTKUI_TILE_MIN       = 50
integer,  parameter :: GTKUI_TILE_MAX       = 1000
real(dp), parameter :: GTKUI_ZOOM_FACTOR    = 1.15_dp

! display panel (select2d.py): histogram size and the pick distance of its limit lines (pixels)
integer, parameter :: GTKUI_HIST_W      = 400
integer, parameter :: GTKUI_HIST_H      = 175
integer, parameter :: GTKUI_HIT_RADIUS  = 7
integer, parameter :: GTKUI_SLIDER_W    = 175

end module simple_defs_gtk
