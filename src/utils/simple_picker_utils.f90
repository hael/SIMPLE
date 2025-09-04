module simple_picker_utils 
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_image,      only: image
use simple_pickseg,    only: pickseg
use simple_pickgau
implicit none

public :: exec_gaupick, exec_segpick
private
#include "simple_local_flags.inc"

real,    parameter :: SMPD_SHRINK1 = 4.0, SMPD_SHRINK2 = 2.0, BOX_EXP_FAC = 1.0
integer, parameter :: OFFSET       = 3
logical, parameter :: L_WRITE      = .false.
logical, parameter :: L_DEBUG      = .false.

contains

    subroutine exec_gaupick( micname, boxfile_out, smpd, nptcls, pickrefs, dir_out, moldiam_opt, append)
        use simple_strings, only: str2real, parsestr
        character(len=*),           intent(in)    :: micname
        character(len=LONGSTRLEN),  intent(out)   :: boxfile_out
        real,                       intent(in)    :: smpd    !< sampling distance in A
        integer,                    intent(out)   :: nptcls
        class(image),     optional, intent(inout) :: pickrefs(:)
        character(len=*), optional, intent(in)    :: dir_out
        real,             optional, intent(out)   :: moldiam_opt
        logical,          optional, intent(in)    :: append
        type(pickgau)             :: gaup, gaup_refine
        real,         allocatable :: moldiams(:)
        character(len=LONGSTRLEN) :: boxfile
        character(len=4), allocatable :: moldiams_str(:)
        real    :: maxdiam, mmoldiam_opt
        integer :: box, istr, num_entries, idiam, num_moldiams
        logical :: l_roi, l_backgr_subtr, l_append
        l_append = .false.
        if(present(append)) l_append = append
        boxfile = basename(fname_new_ext(trim(micname),'box'))
        if( present(dir_out) ) boxfile = trim(dir_out)//'/'//trim(boxfile)
        l_roi          = trim(params_glob%pick_roi).eq.'yes'
        l_backgr_subtr = l_roi .or. (trim(params_glob%backgr_subtr).eq.'yes')
        call read_mic_raw(micname, smpd, subtr_backgr=l_backgr_subtr)
        if( params_glob%nmoldiams > 1 )then
            moldiams = equispaced_vals(params_glob%moldiam, params_glob%moldiam_max, params_glob%nmoldiams)
            call gaupick_multi(params_glob%pcontrast, SMPD_SHRINK1, moldiams, boxfile, offset=OFFSET, moldiam_opt=mmoldiam_opt, append=l_append)
            if( present(moldiam_opt) ) moldiam_opt = mmoldiam_opt
            deallocate(moldiams)
            nptcls      = 1 ! arbitrary, needs be > 0
            boxfile_out = trim(boxfile)
        else if( params_glob%nmoldiams .eq. 1 .and. params_glob%interactive .eq. 'yes' ) then ! for stream non-interactive initial picking
            allocate(moldiams(1))
            moldiams(1) = params_glob%moldiam
            call gaupick_multi(params_glob%pcontrast, SMPD_SHRINK1, moldiams, boxfile, offset=OFFSET, moldiam_opt=mmoldiam_opt, append=l_append)
            if( present(moldiam_opt) ) moldiam_opt = mmoldiam_opt
            deallocate(moldiams)
            nptcls      = 1 ! arbitrary, needs be > 0
            boxfile_out = trim(boxfile)
        else if( .not. (params_glob%multi_moldiams  .eq. '') )then
            ! multiple moldiam pick that uses multiple gaussians, generates .box file outputs
            istr        = 1
            num_entries = 0
            ! find number of molecular diameters
            do 
                if(params_glob%multi_moldiams(istr:istr) .eq. ' ') then
                    num_entries = num_entries + 1
                    exit
                end if
                if(istr > len(params_glob%multi_moldiams)) exit
                if(params_glob%multi_moldiams(istr:istr) .eq. ',') then
                    num_entries = num_entries + 1
                end if
                istr = istr + 1
            end do
            ! parse moldiams from comma-separated string of numbers into real-valued array
            allocate(moldiams_str(num_entries))
            call parsestr(params_glob%multi_moldiams,',',moldiams_str,num_moldiams)
            allocate(moldiams(num_moldiams))
            do idiam = 1, num_moldiams
                moldiams(idiam) = str2real(moldiams_str(idiam))
            end do
            deallocate(moldiams_str)
            ! execute multiple gaussian pick
            call gaup%new_gaupicker_multi(       params_glob%pcontrast, SMPD_SHRINK1, moldiams, offset=OFFSET, roi=l_roi)
            call gaup_refine%new_gaupicker_multi(params_glob%pcontrast, SMPD_SHRINK2, moldiams, offset=1)
            call gaup%gaupick(gaup_refine)
            maxdiam = maxval(moldiams) + maxval(moldiams) * BOX_EXP_FAC
            box     = find_larger_magic_box(round2even(maxdiam / smpd))
            call gaup_refine%report_boxfile(box, smpd, boxfile, nptcls)
            if( nptcls == 0 )then
                boxfile_out = ''
            else
                boxfile_out = simple_abspath(boxfile)
            endif
            call gaup%kill
            call gaup_refine%kill
            deallocate(moldiams)
        else
            ! single moldiam pick
            if( present(pickrefs) )then
                call gaup%new_refpicker(       params_glob%pcontrast, SMPD_SHRINK1, pickrefs, offset=OFFSET, roi=l_roi)
                call gaup_refine%new_refpicker(params_glob%pcontrast, SMPD_SHRINK2, pickrefs, offset=1)
            else
                call gaup%new_gaupicker(       params_glob%pcontrast, SMPD_SHRINK1, params_glob%moldiam, params_glob%moldiam, offset=OFFSET, roi=l_roi)
                call gaup_refine%new_gaupicker(params_glob%pcontrast, SMPD_SHRINK2, params_glob%moldiam, params_glob%moldiam, offset=1)
            endif
            call gaup%gaupick(gaup_refine)
            ! write
            maxdiam = params_glob%moldiam + params_glob%moldiam * BOX_EXP_FAC
            box     = find_larger_magic_box(round2even(maxdiam / smpd))
            call gaup_refine%report_boxfile(box, smpd, boxfile, nptcls)
            if( nptcls == 0 )then
                boxfile_out = ''
            else
                boxfile_out = simple_abspath(boxfile)
            endif
            call gaup%kill
            call gaup_refine%kill
        endif
    end subroutine exec_gaupick

    subroutine exec_segpick( micname, boxfile_out, nptcls, dir_out, moldiam, winsz )
        character(len=*),           intent(in)  :: micname
        character(len=LONGSTRLEN),  intent(out) :: boxfile_out
        integer,                    intent(out) :: nptcls
        character(len=*), optional, intent(in)  :: dir_out
        real,             optional, intent(in)  :: moldiam
        integer,          optional, intent(in)  :: winsz
        character(len=LONGSTRLEN) :: boxfile
        type(pickseg) :: picker
        boxfile = basename(fname_new_ext(trim(micname),'box'))
        if( present(dir_out) ) boxfile = trim(dir_out)//'/'//trim(boxfile)
        if( present(moldiam) )then
            call picker%pick(micname, moldiam=moldiam)
        elseif( present(winsz) )then
            call picker%pick(micname, winsz=winsz)
        else
            call picker%pick(micname)
        endif
        call picker%report_boxfile(boxfile, nptcls)
        if( nptcls == 0 )then
            boxfile_out = ''
        else
            boxfile_out = simple_abspath(boxfile)
        endif
    end subroutine exec_segpick

end module simple_picker_utils
