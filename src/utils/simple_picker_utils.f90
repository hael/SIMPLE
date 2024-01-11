module simple_picker_utils 
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_image,      only: image
use simple_pickgau
implicit none

public :: exec_gaupick
private
#include "simple_local_flags.inc"

real,    parameter :: SMPD_SHRINK1 = 4.0, SMPD_SHRINK2 = 2.0, BOX_EXP_FAC = 0.333
integer, parameter :: OFFSET       = 3
logical, parameter :: L_WRITE      = .false.
logical, parameter :: L_DEBUG      = .false.

contains

    subroutine exec_gaupick( micname, boxfile_out, smpd, nptcls, pickrefs )
        character(len=*),          intent(in)    :: micname
        character(len=LONGSTRLEN), intent(out)   :: boxfile_out
        real,                      intent(in)    :: smpd    !< sampling distance in A
        integer,                   intent(out)   :: nptcls
        class(image), optional,    intent(inout) :: pickrefs(:)
        character(len=LONGSTRLEN) :: boxfile
        type(pickgau) :: gaup, gaup_refine
        integer, allocatable :: pos(:,:)
        integer :: box
        real    :: maxdiam
        boxfile = basename(fname_new_ext(trim(micname),'box'))
        call read_mic_raw(micname, smpd)
        if( present(pickrefs) )then
            call gaup%new_gaupicker(       params_glob%pcontrast, SMPD_SHRINK1, params_glob%moldiam, offset=OFFSET)
            call gaup_refine%new_gaupicker(params_glob%pcontrast, SMPD_SHRINK2, params_glob%moldiam, offset=1)
        else
            call gaup%new_refpicker(       params_glob%pcontrast, SMPD_SHRINK1, params_glob%mskdiam, pickrefs, offset=OFFSET)
            call gaup_refine%new_refpicker(params_glob%pcontrast, SMPD_SHRINK2, params_glob%mskdiam, pickrefs, offset=1)
        endif
        call gaup%gaupick(gaup_refine)
        nptcls = gaup_refine%get_nboxes()
        if( nptcls > 0 )then
            ! write coordinates
            call gaup_refine%get_positions(pos, smpd_new=smpd)
            maxdiam = params_glob%moldiam + params_glob%moldiam * BOX_EXP_FAC
            box     = find_larger_magic_box(round2even(maxdiam / smpd))
            call write_boxfile(nptcls, pos, box, boxfile)
            call make_relativepath(CWD_GLOB, boxfile, boxfile_out) ! returns absolute path
        else
            ! no particles found
            boxfile_out = ''
        endif
        call gaup%kill
        call gaup_refine%kill
    end subroutine exec_gaupick

    subroutine write_boxfile( n, coordinates, box, fname )
        integer,          intent(in) :: n, coordinates(n,2), box
        character(len=*), intent(in) :: fname
        integer :: funit, ibox, iostat
        call fopen(funit, status='REPLACE', action='WRITE', file=trim(adjustl(fname)), iostat=iostat)
        call fileiochk('simple_picker_utils; write_boxfile ', iostat)
        do ibox = 1,n
            write(funit,'(I7,I7,I7,I7,I7)') coordinates(ibox,1), coordinates(ibox,2), box, box, -3
        end do
        call fclose(funit)
    end subroutine write_boxfile

end module simple_picker_utils
