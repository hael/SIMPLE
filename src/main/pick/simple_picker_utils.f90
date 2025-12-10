module simple_picker_utils 
!$ use omp_lib
!$ use omp_lib_kinds
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_image,      only: image
use simple_pickseg,    only: pickseg
use simple_pickref
use simple_picksegdiam
implicit none

public :: exec_refpick, exec_segpick, exec_segdiampick, exec_gaupick
private
#include "simple_local_flags.inc"

real,    parameter :: SMPD_SHRINK1 = 4.0, SMPD_SHRINK2 = 2.0
real,    parameter :: GAUSIG       = 2.5
integer, parameter :: OFFSET       = 3
logical, parameter :: L_WRITE      = .false.
logical, parameter :: L_DEBUG      = .false.

contains

    subroutine exec_refpick( micname, boxfile_out, thumb_den_out, smpd, nptcls, pickrefs, dir_out, nboxes_max )
        use simple_string_utils, only: str2real, parsestr
        class(string),           intent(in)    :: micname
        class(string),           intent(out)   :: boxfile_out, thumb_den_out
        real,                    intent(in)    :: smpd    !< sampling distance in A
        integer,                 intent(out)   :: nptcls
        class(image),  optional, intent(inout) :: pickrefs(:)
        class(string), optional, intent(in)    :: dir_out
        integer,       optional, intent(in)    :: nboxes_max
        type(pickref)                  :: refp, refp_refine
        type(string) :: boxfile, fbody_here, ext, fname_thumb_den
        real    :: maxdiam
        integer :: box
        logical :: l_roi, l_backgr_subtr
        fbody_here      = basename(micname)
        ext             = fname2ext(fbody_here)
        fbody_here      = get_fbody(fbody_here, ext)
        fname_thumb_den = fbody_here//DEN_SUFFIX//JPG_EXT
        boxfile         = basename(fname_new_ext(micname,'box'))
        if( present(dir_out) ) boxfile         = dir_out//'/'//boxfile%to_char()
        if( present(dir_out) ) fname_thumb_den = dir_out//'/'//fname_thumb_den%to_char()
        l_roi           = trim(params_glob%pick_roi).eq.'yes'
        l_backgr_subtr  = l_roi .or. (trim(params_glob%backgr_subtr).eq.'yes')
        call read_mic_raw_pickref(micname, smpd, subtr_backgr=l_backgr_subtr)
        if( present(nboxes_max) )then
            call refp%new(params_glob%pcontrast, SMPD_SHRINK1, pickrefs, offset=OFFSET, roi=l_roi, nboxes_max=params_glob%nboxes_max)
        else
            call refp%new(params_glob%pcontrast, SMPD_SHRINK1, pickrefs, offset=OFFSET, roi=l_roi)
        endif
        call refp_refine%new(params_glob%pcontrast, SMPD_SHRINK2, pickrefs, offset=1)
        call refp%refpick(refp_refine)
        ! write
        maxdiam = refp%get_maxdiam() + refp%get_maxdiam() * BOX_EXP_FAC
        box     = find_larger_magic_box(round2even(maxdiam / smpd))
        call refp_refine%report_boxfile(box, smpd, boxfile, nptcls)
        call refp%report_thumb_den(fname_thumb_den)
        if( nptcls == 0 )then
            boxfile_out   = ''
            thumb_den_out = ''
        else
            boxfile_out   = simple_abspath(boxfile)
            thumb_den_out = simple_abspath(fname_thumb_den)
        endif
        call refp%kill
        call refp_refine%kill
    end subroutine exec_refpick

    subroutine exec_segpick( micname, boxfile_out, nptcls, dir_out, moldiam, winsz )
        class(string),           intent(in)  :: micname
        class(string),           intent(out) :: boxfile_out
        integer,                 intent(out) :: nptcls
        class(string), optional, intent(in)  :: dir_out
        real,          optional, intent(in)  :: moldiam
        integer,       optional, intent(in)  :: winsz
        type(string)  :: boxfile
        type(pickseg) :: picker
        boxfile = basename(fname_new_ext(micname,'box'))
        if( present(dir_out) ) boxfile = dir_out//'/'//boxfile%to_char()
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

    subroutine exec_segdiampick( micname, boxfile_out, smpd, nptcls, moldiam_max, dir_out )
        class(string),           intent(in)  :: micname
        class(string),           intent(out) :: boxfile_out
        real,                    intent(in)  :: smpd    !< sampling distance in A
        integer,                 intent(out) :: nptcls
        real,                    intent(in)  :: moldiam_max
        class(string), optional, intent(in)  :: dir_out
        type(string)      :: boxfile
        type(picksegdiam) :: picker
        boxfile = basename(fname_new_ext(micname,'box'))
        if( present(dir_out) ) boxfile = dir_out//'/'//boxfile%to_char()
        call picker%pick(micname, SMPD, moldiam_max, params_glob%pcontrast)
        call picker%write_pos_and_diams(boxfile, nptcls)
        if( nptcls == 0 )then
            boxfile_out = ''
        else
            boxfile_out = simple_abspath(boxfile)
        endif
    end subroutine exec_segdiampick

    subroutine exec_gaupick( micname, boxfile_out, smpd, nptcls, dir_out, nboxes_max )
        use simple_strategy2D_utils, only: alloc_imgarr, dealloc_imgarr
        class(string),              intent(in)  :: micname
        class(string),              intent(out) :: boxfile_out
        real,                       intent(in)  :: smpd    !< sampling distance in A
        integer,                    intent(out) :: nptcls
        class(string),    optional, intent(in)  :: dir_out
        integer,          optional, intent(in)  :: nboxes_max
        type(image),      allocatable :: pickrefs(:)
        type(pickref),    allocatable :: comprefp_refine(:)
        real,             allocatable :: moldiams(:)
        character(len=4), allocatable :: moldiams_str(:)
        type(pickref) :: refp, refp_refine
        type(string)  :: kind, boxfile
        real          :: maxdiam, sig, r, ang
        integer       :: box, i, istr, num_entries, nmoldiams, sel
        logical       :: l_roi, l_backgr_subtr, l_competitive
        boxfile = basename(fname_new_ext(micname,'box'))
        if( present(dir_out) ) boxfile = dir_out//'/'//boxfile%to_char()
        l_roi          = trim(params_glob%pick_roi).eq.'yes'
        l_backgr_subtr = l_roi .or. (trim(params_glob%backgr_subtr).eq.'yes')
        l_competitive  = .false.
        call read_mic_raw_pickref(micname, smpd, subtr_backgr=l_backgr_subtr)
        ! execution fork
        nmoldiams = 1
        if( params_glob%nmoldiams > 1 )then
            ! array of diameters
            nmoldiams = params_glob%nmoldiams
            box       = round2even((1.0+BOX_EXP_FAC)*params_glob%moldiam_max/ smpd) + 2
            moldiams  = equispaced_vals(params_glob%moldiam, params_glob%moldiam_max, nmoldiams)
            call alloc_imgarr(1, [box, box,1], smpd, pickrefs)
            l_competitive  = .true.
        else if( (params_glob%nmoldiams==1) .and. (params_glob%interactive.eq.'yes') )then
            ! still relevant?
            box = round2even((1.0+BOX_EXP_FAC)*params_glob%moldiam/ smpd) + 2
            sig = params_glob%moldiam / (smpd * GAUSIG)
            call alloc_imgarr(1, [box, box,1], smpd, pickrefs)
            call pickrefs(1)%gauimg2D(sig, sig)
        else if( .not.(params_glob%multi_moldiams.eq.'') )then
            ! multiple moldiam pick that uses multiple gaussians, generates .box file outputs
            istr        = 1
            num_entries = 0
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
            call parsestr(params_glob%multi_moldiams,',',moldiams_str,nmoldiams)
            allocate(moldiams(nmoldiams))
            do i = 1, nmoldiams
                moldiams(i) = str2real(moldiams_str(i))
            end do
            deallocate(moldiams_str)
            box = round2even((1.0+BOX_EXP_FAC)*maxval(moldiams)) + 2
            call alloc_imgarr(1, [box, box,1], smpd, pickrefs)
            l_competitive = .true.
        else
            ! single gaussian and optional shape
            kind = trim(params_glob%pickkind)
            if(params_glob%interactive == "yes") then
                if( params_glob%ring == "yes" ) kind = 'ring'
            endif
            box = round2even((1.0+BOX_EXP_FAC)*params_glob%moldiam/ smpd) + 2
            sig = params_glob%moldiam / (smpd * GAUSIG)
            select case(kind%to_char())
                case('gau')
                    ! single gaussian
                    call alloc_imgarr(1, [box, box,1], smpd, pickrefs)
                    call pickrefs(1)%gauimg2D(sig, sig)
                case('disc')
                    ! small specimen in membrane patch
                    ! 1:gaussian, 2:top view, 3-8:sideviews
                    call alloc_imgarr(8, [box, box,1], smpd, pickrefs)
                    call pickrefs(1)%gauimg2D(sig, sig)
                    r = 0.5*params_glob%moldiam/smpd - 1.0
                    pickrefs(2) = 1.
                    call pickrefs(2)%mask(r, 'soft', backgr=0.)
                    ang = 0.
                    do i = 3,8
                        call pickrefs(i)%disc_sideview([box, box,1], smpd, r)
                        if( ang > 0.5 ) call pickrefs(i)%rtsq(ang,0.,0.)
                        call pickrefs(i)%mask(real(box/2-1),'hard')
                        ang = ang + 30.0
                    enddo
                case('ring')
                    ! TM complex/apoferritin-like shape
                    ! 1:gaussian, 2:near hollow ring
                    call alloc_imgarr(2, [box, box,1], smpd, pickrefs)
                    call pickrefs(1)%gauimg2D(sig, sig)
                    call pickrefs(2)%soft_ring([box,box,1], smpd, 0.5*params_glob%moldiam/smpd)
                case DEFAULT
                    THROW_HARD('Unsupported picking mode: '//kind%to_char())
            end select
        endif
        ! Effective picking
        if( l_competitive )then
            ! competitive picking
            allocate(comprefp_refine(nmoldiams))
            do i = 1,nmoldiams
                ! reference
                sig = moldiams(i) / (smpd * GAUSIG)
                call pickrefs(1)%gauimg2D(sig, sig)
                ! pick
                if( present(nboxes_max) )then
                    call refp%new(params_glob%pcontrast, SMPD_SHRINK1, pickrefs, offset=OFFSET, roi=l_roi, nboxes_max=params_glob%nboxes_max)
                else
                    call refp%new(params_glob%pcontrast, SMPD_SHRINK1, pickrefs, offset=OFFSET, roi=l_roi)
                endif
                call comprefp_refine(i)%new(params_glob%pcontrast, SMPD_SHRINK2, pickrefs, offset=1)
                call refp%refpick(comprefp_refine(i))
            enddo
            ! merge
            call multiref_merge(nmoldiams, comprefp_refine(:), sel)
            ! output
            maxdiam = (1.0 + BOX_EXP_FAC ) * moldiams(sel)
            box     = find_larger_magic_box(round2even(maxdiam / smpd))
            call comprefp_refine(1)%report_boxfile(box, smpd, boxfile, nptcls)
            do i = 1,nmoldiams
                call comprefp_refine(i)%kill
            enddo
            deallocate(comprefp_refine)
        else
            ! Standard picking
            if( present(nboxes_max) )then
                call refp%new(params_glob%pcontrast, SMPD_SHRINK1, pickrefs, offset=OFFSET, roi=l_roi, nboxes_max=params_glob%nboxes_max)
            else
                call refp%new(params_glob%pcontrast, SMPD_SHRINK1, pickrefs, offset=OFFSET, roi=l_roi)
            endif
            call refp_refine%new(params_glob%pcontrast, SMPD_SHRINK2, pickrefs, offset=1)
            call refp%refpick(refp_refine)
            ! Output
            maxdiam = (1.0 + BOX_EXP_FAC ) * refp%get_maxdiam()
            box     = find_larger_magic_box(round2even(maxdiam / smpd))
            call refp_refine%report_boxfile(box, smpd, boxfile, nptcls)
            call refp_refine%kill
        endif
        ! boxfile
        if( nptcls == 0 )then
            boxfile_out = ''
        else
            boxfile_out = simple_abspath(boxfile)
        endif
        ! cleanup
        call refp%kill
        call dealloc_imgarr(pickrefs)
    end subroutine exec_gaupick

end module simple_picker_utils
