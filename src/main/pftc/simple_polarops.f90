module simple_polarops
use simple_core_module_api
use simple_builder,          only: builder, build_glob
use simple_image,            only: image
use simple_parameters,       only: params_glob
use simple_polarft_calc,     only: polarft_calc, pftc_glob
use simple_sp_project,       only: sp_project
use simple_strategy2D_utils, only: calc_cavg_offset
use simple_imgarr_utils,     only: alloc_imgarr, write_imgarr, dealloc_imgarr
use simple_cmdline,          only: cmdline
implicit none

! State
public :: polar_cavger_new
public :: polar_cavger_zero_pft_refs
public :: polar_cavger_set_ref_pftc
public :: polar_cavger_calc_pops
public :: polar_cavger_update_sums
public :: polar_cavger_kill
public :: center_3Dpolar_refs
! Restore
public :: polar_cavger_merge_eos_and_norm2D
public :: polar_cavger_merge_eos_and_norm
public :: polar_cavger_calc_and_write_frcs_and_eoavg
public :: polar_prep2Dref
public :: polar_filterrefs
public :: polar_cavger_gen2Dclassdoc
! I/O
public :: polar_cavger_refs2cartesian
public :: polar_cavger_read
public :: polar_cavger_write
public :: polar_cavger_writeall
public :: polar_cavger_writeall_pftcrefs
public :: polar_cavger_write_cartrefs
public :: polar_cavger_read_all
public :: polar_cavger_readwrite_partial_sums
public :: polar_cavger_assemble_sums_from_parts
public :: polar_cavger_dims_from_header 
private
#include "simple_local_flags.inc"

complex(dp), allocatable :: pfts_even(:,:,:), pfts_odd(:,:,:), pfts_merg(:,:,:)
real(dp),    allocatable :: ctf2_even(:,:,:), ctf2_odd(:,:,:)
integer,     allocatable :: prev_eo_pops(:,:), eo_pops(:,:)
real                     :: smpd       = 0.
integer                  :: ncls       = 0
integer                  :: kfromto(2) = 0
integer                  :: pftsz      = 0
integer                  :: nrots      = 0
logical                  :: l_comlin   = .false.

!-----------------------------
! Interfaces implemented in submodules
!-----------------------------
interface

    ! ===== STATE: simple_polarops_state.f90

    module subroutine polar_cavger_new(pftc, comlin, nrefs)
        class(polarft_calc),     intent(in) :: pftc
        logical,                 intent(in) :: comlin
        integer,       optional, intent(in) :: nrefs
    end subroutine

    module subroutine polar_cavger_zero_pft_refs()
    end subroutine

    module subroutine polar_cavger_set_ref_pftc(icls, which, pftc)
        integer,             intent(in)    :: icls
        character(len=*),    intent(in)    :: which
        class(polarft_calc), intent(inout) :: pftc
    end subroutine

    module subroutine polar_cavger_calc_pops(spproj)
        class(sp_project), target, intent(in) :: spproj
    end subroutine

    module subroutine polar_cavger_update_sums(nptcls, pinds, spproj, pftc, incr_shifts, is3D)
        integer,                         intent(in)    :: nptcls
        integer,                         intent(in)    :: pinds(nptcls)
        class(sp_project),               intent(inout) :: spproj
        class(polarft_calc), target,     intent(inout) :: pftc
        real,                  optional, intent(in)    :: incr_shifts(2,nptcls)
        logical,               optional, intent(in)    :: is3D
    end subroutine

    module subroutine polar_cavger_kill()
    end subroutine

    module subroutine center_3Dpolar_refs(pftc, algndoc, algnrefs)
        class(polarft_calc), intent(inout) :: pftc
        class(oris),         intent(inout) :: algndoc
        class(oris),         intent(in)    :: algnrefs
    end subroutine
    
    ! ===== RESTORE: simple_polarops_restore.f90

    module subroutine polar_cavger_merge_eos_and_norm2D()
    end subroutine

    module subroutine polar_cavger_merge_eos_and_norm(reforis, cl_weight)
        type(oris),           intent(in) :: reforis
        real,       optional, intent(in) :: cl_weight
    end subroutine

    module subroutine polar_cavger_calc_and_write_frcs_and_eoavg(fname, cline)
        class(string), intent(in) :: fname
        type(cmdline), intent(in) :: cline
    end subroutine

    module subroutine polar_prep2Dref(icls, gaufilt)
        integer, intent(in) :: icls
        logical, intent(in) :: gaufilt
    end subroutine

    module subroutine polar_filterrefs( icls, filter )
        integer, intent(in) :: icls
        real,    intent(in) :: filter(:)
    end subroutine

    module subroutine polar_cavger_gen2Dclassdoc(spproj)
        class(sp_project), target, intent(inout) :: spproj
    end subroutine

    ! ===== I/O: simple_polarops_io.f90

    module subroutine polar_cavger_refs2cartesian(pftc, cavgs, which, pfts_in)
        class(polarft_calc),     intent(in)    :: pftc
        type(image),             intent(inout) :: cavgs(ncls)
        character(len=*),        intent(in)    :: which
        complex(dp),   optional, intent(in)    :: pfts_in(1:pftsz,kfromto(1):kfromto(2),1:ncls)
    end subroutine

    module subroutine polar_cavger_read( fname, which )
        class(string),    intent(in) :: fname
        character(len=*), intent(in) :: which
    end subroutine

    module subroutine polar_cavger_write(fname, which)
        class(string),    intent(in) :: fname
        character(len=*), intent(in) :: which
    end subroutine

    module subroutine polar_cavger_writeall(tmpl_fname)
        class(string), intent(in) :: tmpl_fname
    end subroutine

    module subroutine polar_cavger_writeall_pftcrefs(tmpl_fname)
        class(string), intent(in) :: tmpl_fname
    end subroutine

    module subroutine polar_cavger_write_cartrefs(pftc, tmpl_fname, which)
        class(polarft_calc), intent(in) :: pftc
        class(string),       intent(in) :: tmpl_fname
        character(len=*),    intent(in) :: which
    end subroutine

    module subroutine polar_cavger_read_all(fname)
        class(string), intent(in) :: fname
    end subroutine

    module subroutine polar_cavger_readwrite_partial_sums(which)
        character(len=*), intent(in) :: which
    end subroutine

    module subroutine polar_cavger_assemble_sums_from_parts(reforis, clin_anneal)
        type(oris), optional, intent(in) :: reforis
        real,       optional, intent(in) :: clin_anneal
    end subroutine

    module subroutine polar_cavger_dims_from_header(fname, pftsz_here, kfromto_here, ncls_here)
        class(string), intent(in)    :: fname
        integer,       intent(inout) :: pftsz_here, kfromto_here(2), ncls_here
    end subroutine

    module subroutine open_pft_or_ctf2_array_for_write( fname, funit )
        class(string), intent(in)  :: fname
        integer,       intent(out) :: funit
    end subroutine open_pft_or_ctf2_array_for_write

    module subroutine write_pft_array_local( funit, array )
        integer,     intent(in) :: funit
        complex(dp), intent(in) :: array(pftsz,kfromto(1):kfromto(2),ncls)
    end subroutine write_pft_array_local

    module subroutine write_pft_array( array, fname )
        complex(dp),   intent(in) :: array(pftsz,kfromto(1):kfromto(2),ncls)
        class(string), intent(in) :: fname
    end subroutine write_pft_array

    module subroutine write_ctf2_array_local( funit, array )
        integer,  intent(in) :: funit
        real(dp), intent(in) :: array(pftsz,kfromto(1):kfromto(2),ncls)
    end subroutine write_ctf2_array_local

    module subroutine write_ctf2_array( array, fname )
        real(dp),      intent(in) :: array(pftsz,kfromto(1):kfromto(2),ncls)
        class(string), intent(in) :: fname
    end subroutine write_ctf2_array

    module subroutine open_pft_array_for_read( fname, array, funit, dims, buffer )
        class(string),            intent(in)    :: fname
        complex(dp), allocatable, intent(inout) :: array(:,:,:)
        integer,                  intent(out)   :: funit, dims(4)
        complex(sp), allocatable, intent(inout) :: buffer(:,:,:)
    end subroutine open_pft_array_for_read

    module subroutine transfer_pft_array_buffer( array, funit, dims, buffer )
        complex(dp), intent(inout) :: array(pftsz,kfromto(1):kfromto(2),ncls)
        integer,     intent(in)    :: funit, dims(4)
        complex(sp), intent(inout) :: buffer(dims(1),dims(2):dims(3),dims(4))
    end subroutine transfer_pft_array_buffer

    module subroutine read_pft_array( fname, array)
        class(string),            intent(in)    :: fname
        complex(dp), allocatable, intent(inout) :: array(:,:,:)
    end subroutine read_pft_array

    module subroutine open_ctf2_array_for_read( fname, array, funit, dims, buffer )
        class(string),         intent(in)    :: fname
        real(dp), allocatable, intent(inout) :: array(:,:,:)
        integer,               intent(out)   :: funit, dims(4)
        real(sp), allocatable, intent(inout) :: buffer(:,:,:)
    end subroutine open_ctf2_array_for_read

    module subroutine transfer_ctf2_array_buffer( array, funit, dims, buffer )
        real(dp), intent(inout) :: array(pftsz,kfromto(1):kfromto(2),ncls)
        integer,  intent(in)    :: funit, dims(4)
        real(sp), intent(inout) :: buffer(dims(1),dims(2):dims(3),dims(4))
    end subroutine transfer_ctf2_array_buffer

    module subroutine read_ctf2_array( fname, array )
        class(string),         intent(in)    :: fname
        real(dp), allocatable, intent(inout) :: array(:,:,:)
    end subroutine read_ctf2_array

    module subroutine pft2img()
    end subroutine pft2img

end interface

contains

    ! TEST UNIT

    subroutine test_polarops
        use simple_cmdline,    only: cmdline
        use simple_parameters, only: parameters
        integer,     parameter :: N=128
        integer,     parameter :: NIMGS=200
        integer,     parameter :: NCLS=5
        type(image)            :: tmpl_img, img, cavgs(NCLS)
        type(cmdline)          :: cline
        type(polarft_calc)     :: pftc
        type(parameters)       :: p
        type(builder)          :: b
        real    :: ang, shift(2), shifts(2,NIMGS)
        integer :: pinds(NIMGS), i, eo, icls
        ! dummy structure
        call tmpl_img%soft_ring([N,N,1], 1., 8.)
        call tmpl_img%fft
        call tmpl_img%shift2Dserial([ 8.,-16.])
        call img%soft_ring([N,N,1], 1., 12.)
        call img%fft
        call img%shift2Dserial([ 32., 0.])
        call tmpl_img%add(img)
        call img%soft_ring([N,N,1], 1., 16.)
        call img%fft
        call img%shift2Dserial([ -16., 8.])
        call tmpl_img%add(img)
        call img%soft_ring([N,N,1], 1., 32.)
        call img%fft
        call tmpl_img%add(img)
        call tmpl_img%ifft
        call tmpl_img%write(string('template.mrc'))
        ! init of options & parameters
        call cline%set('prg',    'xxx')
        call cline%set('objfun', 'cc')
        call cline%set('smpd',   1.0)
        call cline%set('box',    N)
        call cline%set('ctf',    'no')
        call cline%set('oritype','ptcl2D')
        call cline%set('ncls',    NCLS)
        call cline%set('nptcls',  NIMGs)
        call cline%set('lp',      3.)
        call cline%set('nthr',    8)
        call cline%set('mskdiam', real(N)/2-10.)
        ! Calculators
        call b%init_params_and_build_strategy2D_tbox(cline, p)
        call pftc%new(NCLS, [1,NIMGS], p%kfromto)
        pinds = (/(i,i=1,NIMGS)/)
        call b%img_crop_polarizer%init_polarizer(pftc, p%alpha)
        do i = 1,NIMGS
            shift = 10.*[ran3(), ran3()] - 5.
            ! ang   = 360. * ran3()
            ang   = 0.
            eo    = 0
            if( .not.is_even(i) ) eo = 1
            icls  = ceiling(ran3()*4.)
            call img%copy_fast(tmpl_img)
            call img%fft
            call img%shift2Dserial(-shift)
            call img%ifft
            call img%rtsq(ang, 0.,0.)
            call img%add_gauran(2.)
            call img%write(string('rotimgs.mrc'), i)
            call img%fft
            call b%spproj_field%set_euler(i, [0.,0.,ang])
            call b%spproj_field%set_shift(i, shift)
            call b%spproj_field%set(i,'w',1.0)
            call b%spproj_field%set(i,'state',1)
            call b%spproj_field%set(i,'class', icls)
            call b%spproj_field%set(i,'eo',eo)
            shifts(:,i) = -shift
            call b%img_crop_polarizer%polarize(pftc, img, i, isptcl=.true., iseven=eo==0, mask=b%l_resmsk)
        enddo
        call polar_cavger_new(pftc,.false.)
        call polar_cavger_update_sums(NIMGS, pinds, b%spproj, pftc, shifts)
        call polar_cavger_merge_eos_and_norm2D
        call polar_cavger_calc_and_write_frcs_and_eoavg(string(FRCS_FILE), cline)
        ! write
        call polar_cavger_write(string('cavgs_even.bin'), 'even')
        call polar_cavger_write(string('cavgs_odd.bin'),  'odd')
        call polar_cavger_write(string('cavgs.bin'),      'merged')
        call polar_cavger_refs2cartesian(pftc, cavgs, 'even')
        call write_imgarr(cavgs, string('cavgs_even.mrc'))
        call polar_cavger_refs2cartesian(pftc, cavgs, 'odd')
        call write_imgarr(cavgs, string('cavgs_odd.mrc'))
        call polar_cavger_refs2cartesian(pftc, cavgs, 'merged')
        call write_imgarr(cavgs, string('cavgs_merged.mrc'))
        call polar_cavger_kill
        ! read & write again
        call polar_cavger_new(pftc,.false.)
        call polar_cavger_read(string('cavgs_even.bin'), 'even')
        call polar_cavger_read(string('cavgs_odd.bin'),  'odd')
        call polar_cavger_read(string('cavgs.bin'),      'merged')
        call polar_cavger_refs2cartesian(pftc, cavgs, 'even')
        call write_imgarr(cavgs, string('cavgs2_even.mrc'))
        call polar_cavger_refs2cartesian(pftc, cavgs, 'odd')
        call write_imgarr(cavgs, string('cavgs2_odd.mrc'))
        call polar_cavger_refs2cartesian(pftc, cavgs, 'merged')
        call write_imgarr(cavgs, string('cavgs2_merged.mrc'))
        call polar_cavger_kill
    end subroutine test_polarops

end module simple_polarops
