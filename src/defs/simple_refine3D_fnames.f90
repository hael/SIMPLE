!@descr: filename helpers for refine3D handoff and output artifacts
module simple_refine3D_fnames
use simple_defs_fname, only: BIN_EXT, MRC_EXT, TXT_EXT, &
    &CAVGS_ITER_FBODY, FSC_FBODY, POLAR_REFS_FBODY, STARTVOL_FBODY, VOL_FBODY
use simple_string,       only: string
use simple_string_utils, only: int2str_pad
implicit none

private
public :: refine3D_state_vol_fbody
public :: refine3D_state_vol_fname
public :: refine3D_state_vol_suffix_fname
public :: refine3D_state_halfvol_fname
public :: refine3D_startvol_fbody
public :: refine3D_startvol_fname
public :: refine3D_startvol_half_fname
public :: refine3D_fsc_fbody
public :: refine3D_fsc_fname
public :: refine3D_fsc_plot_fbody
public :: refine3D_resolution_txt_fbody
public :: refine3D_iter_refs_fname
public :: refine3D_iter_vol_fname
public :: refine3D_partial_rec_fbody
public :: refine3D_partial_rec_glob
public :: refine3D_partial_rec_fname
public :: refine3D_partial_rho_fname
public :: refine3D_polar_cavgs_part_fname
public :: refine3D_polar_ctfsqsums_part_fname
public :: refine3D_polar_sums_fname
public :: refine3D_polar_ctf2_fname
public :: refine3D_obsfield_part_fname
public :: refine3D_obsfield_fname
public :: refine3D_polar_refs_fbody
public :: refine3D_polar_refs_fname
public :: refine3D_polar_ref_part_fname
public :: refine3D_bench_fname
public :: refine3D_distr_bench_fname
public :: refine3D_volassemble_bench_fname

contains

    type(string) function state_tag( state ) result(tag)
        integer, intent(in) :: state
        tag = int2str_pad(state, 2)
    end function state_tag

    type(string) function iter_tag( iter ) result(tag)
        integer, intent(in) :: iter
        tag = int2str_pad(iter, 3)
    end function iter_tag

    type(string) function part_tag( part, numlen ) result(tag)
        integer, intent(in) :: part, numlen
        tag = int2str_pad(part, max(1, numlen))
    end function part_tag

    type(string) function half_suffix( half ) result(suffix)
        character(len=*), intent(in) :: half
        select case(trim(half))
            case('even', '_even')
                suffix = '_even'
            case('odd', '_odd')
                suffix = '_odd'
            case default
                suffix = '_'//trim(half)
        end select
    end function half_suffix

    type(string) function refine3D_state_vol_fbody( state ) result(fname)
        integer, intent(in) :: state
        fname = string(VOL_FBODY)//state_tag(state)
    end function refine3D_state_vol_fbody

    type(string) function refine3D_state_vol_fname( state ) result(fname)
        integer, intent(in) :: state
        fname = refine3D_state_vol_fbody(state)
        fname = fname//MRC_EXT
    end function refine3D_state_vol_fname

    type(string) function refine3D_state_vol_suffix_fname( state, suffix ) result(fname)
        integer,          intent(in) :: state
        character(len=*), intent(in) :: suffix
        fname = refine3D_state_vol_fbody(state)//trim(suffix)//MRC_EXT
    end function refine3D_state_vol_suffix_fname

    type(string) function refine3D_state_halfvol_fname( state, half, unfil ) result(fname)
        integer,          intent(in) :: state
        character(len=*), intent(in) :: half
        logical,          intent(in), optional :: unfil
        fname = refine3D_state_vol_fbody(state)//half_suffix(half)
        if( present(unfil) )then
            if( unfil ) fname = fname//'_unfil'
        endif
        fname = fname//MRC_EXT
    end function refine3D_state_halfvol_fname

    type(string) function refine3D_startvol_fbody( state ) result(fname)
        integer, intent(in) :: state
        fname = string(STARTVOL_FBODY)//state_tag(state)
    end function refine3D_startvol_fbody

    type(string) function refine3D_startvol_fname( state ) result(fname)
        integer, intent(in) :: state
        fname = refine3D_startvol_fbody(state)
        fname = fname//MRC_EXT
    end function refine3D_startvol_fname

    type(string) function refine3D_startvol_half_fname( state, half, unfil ) result(fname)
        integer,          intent(in) :: state
        character(len=*), intent(in) :: half
        logical,          intent(in), optional :: unfil
        fname = string(STARTVOL_FBODY)//state_tag(state)//half_suffix(half)
        if( present(unfil) )then
            if( unfil ) fname = fname//'_unfil'
        endif
        fname = fname//MRC_EXT
    end function refine3D_startvol_half_fname

    type(string) function refine3D_fsc_fbody( state ) result(fname)
        integer, intent(in) :: state
        fname = string(FSC_FBODY)//state_tag(state)
    end function refine3D_fsc_fbody

    type(string) function refine3D_fsc_fname( state ) result(fname)
        integer, intent(in) :: state
        fname = refine3D_fsc_fbody(state)//BIN_EXT
    end function refine3D_fsc_fname

    type(string) function refine3D_fsc_plot_fbody( state, iter ) result(fname)
        integer, intent(in) :: state, iter
        fname = refine3D_fsc_fbody(state)//'_iter'//iter_tag(iter)
    end function refine3D_fsc_plot_fbody

    type(string) function refine3D_resolution_txt_fbody( state, iter ) result(fname)
        integer, intent(in) :: state
        integer, intent(in), optional :: iter
        fname = string('RESOLUTION_STATE')//state_tag(state)
        if( present(iter) ) fname = fname//'_ITER'//iter_tag(iter)
    end function refine3D_resolution_txt_fbody

    type(string) function refine3D_iter_refs_fname( iter ) result(fname)
        integer, intent(in) :: iter
        fname = string(CAVGS_ITER_FBODY)//iter_tag(iter)
        fname = fname//MRC_EXT
    end function refine3D_iter_refs_fname

    type(string) function refine3D_iter_vol_fname( state, iter, suffix ) result(fname)
        integer,          intent(in) :: state, iter
        character(len=*), intent(in), optional :: suffix
        fname = refine3D_state_vol_fbody(state)//'_iter'//iter_tag(iter)
        if( present(suffix) ) fname = fname//trim(suffix)
        fname = fname//MRC_EXT
    end function refine3D_iter_vol_fname

    type(string) function refine3D_partial_rec_fbody( state, part, numlen ) result(fname)
        integer, intent(in) :: state, part, numlen
        fname = refine3D_state_vol_fbody(state)//'_part'//part_tag(part, numlen)
    end function refine3D_partial_rec_fbody

    function refine3D_partial_rec_glob( path ) result(pattern)
        character(len=*), intent(in), optional :: path
        character(len=:), allocatable :: pattern
        if( present(path) )then
            pattern = trim(path)//'*'//VOL_FBODY//'*part*'
        else
            pattern = '*'//VOL_FBODY//'*part*'
        endif
    end function refine3D_partial_rec_glob

    type(string) function refine3D_partial_rec_fname( state, part, numlen, half ) result(fname)
        integer,          intent(in) :: state, part, numlen
        character(len=*), intent(in) :: half
        fname = refine3D_partial_rec_fbody(state, part, numlen)//half_suffix(half)
        fname = fname//MRC_EXT
    end function refine3D_partial_rec_fname

    type(string) function refine3D_partial_rho_fname( state, part, numlen, half ) result(fname)
        integer,          intent(in) :: state, part, numlen
        character(len=*), intent(in) :: half
        fname = string('rho_')//refine3D_partial_rec_fname(state, part, numlen, half)
    end function refine3D_partial_rho_fname

    type(string) function refine3D_polar_cavgs_part_fname( half, part, numlen ) result(fname)
        character(len=*), intent(in) :: half
        integer,          intent(in) :: part, numlen
        fname = string('cavgs')//half_suffix(half)//'_part'//part_tag(part, numlen)//BIN_EXT
    end function refine3D_polar_cavgs_part_fname

    type(string) function refine3D_polar_ctfsqsums_part_fname( half, part, numlen ) result(fname)
        character(len=*), intent(in) :: half
        integer,          intent(in) :: part, numlen
        fname = string('ctfsqsums')//half_suffix(half)//'_part'//part_tag(part, numlen)//BIN_EXT
    end function refine3D_polar_ctfsqsums_part_fname

    type(string) function refine3D_polar_sums_fname( half ) result(fname)
        character(len=*), intent(in) :: half
        fname = string('polar_sums')//half_suffix(half)//BIN_EXT
    end function refine3D_polar_sums_fname

    type(string) function refine3D_polar_ctf2_fname( half ) result(fname)
        character(len=*), intent(in) :: half
        fname = string('polar_ctf2')//half_suffix(half)//BIN_EXT
    end function refine3D_polar_ctf2_fname

    type(string) function refine3D_obsfield_part_fname( state, part, numlen ) result(fname)
        integer, intent(in) :: state, part, numlen
        fname = string('obsfield_state')//state_tag(state)//'_part'//part_tag(part, numlen)//BIN_EXT
    end function refine3D_obsfield_part_fname

    type(string) function refine3D_obsfield_fname( state ) result(fname)
        integer, intent(in) :: state
        fname = string('obsfield_state')//state_tag(state)//BIN_EXT
    end function refine3D_obsfield_fname

    type(string) function refine3D_polar_refs_fbody() result(fname)
        fname = string(POLAR_REFS_FBODY)
    end function refine3D_polar_refs_fbody

    type(string) function refine3D_polar_refs_fname( half ) result(fname)
        character(len=*), intent(in), optional :: half
        fname = refine3D_polar_refs_fbody()
        if( present(half) ) fname = fname//half_suffix(half)
        fname = fname//BIN_EXT
    end function refine3D_polar_refs_fname

    type(string) function refine3D_polar_ref_part_fname( state, part, numlen, half ) result(fname)
        integer,          intent(in) :: state, part, numlen
        character(len=*), intent(in) :: half
        fname = refine3D_polar_refs_fbody()//'_s'//state_tag(state)//'_part'//part_tag(part, numlen)//&
            &half_suffix(half)//BIN_EXT
    end function refine3D_polar_ref_part_fname

    type(string) function refine3D_bench_fname( iter ) result(fname)
        integer, intent(in) :: iter
        fname = string('REFINE3D_BENCH_ITER')//iter_tag(iter)//TXT_EXT
    end function refine3D_bench_fname

    type(string) function refine3D_distr_bench_fname( iter ) result(fname)
        integer, intent(in) :: iter
        fname = string('DISTR_REFINE3D_BENCH_ITER')//iter_tag(iter)//TXT_EXT
    end function refine3D_distr_bench_fname

    type(string) function refine3D_volassemble_bench_fname( iter ) result(fname)
        integer, intent(in) :: iter
        fname = string('VOLASSEMBLE_BENCH_ITER')//iter_tag(iter)//TXT_EXT
    end function refine3D_volassemble_bench_fname

end module simple_refine3D_fnames
