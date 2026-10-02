!@descr: filename helpers for refine3D output artifacts
module simple_refine3D_fnames
use simple_defs_fname, only: BIN_EXT, MRC_EXT, TXT_EXT, JPG_EXT,&
    &CAVGS_ITER_FBODY, FSC_FBODY, STARTVOL_FBODY, VOL_FBODY
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
public :: refine3D_partial_rec_fname
public :: refine3D_partial_rho_fname
public :: refine3D_pcg_raw_accum_fname
public :: refine3D_pcg_trail_accum_fname
public :: refine3D_trail_rec_fbody
public :: refine3D_trail_rec_fname
public :: refine3D_trail_rho_fname
public :: refine3D_trail_manifest_fname
public :: refine3D_frozen_context_fname
public :: refine3D_frozen_rec_fbody
public :: refine3D_frozen_manifest_fname
public :: refine3D_frozen_pcg_fname
public :: refine3D_reproj_model_fname
public :: refine3D_cart_refvols_fname
public :: refine3D_bench_fname
public :: refine3D_strategy_bench_fname
public :: refine3D_volassemble_bench_fname
public :: refine3D_oris_heatmap_fname
public :: refine3D_cfar_summary_fname

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

    ! One atomically published worker artifact containing raw full-range B and
    ! real D for exactly one (state,half,part). The master is the only consumer
    ! allowed to fold or finalize these sufficient statistics.
    type(string) function refine3D_pcg_raw_accum_fname( state, part, numlen, half ) result(fname)
        integer,          intent(in) :: state, part, numlen
        character(len=*), intent(in) :: half
        fname = string('pcg_raw_state')//state_tag(state)//half_suffix(half)// &
            &'_part'//part_tag(part, numlen)//BIN_EXT
    end function refine3D_pcg_raw_accum_fname

    ! Persistent, already reduced PCG continuation statistic for one
    ! (state,half). It uses the raw-accumulator file format and is written
    ! atomically only by the shared/master owner after u/f blending.
    type(string) function refine3D_pcg_trail_accum_fname( state, half ) result(fname)
        integer,          intent(in) :: state
        character(len=*), intent(in) :: half
        fname = string('pcg_trail_state')//state_tag(state)//half_suffix(half)//BIN_EXT
    end function refine3D_pcg_trail_accum_fname

    ! Persistent trailing-reconstruction accumulator chain (blended, unregularized
    ! e/o Fourier sums + sampling densities). Deliberately does not contain the
    ! VOL_FBODY stem so partial-reconstruction globs and cleanup never match it.
    type(string) function refine3D_trail_rec_fbody( state ) result(fname)
        integer, intent(in) :: state
        fname = string('trailrec_state')//state_tag(state)
    end function refine3D_trail_rec_fbody

    type(string) function refine3D_trail_rec_fname( state, half ) result(fname)
        integer,          intent(in) :: state
        character(len=*), intent(in) :: half
        fname = refine3D_trail_rec_fbody(state)//half_suffix(half)
        fname = fname//MRC_EXT
    end function refine3D_trail_rec_fname

    type(string) function refine3D_trail_rho_fname( state, half ) result(fname)
        integer,          intent(in) :: state
        character(len=*), intent(in) :: half
        fname = string('rho_')//refine3D_trail_rec_fname(state, half)
    end function refine3D_trail_rho_fname

    !> Manifest tying the four chain files together as one validated artifact
    !! set: provenance (box/sampling/population/state layout), a generation
    !! counter, and the byte size of each component. Written last; readers must
    !! treat a missing or mismatching manifest as an invalid chain.
    type(string) function refine3D_trail_manifest_fname( state ) result(fname)
        integer, intent(in) :: state
        fname = refine3D_trail_rec_fbody(state)//TXT_EXT
    end function refine3D_trail_manifest_fname

    ! abinitio3D_addon frozen contributions. The run-level context ties every
    ! frozen set of one add-on run to its run identifier; the per-state sets
    ! carry their consuming box in the name, so a consumer selects its set by
    ! its own box_crop. None of these contain the VOL_FBODY or trailrec stems,
    ! so partial-reconstruction globs, chain validation and cleanup never
    ! match them.
    type(string) function refine3D_frozen_context_fname() result(fname)
        fname = string('abinitio3D_addon_frozen_context')//TXT_EXT
    end function refine3D_frozen_context_fname

    !> Gridding frozen set: <fbody>_{even,odd}.mrc and rho_<fbody>_{even,odd}.mrc
    type(string) function refine3D_frozen_rec_fbody( state, box ) result(fname)
        integer, intent(in) :: state, box
        fname = string('frozen_state')//state_tag(state)//'_box'//int2str_pad(box, 4)
    end function refine3D_frozen_rec_fbody

    !> Manifest of one gridding frozen set, written last
    type(string) function refine3D_frozen_manifest_fname( state, box ) result(fname)
        integer, intent(in) :: state, box
        fname = refine3D_frozen_rec_fbody(state, box)//TXT_EXT
    end function refine3D_frozen_manifest_fname

    !> PCG frozen raw (B,D) pair for one (state,half) at one box, in the raw
    !! accumulator format with its own provenance tag
    type(string) function refine3D_frozen_pcg_fname( state, box, half ) result(fname)
        integer,          intent(in) :: state, box
        character(len=*), intent(in) :: half
        fname = string('frozen_pcg_state')//state_tag(state)//'_box'//int2str_pad(box, 4)// &
            &half_suffix(half)//BIN_EXT
    end function refine3D_frozen_pcg_fname

    type(string) function refine3D_reproj_model_fname( half ) result(fname)
        character(len=*), intent(in) :: half
        fname = string('reprojection_model')//half_suffix(half)//BIN_EXT
    end function refine3D_reproj_model_fname

    !> The prepared Cartesian reference volumes of every state for one half-set
    !! (cart_refvols_even.bin, cart_refvols_odd.bin; plan section 6.6).
    type(string) function refine3D_cart_refvols_fname( half ) result(fname)
        character(len=*), intent(in) :: half
        fname = string('cart_refvols')//half_suffix(half)//BIN_EXT
    end function refine3D_cart_refvols_fname

    !> per-iteration bench record; with part present the collision-free per-partition
    !! record (REFINE3D_BENCH_ITERnnn_PARTppp.txt), the plain name stays partition 1's legacy file
    type(string) function refine3D_bench_fname( iter, part, numlen ) result(fname)
        integer,           intent(in) :: iter
        integer, optional, intent(in) :: part, numlen
        integer :: nl
        fname = string('REFINE3D_BENCH_ITER')//iter_tag(iter)
        if( present(part) )then
            nl = 1
            if( present(numlen) ) nl = numlen
            fname = fname//'_PART'//part_tag(part, nl)
        endif
        fname = fname//TXT_EXT
    end function refine3D_bench_fname

    type(string) function refine3D_strategy_bench_fname( iter ) result(fname)
        integer, intent(in) :: iter
        fname = string('REFINE3D_STRATEGY_BENCH_ITER')//iter_tag(iter)//TXT_EXT
    end function refine3D_strategy_bench_fname

    type(string) function refine3D_volassemble_bench_fname( iter ) result(fname)
        integer, intent(in) :: iter
        fname = string('VOLASSEMBLE_BENCH_ITER')//iter_tag(iter)//TXT_EXT
    end function refine3D_volassemble_bench_fname

    type(string) function refine3D_oris_heatmap_fname( state ) result(fname)
        integer, intent(in) :: state
        fname = string('orientations_distribution_state')//state_tag(state)//JPG_EXT
    end function refine3D_oris_heatmap_fname

    type(string) function refine3D_cfar_summary_fname( iter ) result(fname)
        integer, intent(in) :: iter
        fname = string('CFAR_SUMMARY_ITER')//iter_tag(iter)//TXT_EXT
    end function refine3D_cfar_summary_fname

end module simple_refine3D_fnames
