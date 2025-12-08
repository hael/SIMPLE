module simple_sp_project
include 'simple_lib.f08'
use simple_cmdline,           only: cmdline
use simple_discrete_stack_io, only: dstack_io
use simple_map_reduce,        only: split_nobjs_even
use json_kinds
use json_module
use simple_gui_utils
use simple_histogram
use simple_image
use simple_rec_list
use simple_stack_io
use simple_starfile
implicit none

public :: sp_project, oritype2segment
private
#include "simple_local_flags.inc"

integer(kind(ENUM_ORISEG)), parameter :: MAXN_OS_SEG = 13


type sp_project
    ! ORIS REPRESENTATIONS OF BINARY FILE SEGMENTS
    type(oris) :: os_mic    ! segment 1
    type(oris) :: os_stk    ! segment 2
    type(oris) :: os_ptcl2D ! segment 3
    type(oris) :: os_cls2D  ! segment 4
    type(oris) :: os_cls3D  ! segment 5
    type(oris) :: os_ptcl3D ! segment 6
    type(oris) :: os_out    ! segment 7
    type(oris) :: os_optics ! segment 8

    ! PROJECT / JOB / ENV
    type(oris) :: projinfo  ! segment 11
    type(oris) :: jobproc   ! segment 12
    type(oris) :: compenv   ! segment 13

    ! binary file-handler
    type(binoris) :: bos
contains

    ! CORE
    procedure          :: new_seg_with_ptr
    procedure          :: kill
    procedure, private :: assign
    generic            :: assignment(=) => assign
    procedure          :: copy
    ! CORE - field updaters
    generic            :: update_projinfo => update_projinfo_1, update_projinfo_2
    procedure, private :: update_projinfo_1
    procedure, private :: update_projinfo_2
    procedure          :: update_compenv
    procedure          :: append_project
    procedure          :: append_job_descr2jobproc
    procedure          :: replace_project
    procedure          :: merge_algndocs
    procedure          :: projrecords2proj
    ! CORE - getters/setters
    procedure          :: ptr2oritype
    procedure          :: get_smpd
    procedure          :: get_ctfflag
    procedure          :: get_ctfflag_type
    procedure          :: get_n_insegment
    procedure          :: get_n_insegment_state
    procedure          :: get_ctfparams
    procedure          :: get_sp_oris
    procedure          :: set_sp_oris

    ! MIC
    procedure          :: add_single_movie
    procedure          :: add_movies
    procedure          :: add_intgs
    procedure          :: report_state2mic
    ! MIC - getters
    procedure          :: has_boxfile
    procedure          :: get_movies_table
    procedure          :: get_mics_table
    procedure          :: get_micparams
    procedure          :: get_micname
    procedure          :: get_mic_kind
    procedure          :: set_boxfile
    procedure          :: get_nmovies
    procedure          :: get_nintgs
    procedure          :: get_nframes
    procedure          :: get_mic2stk_inds

    ! STK
    procedure          :: add_stk
    procedure          :: add_single_stk
    procedure, private :: add_stktab_1
    procedure, private :: add_stktab_2
    generic            :: add_stktab => add_stktab_1, add_stktab_2
    procedure          :: split_stk
    procedure          :: write_substk
    procedure, private :: add_scale_tag
    procedure          :: report_state2stk
    ! STK - getters
    procedure          :: has_phaseplate
    procedure          :: get_box
    procedure          :: get_stkname
    procedure          :: get_stkname_and_ind
    procedure          :: get_nstks

    ! CLS2D/3D
    procedure          :: shape_ranked_cavgs2jpg
    procedure          :: cavgs2jpg
    procedure          :: cavgs2mrc
    ! CLS2D/3D - getters
    procedure          :: get_selected_clsinds
    procedure          :: get_cavgs_stk
    procedure          :: set_cavgs_thumb

    ! PTCL2D/3D
    procedure          :: map_ptcl_ind2stk_ind
    procedure          :: map_cavgs_selection
    procedure          :: map2Dshifts23D
    procedure          :: map2ptcls
    procedure          :: map2ptcls_state
    procedure          :: map_cls2D_flag_to_ptcls
    procedure          :: map_ptcls_state_to_cls
    procedure          :: prune_particles
    procedure          :: scale_projfile
    ! PTCL2D/3D - getters
    procedure          :: is_virgin_field
    procedure          :: count_state_gt_zero
    procedure          :: set_ptcl2D_thumb
    procedure          :: get_nptcls
    procedure          :: get_boxcoords
    procedure          :: set_boxcoords

    ! OUT
    procedure          :: add_cavgs2os_out
    procedure          :: add_frcs2os_out
    procedure          :: add_fsc2os_out
    procedure          :: add_vol2os_out
    procedure          :: add_entry2os_out
    procedure          :: remove_entry_from_osout
    ! OUT - Getters
    procedure          :: isthere_in_osout
    procedure          :: get_all_vols
    procedure          :: get_all_fscs
    procedure          :: get_vol
    procedure          :: get_fsc
    procedure          :: get_frcs
    procedure          :: get_imginfo_from_osout
    procedure          :: get_imgdims_from_osout

    ! OPTICS
    procedure          :: import_optics_map

    ! I/O – Printers
    procedure          :: print_info
    procedure          :: print_info_json
    procedure          :: print_segment
    procedure          :: print_segment_json
    ! I/O – Readers
    procedure          :: read
    procedure          :: read_mic_stk_ptcl2D_segments
    procedure          :: read_non_data_segments
    procedure          :: read_ctfparams_state_eo
    procedure          :: read_segment
    procedure, private :: segreader
    procedure          :: read_segments_info
    procedure          :: read_data_info
    ! I/O – Writers
    procedure          :: write
    procedure          :: write_segment_inside
    procedure          :: write_non_data_segments
    procedure          :: write_segment2txt
    procedure          :: write_mics_star
    procedure          :: write_ptcl2D_star
    procedure          :: write_optics_map
    procedure, private :: segwriter
    procedure          :: segwriter_inside
end type sp_project

interface

    ! CORE

    module subroutine new_seg_with_ptr( self, n, oritype, os_ptr )
        class(sp_project), target, intent(inout) :: self
        integer,                   intent(in)    :: n
        character(len=*),          intent(in)    :: oritype
        class(oris), pointer,      intent(inout) :: os_ptr
    end subroutine new_seg_with_ptr

    module subroutine kill( self )
        class(sp_project), intent(inout) :: self
    end subroutine kill

    module subroutine assign( self_out, self_in )
        class(sp_project), intent(inout) :: self_out
        class(sp_project), intent(in)    :: self_in
    end subroutine assign

    module subroutine copy( self_out, self_in )
        class(sp_project), intent(inout) :: self_out
        class(sp_project), intent(in)    :: self_in
    end subroutine copy

    ! CORE - field updaters

    module function oritype2segment(oritype) result(iseg)
        character(len=*), intent(in) :: oritype
        integer(kind(ENUM_ORISEG))   :: iseg
    end function oritype2segment

    module subroutine update_projinfo_1( self, cline )
        class(sp_project), intent(inout) :: self
        class(cmdline),    intent(in)    :: cline
    end subroutine update_projinfo_1

    module subroutine update_projinfo_2( self, projfile )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: projfile
    end subroutine update_projinfo_2

    module subroutine update_compenv( self, cline )
        class(sp_project), intent(inout) :: self
        class(cmdline),    intent(in)    :: cline
    end subroutine update_compenv

    module subroutine append_project( self1, self2 )
        class(sp_project), intent(inout) :: self1
        class(sp_project), intent(in)    :: self2
    end subroutine append_project

    module subroutine append_job_descr2jobproc( self, exec_dir, job_descr, did_update )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: exec_dir
        class(chash),      intent(inout) :: job_descr
        logical,           intent(out)   :: did_update
    end subroutine append_job_descr2jobproc

    module subroutine replace_project( self, projfile_src, oritype )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: projfile_src
        character(len=*),  intent(in)    :: oritype
    end subroutine replace_project

    !> for merging alignment documents from SIMPLE runs in distributed mode
    module subroutine merge_algndocs( self, nptcls, ndocs, oritype, fbody, numlen_in )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: nptcls, ndocs
        character(len=*),  intent(in)    :: oritype, fbody
        integer, optional, intent(in)    :: numlen_in
    end subroutine merge_algndocs

    module subroutine projrecords2proj( spproj, project_list )
        class(sp_project), intent(inout) :: spproj
        class(rec_list ),  intent(inout) :: project_list
    end subroutine projrecords2proj

    ! CORE - getters/setters

    module subroutine ptr2oritype( self, oritype, os_ptr )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        class(oris),      pointer, intent(inout) :: os_ptr
    end subroutine ptr2oritype

    module real function get_smpd( self )
        class(sp_project), target, intent(inout) :: self
    end function get_smpd

    module function get_ctfflag( self, oritype, iptcl )result( ctfflag )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        integer,         optional, intent(in)    :: iptcl
        character(len=STDLEN) :: ctfflag
    end function get_ctfflag

    module integer(kind(ENUM_CTFFLAG)) function get_ctfflag_type( self, oritype, iptcl )    
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        integer,         optional, intent(in)    :: iptcl
    end function get_ctfflag_type

    module function get_ctfparams( self, oritype, iptcl ) result( ctfvars )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        integer,                   intent(in)    :: iptcl
        class(oris), pointer  :: ptcl_field
        character(len=STDLEN) :: ctfflag, phaseplate
        type(ctfparams)       :: ctfvars
    end function get_ctfparams

    module integer function get_n_insegment( self, oritype )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
    end function get_n_insegment

    module integer function get_n_insegment_state( self, oritype, state )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        integer,                   intent(in)    :: state
    end function get_n_insegment_state

    module subroutine get_sp_oris( self, which_imgkind, os )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: which_imgkind
        class(oris),       intent(inout) :: os
    end subroutine get_sp_oris

    module subroutine set_sp_oris( self, which_imgkind, os )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: which_imgkind
        class(oris),       intent(inout) :: os
    end subroutine set_sp_oris

    ! MIC

    module subroutine add_single_movie( self, moviename, ctfvars )
        class(sp_project), target, intent(inout) :: self
        class(string),             intent(in)    :: moviename
        type(ctfparams),           intent(in)    :: ctfvars
    end subroutine add_single_movie

    module subroutine add_movies( self, movies_array, ctfvars, singleframe, verbose )
        class(sp_project), target, intent(inout) :: self
        type(string),              intent(in)    :: movies_array(:)
        type(ctfparams),           intent(in)    :: ctfvars
        logical,         optional, intent(in)    :: singleframe
        logical,         optional, intent(in)    :: verbose
    end subroutine add_movies

    module subroutine add_intgs( self, intgs_array, os, ctfvars )
        class(sp_project), target, intent(inout) :: self
        class(string),             intent(in)    :: intgs_array(:)
        class(oris),               intent(in)    :: os
        type(ctfparams),           intent(in)    :: ctfvars
    end subroutine add_intgs

    module subroutine report_state2mic( self, states )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: states(:)
    end subroutine report_state2mic

    ! MIC - getters

    module logical function has_boxfile( self )
        class(sp_project), target, intent(in) :: self
    end function has_boxfile

    module subroutine get_movies_table( self, moviestab )
        class(sp_project),         intent(inout) :: self
        type(string), allocatable, intent(out)   :: moviestab(:)
    end subroutine get_movies_table

    module subroutine get_mics_table( self, micstab, orimap)
        class(sp_project),         intent(inout) :: self
        type(string), allocatable, intent(out)   :: micstab(:)
        integer,      allocatable, intent(out)   :: orimap(:)
    end subroutine get_mics_table

    module function get_micparams( self, imic ) result( ctfvars )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: imic
        type(ctfparams) :: ctfvars
    end function get_micparams
    
    module function get_micname( self, iptcl ) result( micname )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: iptcl
        character(len=XLONGSTRLEN) :: micname
    end function get_micname

    module function get_mic_kind( self, imic ) result( mickind )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: imic
        type(string) :: mickind
    end function get_mic_kind

    module subroutine set_boxfile( self, i, boxfname, nptcls )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: i
        class(string),     intent(in)    :: boxfname
        integer, optional, intent(in)    :: nptcls
    end subroutine set_boxfile

    module integer function get_nmovies( self )
        class(sp_project), target, intent(inout) :: self
    end function get_nmovies

    module integer function get_nintgs( self )
        class(sp_project), target, intent(inout) :: self
    end function get_nintgs

    module integer function get_nframes( self )
        class(sp_project), target, intent(inout) :: self
    end function get_nframes

    module subroutine get_mic2stk_inds( self, mic2stk_inds, stk2mic_inds )
        class(sp_project),    intent(inout) :: self
        integer, allocatable, intent(inout) :: mic2stk_inds(:), stk2mic_inds(:)
    end subroutine get_mic2stk_inds

    ! STK

    module subroutine add_stk( self, stk, ctfvars )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: stk
        type(ctfparams),   intent(in)    :: ctfvars
    end subroutine add_stk

    module subroutine add_single_stk( self, stk, ctfvars, os )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: stk
        type(ctfparams),   intent(in)    :: ctfvars
        class(oris),       intent(inout) :: os
    end subroutine add_single_stk

    module subroutine add_stktab_1( self, stkfnames, os )
        class(sp_project), intent(inout) :: self
        class(string),     intent(inout) :: stkfnames(:)
        class(oris),       intent(inout) :: os
    end subroutine add_stktab_1

    module subroutine add_stktab_2( self, stkfnames, ctfvars, os )
        class(sp_project), intent(inout) :: self
        class(string),     intent(inout) :: stkfnames(:)
        type(ctfparams),   intent(in)    :: ctfvars
        class(oris),       intent(inout) :: os
    end subroutine add_stktab_2

    module subroutine split_stk( self, nparts, dir )
        class(sp_project),       intent(inout) :: self
        integer,                 intent(in)    :: nparts
        class(string), optional, intent(in)    :: dir
    end subroutine split_stk

    module subroutine write_substk( self, fromto, stkout )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: fromto(2)
        class(string),     intent(in)    :: stkout
    end subroutine write_substk

    module subroutine add_scale_tag( self, dir )
        class(sp_project),       intent(inout) :: self
        class(string), optional, intent(in)    :: dir
    end subroutine add_scale_tag

    module subroutine report_state2stk( self, states )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: states(:)
    end subroutine report_state2stk

    ! STK - Getters

    module logical function has_phaseplate( self, oritype )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
    end function has_phaseplate

    module integer function get_box( self )
        class(sp_project), target, intent(in) :: self
    end function get_box

    module function get_stkname( self, imic ) result( stkname )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: imic
        type(string) :: stkname
    end function get_stkname

    module subroutine get_stkname_and_ind( self, oritype, iptcl, stkname, ind_in_stk )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        integer,                   intent(in)    :: iptcl
        class(string),             intent(out)   :: stkname
        integer,                   intent(out)   :: ind_in_stk
    end subroutine get_stkname_and_ind

    module function get_nstks( self ) result(nstks)
        class(sp_project), target, intent(in) :: self
        integer                               :: nstks
    end function get_nstks

    ! CLS2D/3D

    module function get_selected_clsinds( self ) result( clsinds )
        class(sp_project), intent(inout) :: self
        integer,             allocatable :: clsinds(:)
    end function get_selected_clsinds

    module subroutine set_cavgs_thumb( self, projfile )
        class(sp_project),  intent(inout) :: self
        class(string),      intent(in)    :: projfile
    end subroutine set_cavgs_thumb

    module subroutine shape_ranked_cavgs2jpg( self, cavg_inds, jpgname, xtiles, ytiles, mskdiam_px )
        class(sp_project),    intent(inout) :: self
        integer, allocatable, intent(inout) :: cavg_inds(:)
        class(string),        intent(in)    :: jpgname
        integer,              intent(out)   :: xtiles, ytiles
        integer, optional,    intent(in)    :: mskdiam_px
    end subroutine shape_ranked_cavgs2jpg

    module subroutine cavgs2jpg( self, cavg_inds, jpgname, xtiles, ytiles )
        class(sp_project),    intent(inout) :: self
        integer, allocatable, intent(inout) :: cavg_inds(:)
        class(string),        intent(in)    :: jpgname
        integer,              intent(out)   :: xtiles, ytiles
    end subroutine cavgs2jpg

    module subroutine cavgs2mrc( self )
        class(sp_project), intent(inout) :: self
    end subroutine cavgs2mrc

    ! PTCL2D/3D

    module subroutine map_ptcl_ind2stk_ind( self, oritype, iptcl, stkind, ind_in_stk )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
        integer,                   intent(in)    :: iptcl
        integer,                   intent(out)   :: stkind
        integer,                   intent(out)   :: ind_in_stk
    end subroutine map_ptcl_ind2stk_ind

    module subroutine map_cavgs_selection( self, states )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: states(:)
    end subroutine map_cavgs_selection

    module subroutine map2Dshifts23D( self )
        class(sp_project), intent(inout) :: self
    end subroutine map2Dshifts23D

    module subroutine map2ptcls( self )
        class(sp_project), intent(inout) :: self
    end subroutine map2ptcls

    module subroutine map2ptcls_state( self, append, maxpop )
        class(sp_project), intent(inout) :: self
        logical, optional, intent(in)    :: append
        integer, optional, intent(in)    :: maxpop
    end subroutine map2ptcls_state

    module subroutine map_cls2D_flag_to_ptcls( self, flag )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: flag
    end subroutine map_cls2D_flag_to_ptcls

    module subroutine map_ptcls_state_to_cls( self )
        class(sp_project), intent(inout) :: self
    end subroutine map_ptcls_state_to_cls

    module subroutine prune_particles( self )
        class(sp_project), target, intent(inout) :: self
    end subroutine prune_particles

    module subroutine scale_projfile( self, smpd_target, new_projfile, cline, cline_scale, dir )
        ! this probably needs an oritype input for dealing with scale class averages
        class(sp_project),       intent(inout) :: self
        real,                    intent(inout) :: smpd_target
        class(string),           intent(out)   :: new_projfile
        class(cmdline),          intent(inout) :: cline
        class(cmdline),          intent(out)   :: cline_scale
        class(string), optional, intent(in)    :: dir
    end subroutine scale_projfile

    ! PTCL2D/3D - getters

    module logical function is_virgin_field( self, oritype )
        class(sp_project), target, intent(inout) :: self
        character(len=*),          intent(in)    :: oritype
    end function is_virgin_field

    module integer function count_state_gt_zero( self )
        class(sp_project), target, intent(inout) :: self
    end function count_state_gt_zero

    module subroutine set_ptcl2D_thumb( self, projfile, indices, boxsize )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: projfile
        integer,           intent(in)    :: indices(:)
        integer,           intent(out)   :: boxsize
    end subroutine set_ptcl2D_thumb

    module integer function get_nptcls( self )
        class(sp_project), target, intent(in) :: self
    end function get_nptcls

    module subroutine get_boxcoords( self, iptcl, coords )
        class(sp_project), target, intent(in)  :: self
        integer,                   intent(in)  :: iptcl
        integer,                   intent(out) :: coords(2)
    end subroutine get_boxcoords

    module subroutine set_boxcoords( self, iptcl, coords )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: iptcl, coords(2)
    end subroutine set_boxcoords

    ! OUT

    module subroutine add_cavgs2os_out( self, stk, smpd, imgkind, clspath )
        class(sp_project),          intent(inout) :: self
        class(string),              intent(in)    :: stk
        real,                       intent(in)    :: smpd
        character(len=*), optional, intent(in)    :: imgkind
        logical,          optional, intent(in)    :: clspath
    end subroutine add_cavgs2os_out

    module subroutine add_frcs2os_out( self, frc, which_imgkind )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: frc
        character(len=*),  intent(in)    :: which_imgkind
    end subroutine add_frcs2os_out

    module subroutine add_fsc2os_out( self, fsc, state, box )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: fsc
        integer,           intent(in)    :: state, box
    end subroutine add_fsc2os_out

    module subroutine add_vol2os_out( self, vol, smpd, state, which_imgkind, box, pop )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: vol
        character(len=*),  intent(in)    :: which_imgkind
        real,              intent(in)    :: smpd
        integer,           intent(in)    :: state
        integer, optional, intent(in)    :: box, pop
    end subroutine add_vol2os_out

    module subroutine add_entry2os_out( self, which_imgkind, ind )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: which_imgkind
        integer,           intent(out)   :: ind
    end subroutine add_entry2os_out

    module subroutine remove_entry_from_osout( self, which_imgkind, state )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: which_imgkind
        integer,           intent(in)    :: state
    end subroutine remove_entry_from_osout

    ! OUT - Getters

    module logical function isthere_in_osout( self, which_imgkind, state )
        class(sp_project), intent(in) :: self
        character(len=*),  intent(in) :: which_imgkind
        integer,           intent(in) :: state
    end function isthere_in_osout

    module subroutine get_cavgs_stk( self, stkname, ncls, smpd, imgkind, fail, out_ind, box )
        class(sp_project),          intent(inout) :: self
        class(string),              intent(inout) :: stkname
        integer,                    intent(out)   :: ncls
        real,                       intent(out)   :: smpd
        character(len=*), optional, intent(in)    :: imgkind
        logical,          optional, intent(in)    :: fail
        integer,          optional, intent(inout) :: out_ind
        real,             optional, intent(inout) :: box
    end subroutine get_cavgs_stk

    module subroutine get_vol( self, imgkind, state, vol_fname, smpd, box )
        class(sp_project), intent(in)    :: self
        character(len=*),  intent(in)    :: imgkind
        integer,           intent(in)    :: state
        type(string),      intent(inout) :: vol_fname
        real,              intent(out)   :: smpd
        integer,           intent(out)   :: box
    end subroutine get_vol

    module subroutine get_all_vols( self, orisout )
        class(sp_project), intent(in)    :: self
        type(oris),        intent(inout) :: orisout
    end subroutine get_all_vols

    module subroutine get_fsc( self, state, fsc_fname, box )
        class(sp_project), intent(in)    :: self
        integer,           intent(in)    :: state
        class(string),     intent(inout) :: fsc_fname
        integer,           intent(out)   :: box
    end subroutine get_fsc

    module subroutine get_all_fscs( self, orisout )
        class(sp_project), intent(in)    :: self
        type(oris),        intent(inout) :: orisout
    end subroutine get_all_fscs

    module subroutine get_frcs( self, frcs, which_imgkind, fail )
        class(sp_project), intent(in)    :: self
        class(string),     intent(inout) :: frcs
        character(len=*),  intent(in)    :: which_imgkind
        logical, optional, intent(in)    :: fail

    end subroutine get_frcs

    module subroutine get_imginfo_from_osout( self, smpd, box, nptcls )
        class(sp_project), intent(inout) :: self
        real,              intent(out)   :: smpd
        integer,           intent(out)   :: box, nptcls
    end subroutine get_imginfo_from_osout

    module subroutine get_imgdims_from_osout( self, iseg, smpd, box )
        class(sp_project), intent(inout) :: self
        integer,           intent(in)    :: iseg
        real,              intent(out)   :: smpd
        integer,           intent(out)   :: box
    end subroutine get_imgdims_from_osout

    ! OPTICS

    module subroutine import_optics_map( self, mapfileprefix )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: mapfileprefix
        type(nrtxtfile)      :: mapfile
        real,    allocatable :: map_entries(:,:)
        integer, allocatable :: mics_optics_map(:)
        real                 :: min_importind, max_importind
        integer              :: il, nl, imic, iptcl, importind, stkind
    end subroutine import_optics_map

    ! I/O - Printers

    module subroutine print_info( self, fname )
        class(sp_project),     intent(inout) :: self
        class(string),         intent(in)    :: fname
    end subroutine print_info

    module subroutine print_info_json( self, fname )
        class(sp_project), intent(inout)   :: self
        class(string),     intent(in)      :: fname
    end subroutine print_info_json

    module subroutine print_segment( self, oritype, fromto )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: oritype
        integer, optional, intent(in)    :: fromto(2)
    end subroutine print_segment

    module subroutine print_segment_json( self, oritype, projfile, fromto, sort_key, sort_asc, hist, nran, boxes, plot_key )
        class(sp_project),           intent(inout) :: self
        character(len=*),            intent(in)    :: oritype
        class(string),               intent(in)    :: projfile
        character(len=*),  optional, intent(in)    :: sort_key, sort_asc, hist, plot_key
        integer,           optional, intent(in)    :: fromto(2), nran
        logical,           optional, intent(in)    :: boxes
    end subroutine print_segment_json

    ! I/O - Readers

    module subroutine read( self, fname, wthreads )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: fname
        logical, optional, intent(in)    :: wthreads
    end subroutine read

    module subroutine read_non_data_segments( self, fname )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: fname
    end subroutine read_non_data_segments

    module subroutine read_ctfparams_state_eo( self, fname )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: fname
    end subroutine read_ctfparams_state_eo

    module subroutine read_mic_stk_ptcl2D_segments( self, fname, wthreads )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: fname
        logical, optional, intent(in)    :: wthreads
    end subroutine read_mic_stk_ptcl2D_segments

    module subroutine read_segment( self, oritype, fname, fromto, wthreads )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: oritype
        class(string),     intent(in)    :: fname
        integer, optional, intent(in)    :: fromto(2)
        logical, optional, intent(in)    :: wthreads
    end subroutine read_segment

    module subroutine segreader( self, isegment, fromto, only_ctfparams_state_eo, wthreads )
        class(sp_project),          intent(inout) :: self
        integer(kind(ENUM_ORISEG)), intent(in)    :: isegment
        integer, optional,          intent(in)    :: fromto(2)
        logical, optional,          intent(in)    :: only_ctfparams_state_eo, wthreads
    end subroutine segreader

    module subroutine read_segments_info( self, fname, seginds, seginfos )
        class(sp_project),                  intent(inout) :: self
        class(string),                      intent(in)  :: fname
        type(binoris_seginfo), allocatable, intent(out) :: seginfos(:)
        integer,               allocatable, intent(out) :: seginds(:)
    end subroutine read_segments_info

    module subroutine read_data_info( self, fname, nmics, nstks, nptcls )
        class(sp_project),   intent(inout) :: self
        class(string),       intent(in)    :: fname
        integer,             intent(out)   :: nmics, nstks, nptcls
    end subroutine read_data_info

    ! I/O - Writers

    module subroutine write( self, fname, fromto, isegment, tempfile )
        class(sp_project),                    intent(inout) :: self
        class(string),              optional, intent(in)    :: fname
        integer,                    optional, intent(in)    :: fromto(2)
        integer(kind(ENUM_ORISEG)), optional, intent(in)    :: isegment
        logical,                    optional, intent(in)    :: tempfile
    end subroutine write

    module subroutine write_segment_inside( self, oritype, fname, fromto )
        class(sp_project),          intent(inout) :: self
        character(len=*),           intent(in)    :: oritype
        class(string),    optional, intent(in)    :: fname
        integer,          optional, intent(in)    :: fromto(2)
    end subroutine write_segment_inside

    module subroutine write_non_data_segments( self, fname )
        class(sp_project), intent(inout) :: self
        class(string),     intent(in)    :: fname
    end subroutine write_non_data_segments

    module subroutine write_segment2txt( self, oritype, fname, fromto )
        class(sp_project), intent(inout) :: self
        character(len=*),  intent(in)    :: oritype
        class(string),     intent(in)    :: fname
        integer, optional, intent(in)    :: fromto(2)
    end subroutine write_segment2txt

    module subroutine segwriter( self, isegment, fromto )
        class(sp_project),          intent(inout) :: self
        integer(kind(ENUM_ORISEG)), intent(in)    :: isegment
        integer, optional,         intent(in)     :: fromto(2)
    end subroutine segwriter

    module subroutine segwriter_inside( self, isegment, fromto )
        class(sp_project),          intent(inout) :: self
        integer(kind(ENUM_ORISEG)), intent(in)    :: isegment
        integer, optional,          intent(in)    :: fromto(2)
    end subroutine segwriter_inside

    module subroutine write_mics_star( self, fname )
        class(sp_project),       intent(inout) :: self
        class(string), optional, intent(in)    :: fname
    end subroutine write_mics_star

    module subroutine write_ptcl2D_star( self, fname )
        class(sp_project),       intent(inout) :: self
        class(string), optional, intent(in)    :: fname
    end subroutine write_ptcl2D_star

    module subroutine write_optics_map( self, fname_prefix )
        class(sp_project),          intent(inout) :: self
        character(len=*),           intent(in)    :: fname_prefix
    end subroutine write_optics_map

end interface

contains

    !------ Non-type-bound helpers ------

    module integer(kind(ENUM_ORISEG)) function oritype2segment( oritype )
        character(len=*), intent(in) :: oritype
        select case(trim(oritype))
            case('mic')
                oritype2segment = MIC_SEG
            case('stk')
                oritype2segment = STK_SEG
            case('ptcl2D')
                oritype2segment = PTCL2D_SEG
            case('cls2D')
                oritype2segment = CLS2D_SEG
            case('cls3D')
                oritype2segment = CLS3D_SEG
            case('ptcl3D')
                oritype2segment = PTCL3D_SEG
            case('out')
                oritype2segment = OUT_SEG
            case('optics')
                oritype2segment = OPTICS_SEG
            case('projinfo')
                oritype2segment = PROJINFO_SEG
            case('jobproc')
                oritype2segment = JOBPROC_SEG
            case('compenv')
                oritype2segment = COMPENV_SEG
            case DEFAULT
                THROW_HARD('unsupported oritype flag; oritype2segment')
            end select
    end function oritype2segment

end module simple_sp_project
