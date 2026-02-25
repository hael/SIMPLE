!@descr: aggregates SIMPLE (non-stream) ui program constructors
module simple_ui_simple_group
use simple_ui_hash,      only: ui_hash
! SIMPLE program constructor interfaces
use simple_ui_project,   only: construct_project_programs,   print_project_programs
use simple_ui_preproc,   only: construct_preproc_programs,   print_preproc_programs
use simple_ui_cluster2D, only: construct_cluster2D_programs, print_cluster2D_programs
use simple_ui_cavgproc,  only: construct_cavgproc_programs,  print_cavgproc_programs
use simple_ui_abinitio3D,only: construct_abinitio3D_programs,print_abinitio3D_programs
use simple_ui_refine3D,  only: construct_refine3D_programs,  print_refine3D_programs
use simple_ui_denoise,   only: construct_denoise_programs,   print_denoise_programs
use simple_ui_filter,    only: construct_filter_programs,    print_filter_programs
use simple_ui_image,     only: construct_image_programs,     print_image_programs
use simple_ui_mask,      only: construct_mask_programs,      print_mask_programs
use simple_ui_ori,       only: construct_ori_programs,       print_ori_programs
use simple_ui_print,     only: construct_print_programs,     print_print_programs
use simple_ui_res,       only: construct_res_programs,       print_res_programs
use simple_ui_sim,       only: construct_sim_programs,       print_sim_programs
use simple_ui_validate,  only: construct_validate_programs,  print_validate_programs
use simple_ui_sym,       only: construct_symmetry_programs,  print_symmetry_programs
use simple_ui_dock,      only: construct_dock_programs,      print_dock_programs
use simple_ui_volume,    only: construct_volume_programs,    print_volume_programs
use simple_ui_other,     only: construct_other_programs,     print_other_programs
implicit none

public :: add_simple_programs, print_simple_programs
private


contains

    subroutine add_simple_programs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call construct_project_programs(prgtab)
        call construct_preproc_programs(prgtab)
        call construct_cluster2D_programs(prgtab)
        call construct_cavgproc_programs(prgtab)
        call construct_abinitio3D_programs(prgtab)
        call construct_refine3D_programs(prgtab)
        call construct_denoise_programs(prgtab)
        call construct_filter_programs(prgtab)
        call construct_image_programs(prgtab)
        call construct_mask_programs(prgtab)
        call construct_ori_programs(prgtab)
        call construct_print_programs(prgtab)
        call construct_res_programs(prgtab)
        call construct_sim_programs(prgtab)
        call construct_validate_programs(prgtab)
        call construct_symmetry_programs(prgtab)
        call construct_dock_programs(prgtab)
        call construct_volume_programs(prgtab)
        call construct_other_programs(prgtab)
    end subroutine add_simple_programs


    subroutine print_simple_programs( logfhandle )
        integer, intent(in) :: logfhandle
        call print_project_programs(logfhandle)
        call print_preproc_programs(logfhandle)
        call print_cluster2D_programs(logfhandle)
        call print_cavgproc_programs(logfhandle)
        call print_abinitio3D_programs(logfhandle)
        call print_refine3D_programs(logfhandle)
        call print_denoise_programs(logfhandle)
        call print_filter_programs(logfhandle)
        call print_image_programs(logfhandle)
        call print_mask_programs(logfhandle)
        call print_ori_programs(logfhandle)
        call print_print_programs(logfhandle)
        call print_res_programs(logfhandle)
        call print_sim_programs(logfhandle)
        call print_validate_programs(logfhandle)
        call print_symmetry_programs(logfhandle)
        call print_dock_programs(logfhandle)
        call print_volume_programs(logfhandle)
        call print_other_programs(logfhandle)
    end subroutine print_simple_programs

end module simple_ui_simple_group
