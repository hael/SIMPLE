!@descr: for all highlevel tests
module simple_commanders_test_highlevel
use simple_commanders_api
use simple_stream_api
use simple_commanders_project_core, only: commander_new_project
use simple_commanders_project_mov,  only: commander_import_movies
use simple_commanders_reproject,    only: commander_reproject
use simple_commanders_pick,         only: commander_pick, commander_extract
use simple_commanders_sim,          only: commander_simulate_particles, commander_simulate_movie
use simple_commanders_preprocess,   only: commander_ctf_estimate, commander_motion_correct
use simple_commanders_solve2D,      only: commander_solve2D
use simple_commanders_solve3D,      only: commander_solve3D
use simple_test_utils,              only: set_fixed_seed
use simple_commanders_validate,     only: commander_mini_stream
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_mini_stream
  contains
    procedure :: execute      => exec_test_mini_stream_quantitative
end type commander_test_mini_stream

type, extends(commander_base) :: commander_test_simulate_particles
  contains
    procedure :: execute      => exec_test_simulate_particles
end type commander_test_simulate_particles

type, extends(commander_base) :: commander_test_solve3D_addon
  contains
    procedure :: execute      => exec_test_solve3D_addon
end type commander_test_solve3D_addon

type, extends(commander_base) :: commander_generate_solve3D_addon_snapshots
    contains
        procedure :: execute      => exec_generate_solve3D_addon_snapshots
end type commander_generate_solve3D_addon_snapshots

type, extends(commander_base) :: commander_test_simulated_workflow
  contains
    procedure :: execute      => exec_test_simulated_workflow
end type commander_test_simulated_workflow

type, extends(commander_base) :: commander_test_pcg_recon
  contains
    procedure :: execute      => exec_test_pcg_recon
end type commander_test_pcg_recon

type, extends(commander_base) :: commander_test_pcg_frac_update
  contains
    procedure :: execute      => exec_test_pcg_frac_update
end type commander_test_pcg_frac_update

type, extends(commander_base) :: commander_test_rec3D_backends
  contains
    procedure :: execute      => exec_test_rec3D_backends
end type commander_test_rec3D_backends

type, extends(commander_base) :: commander_test_cont_refine3D_1jxy
  contains
    procedure :: execute      => exec_test_cont_refine3D_1jxy
end type commander_test_cont_refine3D_1jxy

type, extends(commander_base) :: commander_test_flex_pca_blobs
  contains
    procedure :: execute      => exec_test_flex_pca_blobs
end type commander_test_flex_pca_blobs

interface

    module subroutine exec_test_mini_stream_quantitative( self, cline )
        class(commander_test_mini_stream), intent(inout) :: self
        class(cmdline), intent(inout) :: cline
    end subroutine exec_test_mini_stream_quantitative

    module subroutine exec_test_simulate_particles( self, cline )
        class(commander_test_simulate_particles), intent(inout) :: self
        class(cmdline), intent(inout) :: cline
    end subroutine exec_test_simulate_particles

    module subroutine exec_test_simulated_workflow( self, cline )
        class(commander_test_simulated_workflow), intent(inout) :: self
        class(cmdline), intent(inout) :: cline
    end subroutine exec_test_simulated_workflow

    module subroutine exec_test_pcg_recon( self, cline )
        class(commander_test_pcg_recon), intent(inout) :: self
        class(cmdline), intent(inout) :: cline
    end subroutine exec_test_pcg_recon

    module subroutine exec_test_pcg_frac_update( self, cline )
        class(commander_test_pcg_frac_update), intent(inout) :: self
        class(cmdline), intent(inout) :: cline
    end subroutine exec_test_pcg_frac_update

    module subroutine exec_test_rec3D_backends( self, cline )
        class(commander_test_rec3D_backends), intent(inout) :: self
        class(cmdline), intent(inout) :: cline
    end subroutine exec_test_rec3D_backends

    module subroutine exec_test_solve3D_addon( self, cline )
        class(commander_test_solve3D_addon), intent(inout) :: self
        class(cmdline), intent(inout) :: cline
    end subroutine exec_test_solve3D_addon

    module subroutine exec_generate_solve3D_addon_snapshots( self, cline )
        class(commander_generate_solve3D_addon_snapshots), intent(inout) :: self
        class(cmdline), intent(inout) :: cline
    end subroutine exec_generate_solve3D_addon_snapshots

    module subroutine exec_test_cont_refine3D_1jxy( self, cline )
        class(commander_test_cont_refine3D_1jxy), intent(inout) :: self
        class(cmdline), intent(inout) :: cline
    end subroutine exec_test_cont_refine3D_1jxy

    module subroutine exec_test_flex_pca_blobs( self, cline )
        class(commander_test_flex_pca_blobs), intent(inout) :: self
        class(cmdline), intent(inout) :: cline
    end subroutine exec_test_flex_pca_blobs

end interface

end module simple_commanders_test_highlevel
