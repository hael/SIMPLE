! 3D reconstruction - master module
module simple_rec_master
include 'simple_lib.f08'
use simple_builder,    only: build_glob
use simple_parameters, only: params_glob
use simple_qsys_funs,  only: qsys_job_finished
implicit none

public :: exec_rec, exec_rec_soft
private
#include "simple_local_flags.inc"

contains

    subroutine exec_rec
        character(len=:), allocatable :: fbody
        integer :: s
        if( .not. params_glob%l_distr_exec ) THROW_HARD('eo .ne. no not supported here, use simple_distr_exec!')
        ! rebuild build_glob%vol according to box size (beacuse it is otherwise boxmatch)
        call build_glob%vol%new([params_glob%box,params_glob%box,params_glob%box], params_glob%smpd)
        do s=1,params_glob%nstates
            if( build_glob%spproj_field%get_pop(s, 'state') == 0 ) cycle ! empty state
            fbody = 'recvol_state'
            call build_glob%eorecvol%eorec_distr( build_glob%spproj, build_glob%spproj_field, build_glob%pgrpsyms, s, fbody=fbody)
        end do
        call qsys_job_finished( 'simple_rec_master :: exec_eorec')
        write(logfhandle,'(a,1x,a)') "GENERATED VOLUMES: reconstruct3D*.ext"
    end subroutine exec_rec

    subroutine exec_rec_soft( cline, which_iter )
        use simple_strategy3D_matcher
        use simple_cmdline, only: cmdline
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: which_iter
        call setup_weights_read_o_peaks
        call calc_global_ori_weights
        call calc_proj_weights
        call calc_3Drec( cline, which_iter )
        call qsys_job_finished( 'simple_rec_master :: exec_rec_soft')
    end subroutine exec_rec_soft

end module simple_rec_master
