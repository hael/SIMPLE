! 3D reconstruction - master module
module simple_rec_master
include 'simple_lib.f08'
use simple_builder,    only: build_glob
use simple_parameters, only: params_glob
use simple_qsys_funs,  only: qsys_job_finished
implicit none

public :: exec_rec_master, exec_rec_soft
private
#include "simple_local_flags.inc"

contains

    subroutine exec_rec_master
        select case(params_glob%eo)
            case( 'no' )
                call exec_rec
            case DEFAULT
                call exec_eorec_distr
        end select
    end subroutine exec_rec_master

    subroutine exec_rec
        character(len=:), allocatable :: fbody
        character(len=STDLEN)         :: rho_name
        integer                       :: s
        ! rebuild build_glob%vol according to box size (beacuse it is otherwise boxmatch)
        call build_glob%vol%new([params_glob%box,params_glob%box,params_glob%box], params_glob%smpd)
        do s=1,params_glob%nstates
            if( build_glob%spproj_field%get_pop(s, 'state') == 0 ) cycle ! empty state
            if( params_glob%l_distr_exec )then ! embarrasingly parallel rec
                allocate(fbody, source='recvol_state'//int2str_pad(s,2)//&
                &'_part'//int2str_pad(params_glob%part,params_glob%numlen))
                params_glob%vols(s) = fbody//params_glob%ext
                rho_name      = 'rho_'//fbody//params_glob%ext
                call build_glob%recvol%rec( build_glob%spproj, build_glob%spproj_field, build_glob%pgrpsyms, s, part=params_glob%part)
                call build_glob%recvol%compress_exp
                call build_glob%recvol%write(params_glob%vols(s), del_if_exists=.true.)
                call build_glob%recvol%write_rho(trim(rho_name))
            else ! shared-mem parallel rec
                allocate(fbody, source='recvol_state')
                params_glob%vols(s) = fbody//int2str_pad(s,2)//params_glob%ext
                call build_glob%recvol%rec( build_glob%spproj, build_glob%spproj_field, build_glob%pgrpsyms, s)
                call build_glob%recvol%clip(build_glob%vol)
                call build_glob%vol%write(params_glob%vols(s), del_if_exists=.true.)
            endif
        end do
        write(logfhandle,'(a)') "GENERATED VOLUMES: reconstruct3D*.ext"
        call qsys_job_finished(  'simple_rec_master :: exec_rec')
    end subroutine exec_rec

    subroutine exec_eorec_distr
        character(len=:), allocatable :: fbody
        integer :: s
        if( .not. params_glob%l_distr_exec ) THROW_HARD('eo .ne. no not supported here, use simple_distr_exec!')
        ! rebuild build_glob%vol according to box size (beacuse it is otherwise boxmatch)
        call build_glob%vol%new([params_glob%box,params_glob%box,params_glob%box], params_glob%smpd)
        do s=1,params_glob%nstates
            DebugPrint  'processing state: ', s
            if( build_glob%spproj_field%get_pop(s, 'state') == 0 ) cycle ! empty state
            allocate(fbody, source='recvol_state')
            call build_glob%eorecvol%eorec_distr( build_glob%spproj, build_glob%spproj_field, build_glob%pgrpsyms, s, fbody=fbody)
        end do
        call qsys_job_finished( 'simple_rec_master :: exec_eorec')
        write(logfhandle,'(a,1x,a)') "GENERATED VOLUMES: reconstruct3D*.ext"
    end subroutine exec_eorec_distr

    subroutine exec_rec_soft( cline, which_iter )
        use simple_strategy3D_matcher, only: setup_weights_read_o_peaks, calc_3Drec
        use simple_cmdline,            only: cmdline
        class(cmdline), intent(inout) :: cline
        integer,        intent(in)    :: which_iter
        call setup_weights_read_o_peaks
        call calc_3Drec( cline, which_iter )
    end subroutine exec_rec_soft

end module simple_rec_master
