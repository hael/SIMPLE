module simple_sp_project_handler
#include "simple_lib.f08"
use simple_sp_project, only: sp_project
implicit none

public :: transfer_sp_project_segment, print_sp_project_info
private

contains

    subroutine transfer_sp_project_segment( fname_provider, fname_reciever, which_segment  )
        character(len=*), intent(in) :: fname_provider, fname_reciever, which_segment
        type(sp_project) :: sp_provider, sp_reciever
        select case(trim(which_segment))
            case('stk')
                call sp_provider%read_segment(fname_provider, STK_SEG)
                sp_reciever%os_stk = sp_provider%os_stk
                call sp_reciever%write_segment(fname_reciever, STK_SEG)
            case('ptcl2D')
                call sp_provider%read_segment(fname_provider, PTCL2D_SEG)
                sp_reciever%os_ptcl2D = sp_provider%os_ptcl2D
                call sp_reciever%write_segment(fname_reciever, PTCL2D_SEG)
            case('cls2D')
                call sp_provider%read_segment(fname_provider, CLS2D_SEG)
                sp_reciever%os_cls2D = sp_provider%os_cls2D
                call sp_reciever%write_segment(fname_reciever, CLS2D_SEG)
            case('cls3D')
                call sp_provider%read_segment(fname_provider, CLS3D_SEG)
                sp_reciever%os_cls3D = sp_provider%os_cls3D
                call sp_reciever%write_segment(fname_reciever, CLS3D_SEG)
            case('ptcl3D')
                call sp_provider%read_segment(fname_provider, PTCL3D_SEG)
                sp_reciever%os_ptcl3D = sp_provider%os_ptcl3D
                call sp_reciever%write_segment(fname_reciever, PTCL3D_SEG)
            case('projinfo')
                call sp_provider%read_segment(fname_provider, PROJINFO_SEG)
                sp_reciever%projinfo = sp_provider%projinfo
                call sp_reciever%write_segment(fname_reciever, PROJINFO_SEG)
            case('jobproc')
                call sp_provider%read_segment(fname_provider, JOBPROC_SEG)
                sp_reciever%jobproc = sp_provider%jobproc
                call sp_reciever%write_segment(fname_reciever, JOBPROC_SEG)
            case DEFAULT
                stop 'unsupported which flag; sp_project_handler :: transfer_segment'
        end select
    end subroutine transfer_sp_project_segment

    subroutine print_sp_project_info( fname )
        character(len=*), intent(in) :: fname
        type(sp_project) :: sp_proj
        call sp_proj%read(fname)
        call sp_proj%print_info
        call sp_proj%kill
    end subroutine print_sp_project_info

end module simple_sp_project_handler
