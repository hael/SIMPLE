module simple_sp_project_handler
#include "simple_lib.f08"
use simple_sp_project, only: sp_project
implicit none

public :: transfer_sp_project_segment, print_sp_project_info
private

contains

    subroutine transfer_sp_project_segment( fname_provider, fname_reciever, which_segment )
        character(len=*), intent(in) :: fname_provider, fname_reciever, which_segment
        type(sp_project) :: sp_provider, sp_reciever
        call sp_provider%read_segment(which_segment, fname_provider)
        call sp_reciever%jobproc%write(fname_reciever)
    end subroutine transfer_sp_project_segment

    subroutine print_sp_project_info( fname )
        character(len=*), intent(in) :: fname
        type(sp_project) :: sp_proj
        call sp_proj%read(fname)
        call sp_proj%print_info
        call sp_proj%kill
    end subroutine print_sp_project_info

end module simple_sp_project_handler
