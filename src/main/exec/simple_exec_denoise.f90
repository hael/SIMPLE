!@descr: execution of denoising commanders
module simple_exec_denoise
use simple_cmdline,              only: cmdline
use simple_commanders_resolest,  only: commander_icm2D, commander_icm3D
use simple_commanders_volops,    only: commander_ppca_volvar
use simple_commanders_cluster2D, only: commander_ppca_denoise_classes
use simple_commanders_imgops,    only: commander_ppca_denoise
implicit none

public :: exec_denoise_commander
private

type(commander_icm2D)                :: xicm2D
type(commander_icm3D)                :: xicm3D
type(commander_ppca_denoise)         :: xppca_denoise
type(commander_ppca_denoise_classes) :: xppca_denoise_classes
type(commander_ppca_volvar)          :: xppca_volvar

contains

     subroutine exec_denoise_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(out)   :: l_silent
        logical,             intent(out)   :: l_did_execute
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'icm2D' )
                call xicm2D%execute(cline)
            case( 'icm3D' )
                call xicm3D%execute(cline)
            case( 'ppca_denoise' )
                call xppca_denoise%execute(cline)
            case( 'ppca_denoise_classes' )
                call xppca_denoise_classes%execute(cline)
            case( 'ppca_volvar' )
                call xppca_volvar%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_denoise_commander

end module simple_exec_denoise



