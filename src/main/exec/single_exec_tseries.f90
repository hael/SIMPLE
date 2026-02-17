module single_exec_tseries
use simple_cmdline,               only: cmdline
use single_commanders_trajectory, only: commander_track_particles_distr
use single_commanders_tseries,    only: commander_tseries_import, commander_tseries_make_pickavg,&
&commander_tseries_motion_correct_distr
implicit none

public :: exec_tseries_commander
private

type(commander_track_particles_distr)        :: xtrack
type(commander_tseries_import)               :: xtseries_import
type(commander_tseries_make_pickavg)         :: xtseries_make_pickavg
type(commander_tseries_motion_correct_distr) :: xmcorr

contains
    
    subroutine exec_tseries_commander(which, cline, l_silent, l_did_execute)
        character(len=*), intent(in)    :: which
        class(cmdline),   intent(inout) :: cline
        logical,          intent(inout) :: l_did_execute
        logical,          intent(out)   :: l_silent
        if( l_did_execute )return
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(which))
            case( 'track_particles' )
                call xtrack%execute( cline )
            case( 'tseries_import' )
                call xtseries_import%execute(cline)
            case( 'tseries_make_pickavg' )
                call xtseries_make_pickavg%execute(cline)
            case( 'tseries_motion_correct' )
                call xmcorr%execute( cline )
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_tseries_commander

end module single_exec_tseries
