!@descr: appends the micrograph segments of finished upstream projects to a running micrograph segment
!==============================================================================
! MODULE: simple_mic_import
!
! PURPOSE:
!   The import step every stream stage repeats: read the 'mic' segment of each
!   newly finished upstream project and append its micrographs, in file order,
!   to the stage's global segment. Each file is read once, and the global
!   segment grows once per call. Whether rejected micrographs come along is
!   the caller's explicit choice (the optics stage keeps them; preprocessing
!   keeps only accepted ones on restart).
!
! HOME:
!   In src/main/stream/shared for now; it belongs with the project-file helpers
!   (simple_projfile_utils).
!
! TESTS:
!   simple_mic_import_tester
!==============================================================================
module simple_mic_import
use simple_string,     only: string
use simple_oris,       only: oris
use simple_sp_project, only: sp_project
implicit none

public :: append_mics_from_projects
private

contains

    !> Appends to @p os_mic the micrographs of the project files @p fnames, all of them or,
    !! with @p l_accepted_only, those with state > 0; @p nappended is the number appended.
    subroutine append_mics_from_projects( os_mic, fnames, l_accepted_only, nappended )
        class(oris),   intent(inout) :: os_mic
        class(string), intent(in)    :: fnames(:)
        logical,       intent(in)    :: l_accepted_only
        integer,       intent(out)   :: nappended
        type(sp_project), allocatable :: parts(:)
        integer :: nfiles, ifile, imic, n_old, j
        nappended = 0
        nfiles    = size(fnames)
        if( nfiles == 0 ) return
        allocate(parts(nfiles))
        do ifile = 1,nfiles
            call parts(ifile)%read_segment('mic', fnames(ifile))
            do imic = 1,parts(ifile)%os_mic%get_noris()
                if( selected(parts(ifile)%os_mic, imic) ) nappended = nappended + 1
            enddo
        enddo
        if( nappended > 0 )then
            n_old = os_mic%get_noris()
            if( n_old == 0 )then
                call os_mic%new(nappended, is_ptcl=.false.)
            else
                call os_mic%reallocate(n_old + nappended)
            endif
            j = n_old
            do ifile = 1,nfiles
                do imic = 1,parts(ifile)%os_mic%get_noris()
                    if( .not. selected(parts(ifile)%os_mic, imic) ) cycle
                    j = j + 1
                    call os_mic%transfer_ori(j, parts(ifile)%os_mic, imic)
                enddo
            enddo
        endif
        do ifile = 1,nfiles
            call parts(ifile)%kill
        enddo
        deallocate(parts)

    contains

        logical function selected( os, i )
            class(oris), intent(in) :: os
            integer,     intent(in) :: i
            selected = .true.
            if( l_accepted_only ) selected = os%get_state(i) > 0
        end function selected

    end subroutine append_mics_from_projects

end module simple_mic_import
