module simple_singletons
use simple_defs
!include 'simple_lib.f08'
!use simple_cmdline,          only: cmdline => cmdline_class
!use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_build,            only: b
use simple_params,           only: p
implicit none
private
!public :: pftcc, init_pftcc
public :: b, p
public :: init_params

!type(polarft_corrcalc) :: pftcc
!type(build) ::  b
!type(params_class) :: p
!type(cmdline_class) :: cline
contains

    ! subroutine init_pftcc( nrefs, ptcl_mask, eoarr )
    !     integer,                 intent(in)    :: nrefs
    !     logical, optional,       intent(in)    :: ptcl_mask(:)
    !     integer, optional,       intent(in)    :: eoarr(:)
    !     call pftcc%new(nrefs, ptcl_mask, eoarr)
    ! end subroutine init_pftcc

    subroutine init_params(cline, allow_mix, del_scaled, spproj_a_seg )
        use simple_cmdline,          only: cmdline
        use simple_params
        class(cmdline),     intent(inout) :: cline
        logical, optional, intent(in)    :: allow_mix
        logical, optional, intent(in)    :: del_scaled
        integer, optional, intent(in)    :: spproj_a_seg
        logical :: mix,delscaled
        integer :: seg
        mix=.false.
        delscaled=.false.
        seg=GENERIC_SEG
        if(p%singleton_initiated .eqv. .false.)then
            if(present(allow_mix)) mix = allow_mix
            if(present(del_scaled) ) delscaled = del_scaled
            if (present(spproj_a_seg))then
                call p%new(cline, mix, delscaled, spproj_a_seg)
            else
                call p%new(cline, mix, delscaled)

            endif
            p%singleton_initiated = .true.
        else
            print *, 'WARNING: attempting to re-initialize param singleton '

        endif
    end subroutine init_params

end module simple_singletons
