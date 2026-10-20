!@descr: canonical sigma2 store steps shared by the refine2D and refine3D strategies
module simple_strategy_sigma2
use, intrinsic :: iso_fortran_env, only: int64
use simple_defs,              only: STDLEN
use simple_error,             only: simple_exception
use simple_string,            only: string
use simple_parameters,        only: parameters
use simple_builder,           only: builder
use simple_sigma2_state,      only: sigma2_state_prepare_update
use simple_sigma2_state_file, only: sigma2_state_header, sigma2_state_read_header
implicit none

public :: prepare_canonical_sigma_update, read_sigma2_bench_counts
private
#include "simple_local_flags.inc"

contains

    !> before the parts write the next generation: the committed store validated, nothing left by another
    !! transaction of that generation
    subroutine prepare_canonical_sigma_update( params, build )
        type(parameters), intent(in)    :: params
        type(builder),    intent(inout) :: build
        type(string)          :: state_path
        integer               :: status
        logical               :: found
        character(len=STDLEN) :: message
        call build%spproj%get_sigma2_state_path(state_path, found)
        if( .not. found ) THROW_HARD('particle project has no canonical sigma2 state path')
        call sigma2_state_prepare_update(state_path%to_char(), params%nparts, params%numlen, status, message)
        if( status /= 0 ) THROW_HARD(trim(message))
        call state_path%kill
    end subroutine prepare_canonical_sigma_update

    !> the canonical store's size for a bench record: particles, shells, bytes and the rows updated this iteration;
    !! unchanged when the project has no store
    subroutine read_sigma2_bench_counts( build, nptcls, nshells, nbytes, nrows_updated )
        type(builder),  intent(inout) :: build
        integer,        intent(inout) :: nptcls, nshells, nrows_updated
        integer(int64), intent(inout) :: nbytes
        type(sigma2_state_header) :: header
        type(string)              :: state_path
        integer, allocatable      :: rep_rows(:), smp_rows(:), rec_rows(:)
        integer                   :: status
        logical                   :: found
        character(len=STDLEN)     :: message
        call build%spproj%get_sigma2_state_path(state_path, found)
        if( .not. found ) return
        call sigma2_state_read_header(state_path%to_char(), header, status, message)
        if( status /= 0 ) THROW_HARD(trim(message))
        nptcls  = header%nptcls
        nshells = header%kto - header%kfrom + 1
        nbytes  = header%file_bytes
        call build%spproj_field%get_update_rows(rep_rows, smp_rows, rec_rows)
        nrows_updated = size(smp_rows)
        deallocate(rep_rows, smp_rows, rec_rows)
        call state_path%kill
    end subroutine read_sigma2_bench_counts

end module simple_strategy_sigma2
