!@descr: flex_pca artifact catalog: where part files live and how they are named
!!
!! Every part file (probe/embed parts, state part volumes, PCG raw accumulations) lives under one
!! directory: the run directory by default, or a node-local scratch directory when the master
!! decided so (flex_pca_local_part_dir). Producers and consumers name parts only through
!! flex_pca_part_path / flex_pca_part_fname. The part-file contract shared by every codec is the
!! magic value here; each codec owns its own version, shape header and byte layout.
module simple_flex_pca_artifacts
use simple_core_module_api
use simple_parameters, only: parameters
implicit none
private

public :: flex_pca_part_fname, flex_pca_part_path
public :: flex_pca_set_part_dir, flex_pca_local_part_dir
public :: FLEX_PCA_PART_MAGIC

!> every part file is magic + version + shape header + payload, written to a .tmp and renamed,
!! so a master that finds the final name is guaranteed a complete file
integer, parameter :: FLEX_PCA_PART_MAGIC = 1180053590

!> part-directory prefix ('<dir>/'), unallocated = the run directory (see flex_pca_part_path)
character(len=:), allocatable :: part_dir_prefix

contains

    function flex_pca_part_fname( prefix, part, numlen ) result( fname )
        character(len=*), intent(in) :: prefix
        integer,          intent(in) :: part, numlen
        type(string) :: fname
        fname = flex_pca_part_path('flex_pca_'//prefix//'_part'//int2str_pad(part, numlen)//'.bin')
    end function flex_pca_part_fname

    function flex_pca_part_path( name ) result( path )
        character(len=*), intent(in) :: name
        type(string) :: path
        if( allocated(part_dir_prefix) )then
            path = string(part_dir_prefix)//name
        else
            path = string(name)
        endif
    end function flex_pca_part_path

    subroutine flex_pca_set_part_dir( dir )
        character(len=*), intent(in) :: dir
        if( len_trim(dir) == 0 )then
            if( allocated(part_dir_prefix) ) deallocate(part_dir_prefix)
        else
            part_dir_prefix = trim(dir)//'/'
        endif
    end subroutine flex_pca_set_part_dir

    !> Node-local part directory for the local queue system: parts are written and reduced every
    !! round, so they go to the disk the user already declared local through cache_dir.
    !! Master and workers derive the SAME name from the run directory (no new key travels), which
    !! is safe because local workers run on the master's node in the master's directory. Empty
    !! (= run directory) for any other queue system, when cache_dir is not given, or in shared
    !! memory.
    function flex_pca_local_part_dir( params, nparts ) result( dir )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: nparts
        character(len=:), allocatable :: dir
        type(string) :: cwd
        character(len=:), allocatable :: cwdc
        integer(kind=8) :: h
        integer :: i
        dir = ''
        if( nparts < 2 ) return
        if( trim(params%qsys_name) /= 'local' ) return
        if( params%cache_dir%is_blank() ) return
        call simple_getcwd(cwd)
        cwdc = cwd%to_char()
        h = 1469598103934665603_8            ! FNV-1a over the run directory path
        do i = 1, len_trim(cwdc)
            h = ieor(h, int(ichar(cwdc(i:i)),8))
            h = h * 1099511628211_8
        end do
        h = iand(h, 9223372036854775807_8)
        dir = params%cache_dir%to_char()//'/flex_pca_parts_'//trim(int2str(int(mod(h, 1000000007_8))))
        call cwd%kill
    end function flex_pca_local_part_dir

end module simple_flex_pca_artifacts
