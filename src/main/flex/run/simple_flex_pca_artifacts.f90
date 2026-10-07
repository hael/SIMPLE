!@descr: flex_pca artifact catalog: where part files live and how they are named
!!
!! Every part file (probe/embed parts, state part volumes, PCG raw accumulations) lives under one
!! directory: the run directory by default, or a node-local scratch directory when the master
!! decided so (flex_pca_local_part_dir). One catalog value is owned by the run's rounds object;
!! producers and consumers receive that object explicitly. The part-file contract shared by
!! every codec is the magic value here; each codec owns its own version, shape header and layout.
module simple_flex_pca_artifacts
use simple_core_module_api, only: int2str, int2str_pad, simple_exception, simple_getcwd, string
use simple_parameters, only: parameters
implicit none
private

public :: flex_pca_artifact_catalog, flex_pca_local_part_dir
public :: FLEX_PCA_PART_MAGIC
public :: MERGED_PC_FBODY, MERGED_META, MERGED_EIG_FNAME, PAIRED_MANIFEST

!> every part file is magic + version + shape header + payload, written to a .tmp and renamed,
!! so a master that finds the final name is guaranteed a complete file
integer, parameter :: FLEX_PCA_PART_MAGIC = 1180053590

!> Published products of the paired-fit merge.
character(len=*), parameter :: MERGED_PC_FBODY  = 'flex_pca_merged_pc'
character(len=*), parameter :: MERGED_META      = 'flex_pca_probe_merged.txt'
character(len=*), parameter :: MERGED_EIG_FNAME = 'flex_pca_eigenvalues_merged.txt'
character(len=*), parameter :: PAIRED_MANIFEST  = 'flex_pca_paired.txt'

!> Run-owned part-file namespace. An empty prefix means the run directory.
type :: flex_pca_artifact_catalog
    private
    character(len=:), allocatable :: part_dir_prefix
  contains
    procedure :: new        => artifact_catalog_new
    procedure :: kill       => artifact_catalog_kill
    procedure :: part_fname => artifact_part_fname
    procedure :: part_path  => artifact_part_path
end type flex_pca_artifact_catalog

contains

    function artifact_part_fname( self, prefix, part, numlen ) result( fname )
        class(flex_pca_artifact_catalog), intent(in) :: self
        character(len=*), intent(in) :: prefix
        integer,          intent(in) :: part, numlen
        type(string) :: fname
        fname = self%part_path('flex_pca_'//prefix//'_part'//int2str_pad(part, numlen)//'.bin')
    end function artifact_part_fname

    function artifact_part_path( self, name ) result( path )
        class(flex_pca_artifact_catalog), intent(in) :: self
        character(len=*), intent(in) :: name
        type(string) :: path
        if( allocated(self%part_dir_prefix) )then
            path = string(self%part_dir_prefix)//name
        else
            path = string(name)
        endif
    end function artifact_part_path

    subroutine artifact_catalog_new( self, dir )
        class(flex_pca_artifact_catalog), intent(inout) :: self
        character(len=*), intent(in) :: dir
        call self%kill
        if( len_trim(dir) == 0 )then
            return
        else
            self%part_dir_prefix = trim(dir)//'/'
        endif
    end subroutine artifact_catalog_new

    subroutine artifact_catalog_kill( self )
        class(flex_pca_artifact_catalog), intent(inout) :: self
        if( allocated(self%part_dir_prefix) ) deallocate(self%part_dir_prefix)
    end subroutine artifact_catalog_kill

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
