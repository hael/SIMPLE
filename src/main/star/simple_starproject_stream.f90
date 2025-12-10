module simple_starproject_stream
include 'simple_lib.f08'
!$ use omp_lib
use simple_sp_project, only: sp_project
use simple_cmdline,    only: cmdline
use simple_parameters, only: params_glob
use simple_histogram,  only: histogram
use simple_starfile_wrappers
use simple_math
use CPlot2D_wrapper_module
use simple_rnd
use FoX_dom
implicit none

public :: starproject_stream
private
#include "simple_local_flags.inc"

type starproject_stream
    type(starfile_table_type) :: starfile
    type(string)              :: starfile_name
    type(string)              :: starfile_tmp
    type(string)              :: rootpath
    real                      :: shift_threshold = 0.05
    integer                   :: group_offset = 0
    logical                   :: nicestream   = .false.
    logical                   :: use_beamtilt = .false.
    logical                   :: verbose      = .false.
contains
    ! export
    procedure          :: stream_export_micrographs
    procedure          :: stream_export_optics
    procedure          :: stream_export_particles_2D
    procedure          :: stream_export_pick_diameters
    procedure          :: stream_export_picking_references
    ! starfile
    procedure, private :: starfile_init
    procedure, private :: starfile_deinit
    procedure, private :: starfile_write_table
    procedure, private :: starfile_set_optics_table
    procedure, private :: starfile_set_optics_group_table
    procedure, private :: starfile_set_micrographs_table
    procedure, private :: starfile_set_particles2D_table
    procedure, private :: starfile_set_particles2D_subtable
    procedure, private :: starfile_set_pick_diameters_table
    procedure, private :: starfile_set_clusters2D_table
    ! optics
    procedure, private :: assign_optics_single
    procedure, private :: assign_optics
    procedure          :: copy_optics
    procedure          :: copy_micrographs_optics
end type starproject_stream

type :: starpart
    integer                       :: index
    integer                       :: nstart
    integer                       :: nend
    integer                       :: length
    type(starfile_table_type)     :: startable
    character(len=:), allocatable :: str
end type starpart

contains

    ! starfiles

    subroutine starfile_init( self, fname, outdir, verbose)
        class(starproject_stream), intent(inout) :: self
        class(string),             intent(in)    :: fname
        class(string),             intent(in)    :: outdir
        logical, optional,         intent(in)    :: verbose
        type(string) :: cwd, stem 
        if(present(verbose)) self%verbose = verbose
        self%starfile_name = fname
        self%starfile_tmp  = fname // '.tmp'
        call simple_getcwd(cwd)
        stem = basename(stemname(cwd))
        self%rootpath = stem
        self%nicestream = .true.
        call starfile_table__new(self%starfile)
        call cwd%kill
        call stem%kill
    end subroutine starfile_init

    subroutine starfile_deinit( self )
        class(starproject_stream), intent(inout) :: self
        call starfile_table__delete(self%starfile)
        if(file_exists(self%starfile_tmp)) then
            if(file_exists(self%starfile_name)) call del_file(self%starfile_name)
            call simple_rename(self%starfile_tmp, self%starfile_name)
        end if
    end subroutine starfile_deinit

    subroutine starfile_write_table( self, append )
        class(starproject_stream), intent(inout) :: self
        logical,                   intent(in)    :: append
        integer :: append_int
        append_int = 0
        if(append) append_int = 1
        call starfile_table__open_ofile(self%starfile, self%starfile_tmp%to_char(), append_int)
        call starfile_table__write_ofile(self%starfile)
        call starfile_table__close_ofile(self%starfile)
    end subroutine starfile_write_table

    subroutine starfile_set_optics_table( self, spproj )
        class(starproject_stream),  intent(inout) :: self
        class(sp_project),          intent(inout) :: spproj
        integer                    :: i
        character(len=XLONGSTRLEN) :: str_og
        call starfile_table__clear(self%starfile)
        call starfile_table__new(self%starfile)
        call starfile_table__setIsList(self%starfile, .false.)
        call starfile_table__setname(self%starfile, 'optics')
        do i=1,spproj%os_optics%get_noris()
            if(spproj%os_optics%get(i, 'state') .eq. 0.0 ) cycle
            call starfile_table__addObject(self%starfile)
            ! ints
            if(spproj%os_optics%isthere(i, 'ogid'))   call starfile_table__setValue_int(self%starfile, EMDL_IMAGE_OPTICS_GROUP, spproj%os_optics%get_int(i, 'ogid'))
            if(spproj%os_optics%isthere(i, 'pop'))    call starfile_table__setValue_int(self%starfile, SMPL_OPTICS_POPULATION,  spproj%os_optics%get_int(i, 'pop' ))
            ! doubles
            if(spproj%os_optics%isthere(i, 'kv'))     call starfile_table__setValue_double(self%starfile, EMDL_CTF_VOLTAGE,      real(spproj%os_optics%get(i, 'kv'),    dp))
            if(spproj%os_optics%isthere(i, 'smpd'))   call starfile_table__setValue_double(self%starfile, EMDL_IMAGE_PIXEL_SIZE, real(spproj%os_optics%get(i, 'smpd'),  dp))
            if(spproj%os_optics%isthere(i, 'cs'))     call starfile_table__setValue_double(self%starfile, EMDL_CTF_CS,           real(spproj%os_optics%get(i, 'cs'),    dp))
            if(spproj%os_optics%isthere(i, 'fraca'))  call starfile_table__setValue_double(self%starfile, EMDL_CTF_Q0,           real(spproj%os_optics%get(i, 'fraca'), dp))
            if(spproj%os_optics%isthere(i, 'opcx'))   call starfile_table__setValue_double(self%starfile, SMPL_OPTICS_CENTROIDX, real(spproj%os_optics%get(i, 'opcx'),  dp))
            if(spproj%os_optics%isthere(i, 'opcy'))   call starfile_table__setValue_double(self%starfile, SMPL_OPTICS_CENTROIDY, real(spproj%os_optics%get(i, 'opcy'),  dp))
            ! strings
            call spproj%os_optics%get_static(i, 'ogname', str_og)
            if(spproj%os_optics%isthere(i, 'ogname')) call starfile_table__setValue_string(self%starfile, EMDL_IMAGE_OPTICS_GROUP_NAME, trim(str_og))
        end do
    end subroutine starfile_set_optics_table

    subroutine starfile_set_optics_group_table( self, spproj, ogid )
        class(starproject_stream), intent(inout) :: self
        class(sp_project),         intent(inout) :: spproj
        integer,                   intent(in)    :: ogid
        integer :: i
        call starfile_table__clear(self%starfile)
        call starfile_table__new(self%starfile)
        call starfile_table__setIsList(self%starfile, .false.)
        call starfile_table__setname(self%starfile, 'opticsgroup_' // int2str(ogid))
        do i=1,spproj%os_mic%get_noris()
            if(spproj%os_mic%isthere(i, 'ogid') .and. spproj%os_mic%get_int(i, 'ogid') == ogid) then
                if(spproj%os_mic%get_state(i) .eq. 0 ) cycle
                call starfile_table__addObject(self%starfile)
                ! ints
                call starfile_table__setValue_int(self%starfile, EMDL_IMAGE_OPTICS_GROUP, spproj%os_mic%get_int(i, 'ogid'))
                ! doubles
                if(spproj%os_mic%isthere(i, 'shiftx')) call starfile_table__setValue_double(self%starfile, SMPL_OPTICS_SHIFTX, real(spproj%os_mic%get(i, 'shiftx'), dp))
                if(spproj%os_mic%isthere(i, 'shifty')) call starfile_table__setValue_double(self%starfile, SMPL_OPTICS_SHIFTY, real(spproj%os_mic%get(i, 'shifty'), dp))
            end if
        end do
    end subroutine starfile_set_optics_group_table

    subroutine starfile_set_pick_diameters_table( self, histogram_moldiams )
        class(starproject_stream), intent(inout) :: self
        type(histogram),           intent(in)    :: histogram_moldiams
        integer :: i
        call starfile_table__clear(self%starfile)
        call starfile_table__new(self%starfile)
        call starfile_table__setIsList(self%starfile, .false.)
        call starfile_table__setname(self%starfile, 'pick_diameters')
        do i=1, histogram_moldiams%get_nbins()
            call starfile_table__addObject(self%starfile)
            ! ints
            call starfile_table__setValue_int(self%starfile,    SMPL_PICK_POPULATION, int(histogram_moldiams%get(i)))
            ! doubles
            call starfile_table__setValue_double(self%starfile, SMPL_PICK_DIAMETER,   real(histogram_moldiams%get_x(i), dp))
        end do
    end subroutine starfile_set_pick_diameters_table

    subroutine starfile_set_micrographs_table( self, spproj )
        class(starproject_stream),  intent(inout) :: self
        class(sp_project),          intent(inout) :: spproj
        integer               :: i, pathtrim
        character(len=XLONGSTRLEN) :: str_mov, str_intg, str_mcs, str_boxf, str_ctfj
        pathtrim = 0
        call starfile_table__clear(self%starfile)
        call starfile_table__new(self%starfile)
        call starfile_table__setIsList(self%starfile, .false.)
        call starfile_table__setname(self%starfile, 'micrographs')
        do i=1,spproj%os_mic%get_noris()
            if(spproj%os_mic%get_state(i) .eq. 0 ) cycle
            call starfile_table__addObject(self%starfile)
            ! ints
            if(spproj%os_mic%isthere(i, 'ogid'   )) call starfile_table__setValue_int(self%starfile, EMDL_IMAGE_OPTICS_GROUP, spproj%os_mic%get_int(i, 'ogid'   ))
            if(spproj%os_mic%isthere(i, 'xdim'   )) call starfile_table__setValue_int(self%starfile, EMDL_IMAGE_SIZE_X,       spproj%os_mic%get_int(i, 'xdim'   ))
            if(spproj%os_mic%isthere(i, 'ydim'   )) call starfile_table__setValue_int(self%starfile, EMDL_IMAGE_SIZE_Y,       spproj%os_mic%get_int(i, 'ydim'   ))
            if(spproj%os_mic%isthere(i, 'nframes')) call starfile_table__setValue_int(self%starfile, EMDL_IMAGE_SIZE_Z,       spproj%os_mic%get_int(i, 'nframes'))
            if(spproj%os_mic%isthere(i, 'nptcls' )) call starfile_table__setValue_int(self%starfile, SMPL_N_PTCLS,            spproj%os_mic%get_int(i, 'nptcls' ))
            if(spproj%os_mic%isthere(i, 'nmics'  )) call starfile_table__setValue_int(self%starfile, SMPL_N_MICS,             spproj%os_mic%get_int(i, 'nmics'  ))
            if(spproj%os_mic%isthere(i, 'micid'  )) call starfile_table__setValue_int(self%starfile, SMPL_MIC_ID,             spproj%os_mic%get_int(i, 'micid'  ))
            ! doubles
            if(spproj%os_mic%isthere(i, 'dfx'    )) call starfile_table__setValue_double(self%starfile,  EMDL_CTF_DEFOCUSU,      real(spproj%os_mic%get(i, 'dfx') / 0.0001, dp))
            if(spproj%os_mic%isthere(i, 'dfy'    )) call starfile_table__setValue_double(self%starfile,  EMDL_CTF_DEFOCUSV,      real(spproj%os_mic%get(i, 'dfy') / 0.0001, dp))
            if(spproj%os_mic%isthere(i, 'angast' )) call starfile_table__setValue_double(self%starfile,  EMDL_CTF_DEFOCUS_ANGLE, real(spproj%os_mic%get(i, 'angast'),       dp))
            if(spproj%os_mic%isthere(i, 'phshift')) call starfile_table__setValue_double(self%starfile,  EMDL_CTF_PHASESHIFT,    real(spproj%os_mic%get(i, 'phshift'),      dp))
            if(spproj%os_mic%isthere(i, 'ctfres' )) call starfile_table__setValue_double(self%starfile,  EMDL_CTF_MAXRES,        real(spproj%os_mic%get(i, 'ctfres'),       dp))
            if(spproj%os_mic%isthere(i, 'icefrac')) call starfile_table__setValue_double(self%starfile,  SMPL_ICE_FRAC,          real(spproj%os_mic%get(i, 'icefrac'),      dp))
            if(spproj%os_mic%isthere(i, 'astig'  )) call starfile_table__setValue_double(self%starfile,  SMPL_ASTIGMATISM,       real(spproj%os_mic%get(i, 'astig'),        dp))
            ! strings
            call spproj%os_mic%get_static(i, 'movie',       str_mov)
            call spproj%os_mic%get_static(i, 'intg',        str_intg)
            call spproj%os_mic%get_static(i, 'mc_starfile', str_mcs)
            call spproj%os_mic%get_static(i, 'boxfile',     str_boxf)
            call spproj%os_mic%get_static(i, 'ctfjpg',      str_ctfj)
            str_intg = get_relative_path_here(str_intg)
            str_mcs  = get_relative_path_here(str_mcs)
            str_boxf = get_relative_path_here(str_boxf)
            str_ctfj = get_relative_path_here(str_ctfj)
            if(spproj%os_mic%isthere(i, 'movie'      )) call starfile_table__setValue_string(self%starfile, EMDL_MICROGRAPH_MOVIE_NAME,    trim(str_mov))
            if(spproj%os_mic%isthere(i, 'intg'       )) call starfile_table__setValue_string(self%starfile, EMDL_MICROGRAPH_NAME,          trim(str_intg))
            if(spproj%os_mic%isthere(i, 'mc_starfile')) call starfile_table__setValue_string(self%starfile, EMDL_MICROGRAPH_METADATA_NAME, trim(str_mcs))
            if(spproj%os_mic%isthere(i, 'boxfile'    )) call starfile_table__setValue_string(self%starfile, EMDL_MICROGRAPH_COORDINATES,   trim(str_boxf))
            if(spproj%os_mic%isthere(i, 'ctfjpg'     )) call starfile_table__setValue_string(self%starfile, EMDL_CTF_PSPEC,                trim(str_ctfj))
        end do

        contains

            function get_relative_path_here ( path ) result ( newpath )
                character(len=*), intent(in) :: path
                character(len=XLONGSTRLEN)   :: newpath
                if(pathtrim .eq. 0) pathtrim = index(path, self%rootpath%to_char()) 
                newpath = trim(path(pathtrim:))
            end function get_relative_path_here

    end subroutine starfile_set_micrographs_table

    subroutine starfile_set_particles2D_table( self, spproj )
        class(starproject_stream), intent(inout) :: self
        class(sp_project),         intent(inout) :: spproj
        type(string)               :: stkname
        character(len=XLONGSTRLEN) :: str_stk, str_mic
        integer      :: i, ind_in_stk, stkind, pathtrim, half_boxsize
        pathtrim = 0
        call starfile_table__clear(self%starfile)
        call starfile_table__new(self%starfile)
        call starfile_table__setIsList(self%starfile, .false.)
        call starfile_table__setname(self%starfile, 'particles')
        do i=1,spproj%os_ptcl2d%get_noris()
            if(spproj%os_ptcl2d%get(i, 'state') .eq. 0.0 ) cycle
            call starfile_table__addObject(self%starfile)
            call spproj%get_stkname_and_ind('ptcl2D', i, stkname, ind_in_stk)
            stkind = floor(spproj%os_ptcl2d%get(ind_in_stk, 'stkind'))
            half_boxsize = floor(spproj%os_stk%get(stkind, 'box') / 2.0)
            ! ints
            if(spproj%os_ptcl2d%isthere(i, 'ogid'   )) call starfile_table__setValue_int(self%starfile, EMDL_IMAGE_OPTICS_GROUP, spproj%os_ptcl2d%get_int(i, 'ogid'))
            if(spproj%os_ptcl2d%isthere(i, 'class'  )) call starfile_table__setValue_int(self%starfile, EMDL_PARTICLE_CLASS,     spproj%os_ptcl2d%get_class(i))
            if(spproj%os_ptcl2d%isthere(i, 'gid'    )) call starfile_table__setValue_int(self%starfile, EMDL_MLMODEL_GROUP_NO,   spproj%os_ptcl2d%get_int(i, 'gid'))
            ! doubles
            if(spproj%os_ptcl2d%isthere(i, 'dfx'    )) call starfile_table__setValue_double(self%starfile,  EMDL_CTF_DEFOCUSU,              real(spproj%os_ptcl2d%get(i, 'dfx') / 0.0001,        dp))
            if(spproj%os_ptcl2d%isthere(i, 'dfy'    )) call starfile_table__setValue_double(self%starfile,  EMDL_CTF_DEFOCUSV,              real(spproj%os_ptcl2d%get(i, 'dfy') / 0.0001,        dp))
            if(spproj%os_ptcl2d%isthere(i, 'angast' )) call starfile_table__setValue_double(self%starfile,  EMDL_CTF_DEFOCUS_ANGLE,         real(spproj%os_ptcl2d%get(i, 'angast'),              dp))
            if(spproj%os_ptcl2d%isthere(i, 'phshift')) call starfile_table__setValue_double(self%starfile,  EMDL_CTF_PHASESHIFT,            real(spproj%os_ptcl2d%get(i, 'phshift'),             dp))
            if(spproj%os_ptcl2d%isthere(i, 'e3'     )) call starfile_table__setValue_double(self%starfile,  EMDL_ORIENT_PSI,                real(spproj%os_ptcl2d%get(i, 'e3'),                  dp))
            if(spproj%os_ptcl2d%isthere(i, 'xpos'   )) call starfile_table__setValue_double(self%starfile,  EMDL_IMAGE_COORD_X,             real(spproj%os_ptcl2d%get(i, 'xpos') + half_boxsize, dp))
            if(spproj%os_ptcl2d%isthere(i, 'ypos'   )) call starfile_table__setValue_double(self%starfile,  EMDL_IMAGE_COORD_Y,             real(spproj%os_ptcl2d%get(i, 'ypos') + half_boxsize, dp))
            if(spproj%os_ptcl2d%isthere(i, 'x'      )) call starfile_table__setValue_double(self%starfile,  EMDL_ORIENT_ORIGIN_X_ANGSTROM,  real(spproj%os_ptcl2d%get(i, 'x'),                   dp))
            if(spproj%os_ptcl2d%isthere(i, 'y'      )) call starfile_table__setValue_double(self%starfile,  EMDL_ORIENT_ORIGIN_Y_ANGSTROM,  real(spproj%os_ptcl2d%get(i, 'y'),                   dp))
            ! strings
            if(trim(str_stk) .ne. '' .and. ind_in_stk .gt. 0) then
                call stkname%to_static(str_stk)
                str_stk = get_relative_path_here(str_stk)
                str_mic = get_relative_path_here(spproj%get_micname(i))
                call starfile_table__setValue_string(self%starfile, EMDL_IMAGE_NAME,      int2str(ind_in_stk) // '@' // trim(str_stk))
                call starfile_table__setValue_string(self%starfile, EMDL_MICROGRAPH_NAME, trim(str_mic))
            end if

        end do

        contains

            function get_relative_path_here ( path ) result ( newpath )
                character(len=*), intent(in) :: path
                character(len=XLONGSTRLEN)   :: newpath
                if(pathtrim .eq. 0) pathtrim = index(path, self%rootpath%to_char()) 
                newpath = trim(path(pathtrim:))
            end function get_relative_path_here

    end subroutine starfile_set_particles2D_table

    subroutine starfile_set_particles2D_subtable( self, spproj, part )
        class(starproject_stream),  intent(inout)   :: self
        class(sp_project),          intent(inout)   :: spproj
        class(starpart),            intent(inout)   :: part
        type(string)               :: stkname
        character(len=XLONGSTRLEN) :: str_stk, str_mic
        integer      :: i, ind_in_stk, stkind, pathtrim, half_boxsize
        pathtrim = 0
        call starfile_table__new(part%startable)
        call starfile_table__setIsList(part%startable, .false.)
        call starfile_table__setname(part%startable, 'particles')
        do i=part%nstart, part%nend
            if(spproj%os_ptcl2d%get(i, 'state') .eq. 0.0 ) cycle
            call starfile_table__addObject(part%startable)
            call spproj%get_stkname_and_ind('ptcl2D', i, stkname, ind_in_stk)
            stkind = floor(spproj%os_ptcl2d%get(ind_in_stk, 'stkind'))
            half_boxsize = floor(spproj%os_stk%get(stkind, 'box') / 2.0)
            ! ints
            if(spproj%os_ptcl2d%isthere(i, 'ogid'   )) call starfile_table__setValue_int(part%startable, EMDL_IMAGE_OPTICS_GROUP, spproj%os_ptcl2d%get_int(i, 'ogid'))
            if(spproj%os_ptcl2d%isthere(i, 'class'  )) call starfile_table__setValue_int(part%startable, EMDL_PARTICLE_CLASS,     spproj%os_ptcl2d%get_class(i))
            if(spproj%os_ptcl2d%isthere(i, 'gid'    )) call starfile_table__setValue_int(part%startable, EMDL_MLMODEL_GROUP_NO,   spproj%os_ptcl2d%get_int(i, 'gid'))
            ! doubles
            if(spproj%os_ptcl2d%isthere(i, 'dfx'    )) call starfile_table__setValue_double(part%startable,  EMDL_CTF_DEFOCUSU,              real(spproj%os_ptcl2d%get(i, 'dfx') / 0.0001,        dp))
            if(spproj%os_ptcl2d%isthere(i, 'dfy'    )) call starfile_table__setValue_double(part%startable,  EMDL_CTF_DEFOCUSV,              real(spproj%os_ptcl2d%get(i, 'dfy') / 0.0001,        dp))
            if(spproj%os_ptcl2d%isthere(i, 'angast' )) call starfile_table__setValue_double(part%startable,  EMDL_CTF_DEFOCUS_ANGLE,         real(spproj%os_ptcl2d%get(i, 'angast'),              dp))
            if(spproj%os_ptcl2d%isthere(i, 'phshift')) call starfile_table__setValue_double(part%startable,  EMDL_CTF_PHASESHIFT,            real(spproj%os_ptcl2d%get(i, 'phshift'),             dp))
            if(spproj%os_ptcl2d%isthere(i, 'e3'     )) call starfile_table__setValue_double(part%startable,  EMDL_ORIENT_PSI,                real(spproj%os_ptcl2d%get(i, 'e3'),                  dp))
            if(spproj%os_ptcl2d%isthere(i, 'xpos'   )) call starfile_table__setValue_double(part%startable,  EMDL_IMAGE_COORD_X,             real(spproj%os_ptcl2d%get(i, 'xpos') + half_boxsize, dp))
            if(spproj%os_ptcl2d%isthere(i, 'ypos'   )) call starfile_table__setValue_double(part%startable,  EMDL_IMAGE_COORD_Y,             real(spproj%os_ptcl2d%get(i, 'ypos') + half_boxsize, dp))
            if(spproj%os_ptcl2d%isthere(i, 'x'      )) call starfile_table__setValue_double(part%startable,  EMDL_ORIENT_ORIGIN_X_ANGSTROM,  real(spproj%os_ptcl2d%get(i, 'x'),                   dp))
            if(spproj%os_ptcl2d%isthere(i, 'y'      )) call starfile_table__setValue_double(part%startable,  EMDL_ORIENT_ORIGIN_Y_ANGSTROM,  real(spproj%os_ptcl2d%get(i, 'y'),                   dp))
            ! strings
            if( trim(str_stk) .ne. '' .and. ind_in_stk .gt. 0) then
                !$omp critical
                call stkname%to_static(str_stk)
                str_stk = get_relative_path_here(str_stk)
                str_mic = get_relative_path_here(spproj%get_micname(i))
                call starfile_table__setValue_string(part%startable, EMDL_IMAGE_NAME,      int2str(ind_in_stk) // '@' // trim(str_stk))
                call starfile_table__setValue_string(part%startable, EMDL_MICROGRAPH_NAME, trim(str_mic))

                ! if allocatable strings are used by the above subrouitnes, they need to be changed to static

                !$omp end critical
            end if

        end do

        contains

            function get_relative_path_here ( path ) result ( newpath )
                character(len=*), intent(in) :: path
                character(len=XLONGSTRLEN)   :: newpath
                if(pathtrim .eq. 0) pathtrim = index(path, self%rootpath%to_char()) 
                newpath = trim(path(pathtrim:))
            end function get_relative_path_here
        
    end subroutine starfile_set_particles2D_subtable

    subroutine starfile_set_clusters2D_table( self, spproj )
        class(starproject_stream), intent(inout) :: self
        class(sp_project),         intent(inout) :: spproj
        character(len=XLONGSTRLEN) :: str_stk, stkname
        integer                    :: i, stkind, pathtrim
        pathtrim = 0
        call starfile_table__clear(self%starfile)
        call starfile_table__new(self%starfile)
        call starfile_table__setIsList(self%starfile, .false.)
        call starfile_table__setname(self%starfile, 'clusters')
        do i=1,spproj%os_cls2D%get_noris()
            if(spproj%os_cls2D%get(i, 'state') .eq. 0.0 ) cycle
            call starfile_table__addObject(self%starfile)
            ! ints
            if(spproj%os_cls2D%isthere(i, 'ncls')) call starfile_table__setValue_int(self%starfile, SMPL_N_CLS, spproj%os_cls2D%get_int(i, 'ncls'))
            ! doubles
            if(spproj%os_cls2D%isthere(i, 'res'  )) call starfile_table__setValue_double(self%starfile, EMDL_MLMODEL_ESTIM_RESOL_REF, real(spproj%os_cls2D%get(i, 'res'),   dp))
            if(spproj%os_cls2D%isthere(i, 'pop'  )) call starfile_table__setValue_double(self%starfile, EMDL_MLMODEL_PDF_CLASS,       real(spproj%os_cls2D%get(i, 'pop') ,  dp))
            ! strings
            if(spproj%os_cls2D%isthere(i, 'stkind') .and. spproj%os_cls2D%isthere(i, 'stk')) then
                stkind  = floor(spproj%os_cls2D%get(i, 'stkind'))
                call spproj%os_cls2D%get_static(i, 'stk', stkname)
                str_stk = get_relative_path_here(stkname)
                call starfile_table__setValue_string(self%starfile, EMDL_MLMODEL_REF_IMAGE,  int2str(stkind) // '@' // trim(str_stk))
            end if
        end do

        contains

            function get_relative_path_here ( path ) result ( newpath )
                character(len=*), intent(in) :: path
                character(len=XLONGSTRLEN)   :: newpath
                if(pathtrim .eq. 0) pathtrim = index(path, self%rootpath%to_char()) 
                newpath = trim(path(pathtrim:))
            end function get_relative_path_here

    end subroutine starfile_set_clusters2D_table

    ! export

    subroutine stream_export_micrographs( self, spproj, outdir, optics_set, filename)
        class(starproject_stream), intent(inout) :: self
        class(sp_project),         intent(inout) :: spproj
        class(string),             intent(in)    :: outdir
        class(string), optional,   intent(in)    :: filename
        logical,       optional,   intent(in)    :: optics_set
        logical  :: l_optics_set
        if( spproj%os_mic%get_noris() == 0 ) return
        l_optics_set = .false.
        if(present(optics_set)) l_optics_set = optics_set
        if(.not. l_optics_set) call self%assign_optics_single(spproj)
        if(present(filename)) then
            call self%starfile_init(filename, outdir)
        else
            call self%starfile_init(string('micrographs.star'), outdir)
        endif
        call self%starfile_set_optics_table(spproj)
        call self%starfile_write_table(append = .false.)
        call self%starfile_set_micrographs_table(spproj)
        call self%starfile_write_table(append = .true.)
        call self%starfile_deinit()
    end subroutine stream_export_micrographs
 
    subroutine stream_export_optics( self, spproj, outdir )
        class(starproject_stream), intent(inout) :: self
        class(sp_project),         intent(inout) :: spproj
        class(string),             intent(in)    :: outdir
        integer :: ioptics
        if(params_glob%beamtilt .eq. 'yes') then
            self%use_beamtilt = .true.
        else
            self%use_beamtilt = .false.
        end if
        call self%assign_optics(spproj)
        call self%starfile_init(string('optics.star'), outdir)
        call self%starfile_set_optics_table(spproj)
        call self%starfile_write_table(append = .false.)
        do ioptics = 1, spproj%os_optics%get_noris()
            call self%starfile_set_optics_group_table(spproj, ioptics)
            call self%starfile_write_table(append = .true.)
        end do
        call self%starfile_deinit()
    end subroutine stream_export_optics

    subroutine stream_export_particles_2D( self, spproj, outdir, optics_set, filename, verbose)
        class(starproject_stream), intent(inout) :: self
        class(sp_project),         intent(inout) :: spproj
        class(string),             intent(in)    :: outdir
        class(string), optional,   intent(in)    :: filename
        logical,       optional,   intent(in)    :: optics_set, verbose
        type(starpart), allocatable :: starparts(:)
        type(starpart)          :: newpart
        logical                 :: l_optics_set, l_verbose
        integer                 :: i, nptcls, nbatches, fhandle, ok
        integer,      parameter :: BATCHSIZE=10000
        integer,      parameter :: NTHR=4
        integer(timer_int_kind) :: ms0
        real(timer_int_kind)    :: ms_complete
        if( spproj%os_ptcl2D%get_noris() == 0 ) return
        l_optics_set = .false.
        l_verbose    = .false.
        if(present(optics_set)) l_optics_set = optics_set
        if(present(verbose))    l_verbose    = verbose
        if(.not. l_optics_set) call self%assign_optics_single(spproj)
        if(present(filename)) then
            call self%starfile_init(filename, outdir, verbose=l_verbose)
        else
            call self%starfile_init(string('particles2D.star'), outdir, verbose=l_verbose)
        endif
        if(file_exists(self%starfile_tmp)) call del_file(self%starfile_tmp)
        if(self%verbose) ms0 = tic()
        call self%starfile_set_optics_table(spproj)
        call self%starfile_write_table(append = .false.)
        if(self%verbose) then
            ms_complete = toc(ms0)
            print *,'particle star optics section written in :', ms_complete; call flush(6)
        endif
        if(NTHR .le. 1) then
            if(self%verbose) ms0 = tic()
            call self%starfile_set_particles2D_table(spproj)
            call self%starfile_write_table(append = .true.)
            call self%starfile_deinit()
            if(self%verbose) then
                ms_complete = toc(ms0)
                print *,'particle star written in :', ms_complete, 'using single thread'; call flush(6)
            endif
        else
            nptcls = spproj%os_ptcl2d%get_noris()
            nbatches = ceiling(real(nptcls) / real(BATCHSIZE))
            allocate(starparts(nbatches))
            call omp_set_num_threads(NTHR)
            if(self%verbose) ms0 = tic()
            !$omp parallel do private(i, newpart) default(shared) proc_bind(close)
            do i=1, nbatches
                newpart%index  = i 
                newpart%nstart = 1 + ((i - 1) * BATCHSIZE)
                newpart%nend   = i * BATCHSIZE
                if(newpart%nend .gt. nptcls) newpart%nend = nptcls
                call self%starfile_set_particles2D_subtable(spproj, newpart)
                if(i .eq. 1) then
                    call starfile_table__write_omem(newpart%startable, newpart%str, newpart%length)
                else
                    call starfile_table__write_omem(newpart%startable, newpart%str, newpart%length, ignoreheader=.true.)
                endif
                call starfile_table__delete(newpart%startable)
                starparts(newpart%index) = newpart
                ! if( allocated(newpart%str) ) deallocate(newpart%str)

                ! no string allocations or deallocations in OpenMP sections
                ! if starfile_table__write_omem uses allocatable strings they need to be changed to static

            end do
            !$omp end parallel do
            if(self%verbose) then
                ms_complete = toc(ms0)
                print *,'particle star parts generated in :', ms_complete, 'using', NTHR, 'threads'; call flush(6)
            endif
            if(.not. file_exists(self%starfile_tmp)) THROW_HARD("stream_export_particles_2D: starfile headers not written")
            if(self%verbose) ms0 = tic()
            call fopen(fhandle, file=self%starfile_tmp, position='append', iostat=ok)
            do i=1,nbatches
                write(fhandle, '(a)', advance="no") starparts(i)%str
            end do
            !trailing empty line
            write(fhandle, '(a)')
            call fclose(fhandle)
            if(self%verbose) then
                ms_complete = toc(ms0)
                print *,'particle star parts written in :',  ms_complete; call flush(6)
            endif
            if(file_exists(self%starfile_name)) call del_file(self%starfile_name)
            call simple_rename(self%starfile_tmp, self%starfile_name)
        endif
        if(allocated(starparts)) deallocate(starparts)
    end subroutine stream_export_particles_2D

    subroutine stream_export_pick_diameters( self, outdir, histogram_moldiams, filename)
        class(starproject_stream), intent(inout) :: self
        type(histogram),           intent(inout) :: histogram_moldiams
        class(string),             intent(in)    :: outdir
        class(string), optional,   intent(in)    :: filename
        if(present(filename)) then
            call self%starfile_init(filename, outdir)
        else
            call self%starfile_init(string('pick.star'), outdir)
        endif
        call self%starfile_set_pick_diameters_table(histogram_moldiams)
        call self%starfile_write_table(append = .true.)
        call self%starfile_deinit()
    end subroutine stream_export_pick_diameters

    subroutine stream_export_picking_references( self, spproj, outdir, filename)
        class(starproject_stream), intent(inout) :: self
        class(sp_project),         intent(inout) :: spproj
        class(string),             intent(in)    :: outdir
        class(string),optional,    intent(in)    :: filename
        if(present(filename)) then
            call self%starfile_init(filename, outdir)
        else
            call self%starfile_init(string('pickrefs.star'), outdir)
        endif
        call self%starfile_set_clusters2D_table(spproj)
        call self%starfile_write_table(append = .true.)
        call self%starfile_deinit()
    end subroutine stream_export_picking_references

    ! optics

    subroutine assign_optics_single( self, spproj )
        class(starproject_stream),  intent(inout) :: self
        class(sp_project),          intent(inout) :: spproj
        integer                                   :: i
        call spproj%os_mic%set_all2single('ogid', 1.0)
        call spproj%os_stk%set_all2single('ogid', 1.0)
        call spproj%os_ptcl2D%set_all2single('ogid', 1.0)
        call spproj%os_ptcl3D%set_all2single('ogid', 1.0)
        call spproj%os_optics%new(1, is_ptcl=.false.)
        call spproj%os_optics%set(1, "ogid",   1.0)
        call spproj%os_optics%set(1, "ogname", "opticsgroup1")
        do i=1,spproj%os_mic%get_noris()
            if(spproj%os_mic%get(i, 'state') .gt. 0.0 ) exit      
        end do
        call spproj%os_optics%set(1, "smpd",  spproj%os_mic%get(i, "smpd")   )
        call spproj%os_optics%set(1, "cs",    spproj%os_mic%get(i, "cs")     )
        call spproj%os_optics%set(1, "kv",    spproj%os_mic%get(i, "kv")     )
        call spproj%os_optics%set(1, "fraca", spproj%os_mic%get(i, "fraca")  )
        call spproj%os_optics%set(1, "state", 1.0)
        call spproj%os_optics%set(1, "opcx",  0.0)
        call spproj%os_optics%set(1, "opcy",  0.0)
        call spproj%os_optics%set(1, "pop",   real(spproj%os_mic%get_noris()))
    end subroutine assign_optics_single

    subroutine assign_optics( self, spproj )
        class(starproject_stream),  intent(inout) :: self
        class(sp_project),          intent(inout) :: spproj
        type(ori)                                 :: template_ori
        real                                      :: grp_info(10000, 3) ! max 10000 optics groups ! 1: centroid x, 2: centroid y, 3: population
        integer                                   :: i, ntilt, nshift
        self%shift_threshold = params_glob%tilt_thres
        call assign_tiltgroups()
        call assign_shiftgroups()
        do i=1, spproj%os_mic%get_noris()
            if(spproj%os_mic%get(i, 'state') .gt. 0.0 ) exit      
        end do
        call template_ori%new(.false.)
        call template_ori%set("smpd",  spproj%os_mic%get(i, "smpd")   )
        call template_ori%set("cs",    spproj%os_mic%get(i, "cs")     )
        call template_ori%set("kv",    spproj%os_mic%get(i, "kv")     )
        call template_ori%set("fraca", spproj%os_mic%get(i, "fraca")  )
        call template_ori%set("state", 1.0)
        call template_ori%set("pop",   0.0)
        call template_ori%set("ogid",  0.0)
        call template_ori%set("opcx",  0.0)
        call template_ori%set("opcy",  0.0)
        call template_ori%set("ogname", "opticsgroup")
        call spproj%os_optics%new(nshift, is_ptcl=.false.)
        do i = 1, nshift
            call spproj%os_optics%append(i, template_ori)
            call spproj%os_optics%set(i, "ogid", real(i + self%group_offset))
            call spproj%os_optics%set(i, "pop",  grp_info(i, 3))
            call spproj%os_optics%set(i, "opcx", grp_info(i, 1))
            call spproj%os_optics%set(i, "opcy", grp_info(i, 2))
            call spproj%os_optics%set(i, "ogname", "opticsgroup" // int2str(i + self%group_offset))
        end do
        call template_ori%kill()
        call spproj%write()

        contains

            subroutine assign_tiltgroups()
                real, allocatable :: tiltuniq(:)
                integer           :: itilt, imic
                if(self%use_beamtilt) then
                    call elim_dup(spproj%os_mic%get_all('tiltgrp'), tiltuniq)
                    ntilt = size(tiltuniq)
                    do itilt=1, ntilt 
                        do imic=1, spproj%os_mic%get_noris()
                            if(spproj%os_mic%get(imic, 'tiltgrp') == tiltuniq(itilt)) call spproj%os_mic%set(imic, 'tmpgrp', real(itilt))
                        end do
                    end do
                else
                    ntilt = 1
                    call spproj%os_mic%set_all2single('tmpgrp', 1.0)
                end if
                if(allocated(tiltuniq)) deallocate(tiltuniq)
                write(logfhandle,'(A,I8)') '>>> # TILT GROUPS ASSIGNED : ', ntilt
            end subroutine

            subroutine assign_shiftgroups()
                real,    allocatable :: tiltgrps(:), shiftxs(:), shiftys(:)
                real,    allocatable :: shifts(:,:), centroids(:,:)
                integer, allocatable :: populations(:), labels(:), indices(:)
                integer              :: itilt, imic, ishift, tiltgrppop
                write(logfhandle,'(A,F8.2)') '>>> CLUSTERING TILT GROUPS USING SHIFTS AND THRESHOLD : ', self%shift_threshold
                tiltgrps = spproj%os_mic%get_all('tmpgrp')
                shiftxs  = spproj%os_mic%get_all('shiftx')
                shiftys  = spproj%os_mic%get_all('shifty')
                nshift   = 0
                do itilt = 1, ntilt
                    write(logfhandle,'(A,I8)') '      CLUSTERING TILT GROUP : ', itilt
                    tiltgrppop = count(tiltgrps == itilt)
                    allocate(shifts(tiltgrppop, 2))
                    allocate(labels(tiltgrppop   ))
                    allocate(indices(tiltgrppop  ))
                    ishift = 1
                    do imic = 1, spproj%os_mic%get_noris()
                        if(tiltgrps(imic) == itilt) then
                            shifts(ishift, 1) = shiftxs(imic)
                            shifts(ishift, 2) = shiftys(imic)
                            indices(ishift) = imic 
                            ishift = ishift + 1
                        end if
                    end do
                    call h_clust(shifts, self%shift_threshold, labels, centroids, populations)
                    do ishift = 1, size(labels)
                        call spproj%os_mic%set(indices(ishift), 'ogid', real(labels(ishift) + nshift + self%group_offset))
                    end do
                    do ishift = 1, size(populations)
                        grp_info(ishift + nshift, 1) = centroids(ishift, 1)
                        grp_info(ishift + nshift, 2) = centroids(ishift, 2)
                        grp_info(ishift + nshift, 3) = populations(ishift)
                    end do
                    nshift = nshift + size(populations)
                    deallocate(shifts, labels, indices)
                    write(logfhandle,'(A,I8)') '        # SHIFT GROUPS ASSIGNED : ', size(populations)
                end do
                call spproj%os_mic%delete_entry('tmpgrp')
                if(allocated(populations)) deallocate(populations)
                if(allocated(labels))      deallocate(labels)
                if(allocated(indices))     deallocate(indices)
                if(allocated(shifts))      deallocate(shifts)
                if(allocated(centroids))   deallocate(centroids)
                if(allocated(tiltgrps))    deallocate(tiltgrps)
                if(allocated(shiftxs))     deallocate(shiftxs)
                if(allocated(shiftys))     deallocate(shiftys)
            end subroutine assign_shiftgroups

            ! distance threshold based yerarchical clustering
            ! Source https://www.mathworks.com/help/stats/hierarchical-clustering.html#bq_679x-10
            subroutine h_clust(data_in, thresh, labels, centroids, populations)
                real,                 intent(in)  :: data_in(:,:)   ! input data, point coords
                real,                 intent(in)  :: thresh         ! threshold for class merging
                integer,              intent(out) :: labels(:)      ! labels of the elements in vec
                real,    allocatable, intent(out) :: centroids(:,:) ! centroids of the classes
                integer, allocatable, intent(out) :: populations(:) ! number of elements belonging to each class
                real,    allocatable :: mat(:,:)                    ! pariwise distances matrix
                logical, allocatable :: mask(:), outliers(:)
                integer :: N, i, j, cnt, ncls
                integer :: index(1), loc1(1), loc2(1)
                real    :: d
                if( size(data_in, dim = 2) .ne. 2 )then
                    THROW_HARD('Input data should be two dimensional!; h_clust')
                endif
                N = size(data_in, dim = 1) ! number of data points
                ! 1) calc all the couples of distances, using euclid dist
                allocate(mat(N,N), source = 0.)
                do i = 1, N-1
                    do j = i+1, N
                        mat(i,j) = sqrt((data_in(i,1)-data_in(j,1))**2 + (data_in(i,2)-data_in(j,2))**2) ! pariwise euclidean distance
                        mat(j,i) = mat(i,j)
                    enddo
                enddo
                ! 2) Generate binary clusters
                allocate(mask(N),     source = .true. )
                allocate(outliers(N), source = .false.)
                ncls = 0
                do i = 1, N
                    if( mask(i) )then ! if it's not already been clustered
                        mask(i) = .false.
                        ! find the index of the couple
                        d = minval(mat(i,:), mask)
                        index(:) = minloc(mat(i,:), mask)
                        ncls = ncls + 1
                        ! assign labels
                        labels(i) = ncls
                        if(d <= thresh) then ! if it's not an outlier (it has a couple)
                            labels(index(1)) = labels(i)
                            mask(index(1)) = .false. ! index(1) has already been clustered
                        else
                            outliers(i) = .true.
                        endif
                    endif
                enddo
                ! 3) Calculate centroids
                allocate(centroids(ncls,2), source = 0.)
                mask = .true. ! reset
                do i = 1, ncls
                    ! find the first member of the class
                    loc1(:) = minloc(abs(labels-i))
                    if(.not. outliers(loc1(1))) then
                        mask(loc1(1)) = .false.
                        ! find the second member of the class
                        loc2(:) = minloc(abs(labels-i), mask)
                        mask(loc2(1)) = .false.
                        centroids(i,1) = (data_in(loc1(1),1) + data_in(loc2(1),1))/2.
                        centroids(i,2) = (data_in(loc1(1),2) + data_in(loc2(1),2))/2.
                    else ! the class has just one element
                        loc1(:) = minloc(abs(labels-i))
                        centroids(i,1) = data_in(loc1(1),1)
                        centroids(i,2) = data_in(loc1(1),2)
                        mask(loc1(1)) = .false.
                    endif
                enddo
                mask  = .true. ! reset
                ! 4) merge classes
                do i = 1, ncls-1
                    do j = i+1, ncls
                        if(sqrt((centroids(i,1)-centroids(j,1))**2+(centroids(i,2)-centroids(j,2))**2) <= thresh) then ! merge classes
                            ! change label to class j
                            where(labels == j) labels = i
                        endif
                    enddo
                enddo
                ! 5) Reoder labels
                cnt = 0
                do i = 1, ncls
                    if( any(labels== i) )then !there is a class labelled i
                        cnt = cnt + 1
                        where(labels == i) labels = cnt
                    endif
                enddo
                ! 6) recalculate centroids
                deallocate(centroids)
                ncls = maxval(labels) ! the nr of classes is maxval(labels)
                allocate(centroids(ncls,2), source = 0.)
                allocate(populations(ncls), source = 0 )
                mask = .true. ! reset
                do i = 1, ncls
                    populations(i) = count(labels == i) ! counts the nr of elements in the class
                    ! find the all the cnt member of the class and update the centroids
                    do j = 1, populations(i)
                        loc1(:) = minloc(abs(labels-i), mask)
                        mask(loc1(1)) = .false. ! do not consider this member of the class anymore
                        centroids(i,1) = centroids(i,1)+ data_in(loc1(1),1)
                        centroids(i,2) = centroids(i,2)+ data_in(loc1(1),2)
                    enddo
                    centroids(i,:) = centroids(i,:)/real(populations(i))
                enddo
            end subroutine h_clust

    end subroutine assign_optics

    subroutine copy_optics( self, spproj, spproj_src )
        class(starproject_stream),  intent(inout) :: self
        class(sp_project),          intent(inout) :: spproj, spproj_src
        integer, allocatable :: ogmap(:)
        real    :: min_importind, max_importind
        integer :: i, stkind
        call spproj%os_mic%minmax('importind', min_importind, max_importind)
        call spproj%os_optics%copy(spproj_src%os_optics, is_ptcl=.false.)
        allocate(ogmap(int(max_importind)))
        do i=1, int(max_importind)
            ogmap(i) = 1
        end do
        do i=1, spproj_src%os_mic%get_noris()
            if(spproj_src%os_mic%isthere(i, 'importind') .and. spproj_src%os_mic%isthere(i, 'ogid') .and. .not. spproj_src%os_mic%get_int(i, 'importind') .gt. max_importind) then
               ogmap(spproj_src%os_mic%get_int(i, 'importind')) = spproj_src%os_mic%get_int(i, 'ogid')
            end if
        end do
        do i=1, spproj%os_mic%get_noris()
            if(spproj%os_mic%isthere(i, 'importind')) then
                call spproj%os_mic%set(i, 'ogid', ogmap(spproj%os_mic%get_int(i, 'importind')))
            end if
        end do
        if(spproj%os_ptcl2d%get_noris() .gt. 0) then
            do i=1, spproj%os_ptcl2d%get_noris()
                if(spproj%os_ptcl2d%isthere(i, 'stkind')) then
                    stkind = spproj%os_ptcl2d%get_int(i, 'stkind')
                    call spproj%os_ptcl2d%set(i, 'ogid', spproj%os_mic%get_int(stkind, 'ogid'))
                end if
            end do
        end if
        deallocate(ogmap)
    end subroutine copy_optics  

    subroutine copy_micrographs_optics( self, spproj_dest, write, verbose )
        class(starproject_stream), intent(inout) :: self
        class(sp_project),         intent(inout) :: spproj_dest
        logical,         optional, intent(in)    :: write, verbose
        type(sp_project)        :: spproj_optics
        integer(timer_int_kind) :: ms0
        real(timer_int_kind)    :: ms_copy_optics
        logical                 :: l_verbose, l_write
        l_verbose = .false.
        l_write   = .false.
        if( present(verbose) ) l_verbose = verbose
        if( present(write)   ) l_write   = write
        if( (params_glob%projfile_optics .ne. '') .and.&
           &(file_exists(string('../')//params_glob%projfile_optics)) ) then
            if( l_verbose ) ms0 = tic()
            call spproj_optics%read(string('../')//params_glob%projfile_optics)
            call self%copy_optics(spproj_dest, spproj_optics)
            call spproj_optics%kill()
            if( l_write ) call spproj_dest%write
            if( l_verbose )then
                ms_copy_optics = toc(ms0)
                print *,'ms_copy_optics  : ', ms_copy_optics; call flush(6)
            endif
        end if
    end subroutine copy_micrographs_optics

end module simple_starproject_stream
