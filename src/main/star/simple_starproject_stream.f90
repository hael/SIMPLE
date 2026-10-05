!@descr: for star project management in streaming applications
module simple_starproject_stream
use simple_core_module_api
use simple_sp_project, only: sp_project
use simple_cmdline,    only: cmdline
use simple_parameters, only: parameters
use simple_starfile_wrappers
use CPlot2D_wrapper_module
implicit none

public :: starproject_stream
private
#include "simple_local_flags.inc"

type starproject_stream
    type(string)               :: projfile_optics ! params%projfile_optics of the last export
    type(starfile_table_type)  :: starfile
    type(string)               :: starfile_name
    type(string)               :: starfile_tmp
    type(string)               :: rootpath
    logical                    :: nicestream   = .false.
    logical                    :: verbose      = .false.
contains
    ! export
    procedure          :: stream_export_micrographs
    procedure          :: stream_write_optics
    procedure          :: stream_export_particles_2D
    ! starfile
    procedure, private :: starfile_init
    procedure, private :: starfile_deinit
    procedure, private :: starfile_write_table
    procedure, private :: starfile_set_optics_table
    procedure, private :: starfile_set_optics_group_table
    procedure, private :: starfile_set_micrographs_table
    procedure, private :: starfile_set_particles2D_table
    ! optics
    procedure, private :: assign_optics_single
    procedure          :: copy_optics
    procedure          :: copy_micrographs_optics
end type starproject_stream

contains

    ! starfiles

    ! The STAR file @p fname in @p outdir (the working directory when empty; an absolute @p fname
    ! is taken as it is). Paths in it are relative to the folder above @p outdir, the project root.
    subroutine starfile_init( self, params, fname, outdir, verbose)
        class(starproject_stream), intent(inout) :: self
        class(parameters),         intent(in)    :: params
        class(string),             intent(in)    :: fname
        class(string),             intent(in)    :: outdir
        logical, optional,         intent(in)    :: verbose
        type(string) :: cwd, stem
        self%projfile_optics = params%projfile_optics
        if(present(verbose)) self%verbose = verbose
        if( outdir%strlen() > 0 )then
            cwd = simple_abspath(outdir, check_exists=.false.)
        else
            call simple_getcwd(cwd)
        endif
        if( fname%to_char([1,1]) == '/' )then
            self%starfile_name = fname
        else
            self%starfile_name = cwd//'/'//fname
        endif
        self%starfile_tmp  = self%starfile_name // '.tmp'
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
            call starfile_table__setValue_double(self%starfile, EMDL_CTF_PHASESHIFT, &
                &real(rad2deg(spproj%os_mic%get(i, 'phshift')), dp))
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
                if( pathtrim > 0 ) then
                    newpath = trim(path(pathtrim:))
                else
                    newpath = trim(path)
                end if
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
            stkind       = spproj%os_ptcl2d%get_int(i, 'stkind')
            half_boxsize = spproj%os_stk%get_int(stkind, 'box') / 2
            ! ints
            if(spproj%os_ptcl2d%isthere(i, 'ogid'   )) call starfile_table__setValue_int(self%starfile, EMDL_IMAGE_OPTICS_GROUP, spproj%os_ptcl2d%get_int(i, 'ogid'))
            if(spproj%os_ptcl2d%isthere(i, 'class'  )) call starfile_table__setValue_int(self%starfile, EMDL_PARTICLE_CLASS,     spproj%os_ptcl2d%get_class(i))
            if(spproj%os_ptcl2d%isthere(i, 'gid'    )) call starfile_table__setValue_int(self%starfile, EMDL_MLMODEL_GROUP_NO,   spproj%os_ptcl2d%get_int(i, 'gid'))
            ! doubles
            if(spproj%os_ptcl2d%isthere(i, 'dfx'    )) call starfile_table__setValue_double(self%starfile,  EMDL_CTF_DEFOCUSU,              real(spproj%os_ptcl2d%get(i, 'dfx') / 0.0001,        dp))
            if(spproj%os_ptcl2d%isthere(i, 'dfy'    )) call starfile_table__setValue_double(self%starfile,  EMDL_CTF_DEFOCUSV,              real(spproj%os_ptcl2d%get(i, 'dfy') / 0.0001,        dp))
            if(spproj%os_ptcl2d%isthere(i, 'angast' )) call starfile_table__setValue_double(self%starfile,  EMDL_CTF_DEFOCUS_ANGLE,         real(spproj%os_ptcl2d%get(i, 'angast'),              dp))
            call starfile_table__setValue_double(self%starfile, EMDL_CTF_PHASESHIFT, &
                &real(rad2deg(spproj%os_ptcl2d%get(i, 'phshift')), dp))
            if(spproj%os_ptcl2d%isthere(i, 'e3'     )) call starfile_table__setValue_double(self%starfile,  EMDL_ORIENT_PSI,                real(spproj%os_ptcl2d%get(i, 'e3'),                  dp))
            if(spproj%os_ptcl2d%isthere(i, 'xpos'   )) call starfile_table__setValue_double(self%starfile,  EMDL_IMAGE_COORD_X,             real(spproj%os_ptcl2d%get(i, 'xpos') + half_boxsize, dp))
            if(spproj%os_ptcl2d%isthere(i, 'ypos'   )) call starfile_table__setValue_double(self%starfile,  EMDL_IMAGE_COORD_Y,             real(spproj%os_ptcl2d%get(i, 'ypos') + half_boxsize, dp))
            if(spproj%os_ptcl2d%isthere(i, 'x'      )) call starfile_table__setValue_double(self%starfile,  EMDL_ORIENT_ORIGIN_X_ANGSTROM,  real(spproj%os_ptcl2d%get(i, 'x'),                   dp))
            if(spproj%os_ptcl2d%isthere(i, 'y'      )) call starfile_table__setValue_double(self%starfile,  EMDL_ORIENT_ORIGIN_Y_ANGSTROM,  real(spproj%os_ptcl2d%get(i, 'y'),                   dp))
            ! strings: the particle's image in its stack, and its micrograph
            call spproj%get_stkname_and_ind('ptcl2D', i, stkname, ind_in_stk)
            if( stkname%strlen_trim() > 0 .and. ind_in_stk > 0 )then
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
                if( pathtrim > 0 ) then
                    newpath = trim(path(pathtrim:))
                else
                    newpath = trim(path)
                end if
            end function get_relative_path_here

    end subroutine starfile_set_particles2D_table

    ! export

    subroutine stream_export_micrographs( self, params, spproj, outdir, optics_set, filename)
        class(starproject_stream), intent(inout) :: self
        class(parameters),         intent(in)    :: params
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
            call self%starfile_init(params, filename, outdir)
        else
            call self%starfile_init(params, string('micrographs.star'), outdir)
        endif
        call self%starfile_set_optics_table(spproj)
        call self%starfile_write_table(append = .false.)
        call self%starfile_set_micrographs_table(spproj)
        call self%starfile_write_table(append = .true.)
        call self%starfile_deinit()
    end subroutine stream_export_micrographs

    ! Writes optics.star from the optics groups already in spproj; assigns nothing.
    subroutine stream_write_optics( self, params, spproj, outdir )
        class(starproject_stream), intent(inout) :: self
        class(parameters),         intent(in)    :: params
        class(sp_project),         intent(inout) :: spproj
        class(string),             intent(in)    :: outdir
        integer :: ioptics
        call self%starfile_init(params, string('optics.star'), outdir)
        call self%starfile_set_optics_table(spproj)
        call self%starfile_write_table(append = .false.)
        do ioptics = 1, spproj%os_optics%get_noris()
            call self%starfile_set_optics_group_table(spproj, ioptics)
            call self%starfile_write_table(append = .true.)
        end do
        call self%starfile_deinit()
    end subroutine stream_write_optics

    subroutine stream_export_particles_2D( self, params, spproj, outdir, optics_set, filename, verbose)
        class(starproject_stream), intent(inout) :: self
        class(parameters),         intent(in)    :: params
        class(sp_project),         intent(inout) :: spproj
        class(string),             intent(in)    :: outdir
        class(string), optional,   intent(in)    :: filename
        logical,       optional,   intent(in)    :: optics_set, verbose
        logical                 :: l_optics_set, l_verbose
        integer(timer_int_kind) :: ms0
        real(timer_int_kind)    :: ms_complete
        if( spproj%os_ptcl2D%get_noris() == 0 ) return
        l_optics_set = .false.
        l_verbose    = .false.
        if(present(optics_set)) l_optics_set = optics_set
        if(present(verbose))    l_verbose    = verbose
        if(.not. l_optics_set) call self%assign_optics_single(spproj)
        if(present(filename)) then
            call self%starfile_init(params, filename, outdir, verbose=l_verbose)
        else
            call self%starfile_init(params, string('particles2D.star'), outdir, verbose=l_verbose)
        endif
        if(file_exists(self%starfile_tmp)) call del_file(self%starfile_tmp)
        if(self%verbose) ms0 = tic()
        call self%starfile_set_optics_table(spproj)
        call self%starfile_write_table(append = .false.)
        if(self%verbose) then
            ms_complete = toc(ms0)
            print *,'particle star optics section written in :', ms_complete; call flush(6)
        endif
        if(self%verbose) ms0 = tic()
        call self%starfile_set_particles2D_table(spproj)
        call self%starfile_write_table(append = .true.)
        call self%starfile_deinit()
        if(self%verbose) then
            ms_complete = toc(ms0)
            print *,'particle star written in :', ms_complete, 'using single thread'; call flush(6)
        endif
    end subroutine stream_export_particles_2D

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
        if( (self%projfile_optics .ne. '') .and.&
           &(file_exists(string('../')//self%projfile_optics)) ) then
            if( l_verbose ) ms0 = tic()
            call spproj_optics%read(string('../')//self%projfile_optics)
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
