module simple_starproject
include 'simple_lib.f08'
!$ use omp_lib
use simple_sp_project, only: sp_project
use simple_cmdline,    only: cmdline
use simple_parameters, only: params_glob
use simple_starproject_utils
use CPlot2D_wrapper_module
use simple_rnd
use FoX_dom
implicit none

public :: starproject
private
#include "simple_local_flags.inc"

type starproject
    type(star_file)              :: starfile
    type(tilt_info), allocatable :: tiltinfo(:)
    logical                      :: automode = .false.
    !logical                      :: VERBOSE_OUTPUT =.false.
contains
    ! constructor
    procedure, private :: initialise
    ! import
    procedure          :: import_mics
    procedure          :: import_ptcls2D
    procedure          :: import_cls2D
    procedure          :: import_ptcls3D
    procedure, private :: import_stardata
    procedure, private :: populate_stkmap
    procedure, private :: read_starheaders
    procedure          :: check_stk_params
    ! export
    procedure          :: export_mics
    procedure          :: export_cls2D
    procedure          :: export_iter3D
    procedure          :: export_ptcls2D
    procedure          :: export_ptcls3D
    procedure, private :: export_stardata
    ! tilt
    procedure, private :: assign_initial_tiltgroups
    procedure, private :: assign_xml_tiltinfo
    procedure, private :: cluster_tiltinfo
    ! optics
    procedure          :: assign_optics
    procedure, private :: populate_opticsmap
    procedure, private :: plot_opticsgroups
    procedure, private :: sort_optics_maxpop
    procedure, private :: apply_optics_offset
    procedure, private :: get_image_basename
    !other
    procedure          :: set_verbose
    procedure          :: kill
end type starproject

contains

    ! constructor

    subroutine initialise( self )
        class(starproject),  intent(inout) :: self

        ! assign optics flags
        if (.not. allocated(self%starfile%optics%flags)) allocate(self%starfile%optics%flags(0))
        self%starfile%optics%flags = [self%starfile%optics%flags, star_flag(rlnflag="rlnVoltage", splflag="kv")]
        self%starfile%optics%flags = [self%starfile%optics%flags, star_flag(rlnflag="rlnOpticsGroup", splflag="ogid", int=.true.)]
        self%starfile%optics%flags = [self%starfile%optics%flags, star_flag(rlnflag="rlnOpticsGroupName", splflag="ogname", string=.true.)]
        self%starfile%optics%flags = [self%starfile%optics%flags, star_flag(rlnflag="rlnAmplitudeContrast", splflag="fraca")]
        self%starfile%optics%flags = [self%starfile%optics%flags, star_flag(rlnflag="rlnSphericalAberration", splflag="cs")]
      !  self%starfile%optics%flags = [self%starfile%optics%flags, star_flag(rlnflag="rlnMicrographPixelSize", splflag="smpd")]
        self%starfile%optics%flags = [self%starfile%optics%flags, star_flag(rlnflag="rlnImagePixelSize", splflag="smpd")]
        self%starfile%optics%flags = [self%starfile%optics%flags, star_flag(rlnflag="rlnImageSize", splflag="box", int=.true.)]

        ! assign micrographs flags
        if (.not. allocated(self%starfile%micrographs%flags)) allocate(self%starfile%micrographs%flags(0))
        self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnMicrographMovieName", splflag="movie", string=.true.)]
        self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnMicrographName", splflag="intg", string=.true.)]
        self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnMicrographMetadata", splflag="mc_starfile", string=.true.)]
        self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnDefocusU", splflag="dfx", mult=0.0001)]
        self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnDefocusV", splflag="dfy", mult=0.0001)]
        self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnDefocusAngle", splflag="angast")]
        self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnPhaseShift", splflag="phshift")]
        self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnCtfMaxResolution", splflag="ctfres")]
        self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="splIceFrac", splflag="icefrac")]
        self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnOpticsGroup", splflag="ogid", int=.true.)]
        self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnImageSizeX", splflag="xdim", int=.true.)]
        self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnImageSizeY", splflag="ydim", int=.true.)]
        self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnImageSizeZ", splflag="nframes", int=.true.)]
        self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnCtfPowerSpectrum", splflag="ctfjpg", string=.true.)]
        self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnMicrographCoordinates", splflag="boxfile", string=.true.)]
        self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="splNumberParticles", splflag="nptcls", int=.true.)]

        ! assign stk flags
        if (.not. allocated(self%starfile%stacks%flags)) allocate(self%starfile%stacks%flags(0))
        self%starfile%stacks%flags = [self%starfile%stacks%flags, star_flag(rlnflag="rlnImageName", splflag="stk", string=.true., imagesplit=.true., splflag2="stkind")]
        self%starfile%stacks%flags = [self%starfile%stacks%flags, star_flag(rlnflag="rlnOpticsGroup", splflag="ogid", int=.true.)]

        ! assign particles 2D flags
        if (.not. allocated(self%starfile%particles2D%flags)) allocate(self%starfile%particles2D%flags(0))
        self%starfile%particles2D%flags = [self%starfile%particles2D%flags, star_flag(rlnflag="rlnDefocusU", splflag="dfx", mult=0.0001)]
        self%starfile%particles2D%flags = [self%starfile%particles2D%flags, star_flag(rlnflag="rlnDefocusV", splflag="dfy", mult=0.0001)]
        self%starfile%particles2D%flags = [self%starfile%particles2D%flags, star_flag(rlnflag="rlnDefocusAngle", splflag="angast")]
        self%starfile%particles2D%flags = [self%starfile%particles2D%flags, star_flag(rlnflag="rlnAnglePsi", splflag="e3")]
        self%starfile%particles2D%flags = [self%starfile%particles2D%flags, star_flag(rlnflag="rlnCoordinateX", splflag="xpos")]
        self%starfile%particles2D%flags = [self%starfile%particles2D%flags, star_flag(rlnflag="rlnCoordinateY", splflag="ypos")]
        self%starfile%particles2D%flags = [self%starfile%particles2D%flags, star_flag(rlnflag="rlnOriginXAngst", splflag="x")]
        self%starfile%particles2D%flags = [self%starfile%particles2D%flags, star_flag(rlnflag="rlnOriginYAngst", splflag="y")]
        self%starfile%particles2D%flags = [self%starfile%particles2D%flags, star_flag(rlnflag="rlnPhaseShift", splflag="phshift")]
        self%starfile%particles2D%flags = [self%starfile%particles2D%flags, star_flag(rlnflag="rlnOpticsGroup", splflag="ogid", int=.true.)]
        self%starfile%particles2D%flags = [self%starfile%particles2D%flags, star_flag(rlnflag="rlnClassNumber", splflag="class", int=.true.)]
        self%starfile%particles2D%flags = [self%starfile%particles2D%flags, star_flag(rlnflag="rlnGroupNumber", splflag="gid", int=.true.)]
        self%starfile%particles2D%flags = [self%starfile%particles2D%flags, star_flag(rlnflag="rlnVoltage", splflag="kv")]
        self%starfile%particles2D%flags = [self%starfile%particles2D%flags, star_flag(rlnflag="rlnAmplitudeContrast", splflag="fraca")]
        self%starfile%particles2D%flags = [self%starfile%particles2D%flags, star_flag(rlnflag="rlnSphericalAberration", splflag="cs")]
        self%starfile%particles2D%flags = [self%starfile%particles2D%flags, star_flag(rlnflag="rlnImagePixelSize", splflag="smpd")]
        self%starfile%particles2D%flags = [self%starfile%particles2D%flags, star_flag(rlnflag="rlnDetectorPixelSize", splflag="detpix")]
        self%starfile%particles2D%flags = [self%starfile%particles2D%flags, star_flag(rlnflag="rlnMagnification", splflag="mag")]
        
        ! assign particles 3D flags
        if (.not. allocated(self%starfile%particles3D%flags)) allocate(self%starfile%particles3D%flags(0))
        self%starfile%particles3D%flags = [self%starfile%particles3D%flags, star_flag(rlnflag="rlnDefocusU", splflag="dfx", mult=0.0001)]
        self%starfile%particles3D%flags = [self%starfile%particles3D%flags, star_flag(rlnflag="rlnDefocusV", splflag="dfy", mult=0.0001)]
        self%starfile%particles3D%flags = [self%starfile%particles3D%flags, star_flag(rlnflag="rlnDefocusAngle", splflag="angast")]
        self%starfile%particles3D%flags = [self%starfile%particles3D%flags, star_flag(rlnflag="rlnAnglePsi", splflag="e3")]
        self%starfile%particles3D%flags = [self%starfile%particles3D%flags, star_flag(rlnflag="rlnAngleRot", splflag="e1")]
        self%starfile%particles3D%flags = [self%starfile%particles3D%flags, star_flag(rlnflag="rlnAngleTilt", splflag="e2")]
        self%starfile%particles3D%flags = [self%starfile%particles3D%flags, star_flag(rlnflag="rlnCoordinateX", splflag="xpos")]
        self%starfile%particles3D%flags = [self%starfile%particles3D%flags, star_flag(rlnflag="rlnCoordinateY", splflag="ypos")]
        self%starfile%particles3D%flags = [self%starfile%particles3D%flags, star_flag(rlnflag="rlnOriginXAngst", splflag="x")]
        self%starfile%particles3D%flags = [self%starfile%particles3D%flags, star_flag(rlnflag="rlnOriginYAngst", splflag="y")]
        self%starfile%particles3D%flags = [self%starfile%particles3D%flags, star_flag(rlnflag="rlnPhaseShift", splflag="phshift")]
        self%starfile%particles3D%flags = [self%starfile%particles3D%flags, star_flag(rlnflag="rlnOpticsGroup", splflag="ogid", int=.true.)]
        self%starfile%particles3D%flags = [self%starfile%particles3D%flags, star_flag(rlnflag="rlnClassNumber", splflag="class", int=.true.)]
        self%starfile%particles3D%flags = [self%starfile%particles3D%flags, star_flag(rlnflag="rlnGroupNumber", splflag="gid", int=.true.)]
        self%starfile%particles3D%flags = [self%starfile%particles3D%flags, star_flag(rlnflag="rlnVoltage", splflag="kv")]
        self%starfile%particles3D%flags = [self%starfile%particles3D%flags, star_flag(rlnflag="rlnAmplitudeContrast", splflag="fraca")]
        self%starfile%particles3D%flags = [self%starfile%particles3D%flags, star_flag(rlnflag="rlnSphericalAberration", splflag="cs")]
        self%starfile%particles3D%flags = [self%starfile%particles3D%flags, star_flag(rlnflag="rlnImagePixelSize", splflag="smpd")]
        self%starfile%particles3D%flags = [self%starfile%particles3D%flags, star_flag(rlnflag="rlnDetectorPixelSize", splflag="detpix")]
        self%starfile%particles3D%flags = [self%starfile%particles3D%flags, star_flag(rlnflag="rlnMagnification", splflag="mag")]
        
        ! assign clusters2D flags
        if (.not. allocated(self%starfile%clusters2D%flags)) allocate(self%starfile%clusters2D%flags(0))
        self%starfile%clusters2D%flags = [self%starfile%clusters2D%flags, star_flag(rlnflag="rlnReferenceImage", splflag="stk", string=.true., imagesplit=.true., splflag2="class")]
        self%starfile%clusters2D%flags = [self%starfile%clusters2D%flags, star_flag(rlnflag="rlnEstimatedResolution", splflag="res")]
      !  self%starfile%clusters2D%flags = [self%starfile%clusters2D%flags, star_flag(rlnflag="rlnClassDistribution", splflag="pop", int=.true.)]
        self%starfile%clusters2D%flags = [self%starfile%clusters2D%flags, star_flag(rlnflag="rlnClassDistribution", splflag="pop")]
        
        ! assign class3D flags
        if (.not. allocated(self%starfile%class3D%flags)) allocate(self%starfile%class3D%flags(0))
        self%starfile%class3D%flags = [self%starfile%class3D%flags, star_flag(rlnflag="rlnReferenceImage",      splflag="lpvol",    string=.true.)]
        self%starfile%class3D%flags = [self%starfile%class3D%flags, star_flag(rlnflag="splPprocImage",          splflag="ppvol",    string=.true.)]
        self%starfile%class3D%flags = [self%starfile%class3D%flags, star_flag(rlnflag="rlnSpectralIndex",       splflag="specind",  int=.true.)   ]
        self%starfile%class3D%flags = [self%starfile%class3D%flags, star_flag(rlnflag="rlnResolution",          splflag="specres")                ]
        self%starfile%class3D%flags = [self%starfile%class3D%flags, star_flag(rlnflag="rlnAngstromResolution",  splflag="specares")               ]
        self%starfile%class3D%flags = [self%starfile%class3D%flags, star_flag(rlnflag="rlnGoldStandardFsc",     splflag="specfsc")                ]
        self%starfile%class3D%flags = [self%starfile%class3D%flags, star_flag(rlnflag="rlnEstimatedResolution", splflag="speceres")               ]
        self%starfile%class3D%flags = [self%starfile%class3D%flags, star_flag(rlnflag="splFsc05",               splflag="fsc05"                  )]
        self%starfile%class3D%flags = [self%starfile%class3D%flags, star_flag(rlnflag="splFsc0128",             splflag="fsc0128"                )]
        
    end subroutine initialise

    ! import

    subroutine import_mics(self, cline, spproj, filename)
        class(starproject), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        class(sp_project),  intent(inout) :: spproj
        character(len=*),   intent(in)    :: filename
        integer                           :: i
        if( VERBOSE_OUTPUT )then
            write(logfhandle,*) ''
            write(logfhandle,*) char(9), 'importing ' // filename // " to mics"
            write(logfhandle,*) ''
        endif
        if(.not. self%starfile%initialised) call self%initialise()
        self%starfile%filename = filename
        self%starfile%rootdir = cline%get_carg("import_dir")
        call self%read_starheaders()
        call self%populate_opticsmap(spproj%os_optics)
        call self%import_stardata(self%starfile%micrographs, spproj%os_mic, .false., spproj%os_optics)
        do i =1,spproj%os_mic%get_noris()
            call spproj%os_mic%set(i, "imgkind", "mic")
            call spproj%os_mic%set(i, "ctf", "yes")
        end do
    end subroutine import_mics

    subroutine import_ptcls2D(self, cline, spproj, filename)
        class(starproject), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        class(sp_project),  intent(inout) :: spproj
        character(len=*),   intent(in)    :: filename
        integer                           :: i
        if( VERBOSE_OUTPUT )then
            write(logfhandle,*) ''
            write(logfhandle,*) char(9), 'importing ' // filename // " to ptcls2D"
            write(logfhandle,*)
        endif
        if(.not. self%starfile%initialised) call self%initialise()
        self%starfile%filename = filename
        self%starfile%rootdir = cline%get_carg("import_dir")
        call self%read_starheaders()
        call self%populate_opticsmap(spproj%os_optics)
        call self%populate_stkmap(spproj%os_stk, spproj%os_optics)
        call self%import_stardata(self%starfile%particles2D, spproj%os_ptcl2D, .true.)
        do i=1, spproj%os_stk%get_noris()
            call spproj%os_stk%set(i, "imgkind", "ptcl")
            call spproj%os_stk%set(i, "ctf", "yes")
            call spproj%os_stk%set(i, "stkkind", "split")
        end do
    end subroutine import_ptcls2D

    subroutine import_cls2D(self, cline, spproj, filename)
        class(starproject), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        class(sp_project),  intent(inout) :: spproj
        character(len=*)                  :: filename
        integer                           :: i
        if( VERBOSE_OUTPUT )then
            write(logfhandle,*) ''
            write(logfhandle,*) char(9), 'importing ' // filename // " to cls2D"
            write(logfhandle,*) ''
        endif
        if(.not. self%starfile%initialised) call self%initialise()
        self%starfile%filename = filename
        self%starfile%rootdir  = cline%get_carg("import_dir")
        call self%read_starheaders()
        call self%import_stardata(self%starfile%clusters2D, spproj%os_cls2D, .false.)
    end subroutine import_cls2D

    subroutine import_ptcls3D(self, cline, spproj, filename)
        class(starproject), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        class(sp_project),  intent(inout) :: spproj
        character(len=*),   intent(in)    :: filename
        integer                           :: i
        if( VERBOSE_OUTPUT )then
            write(logfhandle,*) ''
            write(logfhandle,*) char(9), 'importing ' // filename // " to ptcls3D"
            write(logfhandle,*)
        endif
        if(.not. self%starfile%initialised) call self%initialise()
        self%starfile%filename = filename
        self%starfile%rootdir  = cline%get_carg("import_dir")
        call self%read_starheaders()
        call self%populate_opticsmap(spproj%os_optics)
        call self%populate_stkmap(spproj%os_stk, spproj%os_optics)
        call self%import_stardata(self%starfile%particles3D, spproj%os_ptcl3D, .true.)
        do i=1, spproj%os_stk%get_noris()
            call spproj%os_stk%set(i, "imgkind", "ptcl")
            call spproj%os_stk%set(i, "ctf", "yes")
            call spproj%os_stk%set(i, "stkkind", "split")
        end do
    end subroutine import_ptcls3D

    subroutine import_stardata(self, stardata, sporis, isptcl, spoptics)
        class(starproject),    intent(inout) :: self
        class(star_data),      intent(inout) :: stardata
        class(oris),           intent(inout) :: sporis
        class(oris), optional, intent(inout) :: spoptics
        character(len=LEN_LINE), allocatable :: splitline(:)
        character(len=XLONGSTRLEN)           :: cwd
        character(len=LEN_LINE)              :: line, entrystr, splitimage
        character(len=LONGSTRLEN)            :: abspath
        type(ori)                            :: opticsori, spori
        logical                              :: isptcl
        real                                 :: rval
        integer                              :: ios, flagsindex, lineindex, ival, ogid, ogmapid, projindex, fhandle
        call simple_getcwd(cwd)
        allocate(splitline(stardata%flagscount))
        if(isptcl) then
            call sporis%new(self%starfile%stkptclcount, .true.)
            !set all stkind to 1
            do ival = 1, self%starfile%stkptclcount
               call sporis%set(ival, "stkind", real(1))
            end do 
        else
            call sporis%new(stardata%dataend - stardata%datastart, .false.)
        end if
        call fopen(fhandle, file=trim(adjustl(self%starfile%filename)), status='old')
        do lineindex = 1,stardata%datastart - 1
            read(fhandle, '(A)', iostat=ios) line
        end do
        do lineindex = 1, stardata%dataend - stardata%datastart
            if(isptcl) then
                projindex = self%starfile%stkmap(lineindex, 2)
            else
                projindex = lineindex
            end if
            read(fhandle, '(A)', iostat=ios) line
            call split_dataline(line, splitline)
            do flagsindex = 1,size(stardata%flags)
                if(stardata%flags(flagsindex)%present) then
                    if(stardata%flags(flagsindex)%imagesplit) then
                        splitimage = trim(adjustl(splitline(stardata%flags(flagsindex)%ind)))
                        if(.not. stardata%flags(flagsindex)%splflag2 == "") then
                            entrystr = trim(adjustl(splitimage(1:index(splitimage, "@") - 1)))
                            read(entrystr,*) ival
                            call sporis%set(projindex, stardata%flags(flagsindex)%splflag2, real(ival))
                        end if
                        entrystr = trim(adjustl(splitimage(index(splitimage, "@") + 1:)))
                    else
                        entrystr = trim(adjustl(splitline(stardata%flags(flagsindex)%ind)))
                    end if
                    if(stardata%flags(flagsindex)%string) then
                        if(index(entrystr, "/") > 1) then ! handles leading slash of absolute paths
                            if(index(entrystr, ':mrc') > 0) then ! relion ctf names!
                                entrystr = trim(adjustl(entrystr(:index(entrystr, ':mrc') - 1)))
                            end if
                            if(file_exists(trim(adjustl(stemname(self%starfile%rootdir))) // "/" // trim(adjustl(entrystr)))) then
                                call make_relativepath(cwd, trim(adjustl(stemname(self%starfile%rootdir))) // "/" // trim(adjustl(entrystr)), abspath, checkexists=.false.)
                                entrystr = trim(adjustl(abspath))
                            else
                                call make_relativepath(cwd, trim(adjustl(self%starfile%rootdir)) // "/" // trim(adjustl(entrystr)), abspath, checkexists=.false.)
                                entrystr = trim(adjustl(abspath))
                            end if
                        end if
                        call sporis%set(projindex, stardata%flags(flagsindex)%splflag, trim(adjustl(entrystr)))
                    else if(stardata%flags(flagsindex)%int) then
                        read(entrystr,*) ival
                        if(stardata%flags(flagsindex)%mult > 0) then
                            ival = ival * stardata%flags(flagsindex)%mult
                        end if
                        call sporis%set(projindex, stardata%flags(flagsindex)%splflag, real(ival))
                    else
                        read(entrystr,*) rval
                        if(stardata%flags(flagsindex)%mult > 0) then
                            rval = rval * stardata%flags(flagsindex)%mult
                        end if
                        call sporis%set(projindex, stardata%flags(flagsindex)%splflag, rval)
                    end if

                end if
            end do
            if(present(spoptics)) then
                ogid = sporis%get(lineindex, "ogid")
                if(ogid > 0) then
                    ogmapid = self%starfile%opticsmap(ogid)
                    call spoptics%get_ori(ogmapid, opticsori)
                    call sporis%get_ori(lineindex, spori)
                    call spori%append_ori(opticsori)
                    call sporis%set_ori(lineindex, spori)
                end if
            end if
            if(isptcl) then
                call sporis%set(projindex, "stkind", real(self%starfile%stkmap(lineindex, 1)))
            end if
            call sporis%set(projindex, "state", real(1))
        end do
        call fclose(fhandle)
        if(allocated(splitline)) deallocate(splitline)
    end subroutine import_stardata

    subroutine populate_stkmap(self, stkoris, opticsoris)
        class(starproject), intent(inout)    :: self
        class(oris),        intent(inout)    :: stkoris
        class(oris),        intent(inout)    :: opticsoris
        type(oris)                           :: stktmp
        type(ori)                            :: oritmpin, oritmpout
        character(len=LEN_LINE), allocatable :: stknames(:)
        character(len=LEN_LINE)              :: entrystr, searchstr
        integer, allocatable                 :: stkzmax(:), stkoriids(:), seppos(:)
        integer                              :: i, j, top, fromp, sepstart, sepend
        logical                              :: newstack, stkexists
        allocate(stknames(0))
        allocate(stkzmax(0))
        allocate(stkoriids(0))
        sepstart = 1
        sepend   = 1
        call self%import_stardata(self%starfile%stacks, stktmp, .false., opticsoris)
        allocate(self%starfile%stkmap(stktmp%get_noris(), 2))
        do i = 1,stktmp%get_noris()
            newstack = .true.
            do j = 1, size(stknames)
                if(trim(adjustl(stknames(j))) == trim(adjustl(stktmp%get_static(i, "stk")))) then
                    newstack = .false.
                    if(stkzmax(j) < int(stktmp%get(i, "stkind"))) then
                        stkzmax(j) = int(stktmp%get(i, "stkind"))
                    end if
                    self%starfile%stkmap(i, 1) = j
                    self%starfile%stkmap(i, 2) = int(stktmp%get(i, "stkind"))
                    exit
                end if
            end do
            if(newstack) then
                entrystr = trim(adjustl(stktmp%get_static(i, "stk")))
                stknames = [stknames, trim(adjustl(entrystr))]
                stkzmax = [stkzmax, int(stktmp%get(i, "stkind"))]
                stkoriids = [stkoriids, i]
                self%starfile%stkmap(i, 1) = size(stkoriids)
                self%starfile%stkmap(i, 2) = int(stktmp%get(i, "stkind"))
            end if
        end do
        call stkoris%new(size(stknames), .false.)
        allocate(self%starfile%stkstates(size(stknames)))
        do i = 1,size(stknames)
            call stktmp%get_ori(stkoriids(i), oritmpin)
            call stkoris%get_ori(i, oritmpout)
            call oritmpout%append_ori(oritmpin)
            call stkoris%set_ori(i, oritmpout)
            ! test each stk path. Fix wrong. State=0 missing
            stkexists = .false.
            allocate(seppos(0))
            entrystr = trim(adjustl(stkoris%get_static(i, "stk")))
            call find_separators(seppos, entrystr)
            ! test using last sepstart and seppos. If file doesn't exist, re-search
            searchstr = trim(adjustl(entrystr(:seppos(sepstart)))) // trim(adjustl(entrystr(seppos(sepend) + 1:)))
            if(file_exists(trim(adjustl(searchstr)))) then
                entrystr = searchstr
                stkexists = .true.
            else
                outer: do sepstart = 1, size(seppos) - 1
                    do sepend = sepstart + 1, size(seppos)
                        searchstr = trim(adjustl(entrystr(:seppos(sepstart)))) // trim(adjustl(entrystr(seppos(sepend) + 1:)))
                        if(file_exists(trim(adjustl(searchstr)))) then
                            entrystr = searchstr
                            stkexists = .true.
                            exit outer
                        end if
                    end do
                end do outer
            end if
            if(stkexists) then
                call stkoris%set(i, "stk", trim(adjustl(entrystr)))
                call stkoris%set(i, "state",  real(1))
                self%starfile%stkstates(i) = 1
            else
                write(logfhandle,*) char(9), "no stack found for ", trim(adjustl(basename(entrystr))), ". pruning referenced particles"
                call stkoris%set(i, "state",  real(0))
                self%starfile%stkstates(i) = 0
            end if
            if(allocated(seppos)) deallocate(seppos)
        end do
        fromp = 1
        self%starfile%stkptclcount = 0
        do i = 1, size(stknames)
            top = fromp + stkzmax(i)
            call stkoris%set(i, "nptcls", real(stkzmax(i)))
            call stkoris%set(i, "fromp",  real(fromp))
            call stkoris%set(i, "top",    real(top - 1))
            fromp = top
            self%starfile%stkptclcount = self%starfile%stkptclcount + stkzmax(i)
        end do
        do i = 1,stktmp%get_noris()
            self%starfile%stkmap(i, 2) =  self%starfile%stkmap(i, 2) + int(stkoris%get(self%starfile%stkmap(i, 1), "fromp")) - 1
        end do
        if(allocated(stknames))  deallocate(stknames)
        if(allocated(stkzmax))   deallocate(stkzmax)
        if(allocated(stkoriids)) deallocate(stkoriids)
        if(allocated(seppos))    deallocate(seppos)
    end subroutine populate_stkmap

    subroutine read_starheaders(self)
        class(starproject), intent(inout) :: self
        character(len=LEN_LINE)           :: line, blockname, flagname, datastring
        logical                           :: flagsopen, dataopen
        integer                           :: ios, delimindex,flagindex, istart, iend, flagid, lineindex, fhandle
        flagsopen = .false.
        dataopen  = .false.
        lineindex = 1
        call fopen(fhandle, file=trim(adjustl(self%starfile%filename)), status='old')
        do
            read(fhandle, '(A)', iostat=ios) line
            if(ios /= 0) exit
            if(len_trim(line) > 0) then
                line = trim(adjustl(line))
                if (index(line, "data_") > 0) then
                    ! start data block
                    flagindex  = 0
                    flagsopen  = .true.
                    dataopen   = .false.
                    ! is old or new format?
                    if(len_trim(line) > 5) then
                        ! new format
                        delimindex = index(line, "data_") + 5
                        blockname  = line(delimindex:)
                    else
                        blockname = 'particles'
                    end if
                else if(flagsopen .AND. (index(line, "_rln") > 0 .OR. index(line, "_spl") > 0)) then
                    delimindex = index(line, "#")
                    if(delimindex > 0) then
                        flagname = line(2:delimindex - 1)
                    else
                        flagname = trim(adjustl(line(2:)))
                    end if
                    flagindex = flagindex + 1
                    select case (blockname)
                        case ("optics")
                            call enable_rlnflag(flagname, self%starfile%optics%flags, flagindex)
                            self%starfile%optics%flagscount = self%starfile%optics%flagscount + 1
                        case ("micrographs")
                            call enable_rlnflag(flagname, self%starfile%micrographs%flags, flagindex)
                            self%starfile%micrographs%flagscount = self%starfile%micrographs%flagscount + 1
                        case ("particles")
                            call enable_rlnflag(flagname, self%starfile%particles2D%flags, flagindex)
                            self%starfile%particles2D%flagscount = self%starfile%particles2D%flagscount + 1
                            call enable_rlnflag(flagname, self%starfile%particles3D%flags, flagindex)
                            self%starfile%particles3D%flagscount = self%starfile%particles3D%flagscount + 1
                            call enable_rlnflag(flagname, self%starfile%stacks%flags, flagindex)
                            self%starfile%stacks%flagscount = self%starfile%stacks%flagscount + 1
                        case ("model_classes")
                            call enable_rlnflag(flagname, self%starfile%clusters2D%flags, flagindex)
                            self%starfile%clusters2D%flagscount = self%starfile%clusters2D%flagscount + 1
                    end select

                else if(flagsopen .AND. index(line, "loop_") == 0) then
                    dataopen = .true.
                    flagsopen = .false.
                    select case (blockname)
                        case ("optics")
                            self%starfile%optics%datastart      = lineindex
                        case ("micrographs")
                            self%starfile%micrographs%datastart = lineindex
                        case ("particles")
                            self%starfile%particles2D%datastart = lineindex
                            self%starfile%particles3D%datastart = lineindex
                            self%starfile%stacks%datastart      = lineindex
                        case ("model_classes")
                            self%starfile%clusters2D%datastart  = lineindex
                    end select
                end if
                if (dataopen .AND. index(line, "# ") == 0) then
                    select case (blockname)
                        case ("optics")
                            self%starfile%optics%dataend = lineindex + 1
                        case ("micrographs")
                            self%starfile%micrographs%dataend = lineindex + 1
                        case ("particles")
                            self%starfile%particles2D%dataend = lineindex + 1
                            self%starfile%particles3D%dataend = lineindex + 1
                            self%starfile%stacks%dataend = lineindex + 1
                        case ("model_classes")
                            self%starfile%clusters2D%dataend = lineindex + 1
                    end select
                end if
            end if
            lineindex = lineindex + 1
        end do
        print *, ""
        call fclose(fhandle)
    end subroutine read_starheaders
    
    subroutine check_stk_params(self, spproj)
        class(starproject), intent(inout) :: self
        class(sp_project),  intent(inout) :: spproj
        class(ImgHead), allocatable       :: header
        integer                           :: i, fromp, top, funit, ios
        real                              :: foundval, detpix, mag, ogid, stkbox, box
        character(len=LONGSTRLEN)         :: datapath
        character(len=1)                  :: format_descriptor
        character(len=7)                  :: stat_str
        do i = 1, spproj%os_stk%get_noris()
            fromp = nint(spproj%os_stk%get(i, "fromp"))
            top   = nint(spproj%os_stk%get(i, "fromp"))
            if(.not. spproj%os_stk%isthere(i, "smpd")) then
                foundval = get_value_from_ptcls(spproj%os_ptcl2D, fromp, top, "smpd")
                if(foundval > 0) then
                    call spproj%os_stk%set(i, "smpd", foundval)
                else 
                    detpix = get_value_from_ptcls(spproj%os_ptcl2D, fromp, top, "detpix")
                    mag    = get_value_from_ptcls(spproj%os_ptcl2D, fromp, top, "mag")
                    if(detpix > 0 .and. mag > 0) call spproj%os_stk%set(i, "smpd", 10000 * detpix/mag)
                end if
            end if
            if(.not. spproj%os_stk%isthere(i, "kv")) then
                foundval = get_value_from_ptcls(spproj%os_ptcl2D, fromp, top, "kv")
                if(foundval > 0) then
                    call spproj%os_stk%set(i, "kv", foundval)  
                end if
            end if
            if(.not. spproj%os_stk%isthere(i, "fraca")) then
                foundval = get_value_from_ptcls(spproj%os_ptcl2D, fromp, top, "fraca")
                if(foundval > 0) then
                    call spproj%os_stk%set(i, "fraca", foundval)  
                end if
            end if
            if(.not. spproj%os_stk%isthere(i, "cs")) then
                foundval = get_value_from_ptcls(spproj%os_ptcl2D, fromp, top, "cs")
                if(foundval > 0) then
                    call spproj%os_stk%set(i, "cs", foundval)  
                end if
            end if
           if(spproj%os_stk%get(i, "box") == 0.0) then
                datapath = spproj%os_stk%get_static(i, "stk")
                format_descriptor = fname2format(trim(adjustl(datapath)))
                select case(format_descriptor)
                    case ('M')
                        allocate(MrcImgHead :: header)
                        call header%new()
                    case ('S')
                        allocate(SpiImgHead :: header)
                        call header%new()
                    case ('J','L')
                        allocate(TiffImgHead :: header)
                        call header%new()
                    case DEFAULT
                        THROW_HARD('unsupported file format')
                end select
                if(file_exists(trim(adjustl(datapath)))) then
                    funit = 0
                    select case(format_descriptor)
                        case ('M','S')
                            call fopen(funit,access='STREAM',file=datapath,action='READ',status=stat_str,iostat=ios)
                            call fileiochk("imgfile::open_local fopen error: "//trim(datapath),ios)
                            call header%read(funit)
                            call fclose(funit)
                        case('J','L')
                            call header%read_tiff(datapath)
                    end select
                    call spproj%os_stk%set(i, "box", real(header%getDim(1)))
                end if
                if(allocated(header)) deallocate(header)
           end if
           ogid   = spproj%os_stk%get(i, 'ogid')
           stkbox = spproj%os_stk%get(i, 'box')
           if(ogid > 0.0 .and. stkbox > 0.0) then
                box = spproj%os_optics%get(nint(ogid), 'box')
                if(box == 0.0) then
                    call spproj%os_optics%set(nint(ogid), 'box', stkbox)
                end if
           end if
        end do
    end subroutine check_stk_params

    ! export

    subroutine export_mics(self, cline, spproj)
        class(starproject), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        class(sp_project),  intent(inout) :: spproj
        integer                           :: i
        if( L_VERBOSE_GLOB ) VERBOSE_OUTPUT = .true.
        self%starfile%filename = "micrographs.star"
        self%starfile%rootdir = cline%get_carg("import_dir")
        if( VERBOSE_OUTPUT )then
            write(logfhandle,*) ''
            write(logfhandle,*) char(9), 'exporting micrographs to ' // trim(adjustl(self%starfile%filename))
            write(logfhandle,*) ''
        endif
        if( file_exists(self%starfile%filename) ) call del_file(self%starfile%filename)
        if(.not. self%starfile%initialised) call self%initialise()
        call enable_splflags(spproj%os_optics, self%starfile%optics%flags)
        call enable_splflags(spproj%os_mic,    self%starfile%micrographs%flags)
        call self%export_stardata(spproj, self%starfile%optics%flags, spproj%os_optics, "optics", exclude="rlnImageSize")
        call self%export_stardata(spproj, self%starfile%micrographs%flags, spproj%os_mic, "micrographs")
    end subroutine export_mics

    subroutine export_cls2D(self, spproj, iter)
        class(starproject), intent(inout) :: self
        class(sp_project),  intent(inout) :: spproj
        integer, optional,  intent(in)    :: iter
        character(len=:), allocatable     :: stkname, relpath
        integer                           :: i, ncls
        real                              :: smpd
        if( L_VERBOSE_GLOB ) VERBOSE_OUTPUT = .true.
        if(present(iter)) then
            self%starfile%filename ="clusters2D_iter"//int2str_pad(iter,3)//".star"
        else
            self%starfile%filename = "clusters2D.star"
        end if
        if( VERBOSE_OUTPUT )then
            write(logfhandle,*) ''
            write(logfhandle,*) char(9), 'exporting clusters2D to ' // trim(adjustl(self%starfile%filename))
            write(logfhandle,*) ''
        endif
        if( file_exists(self%starfile%filename) ) call del_file(self%starfile%filename)
        if(.not. self%starfile%initialised) call self%initialise()
        if(present(iter)) then
            stkname = basename(CWD_GLOB)//"/"//trim(CAVGS_ITER_FBODY)//int2str_pad(iter,3)//trim(params_glob%ext)
            do i=1, spproj%os_cls2D%get_noris()
                call spproj%os_cls2D%set(i, "stk", trim(adjustl(stkname)))
            end do
        end if
        do i=1, spproj%os_cls2D%get_noris()
            if(.not. spproj%os_cls2D%isthere(i, "stk")) then
                call spproj%get_cavgs_stk(stkname, ncls, smpd)
                call spproj%os_cls2D%set(i, "stk", trim(adjustl(stkname)))
            end if
        end do
        call enable_splflags(spproj%os_cls2D, self%starfile%clusters2D%flags)
        call self%export_stardata(spproj, self%starfile%clusters2D%flags, spproj%os_cls2D, "clusters", mapstks=.false.)
    end subroutine export_cls2D

    subroutine export_iter3D(self, spproj, states, iter)
        class(starproject),                     intent(inout)   :: self
        class(sp_project),                      intent(inout)   :: spproj
        integer,                                intent(in)      :: states
        integer,                                intent(in)      :: iter
        type(sp_project)                                        :: classproj, fscproj
        character(len=:),           allocatable                 :: stkname, relpath
        character(len=STDLEN)                                   :: str_state, str_iter
        character(len=LONGSTRLEN)                               :: volassemble_output, pprocvol, lpvol
        character(len=XLONGSTRLEN)                              :: cwd
        character(len=LEN_LINE)                                 :: line
        character(len=LEN_LINE),    allocatable                 :: splitline(:)
        integer                                                 :: i, j, state, fhandle, ios, ok
        real                                                    :: smpd
        real,                       allocatable                 :: fscs(:,:,:), maxres05(:), maxres0128(:)
        logical                                                 :: ex
        if( L_VERBOSE_GLOB ) VERBOSE_OUTPUT = .true.
        call simple_getcwd(cwd)
        if(.not. self%starfile%initialised) call self%initialise()
        str_iter = int2str_pad(iter,3)
        allocate(fscs(states,2000,2))! assumed max 2000 px box
        allocate(maxres05(states))
        allocate(maxres0128(states))
        fscs       = 0.0
        maxres05   = 0.0
        maxres0128 = 0.0
        call classproj%os_cls3D%new(states, is_ptcl=.false.)
        i = 0
        do state = 1, states
            str_state = int2str_pad(state,2)
            volassemble_output = 'RESOLUTION_STATE'//trim(str_state)//'_ITER'//trim(str_iter)
            if(file_exists(volassemble_output)) then
                i = 0
                call fopen(fhandle, file=trim(adjustl(volassemble_output)), status='old')
                do
                    read(fhandle, '(A)', iostat=ios) line
                    if(ios /= 0) exit
                    if(len_trim(line) > 0) then
                        line = trim(adjustl(line))
                        if (index(line, "CORRELATION") > 0) then
                            i = i + 1
                            allocate(splitline(6))
                            call split_dataline(line, splitline)
                            read(splitline(3),*) fscs(state,i,1)
                            read(splitline(6),*) fscs(state,i,2)
                            deallocate(splitline)
                        else if (index(line, "FSC=0.500") > 0) then
                            allocate(splitline(7))
                            call split_dataline(line, splitline)
                            read(splitline(7),*) maxres05(state)
                            deallocate(splitline)
                        else if (index(line, "FSC=0.143") > 0) then
                            allocate(splitline(7))
                            call split_dataline(line, splitline)
                            read(splitline(7),*) maxres0128(state)
                            deallocate(splitline)
                        end if
                    end if
                end do
                call fclose(fhandle)
            end if
            pprocvol = trim(VOL_FBODY)//trim(str_state)//"_iter"//trim(str_iter)//trim(PPROC_SUFFIX)//".mrc"
            lpvol    = trim(VOL_FBODY)//trim(str_state)//"_iter"//trim(str_iter)//trim(LP_SUFFIX)//".mrc"
            if(file_exists(pprocvol)) then
              call classproj%os_cls3D%set(state, 'ppvol', trim(basename(cwd))//"/"//trim(pprocvol))
            end if
            if(file_exists(lpvol)) then
              call classproj%os_cls3D%set(state, 'lpvol', trim(basename(cwd))//"/"//trim(lpvol))
            end if
            call classproj%os_cls3D%set(state, 'state',    1.0)
            call classproj%os_cls3D%set(state, 'speceres', maxres0128(state))
            call classproj%os_cls3D%set(state, 'fsc05',    maxres05(state))
            call classproj%os_cls3D%set(state, 'fsc0128',  maxres0128(state))
        end do
        if(maxval(maxres0128) .eq. 0.0) GO TO 10
        self%starfile%filename ="class3D_iter"//trim(str_iter)//".star"
        inquire(file=trim(adjustl(self%starfile%filename)), exist=ex)
        if (ex) then
            call fopen(fhandle,file=trim(adjustl(self%starfile%filename)), position='append', iostat=ok)
        else
            call fopen(fhandle,file=trim(adjustl(self%starfile%filename)), status='new', iostat=ok)
        endif
        write(fhandle, *) ""
        write(fhandle, *) "data_model_general"
        write(fhandle, *) ""
        write(fhandle, "(A)")        "_rlnReferenceDimensionality       3"
        write(fhandle, "(A)")        "_rlnDataDimensionality            2"
        write(fhandle, "(A)")        "_rlnNrClasses                     " // states
        write(fhandle, "(A,F12.4)")  "_rlnEstimatedResolution           ", minval(maxres0128)
        call fclose(fhandle)
        call enable_splflags(classproj%os_cls3D, self%starfile%class3D%flags)
        call self%export_stardata(classproj, self%starfile%class3D%flags, classproj%os_cls3D, "model_classes", mapstks=.false.)
        do state = 1, states
            call fscproj%os_cls3D%new(i, is_ptcl=.false.)
            do j = 1, i
                call fscproj%os_cls3D%set(j, "state", 1.0)
                call fscproj%os_cls3D%set(j, "specind", real(j) - 1.0)
                if( j > 1) then
                    call fscproj%os_cls3D%set(j, "specres", 1.0/(real(j) - 1.0))
                else
                    call fscproj%os_cls3D%set(j, "specres", 999.99)
                end if 
                call fscproj%os_cls3D%set(j, "specares", real(fscs(state, j, 1)))
                call fscproj%os_cls3D%set(j, "specfsc",  real(fscs(state, j, 2)))
            end do
            call enable_splflags(fscproj%os_cls3D, self%starfile%class3D%flags)
            call self%export_stardata(fscproj, self%starfile%class3D%flags, fscproj%os_cls3D, "model_class_"//state, mapstks=.false.)
            call fscproj%kill()
        end do
        10 CONTINUE
        if(allocated(fscs))       deallocate(fscs)
        if(allocated(maxres05))   deallocate(maxres05)
        if(allocated(maxres0128)) deallocate(maxres0128)
        if(allocated(stkname))    deallocate(stkname)
        if(allocated(relpath))    deallocate(relpath)
        if(allocated(splitline))  deallocate(splitline)
    end subroutine export_iter3D

    subroutine export_ptcls2D(self, cline, spproj)
        class(starproject), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        class(sp_project),  intent(inout) :: spproj
        integer                           :: i
        if( L_VERBOSE_GLOB ) VERBOSE_OUTPUT = .true.
        self%starfile%filename = "particles2D.star"
        if( VERBOSE_OUTPUT )then
            write(logfhandle,*) ''
            write(logfhandle,*) char(9), 'exporting particles2D to ' // trim(adjustl(self%starfile%filename))
            write(logfhandle,*) ''
        endif
        if( file_exists(self%starfile%filename) ) call del_file(self%starfile%filename)
        if(.not. self%starfile%initialised) call self%initialise()
        call enable_splflags(spproj%os_optics, self%starfile%optics%flags)
        call enable_splflags(spproj%os_ptcl2D, self%starfile%particles2D%flags)
        call center_boxes(spproj, spproj%os_ptcl2D)
        call self%export_stardata(spproj, self%starfile%optics%flags, spproj%os_optics, "optics")
        call self%export_stardata(spproj, self%starfile%particles2D%flags, spproj%os_ptcl2D, "particles", mapstks=.true., exclude="rlnAmplitudeContrast")
    end subroutine export_ptcls2D

    subroutine export_ptcls3D(self, cline, spproj)
        class(starproject), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        class(sp_project)                 :: spproj
        integer                           :: i
        if( L_VERBOSE_GLOB ) VERBOSE_OUTPUT = .true.
        self%starfile%filename = "particles3D.star"
        if( VERBOSE_OUTPUT )then
            write(logfhandle,*) ''
            write(logfhandle,*) char(9), 'exporting particles3D to ' // trim(adjustl(self%starfile%filename))
            write(logfhandle,*) ''
        endif
        if( file_exists(self%starfile%filename) ) call del_file(self%starfile%filename)
        if(.not. self%starfile%initialised) call self%initialise()
        call enable_splflags(spproj%os_optics, self%starfile%optics%flags)
        call enable_splflags(spproj%os_ptcl3D, self%starfile%particles3D%flags)
        call self%export_stardata(spproj, self%starfile%optics%flags, spproj%os_optics, "optics")
        call self%export_stardata(spproj, self%starfile%particles3D%flags, spproj%os_ptcl3D, "particles", mapstks=.true., exclude="rlnAmplitudeContrast")
    end subroutine export_ptcls3D

    subroutine export_stardata(self, spproj, flags, sporis, blockname, mapstks, exclude)
        class(starproject), intent(inout)       :: self
        class(sp_project),  intent(inout)       :: spproj
        type(star_flag),    intent(in)          :: flags(:)
        class(oris),        intent(inout)       :: sporis
        character(len=*),   intent(in)          :: blockname
        character(len=*),optional,   intent(in) :: exclude
        logical, optional,  intent(in)          :: mapstks
        character(len=:), allocatable           :: stkname
        character(len=XLONGSTRLEN)              :: cwd
        character(len=LONGSTRLEN)               :: strtmp
        integer                                 :: iflag, iori, ok, stkindex, fhandle
        real                                    :: rval, stkind
        logical                                 :: ex, excludeflag
        call simple_getcwd(cwd)
        inquire(file=trim(adjustl(self%starfile%filename)), exist=ex)
        if (ex) then
            call fopen(fhandle,file=trim(adjustl(self%starfile%filename)), position='append', iostat=ok)
        else
            call fopen(fhandle,file=trim(adjustl(self%starfile%filename)), status='new', iostat=ok)
        endif
        write(fhandle, *) ""
        write(fhandle, *) "data_" // trim(adjustl(blockname))
        write(fhandle, *) ""
        write(fhandle, *) "loop_"
        do iflag=1, size(flags)
            excludeflag = .false.
            if(present(exclude) .and. index(trim(adjustl(flags(iflag)%rlnflag)), exclude) > 0 ) excludeflag = .true.
            if(flags(iflag)%present .and. .not. excludeflag) then
                write(fhandle, *) "_" // flags(iflag)%rlnflag
            end if
        end do
        if( present(mapstks) )then
            if( mapstks )then
                write(fhandle, *) "_rlnImageName"
                if(spproj%os_mic%get_noris() > 0) write(fhandle, *) "_rlnMicrographName"
            endif
        end if
        do iori=1, sporis%get_noris()
            if(sporis%get_state(iori) > 0) then
                do iflag=1, size(flags)
                    excludeflag = .false.
                    if(present(exclude) .and. index(trim(adjustl(flags(iflag)%rlnflag)), exclude) > 0 ) excludeflag = .true.
                    if(flags(iflag)%present .and. .not. excludeflag) then
                        if(flags(iflag)%imagesplit) then
                            rval = sporis%get(iori, trim(adjustl(flags(iflag)%splflag2)))
                            if(flags(iflag)%mult > 0) then
                                rval = rval / flags(iflag)%mult
                            end if
                            write(fhandle, "(I6,A)", advance="no") int(rval), "@"
                            strtmp = trim(adjustl(sporis%get_static(iori, trim(adjustl(flags(iflag)%splflag)))))
                            if(index(strtmp, '../') == 1) then ! Relative path. Adjust to base directory
                                write(fhandle, "(A)", advance="no")  trim(adjustl(strtmp(4:))) // " "
                            else
                                write(fhandle, "(A)", advance="no") trim(adjustl(sporis%get_static(iori, trim(adjustl(flags(iflag)%splflag))))) // " "
                            end if
                        else
                            if(flags(iflag)%string) then
                                strtmp = trim(adjustl(sporis%get_static(iori, trim(adjustl(flags(iflag)%splflag)))))
                                if(index(strtmp, '../') == 1) then ! Relative path. Adjust to base directory
                                    write(fhandle, "(A)", advance="no")  trim(adjustl(strtmp(4:))) // " "
                                else
                                    write(fhandle, "(A)", advance="no") trim(adjustl(sporis%get_static(iori, trim(adjustl(flags(iflag)%splflag))))) // " "
                                end if
                            elseif(flags(iflag)%int) then
                                rval = sporis%get(iori, trim(adjustl(flags(iflag)%splflag)))
                                if(flags(iflag)%mult > 0) then
                                    rval = rval / flags(iflag)%mult
                                end if
                                write(fhandle, "(I6,A)", advance="no") int(rval), " "
                            else
                                rval = sporis%get(iori, trim(adjustl(flags(iflag)%splflag)))
                                if(flags(iflag)%mult > 0) then
                                    rval = rval / flags(iflag)%mult
                                end if
                                write(fhandle, "(F12.4,A)", advance="no") rval, " "
                            end if
                        end if
                    end if
                end do
                if(present(mapstks) .and. mapstks) then
                    call spproj%get_stkname_and_ind('ptcl2D', iori, stkname, stkindex)
                    if(index(stkname, '../') == 1) then ! Relative path. Adjust to base directory
                        write(fhandle, "(I7,A,A,A)", advance="no") int(stkindex), '@', trim(adjustl(stkname(4:))), ' '
                    else
                        write(fhandle, "(I7,A,A,A)", advance="no") int(stkindex), '@', trim(adjustl(stkname)), ' '
                    end if
                    stkind = spproj%os_ptcl2D%get(iori, "stkind")
                    if(stkind <= spproj%os_mic%get_noris()) then
                        strtmp = trim(adjustl(spproj%os_mic%get_static(int(stkind), "intg")))
                        if(index(strtmp, '../') == 1) then ! Relative path. Adjust to base directory
                            write(fhandle, "(A)", advance="no")  trim(adjustl(strtmp(4:))) // " "
                        else
                            write(fhandle, "(A)", advance="no") trim(adjustl(strtmp)) // " "
                        end if
                    end if
                end if
                write(fhandle,*) ""
            end if
        end do
        call fclose(fhandle)
        if(allocated(stkname)) deallocate(stkname)
    end subroutine export_stardata

    ! tilt

    subroutine assign_initial_tiltgroups(self, beamtilt)
        class(starproject),      intent(inout) :: self
        character(len=LONGSTRLEN), allocatable :: epugroups(:)
        character(len=LONGSTRLEN)              :: eputilt
        integer                                :: i, epugroupid
        logical, intent(in)                    :: beamtilt
        if( VERBOSE_OUTPUT ) write(logfhandle,*) ''
        if(.not. beamtilt) then
            if( VERBOSE_OUTPUT ) write(logfhandle,*) char(9), 'beamtilts not being used. assigning single initial beamtilt group'
        else if(index(self%tiltinfo(1)%basename, 'FoilHole') == 0) then
            if( VERBOSE_OUTPUT ) write(logfhandle,*) char(9), 'no EPU filenames detected. assigning single initial beamtilt group'
            ! already assigned initialtiltgroupid=1 in initialisation
        else
            if( VERBOSE_OUTPUT ) write(logfhandle,"(A,A)", advance="no") char(9), 'EPU filenames detected. using these to assign initial beamtilt groups ... '
            allocate(epugroups(0))
            epugroupid = 1
            do i =1,size(self%tiltinfo)
                eputilt = self%tiltinfo(i)%basename(index(self%tiltinfo(i)%basename,'Data_') + 5:)
                eputilt = eputilt(:index(eputilt,'_') - 1)
                if(any(epugroups == eputilt)) then
                self%tiltinfo(i)%initialtiltgroupid = findloc(epugroups, eputilt, 1)
                else
                    epugroups = [epugroups, eputilt]
                    self%tiltinfo(i)%initialtiltgroupid = epugroupid
                    epugroupid = epugroupid + 1
                end if
            end do
            if( VERBOSE_OUTPUT ) write(logfhandle,"(A,I4,A)") "assigned ", size(epugroups), " initial groups"
            deallocate(epugroups)
        end if
        if( VERBOSE_OUTPUT ) write(logfhandle,*) ''
        if(allocated(epugroups)) deallocate(epugroups)
    end subroutine assign_initial_tiltgroups

    subroutine assign_xml_tiltinfo(self, xmldir)
        class(starproject), target, intent(inout) :: self
        character(len=*),           intent(in)    :: xmldir
        character(len=LONGSTRLEN)                 :: eputilt
        type(Node), pointer                       :: xmldoc, beamtiltnode, beamtiltnodex, beamtiltnodey
        integer                                   :: i, j
        if( VERBOSE_OUTPUT ) write(logfhandle,*) ''
        if( VERBOSE_OUTPUT ) write(logfhandle,*) char(9), "reading tilt info from metadata ... "
        do i = 1,size(self%tiltinfo)
            if(file_exists(xmldir // '/' // trim(adjustl(self%tiltinfo(i)%basename)) // '.xml')) then
                xmldoc => parseFile(xmldir // '/' // trim(adjustl(self%tiltinfo(i)%basename)) // '.xml')
                beamtiltnode => item(getElementsByTagname(xmldoc, "BeamShift"), 0)
                beamtiltnodex => item(getElementsByTagname(beamtiltnode, "a:_x"), 0)
                beamtiltnodey => item(getElementsByTagname(beamtiltnode, "a:_y"), 0)
                self%tiltinfo(i)%tiltx = str2real(getTextContent(beamtiltnodex))
                self%tiltinfo(i)%tilty = str2real(getTextContent(beamtiltnodey))
                call destroy(xmldoc)
            else
                if( VERBOSE_OUTPUT ) write(logfhandle, *) char(9), char(9), xmldir // '/' // trim(adjustl(self%tiltinfo(i)%basename)) // '.xml does not exist. Ignoring'
            end if
        end do
        if( VERBOSE_OUTPUT ) write(logfhandle,*) ''
    end subroutine assign_xml_tiltinfo

    subroutine cluster_tiltinfo(self, threshold)
        class(starproject), intent(inout) :: self
        real,               intent(in)    :: threshold
        integer, allocatable              :: populations(:), labels(:), tiltinfopos(:)
        real,    allocatable              :: tilts(:,:), centroids(:,:)
        integer :: i, j, k, tiltcount, groupcount, matchcount
        if( VERBOSE_OUTPUT ) write(logfhandle,*) ''
        if( VERBOSE_OUTPUT ) write(logfhandle,*) char(9), "clustering initial beamtilt groups using tilt info and a threshold of ", threshold, " ... "
        groupcount = 0
        do i = 1, maxval(self%tiltinfo%initialtiltgroupid)
            if( VERBOSE_OUTPUT ) write(logfhandle,*) char(9), char(9), "clustering initial beamtilt group ", i, " ... "
            matchcount = count(self%tiltinfo%initialtiltgroupid == i)
            allocate(tilts(matchcount, 2))
            allocate(labels(matchcount))
            allocate(tiltinfopos(matchcount))
            tiltcount = 1
            do j = 1, size(self%tiltinfo)
                if(self%tiltinfo(j)%initialtiltgroupid == i) then
                    tilts(tiltcount,1) = self%tiltinfo(j)%tiltx
                    tilts(tiltcount,2) = self%tiltinfo(j)%tilty
                    tiltinfopos(tiltcount) = j
                    tiltcount = tiltcount + 1
                end if
            end do
            call h_clust(tilts, threshold, labels, centroids, populations)
            do j = 1, size(tiltinfopos)
                self%tiltinfo(tiltinfopos(j))%finaltiltgroupid = groupcount + labels(j)
            end do
            groupcount = groupcount + size(populations)
            deallocate(tilts, labels, tiltinfopos)
            if( VERBOSE_OUTPUT ) write(logfhandle,*) ''
        end do
        if(allocated(populations)) deallocate(populations)
        if(allocated(labels))      deallocate(labels)
        if(allocated(tiltinfopos)) deallocate(tiltinfopos)
        if(allocated(tilts))       deallocate(tilts)
        if(allocated(centroids))   deallocate(centroids)
    end subroutine cluster_tiltinfo

    ! optics

    subroutine assign_optics(self, cline, spproj)
        class(starproject),    intent(inout)   :: self
        class(cmdline),        intent(inout)   :: cline
        class(sp_project)                      :: spproj
        integer                                :: i, j, element, noptics
        integer                                :: ptclstkid
        real                                   :: ogid, box, stkbox
        character(len=LONGSTRLEN)              :: ogname
        character(len=XLONGSTRLEN)             :: cwd
        logical                                :: beamtilt
        beamtilt  = .true.
        if(cline%get_carg('beamtilt') .eq. 'no') then
            beamtilt = .false.
        end if
        allocate(self%tiltinfo(0))
        if(spproj%os_mic%get_noris() > 0) then
            call self%get_image_basename(spproj%os_mic, 'intg')
        end if
        if(spproj%os_stk%get_noris() > 0) then
            call self%get_image_basename(spproj%os_stk, 'stk')
        end if
        call self%assign_initial_tiltgroups(beamtilt)
        if(cline%get_carg("xmldir") .ne. '') then
            if(.not. index(cline%get_carg("xmldir"), "/") == 1) then
                call simple_getcwd(cwd)
                call cline%set('xmldir',  trim(adjustl(stemname(cwd))) // "/" // trim(adjustl(cline%get_carg("xmldir"))))
            end if
            call self%assign_xml_tiltinfo(cline%get_carg("xmldir"))
        end if
        call self%cluster_tiltinfo(cline%get_rarg("tilt_thres"))
        if(cline%get_rarg("maxpop") > 0) then
            if( VERBOSE_OUTPUT ) write(logfhandle,*) ''
            if( VERBOSE_OUTPUT ) write(logfhandle,*) char(9), "plotting beamtilt groups (prior to splitting on maxpop) ... "
            call self%plot_opticsgroups("optics_groups_pre_maxpop.eps")
            call self%sort_optics_maxpop(int(cline%get_rarg("maxpop")))
            if( VERBOSE_OUTPUT ) write(logfhandle,*) ''
            if( VERBOSE_OUTPUT ) write(logfhandle,*) char(9), "plotting beamtilt groups (after splitting on maxpop) ... "
            call self%plot_opticsgroups("optics_groups_post_maxpop.eps")
        else
            if( VERBOSE_OUTPUT ) write(logfhandle,*) ''
            if( VERBOSE_OUTPUT ) write(logfhandle,*) char(9), "plotting beamtilt groups ... "
            call self%plot_opticsgroups("optics_groups.eps")
        end if
        if(cline%get_rarg("optics_offset") > 0) then
            call self%apply_optics_offset(int(cline%get_rarg("optics_offset")))
        end if
        if(spproj%os_mic%get_noris() > 0) then
            if( VERBOSE_OUTPUT ) write(logfhandle,*) ''
            if( VERBOSE_OUTPUT ) write(logfhandle,*) char(9), "updating micrographs in project file with updated optics groups ... "
            do i = 1,spproj%os_mic%get_noris()
                element = findloc(self%tiltinfo%basename, trim(adjustl(spproj%os_mic%get_static(i, "bsname"))), 1)
                if(element > 0) then
                    call spproj%os_mic%set(i, 'ogid', real(self%tiltinfo(element)%finaltiltgroupid))
                else
                    THROW_HARD('failed to locate optics info for mic basename' // trim(adjustl(spproj%os_mic%get_static(i, "bsname"))))
                end if
            end do
        end if
        if(spproj%os_stk%get_noris() > 0) then
            if( VERBOSE_OUTPUT ) write(logfhandle,*) ''
            if( VERBOSE_OUTPUT ) write(logfhandle,*) char(9), "updating stacks in project file with updated optics groups ... "
            do i = 1, spproj%os_stk%get_noris()
                element = findloc(self%tiltinfo%basename, trim(adjustl(spproj%os_stk%get_static(i, "bsname"))), 1)
                if(element > 0) then
                    call spproj%os_stk%set(i, 'ogid', real(self%tiltinfo(element)%finaltiltgroupid))
                else
                    THROW_HARD('failed to locate optics info for stk basename' // trim(adjustl(spproj%os_stk%get_static(i, "bsname"))))
                end if
            end do
            if( VERBOSE_OUTPUT ) write(logfhandle,*) ''
            if( VERBOSE_OUTPUT ) write(logfhandle,*) char(9), "updating optics groups in project file with box sizes ... "
            noptics = spproj%os_optics%get_noris()
            do i = 1, spproj%os_stk%get_noris()
                ogid   = spproj%os_stk%get(i, 'ogid')
                stkbox = spproj%os_stk%get(i, 'box')
                if(ogid > 0.0 .and. stkbox > 0.0 .and. ogid <= noptics) then
                    box = spproj%os_optics%get(nint(ogid), 'box')
                    if(box == 0.0) then
                        call spproj%os_optics%set(nint(ogid), 'box', stkbox)
                    end if
                end if
            end do
        end if
        if(spproj%os_ptcl2D%get_noris() > 0) then
            if( VERBOSE_OUTPUT ) write(logfhandle,*) ''
            if( VERBOSE_OUTPUT ) write(logfhandle,*) char(9), "updating particles 2d in project file with updated optics groups ... "
            do i = 1, spproj%os_ptcl2D%get_noris()
                ptclstkid = int(spproj%os_ptcl2D%get(i, 'stkind'))
                call spproj%os_ptcl2D%set(i, 'ogid', spproj%os_stk%get(ptclstkid, 'ogid'))
            end do
        end if
        if(spproj%os_ptcl3D%get_noris() > 0) then
            if( VERBOSE_OUTPUT ) write(logfhandle,*) ''
            if( VERBOSE_OUTPUT ) write(logfhandle,*) char(9), "updating particles 3d in project file with updated optics groups ... "
            do i = 1, spproj%os_ptcl3D%get_noris()
                ptclstkid = int(spproj%os_ptcl2D%get(i, 'stkind'))
                call spproj%os_ptcl3D%set(i, 'ogid', spproj%os_stk%get(ptclstkid, 'ogid'))
            end do
        end if

        if(size(self%tiltinfo) > 0) then
            if( VERBOSE_OUTPUT ) write(logfhandle,*) ''
            if( VERBOSE_OUTPUT ) write(logfhandle,*) char(9), "updating optics groups in project file ... "

            call spproj%os_optics%new(maxval(self%tiltinfo%finaltiltgroupid) - minval(self%tiltinfo%finaltiltgroupid) + 1, is_ptcl=.false.)
            do i = minval(self%tiltinfo%finaltiltgroupid), maxval(self%tiltinfo%finaltiltgroupid)
                element = findloc(self%tiltinfo%finaltiltgroupid, i, 1)
                if(element > 0) then
                    write(ogname,"(I6)") i
                    call spproj%os_optics%set(i - minval(self%tiltinfo%finaltiltgroupid) + 1, "ogid", real(i))
                    call spproj%os_optics%set(i - minval(self%tiltinfo%finaltiltgroupid) + 1, "ogname", "opticsgroup" // trim(adjustl(ogname)))
                    call spproj%os_optics%set(i - minval(self%tiltinfo%finaltiltgroupid) + 1, "smpd", self%tiltinfo(element)%smpd)
                    call spproj%os_optics%set(i - minval(self%tiltinfo%finaltiltgroupid) + 1, "cs", self%tiltinfo(element)%cs)
                    call spproj%os_optics%set(i - minval(self%tiltinfo%finaltiltgroupid) + 1, "kv", self%tiltinfo(element)%kv)
                    call spproj%os_optics%set(i - minval(self%tiltinfo%finaltiltgroupid) + 1, "fraca", self%tiltinfo(element)%fraca)
                    call spproj%os_optics%set(i - minval(self%tiltinfo%finaltiltgroupid) + 1, "box", self%tiltinfo(element)%box)
                    call spproj%os_optics%set(i - minval(self%tiltinfo%finaltiltgroupid) + 1, "state", 1.0)
                end if
            end do
        end if
        if( VERBOSE_OUTPUT ) write(logfhandle,*) ''
    end subroutine assign_optics

    subroutine populate_opticsmap(self, opticsoris)
        class(starproject), intent(inout) :: self
        class(oris),        intent(inout) :: opticsoris
        integer                           :: i
        call self%import_stardata(self%starfile%optics, opticsoris, .false.)
        allocate(self%starfile%opticsmap(opticsoris%get_noris()))
        do i = 1,opticsoris%get_noris()
            self%starfile%opticsmap(int(opticsoris%get(i, "ogid"))) = i
        end do
    end subroutine populate_opticsmap

    subroutine plot_opticsgroups(self, fname_eps)
        class(starproject), intent(inout) :: self
        character(len=*),   intent(in)    :: fname_eps
        type(str4arr)       :: title
        type(CPlot2D_type)  :: plot2D
        type(CDataSet_type) :: dataSet
        integer :: i, j
        call CPlot2D__new(plot2D, 'Optics Groups'//C_NULL_CHAR)
        call CPlot2D__SetXAxisSize(plot2D, 400.0_c_double)
        call CPlot2D__SetYAxisSize(plot2D, 400.0_c_double)
        call CPlot2D__SetDrawXAxisGridLines(plot2D, C_TRUE)
        call CPlot2D__SetDrawYAxisGridLines(plot2D, C_TRUE)
        call CPlot2D__SetDrawLegend(plot2D, C_FALSE)
        call CPlot2D__SetFlipY(plot2D, C_TRUE)
        title%str = 'Beamshift_X'//C_NULL_CHAR
        call CPlot2D__SetXAxisTitle(plot2D, title%str)
        title%str = 'Beamshift_Y'//C_NULL_CHAR
        call CPlot2D__SetYAxisTitle(plot2D, title%str)
        do i=1, maxval(self%tiltinfo%finaltiltgroupid)
            call CDataSet__new(dataSet)
            call CDataSet__SetDrawMarker(dataSet, C_TRUE)
            call CDataSet__SetDrawLine(dataSet, C_FALSE)
            call CDataSet__SetMarkerSize(dataSet, real(3.0, c_double))
            call CDataSet__SetDatasetColor(dataSet, real(ran3(), c_double), real(ran3(), c_double), real(ran3(), c_double))
            do j=1, size(self%tiltinfo)
                if(self%tiltinfo(j)%finaltiltgroupid .eq. i) then
                    call CDataSet_addpoint(dataSet, self%tiltinfo(j)%tiltx, self%tiltinfo(j)%tilty)
                end if
            end do
            call CPlot2D__AddDataSet(plot2D, dataset)
            call CDataSet__delete(dataset)
        end do
        call CPlot2D__OutputPostScriptPlot(plot2D, fname_eps//C_NULL_CHAR)
        call CPlot2D__delete(plot2D)
        if( VERBOSE_OUTPUT ) write(logfhandle,*) char(9), char(9), "wrote ", fname_eps
    end subroutine plot_opticsgroups

    subroutine sort_optics_maxpop(self, maxpop)
        class(starproject), intent(inout) :: self
        integer,            intent(in)    :: maxpop
        integer, allocatable :: nextidpop(:,:)
        integer :: i
        allocate(nextidpop(maxval(self%tiltinfo%finaltiltgroupid) , 2))
        do i = 1, size(nextidpop,1)
            nextidpop(i,1) = i
            nextidpop(i,2) = 0
        end do
        do i = 1, size(self%tiltinfo)
            nextidpop(self%tiltinfo(i)%finaltiltgroupid,2) = nextidpop(self%tiltinfo(i)%finaltiltgroupid,2) + 1
            if(nextidpop(self%tiltinfo(i)%finaltiltgroupid, 2) > maxpop) then
                nextidpop(self%tiltinfo(i)%finaltiltgroupid,1) = nextidpop(self%tiltinfo(i)%finaltiltgroupid,1) + size(nextidpop, 1)
                nextidpop(self%tiltinfo(i)%finaltiltgroupid,2) = 0
            end if
            self%tiltinfo(i)%finaltiltgroupid = nextidpop(self%tiltinfo(i)%finaltiltgroupid,1)
        end do
        deallocate(nextidpop)
    end subroutine sort_optics_maxpop

    subroutine apply_optics_offset(self, offset)
        class(starproject), intent(inout) :: self
        integer,            intent(in)    :: offset
        integer :: i
        do i = 1,size(self%tiltinfo)
            self%tiltinfo(i)%finaltiltgroupid = self%tiltinfo(i)%finaltiltgroupid + offset
        end do
    end subroutine apply_optics_offset

    subroutine get_image_basename(self, sporis, splflag)
        class(starproject), intent(inout) :: self
        class(oris),        intent(inout) :: sporis
        character(len=*),   intent(in)    :: splflag
        type(tilt_info)            :: tiltinfo
        character(len=XLONGSTRLEN) :: path
        character(len=LONGSTRLEN)  :: filename
        integer                    :: i, tiltind
        do i=1, sporis%get_noris()
            path = sporis%get_static(i, splflag)
            filename = basename(path)
            if(index(filename, INTGMOV_SUFFIX) > 0) then
                filename = filename(:index(filename, INTGMOV_SUFFIX) - 1)
            end if
            if(index(filename, EXTRACT_STK_FBODY) > 0) then
                filename = filename(index(filename, EXTRACT_STK_FBODY) + len(EXTRACT_STK_FBODY):)
            end if
            if(index(filename, '_fractions') > 0) then
                filename = filename(:index(filename, '_fractions') - 1)
            end if
            if(index(filename, '_EER') > 0) then
                filename = filename(:index(filename, '_EER') - 1)
            end if
            call sporis%set(i, "bsname", trim(adjustl(filename)))
            tiltind = findloc(self%tiltinfo%basename, filename, 1)
            if(tiltind == 0) then
                tiltinfo%basename = filename
                tiltinfo%smpd     = sporis%get(i, "smpd")
                tiltinfo%kv       = sporis%get(i, "kv")
                tiltinfo%fraca    = sporis%get(i, "fraca")
                tiltinfo%cs       = sporis%get(i, "cs")
                tiltinfo%box      = sporis%get(i, "box")
                if(sporis%isthere(i,'tiltx')) then
                    tiltinfo%tiltx = sporis%get(i, "tiltx")
                else
                    tiltinfo%tiltx = 0.0
                end if
                if(sporis%isthere(i,'tilty')) then
                    tiltinfo%tilty = sporis%get(i, "tilty")
                else
                    tiltinfo%tilty = 0.0
                end if
                self%tiltinfo     = [self%tiltinfo, tiltinfo]
            else
                if(self%tiltinfo(tiltind)%box < sporis%get(i, "box")) then
                    self%tiltinfo(tiltind)%box = sporis%get(i, "box")
                end if
            end if
        end do
    end subroutine get_image_basename

    subroutine set_verbose(self)
        class(starproject), intent(inout) :: self
        VERBOSE_OUTPUT = .true.
    end subroutine set_verbose
    
    subroutine kill(self)
        class(starproject), intent(inout) :: self
        if(allocated(self%tiltinfo))                   deallocate(self%tiltinfo)
        if(allocated(self%starfile%opticsmap))         deallocate(self%starfile%opticsmap)
        if(allocated(self%starfile%stkmap))            deallocate(self%starfile%stkmap)
        if(allocated(self%starfile%stkstates))         deallocate(self%starfile%stkstates)
        if(allocated(self%starfile%optics%flags))      deallocate(self%starfile%optics%flags)
        if(allocated(self%starfile%stacks%flags))      deallocate(self%starfile%stacks%flags)
        if(allocated(self%starfile%micrographs%flags)) deallocate(self%starfile%micrographs%flags)
        if(allocated(self%starfile%particles2D%flags)) deallocate(self%starfile%particles2D%flags)
        if(allocated(self%starfile%particles3D%flags)) deallocate(self%starfile%particles3D%flags)
        if(allocated(self%starfile%clusters2D%flags))  deallocate(self%starfile%clusters2D%flags)
        if(allocated(self%starfile%class3D%flags))     deallocate(self%starfile%class3D%flags)
    end subroutine kill

end module simple_starproject
