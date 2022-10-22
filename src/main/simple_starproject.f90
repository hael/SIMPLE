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
    ! export
    procedure          :: export_mics
    procedure          :: export_cls2D
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

        ! assign clusters2D flags
        if (.not. allocated(self%starfile%clusters2D%flags)) allocate(self%starfile%clusters2D%flags(0))
        self%starfile%clusters2D%flags = [self%starfile%clusters2D%flags, star_flag(rlnflag="rlnReferenceImage", splflag="stk", string=.true., imagesplit=.true., splflag2="class")]
        self%starfile%clusters2D%flags = [self%starfile%clusters2D%flags, star_flag(rlnflag="rlnEstimatedResolution", splflag="res")]
        self%starfile%clusters2D%flags = [self%starfile%clusters2D%flags, star_flag(rlnflag="rlnClassDistribution", splflag="pop", int=.true.)]

    end subroutine initialise

    ! import

    subroutine import_mics(self, cline, spproj, filename)
        class(starproject), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        class(sp_project),  intent(inout) :: spproj
        character(len=*),   intent(in)    :: filename
        integer :: i
        if( VERBOSE_OUTPUT )then
            write(logfhandle,*) ''
            write(logfhandle,*) char(9), 'importing3 ' // filename // " to mics"
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
        deallocate(self%starfile%opticsmap)
    end subroutine import_mics

    subroutine import_ptcls2D(self, cline, spproj, filename)
        class(starproject), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        class(sp_project),  intent(inout) :: spproj
        character(len=*),   intent(in)    :: filename
        integer :: i

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
        deallocate(self%starfile%opticsmap, self%starfile%stkmap)
    end subroutine import_ptcls2D

    subroutine import_cls2D(self, cline, spproj, filename)
        class(starproject), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        class(sp_project),  intent(inout) :: spproj
        character(len=*) :: filename
        integer :: i
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
        integer :: i
        if( VERBOSE_OUTPUT )then
            write(logfhandle,*) ''
            write(logfhandle,*) char(9), 'importing ' // filename // " to ptcls3D"
            write(logfhandle,*)
        endif
        if(.not. self%starfile%initialised) call self%initialise()
        self%starfile%filename = filename
        self%starfile%rootdir = cline%get_carg("import_dir")
        call self%read_starheaders()
        call self%populate_opticsmap(spproj%os_optics)
        call self%populate_stkmap(spproj%os_stk, spproj%os_optics)
        call self%import_stardata(self%starfile%particles3D, spproj%os_ptcl3D, .true.)
        do i=1, spproj%os_stk%get_noris()
            call spproj%os_stk%set(i, "imgkind", "ptcl")
            call spproj%os_stk%set(i, "ctf", "yes")
            call spproj%os_stk%set(i, "stkkind", "split")
        end do
        deallocate(self%starfile%opticsmap)
        deallocate(self%starfile%stkmap)
    end subroutine import_ptcls3D

    subroutine import_stardata(self, stardata, sporis, isptcl, spoptics)
        class(starproject),    intent(inout) :: self
        class(star_data),      intent(inout) :: stardata
        class(oris),           intent(inout) :: sporis
        class(oris), optional, intent(inout) :: spoptics
        character(len=LEN_LINE), allocatable :: splitline(:)
        character(len=XLONGSTRLEN) :: cwd
        character(len=LEN_LINE)    :: line, entrystr, splitimage
        character(len=LONGSTRLEN)  :: abspath
        type(ori) :: opticsori, spori
        logical   :: isptcl
        real      :: rval
        integer   :: ios, flagsindex, lineindex, ival, ogid, ogmapid, projindex, fhandle
        call simple_getcwd(cwd)
        allocate(splitline(stardata%flagscount))
        if(isptcl) then
            call sporis%new(self%starfile%stkptclcount, .true.)
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
                            entrystr = splitimage(1:index(splitimage, "@") - 1)
                            read(entrystr,*) ival
                            call sporis%set(projindex, stardata%flags(flagsindex)%splflag2, real(ival))
                        end if
                        entrystr = splitimage(index(splitimage, "@") + 1:)
                    else
                        entrystr = trim(adjustl(splitline(stardata%flags(flagsindex)%ind)))
                    end if
                    if(stardata%flags(flagsindex)%string) then
                        if(index(entrystr, "/") > 1) then ! handles leading slash of absolute paths
                            if(index(self%starfile%rootdir, "job") > 0) then ! relion
                                if(index(entrystr, ":mrc") > 0) then ! relion ctf names!
                                    entrystr = trim(adjustl(entrystr(:index(entrystr, ":mrc") - 1)))
                                end if
                                call make_relativepath(cwd, stemname(stemname(trim(adjustl(self%starfile%rootdir)))) // "/" // trim(adjustl(entrystr)), abspath, checkexists=.false.)
                                call sporis%set(projindex, stardata%flags(flagsindex)%splflag, trim(adjustl(abspath)))
                            else ! other\
                                call make_relativepath(cwd, stemname(trim(adjustl(self%starfile%rootdir))) // "/" // trim(adjustl(entrystr)), abspath, checkexists=.false.)
                                call sporis%set(projindex, stardata%flags(flagsindex)%splflag, trim(adjustl(abspath)))
                            end if
                        else
                            call sporis%set(projindex, stardata%flags(flagsindex)%splflag, trim(adjustl(entrystr)))
                        end if
                    elseif(stardata%flags(flagsindex)%int) then
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
        deallocate(splitline)
        call fclose(fhandle)
    end subroutine import_stardata

    subroutine populate_stkmap(self, stkoris, opticsoris)
        class(starproject), intent(inout) :: self
        class(oris),        intent(inout) :: stkoris
        class(oris),        intent(inout) :: opticsoris
        type(oris) :: stktmp
        type(ori)  :: oritmpin, oritmpout
        character(len=LEN_LINE), allocatable :: stknames(:)
        integer, allocatable :: stkzmax(:), stkoriids(:)
        integer :: i, j, top, fromp
        logical :: newstack
        allocate(stknames(0))
        allocate(stkzmax(0))
        allocate(stkoriids(0))
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
                stknames = [stknames, trim(adjustl(stktmp%get_static(i, "stk")))]
                stkzmax = [stkzmax, int(stktmp%get(i, "stkind"))]
                stkoriids = [stkoriids, i]
                self%starfile%stkmap(i, 1) = size(stkoriids)
                self%starfile%stkmap(i, 2) = int(stktmp%get(i, "stkind"))
            end if
        end do
        call stkoris%new(size(stknames), .false.)
        do i = 1,size(stknames)
            call stktmp%get_ori(stkoriids(i), oritmpin)
            call stkoris%get_ori(i, oritmpout)
            call oritmpout%append_ori(oritmpin)
            call stkoris%set_ori(i, oritmpout)
        end do
        fromp = 1
        self%starfile%stkptclcount = 0
        do i = 1, size(stknames)
            top = fromp + stkzmax(i)
            call stkoris%set(i, "nptcls", real(stkzmax(i)))
            call stkoris%set(i, "fromp", real(fromp))
            call stkoris%set(i, "top",  real(top - 1))
            call stkoris%set(i, "state",  real(1))
            fromp = top
            self%starfile%stkptclcount = self%starfile%stkptclcount + stkzmax(i)
        end do
        do i = 1,stktmp%get_noris()
            self%starfile%stkmap(i, 2) =  self%starfile%stkmap(i, 2) + int(stkoris%get(self%starfile%stkmap(i, 1), "fromp")) - 1
        end do
        deallocate(stknames)
        deallocate(stkzmax)
        deallocate(stkoriids)
    end subroutine populate_stkmap

    subroutine read_starheaders(self)
        class(starproject), intent(inout) :: self
        character(len=LEN_LINE) :: line, blockname, flagname, datastring
        logical :: flagsopen, dataopen
        integer :: ios, delimindex,flagindex, istart, iend, flagid, lineindex, fhandle
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
                    delimindex = index(line, "data_") + 5
                    blockname  = line(delimindex:)
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
        if(.not. self%starfile%initialised) call self%initialise()
        call enable_splflags(spproj%os_optics, self%starfile%optics%flags)
        call enable_splflags(spproj%os_mic,    self%starfile%micrographs%flags)
        call self%export_stardata(spproj, self%starfile%optics%flags, spproj%os_optics, "optics")
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
        if(.not. self%starfile%initialised) call self%initialise()
        if(present(iter)) then
            stkname = basename(CWD_GLOB)//"/"//trim(CAVGS_ITER_FBODY)//int2str_pad(iter,3)//trim(params_glob%ext)
        else
            call spproj%get_cavgs_stk(stkname, ncls, smpd)
        end if
        do i=1, spproj%os_cls2D%get_noris()
            call spproj%os_cls2D%set(i, "stk", trim(adjustl(stkname)))
        end do
        call enable_splflags(spproj%os_cls2D, self%starfile%clusters2D%flags)
        call self%export_stardata(spproj, self%starfile%clusters2D%flags, spproj%os_cls2D, "clusters", .false.)
    end subroutine export_cls2D

    subroutine export_ptcls2D(self, cline, spproj)
        class(starproject), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        class(sp_project),  intent(inout) :: spproj
        integer :: i
        if( L_VERBOSE_GLOB ) VERBOSE_OUTPUT = .true.
        self%starfile%filename = "particles2D.star"
        if( VERBOSE_OUTPUT )then
            write(logfhandle,*) ''
            write(logfhandle,*) char(9), 'exporting particles2D to ' // trim(adjustl(self%starfile%filename))
            write(logfhandle,*) ''
        endif
        if(.not. self%starfile%initialised) call self%initialise()
        call enable_splflags(spproj%os_optics, self%starfile%optics%flags)
        call enable_splflags(spproj%os_ptcl2D, self%starfile%particles2D%flags)
        call center_boxes(spproj, spproj%os_ptcl2D)
        call self%export_stardata(spproj, self%starfile%optics%flags, spproj%os_optics, "optics")
        call self%export_stardata(spproj, self%starfile%particles2D%flags, spproj%os_ptcl2D, "particles", .true.)
    end subroutine export_ptcls2D

    subroutine export_ptcls3D(self, cline, spproj)
        class(starproject), intent(inout) :: self
        class(cmdline),     intent(inout) :: cline
        class(sp_project) :: spproj
        integer           :: i
        if( L_VERBOSE_GLOB ) VERBOSE_OUTPUT = .true.
        self%starfile%filename = "particles3D.star"
        if( VERBOSE_OUTPUT )then
            write(logfhandle,*) ''
            write(logfhandle,*) char(9), 'exporting particles3D to ' // trim(adjustl(self%starfile%filename))
            write(logfhandle,*) ''
        endif
        if(.not. self%starfile%initialised) call self%initialise()
        call enable_splflags(spproj%os_optics, self%starfile%optics%flags)
        call enable_splflags(spproj%os_ptcl3D, self%starfile%particles3D%flags)
        call self%export_stardata(spproj, self%starfile%optics%flags, spproj%os_optics, "optics")
        call self%export_stardata(spproj, self%starfile%particles3D%flags, spproj%os_ptcl3D, "particles", .true.)
    end subroutine export_ptcls3D

    subroutine export_stardata(self, spproj, flags, sporis, blockname, mapstks)
        class(starproject), intent(inout) :: self
        class(sp_project),  intent(inout) :: spproj
        type(star_flag),    intent(in)    :: flags(:)
        class(oris),        intent(inout) :: sporis
        character(len=*),   intent(in)    :: blockname
        logical, optional,  intent(in)    :: mapstks
        character(len=:), allocatable :: stkname
        character(len=XLONGSTRLEN)    :: cwd
        character(len=LONGSTRLEN)     :: strtmp
        integer                       :: iflag, iori, ok, stkindex, fhandle
        real                          :: rval, stkind
        logical                       :: ex
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
            if(flags(iflag)%present) then
                write(fhandle, *) "_" // flags(iflag)%rlnflag
            end if
        end do
        if( present(mapstks) )then
            if( mapstks )then
                write(fhandle, *) "_rlnImageName"
                write(fhandle, *) "_rlnMicrographName"
            endif
        end if
        do iori=1, sporis%get_noris()
            if(sporis%get_state(iori) > 0) then
                do iflag=1, size(flags)
                    if(flags(iflag)%present) then
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
    end subroutine export_stardata

    ! tilt

    subroutine assign_initial_tiltgroups(self)
        class(starproject),      intent(inout) :: self
        character(len=LONGSTRLEN), allocatable :: epugroups(:)
        character(len=LONGSTRLEN)              :: eputilt
        integer :: i, epugroupid
        if( VERBOSE_OUTPUT ) write(logfhandle,*) ''
        if(index(self%tiltinfo(1)%basename, 'FoilHole') == 0) then
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
    end subroutine assign_initial_tiltgroups

    subroutine assign_xml_tiltinfo(self, xmldir)
        class(starproject), target, intent(inout) :: self
        character(len=*),           intent(in)    :: xmldir
        character(len=LONGSTRLEN) :: eputilt
        type(Node), pointer :: xmldoc, beamtiltnode, beamtiltnodex, beamtiltnodey
        integer :: i, j
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
        integer, allocatable :: populations(:), labels(:), tiltinfopos(:)
        real,    allocatable :: tilts(:,:), centroids(:,:)
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
    end subroutine cluster_tiltinfo

    ! optics

    subroutine assign_optics(self, cline, spproj)
        class(starproject),    intent(inout)   :: self
        class(cmdline),         intent(inout)           :: cline
        class(sp_project)                               :: spproj
        integer                                         :: i, j, element
        integer                                         :: ptclstkid
        character(len=LONGSTRLEN)                       :: ogname
        character(len=XLONGSTRLEN)                      :: cwd
        allocate(self%tiltinfo(0))
        if(spproj%os_mic%get_noris() > 0) then
            call self%get_image_basename(spproj%os_mic, 'intg')
        end if
        if(spproj%os_stk%get_noris() > 0) then
            call self%get_image_basename(spproj%os_stk, 'stk')
        end if
        call self%assign_initial_tiltgroups()
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

        deallocate(self%tiltinfo)
        if( VERBOSE_OUTPUT ) write(logfhandle,*) ''
    end subroutine assign_optics

    subroutine populate_opticsmap(self, opticsoris)
        class(starproject), intent(inout) :: self
        class(oris),        intent(inout) :: opticsoris
        integer :: i
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

end module simple_starproject
