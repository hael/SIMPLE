module simple_starproject

include 'simple_lib.f08'

!$ use omp_lib

use simple_sp_project,          only: sp_project
use simple_cmdline,             only: cmdline
use simple_oris,                only: oris
use simple_ori,                 only: ori
use simple_rnd

use CPlot2D_wrapper_module
use FoX_dom

implicit none

public :: star_project

#include "simple_local_flags.inc"

type stk_map
    
    character(LEN=2056) :: stkpath
    integer             :: stkmax


end type stk_map

type star_flag

    character(LEN=64) :: rlnflag
    character(LEN=64) :: rlnflag2   = ''
    character(LEN=64) :: splflag
    character(LEN=64) :: splflag2   = ''
    integer           :: ind        = 0
    logical           :: present    = .false.
    logical           :: string     = .false.
    logical           :: int        = .false.
    real              :: mult       = 0
    logical           :: imagesplit = .false.
    logical           :: addstk     = .false.

end type star_flag
    
    
type star_data
    
    type(star_flag),      dimension(:), allocatable   :: flags
    integer                                           :: flagscount = 0
    integer                                           :: datastart = 0
    integer                                           :: dataend = 0
        
end type star_data


type star_file

    character(len=1024)                             :: filename
    character(len=2048)                             :: rootdir
    type(star_data)                                 :: optics
    type(star_data)                                 :: stacks
    type(star_data)                                 :: micrographs
    type(star_data)                                 :: particles2D
    type(star_data)                                 :: particles3D
    type(star_data)                                 :: clusters2D
    integer,            dimension(:), allocatable   :: opticsmap
    integer,            dimension(:,:), allocatable :: stkmap ! (stkid : z)
    integer                                         :: stkptclcount
    logical                                         :: initialised = .false.
        
end type star_file


type tilt_info

    character(len=LONGSTRLEN)   :: basename
    integer                     :: initialtiltgroupid = 1
    integer                     :: finaltiltgroupid = 1
    real                        :: tiltx = 0.0
    real                        :: tilty = 0.0
    real                        :: smpd = 0.0
    real                        :: kv = 0.0
    real                        :: fraca = 0.0
    real                        :: cs = 0.0
    real                        :: box = 0.0

end type tilt_info


type star_project 

    type(star_file)                         :: starfile
    type(tilt_info),    allocatable         :: tiltinfo(:)
   ! character(len=LONGSTRLEN), allocatable  :: uniquemics(:)
   ! integer, allocatable                    :: epumap(:)
    
    
contains
  
  procedure :: initialise
  procedure :: enable_rlnflag
  procedure :: enable_splflag
  procedure :: get_rlnflagindex
  procedure :: populate_opticsmap
  procedure :: populate_stkmap
  procedure :: split_dataline
  procedure :: import_stardata
  procedure :: export_stardata
  procedure :: read_starheaders
  procedure :: enable_splflags
  procedure :: import_ptcls2D
  procedure :: import_ptcls3D
  procedure :: import_cls2D
  procedure :: import_mics
  procedure :: export_mics
  procedure :: export_ptcls2D
  procedure :: export_ptcls3D
  procedure :: get_image_basename
  procedure :: assign_initial_tiltgroups
  procedure :: assign_xml_tiltinfo
  procedure :: cluster_tiltinfo
  procedure :: plot_opticsgroups
  procedure :: sort_optics_maxpop
  procedure :: apply_optics_offset
  procedure :: assign_optics
  procedure :: h_clust
  
end type star_project

contains

  subroutine initialise(self)
  
    class(star_project),  intent(inout) :: self
    
    ! assign optics flags
    if (.not. allocated(self%starfile%optics%flags)) allocate(self%starfile%optics%flags(0))
    
    self%starfile%optics%flags = [self%starfile%optics%flags, star_flag(rlnflag="rlnVoltage", splflag="kv")]
    self%starfile%optics%flags = [self%starfile%optics%flags, star_flag(rlnflag="rlnOpticsGroup", splflag="ogid", int=.true.)]
    self%starfile%optics%flags = [self%starfile%optics%flags, star_flag(rlnflag="rlnOpticsGroupName", splflag="ogname", string=.true.)]
    self%starfile%optics%flags = [self%starfile%optics%flags, star_flag(rlnflag="rlnAmplitudeContrast", splflag="fraca")]
    self%starfile%optics%flags = [self%starfile%optics%flags, star_flag(rlnflag="rlnSphericalAberration", splflag="cs")]
    self%starfile%optics%flags = [self%starfile%optics%flags, star_flag(rlnflag="rlnMicrographPixelSize", splflag="smpd")]
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
    self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnCtfImage", splflag="ctfeps", string=.true.)]
    self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnCtfPowerSpectrum", splflag="ctfjpg", string=.true.)]
    self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnMicrographCoordinates", splflag="boxfile", string=.true.)]
    
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
    
  end subroutine initialise

  
  subroutine enable_rlnflag(self, rlnflag, flags, flagindex)
  
    class(star_project),    intent(inout)           :: self
    class(star_flag),       dimension(:)            :: flags
    character(LEN=64)                               :: rlnflag
    integer                                         :: flagindex, i
    
    do i = 1, size(flags)
        if(index(flags(i)%rlnflag, trim(adjustl(rlnflag))) > 0 .AND. len_trim(flags(i)%splflag) > 0) then
            flags(i)%present = .true.
            flags(i)%ind = flagindex
            write(logfhandle,*) char(9), char(9), "mapping ", flags(i)%rlnflag, " => ", trim(adjustl(flags(i)%splflag))
            exit
        end if
         
    end do
  
  end subroutine enable_rlnflag
  
  
  subroutine enable_splflag(self, splflag, flags)
  
    class(star_project),    intent(inout)           :: self
    class(star_flag),       dimension(:)            :: flags
    character(LEN=*)                                :: splflag
    integer                                         :: i
    
    do i = 1, size(flags)
        if((index(flags(i)%splflag, trim(adjustl(splflag))) > 0 .OR. index(flags(i)%splflag2, trim(adjustl(splflag))) > 0) .AND. len_trim(flags(i)%rlnflag) > 0) then
            flags(i)%present = .true.
            write(logfhandle,*) char(9), char(9), "mapping ", flags(i)%splflag, " => ", trim(adjustl(flags(i)%rlnflag))
            exit
        end if
         
    end do
  
  end subroutine enable_splflag
  
  
  subroutine get_rlnflagindex(self, rlnflag, flags, flagindex)
  
    class(star_project),    intent(inout)           :: self
    class(star_flag),       dimension(:)            :: flags
    character(LEN=*)                               :: rlnflag
    integer,                intent(inout)           :: flagindex
    integer                                         :: i
    
    flagindex = 0
    
    do i = 1, size(flags)
      
        if(index(flags(i)%rlnflag, trim(adjustl(rlnflag))) > 0) then
        
            flagindex = flags(i)%ind
            exit
            
        end if
        
    end do
    
  end subroutine get_rlnflagindex
  
  
  subroutine populate_opticsmap(self, opticsoris)
      
    class(star_project),    intent(inout)               :: self
    class(oris),            intent(inout)               :: opticsoris
    integer                                             :: i
    
    call self%import_stardata(self%starfile%optics, opticsoris, .false.)
    
    allocate(self%starfile%opticsmap(opticsoris%get_noris()))

    do i = 1, opticsoris%get_noris()
        
        self%starfile%opticsmap(int(opticsoris%get(i, "ogid"))) = i

    end do
    
    
  end subroutine populate_opticsmap
  
  
  subroutine populate_stkmap(self, stkoris, opticsoris)
      
    class(star_project),    intent(inout)               :: self
    class(oris),            intent(inout)               :: stkoris
    class(oris),            intent(inout)               :: opticsoris
    type(oris)                                          :: stktmp
    type(ori)                                           :: oritmpin, oritmpout
    character(len=2048),    dimension(:),allocatable    :: stknames
    integer,                dimension(:),allocatable    :: stkzmax, stkoriids
    integer                                             :: i, j, top, fromp
    logical                                             :: newstack
    
    allocate(stknames(0))
    allocate(stkzmax(0))
    allocate(stkoriids(0))
    
    call self%import_stardata(self%starfile%stacks, stktmp, .false., opticsoris)
   
    allocate(self%starfile%stkmap(stktmp%get_noris(), 2))
    
    do i = 1, stktmp%get_noris()
        
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
    
    do i = 1, size(stknames)
    
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
    
    do i = 1, stktmp%get_noris()
    
        self%starfile%stkmap(i, 2) =  self%starfile%stkmap(i, 2) + int(stkoris%get(self%starfile%stkmap(i, 1), "fromp")) - 1
        
    end do
    
    deallocate(stknames)
    deallocate(stkzmax)
    deallocate(stkoriids)
    
  end subroutine populate_stkmap
  
  
  subroutine split_dataline(self, line, splitline)
  
    class(star_project),    intent(inout)                           :: self
    character(len=2048)                                 :: line
    character(len=2048),     dimension(:),intent(inout)                           :: splitline
    integer :: iend, istart, flagid
     
    flagid = 1
    iend = 1
    istart = 1

    do while (flagid <= size(splitline))
    
        do while (line(istart:istart + 1) .eq. " ") 
                        
            istart = istart + 1
                          
        end do
                        
        iend = index(line(istart + 1:), ' ') + istart

        splitline(flagid) = trim(adjustl(line(istart:iend)))
        
        istart = iend + 1
        
        flagid = flagid + 1
    end do

  end subroutine split_dataline
  
  
  subroutine import_stardata(self, stardata, sporis, isptcl, spoptics)
  
    class(star_project),    intent(inout)               :: self
    class(star_data),       intent(inout)               :: stardata
    class(oris),            intent(inout)               :: sporis
    class(oris),            intent(inout),optional      :: spoptics
    type(ori)                                           :: opticsori, spori
    logical                                             :: isptcl
    character(len=XLONGSTRLEN)                          :: cwd
    character(len=2048)                                 :: line, entrystr
    character(len=LONGSTRLEN)                           :: abspath
    character(len=2048)                                 :: splitimage
    character(len=2048),    dimension(:), allocatable   :: splitline
    real                                                :: rval
    integer                                             :: ios, flagsindex, lineindex, ival, ogid, ogmapid, projindex
    
    call simple_getcwd(cwd)
     
    allocate(splitline(stardata%flagscount))
    
    if(isptcl) then
    
        call sporis%new(self%starfile%stkptclcount, .true.)
        
    else
    
        call sporis%new(stardata%dataend - stardata%datastart, .false.)
        
    end if
    
    open(1, file = trim(adjustl(self%starfile%filename)), status = 'old') 
    
    do lineindex = 1, stardata%datastart - 1
    
        read(1, '(A)', iostat=ios) line
        
    end do
    
    do lineindex = 1, stardata%dataend - stardata%datastart
    
        if(isptcl) then
        
            projindex = self%starfile%stkmap(lineindex, 2)
        
        else
        
            projindex = lineindex
        
        end if
        
        
        read(1, '(A)', iostat=ios) line
        
        call self%split_dataline(line, splitline)
        
        do flagsindex = 1, size(stardata%flags)
        
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
                            
                            call make_relativepath(cwd, stemname(stemname(trim(adjustl(self%starfile%rootdir)))) // "/" // trim(adjustl(entrystr)), abspath)
                            call sporis%set(projindex, stardata%flags(flagsindex)%splflag, trim(adjustl(abspath)))
                        
                        else ! other
               
                            call make_relativepath(cwd, stemname(trim(adjustl(self%starfile%rootdir))) // "/" // trim(adjustl(entrystr)), abspath)
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
   
   close(1)
  
  end subroutine import_stardata
  
  
  subroutine export_stardata(self, spproj, flags, sporis, blockname, mapstks)
  
    class(star_project),    intent(inout)               :: self
    class(star_flag),       dimension(:)                :: flags
    class(oris),            intent(inout)               :: sporis
    class(sp_project),      intent(inout)               :: spproj
    character(len=*)                                    :: blockname
    character(len=:),       allocatable                 :: stkname
    character(len=XLONGSTRLEN)                          :: cwd
    character(len=LONGSTRLEN)                           :: strtmp
    integer                                             :: iflag, iori, ok, stkindex
    real                                                :: rval
    logical,                optional                    :: mapstks
    logical                                             :: ex
   
    call simple_getcwd(cwd)
    
    inquire(file=trim(adjustl(self%starfile%filename)), exist=ex)

    if (ex) then
        open(1,file=trim(adjustl(self%starfile%filename)),position='append',iostat=ok)
    else
        open(1,file=trim(adjustl(self%starfile%filename)),status='new',iostat=ok)
    endif

    write(1, *) ""
    write(1, *) "# version 30001"
    write(1, *) ""
    write(1, *) "data_" // trim(adjustl(blockname))
    write(1, *) ""
    write(1, *) "loop_"
    
    do iflag=1, size(flags)
    
        if(flags(iflag)%present) then
        
            write(1, *) "_" // flags(iflag)%rlnflag
            
        end if
        
    end do
    
    if(present(mapstks) .and. mapstks) then
        
        write(1, *) "_rlnImageName"
        
    end if
    
    
    do iori=1, sporis%get_noris()
        
        if(sporis%get_state(iori) > 0) then
        
            do iflag=1, size(flags)
        
                if(flags(iflag)%present) then
                        
                    if(flags(iflag)%string) then
                        
                        strtmp = trim(adjustl(sporis%get_static(iori, trim(adjustl(flags(iflag)%splflag)))))
                        
                        if(index(strtmp, '../') == 1) then ! Relative path. Adjust to base directory
                            
                            write(1, "(A)", advance="no")  trim(adjustl(strtmp(4:))) // " "

                        else
                        
                            write(1, "(A)", advance="no") trim(adjustl(sporis%get_static(iori, trim(adjustl(flags(iflag)%splflag))))) // " "
                        
                        end if
                        
                                
                    elseif(flags(iflag)%int) then
                        
                        rval = sporis%get(iori, trim(adjustl(flags(iflag)%splflag)))
                          
                        if(flags(iflag)%mult > 0) then
                                
                            rval = rval / flags(iflag)%mult
                                    
                        end if
                        
                        write(1, "(I6,A)", advance="no") int(rval), " "
                                
                    else

                        rval = sporis%get(iori, trim(adjustl(flags(iflag)%splflag)))
                            
                        if(flags(iflag)%mult > 0) then
                                
                            rval = rval / flags(iflag)%mult
                                    
                        end if
                          
                        write(1, "(F12.4,A)", advance="no") rval, " "
                                
                    end if
                    
                end if
                
            end do
            
            if(present(mapstks) .and. mapstks) then
                
                call spproj%get_stkname_and_ind('ptcl2D', iori, stkname, stkindex)
                
                if(index(stkname, '../') == 1) then ! Relative path. Adjust to base directory
                
                    write(1, "(I7,A,A,A)", advance="no") int(stkindex), '@', trim(adjustl(stkname(4:))), ' '
                    
                else
                
                    write(1, "(I7,A,A,A)", advance="no") int(stkindex), '@', trim(adjustl(stkname)), ' '
                
                end if
                
            end if
            
            write(1,*) ""
        
        end if
        
    end do
    
    close(1)
    
    
  end subroutine export_stardata
 
  
  subroutine read_starheaders(self)
  
    class(star_project),    target, intent(inout)       :: self
    character(len=2048)                                 :: line, blockname, flagname, datastring
    logical                                             :: flagsopen, dataopen
    integer                                             :: ios, delimindex,flagindex, istart, iend, flagid, lineindex
    
    flagsopen = .false.
    dataopen = .false.
    lineindex = 1
    
    open(1, file = trim(adjustl(self%starfile%filename)), status = 'old')  
    
    do
        read(1, '(A)', iostat=ios) line
      
        if(ios /= 0) then
            exit
        end if
      
        if(len_trim(line) > 0) then
            
            line = trim(adjustl(line))
            
            if (index(line, "data_") > 0) then
                ! start data block
                
                flagindex = 0
                flagsopen = .true.
                dataopen = .false.
                
                delimindex = index(line, "data_ ") + 6
                blockname = line(delimindex:)

            else if(flagsopen .AND. index(line, "_rln") > 0) then
            
                delimindex = index(line, "#")
                
                if(delimindex > 0) then
                
                    flagname = line(2:delimindex - 1)
                    
                else
                
                    flagname = trim(adjustl(line(2:)))
                
                end if
                
                flagindex = flagindex + 1
                
                select case (blockname)
                
                    case ("optics")
                    
                        call self%enable_rlnflag(flagname, self%starfile%optics%flags, flagindex)
                        self%starfile%optics%flagscount = self%starfile%optics%flagscount + 1
                        
                    case ("micrographs")
                    
                        call self%enable_rlnflag(flagname, self%starfile%micrographs%flags, flagindex)
                        self%starfile%micrographs%flagscount = self%starfile%micrographs%flagscount + 1
                        
                    case ("particles")
                    
                        call self%enable_rlnflag(flagname, self%starfile%particles2D%flags, flagindex)
                        self%starfile%particles2D%flagscount = self%starfile%particles2D%flagscount + 1
                        
                        call self%enable_rlnflag(flagname, self%starfile%particles3D%flags, flagindex)
                        self%starfile%particles3D%flagscount = self%starfile%particles3D%flagscount + 1
                        
                        call self%enable_rlnflag(flagname, self%starfile%stacks%flags, flagindex)
                        self%starfile%stacks%flagscount = self%starfile%stacks%flagscount + 1
                        
                    case ("model_classes")
                    
                        call self%enable_rlnflag(flagname, self%starfile%clusters2D%flags, flagindex)
                        self%starfile%clusters2D%flagscount = self%starfile%clusters2D%flagscount + 1 
                    
                end select
                
                
            else if(flagsopen .AND. index(line, "loop_") == 0) then
                
                dataopen = .true.
                flagsopen = .false.
                
                select case (blockname)
                
                    case ("optics")
                    
                        self%starfile%optics%datastart = lineindex
                        
                    case ("micrographs")
                    
                        self%starfile%micrographs%datastart = lineindex
                        
                    case ("particles")
                    
                        self%starfile%particles2D%datastart = lineindex
                        self%starfile%particles3D%datastart = lineindex
                        self%starfile%stacks%datastart = lineindex 
                    
                    case ("model_classes")  
                        
                        self%starfile%clusters2D%datastart = lineindex
                        
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
    
    close(1)
    
  end subroutine read_starheaders


  subroutine enable_splflags(self, sporis, flags)
  
    class(star_project),        intent(inout)   :: self
    class(oris),                intent(inout)   :: sporis
    class(star_flag),           dimension(:)    :: flags
    type(ori)                                   :: testori
    character(len=XLONGSTRLEN), allocatable     :: keys(:)
    integer                                     :: iori, ikey
    
    ! find 1st non state 0 ori

    do iori = 1, sporis%get_noris()
       
        if(sporis%get_state(iori) > 0) then
             
            call sporis%get_ori(iori, testori)
            
            keys = testori%get_keys()
            
            do ikey=1, size(keys)
            
                call self%enable_splflag(trim(adjustl(keys(ikey))), flags)
                
            end do
            
            exit
        
        end if
        
    end do
    
  end subroutine enable_splflags
  
  
  subroutine import_ptcls2D(self, cline, spproj, filename)
    
    class(star_project),    target, intent(inout)   :: self
    class(cmdline),         intent(inout)           :: cline
    class(sp_project)                               :: spproj
    character(len=*)                                :: filename
    integer                                         :: i
    
    write(logfhandle,*) ''
    write(logfhandle,*) char(9), 'importing ' // filename // " to ptcls2D"
    write(logfhandle,*) 
    
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
    
    deallocate(self%starfile%opticsmap)
    deallocate(self%starfile%stkmap)
    
  end subroutine import_ptcls2D
  
  
  subroutine import_ptcls3D(self, cline, spproj, filename)
    
    class(star_project),    target, intent(inout)   :: self
    class(cmdline),         intent(inout)           :: cline
    class(sp_project)                               :: spproj
    character(len=*)                                :: filename
    integer                                         :: i
    
    write(logfhandle,*) ''
    write(logfhandle,*) char(9), 'importing ' // filename // " to ptcls3D"
    write(logfhandle,*) 
    
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
  
  
  subroutine import_cls2D(self, cline, spproj, filename)
    
    class(star_project),    target, intent(inout)   :: self
    class(cmdline),         intent(inout)           :: cline
    class(sp_project)                               :: spproj
    character(len=*)                                :: filename
    integer                                         :: i
    
    write(logfhandle,*) ''
    write(logfhandle,*) char(9), 'importing ' // filename // " to cls2D"
    write(logfhandle,*) ''
    
    if(.not. self%starfile%initialised) call self%initialise()
    
    self%starfile%filename = filename
    self%starfile%rootdir = cline%get_carg("import_dir")
    
    call self%read_starheaders()

    call self%import_stardata(self%starfile%clusters2D, spproj%os_cls2D, .false.)
    
  end subroutine import_cls2D
  
  
  subroutine import_mics(self, cline, spproj, filename)
    
    class(star_project),    target, intent(inout)   :: self
    class(cmdline),         intent(inout)           :: cline
    class(sp_project)                               :: spproj
    character(len=*)                                :: filename
    integer                                         :: i
    
    write(logfhandle,*) ''
    write(logfhandle,*) char(9), 'importing3 ' // filename // " to mics"
    write(logfhandle,*) ''
    
    if(.not. self%starfile%initialised) call self%initialise()
    
    self%starfile%filename = filename
    
    self%starfile%rootdir = cline%get_carg("import_dir")
    
    call self%read_starheaders()
    
    call self%populate_opticsmap(spproj%os_optics)

    call self%import_stardata(self%starfile%micrographs, spproj%os_mic, .false., spproj%os_optics)
    
    do i=1, spproj%os_mic%get_noris()
        
        call spproj%os_mic%set(i, "imgkind", "mic")
        call spproj%os_mic%set(i, "ctf", "yes")
     
    end do
    
    deallocate(self%starfile%opticsmap)
    
  end subroutine import_mics
  
  
  subroutine export_mics(self, cline, spproj)
    
    class(star_project),    target, intent(inout)   :: self
    class(cmdline),         intent(inout)           :: cline
    class(sp_project)                               :: spproj
    integer                                         :: i
    
    self%starfile%filename = "micrographs.star"
    self%starfile%rootdir = cline%get_carg("import_dir")
    write(logfhandle,*) ''
    write(logfhandle,*) char(9), 'exporting micrographs to ' // trim(adjustl(self%starfile%filename))
    write(logfhandle,*) ''
    
    if(.not. self%starfile%initialised) call self%initialise()
    
    call self%enable_splflags(spproj%os_optics, self%starfile%optics%flags)
    call self%enable_splflags(spproj%os_mic, self%starfile%micrographs%flags)
    
    call self%export_stardata(spproj, self%starfile%optics%flags, spproj%os_optics, "optics")
    call self%export_stardata(spproj, self%starfile%micrographs%flags, spproj%os_mic, "micrographs")
    
  end subroutine export_mics
  
  
  subroutine export_ptcls2D(self, cline, spproj)
    
    class(star_project),    target, intent(inout)   :: self
    class(cmdline),         intent(inout)           :: cline
    class(sp_project)                               :: spproj
    integer                                         :: i
    
    self%starfile%filename = "particles2D.star"
    !self%starfile%rootdir = cline%get_carg("import_dir")
    
    write(logfhandle,*) ''
    write(logfhandle,*) char(9), 'exporting particles2D to ' // trim(adjustl(self%starfile%filename))
    write(logfhandle,*) ''
    
    if(.not. self%starfile%initialised) call self%initialise()
    
    call self%enable_splflags(spproj%os_optics, self%starfile%optics%flags)
    call self%enable_splflags(spproj%os_ptcl2D, self%starfile%particles2D%flags)
    
    call self%export_stardata(spproj, self%starfile%optics%flags, spproj%os_optics, "optics")
    call self%export_stardata(spproj, self%starfile%particles2D%flags, spproj%os_ptcl2D, "particles", .true.)
    
  end subroutine export_ptcls2D
  
  
  subroutine export_ptcls3D(self, cline, spproj)
    
    class(star_project),    target, intent(inout)   :: self
    class(cmdline),         intent(inout)           :: cline
    class(sp_project)                               :: spproj
    integer                                         :: i
    
    self%starfile%filename = "particles3D.star"
    !self%starfile%rootdir = cline%get_carg("import_dir")
    
    write(logfhandle,*) ''
    write(logfhandle,*) char(9), 'exporting particles3D to ' // trim(adjustl(self%starfile%filename))
    write(logfhandle,*) ''
    
    if(.not. self%starfile%initialised) call self%initialise()
    
    call self%enable_splflags(spproj%os_optics, self%starfile%optics%flags)
    call self%enable_splflags(spproj%os_ptcl3D, self%starfile%particles3D%flags)
    
    call self%export_stardata(spproj, self%starfile%optics%flags, spproj%os_optics, "optics")
    call self%export_stardata(spproj, self%starfile%particles3D%flags, spproj%os_ptcl3D, "particles", .true.)
    
  end subroutine export_ptcls3D
  
  
  subroutine get_image_basename(self, sporis, splflag)
  
    class(star_project),    target, intent(inout)   :: self
    class(oris),                    intent(inout)   :: sporis
    type(tilt_info)                                 :: tiltinfo
    character(len=*)                                :: splflag
    character(len=XLONGSTRLEN)                      :: path
    character(len=LONGSTRLEN)                       :: filename
    integer                                         :: i, tiltind
    
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
            tiltinfo%smpd = sporis%get(i, "smpd")
            tiltinfo%kv = sporis%get(i, "kv")
            tiltinfo%fraca = sporis%get(i, "fraca")
            tiltinfo%cs = sporis%get(i, "cs")
            tiltinfo%box = sporis%get(i, "box")

            self%tiltinfo = [self%tiltinfo, tiltinfo]
        
        else
            
            if(self%tiltinfo(tiltind)%box < sporis%get(i, "box")) then
            
                self%tiltinfo(tiltind)%box = sporis%get(i, "box")
            
            end if
            
        end if
        
    
    end do    

  end subroutine get_image_basename
  
  
  subroutine assign_initial_tiltgroups(self)
  
    class(star_project),    target, intent(inout)   :: self
    character(len=LONGSTRLEN), allocatable          :: epugroups(:)
    character(len=LONGSTRLEN)                       :: eputilt
    integer                                         :: i
    integer                                         :: epugroupid
    
    write(logfhandle,*) ''
    
    if(index(self%tiltinfo(1)%basename, 'FoilHole') == 0) then
     
        write(logfhandle,*) char(9), 'no EPU filenames detected. assigning single initial beamtilt group'
        ! already assigned initialtiltgroupid=1 in initialisation
        
     else
     
         write(logfhandle,"(A,A)", advance="no") char(9), 'EPU filenames detected. using these to assign initial beamtilt groups ... '
         
         allocate(epugroups(0))
         
         epugroupid = 1
         
         do i=1, size(self%tiltinfo)
         
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
         
         write(logfhandle,"(A,I4,A)") "assigned ", size(epugroups), " initial groups"
         
         deallocate(epugroups)
         
     end if
     
     write(logfhandle,*) ''
    
  
  end subroutine assign_initial_tiltgroups
  
  
  subroutine assign_xml_tiltinfo(self, xmldir)
  
    class(star_project),    target, intent(inout)   :: self
    character(len=LONGSTRLEN)                       :: eputilt
    character(len=*)                                :: xmldir
    integer                                         :: i, j
    type(Node),             pointer                 :: xmldoc, beamtiltnode, beamtiltnodex, beamtiltnodey
    
    write(logfhandle,*) ''
    write(logfhandle,*) char(9), "reading tilt info from metadata ... "
    
    do i = 1, size(self%tiltinfo)
        
        if(file_exists(xmldir // '/' // trim(adjustl(self%tiltinfo(i)%basename)) // '.xml')) then
        
            xmldoc => parseFile(xmldir // '/' // trim(adjustl(self%tiltinfo(i)%basename)) // '.xml')
            
            beamtiltnode => item(getElementsByTagname(xmldoc, "BeamShift"), 0)
            beamtiltnodex => item(getElementsByTagname(beamtiltnode, "a:_x"), 0)
            beamtiltnodey => item(getElementsByTagname(beamtiltnode, "a:_y"), 0)
            self%tiltinfo(i)%tiltx = str2real(getTextContent(beamtiltnodex))
            self%tiltinfo(i)%tilty = str2real(getTextContent(beamtiltnodey))
            
            call destroy(xmldoc)  
               
        else
            
            write(logfhandle, *) char(9), char(9), xmldir // '/' // trim(adjustl(self%tiltinfo(i)%basename)) // '.xml does not exist. Ignoring'
                
        end if
    
    end do
    
    write(logfhandle,*) ''

  end subroutine assign_xml_tiltinfo
  
  
  subroutine cluster_tiltinfo(self, threshold)
  
    class(star_project),    target, intent(inout)   :: self
    integer                                         :: i, j, k, tiltcount, groupcount, matchcount
    real,                   allocatable             :: tilts(:,:)
    real,                   allocatable             :: centroids(:,:)
    integer,                allocatable             :: populations(:)
    integer,                allocatable             :: labels(:)
    integer,                allocatable             :: tiltinfopos(:)
    real                                            :: threshold
    
    write(logfhandle,*) ''
    write(logfhandle,*) char(9), "clustering initial beamtilt groups using tilt info and a threshold of ", threshold, " ... "
    
    groupcount = 0
            

    do i = 1, maxval(self%tiltinfo%initialtiltgroupid)
    
        write(logfhandle,*) char(9), char(9), "clustering initial beamtilt group ", i, " ... "

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

        call self%h_clust(tilts, threshold, labels, centroids, populations)

        do j = 1, size(tiltinfopos)
            
            self%tiltinfo(tiltinfopos(j))%finaltiltgroupid = groupcount + labels(j)
            
        end do
        
        groupcount = groupcount + size(populations)
        
        deallocate(tilts)
        deallocate(labels)
        deallocate(tiltinfopos)
        write(logfhandle,*) ''
        
    end do

  end subroutine cluster_tiltinfo
  
  
  subroutine plot_opticsgroups(self, fname_eps)
  
    class(star_project),    target, intent(inout)   :: self
    type(str4arr)                                   :: title
    type(CPlot2D_type)                              :: plot2D
    type(CDataSet_type)                             :: dataSet
    integer                                         :: i, j
    character(len=*)                       :: fname_eps
    
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
    
    write(logfhandle,*) char(9), char(9), "wrote ", fname_eps
    
  end subroutine plot_opticsgroups
  
  
  subroutine sort_optics_maxpop(self, maxpop)
  
    class(star_project),    target, intent(inout)   :: self
    integer                                         :: maxpop
    integer,                allocatable             :: nextidpop(:,:)
    integer                                         :: i
    
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
  
    class(star_project),    target, intent(inout)   :: self
    integer                                         :: offset
    integer                                         :: i
    
    do i = 1, size(self%tiltinfo)
    
        self%tiltinfo(i)%finaltiltgroupid = self%tiltinfo(i)%finaltiltgroupid + offset
    
    end do
    
  
  end subroutine apply_optics_offset
  
  
  
  subroutine assign_optics(self, cline, spproj)
    
    class(star_project),    target, intent(inout)   :: self
    class(cmdline),         intent(inout)           :: cline
    class(sp_project)                               :: spproj
    integer                                         :: i, j, element
    !integer                                         :: tiltinfosz
    integer                                         :: ptclstkid
    character(len=LONGSTRLEN)                       :: ogname
    character(len=XLONGSTRLEN)                      :: cwd
        
    !tiltinfosz = spproj%os_mic%get_noris() + spproj%os_stk%get_noris()
    
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
    
        write(logfhandle,*) ''
        write(logfhandle,*) char(9), "plotting beamtilt groups (prior to splitting on maxpop) ... "
        call self%plot_opticsgroups("optics_groups_pre_maxpop.eps")
        
        call self%sort_optics_maxpop(int(cline%get_rarg("maxpop")))
        
        write(logfhandle,*) ''
        write(logfhandle,*) char(9), "plotting beamtilt groups (after splitting on maxpop) ... "
        call self%plot_opticsgroups("optics_groups_post_maxpop.eps") 
        
    else
    
        write(logfhandle,*) ''
        write(logfhandle,*) char(9), "plotting beamtilt groups ... "
        call self%plot_opticsgroups("optics_groups.eps")
        
    end if
    
    if(cline%get_rarg("optics_offset") > 0) then
        
        call self%apply_optics_offset(int(cline%get_rarg("optics_offset")))
    
    end if
    

    
    if(spproj%os_mic%get_noris() > 0) then
    
        write(logfhandle,*) ''
        write(logfhandle,*) char(9), "updating micrographs in project file with updated optics groups ... "
        
        do i = 1, spproj%os_mic%get_noris()
            
            element = findloc(self%tiltinfo%basename, trim(adjustl(spproj%os_mic%get_static(i, "bsname"))), 1)
            
            if(element > 0) then
            
                call spproj%os_mic%set(i, 'ogid', real(self%tiltinfo(element)%finaltiltgroupid))
                
            else
            
                THROW_HARD('failed to locate optics info for mic basename' // trim(adjustl(spproj%os_mic%get_static(i, "bsname"))))
            
            end if
            
        end do
        
    end if
    
    if(spproj%os_stk%get_noris() > 0) then
    
        write(logfhandle,*) ''
        write(logfhandle,*) char(9), "updating stacks in project file with updated optics groups ... "
        
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
    
        write(logfhandle,*) ''
        write(logfhandle,*) char(9), "updating particles in project file with updated optics groups ... "
        
        do i = 1, spproj%os_ptcl2D%get_noris()
            
            ptclstkid = int(spproj%os_ptcl2D%get(i, 'stkind'))
            call spproj%os_ptcl2D%set(i, 'ogid', spproj%os_stk%get(ptclstkid, 'ogid'))
            call spproj%os_ptcl3D%set(i, 'ogid', spproj%os_stk%get(ptclstkid, 'ogid'))
            
        end do
        
    end if
    
    write(logfhandle,*) ''
    write(logfhandle,*) char(9), "updating optics groups in project file ... "
    
    call spproj%os_optics%new(maxval(self%tiltinfo%finaltiltgroupid) - minval(self%tiltinfo%finaltiltgroupid) + 1, .false.)
    
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
    
    deallocate(self%tiltinfo)
    
    write(logfhandle,*) ''
    
  end subroutine assign_optics
  
  
! distance threshold based yerarchical clustering
! Source https://www.mathworks.com/help/stats/hierarchical-clustering.html#bq_679x-10
  subroutine h_clust(self,data_in,thresh,labels,centroids,populations)
        class(star_project),  intent(inout) :: self ! unused
        real,    intent(in)  :: data_in(:,:)   ! input data, point coords
        real,    intent(in)  :: thresh         ! threshold for class merging
        integer, intent(out) :: labels(:)      ! labels of the elements in vec
        real,    allocatable, intent(out) :: centroids(:,:) ! centroids of the classes
        integer, allocatable, intent(out) :: populations(:) ! number of elements belonging to each class
        real,    allocatable :: mat(:,:)                    ! pariwise distances matrix
        logical, allocatable :: mask(:), outliers(:)
        integer :: N, i, j, cnt, ncls, filnum, io_stat
        integer :: index(1), loc1(1), loc2(1)
        real    :: d
        
        if( size(data_in, dim = 2) .ne. 2 ) THROW_HARD('Input data should be two dimensional!; h_clust')

        
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
        allocate(mask(N),source = .true. )
        allocate(outliers(N), source = .false.)
        ncls = 0
        do i = 1, N
          if(mask(i)) then ! if it's not already been clustered
            mask(i) = .false.
            ! find the index of the couple
            d = minval(mat(i,:), mask)
            index(:) = minloc(mat(i,:), mask)
            ncls = ncls + 1
            ! assign lavels
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
           if(any(labels== i)) then !there is a class labelled i
              cnt = cnt + 1                  !increasin order cc
              where(labels == i) labels = cnt
            endif
        enddo
        ! 6) recalculate centroids
        deallocate(centroids)
        ncls = maxval(labels) ! the nb of classes is maxval(labels)
        allocate(centroids(ncls,2),   source = 0.)
        allocate(populations(ncls), source = 0 )
        mask = .true. ! reset
        do i = 1, ncls
          populations(i) = count(labels == i) ! counts the nb of elements in the class
          ! find the all the cnt member of the class and update the centroids
          do j = 1, populations(i)
            loc1(:) = minloc(abs(labels-i), mask)
            mask(loc1(1)) = .false. ! do not consider this member of the class anymore
            centroids(i,1) = centroids(i,1)+ data_in(loc1(1),1)
            centroids(i,2) = centroids(i,2)+ data_in(loc1(1),2)
          enddo
          centroids(i,:) = centroids(i,:)/real(populations(i))
        enddo
        ! output generation, Joseph please change it as you wish

        do i = 1, ncls
            write(logfhandle, *) char(9), char(9), char(9), 'Class ', i, 'Population ', populations(i)
        enddo
  end subroutine h_clust
  
end module simple_starproject
