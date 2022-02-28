module simple_starproject

include 'simple_lib.f08'

!$ use omp_lib

use simple_sp_project,          only: sp_project
use simple_cmdline,             only: cmdline
use simple_oris,                only: oris
use simple_ori,                 only: ori

implicit none

public :: star_project

type stk_map
    
    character(LEN=2056) :: stkpath
    integer             :: stkmax


end type stk_map

type star_flag

    character(LEN=64) :: rlnflag
    character(LEN=64) :: splflag
    integer           :: ind        = 0
    logical           :: present    = .false.
    logical           :: string     = .false.
    logical           :: int        = .false.
    real              :: mult       = 0
    logical           :: imagesplit = .false.
    logical           :: addstk     = .false.
    logical           :: appendroot = .false.

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
    integer,            dimension(:), allocatable   :: opticsmap
    integer,            dimension(:,:), allocatable :: stkmap ! (stkid : z)
    integer                                         :: stkptclcount
    logical                                         :: initialised = .false.
        
end type star_file


type star_project 

    type(star_file) :: starfile
    
contains
  
  procedure :: initialise
  procedure :: enable_rlnflag
  procedure :: enable_splflag
  procedure :: get_rlnflagindex
  procedure :: populate_opticsmap
  procedure :: populate_stkmap
  procedure :: split_dataline
  procedure :: import_stardata
  procedure :: read_starheaders
  procedure :: import_ptcls2D
  procedure :: import_ptcls3D
  procedure :: import_mics

  
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
    
    self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnMicrographMovieName", splflag="movie", string=.true., appendroot=.true.)]
    self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnMicrographName", splflag="intg", string=.true., appendroot=.true.)]
    self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnMicrographMetadata", splflag="mc_starfile", string=.true., appendroot=.true.)]
    self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnDefocusU", splflag="dfx", mult=0.0001)]
    self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnDefocusV", splflag="dfy", mult=0.0001)]
    self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnDefocusAngle", splflag="angast")]
    self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnPhaseShift", splflag="phshift")]
    self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnCtfMaxResolution", splflag="ctfres")]
    self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnOpticsGroup", splflag="ogid")]
    self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnImageSizeX", splflag="xdim")]
    self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnImageSizeY", splflag="ydim")]
    self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnImageSizeZ", splflag="nframes")]
    self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnCtfImage", splflag="ctfeps", string=.true., appendroot=.true.)]
    self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnCtfPowerSpectrum", splflag="ctfjpg", string=.true., appendroot=.true.)]
    self%starfile%micrographs%flags = [self%starfile%micrographs%flags, star_flag(rlnflag="rlnMicrographCoordinates", splflag="boxfile", string=.true., appendroot=.true.)]
    
    ! assign stk flags
    if (.not. allocated(self%starfile%stacks%flags)) allocate(self%starfile%stacks%flags(0))
    
    self%starfile%stacks%flags = [self%starfile%stacks%flags, star_flag(rlnflag="rlnImageName", splflag="stk", string=.true., imagesplit=.true., appendroot=.true.)]
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
            print *, char(9), "Mapping ", trim(adjustl(flags(i)%rlnflag)), char(9), char(9), " => ", trim(adjustl(flags(i)%splflag))
            exit
        end if
         
    end do
  
  end subroutine enable_rlnflag
  
  
  subroutine enable_splflag(self, splflag, flags, flagindex)
  
    class(star_project),    intent(inout)           :: self
    class(star_flag),       dimension(:)            :: flags
    character(LEN=64)                               :: splflag
    integer                                         :: flagindex, i
    
    do i = 1, size(flags)
        if(index(flags(i)%splflag, trim(adjustl(splflag))) > 0 .AND. len_trim(flags(i)%rlnflag) > 0) then
            flags(i)%present = .true.
            flags(i)%ind = flagindex
            print *, char(9), "Mapping ", trim(adjustl(flags(i)%splflag)), char(9), char(9), " => ", trim(adjustl(flags(i)%rlnflag))
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
    character(len=2048)                                 :: line, entrystr
    character(len=2048)                                 :: splitimage
    character(len=2048),    dimension(:), allocatable   :: splitline
    real                                                :: rval
    integer                                             :: ios, flagsindex, lineindex, ival, ogid, ogmapid, projindex
    
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
                    
                    entrystr = splitimage(1:index(splitimage, "@") - 1)
                    
                    read(entrystr,*) ival   
                    
                    call sporis%set(projindex, "stkind", real(ival))
                    
                    entrystr = splitimage(index(splitimage, "@") + 1:)
                
                else
                    
                    entrystr = trim(adjustl(splitline(stardata%flags(flagsindex)%ind)))
                
                end if 
                
                if(stardata%flags(flagsindex)%string) then
                
                    if(stardata%flags(flagsindex)%appendroot) then
                        
                        call sporis%set(projindex, stardata%flags(flagsindex)%splflag, trim(adjustl(self%starfile%rootdir)) // "/" // rev_basename(rev_basename(trim(adjustl(entrystr)))))
                        
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
                
                print *, ""
                print *, "Reading ", trim(adjustl(blockname)), " ... " 
                print *, ""

            else if(flagsopen .AND. index(line, "_rln") > 0) then
            
                delimindex = index(line, "#")
                flagname = line(2:delimindex - 1)
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
                         
                end select
                
            end if
            
        end if
        
        lineindex = lineindex + 1
        
    end do
    
    print *, ""
    
    close(1)
    
  end subroutine read_starheaders

  
  subroutine import_ptcls2D(self, cline, spproj, filename)
    
    class(star_project),    target, intent(inout)   :: self
    class(cmdline),         intent(inout)           :: cline
    class(sp_project)                               :: spproj
    character(len=*)                                :: filename
    integer                                         :: i
    
    print *, "Importing " // filename
    
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
    
    print *, "Importing " // filename
    
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
  
  
  subroutine import_mics(self, cline, spproj, filename)
    
    class(star_project),    target, intent(inout)   :: self
    class(cmdline),         intent(inout)           :: cline
    class(sp_project)                               :: spproj
    character(len=*)                                :: filename
    integer                                         :: i
    
    print *, "Importing " // filename
    
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

end module simple_starproject
