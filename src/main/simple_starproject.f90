module simple_starproject

include 'simple_lib.f08'

!$ use omp_lib

use simple_sp_project,          only: sp_project
use simple_cmdline,             only: cmdline
use simple_oris,                only: oris
use simple_ori,                 only: ori

implicit none

private :: flag

public :: star_project

type flag

  character(LEN=64)   :: rlnflag
  character(LEN=64)   :: splflag
  logical             :: enabled

end type


type star_data

  character(len=64)                               :: title
  logical                                         :: loop
  type(flag), dimension(:), allocatable           :: flags
  character(len=1024), dimension(:), allocatable  :: datalines

  contains
    procedure :: kill_stardata

end type star_data


type star_file

  character(len=1024)                                :: filename
  type(star_data), dimension(:), allocatable         :: datablocks

  contains
    procedure :: kill_starfile

end type star_file


type star_project

contains

  procedure :: assign_optics_single
  procedure :: insert_flag
  procedure :: process_flags
  procedure :: process_data
  procedure :: insert_datablock
  procedure :: write_starfile
  procedure :: mics2star
  procedure :: ptcls2D2star
  procedure :: clusters2D2star
  procedure :: create


end type star_project

contains

  subroutine assign_optics_single(self, spproj)

    class(star_project),  intent(inout)           :: self
    class(sp_project)                             :: spproj
    class(oris),      pointer                     :: os_ptr
    integer                                       ::orisn

    call spproj%new_seg_with_ptr(1, "optics", os_ptr)

    call os_ptr%new(1, is_ptcl=.false.)
    call os_ptr%set(1, "ogname", "opticsgroup1")
    call os_ptr%set(1, "ogid", real(1))
    call os_ptr%set(1, "state", real(1))
    call os_ptr%set(1, "ogdim", real(2))

    if(spproj%os_stk%get_noris() > 0) then

      call os_ptr%set(1, "ogsize", spproj%os_stk%get(1,'box'))
      call os_ptr%set(1, "smpd", spproj%os_stk%get(1,'smpd'))
      call os_ptr%set(1, "cs", spproj%os_stk%get(1,'cs'))
      call os_ptr%set(1, "kv", spproj%os_stk%get(1,'kv'))
      call os_ptr%set(1, "fraca", spproj%os_stk%get(1,'fraca'))

    else if(spproj%os_mic%get_noris() > 0) then

      call os_ptr%set(1, "smpd", spproj%os_mic%get(1,'smpd'))
      call os_ptr%set(1, "cs", spproj%os_mic%get(1,'cs'))
      call os_ptr%set(1, "kv", spproj%os_mic%get(1,'kv'))
      call os_ptr%set(1, "fraca", spproj%os_mic%get(1,'fraca'))

    end if

    ! Set all mics ogid=1
    !$omp parallel do default(shared) private(orisn) schedule(static) proc_bind(close)
    do orisn = 1, spproj%os_mic%get_noris()

      call spproj%os_mic%set(orisn, 'ogid', real(1))

    end do
    !$omp end parallel do

    ! Set all ptcls ogid=1
    !$omp parallel do default(shared) private(orisn) schedule(static) proc_bind(close)
    do orisn = 1, spproj%os_ptcl2D%get_noris()

      call spproj%os_ptcl2D%set(orisn, 'ogid', real(1))

    end do

    !$omp end parallel do


  end subroutine assign_optics_single



  subroutine kill_stardata(self)

    class(star_data), intent(inout) :: self

    if(allocated(self%flags)) then

      deallocate(self%flags)

    end if

    if(allocated(self%datalines)) then

      deallocate(self%datalines)

    end if

  end subroutine kill_stardata

  subroutine kill_starfile(self)

    class(star_file), intent(inout) :: self
    integer                         :: blockn

    do blockn = 1, size(self%datablocks)

      call self%datablocks(blockn)%kill_stardata()

    end do

    if(allocated(self%datablocks)) then

      deallocate(self%datablocks)

    end if

  end subroutine kill_starfile

  subroutine insert_flag(self, micrographscorrected, blockid, newrlnflag, newsplflag)

    class(star_project),  intent(inout)                   :: self
    type(star_file), intent(inout)                  :: micrographscorrected
    character(*)                                    :: newrlnflag
    character(*)                                    :: newsplflag
    type(flag)                                      :: newflag
    integer                                         :: blockid

    newflag%rlnflag = newrlnflag
    newflag%splflag = newsplflag
    newflag%enabled = .false.

    if (.not. allocated(micrographscorrected%datablocks(blockid)%flags)) then

      allocate(micrographscorrected%datablocks(blockid)%flags(0))

    end if

    micrographscorrected%datablocks(blockid)%flags = [micrographscorrected%datablocks(blockid)%flags, newflag]

  end subroutine insert_flag


  subroutine insert_datablock(self, starfile, title, loop, blockid)

    class(star_project),  intent(inout)           :: self
    type(star_file), intent(inout)                :: starfile
    type(star_data)                               :: newdata
    logical                                       :: loop
    character(*)                                  :: title
    integer, intent(inout)                        :: blockid

    if (.not. allocated(starfile%datablocks)) then

      allocate(starfile%datablocks(0))

    end if

    newdata%title = trim(adjustl(title))
    newdata%loop = loop

    starfile%datablocks = [starfile%datablocks, newdata]

    blockid = size(starfile%datablocks)

  end subroutine insert_datablock


  subroutine process_flags(self, starfile, blockid, section)

    class(star_project),  intent(inout)   :: self
    type(oris), intent(inout)             :: section
    type(star_file),intent(inout)         :: starfile
    type(chash)                           :: keys
    type(ori)                             :: element
    type(flag)                            :: newflag
    integer                               :: blockid, n, keysn, flagsn

    ! find first non state 0 ori
    do n = 1, section%get_noris()

      call section%get_ori(n, element)

      if(element%get_state() > 0) then

        keys =  element%ori2chash()
        exit

      endif

    end do

    ! set enabled to true on existing flags
    do keysn = 1, keys%size_of()
      write(*,*) keys%get_key(keysn)
      do flagsn = 1, size(starfile%datablocks(blockid)%flags)

        if(starfile%datablocks(blockid)%flags(flagsn)%splflag == keys%get_key(keysn)) then

          starfile%datablocks(blockid)%flags(flagsn)%enabled = .true.

        end if

      end do

    end do

  end subroutine process_flags


  subroutine process_data(self, starfile, blockid, section)

    class(star_project),  intent(inout)   :: self
    type(oris), intent(inout)             :: section
    type(star_file),intent(inout)         :: starfile
    type(ori)                             :: element
    character(len=1024)                   :: line
    real                                  :: value
    integer                               :: blockid, orisn, flagsn

    allocate(starfile%datablocks(blockid)%datalines(section%get_noris()))

    write(*,*) "STARDATA SIZE : ", sizeof(starfile%datablocks(blockid)%datalines)

    !$omp parallel do default(shared) private(orisn,flagsn,element,line,value) schedule(static) proc_bind(close)

    do orisn = 1, section%get_noris()

      call section%get_ori(orisn, element)

      if(element%get_state() > 0) then

        line = ""

        do flagsn = 1, size(starfile%datablocks(blockid)%flags)

          if(starfile%datablocks(blockid)%flags(flagsn)%enabled) then

            ! IF LOOP FALSE, DO DIFFERENT HERE!
            if(element%ischar(trim(adjustl(starfile%datablocks(blockid)%flags(flagsn)%splflag)))) then

              write(line, fmt="(A,2X,A)")trim(adjustl(line)), element%get_static(trim(adjustl(starfile%datablocks(blockid)%flags(flagsn)%splflag)))

            else

              value = element%get(trim(adjustl(starfile%datablocks(blockid)%flags(flagsn)%splflag)))

              if(ceiling(value) == value) then

                write(line, fmt="(A,2X,I0)") trim(adjustl(line)), ceiling(value)

              else

                write(line, fmt="(A,2X,F8.2)") trim(adjustl(line)), value

              end if

            end if

          end if

        end do

        starfile%datablocks(blockid)%datalines(orisn) = line

      end if

    end do

    !$omp end parallel do

   ! call starfile%datablocks(blockid)%kill()


  end subroutine process_data


  subroutine write_starfile(self, starfile)

    class(star_project),  intent(inout)   :: self
    type(star_file),intent(inout)         :: starfile
    integer                               :: blockn, flagsn, linen

    open(1, file = trim(adjustl(starfile%filename)), status = 'new')

    do blockn = 1, size(starfile%datablocks)

      write(1,*) ""
      write(1, fmt="(A,A)") "data_", trim(adjustl(starfile%datablocks(blockn)%title))

      if(starfile%datablocks(blockn)%loop) then

        write(1,*) ""
        write(1, fmt="(A)") "loop_"

        do flagsn = 1, size(starfile%datablocks(blockn)%flags)

          if(starfile%datablocks(blockn)%flags(flagsn)%enabled) then

            write(1, fmt="(A,A,A,I0)") "_", trim(adjustl(starfile%datablocks(blockn)%flags(flagsn)%rlnflag)), " #", flagsn

          end if

        end do

        write(1,*) ""

        do linen = 1, size(starfile%datablocks(blockn)%datalines)

          ! only write if length is less that 1024 as these lines are empty
          if(len(trim(adjustl(starfile%datablocks(blockn)%datalines(linen)))) < 1024) then

            write(1, fmt="(A)") trim(adjustl(starfile%datablocks(blockn)%datalines(linen)))

          end if

        end do

      else

      end if

    end do

    close(1)

  end subroutine write_starfile



  subroutine mics2star(self, spproj)

    class(star_project),  intent(inout)           :: self
    type(sp_project)                              :: spproj
    type(star_file)                               :: micrographsstar
    integer                                       :: opticsblock, micsblock

    micrographsstar%filename = "micrographs.star"

    call self%insert_datablock(micrographsstar, "optics", .true., opticsblock)

    call self%insert_flag(micrographsstar, opticsblock, "rlnOpticsGroup", "ogname")
    call self%insert_flag(micrographsstar, opticsblock, "rlnOpticsGroupName", "ogid")
    call self%insert_flag(micrographsstar, opticsblock, "rlnAmplitudeContrast", "fraca")
    call self%insert_flag(micrographsstar, opticsblock, "rlnSphericalAberration", "cs")
    call self%insert_flag(micrographsstar, opticsblock, "rlnVoltage", "kv")
    call self%insert_flag(micrographsstar, opticsblock, "rlnImagePixelSize", "smpd")

    call self%insert_datablock(micrographsstar, "micrographs", .true., micsblock)

    call self%insert_flag(micrographsstar, micsblock, "rlnMicrographMovieName", "movie")
    call self%insert_flag(micrographsstar, micsblock, "rlnMicrographName", "intg")
    call self%insert_flag(micrographsstar, micsblock, "rlnMicrographMetadata", "mc_starfile")
    call self%insert_flag(micrographsstar, micsblock, "rlnDefocusU", "dfx")
    call self%insert_flag(micrographsstar, micsblock, "rlnDefocusV", "dfy")
    call self%insert_flag(micrographsstar, micsblock, "rlnDefocusAngle", "angast")
    call self%insert_flag(micrographsstar, micsblock, "rlnPhaseShift", "phshift")
    call self%insert_flag(micrographsstar, micsblock, "rlnCtfMaxResolution", "ctfres")
    call self%insert_flag(micrographsstar, micsblock, "rlnOpticsGroup", "ogid")
    call self%insert_flag(micrographsstar, micsblock, "rlnImageSizeX", "xdim")
    call self%insert_flag(micrographsstar, micsblock, "rlnImageSizeY", "ydim")
    call self%insert_flag(micrographsstar, micsblock, "rlnImageSizeZ", "nframes")
    call self%insert_flag(micrographsstar, micsblock, "rlnCtfImage", "ctfeps")
    call self%insert_flag(micrographsstar, micsblock, "rlnCtfPowerSpectrum", "ctfjpg")
    call self%insert_flag(micrographsstar, micsblock, "rlnMicrographCoordinates", "boxfile")

    call self%process_flags(micrographsstar, opticsblock, spproj%os_optics)
    call self%process_flags(micrographsstar, micsblock, spproj%os_mic)

    call self%process_data(micrographsstar, opticsblock, spproj%os_optics)
    call self%process_data(micrographsstar, micsblock, spproj%os_mic)

    call self%write_starfile(micrographsstar)

    call micrographsstar%kill_starfile()

  end subroutine mics2star



  subroutine ptcls2D2star(self, spproj)

    class(star_project),  intent(inout)           :: self
    type(sp_project)                              :: spproj
    !type(flag), dimension(:), allocatable         :: ptclflags
    type(star_file)                               :: ptclstar
    character(len=:), allocatable                 :: stkname
    character(len=1024)                           :: stkaddr, stkmic
    integer                                       :: opticsblock, ptclblock, orisn, stkind, stkptclind, micind
    integer, dimension(:), allocatable            :: mic2stk_inds, stk2mic_inds

    call spproj%get_mic2stk_inds(mic2stk_inds, stk2mic_inds)

    ! Assign stknames and ptcl indices

    !$omp parallel do default(shared) private(orisn,stkname,stkind,stkaddr,stkptclind, micind) schedule(static) proc_bind(close)
    do orisn = 1, spproj%os_ptcl2D%get_noris()

      call spproj%get_stkname_and_ind('ptcl2D', orisn, stkname, stkptclind)

      write(stkaddr, fmt="(I0,A,A)") stkptclind, "@", trim(adjustl(stkname))

      call spproj%os_ptcl2D%set(orisn, 'stkaddr', trim(adjustl(stkaddr)))

      stkind = spproj%os_ptcl2D%get(orisn, 'stkind')

      micind = stk2mic_inds(stkind)

      write(stkmic, fmt="(A)") trim(adjustl(spproj%os_mic%get_static(micind, 'intg')))

      call spproj%os_ptcl2D%set(orisn, 'intg', trim(adjustl(stkmic)))

      ! set intgs here for ptcls

    end do
    !$omp end parallel do

    ptclstar%filename = "particles.star"

    call self%insert_datablock(ptclstar, "optics", .true., opticsblock)

    call self%insert_flag(ptclstar, opticsblock, "rlnOpticsGroup", "ogname")
    call self%insert_flag(ptclstar, opticsblock, "rlnOpticsGroupName", "ogid")
    call self%insert_flag(ptclstar, opticsblock, "rlnAmplitudeContrast", "fraca")
    call self%insert_flag(ptclstar, opticsblock, "rlnSphericalAberration", "cs")
    call self%insert_flag(ptclstar, opticsblock, "rlnVoltage", "kv")
    call self%insert_flag(ptclstar, opticsblock, "rlnImagePixelSize", "smpd")
    call self%insert_flag(ptclstar, opticsblock, "rlnImageDimensionality", "ogdim")
    call self%insert_flag(ptclstar, opticsblock, "rlnImageSize", "ogsize")

    call self%insert_datablock(ptclstar, "particles", .true., ptclblock)

    call self%insert_flag(ptclstar, ptclblock, "rlnImageName", "stkaddr")
    call self%insert_flag(ptclstar, ptclblock, "rlnMicrographName", "intg")
    call self%insert_flag(ptclstar, ptclblock, "rlnCoordinateX", "xpos")
    call self%insert_flag(ptclstar, ptclblock, "rlnCoordinateY", "ypos")
    call self%insert_flag(ptclstar, ptclblock, "rlnDefocusU", "dfx")
    call self%insert_flag(ptclstar, ptclblock, "rlnDefocusV", "dfy")
    call self%insert_flag(ptclstar, ptclblock, "rlnDefocusAngle", "angast")
    call self%insert_flag(ptclstar, ptclblock, "rlnOpticsGroup", "ogid")
    call self%insert_flag(ptclstar, ptclblock, "rlnAngleRot", "e1")
    call self%insert_flag(ptclstar, ptclblock, "rlnAngleTilt", "e2")
    call self%insert_flag(ptclstar, ptclblock, "rlnAnglePsi", "e3")

    call self%process_flags(ptclstar, opticsblock, spproj%os_optics)
    call self%process_flags(ptclstar, ptclblock, spproj%os_ptcl2D)

    call self%process_data(ptclstar, opticsblock, spproj%os_optics)
    call self%process_data(ptclstar, ptclblock, spproj%os_ptcl2D)

    call self%write_starfile(ptclstar)

    call ptclstar%kill_starfile()

  end subroutine ptcls2D2star



  subroutine clusters2D2star(self, spproj)

    class(star_project),  intent(inout)           :: self
    type(sp_project)                              :: spproj
    type(star_file)                               :: clustersstar
    integer                                       :: generalblock, clustersblock, orisn, stkind
    character(len=:), allocatable                 :: stkname
    character(len=1024)                           :: stkaddr

    stkname = "cluster2D.mrc"

    do orisn = 1, spproj%os_cls2D%get_noris()

      stkind = spproj%os_cls2D%get(orisn, 'class')

      write(stkaddr, fmt="(I0,A,A)") stkind, "@", trim(adjustl(stkname))

      call spproj%os_cls2D%set(orisn, 'stkaddr', trim(adjustl(stkaddr)))

    end do

    clustersstar%filename = "clusters2D.star"

    call self%insert_datablock(clustersstar, "general", .false., generalblock)

    call self%insert_flag(clustersstar, generalblock, "rlnReferenceDimensionality", "")
    call self%insert_flag(clustersstar, generalblock, "rlnDataDimensionality", "")
    call self%insert_flag(clustersstar, generalblock, "rlnOriginalImageSize", "")
    call self%insert_flag(clustersstar, generalblock, "rlnCurrentResolution", "")
    call self%insert_flag(clustersstar, generalblock, "rlnCurrentImageSize", "")
    call self%insert_flag(clustersstar, generalblock, "rlnNrClasses", "")

    call self%insert_datablock(clustersstar, "clusters", .true., clustersblock)

    call self%insert_flag(clustersstar, clustersblock, "rlnReferenceImage", "stkaddr")
    call self%insert_flag(clustersstar, clustersblock, "rlnClassDistribution", "pop")
    call self%insert_flag(clustersstar, clustersblock, "rlnEstimatedResolution", "res")

    call self%process_flags(clustersstar, generalblock, spproj%os_cls2D)
    call self%process_flags(clustersstar, clustersblock, spproj%os_cls2D)

    call self%process_data(clustersstar, generalblock, spproj%os_cls2D)
    call self%process_data(clustersstar, clustersblock, spproj%os_cls2D)

    call self%write_starfile(clustersstar)

    call clustersstar%kill_starfile()

  end subroutine clusters2D2star


  subroutine create(self, cline, spproj)

    class(star_project),  intent(inout)   :: self
    class(cmdline), intent(inout)         :: cline
    type(sp_project)                      :: spproj

    ! Assign single optic group if unset
    if(spproj%os_optics%get_noris() == 0) then
      call self%assign_optics_single(spproj)
    endif


    ! Process micrographs
    if(spproj%os_mic%get_noris() > 0) then
      !Need to assign optics groups else assume 1
      call self%mics2star(spproj)

    endif

    if(spproj%os_ptcl2D%get_noris() > 0) then
      !Need to assign optics groups else assume 1
      call self%ptcls2D2star(spproj)

    endif

    if(spproj%os_cls2D%get_noris() > 0) then
      !Need to assign optics groups else assume 1
      call self%clusters2D2star(spproj)

    endif

   ! call spproj%os_mic%pparms2str(1)
   ! call spproj%os_mic%ori2str(1)

  end subroutine create



end module simple_starproject
