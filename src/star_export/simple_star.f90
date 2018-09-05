!! Importing and exporting Relion Star-formatted files to/from SIMPLE
module simple_star
include 'simple_lib.f08'
use simple_stardoc
use simple_sp_project, only: sp_project
use simple_cmdline,    only: cmdline
use simple_parameters
use simple_binoris_io
use simple_oris, only: oris
implicit none
private
public :: star_project
#include "simple_local_flags.inc"


type star_project
    type(stardoc) :: doc
contains
    procedure :: new
    procedure :: prepareimport
    procedure :: readfile
    procedure :: get_ndatalines
    procedure :: get_nrecs_per_line
    procedure :: check_temp_files
    procedure :: read
    procedure :: print_info
    procedure :: print_valid_import_startypes
    !! Explicit import/export methods
    procedure :: export_micrographs
    procedure :: import_micrographs
    procedure :: export_motion_corrected_micrographs
    procedure :: import_motion_corrected_micrographs
    procedure :: export_ctf_estimation
    procedure :: import_ctf_estimation
    procedure :: export_autopick
    procedure :: import_autopick
    procedure :: export_extract_doseweightedptcls
    procedure :: import_extract_doseweightedptcls
    procedure :: export_class2D
    procedure :: import_class2D
    procedure :: export_class2D_select
    procedure :: import_class2D_select
    procedure :: export_init3Dmodel
    procedure :: import_init3Dmodel
    procedure :: export_class3D
    procedure :: import_class3D
    procedure :: export_refine3D
    procedure :: import_refine3D
    procedure :: export_shiny3D
    procedure :: import_shiny3D
    procedure :: export_all
    procedure :: import_all
    procedure :: import_particles
    procedure :: import_cavgs
    procedure :: kill
end type star_project

interface star_project
    module procedure constructor
end interface star_project

contains

    !>  \brief  is an abstract constructor
    function constructor( filename ) result( self )
        character(len=*),intent(inout),optional :: filename
        type(star_project) :: self
        call self%new(filename)
    end function constructor

    subroutine new( self, filename )
        class(star_project), intent(inout) :: self
        character(len=*),intent(inout),optional :: filename
        call self%kill()

        if(present(filename)) then
            if( .not. file_exists(trim(filename)) )then
                !! import mode
                call self%doc%open4import(filename)
            else
                !! export mode
            endif
        endif
    end subroutine new

    subroutine prepareimport(self, sp, p, filename)
        class(star_project), intent(inout) :: self
        class(sp_project), intent(inout)   :: sp
        class(parameters), intent(inout) :: p
        character(len=*),intent(inout) :: filename
        if( .not. file_exists(trim(filename)) )then
            THROW_HARD('file: '//trim(filename)//' not in cwd; prepare')
        endif
        !! import mode
        call self%doc%open4import(filename)
        call self%doc%close()
    end subroutine prepareimport

    subroutine readfile(self, sp, filename)
        class(star_project), intent(inout) :: self
        class(sp_project), intent(inout)   :: sp
        character(len=*), intent(inout)       :: filename
    end subroutine readfile

    function get_ndatalines(self) result(n)
        class(star_project), intent(inout) :: self
        integer :: n
        n=0
        if(.not. self%doc%existence) THROW_HARD('get_ndatalines doc unopened')
        if( self%doc%num_data_lines == 0) then
            print *," simple_star :: get_ndatalines no data entries"
        else
            n= self%doc%num_data_lines
        endif
    end function get_ndatalines

    function get_nrecs_per_line(self)result(nrecs)
        class(star_project), intent(inout) :: self
        integer :: nrecs
        nrecs=0
        if(.not. self%doc%existence) THROW_HARD('get_ndatalines doc unopened')
        if( self%doc%num_data_elements == 0) then
            print *," simple_star :: get_ndatalines no data entries"
        else
            nrecs = self%doc%num_data_elements
        endif
    end function get_nrecs_per_line

    subroutine read( self, fname )
        class(star_project), intent(inout)     :: self
        character(len=*), optional, intent(inout) :: fname
        call self%doc%setdoprint()
        call self%doc%open4import(fname)

    end subroutine read


    subroutine print_info( self )
        class(star_project), intent(inout)     :: self
        call self%doc%print()
    end subroutine print_info

    subroutine check_temp_files(self,msg)
        class(star_project), intent(inout) :: self
        character(len=*), intent(in)       :: msg
        integer :: nl_filetab, nl_deftab

        if (.not. file_exists('oritab-stardoc.txt')) then
            print *,"star_project; check_temp_files called from "//trim(msg)
            THROW_HARD( 'star_project;  oritab-stardoc.txt not found' )
        endif
        nl_deftab = nlines('oritab-stardoc.txt')
        if(nl_deftab /= self%doc%num_data_lines)then
            print *,"star_project; check_temp_files called from "//trim(msg)
            THROW_HARD('star_project;  oritab-stardoc.txt lines not same as num_data_lines')
        endif
        if (.not. file_exists('filetab-stardoc.txt')) then
            print *,"star_project; check_temp_files called from "//trim(msg)
            THROW_HARD('star_project;  filetab-stardoc.txt not found')
        endif
        nl_filetab= nlines('filetab-stardoc.txt')
        if(nl_deftab /= self%doc%num_data_lines)then
            print *,"star_project; check_temp_files called from "//trim(msg)
            THROW_HARD('star_project;  oritab-stardoc.txt lines not same as num_data_lines')
        endif
        if(nl_filetab /= nl_deftab)then
            print *,"star_project; check_temp_files called from "//trim(msg)
            THROW_HARD('star_project;  oritab-stardoc.txt lines not same as num_data_lines')
        endif

    end subroutine check_temp_files

    subroutine export_micrographs (self, sp, filename)
        class(star_project), intent(inout) :: self
        class(sp_project), intent(inout)   :: sp
        character(len=*), intent(inout)       :: filename
        character(len=KEYLEN),allocatable      :: labels(:)
        labels=(/ 'MicrographNameNoDW' /)
        call self%doc%write(filename, sp, labels)
    end subroutine export_micrographs
    subroutine import_micrographs (self, sp, params, cline, filename)
        class(star_project), intent(inout) :: self
        class(sp_project),   intent(inout) :: sp
        class(parameters),   intent(inout) :: params
        class(cmdline),      intent(inout) :: cline
        character(len=*),    intent(inout) :: filename

    end subroutine import_micrographs

    !! Motion Correct: Simple/Unblur/MotionCorr2
    !! Format of Motion Correct Star_Project File:
    !! data_
    !! loop_
    !! _rlnMicrographNameNoDW #1
    !! _rlnMicrographName #2
    !! [#1 MotionCorr/job026/Micrographs3/FoilHole_24003709_Data_23978423_23978424_20180225_0629-1729_noDW.mrc] [#2 MotionCorr/job026/ Micrographs3/FoilHole_24003709_Data_23978423_23978424_20180225_0629-1729.mrc]
    !! ...
    subroutine export_motion_corrected_micrographs (self, sp, filename)
        class(star_project), intent(inout) :: self
        class(sp_project), intent(inout)   :: sp
        character(len=*), intent(inout) :: filename
        character(len=KEYLEN),allocatable:: labels(:)
        labels=(/  &
            'MicrographNameNoDW',&
            'MicrographName    '/)
        call self%doc%write(filename, sp, labels)
    end subroutine export_motion_corrected_micrographs
    subroutine import_motion_corrected_micrographs (self, sp, params, cline, filename)
        class(star_project), intent(inout) :: self
        class(sp_project),   intent(inout) :: sp
        class(parameters),   intent(inout) :: params
        class(cmdline),      intent(inout) :: cline
        character(len=*),    intent(inout) :: filename
    end subroutine import_motion_corrected_micrographs

    !! CTF Estimation: Simple/GCTF/CTFFIND4
    !! Format of CTF estimation Star file:
    !! data_
    !! loop_
    !! _rlnMicrographNameNoDW #1            <-- unneeded
    !! _rlnMicrographName #2
    !! _rlnCtfImage #3
    !! _rlnDefocusU #4
    !! _rlnDefocusV #5
    !! _rlnDefocusAngle #6
    !! _rlnVoltage #7
    !! _rlnSphericalAberration #8
    !! _rlnAmplitudeContrast #9
    !! _rlnMagnification #10
    !! _rlnDetectorPixelSize #11
    !! _rlnCtfFigureOfMerit #12
    !! [noDW filename #1] [ctf filename #2] ... [#12]
    !! others : CtfMaxResolution
    subroutine export_ctf_estimation (self, sp, filename)
        class(star_project), intent(inout) :: self
        class(sp_project), intent(inout)   :: sp
        character(len=*), intent(inout) :: filename
        character(len=KEYLEN),allocatable:: labels(:)
        labels=(/  'MicrographName     ',&
            'CtfImage           ',&
            'DefocusU           ',&
            'DefocusV           ',&
            'DefocusAngle       ',&
            'Voltage            ',&
            'SphericalAberration',&
            'AmplitudeContrast  ',&
            'Magnification      ',&
            'DetectorPixelSize  ',&
            'CtfFigureOfMerit   ',&
            'CtfMaxResolution   ' /)
        call self%doc%write(filename, sp, labels)
    end subroutine export_ctf_estimation
    subroutine import_ctf_estimation (self, spproj, params, cline, filename)
        use simple_sp_project, only: oritype2segment
        class(star_project), intent(inout) :: self
        class(sp_project),   intent(inout) :: spproj
        class(parameters),   intent(inout) :: params
        class(cmdline),      intent(inout) :: cline
        character(len=*),    intent(inout) :: filename
        integer         :: i, nheaders, nlines
        type(oris)      :: os
        integer         :: ndatlines,noris, n_ori_inputs
        integer(kind(ENUM_ORISEG)) :: iseg
        type(ctfparams) :: ctfvars
        logical         :: inputted_oritab, inputted_plaintexttab, inputted_deftab
        character(len=:), allocatable ::  oritype
        class(oris),      pointer     :: os_ptr
        call cline%set('oritype','mic')
        !! make sure starfile has been parsed and temporary files are in current dir
        call self%check_temp_files('import_ctf_estimation')
        !! set deftab
        call cline%set('deftab', 'oritab-stardoc.txt')
        if(.not. (cline%defined('filetab') .or. cline%defined('stktab')) ) then
            call cline%set('filetab', 'filetab-stardoc.txt')
        endif

        call params%new(cline)

        ndatlines = binread_nlines(params%oritab)
        if(ndatlines /= self%doc%num_data_lines)then
            THROW_HARD('star_project; import_ctf_estimation binread_nlines does not match num_data_lines from starfile')
        endif
        call os%new(ndatlines)
        call binread_oritab(params%oritab, spproj, os, [1,ndatlines])
        !!       call spproj%kill ! for safety

         !! load ctf information into sp%os_stk
         noris = self%get_ndatalines()
         iseg  = oritype2segment(params%oritype)

!         call os%new(noris)
         !! Convert data
         do i = 1, noris

        !     ! call self%o(i)%read(fnr)

        !     ! call sauron_line_parser( line, self%htab, self%chtab )
        !     ! isthere(1) = self%htab%isthere('e1')
        !     ! isthere(2) = self%htab%isthere('e2')
        !     ! isthere(3) = self%htab%isthere('e3')
        !     ! if( any(isthere) )&
        !     !     &call self%set_euler([self%htab%get('e1'),self%htab%get('e2'),self%htab%get('e3')])
        !     ! if( present(nst) )then
        !     !     state = self%o(i)%get_state()
        !     !     nst   = max(1,max(state,nst))
        !     ! endif
         end do

         !!Insert oris into sp project
         call spproj%set_sp_oris(params%oritype, os)

        ! call params%new(cline)
        ! ! parameter input management
        ! inputted_boxtab = cline%defined('boxtab')
        ! ! project file management
        ! if( .not. file_exists(trim(params%projfile)) )then
        !      ! raise exception
        ! endif
        ! call spproj%read(params%projfile)
        ! ! CTF
        ! if( cline%defined('phaseplate') )then
        !     phaseplate = cline%get_carg('phaseplate')
        ! else
        !     allocate(phaseplate, source='no')
        ! endif
        ! ctfvars%smpd  = params%smpd
        ! ctfvars%kv    = params%kv
        ! ctfvars%cs    = params%cs
        ! ctfvars%fraca = params%fraca
        ! select case(params%ctf)
        !     case('yes')
        !         ctfvars%ctfflag = 1
        !     case('no')
        !         ctfvars%ctfflag = 0
        !     case('flip')
        !         ctfvars%ctfflag = 2
        !     case DEFAULT
        !         ! raise exception
        ! end select
        ! ctfvars%l_phaseplate = .false.
        ! if( trim(params%phaseplate) .eq. 'yes' ) ctfvars%l_phaseplate = .true.
        ! ! update project info
        ! call spproj%update_projinfo( cline )

        ! ! updates segment
        ! call spproj%add_movies(params%filetab, ctfvars)
        ! do i = 1, self%doc%num_data_lines
        !     call simple_abspath(, fname, 'simple_sp_project::add_single_movie')

        !     call sp%add_single_movie
        ! end do

        ! ! add boxtab
        ! ! if( inputted_boxtab )then
        ! !     call read_filetable(params%boxtab, boxfnames)
        ! !     nboxf = size(boxfnames)
        ! !     nmovf = nlines(params%filetab)
        ! !     if( nboxf /= nmovf )then
        ! !         ! raise exception
        ! !     endif
        ! !     do i=1,nmovf
        ! !         call simple_abspath(trim(boxfnames(i)), boxf_abspath, 'commander_project :: exec_import_movies')
        ! !         call spproj%os_mic%set(i, 'boxfile', boxf_abspath)
        ! !     end do
        ! ! endif

    end subroutine import_ctf_estimation


    !! Autopick  :  Simple/Relion/Gautomatch/EMAN
    !! Star format:
    !! data_
    !! loop_
    !! _rlnCoordinateX #1
    !! _rlnCoordinateY #2
    !! _rlnClassNumber #3
    !! _rlnAutopickFigureOfMerit #4
    !! _rlnAnglePsi #5
    !!  1219.364633   233.878135            1     1.523455     0.000000
    subroutine export_autopick (self, sp, filename)
        class(star_project), intent(inout) :: self
        class(sp_project), intent(inout)   :: sp
        character(len=*), intent(inout) :: filename
        character(len=KEYLEN), allocatable:: labels(:)
        labels=(/ &
            'CoordinateX          ',&
            'CoordinateY          ',&
            'ClassNumber          ',&
            'AutopickFigureOfMerit',&
            'AnglePsi             ' /)
        call self%doc%write(filename, sp, labels)

    end subroutine export_autopick
    subroutine import_autopick (self, spproj, params, cline, filename)
        class(star_project), intent(inout) :: self
        class(sp_project),   intent(inout) :: spproj
        class(parameters),   intent(inout) :: params
        class(cmdline),      intent(inout) :: cline
        character(len=*),    intent(inout) :: filename
    end subroutine import_autopick


    !! Extract_Doseweighted_Ptcls  :  Simple/Relion/Gautomatch/EMAN
    !! Star format:
    !! data_
    !! loop_
    !! _rlnCoordinateX #1
    !! _rlnCoordinateY #2
    !! _rlnClassNumber #3
    !! _rlnAutopickFigureOfMerit #4
    !! _rlnAnglePsi #5
    !! _rlnImageName #6
    !! _rlnMicrographName #7
    !! _rlnVoltage #8
    !! _rlnDefocusU #9
    !! _rlnDefocusV #10
    !! _rlnDefocusAngle #11
    !! _rlnSphericalAberration #12
    !! _rlnCtfBfactor #13
    !! _rlnCtfScalefactor #14
    !! _rlnPhaseShift #15
    !! _rlnAmplitudeContrast #16
    !! _rlnMagnification #17
    !! _rlnDetectorPixelSize #18
    !! _rlnCtfFigureOfMerit #19
    !! _rlnGroupNumber #20
    !! _rlnAngleRot #21
    !! _rlnAngleTilt #22
    !! _rlnOriginX #23
    !! _rlnOriginY #24
    !! _rlnNormCorrection #25
    !! _rlnLogLikeliContribution #26
    !! _rlnMaxValueProbDistribution #27
    !! _rlnNrOfSignificantSamples #28
    !! _rlnGroupName #29
    !!  [#1] ... [#6 filename] [#7 filename] ... [#29]
    subroutine export_extract_doseweightedptcls (self, sp, filename)
        class(star_project)   , intent(inout) :: self
        class(sp_project)     , intent(inout) :: sp
        character(len=*)      , intent(inout) :: filename
        character(len=KEYLEN) , allocatable   :: labels(:)
          labels=(/ &
            'CoordinateX              ',&
            'CoordinateY              ',&
            'ClassNumber              ',&
            'AutopickFigureOfMerit    ',&
            'AnglePsi                 ',&
            'ImageName                ',&
            'MicrographName           ',&
            'Voltage                  ',&
            'DefocusU                 ',&
            'DefocusV                 ',&
            'DefocusAngle             ',&
            'SphericalAberration      ',&
            'CtfBfactor               ',&
            'CtfScalefactor           ',&
            'PhaseShift               ',&
            'AmplitudeContrast        ',&
            'Magnification            ',&
            'DetectorPixelSize        ',&
            'CtfFigureOfMerit         ',&
            'GroupNumber              ',&
            'AngleRot                 ',&
            'AngleTilt                ',&
            'OriginX                  ',&
            'OriginY                  ',&
            'NormCorrection           ',&
            'LogLikeliContribution    ',&
            'MaxValueProbDistribution ',&
            'NrOfSignificantSamples   ',&
            'GroupName                ' /)
        call self%doc%write(filename, sp, labels)

    end subroutine export_extract_doseweightedptcls
    subroutine import_extract_doseweightedptcls (self, spproj, params, cline, filename)
        class(star_project), intent(inout) :: self
        class(sp_project),   intent(inout) :: spproj
        class(parameters),   intent(inout) :: params
        class(cmdline),      intent(inout) :: cline
        character(len=*),    intent(inout) :: filename
    end subroutine import_extract_doseweightedptcls


    !! Class2D classification : SIMPLE/Relion/CisTEM/CryoSparc
    !! data_
    !!
    !! loop_
    !! _rlnMicrographName #1
    !! _rlnCoordinateX #2
    !! _rlnCoordinateY #3
    !! _rlnVoltage #4
    !! _rlnDefocusU #5
    !! _rlnDefocusV #6
    !! _rlnDefocusAngle #7
    !! _rlnSphericalAberration #8
    !! _rlnDetectorPixelSize #9
    !! _rlnCtfFigureOfMerit #10
    !! _rlnMagnification #11
    !! _rlnAmplitudeContrast #12
    !! _rlnAutopickFigureOfMerit #13
    !! _rlnPhaseShift #14
    !! _rlnImageName #15
    !! _rlnGroupNumber #16
    !! _rlnAngleRot #17
    !! _rlnAngleTilt #18
    !! _rlnAnglePsi #19
    !! _rlnOriginX #20
    !! _rlnOriginY #21
    !! _rlnClassNumber #22
    !! _rlnNormCorrection #23
    !! _rlnLogLikeliContribution #24
    !! _rlnMaxValueProbDistribution #25
    !! _rlnNrOfSignificantSamples #26
    !! _rlnGroupName #27
    !! _rlnRandomSubset #28
    !! _rlnOriginalParticleName #29
    !! _rlnNrOfFrames #30
    !! _rlnAverageNrOfFrames #31
    !! _rlnMovieFramesRunningAverage #32
    !! [#1 MC filename] ... [#15 Polish filename] .. [#29 Extract filename] ..
    subroutine export_class2D (self, sp, filename)
        class(star_project), intent(inout) :: self
        class(sp_project), intent(inout)   :: sp
        character(len=*), intent(inout) :: filename
        character(len=KEYLEN), allocatable:: labels(:)
        labels=(/ &
            'MicrographName           ',&
            'CoordinateX              ',&
            'CoordinateY              ',&
            'Voltage                  ',&
            'DefocusU                 ',&
            'DefocusV                 ',&
            'DefocusAngle             ',&
            'SphericalAberration      ',&
            'DetectorPixelSize        ',&
            'CtfFigureOfMerit         ',&
            'Magnification            ',&
            'AmplitudeContrast        ',&
            'AutopickFigureOfMerit    ',&
            'PhaseShift               ',&
            'ImageName                ',&
            'GroupNumber              ',&
            'AngleRot                 ',&
            'AngleTilt                ',&
            'AnglePsi                 ',&
            'OriginX                  ',&
            'OriginY                  ',&
            'ClassNumber              ',&
            'NormCorrection           ',&
            'LogLikeliContribution    ',&
            'MaxValueProbDistribution ',&
            'NrOfSignificantSamples   ',&
            'GroupName                ',&
            'RandomSubset             ',&
            'OriginalParticleName     ',&
            'NrOfFrames               ',&
            'AverageNrOfFrames        ',&
            'MovieFramesRunningAverage' /)
        call self%doc%write(filename, sp, labels)

    end subroutine export_class2D
    subroutine import_class2D (self, spproj, params, cline, filename)
        class(star_project), intent(inout) :: self
        class(sp_project),   intent(inout) :: spproj
        class(parameters),   intent(inout) :: params
        class(cmdline),      intent(inout) :: cline
        character(len=*),    intent(inout) :: filename
    end subroutine import_class2D

    !! Class 2D select
    !! data_
    !!
    !! loop_
    !! _rlnMicrographName #1
    !! _rlnCoordinateX #2
    !! _rlnCoordinateY #3
    !! _rlnVoltage #4
    !! _rlnDefocusU #5
    !! _rlnDefocusV #6
    !! _rlnDefocusAngle #7
    !! _rlnSphericalAberration #8
    !! _rlnDetectorPixelSize #9
    !! _rlnCtfFigureOfMerit #10
    !! _rlnMagnification #11
    !! _rlnAmplitudeContrast #12
    !! _rlnAutopickFigureOfMerit #13
    !! _rlnPhaseShift #14
    !! _rlnImageName #15
    !! _rlnGroupNumber #16
    !! _rlnAngleRot #17
    !! _rlnAngleTilt #18
    !! _rlnAnglePsi #19
    !! _rlnOriginX #20
    !! _rlnOriginY #21
    !! _rlnClassNumber #22
    !! _rlnNormCorrection #23
    !! _rlnLogLikeliContribution #24
    !! _rlnMaxValueProbDistribution #25
    !! _rlnNrOfSignificantSamples #26
    !! _rlnGroupName #28                                    <-- Possible index error
    !! [#1 MC filename] ... [#15 Polish filename] .. [#27]
    subroutine export_class2D_select (self, sp, filename)
        class(star_project), intent(inout) :: self
        class(sp_project), intent(inout)   :: sp
        character(len=*), intent(inout) :: filename
        character(len=KEYLEN), allocatable:: labels(:)
        labels=(/ &
            'MicrographName           ',&
            'CoordinateX              ',&
            'CoordinateY              ',&
            'Voltage                  ',&
            'DefocusU                 ',&
            'DefocusV                 ',&
            'DefocusAngle             ',&
            'SphericalAberration      ',&
            'DetectorPixelSize        ',&
            'CtfFigureOfMerit         ',&
            'Magnification            ',&
            'AmplitudeContrast        ',&
            'AutopickFigureOfMerit    ',&
            'PhaseShift               ',&
            'ImageName                ',&
            'GroupNumber              ',&
            'AngleRot                 ',&
            'AngleTilt                ',&
            'AnglePsi                 ',&
            'OriginX                  ',&
            'OriginY                  ',&
            'ClassNumber              ',&
            'NormCorrection           ',&
            'LogLikeliContribution    ',&
            'MaxValueProbDistribution ',&
            'NrOfSignificantSamples   ',&
            'GroupName                ' /)
        call self%doc%write(filename, sp, labels)

    end subroutine export_class2D_select
    subroutine import_class2D_select (self, spproj, params, cline, filename)
        class(star_project), intent(inout) :: self
        class(sp_project),   intent(inout) :: spproj
        class(parameters),   intent(inout) :: params
        class(cmdline),      intent(inout) :: cline
        character(len=*),    intent(inout) :: filename
    end subroutine import_class2D_select

    !! Class3D
    !! Possible apps: SIMPLE/EMAN/Relion/CryoSparc or PreviousData or GCTF for local ctf estimate
    !! Star format for Refine3D:
    !! data_model_general
    !!
    !! _rlnReferenceDimensionality                        3
    !! _rlnDataDimensionality                             2
    !! _rlnOriginalImageSize                            128
    !! _rlnCurrentResolution                       9.307315
    !! _rlnCurrentImageSize                             102
    !! _rlnPaddingFactor                           2.000000
    !! _rlnIsHelix                                        0
    !! _rlnFourierSpaceInterpolator                       1
    !! _rlnMinRadiusNnInterpolation                      10
    !! _rlnPixelSize                               2.981249
    !! _rlnNrClasses                                      4
    !! _rlnNrBodies                                       1
    !! _rlnNrGroups                                      39
    !! _rlnTau2FudgeFactor                         4.000000
    !! _rlnNormCorrectionAverage                   0.738236
    !! _rlnSigmaOffsets                            9.181616
    !! _rlnOrientationalPriorMode                         0
    !! _rlnSigmaPriorRotAngle                      0.000000
    !! _rlnSigmaPriorTiltAngle                     0.000000
    !! _rlnSigmaPriorPsiAngle                      0.000000
    !! _rlnLogLikelihood                       2.020615e+09
    !! _rlnAveragePmax                             0.975562
    !!
    !! data_model_classes
    !!
    !! loop_
    !! _rlnReferenceImage #1
    !! _rlnClassDistribution #2
    !! _rlnAccuracyRotations #3
    !! _rlnAccuracyTranslations #4
    !! _rlnEstimatedResolution #5
    !! _rlnOverallFourierCompleteness #6
    !! Class3D/job067/run_it025_class001.mrc     0.280179     1.030000     0.301000     9.307315     0.983456
    !! "                                   "          " "          " "          " "          " "          " " x4
    !!
    !! data_model_class_1                                                                                         \
    !!                                                                                                            |
    !! loop_                                                                                                      |
    !! _rlnSpectralIndex #1                                                                                       |
    !! _rlnResolution #2                                                                                          |
    !! _rlnAngstromResolution #3                                                                                  |
    !! _rlnSsnrMap #4                                                                                             |
    !! _rlnGoldStandardFsc #5                                                                                     |
    !! _rlnFourierCompleteness #6                                                                                 |
    !! _rlnReferenceSigma2 #7                                                                                     |
    !! _rlnReferenceTau2 #8                                                                                       |
    !! _rlnSpectralOrientabilityContribution #9                                                                   |
    !! 0     0 .000000   999.000000 5.066702e+06     0.000000     1.000000 1.888607e-07     0.239225     0.000000 |
    !! "                                   "          " "          " "          " "          " "          " " x64 |
    !!                                                                                                            /x4
    !! data_model_groups
    !!
    !! loop_
    !! _rlnGroupNumber #1
    !! _rlnGroupName #2
    !! _rlnGroupNrParticles #3
    !! _rlnGroupScaleCorrection #4
    !!       1 group_24         1336     1.023934
    !! "                                   "          " "          " "          " "          " "          " " x39
    !!
    !! data_model_group_1                                                                                         \
    !!                                                                                                            |
    !! loop_                                                                                                      |
    !! _rlnSpectralIndex #1                                                                                       |
    !! _rlnResolution #2                                                                                          |
    !! _rlnSigma2Noise #3                                                                                         |
    !! 0     0 .000000 3.814716e-04                                                                               |
    !! " "   "       " "          "x64                                                                            |
    !!                                                                                                            /x39
    !!
    !! data_model_pdf_orient_class_1                                                                              \
    !!                                                                                                            |
    !! loop_                                                                                                      |
    !! _rlnOrientationDistribution #1                                                                             |
    !! 1 .326604e-04                                                                                              |             " "  x~1000                                                                                             |
    !!                                                                                                            /x4
    !!
    !!
    subroutine export_class3D (self, sp, filename)
        class(star_project), intent(inout) :: self
        class(sp_project),   intent(inout) :: sp
        character(len=*),    intent(inout) :: filename
    end subroutine export_class3D
    subroutine import_class3D (self, spproj, params, cline, filename)
        class(star_project), intent(inout) :: self
        class(sp_project),   intent(inout) :: spproj
        class(parameters),   intent(inout) :: params
        class(cmdline),      intent(inout) :: cline
        character(len=*),    intent(inout) :: filename
    end subroutine import_class3D

    subroutine export_init3Dmodel (self, sp, filename)
        class(star_project), intent(inout) :: self
        class(sp_project),   intent(inout) :: sp
        character(len=*),    intent(inout) :: filename
        character(len=KEYLEN), allocatable :: labels(:)
        labels=(/ &
            'MicrographName           ',&
            'CoordinateX              ',&
            'CoordinateY              ',&
            'Voltage                  ',&
            'DefocusU                 ',&
            'DefocusV                 ',&
            'DefocusAngle             ',&
            'SphericalAberration      ',&
            'DetectorPixelSize        ',&
            'CtfFigureOfMerit         ',&
            'Magnification            ',&
            'AmplitudeContrast        ',&
            'AutopickFigureOfMerit    ',&
            'PhaseShift               ',&
            'ImageName                ',&
            'GroupNumber              ',&
            'AngleRot                 ',&
            'AngleTilt                ',&
            'AnglePsi                 ',&
            'OriginX                  ',&
            'OriginY                  ',&
            'ClassNumber              ',&
            'NormCorrection           ',&
            'LogLikeliContribution    ',&
            'MaxValueProbDistribution ',&
            'NrOfSignificantSamples   ',&
            'GroupName                ',&
            'RandomSubset             ',&
            'OriginalParticleName     ',&
            'NrOfFrames               ',&
            'AverageNrOfFrames        ',&
            'MovieFramesRunningAverage',&
            'RandomSubset             ' /)
        call self%doc%write(filename, sp, labels)
    end subroutine export_init3Dmodel
    subroutine import_init3Dmodel (self, spproj, params, cline, filename)
        class(star_project), intent(inout) :: self
        class(sp_project),   intent(inout) :: spproj
        class(parameters),   intent(inout) :: params
        class(cmdline),      intent(inout) :: cline
        character(len=*),    intent(inout) :: filename
    end subroutine import_init3Dmodel


    !! Refine3D SIMPLE
    !! Possible initial model apps: SIMPLE/EMAN/Relion/CryoSparc or PreviousData
    !! Refine and polish apps:  SIMPLE/CisTEM/Relion/CryoSparc
    !! Star format for Refine3D:
    !! data_
    !!
    !! loop_
    !! _rlnMicrographName #1
    !! _rlnCoordinateX #2
    !! _rlnCoordinateY #3
    !! _rlnVoltage #4
    !! _rlnDefocusU #5
    !! _rlnDefocusV #6
    !! _rlnDefocusAngle #7
    !! _rlnSphericalAberration #8
    !! _rlnDetectorPixelSize #9
    !! _rlnCtfFigureOfMerit #10
    !! _rlnMagnification #11
    !! _rlnAmplitudeContrast #12
    !! _rlnAutopickFigureOfMerit #13
    !! _rlnPhaseShift #14
    !! _rlnImageName #15
    !! _rlnGroupNumber #16
    !! _rlnAngleRot #17
    !! _rlnAngleTilt #18
    !! _rlnAnglePsi #19
    !! _rlnOriginX #20
    !! _rlnOriginY #21
    !! _rlnClassNumber #22
    !! _rlnNormCorrection #23
    !! _rlnLogLikeliContribution #24
    !! _rlnMaxValueProbDistribution #25
    !! _rlnNrOfSignificantSamples #26
    !! _rlnGroupName #27
    !! _rlnRandomSubset #28
    !! [#1 MC filename] ... [#15 Extract filename] .. [#28]
    subroutine export_refine3D (self, sp, filename)
        class(star_project), intent(inout) :: self
        class(sp_project), intent(inout)   :: sp
        character(len=*), intent(inout) :: filename
        character(len=KEYLEN), allocatable :: labels(:)
        labels=(/ &
            'MicrographName           ',&
            'CoordinateX              ',&
            'CoordinateY              ',&
            'Voltage                  ',&
            'DefocusU                 ',&
            'DefocusV                 ',&
            'DefocusAngle             ',&
            'SphericalAberration      ',&
            'DetectorPixelSize        ',&
            'CtfFigureOfMerit         ',&
            'Magnification            ',&
            'AmplitudeContrast        ',&
            'AutopickFigureOfMerit    ',&
            'PhaseShift               ',&
            'ImageName                ',&
            'GroupNumber              ',&
            'AngleRot                 ',&
            'AngleTilt                ',&
            'AnglePsi                 ',&
            'OriginX                  ',&
            'OriginY                  ',&
            'ClassNumber              ',&
            'NormCorrection           ',&
            'LogLikeliContribution    ',&
            'MaxValueProbDistribution ',&
            'NrOfSignificantSamples   ',&
            'GroupName                ',&
            'RandomSubset             ',&
            'OriginalParticleName     ',&
            'NrOfFrames               ',&
            'AverageNrOfFrames        ',&
            'MovieFramesRunningAverage',&
            'RandomSubset             ' /)
        call self%doc%write(filename, sp, labels)

    end subroutine export_refine3D
    subroutine import_refine3D (self, spproj, params, cline, filename)
        class(star_project), intent(inout) :: self
        class(sp_project),   intent(inout) :: spproj
        class(parameters),   intent(inout) :: params
        class(cmdline),      intent(inout) :: cline
        character(len=*),    intent(inout) :: filename
    end subroutine import_refine3D


    !! Shiny -- Particle polish, 3D refinement and local masking/resolution enhancements
    !! Refine and polish apps:  SIMPLE/CisTEM/Relion/CryoSparc
    !! Star format:
    !! data_
    !!
    !! loop_
    !! _rlnMicrographName #1
    !! _rlnCoordinateX #2
    !! _rlnCoordinateY #3
    !! _rlnVoltage #4
    !! _rlnDefocusU #5
    !! _rlnDefocusV #6
    !! _rlnDefocusAngle #7
    !! _rlnSphericalAberration #8
    !! _rlnDetectorPixelSize #9
    !! _rlnCtfFigureOfMerit #10
    !! _rlnMagnification #11
    !! _rlnAmplitudeContrast #12
    !! _rlnAutopickFigureOfMerit #13
    !! _rlnPhaseShift #14
    !! _rlnImageName #15
    !! _rlnGroupNumber #16
    !! _rlnAngleRot #17
    !! _rlnAngleTilt #18
    !! _rlnAnglePsi #19
    !! _rlnOriginX #20
    !! _rlnOriginY #21
    !! _rlnClassNumber #22
    !! _rlnNormCorrection #23
    !! _rlnLogLikeliContribution #24
    !! _rlnMaxValueProbDistribution #25
    !! _rlnNrOfSignificantSamples #26
    !! _rlnGroupName #27
    !! _rlnRandomSubset #28
    !! _rlnOriginalParticleName #29
    !! _rlnNrOfFrames #30
    !! _rlnAverageNrOfFrames #31
    !! _rlnMovieFramesRunningAverage #32
    !! [#1 MC filename] ... [#15 Polish filename] ..[#29 Extract filename] [#32]
    subroutine export_shiny3D (self, sp, filename)
        class(star_project), intent(inout) :: self
        class(sp_project), intent(inout)   :: sp
        character(len=*), intent(inout) :: filename
        character(len=KEYLEN), allocatable :: labels(:)
        labels=(/ &
            'MicrographName           ',&
            'CoordinateX              ',&
            'CoordinateY              ',&
            'Voltage                  ',&
            'DefocusU                 ',&
            'DefocusV                 ',&
            'DefocusAngle             ',&
            'SphericalAberration      ',&
            'DetectorPixelSize        ',&
            'CtfFigureOfMerit         ',&
            'Magnification            ',&
            'AmplitudeContrast        ',&
            'AutopickFigureOfMerit    ',&
            'PhaseShift               ',&
            'ImageName                ',&
            'GroupNumber              ',&
            'AngleRot                 ',&
            'AngleTilt                ',&
            'AnglePsi                 ',&
            'OriginX                  ',&
            'OriginY                  ',&
            'ClassNumber              ',&
            'NormCorrection           ',&
            'LogLikeliContribution    ',&
            'MaxValueProbDistribution ',&
            'NrOfSignificantSamples   ',&
            'GroupName                ',&
            'RandomSubset             ',&
            'OriginalParticleName     ',&
            'NrOfFrames               ',&
            'AverageNrOfFrames        ',&
            'MovieFramesRunningAverage' /)
        call self%doc%write(filename, sp, labels)

    end subroutine export_shiny3D
    subroutine import_shiny3D (self, spproj, params, cline, filename)
        class(star_project), intent(inout) :: self
        class(sp_project),   intent(inout) :: spproj
        class(parameters),   intent(inout) :: params
        class(cmdline),      intent(inout) :: cline
        character(len=*),    intent(inout) :: filename
    end subroutine import_shiny3D


    subroutine export_all (self, sp, filename)
        class(star_project), intent(inout) :: self
        class(sp_project),   intent(inout) :: sp
        character(len=*),    intent(inout) :: filename
    end subroutine export_all


    subroutine import_all (self, spproj, params, cline, filename)
        class(star_project), intent(inout) :: self
        class(sp_project),   intent(inout) :: spproj
        class(parameters),   intent(inout) :: params
        class(cmdline),      intent(inout) :: cline
        character(len=*),    intent(inout) :: filename

    end subroutine import_all

    !> Import class averages
    subroutine import_cavgs (self, spproj, params, cline, filename)
        class(star_project), intent(inout) :: self
        class(sp_project),   intent(inout) :: spproj
        class(parameters),   intent(inout) :: params
        class(cmdline),      intent(inout) :: cline
        character(len=*),    intent(inout) :: filename
        character(len=:), allocatable :: stk, stk_abspath
        integer :: ldim(3), nptcls, ind

        call cline%set('oritype','cavgs')
        !! make sure starfile has been parsed and temporary files are in current dir
        call self%check_temp_files('import_cavgs')
        stk='filetab-stardoc.txt'
!!        call spproj%add_cavgs2os_out(params%stk, params%smpd)
        ! full path and existence check
        stk_abspath = simple_abspath(stk,errmsg='star_project :: import_cavgs')
        ! find dimension of inputted stack
        call find_ldim_nptcls(stk_abspath, ldim, nptcls)
        ! add os_out entry
        call spproj%add_entry2os_out('cavg', ind)
        ! fill-in field
        call spproj%os_out%set(ind, 'stk',     trim(stk_abspath))
        call spproj%os_out%set(ind, 'box',     real(ldim(1)))
        call spproj%os_out%set(ind, 'nptcls',  real(nptcls))
        call spproj%os_out%set(ind, 'fromp',   1.0)
        call spproj%os_out%set(ind, 'top',     real(nptcls))
        call spproj%os_out%set(ind, 'smpd',    real(params%smpd))
        call spproj%os_out%set(ind, 'stkkind', 'single')
        call spproj%os_out%set(ind, 'imgkind', 'cavg')
        call spproj%os_out%set(ind, 'ctf',     'no')
        ! add congruent os_cls2D
        call spproj%os_cls2D%new(nptcls)
        call spproj%os_cls2D%set_all2single('state',1.)
    !---------------^^^^ add_cavgs2os_out ----------
        ! update project info
        call spproj%update_projinfo( cline )
        ! update computer environment
        call spproj%update_compenv( cline )

        !! write spproj in commander_star; import_starproject
    end subroutine import_cavgs

    subroutine import_particles (self, spproj, params, cline, filename)
        use simple_sp_project, only: oritype2segment
        use simple_oris,      only: oris
        use simple_nrtxtfile, only: nrtxtfile
        use simple_binoris_io ! use all in there
        class(star_project), intent(inout) :: self
        class(sp_project),   intent(inout) :: spproj
        class(parameters),   intent(inout) :: params
        class(cmdline),      intent(inout) :: cline
        character(len=*),    intent(inout) :: filename
        character(len=:), allocatable :: phaseplate, ctfstr
        real,             allocatable :: line(:)
        type(oris)       :: os
        type(nrtxtfile)  :: paramfile
        logical          :: inputted_oritab, inputted_plaintexttab, inputted_deftab
        integer          :: i, ndatlines, nrecs, n_ori_inputs, lfoo(3)
        type(ctfparams)  :: ctfvars

        call cline%set('oritype','ptcl2D')
        !! make sure starfile has been parsed and temporary files are in current dir
        call self%check_temp_files('import_particles')

        call params%new(cline)
        ! PARAMETER INPUT MANAGEMENT
        ! parameter input flags
        inputted_oritab       = cline%defined('oritab')
        inputted_deftab       = cline%defined('deftab')
        inputted_plaintexttab = cline%defined('plaintexttab')
        n_ori_inputs          = count([inputted_oritab,inputted_deftab,inputted_plaintexttab])

        ! exceptions
        if( n_ori_inputs > 1 )then
            THROW_HARD('multiple parameter sources inputted, please use (oritab|deftab|plaintexttab); exec_import_particles')
        endif
        if( cline%defined('stk') .and. cline%defined('stktab') )then
            THROW_HARD('stk and stktab are both defined on command line, use either or; exec_import_particles')
        endif
        if( cline%defined('stk') .or. cline%defined('stktab') )then
            if( trim(params%ctf) .ne. 'no' )then
                ! there needs to be associated parameters of some form
                if( n_ori_inputs < 1 )then
                    THROW_HARD('stk or stktab input requires associated parameter input when ctf .ne. no (oritab|deftab|plaintexttab)')
                endif
            endif
        else
            THROW_HARD('either stk or stktab needed on command line; exec_import_particles')
        endif


        ! oris input
        if( inputted_oritab )then
            ndatlines = binread_nlines(params%oritab)
            call os%new(ndatlines)
            call binread_oritab(params%oritab, spproj, os, [1,ndatlines])
            call spproj%kill ! for safety
        endif
        if( inputted_deftab )then
            ndatlines = binread_nlines(params%deftab)
            call os%new(ndatlines)
            call binread_ctfparams_state_eo(params%deftab, spproj, os, [1,ndatlines])
            call spproj%kill ! for safety
        endif
        if( cline%defined('stk') )then
            ctfvars%smpd = params%smpd
            select case(trim(params%ctf))
                case('yes')
                    ctfvars%ctfflag = CTFFLAG_YES
                case('no')
                    ctfvars%ctfflag = CTFFLAG_NO
                case('flip')
                    ctfvars%ctfflag = CTFFLAG_FLIP
                case DEFAULT
                    write(*,*)
                    THROW_HARD('unsupported ctf flag: '//trim(params%ctf)//'; exec_extract_ptcls')
            end select
            if( ctfvars%ctfflag .ne. CTFFLAG_NO )then
                ! if importing single stack of extracted particles, these are hard requirements
                if( .not. cline%defined('kv')    ) THROW_HARD('kv (acceleration voltage in kV{300}) input required when importing movies; exec_import_particles')
                if( .not. cline%defined('cs')    ) THROW_HARD('cs (spherical aberration constant in mm{2.7}) input required when importing movies; exec_import_particles')
                if( .not. cline%defined('fraca') ) THROW_HARD('fraca (fraction of amplitude contrast{0.1}) input required when importing movies; exec_import_particles')
                if( cline%defined('phaseplate') )then
                    phaseplate = cline%get_carg('phaseplate')
                else
                    allocate(phaseplate, source='no')
                endif
                ctfvars%kv           = params%kv
                ctfvars%cs           = params%cs
                ctfvars%fraca        = params%fraca
                ctfvars%l_phaseplate = phaseplate .eq. 'yes'
            endif
        else
            ! importing from stktab
            if( n_ori_inputs == 1 )then
                ! sampling distance
                call os%set_all2single('smpd', params%smpd)
                ! acceleration voltage
                if( cline%defined('kv') )then
                    call os%set_all2single('kv', params%kv)
                else
                    do i=1,ndatlines
                        if( .not. os%isthere(i, 'kv') )then
                            write(*,*) 'os entry: ', i, ' lacks acceleration volatage (kv)'
                            THROW_HARD('provide kv on command line or update input document; exec_extract_ptcls')
                        endif
                    end do
                endif
                ! spherical aberration
                if( cline%defined('cs') )then
                    call os%set_all2single('cs', params%cs)
                else
                    do i=1,ndatlines
                        if( .not. os%isthere(i, 'cs') )then
                            write(*,*) 'os entry: ', i, ' lacks spherical aberration constant (cs)'
                            THROW_HARD('provide cs on command line or update input document; exec_extract_ptcls')
                        endif
                    end do
                endif
                ! fraction of amplitude contrast
                if( cline%defined('fraca') )then
                    call os%set_all2single('fraca', params%fraca)
                else
                    do i=1,ndatlines
                        if( .not. os%isthere(i, 'fraca') )then
                            write(*,*) 'os entry: ', i, ' lacks fraction of amplitude contrast (fraca)'
                            THROW_HARD('provide fraca on command line or update input document; exec_extract_ptcls')
                        endif
                    end do
                endif
                ! phase-plate
                if( cline%defined('phaseplate') )then
                    call os%set_all2single('phaseplate', trim(params%phaseplate))
                else
                    do i=1,ndatlines
                        if( .not. os%isthere(i, 'phaseplate') )then
                            call os%set(i, 'phaseplate', 'no')
                        endif
                    end do
                endif
                call os%getter(1, 'phaseplate', phaseplate)
                if( trim(phaseplate) .eq. 'yes' )then
                    if( .not. os%isthere(1,'phshift') )then
                        THROW_HARD('phaseplate .eq. yes requires phshift input, currently lacking; exec_import_particles')
                    endif
                endif
                ! ctf flag
                if( cline%defined('ctf') )then
                    call os%set_all2single('ctf', trim(params%ctf))
                else
                    do i=1,ndatlines
                        if( .not. os%isthere(i, 'ctf') )then
                            call os%set(i, 'ctf', 'yes')
                        endif
                    end do
                endif
                call os%getter(1, 'ctf', ctfstr)
                if( trim(ctfstr) .ne. 'no' )then
                    if( .not. os%isthere(1,'dfx') )then
                        THROW_HARD('ctf .ne. no requires dfx input, currently lacking; exec_import_particles')
                    endif
                endif
            endif
        endif

        ! PROJECT FILE MANAGEMENT
! already read in         call spproj%read(params%projfile)

        ! UPDATE FIELDS
        ! add stack if present
        ! if( cline%defined('stk') )then
        !     if( n_ori_inputs == 0 .and. trim(params%ctf) .eq. 'no' )then
        !         ! get number of particles from stack
        !         call find_ldim_nptcls(params%stk, lfoo, params%nptcls)
        !         call os%new(params%nptcls)
        !         call os%set_all2single('state', 1.0)
        !     endif
            call spproj%add_single_stk(params%stk, ctfvars, os)
!        endif
        ! add list of stacks (stktab) if present
!        if( cline%defined('stktab') ) call spproj%add_stktab(params%stktab, os)



    end subroutine import_particles



    function exporttype2star( exporttype ) result(stype)
        character(len=*),  intent(in) :: exporttype
        integer(kind(ENUM_STARTYPE))  :: stype
        !select case(trim(exporttype))
        ! match oritype equivalent for star export type
        !case('movies') ! movies
        if(  index(trim(exporttype), 'movie')/=0 )then
            stype = MOV_STAR
            !case('mic':'micrographs': 'mcmicrographs':'ctf_estimation') ! movies
        else if( index(trim(exporttype), 'mic')/=0 .or. &
            & index(trim(exporttype), 'ctf')/=0 ) then
            stype = MIC_STAR
            !         case('stk':'select')           ! micrographs with selected boxes
        else if(index(trim(exporttype), 'select')/=0) then
            stype = STK_STAR
            !         case('ptcl2D':'extract')
        else if(index(trim(exporttype), 'extract')/=0) then
            stype = PTCL_STAR
            !         case('cls2D': 'class2d')             ! class averages
        else if(index(trim(exporttype), 'class2D')/=0) then
            stype = CLSAVG_STAR
            !         case('cls3D': 'init3dmodel')
        else if(index(trim(exporttype), 'init3d')/=0) then
            stype = CLSAVG_STAR
            !         case('ptcl3D':'refine3d')            ! 3D particles
        else if(index(trim(exporttype), 'refine3')/=0) then
            stype = PTCL_STAR
            !         case('out': 'post')
        else if(index(trim(exporttype), 'post')/=0) then
            stype = OUT_STAR
            !         case('projinfo': 'all')
        else if(index(trim(exporttype), 'all')/=0 .or. &
            index(trim(exporttype), 'projinfo')/=0) then
            stype = PROJINFO_STAR
            !         case DEFAULT
        else
            THROW_HARD('unsupported exporttype flag; star project :: exporttype2star')
            !        end select
        end if
    end function exporttype2star


    subroutine print_valid_import_startypes(self)
        class(star_project) :: self
        write (*,*) " import_starproject valid startypes:                                      "
        write (*,*) " Accepted values are m|movies|micrographs ==> import micrographs          "
        write (*,*) "                     ctf|ctf_estimation   ==> import mic + ctf params     "
        write (*,*) "                     p|ptcl|particles|stack   ==> import particles        "
        write (*,*) "                     cavgs|class2D|class3D    ==> import class averages   "
    end subroutine print_valid_import_startypes



    subroutine kill(self,keeptabs)
        class(star_project), intent(inout):: self
        logical, optional :: keeptabs
        call self%doc%kill(keeptabs)
    end subroutine kill

end module simple_star
