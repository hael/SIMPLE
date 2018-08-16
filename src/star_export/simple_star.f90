!! Importing and exporting Relion Star-formatted files to/from SIMPLE
module simple_star
include 'simple_lib.f08'
use simple_stardoc
use simple_sp_project, only: sp_project
use simple_cmdline,        only: cmdline
use simple_parameters
use simple_binoris_io
use simple_oris, only: oris
implicit none
private
public :: star_project


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
    procedure :: kill
end type star_project

interface star_project
    module procedure constructor
end interface star_project
#include "simple_local_flags.inc"
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
            write(*,*) 'file: ', trim(filename)
            stop 'file does not exist in cwd;  simple_star :: prepare '
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
        if(.not. self%doc%existence) &
            stop ' simple_star ::  get_ndatalines doc unopened '

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
        if(.not. self%doc%existence) &
            stop ' simple_star ::  get_ndatalines doc unopened '

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
            HALT_NOW( 'star_project;  oritab-stardoc.txt not found' )
        endif
        nl_deftab = nlines('oritab-stardoc.txt')
        if(nl_deftab /= self%doc%num_data_lines)then
            print *,"star_project; check_temp_files called from "//trim(msg)
            HALT_NOW('star_project;  oritab-stardoc.txt lines not same as num_data_lines')
        endif
        if (.not. file_exists('filetab-stardoc.txt')) then
            print *,"star_project; check_temp_files called from "//trim(msg)
            HALT_NOW('star_project;  filetab-stardoc.txt not found')
        endif
        nl_filetab= nlines('filetab-stardoc.txt')
        if(nl_deftab /= self%doc%num_data_lines)then
            print *,"star_project; check_temp_files called from "//trim(msg)
            HALT_NOW('star_project;  oritab-stardoc.txt lines not same as num_data_lines')
        endif
        if(nl_filetab /= nl_deftab)then
            print *,"star_project; check_temp_files called from "//trim(msg)
            HALT_NOW('star_project;  oritab-stardoc.txt lines not same as num_data_lines')
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
        integer         :: ndatlines,noris,iseg, n_ori_inputs
        type(ctfparams) :: ctfvars
        logical         :: inputted_oritab, inputted_plaintexttab, inputted_deftab
        character(len=:), allocatable ::  oritype
        class(oris),      pointer     :: os_ptr
        call cline%set('oritype','mic')
        !! make sure starfile has been parsed and temporary files are in current dir        
        call self%check_temp_files('import_ctf_estimation')
        !! set deftab
        call cline%set('deftab', 'oritab-stardoc.txt')
        if(.not. (cline%defined('stk') .or. cline%defined('stktab')) ) then
            call cline%set('filetab', 'filetab-stardoc.txt')
        endif

  
        ! if( n_ori_inputs > 1 )then
        !     write(*,*) 'ERROR, multiple parameter sources inputted, please use (oritab|deftab|plaintexttab)'
        !     stop 'commander_project :: exec_import_particles'
        ! endif

        call params%new(cline)

        ndatlines = binread_nlines(params%oritab)
        if(ndatlines /= self%doc%num_data_lines)then
            HALT_NOW('star_project; import_ctf_estimation binread_nlines does not match num_data_lines from starfile')
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
        !     write(*,*) 'Project file: ', trim(params%projfile), ' does not exists!'
        !     stop 'ERROR! simple_commander_project :: exec_import_movies'
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
        !         write(*,*) 'ctf flag params%ctf: ', params%ctf
        !         stop 'ERROR! ctf flag not supported; commander_project :: import_movies'
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
        ! !         write(*,*) '# boxfiles: ', nboxf
        ! !         write(*,*) '# movies  : ', nmovf
        ! !         stop 'ERROR! # boxfiles .ne. # movies; commander_project :: exec_import_movies'
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


    function exporttype2star( exporttype ) result(stype)
        character(len=*),  intent(in) :: exporttype
        integer(kind=kind(MIC_STAR))  :: stype
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
            stop 'unsupported exporttype flag; star project :: exporttype2star'
            !        end select
        end if
    end function exporttype2star


    subroutine kill(self,keeptabs)
        class(star_project), intent(inout):: self
        logical, optional :: keeptabs
        call self%doc%kill(keeptabs)
    end subroutine kill

end module simple_star
