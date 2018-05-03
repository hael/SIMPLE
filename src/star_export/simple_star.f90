!! Importing and exporting Relion Star-formatted files to/from SIMPLE
module simple_star
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_sp_project,     only: sp_project
implicit none

type starfile
contains
    procedure :: export_motion_correct
    procedure :: import_motion_correct
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
    procedure :: export_class3D
    procedure :: import_class3D
    procedure :: export_shiny3D
    procedure :: import_shiny3D
end type starfile


contains

!! Motion Correct: Simple/Unblur/MotionCorr2
!! Format of Motion Correct Starfile File:
!! data_
!! loop_
!! _rlnMicrographNameNoDW #1
!! _rlnMicrographName #2
!! [#1 MotionCorr/job026/Micrographs3/FoilHole_24003709_Data_23978423_23978424_20180225_0629-1729_noDW.mrc] [#2 MotionCorr/job026/ Micrographs3/FoilHole_24003709_Data_23978423_23978424_20180225_0629-1729.mrc]
!! ...
subroutine export_motion_correct (self, arg)
    class(starfile), intent(inout) :: self
    real, intent(in) :: arg
end subroutine export_motion_correct
subroutine import_motion_correct (self, arg)
    class(starfile), intent(inout) :: self
    real, intent(in) :: arg
end subroutine import_motion_correct

!! CTF Estimation: Simple/GCTF/CTFFIND4
!! Format of CTF estimation Star file:
!! data_
!! loop_
!! _rlnMicrographNameNoDW #1
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
subroutine export_ctf_estimation (self, arg)
    class(starfile), intent(inout) :: self
    real, intent(in) :: arg
end subroutine export_ctf_estimation
subroutine import_ctf_estimation (self, arg)
    class(starfile), intent(inout) :: self
    real, intent(in) :: arg
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
subroutine export_autopick (self, arg)
    class(starfile), intent(inout) :: self
    real, intent(in) :: arg
end subroutine export_autopick
subroutine import_autopick (self, arg)
    class(starfile), intent(inout) :: self
    real, intent(in) :: arg
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
subroutine export_extract_doseweightedptcls (self, arg)
    class(starfile), intent(inout) :: self
    real, intent(in) :: arg
end subroutine export_extract_doseweightedptcls
subroutine import_extract_doseweightedptcls (self, arg)
    class(starfile), intent(inout) :: self
    real, intent(in) :: arg
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
subroutine export_class2D (self, arg)
    class(starfile), intent(inout) :: self
    real, intent(in) :: arg
end subroutine export_class2D
subroutine import_class2D (self, arg)
    class(starfile), intent(inout) :: self
    real, intent(in) :: arg
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
subroutine export_class2D_select (self, arg)
    class(starfile), intent(inout) :: self
    real, intent(in) :: arg
end subroutine export_class2D_select
subroutine import_class2D_select (self, arg)
    class(starfile), intent(inout) :: self
    real, intent(in) :: arg
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
subroutine export_class3D (self, arg)
    class(starfile), intent(inout) :: self
    real, intent(in) :: arg
end subroutine export_class3D
subroutine import_class3D (self, arg)
    class(starfile), intent(inout) :: self
    real, intent(in) :: arg
end subroutine import_class3D


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
subroutine export_refine3D (self, arg)
    class(starfile), intent(inout) :: self
    real, intent(in) :: arg
end subroutine export_refine3D
subroutine import_refine3D (self, arg)
    class(starfile), intent(inout) :: self
    real, intent(in) :: arg
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
subroutine export_shiny3D (self, arg)
    class(starfile), intent(inout) :: self
    real, intent(in) :: arg
end subroutine export_shiny3D
subroutine import_shiny3D (self, arg)
    class(starfile), intent(inout) :: self
    real, intent(in) :: arg
end subroutine import_shiny3D


end module simple_star
