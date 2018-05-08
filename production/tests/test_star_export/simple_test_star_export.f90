program simple_test_star_export
    include 'simple_lib.f08'
    use simple_star

use simple_cmdline,        only: cmdline
use simple_params,         only: params
use simple_build,          only: build
use simple_sp_project,     only: sp_project
    implicit none

    type(starfile) :: s
    type(params) :: p

    call s%export_motion_corrected_micrographs (p, trim('tmp_mc.star'))
    call exec_cmdline( 'relion_star_loopheader rlnMicrographNameNoDW rlnMicrographName > tmp_mc.star')

! Generate STAR files from separate stacks for each micrograph
! If the input images are in a separate stack for each micrograph, then one could use the following commands to generate the input STAR file:
call exec_cmdline( ' relion_star_loopheader rlnImageName rlnMicrographName rlnDefocusU rlnDefocusV rlnDefocusAngle rlnVoltage rlnSphericalAberration rlnAmplitudeContrast > my_images.star')
call exec_cmdline( ' relion_star_datablock_stack 4 mic1.mrcs mic1.mrcs 10000 10500 30 200 2 0.1  >> my_images.star')
call exec_cmdline( 'relion_star_datablock_stack 3 mic2.mrcs mic2.mrcs 21000 20500 25 200 2 0.1  >> my_images.star')
call exec_cmdline( ' relion_star_datablock_stack 2 mic3.mrcs mic3.mrcs 16000 15000 35 200 2 0.1  >> my_images.star')



! Generate STAR files from particles in single-file format
!call exec_cmdline( ' relion_star_loopheader rlnImageName rlnMicrographName rlnDefocusU rlnDefocusV rlnDefocusAngle rlnVoltage rlnSphericalAberration rlnAmplitudeContrast > my_images.star
! call exec_cmdline( 'relion_star_datablock_singlefiles "mic1/*.spi" mic1 16000 15000 35 200 2 0.1  >> my_images.star
!call exec_cmdline( ' relion_star_datablock_singlefiles "mic2/*.spi" mic2 16000 15000 35 200 2 0.1  >> my_images.star
!call exec_cmdline( ' relion_star_datablock_singlefiles "mic3/*.spi" mic3 16000 15000 35 200 2 0.1  >> my_images.star

! Generate STAR files from XMIPP-style CTFDAT files
! To generate a STAR file from an XMIPP-style ctfdat file, one could use:

!call exec_cmdline( ' relion_star_loopheader rlnImageName rlnMicrographName rlnDefocusU rlnDefocusV rlnDefocusAngle rlnVoltage rlnSphericalAberration rlnAmplitudeContrast > all_images.star
!call exec_cmdline( ' relion_star_datablock_ctfdat all_images.ctfdat>>  all_images.star

! Generate STAR files from FREALIGN-style .par files
! To generate a STAR file from a FREALIGN-style .par file, one could use:

! call exec_cmdline( 'relion_star_loopheader rlnImageName rlnMicrographName rlnDefocusU rlnDefocusV rlnDefocusAngle rlnVoltage rlnSphericalAberration rlnAmplitudeContrast > all_images.star
! call exec_cmdline( 'awk '{if ($1!="C") {print $1"@./my/abs/path/bigstack.mrcs", $8, $9, $10, $11, " 80 2.0 0.1"}  }' < frealign.par >> all_images.star ')
! Assuming the voltage is 80kV, the spherical aberration is 2.0 and the amplitude contrast is 0.1. Also, a single stack is assumed called: /my/abs/path/bigstack.mrcs.



end program simple_test_star_export
