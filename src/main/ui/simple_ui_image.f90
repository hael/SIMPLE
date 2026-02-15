!@descr: "image" UI api (concrete implementation)
module simple_ui_image
use simple_ui_modules
implicit none

type(ui_program), target :: binarize
type(ui_program), target :: convert
type(ui_program), target :: ctf_phaseflip
type(ui_program), target :: ctfops
type(ui_program), target :: normalize_
type(ui_program), target :: scale
type(ui_program), target :: stack
type(ui_program), target :: stackops

contains

    subroutine construct_image_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_binarize(prgtab)
        call new_convert(prgtab)
        call new_ctf_phaseflip(prgtab)
        call new_ctfops(prgtab)
        call new_normalize(prgtab)
        call new_scale(prgtab)
        call new_stack(prgtab)
        call new_stackops(prgtab)
    end subroutine construct_image_programs

    subroutine print_image_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('GENERAL IMAGE PROCESSING:', C_UNDERLINED)
        write(logfhandle,'(A)') binarize%name%to_char()
        write(logfhandle,'(A)') convert%name%to_char()
        write(logfhandle,'(A)') ctf_phaseflip%name%to_char()
        write(logfhandle,'(A)') ctfops%name%to_char()
        write(logfhandle,'(A)') normalize_%name%to_char()
        write(logfhandle,'(A)') scale%name%to_char()
        write(logfhandle,'(A)') stack%name%to_char()
        write(logfhandle,'(A)') stackops%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_image_programs

    subroutine new_binarize( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call binarize%new(&
        &'binarize',&                                     ! name
        &'Binarization routines for volumes and stacks',& ! descr_long
        &'Binarization routines for volumes and stacks',& ! descr_long
        &'simple_exec',&                                  ! executable
        &.false.)                                         ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call binarize%add_input(UI_PARM, 'fill_holes', 'binary', 'Fill holes', 'Fill holes(yes|no){no}', '(yes|no){no}', .false., 'no')
        call binarize%add_input(UI_PARM, 'ndev', 'num', 'Binarization threshold', 'Binarization threshold in # sigmas', '# sigmas', .false., 0.)
        call binarize%add_input(UI_PARM, 'winsz', 'num', 'Half-window size', 'Half-window size(in pixels)', 'winsz in pixels', .false., 15.0)
        ! alternative inputs
        call binarize%add_input(UI_ALT, 'vol1', 'file', 'Volume', 'Volume to binarize',&
        & 'input volume e.g. vol.mrc', .false., '')
        call binarize%add_input(UI_ALT, 'stk', 'file', 'Stack', 'Stack to binarize',&
        & 'input stack e.g. imgs.mrcs', .false., '')
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
        call binarize%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('binarize', binarize, prgtab)
    end subroutine new_binarize

    subroutine new_convert( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call convert%new(&
        &'convert',&                                                    ! name
        &'Convert between SPIDER and MRC formats',&                     ! descr_short
        &'is a program for converting between SPIDER and MRC formats',& ! descr_long
        &'simple_exec',&                                                ! executable
        &.false.)                                                       ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call convert%add_input(UI_IMG, outvol)
        call convert%add_input(UI_IMG, outstk)
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        call convert%add_input(UI_ALT, 'vol1', 'file', 'Volume', 'Volume to convert', &
        & 'input volume e.g. vol.spi', .false., '')
        call convert%add_input(UI_ALT, 'stk', 'file', 'Stack', 'Stack to convert',&
        & 'input stack e.g. imgs.spi', .false., '')
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
        ! add to ui_hash
        call add_ui_program('convert', convert, prgtab)
    end subroutine new_convert

    subroutine new_ctf_phaseflip( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call ctf_phaseflip%new(&
        &'ctf_phaseflip', &                                        ! name
        &'CTF phaseflip images in project',&                       ! descr_short
        &'is a program for CTF phase flipping images in project',& ! descr long
        &'simple_exec',&                                           ! executable
        &.true.)                                                   ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call ctf_phaseflip%add_input(UI_PARM, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
        ! add to ui_hash
        call add_ui_program('ctf_phaseflip', ctf_phaseflip, prgtab)
    end subroutine new_ctf_phaseflip

    subroutine new_ctfops( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call ctfops%new(&
        &'ctfops', &                                         ! name
        &'Apply CTF to stacked images',&                     ! descr_short
        &'is a program for applying CTF to stacked images',& ! descr long
        &'simple_exec',&                                     ! executable
        &.false.)                                            ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call ctfops%add_input(UI_IMG, stk)
        call ctfops%add_input(UI_IMG, outstk)
        ! parameter input/output
        call ctfops%add_input(UI_PARM, smpd)
        call ctfops%add_input(UI_PARM, neg)
        call ctfops%add_input(UI_PARM, oritab)
        call ctfops%add_input(UI_PARM, deftab)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call ctfops%add_input(UI_FILT, ctf)
        call ctfops%add_input(UI_FILT, 'bfac', 'num', 'CTF B-factor','B-factor of CTF in Angstroms^2', &
            'B-factor in Angstroms^2(>0.0){0}', .false., 0.)
        ! mask controls
        ! <empty>
        ! computer controls
        call ctfops%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('ctfops', ctfops, prgtab)
    end subroutine new_ctfops

    subroutine new_normalize( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call normalize_%new(&
        &'normalize',&                         ! name
        &'Normalize volume/stack',&            ! descr_short
        &'is a program for normalization of MRC or SPIDER stacks and volumes',&
        &'simple_exec',&                       ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call normalize_%add_input(UI_PARM, smpd)
        call normalize_%add_input(UI_PARM, 'norm',       'binary', 'Normalize',       'Statistical normalization: avg=zero, var=1(yes|no){no}',    '(yes|no){no}', .false., 'no')
        call normalize_%add_input(UI_PARM, 'noise_norm', 'binary', 'Noise normalize', 'Statistical normalization based on background(yes|no){no}', '(yes|no){no}', .false., 'no')
        ! alternative inputs
        call normalize_%add_input(UI_ALT, 'stk',  'file', 'Stack to normalize',  'Stack of images to normalize', 'e.g. imgs.mrc', .false., '')
        call normalize_%add_input(UI_ALT, 'vol1', 'file', 'Volume to normalize', 'Volume to normalize',          'e.g. vol.mrc',  .false., '')
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        call normalize_%add_input(UI_MASK, mskdiam)
        ! computer controls
        call normalize_%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('normalize', normalize_, prgtab)
    end subroutine new_normalize

    subroutine new_scale( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call scale%new(&
        &'scale', &                                                                               ! name
        &'Re-scaling MRC and SPIDER stacks and volumes',&                                         ! descr_short
        &'is a program for re-scaling, clipping and padding MRC and SPIDER stacks and volumes',&  ! descr_long
        &'simple_exec',&                                                                          ! executable
        &.false.)                                                                                 ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call scale%add_input(UI_PARM, smpd)
        call scale%add_input(UI_PARM, 'newbox', 'num', 'Scaled box size', 'Target for scaled box size in pixels', 'new box in pixels', .false., 0.)
        call scale%add_input(UI_PARM, 'scale', 'num', 'Scaling ratio', 'Target box ratio for scaling(0-1+)', '(0-1+)', .false., 1.)
        call scale%add_input(UI_PARM, clip)
        call scale%add_input(UI_PARM, outvol)
        call scale%add_input(UI_PARM, outstk)
        ! alternative inputs
        call scale%add_input(UI_ALT, stk, required_override=.false.)
        call scale%add_input(UI_ALT, 'vol1', 'file', 'Input volume', 'Input volume to re-scale',&
        &'input volume e.g. vol.mrc', .false., '')
        call scale%add_input(UI_ALT, 'filetab', 'file', 'Stacks list',&
        &'List of stacks of images to rescale', 'list input e.g. stktab.txt', .false., '')
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call scale%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('scale', scale, prgtab)
    end subroutine new_scale

    subroutine new_stack( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call stack%new(&
        &'stack',&                     ! name
        &'Stack images',&              ! descr_short
        &'is a program for stacking individual images (list) or multiple stacks into one',& ! descr_long
        &'simple_exec',&               ! executable
        &.false.)                      ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call stack%add_input(UI_IMG, 'filetab', 'file', 'Stacks list',&
        &'List of stacks of images to stack into one', 'list input e.g. stktab.txt', .true., '')
        call stack%add_input(UI_IMG, outstk)
        ! parameter input/output
        call stack%add_input(UI_PARM, smpd)
        call stack%add_input(UI_PARM, clip)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
        ! add to ui_hash
        call add_ui_program('stack', stack, prgtab)
    end subroutine new_stack

    subroutine new_stackops( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call stackops%new(&
        &'stackops',&                                ! name
        &'Standard stack editing',&                  ! descr_short
        &'is a program that provides standard single-particle image processing routines for MRC or SPIDER&
        & stacks. To extract a particular state, give oritab and set state.&
        & To select the fraction of best particles, give oritab&
        & and set frac. State and frac options can be combined.&
        & To apply noise, give the desired signal-to-noise ratio via snr. To calculate the autocorrelation&
        & function, set acf=yes. To extract a contiguous subset of images from the stack, set&
        & fromp and top. To extract a number of particle images at random, set nran to the desired number.&
        & With avg=yes the global average of the stack is calculated.&
        & If nframesgrp is set to some integer number >1, averages with chunk sizes of nframesgrp are produced,&
        & which may be useful for analysis of dose-fractionated image series. neg inverts the contrast of the images',& ! descr_long
        &'simple_exec',&                             ! executable
        &.false.)                                    ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call stackops%add_input(UI_IMG, stk, required_override=.true.)
        call stackops%add_input(UI_IMG, outstk)
        ! parameter input/output
        call stackops%add_input(UI_PARM, smpd)
        call stackops%add_input(UI_PARM, oritab)
        call stackops%add_input(UI_PARM, mirr)
        call stackops%add_input(UI_PARM, nran)
        call stackops%add_input(UI_PARM, frac)
        call stackops%add_input(UI_PARM, 'state', 'num', 'State index', 'Index of state to extract', 'give state index', .false., 1.)
        call stackops%add_input(UI_PARM, 'class', 'num', 'Class index', 'Index of class to extract', 'give class index', .false., 1.)
        call stackops%add_input(UI_PARM, neg)
        call stackops%add_input(UI_PARM, 'acf',   'binary', 'Autocorrelation, A * conjg(A)', 'Generate autocorrelation function: A * conjg(A)(yes|no){no}', '(yes|no){no}', .false., 'no')
        call stackops%add_input(UI_PARM, 'avg',   'binary', 'Average stack', 'Generate global stack average(yes|no){no}', '(yes|no){no}', .false., 'no')
        call stackops%add_input(UI_PARM, 'nframesgrp', 'num', 'Number of stack entries to group & average', 'Number of stack entries to group & average{0}', '# frames', .false., 0.)
        call stackops%add_input(UI_PARM, 'vis',   'binary', 'Visualize stack images', 'Visualize stack images with gnuplot(yes|no){no}', '(yes|no){no}', .false., 'no')
        call stackops%add_input(UI_PARM, 'snr',   'num', 'Apply noise to give SNR', 'Apply noise to give this signal-to-noise ratio of output', 'signal-to-noise ratio(0.)', .false., 0.)
        call stackops%add_input(UI_PARM, 'fromp', 'num', 'From particle index', 'Start index for stack copy', 'start index', .false., 1.0)
        call stackops%add_input(UI_PARM, 'top',   'num', 'To particle index', 'Stop index for stack copy', 'stop index', .false., 1.0)
        call stackops%add_input(UI_PARM, outfile)
        call stackops%add_input(UI_PARM, 'stats', 'binary', 'Provide statistics', 'Provide statistics about images in stack(yes|no){no}', '(yes|no){no}', .false., 'no')
        call stackops%add_input(UI_PARM, 'roavg', 'binary', 'Rotationally average', 'Rotationally average images in stack(yes|no){no}', '(yes|no){no}', .false., 'no')
        call stackops%add_input(UI_PARM, 'angstep', 'num', 'Angular stepsize', 'Angular stepsize for rotational averaging(in degrees)', 'give degrees', .false., 5.)
        call stackops%add_input(UI_PARM, 'makemovie', 'binary', 'Whether to make a movie', 'Generates images and script to make a movie with FFmpeg(yes|no){no}', '(yes|no){no}', .false., 'no')
        call stackops%add_input(UI_PARM, 'winsz', 'num', 'Window size for local sdev estimation', 'Window size for local sdev estimation(in pixels)', 'winsz in pixels', .false., 5.0)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call stackops%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('stackops', stackops, prgtab)
    end subroutine new_stackops

end module simple_ui_image
