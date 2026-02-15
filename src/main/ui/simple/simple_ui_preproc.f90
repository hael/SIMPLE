!@descr: module defining the user interfaces for pre-processing programs in the simple_exec suite
module simple_ui_preproc
use simple_ui_modules
implicit none

type(ui_program), target :: assign_optics_groups
type(ui_program), target :: ctf_estimate
type(ui_program), target :: extract
type(ui_program), target :: gen_pspecs_and_thumbs
type(ui_program), target :: motion_correct
type(ui_program), target :: pick
type(ui_program), target :: preprocess
type(ui_program), target :: reextract

contains

    subroutine construct_preproc_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_assign_optics_groups(prgtab)
        call new_ctf_estimate(prgtab)
        call new_extract(prgtab)
        call new_gen_pspecs_and_thumbs(prgtab)
        call new_motion_correct(prgtab)
        call new_pick(prgtab)
        call new_preprocess(prgtab)
        call new_reextract(prgtab)
    end subroutine construct_preproc_programs

    subroutine print_preproc_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('PRE-PROCESSING:', C_UNDERLINED)
        write(logfhandle,'(A)') assign_optics_groups%name%to_char()
        write(logfhandle,'(A)') ctf_estimate%name%to_char()
        write(logfhandle,'(A)') extract%name%to_char()
        write(logfhandle,'(A)') gen_pspecs_and_thumbs%name%to_char()
        write(logfhandle,'(A)') motion_correct%name%to_char()
        write(logfhandle,'(A)') pick%name%to_char()
        write(logfhandle,'(A)') preprocess%name%to_char()
        write(logfhandle,'(A)') reextract%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_preproc_programs

    subroutine new_assign_optics_groups( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call assign_optics_groups%new(&
        &'assign_optics_groups', &                                              ! name
        &'Assign optics groups',&                                               ! descr_short
        &'is a program to assign optics groups',&                               ! descr long
        &'simple_exec',&                                                        ! executable
        &.true.)                                                                ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! parameter input/output
        call assign_optics_groups%add_input(UI_PARM, 'xmldir', 'dir', 'Directory containing per movie EPU XML files',&
        & 'Directory containing per movie EPY XML files', 'e.g. /data/datasetid/xml', .false., '')
        call assign_optics_groups%add_input(UI_PARM, 'maxpop', 'num', 'Maximum number of movies/micrographs/stacks in each optics group',&
        & 'Maximum number of movies/micrographs/stacks in each optics group', 'e.g. 100', .false., '')
        call assign_optics_groups%add_input(UI_PARM, 'optics_offset', 'num', 'Numbering offset to apply to optics groups',&
        & 'Numbering offset to apply to optics groups. Aids with combining datasets', 'e.g. 10', .false., '')
        call assign_optics_groups%add_input(UI_PARM, 'tilt_thres', 'num', 'Threshold for hierarchical clustering of beamtilts',&
        & 'Threshold for hierarchical clustering of beamtilts', 'e.g 0.05', .false., 0.05)
        call assign_optics_groups%add_input(UI_PARM, 'beamtilt', 'binary', 'Use beamtilts in optics group assignment',&
        &'Use beamtilt values (if found in EPU filenames) during optics group assignment(yes|no){yes}', 'beamtilt(yes|no){yes}', .false., 'yes')
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
        call add_ui_program('assign_optics_groups', assign_optics_groups, prgtab)
    end subroutine new_assign_optics_groups

    subroutine new_ctf_estimate( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call ctf_estimate%new(&
        &'ctf_estimate', &                                                  ! name
        &'CTF parameter fitting',&                                          ! descr_short
        &'is a distributed SIMPLE workflow for CTF parameter fitting',&     ! descr_long
        &'simple_exec',&                                                    ! executable
        &.true.,&                                                           ! requires sp_project
        &gui_advanced=.false., gui_submenu_list = "CTF estimation,compute") ! GUI            
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call ctf_estimate%add_input(UI_PARM, pspecsz, gui_submenu="CTF estimation")
        call ctf_estimate%add_input(UI_PARM, ctfpatch, gui_submenu="CTF estimation")
        ! alternative inputs
        ! <empty>
        ! search controls
        call ctf_estimate%add_input(UI_SRCH, dfmin, gui_submenu="CTF estimation")
        call ctf_estimate%add_input(UI_SRCH, dfmax, gui_submenu="CTF estimation")
        call ctf_estimate%add_input(UI_SRCH, astigtol, gui_submenu="CTF estimation")
        ! filter controls
        call ctf_estimate%add_input(UI_FILT, lp, required_override=.false., gui_submenu="CTF estimation")
        call ctf_estimate%add_input(UI_FILT, hp, required_override=.false., gui_submenu="CTF estimation")
        ! mask controls
        ! <empty>
        ! computer controls
        call ctf_estimate%add_input(UI_COMP, nparts, gui_submenu="compute")
        call ctf_estimate%add_input(UI_COMP, nthr, gui_submenu="compute")
        ! add to ui_hash
        call add_ui_program('ctf_estimate', ctf_estimate, prgtab)
    end subroutine new_ctf_estimate

    subroutine new_extract( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call extract%new(&
        &'extract', &                                                           ! name
        &'Extract particle images from integrated movies',&                     ! descr_short
        &'is a program for extracting particle images from integrated movies',& ! descr long
        &'simple_exec',&                                                        ! executable
        &.true.,&                                                               ! requires sp_project
        &gui_advanced=.false., gui_submenu_list = "extract,compute")            ! GUI      
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call extract%add_input(UI_IMG, 'dir_box', 'dir', 'Box files directory', 'Directory to read the box files from', 'e.g. boxes/', .false., '',&
        &gui_submenu="extract")
        ! parameter input/output
        call extract%add_input(UI_PARM, box, required_override=.false., gui_submenu="extract", gui_advanced=.false.)
        call extract%add_input(UI_PARM, pcontrast, gui_submenu="extract", gui_advanced=.false.)
        call extract%add_input(UI_PARM, outside, gui_submenu="extract")
        call extract%add_input(UI_PARM, backgr_subtr, gui_submenu="extract")
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call extract%add_input(UI_COMP, nparts, gui_submenu="compute", gui_advanced=.false.)
        call extract%add_input(UI_COMP, nthr, gui_submenu="compute", gui_advanced=.false.)
        ! add to ui_hash
        call add_ui_program('extract', extract, prgtab)
    end subroutine new_extract

    subroutine new_gen_pspecs_and_thumbs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call gen_pspecs_and_thumbs%new(&
        &'gen_pspecs_and_thumbs', &                                              ! name
        &'Motion correction of movies',&                                         ! descr_short
        &'is a distributed workflow for generating power spectra and thumbnails&
        & for imported integrated movies',&                                      ! descr_long
        &'simple_exec',&                                                         ! executable
        &.true.)                                                                 ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call gen_pspecs_and_thumbs%add_input(UI_PARM, pspecsz)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call gen_pspecs_and_thumbs%add_input(UI_COMP, nparts)
        call gen_pspecs_and_thumbs%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('gen_pspecs_and_thumbs', gen_pspecs_and_thumbs, prgtab)
    end subroutine new_gen_pspecs_and_thumbs

    subroutine new_motion_correct( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call motion_correct%new(&
        &'motion_correct', &                                                            ! name
        &'Anisotropic motion correction of movies',&                                    ! descr_short
        &'is a distributed workflow for anisotropic motion correction of movies.&
        & If then total dose is given the individual frames will be filtered accordingly&
        & (dose-weighting strategy). If scale is given, the movie will be Fourier cropped according to&
        & the down-scaling factor (for super-resolution movies).',&                     ! descr_long
        &'simple_exec',&                                                                ! executable
        &.true.,&                                                                       ! requires sp_project
        &gui_advanced=.false., gui_submenu_list = "data,motion correction,compute")     ! GUI                    
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call motion_correct%add_input(UI_IMG, gainref, gui_submenu="data", gui_advanced=.false.)
        ! parameter input/output
        call motion_correct%add_input(UI_PARM, total_dose, gui_submenu="data", gui_advanced=.false.)
        call motion_correct%add_input(UI_PARM, fraction_dose_target, gui_submenu="data")
        call motion_correct%add_input(UI_PARM, max_dose, gui_submenu="data")
        call motion_correct%add_input(UI_PARM, smpd_downscale, gui_submenu="data")
        call motion_correct%add_input(UI_PARM, 'fbody', 'string', 'Template output micrograph name',&
        &'Template output integrated movie name', 'e.g. mic_', .false., '', gui_submenu="data")
        call motion_correct%add_input(UI_PARM, pspecsz, gui_submenu="motion correction")
        call motion_correct%add_input(UI_PARM, eer_fraction, gui_submenu="data")
        call motion_correct%add_input(UI_PARM, flipgain, gui_submenu="motion correction")
        ! alternative inputs
        ! <empty>
        ! search controls
        call motion_correct%add_input(UI_SRCH, trs_mc, gui_submenu ="motion correction")
        call motion_correct%add_input(UI_SRCH, 'bfac', 'num', 'B-factor applied to frames', 'B-factor applied to frames (in Angstroms^2)', &
        &'in Angstroms^2{50}', .false., 50., gui_submenu="motion correction")
        call motion_correct%add_input(UI_SRCH, mcpatch, gui_submenu="motion correction")
        call motion_correct%add_input(UI_SRCH, nxpatch, gui_submenu="motion correction")
        call motion_correct%add_input(UI_SRCH, nypatch, gui_submenu="motion correction")
        call motion_correct%add_input(UI_SRCH, mcconvention, gui_submenu="motion correction")
        call motion_correct%add_input(UI_SRCH, algorithm, gui_submenu="motion correction")
        call motion_correct%add_input(UI_SRCH, mcpatch_thres, gui_submenu="motion correction")
        ! filter controls
        call motion_correct%add_input(UI_FILT, 'lpstart', 'num', 'Initial low-pass limit', 'Low-pass limit to be applied in the first &
        &iterations of movie alignment (in Angstroms){8}', 'in Angstroms{8}', .false., 8., gui_submenu="motion correction")
        call motion_correct%add_input(UI_FILT, 'lpstop', 'num', 'Final low-pass limit', 'Low-pass limit to be applied in the last &
        &iterations of movie alignment (in Angstroms){5}', 'in Angstroms{5}', .false., 5., gui_submenu="motion correction")
        call motion_correct%add_input(UI_FILT, wcrit, gui_submenu="motion correction")
        ! mask controls
        ! <empty>
        ! computer controls
        call motion_correct%add_input(UI_COMP, nparts, gui_submenu="compute")
        call motion_correct%add_input(UI_COMP, nthr, gui_submenu="compute")
        ! add to ui_hash
        call add_ui_program('motion_correct', motion_correct, prgtab)
    end subroutine new_motion_correct

    subroutine new_pick( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call pick%new(&
        &'pick', &                                                         ! name
        &'Template-based particle picking',&                               ! descr_short
        &'is a distributed workflow for template-based particle picking',& ! descr_long
        &'simple_exec',&                                                   ! executable
        &.true.,&                                                          ! requires sp_project
        &gui_advanced=.false., gui_submenu_list = "picking,compute")       ! GUI         
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call pick%add_input(UI_IMG, pickrefs, gui_submenu="picking", gui_advanced=.false.)
        ! parameter input/output
        call pick%add_input(UI_PARM, 'dir', 'dir', 'Output directory', 'Output directory', 'e.g. pick/', .false., 'pick', gui_submenu="picking")
        call pick%add_input(UI_PARM, pcontrast, gui_submenu="picking")
        call pick%add_input(UI_PARM, moldiam, gui_submenu="picking")
        call pick%add_input(UI_PARM, picker, gui_submenu="picking")
        call pick%add_input(UI_PARM, 'nmoldiams', 'num', 'Number of molecular diameters to investigate', 'Number of molecular diameters tested', 'e.g. 5',&
        &.false., 5., gui_submenu="picking")
        call pick%add_input(UI_PARM, moldiam_max, gui_submenu="picking")
        call pick%add_input(UI_PARM, 'multi_moldiams', 'str', 'Comma-separated molecular diameters with which to execute multiple gaussian pick ',&
        &'Molecular diameters with which to execulte multiple gaussian pick', 'e.g. 100,150', .false., '', gui_submenu="picking")
        ! alternative inputs
        ! <empty>
        ! search controls
        call pick%add_input(UI_SRCH, 'ndev', 'num', '# of sigmas for outlier detection', '# of standard deviations threshold for outlier detection{2.5}',&
        &'{2.5}', .false., 2.5, gui_submenu="picking", gui_advanced=.false.)
        call pick%add_input(UI_SRCH, pick_roi, gui_submenu="picking")
        call pick%add_input(UI_SRCH, backgr_subtr, gui_submenu="picking")
        call pick%add_input(UI_SRCH, particle_density, gui_submenu="picking")
        call pick%add_input(UI_SRCH, 'winsz', 'num', 'Window size for sauvola', 'Window size for local sauvola binarisation', 'winsz in pixels ', .false., 32.)
        call pick%add_input(UI_SRCH, nboxes_max, gui_submenu="picking")
        ! filter controls
        call pick%add_input(UI_FILT, lp, gui_submenu="picking")
        ! mask controls
        ! <empty>
        ! computer controls
        call pick%add_input(UI_COMP, nparts, gui_submenu="compute", gui_advanced=.false.)
        call pick%add_input(UI_COMP, nthr,   gui_submenu="compute", gui_advanced=.false.)
        ! add to ui_hash
        call add_ui_program('pick', pick, prgtab)
    end subroutine new_pick

    subroutine new_preprocess( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call preprocess%new(&
        &'preprocess', &                                                                    ! name
        &'Preprocessing',&                                                                  ! descr_short
        &'is a distributed workflow that executes motion_correct, ctf_estimate and pick'//& ! descr_long
        &' in sequence',&
        &'simple_exec',&                                                                    ! executable
        &.true.)                                                                            ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call preprocess%add_input(UI_IMG, gainref)
        ! parameter input/output
        call preprocess%add_input(UI_PARM, total_dose)
        call preprocess%add_input(UI_PARM, fraction_dose_target)
        call preprocess%add_input(UI_PARM, max_dose)
        call preprocess%add_input(UI_PARM, smpd_downscale)
        call preprocess%add_input(UI_PARM, eer_fraction)
        call preprocess%add_input(UI_PARM, 'fbody', 'string', 'Template output micrograph name',&
        &'Template output integrated movie name', 'e.g. mic_', .false., 'mic_')
        call preprocess%add_input(UI_PARM, pspecsz)
        call preprocess%add_input(UI_PARM, numlen)
        call preprocess%add_input(UI_PARM, ctfpatch, required_override=.false.)
        call preprocess%add_input(UI_PARM, flipgain)
        ! alternative inputs
        ! <empty>
        ! search controls
        call preprocess%add_input(UI_SRCH, trs_mc)
        call preprocess%add_input(UI_SRCH, dfmin)
        call preprocess%add_input(UI_SRCH, dfmax)
        call preprocess%add_input(UI_SRCH, astigtol)
        call preprocess%add_input(UI_SRCH, 'bfac', 'num', 'B-factor applied to frames', 'B-factor applied to frames (in Angstroms^2)', 'in Angstroms^2{50}', .false., 50.)
        call preprocess%add_input(UI_SRCH, mcpatch)
        call preprocess%add_input(UI_SRCH, nxpatch)
        call preprocess%add_input(UI_SRCH, nypatch)
        call preprocess%add_input(UI_SRCH, mcconvention)
        call preprocess%add_input(UI_SRCH, algorithm)
        call preprocess%add_input(UI_SRCH, mcpatch_thres)
        ! filter controls
        call preprocess%add_input(UI_FILT, 'lpstart', 'num', 'Initial low-pass limit for movie alignment', 'Low-pass limit to be applied in the first &
        &iterations of movie alignment(in Angstroms){8}', 'in Angstroms{8}', .false., 8.)
        call preprocess%add_input(UI_FILT, 'lpstop', 'num', 'Final low-pass limit for movie alignment', 'Low-pass limit to be applied in the last &
        &iterations of movie alignment(in Angstroms){5}', 'in Angstroms{5}', .false., 5.)
        call preprocess%add_input(UI_FILT, 'lp_ctf_estimate', 'num', 'Low-pass limit for CTF parameter estimation',&
        & 'Low-pass limit for CTF parameter estimation in Angstroms{5}', 'in Angstroms{5}', .false., 5.)
        call preprocess%add_input(UI_FILT, 'hp_ctf_estimate', 'num', 'High-pass limit for CTF parameter estimation',&
        & 'High-pass limit for CTF parameter estimation  in Angstroms{30}', 'in Angstroms{30}', .false., 30.)
        ! mask controls
        ! <empty>
        ! computer controls
        call preprocess%add_input(UI_COMP, nparts)
        call preprocess%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('preprocess', preprocess, prgtab)
    end subroutine new_preprocess

    subroutine new_reextract( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call reextract%new(&
        &'reextract', &                                                         ! name
        &'Re-extract particle images from integrated movies',&                  ! descr_short
        &'is a program for re-extracting particle images from integrated movies based on determined 2D/3D shifts',& ! descr long
        &'simple_exec',&                                                        ! executable
        &.true.)                                                                ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call reextract%add_input(UI_PARM, box,     required_override=.false.)
        call reextract%add_input(UI_PARM, oritype, descr_placeholder_override = '(ptcl2D|ptcl3D){ptcl3D}')
        call reextract%add_input(UI_PARM, pcontrast)
        call reextract%add_input(UI_PARM, backgr_subtr)
        call reextract%add_input(UI_PARM, outside)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call reextract%add_input(UI_COMP, nparts)
        call reextract%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('reextract', reextract, prgtab)
    end subroutine new_reextract

end module simple_ui_preproc
