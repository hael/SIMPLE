!@descr: core helpers for SIMPLE parameter defaults and utility procedures
submodule(simple_parameters) simple_parameters_core
implicit none
#include "simple_local_flags.inc"

contains

    module subroutine init_dynamic_defaults(self)
        class(parameters), intent(inout) :: self
        self%boxfile=''           !< file with EMAN particle coordinates(.txt)
        self%boxtab=''            !< table (text file) of files with EMAN particle coordinates(.txt)
        self%ciffile=''           !< PDBx/mmCIF file
        self%class_assignment=''  !< text file listing class ids assigned to a worker
        self%classdoc=''          !< doc with per-class stats(.txt)
        self%cwd=''
        self%deftab=''            !< file with CTF info(.txt|.simple)
        self%deselfile=''         !< file with indices to be deselected(.txt)
        self%dir=''               !< directory
        self%dir_box=''
        self%dir_exec=''          !< name of execution directory
        self%dir_meta=''          !< grab xml files from here
        self%dir_movies=''        !< grab mrc mrcs files from here
        self%dir_preprocess=''    !< grab preprocessed files from here
        self%dir_prev=''          !< grab previous projects for streaming
        self%dir_refine=''        !< refinement directory
        self%dir_reject='rejected'!< move rejected files to here{rejected}
        self%dir_select='selected'!< move selected files to here{selected}
        self%dir_target=''        !< put output here
        self%exec_dir='./'        !< auto-named execution directory
        self%executable=''        !< name of executable
        self%ext=MRC_EXT          !< file extension{.mrc}
        self%fbody=''             !< file body
        self%filetab=''           !< list of files(.txt)
        self%fname=''             !< file name
        self%frcs=trim(FRCS_FILE) !< binary file with per-class/proj Fourier Ring Correlations(.bin)
        self%fsc='fsc_state01.bin'!< binary file with FSC info{fsc_state01.bin}
        self%gainref=''           !< gain reference for movie alignment
        self%import_dir=''        !< dir to import .star files from for import_starproject
        self%infile2=''           !< file with inputs(.txt)
        self%infile=''            !< file with inputs(.txt)
        self%last_prev_dir=''     !< last previous execution directory
        self%msklist=''           !< table (text file) of mask volume files(.txt)
        self%mskvols(MAXS)=''
        self%niceserver=''        !< address and port of nice server for comms
        self%optics_dir=''        !< directory containing stream optics
        self%oritab2=''           !< 2nd table of orientations(.txt|.simple)
        self%oritab=''            !< table  of orientations(.txt|.simple)
        self%outdir=''            !< manually set output directory name
        self%outfile='outfile'//trim(METADATA_EXT) !< output document
        self%outstk=''            !< output image stack
        self%outvol=''            !< output volume{outvol.ext}
        self%pdbfile2=''          !< PDB file, another one
        self%pdbfile=''           !< PDB file
        self%pdbfiles=''          !< list of PDB files
        self%pdbout=''            !< PDB output file
        self%pdfile='pdfile.bin'
        self%pickrefs=''          !< picking references
        self%plaintexttab=''      !< plain text file of input parameters
        self%prg=''               !< SIMPLE program being executed
        self%projfile=''          !< SIMPLE *.simple project file
        self%projfile_orig=''     !< original SIMPLE *.simple project file
        self%projfile_merged=''   !< merged SIMPLE *.simple project file output
        self%projfile_optics=''   !< SIMPLE *.simple project file containing optics group definitions
        self%projfile_ref=''      !< SIMPLE project containing reference assignments
        self%projfile_target=''   !< another SIMPLE *.simple project file
        self%projname=''          !< SIMPLE  project name
        self%projtab=''           !< table of SIMPLE *.simple project files
        self%refs=''              !< initial2Dreferences.ext
        self%refs_even=''
        self%refs_odd=''
        self%snapshot=''          !< path to write snapshot project file to
        self%star_datadir=''      !< STAR-generated data directory
        self%star_mic=''          !< STAR-formatted EM file (micrographs.star)
        self%star_model=''        !< STAR-formatted EM file (model.star)
        self%star_ptcl=''         !< STAR-formatted EM file (data.star)
        self%starfile=''          !< STAR-formatted EM file (proj.star)
        self%stk_den=''           !< denoised particle stack paired with stk
        self%stk2=''              !< 2nd stack(in selection map: selected(cavgs).ext)
        self%stk3=''              !< 3d stack (in selection map (cavgs)2selectfrom.ext)
        self%stk=''               !< particle stack with all images(ptcls.ext)
        self%stk_backgr=''        !< stack with image for background subtraction
        self%stktab_den=''        !< list of denoised per-micrograph stacks paired with stktab
        self%stktab=''            !< list of per-micrograph stacks
        self%subprojname=''       !< SIMPLE  subproject name
        self%vol=''
        self%vol_even=''          !< even reference volume
        self%vol_odd=''           !< odd  reference volume
        self%vols(MAXS)=''
        self%vols_even(MAXS)=''
        self%vols_odd(MAXS)=''
        self%worker_priority=''   !< priority to submit jobs with
        self%worker_server=''     !< address and port of worker server for job submission
        self%xmldir=''
        self%xmlloc=''
    end subroutine init_dynamic_defaults

    module function is_final_planned_iter(self) result(final_planned)
        class(parameters), intent(in) :: self
        logical :: final_planned
        final_planned = (self%which_iter - self%startit + 1) >= self%maxits
    end function is_final_planned_iter

    module subroutine set_img_format(self, ext)
        class(parameters), intent(inout) :: self
        character(len=*),  intent(in)    :: ext
        select case(trim(ext))
            case('M','D','B')
                self%ext = MRC_EXT
            case('S')
                self%ext = '.spi'
            case('J','K','L')
                self%ext = MRC_EXT
            case DEFAULT
                write(logfhandle,*)'format: ', trim(ext)
                THROW_HARD('This file format is not supported by SIMPLE')
        end select
    end subroutine set_img_format

end submodule simple_parameters_core
