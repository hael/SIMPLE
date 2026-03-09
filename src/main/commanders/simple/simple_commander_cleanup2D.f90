!@descr: cleanup for molecule with large heterogeneity 
module simple_commanders_cleanup2D
use simple_commanders_api
use simple_commanders_project_core, only: commander_extract_subproj
use simple_commanders_abinitio2D,   only: commander_abinitio2D
use simple_commanders_cavgs,        only: commander_cluster_cavgs
use simple_projfile_utils,          only: merge_chunk_projfiles
use simple_imgarr_utils,            only: read_cavgs_into_imgarr, dealloc_imgarr
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_cleanup2D
    contains
        procedure :: execute      => exec_cleanup2D
end type commander_cleanup2D

contains 
    subroutine exec_cleanup2D( self, cline )
        class(commander_cleanup2D), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        integer, parameter :: NCLS_COARSE   = 30
        integer, parameter :: NCLS_FINE_MAX = 300, NCLS_FINE_MIN = 6
        real,    parameter :: LP_COARSE     = 6.
        type(cmdline)                   :: cline_coarse_abinitio2D, cline_subproj_abinitio2D
        type(cmdline)                   :: cline_new_subproject
        type(commander_abinitio2D)      :: xabinitio2D
        ! type(commander_cluster_cavgs)   :: xcluster_cavgs
        type(image),  allocatable       :: cavgs_coarse(:)
        real,         allocatable       :: diams(:), shifts(:,:), state_labels(:)
        integer,      allocatable       :: class_labels(:)
        type(sp_project)     :: spproj, spproj_sub
        type(parameters)     :: params
        type(string)         :: projname_sub, projfile_sub
        integer :: icls_coarse, ncls_sub, nptcls_sub, box_coarse
        real    :: mskdiam_sub, moldiam_sub
        if( .not. cline%defined('nptcls_per_cls') ) call cline%set('nptcls_per_cls', 500)
        if( .not. cline%defined('nparts')         ) THROW_HARD('NPARTS need to be defined')
        if( .not. cline%defined('nthr')           ) THROW_HARD('NTHR need to be defined')
        if( .not. cline%defined('mkdir')          ) call cline%set('mkdir', 'yes')
        ! master parameters
        call params%new(cline)
        call cline%set('mkdir', 'no')
        call spproj%read(params%projfile)
        ! do first pass of coarse ab initio 2D
        cline_coarse_abinitio2D = cline
        call cline_coarse_abinitio2D%set('mskdiam',          0.)
        call cline_coarse_abinitio2D%set('ncls',    NCLS_COARSE)
        call cline_coarse_abinitio2D%set('lpstop',    LP_COARSE)
        call cline_coarse_abinitio2D%set('center',        'yes')
        call cline_coarse_abinitio2D%set('cls_init',     'rand')
        call cline_coarse_abinitio2D%printline
        call xabinitio2D%execute(cline_coarse_abinitio2D)
        call spproj%read(params%projfile)

        print *, 'read projfile: ', params%projfile%to_char()

        ! estimate per coarse-class mask diameters
        cavgs_coarse = read_cavgs_into_imgarr(spproj)
        call automask2D(params, cavgs_coarse, params%ngrow, nint(params%winsz), params%edge, diams, shifts)
        call dealloc_imgarr(cavgs_coarse)
        ! loop to:
        ! (1) estimate subproject mask diameter from coarse class average
        ! (2) extract subproject from coarse clustering (class index)
        ! (3) execute distributed ab initio 2D on subproject with dynamically estimated ncls
        do icls_coarse = 1, NCLS_COARSE
            ! (1) estimate subproject mask diameter from coarse class average
            box_coarse  = round2even(diams(icls_coarse) / params%smpd + 2. * COSMSKHALFWIDTH)
            moldiam_sub = params%smpd * real(box_coarse) 
            mskdiam_sub = moldiam_sub * MSK_EXP_FAC

            print *, 'coarse class index: ', icls_coarse, ' estimated mask diameter: ', mskdiam_sub
            
            ! (2) extract subproject from coarse clustering (class index)
            projname_sub = 'subproj'//int2str_pad(icls_coarse,2)
            projfile_sub = 'subproj'//int2str_pad(icls_coarse,2)//'.simple'
            ! get class labels
            class_labels = spproj%os_ptcl2D%get_all_asint('class')
            ! update project info
            call cline_new_subproject%set('projname', projname_sub) 
            call spproj_sub%update_projinfo(cline_new_subproject)
            ! update computer environment
            call spproj_sub%update_compenv(cline_new_subproject)
            ! copy relevant project fields
            spproj_sub%os_stk    = spproj%os_stk
            spproj_sub%os_ptcl2D = spproj%os_ptcl2D
            spproj_sub%os_ptcl3D = spproj%os_ptcl3D            
            ! compress ptcl 2D & 3D fields
            allocate(state_labels(size(class_labels)))
            where(class_labels == icls_coarse)
                state_labels = 1.0
            else where
                state_labels = 0.0
            endwhere
            call spproj_sub%os_ptcl2D%set_all('state', state_labels)
            call spproj_sub%os_ptcl3D%set_all('state', state_labels)
            call spproj_sub%prune_particles
            deallocate(state_labels)
            ! write subproject file
            call spproj_sub%write 



            ! (3) execute distributed ab initio 2D on subproject with dynamically estimated ncls

            nptcls_sub = spproj_sub%get_nptcls()

            print *, 'nptcls_sub: ', nptcls_sub

            ncls_sub   = max(NCLS_FINE_MIN, min(NCLS_FINE_MAX, nint(real(nptcls_sub) / real(params%nptcls_per_cls))))
            
            print *, 'ncls_sub: ', ncls_sub

            ! setup command line
            cline_subproj_abinitio2D = cline
            call cline_subproj_abinitio2D%set('prg',      'abinitio2D')
            call cline_subproj_abinitio2D%set('projfile', projfile_sub)
            call cline_subproj_abinitio2D%set('ncls',         ncls_sub)
            call cline_subproj_abinitio2D%set('mskdiam',   mskdiam_sub)
            call cline_subproj_abinitio2D%set('center',          'yes')
            call cline_subproj_abinitio2D%set('cls_init',       'rand')
            call cline_subproj_abinitio2D%set('mkdir',           'yes')

            call cline_subproj_abinitio2D%printline

            call xabinitio2D%execute(cline_subproj_abinitio2D)
            call spproj_sub%kill
            call simple_chdir('../')
        end do

        call simple_end('**** SIMPLE_CLEANUP2D NORMAL STOP ****', print_simple = .false.)
    end subroutine exec_cleanup2D

end module simple_commanders_cleanup2D
    
