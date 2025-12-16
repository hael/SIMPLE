module simple_commanders_project_cls
use simple_commander_module_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_import_cavgs
  contains
    procedure :: execute      => exec_import_cavgs
end type commander_import_cavgs

type, extends(commander_base) :: commander_export_cavgs
  contains
    procedure :: execute      => exec_export_cavgs
end type commander_export_cavgs

type, extends(commander_base) :: commander_sample_classes
  contains
    procedure :: execute      => exec_sample_classes
end type commander_sample_classes

contains

    subroutine exec_import_cavgs( self, cline )
        class(commander_import_cavgs), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(simple_nice_communicator) :: nice_communicator
        type(parameters)               :: params
        type(sp_project)               :: spproj
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        ! nice communicator init
        call nice_communicator%init(params%niceprocid, params%niceserver)
        call nice_communicator%cycle()
        if( file_exists(params%projfile) ) call spproj%read(params%projfile)
        call spproj%add_cavgs2os_out(params%stk, params%smpd)
        if( cline%defined('frcs') ) call spproj%add_frcs2os_out(params%frcs,'frc2D')
        call spproj%os_cls2D%set_all2single('state',1.)
        ! update project info
        call spproj%update_projinfo( cline )
        ! update computer environment
        call spproj%update_compenv( cline )
        ! WRITE PROJECT FILE
        call spproj%write ! full write since this is guaranteed to be the first import
        call nice_communicator%terminate()
        call simple_end('**** IMPORT_CAVGS NORMAL STOP ****')
    end subroutine exec_import_cavgs

    subroutine exec_export_cavgs( self, cline )
        class(commander_export_cavgs),   intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(simple_nice_communicator) :: nice_communicator
        type(parameters)               :: params
        type(sp_project)               :: spproj
        type(image)                    :: img
        type(string)                   :: cavgs_fname
        logical,           allocatable :: lstates(:)
        integer :: ldim(3), icls, ncls, ncavgs, cnt
        real    :: smpd, smpd_phys
        call cline%set('oritype', 'cls2D')
        call params%new(cline)
        ! nice communicator init
        call nice_communicator%init(params%niceprocid, params%niceserver)
        call nice_communicator%cycle()
        ! read files and sanity checks
        if( .not.file_exists(params%projfile) ) THROW_HARD('Project file does not exist!')
        call spproj%read_segment(params%oritype,params%projfile)
        if( spproj%os_cls2D%get_noris() == 0 ) THROW_HARD('Absent cls2D field!')
        call spproj%read_segment('out',params%projfile)
        lstates = (spproj%os_cls2D%get_all('state') > 0.5)
        if( count(lstates) == 0 ) THROW_HARD('All class averages are deselected')
        call spproj%get_cavgs_stk(cavgs_fname, ncls, smpd)
        if( spproj%os_cls2D%get_noris() /= ncls ) THROW_HARD('Inconsistent # of entries cls2D/out!')
        call find_ldim_nptcls(cavgs_fname, ldim, ncavgs, smpd=smpd_phys)
        if(ncavgs /= ncls)    THROW_HARD('Inconsistent # of cls2D cavgs & physical cavgs!')
        if( abs(smpd-smpd_phys) > 0.001 ) THROW_HARD('Inconsistent sampling distancs in project & physical cavgs!')
        ! copy selected cavgs
        cnt     = 0
        ldim(3) = 1
        call img%new(ldim, smpd)
        do icls = 1,ncls
            if( lstates(icls) )then
                call img%read(cavgs_fname,icls)
                cnt = cnt + 1
                call img%write(params%outstk,cnt)
            endif
        enddo
        ! the end
        call nice_communicator%terminate()
        call simple_end('**** EXPORT_CAVGS NORMAL STOP ****')
    end subroutine exec_export_cavgs

    subroutine exec_sample_classes( self, cline )
        class(commander_sample_classes), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(simple_nice_communicator)  :: nice_communicator
        type(parameters)                :: params
        type(sp_project)                :: spproj, spproj_part
        integer,            allocatable :: states(:), tmpinds(:), clsinds(:), states_map(:), clustinds(:)
        real,               allocatable :: rstates(:)
        logical,            allocatable :: clustmask(:)
        type(class_sample), allocatable :: clssmp(:)
        type(string)                    :: projfname, fbody
        integer                         :: noris, ncls, icls, iclust, nclust
        integer                         :: nptcls, ipart
        call cline%set('oritype', 'cls2D')
        if( .not. cline%defined('mkdir')           ) call cline%set('mkdir',           'yes')
        if( .not. cline%defined('greedy_sampling') ) call cline%set('greedy_sampling', 'yes')
        if( .not. cline%defined('ranked_parts')    ) call cline%set('ranked_parts',    'yes')
        call params%new(cline, silent=.true.)
         ! nice communicator init
        call nice_communicator%init(params%niceprocid, params%niceserver)
        call nice_communicator%cycle()
        ! read project (almost all or largest segments are updated)
        call spproj%read(params%projfile)
        ! check number of oris in field
        noris   = spproj%get_n_insegment(params%oritype)
        ncls    = spproj%os_cls2D%get_noris()
        nptcls  = spproj%os_ptcl2D%get_noris()
        tmpinds = (/(icls,icls=1,ncls)/)
        rstates = spproj%os_cls2D%get_all('state')
        clsinds = pack(tmpinds, mask=rstates > 0.5)
        allocate(states(nptcls), states_map(nptcls), source=0)
        if( trim(params%partition).eq.'yes' )then
            clustinds = spproj%os_cls2D%get_all_asint('cluster')
            clustinds = pack(clustinds, mask=rstates > 0.5)
            nclust    = maxval(clustinds)
            allocate(clustmask(nclust), source=.false.)
            do iclust = 1, nclust
                if( count(clustinds == iclust) > 0 ) clustmask(iclust) = .true.
            end do
            tmpinds = (/(iclust,iclust=1,nclust)/)
            clsinds = pack(tmpinds, mask=clustmask)
            call spproj%os_ptcl2D%get_class_sample_stats(clsinds, clssmp, label='cluster')
        else
            ncls    = spproj%os_cls2D%get_noris()
            tmpinds = (/(icls,icls=1,ncls)/)
            clsinds = pack(tmpinds, mask=rstates > 0.5)
            call spproj%os_ptcl2D%get_class_sample_stats(clsinds, clssmp)
        endif
        if( cline%defined('nparts') )then
            if( trim(params%ranked_parts).eq.'yes' )then
                if( cline%defined('nptcls_per_part') )then
                    call spproj%os_ptcl2D%sample_ranked_parts(clssmp, params%nparts, states, params%nptcls_per_part)
                else
                    call spproj%os_ptcl2D%sample_ranked_parts(clssmp, params%nparts, states)
                endif
                fbody = RANKPROJPARTFBODY
            else
                if( cline%defined('nptcls_per_part') )then
                    call spproj%os_ptcl2D%sample_balanced_parts(clssmp, params%nparts, states, params%nptcls_per_part)
                else
                    call spproj%os_ptcl2D%sample_balanced_parts(clssmp, params%nparts, states)
                endif
                fbody = BALPROJPARTFBODY
            endif
            do ipart = 1, params%nparts
                ! count # particles in part
                write(logfhandle,*) '# particles in part '//int2str(ipart)//': ', count(states == ipart)
                ! copy project
                projfname = fbody//int2str(ipart)//'.simple'
                call simple_copy_file(params%projfile, projfname)
                call spproj_part%read(projfname)
                ! create state mapping
                where(states == ipart)
                    states_map = 1
                elsewhere
                    states_map = 0
                endwhere
                ! communicate state mapping to copied project
                call spproj_part%os_ptcl2D%set_all('state', real(states_map))
                call spproj_part%os_ptcl3D%set_all('state', real(states_map))
                ! map ptcl states to classes
                call spproj_part%map_ptcls_state_to_cls
                ! write project
                call spproj_part%write(projfname)
                ! destruct
                call spproj_part%kill
            end do
        else
            if( .not. cline%defined('nsample') ) THROW_HARD('Need NSAMPLE at the command line for this mode of sampling')
            if( cline%defined('frac_best') )then
                call spproj%os_ptcl2D%sample_balanced(    clssmp, params%nsample, params%frac_best,     states)
            else if( cline%defined('frac_worst') )then
                call spproj%os_ptcl2D%sample_balanced_inv(clssmp, params%nsample, params%frac_worst,    states)
            else
                call spproj%os_ptcl2D%sample_balanced(    clssmp, params%nsample, params%l_greedy_smpl, states)
            endif
            call spproj%os_ptcl2D%set_all('state', real(states))
            call spproj%os_ptcl3D%set_all('state', real(states))
            call spproj%map_ptcls_state_to_cls
        endif
        ! final full write
        call spproj%write(params%projfile)
        call nice_communicator%terminate()
        call simple_end('**** SAMPLE_CLASSES NORMAL STOP ****')
    end subroutine exec_sample_classes

end module simple_commanders_project_cls
