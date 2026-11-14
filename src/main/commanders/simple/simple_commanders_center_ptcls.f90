!@descr: center_ptcls, centring particles at scale: a sample's few low-resolution class averages, then one pass over every particle against them under the polar shift map's wide prior; only the shifts and their widths are written
! doc/algorithms/center_ptcls.md is the design; the pass leaves every particle's shift, its posterior
! width (shsig) and the move (shincr), so that a later run's dynamic centring searches the window, not the disk.
module simple_commanders_center_ptcls
use simple_commander_base,      only: commander_base
use simple_error,               only: simple_exception
use simple_cmdline,             only: cmdline
use simple_parameters,          only: parameters
use simple_sp_project,          only: sp_project
use simple_string,              only: string
use simple_commanders_refine2D, only: commander_refine2D
implicit none

public :: commander_center_ptcls
private
#include "simple_local_flags.inc"

integer, parameter :: CENTER_NSAMPLE     = 10000 !< the references' sample: enough particles for ten low-resolution classes
integer, parameter :: CENTER_NCLS        = 10    !< ten references centred twice as many trpv1 picks within 3 A as the rotational average (Phase 1)
integer, parameter :: CENTER_REF_MAXITS  = 15    !< the references' iterations: refine2D stops earlier on the sample (7 to 13), more changed nothing
real,    parameter :: CENTER_LP          = 20.   !< the references' band (A)
real,    parameter :: CENTER_LPSTOP      = 12.   !< the pass's band (A): finer than the references', which centred more particles
real,    parameter :: CENTER_SHIFT_SIGMA = 10.   !< the shift prior (A): a 3-sigma disk of 30 A covers the picking error's 95th percentile
character(len=*), parameter :: SAMPLE_PROJ = 'center_sample.simple'
character(len=*), parameter :: PASS_PROJ   = 'center_pass.simple'
character(len=*), parameter :: PASS_DIR    = 'center_pass'        !< the pass runs as a first iteration (every particle not yet centred) in its own folder

type, extends(commander_base) :: commander_center_ptcls
  contains
    procedure :: execute => exec_center_ptcls
end type commander_center_ptcls

contains

    !> The project's active particles centred: references from a random sample (refine2D at lp on the library's crop
    !! for that band, ncls classes, centred), then one refine2D pass of every particle against them at lpstop under
    !! the polar shift map; the shifts and their widths go into ptcl2D and ptcl3D, nothing else changes
    subroutine exec_center_ptcls( self, cline )
        use simple_estimate_ssnr,   only: lpstages_setlims
        use simple_type_defs,       only: lp_crop_inf
        use simple_refine2D_fnames, only: refine2D_cavgs_fname
        use simple_ran_tabu,        only: ran_tabu
        use simple_defs,            only: logfhandle
        use simple_string_utils,    only: int2str
        use simple_jiffys,          only: simple_end
        use simple_syslib,          only: simple_mkdir, simple_chdir, simple_abspath, simple_getcwd
        class(commander_center_ptcls), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters), allocatable :: params
        type(sp_project), allocatable :: spproj, work
        type(commander_refine2D) :: xrefine2D
        type(cmdline)            :: cline_ref
        type(lp_crop_inf)        :: lpinfo(1)
        type(ran_tabu)           :: rt
        type(string)             :: refs, projabs, cwd
        integer, allocatable     :: pinds(:), pinds_sample(:)
        integer :: nptcls, nsample, box_crop, endit, i, iptcl
        real    :: sh(2)
        ! this program's defaults, before the single params%new
        if( .not. cline%defined('ncls')        ) call cline%set('ncls',        CENTER_NCLS)
        if( .not. cline%defined('nsample')     ) call cline%set('nsample',     CENTER_NSAMPLE)
        if( .not. cline%defined('lp')          ) call cline%set('lp',          CENTER_LP)
        if( .not. cline%defined('lpstop')      ) call cline%set('lpstop',      CENTER_LPSTOP)
        if( .not. cline%defined('shift_sigma') ) call cline%set('shift_sigma', CENTER_SHIFT_SIGMA)
        if( .not. cline%defined('mkdir')       ) call cline%set('mkdir',       'yes')
        call cline%set('oritype', 'ptcl2D')
        call cline%set('objfun',  'euclid')
        if( cline%get_rarg('lpstop') >= cline%get_rarg('lp') ) THROW_HARD('the pass band (lpstop) must be finer than the references'' (lp); center_ptcls')
        if( .not. (cline%get_rarg('shift_sigma') > 0.) ) THROW_HARD('center_ptcls centres under the polar shift map: shift_sigma > 0')
        allocate(params, spproj, work)
        call params%new(cline)
        call cline%set('mkdir', 'no')
        call spproj%read(params%projfile)
        call spproj%os_ptcl2D%get_pinds(1, 'state', pinds)
        if( .not. allocated(pinds) ) THROW_HARD('no active particle; center_ptcls')
        nptcls = size(pinds)
        if( nptcls < params%ncls ) THROW_HARD('fewer active particles than classes; center_ptcls')
        ! the library's crop for the references' band (solve3D's rule), never above the native box
        call lpstages_setlims(params%box, 1, params%smpd, params%lp, params%lp, lpinfo, verbose=.false.)
        box_crop = min(params%box, lpinfo(1)%box_crop)
        write(logfhandle,'(A,I0,A,I0,A,F6.2,A)') '>>> CENTER_PTCLS: ', nptcls, ' active particles; working box ', box_crop, &
            &' (', real(params%box) / real(box_crop) * params%smpd, ' A per pixel)'
        ! the sample's project: a random draw of nsample active particles (all when fewer), no 2D clustering
        nsample = min(params%nsample, nptcls)
        pinds_sample = pinds
        if( nsample < nptcls )then
            rt = ran_tabu(nptcls)
            call rt%shuffle(pinds_sample)
            call rt%kill
        endif
        call work%read(params%projfile)
        call work%os_ptcl2D%delete_2Dclustering(keepshifts=.true.)
        do i = 1,nptcls
            call work%os_ptcl2D%set_state(pinds(i), 0)
        enddo
        do i = 1,nsample
            call work%os_ptcl2D%set_state(pinds_sample(i), 1)
        enddo
        call work%update_projinfo(string(SAMPLE_PROJ))
        call work%write(string(SAMPLE_PROJ))
        call work%kill
        write(logfhandle,'(A,I0,A,I0,A,F6.2,A)') '>>> CENTER_PTCLS: references from ', nsample, ' particles, ', params%ncls, &
            &' classes at ', params%lp, ' A'
        call set_child(cline_ref, SAMPLE_PROJ, params%lp)
        call cline_ref%set('cls_init', 'rand')
        call cline_ref%set('maxits',   CENTER_REF_MAXITS)
        call cline_ref%set('center',   'yes')
        call xrefine2D%execute(cline_ref)
        endit = cline_ref%get_iarg('endit')
        refs  = simple_abspath(refine2D_cavgs_fname(endit, params%ext%to_char()))
        call cline_ref%kill
        ! the pass: every active particle of a copy of the project, as the first iteration of its own run in its own
        ! folder, so that every particle starts not yet centred and searches the disk (a later iteration would take
        ! the stored width, zero for raw picks, and search only the window)
        projabs = simple_abspath(params%projfile)
        call simple_getcwd(cwd)
        call simple_mkdir(PASS_DIR)
        call simple_chdir(PASS_DIR)
        call work%read(projabs)
        call work%os_ptcl2D%delete_2Dclustering(keepshifts=.true.)
        call work%update_projinfo(string(PASS_PROJ))
        call work%write(string(PASS_PROJ))
        call work%kill
        write(logfhandle,'(A,I0,A,F6.2,A)') '>>> CENTER_PTCLS: one pass of ', nptcls, ' particles at ', params%lpstop, ' A'
        call set_child(cline_ref, PASS_PROJ, params%lpstop)
        call cline_ref%set('refs',   refs%to_char())
        call cline_ref%set('maxits', 1)
        call cline_ref%set('center', 'no')
        call xrefine2D%execute(cline_ref)
        call cline_ref%kill
        ! the shifts, their widths and moves into the owner's project, in both particle fields
        call work%read_segment('ptcl2D', string(PASS_PROJ))
        call simple_chdir(cwd)
        do i = 1,nptcls
            iptcl = pinds(i)
            sh    = work%os_ptcl2D%get_2Dshift(iptcl)
            call spproj%os_ptcl2D%set_shift(iptcl, sh)
            call spproj%os_ptcl2D%set(iptcl, 'shsig',  work%os_ptcl2D%get(iptcl, 'shsig'))
            call spproj%os_ptcl2D%set(iptcl, 'shincr', work%os_ptcl2D%get(iptcl, 'shincr'))
            if( spproj%os_ptcl3D%get_noris() >= iptcl )then
                call spproj%os_ptcl3D%set_shift(iptcl, sh)
                call spproj%os_ptcl3D%set(iptcl, 'shsig',  work%os_ptcl2D%get(iptcl, 'shsig'))
                call spproj%os_ptcl3D%set(iptcl, 'shincr', work%os_ptcl2D%get(iptcl, 'shincr'))
            endif
        enddo
        call work%kill
        call spproj%write_segment_inside('ptcl2D', params%projfile)
        if( spproj%os_ptcl3D%get_noris() > 0 ) call spproj%write_segment_inside('ptcl3D', params%projfile)
        write(logfhandle,'(A,I0,A)') '>>> CENTER_PTCLS: ', nptcls, ' particles centred'
        call spproj%kill
        deallocate(params, spproj, work)
        call simple_end('**** SIMPLE_CENTER_PTCLS NORMAL STOP ****')

      contains

        !> a refine2D child on a work project at a band: the probabilistic dense table under the polar shift map, on
        !! the library's crop; the caller adds what differs between the references and the pass
        subroutine set_child( cl, projfile, lp )
            type(cmdline),    intent(inout) :: cl
            character(len=*), intent(in)    :: projfile
            real,             intent(in)    :: lp
            call cl%set('prg',         'refine2D')
            call cl%set('projfile',    projfile)
            call cl%set('oritype',     'ptcl2D')
            call cl%set('ncls',        params%ncls)
            call cl%set('refine',      'prob')
            call cl%set('objfun',      'euclid')
            call cl%set('lp',          lp)
            call cl%set('box_crop',    box_crop)
            call cl%set('shift_sigma', params%shift_sigma)
            call cl%set('mskdiam',     params%mskdiam)
            call cl%set('nthr',        params%nthr)
            if( cline%defined('nparts') ) call cl%set('nparts', params%nparts)
            call cl%set('mkdir',       'no')
        end subroutine set_child

    end subroutine exec_center_ptcls

end module simple_commanders_center_ptcls
