!@descr: fractional-update sampling units (balance=class|cavg): one unit per selected 2D class, cavg groups by class-average correlation
! Under cavg the selected class averages are clustered into nclust groups by average linkage on their aligned
! correlation; with no more selected classes than nclust every class is its own group. The sampler never uses
! 3D maps, poses or projection directions.
module simple_view_partition_sampling
use simple_pftc_srch_api
use simple_clustering_utils, only: cluster_dmat
use simple_strategy2D_utils, only: prep_cavgs4clust, calc_cluster_cavgs_dmat
use simple_oris,             only: class_sample_quotas, class_sample_sweep
implicit none

public :: make_class_samples, report_class_sample_coverage, has_selected_cavgs, VIEW_PARTITION_FILE
private
#include "simple_local_flags.inc"

character(len=*), parameter :: VIEW_PARTITION_FILE  = 'view_partition.txt'
character(len=*), parameter :: VIEW_PARTITION_FBODY = 'view_partition'
character(len=*), parameter :: VIEW_PARTITION_CRIT  = 'cc' !< in-plane, shift and mirror invariant correlation
real,             parameter :: VIEW_PARTITION_LP    = 6.  !< low-pass limit of the class-average alignment (as cluster_cavgs)
real,             parameter :: VIEW_PARTITION_TRS   = 10. !< shift search range of the class-average alignment (as cluster_cavgs)
real,             parameter :: VISITS_WARN_FAC      = 10. !< warn when a unit's particles are visited this many times the target

contains

    !> whether the project carries selected class averages (cls2D state > 0 and a class-average stack)
    logical function has_selected_cavgs( spproj ) result( l_has )
        class(sp_project), intent(inout) :: spproj
        type(string) :: stkname
        integer      :: ncls
        real         :: smpd
        l_has = .false.
        if( spproj%os_cls2D%get_noris() < 1 ) return
        if( .not. any(spproj%os_cls2D%get_all_asint('state') > 0) ) return
        call spproj%get_cavgs_stk(stkname, ncls, smpd, fail=.false.)
        l_has = ncls > 0 .and. file_exists(stkname)
        call stkname%kill
    end function has_selected_cavgs

    !> One class_sample per selected 2D class, best ptcl2D score first, for balance (class|cavg): group 0
    !! under class, the class-average group under cavg. l_drop_inactive3D leaves out rows inactive in
    !! ptcl3D (a workflow whose ptcl3D states are its own). The cavg clustering runs on nthr threads.
    subroutine make_class_samples( params, spproj, balance, nthr, l_drop_inactive3D, clssmp )
        class(parameters),               intent(in)    :: params
        class(sp_project),               intent(inout) :: spproj
        character(len=*),                intent(in)    :: balance
        integer,                         intent(in)    :: nthr
        logical,                         intent(in)    :: l_drop_inactive3D
        type(class_sample), allocatable, intent(inout) :: clssmp(:)
        integer, allocatable :: clsinds(:), labels(:)
        integer :: i
        select case(trim(balance))
            case('class', 'cavg')
            case DEFAULT
                THROW_HARD('class sampling units need balance=class|cavg, not '//trim(balance))
        end select
        if( spproj%os_cls2D%get_noris() < 1 ) THROW_HARD('balance='//trim(balance)//' needs a cls2D segment')
        clsinds = pack([(i, i=1,spproj%os_cls2D%get_noris())], mask=spproj%os_cls2D%get_all_asint('state') > 0)
        if( size(clsinds) < 1 ) THROW_HARD('balance='//trim(balance)//' needs selected classes (cls2D state > 0)')
        if( allocated(clssmp) ) call deallocate_class_samples(clssmp)
        call spproj%os_ptcl2D%get_class_sample_stats(clsinds, clssmp)
        if( trim(balance) == 'cavg' )then
            call cluster_selected_cavgs(params, spproj, nthr, clsinds, labels)
            do i = 1, size(clssmp)
                clssmp(i)%group = labels(i)
            end do
        endif
        if( l_drop_inactive3D ) call drop_inactive_rows(spproj, clssmp)
    end subroutine make_class_samples

    !> class-average group of every selected class (clsinds order); every class is its own group when
    !! there are no more selected classes than nclust
    subroutine cluster_selected_cavgs( params, spproj, nthr, clsinds, labels )
        class(parameters),    intent(in)    :: params
        class(sp_project),    intent(inout) :: spproj
        integer,              intent(in)    :: nthr
        integer,              intent(in)    :: clsinds(:)
        integer, allocatable, intent(inout) :: labels(:)
        type(cmdline)            :: cline_clust
        type(parameters)         :: params_clust
        type(image), allocatable :: cavg_imgs(:)
        real,        allocatable :: mm(:,:), dmat(:,:)
        integer,     allocatable :: states(:), clspops(:), clsinds_clust(:), i_medoids(:)
        logical,     allocatable :: l_sel(:)
        integer(timer_int_kind) :: t_start
        integer :: nsel, nclust, ldim(3), i
        t_start = tic()
        if( params%nclust < 1 ) THROW_HARD('balance=cavg needs nclust >= 1')
        if( .not. has_selected_cavgs(spproj) ) THROW_HARD('balance=cavg needs selected class averages')
        states = spproj%os_cls2D%get_all_asint('state')
        l_sel  = states > 0
        nsel   = count(l_sel)
        if( nsel <= params%nclust )then
            labels = [(i, i=1,nsel)]
            write(logfhandle,'(A,I0,A,I0,A)') '>>> CLASS-AVERAGE GROUPS: one group per class (', nsel, &
                &' selected classes <= nclust=', params%nclust, ')'
        else
            ! the class averages are aligned and compared under their own parameters: objfun=cc, no CTF
            call cline_clust%set('prg',      'cluster_cavgs')
            call cline_clust%set('projfile', params%projfile)
            call cline_clust%set('mskdiam',  params%mskdiam)
            call cline_clust%set('nthr',     max(1, nthr))
            call cline_clust%set('oritype',  'cls2D')
            call cline_clust%set('ctf',      'no')
            call cline_clust%set('objfun',   'cc')
            call cline_clust%set('mkdir',    'no')
            call cline_clust%set('lp',       VIEW_PARTITION_LP)
            call cline_clust%set('trs',      VIEW_PARTITION_TRS)
            ! sets the OpenMP team and nthr_glob to nthr until restored below
            call params_clust%new(cline_clust, silent=.true.)
            call prep_cavgs4clust(spproj, cavg_imgs, params_clust%mskdiam, clspops, clsinds_clust, l_sel, mm)
            ldim              = cavg_imgs(1)%get_ldim()
            params_clust%smpd = cavg_imgs(1)%get_smpd()
            params_clust%box  = ldim(1)
            params_clust%msk  = min(real(params_clust%box/2) - COSMSKHALFWIDTH - 1., &
                &0.5 * params_clust%mskdiam / params_clust%smpd)
            dmat = calc_cluster_cavgs_dmat(params_clust, cavg_imgs, [minval(mm(:,1)), maxval(mm(:,2))], VIEW_PARTITION_CRIT)
            nclust = params%nclust
            call cluster_dmat(dmat, 'avglink', nclust, i_medoids, labels)
            call dealloc_imgarr(cavg_imgs)
            call cline_clust%kill
            !$ call omp_set_num_threads(params%nthr)
            nthr_glob = params%nthr
            if( any(clsinds_clust /= clsinds) ) THROW_HARD('class-average clustering changed the class order')
            write(logfhandle,'(A,I0,A,I0,A,I0,A,F8.1,A)') '>>> CLASS-AVERAGE GROUPS: ', nsel, ' class averages in ', &
                &nclust, ' groups by average linkage on correlation (', max(1, nthr), ' threads, ', toc(t_start), ' s)'
        endif
        ! inspection output
        cavg_imgs = read_cavgs_into_imgarr(spproj, mask=l_sel)
        call write_imgarr(nsel, cavg_imgs, labels, VIEW_PARTITION_FBODY, params%ext%to_char())
        call dealloc_imgarr(cavg_imgs)
        call write_partition_table
        if( allocated(clsinds_clust) ) deallocate(clsinds_clust)

    contains

        subroutine write_partition_table
            integer :: funit, io_stat, k
            call fopen(funit, file=string(VIEW_PARTITION_FILE), status='REPLACE', action='WRITE', iostat=io_stat)
            call fileiochk('cluster_selected_cavgs; '//VIEW_PARTITION_FILE, io_stat)
            write(funit,'(A)') '# class group'
            do k = 1, nsel
                write(funit,'(I6,1X,I4)') clsinds(k), labels(k)
            enddo
            call fclose(funit)
        end subroutine write_partition_table

    end subroutine cluster_selected_cavgs

    !> only active particles are sampled: rows with ptcl3D state 0 leave their unit
    subroutine drop_inactive_rows( spproj, clssmp )
        class(sp_project),  intent(inout) :: spproj
        type(class_sample), intent(inout) :: clssmp(:)
        integer, allocatable :: states(:)
        logical, allocatable :: l_keep(:)
        integer :: i
        if( spproj%os_ptcl3D%get_noris() /= spproj%os_ptcl2D%get_noris() ) return
        states = spproj%os_ptcl3D%get_all_asint('state')
        do i = 1, size(clssmp)
            if( .not. allocated(clssmp(i)%pinds) ) cycle
            l_keep = states(clssmp(i)%pinds) > 0
            if( all(l_keep) ) cycle
            clssmp(i)%pinds = pack(clssmp(i)%pinds, mask=l_keep)
            clssmp(i)%ccs   = pack(clssmp(i)%ccs,   mask=l_keep)
            clssmp(i)%pop   = size(clssmp(i)%pinds)
        end do
    end subroutine drop_inactive_rows

    !> The unit table, printed once before the first stage: per group and unit the population, the
    !! per-draw quota and the expected visits per particle over nplanned draws (iterations, or frequency
    !! blocks under cohorts), the draws one sweep needs, and a warning on short or excessive coverage.
    subroutine report_class_sample_coverage( clssmp, balance, ntarget, nplanned, draw_label )
        type(class_sample), intent(in) :: clssmp(:)
        character(len=*),   intent(in) :: balance, draw_label
        integer,            intent(in) :: ntarget, nplanned
        real,    allocatable :: quotas(:), visits(:), keys(:)
        integer, allocatable :: order(:)
        type(string) :: msg
        integer :: nunits, ngroups, nactive, sweep, i, k
        real    :: target, vmin, vmax
        nunits  = size(clssmp)
        if( nunits < 1 ) return
        nactive = sum(clssmp(:)%pop)
        allocate(quotas(nunits), visits(nunits), keys(nunits))
        call class_sample_quotas(clssmp, ntarget, quotas)
        sweep  = class_sample_sweep(clssmp, ntarget)
        visits = 0.
        where( clssmp(:)%pop > 0 ) visits = real(nplanned) * quotas / real(max(1, clssmp(:)%pop))
        target = real(nplanned) * real(min(ntarget, nactive)) / real(max(1, nactive))
        ngroups = count(clssmp(:)%group == 0)
        do i = 1, nunits
            if( clssmp(i)%group > 0 )then
                if( findloc(clssmp(:i-1)%group, clssmp(i)%group, dim=1) == 0 ) ngroups = ngroups + 1
            endif
        end do
        ! group order: grouped units by group index, then the units that are their own group
        order = [(i, i=1,nunits)]
        do i = 1, nunits
            if( clssmp(i)%group > 0 )then
                keys(i) = real(clssmp(i)%group) * real(nunits + 1) + real(i)
            else
                keys(i) = real(maxval(clssmp(:)%group) + 1) * real(nunits + 1) + real(i)
            endif
        end do
        call hpsort(keys, order)
        vmin = minval(visits, mask=clssmp(:)%pop > 0)
        vmax = maxval(visits, mask=clssmp(:)%pop > 0)
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A,A,A)')       '>>> FRACTIONAL-UPDATE SAMPLING UNITS (balance=', trim(balance), ')'
        write(logfhandle,'(A,I10)')       '    Active particles in units      : ', nactive
        write(logfhandle,'(A,I10)')       '    Particles drawn per draw       : ', ntarget
        write(logfhandle,'(A,I10,A,I0)')  '    Groups / units                 : ', ngroups, ' / ', nunits
        write(logfhandle,'(A,I10,A)')     '    Draws for one full sweep       : ', sweep, ' ('//draw_label//'s)'
        write(logfhandle,'(A,I10,A)')     '    Draws planned                  : ', nplanned, ' ('//draw_label//'s)'
        write(logfhandle,'(A,F10.2)')     '    Target visits per particle     : ', target
        write(logfhandle,'(A,2F10.2)')    '    Visits per particle min/max    : ', vmin, vmax
        write(logfhandle,'(A)')           '    Group   Class   Particles     Quota    Visits'
        do k = 1, nunits
            i = order(k)
            write(logfhandle,'(4X,I5,2X,I6,2X,I10,2X,F8.2,2X,F8.2)') &
                &clssmp(i)%group, clssmp(i)%clsind, clssmp(i)%pop, quotas(i), visits(i)
        end do
        write(logfhandle,'(A)') ''
        if( vmin < 1. )then
            msg = 'sampling coverage short: the least-visited unit is visited '//real2str(vmin)//&
                &' times per particle over the planned '//draw_label//'s; one sweep needs '//int2str(sweep)
            THROW_WARN(msg%to_char())
        endif
        if( vmax > VISITS_WARN_FAC * target )then
            msg = 'sampling coverage uneven: the most-visited unit is visited '//real2str(vmax)//&
                &' times per particle, above '//real2str(VISITS_WARN_FAC)//' times the target '//real2str(target)
            THROW_WARN(msg%to_char())
        endif
        call msg%kill
    end subroutine report_class_sample_coverage

end module simple_view_partition_sampling
