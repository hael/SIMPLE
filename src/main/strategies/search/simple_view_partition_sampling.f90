!@descr: view-cluster sampling groups for partition=yes: selected class averages clustered by aligned correlation, one group per cluster
! Within a group the particles are ordered by rank fraction inside their own class; ccs holds
! 1 - that fraction (an order key, not a correlation).
module simple_view_partition_sampling
use simple_pftc_srch_api
use simple_clustering_utils, only: cluster_dmat
use simple_strategy2D_utils, only: prep_cavgs4clust, calc_cluster_cavgs_dmat
implicit none

public :: make_view_partition_class_samples, VIEW_PARTITION_FILE
private
#include "simple_local_flags.inc"

character(len=*), parameter :: VIEW_PARTITION_FILE  = 'view_partition.txt'
character(len=*), parameter :: VIEW_PARTITION_FBODY = 'view_partition'
real,             parameter :: VIEW_PARTITION_LP    = 6.  !< low-pass limit of the class-average alignment (as cluster_cavgs)
real,             parameter :: VIEW_PARTITION_TRS   = 10. !< shift search range of the class-average alignment (as cluster_cavgs)

contains

    !> one class_sample per view cluster of the selected class averages; ntarget (particles sampled
    !! per iteration) only feeds the report; clustering runs on nthr threads
    subroutine make_view_partition_class_samples( params, spproj, ntarget, nthr, clssmp )
        class(parameters),               intent(in)    :: params
        class(sp_project),               intent(inout) :: spproj
        integer,                         intent(in)    :: ntarget
        integer,                         intent(in)    :: nthr
        type(class_sample), allocatable, intent(inout) :: clssmp(:)
        type(cmdline)            :: cline_clust
        type(parameters)         :: params_clust
        type(image), allocatable :: cavg_imgs(:)
        real,        allocatable :: mm(:,:), dmat(:,:)
        integer,     allocatable :: states(:), clspops(:), clsinds(:), labels(:), i_medoids(:)
        logical,     allocatable :: l_sel(:)
        integer(timer_int_kind) :: t_start, t_dmat
        real(timer_int_kind)    :: rt_dmat
        integer :: nsel, nclust, ldim(3), i
        logical :: l_clustered
        t_start = tic()
        rt_dmat = 0.
        if( spproj%os_cls2D%get_noris() < 1 ) THROW_HARD('partition=yes needs class averages (an empty cls2D segment)')
        states = spproj%os_cls2D%get_all_asint('state')
        l_sel  = states > 0
        nsel   = count(l_sel)
        if( nsel < 1 )          THROW_HARD('partition=yes needs selected class averages (cls2D state > 0)')
        if( params%nclust < 1 ) THROW_HARD('partition=yes needs nclust >= 1')
        select case(trim(params%clust_crit))
            case('cc', 'sig', 'res', 'hybrid')
            case DEFAULT
                THROW_HARD('partition=yes supports clust_crit=cc|sig|res|hybrid, not '//trim(params%clust_crit))
        end select
        l_clustered = nsel > params%nclust
        if( .not. l_clustered )then
            nclust  = nsel
            clsinds = pack([(i, i=1,size(states))], mask=l_sel)
            clspops = pack(spproj%os_cls2D%get_all_asint('pop'), mask=l_sel)
            labels  = [(i, i=1,nsel)]
        else
            nclust = params%nclust
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
            call prep_cavgs4clust(spproj, cavg_imgs, params_clust%mskdiam, clspops, clsinds, l_sel, mm)
            ldim              = cavg_imgs(1)%get_ldim()
            params_clust%smpd = cavg_imgs(1)%get_smpd()
            params_clust%box  = ldim(1)
            params_clust%msk  = min(real(params_clust%box/2) - COSMSKHALFWIDTH - 1., &
                &0.5 * params_clust%mskdiam / params_clust%smpd)
            write(logfhandle,'(A,I0,A,I0,A,A,A,I0,A)') '>>> VIEW PARTITION: clustering ', nsel, ' class averages into ', &
                &nclust, ' groups by average linkage on clust_crit=', trim(params%clust_crit), ' (', params_clust%nthr, ' threads)'
            t_dmat  = tic()
            dmat    = calc_cluster_cavgs_dmat(params_clust, cavg_imgs, [minval(mm(:,1)), maxval(mm(:,2))], trim(params%clust_crit))
            rt_dmat = toc(t_dmat)
            call cluster_dmat(dmat, 'avglink', nclust, i_medoids, labels)
            call dealloc_imgarr(cavg_imgs)
            call cline_clust%kill
            !$ call omp_set_num_threads(params%nthr)
            nthr_glob = params%nthr
        endif
        ! inspection output
        cavg_imgs = read_cavgs_into_imgarr(spproj, mask=l_sel)
        call write_imgarr(nsel, cavg_imgs, labels, VIEW_PARTITION_FBODY, params%ext%to_char())
        call dealloc_imgarr(cavg_imgs)
        call write_partition_table
        call make_groups
        call print_report

    contains

        !> one class_sample per group, particles ordered by rank fraction within their class
        subroutine make_groups
            integer, allocatable :: members(:), pinds_cls(:), pinds_grp(:)
            real,    allocatable :: corrs_cls(:), keys(:)
            integer :: ic, im, j, npop
            if( allocated(clssmp) ) call deallocate_class_samples(clssmp)
            allocate(clssmp(nclust))
            do ic = 1, nclust
                members = pack(clsinds, mask=labels == ic)
                allocate(pinds_grp(0), keys(0))
                do im = 1, size(members)
                    call spproj%os_ptcl2D%get_pinds(members(im), 'class', pinds_cls)
                    if( .not. allocated(pinds_cls) ) cycle
                    npop = size(pinds_cls)
                    if( npop > 0 )then
                        allocate(corrs_cls(npop))
                        do j = 1, npop
                            corrs_cls(j) = spproj%os_ptcl2D%get(pinds_cls(j), 'corr')
                        enddo
                        call hpsort(corrs_cls, pinds_cls)
                        call reverse(pinds_cls)
                        pinds_grp = [pinds_grp, pinds_cls]
                        keys      = [keys, [((real(j) - 0.5) / real(npop), j=1,npop)]]
                        deallocate(corrs_cls)
                    endif
                    deallocate(pinds_cls)
                enddo
                if( size(pinds_grp) > 1 ) call hpsort(keys, pinds_grp)
                clssmp(ic)%clsind = ic
                clssmp(ic)%pop    = size(pinds_grp)
                allocate(clssmp(ic)%pinds(clssmp(ic)%pop), source=pinds_grp)
                allocate(clssmp(ic)%ccs(clssmp(ic)%pop),   source=1. - keys)
                deallocate(pinds_grp, keys, members)
            enddo
        end subroutine make_groups

        !> per-group populations and per-iteration samples (water-filling, as sample4update_class)
        subroutine print_report
            integer :: pops(nclust), nsmp(nclust), ncls_grp(nclust), ic, ntot, nsmp_tot
            real    :: pct_ptcls(nclust), pct_smp(nclust)
            pops = clssmp(:)%pop
            do ic = 1, nclust
                ncls_grp(ic) = count(labels == ic)
            enddo
            nsmp = 0
            do while( sum(nsmp) < ntarget .and. any(nsmp < pops) )
                where( nsmp < pops ) nsmp = nsmp + 1
            enddo
            ntot      = sum(pops)
            nsmp_tot  = sum(nsmp)
            pct_ptcls = 100. * real(pops) / real(max(1, ntot))
            pct_smp   = 100. * real(nsmp) / real(max(1, nsmp_tot))
            write(logfhandle,'(A)') ''
            write(logfhandle,'(A)')               '>>> VIEW PARTITION SAMPLING GROUPS'
            if( l_clustered )then
                write(logfhandle,'(A,A)')         '    Grouping                      :   average linkage on clust_crit=', &
                    &trim(params%clust_crit)
                write(logfhandle,'(A,I12)')       '    Clustering threads            : ', max(1, nthr)
                write(logfhandle,'(A,F12.1)')     '    Distance matrix time (s)      : ', rt_dmat
            else
                write(logfhandle,'(A,I0,A)')      '    Grouping                      :   one group per class '//&
                    &'(selected classes <= nclust=', params%nclust, ')'
            endif
            write(logfhandle,'(A,F12.1)')         '    Total time (s)                : ', toc(t_start)
            write(logfhandle,'(A,I12)')           '    Selected classes              : ', nsel
            write(logfhandle,'(A,I12)')           '    Groups                        : ', nclust
            write(logfhandle,'(A,I12)')           '    Active particles              : ', ntot
            write(logfhandle,'(A,I12,F9.2,A)')    '    Sampled per iteration         : ', nsmp_tot, &
                &100. * real(nsmp_tot) / real(max(1, ntot)), '%'
            if( count(pops > 0) > 0 )then
                write(logfhandle,'(A,F12.2)')     '    Group share max/min, particles: ', &
                    &maxval(pct_ptcls) / minval(pct_ptcls, mask=pops > 0)
                write(logfhandle,'(A,F12.2)')     '    Group share max/min, sampled  : ', &
                    &maxval(pct_smp) / max(TINY, minval(pct_smp, mask=pops > 0))
            endif
            write(logfhandle,'(A)')               '    Group  Classes     Particles  Pct ptcls       Sampled  Pct sample'
            do ic = 1, nclust
                write(logfhandle,'(4X,I5,2X,I7,2X,I12,2X,F9.2,2X,I12,2X,F10.2)') &
                    &ic, ncls_grp(ic), pops(ic), pct_ptcls(ic), nsmp(ic), pct_smp(ic)
            enddo
            write(logfhandle,'(A)') ''
        end subroutine print_report

        subroutine write_partition_table
            integer :: funit, io_stat, k
            call fopen(funit, file=string(VIEW_PARTITION_FILE), status='REPLACE', action='WRITE', iostat=io_stat)
            call fileiochk('make_view_partition_class_samples; '//VIEW_PARTITION_FILE, io_stat)
            write(funit,'(A)') '# class group class_population'
            do k = 1, nsel
                write(funit,'(I6,1X,I4,1X,I8)') clsinds(k), labels(k), clspops(k)
            enddo
            call fclose(funit)
        end subroutine write_partition_table

    end subroutine make_view_partition_class_samples

end module simple_view_partition_sampling
