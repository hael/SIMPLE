module simple_cluster_cavgs
include 'simple_lib.f08'
use simple_strategy2D3D_common
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_aff_prop,         only: aff_prop
use simple_oris,             only: oris
use simple_builder,          only: build_glob
use simple_parameters,       only: params_glob
implicit none

public :: cluster_cavgs_exec
private

type(polarft_corrcalc) :: pftcc
real, allocatable      :: corrmat(:,:)
type(aff_prop)         :: aprop
integer, allocatable   :: centers(:), labels(:)

contains

    subroutine cluster_cavgs_exec
        character(len=:), allocatable :: fname
        real,             allocatable :: tmp(:), clspops(:), clsres(:), pops(:), mempops(:), memres(:)
        logical,          allocatable :: mask(:), included(:)
        integer,          allocatable :: inds(:), meminds(:), meminds_native(:)
        integer,          parameter   :: hlen=50
        logical    :: err
        real       :: simsum, corr_within_cls_avg, corr_within_cls_min, corr_within_cls_max, corr_med
        real       :: corr_between_cls_avg, corr_between_cls_min, corr_between_cls_max, pop, popmin
        real       :: popmax, popmed, popave, popsdev, popvar, szmax, scale
        integer    :: iptcl, jptcl, icen, jcen, cnt, npairs, i, j, nmems, curr_pop
        type(oris) :: ceninfo, clsdoc
        ! read class doc
        call clsdoc%new(params_glob%nptcls)
        call clsdoc%read(params_glob%classdoc)
        ! xtract class populations
        clspops = clsdoc%get_all('pop')
        ! xtract class resolutions
        clsres  = clsdoc%get_all('res')
        ! set bp range
        params_glob%kfromto(1) = max(2, calc_fourier_index(params_glob%hp, params_glob%boxmatch, params_glob%smpd))
        params_glob%kfromto(2) = calc_fourier_index(params_glob%lp, params_glob%boxmatch, params_glob%smpd)
        params_glob%kstop      = params_glob%kfromto(2)
        ! prep pftcc
        call preppftcc4cluster
        ! memoize FFTs for improved performance
        call pftcc%memoize_ffts
        ! calculate similarity matrix in parallel
        corrmat = pftcc%calc_roinv_corrmat()
        npairs = (params_glob%nptcls * (params_glob%nptcls - 1)) / 2
        allocate(tmp(npairs), mask(params_glob%nptcls), inds(params_glob%nptcls), included(params_glob%nptcls))
        cnt = 0
        do iptcl=1,params_glob%nptcls - 1
            do jptcl=iptcl + 1,params_glob%nptcls
                cnt = cnt + 1
                tmp(cnt) = corrmat(iptcl,jptcl)
            end do
        end do
        ! calculate similarity matrix in parallel of mirrored even references
        call pftcc%swap_ptclsevenodd
        corrmat = pftcc%calc_roinv_corrmat()
        cnt = 0
        do iptcl=1,params_glob%nptcls - 1
            do jptcl=iptcl + 1,params_glob%nptcls
                cnt = cnt + 1
                tmp(cnt) = max(corrmat(iptcl,jptcl), tmp(cnt))
            end do
        end do
        corr_med = median_nocopy(tmp)
        ! perform clustering with affinity propagation
        call aprop%new(params_glob%nptcls, corrmat, pref=corr_med/2.)
        call aprop%propagate(centers, labels, simsum)
        ! report clustering solution
        params_glob%ncls = size(centers)
        call clsdoc%new(params_glob%nptcls)
        do iptcl=1,params_glob%nptcls
            call clsdoc%set(iptcl, 'class', real(labels(iptcl)))
        end do
        call clsdoc%write('aff_prop_clustering'//trim(TXT_EXT), [1,params_glob%nptcls])
        ! calculate within cluster correlations
        corr_within_cls_avg =  0.
        corr_within_cls_min =  1.0
        corr_within_cls_max = -1.0
        cnt                 = 0
        do iptcl=1,params_glob%nptcls - 1
            do jptcl=iptcl + 1,params_glob%nptcls
                if( labels(iptcl) == labels(jptcl) )then
                    corr_within_cls_avg = corr_within_cls_avg + corrmat(iptcl,jptcl)
                    if( corrmat(iptcl,jptcl) < corr_within_cls_min ) corr_within_cls_min = corrmat(iptcl,jptcl)
                    if( corrmat(iptcl,jptcl) > corr_within_cls_max ) corr_within_cls_max = corrmat(iptcl,jptcl)
                    cnt = cnt + 1
                endif
            end do
        end do
        corr_within_cls_avg = corr_within_cls_avg / real(cnt)
        ! calculate between cluster correlations
        corr_between_cls_avg =  0.
        corr_between_cls_min =  1.0
        corr_between_cls_max = -1.0
        do icen=1,params_glob%ncls - 1
            do jcen=icen + 1, params_glob%ncls
                corr_between_cls_avg = corr_between_cls_avg + corrmat(centers(icen),centers(jcen))
                if( corrmat(centers(icen),centers(jcen)) < corr_between_cls_min )&
                &corr_between_cls_min = corrmat(centers(icen),centers(jcen))
                if( corrmat(centers(icen),centers(jcen)) > corr_between_cls_max )&
                &corr_between_cls_max = corrmat(centers(icen),centers(jcen))
            end do
        end do
        corr_between_cls_avg = corr_between_cls_avg / real((params_glob%ncls * (params_glob%ncls -1 )) / 2)
        write(*,'(a)') '>>> CLUSTERING STATISTICS'
        write(*,'(a,1x,f8.4)') '>>> # CLUSTERS FOUND        ', real(params_glob%ncls)
        write(*,'(a,1x,f8.4)') '>>> WITHIN  CLUSTER CORR    ', corr_within_cls_avg
        write(*,'(a,1x,f8.4)') '>>> BETWEEN CLUSTER CORR    ', corr_between_cls_avg
        call ceninfo%new(params_glob%ncls)
        do icen=1,params_glob%ncls
            where( labels == icen )
                mask = .true.
            else where
                mask = .false.
            end where
            pop = sum(clspops, mask)
            call ceninfo%set(icen, 'class',  real(icen))
            call ceninfo%set(icen, 'center', real(centers(icen)))
            call ceninfo%set(icen, 'pop',    pop)
            allocate(fname, source='cluster'//int2str_pad(icen,3)//'imgs'//params_glob%ext)
            cnt = 0
            do iptcl=1,params_glob%nptcls
                if( labels(iptcl) == icen )then
                    cnt = cnt + 1
                    call read_img( iptcl )
                    call build_glob%img%write(fname, cnt)
                endif
            enddo
            deallocate(fname)
        end do
        pops   = ceninfo%get_all('pop')
        popmin = minval(pops)
        popmax = maxval(pops)
        popmed = median(pops)
        call moment(pops, popave, popsdev, popvar, err)
        write(*,'(a,1x,f8.2)') '>>> MINIMUM POPULATION    ', popmin
        write(*,'(a,1x,f8.2)') '>>> MAXIMUM POPULATION    ', popmax
        write(*,'(a,1x,f8.2)') '>>> MEDIAN  POPULATION    ', popmed
        write(*,'(a,1x,f8.2)') '>>> AVERAGE POPULATION    ', popave
        write(*,'(a,1x,f8.2)') '>>> SDEV OF POPULATION    ', popsdev
        if( params_glob%balance > 0 )then
            inds     = (/(i,i=1,params_glob%nptcls)/)
            included = .true.
            do icen=1,params_glob%ncls
                if( pops(icen) > params_glob%balance )then
                    ! identify the members of the cluster
                    where( labels == icen )
                        mask = .true.
                    else where
                        mask = .false.
                    end where
                    ! # members
                    nmems = count(mask)
                    ! indices of members
                    allocate( meminds(nmems) )
                    meminds = (/(i,i=1,nmems)/)
                    meminds_native = pack(inds, mask=mask)
                    ! resolutions of members
                    memres  = pack(clsres, mask=mask)
                    ! populations of members
                    mempops = pack(clspops, mask=mask)
                    ! rank according to resolution
                    call hpsort(memres, meminds)
                    curr_pop = 0
                    do i=1,nmems
                        curr_pop = curr_pop + nint(mempops(meminds(i)))
                        if( curr_pop <= params_glob%balance )then
                            ! update pops
                            pops(icen) = curr_pop
                        else
                            ! reject this class
                            included(meminds_native(i)) = .false.
                            clspops(meminds_native(i))  = 0.

                        endif
                    end do
                    deallocate(meminds, memres, mempops)
                endif
            end do
            write(*,'(a)') '>>> CLUSTERING STATISTICS AFTER BALANCING'
            popmin = minval(pops)
               popmax = maxval(pops)
            popmed = median_nocopy(pops)
            call moment(pops, popave, popsdev, popvar, err)
            write(*,'(a,1x,f8.2)') '>>> MINIMUM POPULATION    ', popmin
            write(*,'(a,1x,f8.2)') '>>> MAXIMUM POPULATION    ', popmax
            write(*,'(a,1x,f8.2)') '>>> MEDIAN  POPULATION    ', popmed
            write(*,'(a,1x,f8.2)') '>>> AVERAGE POPULATION    ', popave
            write(*,'(a,1x,f8.2)') '>>> SDEV OF POPULATION    ', popsdev
        endif
        ! produce a histogram of class populations
        szmax = maxval(pops)
        ! scale to max 50 *:s
        scale = 1.0
        do while( nint(scale*szmax) > hlen )
            scale = scale - 0.001
        end do
        write(*,'(a)') '>>> HISTOGRAM OF CLUSTER POPULATIONS'
        call hpsort(pops)
        do icen=1,params_glob%ncls
            write(*,*) nint(pops(icen)),"|",('*', j=1,nint(pops(icen)*scale))
        end do
        call ceninfo%write('aff_prop_ceninfo'//trim(TXT_EXT), [1,params_glob%ncls])
        ! write selected class averages
        if( params_glob%balance > 0 )then
            cnt = 0
            do iptcl=1,params_glob%nptcls
                if( included(iptcl) )then
                    cnt = cnt + 1
                    call read_img( iptcl )
                    call build_glob%img%write('aff_prop_selected_imgs'//params_glob%ext, cnt)
                endif
            enddo
            write(*,'(i5,1x,a)') cnt, 'SELECTED CLASS AVERAGES WRITTEN TO FILE: aff_prop_selected_imgs.ext'
        endif
        call ceninfo%kill
        call aprop%kill
    end subroutine cluster_cavgs_exec

    !>  \brief  prepares the polarft corrcalc object for clustering
    subroutine preppftcc4cluster
        use simple_polarizer, only: polarizer
        type(polarizer), allocatable :: match_imgs(:), mirr_match_imgs(:)
        integer :: iptcl
        ! create the polarft_corrcalc object
        call pftcc%new(params_glob%nptcls, [1,params_glob%nptcls])
        ! prepare the polarizer images
        call build_glob%img_match%init_polarizer(pftcc, params_glob%alpha)
        allocate(match_imgs(params_glob%nptcls), mirr_match_imgs(params_glob%nptcls))
        do iptcl=1,params_glob%nptcls
            call match_imgs(iptcl)%new([params_glob%boxmatch, params_glob%boxmatch, 1], params_glob%smpd, wthreads=.false.)
            call match_imgs(iptcl)%copy_polarizer(build_glob%img_match)
            call mirr_match_imgs(iptcl)%new([params_glob%boxmatch, params_glob%boxmatch, 1], params_glob%smpd, wthreads=.false.)
            call mirr_match_imgs(iptcl)%copy_polarizer(build_glob%img_match)
        end do
        ! prepare image batch
        call prepimgbatch( params_glob%nptcls)
        call read_imgbatch( [1,params_glob%nptcls])
        ! PREPARATION OF PFTS IN PFTCC
        ! read references and transform into polar coordinates
        !$omp parallel do default(shared) private(iptcl)&
        !$omp schedule(static) proc_bind(close)
        do iptcl=1,params_glob%nptcls
            ! clip image if needed
            call build_glob%imgbatch(iptcl)%clip(match_imgs(iptcl))
            ! apply mask
            if( params_glob%l_innermsk )then
                call match_imgs(iptcl)%mask(params_glob%msk, 'soft', inner=params_glob%inner, width=params_glob%width)
            else
                call match_imgs(iptcl)%mask(params_glob%msk, 'soft')
            endif
            ! mirror
            call mirr_match_imgs(iptcl)%set_rmat( match_imgs(iptcl)%get_rmat() )
            call mirr_match_imgs(iptcl)%mirror('x')
            ! move to Fourier space
            call match_imgs(iptcl)%fft()
            call mirr_match_imgs(iptcl)%fft()
            ! transfer to polar coordinates in even
            call match_imgs(iptcl)%polarize(pftcc, iptcl, isptcl=.false., iseven=.true. )
            ! put mirror image in odd
            call mirr_match_imgs(iptcl)%polarize(pftcc, iptcl, isptcl=.false., iseven=.false. )
            ! copy reference to particle
            call pftcc%cp_even_ref2ptcl(iptcl, iptcl)
        end do
        !$omp end parallel do

        ! DESTRUCT
        do iptcl=1,params_glob%nptcls
            call mirr_match_imgs(iptcl)%kill_polarizer
            call mirr_match_imgs(iptcl)%kill
            call match_imgs(iptcl)%kill_polarizer
            call match_imgs(iptcl)%kill
            call build_glob%imgbatch(iptcl)%kill
        end do
        deallocate(match_imgs, mirr_match_imgs, build_glob%imgbatch)

    end subroutine preppftcc4cluster

end module simple_cluster_cavgs
