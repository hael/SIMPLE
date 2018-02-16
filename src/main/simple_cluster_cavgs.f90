module simple_cluster_cavgs
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_build,            only: build
use simple_params,           only: params
use simple_corrmat,          only: calc_roinv_corrmat
use simple_aff_prop,         only: aff_prop
use simple_hadamard_common   ! use all in there
use simple_defs              ! use all in there
use simple_binoris_io        ! use all in there
implicit none

public :: cluster_cavgs_exec
private

type(polarft_corrcalc) :: pftcc
real, allocatable      :: corrmat(:,:)
type(aff_prop)         :: aprop
integer, allocatable   :: centers(:), labels(:)

contains

	subroutine cluster_cavgs_exec( b, p )
		use simple_oris,    only: oris
		use simple_math,    only: calc_fourier_index, median_nocopy, hpsort, median
		use simple_strings, only: int2str_pad
		use simple_stat,    only: moment
		class(build),  intent(inout) :: b
        class(params), intent(inout) :: p
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
        call clsdoc%new_clean(p%nptcls)
        call clsdoc%read(p%classdoc)
        ! xtract class populations
        clspops = clsdoc%get_all('pop')
        ! xtract class resolutions
        clsres  = clsdoc%get_all('res')
        ! set bp range
        p%kfromto(1) = max(2, calc_fourier_index(p%hp, p%boxmatch, p%smpd))
    	p%kfromto(2) = calc_fourier_index(p%lp, p%boxmatch, p%smpd)
    	p%lp_dyn     = p%lp
        ! prep pftcc
        call preppftcc4cluster( b, p )
        ! memoize FFTs for improved performance
        call pftcc%memoize_ffts
		! calculate similarity matrix in parallel
		call calc_roinv_corrmat( pftcc, corrmat )
		npairs = (p%nptcls * (p%nptcls - 1)) / 2
		allocate(tmp(npairs), mask(p%nptcls), inds(p%nptcls), included(p%nptcls))
		cnt = 0
		do iptcl=1,p%nptcls - 1
			do jptcl=iptcl + 1,p%nptcls
				cnt = cnt + 1
				tmp(cnt) = corrmat(iptcl,jptcl)
			end do
		end do
		corr_med = median_nocopy(tmp)
		! perform clustering with affinity propagation
		call aprop%new(p%nptcls, corrmat, pref=corr_med)
		call aprop%propagate(centers, labels, simsum)
		! report clustering solution
		p%ncls = size(centers)
		call clsdoc%new_clean(p%nptcls)
		do iptcl=1,p%nptcls
			call clsdoc%set(iptcl, 'class', real(labels(iptcl)))
		end do
		call binwrite_oritab('aff_prop_clustering'//trim(METADATA_EXT), clsdoc, [1,p%nptcls])
		! calculate within cluster correlations
		corr_within_cls_avg =  0.
		corr_within_cls_min =  1.0
		corr_within_cls_max = -1.0
		cnt                 = 0
		do iptcl=1,p%nptcls - 1
			do jptcl=iptcl + 1,p%nptcls
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
		do icen=1,p%ncls - 1
			do jcen=icen + 1, p%ncls
				corr_between_cls_avg = corr_between_cls_avg + corrmat(centers(icen),centers(jcen))
				if( corrmat(centers(icen),centers(jcen)) < corr_between_cls_min )&
				&corr_between_cls_min = corrmat(centers(icen),centers(jcen))
				if( corrmat(centers(icen),centers(jcen)) > corr_between_cls_max )&
				&corr_between_cls_max = corrmat(centers(icen),centers(jcen))
			end do
		end do
		corr_between_cls_avg = corr_between_cls_avg / real((p%ncls * (p%ncls -1 )) / 2)
		write(*,'(a)') '>>> CLUSTERING STATISTICS'
		write(*,'(a,1x,f8.4)') '>>> # CLUSTERS FOUND        ', real(p%ncls)
		write(*,'(a,1x,f8.4)') '>>> WITHIN  CLUSTER CORR    ', corr_within_cls_avg
		write(*,'(a,1x,f8.4)') '>>> BETWEEN CLUSTER CORR    ', corr_between_cls_avg
		call ceninfo%new_clean(p%ncls)
		do icen=1,p%ncls
			where( labels == icen )
				mask = .true.
			else where
				mask = .false.
			end where
			pop = sum(clspops, mask)
			call ceninfo%set(icen, 'class',  real(icen))
			call ceninfo%set(icen, 'center', real(centers(icen)))
			call ceninfo%set(icen, 'pop',    pop)
			allocate(fname, source='cluster'//int2str_pad(icen,3)//'imgs'//p%ext)
			cnt = 0
			do iptcl=1,p%nptcls
				if( labels(iptcl) == icen )then
					cnt = cnt + 1
					call read_img( b, p, iptcl )
					call b%img%write(fname, cnt)
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
        if( p%balance > 0 )then
        	inds     = (/(i,i=1,p%nptcls)/)
        	included = .true.
        	do icen=1,p%ncls
        		if( pops(icen) > p%balance )then
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
						if( curr_pop <= p%balance )then
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
        do icen=1,p%ncls
            write(*,*) nint(pops(icen)),"|",('*', j=1,nint(pops(icen)*scale))
        end do
		call binwrite_oritab('aff_prop_ceninfo'//trim(METADATA_EXT), ceninfo, [1,p%ncls])
		! write selected class averages
		if( p%balance > 0 )then
			cnt = 0
			do iptcl=1,p%nptcls
				if( included(iptcl) )then
					cnt = cnt + 1
					call read_img( b, p, iptcl )
					call b%img%write('aff_prop_selected_imgs'//p%ext, cnt)
				endif
			enddo
			write(*,'(i5,1x,a)') cnt, 'SELECTED CLASS AVERAGES WRITTEN TO FILE: aff_prop_selected_imgs.ext'
		endif
		call ceninfo%kill
		call aprop%kill
	end subroutine cluster_cavgs_exec

	!>  \brief  prepares the polarft corrcalc object for clustering
    subroutine preppftcc4cluster( b, p )
    	use simple_polarizer, only: polarizer
        class(build),  intent(inout) :: b
        class(params), intent(inout) :: p
    	type(polarizer), allocatable :: match_imgs(:)
    	integer :: iptcl
    	logical :: do_center
		! create the polarft_corrcalc object
		call pftcc%new(p%nptcls, p)
		! prepare the polarizer images
		call b%img_match%init_polarizer(pftcc, p%alpha)
		allocate(match_imgs(p%nptcls))
        do iptcl=1,p%nptcls
            call match_imgs(iptcl)%new([p%boxmatch, p%boxmatch, 1], p%smpd)
            call match_imgs(iptcl)%copy_polarizer(b%img_match)
        end do
        ! prepare image batch
        call prepimgbatch(b, p, p%nptcls)
        call read_imgbatch( b, p, [1,p%nptcls])

        ! PREPARATION OF PFTS IN PFTCC
        ! read references and transform into polar coordinates
        !$omp parallel do default(shared) private(iptcl)&
        !$omp schedule(static) proc_bind(close)
        do iptcl=1,p%nptcls
            ! clip image if needed
	        call b%imgbatch(iptcl)%clip(match_imgs(iptcl))
	        ! apply mask
	        if( p%l_innermsk )then
	            call match_imgs(iptcl)%mask(p%msk, 'soft', inner=p%inner, width=p%width)
	        else
	            call match_imgs(iptcl)%mask(p%msk, 'soft')
	        endif
	        ! move to Fourier space
	        call match_imgs(iptcl)%fwd_ft
            ! transfer to polar coordinates
            call match_imgs(iptcl)%polarize(pftcc, iptcl, isptcl=.false., iseven=.true. )
            ! put in both even and odd positions
	        call pftcc%cp_even2odd_ref(iptcl)
	        ! copy reference to particle
	        call pftcc%cp_even_ref2ptcl(iptcl, iptcl)
        end do
        !$omp end parallel do

		! DESTRUCT
        do iptcl=1,p%nptcls
            call match_imgs(iptcl)%kill_polarizer
            call match_imgs(iptcl)%kill
            call b%imgbatch(iptcl)%kill
        end do
        deallocate(match_imgs, b%imgbatch)

	end subroutine preppftcc4cluster

end module simple_cluster_cavgs
