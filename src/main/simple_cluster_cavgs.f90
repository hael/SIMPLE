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
		use simple_oris, only: oris
		use simple_math, only: calc_fourier_index
		class(build),  intent(inout) :: b
        class(params), intent(inout) :: p
        real       :: simsum
        integer    :: iptcl, icen
        type(oris) :: ceninfo
        ! set bp range
        p%kfromto(1) = max(2, calc_fourier_index(p%hp, p%boxmatch, p%smpd))
    	p%kfromto(2) = calc_fourier_index(p%lp, p%boxmatch, p%smpd)
    	p%lp_dyn     = p%lp
    	call b%a%set_all2single('lp',p%lp)
        ! prep pftcc
        call preppftcc4cluster( b, p )
        ! memoize FFTs for improved performance
        call pftcc%memoize_ffts
		! calculate similarity matrix in parallel
		call calc_roinv_corrmat( pftcc, p%trs, corrmat )
		! perform clustering with affinity propagation
		call aprop%new(p%nptcls, corrmat)
		call aprop%propagate(centers, labels, simsum)
		! report clustering solution
		p%ncls = size(centers)
		write(*,'(a,1x,i9)') '# CLUSTERS FOUND:', p%ncls
		do iptcl=1,p%nptcls
			call b%a%set(iptcl, 'class', real(labels(iptcl)))
		end do
		call binwrite_oritab('aff_prop_clustering'//METADATEXT, b%a, [1,p%nptcls])
		call ceninfo%new_clean(p%ncls)
		do icen=1,p%ncls 
			call ceninfo%set(icen, 'class',  real(icen))
			call ceninfo%set(icen, 'center', real(centers(icen)))
		end do
		call binwrite_oritab('aff_prop_ceninfo'//METADATEXT, ceninfo, [1,p%ncls])
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
