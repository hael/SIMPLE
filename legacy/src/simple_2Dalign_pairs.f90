module simple_2Dalign_pairs
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_hadamard_common   ! use all in there
use simple_pftcc_shsrch      ! use all in there
use simple_magic_boxes       ! use all in there
use simple_build,            only: build
use simple_params,           only: params
use simple_cmdline,          only: cmdline
implicit none

public :: init_2Dalign_pairs, simmat_2Dalign_pairs
private

type(polarft_corrcalc) :: pftcc
real :: lims(2,2)

contains

    subroutine init_2Dalign_pairs( p, b, cline )
        class(build),   intent(inout) :: b
        class(params),  intent(inout) :: p
        class(cmdline), intent(inout) :: cline   
        integer :: iptcl, box_new
        real    :: msk_new, smpd_new
        ! read images into array b%imgs
        call read_imgs_from_stk( b, p )
        ! scale
        call autoscale(p%box, p%msk, p%smpd, box_new, msk_new, smpd_new )
        p%box  = box_new
        p%msk  = msk_new
        p%smpd = smpd_new
        do iptcl=1,p%nptcls
            call b%imgs(iptcl)%fwd_ft
            call b%imgs(iptcl)%clip_inplace([p%box,p%box,1])
            call b%imgs(iptcl)%bwd_ft
        end do
        ! set Fourier index range
        call set_bp_range(b, p, cline)
        ! prepare corr calculator
        call pftcc%new(p%nptcls, [1,p%nptcls], [p%box,p%box,1], p%kfromto, p%ring2, 'no')
        ! prepare the polarizer
        call b%img%init_imgpolarizer(pftcc)
        do iptcl=1,p%nptcls
            b%img = b%imgs(iptcl) ! put the original image back
            ! prepare as we normally prepare the references
            ! assumes that b%a is setup accordingly (class=1 .. nptcls for entries 1 .. nptcls)
            call prep2Dref(p, b%img, b%a, iptcl)
            ! but transfer to polar coordinates as particle (full rotation range)
            call b%img%imgpolarizer(pftcc, iptcl)
        end do
        ! set limits for shift search
        lims(:,1) = -p%trs
        lims(:,2) =  p%trs 
    end subroutine init_2Dalign_pairs

    subroutine simmat_2Dalign_pairs( p, simmat )
        use simple_jiffys, only: alloc_err
        class(params),     intent(inout) :: p
        real, allocatable, intent(out)   :: simmat(:,:)
        real    :: cc(pftcc%get_nrots())
        integer :: alloc_stat, iptcl, jptcl
        if( allocated(simmat) ) deallocate(simmat)
        allocate(simmat(p%nptcls,p%nptcls), stat=alloc_stat)
        call alloc_err("In: simple_2Dalign_pairs :: simmat_2Dalign_pairs", alloc_stat)
        forall(iptcl=1:p%nptcls) simmat(iptcl,iptcl) = 1.0
        !$omp parallel default(shared) private(jptcl)
        do iptcl=1,p%nptcls-1
            !$omp single
            call pftcc%cp_ptcl2ref(iptcl, 1)
            !$omp end single nowait
            !$omp do schedule(auto)
            do jptcl=iptcl+1,p%nptcls
                cc = pftcc%gencorrs_serial(1, jptcl)
                simmat(iptcl,jptcl) = maxval(cc)
                simmat(jptcl,iptcl) = simmat(iptcl,jptcl) 
            end do
            !$omp end do
        end do
        !$omp end parallel  
    end subroutine simmat_2Dalign_pairs

    subroutine kill_2Dalign_pairs
        call pftcc%kill
    end subroutine kill_2Dalign_pairs

end module simple_2Dalign_pairs
