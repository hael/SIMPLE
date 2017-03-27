program simple_profile_corr
!$ use omp_lib
!$ use omp_lib_kinds
use simple_math,             only: csq
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_oris,             only: oris
implicit none
integer, parameter :: MODE=3
integer, parameter :: BOX=240, NK=50, RING2=100, NTHR=8, NREFS=10, NPTCLS=1000
real,    parameter :: SMPD=1.77
type(polarft_corrcalc) :: pftcc
type(oris)             :: a
integer                :: nrots, iptcl, iref, i
real                   :: corr
real, allocatable      :: corrmat3d(:,:,:)
!$ call omp_set_num_threads(NTHR)   
call a%new(NPTCLS)
do i=1,NPTCLS
    call a%set(i, 'kv',   300.)
    call a%set(i, 'cs',    2.7)
    call a%set(i, 'fraca', 0.1)
    call a%set(i, 'dfx',   2.0)
    call a%set(i, 'dfy',   2.0)
    call a%set(i, 'angast', 0.)
    call a%set(i, 'kv',   300.)
end do 
call pftcc%new(NREFS, [1,NPTCLS], [BOX,BOX,1], [2,NK], RING2, NTHR, 'yes')
call pftcc%create_polar_ctfmats(SMPD, a)
nrots = pftcc%get_nrots()
allocate(corrmat3d(NPTCLS,NREFS,nrots))


select case(MODE)
    ! 2.6 s -> 1.3 s
    case(1)
        
        do iptcl=1,NPTCLS
            !$omp parallel do schedule(auto) default(shared) private(iref)
            do iref=1,NREFS
                call pftcc%apply_ctf_single(iptcl, iref)
                corrmat3d(iptcl,iref,:) = pftcc%gencorrs_serial(iref,iptcl)
            end do
            !$omp end parallel do
        end do
        
    case(2)
    ! 5 s
        do iptcl=1,NPTCLS
            do iref=1,NREFS
                call pftcc%apply_ctf_single(iptcl, iref)
                corrmat3d(iptcl,iref,:) = pftcc%gencorrs(iref,iptcl)
            end do 
        end do
    case(3)
    ! 1.3 s
        call pftcc%gencorrs_all_cpu(corrmat3d)
end select

    




end program simple_profile_corr
