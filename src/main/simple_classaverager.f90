module simple_classaverager
include 'simple_lib.f08'
!$ use omp_lib
use simple_builder,    only: build_glob
use simple_parameters, only: params_glob
use simple_ctf,        only: ctf
use simple_image,      only: image, image_ptr
use simple_stack_io,   only: stack_io
implicit none

public :: cavger_new, cavger_transf_oridat, cavger_gen2Dclassdoc, cavger_assemble_sums,&
cavger_merge_eos_and_norm, cavger_calc_and_write_frcs_and_eoavg, cavger_write, cavger_read,&
cavger_readwrite_partial_sums, cavger_assemble_sums_from_parts, cavger_kill, cavgs_even, cavgs_odd, cavgs_merged
private
#include "simple_local_flags.inc"

type ptcl_record
    type(ctf) :: tfun          !< transfer function
    real      :: pw      = 0.0 !< particle weight
    real      :: dfx     = 0.0 !< defocus in x (microns)
    real      :: dfy     = 0.0 !< defocus in y (microns)
    real      :: angast  = 0.0 !< angle of astigmatism (in degrees)
    real      :: phshift = 0.0 !< additional phase shift from the Volta
    real      :: e3            !< in-plane rotations
    real      :: shift(2)      !< rotational origin shift
    integer   :: pind    = 0   !< particle index in stack
    integer   :: eo      = -1  !< even is 0, odd is 1, default is -1
    integer   :: class         !< class assignment
end type ptcl_record

integer                        :: ctfflag                  !< ctf flag <yes=1|no=0|flip=2>
integer                        :: istart     = 0, iend = 0 !< particle index range
integer                        :: partsz     = 0           !< size of partition
integer                        :: ncls       = 0           !< # classes
integer                        :: filtsz     = 0           !< size of filter function or FSC
integer                        :: ldim(3)    = [0,0,0]     !< logical dimension of image
integer                        :: ldim_pd(3) = [0,0,0]     !< logical dimension of image, padded
real                           :: smpd       = 0.          !< sampling distance
type(ptcl_record), allocatable :: precs(:)                 !< particle records
type(image),       allocatable :: cavgs_even(:)            !< class averages
type(image),       allocatable :: cavgs_odd(:)             !< -"-
type(image),       allocatable :: cavgs_merged(:)          !< -"-
type(image),       allocatable :: ctfsqsums_even(:)        !< CTF**2 sums for Wiener normalisation
type(image),       allocatable :: ctfsqsums_odd(:)         !< -"-
type(image),       allocatable :: ctfsqsums_merged(:)      !< -"-
integer,           allocatable :: prev_eo_pops(:,:)
logical,           allocatable :: pptcl_mask(:)
logical                        :: phaseplate = .false.     !< Volta phaseplate images or not
logical                        :: l_bilinear = .true.      !< whether to use bilinear or convolution interpolation
logical                        :: exists     = .false.     !< to flag instance existence

integer(timer_int_kind) :: t_class_loop,t_batch_loop, t_gridding, t_init, t_tot
real(timer_int_kind)    :: rt_class_loop,rt_batch_loop, rt_gridding, rt_init, rt_tot
character(len=STDLEN)   :: benchfname

contains

    !>  \brief  is a constructor
    subroutine cavger_new( ptcl_mask )
        logical, optional, intent(in) :: ptcl_mask(params_glob%fromp:params_glob%top)
        integer :: icls
        ! destruct possibly pre-existing instance
        call cavger_kill
        if( present(ptcl_mask) )then
            allocate(pptcl_mask(params_glob%fromp:params_glob%top), source=ptcl_mask)
        else
            allocate(pptcl_mask(params_glob%fromp:params_glob%top), source=.true.)
        endif
        ncls       = params_glob%ncls
        ! work out range and partsz
        if( params_glob%l_distr_exec )then
            istart = params_glob%fromp
            iend   = params_glob%top
        else
            istart = 1
            iend   = params_glob%nptcls
        endif
        partsz     = size(pptcl_mask)
        ! CTF logics
        ctfflag    = build_glob%spproj%get_ctfflag_type('ptcl2D',iptcl=params_glob%fromp)
        ! set phaseplate flag
        phaseplate = build_glob%spproj%has_phaseplate('ptcl2D')
        ! smpd
        smpd       = params_glob%smpd
        ! set ldims
        ldim       = [params_glob%box,params_glob%box,1]
        ldim_pd    = [params_glob%boxpd,params_glob%boxpd,1]
        ldim_pd(3) = 1
        filtsz     = build_glob%img%get_filtsz()
        ! build arrays
        allocate(precs(partsz), cavgs_even(ncls), cavgs_odd(ncls),&
        &cavgs_merged(ncls), ctfsqsums_even(ncls),&
        &ctfsqsums_odd(ncls), ctfsqsums_merged(ncls), prev_eo_pops(ncls,2))
        prev_eo_pops = 0
        !$omp parallel do default(shared) private(icls) schedule(static) proc_bind(close)
        do icls=1,ncls
            call cavgs_even(icls)%new(ldim_pd,params_glob%smpd,wthreads=.false.)
            call cavgs_odd(icls)%new(ldim_pd,params_glob%smpd,wthreads=.false.)
            call cavgs_merged(icls)%new(ldim_pd,params_glob%smpd,wthreads=.false.)
            call ctfsqsums_even(icls)%new(ldim_pd,params_glob%smpd,wthreads=.false.)
            call ctfsqsums_odd(icls)%new(ldim_pd,params_glob%smpd,wthreads=.false.)
            call ctfsqsums_merged(icls)%new(ldim_pd,params_glob%smpd,wthreads=.false.)
        end do
        !$omp end parallel do
        ! flag existence
        exists = .true.
    end subroutine cavger_new

    ! setters/getters

    !>  \brief  transfers orientation data to the instance
    subroutine cavger_transf_oridat( spproj )
        use simple_sp_project, only: sp_project
        class(sp_project), intent(inout) :: spproj
        type(ctfparams)   :: ctfvars(nthr_glob)
        integer           :: i, icls, cnt, iptcl, ithr
        ! build index map
        cnt = 0
        do iptcl=istart,iend
            cnt = cnt + 1
            ! exclusion
            precs(cnt)%pind = 0
            if( spproj%os_ptcl2D%get_state(iptcl) == 0 ) cycle
            if( spproj%os_ptcl2D%get(iptcl,'w') < TINY ) cycle
            if( .not.pptcl_mask(iptcl)                 ) cycle
            precs(cnt)%pind = iptcl
        enddo
        ! fetch data from project
        !$omp parallel do default(shared) private(cnt,iptcl,ithr) schedule(static) proc_bind(close)
        do cnt = 1,partsz
            iptcl              = precs(cnt)%pind
            if( iptcl == 0 ) cycle
            ithr               = omp_get_thread_num() + 1
            precs(cnt)%eo      = nint(spproj%os_ptcl2D%get(iptcl,'eo'))
            precs(cnt)%pw      = spproj%os_ptcl2D%get(iptcl,'w')
            ctfvars(ithr)      = spproj%get_ctfparams('ptcl2D',iptcl)
            precs(cnt)%tfun    = ctf(params_glob%smpd, ctfvars(ithr)%kv, ctfvars(ithr)%cs, ctfvars(ithr)%fraca)
            precs(cnt)%dfx     = ctfvars(ithr)%dfx
            precs(cnt)%dfy     = ctfvars(ithr)%dfy
            precs(cnt)%angast  = ctfvars(ithr)%angast
            precs(cnt)%phshift = 0.
            if( phaseplate ) precs(cnt)%phshift = ctfvars(ithr)%phshift
            precs(cnt)%class   = nint(spproj%os_ptcl2D%get(iptcl, 'class'))
            precs(cnt)%e3      = spproj%os_ptcl2D%e3get(iptcl)
            precs(cnt)%shift   = spproj%os_ptcl2D%get_2Dshift(iptcl)
        end do
        !$omp end parallel do
        prev_eo_pops = 0
        if( trim(params_glob%stream).eq.'yes' .and. spproj%os_cls2D%get_noris() == ncls )then
            do i = 1,ncls
                icls = nint(spproj%os_cls2D%get(i,'class'))
                if( .not.spproj%os_cls2D%isthere(i,'prev_pop_even') ) cycle
                prev_eo_pops(icls,1) = nint(spproj%os_cls2D%get(i,'prev_pop_even'))
                prev_eo_pops(icls,2) = nint(spproj%os_cls2D%get(i,'prev_pop_odd'))
            enddo
        endif
    end subroutine cavger_transf_oridat

    !>  \brief prepares a 2D class document with class index, resolution,
    !!         poulation, average correlation and weight
    subroutine cavger_gen2Dclassdoc( spproj )
        use simple_sp_project, only: sp_project
        class(sp_project), intent(inout) :: spproj
        integer  :: pops(params_glob%ncls)
        real(dp) :: corrs(params_glob%ncls), ws(params_glob%ncls), specscores(params_glob%ncls)
        real     :: frc05, frc0143, rstate
        integer  :: i, iptcl, icls, pop, nptcls
        nptcls     = spproj%os_ptcl2D%get_noris()
        pops       = 0
        corrs      = 0.d0
        ws         = 0.d0
        specscores = 0.d0
        !$omp parallel do default(shared) private(iptcl,rstate,icls) schedule(static)&
        !$omp proc_bind(close) reduction(+:pops,corrs,ws,specscores)
        do iptcl=1,nptcls
            rstate = spproj%os_ptcl2D%get(iptcl,'state')
            if( rstate < 0.5 )cycle
            icls = nint(spproj%os_ptcl2D%get(iptcl,'class'))
            if( icls<1 .or. icls>params_glob%ncls )cycle
            pops(icls)       = pops(icls)      + 1
            corrs(icls)      = corrs(icls)     + spproj%os_ptcl2D%get(iptcl,'corr')
            ws(icls)         = ws(icls)        + spproj%os_ptcl2D%get(iptcl,'w')
            specscores(icls) = specscores(icls)+ spproj%os_ptcl2D%get(iptcl,'specscore')
        enddo
        !$omp end parallel do
        if( trim(params_glob%stream).eq.'yes'  .and.&
            &spproj%os_cls2D%get_noris()==ncls .and. params_glob%update_frac<.99 )then
            do i = 1,ncls
                icls = nint(spproj%os_cls2D%get(i,'class'))
                if( .not.spproj%os_cls2D%isthere(i,'prev_pop_even') ) cycle
                prev_eo_pops(icls,1) = nint(spproj%os_cls2D%get(i,'prev_pop_even'))
                prev_eo_pops(icls,2) = nint(spproj%os_cls2D%get(i,'prev_pop_odd'))
                pop = sum(prev_eo_pops(icls,:))
                if( pop == 0 ) cycle
                corrs(icls)      = corrs(icls)      + real(pop) * spproj%os_cls2D%get(i,'corr')
                ws(icls)         = ws(icls)         + real(pop) * spproj%os_cls2D%get(i,'w')
                specscores(icls) = specscores(icls) + real(pop) * spproj%os_cls2D%get(i,'specscore')
                pops(icls)       = pops(icls) + pop
            enddo
        endif
        where(pops>1)
            corrs      = corrs / real(pops)
            ws         = ws / real(pops)
            specscores = specscores / real(pops)
        elsewhere
            corrs      = -1.
            ws         = 0.
            specscores = 0.
        end where
        call spproj%os_cls2D%new(params_glob%ncls, is_ptcl=.false.)
        do icls=1,params_glob%ncls
            pop = pops(icls)
            call build_glob%clsfrcs%estimate_res(icls, frc05, frc0143)
            call spproj%os_cls2D%set(icls, 'class',     real(icls))
            call spproj%os_cls2D%set(icls, 'pop',       real(pop))
            call spproj%os_cls2D%set(icls, 'res',       frc0143)
            call spproj%os_cls2D%set(icls, 'corr',      real(corrs(icls)))
            call spproj%os_cls2D%set(icls, 'w',         real(ws(icls)))
            call spproj%os_cls2D%set(icls, 'specscore', real(specscores(icls)))
            if( pop > 0 )then
                call spproj%os_cls2D%set(icls, 'state', 1.0) ! needs to be default val if no selection has been done
            else
                call spproj%os_cls2D%set(icls, 'state', 0.0) ! exclusion
            endif
        end do
    end subroutine cavger_gen2Dclassdoc

    !>  \brief  is for initialization of the sums
    subroutine init_cavgs_sums
        integer :: icls
        !$omp parallel do default(shared) private(icls) schedule(static) proc_bind(close)
        do icls=1,ncls
            call cavgs_even(icls)%new(ldim_pd,smpd,wthreads=.false.)
            call cavgs_odd(icls)%new(ldim_pd,smpd,wthreads=.false.)
            call cavgs_merged(icls)%new(ldim_pd,smpd,wthreads=.false.)
            call cavgs_even(icls)%zero_and_flag_ft
            call cavgs_odd(icls)%zero_and_flag_ft
            call cavgs_merged(icls)%zero_and_flag_ft
            call ctfsqsums_even(icls)%zero_and_flag_ft
            call ctfsqsums_odd(icls)%zero_and_flag_ft
            call ctfsqsums_merged(icls)%zero_and_flag_ft
        end do
        !$omp end parallel do
    end subroutine init_cavgs_sums

    !>  \brief  is for calculating class population
    integer function class_pop( class )
        integer, intent(in) :: class
        class_pop = sum(eo_class_pop(class))
    end function class_pop

    !>  \brief  is for calculating even/odd class population
    function eo_class_pop( class ) result( pops )
        integer, intent(in) :: class
        integer :: pops(2), iprec
        pops = 0
        do iprec=1,partsz
            if( precs(iprec)%pind > 0 .and. precs(iprec)%class .eq. class )then
                if( precs(iprec)%eo == 1 )then
                    pops(2) = pops(2) + 1
                else
                    pops(1) = pops(1) + 1
                endif
            endif
        end do
    end function eo_class_pop

    ! calculators

    !>  \brief  is for assembling the sums in distributed/non-distributed mode
    !!          using gridding interpolation in Fourier space
    subroutine cavger_assemble_sums( do_frac_update )
        use simple_kbinterpol, only: kbinterpol
        logical, intent(in)           :: do_frac_update
        integer, parameter            :: READBUFFSZ = 1024
        complex, parameter            :: zero = cmplx(0.,0.)
        type(kbinterpol)              :: kbwin
        type(stack_io)                :: stkio_r
        type(image_ptr)               :: pcmat(params_glob%nthr), prhomat(params_glob%nthr)
        type(image),      allocatable :: cgrid_imgs(:), read_imgs(:)
        character(len=:), allocatable :: stk_fname
        complex,          allocatable :: cmats(:,:,:)
        real,             allocatable :: rhos(:,:,:), w(:,:)
        integer,          allocatable :: cyc1(:), cyc2(:)
        complex :: fcomp, cswap
        real    :: loc(2), mat(2,2), dist(2), pw, add_phshift
        integer :: fdims(3), lims(3,2), phys_cmat(2), win_corner(2), cyc_limsR(2,2),cyc_lims(3,2)
        integer :: iprec, i, j, sh, iwinsz, nyq, ind_in_stk, foffset, ok
        integer :: wdim, h, k, l, m, ll, mm, incr, icls, iptcl, interp_shlim, interp_shlim_sq
        integer :: first_iprec, first_stkind, fromp, top, istk, nptcls_in_stk, nstks, last_stkind, stkind
        integer :: ibatch, nbatches, istart, iend, ithr, nptcls_in_batch, first_pind, last_pind
        if( .not. params_glob%l_distr_exec ) write(logfhandle,'(a)') '>>> ASSEMBLING CLASS SUMS'
        ! init cavgs
        call init_cavgs_sums
        if( do_frac_update )then
            call cavger_readwrite_partial_sums( 'read' )
            call cavger_apply_weights( 1. - params_glob%update_frac )
        endif
        kbwin  = kbinterpol(KBWINSZ, params_glob%alpha)
        wdim   = kbwin%get_wdim()
        iwinsz = ceiling(kbwin%get_winsz() - 0.5)
        ! Number stacks
        first_pind = 0
        do i = 1,partsz
            if( precs(i)%pind > 0 )then
                first_pind = precs(i)%pind
                exit
            endif
        enddo
        call build_glob%spproj%map_ptcl_ind2stk_ind(params_glob%oritype, first_pind, first_stkind, ind_in_stk)
        last_pind = 0
        do i = partsz,1,-1
            if( precs(i)%pind > 0 )then
                last_pind = precs(i)%pind
                exit
            endif
        enddo
        call build_glob%spproj%map_ptcl_ind2stk_ind(params_glob%oritype, last_pind, last_stkind,  ind_in_stk)
        nstks = last_stkind - first_stkind + 1
        ! Objects allocations
        allocate(read_imgs(READBUFFSZ), cgrid_imgs(params_glob%nthr), cyc1(wdim), cyc2(wdim), w(wdim, wdim))
        !$omp parallel default(shared) proc_bind(close) private(i)
        !$omp do schedule(static)
        do i = 1,READBUFFSZ
            call read_imgs(i)%new(ldim, params_glob%smpd, wthreads=.false.)
        enddo
        !$omp end do nowait
        !$omp do schedule(static)
        do i = 1,params_glob%nthr
            call cgrid_imgs(i)%new(ldim_pd, params_glob%smpd, wthreads=.false.)
        enddo
        !$omp end do nowait
        !$omp end parallel
        lims            = cgrid_imgs(1)%loop_lims(2)
        cyc_lims        = cgrid_imgs(1)%loop_lims(3)
        nyq             = cgrid_imgs(1)%get_lfny(1)
        fdims           = cgrid_imgs(1)%get_array_shape()
        foffset         = fdims(2) / 2
        interp_shlim    = nyq + 1
        interp_shlim_sq = interp_shlim**2
        cyc_limsR(:,1)  = cyc_lims(1,:)  ! limits in fortran layered format
        cyc_limsR(:,2)  = cyc_lims(2,:)  ! to avoid copy on cyci_1d call
        allocate(cmats(fdims(1),fdims(2),READBUFFSZ), rhos(fdims(1),fdims(2),READBUFFSZ))
        ! Main loop
        first_iprec  = 1 ! first particle record in current stack
        do istk = first_stkind,last_stkind
            ! Particles range
            stkind = first_stkind + istk - 1
            fromp  = max(precs(1)%pind,      nint(build_glob%spproj%os_stk%get(istk,'fromp')))
            top    = min(precs(partsz)%pind, nint(build_glob%spproj%os_stk%get(istk,'top')))
            nptcls_in_stk = top - fromp + 1 ! # of particles in stack
            if( nptcls_in_stk == 0 )cycle
            call build_glob%spproj%get_stkname_and_ind(params_glob%oritype, fromp, stk_fname, ind_in_stk)
            ind_in_stk = ind_in_stk - 1 ! because incremented in the loop below
            ! open stack
            if( nptcls_in_stk < READBUFFSZ )then
                call stkio_r%open(stk_fname, smpd, 'read', bufsz=nptcls_in_stk)
            else
                call stkio_r%open(stk_fname, smpd, 'read', bufsz=READBUFFSZ)
            endif
            ! Read batches loop
            nbatches = ceiling(real(nptcls_in_stk)/real(READBUFFSZ)) ! will be 1 most of the tme
            do ibatch = 1,nbatches
                istart = (ibatch - 1)              * READBUFFSZ + 1  ! first index in current batch, will be 1      most of the time
                iend   = min(nptcls_in_stk, istart + READBUFFSZ - 1) ! last  index in current batch, will be nptcls most of the time
                nptcls_in_batch = iend - istart + 1                  ! nptcls_in_batch will be nptcls_in_stk most of the rime
                j      = 0                                           ! index in batch
                do i = istart,iend
                    iprec      = first_iprec + i - 1
                    ind_in_stk = ind_in_stk + 1
                    j          = j + 1
                    if( precs(iprec)%pind == 0 ) cycle
                    call stkio_r%read(ind_in_stk, read_imgs(j))
                enddo
                ! Interpolation loop
                !$omp parallel default(shared) proc_bind(close)&
                !$omp private(i,iprec,iptcl,fcomp,win_corner,add_phshift,mat,ithr,h,k,l,m,ll,mm,dist,loc,sh,phys_cmat,cyc1,cyc2,w,incr,pw)
                !$omp do schedule(static)
                do i = 1,nptcls_in_batch
                    iprec        = first_iprec + istart + i - 2
                    if( precs(iprec)%pind == 0 ) cycle
                    iptcl        = precs(iprec)%pind
                    ithr         = omp_get_thread_num()+1
                    cmats(:,:,i) = zero
                    rhos(:,:,i)  = 0.
                    ! normalize & pad & FFT
                    call read_imgs(i)%noise_norm_pad_fft(build_glob%lmsk, cgrid_imgs(ithr))
                    ! apply CTF, shift
                    add_phshift = 0.
                    if( phaseplate ) add_phshift = precs(iprec)%phshift
                    if( ctfflag /= CTFFLAG_NO )then
                        if( ctfflag == CTFFLAG_FLIP )then
                            call precs(iprec)%tfun%apply_and_shift(cgrid_imgs(ithr), 1, lims, rhos(:,:,i), -precs(iprec)%shift(1),&
                            &-precs(iprec)%shift(2), precs(iprec)%dfx, precs(iprec)%dfy, precs(iprec)%angast, add_phshift)
                        else
                            call precs(iprec)%tfun%apply_and_shift(cgrid_imgs(ithr), 2, lims, rhos(:,:,i), -precs(iprec)%shift(1),&
                            &-precs(iprec)%shift(2), precs(iprec)%dfx, precs(iprec)%dfy, precs(iprec)%angast, add_phshift)
                        endif
                    else
                        call precs(iprec)%tfun%apply_and_shift(cgrid_imgs(ithr), 3, lims, rhos(:,:,i), -precs(iprec)%shift(1),&
                        &-precs(iprec)%shift(2), precs(iprec)%dfx, precs(iprec)%dfy, precs(iprec)%angast, add_phshift)
                    endif
                    ! Rotation matrix
                    call rotmat2d(-precs(iprec)%e3, mat)
                    ! Interpolation
                    if( l_bilinear )then
                        ! bi-linear interpolation
                        do h=lims(1,1),lims(1,2)
                            do k=lims(2,1),lims(2,2)
                                sh = nint(hyp(real(h),real(k)))
                                if( sh > interp_shlim )cycle
                                loc = matmul(real([h,k]),mat)
                                ! interpolation
                                win_corner = floor(loc) ! bottom left corner
                                dist  = loc - real(win_corner)
                                l     = cyci_1d(cyc_limsR(:,1), win_corner(1))
                                ll    = cyci_1d(cyc_limsR(:,1), win_corner(1)+1)
                                m     = cyci_1d(cyc_limsR(:,2), win_corner(2))
                                mm    = cyci_1d(cyc_limsR(:,2), win_corner(2)+1)
                                fcomp =         (1.-dist(1))*(1.-dist(2)) * cgrid_imgs(ithr)%get_fcomp2D(l, m)  ! bottom left corner
                                fcomp = fcomp + (1.-dist(1))*dist(2)      * cgrid_imgs(ithr)%get_fcomp2D(l, mm) ! bottom right corner
                                fcomp = fcomp + dist(1)*(1.-dist(2))      * cgrid_imgs(ithr)%get_fcomp2D(ll,m)  ! upper left corner
                                fcomp = fcomp + dist(1)*dist(2)           * cgrid_imgs(ithr)%get_fcomp2D(ll,mm) ! upper right corner
                                ! set
                                phys_cmat = cgrid_imgs(ithr)%comp_addr_phys(h,k)
                                cmats(phys_cmat(1),phys_cmat(2),i) = fcomp
                            end do
                        end do
                    else
                        ! Kaiser-Bessel
                        do h=lims(1,1),lims(1,2)
                            do k=lims(2,1),lims(2,2)
                                sh = nint(hyp(real(h),real(k)))
                                if( sh > interp_shlim )cycle
                                loc = matmul(real([h,k]),mat)
                                win_corner = nint(loc) - iwinsz
                                ! weights kernel
                                w = 1.
                                do l=1,wdim
                                    incr = l - 1
                                    ! circular addresses
                                    cyc1(l) = cyci_1d(cyc_limsR(:,1), win_corner(1) + incr)
                                    cyc2(l) = cyci_1d(cyc_limsR(:,2), win_corner(2) + incr)
                                    ! interpolation kernel matrix
                                    w(l,:) = w(l,:) * kbwin%apod( real(win_corner(1) + incr) - loc(1) )
                                    w(:,l) = w(:,l) * kbwin%apod( real(win_corner(2) + incr) - loc(2) )
                                enddo
                                w = w / sum(w)
                                ! interpolation
                                fcomp = zero
                                do l=1,wdim
                                    do m=1,wdim
                                        if( w(l,m) < TINY ) cycle
                                        fcomp = fcomp + cgrid_imgs(ithr)%get_fcomp2D(cyc1(l),cyc2(m)) * w(l,m)
                                    end do
                                end do
                                ! set
                                phys_cmat = cgrid_imgs(ithr)%comp_addr_phys(h,k)
                                cmats(phys_cmat(1),phys_cmat(2),i) = fcomp
                            end do
                        end do
                    endif
                enddo
                !$omp end do nowait
                ! Sum over classes
                !$omp do schedule(static)
                do icls = 1,ncls
                    ithr = omp_get_thread_num() + 1
                    do i = 1,nptcls_in_batch
                        iprec = first_iprec + istart + i - 2
                        if( precs(iprec)%pind == 0 ) cycle
                        if( precs(iprec)%class == icls )then
                            pw = precs(iprec)%pw
                            select case(precs(iprec)%eo)
                                case(0,-1)
                                    call cavgs_even(icls)%get_cmat_ptr(pcmat(ithr)%cmat)
                                    call ctfsqsums_even(icls)%get_cmat_ptr(prhomat(ithr)%cmat)
                                case(1)
                                    call cavgs_odd(icls)%get_cmat_ptr(pcmat(ithr)%cmat)
                                    call ctfsqsums_odd(icls)%get_cmat_ptr(prhomat(ithr)%cmat)
                            end select
                            pcmat(ithr)%cmat(:,:,1)   = pcmat(ithr)%cmat(:,:,1)   + pw * cmats(:,:,i)
                            prhomat(ithr)%cmat(:,:,1) = prhomat(ithr)%cmat(:,:,1) + pw * cmplx(rhos(:,:,i),0.0)
                        endif
                    enddo
                enddo
                !$omp end do nowait
                !$omp end parallel
            enddo ! end read batches loop
            ! close stack
            call stkio_r%close
            ! Keeeping track of particles records
            first_iprec = first_iprec + nptcls_in_stk
        enddo
        ! performs final quadrant swap on e/o rhos
        !$omp parallel do schedule(static) private(icls,ithr,k,ok,h,cswap) default(shared) proc_bind(close)
        do icls = 1,ncls
            ithr = omp_get_thread_num() + 1
            call ctfsqsums_even(icls)%get_cmat_ptr(prhomat(ithr)%cmat)
            do k = 1,foffset
                ok = k + foffset
                do h = 1,fdims(1)
                    cswap = prhomat(ithr)%cmat(h,k,1)
                    prhomat(ithr)%cmat(h, k,1) = prhomat(ithr)%cmat(h,ok,1)
                    prhomat(ithr)%cmat(h,ok,1) = cswap
                end do
            end do
            call ctfsqsums_odd(icls)%get_cmat_ptr(prhomat(ithr)%cmat)
            do k = 1,foffset
                ok = k + foffset
                do h = 1,fdims(1)
                    cswap = prhomat(ithr)%cmat(h,k,1)
                    prhomat(ithr)%cmat(h, k,1) = prhomat(ithr)%cmat(h,ok,1)
                    prhomat(ithr)%cmat(h,ok,1) = cswap
                end do
            end do
        enddo
        !$omp end parallel do
        ! Cleanup
        do i = 1,READBUFFSZ
            call read_imgs(i)%kill
        enddo
        do i = 1,params_glob%nthr
            call cgrid_imgs(i)%kill
        enddo
        deallocate(read_imgs,cgrid_imgs)
        if( .not. params_glob%l_distr_exec ) call cavger_merge_eos_and_norm
    end subroutine cavger_assemble_sums

    !>  \brief  merges the even/odd pairs and normalises the sums
    subroutine cavger_merge_eos_and_norm
        type(image) :: gridcorrection_img
        integer     :: icls, eo_pop(2), pop
        call cavger_prep_gridding_correction(gridcorrection_img)
        !$omp parallel do default(shared) private(icls,eo_pop,pop) schedule(static) proc_bind(close)
        do icls=1,ncls
            eo_pop = prev_eo_pops(icls,:) + eo_class_pop(icls)
            pop    = sum(eo_pop)
            if(pop == 0)then
                call cavgs_merged(icls)%zero_and_unflag_ft
                call cavgs_even(icls)%zero_and_unflag_ft
                call cavgs_odd(icls)%zero_and_unflag_ft
                call ctfsqsums_merged(icls)%zero_and_flag_ft
            else
                call cavgs_merged(icls)%zero_and_flag_ft
                call cavgs_merged(icls)%add(cavgs_even(icls))
                call cavgs_merged(icls)%add(cavgs_odd(icls))
                call ctfsqsums_merged(icls)%zero_and_flag_ft
                call ctfsqsums_merged(icls)%add(ctfsqsums_even(icls))
                call ctfsqsums_merged(icls)%add(ctfsqsums_odd(icls))
                ! (w*CTF)**2 density correction
                if(eo_pop(1) > 1) call cavgs_even(icls)%ctf_dens_correct(ctfsqsums_even(icls))
                if(eo_pop(2) > 1) call cavgs_odd(icls)%ctf_dens_correct(ctfsqsums_odd(icls))
                if(pop > 1)       call cavgs_merged(icls)%ctf_dens_correct(ctfsqsums_merged(icls))
                call cavgs_even(icls)%ifft()
                call cavgs_odd(icls)%ifft()
                call cavgs_merged(icls)%ifft()
            endif
            call cavgs_even(icls)%clip_inplace(ldim)
            call cavgs_odd(icls)%clip_inplace(ldim)
            call cavgs_merged(icls)%clip_inplace(ldim)
            ! gridding correction
            call cavgs_even(icls)%div(gridcorrection_img)
            call cavgs_odd(icls)%div(gridcorrection_img)
            call cavgs_merged(icls)%div(gridcorrection_img)
        end do
        !$omp end parallel do
        call gridcorrection_img%kill
    end subroutine cavger_merge_eos_and_norm

    !>  \brief  calculates Fourier ring correlations
    subroutine cavger_calc_and_write_frcs_and_eoavg( fname )
        character(len=*), intent(in) :: fname
        type(image), allocatable     :: even_imgs(:), odd_imgs(:)
        real,        allocatable     :: frc(:)
        integer :: icls, find, find_plate
        allocate(even_imgs(ncls), odd_imgs(ncls), frc(filtsz))
        do icls=1,ncls
            call even_imgs(icls)%copy(cavgs_even(icls))
            call odd_imgs(icls)%copy(cavgs_odd(icls))
        end do
        !$omp parallel do default(shared) private(icls,frc,find,find_plate) schedule(static) proc_bind(close)
        do icls=1,ncls
            call even_imgs(icls)%mask(params_glob%msk, 'soft')
            call odd_imgs(icls)%mask(params_glob%msk, 'soft')
            call even_imgs(icls)%fft()
            call odd_imgs(icls)%fft()
            call even_imgs(icls)%fsc(odd_imgs(icls), frc)
            find_plate = 0
            if( phaseplate ) call phaseplate_correct_fsc(frc, find_plate)
            call build_glob%clsfrcs%set_frc(icls, frc, 1)
            ! average low-resolution info between eo pairs to keep things in register
            find = build_glob%clsfrcs%estimate_find_for_eoavg(icls, 1)
            find = max(find, find_plate)
            call cavgs_merged(icls)%fft()
            call cavgs_even(icls)%fft()
            call cavgs_odd(icls)%fft()
            call cavgs_even(icls)%insert_lowres_serial(cavgs_merged(icls), find)
            call cavgs_odd(icls)%insert_lowres_serial(cavgs_merged(icls), find)
            call cavgs_merged(icls)%ifft()
            call cavgs_even(icls)%ifft()
            call cavgs_odd(icls)%ifft()
        end do
        !$omp end parallel do
        ! write FRCs
        call build_glob%clsfrcs%write(fname)
        ! destruct
        do icls=1,ncls
            call even_imgs(icls)%kill
            call odd_imgs(icls)%kill
        end do
        deallocate(even_imgs, odd_imgs, frc)
    end subroutine cavger_calc_and_write_frcs_and_eoavg

    ! I/O

    !>  \brief  writes class averages to disk
    subroutine cavger_write( fname, which )
        character(len=*),  intent(in) :: fname, which
        integer :: icls
        select case(which)
            case('even')
                do icls=1,ncls
                    call cavgs_even(icls)%write(fname, icls)
                end do
            case('odd')
                do icls=1,ncls
                    call cavgs_odd(icls)%write(fname, icls)
                end do
            case('merged')
                    do icls=1,ncls
                    call cavgs_merged(icls)%write(fname, icls)
                end do
            case DEFAULT
                THROW_HARD('unsupported which flag')
        end select
        call update_stats

    contains

        subroutine update_stats
            integer :: icls, cnt
            real    :: stats(4),minv,maxv,meanv,stdevv
            logical :: l_err
            stats(1)   = huge(stats(1))
            stats(2)   = -stats(1)
            stats(3:4) = 0.
            cnt        = 1
            do icls = 1,ncls
                select case(which)
                    case('even')
                        call cavgs_even(icls)%stats(meanv, stdevv, maxv, minv, errout=l_err)
                    case('odd')
                        call cavgs_odd(icls)%stats(meanv, stdevv, maxv, minv, errout=l_err)
                    case('merged')
                        call cavgs_merged(icls)%stats(meanv, stdevv, maxv, minv, errout=l_err)
                end select
                if( .not.l_err )then
                    cnt = cnt + 1
                    stats(1) = min(stats(1),minv)
                    stats(2) = max(stats(2),maxv)
                    stats(3) = stats(3) + meanv
                    stats(4) = stats(4) + stdevv**2.
                endif
            enddo
            if( cnt > 1 )then
                ! updates header, size, stack & removes box file
                stats(3) = stats(3) / real(cnt)
                stats(4) = sqrt(stats(4) / real(cnt))
                select case(which)
                    case('even')
                        call cavgs_even(1)%update_header_stats(fname,stats)
                    case('odd')
                        call cavgs_odd(1)%update_header_stats(fname,stats)
                    case('merged')
                        call cavgs_merged(1)%update_header_stats(fname,stats)
                end select
            endif
        end subroutine update_stats

    end subroutine cavger_write

    !>  \brief  reads class averages from disk
    subroutine cavger_read( fname, which )
        character(len=*),  intent(in) :: fname, which
        type(stack_io) :: stkio_r
        integer        :: icls
        select case(which)
            case('even')
                call stkio_r%open(trim(fname), smpd, 'read', bufsz=ncls)
                do icls=1,ncls
                    call cavgs_even(icls)%new(ldim,smpd,wthreads=.false.)
                    call stkio_r%read(icls, cavgs_even(icls))
                end do
                call stkio_r%close
            case('odd')
                call stkio_r%open(trim(fname), smpd, 'read', bufsz=ncls)
                do icls=1,ncls
                    call cavgs_odd(icls)%new(ldim,smpd,wthreads=.false.)
                    call stkio_r%read(icls, cavgs_odd(icls))
                end do
                call stkio_r%close
            case('merged')
                call stkio_r%open(trim(fname), smpd, 'read', bufsz=ncls)
                do icls=1,ncls
                    call cavgs_merged(icls)%new(ldim,smpd,wthreads=.false.)
                    call stkio_r%read(icls, cavgs_merged(icls))
                end do
                call stkio_r%close
            case DEFAULT
                THROW_HARD('unsupported which flag')
        end select
    end subroutine cavger_read

    !>  \brief  writes partial class averages to disk (distributed execution)
    subroutine cavger_readwrite_partial_sums( which )
        character(len=*), intent(in)  :: which
        integer                       :: icls, ldim_here(3)
        character(len=:), allocatable :: cae, cao, cte, cto
        type(stack_io)                :: stkio(4)
        logical                       :: is_ft
        allocate(cae, source='cavgs_even_part'//int2str_pad(params_glob%part,params_glob%numlen)//params_glob%ext)
        allocate(cao, source='cavgs_odd_part'//int2str_pad(params_glob%part,params_glob%numlen)//params_glob%ext)
        allocate(cte, source='ctfsqsums_even_part'//int2str_pad(params_glob%part,params_glob%numlen)//params_glob%ext)
        allocate(cto, source='ctfsqsums_odd_part'//int2str_pad(params_glob%part,params_glob%numlen)//params_glob%ext)
        select case(trim(which))
            case('read')
                call stkio(1)%open(cae, smpd, 'read', bufsz=ncls)
                call stkio(2)%open(cao, smpd, 'read', bufsz=ncls)
                call stkio(3)%open(cte, smpd, 'read', bufsz=ncls)
                call stkio(4)%open(cto, smpd, 'read', bufsz=ncls)
                do icls=1,ncls
                    call stkio(1)%read(icls, cavgs_even(icls))
                    call stkio(2)%read(icls, cavgs_odd(icls))
                    call stkio(3)%read(icls, ctfsqsums_even(icls))
                    call stkio(4)%read(icls, ctfsqsums_odd(icls))
                end do
                call stkio(1)%close
                call stkio(2)%close
                call stkio(3)%close
                call stkio(4)%close
            case('write')
                is_ft = cavgs_even(1)%is_ft()
                ldim_here  = cavgs_even(1)%get_ldim()
                call stkio(1)%open(cae, smpd, 'write', bufsz=ncls, is_ft=is_ft, box=ldim_here(1))
                call stkio(2)%open(cao, smpd, 'write', bufsz=ncls, is_ft=is_ft, box=ldim_here(1))
                call stkio(3)%open(cte, smpd, 'write', bufsz=ncls, is_ft=is_ft, box=ldim_here(1))
                call stkio(4)%open(cto, smpd, 'write', bufsz=ncls, is_ft=is_ft, box=ldim_here(1))
                do icls=1,ncls
                    call stkio(1)%write(icls, cavgs_even(icls))
                    call stkio(2)%write(icls, cavgs_odd(icls))
                    call stkio(3)%write(icls, ctfsqsums_even(icls))
                    call stkio(4)%write(icls, ctfsqsums_odd(icls))
                end do
                call stkio(1)%close
                call stkio(2)%close
                call stkio(3)%close
                call stkio(4)%close
            case DEFAULT
                THROW_HARD('unknown which flag; only read & write supported; cavger_readwrite_partial_sums')
        end select
        deallocate(cae, cao, cte, cto)
    end subroutine cavger_readwrite_partial_sums

    subroutine cavger_apply_weights( w )
        real, intent(in) :: w
        integer :: icls
        !$omp parallel do default(shared) private(icls) schedule(static) proc_bind(close)
        do icls=1,ncls
            call cavgs_even(icls)%mul(w)
            call ctfsqsums_even(icls)%mul(w)
            call cavgs_odd(icls)%mul(w)
            call ctfsqsums_odd(icls)%mul(w)
        end do
        !$omp end parallel do
    end subroutine cavger_apply_weights

    !>  \brief  re-generates the object after distributed execution
    subroutine cavger_assemble_sums_from_parts
        complex(kind=c_float_complex), pointer :: cmat_ptr1(:,:,:) => null()
        complex(kind=c_float_complex), pointer :: cmat_ptr2(:,:,:) => null()
        complex(kind=c_float_complex), pointer :: cmat_ptr3(:,:,:) => null()
        complex(kind=c_float_complex), pointer :: cmat_ptr4(:,:,:) => null()
        integer(timer_int_kind)       ::  t_init,  t_io,  t_workshare_sum,  t_set_sums,  t_merge_eos_and_norm,  t_tot
        real(timer_int_kind)          :: rt_init, rt_io, rt_workshare_sum, rt_set_sums, rt_merge_eos_and_norm, rt_tot
        complex,          allocatable :: csums(:,:,:,:)
        character(len=:), allocatable :: cae, cao, cte, cto
        character(len=STDLEN)         :: benchfname
        type(image) :: imgs4read(4)
        integer     :: ipart, icls, array_shape(3), ldim_here(3), fnr
        if( L_BENCH_GLOB )then
            t_init = tic()
            t_tot  = t_init
        endif
        call init_cavgs_sums
        ! construct image objs for read/sum
        ldim_here    = ldim_pd
        ldim_here(3) = ncls
        call imgs4read(1)%new(ldim_here, smpd)
        call imgs4read(2)%new(ldim_here, smpd)
        call imgs4read(3)%new(ldim_here, smpd)
        call imgs4read(4)%new(ldim_here, smpd)
        call imgs4read(1)%set_ft(.true.)
        call imgs4read(2)%set_ft(.true.)
        call imgs4read(3)%set_ft(.true.)
        call imgs4read(4)%set_ft(.true.)
        call imgs4read(1)%get_cmat_ptr(cmat_ptr1)
        call imgs4read(2)%get_cmat_ptr(cmat_ptr2)
        call imgs4read(3)%get_cmat_ptr(cmat_ptr3)
        call imgs4read(4)%get_cmat_ptr(cmat_ptr4)
        ! construct complex matrices for parallel summation
        array_shape = imgs4read(1)%get_array_shape()
        allocate(csums(4,array_shape(1),array_shape(2),array_shape(3)), source=cmplx(0.,0.))
        if( L_BENCH_GLOB )then
            ! end of init
            rt_init = toc(t_init)
            ! initialise incremental timers before loop
            rt_io            = 0.
            rt_workshare_sum = 0.
        endif
        do ipart=1,params_glob%nparts
            if( L_BENCH_GLOB ) t_io = tic()
            ! look for files
            allocate(cae, source='cavgs_even_part'    //int2str_pad(ipart,params_glob%numlen)//params_glob%ext)
            allocate(cao, source='cavgs_odd_part'     //int2str_pad(ipart,params_glob%numlen)//params_glob%ext)
            allocate(cte, source='ctfsqsums_even_part'//int2str_pad(ipart,params_glob%numlen)//params_glob%ext)
            allocate(cto, source='ctfsqsums_odd_part' //int2str_pad(ipart,params_glob%numlen)//params_glob%ext)
            ! serial read
            call imgs4read(1)%read(cae)
            call imgs4read(2)%read(cao)
            call imgs4read(3)%read(cte)
            call imgs4read(4)%read(cto)
            deallocate(cae, cao, cte, cto)
            if( L_BENCH_GLOB )then
                rt_io = rt_io + toc(t_io)
                t_workshare_sum = tic()
            endif
            ! parallel summation
            !$omp parallel workshare proc_bind(close)
            csums(1,:,:,:) = csums(1,:,:,:) + cmat_ptr1(:,:,:)
            csums(2,:,:,:) = csums(2,:,:,:) + cmat_ptr2(:,:,:)
            csums(3,:,:,:) = csums(3,:,:,:) + cmat_ptr3(:,:,:)
            csums(4,:,:,:) = csums(4,:,:,:) + cmat_ptr4(:,:,:)
            !$omp end parallel workshare
            if( L_BENCH_GLOB ) rt_workshare_sum = rt_workshare_sum + toc(t_workshare_sum)
        end do
        if( L_BENCH_GLOB ) t_set_sums = tic()
        ! update image objects in parallel
        !$omp parallel do default(shared) private(icls) schedule(static) proc_bind(close)
        do icls=1,ncls
            call cavgs_even(icls)    %set_cmat(csums(1,:,:,icls))
            call cavgs_odd(icls)     %set_cmat(csums(2,:,:,icls))
            call ctfsqsums_even(icls)%set_cmat(csums(3,:,:,icls))
            call ctfsqsums_odd(icls) %set_cmat(csums(4,:,:,icls))
        end do
        !$omp end parallel do
        if( L_BENCH_GLOB ) rt_set_sums = rt_set_sums + toc(t_set_sums)
        ! destruct
        call imgs4read(1)%kill
        call imgs4read(2)%kill
        call imgs4read(3)%kill
        call imgs4read(4)%kill
        deallocate(csums)
        ! merge eo-pairs and normalize
        if( L_BENCH_GLOB ) t_merge_eos_and_norm = tic()
        call cavger_merge_eos_and_norm()
        if( L_BENCH_GLOB )then
            rt_merge_eos_and_norm = toc(t_merge_eos_and_norm)
            rt_tot                = toc(t_tot)
            benchfname = 'CAVGASSEMBLE_BENCH.txt'
            call fopen(fnr, FILE=trim(benchfname), STATUS='REPLACE', action='WRITE')
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,1x,f9.2)') 'initialisation       : ', rt_init
            write(fnr,'(a,1x,f9.2)') 'I/O                  : ', rt_io
            write(fnr,'(a,1x,f9.2)') 'workshare sum        : ', rt_workshare_sum
            write(fnr,'(a,1x,f9.2)') 'set sums             : ', rt_set_sums
            write(fnr,'(a,1x,f9.2)') 'merge eo-pairs & norm: ', rt_merge_eos_and_norm
            write(fnr,'(a,1x,f9.2)') 'total time           : ', rt_tot
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
            write(fnr,'(a,1x,f9.2)') 'initialisation        : ', (rt_init/rt_tot)               * 100.
            write(fnr,'(a,1x,f9.2)') 'I/O                   : ', (rt_io/rt_tot)                 * 100.
            write(fnr,'(a,1x,f9.2)') 'workshare sum         : ', (rt_workshare_sum/rt_tot)      * 100.
            write(fnr,'(a,1x,f9.2)') 'set sums              : ', (rt_set_sums/rt_tot)           * 100.
            write(fnr,'(a,1x,f9.2)') 'merge eo-pairs & norm : ', (rt_merge_eos_and_norm/rt_tot) * 100.
            write(fnr,'(a,1x,f9.2)') '% accounted for       : ',&
            &((rt_init+rt_io+rt_workshare_sum+rt_set_sums+rt_merge_eos_and_norm)/rt_tot) * 100.
            call fclose(fnr)
        endif
    end subroutine cavger_assemble_sums_from_parts

    !>  \brief  corrects for Fourier domain bilinear interpolation
    subroutine cavger_prep_gridding_correction( img )
        class(image), intent(inout) :: img
        real    :: center(3),dist(2),pid,sinc,pad_sc
        integer :: i,j
        call img%new(ldim,smpd)
        if( l_bilinear )then
            center = real(ldim/2 + 1)
            pad_sc = 1. / real(ldim_pd(1))
            do i = 1, ldim(1)
                dist(1) = pad_sc * (real(i) - center(1))
                do j = 1, ldim(2)
                    dist(2) = pad_sc * (real(j) - center(2))
                    pid     = PI * sqrt(sum(dist**2.))
                    if( pid < TINY )then
                        sinc = 1.
                    else
                        sinc = sin(pid) / pid
                    endif
                    call img%set([i,j,1], sinc*sinc)
                enddo
            enddo
        else
            img = 1.
        endif
    end subroutine cavger_prep_gridding_correction

    ! destructor

    !>  \brief  is a destructor
    subroutine cavger_kill
        integer :: icls
        if( exists )then
            do icls=1,ncls
                call cavgs_even(icls)%kill
                call cavgs_odd(icls)%kill
                call cavgs_merged(icls)%kill
                call ctfsqsums_even(icls)%kill
                call ctfsqsums_odd(icls)%kill
                call ctfsqsums_merged(icls)%kill
            end do
            deallocate( cavgs_even, cavgs_odd, cavgs_merged, ctfsqsums_even,&
            &ctfsqsums_odd, ctfsqsums_merged, pptcl_mask, prev_eo_pops)
            deallocate(precs)
            istart = 0
            iend   = 0
            partsz = 0
            ncls   = 0
            exists = .false.
        endif
    end subroutine cavger_kill

end module simple_classaverager
