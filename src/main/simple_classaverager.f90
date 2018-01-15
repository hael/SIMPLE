module simple_classaverager
#include "simple_lib.f08"
use simple_ctf,     only: ctf
use simple_build,   only: build
use simple_params,  only: params
use simple_image,   only: image
use simple_timer    ! use all in there
implicit none

public :: cavger_new, cavger_transf_oridat, cavger_get_cavg, cavger_set_cavg, cavger_assemble_sums,&
cavger_merge_eos_and_norm, cavger_calc_and_write_frcs_and_eoavg, cavger_write, cavger_read,&
cavger_write_partial_sums, cavger_assemble_sums_from_parts, cavger_kill, cavgs_even, cavgs_odd, cavgs_merged
private

type ptcl_record
    type(ctf)            :: tfun                                !< transfer function
    integer              :: pind    = 0                         !< particle index in stack
    integer              :: eo      = -1                        !< even is 0, odd is 1, default is -1
    real                 :: pw      = 0.0                       !< particle weight
    real                 :: dfx     = 0.0                       !< defocus in x (microns)
    real                 :: dfy     = 0.0                       !< defocus in y (microns)
    real                 :: angast  = 0.0                       !< angle of astigmatism (in degrees)
    real                 :: phshift = 0.0                       !< additional phase shift from the Volta
    integer, allocatable :: classes(:)                          !< class assignments (can be many per particle in 3D case, hence the array)
    integer, allocatable :: states(:)                           !< state assignments
    integer, allocatable :: inpl_inds(:)                        !< in-plane rotation indices
    real,    allocatable :: ows(:)                              !< orientation weights
    real,    allocatable :: e3s(:)                              !< in-plane rotations
    real,    allocatable :: shifts(:,:)                         !< rotational origin shifts
end type ptcl_record

class(build),      pointer     :: bp => null()                  !< pointer to build
class(params),     pointer     :: pp => null()                  !< pointer to params
type(CTFFLAGTYPE)              :: ctfflag                       !< ctf flag <yes|no|mul|flip>
integer                        :: istart          = 0, iend = 0 !< particle index range
integer                        :: partsz          = 0           !< size of partition
integer                        :: ncls            = 0           !< # classes
integer                        :: nstates         = 0           !< # states
integer                        :: filtsz          = 0           !< size of filter function or FSC
integer                        :: ldim(3)         = [0,0,0]     !< logical dimension of image
integer                        :: ldim_pd(3)      = [0,0,0]     !< logical dimension of image, padded
real                           :: smpd            = 0.          !< sampling distance
type(ptcl_record), allocatable :: precs(:)                      !< particle records
type(image),       allocatable :: cavgs_even(:,:)               !< class averages
type(image),       allocatable :: cavgs_odd(:,:)                !< -"-
type(image),       allocatable :: cavgs_merged(:,:)             !< -"-
type(image),       allocatable :: ctfsqsums_even(:,:)           !< CTF**2 sums for Wiener normalisation
type(image),       allocatable :: ctfsqsums_odd(:,:)            !< -"-
type(image),       allocatable :: ctfsqsums_merged(:,:)         !< -"-
real,              allocatable :: inpl_rots(:)                  !< in-plane rotations (sign shifted)
logical                        :: phaseplate    = .false.       !< Volta phaseplate images or not
logical                        :: l_is_class    = .true.        !< for prime2D or not
logical                        :: l_hard_assign = .true.        !< npeaks == 1 or not
logical                        :: exists        = .false.       !< to flag instance existence

integer, parameter      :: BATCHTHRSZ = 50
logical, parameter      :: L_BENCH    = .false.
integer(timer_int_kind) :: t_batch_loop, t_gridding, t_tot
real(timer_int_kind)    :: rt_batch_loop, rt_gridding, rt_tot
character(len=STDLEN)   :: benchfname

contains

    !>  \brief  is a constructor
    !!          data is now managed so that all exclusions are taken care of here
    !!          which means properly balanced batches can be produced for both soft
    !!          and hard clustering solutions
    subroutine cavger_new( b, p, which )
        class(build),  target, intent(inout) :: b     !< builder
        class(params), target, intent(inout) :: p     !< params
        character(len=*),      intent(in)    :: which !< class/proj
        integer :: alloc_stat, istate, icls
        ! destruct possibly pre-existing instance
        call cavger_kill
        ! set pointers
        bp => b
        pp => p
        ! set nstates
        nstates = p%nstates
        ! class or proj
        select case(which)
            case('class')
                l_is_class = .true.
                ncls       = p%ncls
            case('proj')
                l_is_class = .false.
                ! possible reduction of # projection directions used
                ! for the class average representation
                ncls = min(NSPACE_BALANCE,p%nspace)
            case DEFAULT
                stop 'unsupported which flag; simple_classaverager :: cavger_new'
        end select
        ! work out range and partsz
        if( p%l_distr_exec )then
            istart = p%fromp
            iend   = p%top
            partsz = iend - istart + 1
        else
            istart = 1
            iend   = p%nptcls
            partsz = p%nptcls
        endif
        ! CTF logics
        select case(p%ctf)
            case('no')
                ctfflag%flag = CTFFLAG_NO
            case('yes')
                ctfflag%flag = CTFFLAG_YES
            case('mul')
                stop 'ERROR ctf=mul deprecated; simple_classaverager_dev :: cavger_new'
            case('flip')
                ctfflag%flag = CTFFLAG_FLIP
        end select
        phaseplate = p%tfplan%l_phaseplate
        ! smpd
        smpd       = p%smpd
        ! set ldims
        ldim       = [pp%box,pp%box,1]
        ldim_pd    = [pp%boxpd,pp%boxpd,1]
        ldim_pd(3) = 1
        filtsz     = b%img%get_filtsz()
        ! build arrays
        allocate(precs(partsz), cavgs_even(p%nstates,ncls), cavgs_odd(p%nstates,ncls),&
        &cavgs_merged(p%nstates,ncls), ctfsqsums_even(p%nstates,ncls),&
        &ctfsqsums_odd(p%nstates,ncls), ctfsqsums_merged(p%nstates,ncls), stat=alloc_stat)
        call alloc_errchk('cavger_new; simple_classaverager_dev', alloc_stat)
        do istate=1,p%nstates
            do icls=1,ncls
                call cavgs_even(istate,icls)%new(ldim,p%smpd,wthreads=.false.)
                call cavgs_odd(istate,icls)%new(ldim,p%smpd,wthreads=.false.)
                call cavgs_merged(istate,icls)%new(ldim,p%smpd,wthreads=.false.)
                call ctfsqsums_even(istate,icls)%new(ldim,p%smpd,wthreads=.false.)
                call ctfsqsums_odd(istate,icls)%new(ldim,p%smpd,wthreads=.false.)
                call ctfsqsums_merged(istate,icls)%new(ldim,p%smpd,wthreads=.false.)
            end do
        end do
        ! flag existence
        exists = .true.
    end subroutine cavger_new

    ! setters/getters

    !>  \brief  transfers orientation data to the instance
    subroutine cavger_transf_oridat( a )
        use simple_ori,  only: ori
        use simple_oris, only: oris
        class(oris), intent(in)    :: a
        real, allocatable :: ori_weights(:)
        type(ori)         :: orientation
        type(oris)        :: a_here
        integer           :: alloc_stat, cnt, n_incl, iori
        integer           :: cnt_ori, istate, icls, iptcl
        logical           :: l_reduce_projs
        ! create a copy of a that can be modified
        a_here = a
        cnt    = 0
        ! fetch data from a_here
        do iptcl=istart,iend
            cnt = cnt + 1
            ! exclusion condition
            if( nint(a_here%get(iptcl,'state')) == 0 .or.&
               &nint(a_here%get(iptcl,'state_balance')) == 0 .or.&
               &a_here%get(iptcl,'w') < TINY )then
                precs(cnt)%pind  = 0
                cycle
            endif
            ! parameter transfer
            precs(cnt)%pind  = iptcl
            precs(cnt)%eo    = nint(a_here%get(iptcl,'eo'))
            precs(cnt)%pw    = a_here%get(iptcl,'w')
            precs(cnt)%tfun  = ctf(pp%smpd, a_here%get(iptcl,'kv'), a_here%get(iptcl,'cs'), a_here%get(iptcl,'fraca'))
            select case(pp%tfplan%mode)
                case('astig') ! astigmatic CTF
                    precs(cnt)%dfx    = a_here%get(iptcl,'dfx')
                    precs(cnt)%dfy    = a_here%get(iptcl,'dfy')
                    precs(cnt)%angast = a_here%get(iptcl,'angast')
                case('noastig') ! non-astigmatic CTF
                    precs(cnt)%dfx    = a_here%get(iptcl,'dfx')
                    precs(cnt)%dfy    = precs(cnt)%dfx
                    precs(cnt)%angast = 0.
            end select
            precs(cnt)%phshift = 0.
            if( phaseplate ) precs(cnt)%phshift = a_here%get(iptcl,'phshift')
        end do
        l_hard_assign = .true.
        cnt = 0
        do iptcl=istart,iend
            cnt = cnt + 1
            ! inclusion condition
            if( precs(cnt)%pind > 0 )then
                ! allocate & set info in record
                if( allocated(precs(cnt)%classes) )  deallocate(precs(cnt)%classes)
                if( allocated(precs(cnt)%inpl_inds)) deallocate(precs(cnt)%inpl_inds)
                if( allocated(precs(cnt)%states)  )  deallocate(precs(cnt)%states)
                if( allocated(precs(cnt)%ows)     )  deallocate(precs(cnt)%ows)
                if( allocated(precs(cnt)%e3s)     )  deallocate(precs(cnt)%e3s)
                if( allocated(precs(cnt)%shifts)  )  deallocate(precs(cnt)%shifts)
                allocate( precs(cnt)%classes(1),  precs(cnt)%states(1),&
                          precs(cnt)%ows(1),      precs(cnt)%e3s(1),&
                          precs(cnt)%shifts(1,2), precs(cnt)%inpl_inds(1), stat=alloc_stat )
                call alloc_errchk('cavger_new; simple_classaverager, record arrays', alloc_stat)
                precs(cnt)%classes(1)   = nint(a_here%get(iptcl, 'class'))
                precs(cnt)%inpl_inds(1) = nint(a_here%get(iptcl, 'inpl'))
                precs(cnt)%states(1)    = nint(a_here%get(iptcl, 'state'))
                precs(cnt)%ows(1)       = a_here%get(iptcl, 'w')
                precs(cnt)%e3s(1)       = a_here%e3get(iptcl)
                precs(cnt)%shifts(1,1)  = a_here%get(iptcl, 'x')
                precs(cnt)%shifts(1,2)  = a_here%get(iptcl, 'y')
            endif
        end do
        ! endif
        call a_here%kill
    end subroutine cavger_transf_oridat

    !>  \brief  is for initialization of the sums
    subroutine init_cavgs_sums
        integer :: istate, icls
        do istate=1,nstates
            do icls=1,ncls
                call cavgs_even(istate,icls)%new(ldim,smpd,wthreads=.false.)
                call cavgs_odd(istate,icls)%new(ldim,smpd,wthreads=.false.)
                call cavgs_merged(istate,icls)%new(ldim,smpd,wthreads=.false.)
                call cavgs_even(istate,icls)%zero_and_flag_ft
                call cavgs_odd(istate,icls)%zero_and_flag_ft
                call cavgs_merged(istate,icls)%zero_and_flag_ft
                call ctfsqsums_even(istate,icls)%zero_and_flag_ft
                call ctfsqsums_odd(istate,icls)%zero_and_flag_ft
                call ctfsqsums_merged(istate,icls)%zero_and_flag_ft
            end do
        end do
    end subroutine init_cavgs_sums

    !>  \brief  is for getting allocatable arrays with particle/record/ori indices
    subroutine get_indices( state, class, pinds, iprecs, ioris )
        integer,              intent(in)  :: state, class
        integer, allocatable, intent(out) :: pinds(:)
        integer, allocatable, intent(out) :: iprecs(:)
        integer, allocatable, intent(out) :: ioris(:)
        integer :: pop, alloc_stat, i, sz, iprec, cnt
        logical, allocatable :: l_state_class(:)
        pop = class_pop(state, class)
        if( allocated(pinds) )  deallocate(pinds)
        if( allocated(iprecs) ) deallocate(iprecs)
        if( allocated(ioris)  ) deallocate(ioris)
        allocate(pinds(pop), iprecs(pop), ioris(pop), stat=alloc_stat)
        call alloc_errchk('get_iprecs_ioris; simple_classaverager', alloc_stat)
        cnt = 0
        do iprec=1,partsz
            if( allocated(precs(iprec)%classes) )then
                sz = size(precs(iprec)%classes)
                allocate(l_state_class(sz))
                where( precs(iprec)%states .eq. state .and. precs(iprec)%classes .eq. class )
                    l_state_class = .true.
                else where
                    l_state_class = .false.
                endwhere
                if( any(l_state_class) )then
                    do i=1,sz
                        if( l_state_class(i) )then
                            cnt = cnt + 1
                            pinds(cnt)  = precs(iprec)%pind
                            iprecs(cnt) = iprec
                            ioris(cnt)  = i
                        endif
                    enddo
                endif
                deallocate(l_state_class)
            endif
        end do
    end subroutine get_indices

    !>  \brief  is for calculating class population
    function class_pop( state, class ) result( pop )
        integer, intent(in) :: state, class
        integer :: pop, iprec, sz
        logical, allocatable :: l_state_class(:)
        pop = 0
        do iprec=1,partsz
            if( allocated(precs(iprec)%classes) )then
                sz = size(precs(iprec)%classes)
                allocate(l_state_class(sz))
                where( precs(iprec)%states .eq. state .and. precs(iprec)%classes .eq. class )
                    l_state_class = .true.
                else where
                    l_state_class = .false.
                endwhere
                pop = pop + count(l_state_class)
                deallocate(l_state_class)
            endif
        end do
    end function class_pop

    !>  \brief  is for getting a class average
    subroutine cavger_get_cavg( class, which, img, state )
        integer,              intent(in)    :: class
        character(len=*),     intent(in)    :: which
        class(image),         intent(inout) :: img
        integer, optional,    intent(in)    :: state
        integer :: sstate
        sstate = 1
        if( present(state) ) sstate = state
        select case(which)
            case('even')
                call img%copy(cavgs_even(sstate,class))
            case('odd')
                call img%copy(cavgs_odd(sstate,class))
            case('merged')
                call img%copy(cavgs_merged(sstate,class))
            case DEFAULT
                stop 'unsupported which flag; simple_classaverager :: cavger_get_cavg'
        end select
    end subroutine cavger_get_cavg

    !>  \brief  is for setting a class average
    subroutine cavger_set_cavg( class, which, img, state )
        integer,              intent(in)    :: class
        character(len=*),     intent(in)    :: which
        class(image),         intent(in)    :: img
        integer, optional,    intent(in)    :: state
        integer :: sstate
        sstate = 1
        if( present(state) ) sstate = state
        select case(which)
            case('even')
                call cavgs_even(sstate,class)%copy(img)
            case('odd')
                call cavgs_odd(sstate,class)%copy(img)
            case('merged')
                call cavgs_merged(sstate,class)%copy(img)
            case DEFAULT
                stop 'unsupported which flag; simple_classaverager :: cavger_set_cavg'
        end select
    end subroutine cavger_set_cavg

    ! calculators

    !>  \brief  is for assembling the sums in distributed/non-distributed mode
    !!          using gridding interpolation in Fourier space or quadratic
    !!          interpolation in real-space
    subroutine cavger_assemble_sums
        use simple_kbinterpol,      only: kbinterpol
        use simple_prep4cgrid,      only: prep4cgrid
        use simple_map_reduce,      only: split_nobjs_even
        use simple_hadamard_common, only: read_img
        type(kbinterpol)         :: kbwin
        type(prep4cgrid)         :: gridprep
        type(image)              :: cls_imgsum_even, cls_imgsum_odd
        type(image), allocatable :: batch_imgs(:), cgrid_imgs(:)
        complex,     allocatable :: cmat_even(:,:,:), cmat_odd(:,:,:)
        real,        allocatable :: rho(:,:), rho_even(:,:), rho_odd(:,:), w(:,:)
        integer,     allocatable :: ptcls_inds(:), batches(:,:), iprecs(:)
        integer,     allocatable :: ioris(:), cyc1(:), cyc2(:)
        complex   :: zero
        real      :: loc(2), mat(2,2), pw, add_phshift
        integer   :: cnt_progress, nbatches, batch, icls_pop, iprec, iori, i, batchsz, fnr, sh
        integer   :: lims(3,2), istate, nyq, logi(3), phys(3), win(2,2), lims_small(3,2), phys_cmat(3)
        integer   :: cyc_lims(3,2), alloc_stat, wdim, h, k, l, m, incr, icls, iptcl, batchsz_max
        if( .not. pp%l_distr_exec ) write(*,'(a)') '>>> ASSEMBLING CLASS SUMS'
        ! init
        call init_cavgs_sums
        kbwin = kbinterpol(KBWINSZ, pp%alpha)
        zero  = cmplx(0.,0.)
        wdim  = kbwin%get_wdim()
        ! determines max batch size
        batchsz_max = 0
        do istate=1,nstates
            ! class loop
            do icls=1,ncls
                ! batch planning
                icls_pop = class_pop(istate, icls)
                if( icls_pop < 2 ) cycle
                nbatches = ceiling(real(icls_pop)/real(pp%nthr*BATCHTHRSZ))
                batches  = split_nobjs_even(icls_pop, nbatches)
                ! batch loop
                do batch=1,nbatches
                    ! prep batch
                    batchsz = batches(batch,2) - batches(batch,1) + 1
                    if( batchsz > batchsz_max ) batchsz_max = batchsz
                end do
            end do
        end do
        if( allocated(batches) ) deallocate(batches)
        ! pre-allocations
        allocate(batch_imgs(batchsz_max), cgrid_imgs(batchsz_max),&
                &cyc1(wdim), cyc2(wdim), w(wdim, wdim))
        do i=1,batchsz_max
            call batch_imgs(i)%new(ldim, pp%smpd,    wthreads=.false.)
            call cgrid_imgs(i)%new(ldim_pd, pp%smpd, wthreads=.false.)
        end do
        lims_small = batch_imgs(1)%loop_lims(2)
        lims       = cgrid_imgs(1)%loop_lims(2)
        cyc_lims   = cgrid_imgs(1)%loop_lims(3)
        cmat_even  = cgrid_imgs(1)%get_cmat()
        cmat_odd   = cgrid_imgs(1)%get_cmat()
        nyq        = cgrid_imgs(1)%get_lfny(1)
        allocate( rho(lims_small(1,1):lims_small(1,2),lims_small(2,1):lims_small(2,2)),&
                  rho_even(lims_small(1,1):lims_small(1,2),lims_small(2,1):lims_small(2,2)),&
                 &rho_odd( lims_small(1,1):lims_small(1,2),lims_small(2,1):lims_small(2,2)), stat=alloc_stat)
        call gridprep%new(bp%img, kbwin, ldim_pd)
        if( L_BENCH )then
            rt_batch_loop = 0.
            rt_gridding   = 0.
            rt_tot        = 0.
            t_tot         = tic()
        endif
        ! state loop
        cnt_progress = 0
        do istate=1,nstates
            ! class loop
            do icls=1,ncls
                cnt_progress = cnt_progress + 1
                call progress(cnt_progress, nstates * ncls)
                icls_pop = class_pop(istate, icls)
                if( icls_pop == 0 ) cycle
                call get_indices(istate, icls, ptcls_inds, iprecs, ioris)
                ! class temporary matrices
                cmat_even = zero
                cmat_odd  = zero
                rho       = 0.
                rho_even  = 0.
                rho_odd   = 0.
                ! batch planning
                nbatches = ceiling(real(icls_pop)/real(pp%nthr*BATCHTHRSZ))
                batches  = split_nobjs_even(icls_pop, nbatches)
                ! batch loop, prep
                do batch=1,nbatches
                    ! prep batch
                    batchsz = batches(batch,2) - batches(batch,1) + 1
                    ! read images
                    if( L_BENCH ) t_batch_loop = tic()
                    do i=1,batchsz
                        iptcl = ptcls_inds(batches(batch,1) + i - 1)
                        call read_img( bp, pp, iptcl )
                        batch_imgs(i) = bp%img
                    enddo
                    ! batch particles loop
                    if( L_BENCH ) rt_batch_loop = rt_batch_loop + toc(t_batch_loop)
                    if( L_BENCH ) t_gridding = tic()
                    !$omp parallel do default(shared) schedule(static) reduction(+:cmat_even,cmat_odd,rho_even,rho_odd) proc_bind(close)&
                    !$omp private(i,iprec,iori,add_phshift,rho,pw,mat,h,k,l,m,loc,sh,win,logi,phys,phys_cmat,cyc1,cyc2,w,incr)
                    ! batch loop, direct Fourier interpolation
                    do i=1,batchsz
                        iprec = iprecs(batches(batch,1) + i - 1)
                        iori  = ioris(batches(batch,1)  + i - 1)
                        ! normalise and FFT
                        call batch_imgs(i)%norm
                        call batch_imgs(i)%fwd_ft
                        ! apply CTF and shift
                        if( phaseplate )then
                            add_phshift = precs(iprec)%phshift
                        else
                            add_phshift = 0.
                        endif
                        if( ctfflag%flag /= CTFFLAG_NO )then
                            if( ctfflag%flag == CTFFLAG_FLIP )then
                                call precs(iprec)%tfun%apply_and_shift(batch_imgs(i), 1, lims_small, rho, -precs(iprec)%shifts(iori,1),&
                                    &-precs(iprec)%shifts(iori,2), precs(iprec)%dfx, precs(iprec)%dfy,&
                                    &precs(iprec)%angast, add_phshift)
                            else
                                call precs(iprec)%tfun%apply_and_shift(batch_imgs(i), 2, lims_small, rho, -precs(iprec)%shifts(iori,1),&
                                    &-precs(iprec)%shifts(iori,2), precs(iprec)%dfx, precs(iprec)%dfy,&
                                    &precs(iprec)%angast, add_phshift)
                            endif
                        else
                            call precs(iprec)%tfun%apply_and_shift(batch_imgs(i), 3, lims_small, rho, -precs(iprec)%shifts(iori,1),&
                                &-precs(iprec)%shifts(iori,2), precs(iprec)%dfx, precs(iprec)%dfy,&
                                &precs(iprec)%angast, add_phshift)
                        endif
                        ! prep weight
                        if( l_hard_assign )then
                            pw = precs(iprec)%pw
                        else
                            pw = precs(iprec)%pw * precs(iprec)%ows(iori)
                        endif
                        ! sampling density update
                        select case(precs(iprec)%eo)
                            case(0,-1)
                                rho_even = rho_even + pw * rho
                            case(1)
                                rho_odd  = rho_odd + pw * rho
                        end select
                        ! reverse FFT and prepare for gridding
                        call batch_imgs(i)%bwd_ft
                        call gridprep%prep_serial(batch_imgs(i), cgrid_imgs(i))
                        ! rotation
                        mat = rotmat2d( -precs(iprec)%e3s(iori) )
                        ! Fourier components loop
                        do h=lims(1,1),lims(1,2)
                            do k=lims(2,1),lims(2,2)
                                sh = nint(hyp(real(h),real(k)))
                                if( sh > nyq + 1 )cycle
                                loc = matmul(real([h,k]),mat)
                                call sqwin_2d(loc(1),loc(2), KBWINSZ, win)
                                w = pw
                                do l=1,wdim
                                    incr = l - 1
                                    ! circular addresses
                                    cyc1(l) = cyci_1d(cyc_lims(1,:), win(1,1) + incr)
                                    cyc2(l) = cyci_1d(cyc_lims(2,:), win(2,1) + incr)
                                    ! interpolation kernel matrix
                                    w(l,:) = w(l,:) * kbwin%apod( real(win(1,1) + incr) - loc(1) )
                                    w(:,l) = w(:,l) * kbwin%apod( real(win(2,1) + incr) - loc(2) )
                                enddo
                                ! point of addition
                                phys_cmat = cgrid_imgs(i)%comp_addr_phys([h,k,0])
                                select case(precs(iprec)%eo)
                                    case(0,-1)
                                        ! interpolation
                                        do l=1,wdim
                                            do m=1,wdim
                                                if( w(l,m) == 0. ) cycle
                                                logi       = [cyc1(l), cyc2(m), 0]
                                                phys       = cgrid_imgs(i)%comp_addr_phys(logi)
                                                cmat_even(phys_cmat(1),phys_cmat(2),phys_cmat(3)) = cmat_even(phys_cmat(1),phys_cmat(2),phys_cmat(3)) +&
                                                    &cgrid_imgs(i)%get_fcomp(logi, phys) * w(l,m)
                                            end do
                                        end do
                                    case(1)
                                        ! interpolation
                                        do l=1,wdim
                                            do m=1,wdim
                                                if( w(l,m) == 0. ) cycle
                                                logi       = [cyc1(l), cyc2(m), 0]
                                                phys       = cgrid_imgs(i)%comp_addr_phys(logi)
                                                cmat_odd(phys_cmat(1),phys_cmat(2),phys_cmat(3)) = cmat_odd(phys_cmat(1),phys_cmat(2),phys_cmat(3)) +&
                                                    &cgrid_imgs(i)%get_fcomp(logi, phys) * w(l,m)
                                            end do
                                        end do
                                end select
                            end do
                        end do
                    enddo
                    !$omp end parallel do
                    if( L_BENCH ) rt_gridding = rt_gridding + toc(t_gridding)
                enddo ! batch loop
                ! put back cmats
                call cls_imgsum_even%new(ldim_pd, pp%smpd)
                call cls_imgsum_odd%new(ldim_pd, pp%smpd)
                call cls_imgsum_even%set_cmat(cmat_even)
                call cls_imgsum_odd%set_cmat(cmat_odd)
                ! real space & clipping
                call cls_imgsum_even%bwd_ft
                call cls_imgsum_odd%bwd_ft
                call cls_imgsum_even%clip_inplace(ldim)
                call cls_imgsum_odd%clip_inplace(ldim)
                ! back to Fourier space
                call cls_imgsum_even%fwd_ft
                call cls_imgsum_odd%fwd_ft
                ! updates cavgs & rhos
                !$omp parallel do collapse(2) default(shared) private(h,k,logi,phys)&
                !$omp schedule(static) proc_bind(close)
                do h=lims_small(1,1),lims_small(1,2)
                    do k=lims_small(2,1),lims_small(2,2)
                        logi = [h,k,0]
                        phys = cls_imgsum_even%comp_addr_phys([h,k,0])
                        call cavgs_even(istate,icls)%set_fcomp(logi, phys, cls_imgsum_even%get_fcomp(logi, phys))
                        call cavgs_odd(istate,icls)%set_fcomp( logi, phys, cls_imgsum_odd%get_fcomp(logi, phys))
                        call ctfsqsums_even(istate,icls)%set_fcomp(logi, phys, cmplx(rho_even(h,k),0.))
                        call ctfsqsums_odd(istate,icls)%set_fcomp( logi, phys, cmplx(rho_odd(h,k), 0.))
                    enddo
                enddo
                !$omp end parallel do
                deallocate(ptcls_inds, batches, iprecs, ioris)
            enddo ! class loop
        enddo
        ! batch cleanup
        call cls_imgsum_even%kill
        call cls_imgsum_odd%kill
        call gridprep%kill
        do i=1,batchsz_max
            call batch_imgs(i)%kill
            call cgrid_imgs(i)%kill
        enddo
        if( allocated(cmat_even) ) deallocate(cmat_even)
        if( allocated(cmat_odd)  ) deallocate(cmat_odd)
        deallocate(rho, rho_even, rho_odd, batch_imgs, cgrid_imgs, cyc1, cyc2, w)
        if( .not. pp%l_distr_exec ) call cavger_merge_eos_and_norm
        if( L_BENCH )then
            rt_tot = rt_tot + toc(t_tot)
            benchfname = 'CLASSAVERAGER_BENCH.txt'
            call fopen(fnr, FILE=trim(benchfname), STATUS='REPLACE', action='WRITE')
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,1x,f9.2)') 'batch loop : ', rt_batch_loop
            write(fnr,'(a,1x,f9.2)') 'gridding   : ', rt_gridding
            write(fnr,'(a,1x,f9.2)') 'total time : ', rt_tot
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
            write(fnr,'(a,1x,f9.2)') 'batch loop : ', (rt_batch_loop/rt_tot) * 100.
            write(fnr,'(a,1x,f9.2)') 'gridding   : ', (rt_gridding/rt_tot)   * 100.
            call fclose(fnr)
        endif
    end subroutine cavger_assemble_sums

    !>  \brief  merges the even/odd pairs and normalises the sums
    subroutine cavger_merge_eos_and_norm
        integer :: istate, icls
        !$omp parallel do default(shared) private(istate,icls) schedule(static) proc_bind(close)
        do icls=1,ncls
            do istate=1,nstates
                call cavgs_merged(istate,icls)%zero_and_flag_ft
                call cavgs_merged(istate,icls)%add(cavgs_even(istate,icls))
                call cavgs_merged(istate,icls)%add(cavgs_odd(istate,icls))
                call ctfsqsums_merged(istate,icls)%zero_and_flag_ft
                call ctfsqsums_merged(istate,icls)%add(ctfsqsums_even(istate,icls))
                call ctfsqsums_merged(istate,icls)%add(ctfsqsums_odd(istate,icls))
                ! (w*CTF)**2 density correction
                call cavgs_even(istate,icls)%ctf_dens_correct(ctfsqsums_even(istate,icls))
                call cavgs_even(istate,icls)%bwd_ft
                call cavgs_even(istate,icls)%norm
                call cavgs_odd(istate,icls)%ctf_dens_correct(ctfsqsums_odd(istate,icls))
                call cavgs_odd(istate,icls)%bwd_ft
                call cavgs_odd(istate,icls)%norm
                call cavgs_merged(istate,icls)%ctf_dens_correct(ctfsqsums_merged(istate,icls))
                call cavgs_merged(istate,icls)%bwd_ft
                call cavgs_merged(istate,icls)%norm
            end do
        end do
        !$omp end parallel do
    end subroutine cavger_merge_eos_and_norm

    !>  \brief  calculates Fourier ring correlations
    subroutine cavger_calc_and_write_frcs_and_eoavg( fname )
        character(len=*), intent(in) :: fname
        type(image), allocatable     :: even_imgs(:,:), odd_imgs(:,:)
        real,        allocatable     :: frc(:)
        integer :: istate, icls, find, find_plate
        ! serial code for allocation/copy
        allocate(even_imgs(nstates,ncls), odd_imgs(nstates,ncls), frc(filtsz))
        do istate=1,nstates
            do icls=1,ncls
                call even_imgs(istate,icls)%copy(cavgs_even(istate,icls))
                call odd_imgs(istate,icls)%copy(cavgs_odd(istate,icls))
            end do
        end do
        ! parallel loop to do the job
        !$omp parallel do default(shared) private(istate,icls,frc,find,find_plate) schedule(static) proc_bind(close)
        do icls=1,ncls
            do istate=1,nstates
                call even_imgs(istate,icls)%norm
                call odd_imgs(istate,icls)%norm
                if( pp%l_innermsk )then
                    call even_imgs(istate,icls)%mask(pp%msk, 'soft', inner=pp%inner, width=pp%width)
                    call odd_imgs(istate,icls)%mask(pp%msk, 'soft', inner=pp%inner, width=pp%width)
                else
                    call even_imgs(istate,icls)%mask(pp%msk, 'soft')
                    call odd_imgs(istate,icls)%mask(pp%msk, 'soft')
                endif
                call even_imgs(istate,icls)%fwd_ft
                call odd_imgs(istate,icls)%fwd_ft
                call even_imgs(istate,icls)%fsc(odd_imgs(istate,icls), frc)
                find_plate = 0
                if( phaseplate ) call phaseplate_correct_fsc(frc, find_plate)
                call bp%projfrcs%set_frc(icls, frc, istate)
                ! average low-resolution info between eo pairs to keep things in register
                find = bp%projfrcs%estimate_find_for_eoavg(icls, istate)
                find = max(find, find_plate)
                call cavgs_merged(istate,icls)%fwd_ft
                call cavgs_even(istate,icls)%fwd_ft
                call cavgs_odd(istate,icls)%fwd_ft
                call cavgs_even(istate,icls)%insert_lowres_serial(cavgs_merged(istate,icls), find)
                call cavgs_odd(istate,icls)%insert_lowres_serial(cavgs_merged(istate,icls), find)
                call cavgs_merged(istate,icls)%bwd_ft
                call cavgs_even(istate,icls)%bwd_ft
                call cavgs_odd(istate,icls)%bwd_ft
                ! destruct
                call even_imgs(istate,icls)%kill
                call odd_imgs(istate,icls)%kill
            end do
        end do
        !$omp end parallel do
        ! write FRCs
        call bp%projfrcs%write(fname)
        ! destruct
        deallocate(even_imgs, odd_imgs, frc)
    end subroutine cavger_calc_and_write_frcs_and_eoavg

    ! I/O

    !>  \brief  writes class averages to disk
    subroutine cavger_write( fname, which, state )
        character(len=*),  intent(in)    :: fname, which
        integer, optional, intent(in)    :: state
        integer :: sstate, icls
        sstate = 1
        if( present(state) ) sstate = state
        select case(which)
            case('even')
                do icls=1,ncls
                    call cavgs_even(sstate, icls)%write(fname, icls)
                end do
            case('odd')
                do icls=1,ncls
                    call cavgs_odd(sstate, icls)%write(fname, icls)
                end do
            case('merged')
                 do icls=1,ncls
                    call cavgs_merged(sstate, icls)%write(fname, icls)
                end do
            case DEFAULT
                stop 'unsupported which flag; simple_classaverager :: cavger_get_cavg'
        end select
    end subroutine cavger_write

    !>  \brief  reads class averages from disk
    subroutine cavger_read( fname, which, state )
        character(len=*),  intent(in) :: fname, which
        integer, optional, intent(in) :: state
        integer :: sstate, icls
        if( .not. file_exists(fname) )then
            write(*,*) 'file does not exist in cwd: ', trim(fname)
            stop 'simple_classaverager :: read'
        endif
        sstate = 1
        if( present(state) ) sstate = state
        select case(which)
            case('even')
                do icls=1,ncls
                    call cavgs_even(sstate,icls)%new(ldim,smpd)
                    call cavgs_even(sstate, icls)%read(fname, icls)
                end do
            case('odd')
                do icls=1,ncls
                    call cavgs_odd(sstate,icls)%new(ldim,smpd)
                    call cavgs_odd(sstate, icls)%read(fname, icls)
                end do
            case('merged')
                 do icls=1,ncls
                    call cavgs_merged(sstate,icls)%new(ldim,smpd)
                    call cavgs_merged(sstate, icls)%read(fname, icls)
                end do
            case DEFAULT
                stop 'unsupported which flag; simple_classaverager :: cavger_read'
        end select
    end subroutine cavger_read

    !>  \brief  writes partial class averages to disk (distributed execution)
    subroutine cavger_write_partial_sums
        integer                       :: istate, icls
        character(len=:), allocatable :: cae, cao, cte, cto
        do istate=1,nstates
            if( nstates > 1 )then
                allocate(cae, source='cavgs'//'_state'//int2str_pad(istate,2)//'_even_part'//int2str_pad(pp%part,pp%numlen)//pp%ext)
                allocate(cao, source='cavgs'//'_state'//int2str_pad(istate,2)//'_odd_part'//int2str_pad(pp%part,pp%numlen)//pp%ext)
                allocate(cte, source='ctfsqsums'//'_state'//int2str_pad(istate,2)//'_even_part'//int2str_pad(pp%part,pp%numlen)//pp%ext)
                allocate(cto, source='ctfsqsums'//'_state'//int2str_pad(istate,2)//'_odd_part'//int2str_pad(pp%part,pp%numlen)//pp%ext)
            else
                allocate(cae, source='cavgs_even_part'//int2str_pad(pp%part,pp%numlen)//pp%ext)
                allocate(cao, source='cavgs_odd_part'//int2str_pad(pp%part,pp%numlen)//pp%ext)
                allocate(cte, source='ctfsqsums_even_part'//int2str_pad(pp%part,pp%numlen)//pp%ext)
                allocate(cto, source='ctfsqsums_odd_part'//int2str_pad(pp%part,pp%numlen)//pp%ext)
            endif
            do icls=1,ncls
                call cavgs_even(istate, icls)%write(cae, icls)
                call cavgs_odd(istate, icls)%write(cao, icls)
                call ctfsqsums_even(istate, icls)%write(cte, icls)
                call ctfsqsums_odd(istate, icls)%write(cto, icls)
            end do
            deallocate(cae, cao, cte, cto)
        end do
    end subroutine cavger_write_partial_sums

    !>  \brief  re-generates the object after distributed execution
    subroutine cavger_assemble_sums_from_parts
        type(image), allocatable :: imgs4read(:)
        character(len=:), allocatable :: cae, cao, cte, cto
        integer :: ipart, istate, icls
        call init_cavgs_sums
        allocate(imgs4read(4))
        call imgs4read(1)%new(ldim, smpd)
        call imgs4read(1)%set_ft(.true.)
        call imgs4read(2)%new(ldim, smpd)
        call imgs4read(2)%set_ft(.true.)
        call imgs4read(3)%new(ldim, smpd)
        call imgs4read(3)%set_ft(.true.)
        call imgs4read(4)%new(ldim, smpd)
        call imgs4read(4)%set_ft(.true.)
        do istate=1,nstates
            do ipart=1,pp%nparts
                if( nstates > 1 )then
                    allocate(cae, source='cavgs'//'_state'//int2str_pad(istate,2)//'_even_part'//int2str_pad(ipart,pp%numlen)//pp%ext)
                    allocate(cao, source='cavgs'//'_state'//int2str_pad(istate,2)//'_odd_part'//int2str_pad(ipart,pp%numlen)//pp%ext)
                    allocate(cte, source='ctfsqsums'//'_state'//int2str_pad(istate,2)//'_even_part'//int2str_pad(ipart,pp%numlen)//pp%ext)
                    allocate(cto, source='ctfsqsums'//'_state'//int2str_pad(istate,2)//'_odd_part'//int2str_pad(ipart,pp%numlen)//pp%ext)
                else
                    allocate(cae, source='cavgs_even_part'//int2str_pad(ipart,pp%numlen)//pp%ext)
                    allocate(cao, source='cavgs_odd_part'//int2str_pad(ipart,pp%numlen)//pp%ext)
                    allocate(cte, source='ctfsqsums_even_part'//int2str_pad(ipart,pp%numlen)//pp%ext)
                    allocate(cto, source='ctfsqsums_odd_part'//int2str_pad(ipart,pp%numlen)//pp%ext)
                endif
                if( .not. file_exists(cae) )then
                    write(*,*) 'File does not exists: ', trim(cae)
                    stop 'In: simple_classaverager :: cavger_assemble_sums_from_parts'
                endif
                if( .not. file_exists(cao) )then
                    write(*,*) 'File does not exists: ', trim(cao)
                    stop 'In: simple_classaverager :: cavger_assemble_sums_from_parts'
                endif
                if( .not. file_exists(cte) )then
                    write(*,*) 'File does not exists: ', trim(cte)
                    stop 'In: simple_classaverager :: cavger_assemble_sums_from_parts'
                endif
                if( .not. file_exists(cto) )then
                    write(*,*) 'File does not exists: ', trim(cto)
                    stop 'In: simple_classaverager :: cavger_assemble_sums_from_parts'
                endif
                do icls=1,ncls
                    call imgs4read(1)%read(cae, icls)
                    call imgs4read(2)%read(cao, icls)
                    call imgs4read(3)%read(cte, icls)
                    call imgs4read(4)%read(cto, icls)
                    call cavgs_even(istate,icls)%add_workshare(imgs4read(1), cavgs_odd(istate,icls),imgs4read(2),&
                        &ctfsqsums_even(istate,icls), imgs4read(3), ctfsqsums_odd(istate,icls), imgs4read(4))
                end do
                deallocate(cae, cao, cte, cto)
            end do
        end do
        call imgs4read(1)%kill
        call imgs4read(2)%kill
        call imgs4read(3)%kill
        call imgs4read(4)%kill
        deallocate(imgs4read)
        call cavger_merge_eos_and_norm()
    end subroutine cavger_assemble_sums_from_parts

    ! destructor

    !>  \brief  is a destructor
    subroutine cavger_kill
        integer :: istate, icls, iprec
        if( exists )then
            bp => null()
            pp => null()
            do istate=1,nstates
                do icls=1,ncls
                    call cavgs_even(istate,icls)%kill
                    call cavgs_odd(istate,icls)%kill
                    call cavgs_merged(istate,icls)%kill
                    call ctfsqsums_even(istate,icls)%kill
                    call ctfsqsums_odd(istate,icls)%kill
                    call ctfsqsums_merged(istate,icls)%kill
                end do
            end do
            deallocate( cavgs_even, cavgs_odd, cavgs_merged,&
            &ctfsqsums_even, ctfsqsums_odd, ctfsqsums_merged)
            do iprec=1,partsz
                if( allocated(precs(iprec)%classes) ) deallocate(precs(iprec)%classes)
                if( allocated(precs(iprec)%states)  ) deallocate(precs(iprec)%states)
                if( allocated(precs(iprec)%ows)     ) deallocate(precs(iprec)%ows)
                if( allocated(precs(iprec)%e3s)     ) deallocate(precs(iprec)%e3s)
                if( allocated(precs(iprec)%shifts)  ) deallocate(precs(iprec)%shifts)
            end do
            deallocate(precs)
            if( allocated(inpl_rots) ) deallocate(inpl_rots)
            istart        = 0
            iend          = 0
            partsz        = 0
            ncls          = 0
            nstates       = 0
            l_is_class    = .true.
            l_hard_assign = .true.
            exists        = .false.
        endif
    end subroutine cavger_kill

end module simple_classaverager
