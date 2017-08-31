module simple_classaverager


implicit none

public :: classaverager
private

type cavger_ptcl_record
    integer              :: pind      = 0
    integer              :: state_bal = 0
    integer              :: eo        = 1 ! even is 0, odd is 1, default is -1
    real                 :: pw        = 0.0
    integer              :: class     = 0
    integer, allocatable :: states(:)
    real,    allocatable :: ows(:)
    real,    allocatable :: e3s(:)
    real,    allocatable :: shifts(:,:)
end type cavger_ptcl_record


type classaverager
    private
    integer :: istart, iend
    type(cavger_ptcl_record), allocatable :: pinfo(:)
    type(image), allocatable :: cavgs_even(:)       !< class averages
    type(image), allocatable :: cavgs_odd(:)        ! -"-
    type(image), allocatable :: cavgs_merged(:)     ! -"-
    type(image), allocatable :: ctfsqsums_even(:)   !< CTF**2 sums for Wiener normalisation
    type(image), allocatable :: ctfsqsums_odd(:)    !< -"-
    type(image), allocatable :: ctfsqsums_merged(:) !< -"-
    logical :: exists = .false.
contains
    procedure :: new


end type classaverager

contains

    subroutine new( self, b, p, grid, prime3Dsrchobj )

    end subroutine new



    subroutine assemble_sums( b, p, grid )
        use simple_projector_hlev, only: rot_imgbatch
        use simple_map_reduce,     only: split_nobjs_even
        use simple_oris,           only: oris
        use simple_ctf,            only: ctf
        class(build),      intent(inout) :: b
        class(params),     intent(inout) :: p
        logical, optional, intent(in)    :: grid
        type(oris)               :: a_here, batch_oris
        type(ori)                :: orientation
        type(image)              :: batch_imgsum, cls_imgsum
        type(image), allocatable :: batch_imgs(:) 
        integer,     allocatable :: ptcls_inds(:), batches(:,:)
        logical,     allocatable :: batch_mask(:)
        real      :: w
        integer   :: icls, iptcl, istart, iend, inptcls, icls_pop
        integer   :: i, nbatches, batch, batchsz, cnt
        logical   :: l_grid
        integer, parameter :: BATCHTHRSZ = 20
        l_grid = .true.
        if( present(grid) ) l_grid = grid
        if( .not. p%l_distr_exec )then
            write(*,'(a)') '>>> ASSEMBLING CLASS SUMS'
        endif


        ! init
         do icls=1,p%ncls
            self%cavgs_even(icls)        = 0.
            self%cavgs_odd(icls)         = 0.
            self%cavgs_merged(icls)      = 0.
            self%ctfsqsums_even(icls)    = cmplx(0.,0.)
            self%ctfsqsums_odd(icls)     = cmplx(0.,0.)
            self%ctfsqsums_merged(icls)  = cmplx(0.,0.)
        end do

        ! work out range and oris
        if( p%l_distr_exec )then
            istart  = p%fromp
            iend    = p%top
            inptcls = iend - istart +1
            call a_here%new(inptcls)
            cnt = 0
            do iptcl = istart, iend
                cnt = cnt + 1
                call a_here%set_ori(cnt, b%a%get_ori(iptcl))
            enddo
        else
            istart  = 1
            iend    = p%nptcls
            inptcls = p%nptcls
            a_here  = b%a
        endif



        ! cluster loop
        do icls = 1, p%ncls
            call progress(icls,p%ncls)
            icls_pop = a_here%get_pop( icls, 'class' )
            if(icls_pop == 0)cycle
            call cls_imgsum%new([p%box, p%box, 1], p%smpd)
            ptcls_inds = a_here%get_pinds( icls, 'class' )
            ! batch planning
            nbatches = ceiling(real(icls_pop)/real(p%nthr*BATCHTHRSZ))
            batches  = split_nobjs_even(icls_pop, nbatches)
            ! batch loop
            do batch = 1, nbatches
                ! prep batch
                batchsz = batches(batch,2) - batches(batch,1) + 1
                allocate(batch_imgs(batchsz), batch_mask(batchsz))
                batch_mask = .true.
                call batch_oris%new(batchsz)
                ! batch particles loop
                do i = 1,batchsz
                    iptcl       = istart - 1 + ptcls_inds(batches(batch,1)+i-1)
                    orientation = b%a%get_ori(iptcl)
                    call batch_oris%set_ori(i, orientation)
                    ! stash images (this goes here or suffer bugs)
                    call read_img_from_stk( b, p, iptcl )
                    batch_imgs(i) = b%img
                    ! enforce state, balancing and weight exclusions
                    if( nint(orientation%get('state')) == 0 .or.&
                        &nint(orientation%get('state_balance')) == 0 .or.&
                        &orientation%get('w') < TINY )then
                        batch_mask(i) = .false.
                        cycle
                    endif
                    ! CTF square sum & shift
                    call apply_ctf_and_shift(batch_imgs(i), orientation)
                enddo
                if( l_grid )then
                    ! rotate batch by gridding
                    call rot_imgbatch(batch_imgs, batch_oris, batch_imgsum, p%msk, batch_mask)
                else
                    ! real space rotation
                    call batch_imgsum%new([p%box, p%box, 1], p%smpd)
                    do i = 1,batchsz
                        if( .not. batch_mask(i) ) cycle
                        iptcl = istart - 1 + ptcls_inds(batches(batch,1)+i-1)
                        orientation = b%a%get_ori(iptcl)
                        w = orientation%get('w')
                        call batch_imgs(i)%rtsq( -orientation%e3get(), 0., 0. )
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !call batch_imgsum%add(batch_imgs(i), w)
                    enddo
                endif
                ! batch summation
                call cls_imgsum%add( batch_imgsum )
                ! batch cleanup
                do i = 1, batchsz
                    call batch_imgs(i)%kill
                enddo
                deallocate(batch_imgs, batch_mask)
            enddo
            ! set class
            b%cavgs(icls) = cls_imgsum
            ! class cleanup
            deallocate(ptcls_inds)
        enddo
        if( .not.p%l_distr_exec ) call prime2D_norm_sums( b, p )

        contains

            !> image is shifted and Fted on exit and the class CTF square sum updated
            subroutine apply_ctf_and_shift( img, o )
                class(image), intent(inout) :: img
                class(ori),   intent(inout) :: o
                type(image) :: ctfsq
                type(ctf)   :: tfun
                real        :: dfx, dfy, angast, x, y, pw
                call ctfsq%new(img%get_ldim(), p%smpd)
                call ctfsq%set_ft(.true.)
                tfun = ctf(img%get_smpd(), o%get('kv'), o%get('cs'), o%get('fraca'))
                ! set CTF and shift parameters
                select case(p%tfplan%mode)
                    case('astig') ! astigmatic CTF
                        dfx    = o%get('dfx')
                        dfy    = o%get('dfy')
                        angast = o%get('angast')
                    case('noastig') ! non-astigmatic CTF
                        dfx    = o%get('dfx')
                        dfy    = dfx
                        angast = 0.
                end select
                x  = -o%get('x')
                y  = -o%get('y')
                pw = o%get('w')
                ! apply
                call img%fwd_ft
                ! take care of the nominator
                select case(p%tfplan%flag)
                    case('yes')  ! multiply with CTF
                        call tfun%apply_and_shift(img, ctfsq, x, y, dfx, 'ctf', dfy, angast)
                    case('flip') ! multiply with abs(CTF)
                        call tfun%apply_and_shift(img, ctfsq, x, y, dfx, 'abs', dfy, angast)
                    case('mul','no')
                        call tfun%apply_and_shift(img, ctfsq, x, y, dfx, '', dfy, angast)
                end select
                ! add to sum
                call  b%ctfsqsums(icls)%add(ctfsq, pw)
            end subroutine apply_ctf_and_shift

    end subroutine assemble_sums





end module simple_classaverager
