!>  \brief  SIMPLE projector class
module simple_projector
use simple_winfuns   ! singleton
use simple_defs      ! singleton
use simple_image,    only: image
use simple_ori,      only: ori
use simple_oris,     only: oris
use simple_params,   only: params
use simple_gridding, only: prep4cgrid
use simple_jiffys,   only: alloc_err, progress
use simple_math,     only: recwin_3d, euclid
implicit none

public :: projector
private

logical, parameter :: debug = .false.
logical, parameter :: warn  = .false.

type :: projector
    private
    type(winfuns)         :: wfuns         !< window functions object
    character(len=STDLEN) :: wfun_str='kb' !< window function, string descriptor
    character(len=STDLEN) :: imgkind='em'  !< image kind descriptor
    real                  :: winsz=1.5     !< window half-width
    real                  :: alpha=2.      !< oversampling ratio
    real                  :: harwin=2.     !< rounded window half-width
  contains
    ! GETTER
    procedure :: get_wfuns
    ! FOURIER PROJECTORS
    procedure          :: fproject
    procedure, private :: fproject_polar_1
    procedure, private :: fproject_polar_2
    generic            :: fproject_polar => fproject_polar_1, fproject_polar_2
    procedure, private :: img2polarft_1
    procedure, private :: img2polarft_2
    generic            :: img2polarft => img2polarft_1, img2polarft_2
    procedure          :: extr_gridfcomp
    ! REAL-SPACE PROJECTOR
    procedure          :: rproject
    ! HIGH-LEVEL PROJECTORS
    procedure          :: projvol
    procedure          :: rotvol
    procedure          :: rotimg
    procedure          :: fprojvol_polar
end type projector

interface projector
    module procedure :: constructor
end interface

contains

    ! CONSTRUCTOR

    !> \brief  is a constructor
    function constructor( wfun, winsz, alpha, imgkind ) result( self )
        use simple_params, only: params
        character(len=*), optional, intent(in) :: wfun    !< window function, string descriptor
        real,             optional, intent(in) :: winsz   !< window half-width
        real,             optional, intent(in) :: alpha   !< oversampling ratio
        character(len=*), optional, intent(in) :: imgkind !< image kind (em/xfel) 
        type(projector)                        :: self    !< instance
        if( present(wfun) )then
            self%wfun_str = wfun
            if( self%wfun_str .eq. 'nn' )then
                self%winsz  = 0.5
                self%harwin = real(ceiling(self%winsz))
            endif
        endif
        if( present(winsz) )then
            self%winsz = winsz
            self%harwin = real(ceiling(winsz))
        endif
        if( present(alpha) )then
            self%alpha = alpha
        endif
        if( present(imgkind) )then
            self%imgkind = imgkind
        endif
        self%wfuns = winfuns(self%wfun_str,self%winsz,self%alpha)
    end function constructor

    ! GETTER
    
    !>  \brief  return the window functions used by reconstructor
    function get_wfuns( self ) result( wfs )
        class(projector), intent(inout) :: self
        type(winfuns) :: wfs
        wfs = self%wfuns
    end function get_wfuns
    
    ! FOURIER PROJECTORS
    
    !> \brief  extracts a Fourier plane from a volume
    subroutine fproject( self, fvol, e, fplane, lp_dyn )
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(projector), intent(inout) :: self
        class(image),     intent(inout) :: fvol
        class(ori),       intent(in)    :: e
        class(image),     intent(inout) :: fplane
        real, optional,   intent(in)    :: lp_dyn
        complex :: comp
        integer :: h, hh, k, kk, lims(3,2), ldim(3), ldim_ptcl(3), sqarg, sqlp
        real    :: vec(3), loc(3)
        if( warn)then
            if( .not. fvol%exists() )   stop 'fvol needs to exists before call; fproject_1; simple_projector'
            if( .not. fvol%is_ft() )    stop 'volume needs to be FTed before call; fproject_1; simple_projector'
            if( .not. fplane%exists() ) stop 'fplane needs to exists before call; fproject_1; simple_projector'
            ldim = fvol%get_ldim()
            ldim_ptcl = fplane%get_ldim()
            if( ldim(3) == 1 )          stop 'only for interpolation from 3D mages; fproject_1; simple_projector'
            if( ldim(1) == ldim_ptcl(1) .and. ldim(2) == ldim_ptcl(2) )then
            else
                print *, 'ldim1 vol/ptcl:', ldim(1), ldim_ptcl(1)
                print *, 'ldim2 vol/ptcl:', ldim(2), ldim_ptcl(2)
                stop 'nonconformable dims btw vol & ptcl; fproject_1; simple_projector'
            endif
        endif
        call fplane%set_ft(.true.)
        if( present(lp_dyn) )then
            lims = fvol%loop_lims(1,lp_dyn)
        else
            lims = fvol%loop_lims(2) ! Nyqvist default low-pass limit
        endif
        sqlp = (maxval(lims(:,2)))**2
        !$omp parallel do schedule(auto) shared(lims,sqlp) private(h,k,hh,kk,sqarg,vec,loc,comp)
        do h=lims(1,1),lims(1,2)
            hh = h*h
            do k=lims(2,1),lims(2,2)
                kk = k*k
                sqarg = hh+kk
                if( sqarg <= sqlp )then
                    vec(1) = real(h)
                    vec(2) = real(k)
                    vec(3) = 0.
                    loc = matmul(vec,e%get_mat())
                    comp = self%extr_gridfcomp(fvol,loc)
                else
                    comp = cmplx(0.,0.)
                endif
                call fplane%set_fcomp([h,k,0],comp)
            end do
        end do
        !$omp end parallel do
    end subroutine fproject

    !> \brief  extracts a polar FT from a volume
    subroutine fproject_polar_1( self, fvol, e, pimg, memoize )
        use simple_polarft, only: polarft
        use simple_math,    only: deg2rad
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(projector),  intent(inout) :: self    !< projector object
        class(image),      intent(inout) :: fvol    !< Fourier volume
        class(ori),        intent(inout) :: e       !< orientation
        class(polarft),    intent(inout) :: pimg    !< polar image
        logical, optional, intent(in)    :: memoize !< memoize or not
        integer :: i, k, ldim(3), pdims(3), ldim_polft(3)
        real    :: vec(3), loc(3)
        complex :: comp
        logical :: mmemoize
        ldim = fvol%get_ldim()
        if( ldim(3) == 1 )        stop 'only for interpolation from 3D images; fproject_polar_1; simple_projector'
        if( .not. fvol%is_ft() )  stop 'volume needs to be FTed before call; fproject_polar_1; simple_projector'
        if( .not. pimg%exists() ) stop 'polarft object needs to be created; fproject_polar_1; simple_projector'
        mmemoize = .true.
        if( present(memoize) ) mmemoize = memoize
        ldim_polft(1:2) = ldim(1:2)
        ldim_polft(3)   = 1
        call pimg%set_ldim(ldim_polft)
        pdims = pimg%get_dims()
        !$omp parallel do schedule(auto) default(shared) private(i,k,vec,loc,comp)
        do i=1,pdims(1)
            do k=pdims(2),pdims(3)
                vec(:2) = pimg%get_coord(i,k)
                vec(3)  = 0.
                loc     = matmul(vec,e%get_mat())
                comp    = self%extr_gridfcomp(fvol,loc)
                call pimg%set_fcomp(i, k, comp)
            end do
        end do
        !$omp end parallel do
        ! calculate the square sums required for correlation calculation
        if( mmemoize ) call pimg%memoize_sqsums
    end subroutine fproject_polar_1

    !> \brief  extracts a polar FT from a volume
    subroutine fproject_polar_2( self, iref, fvol, e, p, pftcc )
        use simple_polarft, only: polarft
        use simple_polarft_corrcalc, only: polarft_corrcalc
        use simple_params, only: params
        class(projector),        intent(inout) :: self  !< projector instance
        integer,                 intent(in)    :: iref  !< logical reference index [1,nrefs]
        class(image),            intent(inout) :: fvol  !< Fourier volume
        class(ori),              intent(inout) :: e     !< orientation
        class(params),           intent(in)    :: p     !< parameters class
        class(polarft_corrcalc), intent(inout) :: pftcc !< polarft_corrcalc object to be filled
        type(polarft) :: pimg
        complex :: comp
        integer :: irot, k, ldim_fvol(3), ldim_pft(3), refsz
        ldim_fvol = fvol%get_ldim()
        ldim_pft = pftcc%get_ldim()
        if( .not. all(ldim_fvol(1:2) == ldim_pft(1:2)) ) stop 'nonconforming logical dimensions; fproject_polar_2; simple_projector'
        call pimg%new([p%kfromto(1),p%kfromto(2)],p%ring2,ptcl=.false.)
        call self%fproject_polar_1(fvol, e, pimg, memoize=.false.)
        refsz = pftcc%get_refsz()
        do irot=1,refsz
            do k=p%kfromto(1),p%kfromto(2)
                comp = pimg%get_fcomp(irot,k)
                call pftcc%set_ref_fcomp(iref, irot, k, comp)
            end do
        end do
        call pftcc%memoize_sqsum_ref(iref)
        call pimg%kill
    end subroutine fproject_polar_2
    
    !> \brief  transfers a 2D Cartesian image to polar coordinates Fourier coordinatess
    subroutine img2polarft_1( self, img, pimg )
        use simple_polarft,  only: polarft
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(projector), intent(inout) :: self !< instance
        class(image),     intent(inout) :: img  !< image to be transformed
        class(polarft),   intent(inout) :: pimg !< polarft ouput
        integer :: i, k, ldim(3), pdims(3)
        real    :: vec(3)
        complex :: comp
        ldim = img%get_ldim()
        call pimg%set_ldim(ldim)
        if( ldim(3) > 1 )         stop 'only for interpolation from 2D images; img2polarft_1; simple_projector'
        if( .not. img%is_ft() )   stop 'image needs to FTed before this operation; simple_projector::img2polarft_1'
        if( .not. pimg%exists() ) stop 'polarft object needs to be created before call; img2polarft_1; simple_projector'
        pdims = pimg%get_dims()
        !$omp parallel do schedule(auto) default(shared) private(i,k,vec,comp)
        do i=1,pdims(1)
            do k=pdims(2),pdims(3)
                vec(:2) = pimg%get_coord(i,k)
                vec(3)  = 0.
                comp    = self%extr_gridfcomp(img,vec)
                call pimg%set_fcomp(i, k, comp)
            end do
        end do
        !$omp end parallel do
        ! calculate the square sums required for correlation calculation
        call pimg%memoize_sqsums
    end subroutine img2polarft_1
    
    !> \brief  transfers a 2D Cartesian image to polar coordinates Fourier coordinatess
    subroutine img2polarft_2( self, ind, img, pftcc, isptcl )
        use simple_polarft_corrcalc, only: polarft_corrcalc
        use gnufor2
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(projector),        intent(inout) :: self   !< instance
        integer,                 intent(in)    :: ind    !< logical particle index [fromp,top]
        class(image),            intent(inout) :: img    !< Cartesian image
        class(polarft_corrcalc), intent(inout) :: pftcc  !< polarft_corrcalc object to be filled
        logical, optional,       intent(in)    :: isptcl !< to indicate whether particle or reference
        complex, allocatable :: pft(:,:)
        integer :: i, k, ldim_img(3), ldim_pft(3), pdims(3), alloc_stat
        real    :: vec(3)
        logical :: iisptcl
        iisptcl = .true.
        if( present(isptcl) )iisptcl = isptcl
        ldim_img    = img%get_ldim()
        if( ldim_img(3) > 1 )      stop 'only for interpolation from 2D images; img2polarft_2; simple_projector'
        if( .not. pftcc%exists() ) stop 'polarft_corrcalc object needs to be created; img2polarft_2; simple_projector'
        ldim_pft = pftcc%get_ldim()
        if( .not. all(ldim_img == ldim_pft) )then
            print *, 'ldim_img: ', ldim_img
            print *, 'ldim_pft: ', ldim_pft
            stop 'logical dimensions do not match; img2polarft_2; simple_projector'
        endif
        if( .not. img%is_ft() ) stop 'image needs to FTed before this operation; simple_projector::img2polarft_2'
        pdims(1)   = pftcc%get_nrots()
        pdims(2:3) = pftcc%get_kfromto()
        allocate( pft(pdims(1),pdims(2):pdims(3)), stat=alloc_stat )
        call alloc_err("In: img2polarft_2; simple_projector", alloc_stat)
        !$omp parallel do schedule(auto) default(shared) private(i,k,vec)
        do i=1,pdims(1)
            do k=pdims(2),pdims(3)
                vec(:2)  = pftcc%get_coord(i,k)
                vec(3)   = 0.
                pft(i,k) = self%extr_gridfcomp(img,vec)
            end do
        end do
        !$omp end parallel do
        if( iisptcl )then
            call pftcc%set_ptcl_pft(ind, pft)
        else
            call pftcc%set_ref_pft(ind, pft)
        endif
        ! kill the remains
        deallocate(pft)
    end subroutine img2polarft_2
    
    !> \brief  extracts a Fourier component from a transform by gridding
    function extr_gridfcomp( self, fimg, loc ) result( comp_sum )
        use simple_math, only: cyci_1d
        class(projector), intent(inout)       :: self
        class(image), intent(inout)           :: fimg
        real, intent(in)                      :: loc(3)
        real,    allocatable :: w1(:), w2(:), w3(:)
        integer, allocatable :: cyc1(:), cyc2(:), cyc3(:)
        integer           :: alloc_stat, i, j, m
        integer           :: lims(3,2), win(3,2)
        complex           :: comp, comp_sum, zero
        real              :: harwin_here
        harwin_here = 2.
        zero        = cmplx(0.,0.)
        comp_sum    = zero
        win         = recwin_3d(loc(1),loc(2),loc(3),harwin_here)
        lims        = fimg%loop_lims(3)
        allocate( w1(win(1,1):win(1,2)), w2(win(2,1):win(2,2)), &
            & cyc1(win(1,1):win(1,2)), cyc2(win(2,1):win(2,2)), stat=alloc_stat)
        call alloc_err("In: extr_gridfcomp; simple_projector, 1", alloc_stat)
        do i=win(1,1),win(1,2)
            w1(i)   = self%wfuns%eval_apod(real(i)-loc(1))
            cyc1(i) = cyci_1d( lims(1,:),i )
        end do
        do j=win(2,1),win(2,2)
            w2(j)   = self%wfuns%eval_apod(real(j)-loc(2))
            cyc2(j) = cyci_1d( lims(2,:),j )
        end do
        if( fimg%is_2d() )then
            ! 2D
            do i=win(1,1),win(1,2)
                if( w1(i) == 0. ) cycle
                do j=win(2,1),win(2,2)
                    if( w2(j) == 0. ) cycle
                    comp = fimg%get_fcomp( [cyc1(i),cyc2(j),0] )
                    if( comp .eq. zero ) cycle
                    comp_sum = comp_sum+comp*w1(i)*w2(j)
                end do
            end do
            deallocate( w1,w2,cyc1,cyc2 )
        else
            ! 3D
            allocate( w3(win(3,1):win(3,2)), cyc3(win(3,1):win(3,2)), stat=alloc_stat)
            call alloc_err("In: extr_gridfcomp; simple_projector, 2", alloc_stat)
            do m=win(3,1),win(3,2)
                w3(m)   = self%wfuns%eval_apod(real(m)-loc(3))
                cyc3(m) = cyci_1d( lims(3,:),m )
            end do
            do i=win(1,1),win(1,2)
                if( w1(i) == 0. ) cycle
                do j=win(2,1),win(2,2)
                    if( w2(j) == 0. ) cycle
                    do m=win(3,1),win(3,2)
                        if( w3(m) == 0. ) cycle
                        comp = fimg%get_fcomp( [cyc1(i),cyc2(j),cyc3(m)] )
                        if( comp .eq. zero ) cycle
                        comp_sum = comp_sum+comp*w1(i)*w2(j)*w3(m)
                    end do
                end do
            end do
            deallocate(w1,w2,w3,cyc1,cyc2,cyc3)
        endif
    end function extr_gridfcomp
    
    ! REAL-SPACE PROJECTOR
    
    !>  \brief  project a 3d volume onto a 2d slice (real-space projection from Frealix)
    subroutine rproject(self, vol, e, img, maxrad )
        use simple_math,  only: recwin_3d, hyp
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(projector), intent(inout) :: self   !< projector instance
        type(image),      intent(inout) :: vol    !< volume to be projected
        class(ori),       intent(in)    :: e      !< Euler angle
        type(image),      intent(inout) :: img    !< resulting projection image
        real,             intent(in)    :: maxrad !< project inside this radius
        real, parameter :: harwin_here=1.5
        integer :: i, j, k, orig(3), win(3,2), s, t, u, ldim(3), sosh, tosh, uosh
        real    :: inp_coos(3), out_coos(3), rcomp, rcomp_sum
        real    :: voxval, w, mmaxrad, sqmaxrad, rad(3)
        logical :: didft
        if( self%imgkind .eq. 'xfel' ) stop 'routine not intended for xfel-kind images; simple_projector::rproject'
        ldim = vol%get_ldim()
        if( ldim(3) == 1 )             stop 'only for interpolation from 3D mages; project; simple_projector'
        if( .not. vol%even_dims() )    stop 'even dimensions assumed; project; simple_projector'
        didft = .false.
        if( vol%is_ft() )then
            call vol%bwd_ft
            didft = .true.
        endif
        call img%new([ldim(1),ldim(2),1],vol%get_smpd())
        ! initialize
        orig(1) = ldim(1)/2+1
        orig(2) = ldim(2)/2+1
        orig(3) = ldim(3)/2+1
        inp_coos = 0.
        out_coos = 0.
        mmaxrad = min(maxrad,real(ldim(1))/2-1-harwin_here)
        sqmaxrad = mmaxrad**2
        rad = 0.
        !$omp parallel default(shared) private(j,out_coos,rad,i,k,inp_coos,voxval,win,rcomp_sum,s,sosh,t,tosh,u,uosh,rcomp,w)
        !$omp do schedule(auto)
        do j=1,ldim(2)
            out_coos(2) = real(j-orig(2))
            rad(2) = out_coos(2)**2.
            do i=1,ldim(1)
                out_coos(1) = real(i-orig(1))
                rad(1) = rad(2)+out_coos(1)**2
                ! check whether we are within the radius
                if( rad(1) .le. sqmaxrad )then
                    ! iterate over all the voxels in the ray which lands on the current pixel
                    do k=1,ldim(3)
                        out_coos(3) = real(k-orig(3))
                        rad(3) = rad(1)+out_coos(3)**2
                        ! check that we are within the radius
                        if(rad(3) .le. sqmaxrad) then
                            inp_coos = matmul(out_coos,e%get_mat())
                            ! get the value at the present point
                            voxval = vol%get([nint(inp_coos(1)+orig(1)),&
                            nint(inp_coos(2)+orig(2)),nint(inp_coos(3)+orig(3))])
                            ! check whether the we are at a 0. voxel 
                            ! (this introduces slight artefacts near boundary of mask)
                            if( voxval .ne. 0. )then
                                win = recwin_3d(inp_coos(1), inp_coos(2), inp_coos(3), harwin_here)
                                rcomp_sum = 0.
                                do s=win(1,1),win(1,2)
                                    sosh = s+orig(1)
                                    do t=win(2,1),win(2,2)
                                        tosh = t+orig(2)
                                        do u=win(3,1),win(3,2)
                                            uosh = u+orig(3)
                                            rcomp = vol%get([sosh,tosh,uosh])
                                            w = wfun_here(hyp(inp_coos(1)-real(s),inp_coos(2)-&
                                            real(t),inp_coos(3)-real(u)))
                                            rcomp_sum = rcomp_sum+rcomp*w
                                        end do
                                    end do
                                end do
                                call img%add(rcomp_sum**2.,i,j,1) 
                            endif                    
                        endif
                    enddo
                endif
            enddo
        enddo
        !$omp end do nowait
        !$omp end parallel
        if( didft )then
            call vol%fwd_ft
        endif
        
        contains
        
            !>  \brief  linear interpolation with window halfwidth=1
            function wfun_here( x ) result( r )
                real, intent(in) :: x
                real ::  ax, r
                ax = abs(x)
                if( ax .gt. 1. )then
                    r = 0.
                else
                    r = 1.-ax/harwin_here
                endif
            end function
            
    end subroutine rproject
    
    ! HIGH-LEVEL PROJECTORS
    
    !>  \brief  generates an array of projection images of volume vol in orientations o
    function projvol( self, vol, o, p, top ) result( imgs )
        use simple_image,    only: image
        use simple_oris,     only: oris
        use simple_params,   only: params
        use simple_gridding, only: prep4cgrid
        class(projector),  intent(inout) :: self    !< projector instance
        class(image),      intent(inout) :: vol     !< volume to project
        class(oris),       intent(inout) :: o       !< orientations
        class(params),     intent(in)    :: p       !< parameters
        integer, optional, intent(in)    :: top     !< stop index
        type(image), allocatable         :: imgs(:) !< resulting images
        type(image)                      :: vol_pad, img_pad
        integer                          :: n, i, alloc_stat
        call vol_pad%new([p%boxpd,p%boxpd,p%boxpd], p%smpd, self%imgkind)
        if( self%imgkind .eq. 'xfel' )then
            call vol%pad(vol_pad)
        else
            call prep4cgrid(vol, vol_pad, p%msk, wfuns=self%wfuns)
        endif
        call img_pad%new([p%boxpd,p%boxpd,1], p%smpd, self%imgkind)
        if( present(top) )then
            n = top
        else
            n = o%get_noris()
        endif
        allocate( imgs(n), stat=alloc_stat )
        call alloc_err('projvol; simple_projector', alloc_stat)
        write(*,'(A)') '>>> GENERATES PROJECTIONS' 
        do i=1,n
            call progress(i, n)
            call imgs(i)%new([p%box,p%box,1], p%smpd, self%imgkind)
            call self%fproject(vol_pad, o%get_ori(i), img_pad)
            if( self%imgkind .eq. 'xfel' )then
                ! no back transformation or normalisation, just clipping
                call img_pad%clip(imgs(i))
            else
                call img_pad%bwd_ft
                call img_pad%clip(imgs(i))
                call imgs(i)%norm
            endif
        end do
        call vol_pad%kill
        call img_pad%kill
    end function projvol
    
    !>  \brief  rotates a volume by Euler angle o using Fourier gridding
    function rotvol( self, vol, o, p ) result( rovol )
        use simple_image,    only: image
        use simple_ori,      only: ori
        use simple_params,   only: params
        use simple_gridding, only: prep4cgrid
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(projector), intent(inout) :: self !< projector instance
        class(image),     intent(inout) :: vol  !< volume to project
        class(ori),       intent(inout) :: o    !< orientation
        class(params),    intent(in)    :: p    !< parameters
        type(image)                     :: vol_pad, rovol_pad, rovol  
        integer                         :: h,k,l,lims(3,2)
        real                            :: loc(3)
         if( self%imgkind .eq. 'xfel' ) stop 'routine not intended for xfel-kind images; simple_projector::rotvol'
        call vol_pad%new([p%boxpd,p%boxpd,p%boxpd], p%smpd)
        rovol_pad = vol_pad
        call rovol_pad%set_ft(.true.)
        call rovol%new([p%box,p%box,p%box], p%smpd)
        call prep4cgrid(vol, vol_pad, p%msk, wfuns=self%wfuns)
        lims = vol_pad%loop_lims(2)
        write(*,'(A)') '>>> ROTATING VOLUME'
        !$omp parallel do default(shared) private(h,k,l,loc) schedule(auto)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    loc  = matmul(real([h,k,l]),o%get_mat())
                    call rovol_pad%set_fcomp([h,k,l], self%extr_gridfcomp(vol_pad, loc))
                end do 
            end do
        end do
        !$omp end parallel do
        call rovol_pad%bwd_ft
        call rovol_pad%clip(rovol)     
        call rovol%norm
        call vol_pad%kill
        call rovol_pad%kill
    end function rotvol
    
    !>  \brief  rotates an image by angle ang using Fourier gridding
    subroutine rotimg( self, img, ang, msk, roimg )
        use simple_image,    only: image
        use simple_params,   only: params
        use simple_gridding, only: prep4cgrid
        use simple_math,     only: rotmat2d
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(projector), intent(inout) :: self  !< projector instance
        class(image),     intent(inout) :: img   !< image to rotate
        real,             intent(in)    :: ang   !< angle of rotation
        real,             intent(in)    :: msk   !< mask radius (in pixels)
        class(image),     intent(out)   :: roimg !< rotated image
        type(image)                     :: img_pad, roimg_pad, img_copy
        integer                         :: h,k,lims(3,2),ldim(3),ldim_pd(3) 
        real                            :: loc(3),mat(2,2),smpd
        if( self%imgkind .eq. 'xfel' ) stop 'routine not intended for xfel-kind images; simple_projector::rotimg'
        ldim       = img%get_ldim()
        ldim_pd    = nint(self%alpha)*ldim
        ldim_pd(3) = 1
        smpd       = img%get_smpd()
        call roimg%new(ldim, smpd)
        call img_pad%new(ldim_pd, smpd)
        call roimg_pad%new(ldim_pd, smpd)
        roimg_pad = cmplx(0.,0.)
        img_copy  = img
        call img_copy%mask(msk, 'soft')
        call prep4cgrid(img, img_pad, msk, wfuns=self%wfuns)
        lims = img_pad%loop_lims(2)
        mat = rotmat2d(ang)
        !$omp parallel do default(shared) private(h,k,loc) schedule(auto)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)                
                loc(:2) = matmul(real([h,k]),mat)
                loc(3)  = 0.
                call roimg_pad%set_fcomp([h,k,0], self%extr_gridfcomp(img_pad, loc))
            end do
        end do
        !$omp end parallel do
        call roimg_pad%bwd_ft
        call roimg_pad%clip(roimg)
        call img_pad%kill
        call roimg_pad%kill
        call img_copy%kill
    end subroutine rotimg
    
    !>  \brief  generates an array of polar Fourier transforms of volume vol in orientations o
    subroutine fprojvol_polar( self, vol, oset, p, pimgs, s )
        use simple_image,    only: image
        use simple_oris,     only: oris
        use simple_params,   only: params
        use simple_gridding, only: prep4cgrid
        use simple_polarft,  only: polarft
        use simple_ori,      only: ori
        class(projector), intent(inout) :: self                      !< projector instance
        class(image),     intent(inout) :: vol                       !< volume to project
        class(oris),      intent(inout) :: oset                      !< orientations
        class(params),    intent(in)    :: p                         !< parameters
        class(polarft),   intent(inout) :: pimgs(p%nstates,p%nspace) !< resulting polar FTs
        integer,          intent(in)    :: s                         !< state
        type(image) :: vol4grid
        type(ori)   :: o
        integer     :: n, i, ldim(3)
        ldim = vol%get_ldim()
        call vol4grid%new(ldim, p%smpd, self%imgkind)
        if( self%imgkind .eq. 'xfel' )then
            call vol%pad(vol4grid)
        else
            call prep4cgrid(vol, vol4grid, p%msk, wfuns=self%wfuns)
        endif
        n = oset%get_noris()
        write(*,'(A)') '>>> GENERATES PROJECTIONS' 
        do i=1,n
            call progress(i, n)
            call pimgs(s,i)%new([p%kfromto(1),p%kfromto(2)],p%ring2,ptcl=.false.)
            o = oset%get_ori(i)
            call self%fproject_polar(vol4grid, o, pimgs(s,i))
        end do
        call vol4grid%kill
    end subroutine fprojvol_polar

end module simple_projector
