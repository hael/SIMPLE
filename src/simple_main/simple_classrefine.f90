module simple_classrefine
use simple_params,           only: params
use simple_build,            only: build
use simple_prime2D_srch,     only: prime2D_srch
use simple_ori,              only: ori
use simple_oris,             only: oris
use simple_image,            only: image
use simple_polarft_corrcalc, only: polarft_corrcalc
use simple_convergence,      only: convergence
use simple_prime_srch,       only: prime_srch
use simple_hadamard_common
use simple_filterer
use simple_jiffys
use simple_stat
use simple_defs
implicit none

type classrefine
    private
    class(build),  pointer     :: bp => null()      !< pointer to builder
    class(params), pointer     :: pp => null()      !< pointer to params
    type(image),   allocatable :: imgs(:)           !< experimental images
    integer,       allocatable :: pinds(:)          !< particle indices
    real,          allocatable :: wmat(:,:)         !< shellweights
    type(prime_srch)           :: srch_common       !< functionalities common to primesrch2D/3D
    type(oris)                 :: a                 !< alignment params
    type(oris)                 :: a_modified        !< alignment params, modified
    type(prime2D_srch)         :: srch2D            !< search object
    type(polarft_corrcalc)     :: pftcc_valid       !< PFTCC object (static lowres lim)
    type(polarft_corrcalc)     :: pftcc_dynres      !< PFTCC object (dynamic lowres lim)
    type(image)                :: avg               !< local cavg
    type(image)                :: ref               !< local reference
    type(image)                :: imgsum            !< Wiener nominator
    type(image)                :: ctfsqsum          !< Wiener denominator
    type(convergence)          :: conv              !< convergence object
    type(ctfplan)              :: tfplan            !< CTF plan     
    real                       :: smpd              !< sampling distance
    real                       :: msk               !< mask radius
    real                       :: lp_opt            !< optimal low-pass limit found
    real                       :: validlp           !< low-pass limit for validation
    integer                    :: startlpstage      !< start index for lp stepping
    integer                    :: class             !< class being refined
    integer                    :: ring2             !< rounded integer mask radius
    integer                    :: kfromto_valid(2)  !< static Fourier index range
    integer                    :: ldim(3)           !< logical image dimension
    integer                    :: nimgs  = 0        !< nr of imgs in class
    integer                    :: filtsz = 0        !< filter size
    logical                    :: exists = .false.  !< to indicate existence
  contains
    ! CONSTRUCTOR
    procedure          :: new
    ! MASTER REFINEMENT METHOD
    procedure          :: refine_master
    ! GETTERS
    procedure          :: update_parent_oritab
    procedure          :: get_avg
    procedure          :: get_nimgs
    ! I/O
    procedure          :: write_shellweights
    ! SLAVES
    procedure, private :: refine
    procedure, private :: validate
    procedure, private :: calc_wiener_sums
    procedure, private :: add_image_to_wiener_sums
    procedure, private :: subtr_image_from_wiener_sums
    procedure, private :: prep_ref
    procedure, private :: prep_avg
    procedure, private :: center_shifts
    procedure, private :: calc_shellweights
    ! DESTRUCTOR
    procedure          :: kill
end type classrefine

integer, parameter :: NDLPS=3
real, parameter    :: DYNLPS(NDLPS)=[20.0,10.0,6.0]
integer, parameter :: MAXITS=15
logical            :: DEBUG=.false.

contains

    ! CONSTRUCTOR
    
    !>  \brief  is a constructor
    subroutine new( self, b, p, cline, stk, class, msk, ctf )
        use simple_math,    only: round2even
        use simple_cmdline, only: cmdline
        class(classrefine),     intent(inout) :: self
        class(build),   target, intent(inout) :: b
        class(params),  target, intent(inout) :: p
        class(cmdline),         intent(in)    :: cline
        character(len=*),       intent(in)    :: stk
        integer,                intent(in)    :: class
        real,                   intent(in)    :: msk
        character(len=*),       intent(in)    :: ctf
        type(ori) :: o
        integer   :: iimg, alloc_stat
        ! set pointers
        self%bp => b
        self%pp => p
        ! retrieve inputs
        self%ring2 = round2even(msk)
        self%msk   = msk
        self%class = class
        ! make CTF plan
        self%tfplan%flag = ctf
        if( b%a%isthere('dfx') .and. b%a%isthere('dfy'))then
            self%tfplan%mode = 'astig'
        else if( b%a%isthere('dfx') )then
            self%tfplan%mode = 'noastig'
        else
            self%tfplan%mode = 'no'
        endif
        ! build images
        self%ldim = b%img%get_ldim()
        self%smpd = b%img%get_smpd()
        call self%avg%new(self%ldim, self%smpd)
        call self%ref%new(self%ldim, self%smpd)
        call self%imgsum%new(self%ldim, self%smpd)
        call self%ctfsqsum%new(self%ldim, self%smpd)
        ! read in images and orientations
        self%pinds = b%a%get_cls_pinds(class)
        self%nimgs = size(self%pinds)
        call self%a%new(self%nimgs)
        allocate(self%imgs(self%nimgs), stat=alloc_stat)
        call alloc_err("In: simple_classrefine :: new, 1", alloc_stat)
        do iimg=1,self%nimgs
            call b%img%read(stk, self%pinds(iimg))
            self%imgs(iimg) = b%img
            o = b%a%get_ori(self%pinds(iimg))
            call o%set('class', 1.0)
            call self%a%set_ori(iimg, o)
        end do
        ! allocate weight matrix
        self%filtsz = self%imgs(1)%get_filtsz()
        allocate( self%wmat(self%nimgs,self%filtsz), stat=alloc_stat)
        call alloc_err("In: simple_classrefine :: new, 2", alloc_stat)
        self%wmat  = 1.0
        ! create srch object
        call self%srch2D%new(self%pp, ncls=1)
        ! create convergence object
        self%conv = convergence(self%a_modified, self%pp, cline)
        ! create srch_common object
        self%srch_common = prime_srch(self%pp, 1, round2even(twopi*real(self%ring2)))
        ! indicate existence
        self%exists = .true.
        if( DEBUG ) print *, 'constructed new classrefine object; simple_classrefine :: new'
    end subroutine new

    ! MASTER REFINEMENT METHOD

    !>  \brief  does the multi-resolution refinement
    subroutine refine_master( self, lplim )
        use simple_math, only: find
        class(classrefine), intent(inout) :: self
        real,               intent(in)    :: lplim
        type(oris) :: a_stash(NDLPS)
        integer    :: loc(1), lpstage
        real       :: validcorrs(NDLPS), dist
        logical    :: available(NDLPS)
        ! set resolution range
        call find(DYNLPS, NDLPS, lplim, self%startlpstage, dist)
        self%validlp = DYNLPS(self%startlpstage)
        self%kfromto_valid(1) = max(2,self%bp%img%get_find(self%pp%hp))
        self%kfromto_valid(2) = self%bp%img%get_find(self%validlp)
        ! generate initial solution
        call self%center_shifts(self%a)
        call self%calc_wiener_sums(self%a)
        ! start the process
        available = .false.
        do lpstage=self%startlpstage,NDLPS
            if( DYNLPS(lpstage) >= 2.0*self%smpd )then
                available(lpstage) = .true.
                call self%refine(lpstage)
                if( DEBUG ) call self%avg%write('classrefine_avg_stage'//int2str(lpstage)//'.mrc', 1)
                validcorrs(lpstage) = self%validate(self%a_modified)
                if( DEBUG ) print *, 'validcorr: ', validcorrs(lpstage)
                if( lpstage ==  self%startlpstage )then
                    ! update self%a
                    self%a = self%a_modified
                else
                    ! update self%a if improvement
                    if( validcorrs(lpstage) >= validcorrs(lpstage-1) ) self%a = self%a_modified
                endif
                ! stash modified parameters for selection before return
                a_stash(lpstage) = self%a_modified
            else
                available(lpstage) = .false.
            endif
        end do
        loc = maxloc(validcorrs, available)
        self%lp_opt = DYNLPS(loc(1))
        ! set a to the best solution identified
        self%a = a_stash(loc(1))
        ! update the average
        call self%calc_wiener_sums(self%a)
        call self%prep_avg
        if( DEBUG ) print *, 'completed; simple_classrefine :: refine_master'
    end subroutine refine_master

    ! GETTERS

    !>  \brief  updates the oritab in the builder with the new alignment parameters
    subroutine update_parent_oritab( self )
        class(classrefine), intent(inout) :: self
        integer   :: iimg
        type(ori) :: o
        do iimg=1,self%nimgs
            o = self%a%get_ori(iimg)
            call o%set('class', real(self%class))
            call self%bp%a%set_ori(self%pinds(iimg), o)
            call self%bp%a%set(self%pinds(iimg), 'lp', self%lp_opt)
        end do
    end subroutine update_parent_oritab

    !>  \brief  for getting the final refined average
    function get_avg( self ) result( avgout )
        class(classrefine), intent(inout) :: self
        type(image) :: avgout
        avgout = self%avg
    end function get_avg

    !>  \brief  for getting the number of images in the class being refined
    integer function get_nimgs( self )
        class(classrefine), intent(in) :: self
        get_nimgs = self%nimgs
    end function get_nimgs

    ! I/O

    !>  \brief  writes the shellweights to file
    subroutine write_shellweights( self, fnam )
        class(classrefine), intent(in) :: self
        character(len=*),   intent(in) :: fnam
        call arr2D2file(self%wmat, fnam)
    end subroutine write_shellweights

    ! SLAVES

    !>  \brief  refines the average in the given resolution range (lpstage)
    !!          assumes that the Wiener sums have been assembled beforehand
    subroutine refine( self, lpstage )
        class(classrefine), intent(inout) :: self
        integer,            intent(in)    :: lpstage
        type(ori)  :: o
        type(oris) :: a_best
        integer    :: i, iimg, kfromto_dynamic(2)
        real       :: dfx, dfy, angast, corr, corr_best
        logical    :: converged
        self%a_modified = self%a
        a_best          = self%a
        ! make the dynamic correlator
        kfromto_dynamic(1) = self%kfromto_valid(1)
        kfromto_dynamic(2) = self%bp%img%get_find(DYNLPS(lpstage))
        call self%pftcc_dynres%new(1, [1,1], self%ldim,&
        kfromto_dynamic, self%ring2, self%tfplan%flag)
        write(*,'(a,1x,i1,1x,a,1x,i1)') '>>> REFINEMENT STAGE', lpstage, '/', NDLPS
        corr_best = -1.
        do i=1,MAXITS
            do iimg=1,self%nimgs
                ! prepare reference
                o = self%a_modified%get_ori(iimg)
                if( self%tfplan%mode .eq. 'astig' )then ! astigmatic CTF
                    dfx    = o%get('dfx')
                    dfy    = o%get('dfy')
                    angast = o%get('angast')
                else if( self%tfplan%mode .eq. 'noastig' )then
                    dfx    = o%get('dfx')
                    dfy    = dfx
                    angast = 0.
                else
                    stop 'Unsupported ctf mode; simple_classrefine :: refine'
                endif
                call self%subtr_image_from_wiener_sums(iimg, self%a_modified)
                call self%prep_ref(iimg)
                call self%bp%proj%img2polarft(1, self%ref, self%pftcc_dynres, isptcl=.false.)
                if( self%tfplan%flag .ne. 'no' )then
                    call self%pftcc_dynres%apply_ctf(self%smpd, self%bp%tfun, dfx, dfy, angast)
                endif
                ! prepare particle
                self%bp%img = self%imgs(iimg)
                call prepimg4align(self%bp, self%pp, o)
                call self%bp%proj%img2polarft(1, self%bp%img, self%pftcc_dynres)
                ! search
                o = self%a_modified%get_ori(iimg)
                call self%srch2D%prep4srch(self%pftcc_dynres, 1, DYNLPS(lpstage), o)
                call self%srch2D%stochastic_srch_shc(self%pftcc_dynres, 1)
                call self%srch2D%get_cls(o)
                ! update parameters 
                call self%a_modified%set_ori(iimg,o)
                ! add image back to Wiener sums
                call self%add_image_to_wiener_sums(iimg, self%a_modified)
            end do
            if( i < 3 )then ! minimum of three iterations
                converged = .false.
            else
                converged = self%conv%check_conv2D(ncls=1)
            endif
            corr = self%a_modified%median('corr')
            if( corr >= corr_best )then
                write(*,'(a,1x,i2,1x,a,1x,f7.3)') 'ITERATION:', i, 'CORRELATION:', corr
                call self%center_shifts(self%a_modified)
                a_best = self%a_modified
                corr_best = corr
            else
                self%a_modified = a_best
            endif
            call self%calc_wiener_sums(self%a_modified)
            call self%prep_avg
            call self%calc_shellweights(self%a_modified)
            if( converged ) exit
        end do
        if( DEBUG ) print *, 'completed; simple_classrefine :: refine'
    end subroutine refine

    !>  \brief  prepares the correlators for search
    function validate( self, os ) result( corrvalid )
        use simple_math, only: median_nocopy
        class(classrefine), intent(inout) :: self
        class(oris),        intent(inout) :: os
        integer     :: iimg, rot
        type(image) :: ref_local
        type(ori)   :: o
        real        :: corrs(self%nimgs), corrvalid, dfx, dfy, angast
        
        call self%pftcc_valid%new(1, [1,1], self%ldim,&
        self%kfromto_valid, self%ring2, self%tfplan%flag)
        ! do the reference preparation
        ref_local = self%avg
        call prep2Dref(self%pp, ref_local)
        ! transfer the reference to polar coordiates
        call self%bp%proj%img2polarft(1, ref_local, self%pftcc_valid, isptcl=.false.)
        do iimg=1,self%nimgs
            ! prepare reference
            o = os%get_ori(iimg)
            if( self%tfplan%mode .eq. 'astig' )then ! astigmatic CTF
                dfx    = o%get('dfx')
                dfy    = o%get('dfy')
                angast = o%get('angast')
            else if( self%tfplan%mode .eq. 'noastig' )then
                dfx    = o%get('dfx')
                dfy    = dfx
                angast = 0.
            else
                stop 'Unsupported ctf mode; simple_classrefine :: refine'
            endif
            if( self%tfplan%flag .ne. 'no' )then
                call self%pftcc_valid%apply_ctf(self%smpd, self%bp%tfun, dfx, dfy, angast)
            endif
            ! transfer the particle to polar coordiates
            self%bp%img = self%imgs(iimg)
            rot         = self%srch_common%roind(360.-o%e3get())  ! in-plane angle index
            call prepimg4align(self%bp, self%pp, o)
            call self%bp%proj%img2polarft(1, self%bp%img, self%pftcc_valid)
            ! calculate corr
            corrs(iimg) = self%pftcc_valid%corr(1, 1, rot)
        end do
        ! calculate a joint (median) correlation coefficient
        corrvalid = median_nocopy(corrs)
        call ref_local%kill
        if( DEBUG ) print *, 'completed; simple_classrefine :: validate'
    end function validate

    !>  \brief  assemblage of the Wiener sums
    subroutine calc_wiener_sums( self, os )
        class(classrefine), intent(inout) :: self
        class(oris),        intent(inout) :: os
        integer   :: iimg
        type(ori) :: o
        self%imgsum   = 0.
        self%ctfsqsum = cmplx(0.,0.)
        do iimg=1,self%nimgs
            call self%add_image_to_wiener_sums(iimg, os)
        end do
        if( DEBUG ) print *, 'calculated wiener sums; simple_classrefine :: calc_wiener_sums'
    end subroutine calc_wiener_sums

    !>  \brief  adds an image to the Wiener sums
    subroutine add_image_to_wiener_sums( self, iimg, os )
        class(classrefine), intent(inout) :: self
        integer,            intent(in)    :: iimg
        class(oris),        intent(inout) :: os
        type(ori) :: o
        o = os%get_ori(iimg)
        self%bp%img = self%imgs(iimg)
        call wiener_restore2D_online(self%bp%img, o, self%bp%tfun,&
        self%tfplan, self%imgsum, self%msk, self%wmat(iimg,:))
        call assemble_ctfsqsum_online(self%bp%img, o, self%bp%tfun,&
        self%tfplan, self%ctfsqsum, self%wmat(iimg,:))
    end subroutine add_image_to_wiener_sums

    !>  \brief  subtracts an image from the Wiener sums
    subroutine subtr_image_from_wiener_sums( self, iimg, os )
        class(classrefine), intent(inout) :: self
        integer,            intent(in)    :: iimg
        class(oris),        intent(inout) :: os
        type(ori) :: o
        o = os%get_ori(iimg)
        self%bp%img = self%imgs(iimg)
        call wiener_restore2D_online(self%bp%img, o, self%bp%tfun,&
        self%tfplan, self%imgsum, self%msk, self%wmat(iimg,:), add=.false.)
        call assemble_ctfsqsum_online(self%bp%img, o, self%bp%tfun,&
        self%tfplan, self%ctfsqsum, self%wmat(iimg,:), add=.false.)
    end subroutine subtr_image_from_wiener_sums

    !>  \brief  density correction and masking of the reference 2D template
    subroutine prep_ref( self, iimg )
        class(classrefine), intent(inout) :: self
        integer,            intent(in)    :: iimg
        type(ori) :: o
        ! do the density correction
        self%ref = self%imgsum
        call self%ref%fwd_ft
        call self%ref%ctf_dens_correct(self%ctfsqsum)
        call self%ref%bwd_ft
        ! do the reference preparation
        call prep2Dref(self%pp, self%ref)
    end subroutine prep_ref

    !>  \brief  density correction of the refined average
    subroutine prep_avg( self )
        class(classrefine), intent(inout) :: self
        ! do the density correction
        self%avg = self%imgsum
        call self%avg%fwd_ft
        call self%avg%ctf_dens_correct(self%ctfsqsum)
        call self%avg%bwd_ft
    end subroutine prep_avg

    !>  \brief  centers the shifts to avoid drifting
    subroutine center_shifts( self, os )
        class(classrefine), intent(inout) :: self
        class(oris),        intent(inout) :: os
        real, allocatable :: xshifts(:), yshifts(:)
        real              :: xavg, yavg
        xshifts = os%get_all('x')
        yshifts = os%get_all('y')
        xavg    = sum(xshifts)/real(self%nimgs)
        yavg    = sum(yshifts)/real(self%nimgs)
        xshifts = xshifts - xavg
        yshifts = yshifts - yavg
        call os%set_all('x', xshifts)
        call os%set_all('y', yshifts)
    end subroutine center_shifts

    !>  \brief  calculates the shellweights used in the Wiener restoration
    subroutine calc_shellweights( self, os )
        use simple_estimate_ssnr, only: ssnr2optlp
        class(classrefine), intent(inout) :: self
        class(oris),        intent(inout) :: os
        integer           :: iimg, ishell
        real, allocatable :: res(:), corrs(:), weights_tmp(:)
        type(image)       :: ref_local, ref_tmp, fdiff
        type(ori)         :: o
        real              :: dfx, dfy, angast, wsums(self%nimgs)
        integer           :: j
        ! do the reference preparation
        ref_tmp = self%avg
        call prep2Dref(self%pp, ref_tmp)
        ! calculate FRCs
        do iimg=1,self%nimgs
            ! prepare reference
            ref_local = ref_tmp
            o = os%get_ori(iimg)
            if( self%tfplan%mode .eq. 'astig' )then ! astigmatic CTF
                dfx    = o%get('dfx')
                dfy    = o%get('dfy')
                angast = o%get('angast')
            else if( self%tfplan%mode .eq. 'noastig' )then
                dfx    = o%get('dfx')
                dfy    = dfx
                angast = 0.
            else
                stop 'Unsupported ctf mode; simple_classrefine :: calc_shellweights'
            endif
            call self%bp%tfun%apply(ref_local, dfx, 'ctf', dfy, angast)
            ! prepare particle
            self%bp%img = self%imgs(iimg)
            call prepimg4align(self%bp, self%pp, o)
            ! calc FRC
            call ref_local%fsc(self%bp%img, res, corrs)
            self%wmat(iimg,:) = corrs
        end do
        ! create weights, shell-by-shell
        do ishell=1,self%filtsz
            weights_tmp = corrs2weights(self%wmat(:,ishell))*real(self%nimgs)
            self%wmat(:,ishell) = weights_tmp
            deallocate(weights_tmp)
        end do
        ! set wsums (integrated shellweights)
        do iimg=1,self%nimgs
            wsums(iimg) = sum(self%wmat(iimg,:))
        end do
        call normalize_sigm(wsums)
        ! modify the shellweights accordingly
        do iimg=1,self%nimgs
            self%wmat(iimg,:) = self%wmat(iimg,:)*wsums(iimg)
        end do
        call ref_local%kill
    end subroutine calc_shellweights

    ! DESTRUCTOR

    !>  \brief  is a destructor
    subroutine kill( self )
        class(classrefine), intent(inout) :: self
        integer :: iimg
        if( self%exists )then
            self%bp => null()
            self%pp => null()
            do iimg=1,self%nimgs
                call self%imgs(iimg)%kill
            end do
            deallocate(self%imgs, self%pinds, self%wmat)
            call self%srch_common%kill
            call self%a%kill
            call self%a_modified%kill
            call self%srch2D%kill
            call self%pftcc_valid%kill
            call self%pftcc_dynres%kill
            call self%avg%kill
            call self%ref%kill
            call self%imgsum%kill
            call self%ctfsqsum%kill
            self%exists = .false.
        endif
    end subroutine kill

end module simple_classrefine