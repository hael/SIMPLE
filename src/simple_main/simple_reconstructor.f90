module simple_reconstructor
use simple_fftw3
use simple_defs
use simple_image,   only: image
use simple_winfuns, only: winfuns
use simple_ctf,     only: ctf
implicit none

public :: reconstructor
private

type, extends(image) :: reconstructor
    private
    type(winfuns)               :: wfuns                 !< window functions object
    type(c_ptr)                 :: kp                    !< c pointer for fftw allocation
    type(ctf)                   :: tfun                  !< CTF object
    real(kind=c_float), pointer :: rho(:,:,:)=>null()    !< sampling+CTF**2 density
    integer                     :: rhosz                 !< size of the sampling density matrix
    character(len=STDLEN)       :: wfun_str='kb'         !< window function, string descriptor
    real                        :: winsz=1.              !< window half-width
    real                        :: alpha=2.              !< oversampling ratio
    real                        :: dens_const=1.0        !< density estimation constant, old val: 1/nptcls
    integer                     :: lfny=0                !< Nyqvist Fourier index
    character(len=STDLEN)       :: ctfflag=''            !< ctf flag <yes|no|mul|flip>
    character(len=STDLEN)       :: ikind                 !< image kind (em/xfel)
    logical                     :: tfneg=.false.         !< invert contrast or not
    logical                     :: tfastig=.false.       !< astigmatic CTF or not
    logical                     :: rho_allocated=.false. !< existence of rho matrix
  contains
    ! CONSTRUCTOR
    procedure          :: alloc_rho
    ! SETTERS
    procedure          :: reset
    ! GETTER
    procedure          :: get_wfuns
    ! CHECKERS
    procedure          :: rho_contains_nans
    ! I/O
    procedure          :: write_rho
    procedure          :: read_rho
    ! INTERPOLATION
    procedure, private :: inout_fcomp
    procedure, private :: calc_tfun_vals
    procedure          :: inout_fplane
    procedure          :: sampl_dens_correct
    ! SUMMATION
    procedure          :: sum
    ! RECONSTRUCTION
    procedure          :: rec
    ! DESTRUCTOR
    procedure          :: dealloc_rho
end type reconstructor

real            :: dfx=0., dfy=0., angast=0.
logical         :: debug=.false.
real, parameter :: SHTHRESH=0.0001

contains

    ! CONSTRUCTOR

    !>  \brief  allocates the sampling density matrix
    subroutine alloc_rho( self, p )
        use simple_params, only: params
        use simple_math,   only: fdim
        class(reconstructor),         intent(inout) :: self  !< instance
        class(params),                intent(in)    :: p     !< parameters
        integer                       :: ld_here(3)          !< logical dimension here
        integer                       :: rho_shape(3)        !< shape of the rho matrix
        integer                       :: rho_lims(3,2)       !< bounds of the rho matrix (xfel-kind images)
        character(len=:), allocatable :: ikind_tmp
        if( .not. self%exists() ) stop 'construct image before allocating rho; alloc_rho; simple_reconstructor'
        ld_here = self%get_ldim()
        if( ld_here(3) < 2 ) stop 'reconstructor need to be 3D 4 now; alloc_rho; simple_reconstructor'
        self%dens_const = 1./real(p%nptcls)
        self%wfun_str   = p%wfun
        self%winsz      = p%winsz
        self%alpha      = p%alpha
        self%ctfflag    = p%ctf
        self%tfneg = .false.
        if( p%neg .eq. 'yes' ) self%tfneg = .true.
        self%tfastig = .false.
        if( p%tfplan%mode .eq. 'astig' ) self%tfastig = .true.
        self%wfuns = winfuns(self%wfun_str,self%winsz,self%alpha)
        ikind_tmp  = self%get_imgkind()
        self%ikind = ikind_tmp
        if( trim(self%ikind) .eq. 'xfel' )then
            rho_lims = self%get_cmat_lims()
            allocate(self%rho(rho_lims(1,1):rho_lims(1,2),rho_lims(2,1):rho_lims(2,2),rho_lims(3,1):rho_lims(3,2)))
        else
            ! Work out dimensions of the rho array
            rho_shape(1)   = fdim(ld_here(1))
            rho_shape(2:3) = ld_here(2:3)
            ! Letting FFTW do the allocation in C ensures that we will be using aligned memory
            self%kp = fftwf_alloc_real(int(product(rho_shape),c_size_t))
            ! Set up the rho array which will point at the allocated memory
            call c_f_pointer(self%kp,self%rho,rho_shape)
        endif
        ! Set the record size of stack entry
        inquire(iolength=self%rhosz) self%rho
        self%rho_allocated = .true.
        call self%reset
    end subroutine alloc_rho
    
    ! SETTERS
    
    !>  \brief  resets the reconstructor object before reconstruction
    subroutine reset( self )
        class(reconstructor), intent(inout) :: self
        self     = cmplx(0.,0.)
        self%rho = 0.
    end subroutine reset
    
    ! GETTERS
    
    !>  \brief  return the window functions used by reconstructor
    function get_wfuns( self ) result( wfs )
        class(reconstructor), intent(inout) :: self
        type(winfuns) :: wfs
        wfs = self%wfuns
    end function get_wfuns
    
    ! CHECKERS
    
    !>  \brief  is for checking the numerical soundness of rho
    logical function rho_contains_nans( self )
        use ieee_arithmetic
        class(reconstructor), intent(in) :: self
        integer :: i, j ,k
        rho_contains_nans = .false.
        do i=1,size(self%rho,1)
            do j=1,size(self%rho,2)
                do k=1,size(self%rho,3)
                    if( ieee_is_nan(self%rho(i,j,k)) )then
                        rho_contains_nans = .true.
                        return
                    endif
                end do
            end do
        end do
    end function rho_contains_nans
    
    ! I/O

    !>  \brief  is for writing the sampling density (rho)
    subroutine write_rho( self, kernam )
        use simple_filehandling, only: get_fileunit, fopen_err
        class(reconstructor), intent(in) :: self
        character(len=*),     intent(in) :: kernam
        integer                          :: filnum, ier
        filnum = get_fileunit( )
        open(unit=filnum, status='replace', action='write', file=kernam,&
        access='direct', form='unformatted', recl=self%rhosz, iostat=ier)
        call fopen_err('write_rho; simple_reconstructor', ier)
        write(filnum, rec=1) self%rho
        close(unit=filnum)
    end subroutine write_rho
    
    !>  \brief  is for reading the sampling density (rho)
    subroutine read_rho( self, kernam )
        use simple_filehandling, only: get_fileunit, fopen_err
        class(reconstructor), intent(inout) :: self
        character(len=*),     intent(in)    :: kernam
        integer                             :: filnum, ier
        filnum = get_fileunit( )
        open(unit=filnum, status='old', action='read', file=kernam,&
        access='direct', form='unformatted', recl=self%rhosz, iostat=ier)
        call fopen_err('read_rho; simple_reconstructor', ier)
        read(filnum, rec=1) self%rho
        close(unit=filnum)
    end subroutine read_rho
    
    ! INTERPOLATION
    
    ! !> \brief  inserts or uninserts a Fourier plane component to the Fourier volume
    subroutine inout_fcomp( self, h, k, e, inoutmode, comp, oshift, pwght )
        use simple_ori,    only: ori
        use simple_math,   only: recwin_3d, euclid, hyp, cyci_1d
        use simple_jiffys, only: alloc_err
        class(reconstructor), intent(inout) :: self      !< the objetc
        integer,              intent(in)    :: h, k      !< Fourier indices
        class(ori),           intent(inout) :: e         !< orientation
        logical,              intent(in)    :: inoutmode !< add = .true., subtract = .false.
        complex,              intent(in)    :: comp      !< input component, if not given only sampling density calculation
        complex,              intent(in)    :: oshift    !< origin shift
        real, optional,       intent(in)    :: pwght     !< external particle weight (affects both fplane and rho)
        integer              :: i,j,m,kwin(3,2),phys(3),alloc_stat,inds(3),nn(3),lims(3,2), sh
        real                 :: w,vec(3),loc(3),tval,tvalsq,x,dist,mindist
        real, allocatable    :: kw1(:),kw2(:),kw3(:)
        integer, allocatable :: cyci1(:),cyci2(:),cyci3(:)
        complex              :: mod_comp
        if( comp == cmplx(0.,0.) ) return
        lims = self%loop_lims(3)
        ! calculate nonuniform sampling location
        vec(1) = real(h)
        vec(2) = real(k)
        vec(3) = 0.
        loc  = matmul(vec,e%get_mat())
        ! evaluate the transfer function
        call self%calc_tfun_vals(vec, tval, tvalsq)
        ! calculate kernel values
        kwin = recwin_3d(loc(1), loc(2), loc(3), self%winsz)
        allocate( kw1(kwin(1,1):kwin(1,2)), kw2(kwin(2,1):kwin(2,2)), kw3(kwin(3,1):kwin(3,2)), &
            & cyci1(kwin(1,1):kwin(1,2)), cyci2(kwin(2,1):kwin(2,2)), cyci3(kwin(3,1):kwin(3,2)), &
            & stat=alloc_stat )
        call alloc_err("In: inout_fcomp; simple_reconstructor", alloc_stat)
        do i=kwin(1,1),kwin(1,2)
            kw1(i) = self%wfuns%eval_apod(real(i)-loc(1))
            cyci1(i) = cyci_1d( lims(1,:),i )
        end do
        do j=kwin(2,1),kwin(2,2)
            kw2(j) = self%wfuns%eval_apod(real(j)-loc(2))
            cyci2(j) = cyci_1d( lims(2,:),j )
        end do
        do m=kwin(3,1),kwin(3,2)
            kw3(m) = self%wfuns%eval_apod(real(m)-loc(3))
            cyci3(m) = cyci_1d( lims(3,:),m )
        end do
        mindist = huge(x)   
        ! convolution interpolation
        do i=kwin(1,1),kwin(1,2)
            if( kw1(i) == 0. ) cycle
            inds(1) = cyci1(i)
            do j=kwin(2,1),kwin(2,2)
                if( kw2(j) == 0. ) cycle
                inds(2) = cyci2(j)
                do m=kwin(3,1),kwin(3,2)
                    if( kw3(m) == 0. ) cycle
                    ! calculate cyclic indices
                    inds(3) = cyci3(m)
                    ! calculate physical indices
                    phys = self%comp_addr_phys(inds)
                    ! calculate weight
                    w = kw1(i)*kw2(j)*kw3(m)*self%dens_const
                    ! keep track of the nearest neighbor
                    dist = euclid(loc,real([i,j,m]))
                    if( dist < mindist )then
                        mindist = dist
                        nn = [i,j,m]
                    endif
                    if( present(pwght) ) w = w*pwght
                    mod_comp = (comp*tval*w)*oshift ! CTF and w modulates the component before origin shift
                    if( inoutmode )then ! add
                        call self%add(inds, mod_comp, phys_in=phys) 
                        ! CTF**2 modulates the sampling density
                        self%rho(phys(1),phys(2),phys(3)) = self%rho(phys(1),phys(2),phys(3))+tvalsq*w
                    else ! subtract
                        call self%subtr(inds, mod_comp, phys_in=phys)
                        ! CTF**2 modulates the sampling density
                        self%rho(phys(1),phys(2),phys(3)) = self%rho(phys(1),phys(2),phys(3))-tvalsq*w 
                    endif
                end do
            end do
        end do
        deallocate(kw1,kw2,kw3,cyci1,cyci2,cyci3)   
    end subroutine inout_fcomp
   
    !> \brief  for evaluating the transfer function
    subroutine calc_tfun_vals( self, vec, tval, tvalsq )
        use simple_ori,  only: ori
        use simple_math, only: hyp
        use simple_ctf,  only: ctf
        class(reconstructor), intent(inout) :: self         !< instance
        real,                 intent(in)    :: vec(3)       !< nonuniform sampling location
        real,                 intent(out)   :: tval, tvalsq !< CTF and CTF**2.
        real    :: sqSpatFreq,ang,inv1,inv2
        integer :: ldim(3)
        if( self%ctfflag .ne. 'no' )then
            ldim       = self%get_ldim()
            inv1       = vec(1)*(1./real(ldim(1)))
            inv2       = vec(2)*(1./real(ldim(2)))
            sqSpatFreq = inv1*inv1+inv2*inv2
            ang        = atan2(vec(2), vec(1))
            ! calculate CTF and CTF**2 values 
            ! dfx, dfy, angast are class vars that are set by inout fplane
            tval   = self%tfun%eval(sqSpatFreq, dfx, dfy, angast, ang) ! no bfactor 4 now
            tvalsq = min(1.,max(tval**2.,0.001))
            if( self%ctfflag .eq. 'flip' ) tval = abs(tval)
            if( self%ctfflag .eq. 'mul'  ) tval = 1.
            if( self%tfneg ) tval = -tval
        else
            tval   = 1.
            tvalsq = tval
            if( self%tfneg ) tval = -tval
        endif
    end subroutine calc_tfun_vals
    
    !> \brief  for gridding or ungridding a Fourier plane
    subroutine inout_fplane( self, o, inoutmode, fpl, pwght, mul, shellweights )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_math, only: deg2rad, hyp
        use simple_ori,  only: ori
        class(reconstructor), intent(inout) :: self      !< instance
        class(ori),           intent(inout) :: o         !< orientation
        logical,              intent(in)    :: inoutmode !< add = .true., subtract = .false.
        class(image),         intent(inout) :: fpl       !< Fourier plane
        real, optional,       intent(in)    :: pwght     !< external particle weight (affects both fplane and rho)
        real, optional,       intent(in)    :: mul
        real, optional,       intent(in)    :: shellweights(:)
        integer                             :: h, k, lims(3,2), sh, lfny
        complex                             :: oshift=cmplx(1.,0.)
        real                                :: x=0., y=0., xtmp, ytmp, pw, shw
        logical                             :: pwght_present
        if( .not. fpl%is_ft() )       stop 'image need to be FTed; inout_fplane; simple_reconstructor'
        if( .not. (self.eqsmpd.fpl) ) stop 'scaling not yet implemented; inout_fplane; simple_reconstructor'
        pwght_present = present(pwght)
        if( self%ctfflag .ne. 'no' )then ! make CTF object & get CTF info
            self%tfun = ctf(self%get_smpd(), o%get('kv'), o%get('cs'), o%get('fraca'))
            dfx  = o%get('dfx')
            if( self%tfastig )then ! astigmatic CTF model
                dfy = o%get('dfy')
                angast = o%get('angast')
            else ! non-astigmatic CTF model
                dfy = dfx
                angast = 0.
            endif
        endif
        lims = self%loop_lims(2)
        x    = o%get('x')
        y    = o%get('y')
        xtmp = 0.
        ytmp = 0.
        if( abs(x) > SHTHRESH .or. abs(y) > SHTHRESH )then ! shift the image prior to insertion
            if( present(mul) )then
                x = x*mul
                y = y*mul
            endif
            xtmp = x
            ytmp = y
            if( abs(x) < 1e-6 ) xtmp = 0.
            if( abs(y) < 1e-6 ) ytmp = 0.
        endif
        if( present(shellweights) )then
            lfny = size(shellweights)
            !$omp parallel do default(shared) private(h,k,oshift,sh,pw) schedule(auto)
            do h=lims(1,1),lims(1,2)
                do k=lims(1,1),lims(1,2)
                    sh     = min(max(1,nint(hyp(real(h),real(k)))),lfny)
                    oshift = fpl%oshift([h,k,0], [-xtmp,-ytmp,0.], ldim=2)
                    if( pwght_present )then
                        pw = pwght*shellweights(sh)
                    else
                        pw = shellweights(sh)
                    endif
                    call self%inout_fcomp(h,k,o,inoutmode,fpl%get_fcomp([h,k,0]),oshift,pw)
                end do
            end do
            !$omp end parallel do
        else
            !$omp parallel do default(shared) private(h,k,oshift) schedule(auto)
            do h=lims(1,1),lims(1,2)
                do k=lims(1,1),lims(1,2)
                    oshift = fpl%oshift([h,k,0], [-xtmp,-ytmp,0.], ldim=2)
                    call self%inout_fcomp(h,k,o,inoutmode,fpl%get_fcomp([h,k,0]),oshift,pwght)
                end do
            end do
            !$omp end parallel do
        endif
    end subroutine inout_fplane
    
    !> \brief  for sampling density compensation & Wiener normalization
    subroutine sampl_dens_correct( self, self_out )
        class(reconstructor),   intent(inout) :: self
        class(image), optional, intent(inout) :: self_out
        integer                               :: h, k, l, lims(3,2), phys(3)
        ! set constants
        lims = self%loop_lims(2)
        if( present(self_out) ) self_out = self
        !$omp parallel do default(shared) private(h,k,l,phys) schedule(auto)
        do h=lims(1,1),lims(1,2)
            do k=lims(2,1),lims(2,2)
                do l=lims(3,1),lims(3,2)
                    phys = self%comp_addr_phys([h,k,l])
                    if( present(self_out) )then
                        call self_out%div([h,k,l],self%rho(phys(1),phys(2),phys(3)),phys_in=phys)
                    else
                        call self%div([h,k,l],self%rho(phys(1),phys(2),phys(3)),phys_in=phys) 
                    endif                    
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine sampl_dens_correct
    
    ! SUMMATION
    
    !> \brief  for summing reconstructors generated by parallel execution
    subroutine sum( self, self_in )
         class(reconstructor), intent(inout) :: self
         class(reconstructor), intent(in)    :: self_in
         call self%add(self_in)
         !$omp parallel workshare
         self%rho = self%rho+self_in%rho
         !$omp end parallel workshare
    end subroutine sum
    
    ! RECONSTRUCTION
    
    !> \brief  for reconstructing Fourier volumes according to the orientations 
    !!         and states in o, assumes that stack is open   
    subroutine rec( self, fname, p, o, se, state, mul, eo, part, wmat )
        use simple_oris,     only: oris
        use simple_sym,      only: sym
        use simple_params,   only: params
        use simple_ran_tabu, only: ran_tabu
        use simple_filterer, only: resample_filter
        use simple_jiffys,   only: find_ldim_nptcls, progress
        use simple_gridding  ! use all in there
        class(reconstructor), intent(inout) :: self      !< object
        character(len=*),     intent(inout) :: fname     !< spider/MRC stack filename
        class(params),        intent(in)    :: p         !< parameters
        class(oris),          intent(inout) :: o         !< orientations
        class(sym),           intent(inout) :: se        !< symmetry element
        integer,              intent(in)    :: state     !< state to reconstruct
        real,    optional,    intent(in)    :: mul       !< shift multiplication factor
        integer, optional,    intent(in)    :: eo        !< even(2) or odd(1)
        integer, optional,    intent(in)    :: part      !< partition (4 parallel rec)
        real,    optional,    intent(in)    :: wmat(:,:) !< shellweights
        type(image)       :: img, img_pd
        type(ran_tabu)    :: rt
        integer           :: statecnt(p%nstates), i, cnt, n, ldim(3), filtsz, filtsz_pad
        real, allocatable :: res(:), res_pad(:), wresamp(:)
        logical           :: doshellweight
        call find_ldim_nptcls(fname, ldim, n)
        if( n /= o%get_noris() ) stop 'inconsistent nr entries; rec; simple_reconstructor'
        doshellweight = present(wmat)
        ! make random number generator
        rt = ran_tabu(n)
        ! make the images
        call img_pd%new([p%boxpd,p%boxpd,1],self%get_smpd(),p%imgkind)
        call img%new([p%box,p%box,1],self%get_smpd(),p%imgkind)
        ! make the stuff needed for shellweighting
        res        = img%get_res()
        res_pad    = img_pd%get_res()
        filtsz_pad = img_pd%get_filtsz()
        ! calculate particle weights
        if( p%frac < 0.99 ) call o%calc_hard_ptcl_weights(p%frac, bystate=.true.)
        ! zero the Fourier volume and rho
        call self%reset
        write(*,'(A)') '>>> KAISER-BESSEL INTERPOLATION'
        statecnt = 0
        cnt = 0
        do i=1,p%nptcls
            call progress(i, p%nptcls) 
            if( i <= p%top .and. i >= p%fromp )then
                cnt = cnt+1
                if( nint(o%get(i,'state')) == state )then
                    statecnt(state) = statecnt(state)+1
                    if( present(eo) )then
                        if( mod(cnt,2) == 0 .and. eo == 2 )then
                            call rec_dens
                        else if( mod(cnt,2) /= 0 .and. eo == 1 )then
                            call rec_dens
                        endif
                    else
                        call rec_dens
                    endif
                endif
            endif
        end do
        if( present(part) )then
            return
        else
            write(*,'(A)') '>>> SAMPLING DENSITY (RHO) CORRECTION (JACKSON) & WIENER NORMALIZATION'
            call self%sampl_dens_correct
        endif
        if( p%l_xfel )then
            ! no bwd_ft or normalisation
        else
            call self%bwd_ft
            ! HAD TO TAKE OUT BECAUSE PGI COMPILER BAILS
            ! call self%norm
        endif
        call img%kill
        call img_pd%kill
        call rt%kill
        ! report how many particles were used to reconstruct each state
        if( p%nstates > 1 )then
            write(*,'(a,1x,i3,1x,a,1x,i6)') '>>> NR OF PARTICLES INCLUDED IN STATE:', state, 'WAS:', statecnt(state)
        endif
        
        contains
        
            !> \brief  the densty reconstruction functionality
            subroutine rec_dens
                use simple_ori, only: ori
                type(ori) :: orientation, o_sym
                integer   :: j, state
                real      :: pw
                state = nint(o%get(i, 'state'))
                if( state == 0 )return
                pw = 1.
                if( p%frac < 0.99 ) pw = o%get(i, 'w')
                if( pw > 0. )then
                    orientation = o%get_ori(i)
                    if( p%l_distr_exec )then
                        call img%read(p%stk_part, cnt, p%l_xfel)
                    else
                        call img%read(fname, i, p%l_xfel)
                    endif
                    if( p%l_xfel )then
                        call img%pad(img_pd)
                    else
                        call prep4cgrid(img, img_pd, p%msk, wfuns=self%wfuns)
                    endif
                    if( doshellweight )then
                        wresamp = resample_filter(wmat(i,:), res, res_pad)
                    else
                        allocate(wresamp(filtsz_pad))
                        wresamp = 1.0
                    endif
                    if( p%pgrp == 'c1' )then
                        call self%inout_fplane(orientation, .true., img_pd, pwght=pw, mul=mul, shellweights=wresamp)
                    else
                        do j=1,se%get_nsym()
                            o_sym = se%apply(orientation, j)
                            call self%inout_fplane(o_sym, .true., img_pd, pwght=pw, mul=mul, shellweights=wresamp)
                        end do
                    endif
                    deallocate(wresamp)
                endif
            end subroutine rec_dens
            
    end subroutine rec
    
    !>  \brief  is a destructor
    subroutine dealloc_rho( self )
        class(reconstructor), intent(inout) :: self  
        if( self%rho_allocated )then
            call fftwf_free(self%kp)
            self%rho  => null()
        endif
    end subroutine dealloc_rho
    
end module simple_reconstructor