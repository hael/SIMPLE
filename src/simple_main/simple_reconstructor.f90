module simple_reconstructor
use simple_fftw3
use simple_defs
use simple_image,   only: image
use simple_winfuns, only: winfuns
use simple_ctf,     only: ctf
use simple_jiffys,  only: alloc_err
use simple_kbinterpol ! use all in there
implicit none

public :: reconstructor
private

type, extends(image) :: reconstructor
    private
    type(winfuns)               :: wfuns                        !< window functions object
    type(c_ptr)                 :: kp                           !< c pointer for fftw allocation
    type(ctf)                   :: tfun                         !< CTF object
    real(kind=c_float), pointer :: rho(:,:,:)=>null()           !< sampling+CTF**2 density
    complex,        allocatable :: cmat_exp(:,:,:)              !< Fourier components of expanded reconstructor
    real,           allocatable :: rho_exp(:,:,:)               !< sampling+CTF**2 density of expanded reconstructor
    character(len=STDLEN)       :: wfun_str      = 'kb'         !< window function, string descriptor
    real                        :: winsz         = 1.           !< window half-width
    real                        :: alpha         = 2.           !< oversampling ratio
    real                        :: dens_const    = 1.           !< density estimation constant, old val: 1/nptcls
    integer                     :: lfny          = 0            !< Nyqvist Fourier index
    integer                     :: ldim_exp(3,2) = 0            !< dimensions of the Fourier components expanded matrix
    character(len=STDLEN)       :: ctfflag       = ''           !< ctf flag <yes|no|mul|flip>
    character(len=STDLEN)       :: ikind                        !< image kind (em/xfel)
    logical                     :: tfneg              = .false. !< invert contrast or not
    logical                     :: tfastig            = .false. !< astigmatic CTF or not
    logical                     :: rho_allocated      = .false. !< existence of rho matrix
    logical                     :: rho_exp_allocated  = .false. !< existence of rho expanded matrix
    logical                     :: cmat_exp_allocated = .false. !< existence of Fourien components expanded matrix
  contains
    ! CONSTRUCTORS
    procedure          :: alloc_rho
    ! SETTERS
    procedure          :: reset
    ! GETTER
    procedure          :: get_wfuns
    ! I/O
    procedure          :: write_rho
    procedure          :: read_rho
    ! INTERPOLATION
    procedure, private :: inout_fcomp
    procedure, private :: calc_tfun_vals
    procedure          :: inout_fplane
    procedure          :: sampl_dens_correct
    procedure          :: compress_exp
    ! SUMMATION
    procedure          :: sum
    ! RECONSTRUCTION
    procedure          :: rec
    ! DESTRUCTORS
    procedure          :: dealloc_exp
    procedure          :: dealloc_rho
end type reconstructor

real            :: dfx=0., dfy=0., angast=0.
logical         :: debug=.false.
real, parameter :: SHTHRESH=0.0001

contains

    ! CONSTRUCTOR

    !>  \brief  allocates the sampling density matrix
    subroutine alloc_rho( self, p, expand )
        use simple_params, only: params
        use simple_math,   only: fdim
        class(reconstructor),         intent(inout) :: self  !< instance
        class(params),                intent(in)    :: p     !< parameters
        logical,            optional, intent(in)    :: expand
        integer :: ld_here(3), rho_shape(3), rho_lims(3,2), lims(3,2), alloc_stat, dim
        character(len=:), allocatable :: ikind_tmp
        logical :: l_expand = .true.
        if(.not. self%exists() ) stop 'construct image before allocating rho; alloc_rho; simple_reconstructor'
        if(      self%is_2d()  )stop 'only for volumes; alloc_rho; simple_reconstructor'
        if( present(expand) )l_expand = expand
        ld_here = self%get_ldim()
        if( ld_here(3) < 2 ) stop 'reconstructor need to be 3D 4 now; alloc_rho; simple_reconstructor'
        call self%dealloc_rho
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
        self%rho_allocated = .true.
        if( l_expand )then
            ! setup expanded matrices
            lims = self%loop_lims(2)
            dim  = ceiling(sqrt(2.)*maxval(abs(lims))) + ceiling(self%winsz)
            self%ldim_exp(1,:) = [-dim, dim]
            self%ldim_exp(2,:) = [-dim, dim]
            self%ldim_exp(3,:) = [-dim, dim]
            allocate(self%cmat_exp( self%ldim_exp(1,1):self%ldim_exp(1,2),&
                                    &self%ldim_exp(2,1):self%ldim_exp(2,2),&
                                    &self%ldim_exp(3,1):self%ldim_exp(3,2)), stat=alloc_stat)
            call alloc_err("In: alloc_rho; simple_reconstructor 1", alloc_stat)
            allocate(self%rho_exp( self%ldim_exp(1,1):self%ldim_exp(1,2),&
                                    &self%ldim_exp(2,1):self%ldim_exp(2,2),&
                                    &self%ldim_exp(3,1):self%ldim_exp(3,2)), stat=alloc_stat)
            call alloc_err("In: alloc_rho; simple_reconstructor 2", alloc_stat)
            self%cmat_exp           = cmplx(0.,0.)
            self%rho_exp            = 0.
            self%cmat_exp_allocated = .true.
            self%rho_exp_allocated  = .true.
        endif
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
    
    ! I/O

    !>  \brief  is for writing the sampling density (rho)
    subroutine write_rho( self, kernam )
        use simple_filehandling, only: get_fileunit, fopen_err, del_file
        class(reconstructor), intent(in) :: self
        character(len=*),     intent(in) :: kernam
        character(len=100) :: io_message
        integer :: filnum, ier, io_stat
        call del_file(trim(kernam))
        filnum = get_fileunit( )
        open(unit=filnum, status='NEW', action='WRITE', file=trim(kernam), access='STREAM', iostat=ier)
        call fopen_err('write_rho; simple_reconstructor', ier)
        write(filnum, pos=1, iostat=io_stat, iomsg=io_message) self%rho    
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(write_rho): I/O error ', io_stat, ' when writing to: ', trim(kernam)
            write(*,'(2a)') 'IO error message was: ', io_message
            stop 'I/O error; write_rho; simple_reconstructor'
        endif
        close(unit=filnum)
    end subroutine write_rho
    
    !>  \brief  is for reading the sampling density (rho)
    subroutine read_rho( self, kernam )
        use simple_filehandling, only: get_fileunit, fopen_err
        class(reconstructor), intent(inout) :: self
        character(len=*),     intent(in)    :: kernam
        character(len=100) :: io_message
        integer :: filnum, ier, io_stat
        filnum = get_fileunit( )
        open(unit=filnum, status='OLD', action='READ', file=kernam, access='STREAM', iostat=ier)
        call fopen_err('read_rho; simple_reconstructor', ier)
        read(filnum, pos=1, iostat=io_stat, iomsg=io_message) self%rho
        if( io_stat .ne. 0 )then
            write(*,'(a,i0,2a)') '**ERROR(read_rho): I/O error ', io_stat, ' when reading from: ', trim(kernam)
            write(*,'(2a)') 'IO error message was: ', io_message
            stop 'I/O error; read_rho; simple_reconstructor'
        endif
        close(unit=filnum)
    end subroutine read_rho
    
    ! INTERPOLATION

    !> \brief  inserts or uninserts a Fourier plane component to the expanded Fourier matrix
    subroutine inout_fcomp( self, h, k, e, inoutmode, comp, oshift, pwght)
        use simple_ori,    only: ori
        use simple_math,   only: sqwin_3d, euclid, hyp, cyci_1d
        use simple_jiffys, only: alloc_err
        class(reconstructor), intent(inout) :: self      !< the objetc
        integer,              intent(in)    :: h, k      !< Fourier indices
        class(ori),           intent(inout) :: e         !< orientation
        logical,              intent(in)    :: inoutmode !< add = .true., subtract = .false.
        complex,              intent(in)    :: comp      !< input component, if not given only sampling density calculation
        complex,              intent(in)    :: oshift    !< origin shift
        real, optional,       intent(in)    :: pwght     !< external particle weight (affects both fplane and rho)
        real, allocatable :: w(:,:,:)
        real              :: vec(3), loc(3), tval, tvalsq
        integer           :: i, win(3,2), wdim
        if(comp == cmplx(0.,0.)) return
        ! calculate non-uniform sampling location
        vec = [real(h), real(k), 0.]
        loc = matmul(vec, e%get_mat())
        ! evaluate the transfer function
        call self%calc_tfun_vals(vec, tval, tvalsq)
        ! initiate kernel matrix
        win  = sqwin_3d(loc(1), loc(2), loc(3), self%winsz)
        wdim = 2*ceiling(self%winsz) + 1
        allocate(w(wdim, wdim, wdim))
        ! (weighted) kernel values
        if(present(pwght))then
            w = self%dens_const*pwght
        else
            w = self%dens_const
        endif
        do i=1,wdim
            where(w(i,:,:)>0.)w(i,:,:) = w(i,:,:) * self%wfuns%eval_apod(real(win(1,1)+i-1)-loc(1))
            where(w(:,i,:)>0.)w(:,i,:) = w(:,i,:) * self%wfuns%eval_apod(real(win(2,1)+i-1)-loc(2))
            where(w(:,:,i)>0.)w(:,:,i) = w(:,:,i) * self%wfuns%eval_apod(real(win(3,1)+i-1)-loc(3))
        enddo
        ! expanded matrices update
        if( inoutmode )then
            ! addition
            ! CTF and w modulates the component before origin shift
            self%cmat_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) =&
            &self%cmat_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) + (comp*tval*w)*oshift
            self%rho_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) =&
            &self%rho_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) + tvalsq*w
        else
            ! substraction
            ! CTF and w modulates the component before origin shift
            self%cmat_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) =&
            &self%cmat_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) - (comp*tval*w)*oshift
            self%rho_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) =&
            &self%rho_exp(win(1,1):win(1,2), win(2,1):win(2,2), win(3,1):win(3,2)) - tvalsq*w
        endif
        deallocate(w)
    end subroutine inout_fcomp
   
    !> \brief  for evaluating the transfer function
    subroutine calc_tfun_vals( self, vec, tval, tvalsq )
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
        real,    optional,    intent(in)    :: pwght     !< external particle weight (affects both fplane and rho)
        real,    optional,    intent(in)    :: mul
        real,    optional,    intent(in)    :: shellweights(:)
        integer :: h, k, lims(3,2), sh, lfny, logi(3), phys(3)
        complex :: oshift=cmplx(1.,0.)
        real    :: x=0., y=0., xtmp, ytmp, pw
        logical :: pwght_present
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
        if(present(shellweights))then
            lfny = size(shellweights)
            !$omp parallel do collapse(2) default(shared) schedule(auto)&
            !$omp private(h,k,oshift,sh,pw,logi,phys)
            do h=lims(1,1),lims(1,2)
                do k=lims(1,1),lims(1,2)
                    logi   = [h,k,0]
                    sh     = min(max(1,nint(hyp(real(h),real(k)))),lfny)
                    oshift = fpl%oshift(logi, [-xtmp,-ytmp,0.], ldim=2)
                    if( pwght_present )then
                        pw = pwght*shellweights(sh)
                    else
                        pw = shellweights(sh)
                    endif
                    phys = fpl%comp_addr_phys(logi)
                    call self%inout_fcomp(h,k,o,inoutmode,fpl%get_fcomp(logi,phys),oshift,pw)
                end do
            end do
            !$omp end parallel do
        else
            !$omp parallel do collapse(2) default(shared) private(h,k,oshift,logi,phys)&
            !$omp schedule(auto)
            do h=lims(1,1),lims(1,2)
                do k=lims(1,1),lims(1,2)
                    logi   = [h,k,0]
                    oshift = fpl%oshift(logi, [-xtmp,-ytmp,0.], ldim=2)
                    phys   = fpl%comp_addr_phys(logi)
                    call self%inout_fcomp(h,k,o,inoutmode,fpl%get_fcomp(logi,phys),oshift,pwght)
                end do
            end do
            !$omp end parallel do
        endif
    end subroutine inout_fplane
    
    !> \brief  for sampling density compensation & Wiener normalization
    subroutine sampl_dens_correct( self, self_out )
        class(reconstructor),   intent(inout) :: self
        class(image), optional, intent(inout) :: self_out
        integer :: h, k, l, lims(3,2), phys(3)
        ! set constants
        lims = self%loop_lims(2)
        if( present(self_out) ) self_out = self
        !$omp parallel do collapse(3) default(shared) private(h,k,l,phys) schedule(auto)
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

    !>  \brief converts the expanded matrix to standard imaginary representation
    subroutine compress_exp( self )
        class(reconstructor), intent(inout) :: self
        complex :: comp, zero
        integer :: lims(3,2), logi(3), phys(3), h, k, m
        if(.not. self%cmat_exp_allocated .or. .not.self%rho_allocated)then
            stop 'expanded complex or rho matrices do not exist; simple_reconstructor::compress_exp'
        endif
        call self%reset
        lims = self%loop_lims(3)
        zero = cmplx(0.,0.)
        ! Fourier components & rho matrices compression
        ! Can't be threaded because of add_fcomp
        do h = lims(1,1),lims(1,2)
            do k = lims(2,1),lims(2,2)
                if(any(self%cmat_exp(h,k,:).ne.zero))then
                    do m = lims(3,1),lims(3,2)
                        comp = self%cmat_exp(h,k,m)
                        if(comp .eq. zero)cycle
                        logi = [h, k, m]
                        phys = self%comp_addr_phys(logi)
                        ! addition because FC and its Friedel mate must be summed
                        ! as the expansion updates one or the other
                        call self%add_fcomp(logi, phys, comp)
                        self%rho(phys(1),phys(2),phys(3)) = &
                            &self%rho(phys(1),phys(2),phys(3)) + self%rho_exp(h,k,m)
                    end do
                endif
            end do 
        end do
    end subroutine compress_exp
    
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
        type(image) :: img, img_pd
        integer     :: statecnt(p%nstates), i, cnt, n, ldim(3)
        integer     :: state_here, state_glob
        logical     :: doshellweight
        call find_ldim_nptcls(fname, ldim, n)
        if( n /= o%get_noris() ) stop 'inconsistent nr entries; rec; simple_reconstructor'
        doshellweight = present(wmat)
        ! stash global state index
        state_glob = state 
        ! make random number generator
        ! make the images
        call img_pd%new([p%boxpd,p%boxpd,1],self%get_smpd(),p%imgkind)
        call img%new([p%box,p%box,1],self%get_smpd(),p%imgkind)
        ! calculate particle weights
        call o%calc_spectral_weights(p%frac, bystate=.true.)
        ! zero the Fourier volume and rho
        call self%reset
        write(*,'(A)') '>>> KAISER-BESSEL INTERPOLATION'
        statecnt = 0
        cnt      = 0
        do i=1,p%nptcls
            call progress(i, p%nptcls) 
            if( i <= p%top .and. i >= p%fromp )then
                cnt = cnt+1
                state_here = nint(o%get(i,'state')) 
                if( state_here > 0 .and. (state_here == state) )then
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
            call self%compress_exp
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
                if( state == 0 ) return
                pw = 1.
                if( p%frac < 0.99 ) pw = o%get(i, 'w')
                if( pw > 0. )then
                    orientation = o%get_ori(i)
                    if( p%l_distr_exec )then
                        call img%read(p%stk_part, cnt, isxfel=p%l_xfel)
                    else
                        call img%read(fname, i, isxfel=p%l_xfel)
                    endif
                    if( p%l_xfel )then
                        call img%pad(img_pd)
                    else
                        call prep4cgrid(img, img_pd, p%msk)
                    endif
                    if( p%pgrp == 'c1' )then
                        if( doshellweight )then
                            call self%inout_fplane(orientation, .true., img_pd, pwght=pw, mul=mul,&
                                &shellweights=wmat(i,:))
                        else
                            call self%inout_fplane(orientation, .true., img_pd, pwght=pw, mul=mul)
                        endif
                    else
                        do j=1,se%get_nsym()
                            o_sym = se%apply(orientation, j)
                            if( doshellweight )then
                                call self%inout_fplane(o_sym, .true., img_pd, pwght=pw, mul=mul,&
                                    &shellweights=wmat(i,:))
                            else
                                call self%inout_fplane(o_sym, .true., img_pd, pwght=pw, mul=mul)
                            endif
                        end do
                    endif
                endif
            end subroutine rec_dens
    end subroutine rec

    ! DESTRUCTORS

    !>  \brief  is the expanded destructor
    subroutine dealloc_exp( self )
        class(reconstructor), intent(inout) :: self
        if(allocated(self%rho_exp)) deallocate(self%rho_exp)
        if(allocated(self%cmat_exp))deallocate(self%cmat_exp)
        self%rho_exp_allocated  = .false.
        self%cmat_exp_allocated = .false.
    end subroutine dealloc_exp

    !>  \brief  is a destructor
    subroutine dealloc_rho( self )
        class(reconstructor), intent(inout) :: self
        call self%dealloc_exp
        if( self%rho_allocated )then
            call fftwf_free(self%kp)
            self%rho => null()
            self%rho_allocated = .false.
        endif
    end subroutine dealloc_rho

end module simple_reconstructor