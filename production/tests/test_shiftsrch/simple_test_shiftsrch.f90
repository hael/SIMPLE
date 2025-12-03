program simple_test_shiftsrch
! include 'simple_lib.f08'
! use simple_polarft_calc,  only: polarft_calc
! use simple_cmdline,           only: cmdline
! use simple_builder,           only: builder
! use simple_image,             only: image
! use simple_parameters,        only: parameters, params_glob
! use simple_polarizer,         only: polarizer
! use simple_pftc_shsrch_grad, only: pftc_shsrch_grad  ! gradient-based in-plane angle and shift search
! use simple_commanders_volops,  only: commander_reproject
! implicit none
! type(cmdline)                 :: cline, cline_projection
! type(builder)                 :: b
! type(parameters)              :: p
! type(polarft_calc)        :: pftc
! type(polarizer)               :: img_copy
! type(pftc_shsrch_grad)       :: grad_shsrch_obj           !< origin shift search object, L-BFGS with gradient
! type(commander_reproject)     :: xreproject
! character(len=:), allocatable :: cmd
! logical                :: be_verbose=.false.
! real,    parameter     :: SHMAG=1.0
! integer, parameter     :: N_PTCLS = 9
! real,    allocatable   :: corrs(:), norm_const(:, :)
! real                   :: corrmax, corr, cxy(3), lims(2,2), sh(2)
! integer                :: xsh, ysh, xbest, ybest, i, irot, rc
! real, allocatable      :: sigma2_noise(:,:)      !< the sigmas for alignment & reconstruction (from groups)
! logical                :: mrc_exists
! if( command_argument_count() < 3 )then
!     write(logfhandle,'(a)',advance='no') 'ERROR! Usage: simple_test_shiftsrch stk=<particles.ext> mskdiam=<mask radius(in pixels)>'
!     write(logfhandle,'(a)') ' smpd=<sampling distance(in A)> [nthr=<number of threads{1}>] [verbose=<yes|no{no}>]'
!     write(logfhandle,'(a)') 'Example: https://www.rcsb.org/structure/1jyx with smpd=1. mskdiam=180'
!     write(logfhandle,'(a)') 'DEFAULT TEST (example above) is running now...'
!     inquire(file="1JYX.mrc", exist=mrc_exists)
!     if( .not. mrc_exists )then
!         write(*, *) 'Downloading the example dataset...'
!         cmd = 'curl -s -o 1JYX.pdb https://files.rcsb.org/download/1JYX.pdb'
!         call execute_command_line(cmd, exitstat=rc)
!         write(*, *) 'Converting .pdb to .mrc...'
!         cmd = 'e2pdb2mrc.py 1JYX.pdb 1JYX.mrc'
!         call execute_command_line(cmd, exitstat=rc)
!         cmd = 'rm 1JYX.pdb'
!         call execute_command_line(cmd, exitstat=rc)
!         write(*, *) 'Projecting 1JYX.mrc...'
!         call cline_projection%set('vol1'      , '1JYX.mrc')
!         call cline_projection%set('smpd'      , 1.)
!         call cline_projection%set('pgrp'      , 'c1')
!         call cline_projection%set('mskdiam'   , 180.)
!         call cline_projection%set('nspace'    , 6.)
!         call cline_projection%set('nthr'      , 16.)
!         call xreproject%execute(cline_projection)
!         call cline%set('stk'    , 'reprojs.mrcs')
!         call cline%set('smpd'   , 1.)
!         call cline%set('nthr'   , 16.)
!         call cline%set('stk'    , 'reprojs.mrcs')
!         call cline%set('mskdiam', 180.)
!     endif
! endif
! call cline%parse_oldschool
! call cline%checkvar('stk',      1)
! call cline%checkvar('mskdiam',  2)
! call cline%checkvar('smpd',     3)
! call cline%check
! be_verbose = .false.
! if( cline%defined('verbose') )then
!     if( trim(cline%get_carg('verbose')) .eq. 'yes' )then
!         be_verbose = .true.
!     endif
! endif
! call p%new(cline)
! p%kfromto(1) = 2
! p%kfromto(2) = 40
! allocate( sigma2_noise(p%kfromto(1):p%kfromto(2), 1:N_PTCLS), source=1. )
! call b%build_general_tbox(p, cline)
! call pftc%new(N_PTCLS, [1,N_PTCLS], p%kfromto)
! call pftc%assign_sigma2_noise(sigma2_noise)
! allocate(corrs(pftc%get_nrots()), norm_const(pftc%get_nrots(), 2))
! call img_copy%new([p%box_crop,p%box_crop,1],p%smpd_crop)
! call img_copy%init_polarizer(pftc, p%alpha)
! call b%img%read(p%stk, 1)
! call b%img%norm
! call b%img%fft
! call b%img%clip_inplace([p%box_crop,p%box_crop,1])
! call img_copy%polarize(pftc, b%img, 1, isptcl=.false., iseven=.true., mask=b%l_resmsk)
! call img_copy%polarize(pftc, b%img, 1, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
! call pftc%shift_ptcl(1, [SHMAG,0.,0.]) ! left
! call img_copy%polarize(pftc, b%img, 2, isptcl=.false., iseven=.true., mask=b%l_resmsk)
! call img_copy%polarize(pftc, b%img, 2, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
! call pftc%shift_ptcl(2, [0.,SHMAG,0.]) ! down
! call img_copy%polarize(pftc, b%img, 3, isptcl=.false., iseven=.true., mask=b%l_resmsk)
! call img_copy%polarize(pftc, b%img, 3, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
! call pftc%shift_ptcl(3, [-SHMAG,0.,0.]) ! right
! call img_copy%polarize(pftc, b%img, 4, isptcl=.false., iseven=.true., mask=b%l_resmsk)
! call img_copy%polarize(pftc, b%img, 4, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
! call pftc%shift_ptcl(4, [0.,-SHMAG,0.]) ! up
! call img_copy%polarize(pftc, b%img, 5, isptcl=.false., iseven=.true., mask=b%l_resmsk)
! call img_copy%polarize(pftc, b%img, 5, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
! call pftc%gen_sigma_contrib(5,5,[SHMAG,SHMAG],1,sigma2_noise(:,5))
! call pftc%assign_sigma2_noise(sigma2_noise)
! call pftc%shift_ptcl(5, [SHMAG,SHMAG,0.]) ! left + down
! call img_copy%polarize(pftc, b%img, 6, isptcl=.false., iseven=.true., mask=b%l_resmsk)
! call img_copy%polarize(pftc, b%img, 6, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
! call pftc%shift_ptcl(6, [-SHMAG,-SHMAG,0.]) ! right + up
! call img_copy%polarize(pftc, b%img, 7, isptcl=.false., iseven=.true., mask=b%l_resmsk)
! call img_copy%polarize(pftc, b%img, 7, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
! call pftc%shift_ptcl(7, [-SHMAG,SHMAG,0.]) ! right + down
! call img_copy%polarize(pftc, b%img, 8, isptcl=.false., iseven=.true., mask=b%l_resmsk)
! call img_copy%polarize(pftc, b%img, 8, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
! call pftc%shift_ptcl(8, [SHMAG,-SHMAG,0.]) ! left + up
! call img_copy%polarize(pftc, b%img, 9, isptcl=.false., iseven=.true., mask=b%l_resmsk)
! call img_copy%polarize(pftc, b%img, 9, isptcl=.true.,  iseven=.true., mask=b%l_resmsk)
! call pftc%shift_ptcl(9, [0.,0.,0.]) ! no shift
! call pftc%set_with_ctf(.false.)
! call b%img%ifft
! call b%img%read(p%stk, 5)
! call b%img%norm
! call b%img%fft
! call b%img%clip_inplace([p%box_crop,p%box_crop,1])
! call img_copy%fft
! call img_copy%polarize(pftc, b%img, 1, isptcl=.false., iseven=.true., mask=b%l_resmsk)
! call pftc%memoize_refs
! do i = 1, N_PTCLS
!     call pftc%memoize_sqsum_ptcl(i)
! enddo
! call pftc%memoize_ptcls
! lims(1,1) = -6.
! lims(1,2) =  6.
! lims(2,1) = -6.
! lims(2,2) =  6.
! call grad_shsrch_obj%new(lims, opt_angle=.false.)
! call grad_shsrch_obj%set_indices(5, 5)
! irot = 1
! cxy  = grad_shsrch_obj%minimize(irot)
! print *, cxy(1), cxy(2:3), irot
! params_glob%nstates = 2
! do i=5,5
!     call pftc%gen_corrs(i, i, corrs)
!     print *, 'corr: ', maxval(corrs)
!     corrmax = 0.
!     do xsh=-2,2
!         do ysh=-2,2
!             call pftc%gen_corrs(i, i, real([xsh,ysh]), corrs)
!             corr  = maxval(corrs)

!             print *, 'corr: ', corr, xsh, ysh

!             if( corr > corrmax )then
!                 corrmax = corr
!                 xbest   = xsh
!                 ybest   = ysh
!             endif
!         enddo
!     enddo
!     print *, xbest, ybest, corrmax
! enddo
! call pftc%calc_shift(5, 5, sh, rot_in=1)
! print *, 'calculated shift = ', sh

! contains

!     subroutine calc_shift( self, iref, iptcl, sh, rot_in )
!         class(polarft_calc), intent(inout) :: self
!         integer,                 intent(in)    :: iref
!         integer,                 intent(in)    :: iptcl
!         real,                    intent(inout) :: sh(2)
!         integer,  optional,      intent(in)    :: rot_in
!         complex,  pointer   :: pft_ref(:,:), pft_ref_tmp(:,:)
!         real(dp), pointer   :: args1(:,:), args2(:,:)
!         integer,  parameter :: NPI   = 5,&       ! number of trigonometry periods
!                               &NLINS = 2,&       ! number of 2x2 linear systems
!                               &NEQS  = 2*NLINS
!         integer :: ind, ithr, i, j, k, irot, cand1_cnt, cand2_cnt, eq_cnt, k_cands(NEQS), r_cands(NEQs)
!         logical :: rots(self%pftsz), ks(self%kfromto(1):self%kfromto(2))
!         complex :: AB
!         real    :: C_1st, T_1st, C_2nd, T_2nd, RHS_1st(NPI), RHS_2nd(NPI), x_comp, y_comp, trs, abspft,&
!                 &cand1_x(NPI*NPI), cand1_y(NPI*NPI), cand2_x(NPI*NPI), cand2_y(NPI*NPI), minval, minx, miny, val
!         ind  = self%pinds(iptcl)
!         ithr = omp_get_thread_num() + 1
!         pft_ref     => self%heap_vars(ithr)%pft_ref
!         pft_ref_tmp => self%heap_vars(ithr)%pft_ref_tmp
!         args1       => self%heap_vars(ithr)%pft_r1_8
!         args2       => self%heap_vars(ithr)%pft_r2_8
!         if( present(rot_in) )then
!             if( self%iseven(ind) )then
!                 pft_ref_tmp = self%pfts_refs_even(:,:,iref)
!             else
!                 pft_ref_tmp = self%pfts_refs_odd(:,:,iref)
!             endif
!             call self%rotate_pft(pft_ref_tmp, rot_in, pft_ref)
!             call self%rotate_pft(self%argtransf(           1:self%pftsz,:), rot_in, args1)
!             call self%rotate_pft(self%argtransf(self%pftsz+1:          ,:), rot_in, args2)
!         else
!             if( self%iseven(ind) )then
!                 pft_ref = self%pfts_refs_even(:,:,iref)
!             else
!                 pft_ref = self%pfts_refs_odd(:,:,iref)
!             endif
!             args1 = self%argtransf(           1:self%pftsz,:)
!             args2 = self%argtransf(self%pftsz+1:          ,:)
!         endif
!         if( self%with_ctf ) pft_ref = pft_ref * self%ctfmats(:,:,ind)
!         trs    = params_glob%trs
!         eq_cnt = 0
!         ! generating candidates for equations
!         rots = .false.
!         ks   = .false.
!         do k = self%kfromto(1), self%kfromto(2)
!             if( ks(k) )cycle
!             do irot = 1, self%pftsz
!                 if( rots(irot) .or. ks(k) )cycle
!                 if( eq_cnt > NEQS )cycle
!                 abspft = real(self%pfts_ptcls(irot,k,ind) * conjg(self%pfts_ptcls(irot,k,ind)))
!                 if( abspft < TINY )cycle
!                 abspft = real(pft_ref(irot,k) * conjg(pft_ref(irot,k)))
!                 if( abspft < TINY )cycle
!                 eq_cnt          = eq_cnt + 1
!                 k_cands(eq_cnt) = k
!                 r_cands(eq_cnt) = irot
!                 rots(irot)      = .true.
!                 ks(k)           = .true.
!             enddo
!         enddo
!         ! not enough candidates to calculate the shifts
!         if( eq_cnt < NEQS )then
!             sh = 0.
!             return
!         endif
!         ! first pair
!         k       = k_cands(1)
!         irot    = r_cands(1)
!         AB      = self%pfts_ptcls(irot,k,ind)/pft_ref(irot,k)
!         RHS_1st = atan(aimag(AB)/real(AB))
!         C_1st   = args1(irot,k)
!         T_1st   = args2(irot,k)
!         do j = 1, NPI
!             RHS_1st(j) = RHS_1st(j) + (j-1-NPI/2)*PI
!         enddo
!         k       = k_cands(2)
!         irot    = r_cands(2)
!         AB      = self%pfts_ptcls(irot,k,ind)/pft_ref(irot,k)
!         RHS_2nd = atan(aimag(AB)/real(AB))
!         C_2nd   = args1(irot,k)
!         T_2nd   = args2(irot,k)
!         do j = 1, NPI
!             RHS_2nd(j) = RHS_2nd(j) + (j-1-NPI/2)*PI
!         enddo
!         cand1_cnt = 0
!         do i = 1, NPI
!             do j = 1, NPI
!                 x_comp = -(RHS_1st(i)*T_2nd - RHS_2nd(j)*T_1st)/(T_1st*C_2nd - T_2nd*C_1st)
!                 y_comp =  (RHS_1st(i)*C_2nd - RHS_2nd(j)*C_1st)/(T_1st*C_2nd - T_2nd*C_1st)
!                 if( x_comp < -trs .or. x_comp > trs .or. y_comp < -trs .or. y_comp > trs )cycle
!                 cand1_cnt          = cand1_cnt + 1
!                 cand1_x(cand1_cnt) = x_comp
!                 cand1_y(cand1_cnt) = y_comp
!             enddo
!         enddo
!         ! another pair
!         k       = k_cands(3)
!         irot    = r_cands(3)
!         AB      = self%pfts_ptcls(irot,k,ind)/pft_ref(irot,k)
!         RHS_1st = atan(aimag(AB)/real(AB))
!         C_1st   = args1(irot,k)
!         T_1st   = args2(irot,k)
!         do j = 1, NPI
!             RHS_1st(j) = RHS_1st(j) + (j-1-NPI/2)*PI
!         enddo
!         k       = k_cands(4)
!         irot    = r_cands(4)
!         AB      = self%pfts_ptcls(irot,k,ind)/pft_ref(irot,k)
!         RHS_2nd = atan(aimag(AB)/real(AB))
!         C_2nd   = args1(irot,k)
!         T_2nd   = args2(irot,k)
!         do j = 1, NPI
!             RHS_2nd(j) = RHS_2nd(j) + (j-1-NPI/2)*PI
!         enddo
!         cand2_cnt = 0
!         do i = 1, NPI
!             do j = 1, NPI
!                 x_comp = -(RHS_1st(i)*T_2nd - RHS_2nd(j)*T_1st)/(T_1st*C_2nd - T_2nd*C_1st)
!                 y_comp =  (RHS_1st(i)*C_2nd - RHS_2nd(j)*C_1st)/(T_1st*C_2nd - T_2nd*C_1st)
!                 if( x_comp < -trs .or. x_comp > trs .or. y_comp < -trs .or. y_comp > trs )cycle
!                 cand2_cnt          = cand2_cnt + 1
!                 cand2_x(cand2_cnt) = x_comp
!                 cand2_y(cand2_cnt) = y_comp
!             enddo
!         enddo
!         ! finding the min between cand1 and cand2
!         minval = huge(minval)
!         do i = 1, cand1_cnt
!             do j = 1, cand2_cnt
!                 val = sqrt((cand1_x(i)-cand2_x(j))**2 + (cand1_y(i)-cand2_y(j))**2)
!                 if( val < minval )then
!                     minval = val
!                     minx   = (cand1_x(i)+cand2_x(j)) / 2.
!                     miny   = (cand1_y(i)+cand2_y(j)) / 2.
!                 endif
!             enddo
!         enddo
!         sh = [minx, miny]
!     end subroutine calc_shift

end program simple_test_shiftsrch
