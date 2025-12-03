program simple_test_shift
! !$ use omp_lib
! !$ use omp_lib_kinds
! include 'simple_lib.f08'
! use simple_polarft_calc,  only: polarft_calc
! use simple_cmdline,           only: cmdline
! use simple_builder,           only: builder
! use simple_image,             only: image
! use simple_parameters,        only: parameters
! use simple_polarizer,         only: polarizer
! use simple_pftc_shsrch_grad, only: pftc_shsrch_grad  ! gradient-based in-plane angle and shift search
! use simple_strategy2D3D_common
! use simple_simulator
! use simple_ctf
! use simple_ori
! use simple_classaverager
! use simple_euclid_sigma2
! implicit none
! type(cmdline)                 :: cline
! type(builder)                 :: b
! type(parameters)              :: p
! type(polarft_calc)        :: pftc
! type(pftc_shsrch_grad)       :: grad_shsrch_obj           !< origin shift search object, L-BFGS with gradient
! type(ctf)                     :: tfun
! type(ori)                     :: o
! type(oris)                    :: os
! type(image), allocatable      :: match_imgs(:), ptcl_match_imgs(:)
! ! character(len=:), allocatable :: cmd
! logical                :: be_verbose=.false.
! real,    parameter     :: SHMAG = 3.0
! real,    parameter     :: SNR   = 0.1
! real,    parameter     :: BFAC  = 10.
! integer, parameter     :: SH_ITERS = 5, N_SH = 5, N_PTCLS = N_SH**2
! integer, allocatable   :: pinds(:)
! type(ctfparams)        :: ctfparms
! type(euclid_sigma2)    :: eucl
! real, allocatable      :: sigma2_group(:,:,:), truth_sh(:,:)
! real                   :: cxy(3), lims(2,2), lims_init(2,2), rnd_shifts(2,N_SH), sh, shifts_cnt(N_SH), shifts_dist(N_SH),&
!                          &min_truth, correct_cnt, truth_cnt(N_SH)
! integer                :: xsh, ysh, xbest, ybest, i, irot, ne, no, iptcl, nptcls2update, ithr, iter
! logical                :: mrc_exists
! if( command_argument_count() < 4 )then
!     write(logfhandle,'(a)',advance='no') 'ERROR! required arguments: '

! endif
! call cline%parse_oldschool
! call cline%checkvar('stk',      1)
! call cline%checkvar('mskdiam',  2)
! call cline%checkvar('smpd',     3)
! call cline%checkvar('lp',       4)
! call cline%check
! call cline%set('oritype','ptcl2D')
! if( .not.cline%defined('objfun') ) call cline%set('objfun', 'euclid')
! call cline%set('ml_reg', 'no')
! call cline%set('ncls',   1.)
! call cline%set('kv',     300)
! call cline%set('cs',     2.7)
! call cline%set('fraca',  0.1)
! call cline%set('nptcls',  N_PTCLS)
! be_verbose = .false.

! ! general input
! call b%init_params_and_build_strategy2D_tbox(cline, p)
! call b%spproj%projinfo%new(1,is_ptcl=.false.)
! call b%spproj%projinfo%set(1,'projname', 'test')
! call b%spproj%os_ptcl2D%kill
! call b%spproj%os_ptcl3D%kill
! p%fromp  = 1
! p%top    = p%nptcls
! p%frcs   = trim(FRCS_FILE)
! ctfparms%smpd  = p%smpd
! ctfparms%kv    = p%kv
! ctfparms%cs    = p%cs
! if( trim(p%ctf) .eq. 'no' )then
!     ctfparms%ctfflag = CTFFLAG_NO
! endif
! ctfparms%fraca = p%fraca
! tfun = ctf(p%smpd, p%kv, p%cs, p%fraca)

! ! generate particles
! call b%img%read(p%stk, p%iptcl)
! call prepimgbatch(N_PTCLS)
! call os%new(p%nptcls,is_ptcl=.true.)
! if( trim(p%ctf) .eq. 'no' )then
! else
!     call os%rnd_ctf(p%kv, p%cs, p%fraca, 2.5, 1.5, 0.001)
! endif
! call seed_rnd
! do i = 1, N_SH
!     sh = SHMAG * real(i - 1) / real(N_SH - 1)
!     if( ran3() < 0.5 ) sh = -sh
!     rnd_shifts(1,i) = sh
!     if( ran3() < 0.5 ) sh = -sh
!     rnd_shifts(2,i) = sh
! enddo
! allocate(truth_sh(p%fromp:p%top,2))
! do iptcl = p%fromp,p%top
!     call os%set(iptcl,'state',1.)
!     call os%set(iptcl,'w',    1.)
!     call os%set(iptcl,'class',1.)
!     call b%img%read(p%stk, 1)
!     truth_sh(iptcl,:) = rnd_shifts(:,floor(ran3() * N_SH) + 1)
!     call os%set(iptcl,'x', truth_sh(iptcl,1))
!     call os%set(iptcl,'y', truth_sh(iptcl,2))
!     call os%get_ori(iptcl, o)
!     call b%img%pad(b%img_pad)
!     call b%img_pad%fft
!     call b%img_pad%shift2Dserial(truth_sh(iptcl,:) )
!     call simimg(b%img_pad, o, tfun, p%ctf, SNR, bfac=BFAC)
!     call b%img_pad%clip(b%imgbatch(iptcl))
!     call b%imgbatch(iptcl)%write('particles.mrc',iptcl)
! enddo
! do iptcl = p%fromp,p%top
!     call os%set(iptcl,'x', 0.)
!     call os%set(iptcl,'y', 0.)
! enddo
! call b%spproj%add_single_stk('particles.mrc', ctfparms, os)
! call b%spproj_field%partition_eo
! call b%spproj_field%sample4update_all([p%fromp,p%top],nptcls2update, pinds, .true.)

! ! pftc
! call pftc%new(p%nptcls, [1,p%nptcls], p%kfromto)
! call eucl%new('dummy.dat', p%box)
! call eucl%allocate_ptcls
! allocate(match_imgs(p%ncls),ptcl_match_imgs(nthr_glob))
! call pftc%reallocate_ptcls(p%nptcls, pinds)
! do ithr = 1,nthr_glob
!     call ptcl_match_imgs(ithr)%new([p%box_crop, p%box_crop, 1], p%smpd_crop, wthreads=.false.)
! enddo

! ! references
! call restore_read_polarize_cavgs(0)

! ! particles
! !$omp parallel do default(shared) private(iptcl,ithr)&
! !$omp schedule(static) proc_bind(close)
! do iptcl = 1,p%nptcls
!     ithr  = omp_get_thread_num() + 1
!     call prepimg4align(iptcl, b%imgbatch(iptcl), ptcl_match_imgs(ithr))
!     call b%imgbatch(iptcl)%ifft
!     call b%img_crop_polarizer%polarize(pftc, ptcl_match_imgs(ithr), iptcl, .true., .true.)
!     call pftc%set_eo(iptcl, nint(b%spproj_field%get(iptcl,'eo'))<=0 )
! end do
! !$omp end parallel do
! call pftc%create_polar_absctfmats(b%spproj, 'ptcl2D')

! ! prep for shift search
! lims(:,1)       = -p%trs
! lims(:,2)       =  p%trs
! lims_init(:,1)  = -SHC_INPL_TRSHWDTH
! lims_init(:,2)  =  SHC_INPL_TRSHWDTH
! call grad_shsrch_obj%new(lims, lims_init=lims_init, maxits=p%maxits_sh, opt_angle=.true., coarse_init=.false.)

! ! initial sigma2
! allocate( sigma2_group(2,1,1:fdim(p%box)-1), source=0. )
! shifts_cnt = 0
! min_truth  = 0.
! do iter = 1, SH_ITERS
!     ne = 0
!     no = 0
!     do iptcl = p%fromp,p%top
!         call b%spproj_field%get_ori(iptcl, o)
!         call eucl%calc_sigma2(pftc, iptcl, o, 'class')
!         if( o%get_eo() == 0 )then
!             ne = ne+1
!             sigma2_group(1,1,:) = sigma2_group(1,1,:) + eucl%sigma2_part(:,iptcl)
!         else
!             no = no+1
!             sigma2_group(2,1,:) = sigma2_group(2,1,:) + eucl%sigma2_part(:,iptcl)
!         endif
!     enddo
!     sigma2_group(1,:,:) = sigma2_group(1,:,:) / real(ne)
!     sigma2_group(2,:,:) = sigma2_group(2,:,:) / real(no)
!     call write_groups_starfile(sigma2_star_from_iter(iter-1), sigma2_group, 1)
!     call eucl%read_groups(b%spproj_field)
!     do iptcl = p%fromp,p%top
!         call pftc%memoize_sqsum_ptcl(iptcl)
!     enddo
!     call pftc%memoize_ptcls

!     ! shift search
!     do iptcl = p%fromp,p%top
!         call grad_shsrch_obj%set_indices(1, iptcl)
!         irot = 1 ! zero angle
!         cxy  = grad_shsrch_obj%minimize(irot=irot)
!         if( iptcl == 4 )then
!             print *,'irot  ',     irot
!             print *,'score ',     cxy(1)
!             print *,'cur shift ', b%spproj_field%get_2Dshift(iptcl)
!             print *,'opt shift ', cxy(2:3)
!             print *,'truth ',     truth_sh(iptcl,:)
!         endif
!         call b%spproj_field%set_shift(iptcl, cxy(2:3)) !!
!         ! correct shifts statistics
!         if( iter == SH_ITERS )then
!             do i = 1, N_SH
!                 shifts_dist(i) = euclid(cxy(2:3), rnd_shifts(:,i))
!             enddo
!             min_truth                         = min_truth + minval(shifts_dist)
!             shifts_cnt(minloc(shifts_dist,1)) = shifts_cnt(minloc(shifts_dist,1)) + 1
!         endif
!     enddo

!     call restore_read_polarize_cavgs(iter)
! enddo

! print *, 'SHIFT COUNTS with total/total_min_dist = ', sum(shifts_cnt), min_truth
! do i = 1, N_SH
!     print *, shifts_cnt(i)
! enddo

! truth_cnt = 0
! min_truth = 0.
! do iptcl = p%fromp,p%top
!     call b%spproj_field%set_shift(iptcl, truth_sh(iptcl,:))
!     do i = 1, N_SH
!         shifts_dist(i) = euclid(truth_sh(iptcl,:), rnd_shifts(:,i))
!     enddo
!     min_truth                        = min_truth + minval(shifts_dist)
!     truth_cnt(minloc(shifts_dist,1)) = truth_cnt(minloc(shifts_dist,1)) + 1
! enddo

! correct_cnt = 0
! print *, 'TRUTH SHIFT COUNTS with total/total_min_dist = ', sum(truth_cnt), min_truth
! do i = 1, N_SH
!     print *, truth_cnt(i)
!     correct_cnt = correct_cnt + min(truth_cnt(i), shifts_cnt(i))
! enddo

! print *, 'correct shifts cnt / total = ', correct_cnt, ' / ', sum(truth_cnt)

! call restore_read_polarize_cavgs(iter)


! contains

!     subroutine restore_read_polarize_cavgs( iter )
!         integer, intent(in) :: iter
!         p%which_iter = iter
!         call cavger_kill()
!         call cavger_new
!         call cavger_transf_oridat( b%spproj )
!         call cavger_assemble_sums( .false. )
!         call cavger_merge_eos_and_norm
!         call cavger_calc_and_write_frcs_and_eoavg(p%frcs, p%which_iter)
!         p%refs      = trim(CAVGS_ITER_FBODY)//int2str_pad(p%which_iter,3)//p%ext
!         p%refs_even = trim(CAVGS_ITER_FBODY)//int2str_pad(p%which_iter,3)//'_even'//p%ext
!         p%refs_odd  = trim(CAVGS_ITER_FBODY)//int2str_pad(p%which_iter,3)//'_odd'//p%ext
!         call cavger_write(trim(p%refs),      'merged')
!         call cavger_write(trim(p%refs_even), 'even'  )
!         call cavger_write(trim(p%refs_odd),  'odd'   )
!         call b%clsfrcs%read(FRCS_FILE)
!         call cavger_read(trim(p%refs_even), 'even' )
!         call cavger_read(trim(p%refs_odd),  'odd' )
!         call b%img_crop_polarizer%init_polarizer(pftc, p%alpha)
!         call match_imgs(1)%new([p%box_crop, p%box_crop, 1], p%smpd_crop, wthreads=.false.)
!         call prep2Dref(cavgs_even(1), match_imgs(1), 1, center=.false.)
!         call b%img_crop_polarizer%polarize(pftc, match_imgs(1), 1, isptcl=.false., iseven=.true.)
!         call prep2Dref(cavgs_odd(1), match_imgs(1), 1, center=.false.)
!         call b%img_crop_polarizer%polarize(pftc, match_imgs(1), 1, isptcl=.false., iseven=.false.)
!         call pftc%memoize_refs
!     end subroutine restore_read_polarize_cavgs

end program simple_test_shift