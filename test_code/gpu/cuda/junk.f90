
  write(*,*) "before the job_descr%set"
      write(*,*) "qsys_cleanup_iter"
       write(*,*) "before the job_descr%set"
       write(*,*) "before the qscripts%generate_scripts"
       write(*,*) "algnfbody: ", algnfbody
       write(*,*) "before the merge align doc setter"
       write(*,*) "before the merge align doc execute"
       write(*,*) "before the qscripts%schedule_jobs"
       write(*,*) "before the merge align doc execute"
       write(*,*) "before the xprime3D_init"
       write(*,*) "after the xprime3D_init"


       
  job_descr%chash2str() 


  !removed from the commnad line class
  funit = get_fileunit( )
  open(unit=funit, file='cmdline.txt', status='replace',&
  iostat=file_stat, action='write', form='formatted')

  write(funit,*) trim(self%cmds(i)%key), ' ', self%cmds(i)%rarg
  write(funit,*) trim(self%cmds(i)%key), ' ', trim(self%cmds(i)%carg)
  close(funit)
 
  write(*,*)'************************************************************'
  write(*,*)' Device fills in object(devQ) and the data structure(devD)  '
  write(*,*)'************************************************************'
  err  = get_dev_count_c(devD)
  ndev = devD%ndev
  allocate(a_devD(0:devD%ndev-1))
  do idev = 0, ndev-1
     call devQ%new_deviceQuery_gpu(devD,idev)
     !call Sanity_check_gpu(devQ, devD)
     !mapping the data structures into an array
     a_devD(idev) = devD
  end do

  integer                         :: get_polarft_corr_gpu_c_
  integer                         :: get_dev_count_c


  !logic flippers
    procedure :: deviceQuery_logic

    subroutine deviceQuery_logic(self,ncnt,ndev,a_devD)
    class(build), intent(in) :: self
    integer :: ncnt
    integer :: ndev
    type(deviceDetails)    :: a_devD(0:ndev-1)
    !local variables
    !indexers
    integer :: idev
    !counters
    integer :: counter
    !start of the execution commands
    counter = 0
    do idev = 0,ndev-1
       if (self%a_devD(idev)%is_SMsuitable == 1 ) then
          counter = counter + 1
          write(*,'(5x,i2,20x,i1,20x,i1)') idev,&
               self%a_devD(idev)%is_SMsuitable,&
               self%a_devD(idev)%is_ecc
          has_gpu = .true.
       end if
    end do
    ncnt = counter
    return
  end subroutine deviceQuery_logic




  write(*,*) "in calc_corrs_on_gpu: at start "
        write(*,*) "in calc_corrs_on_gpu: before expand dimn"
           write(*,*) "in calc_corrs_on_gpu: before gencorr_all"
                   write(*,*) "in before calc_corrs_on_gpu:"
                   write(*,*) "in after calc_corrs_on_gpu:"
            write(*,*) "in calc_corrs_on_gpu: after gencorr_all"

  write(*,*) "has_gpu",has_gpu
  write(*,*) "use_gpu",use_gpu
  write(*,*) "has_multi_gpu",has_multi_gpu
  
  write(*,*) "shc"
  write(*,*) "in calc_corrs_on_gpu just before self%srch_common%calc_corrs_on_gpu" 
  write(*,*) "else"
  write(*,*) "in calc_corrs_on_gpu just before self%srch_common%calc_corrs_on_gpu" 


  interface exec_mask
     module procedure mask_soft, ,mask_hard 
  end interface exec_mask

  subroutine mask_soft(self, cline , flag )
    include 'exec_mask-inc.f90'

!in
    
    
    return
  end subroutine mask_soft
  


  read(unit,'(3(a,x,a,x),a,x,i1,x,f7.4,x,f8.3,x,i5,x,i5,x,i5,x,i3,a,x,f12.4,x,a,x,a)') &

  write(*,*)'******************************************************************'
  write(*,*)'   Device fills in the object(devQ) and the data structure(devD)  '
  write(*,*)'******************************************************************'
  rc = get_dev_count_c(devD)
  ndev = devD%ndev
  allocate(a_devD(0:devD%ndev-1))
  allocate(devname(0:devD%ndev-1,0:DEVNAME_STRING_LENGTH))
  do idev = 0, ndev-1
     call devQ%new_deviceQuery_gpu(devD,idev)
     devname(idev,:) = devQ%get_devname(0)
     devname(idev,:) = devname(idev,0:strlen(devname(idev,:)))
     !call Sanity_check_gpu(devQ, devD)
     !mapping the data structures into an array
     a_devD(idev) = devD
  end do


  interface error_handling
     subroutine file_err_define(err_name,err_msg, err_id,  err_action)
       character(len=*), intent(in) :: err_name !name of the error
       character(len=*), intent(in) :: err_msg  !error message
       integer, intent(out) :: err_id           !code of the error
       character(len=*), intent(in), optional :: err_action !Yet to be tested
     end subroutine file_err_define
  end interface error_handling

  interface error_handling

     subroutine file_err_define(err_name,err_msg, err_id,  err_action)
       character(len=*), intent(in) :: err_name !name of the error
       character(len=*), intent(in) :: err_msg  !error message
       integer, intent(out) :: err_id           !code of the error
       character(len=*), intent(in), optional :: err_action !Yet to be tested
     end subroutine file_err_define
     
     subroutine file_err_throw(err_msg, err_id, err_name)
       use simple_yaml_strings
       integer, intent(in), optional :: err_id      !error code to be raised.
                                                 ! already defined f_err_define
       character(len=*), intent(in), optional :: err_name  ! error name
       character(len=*), intent(in), optional :: err_msg   ! error message
     end subroutine file_err_throw

  end interface error_handling

  
  write(*,*) unt," ",trim(file)," ", f_status," ", f_form," ",f_position," ",f_action

  write(*,*) "after file_unit: ",file,unit,unt


  write(*,*) "tmr_name: ",tmr_name
  !err = convert_int2char_pos_c(char_out,i)
   
  
  
  write(*,*) "ibench = ",ibench
  write(*,*) "ibench_write = ",ibench_write
  write(*,*) "b%ibench = ",b%s_bench%bench_i
  write(*,*) "b%ibench_write = ",b%s_bench%bench_write_i




  
  double precision              :: elps_prep_cpu
  double precision,dimension(2) :: st_prep_cpu, et_prep_cpu


  call gettimeofday_c(st_prep_cpu)

  call gettimeofday_c(et_prep_cpu)
  call elapsed_time_c(st_prep_cpu,et_prep_cpu,elps_prep_cpu)

  write(*,*) "elapsed_time preppftcc4align(secs): ",elps_prep_cpu



  character(len=80),allocatable :: char_primeExec_tmr(:)

   allocate(char_primeExec_tmr(p%maxits))
   write(*,*) "length: ",length, "p%maxits: ",p%maxits

      write(*,*) char_out

      write(*,*)"tmr_name: ", tmr_name
      char_primeExec_tmr(i) = tmr_name
      write(*,*) "char_primeExec_tmr(",i,"):",char_primeExec_tmr(i)

deallocate(char_primeExec_tmr)




  write(*,*) " In stochastic_srch_gpu in line 808 inpl_ind_new: ",&
       self%proj_space(self%proj_space_inds(self%nrefs))%inpl_ind

  write(*,*) " In stochastic_srch_gpu in line 808 inpl_ind_old: ",self%proj_space(self%nrefs)%inpl_ind

  write(*,*) "line 779", cnt_glob,ref, inpl


!!! in polarft_corrcacl

      write(*,*) "writing corrmat3dout to fort.7500"

      if ( debug_cpu .eqv. .true. ) then
         do i=1,self%nptcls
            do irot=1,self%nptcls
               do ik=1,self%nk
                  write(7500,*)i, irot, ik, corrmat3dout(i,irot,ik)
               end do
            end do
         end do
      end if

!      write(*,*) "stop at line 455"
!      stop 



      

!!! in project_corrmat3D_greedy

      write(*,*) "writing corrmat3d to fort.1100"
      write(1100,*) corrmat3d


      write(*,*) "writing inplmat to fort.2100"
      write(2100,*) inplmat
        write(*,*) "writing corrmat2d to fort.3100"
        write(3100,*) corrmat2d


!!! in project_corrmat3D_greedy

      write(*,*) "lda,ldb,ldc: ",lda,ldb,ldc

      write(1500,*) corrmat3dout
      write(2500,*) sumb_vec

        write(*,*) "get_polarft_gencorrall_gpu_c(...,Z,N,...) line 425"

        write(3500,*)self%pfts_refs
        write(4500,*)self%pfts_ptcls
        write(5500,*)self%sqsums_refs
        write(6500,*)self%sqsums_ptcls




write(*,*) " in this routine ......."



  
     do ivx=1,vx - (vx-2)
        do ivy=1,vy - (vy-2)
           do ivz=1,vz - (vz-1)
              write(*,*)ivx,ivy,ivz,cmat2sh(ivx,ivy,ivz)
           end do
        end do
     end do

     write(*,*)'                           C    N                              '
     r_gpu = 1.0
     call start_timer_cpu("recastcn_gpu")
     call gettimeofday_c(st_C_r)
     !Kernel implementation has changed the argument list has changed
     !TODO: need to fix the argument list of the function call
     err = get_carte2d_ftext_corr_gpu_c(a_devD,                     &
                                        s_carte,"C","N",            &
                                        r_gpu,                      &        
                                        cmat1,cmat2,                &
                                        vx,vy,vz,                   &
                                        lda,ldb,ldc,alpha,          &
                                        s_bench, s_debug_gpu) 
     call gettimeofday_c(et_C_r)
     call elapsed_time_c(st_C_r,et_C_r,elps_C)
     call stop_timer_cpu("recastcn_gpu")
     if (err /= RC_SUCCESS ) write(*,*) "get_carte2d_ftext_corr_gpu_c=",err

     r_gpu = calc_corr(s_carte%r_polar,s_carte%sumasq_polar*s_carte%sumbsq_polar)


     speedup_RC = elps_R/elps_C

     write(*,*) " Speed up from recats to GPU(C) : ",speedup_RC
     if (speedup_RC<1.0) write(*,*)"speedup < 1, try bigger Volume"


!        open(5,file='2Dcarte_corr_OR.asc',status='unknown',position='append')
!        open(6,file='2Dcarte_corr_RF.asc',status='unknown',position='append')
!        open(7,file='2Dcarte_corr_OF.asc',status='unknown',position='append')

!        write(5,'(1x,i5,2x,f15.8)')ipart, speedup_OR
!        write(6,'(1x,i5,2x,f15.8)')ipart, speedup_RF
!        write(7,'(1x,i5,2x,f15.8)')ipart, speedup_OF

!        close(5)
!        close(6)
!        close(7)

  
#if defined (BENCH)
   if( ibench_write ) then
      open(4,file='2Dcarte_corr_Al.asc',status='unknown',position='append')
      write(4,'(1x,a,7x,a,7x,a,4x,a,3x,a,7x,a,7x,a)')&
           "ipart","Old(secs)","Recast(secs)","Recast-GPU(secs)", &
           "speedup_OR","speedup_RG","speedup_OG"
   end if
#endif

  !closing files

#if defined (BENCH)
  if( ibench_write ) then
  end if
#endif


!>  \brief is calculating physical indicies and logical mask to prepare for corr calc
    subroutine fcorr_phys_lmsk( self, phys, lmsk, lims, lp_dyn, hp_dyn )
        class(image), intent(in)          :: self
        integer, allocatable, intent(out) :: phys(:,:,:,:)
        logical, allocatable, intent(out) :: lmsk(:,:,:)
        integer, intent(out)              :: lims(3,2)
        real, intent(in), optional        :: lp_dyn, hp_dyn
        real                              :: r, sumasq, sumbsq
        integer                           :: h, hh, k, kk, l, ll, sqarg, sqlp, sqhp
        if( allocated(phys) ) deallocate(phys)
        if( allocated(lmsk) ) deallocate(lmsk)
        if( present(lp_dyn) )then
            lims = self%fit%loop_lims(1,lp_dyn)
        else
            lims = self%fit%loop_lims(2) ! Nyqvist default low-pass limit
         endif
        allocate( phys(lims(1,1):lims(1,2),lims(2,1):lims(2,2),lims(3,1):lims(3,2),3),&
                  lmsk(self%ldim(1),self%ldim(2),self%ldim(3)) )
        phys = 0
        lmsk = .false.
        sqlp = (maxval(lims(:,2)))**2
        if( present(hp_dyn) )then
            sqhp = max(2,self%get_find(hp_dyn))**2
        else
            sqhp = 2 ! index 2 default high-pass limit
         endif
        !$omp parallel do default(shared) private(h,hh,k,kk,l,ll,sqarg) schedule(auto)
        do h=lims(1,1),lims(1,2)
            hh = h*h
            do k=lims(2,1),lims(2,2)
                kk = k*k
                do l=lims(3,1),lims(3,2)
                    ll = l*l
                    sqarg = hh+kk+ll
                    if( sqarg <= sqlp .and. sqarg >= sqhp  )then
                        phys(h,k,l,:) = self%fit%comp_addr_phys([h,k,l])
                        lmsk(h,k,l)   = .true.
                     endif
                end do
            end do
        end do
        !$omp end parallel do
         write(*,*) "after the omp loop"
    end subroutine fcorr_phys_lmsk


  
!$omp parallel default(shared) private(iptcl)
!$omp workshare 
corrmat2d_tmp = maxval(corrmat3d_new, dim=3)
inplmat2d_tmp = maxloc(corrmat3d_new, dim=3)
corrmat2d_old = maxval(corrmat3d_old, dim=3)
inplmat2d_old = maxloc(corrmat3d_old, dim=3)
!$omp end workshare
!!!$omp do schedule(auto) 
!!do iptcl=1,p%nptcls
!!    indices(iptcl) = iptcl ! generate indices for cshifting
!!end do
!!!$omp end do nowait

! initialise validation variables
! corr      = 0.
! sdiff     = 0.
! rel_sdiff = 0.
! dist      = 0.

!$omp end parallel

!!cnt       = 0
!!#if defined (BENCH)
!!call gettimeofday_c(st_L_r)
!!#endif
!!do iref=1,p%nptcls
!!   if( iref /= 1 )then
!!      indices = cshift(indices, shift=1)
!!   endif
!!   !$omp parallel do schedule(auto) private(iptcl)
!!   do iptcl=1,p%nptcls
!!       corrmat2d_new(indices(iptcl),iref) = corrmat2d_tmp(iref,iptcl)
!!       inplmat2d_new(indices(iptcl),iref) = inplmat2d_tmp(iref,iptcl)
!!!       call metrics_pft_corrs(corrmat3d_new(indices(iptcl),iref,:),&
!!!       corrmat3d_old(iref,iptcl,:), corr, sdiff, rel_sdiff, dist)
!!      cnt = cnt+1
!!   end do
!!   !$omp end parallel do
!!end do
!!#if defined (BENCH)
!!call gettimeofday_c(et_L_r)
!!call elapsed_time_c(st_L_r,et_L_r,elps_L)
!!write(*,*)"---------------------------------------------------"
!!write(*,'(a,x,f15.8)')"Matrix mapping (secs): ", elps_L
!!write(*,*)"---------------------------------------------------"
!!#endif



  interface gettimeofday
           subroutine gettimeofday_c_(time)
             implicit none
             double precision, dimension(2) :: time
           end subroutine gettimeofday_c_
        end interface gettimeofday


        write(*,*)" hello subroutine getBlockInverseMatrix_gpu"


  write(*,*)'(x,a,x,f15.5,x,a,x,f15.5)') &
       " Time CPU (secs)          : ",elps_mtml_cpu," Time GPU (secs): ", elps_mtml_gpu

  interface
#if defined (LINUX)
     function print_s_devD_struct(s_devD)
       import :: deviceDetails
       type(deviceDetails) :: s_devD(:)
       integer :: print_s_bench_struct
     end function print_s_devD_struct
#endif
  end interface
  
#if defined (CUDA)
#else
  write(*,*)"**************************WARNING******************************"
  write(*,*)"You need to compile with -DCUDA                                "
  write(*,*)"to acces the CUDA environment computation using GPU            "
  write(*,*)"switching back to the CPU version of corr function             "
  write(*,*)"***************************************************************"
#endif



 !type(c_ptr)          :: dev_name
  do idev=1,ndev
     call c_f_pointer(a_devD(idev)%dev_name, devname, [DEVNAME_STRING_LENGTH])
     write(*,*) devname
     call c_f_pointer(a_devD(idev)%dn, devname, [DEVNAME_STRING_LENGTH])
     devname
  end do

  write(*,*)"dn               : ",a_devD(1:ndev)%dn
  

  real(sp),dimension(0:1) :: totGlobalMem_MB
  character,pointer,dimension(ndev) :: devname(:)
  
  totGlobalMem_MB(idev) = devQ%get_totGlobalMem_MB(idev)
  call c_f_pointer(devD%dev_name, devname, [DEVNAME_STRING_LENGTH])

  totGlobalMem_MB(0:ndev-1) =  a_devD(0:ndev-1)%tot_global_mem_MB
  write(*,*) totGlobalMem_MB, a_devD%tot_global_mem_MB

  write(*,*) totGlobalMem_MB, devname

  write(*,*) s_debug_gpu%debug, s_debug_gpu%debug_cpu, s_debug_gpu%debug_high, s_debug_gpu%debug_write, s_debug_gpu%debug_write_C

  !Intermediate variable for the corrmat3d=corrmat3d_1 x corrmat3d_2
  real(sp), allocatable         ::    corrmat3d_1(:,:)!cor matrix 3d _1
  real(sp), allocatable         ::    corrmat3d_2(:,:)!cor matrix 3d _2
  real(sp), allocatable         ::    corrmat3d_s(:,:,:)!cor matrix 3d stitched

  !intermdiate variables allocations for the corrmat3d_{1,2}
  allocate(  corrmat3d_1(ipart, ipart       ))
  allocate(  corrmat3d_2(ipart, nrot        ))
  allocate(  corrmat3d_s(ipart, ipart, nrot ))

  corrmat3d_1       = 0.0
  corrmat3d_2       = 0.0
  corrmat3d_s       = 0.0

  deallocate(corrmat3d_1)
  deallocate(corrmat3d_2)
  deallocate(corrmat3d_s)


  real(sp)                      :: sum_diff_crmt3s
  sum_diff_crmt3s = sum(corrmat3d-corrmat3d_s)
  write(*,*) "Differences sum(corrmat3d-corrmat3d_s)  : ",sum_diff_crmt3s

subroutine get_polarft_gencorrall_cpu(corrmat3d, corrmat3d_1, corrmat3d_2, &...)

  !Intermediate variables for the cormat 3d variables
  real(sp)              ::     corrmat3d_1(ipart,*)
  real(sp)              ::     corrmat3d_2(ipart,*)
  
  
  !$omp do schedule(auto)
  do iptcl=1,ipart
     corrmat3d(iptcl,jpart,irot) = sum(hadamard_prod(iptcl,:nrot/2,:nk))
     !corrmat3d_1(iptcl,jpart) = sum(hadamard_prod(iptcl,:nrot/2,:nk))
     !corrmat3d_2(iptcl,irot) = sum(hadamard_prod(iptcl,:nrot/2,:nk))
  end do
  !$omp end do

  !$omp workshare
  corrmat3d(:,jpart,irot) = corrmat3d(:,jpart,irot)/&
       sqrt( sqsums_refs(:ipart) * sqsums_ptcls(d1lim(1):d1lim(2)) )

  corrmat3d_1(:,jpart) = sum(hadamard_prod(iptcl,:nrot/2,:nk))/&
       sqrt( sqsums_refs(:ipart) * sqsums_ptcls(d1lim(1):d1lim(2)) )
  
  corrmat3d_2(:,irot) = sum(hadamard_prod(iptcl,:nrot/2,:nk))/&
       sqrt( sqsums_refs(:ipart) * sqsums_ptcls(d1lim(1):d1lim(2)) )
  
  !$omp end workshare nowait

  d1lim(2) = 0
    d2lim(2) = 0
    do jpart = 1,ipart
       d1lim(1) = jpart
       d1lim(2) = jpart + ipart - 1
       do irot = 1,nrot
          d2lim(1) = irot 
          d2lim(2) = irot + refsz - 1

          corrmat3d_s(:,jpart,irot) = corrmat3d_1(:,jpart)

          corrmat3d_s(:,jpart,irot) = corrmat3d_2(:,irot)
          
       end do
    end do
     
    !     do irot = 1,nrot
    !        corrmat3d_s(:,:,irot) = corrmat3d_1(:,:)
    !     end do
    !     do jpart = 1, ipart
    !        corrmat3d_s(:,jpart,:) = corrmat3d_2(:,:)
    !     end do

    where ( corrmat3d_s ==  infinity ) corrmat3d_s = 0.0
    where ( corrmat3d_s == -infinity ) corrmat3d_s = 0.0

         write(*,'(10x,a,5x,a,10x,a,4x,a,5x,a,1x,a)') &
          "ipart","irot","ik","(corrmat3d)","(corrmat3d_gpu)","(corrmat3d_s)"
     do i=1,ipart - (ipart-2)
        do irot=1,nrot - (nrot-2)
           do ik=1,nk - (nk-2)
              write(*,*)i, irot, ik, corrmat3d(i,irot,ik), corrmat3d_gpu(i,irot,ik),corrmat3d_s(i,irot,ik)
           end do
        end do
     end do

    
  do i=1,ipart
     do j=1,nrot/2
        do ik=1,nk
           write(5000,*)i,j , ik, hadamard_prod(i,j,ik)
        end do
     end do
  end do
  
  do i=1,ipart
     do j=1,ipart
        write(20,*)i, j, ik, corrmat2d(i,j), corrmat2d_gpu(i,j)
        
        do irot=1,nrot
           write(30,*)i, j, ik, corrmat3d(i,j,irot), corrmat3d_gpu(i,j,irot)
        end do
     end do
  end do

write(*,*) "N       corrmat2d(npart,npart): ",size_crmt2d
write(*,*) "N  corrmat3d(npart,npart,nrot): ",size_crmt3d

  
  write(*,*) "Relative Differences sum((corrmat3d-corrmat3d_gpu)/corrmat3d_gpu): ",sum((corrmat3d-corrmat3d_gpu)/corrmat3d_gpu)

  program infinity
    implicit none
    integer :: inf
    real :: infi
    equivalence (inf,infi) !Stores two variable at the same address
    data inf/z'7f800000'/ !hex for +infinity
    write(*,*)infi
  end program infinity

 do i=1,ipart! - (ipart-3)
       do j=1,nrot/2! - (nrot-3)
          do ik=1,nk! - (nk-3)
             write(5000,*)i,j , ik, hadamard_prod(i,j,ik)
          end do
       end do
    end do


    write(3000,*) corrmat3d
    write(4000,*) corrmat3d_gpu



  
     d1lim(2) = 0
     d2lim(2) = 0
     do jpart = 1,ipart
        d1lim(1) = jpart
        d1lim(2) = jpart + ipart - 1
        do irot = 1,nrot
           d2lim(1) = irot 
           d2lim(2) = irot + refsz - 1

           !$omp parallel default(shared) private(iptcl)
           !$omp workshare
           hadamard_prod = real(pft1(:,:,:)*&
                conjg(pft2(d1lim(1):d1lim(2),d2lim(1):d2lim(2),:)))
           !$omp end workshare

           ! correlation normalization

           !$omp do schedule(auto)
           do iptcl=1,ipart
              corrmat3d(iptcl,jpart,irot) = sum(hadamard_prod(iptcl,:,:))
           end do
           !$omp end do

           !$omp workshare
           corrmat3d(:,jpart,irot) = corrmat3d(:,jpart,irot)/&
                sqrt( sqsums_refs(:) * sqsums_ptcls(d1lim(1):d1lim(2)) )
           !$omp end workshare nowait

           !$omp end parallel

        end do
     end do



  
     write(*,*)"size(pft(npart,nrot,nk)) : ",ipart*nrot*nk/1.e6,"(Millions)"
     write(*,*)"Full sum of r is         : ",s_polar%r_polar
     write(*,*)"Full sum of sumasq is    : ",s_polar%sumasq_polar
     write(*,*)"Full sum of sumbsq is    : ",s_polar%sumbsq_polar
     write(*,*)"the correlator"
     r_gpu = calc_corr(s_polar%r_polar,s_polar%sumasq_polar*s_polar%sumbsq_polar)
     write(*,*)"after calc_corr r_gpu    : ",r_gpu
     write(*,*)"Elapsed time for corr_gpu: ",real(elps_corr_NN_gpu),"(seconds)"


     !allocating the complex matrices
     allocate(sumX_vec(ipart))

     call getCRandomGaussianDistr_3D(ipart,nrot,nk,pft1)
     call getCRandomGaussianDistr_3D(ipart,nrot,nk,pft2)

     write(*,'(10x,a,5x,a,10x,a,16x,a,26x,a)')"ipart","irot","ik","(pft1)","(pft2)"
     do i=1,ipart - (ipart-2)
        do irot=1,nrot - (nrot-2)
           do ik=1,nk - (nk-2)
              write(*,*)i, irot,ik,pft1(i,irot,ik), pft2(i,irot,ik)
           end do
        end do
     end do
     write(*,'(50x,a,12x,a)')"conjg(pft2)","real(pft1(i1,irot,ik)*conjg(pft2(i2,irot,ik)))"
     do i=1,ipart - (ipart-2)
        do irot=1,nrot - (nrot-2)
           do ik=1,nk - (nk-2)
              write(*,*)i, irot, ik, conjg(pft2(i,irot,ik)) , &
                   real( pft1(i,irot,ik)*conjg(pft2(i,irot,ik)) ), &
                   real(pft1(i,irot,ik)) * real(pft2(i,irot,ik)) + &
                   imag(pft1(i,irot,ik)) * imag(pft2(i,irot,ik)) 
           end do
        end do
     end do


     i1     = 1
     i2     = 1
     sumasq = 0.
     sumbsq = 0.
     r      = 0.

     open(1,file='AAstar_cpu_GNUf.log',status='unknown',action='write')

     do jpart=1,ipart
        do irot=1,nrot
           do ik=1,nk
              r = r+real(pft1(i1,irot,ik)*conjg(pft2(i2,irot,ik)))
              sumasq = sumasq+csq(pft1(i1,irot,ik))
              sumbsq = sumbsq+csq(pft2(i2,irot,ik))
!              write(1,'(i5,x,i5,x,i5,x,f20.8)')i1, irot, ik,csq(pft1(i1,irot,ik))
           end do
        end do
        i1 = i1+1
        if( i1 > ipart ) i1 = 1
        i2 = i2+1
        if( i2 > ipart ) i2 = 1
        !write(*,*)i1,i2,i,j,r
     end do

     close(1)


     write(*,*)"size(pft(npart,nrot,nk)) : ",ipart*nrot*nk/1.e6,"(Millions)"
     !write(*,*)"Full sum of real(pft1) is: ",sum(real(pft1))
     write(*,*)"Full sum of r is         : ",r
     write(*,*)"Full sum of sumasq is    : ",sumasq
     write(*,*)"Full sum of sumbsq is    : ",sumbsq
     write(*,*)"the correlator"
     r = calc_corr(r,sumasq*sumbsq)
     write(*,*)"after calc_corr r        : ",r
     write(*,*)"Elapsed time for corr_cpu: ",real(elps_corr_cpu),"(seconds)"

     write(*,*)'                                                               '
     write(*,*)'***************GPU corr****************************************'
     write(*,*)'                                                               '

     write(*,*)'                           X    N                              '

     call start_timer_cpu("XN_gpu")
     call gettimeofday_c(st_XN_r)
     err = get_polarft_corr_gpu_c_(s_polar,"X","N",            &
                                   sumX_vec,                   &
                                   pft1,pft2,                  &
                                   ipart,nrot,nk,              &
                                   lda,ldb,ldc,dble_alpha)
     call gettimeofday_c(et_XN_r)
     call elapsed_time_c(st_XN_r,et_XN_r,elps_corr_XN_gpu)
     call stop_timer_cpu("XN_gpu")

     write(*,*)"size(pft(npart,nrot,nk)) : ",ipart*nrot*nk/1.e6,"(Millions)"
     write(*,*)"Elapsed time for corr_gpu: ",real(elps_corr_XN_gpu),"(seconds)"

     deallocate(sumX_vec)




  end do



  write(*,*) "writting the sqsums_ref(:)"
  write(*,*) sqsums_refs
  write(*,*) "writting the sqsums_ptcls(:)"
  write(*,*) sqsums_ptcls

  !!        complex, allocatable :: pfts_ptcls_tmp(:,:,:)
!!        real, allocatable    :: sqsums_ptcls_tmp(:)

!!        complex, allocatable :: pfts_ptcls(:,:,:)

!!        real, allocatable    :: sqsums_ptcls(:)

!!        complex(sp), allocatable :: pfts_refs(:,:,:)

!!        integer :: alloc_stat ! shapearr(3), 


!        if( size(self%pfts_ptcls,dim=1) == 2*self%nptcls )then
            ! the arrays are already expanded, all we need to do is to copy data over
!            !$omp parallel workshare default(shared)
!            self%pfts_ptcls(self%nptcls+1:,:,:) = self%pfts_ptcls(:self%nptcls,:,:)
!            self%sqsums_ptcls(self%nptcls+1:)   = self%sqsums_ptcls(:self%nptcls)
!            !$omp end parallel workshare
!        else if( size(self%pfts_ptcls,dim=1) == self%nptcls )then
!!        allocate( pfts_ptcls_tmp(self%nptcls,self%ptclsz,self%nk)) !,& source=self%pfts_ptcls(:self%nptcls,:,:), stat=alloc_stat )

!!        pfts_ptcls_tmp = self%pfts_ptcls(:self%nptcls,:,:)
             
!        call alloc_err("In: gencorrs_all, 1", alloc_stat)

!!        allocate( sqsums_ptcls_tmp(self%nptcls))!
!!        sqsums_ptcls_tmp = self%sqsums_ptcls(:self%nptcls) !, stat=alloc_stat )

!        call alloc_err("In: gencorrs_all, 2", alloc_stat)
!        deallocate(self%pfts_ptcls)
!        deallocate(self%sqsums_ptcls)

!!        allocate( pfts_ptcls(2*self%nptcls,self%ptclsz,self%nk) ) !, stat=alloc_stat )
!        call alloc_err("In: gencorrs_all, 3", alloc_stat)

!!        allocate( sqsums_ptcls(2*self%nptcls) ) !, stat=alloc_stat )
!        call alloc_err("In: gencorrs_all, 4", alloc_stat)

!        allocate( pfts_refs(self%nptcls,self%refsz,self%nk) ) !,& source=self%pfts_refs, stat=alloc_stat )
!!        allocate( pfts_refs(self%nptcls,self%refsz,self%nk) ) !,& source=self%pfts_refs, stat=alloc_stat )

!        pfts_refs = 1.23456789 !self%pfts_refs
!!        pfts_refs = self%pfts_refs

        
!!        !$omp parallel workshare default(shared)
!!        pfts_ptcls(:self%nptcls,:,:)   = pfts_ptcls_tmp
!!        pfts_ptcls(self%nptcls+1:,:,:) = pfts_ptcls_tmp
!!        sqsums_ptcls(:self%nptcls)     = sqsums_ptcls_tmp
!!        sqsums_ptcls(self%nptcls+1:)   = sqsums_ptcls_tmp
!!        !$omp end parallel workshare

!        else
!            stop 'nonconforming dimension (1) of self%pfts_ptcls; expand_dim; simple_polarft_corrcalc'
!        endif
        
!                write(*,*) "#if defined (CUDA) && defined (MAGMA)", &
!                     d1lim(2),d1lim(1),d2lim(2),d2lim(1),self%nk, &
!                     d1lim(2)-d1lim(1),d2lim(2)-d2lim(1),self%nk
                ! (1) do the Hadamard product on the GPU

!                !$omp parallel default(shared) private(iptcl)
!                !$omp workshare
!                hadamard_prod = real( self%pfts_refs(:,:,:) * &
!                     conjg(self%pfts_ptcls(d1lim(1):d1lim(2),d2lim(1):d2lim(2),:)) )
!                !$omp end workshare

                ! (2) do the double sum on GPU or use the OpenMP construct below

                ! correlation normalization
!                !$omp do schedule(auto)
!                do iptcl=1,self%nptcls
!                    corrmat3d(iptcl,ptcl_ind,irot) = sum(hadamard_prod(iptcl,:,:))
!                end do
!                !$omp end do


                !pfts_refs = 1.1234567 !irot+self%refsz-1
!                do i=1,self%nptcls
!                   do j=1,self%nrots
!                      do k=1,self%nk
!                         write(3000,*)                                                    &
!                              i,j,k, real(pfts_refs(i,j,k)), aimag(pfts_refs(i,j,k)), &
!                              i,j,k, real(temp(i,j,k)), aimag(temp(i,j,k)),                     &
!                              i,j,k, real(conjg(temp(i,j,k))), aimag(conjg(temp(i,j,k)))
!                      end do
!                   end do
!                end do

 
!                do iptcl=1,self%nptcls
!                   write(*,*) sum(hadamard_prod(iptcl,:,:)),sumb_vec(iptcl)
!                end do

  

  do i=1,self%nptcls
     do j=1,self%nrots
        do k=1,self%nk
           write(3000,'(x,i,x,i,x,i,x,15.8f,x,15.8f,x,i,x,i,x,i,x,15.8f,x,15.8f,x,i,x,i,x,i,x,15.8f,x,15.8f)',
           i,j,k, real(A(i,j,k)), cuCimagf(A(i,j,k)),
           i,j,k, real(B(i,j,k)), cuCimagf(B(i,j,k)),
           i,j,k, real(conjg(B(i,j,k))), cuCimagf(conjg(B(i,j,k)))
        end do
     end do
  end do

  

     write(*,'(x,a,12x,a,12x,a)')"i","r1_vec","r2_vec"
     do i=1,ipart - (ipart-5)
        write(*,*)i, r1_vec(i), r2_vec(i)
     end do
     write(*,*) "the difference between sum(r1_vec-r2_vec): ",sum(r1_vec-r2_vec)

     do jpart=1,ipart-(ipart-16)
        write(*,'(x,i5,x,f15.8)') jpart,hadamard_prod_1D(jpart)
     end do

     do irot=1,nrot-(nrot-4)
        do ik=1,nk-(nk-4)
           write(*,'(x,i5,x,i5,f15.8)') irot, ik, hadamard_prod_23D(irot,ik)
        end do
     end do

     do jpart=1,ipart-(ipart-4)
        do irot=1,nrot-(nrot-4)
           do ik=1,nk-(nk-4)
              write(*,'(x,i5,x,i5,x,i5,f15.8)') jpart,irot, ik, hadamard_prod(jpart,irot,ik)
           end do
        end do
     end do






  real(sp), allocatable         :: pft1_cn_r(:,:)
  real(sp), allocatable         :: pft1_cn_i(:,:)

  real(sp), allocatable         :: pft1_cr_r(:,:)
  real(sp), allocatable         :: pft1_cr_i(:,:)

     allocate(pft1_cn_r(nrot,nk))
     allocate(pft1_cn_i(nrot,nk))

     allocate(pft1_cr_r(ipart,nk))
     allocate(pft1_cr_i(ipart,nk))



     write(*,*)'                                                               '
     write(*,*)'***************Cshift test on one of the PFT(1)****************'
     write(*,*)'                                                               '

     pft1_cn_r = 0.0
     pft1_cn_i = 0.0

     call start_timer_cpu("cshift_cpu")
     call gettimeofday_c(st_cshift_cpu)

     chunk = (ipart * nrot * nk ) / (24.0*10000.0)

!!$OMP parallel shared(pft1_cn_r,pft1_cn_i,pft1,chunk) private(jpart)
!!$OMP do schedule (dynamic ,chunk)
!     do jpart = 1,ipart
!        pft1_cn_r(:,:) = real(pft1(jpart,:,:))
!        pft1_cr_i(:,:) = imag(pft1(jpart,:,:))

!        pft1_cn_r(:,:) = pft1_cn_r(:,:) * pft1_cn_r(:,:) + pft1_cr_i(:,:) * pft1_cr_i(:,:) 
!        pft1_cn_i(:,:) = pft1_cn_r(:,:) * pft1_cn_r(:,:) - pft1_cr_i(:,:) * pft1_cr_i(:,:) 

!        pft1(jpart,:,:) = cmplx(pft1_cn_r,pft1_cn_i,sp)
!     end do

!!$OMP end do nowait
!!$OMP end parallel

!!$OMP parallel shared(pft1_cn_r,pft1_cn_i,pft1_cr_r,pft1_cr_i,pft1,chunk) private(jpart,irot)
!!$OMP do schedule (dynamic ,chunk)

     do jpart = 1,ipart
        pft1_cn_r(:,:) = cshift(real(pft1(jpart,:,:)),dim=1,shift=1)
        pft1_cn_i(:,:) = cshift(imag(pft1(jpart,:,:)),dim=1,shift=1)

        do irot=1,nrot
           pft1_cr_r(:,:) = cshift(real(pft1(:,irot,:)),dim=1,shift=1)
           pft1_cr_i(:,:) = cshift(imag(pft1(:,irot,:)),dim=1,shift=1)

           pft1(:,irot,:) = cmplx(pft1_cr_r,pft1_cr_i,sp)

        end do

        pft1(jpart,:,:) = cmplx(pft1_cn_r,pft1_cn_i,sp)
     end do
     
!!$OMP end do nowait
!!$OMP end parallel

     call gettimeofday_c(et_cshift_cpu)
     call elapsed_time_c(st_cshift_cpu,et_cshift_cpu,elps_cshift_cpu)
     call stop_timer_cpu("cshift_cpu")
     write(*,*)"Elapsed time for cshift_cpu: ",real(elps_cshift_cpu),"(seconds)"


     deallocate(pft1_cn_r)  
     deallocate(pft1_cn_i)

     deallocate(pft1_cr_r)
     deallocate(pft1_cr_i)







     do i=1,ipart
        do irot=1,nrot
           do ik=1,nk
              write(1000,'(1x,i3,1x,i3,1x,i3,1x,f15.8,1x,f15.8,1x,f15.8,1x,f15.8)')i,irot,ik, &
                   real(pft1(i,irot,ik)), imag(pft1(i,irot,ik)), &
                   real(pft2(i,irot,ik)), imag(pft2(i,irot,ik))
           end do
        end do
     end do





#if defined (LINUX)
     function get_polarft_corr_gpu(TRANSA,TRANSB,              &
                                   r_gpu,sumasq_gpu,sumbsq_gpu,&
                                   pft1,pft2,                  &
                                   npart,nrot,nk,              &
                                   lda,ldb,ldc,dble_alpha,dble_beta)
       character   :: TRANSA,TRANSB
       real(dp)    :: r_gpu
       real(dp)    :: sumasq_gpu, sumbsq_gpu
       complex(dp) :: pft1(npart,nrot,*)
       complex(dp) :: pft2(npart,nrot,*)
       integer     :: lda
       integer     :: ldb,ldc
       real(dp)    :: dble_alpha,dble_beta
       integer :: get_polarft_corr_gpu
     end function get_polarft_corr_gpu
#elif defined (MACOSX)
  interface
     function get_polarft_corr_gpu_c_(TRANSA,TRANSB,              &
                                      r_gpu,sumasq_gpu,sumbsq_gpu,&
                                      pft1,pft2,                  &
                                      npart,nrot,nk,              &
                                      lda,ldb,ldc,dble_alpha,dble_beta)
       use simple_defs
      character   :: TRANSA,TRANSB
       real(dp)    :: r_gpu
       real(dp)    :: sumasq_gpu, sumbsq_gpu
       complex(dp) :: pft1(npart,nrot,*)
       complex(dp) :: pft2(npart,nrot,*)
       integer     :: lda
       integer     :: ldb,ldc
       real(dp)    :: dble_alpha,dble_beta
       integer :: get_polarft_corr_gpu_c_
     end function get_polarft_corr_gpu_c_
#endif
  end interface






#if defined (MACOSX) && defined (CUDA)

#elif defined (LINUX) && defined (CUDA)
     err = get_polarft_corr_gpu("N","N",                       &
                                r_gpu,sumasq_gpu,sumbsq_gpu,   &
                                pft1,pft2,                     &
                                ipart,nrot,nk,                 &
                                lda,ldb,ldc,dble_alpha,dble_beta)
#endif


     call zz2dgemm_ElmtWs_tesla_gpu("N", "N", &
                                    inmx, inmx, inmx,  &
                                    dble_alpha,  &
                                    devPtrA_pft1, lda,  &
                                    devPtrA_pft2, ldb,  &
                                    dble_beta,  &
                                    devPtrA_D1, ldc)



  !timer variables

  double precision              :: elps_sumasq
  double precision,dimension(2) :: st_sumasq, et_sumasq

  double precision              :: elps_sumbsq
  double precision,dimension(2) :: st_sumbsq, et_sumbsq

  double precision              :: elps_r
  double precision,dimension(2) :: st_r, et_r
  double precision              :: elps_sum
  double precision,dimension(2) :: st_sum, et_sum




     allocate(pft1(inmx,inmx))
     allocate(pft2(inmx,inmx))


     call get1to9RowMajZdpMat_cpu(inmx,lda,pft1)
     call get1to9ColMajZdpMat_cpu(inmx,lda,pft2)
    



     write(*,'(34x,a,32x,a)')"(pft1) row major","(pft2) column major"
     do i=1,inmx - (inmx-2)
        do j=1,inmx - (inmx-2)
           write(*,*)i, j,pft1(i,j), pft2(i,j)
        end do
     end do
     write(*,'(34x,a)')"conjg(pft2)"
     do i=1,inmx - (inmx-2)
        do j=1,inmx - (inmx-2)
           write(*,*)i, j,conjg(pft2(i,j))
        end do
     end do




     i1     = rot1
     i2     = rot2
     sumasq = 0.
     sumbsq = 0.
     r      = 0.
     call start_timer_cpu("corr_cpu")
     write(*,*)"the correlator"
     !  do i=1,nradial1/2
     do i=1,nradial1
        do j=khp,klp
           r = r+real(pft1(i1,j)*conjg(pft2(i2,j)))
           sumasq = sumasq+csq_dble(pft1(i1,j))
           sumbsq = sumbsq+csq_dble(pft2(i2,j))
           !write(*,*)i, j, r
        end do
        i1 = i1+1
        if( i1 > nradial1 ) i1 = 1
        i2 = i2+1
        if( i2 > nradial2 ) i2 = 1
        !write(*,*)i1,i2,i,j,r
     end do
     call stop_timer_cpu("corr_cpu")

     !write(*,'(1x,a,1x,2f20.8)')"Full sum of the Double  precision matrix r is: ",r
     r = calc_corr_dble(r,sumasq*sumbsq)
     write(*,*)"after calc_corr_dble r= ",r!," sumasq= ",sumasq," sumbsq= ",sumbsq









     !Now calluating the r value
     allocate(DA_gpu(inmx,inmx))
     write(*,*) 'before the D1 alloc'
     write(*,*) 'after the D1 alloc'
     write(*,*) 'before the D1 set'

     call start_timer_cpu("r_gpu")
     call gettimeofday_c(st_r)
     write(*,*)"r_gpu"
     call gettimeofday_c(et_r)
     call elapsed_time_c(st_r,et_r,elps_r)
     write(*,'(x,a,x,f20.8)') "the elapse time for the tesla_r(s)= ",elps_r
     call stop_timer_cpu("r_gpu")

     DA_gpu = 0.0

     r_gpu = sum(DA_gpu)
     !freein the ressources on device for the first r calculation.
     deallocate(DA_gpu)

     !***********Now calculate the sumsq(a and b)*********************
     !******sumasq*****
     allocate(DA_gpu(inmx,inmx))

     call gettimeofday_c(st_sumasq)
     write(*,*)"sumsq"
     call gettimeofday_c(et_sumasq)
     call elapsed_time_c(st_sumasq,et_sumasq,elps_sumasq)
     write(*,'(x,a,x,f20.8)') "the elapse time for the tesla_sumasq(s)= ",elps_sumasq

     sumasq_gpu = 0.0d0
     sumasq_gpu = sum(DA_gpu)
     err = cublas_free(devPtrA_D2)
     deallocate(DA_gpu)

     !******sumbsq*****
     allocate(DA_gpu(inmx,inmx))

     call gettimeofday_c(st_sumbsq)
     write (*,*)"hello"
     call gettimeofday_c(et_sumbsq)
     call elapsed_time_c(st_sumbsq,et_sumbsq,elps_sumbsq)
     write(*,'(x,a,x,f20.8)') "the elapse time for the tesla_sumbsq(s)= ",elps_sumbsq

     DA_gpu = 0.0

     sumbsq_gpu = 0.0d0
     sumbsq_gpu = sum(DA_gpu)
     err = cublas_free(devPtrA_D3)
     deallocate(DA_gpu)

     !******Now calculating the correlator*****
     r = 0.0d0
     r = calc_corr_dble(r_gpu,sumasq_gpu*sumbsq_gpu)
     write(*,*)"after calc_corr_dble r= ",r!," sumasq= ",sumasq," sumbsq= ",sumbsq






    if ( present(nx) ) then
       call c_devFFT_cpu%new_fft_1D_DZ_cpu(nx)
    else if ( present(nx) .and. present(ny) ) then
       call c_devFFT_cpu%new_fft_2D_DZ_cpu(nx,ny)
    else if ( present(nx) .and. present(ny) .and. present(nz) ) then
       call c_devFFT_cpu%new_fft_3D_DZ_cpu(nx,ny,nz)
    end if

  ui = cmplx(-10.0d0,-10.0d0 )
  uf = cmplx( 10.0d0, 10.0d0 )

  allocate(fxy(nx,ny))
  allocate(x(nx))
  allocate(y(ny))

  x = 0.0d0
  y = 0.0d0
  fxy = 0.0d0

  call atom_like_2D_double(ui,uf,nx,ny,a,b,c,x,y,fxy,n)

  deallocate(x)
  deallocate(y)
  deallocate(fxy)


  do i=1,nx
     do j=1,ny
        write(2000,*)real(u(i)),real(v(j)),real(fuv(i,j))
        write(3000,*)imag(u(i)),imag(v(j)),imag(fuv(i,j))
     end do
  end do
  


  write(*,*)i,j,u(i),v(j),fxy(i,j),fhp(i,j)!,fuv_gpu(i,j)/(nx*ny)
  write(*,'(2x,i1,2x,i1,3(x,f25.15),x,a,(x,f25.15,f25.15),x,a,x,f25.15)')&
       i,j,x(i),y(j),fxy(i,j),"(",fhp(i,j),")",fxy_gpu(i,j)/(nx*ny)

    if ( nx >= 5 ) then
       write(*,*) "the first 5 entries of:"
       write(*,*) "data:f(u,v)     Fourier: f(h,p)"
       do i=1,nx-(nx-3)
          do j=1,ny-(ny-3)
             write(*,*)fuv(i,j), fhp(i,j)/(nx*ny)
          end do
       end do
    end if

  do j=1,30
     do i=1,nx
        write(2000+j,*)x(i), j * fx(i)
     end do
  end do


    interface external_c_function_gather_fft_2D_gpu_c
       function gather_fft_2D_gpu_c(nx, ny,fuv,fhp,sign)
         use simple_defs
         integer  :: nx,ny
         integer  :: sign
         complex(dp) :: fuv(nx,*)
         complex(dp) :: fhp(nx,*)
         integer :: gather_fft_2D_gpu_c
       end function gather_fft_2D_gpu_c
    end interface external_c_function_gather_fft_2D_gpu_c


#define devptr_t integer*8
    devptr_t              :: devPtrA
    integer                         :: size_of_elt
    
    rc = cublas_get_matrix ( nx, ny, size_of_elt, devPtrA, ny, fhp(1:nx,1:ny), ny )

    rc = simple_cublas_free(devPtrA)

    interface external_c_function_gather_fft_2D_gpu_c
       function gather_fft_2D_gpu_c(nx, ny, devPtrA, fhp, sign)
         use simple_defs
         integer  :: nx,ny
         integer  :: sign
         devptr_t              :: devPtrA
         complex(dp) :: fhp(nx,*)
         integer :: gather_fft_2D_gpu_c
       end function gather_fft_2D_gpu_c
    end interface external_c_function_gather_fft_2D_gpu_c



    size_of_elt = 2 * size_of_double_complex
    write(*,*) size_of_elt, size_of_double_complex

    rc = cublas_alloc(nx*ny, size_of_elt, devPtrA)
!    if (rc .ne. 0 ) call simple_cudblas_stat_return(rc)

    rc = cublas_set_matrix (nx, ny, size_of_elt, fuv(1:nx,1:ny), nx, devPtrA, ny )
!    if ( rc .ne. 0 ) call simple_cudblas_stat_return(rc)

    !fhp(1:nx,1:ny)=0.0d0 !initializing the inverted matrix


    rc = gather_fft_2D_gpu_c(nx, ny, devPtrA, fhp, sign)




    if ( nx >= 5 ) then
       write(*,*) "the first 5 entries of:"
       write(*,*) "data:f(u,v)     Fourier: f(h,p)"
       do i=1,nx-(nx-3)
          do j=1,ny-(ny-3)
             write(*,*)fuv(i,j), fhp(i,j)
          end do
       end do
    end if



  if ( nx >= 5 ) then
     write(*,*) "the first 3x(i,j) = 9 entries of:"
     write(*,'(2x,a,5(10x,a))')"u","v","w","data:f(u,v,w)","Fourier: f(h,p,q)","Inverse Fourier"
     do i=1,nx-(nx-3)
        do j=1,ny-(ny-3)
           do k=1,nz-(nz-3)
              write(*,*)i,j,k,u(i),v(j),w(k),fuvw(i,j,k),fhpq(i,j,k),fuvw_gpu(i,j,k)/(nx*ny*nz)
           end do
        end do
     end do
  end if



  if ( nx >= 5 ) then
     write(*,*) "the first 5 entries of:"
     write(*,*) "u         data:f(u)     Fourier: f(h)      Inverse Fourier"
     do i=1,nx-(nx-5)
        write(*,*)u(i), fu(i), fh(i), fu_gpu(i)/nx
     end do
  end if

  do i=1,nx
     do j=1,ny
        write(2000,*)real(u(i)),real(v(j)),real(fuv(i,j))
        write(3000,*)imag(u(i)),imag(v(j)),imag(fuv(i,j))
     end do
  end do

  do i=1,nx
     do j=1,ny
        write(2000,*) x(i),y(j),fxy(i,j)
     end do
  end do


  complex(dp),allocatable :: fh_cpu(:)
  complex(dp),allocatable :: fz_cpu(:)
  complex(dp),allocatable :: fh_gpu(:)
  complex(dp),allocatable :: fz_gpu(:)

  allocate(fh_cpu(nstep))
  allocate(fz_cpu(nstep))

  allocate(fh_gpu(nstep))
  allocate(fz_gpu(nstep))

  fh_gpu = 0.0d0
  call t_cuFFT_gpu%gather_fft_1D_gpu(nstep,fz,fh_cpu,FFTW_FORWARD)

  do i=1,nstep
     write(2000,*)z(i), fz(i), fh_cpu(i), fz_cpu(i)/nstep
  end do

  deallocate(fh_cpu)
  deallocate(fz_cpu)
  deallocate(fh_gpu)
  deallocate(fz_gpu)

!
!---------------------Start of the Removed from the testing codes---------------
!
  interface
     function atom_like_real(xi,xf,nstep,a,b,c,x,n) result(fx)
       use simple_defs
       real(dp) :: xi,xf
       integer  :: nstep
       real(dp) :: a,b,c,n
       real(dp),dimension(nstep) :: x,fx
     end function atom_like_real
     function atom_like_complex(zi,zf,nstep,a,b,c,z,n) result(fz)
       use simple_defs
       integer  :: nstep
       real(dp) :: a,b,c,n
       complex(dp) :: zi,zf
       complex(dp),dimension(nstep) :: z,fz
     end function atom_like_complex
  end interface

!Generating fuinctional for testign purposes.
function atom_like_complex(zi,zf,nstep,a,b,c,z,n) result(fz)
  use simple_defs
  integer  :: nstep
  real(dp) :: a,b,c,n
  complex(dp) :: zi,zf
  complex(dp),dimension(nstep) :: z,fz
  !local variables
  real(dp) :: re_delta,im_delta
  real(dp) :: re_z,im_z
  complex(dp) :: z_i
  !counters
  integer :: i

  re_delta = ( real(zf) - real(zi) ) / (nstep-1)
  im_delta = ( imag(zf) - imag(zi) ) / (nstep-1)
  z_i = zi
  do i=1,nstep
     re_z = real(zi) + (i-1) * re_delta
     im_z = imag(zi) + (i-1) * im_delta
     z(i) = cmplx(re_z , im_z, dp)
  end do

  fz(:) = cmplx( &
       cos(real(z(:))**b) / ( (real(z(:))**a) + c )**n, &
       cos(imag(z(:))**b) / ( (imag(z(:))**a) + c )**n, dp)

end function atom_like_complex

!Generating fuinctional for testign purposes.
function atom_like_real(xi,xf,nstep,a,b,c,x,n) result(fx)
  use simple_defs
  integer  :: nstep
  real(dp) :: xi,xf
  real(dp) :: a,b,c,n
  real(dp) :: delta
  real(dp),dimension(nstep) :: x,fx
  real(dp) :: x_i
  !counters
  integer :: i

  delta = ( xf - xi ) / (nstep-1)
  x_i = xi
  do i=1,nstep
     x(i) = xi + (i-1) * delta
  end do

  fx(:) = cos(x(:)**b) / ( (x(:)**a) + c )**n

end function atom_like_real


!
!---------------------End of the Removed from the testing codes-----------------
!



!-------------------------------------------------------------------------------

!Generating fuinctional for testign purposes.
function atom_like_complex(zi,zf,nstep,a,b,c,z,n) result(fz)
  use simple_defs
  integer  :: nstep
  real(dp) :: a,b,c,n
  complex(dp) :: zi,zf
  complex(dp),allocatable :: z(:),fz(:)
  !local variables
  real(dp) :: re_delta,im_delta
  real(dp) :: re_z,im_z
  complex(dp) :: z_i
  !counters
  integer :: i

  !allocating the resources 
  allocate(fz(nstep))

  re_delta = ( real(zf) - real(zi) ) / (nstep-1)
  im_delta = ( imag(zf) - imag(zi) ) / (nstep-1)
  z_i = zi
  do i=1,nstep
     re_z = real(zi) + (i-1) * re_delta
     im_z = imag(zi) + (i-1) * im_delta
     z(i) = cmplx(re_z , im_z, dp)
  end do

  fz(:) = cmplx( &
       cos(real(z(:))**b) / ( (real(z(:))**a) + c )**n, &
       cos(imag(z(:))**b) / ( (imag(z(:))**a) + c )**n, &
       dp)

end function atom_like_complex

!Generating fuinctional for testign purposes.
function atom_like_real(xi,xf,nstep,a,b,c,x,n) result(fx)
  use simple_defs
  integer  :: nstep
  real(dp) :: xi,xf
  real(dp) :: a,b,c,n
  real(dp) :: delta
  real(dp),allocatable :: x(:),fx(:)
  real(dp) :: x_i
  !counters
  integer :: i

  !allocating the resources 
  allocate(fx(nstep))

  delta = ( xf - xi ) / (nstep-1)
  x_i = xi
  do i=1,nstep
     x(i) = xi + (i-1) * delta
  end do

  fx(:) = cos(x(:)**b) / ( (x(:)**a) + c )**n

end function atom_like_real

!-------------------------------------------------------------------------------


!Generating fuinctional for testign purposes.
function atom_like_complex(zi,zf,nstep,a,b,c,z,n) result(fz)
  use simple_defs
  integer  :: nstep
  real(dp) :: a,b,c,n
  complex(dp) :: zi,zf
  complex(dp),allocatable :: z(:),fz(:)
  !local variables
  real(dp) :: re_delta,im_delta
  real(dp) :: re_z,im_z
  complex(dp) :: z_i
  !counters
  integer :: i

  !allocating the resources 
  allocate(fz(nstep))

  re_delta = ( real(zf) - real(zi) ) / (nstep-1)
  im_delta = ( imag(zf) - imag(zi) ) / (nstep-1)
  z_i = zi
  do i=1,nstep
     re_z = real(zi) + (i-1) * re_delta
     im_z = imag(zi) + (i-1) * im_delta
     z(i) = cmplx(re_z , im_z, dp)
  end do

  fz(:) = cmplx( &
       cos(real(z(:))**b) / ( (real(z(:))**a) + c )**n, &
       cos(imag(z(:))**b) / ( (imag(z(:))**a) + c )**n, &
       dp)

end function atom_like_complex

!Generating fuinctional for testign purposes.
function atom_like_real(xi,xf,nstep,a,b,c,x,n) result(fx)
  use simple_defs
  integer  :: nstep
  real(dp) :: xi,xf
  real(dp) :: a,b,c,n
  real(dp) :: delta
  real(dp),allocatable :: x(:),fx(:)
  real(dp) :: x_i
  !counters
  integer :: i

  !allocating the resources 
  allocate(fx(nstep))

  delta = ( xf - xi ) / (nstep-1)
  x_i = xi
  do i=1,nstep
     x(i) = xi + (i-1) * delta
  end do

  fx(:) = cos(x(:)**b) / ( (x(:)**a) + c )**n

end function atom_like_real


  use simple_fftw3


   !CPU
   function gather_fft_1D_cpu_cpp(nstep,fz,fh,sign)
     use simple_defs
     integer  :: nstep
     integer  :: sign
     complex(dp) :: fz(*)
     complex(dp) :: fh(*)
     integer :: gather_fft_1D_cpu_cpp
   end function gather_fft_1D_cpu_cpp


   procedure :: gather_fft_1D_cpu



  !getting the fourier transform on CPU for 1D data
  subroutine gather_fft_1D_cpu(c_devFFT_gpu,nstep,fz,fh,sign)
    use simple_defs
    class(cuFFT_gpu), intent(inout) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer  :: nstep,sign
    complex(dp) :: fz(*)
    complex(dp) :: fh(:)
    integer :: rc !return code
    interface external_c_function_gather_fft_1D_cpu_c
       function gather_fft_1D_cpu_c(nstep,fz,fh,sign)
         use simple_defs
         integer  :: nstep
         integer  :: sign
         complex(dp) :: fz(*)
         complex(dp) :: fh(:)
         integer :: gather_fft_1D_cpu_c
       end function gather_fft_1D_cpu_c
    end interface external_c_function_gather_fft_1D_cpu_c

#if defined (MACOSX) && defined (CUDA)
    rc = gather_fft_1D_cpu_c(nstep,fz,fh,sign)
#elif defined (LINUX) && defined (CUDA)
    rc = gather_fft_1D_cpu_cpp(nstep,fz,fh,sign)
#else
    call c_devFFT_gpu%terminate()
#endif

    return
  end subroutine gather_fft_1D_cpu



testfunc


  write(*,*) zi,zf

  write(*,*) FFTW_FORWARD, FFTW_BACKWARD 

  do i=1,nstep
     write(2000,*)z(i), fz(i)
  end do

  do i=1,nstep
     write(2000,*)z(i), fz(i), fh_cpu(i), fz_cpu(i)/nstep
  end do

  subroutine get_fft_1D_cpu(xi,xf,nstep,x,fx,fh)
    use simple_defs
    use simple_fftw3
    integer  :: nstep
    real(dp) :: xi,xf
    real(dp),allocatable :: x(:),fx(:)
    complex(dp),allocatable :: fh(:)
  end subroutine get_fft_1D_cpu

!getting the fourier transform on CPU for comparision checks with GPU
subroutine get_fft_1D_cpu(xi,xf,nstep,x,fx,fh)
  use simple_defs
  use simple_fftw3
  integer  :: nstep
  real(dp) :: xi,xf
  real(dp),allocatable :: x(:),fx(:)
  complex(dp),allocatable :: fh(:)
  !fast fourier transform
  integer*8 :: nstep_c
  type(c_ptr) :: p  !c pointer to the fast fourier fftw allocation
  type(c_ptr) :: plan_fwd
  integer(c_int) :: flags
  real(c_float),pointer :: in(:)=>null()
  complex(c_float_complex),pointer :: out(:)=>null()

  !converting to integer 8 for c_binding 
  nstep_c = nstep
  !Allocating the memory
  p = fftwf_alloc_real(nstep_c)
!  in => real(fx)
!  out => real(fh)
!  call c_f_pointer(p,out,[nstep])
  plan_fwd = fftwf_plan_dft_r2c_1d(nstep,in,out,fftw_estimate)

  !write()
  
  !freeing the ressources
  call fftwf_free(p)

  return
end subroutine get_fft_1D_cpu



  integer :: i
  integer :: j

  do j=1,30
     do i=1,nstep
        write(2000+j,*)x(i), j * fx(i), j*fh(i)
     end do
  end do

  integer :: i
  integer :: j

  do j=1,30
     do i=1,nstep
        write(2000+j,*)x(i), j * fx(i)
     end do
  end do

  integer :: i
  integer :: j

  do j=1,300
     do i=1,nstep
        write(2000+j,*)c_devFFT_gpu%x(i), j * c_devFFT_gpu%fx(i)
     end do
  end do

  do i=1,nstep
     write(1000,*)x(i), fx(i)
  end do

!*******************************************************************************
!     now testign the corr function 
!
!*******************************************************************************

     write(*,*)'                                                               '
     write(*,*)'***************************************************************'
     write(*,*)'     now testign the corr function                             '
     write(*,*)'***************************************************************'

     write(*,*)'                                                               '
     write(*,*)'***************CPU corr****************************************'
     write(*,*)'                                                               '

     write(*,*)'                                                               '
     write(*,*)'***************GPU corr****************************************'
     write(*,*)'                                                               '

     !***********Now calculate the sumsq(a and b)*********************
     !******sumasq*****
     !******sumbsq*****
     !******Now calculating the correlator*****

     write(*,*)'                                                               '
     write(*,*)'***************************************************************'
     write(*,*)'                                                               '


     call gettimeofday_c(st_sum)
     r_gpu = sum(DA_gpu)
     call gettimeofday_c(et_sum)
     call elapsed_time_c(st_sum,et_sum,elps_sum)
     write(*,'(x,a,x,f20.8)') "the elapse time for the sum r(s)= ",elps_sum

     call gettimeofday_c(st_sum)
     sumasq_gpu = sum(DA_gpu)
     call gettimeofday_c(et_sum)
     call elapsed_time_c(st_sum,et_sum,elps_sum)
     write(*,'(x,a,x,f20.8)') "the elapse time for the sum r(s)= ",elps_sum

     call gettimeofday_c(st_sum)
     sumbsq_gpu = sum(DA_gpu)
     call gettimeofday_c(et_sum)
     call elapsed_time_c(st_sum,et_sum,elps_sum)
     write(*,'(x,a,x,f20.8)') "the elapse time for the sum bsq(s)= ",elps_sum



     !checking the matrix DA_gpu after CUDA experience
     write(*,*)
     write(*,*) "The DA_gpu(ilda,jlda))"
     write(*,*)
     do ilda = 1,LDA - (inmx-2)
        do jlda = 1,LDA - (inmx-2)
           write(*,*) ilda, jlda, DA_gpu(ilda,jlda)
        end do
     end do

  !pft device pointers to be deleted after the test
  devptr_t                      ::  devPtrA_pft11
  devptr_t                      ::  devPtrA_pft12

  devptr_t                      ::  devPtrA_pft21
  devptr_t                      ::  devPtrA_pft22




     err = cublas_free(devPtrA_pft1)
     err = cublas_free(devPtrA_pft2)



     err = cublas_alloc(inmx*inmx, size_of_double_complex, devPtrA_pft21)
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)
     err = cublas_alloc(inmx*inmx, size_of_double_complex, devPtrA_pft22)
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)

     !setting up the pft 1 and 2 matrix on device
     err = cublas_set_matrix (inmx, inmx, size_of_double_complex, pft1, inmx, devPtrA_pft21, inmx )
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)
     err = cublas_set_matrix (inmx, inmx, size_of_double_complex, pft2, inmx, devPtrA_pft22, inmx )
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)

     err = cublas_free(devPtrA_pft21)
     err = cublas_free(devPtrA_pft22)


     err = cublas_alloc(inmx*inmx, size_of_double_complex, devPtrA_pft11)
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)
     err = cublas_alloc(inmx*inmx, size_of_double_complex, devPtrA_pft12)
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)

     !setting up the pft 1 and 2 matrix on device
     err = cublas_set_matrix (inmx, inmx, size_of_double_complex, pft1, inmx, devPtrA_pft11, inmx )
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)
     err = cublas_set_matrix (inmx, inmx, size_of_double_complex, pft2, inmx, devPtrA_pft12, inmx )
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)

     err = cublas_free(devPtrA_pft11)
     err = cublas_free(devPtrA_pft12)



!*****************first part of the testing_corr_gpu.f90************************

     allocate(matSymZE_gpu(inmx,inmx))
     allocate(invSymZA(inmx,inmx))

     matSymZE_gpu = 0.0d0

     allocate(matSymZF_gpu(inmx,inmx))

     call random_seed(SIZE=k)       !initializing the seed of the random generator

     call getSymDRandomMat_cpu(inmx,lda,matSymZE_gpu)

     !checking the matrix matSymZE_gpu after initialisation
     write(*,*)
     write(*,*) "The product of matrix matSymZE_gpu"
     write(*,*)
     do ilda = 1,LDA - (inmx-2)
        do jlda = 1,LDA - (inmx-2)
           write(*,'(x,i4,x,i4,4x,f16.8,4x,f16.8)') ilda, jlda, matSymZE_gpu(ilda,jlda)
        end do
     end do

    ! now using the magma to get iunverse matrix

     size_of_elt = kind(matSymZE_gpu)
      
     err = cublas_alloc(inmx*inmx, size_of_elt, devPtrA_new)
     if (err .ne. 0 ) call simple_cudblas_stat_return(err)

     err = cublas_set_matrix (inmx, inmx, size_of_elt ,matSymZE_gpu, inmx, devPtrA_new, inmx )
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)

     invSymZA = 0.0d0 !initializing the inverted matrix
     call getBlockInverseDMatrix_gpu(inmx,inmx,lda,size_of_elt,inmx,devPtrA_new,invSymZA)
 
     !checking the matrix matSymZE_gpu after CUDA experience
     write(*,*)
     write(*,*) "The Real(matSymZE_gpu(ilda,jlda))," , ",aimag(matSymZE_gpu(ilda,jlda))"
     write(*,*)
     do ilda = 1,LDA - (inmx-2)
        do jlda = 1,LDA - (inmx-2)
           write(*,*) ilda, jlda, matSymZE_gpu(ilda,jlda)
        end do
     end do
     !freein memory on the device before doing the matrix multiplication test
     err = cublas_free(devPtrA_new)
     
     !matSymZF_gpu = matmul(matSymZE_gpu,invSymZA)   !this is the CPU slow way
     !this is the cublas fast way for matrix multiplication
     alpha = (1.0d0,0.0d0)
     beta = (0.0d0,0.0d0)
     ldb = lda
     ldc = lda
     matSymZF_gpu = 0.0d0
     !allocating memory
     err = cublas_alloc(inmx*inmx, size_of_elt, devPtrA_invSymZA)
     err = cublas_alloc(inmx*inmx, size_of_elt, devPtrA_SymZE_gpu)
     err = cublas_alloc(inmx*inmx, size_of_elt, devPtrA_SymZF_gpu)
     !setting up the matrix on device
     err = cublas_set_matrix (inmx, inmx, size_of_elt ,invSymZA, inmx, devPtrA_invSymZA, inmx )
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)
     err = cublas_set_matrix (inmx, inmx, size_of_elt ,matSymZE_gpu, inmx, devPtrA_SymZE_gpu, inmx )
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)
     !performing the matrix multiplication
!     call start_timer_cpu("cublas_ZGEMM")
     call cublas_DGEMM("N", "N", inmx, inmx, inmx, alpha, devPtrA_SymZE_gpu, &
          lda, devPtrA_invSymZA, ldb, beta, devPtrA_SymZF_gpu, ldc)
!     call stop_timer_cpu("cublas_ZGEMM")
     err = cublas_get_matrix ( inmx, inmx, size_of_elt, devPtrA_SymZF_gpu, inmx, matSymZF_gpu, inmx )

     write(*,'(1x,a,1x,2f20.8)')"Full sum of the Double Complex precision matrix matSymZF_gpu is: ",sum(matSymZF_gpu)

     !freeing the ressources on GPU
     err = cublas_free(devPtrA_invSymZA)
     err = cublas_free(devPtrA_SymZE_gpu)
     err = cublas_free(devPtrA_SymZF_gpu)

!****************************************

     DA_gpu(1:inmx,1:inmx) = Real(invSymZA(1:inmx,1:inmx))


     deallocate(matSymZE_gpu)
     deallocate(matSymZF_gpu)
     deallocate(invSymZA)

!*******************************************************************************

  real(sp)                        :: devPtr_kfromto
  real(sp)                        :: devPtr_coords
  real(sp)                        :: devptr_angtab

  devPtr_kfromto = 5.0
  devPtr_coords = 6.0
  devPtr_angtab = 7.0
    
  write(*,*) 'The size of the matrix exceeds n: ',inmx
        write(*,*) 'Performing a vector decomposition for checking inversion'

        !disecting matrix invSymZA and matSymZE_gpu into the real and imaginary part
        !into vectors

        !allocating the ressources for the computation
        allocate(vecZEr(inmx))
        allocate(vecZEi(inmx))
        allocate(vecZAr(inmx))
        allocate(vecZAi(inmx))

        vecZEr(1:inmx) = Real(matSymZE_gpu(1,1:inmx))
        vecZEi(1:inmx) = aimag(matSymZE_gpu(1,1:inmx))

        vecZAr(1:inmx) = Real(invSymZA(1:inmx,1))
        vecZAi(1:inmx) = aimag(invSymZA(1:inmx,1))

        !allocating space on the device
        err = cublas_alloc(inmx, size_of_double, devPtrA_vecZEr)
        if (err .ne. 0 ) call simple_cudblas_stat_return(err)
        err = cublas_alloc(inmx, size_of_double, devPtrA_vecZEi)
        if (err .ne. 0 ) call simple_cudblas_stat_return(err)
        err = cublas_alloc(inmx, size_of_double, devPtrA_vecZAr)
        if (err .ne. 0 ) call simple_cudblas_stat_return(err)
        err = cublas_alloc(inmx, size_of_double, devPtrA_vecZAi)
        if (err .ne. 0 ) call simple_cudblas_stat_return(err)

        !setting the vectors on device
        err = cublas_set_vector(inmx, size_of_double ,vecZEr, 1, devPtrA_vecZEr, 1)
        if (err .ne. 0 ) call simple_cudblas_stat_return(err)
        err = cublas_set_vector(inmx, size_of_double ,vecZEi, 1, devPtrA_vecZEi, 1)
        if (err .ne. 0 ) call simple_cudblas_stat_return(err)
        err = cublas_set_vector(inmx, size_of_double ,vecZAr, 1, devPtrA_vecZAr, 1)
        if (err .ne. 0 ) call simple_cudblas_stat_return(err)
        err = cublas_set_vector(inmx, size_of_double ,vecZAi, 1, devPtrA_vecZAi, 1)
        if (err .ne. 0 ) call simple_cudblas_stat_return(err)

        err = cublas_alloc(1, size_of_double, devPtrA_result_mat11)

        cublas_ddot(inmx, devPtrA_vecZEr, 1, devPtrA_vecZAr, 1)

        !Real(result_mat(1,1)) = cublas_ddot()


        !dealocating the ressources on device
        err = simple_cublas_free(devPtrA_vecZEr)
        err = simple_cublas_free(devPtrA_vecZEi)
        err = simple_cublas_free(devPtrA_vecZAr)
        err = simple_cublas_free(devPtrA_vecZAi)

        !dealocating the ressources on host
        deallocate(vecZEr)
        deallocate(vecZEi)
        deallocate(vecZAr)
        deallocate(vecZAi)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to fill a complex 2d array of size (n)x(m) with random numbers.
! using the random generator from the GPU magma routines.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! SOURCE
  subroutine getBlockInverseSMatrix_gpu(n,m,lda,size_of_elt,lwork,devPtrA,invSA_gpu)
    implicit none

    !global variables

    integer                         :: n,m
    integer                         :: lda
    integer                         :: size_of_elt
    integer                         :: lwork

    real(sp)                        :: invSA_gpu(LDA,*)

    devptr_t                        :: devPtrA

    !local variable

    integer                         :: err
    integer,allocatable             :: ipiv(:)
    real(sp),allocatable            :: work(:)

    !timer variables

    double precision                :: elps_magma_block_sgetrf, elps_magma_block_sgetri
    double precision,dimension(2)   :: st_magma_block_sgetrf, et_magma_block_sgetrf
    double precision,dimension(2)   :: st_magma_block_sgetri, et_magma_block_sgetri

    !start of the execution commands
    
#if defined (MAGMA)

    write(*,*)" hello subroutine getBlockInverseSMatrix_gpu"

    allocate(ipiv(n))
    allocate(work(lwork))

    call magmaf_init()

    call gettimeofday_c(st_magma_block_sgetrf)
    call start_timer_cpu("magma_sgetrf_gpu")
!    call magmaf_sgetrf_gpu(n, m, devPtrA, n, ipiv, err ) !The magma call with the f control
    call magma_sgetrf_gpu(n, m, devPtrA, n, ipiv, err )  !my wrapper call to magma
    call stop_timer_cpu("magma_sgetrf_gpu")
    if (err .ne. 0) write(*,*) "the LU decomposition in magma_sgetrf_gpu has failed"
    call gettimeofday_c(et_magma_block_sgetrf)
    call elapsed_time_c(st_magma_block_sgetrf,et_magma_block_sgetrf,elps_magma_block_sgetrf)


    call gettimeofday_c(st_magma_block_sgetri)
    call start_timer_cpu("magma_sgetri_gpu")
!    call magma_sgetri_block_gpu_v2( n, devPtrA, lda, ipiv, work, lwork, err )
    call stop_timer_cpu("magma_sgetri_gpu")
    call gettimeofday_c(et_magma_block_sgetri)
    call elapsed_time_c(st_magma_block_sgetri,et_magma_block_sgetri,elps_magma_block_sgetri)

    err = cublas_get_matrix ( n, m, size_of_elt, devPtrA, m, invSA_gpu, m )

    deallocate(ipiv)
    deallocate(work)

    !benchmark output.
#ifdef BENCH
    write(2500,'(i5,a,f10.5,a,f15.8,a,f15.8,a,f15.8,a,f15.8,a,f15.8,a,f15.8)') &
         n, ", ", &
         log10(dble(n)), ", ", &
         elps_magma_block_sgetrf, ", ", &
         elps_magma_block_sgetri, ", ", &
         (elps_magma_block_sgetrf + elps_magma_block_sgetri), ", ", &
         Log10(elps_magma_block_sgetrf), ", ", &
         Log10(elps_magma_block_sgetri), ", ", &
         Log10(elps_magma_block_sgetrf+elps_magma_block_sgetri)
#endif
    !end of benchmark output

    !finalising the magma environment
    call magmaf_finalize()

    write(*,*)" by subroutine getBlockInverseSMatrix_gpu"

#endif

    return
  end subroutine getBlockInverseSMatrix_gpu






    complex(dp),allocatable         :: matZA_gpu(:,:)

    devptr_t                        :: devPtrA_in


    allocate(matZA_gpu(n,m))

    matZA_gpu(1:n,1:m) = invZA_gpu(1:n,1:m)

    size_of_elt = 2 * kind(matZA_gpu)
    
    err = cublas_alloc(n*m, size_of_elt, devPtrA_in)
    if (err .ne. 0 ) call simple_cudblas_stat_return(err)
    
    err = cublas_set_matrix (n, m, size_of_elt ,matZA_gpu, m, devPtrA_in, m )
    if ( err .ne. 0 ) call simple_cudblas_stat_return(err)

!  integer,external :: cublas_get_error      !function to get the error from the GPU
!  integer,external :: cublas_shutdown       !shutting down the CUBLAS device
!  integer,external :: cublas_set_matrix     !setting the matrix on the GPU
!  integer,external :: cublas_get_matrix     !getting the matrix from the GPU
!  integer,external :: cublas_set_vector     !setting the matrix on the GPU
!  integer,external :: cublas_get_vector     !setting the matrix on the GPU


!  integer,external :: cublas_init           !intializing the memory for the GPU
!  integer,external :: cublas_alloc          !allocatign the memory for the GPU
!  integer,external :: cublas_free           !freeing the memory from the GPU


  write(2000,'(1x,i5,2x,f10.5,2x,f15.8,2x,f15.8,10x,f15.8,10x,f15.8,12x,f15.8,12x,f15.8)') &
       n, &
       log10(dble(n)), &
       elps_magma_block_zgetrf, &
       elps_magma_block_zgetri, &
       (elps_magma_block_zgetrf + elps_magma_block_zgetri), &
       Log10(elps_magma_block_zgetrf), &
       Log10(elps_magma_block_zgetri), &
       Log10(elps_magma_block_zgetrf+elps_magma_block_zgetri)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to fill a complex 2d array of size (n)x(m) with random numbers.
! and construct symmetric matrix A such that A^{T} = A
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE

  subroutine getSymDRandomMat_cpu(n,lda,matSymDA)
    implicit none
  
    ! global varaibles

    integer                                    :: LDA,n
    real(dp)                                   :: matSymDA(LDA,*)

    ! local variables

    double precision,dimension(n*n)           :: harvestr
    double precision,dimension(n,n)           :: re_matSymDA

    integer                                   :: i,j

    ! start of the execution commands

    matSymDA(1:n,1:n) = 0.0d0         !initializing the matrix

    call random_number( harvestr )
    where( harvestr == 0.0d0 ) harvestr = 1.0d0

    re_matSymDA(1:n,1:n) = 0.0d0 

    do i=1,n
       do j=i,n
          re_matSymDA(i,j) = harvestr( (j - 1)*lda + (i-1) + 1 )
          re_matSymDA(j,i) = re_matSymDA(i,j)
       end do
    end do

    matSymDA(1:n,1:n) = re_matSymDA(1:n,1:n)

    return
  end subroutine getSymDRandomMat_cpu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to fill a complex 2d array of size (n)x(m) with random numbers.
! and construct symmetric matrix A such that A^{T} = A
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE

  subroutine getSymZRandomMat_cpu(n,lda,matSymZA)
    implicit none
  
    ! global varaibles

    integer                                    :: LDA,n
    complex(dp)                                :: matSymZA(LDA,*)

    ! local variables

    double precision,dimension(n*n)           :: harvestr
    double precision,dimension(n*n)           :: harvesti

    double precision,dimension(n,n)           :: re_matSymZA
    double precision,dimension(n,n)           :: im_matSymZA

    integer                                   :: i,j

    ! start of the execution commands

    matSymZA(1:n,1:n) = 0.0d0         !initializing the matrix

    call random_number( harvestr )
    where( harvestr == 0.0d0 ) harvestr = 1.0d0
    call random_number( harvesti )
    where( harvesti == 0.0d0 ) harvesti = 1.0d0

    re_matSymZA(1:n,1:n) = 0.0d0 
    im_matSymZA(1:n,1:n) = 0.0d0 

    do i=1,n
       do j=i,n
          re_matSymZA(i,j) = harvestr( (j - 1)*lda + (i-1) + 1 )
          im_matSymZA(i,j) = harvesti( (j - 1)*lda + (i-1) + 1 )
          re_matSymZA(j,i) = re_matSymZA(i,j)
          im_matSymZA(j,i) = im_matSymZA(i,j)
       end do
    end do

    matSymZA(1:n,1:n) = cmplx(re_matSymZA(1:n,1:n) , im_matSymZA(1:n,1:n) , dp)

    return
  end subroutine getSymZRandomMat_cpu
type(image_smat_commander)         :: ximage_smat
type(iminfo_commander)             :: ximinfo
type(integrate_movies_commander)   :: xintegrate_movies
type(map2ptcls_commander)          :: xmap2ptcls
type(merge_algndocs_commander)     :: xmerge_algndocs
type(merge_similarities_commander) :: xmerge_similarities
type(multiptcl_init_commander)     :: xsimple_multiptcl_init
type(npeaks_commander)             :: xnpeaks
case('simple_npeaks')
!==Program simple_npeaks
!
! <npeaks/begin> is a program for checking the number of nonzero orientation weights (number of correlation peaks included in the weighted reconstruction) <npeaks/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
   !set required keys
   keys_required(1)  = 'smpd'
   keys_required(2)  = 'box'
   keys_required(3)  = 'lp'
   !set optionmal keys
   keys_optional(1)  = 'nspace'
   keys_optional(2)  = 'moldiam'
   keys_optional(3)  = 'pgrp'

   ! parse command line
   call cline%parse(keys_required(:3), keys_optional(:3))

   !execute
   call xnpeaks%execute(cline)

case('simple_multiptcl_init')

   !set required keys
   keys_required(1)  = 'stk'
   keys_required(2)  = 'smpd'
   keys_required(3)  = 'oritab'
   !set optionmal keys
   keys_optional(1)  = 'msk'
   keys_optional(2)  = 'nstates'
   keys_optional(3)  = 'lp'
   keys_optional(4)  = 'eo'
   keys_optional(5)  = 'frac'
   keys_optional(6)  = 'nthr'
   keys_optional(7)  = 'norec'
   keys_optional(8)  = 'state2split'
   !less optional
   keys_optional(9)  = 'mul'
   keys_optional(10)  = 'ctf'
   keys_optional(11)  = 'kv'
   keys_optional(12)  = 'fraca'
   keys_optional(13)  = 'cs'
   keys_optional(14)  = 'deftab'
   keys_optional(15)  = 'errify'
   keys_optional(16)  = 'inner'
   keys_optional(17)  = 'width'
   keys_optional(18)  = 'zero'

   ! parse command line
   call cline%parse(keys_required(:3), keys_optional(:18))

   !execute
   call xsimple_multiptcl_init%execute(cline)

case('simple_merge_similarities')
  !==Program simple_merge_similarities
!
! <simple\_merge_similarities/begin> is a program for splitting calculations between pairs of objects 
! into balanced partitions. <simple\_merge_similarities/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2016
!

   !set required keys
   keys_required(1)  = 'nptcls'
   !set optionmal keys
   keys_optional(1)  = 'npart'

   ! parse command line
   call cline%parse(keys_required(:1), keys_optional(:1))
 
   !execute
   call xmerge_similarities%execute(cline)

case('simple_merge_algndocs')
!==Program simple_merge_algndocs
!
! <merge_algndocs/begin> is a program for merging alignment documents produced by PRIME2D/3D 
! when run in distributed mode using \texttt{distr\_simple.pl}. <merge_algndocs/end>
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund 2015
!
   !set required keys
   keys_required(1)  = 'fbody'
   keys_required(2)  = 'ndocs'
   keys_required(3)  = 'outfile'
   !set optionmal keys
   keys_optional(1)  = 'oritab'
   
   ! parse command line
   call cline%parse(keys_required(:3), keys_optional(:1))
 
   !execute
   call xmerge_algndocs%execute(cline)

case('simple_map2ptcls')
!==Program simple_map2ptcls
!
! <map2ptcls/begin> is a program for mapping parameters that have been obtained using class averages to 
! the individual particle images. There are many functionalities present that will become critical in
! future releases. Right now we recommend using this program exclusively to exclude the particles
! corresponding to deselected class averages. See the workflows section of the manual for further info.
! <map2ptcls/end> 
!
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_.
! Redistribution and modification is regulated by the GNU General Public License.
! *Author:* Hans Elmlund X-mas 2015
!
   !set required keys
   keys_required(1)  = 'stk'
   keys_required(2)  = 'stk2'
   keys_required(3)  = 'stk3'
   keys_required(4)  = 'oritab'
   !set optionmal keys
   keys_optional(1)  = 'oritab2'
   keys_optional(2)  = 'comlindoc'
   keys_optional(3)  = 'doclist'
   keys_optional(4)  = 'deftab'
   keys_optional(5)  = 'outfile'
   keys_optional(6)  = 'mul'
   keys_optional(7)  = 'nthr'

   ! parse command line
   call cline%parse(keys_required(:4), keys_optional(:7))
 
   !execute
   call xmap2ptcls%execute(cline)


   public :: image_smat_commander
public :: iminfo_commander
public :: integrate_movies_commander
public :: map2ptcls_commander
public :: merge_algndocs_commander
public :: merge_similarities_commander
public :: multiptcl_init_commander
public :: npeaks_commander


type, extends(commander_base) :: image_smat_commander
 contains
   procedure :: execute      => exec_image_smat
end type image_smat_commander

type, extends(commander_base) :: iminfo_commander
 contains
   procedure :: execute      => exec_iminfo
end type iminfo_commander

type, extends(commander_base) :: integrate_movies_commander
 contains
   procedure :: execute      => exec_integrate_movies
end type integrate_movies_commander

type, extends(commander_base) :: map2ptcls_commander
 contains
   procedure :: execute      => exec_map2ptcls
end type map2ptcls_commander

type, extends(commander_base) :: merge_algndocs_commander
 contains
   procedure :: execute      => exec_merge_algndocs
end type merge_algndocs_commander

type, extends(commander_base) :: merge_similarities_commander
 contains
   procedure :: execute      => exec_merge_similarities
end type merge_similarities_commander

type, extends(commander_base) :: multiptcl_init_commander
 contains
   procedure :: execute      => exec_multiptcl_init
end type multiptcl_init_commander

type, extends(commander_base) :: npeaks_commander
 contains
   procedure :: execute      => exec_npeaks
end type npeaks_commander




  subroutine exec_npeaks(self,cline)
    use simple_oris,    only: oris
    class(npeaks_commander), intent(inout) :: self
    class(cmdline),           intent(inout) :: cline
    !local variables
    type(build)   :: b
    type(params)  :: p
    if( .not. cline%defined('nspace') ) call cline%set('nspace', 1000.)
    if( .not. cline%defined('lp') )     call cline%set('lp',       20.)
    p = params(cline) ! parameters generated
    call b%build_general_tbox(p, cline)
    p%npeaks = min(10,b%e%find_npeaks(p%lp, p%moldiam))
    write(*,'(A,1X,I4)') '>>> NPEAKS:', p%npeaks
    ! END GRACEFULLY
    call simple_end('**** SIMPLE_NPEAKS NORMAL STOP ****')
    return
  end subroutine exec_npeaks
  
  subroutine exec_multiptcl_init(self,cline)
    use simple_rec_master, only: exec_rec, exec_eorec
    class(multiptcl_init_commander), intent(inout) :: self
    class(cmdline),           intent(inout) :: cline
    !local variables
    type(params)      :: p
    type(build)       :: b
    real, allocatable :: optlp

    !
    !TODO: needs to verify the setter here
    !
    call cline%set('trs', 3.) ! to assure that shifts are being used

    p = params(cline)         ! constants & derived constants produced
    call b%build_general_tbox(p, cline)
    if( cline%defined('state2split') )then
       call b%a%split_state(p%state2split)
       p%nstates = p%nstates+1
    else
       call b%a%rnd_states(p%nstates)
    endif
    if( p%errify.eq.'yes' )call errify_oris
    if( p%norec .eq. 'no' )then
       if( cline%defined('lp') )then
          call b%build_rec_tbox(p)
          call exec_rec(b, p, cline, 'startvol')
       else
          call b%build_eo_rec_tbox(p)
          call exec_eorec(b, p, cline, 'startvol')
       endif
    endif
    if( p%zero.eq.'yes' )call b%a%set_all('corr',0.)
    call b%a%write('multiptcl_startdoc.txt')
    call simple_end('**** SIMPLE_MULTIPTCL_INIT NORMAL STOP ****')    

  contains
    
    subroutine errify_oris()
      real :: sherr,angerr
      if( cline%defined('trs') )then
         sherr = p%trs
      else
         sherr = 3.
      endif
      angerr = 15.
      write(*,'(A,F6.2,A)')'>>> IN-PLANE SHIFT   ERROR INTRODUCED: ',sherr,' PIXELS'
      write(*,'(A,F6.2,A)')'>>> IN-PLANE ANGULAR ERROR INTRODUCED: ',angerr,' DEGREES'
      call b%a%introd_alig_err( angerr,sherr )
    end subroutine errify_oris

  end subroutine exec_multiptcl_init
    
  subroutine exec_merge_similarities(self,cline)
    use simple_map_reduce, only: merge_similarities_from_parts
    class(merge_similarities_commander), intent(inout) :: self
    class(cmdline),           intent(inout) :: cline
    !local variables
    type(params)      :: p
    real, allocatable :: simmat(:,:)
    integer           :: filnum, io_stat

    p = params(cline) ! parameters generated
    simmat = merge_similarities_from_parts(p%nptcls, p%npart)
    filnum = get_fileunit()
    open(unit=filnum, status='REPLACE', action='WRITE', file='smat.bin', access='STREAM')
    write(unit=filnum,pos=1,iostat=io_stat) simmat
    if( io_stat .ne. 0 )then
       write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when writing to smat.bin'
       stop 'I/O error; simple_merge_similarities'
    endif
    close(filnum)
    call simple_end('**** SIMPLE_MERGE_SIMILARITIES NORMAL STOP ****')
    return
  end subroutine exec_merge_similarities

  subroutine exec_merge_algndocs(self,cline)
    use simple_oris,    only: oris
    class(merge_algndocs_commander), intent(inout) :: self
    class(cmdline),           intent(inout) :: cline
    !local variables
    type(params)          :: p
    type(oris)            :: o, o_read
    integer               :: i, j, istart, istop, ptcls_per_part, nentries, leftover, cnt, nentries_all, numlen
    character(len=STDLEN) :: fname
    logical               :: here, useoritab

    p = params(cline) ! parameters generated
    useoritab = .false.
    if( cline%defined('oritab') )then
       if( file_exists(p%oritab) ) useoritab = .true.
    endif
    if( useoritab )then
       if( nlines(p%oritab) /= p%nptcls )then
          stop 'the inputted nptcls is not consistent with the nptcls in oritab!'
       endif
       ! create object for orientations
       o = oris(p%nptcls)
       ! read previous orientations
       call o%read(p%oritab)
    endif
    ptcls_per_part = p%nptcls/p%ndocs
    leftover       = p%nptcls-ptcls_per_part*p%ndocs
    istop          = 0
    numlen         = len(int2str(p%ndocs))
    do i=1,p%ndocs
       fname  = trim(adjustl(p%fbody))//int2str_pad(i,numlen)//'.txt'
       if( i == p%ndocs )then
          istart = istop+1
          istop  = p%nptcls
       else
          if( leftover == 0 )then
             istart = istop+1;
             istop  = istart+ptcls_per_part-1;
          else
             istop  = i*(ptcls_per_part+1)
             istart = istop-(ptcls_per_part+1)+1
             leftover = leftover-1
          endif
       endif
       ! calculate the number of all entries
       nentries_all = istop-istart+1 
       ! calculate the actual number of entries
       inquire(FILE=fname, EXIST=here)
       if( here )then
          nentries = nlines(fname)
       else
          nentries = 0
       endif
       ! check if oritab is there to fill-in blanks
       if( nentries < nentries_all )then
          if( .not. useoritab )then
             stop 'need previous oritab to fill-in blanks; simple_merge_algndocs'
          endif
       endif
       ! print partition info
       write(*,'(a,1x,i3,1x,a,1x,i6,1x,i6)') 'partition:', i, 'from/to:', istart, istop
       if( nentries > 0 )then
          o_read = oris(nentries)
          call o_read%read(fname)
       endif
       ! read
       if( useoritab )then ! read and fill-in from oritab
          cnt = 0
          do j=istart,istop
             cnt = cnt+1
             if( cnt <= nentries )then
                call o%set_ori(j,o_read%get_ori(cnt))
             else
                exit
             endif
          end do
       else                                ! just merge (all ptcls is there)
          call o%merge(o_read)
       endif
    end do
    call o%write(p%outfile)
    call simple_end('**** SIMPLE_MERGE_ALGNDOCS NORMAL STOP ****')

    return
  end subroutine exec_merge_algndocs
  
  subroutine exec_map2ptcls(self,cline)
    use simple_oris,    only: oris
    use simple_ori,     only: ori
    use simple_image,   only: image
    use simple_corrmat
    class(map2ptcls_commander), intent(inout) :: self
    class(cmdline),           intent(inout) :: cline
    !local variables
    type state_organiser
       integer, allocatable :: particles(:)
       integer   :: cls_orig = 0 
       integer   :: cls_sel  = 0
       integer   :: istate   = 0
       type(ori) :: ori3d
    end type state_organiser
    type(params)                       :: p
    type(build)                        :: b
    type(state_organiser), allocatable :: labeler(:)
    type(image), allocatable           :: imgs_sel(:), imgs_cls(:)
    type(oris)                         :: o_comlindoc, o_state, a_copy, o_oritab2
    type(ori)                          :: ori2d, ori_comp
    real, allocatable                  :: correlations(:,:)
    integer, allocatable               :: statepops(:), state_particles(:), rejected_particles(:)
    integer                            :: isel, nsel, loc(1), iptcl, pind, icls
    integer                            :: nlines_oritab, nlines_oritab2
    integer                            :: nlines_comlindoc, nlines_deftab
    integer                            :: cnt, istate, funit, iline, nls, lfoo(3)
    logical, allocatable               :: statedoc_exists(:), selected(:)
    character(len=STDLEN)              :: statedoc
    real                               :: corr
    logical, parameter                 :: debug=.false.
    
    p = params(cline)                   ! parameters generated
    call b%build_general_tbox(p, cline) ! general objects built
    ! find number of selected cavgs
    call find_ldim_nptcls(p%stk2, lfoo, nsel)
    if( debug ) print *, 'nsel: ', nsel
    ! find number of original cavgs
    call find_ldim_nptcls(p%stk3, lfoo, p%ncls)
    if( debug ) print *, 'ncls: ', p%ncls
    if( p%ncls < nsel )then
       stop 'nr of original clusters cannot be less than the number of selected ones'
    endif
    ! find number of lines in input document
    nlines_oritab = nlines(p%oritab)
    if( debug ) print *, 'nlines_oritab: ', nlines_oritab
    if( nlines_oritab /= p%nptcls )then 
       stop 'nr lines in oritab .ne. nr images in particle stack; must be congruent!'
    endif
    if( cline%defined('deftab') )then
       nlines_deftab = nlines(p%deftab)
       if( nlines_oritab /= nlines_deftab )then
          stop 'nr lines in oritab .ne. nr lines in deftab; must be congruent!'
       endif
    endif
    if( cline%defined('doclist') )then
       if( .not. cline%defined('comlindoc') )then
          if( nlines(p%doclist) /= 1 )then
            stop 'need a comlindoc together with statelist'
         endif
      endif
   endif
   if( cline%defined('comlindoc') .and. cline%defined('oritab2') )then
      stop 'either comlindoc or oritab2 can be inputted, not both'
   endif
   if( cline%defined('doclist') .and. cline%defined('oritab2') )then
      stop 'either doclist or oritab2 can be inputted, not both'
   endif
   allocate(imgs_sel(nsel), imgs_cls(p%ncls))
   ! read images
   do isel=1,nsel
      call imgs_sel(isel)%new([p%box,p%box,1], p%smpd)
      call imgs_sel(isel)%read(p%stk2, isel)
   end do
   do icls=1,p%ncls
      call imgs_cls(icls)%new([p%box,p%box,1], p%smpd)
      call imgs_cls(icls)%read(p%stk3, icls)
   end do
   write(*,'(a)') '>>> CALCULATING CORRELATIONS'
   call calc_cartesian_corrmat(imgs_sel, imgs_cls, correlations)
   ! find selected clusters & map selected to original clusters & extract the particle indices
   allocate(labeler(nsel), selected(p%ncls))
   ! initialise selection array
   selected = .false.
   write(*,'(a)') '>>> MAPPING SELECTED TO ORIGINAL CLUSTERS'
   do isel=1,nsel
      loc = maxloc(correlations(isel,:))
      labeler(isel)%cls_orig  = loc(1)
      selected(loc(1)) = .true.
      if( debug ) print *, 'found orig clsind: ', labeler(isel)%cls_orig
      labeler(isel)%cls_sel   = isel
      if( debug ) print *, 'selected class index: ', labeler(isel)%cls_sel
      labeler(isel)%particles = b%a%get_cls(labeler(isel)%cls_orig)
      if( debug ) print *, 'got this number of partices: ', size(labeler(isel)%particles)
   end do
   ! erase deselected (by setting their state to zero)
   do icls=1,p%ncls
      if( selected(icls) ) cycle
      if( b%a%get_clspop(icls) > 0 )then
         rejected_particles = b%a%get_cls(icls)
         do iptcl=1,size(rejected_particles)
            call b%a%set(rejected_particles(iptcl), 'state', 0.)
         end do
         deallocate(rejected_particles)
      endif
   end do
   ! parse state info
   if( cline%defined('comlindoc') )then
      write(*,'(a)') '>>> PROCESSING COMLIN STATE ASSIGNMENT DOCUMENT'
      nlines_comlindoc = nlines(p%comlindoc)
      if( nsel /= nlines_comlindoc )then
         stop 'nr lines in comlindoc .ne. nr of selected clusters; must be congruent!'
      endif
      ! make a new oris object and read in the comlin clustering (state) info
      o_comlindoc = oris(nsel)
      call o_comlindoc%read(p%comlindoc)
      if( .not. o_comlindoc%isthere('state') )then
         write(*,*) 'no state labeling in comlindoc, perhaps you clustered with label=class'
         stop 'please, re-cluster with label=state'
      endif
      ! set the number of states
      p%nstates = o_comlindoc%get_nstates()
      ! map states to selected cavgs
      do isel=1,nsel
         labeler(isel)%istate = nint(o_comlindoc%get(isel, 'state'))
      end do
      ! extract the state populations
      allocate( statepops(p%nstates) )
      do istate=1,p%nstates
         statepops(istate) = o_comlindoc%get_statepop(istate)
      end do
   else
      ! set default values for nstates, statepop and state
      p%nstates         = 1
      allocate( statepops(1) )
      statepops(1)      = nsel
      labeler(:)%istate = 1
   endif
   if( cline%defined('oritab2') )then
      if( .not. file_exists(p%oritab2) ) stop 'Inputted oritab2 does not exist in the cwd'
      nlines_oritab2 = nlines(p%oritab2)
      if( nlines_oritab2 /= nsel ) stop 'Nr lines in oritab2 /= nr of selected cavgs'
      o_oritab2 = oris(nsel)
      call o_oritab2%read(p%oritab2)
      ! compose orientations and set states
      do isel=1,nsel
         ! get 3d ori
         labeler(isel)%ori3d  = o_oritab2%get_ori(isel)
         labeler(isel)%istate = nint(labeler(isel)%ori3d%get('state'))
         corr                 = labeler(isel)%ori3d%get('corr')
         do iptcl=1,size(labeler(isel)%particles)
            ! get particle index 
            pind = labeler(isel)%particles(iptcl)
            ! get 2d ori
            ori2d = b%a%get_ori(pind)
            if( cline%defined('mul') )then
               call ori2d%set('x', p%mul*ori2d%get('x'))
               call ori2d%set('y', p%mul*ori2d%get('y'))
            endif
            ! transfer original parameters in b%a 
            ori_comp = b%a%get_ori(pind)
            ! compose ori3d and ori2d
            call labeler(isel)%ori3d%compose3d2d(ori2d, ori_comp)
            ! set parameters in b%a
            call b%a%set_ori(pind,ori_comp)
            call b%a%set(pind, 'corr', corr)
         end do
      end do
   endif
   ! map states to particles
   do isel=1,nsel
      do iptcl=1,size(labeler(isel)%particles)
         ! get particle index
         pind = labeler(isel)%particles(iptcl)
         call b%a%set(pind, 'state', real(labeler(isel)%istate))
      end do
   end do
   ! parse ori info
   if( cline%defined('doclist') )then
      write(*,'(a)') '>>> COMBINING 3D ORIS (CAVGS) WITH 2D ALIGNMENT (PARTICLES)'
      if( nlines(p%doclist) /= p%nstates )then
         stop 'the number of lines in doclist does not match the number of states in comlindoc'
      endif
      allocate(statedoc_exists(p%nstates))
      ! read in 3d orientations
      funit = get_fileunit()
      open(unit=funit, status='old', file=p%doclist)
      do istate=1,p%nstates
         ! read the relevant statedoc
         read(funit,'(a256)') statedoc
         statedoc_exists(istate) = file_exists(statedoc)
         if( statedoc_exists(istate) )then
            nls = nlines(statedoc)
            if( nls /= statepops(istate) )then
               write(*,*) 'the nr of lines in statedoc: ', trim(statedoc),&
                    'does not match pop size: ', statepops(istate), 'in comlindoc'
               stop
            endif
            o_state = oris(nls)
            call o_state%read(statedoc)
         else
            ! make a fake o_state
            o_state = oris(statepops(istate))
            do iline=1,statepops(istate)
               call o_state%set(iline, 'state', 0.)
            end do
            statepops(istate) = 0
         endif
         cnt = 0
         do isel=1,nsel
            if( labeler(isel)%istate == istate )then
               cnt = cnt+1
               labeler(isel)%ori3d = o_state%get_ori(cnt)
            endif
         end do
      end do
      close(funit)
      ! wipe out the states for which no docs are provided
      do isel=1,nsel
         do iptcl=1,size(labeler(isel)%particles)
            ! get particle index
            pind = labeler(isel)%particles(iptcl)
            ! get state index
            istate = nint(b%a%get(pind, 'state'))
            if( .not. statedoc_exists(istate) )then
               call b%a%set(pind, 'state', 0.)
            endif
         end do
      end do
      ! compose orientations
      do isel=1,nsel
         do iptcl=1,size(labeler(isel)%particles)
            ! get particle index
            pind = labeler(isel)%particles(iptcl)
            ! get 2d ori
            ori2d = b%a%get_ori(pind)
            if( cline%defined('mul') )then
               call ori2d%set('x', p%mul*ori2d%get('x'))
               call ori2d%set('y', p%mul*ori2d%get('y'))
            endif
            ! transfer original parameters in b%a
            ori_comp = b%a%get_ori(pind)
            ! compose ori3d and ori2d
            call labeler(isel)%ori3d%compose3d2d(ori2d, ori_comp)
            ! set parameters in b%a
            call b%a%set_ori(pind,ori_comp)
            call b%a%set(pind, 'corr',  labeler(isel)%ori3d%get('corr'))
         end do
      end do
      ! relabel states in consequtive order
      if( any(statepops == 0) )then
         a_copy = b%a
         cnt    = 0
         do istate=1,p%nstates
            if( statepops(istate) > 0 )then
               cnt = cnt+1
               write(*,'(a,1x,i3,1x,a,1x,i3)') '>>> THE STATE THAT WAS FORMERLY:', istate, 'IS NOW:', cnt
               state_particles = nint(a_copy%get_ptcls_in_state(istate))
               do iptcl=1,size(state_particles)
                  call b%a%set(state_particles(iptcl),'state',real(cnt))
               end do
            else
               write(*,'(a,1x,i3,1x,a)') '>>> THE STATE THAT WAS FORMERLY:', istate, 'HAS BEEN EXCLUDED'
            endif
         end do
      endif
   endif
   call b%a%write(p%outfile)
   call simple_end('**** SIMPLE_MAP2PTCLS NORMAL STOP ****')

   return
 end subroutine exec_map2ptcls
 
  subroutine exec_integrate_movies(self,cline)
    use simple_imgfile,   only: imgfile
    use simple_image,     only: image

    class(integrate_movies_commander), intent(inout) :: self
    class(cmdline),           intent(inout) :: cline
    !local variables
    type(params)                       :: p
    type(build)                        :: b
    integer                            :: nmovies, nframes, frame, alloc_stat, lfoo(3)
    integer                            :: numlen, ldim(3), fromto(2), movie, ifoo
    character(len=STDLEN), allocatable :: movienames(:)
    character(len=:), allocatable      :: cpcmd, new_name
    real                               :: x, y
    type(image), allocatable           :: img_frames(:)
    type(image)                        :: img_sum, pspec
    logical, parameter                 :: debug = .false.
    p = params(cline,checkdistr=.false.) ! constants & derived constants produced
    call b%build_general_tbox(p,cline,do3d=.false.)
    call read_filetable(p%filetab, movienames)
    nmovies = size(movienames)
    if( debug ) write(*,*) 'read the movie filenames'
    ! FIND LDIM AND NUMLEN (LENGTH OF NUMBER STRING)
    call find_ldim_nptcls(movienames(1), ldim, ifoo)
    if( debug ) write(*,*) 'logical dimension: ', ldim
    ldim(3) = 1 ! to correct for the stupide 3:d dim of mrc stacks
    numlen = len(int2str(nmovies))
    if( debug ) write(*,*) 'length of number string: ', numlen
    ! DETERMINE LOOP RANGE
    if( cline%defined('part') )then
       if( cline%defined('fromp') .and. cline%defined('top') )then
          fromto(1) = p%fromp
          fromto(2) = p%top
       else
          stop 'fromp & top args need to be defined in parallel execution; simple_integrate_movies'
       endif
    else
       fromto(1) = 1
       fromto(2) = nmovies
    endif
    if( debug ) write(*,*) 'fromto: ', fromto(1), fromto(2)
    ! CREATE SUM
    call img_sum%new([ldim(1),ldim(2),1], p%smpd)
    ! LOOP OVER EXPOSURES (MOVIES)
    do movie=fromto(1),fromto(2)
       if( .not. file_exists(movienames(movie)) )then
          write(*,*) 'inputted movie stack does not exist: ', trim(adjustl(movienames(movie)))
       endif
       ! GET NUMBER OF FRAMES FROM STACK
       call find_ldim_nptcls(movienames(movie), lfoo, nframes)
       if( debug ) write(*,*) 'number of frames: ', nframes
       ! CREATE FRAMES & READ
       allocate(img_frames(nframes), stat=alloc_stat)
       call alloc_err('In: simple_integrate_movies', alloc_stat)
       img_sum = 0.
       do frame=1,nframes
          call img_frames(frame)%new([ldim(1),ldim(2),1], p%smpd)
          call img_frames(frame)%read(movienames(movie),frame)
          if( cline%defined('oritab') )then
             ! SHIFT FRAME ACCORDING TO GLOBAL SHIFT (DRIFT CORRECTION)
             if( b%a%isthere(movie, 'x'//int2str(frame)) .and. b%a%isthere(movie, 'y'//int2str(frame)) )then
                call img_frames(frame)%fwd_ft
                x = b%a%get(movie, 'x'//int2str(frame))
                y = b%a%get(movie, 'y'//int2str(frame))
                if( debug ) print *, 'shifting frame: ', x, y
                call img_frames(frame)%shift(-x,-y)
                call img_frames(frame)%bwd_ft
                call img_sum%add(img_frames(frame))
            else
               stop 'no shift parameters available for alignment in oritab'
            endif
         else
            call img_sum%add(img_frames(frame))
         endif
      end do
      ! RENAME MOVIE
      allocate(new_name, source=trim(adjustl(p%fbody))//int2str_pad(movie, numlen)//p%ext)
      if( .not. file_exists(new_name))then
         allocate(cpcmd, source='cp '//trim(adjustl(movienames(movie)))//' ./'//new_name)
         call system(cpcmd)
         deallocate(cpcmd)
      endif
      ! DESTROY OBJECTS AND DEALLOCATE
      do frame=1,nframes
         call img_frames(frame)%kill
      end do
      deallocate(new_name,img_frames)
      call img_sum%write(trim(adjustl(p%fbody))//'_intg'//int2str_pad(movie, numlen)//p%ext)
      pspec = img_sum%mic2spec(p%pspecsz, trim(adjustl(p%speckind)))
      call pspec%write(trim(adjustl(p%fbody))//'_pspec'//int2str_pad(movie, numlen)//p%ext)
   end do
   call simple_end('**** SIMPLE_INTEGRATE_MOVIES NORMAL STOP ****')

   return
 end subroutine exec_integrate_movies

  subroutine exec_iminfo(self,cline)
    use simple_image,   only: image
    use simple_imgfile
    class(iminfo_commander), intent(inout) :: self
    class(cmdline),           intent(inout) :: cline
    !local variables
    type(params)      :: p
    type(image)       :: img
    character(len=20) :: conv
    character(len=1)  :: form
    integer           :: ldim(3), maxim, i, iform, n_nans, mode
    real              :: smpd, sdev, ave, minv, maxv

    p = params(cline) ! constants & derived constants produced
    if( cline%defined('fname') )then
       call find_ldim_nptcls(p%fname, ldim, maxim, doprint=.true.)
    endif
    p%box  = ldim(1)
    p%smpd = smpd
    call img%new([ldim(1),ldim(2),1],p%smpd)
    if( p%vis .eq. 'yes' .or. p%stats .ne. 'no' )then
       do i=1,maxim
          call img%read(p%fname, i)
          if( p%stats .ne. 'no' )then
             call img%cure(maxv, minv, ave, sdev, n_nans)
             if( p%stats .eq. 'print' .or. n_nans > 0 )then
                write(*,*) '*********IMAGE********', i, '*******'
                write(*,*) 'maxv = ',   maxv
                write(*,*) 'minv = ',   minv
                write(*,*) 'ave = ',    ave
                write(*,*) 'sdev = ',   sdev
                write(*,*) 'n_nans = ', n_nans
             endif
          endif
          if( p%vis .eq. 'yes' ) call img%vis
       end do
    endif
    call simple_end('**** SIMPLE_IMINFO NORMAL STOP ****')
 
    return
  end subroutine exec_iminfo
  
  subroutine exec_image_smat(self, cline)
    use simple_corrmat  ! singleton
    use simple_ori,     only: ori
    use simple_imgfile, only: imgfile
    use simple_image,   only: image
    class(image_smat_commander), intent(inout) :: self
    class(cmdline),           intent(inout) :: cline
    !local variables
    type(params)         :: p
    type(build)          :: b
    integer              :: iptcl, alloc_stat, funit, io_stat
    real, allocatable    :: corrmat(:,:)
    logical              :: debug=.false.
    p = params(cline, .false.)                    ! constants & derived constants produced
    call b%build_general_tbox(p, cline, .false., .true.) ! general objects built (no oritab reading)
    allocate(b%imgs_sym(p%nptcls), stat=alloc_stat)
    call alloc_err('In: simple_image_smat, 1', alloc_stat)
    do iptcl=1,p%nptcls
       call b%imgs_sym(iptcl)%new([p%box,p%box,1], p%smpd, p%imgkind)
       call b%imgs_sym(iptcl)%read(p%stk, iptcl)
    end do
    write(*,'(a)') '>>> CALCULATING CORRELATIONS'
    if( cline%defined('lp') )then
       if( .not. cline%defined('msk') ) stop 'need mask radius (msk) 4 Fourier corr calc!'
       call calc_cartesian_corrmat(b%imgs_sym, corrmat, p%msk, p%lp)
    else
       if( cline%defined('msk') )then
          call calc_cartesian_corrmat(b%imgs_sym, corrmat, p%msk)
       else
          call calc_cartesian_corrmat(b%imgs_sym, corrmat)
       endif
    endif
    funit = get_fileunit()
    open(unit=funit, status='REPLACE', action='WRITE', file='img_smat.bin', access='STREAM')
    write(unit=funit,pos=1,iostat=io_stat) corrmat
    if( io_stat .ne. 0 )then
       write(*,'(a,i0,a)') 'I/O error ', io_stat, ' when writing to image_smat.bin'
       stop 'I/O error; simple_image_smat'
    endif
    close(funit)
    call simple_end('**** SIMPLE_IMAGE_SMAT NORMAL STOP ****')

    return
  end subroutine exec_image_smat

public :: recvol_commander

  type, extends(commander_base) :: recvol_commander
 contains
   procedure :: execute      => exec_recvol
end type recvol_commander


  subroutine exec_recvol(self,cline)
use simple_rec_master, only: exec_rec
class(recvol_commander), intent(inout) :: self
    class(cmdline),           intent(inout) :: cline
    !local variables
type(params)  :: p
type(build)   :: b

call cline%set('trs', 5.) ! to assure that shifts are being used
p = params(cline)         ! constants & derived constants produced
call b%build_general_tbox(p, cline)
call b%build_rec_tbox(p)
call exec_rec(b, p, cline)
call simple_end('**** SIMPLE_RECVOL NORMAL STOP ****')    


    return
  end subroutine exec_recvol
