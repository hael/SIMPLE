    !>  \brief  is for generating all rotational correlations for all nptcls particles and all nrefs
    !           references (nptcls=nrefs) (assumes that the first and second dimensions of the particle
    !           matrix pfts_ptcls have been expanded)
    subroutine gencorrs_all( self, corrmat2dout, inplmat2dout )
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_jiffys, only: progress
        use simple_jiffys, only: alloc_err
        use simple_deviceQuery_gpu
        class(polarft_corrcalc), intent(inout) :: self                         !< instance

        real, intent(out)             :: corrmat2dout(self%nptcls,self%nptcls) !< output correlation matrix
        integer, intent(out)          :: inplmat2dout(self%nptcls,self%nptcls) !< output inplane rot index matrix

        real    ::     corrmat2d(self%nptcls,self%nptcls)
        real    ::     corrmat3d(self%nptcls,self%nptcls,self%nrots)
        real    :: hadamard_prod(self%nptcls,self%refsz,self%nk)
        real(sp), allocatable         :: sumb_vec(:)
        complex(sp), allocatable      :: temp(:,:,:)

        integer :: ptcl_ind, iptcl, iref, irot, d1lim(2), d2lim(2), indices(self%nptcls)
        integer :: inplmat2d(self%nptcls,self%nptcls)

        !kernel variables
        integer                       :: err
        integer                       :: lda,ldb,ldc
        real(sp)                      :: alpha

        !deviceDetails
        type(deviceQuery_gpu)         :: devQ
        type(deviceDetails)           :: devD
        !type(deviceDetails),pointer, dimension(:) :: a_devD_corr=>NULL()
        type(deviceDetails),allocatable :: a_devD(:)
        integer                       :: ndev
        !timming variables
        double precision              :: elps_L
        double precision,dimension(2) :: st_L_r, et_L_r
        ! data structure
        logical, parameter :: ibench=.true.        !< benchmark result or not
        logical, parameter :: ibench_write=.true.  !< write benchmark result or not
        
        logical, parameter :: debug=.false.          !< debug indicator
        logical, parameter :: debug_cpu=.false.     !false/true
        logical, parameter :: debug_high=.false.     !true/false
        logical, parameter :: debug_write=.false.   !false/true
        logical, parameter :: debug_write_C=.false. !false/true

        type(t_debug_gpu)             :: s_debug_gpu
        type(t_bench)                 :: s_bench

        !temporray index variables
        integer :: i,j,k
        integer :: idev
        
        !function calls
        integer :: get_polarft_corr_gpu_c_
        integer :: get_dev_count_c
        integer :: get_polarft_gencorrall_gpu_c_

        if( allocated(self%pfts_refs_ctf) )then
            stop 'CTF modulation of the references not possible with simple_polarft_corrcalc :: gencorrs_all'
        endif

        allocate(sumb_vec(self%nptcls))
        allocate(temp(self%nptcls,self%refsz,self%nk))
        
        corrmat3d     = 0.
        corrmat2d     = 0.
        hadamard_prod = 0.
        sumb_vec      = 0.0
        temp          = 0.0
        inplmat2d     = 0
        self%s_polar%r_polar         = R_POLAR_D
        self%s_polar%sumasq_polar    = SUMASQ_POLAR_D
        self%s_polar%sumbsq_polar    = SUMBSQ_POLAR_D
        self%s_polar%ikrnl           = KERNEL_D
        self%s_polar%threadsPerBlock = THREADSPERBLOCK_D
        self%s_polar%nx              = NX_3D_D
        self%s_polar%ny              = NY_3D_D
        self%s_polar%nz              = NZ_3D_D
        alpha         = 1.0

        !first setup the debugging options first using 0 and 1.
        if (        debug) s_debug_gpu%debug_i         = 1 !else 0
        if (    debug_cpu) s_debug_gpu%debug_cpu_i     = 1
        if (   debug_high) s_debug_gpu%debug_high_i    = 1
        if (  debug_write) s_debug_gpu%debug_write_i   = 1
        if (debug_write_C) s_debug_gpu%debug_write_C_i = 1
#if defined (BENCH)
        if (ibench) s_bench%bench_i = 1
        if (ibench_write) s_bench%bench_write_i = 1
#endif

        write(*,*) "In simple_polarft_corrcalc.f90 use_gpu: ",use_gpu
        write(*,*) "In simple_polarft_corrcalc.f90 has_gpu: ",has_gpu
        
        !New GPU Z N integration starts here
#if defined (CUDA) && defined (MAGMA)

        write(*,*)'******************************************************************'
        write(*,*)'   Device fills in the object(devQ) and the data structure(devD)  '
        write(*,*)'******************************************************************'
        err = get_dev_count_c(devD)
        ndev = devD%ndev
        allocate(a_devD(0:devD%ndev-1))
        do idev = 0, ndev-1
           call devQ%new_deviceQuery_gpu(devD,idev)
           !call Sanity_check_gpu(devQ, devD)
           !mapping the data structures into an array
           a_devD(idev) = devD
        end do
!!#endif
        lda           = self%nptcls
        ldb           = self%nrots
        ldc           = self%nptcls

        err = get_polarft_gencorrall_gpu_c_(a_devD,                 &
                                            self%s_polar,"Z","N",   &
                                            sumb_vec,corrmat3d,     &
                                            self%pfts_refs,         &!pft1
                                            self%pfts_ptcls,        &!pft2
                                            self%sqsums_refs,       &!sqsums_refs
                                            self%sqsums_ptcls,      &!sqsums_ptcls
                                            self%nptcls,self%nrots, &!ipart,nrot
                                            self%nk,                &!nk
                                            lda,ldb,ldc,alpha,      &
                                            s_bench, s_debug_gpu) 

        deallocate(a_devD)
#else
        !New GPU Z N integration ends here

!!!        lda           = self%nptcls
!!!        ldb           = self%refsz
!!!        ldc           = self%nptcls
        
!!!        do ptcl_ind=1,self%nptcls
!!!            d1lim(1) = ptcl_ind
!!!            d1lim(2) = ptcl_ind+self%nptcls-1
!!!            do irot =1,self%nrots
!!!                d2lim(1) = irot
!!!                d2lim(2) = irot+self%refsz-1
!!#if defined (CUDA) && defined (MAGMA)

! GPU starts here
!!                temp = 0.0
!!                temp = self%pfts_ptcls(d1lim(1):d1lim(2),d2lim(1):d2lim(2),:)

!!                err = get_polarft_corr_gpu_c_(a_devD,                          &
!!                                              self%s_polar,"X","N",            &
!!                                              sumb_vec,                        &
!!                                              self%pfts_refs,                  &
!!                                              temp,                            &
!!                                              self%nptcls,self%refsz,self%nk,  &
!!                                              lda,ldb,ldc,alpha,               &
!!                                              s_bench, s_debug_gpu)
! GPU ends here

!!                !$omp parallel default(shared) private(iptcl)
                
!!                !$omp workshare
!!                corrmat3d(:,ptcl_ind,irot) = sumb_vec(:)/&
!!                sqrt(self%sqsums_refs(:)*self%sqsums_ptcls(d1lim(1):d1lim(2)))
!!                !$omp end workshare nowait

!!                !$omp end parallel
!!#else
                !write(*,*) "#elseif defined (CUDA) && defined (MAGMA)"
!!!                !$omp parallel default(shared) private(iptcl)
!!!                !$omp workshare
!!!                hadamard_prod = real(self%pfts_refs(:,:,:)*&
!!!                conjg(self%pfts_ptcls(d1lim(1):d1lim(2),d2lim(1):d2lim(2),:)))       
!!!                !$omp end workshare

                ! correlation normalization

!!!                !$omp do schedule(auto)
!!!                do iptcl=1,self%nptcls
!!!                    corrmat3d(iptcl,ptcl_ind,irot) = sum(hadamard_prod(iptcl,:,:))
!!!                end do
!!!                !$omp end do
                
!!!                !$omp workshare
!!!                corrmat3d(:,ptcl_ind,irot) = corrmat3d(:,ptcl_ind,irot)/&
!!!                sqrt(self%sqsums_refs(:)*self%sqsums_ptcls(d1lim(1):d1lim(2)))
!!!                !$omp end workshare nowait
!!!                !$omp end parallel      
!!#endif        
!!!            end do
!!!        end do


        write(*,*) "**************************** Need to -DCUDA*******"
#endif

        deallocate(sumb_vec)
        deallocate(temp)
        
        !$omp parallel default(shared) private(iptcl)
        !$omp workshare 
        corrmat2d = maxval(corrmat3d, dim=3)
        inplmat2d = maxloc(corrmat3d, dim=3)
        !$omp end workshare
        !$omp do schedule(auto) 
        do iptcl=1,self%nptcls
            indices(iptcl) = iptcl
        end do
        !$omp end do nowait

#if defined (BENCH)
        call gettimeofday_c(st_L_r)
#endif
        !$omp end parallel
        do iref=1,self%nptcls
           if( iref /= 1 )then
              indices = cshift(indices, shift=1)
           endif
           !$omp parallel do schedule(auto) private(iptcl)
           do iptcl=1,self%nptcls
              corrmat2dout(indices(iptcl),iref) = corrmat2d(iref,iptcl)
              inplmat2dout(indices(iptcl),iref) = inplmat2d(iref,iptcl)
           end do
           !$omp end parallel do
        end do
#if defined (BENCH)
        call gettimeofday_c(et_L_r)
        call elapsed_time_c(st_L_r,et_L_r,elps_L)
        write(*,*)"Matrix mapping"
        write(*,*)"---------------------------------------------------"
        write(*,'(f15.8)') elps_L
#endif

    end subroutine gencorrs_all