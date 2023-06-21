! This programs the surface-core NanoParticle (NP) structure from 2D NP class averages
! and 2D reprojections of the 3D volume. A stack containing the differences between
! the class averages and reprojections is output. Average intensities over shells 
! starting from the NP center out to the NP surface are reported in a .csv file. 
! If the surface-core model is correct, the fractional differences should be higher 
! towards the surface, where the structure is more flexible and dynamic.
program simple_test_2D_core_finder

    include 'simple_lib.f08'
    use simple_image,        only: image
    implicit none
    
    real, parameter     :: smpd=0.358, nprad_A=14.0, shell_size_A = 2.0, res=1.0, maskdiam_A=20.0
    logical, parameter  :: normalize_real_space = .false.

    type(image)                     :: cavg, reproj, diff
    type(image)                     :: temp
    real(kind=c_float), pointer     :: rmat(:,:,:)=>null(), rmat_cavg(:,:,:)=>null(), rmat_reproj(:,:,:)=>null()
    real                :: nprad_vox=nprad_A/smpd, shell_size_vox=shell_size_A/smpd, mskdiam_vox=maskdiam_A/smpd
    real                :: r, mean, center(2), minmax(2)
    integer             :: iref, icavg, ldim1(3), ldim2(3), ldim_refs(3), ifoo, nshells, i, j, n, funit
    logical, allocatable       :: is_populated(:), mask(:,:,:)
    character(len=256)         :: fn_cavgs, fn_reprojs
    character(*), parameter    :: mask_type="soft"
    character(*), parameter    :: fn_diff='cavgs_minus_reprojections.mrc', fn_shells='shells.mrc'
    character(*), parameter    :: fn_cavgs_out='cavgs_out.mrc', fn_reprojs_out='reprojs_out.mrc'
    character(*), parameter    :: fn_results='results.csv'
    character(*), parameter    :: results_header='INDEX'//CSV_DELIM//'RMIN'//CSV_DELIM//'RMAX'//CSV_DELIM//&
                                  &'COUNT'//CSV_DELIM//'MEAN_DIFF'//CSV_DELIM//'MIN_DIFF'//CSV_DELIM//'MAX_DIFF'

    ! Handle command-line input
    if( command_argument_count() /= 2 )then
        write(logfhandle,'(a)')  'Usage: simple_test_2D_core_finder cavgs.mrc reprojs.mrc'
        write(logfhandle,'(a)')  'cavgs.mrc: stack of 2D class averages'
        write(logfhandle,'(2a)') 'reprojs.mrc: stack of 2D reprojections of 3D volume along '//&
                                &'same orientations as cavgs.mrc', NEW_LINE('a')
        stop
    else
        call get_command_argument(1, fn_cavgs)
        call get_command_argument(2, fn_reprojs)
    end if

    ! Read in stacks
    call find_ldim_nptcls(fn_cavgs, ldim1, ifoo)
    call find_ldim_nptcls(fn_reprojs, ldim2, ifoo)
    if (ldim1(1) /= ldim2(1) .or. ldim1(2) /= ldim2(2)) then
        print *, "CAVGS AND REPROJS OF DIFFERENT DIMENSIONS"
        stop
    else if (ldim1(3) /= ldim2(3)) then
        print *, "CAVGS AND REPROJS OF DIFFERENT NUMBER"
        stop
    end if
    
    ! Allocations and parameter calculations
    ldim_refs = [ldim1(1), ldim1(2), 1]
    nshells = nprad_vox / shell_size_vox
    center = ldim_refs(1:2) / 2 + 0.5
    allocate(mask(ldim_refs(1), ldim_refs(2), ldim_refs(3)), source=.false.)
    call cavg%new(ldim_refs, smpd)
    call reproj%new(ldim_refs, smpd)
    call diff%new(ldim_refs, smpd)

    ! Results CSV File
    call fopen(funit, FILE=trim(fn_results), STATUS='REPLACE', action='WRITE')
    write(funit, '(A)') results_header
    601 format(F12.4,A2)
    602 format(F12.4)

    iref = 1 ! Index counting only populated class averages
    do icavg=1,ldim1(3)
        call cavg%read(fn_cavgs, icavg)
        ! If cavg is empty then ignore
        minmax = cavg%minmax()
        if (minmax(1) == 0. .and. minmax(2) == 0.) then
            call cavg%zero()
            call reproj%zero()
            cycle
        end if
        call reproj%read(fn_reprojs, icavg)

        ! Apply circular mask prior to normalization
        call cavg%mask(mskdiam_vox, mask_type)
        call reproj%mask(mskdiam_vox, mask_type)

        ! Normalize for comparison
        if (normalize_real_space) then
            ! Low-pass filter at resolution in reconstruction (1Å)
            call cavg%fft()
            call reproj%fft()
            call cavg%lp(calc_fourier_index(res, ldim_refs(1), smpd))
            call reproj%lp(calc_fourier_index(res, ldim_refs(1), smpd))
            call cavg%ifft()
            call reproj%ifft()
            call cavg%norm()
            call reproj%norm()
        else
            call match_fourier_variance(cavg, reproj)
        end if
        call reproj%mask(mskdiam_vox, mask_type)

        ! Signal Subtraction
        diff = cavg - reproj

        ! Analysis should inlcude low-pass filter if normalized in Fourier space
        if (.not. normalize_real_space) then
            call diff%fft()
            call diff%lp(calc_fourier_index(res, ldim_refs(1), smpd))
            call diff%ifft()
        end if

        ! Take averages of shells out to the NP radius
        call diff%get_rmat_ptr(rmat)
        call cavg%get_rmat_ptr(rmat_cavg)
        call reproj%get_rmat_ptr(rmat_reproj)
        mean = 0. 
        do n=0, nshells
            mask = .false. ! For now use mask in case we want to calculate other stats in the future
            do i=1, ldim_refs(1)
                do j=1, ldim_refs(2)
                    r = sqrt(real((i - center(1))**2 + (j - center(2))**2))
                    if (r > n*shell_size_vox .and. r < (n+1)*shell_size_vox) then
                        mask(i,j,1) = .true.
                    end if
                end do
            end do
            !mean = sum(abs(rmat), mask=mask) / sum(abs(reproj%get_rmat()), mask=mask)
            ! Mean fractional difference since the signal is higher at the core than at the surface
            mean = 2 * sum(abs(rmat), mask=mask) / (sum(abs(rmat_reproj), mask=mask) + sum(abs(rmat_cavg), mask=mask))
            ! Write csv file containing statistics
            write(funit,601,advance='no') real(iref-1),                CSV_DELIM ! INDEX
            write(funit,601,advance='no') n*shell_size_A,              CSV_DELIM ! RMIN
            write(funit,601,advance='no') (n+1)*shell_size_A,          CSV_DELIM ! RMAX
            write(funit,601,advance='no') real(count(mask)),           CSV_DELIM ! COUNT
            write(funit,601,advance='no') mean,                        CSV_DELIM ! MEAN_DIFF
            write(funit,601,advance='no') minval(rmat, mask=mask),     CSV_DELIM ! MIN_DIFF
            write(funit,602)              maxval(rmat, mask=mask)                ! MAX_DIFF
            ! Write images
            call cavg%write(fn_cavgs_out, iref)
            call reproj%write(fn_reprojs_out, iref)
            call diff%write(fn_diff, iref)
        end do
        ! Cleanup
        call cavg%zero()
        call reproj%zero()
        call diff%zero()
        iref = iref + 1
    end do
    call fclose(funit)
    call cavg%kill()
    call reproj%kill()
    call diff%kill()
    deallocate(mask)

contains

    ! Normalizes img2 to match the normalizaiton of img1 in preparation for subtraction
    ! by matching the energy of img2 to img1 over rings in Fourier space.
    ! Adapted from E. Fernandez-Gimenez, et. al., 2021, doi: https://doi.org/10.1016/j.jsb.2021.107780
    subroutine match_fourier_variance(img1, img2)
        type(image), intent(in)                :: img1
        type(image), intent(inout)             :: img2
        type(image)                            :: img1_copy, img2_orig_ft
        complex(kind=c_float_complex), pointer :: cmat1(:,:,:)=>null(), cmat2(:,:,:)=>null(), cmat2_orig(:,:,:)=>null()
        complex                                :: phase
        real(kind=c_float),            pointer :: rmat2(:,:,:)=>null()
        real, allocatable                      :: img1_shavgs(:), img2_shavgs(:)
        real                                   :: minmax(2)
        integer, allocatable        :: sh_counts(:)
        integer, parameter          :: niter=1
        integer                     :: iter, i, j, h, k, sh, lims(3,2), phys(3), ldim(3), array_shape(3), nshells

        if (.not. img1%exists()) then
            write(logfhandle,'(a)')  'match_fourier_variance: img1 does not exist!'
            stop
        end if
        if (.not. img2%exists()) then
            write(logfhandle,'(a)')  'match_fourier_variance: img2 does not exist!'
            stop
        end if
        ldim = img1%get_ldim()

        call img1_copy%copy(img1)    ! Need to preserve img1
        if (.not. img1_copy%is_ft()) call img1_copy%fft()
        if (.not. img2%is_ft()) call img2%fft()
        call img2_orig_ft%copy(img2) ! Save phase information for later

        ! Get matrix pointers
        call img1_copy%get_cmat_ptr(cmat1)
        call img2%get_cmat_ptr(cmat2)
        call img2%get_rmat_ptr(rmat2)
        call img2_orig_ft%get_cmat_ptr(cmat2_orig)
        
        do iter=1,niter

            ! Find average energy of shells in Fourier space
            lims = img1%loop_lims(2)
            nshells = 1 + nint( hyp( max( abs(lims(1,1)),lims(1,2) ), max( abs(lims(2,1)),lims(2,2) ) ) )
            allocate(sh_counts(nshells), source=0)
            allocate(img1_shavgs(nshells), img2_shavgs(nshells), source=0.)
            do h=lims(1,1), lims(1,2)
                do k=lims(2,1), lims(2,2)
                    phys = img2%comp_addr_phys(h,k,0)  ! compute physical address
                    sh = nint(hyp(h,k)) + 1 ! find shell (+1 so that sh is never 0)
                    img1_shavgs(sh) = ( img1_shavgs(sh)*sh_counts(sh) + abs(cmat1(phys(1),phys(2),phys(3))) ) / (sh_counts(sh) + 1)
                    img2_shavgs(sh) = ( img2_shavgs(sh)*sh_counts(sh) + abs(cmat2(phys(1),phys(2),phys(3))) ) / (sh_counts(sh) + 1)
                    sh_counts(sh) = sh_counts(sh) + 1
                end do
            end do
            
            ! Match energy by shell
            do h=lims(1,1), lims(1,2)
                do k=lims(2,1), lims(2,2)
                    phys = img2%comp_addr_phys(h,k,0)
                    sh = nint(hyp(h,k)) + 1
                    cmat2(phys(1),phys(2),phys(3)) = img1_shavgs(sh)/img2_shavgs(sh)*cmat2(phys(1),phys(2),phys(3))
                    sh_counts(sh) = sh_counts(sh) + 1
                end do
            end do
            deallocate(sh_counts, img1_shavgs, img2_shavgs)

            ! ! Match greyscale range of img2 to img1
            ! call img2%ifft()
            ! minmax = img1%minmax()
            ! do i=1, ldim(1)
            !     do j=1, ldim(2)
            !         if (rmat2(i,j,1) > minmax(2)) then
            !             rmat2(i,j,1) = minmax(2)
            !         else if (rmat2(i,j,1) < minmax(1)) then
            !             rmat2(i,j,1) = minmax(1)
            !         end if
            !     end do
            ! end do
            ! call img2%fft()

            ! Preserve phase information of img2
            array_shape = img2%get_array_shape()
            do i=1, array_shape(1)
                do j=1, array_shape(2)
                    phase = cmat2_orig(i,j,1) / abs(cmat2_orig(i,j,1))
                    cmat2(i,j,1) = abs(cmat2(i,j,1))*phase
                end do
            end do
        end do

        call img2%ifft()
        ! Cleanup
        call img1_copy%kill()
        call img2_orig_ft%kill()
    end subroutine match_fourier_variance

end program simple_test_2D_core_finder

