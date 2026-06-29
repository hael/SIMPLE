!@descr: shared types for class-average compatibility analysis
module simple_cavg_compatibility_analysis
use unix,                only: c_float
use simple_defs,         only: logfhandle, nthr_glob, SHC_INPL_TRSHWDTH
use simple_string,       only: string
use simple_image,        only: image
use simple_error,        only: simple_exception
use simple_image_bin,    only: image_bin
use simple_imgarr_utils, only: read_cavgs_into_imgarr, write_imgarr, dealloc_imgarr
use simple_segmentation, only: otsu_img
use simple_math,         only: otsu
use simple_math_ft,      only: calc_fourier_index
use simple_srch_sort_loc, only: hpsort
use simple_parameters,   only: parameters
use simple_cmdline,      only: cmdline
use simple_type_defs,    only: OBJFUN_CC
use simple_builder,      only: builder
use simple_sp_project,   only: sp_project
use simple_defs_fname,   only: METADATA_EXT
use simple_commanders_cavgs, only: commander_cluster_cavgs

implicit none

public :: cavg_compatibility_analysis, run_cluster_overfitting_rejection

private
#include "simple_local_flags.inc"

integer, parameter :: ANALYSIS_BOXSIZE    = 128
integer, parameter :: ANALYSIS_MORPH_SIZE = 5
logical, parameter :: ANALYSIS_AUTOTUNE_SIZE_MODEL = .true.

type :: image_pointset
    type(image)          :: img
    type(image_bin)      :: mask
    integer, allocatable :: pts(:,:)
    integer              :: npts             = 0
    integer              :: rejection_reason = 0
    real                 :: thr              = 0.0
    real                 :: variance         = 0.0
    real                 :: sum_in_mask      = 0.0
    real                 :: local_var_in_mask  = 0.0
    real                 :: local_var_out_mask = 0.0
    real                 :: mean_hausdorff   = 0.0
    logical              :: is_rejected      = .false.
    logical              :: is_compatible    = .false.
end type image_pointset

type :: cavg_compatibility_analysis
    private
    type(sp_project)                  :: spproj
    type(image_pointset), allocatable :: pointsets(:)
    integer,              allocatable :: hausdorff_tbl(:,:)
    type(string)                      :: input_stkname
    integer                           :: npointsets = 0
  contains
    procedure :: new
    procedure :: kill
    procedure :: analyse
    procedure :: print_rejection_reasons
    procedure :: generate_pointset
    procedure :: infer_compatible_size_subset
    procedure :: write_selected_rejected_stacks
    procedure :: run_cluster_overfitting_rejection
    
end type cavg_compatibility_analysis

contains

    subroutine new( self, spproj )
        class(cavg_compatibility_analysis), intent(inout) :: self
        type(sp_project),                      intent(in) :: spproj
        type(image),                          allocatable :: imgs(:), img_out(:), mask_out(:)
        type(string)                                      :: out_imgs, out_masks
        integer                                           :: iimg, ncls
        real                                              :: smpd
        call self%kill()
        self%spproj = spproj
        call self%spproj%get_cavgs_stk(self%input_stkname, ncls, smpd)
        imgs            = read_cavgs_into_imgarr(self%spproj)
        self%npointsets = size(imgs)
        if( self%npointsets < 1 ) THROW_HARD('simple_test_hausdorff: no images found in input stack')
        allocate(self%pointsets(self%npointsets), img_out(self%npointsets), mask_out(self%npointsets))
        do iimg = 1, self%npointsets
            call imgs(iimg)%set_smpd(smpd)
            call self%generate_pointset(iimg, imgs(iimg))
            call img_out(iimg)%copy(self%pointsets(iimg)%img)
            if( self%pointsets(iimg)%is_rejected )then
                call mask_out(iimg)%copy(self%pointsets(iimg)%img)
                call mask_out(iimg)%zero()
            else
                call mask_out(iimg)%copy(self%pointsets(iimg)%mask)
            end if
        end do

        out_imgs  = self%input_stkname//'.cavg_compatibility_imgs.mrcs'
        out_masks = self%input_stkname//'.cavg_compatibility_masks.mrcs'
        call write_imgarr(img_out,   out_imgs)
        call write_imgarr(mask_out, out_masks)

        call dealloc_imgarr(img_out)
        call dealloc_imgarr(mask_out)
        call dealloc_imgarr(imgs)
    end subroutine new

    subroutine kill( self )
        class(cavg_compatibility_analysis), intent(inout) :: self
        integer :: i
        if( allocated(self%pointsets) )then
            do i = 1, size(self%pointsets)
                call self%pointsets(i)%img%kill()
                call self%pointsets(i)%mask%kill_bimg()
                if( allocated(self%pointsets(i)%pts) ) deallocate(self%pointsets(i)%pts)
            end do
            deallocate(self%pointsets)
        end if
        if( allocated(self%hausdorff_tbl) ) deallocate(self%hausdorff_tbl)
        call self%input_stkname%kill()
        call self%spproj%kill()
        self%npointsets = 0
    end subroutine kill

    subroutine analyse( self )
        class(cavg_compatibility_analysis), intent(inout) :: self
        real, allocatable :: similarity_mat(:,:)
        integer, allocatable :: sim_idx_map(:)
        real :: hp_use, lp_use, trs_use
        integer :: i, nrej, ncomp

        write(logfhandle,'(A,I0)') 'cavg_compatibility analyse: npointsets=', self%npointsets

        call self%run_cluster_overfitting_rejection(5000.0)
        nrej = 0
        do i = 1, self%npointsets
            if( self%pointsets(i)%is_rejected ) nrej = nrej + 1
        end do
        write(logfhandle,'(A,I0)') 'after cluster overfitting rejection: nrejected=', nrej

        call self%infer_compatible_size_subset()
        nrej = 0
        do i = 1, self%npointsets
            if( self%pointsets(i)%is_rejected ) nrej = nrej + 1
        end do
        write(logfhandle,'(A,I0)') 'after infer_compatible_size_subset: nrejected=', nrej

        call self%write_selected_rejected_stacks()

        call self%print_rejection_reasons()
        write(logfhandle,'(A)') 'rejection reasons written to rejection_reasons.txt'

        if( allocated(similarity_mat) ) deallocate(similarity_mat)
        if( allocated(sim_idx_map) ) deallocate(sim_idx_map)
        

    end subroutine analyse

    subroutine write_selected_rejected_stacks( self )
        class(cavg_compatibility_analysis), intent(inout) :: self
        type(image), allocatable :: selected_imgs(:), rejected_imgs(:)
        type(string)             :: selected_out, rejected_out
        integer                  :: i, nsel, nrej, isel, irej

        nsel = 0
        nrej = 0
        do i = 1, self%npointsets
            if( self%pointsets(i)%is_compatible ) nsel = nsel + 1
            if( self%pointsets(i)%is_rejected ) nrej = nrej + 1
        end do

        selected_out = self%input_stkname//'.selected_cavgs.mrcs'
        rejected_out = self%input_stkname//'.rejected_cavgs.mrcs'

        if( nsel > 0 )then
            allocate(selected_imgs(nsel))
            isel = 0
            do i = 1, self%npointsets
                if( .not. self%pointsets(i)%is_compatible ) cycle
                isel = isel + 1
                call selected_imgs(isel)%copy(self%pointsets(i)%img)
            end do
            call write_imgarr(selected_imgs, selected_out)
            call dealloc_imgarr(selected_imgs)
        end if

        if( nrej > 0 )then
            allocate(rejected_imgs(nrej))
            irej = 0
            do i = 1, self%npointsets
                if( .not. self%pointsets(i)%is_rejected ) cycle
                irej = irej + 1
                call rejected_imgs(irej)%copy(self%pointsets(i)%img)
            end do
            call write_imgarr(rejected_imgs, rejected_out)
            call dealloc_imgarr(rejected_imgs)
        end if

        write(logfhandle,'(A,I0,A,A)') 'selected cavg stack count=', nsel, ' file=', selected_out%to_char()
        write(logfhandle,'(A,I0,A,A)') 'rejected cavg stack count=', nrej, ' file=', rejected_out%to_char()
    end subroutine write_selected_rejected_stacks
   
    subroutine infer_compatible_size_subset( self )
        class(cavg_compatibility_analysis), intent(inout) :: self
        logical, allocatable :: allowed(:), compatible(:)
        real,    allocatable :: min_dim(:), max_dim(:)
        integer :: n, ii, nallowed, ncomp, nrejected_here
        integer :: ldim(3)
        real    :: dmax, dmin, axis_c, axis_b, axis_a, used_relax, used_qlo, used_qhi
        logical :: do_autotune
        real(kind=c_float), pointer :: rmat_mask(:,:,:) => null()

        n = self%npointsets
        if( n < 1 ) return

        do ii = 1, n
            self%pointsets(ii)%is_compatible = .false.
        end do

        allocate(allowed(n), source=.false.)
        allocate(min_dim(n), max_dim(n), source=0.0)
        nallowed = 0
        do ii = 1, n
            allowed(ii) = .not. self%pointsets(ii)%is_rejected
            if( .not. allowed(ii) ) cycle
            nallowed = nallowed + 1

            ldim = self%pointsets(ii)%mask%get_ldim()
            call self%pointsets(ii)%mask%get_rmat_ptr(rmat_mask)

            dmax = estimate_max_feret(rmat_mask, ldim, self%pointsets(ii)%mask%get_smpd())
            dmin = estimate_min_feret(rmat_mask, ldim, self%pointsets(ii)%mask%get_smpd())
            nullify(rmat_mask)

            max_dim(ii) = dmax
            min_dim(ii) = dmin
            write(logfhandle,'(A,I0,A,ES14.6,A,ES14.6)') 'infer_compatible_size_subset dims: idx=', ii, &
                ' min_dim=', min_dim(ii), ' max_dim=', max_dim(ii)
        end do
        write(logfhandle,'(A,I0)') 'infer_compatible_size_subset: allowed candidates=', nallowed

        if( nallowed <= 2 )then
            do ii = 1, n
                self%pointsets(ii)%is_compatible = allowed(ii)
            end do
            deallocate(allowed, min_dim, max_dim)
            return
        end if
        do_autotune = ANALYSIS_AUTOTUNE_SIZE_MODEL
        call find_reprojection_support(min_dim, max_dim, allowed, axis_c, axis_b, axis_a, compatible, do_autotune, used_relax, used_qlo, used_qhi)
        ncomp = count(compatible)
        write(logfhandle,'(A,ES14.6,A,ES14.6,A,ES14.6,A,I0)') 'infer_compatible_size_subset: latent_axes c/b/a=', &
            axis_c, '/', axis_b, '/', axis_a, ' compatible_count=', ncomp
        write(logfhandle,'(A,L1,A,ES14.6,A,ES14.6,A,ES14.6)') 'infer_compatible_size_subset: autotune=', do_autotune, &
            ' relax=', used_relax, ' qlo=', used_qlo, ' qhi=', used_qhi

        nrejected_here = 0
        do ii = 1, n
            self%pointsets(ii)%is_compatible = allowed(ii) .and. compatible(ii)
            if( allowed(ii) .and. .not. compatible(ii) )then
                self%pointsets(ii)%is_rejected = .true.
                self%pointsets(ii)%rejection_reason = 11
                self%pointsets(ii)%is_compatible = .false.
                nrejected_here = nrejected_here + 1
                write(logfhandle,'(A,I0,A)') 'infer_compatible_size_subset: rejecting idx=', ii, ' reason=size_incompatible_subset'
            end if
        end do
        write(logfhandle,'(A,I0,A,I0)') 'infer_compatible_size_subset: compatible count=', ncomp, ' newly_rejected=', nrejected_here

        deallocate(allowed, compatible, min_dim, max_dim)

    contains

        real function estimate_max_feret(mask_rmat, ldim_local, smpd_local) result(max_feret)
            real(kind=c_float), intent(in) :: mask_rmat(:,:,:)
            integer,            intent(in) :: ldim_local(3)
            real,               intent(in) :: smpd_local
            integer, parameter :: NFERET = 180
            integer :: npix, x, y, ip, itheta
            real, allocatable :: xpts(:), ypts(:)
            real :: theta, cth, sth, proj, pmin, pmax, width, pi_v

            npix = count(mask_rmat(1:ldim_local(1), 1:ldim_local(2), 1) > 0.5_c_float)
            if( npix <= 0 )then
                max_feret = 0.0
                return
            end if

            allocate(xpts(npix), ypts(npix))
            ip = 0
            do y = 1, ldim_local(2)
                do x = 1, ldim_local(1)
                    if( mask_rmat(x,y,1) > 0.5_c_float )then
                        ip = ip + 1
                        xpts(ip) = real(x) * smpd_local
                        ypts(ip) = real(y) * smpd_local
                    end if
                end do
            end do

            pi_v = acos(-1.0)
            max_feret = 0.0
            do itheta = 0, NFERET - 1
                theta = pi_v * real(itheta) / real(NFERET)
                cth = cos(theta)
                sth = sin(theta)
                pmin = huge(1.0)
                pmax = -huge(1.0)
                do ip = 1, npix
                    proj = xpts(ip) * cth + ypts(ip) * sth
                    if( proj < pmin ) pmin = proj
                    if( proj > pmax ) pmax = proj
                end do
                width = (pmax - pmin) + smpd_local
                if( width > max_feret ) max_feret = width
            end do

            deallocate(xpts, ypts)
        end function estimate_max_feret

        real function estimate_min_feret(mask_rmat, ldim_local, smpd_local) result(min_feret)
            real(kind=c_float), intent(in) :: mask_rmat(:,:,:)
            integer,            intent(in) :: ldim_local(3)
            real,               intent(in) :: smpd_local
            integer, parameter :: NFERET = 180
            integer :: npix, x, y, ip, itheta
            real, allocatable :: xpts(:), ypts(:)
            real :: theta, cth, sth, proj, pmin, pmax, width, pi_v

            npix = count(mask_rmat(1:ldim_local(1), 1:ldim_local(2), 1) > 0.5_c_float)
            if( npix <= 0 )then
                min_feret = 0.0
                return
            end if

            allocate(xpts(npix), ypts(npix))
            ip = 0
            do y = 1, ldim_local(2)
                do x = 1, ldim_local(1)
                    if( mask_rmat(x,y,1) > 0.5_c_float )then
                        ip = ip + 1
                        xpts(ip) = real(x) * smpd_local
                        ypts(ip) = real(y) * smpd_local
                    end if
                end do
            end do

            pi_v = acos(-1.0)
            min_feret = huge(1.0)
            do itheta = 0, NFERET - 1
                theta = pi_v * real(itheta) / real(NFERET)
                cth = cos(theta)
                sth = sin(theta)
                pmin = huge(1.0)
                pmax = -huge(1.0)
                do ip = 1, npix
                    proj = xpts(ip) * cth + ypts(ip) * sth
                    if( proj < pmin ) pmin = proj
                    if( proj > pmax ) pmax = proj
                end do
                width = (pmax - pmin) + smpd_local
                if( width < min_feret ) min_feret = width
            end do

            deallocate(xpts, ypts)
        end function estimate_min_feret

        subroutine find_reprojection_support(mins, maxs, allowed_mask, c_out, b_out, a_out, support_mask, autotune, relax_out, qlo_out, qhi_out)
            real,               intent(in)  :: mins(:), maxs(:)
            logical,            intent(in)  :: allowed_mask(:)
            real,               intent(out) :: c_out, b_out, a_out
            logical, allocatable, intent(out) :: support_mask(:)
            logical,            intent(in)  :: autotune
            real,               intent(out) :: relax_out, qlo_out, qhi_out
            integer, parameter :: NRELAX = 5, NQ = 3
            real, parameter :: RELAX_GRID(NRELAX) = [0.05, 0.07, 0.10, 0.13, 0.18]
            real, parameter :: QLOW_GRID(NQ) = [0.03, 0.08, 0.12]
            real, parameter :: QHIGH_GRID(NQ) = [0.88, 0.93, 0.97]
            integer :: nn, i, j, k, nvalid, best_count
            integer :: ir, iq, jq, final_count, best_final_count
            real, allocatable :: cands(:), mins_sup(:), maxs_sup(:)
            real    :: b_cand, left, right
            real    :: spread, best_spread, score, best_score
            real    :: relax_cur, qlo_cur, qhi_cur
            real    :: c_cur, b_cur, a_cur
            real    :: c_best, b_best, a_best
            logical, allocatable :: tmp_support(:), support_cur(:), best_support(:)

            nn = size(mins)
            allocate(support_mask(nn), source=.false.)
            allocate(tmp_support(nn), source=.false.)
            allocate(support_cur(nn), source=.false.)
            allocate(best_support(nn), source=.false.)
            c_out = 0.0
            b_out = 0.0
            a_out = 0.0
            relax_out = 0.07
            qlo_out = 0.08
            qhi_out = 0.93

            nvalid = 0
            do i = 1, nn
                if( allowed_mask(i) ) nvalid = nvalid + 1
            end do
            if( nvalid <= 0 )then
                deallocate(tmp_support, support_cur, best_support)
                return
            end if

            allocate(cands(2*nvalid))
            k = 0
            do i = 1, nn
                if( .not. allowed_mask(i) ) cycle
                k = k + 1
                cands(k) = mins(i)
                k = k + 1
                cands(k) = maxs(i)
            end do

            best_count = -1
            best_spread = huge(1.0)
            best_final_count = -1
            best_score = -huge(1.0)
            c_best = 0.0
            b_best = 0.0
            a_best = 0.0

            do ir = 1, NRELAX
                if( autotune )then
                    relax_cur = RELAX_GRID(ir)
                else
                    if( ir > 1 ) cycle
                    relax_cur = 0.07
                end if

                do iq = 1, NQ
                    if( autotune )then
                        qlo_cur = QLOW_GRID(iq)
                    else
                        if( iq > 1 ) cycle
                        qlo_cur = 0.08
                    end if

                    do jq = 1, NQ
                        if( autotune )then
                            qhi_cur = QHIGH_GRID(jq)
                        else
                            if( jq > 1 ) cycle
                            qhi_cur = 0.93
                        end if
                        if( qhi_cur <= qlo_cur ) cycle

                        best_count = -1
                        best_spread = huge(1.0)
                        b_cur = 0.0
                        support_cur = .false.

                        do i = 1, k
                            b_cand = cands(i)
                            final_count = 0
                            spread = 0.0
                            tmp_support = .false.

                            do j = 1, nn
                                if( .not. allowed_mask(j) ) cycle
                                left = mins(j) * (1.0 - relax_cur)
                                right = maxs(j) * (1.0 + relax_cur)
                                if( b_cand >= left .and. b_cand <= right )then
                                    final_count = final_count + 1
                                    tmp_support(j) = .true.
                                    spread = spread + abs(b_cand - 0.5 * (mins(j) + maxs(j))) / max(maxs(j)-mins(j), 1.0e-6)
                                end if
                            end do

                            if( final_count > best_count .or. (final_count == best_count .and. spread < best_spread) )then
                                best_count = final_count
                                best_spread = spread
                                b_cur = b_cand
                                support_cur = tmp_support
                            end if
                        end do

                        if( best_count <= 0 ) cycle

                        allocate(mins_sup(best_count), maxs_sup(best_count))
                        j = 0
                        do i = 1, nn
                            if( .not. support_cur(i) ) cycle
                            j = j + 1
                            mins_sup(j) = mins(i)
                            maxs_sup(j) = maxs(i)
                        end do

                        c_cur = percentile_real(mins_sup, qlo_cur)
                        a_cur = percentile_real(maxs_sup, qhi_cur)
                        if( c_cur > b_cur ) c_cur = minval(mins_sup)
                        if( a_cur < b_cur ) a_cur = maxval(maxs_sup)

                        final_count = 0
                        spread = 0.0
                        tmp_support = .false.
                        do i = 1, nn
                            if( .not. allowed_mask(i) ) cycle
                            left = mins(i) * (1.0 - relax_cur)
                            right = maxs(i) * (1.0 + relax_cur)
                            if( b_cur >= left .and. b_cur <= right .and. mins(i) >= c_cur * (1.0 - relax_cur) .and. maxs(i) <= a_cur * (1.0 + relax_cur) )then
                                final_count = final_count + 1
                                tmp_support(i) = .true.
                                spread = spread + abs(b_cur - 0.5 * (mins(i) + maxs(i))) / max(maxs(i)-mins(i), 1.0e-6)
                            end if
                        end do

                        score = real(final_count) - 0.01 * spread - 0.10 * relax_cur
                        if( score > best_score .or. (abs(score - best_score) <= 1.0e-6 .and. final_count > best_final_count) )then
                            best_score = score
                            best_final_count = final_count
                            c_best = c_cur
                            b_best = b_cur
                            a_best = a_cur
                            relax_out = relax_cur
                            qlo_out = qlo_cur
                            qhi_out = qhi_cur
                            best_support = tmp_support
                        end if

                        deallocate(mins_sup, maxs_sup)
                    end do
                end do
            end do

            if( best_final_count <= 0 )then
                deallocate(cands, tmp_support, support_cur, best_support)
                return
            end if

            c_out = c_best
            b_out = b_best
            a_out = a_best
            support_mask = best_support

            deallocate(cands, tmp_support, support_cur, best_support)
        end subroutine find_reprojection_support

        real function percentile_real(arr, q) result(pval)
            real, intent(in) :: arr(:)
            real, intent(in) :: q
            real, allocatable :: sorted(:)
            integer :: nvals, idx

            nvals = size(arr)
            if( nvals <= 0 )then
                pval = 0.0
                return
            end if

            allocate(sorted(nvals))
            sorted = arr
            call hpsort(sorted)

            idx = 1 + int(q * real(max(nvals - 1, 0)))
            idx = max(1, min(nvals, idx))
            pval = sorted(idx)
            deallocate(sorted)
        end function percentile_real

    end subroutine infer_compatible_size_subset

    subroutine print_rejection_reasons( self )
        class(cavg_compatibility_analysis), intent(inout) :: self
        integer :: i, funit
        character(len=64) :: reason_txt

        open(newunit=funit, file='rejection_reasons.txt', status='replace', action='write')
        write(funit,'(A)') '# class_index reason_code reason_text'
        write(logfhandle,'(A)') 'Rejected class averages:'

        do i = 1, self%npointsets
            if( .not. self%pointsets(i)%is_rejected ) cycle
            select case(self%pointsets(i)%rejection_reason)
            case(1)
                reason_txt = 'sum_in_mask_low_outlier'
            case(2)
                reason_txt = 'sum_in_mask_high_outlier'
            case(3)
                reason_txt = 'zero_variance'
            case(4)
                reason_txt = 'mask_outside_circular_support'
            case(5)
                reason_txt = 'high_hausdorff_outlier'
            case(6)
                reason_txt = 'not_in_compatible_subset'
            case(7)
                reason_txt = 'low_in_local_variance'
            case(8)
                reason_txt = 'low_out_local_variance'
            case(9)
                reason_txt = 'low_dynamic_range_otsu'
            case(10)
                reason_txt = 'dynamic_range_outlier_pre_otsu'
            case(11)
                reason_txt = 'size_incompatible_subset'
            case(12)
                reason_txt = 'haralick_texture_outlier'
            case(13)
                reason_txt = 'suspected_overfitting'
            case default
                reason_txt = 'unknown'
            end select
            write(funit,'(I0,1X,I0,1X,A)') i, self%pointsets(i)%rejection_reason, trim(reason_txt)
            write(logfhandle,'(A,I0,A,I0,A,A)') '  idx=', i, ' code=', self%pointsets(i)%rejection_reason, ' reason=', trim(reason_txt)
        end do

        close(funit)
    end subroutine print_rejection_reasons

    subroutine generate_pointset( self, i, img )
        class(cavg_compatibility_analysis), intent(inout) :: self
        type(image),                        intent(inout) :: img
        integer,                            intent(in)    :: i
        logical,                            allocatable   :: l_circ(:,:,:)
        real(kind=c_float),                 pointer       :: rmat_mask(:,:,:) => null(), rmat_src(:,:,:) => null()
        type(image) :: circ_mask
        integer     :: ldim(3)
        integer     :: ldim_target(3)
        integer     :: imorph, k, x, y

        ! Reset per-image outputs so repeated calls on the same slot are safe.
        if( allocated(self%pointsets(i)%pts) ) deallocate(self%pointsets(i)%pts)
        self%pointsets(i)%npts        = 0
        self%pointsets(i)%rejection_reason = 0
        self%pointsets(i)%sum_in_mask = 0.0
        self%pointsets(i)%local_var_in_mask = 0.0
        self%pointsets(i)%local_var_out_mask = 0.0
        self%pointsets(i)%mean_hausdorff = 0.0
        self%pointsets(i)%thr         = 0.0

        ! Preprocess the input image: zero edge average and bandpass filter.
        call img%zero_edgeavg()
        call img%bp(0., 10.)

        ! Frequency-domain resize to 128x128: FFT -> crop/pad -> inverse FFT.
        call img%fft()
        ldim        = img%get_ldim()
        ldim_target = [ANALYSIS_BOXSIZE, ANALYSIS_BOXSIZE, ldim(3)]
        if( ldim(1) > ldim_target(1) .or. ldim(2) > ldim_target(2) )then
            call img%clip_inplace(ldim_target)
        else if( ldim(1) < ldim_target(1) .or. ldim(2) < ldim_target(2) )then
            call img%pad_inplace(ldim_target)
        end if
        call img%ifft()
        call img%set_smpd(img%get_smpd() * real(ldim(1)) / real(ldim_target(1)))

        ! Keep a processed copy of the resized image for downstream analysis/output.
        self%pointsets(i)%img           = img
        self%pointsets(i)%variance      = img%variance()
        self%pointsets(i)%is_rejected   = self%pointsets(i)%variance == 0.0
        self%pointsets(i)%is_compatible = .false.
        if( self%pointsets(i)%is_rejected ) then
            self%pointsets(i)%rejection_reason = 3
            return
        end if
        ! Compute Otsu threshold and binary mask.
        call otsu_img(img, thresh=self%pointsets(i)%thr)
        ! Apply 5 px morphological closing to the binary Otsu mask.
        call self%pointsets(i)%mask%transfer2bimg(img)
        do imorph = 1, ANALYSIS_MORPH_SIZE
            call self%pointsets(i)%mask%dilate()
        end do
        do imorph = 1, ANALYSIS_MORPH_SIZE
            call self%pointsets(i)%mask%erode()
        end do
        ! Build a hard circular support mask with diameter equal to the box size.
        ldim = self%pointsets(i)%mask%get_ldim()
        call circ_mask%disc(ldim, self%pointsets(i)%mask%get_smpd(), 0.5 * real(min(ldim(1), ldim(2))), l_circ)
        call self%pointsets(i)%mask%get_rmat_ptr(rmat_mask)
        ! Any foreground outside support marks this class average as rejected.
        if( any(rmat_mask(1:ldim(1),1:ldim(2),1:ldim(3)) > 0.5 .and. .not. l_circ(1:ldim(1),1:ldim(2),1:ldim(3))) ) then
            self%pointsets(i)%is_rejected = .true.
            self%pointsets(i)%rejection_reason = 4
            ! Early exit: free temporary support resources before returning.
            nullify(rmat_mask)
            if( allocated(l_circ) ) deallocate(l_circ)
            call circ_mask%kill()
            return
        end if
        if( allocated(l_circ) ) deallocate(l_circ)
        call circ_mask%kill()

        ! Extract foreground coordinates from the final binary mask.
        self%pointsets(i)%npts = count(rmat_mask(1:ldim(1), 1:ldim(2), 1) > 0.5)
        if( self%pointsets(i)%npts == 0 )then
            allocate(self%pointsets(i)%pts(2,1), source=0)
        else
            allocate(self%pointsets(i)%pts(2,self%pointsets(i)%npts))
            k = 0
            do y = 1, ldim(2)
                do x = 1, ldim(1)
                    if( rmat_mask(x,y,1) > 0.5 )then
                        k = k + 1
                        self%pointsets(i)%pts(1,k) = x
                        self%pointsets(i)%pts(2,k) = y
                    end if
                end do
            end do
        end if
        call self%pointsets(i)%img%get_rmat_ptr(rmat_src)
        ! Compute source-image intensity sum under the binary foreground mask.
        self%pointsets(i)%sum_in_mask = sum( &
            real(rmat_src(1:ldim(1), 1:ldim(2), 1), kind=kind(self%pointsets(i)%sum_in_mask)),&
            mask=rmat_mask(1:ldim(1), 1:ldim(2), 1) > 0.5&
        )

        call self%pointsets(i)%img%loc_var_masked(rmat_mask(1:ldim(1), 1:ldim(2), 1), 10, &
            self%pointsets(i)%local_var_in_mask, self%pointsets(i)%local_var_out_mask)

        write(logfhandle,'(A,I0,A,ES14.6,A,ES14.6)') 'generate_pointset: idx=', i, ' var_in=', self%pointsets(i)%local_var_in_mask, &
            ' var_out=', self%pointsets(i)%local_var_out_mask
        nullify(rmat_mask)
        nullify(rmat_src)
        
    end subroutine generate_pointset

    subroutine run_cluster_overfitting_rejection( self, mskdiam, nclust_in )
        class(cavg_compatibility_analysis), intent(inout) :: self
        real,                     parameter  :: LOWVAR_CLUSTER_REJECT_FRAC = 0.60
        integer,                  parameter  :: REJECT_REASON_SUSPECTED_OVERFITTING = 13
        real,                                  intent(in) :: mskdiam
        integer,                     optional, intent(in) :: nclust_in
        integer,                              allocatable :: states(:), clusters(:)
        integer,                              allocatable :: per_cluster_lowvar(:), per_cluster_total(:)
        logical,                              allocatable :: reject_cluster(:)
        real,                                 allocatable :: vals_in(:), vals_out(:), absdev(:)
        type(string)                                      :: projname
        type(cmdline)                                     :: cline
        type(sp_project)                                  :: spproj
        type(commander_cluster_cavgs)                     :: commander
        integer                                           :: i, nclustered, nclusters, k, iclust
        real                                              :: thr_in, thr_out, med_in_all, med_out_all, mad_in_all, mad_out_all
        real                                              :: vin, vout, frac_lowvar
        projname = 'cluster_cavgs'//METADATA_EXT

        allocate(states(self%npointsets), source=1)
        do i = 1, self%npointsets
            if( self%pointsets(i)%is_rejected ) states(i) = 0
        end do
        call self%spproj%map_cavgs_selection(states)
        deallocate(states)
        call self%spproj%write(projname)
        call cline%set('projfile', projname)
        call cline%set('mskdiam', mskdiam)
        if( present(nclust_in) ) call cline%set('ncls', nclust_in)
        call commander%execute(cline)
        call spproj%read(projname)
        clusters = spproj%os_cls2D%get_all_asint('cluster')
        if( size(clusters) == self%npointsets )then
            nclustered = count(clusters > 0)
            nclusters  = maxval(clusters)
            write(logfhandle,'(A,I0,A,I0)') 'run_exec_cluster_cavgs: clustered classes=', nclustered, ' total=', self%npointsets

            k = 0
            do i = 1, self%npointsets
                if( self%pointsets(i)%is_rejected ) cycle
                k = k + 1
            end do
            if( k > 0 .and. nclusters > 0 )then
                allocate(vals_in(k), vals_out(k))
                k = 0
                do i = 1, self%npointsets
                    if( self%pointsets(i)%is_rejected ) cycle
                    k = k + 1
                    vals_in(k)  = log(max(self%pointsets(i)%local_var_in_mask,  1.0e-8))
                    vals_out(k) = log(max(self%pointsets(i)%local_var_out_mask, 1.0e-8))
                end do

                med_in_all = median_real_local(vals_in)
                med_out_all = median_real_local(vals_out)
                allocate(absdev(size(vals_in)))
                absdev = abs(vals_in - med_in_all)
                mad_in_all = median_real_local(absdev)
                deallocate(absdev)
                allocate(absdev(size(vals_out)))
                absdev = abs(vals_out - med_out_all)
                mad_out_all = median_real_local(absdev)
                deallocate(absdev)

                thr_in  = med_in_all  - 0.5 * max(mad_in_all,  1.0e-3)
                thr_out = med_out_all - 0.5 * max(mad_out_all, 1.0e-3)

                allocate(per_cluster_lowvar(nclusters), source=0)
                allocate(per_cluster_total(nclusters), source=0)
                allocate(reject_cluster(nclusters), source=.false.)
                do i = 1, self%npointsets
                    if( clusters(i) < 1 .or. clusters(i) > nclusters ) cycle
                    per_cluster_total(clusters(i)) = per_cluster_total(clusters(i)) + 1
                    if( self%pointsets(i)%is_rejected ) cycle
                    vin = log(max(self%pointsets(i)%local_var_in_mask,  1.0e-8))
                    vout = log(max(self%pointsets(i)%local_var_out_mask, 1.0e-8))
                    if( vin <= thr_in .and. vout <= thr_out )then
                        per_cluster_lowvar(clusters(i)) = per_cluster_lowvar(clusters(i)) + 1
                    end if
                end do

                write(logfhandle,'(A,ES14.6,A,ES14.6)') 'run_exec_cluster_cavgs: localvar thresholds in/out=', thr_in, ' / ', thr_out
                do iclust = 1, nclusters
                    if( per_cluster_total(iclust) > 0 )then
                        frac_lowvar = real(per_cluster_lowvar(iclust)) / real(per_cluster_total(iclust))
                    else
                        frac_lowvar = 0.0
                    end if
                    if( frac_lowvar > LOWVAR_CLUSTER_REJECT_FRAC ) reject_cluster(iclust) = .true.
                    write(logfhandle,'(A,I0,A,I0,A,I0,A,F6.2,A,L1)') 'run_exec_cluster_cavgs: cluster=', iclust, &
                        ' lowvar_in_and_out=', per_cluster_lowvar(iclust), ' total=', per_cluster_total(iclust), &
                        ' lowvar_pct=', 100.0 * frac_lowvar, ' reject_cluster=', reject_cluster(iclust)
                end do

                allocate(states(self%npointsets), source=1)
                do i = 1, self%npointsets
                    if( clusters(i) >= 1 .and. clusters(i) <= nclusters )then
                        if( reject_cluster(clusters(i)) )then
                            self%pointsets(i)%is_rejected = .true.
                            if( self%pointsets(i)%rejection_reason == 0 ) then
                                self%pointsets(i)%rejection_reason = REJECT_REASON_SUSPECTED_OVERFITTING
                            end if
                        end if
                    end if
                    if( self%pointsets(i)%is_rejected ) states(i) = 0
                end do

                deallocate(states)

                deallocate(reject_cluster, per_cluster_lowvar, per_cluster_total, vals_in, vals_out)
            end if
            
        else
            write(logfhandle,'(A,I0,A,I0)') 'run_exec_cluster_cavgs: cluster size mismatch cls2D=', size(clusters), ' pointsets=', self%npointsets
        end if
        if( allocated(clusters) ) deallocate(clusters)
        
        call spproj%kill()
        call cline%kill()

    contains

        function median_real_local(arr) result(med_val)
            real, intent(in) :: arr(:)
            real :: med_val
            real, allocatable :: sorted(:)
            integer :: n, mid
            if( size(arr) == 0 )then
                med_val = 0.0
                return
            end if
            n = size(arr)
            mid = (n + 1) / 2
            allocate(sorted(n))
            sorted = arr
            call hpsort(sorted)
            if( mod(n,2) == 1 )then
                med_val = sorted(mid)
            else
                med_val = 0.5 * (sorted(mid) + sorted(mid+1))
            end if
            deallocate(sorted)
        end function median_real_local

    end subroutine run_cluster_overfitting_rejection

end module simple_cavg_compatibility_analysis
