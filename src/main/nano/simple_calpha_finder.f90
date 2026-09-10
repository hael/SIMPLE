!@descr: Buccaneer-inspired oriented target detection of alpha carbons in cryo-EM maps
module simple_calpha_finder
use simple_core_module_api
use simple_image, only: image
use simple_atoms, only: atoms
implicit none

public :: calpha_finder
private
#include "simple_local_flags.inc"

real, parameter :: PI_LOCAL         = 3.14159265358979323846
! Ideal peptide-backbone geometry in Angstroms.
real, parameter :: N_CA_BOND        = 1.458
real, parameter :: CA_C_BOND        = 1.525
real, parameter :: N_CA_C_ANGLE     = 111.2
real, parameter :: MIN_TARGET_SIGMA = 0.85

type :: calpha_finder
    private
    real                            :: smpd      = 0.
    real                            :: radius    = 0.
    integer                         :: halfwidth = 0
    integer                         :: nside     = 0
    real(kind=c_float), allocatable :: target_mean(:,:,:)
    real(kind=c_float), allocatable :: target_weight(:,:,:)
    logical                         :: existence = .false.
contains
    procedure :: new
    procedure :: search
    procedure :: kill
end type calpha_finder

contains

    subroutine new(self, smpd, radius)
        class(calpha_finder), intent(inout) :: self
        real,                 intent(in)    :: smpd, radius
        real    :: atom_sites(3,3), amplitudes(3), offset(3), delta(3)
        real    :: bond_angle, distance, sigma, value
        integer :: iatom, ix, iy, iz

        if(smpd   <= 0.) THROW_HARD('Target sampling distance must be positive; calpha_finder%new')
        if(radius <= 0.) THROW_HARD('Target radius must be positive; calpha_finder%new')

        call self%kill()
        self%smpd      = smpd
        self%radius    = radius
        self%halfwidth = ceiling(radius / self%smpd)
        self%nside     = 2 * self%halfwidth + 1
        allocate(self%target_mean(self%nside,self%nside,self%nside), source=0.)
        allocate(self%target_weight(self%nside,self%nside,self%nside), source=0.)

        ! Canonical frame: CA at the origin, CA-C on +x, and N in the xy plane.
        bond_angle = N_CA_C_ANGLE * PI_LOCAL / 180.
        atom_sites(:,1) = [0.,0.,0.]
        atom_sites(:,2) = [N_CA_BOND*cos(bond_angle), N_CA_BOND*sin(bond_angle), 0.]
        atom_sites(:,3) = [CA_C_BOND,0.,0.]
        amplitudes      = [1.25,1.0,1.0]
        ! Keep the expected atom density resolved on coarser input grids.
        sigma = max(MIN_TARGET_SIGMA, 0.75 * smpd)
        do iz = 1, self%nside
            do iy = 1, self%nside
                do ix = 1, self%nside
                    offset   = real([ix,iy,iz] - (self%halfwidth + 1)) * self%smpd
                    distance = sqrt(sum(offset * offset))
                    if(distance > radius) cycle
                    value = 0.
                    do iatom = 1, size(atom_sites,2)
                        delta = offset - atom_sites(:,iatom)
                        value = value + amplitudes(iatom) * &
                            exp(-0.5 * sum(delta * delta) / (sigma * sigma))
                    enddo
                    self%target_mean(ix,iy,iz) = value
                    self%target_weight(ix,iy,iz) = 0.5 * &
                        (1. + cos(PI_LOCAL * distance / radius))
                enddo
            enddo
        enddo
        self%existence = .true.
        write(logfhandle,'(A,I0,A,F5.2,A)') 'Constructed analytic C-alpha target on a ', &
            self%nside, '-voxel grid with sigma ', sigma, ' A'
    end subroutine new

    subroutine search( self, workvol, angstep, npeaks, score_threshold, pdbout, scorevol )
        class(calpha_finder), intent(in)    :: self
        class(image),         intent(inout) :: workvol
        real,                 intent(in)    :: angstep
        integer,              intent(in)    :: npeaks
        real,                 intent(in)    :: score_threshold
        class(string),        intent(in)    :: pdbout, scorevol
        type(image) :: work_ft, work_sq_ft, weighted_template, support_template
        type(image) :: data_target_corr, data_sum_corr, data_sq_sum_corr, best_scores
        real(kind=c_float), pointer :: workptr(:,:,:), sqptr(:,:,:), bestptr(:,:,:)
        real,    allocatable :: rotations(:,:,:)
        integer, allocatable :: best_rotation(:,:,:)
        real(dp) :: sum_w, sum_t, sum_tt
        integer  :: ldim(3), irot, nrot

        if(.not.self%existence)  THROW_HARD('C-alpha finder is not initialized; search')
        if(.not.workvol%is_3d()) THROW_HARD('Search input must be a 3D map; search')
        if(workvol%is_ft())      THROW_HARD('Search input must be in real space; search')
        if(abs(workvol%get_smpd() - self%smpd) > 1.e-4) &
            THROW_HARD('Target and search map must have the same sampling distance')
        if(any(mod(workvol%get_ldim(),2) /= 0)) &
            THROW_HARD('C-alpha FFT search currently requires even map dimensions')
        if(angstep <= 0. .or. angstep > 180.) THROW_HARD('angstep must be in (0,180]')
        if(npeaks <= 0) THROW_HARD('npeaks must be positive')

        ldim = workvol%get_ldim()
        if(any(ldim <= 2 * (ceiling(self%radius / self%smpd) + 1))) &
            THROW_HARD('Search map is too small for the C-alpha target')
        call build_rotation_grid(angstep, rotations)
        nrot = size(rotations,3)
        allocate(best_rotation(ldim(1),ldim(2),ldim(3)), source=0)
        call best_scores%new(ldim, self%smpd)
        call best_scores%get_rmat_ptr(bestptr)
        bestptr(1:ldim(1),1:ldim(2),1:ldim(3)) = -1.

        work_ft    = workvol
        call work_ft%fft()
        work_sq_ft = workvol
        call work_sq_ft%get_rmat_ptr(sqptr)
        call workvol%get_rmat_ptr(workptr)
        sqptr(1:ldim(1),1:ldim(2),1:ldim(3)) = &
            workptr(1:ldim(1),1:ldim(2),1:ldim(3))**2
        call work_sq_ft%fft()
        call weighted_template%new(ldim, self%smpd)
        call support_template%new(ldim, self%smpd)
        call data_target_corr%new(ldim, self%smpd)
        call data_sum_corr%new(ldim, self%smpd)
        call data_sq_sum_corr%new(ldim, self%smpd)

        write(logfhandle,'(A,I0,A,F6.1,A)') 'Searching ', nrot, &
            ' coarse orientations at ', angstep, ' degree spacing'
        do irot = 1, nrot
            call fill_rotated_target(self, rotations(:,:,irot), weighted_template, &
                support_template, sum_w, sum_t, sum_tt)
            call work_ft%ccf_into(weighted_template, data_target_corr)
            call work_ft%ccf_into(support_template, data_sum_corr)
            call work_sq_ft%ccf_into(support_template, data_sq_sum_corr)
            call update_scores(ldim, sum_w, sum_t, sum_tt, data_target_corr, &
                data_sum_corr, data_sq_sum_corr, irot, best_scores, best_rotation)
            if(mod(irot, max(1,nrot/10)) == 0 .or. irot == nrot) &
                write(logfhandle,'(A,I0,A,I0)') 'C-alpha orientations: ', irot, '/', nrot
        enddo

        call zero_score_border(self, best_scores)
        call best_scores%write(scorevol)
        call write_candidates(self, best_scores, best_rotation, rotations, npeaks, &
            score_threshold, pdbout)
        call weighted_template%kill()
        call support_template%kill()
        call data_target_corr%kill()
        call data_sum_corr%kill()
        call data_sq_sum_corr%kill()
        call work_ft%kill()
        call work_sq_ft%kill()
        call best_scores%kill()
        deallocate(rotations, best_rotation)
    end subroutine search

    subroutine interpolation_cell(pos, ldim, base, frac, valid)
        real,    intent(in)  :: pos(3)
        integer, intent(in)  :: ldim(3)
        integer, intent(out) :: base(3)
        real,    intent(out) :: frac(3)
        logical, intent(out) :: valid
        integer :: idim
        valid = .true.
        do idim = 1, 3
            if(pos(idim) < 1. .or. pos(idim) > real(ldim(idim)))then
                valid = .false.
                return
            endif
            if(pos(idim) >= real(ldim(idim)))then
                base(idim) = ldim(idim) - 1
                frac(idim) = 1.
            else
                base(idim) = floor(pos(idim))
                frac(idim) = pos(idim) - real(base(idim))
            endif
        enddo
    end subroutine interpolation_cell

    real(dp) function trilinear_value(array, base, frac)
        real(kind=c_float), intent(in) :: array(:,:,:)
        integer,            intent(in) :: base(3)
        real,               intent(in) :: frac(3)
        integer  :: dx, dy, dz
        real(dp) :: wx, wy, wz
        trilinear_value = 0._dp
        do dz = 0, 1
            wz = merge(real(frac(3),dp), 1._dp-real(frac(3),dp), dz == 1)
            do dy = 0, 1
                wy = merge(real(frac(2),dp), 1._dp-real(frac(2),dp), dy == 1)
                do dx = 0, 1
                    wx = merge(real(frac(1),dp), 1._dp-real(frac(1),dp), dx == 1)
                    trilinear_value = trilinear_value + wx * wy * wz * &
                        real(array(base(1)+dx,base(2)+dy,base(3)+dz),dp)
                enddo
            enddo
        enddo
    end function trilinear_value

    subroutine build_rotation_grid(angstep, rotations)
        real,              intent(in)  :: angstep
        real, allocatable, intent(out) :: rotations(:,:,:)
        real    :: delta, zaxis(3), uaxis(3), vaxis(3), xaxis(3), yaxis(3)
        real    :: z, phi, roll, radial
        integer :: naxes, nroll, iaxis, iroll, irot
        delta = angstep * PI_LOCAL / 180.
        naxes = max(1, ceiling(4. * PI_LOCAL / (delta * delta)))
        nroll = max(1, ceiling(2. * PI_LOCAL / delta))
        allocate(rotations(3,3,naxes*nroll))
        irot = 0
        do iaxis = 1, naxes
            z = 1. - 2. * (real(iaxis) - 0.5) / real(naxes)
            radial = sqrt(max(0., 1. - z*z))
            phi = PI_LOCAL * (3. - sqrt(5.)) * real(iaxis - 1)
            zaxis = [radial*cos(phi), radial*sin(phi), z]
            if(abs(zaxis(3)) < 0.9)then
                uaxis = cross([0.,0.,1.], zaxis)
            else
                uaxis = cross([0.,1.,0.], zaxis)
            endif
            uaxis = uaxis / sqrt(sum(uaxis*uaxis))
            vaxis = cross(zaxis, uaxis)
            do iroll = 1, nroll
                roll  = 2. * PI_LOCAL * real(iroll - 1) / real(nroll)
                xaxis               = cos(roll) * uaxis + sin(roll) * vaxis
                yaxis               = cross(zaxis, xaxis)
                irot                = irot + 1
                rotations(:,1,irot) = xaxis
                rotations(:,2,irot) = yaxis
                rotations(:,3,irot) = zaxis
            enddo
        enddo
    end subroutine build_rotation_grid

    subroutine fill_rotated_target( self, rotation, weighted_template, support_template, &
        sum_w, sum_t, sum_tt )
        class(calpha_finder), intent(in)    :: self
        real,                 intent(in)    :: rotation(3,3)
        class(image),         intent(inout) :: weighted_template, support_template
        real(dp),             intent(out)   :: sum_w, sum_t, sum_tt
        real(kind=c_float), pointer :: weighted(:,:,:), support(:,:,:)
        integer :: ldim(3), center(3), reach, ix, iy, iz
        real    :: map_offset(3), target_offset(3), target_value, weight
        logical :: valid

        call weighted_template%zero_and_unflag_ft()
        call support_template%zero_and_unflag_ft()
        call weighted_template%get_rmat_ptr(weighted)
        call support_template%get_rmat_ptr(support)
        ldim   = weighted_template%get_ldim()
        center = ldim / 2 + 1
        reach  = ceiling(self%radius / self%smpd) + 1
        sum_w  = 0._dp
        sum_t  = 0._dp
        sum_tt = 0._dp
        do iz = max(1,center(3)-reach), min(ldim(3),center(3)+reach)
            do iy = max(1,center(2)-reach), min(ldim(2),center(2)+reach)
                do ix = max(1,center(1)-reach), min(ldim(1),center(1)+reach)
                    map_offset         = real([ix,iy,iz] - center) * self%smpd
                    if(sum(map_offset * map_offset) > self%radius * self%radius) cycle
                    target_offset      = matmul(transpose(rotation), map_offset)
                    call sample_target(self, target_offset, target_value, weight, valid)
                    if(.not.valid .or. weight <= 0.) cycle
                    support(ix,iy,iz)  = weight
                    weighted(ix,iy,iz) = weight * target_value
                    sum_w              = sum_w  + real(weight,dp)
                    sum_t              = sum_t  + real(weight * target_value,dp)
                    sum_tt             = sum_tt + real(weight * target_value * target_value,dp)
                enddo
            enddo
        enddo
    end subroutine fill_rotated_target

    subroutine sample_target( self, xyz, value, weight, valid )
        class(calpha_finder), intent(in)  :: self
        real,                 intent(in)  :: xyz(3)
        real,                 intent(out) :: value, weight
        logical,              intent(out) :: valid
        real    :: pos(3), frac(3)
        integer :: base(3), ldim(3)
        ldim = [self%nside,self%nside,self%nside]
        pos  = xyz / self%smpd + real(self%halfwidth + 1)
        call interpolation_cell(pos, ldim, base, frac, valid)
        if(.not.valid)then
            value  = 0.
            weight = 0.
            return
        endif
        value  = real(trilinear_value(self%target_mean, base, frac))
        weight = real(trilinear_value(self%target_weight, base, frac))
    end subroutine sample_target

    subroutine update_scores( ldim, sum_w, sum_t, sum_tt, data_target_corr, data_sum_corr, &
        data_sq_sum_corr, irot, best_scores, best_rotation )
        integer,        intent(in)    :: ldim(3), irot
        real(dp),       intent(in)    :: sum_w, sum_t, sum_tt
        class(image),   intent(inout) :: data_target_corr, data_sum_corr, data_sq_sum_corr
        class(image),   intent(inout) :: best_scores
        integer,        intent(inout) :: best_rotation(:,:,:)
        real(kind=c_float), pointer :: td(:,:,:), d(:,:,:), dd(:,:,:), best(:,:,:)
        real(dp) :: scale, target_var, data_var, numerator, score, data_sum, data_sq_sum
        integer  :: ix, iy, iz

        target_var = sum_tt - sum_t * sum_t / sum_w
        if(target_var <= epsilon(1._dp)) return
        scale = real(product(ldim),dp)
        call data_target_corr%get_rmat_ptr(td)
        call data_sum_corr%get_rmat_ptr(d)
        call data_sq_sum_corr%get_rmat_ptr(dd)
        call best_scores%get_rmat_ptr(best)
        !$omp parallel do collapse(3) default(shared) private(ix,iy,iz,data_sum,data_sq_sum,data_var,numerator,score)
        do iz = 1, ldim(3)
            do iy = 1, ldim(2)
                do ix = 1, ldim(1)
                    data_sum    = scale * real(d(ix,iy,iz),dp)
                    data_sq_sum = scale * real(dd(ix,iy,iz),dp)
                    data_var    = data_sq_sum - data_sum * data_sum / sum_w
                    if(data_var <= epsilon(1._dp)) cycle
                    numerator   = scale * real(td(ix,iy,iz),dp) - sum_t * data_sum / sum_w
                    score       = numerator / sqrt(target_var * data_var)
                    if(score > real(best(ix,iy,iz),dp))then
                        best(ix,iy,iz)          = real(score)
                        best_rotation(ix,iy,iz) = irot
                    endif
                enddo
            enddo
        enddo
        !$omp end parallel do
    end subroutine update_scores

    subroutine zero_score_border( self, scores )
        class(calpha_finder), intent(in)    :: self
        class(image),         intent(inout) :: scores
        real(kind=c_float), pointer :: values(:,:,:)
        integer :: ldim(3), border
        ldim   = scores%get_ldim()
        border = ceiling(self%radius / self%smpd) + 1
        call scores%get_rmat_ptr(values)
        values(1:border,:,:)                 = -1.
        values(ldim(1)-border+1:ldim(1),:,:) = -1.
        values(:,1:border,:)                 = -1.
        values(:,ldim(2)-border+1:ldim(2),:) = -1.
        values(:,:,1:border)                 = -1.
        values(:,:,ldim(3)-border+1:ldim(3)) = -1.
    end subroutine zero_score_border

    subroutine write_candidates( self, scores, best_rotation, rotations, npeaks, threshold, pdbout )
        class(calpha_finder), intent(in)    :: self
        class(image),         intent(inout) :: scores
        integer,              intent(in)    :: best_rotation(:,:,:), npeaks
        real,                 intent(in)    :: rotations(:,:,:), threshold
        class(string),        intent(in)    :: pdbout
        real(kind=c_float), pointer :: scoremat(:,:,:)
        real, allocatable    :: selection(:,:,:), eligible_scores(:), xyz(:,:), found_scores(:)
        real, allocatable    :: found_rotations(:,:,:)
        integer, allocatable :: eligible_indices(:)
        integer :: ldim(3), location(3), neligible, icandidate, nfound, irot
        integer :: ix, iy, iz, linear_index, plane_index

        ldim = scores%get_ldim()
        allocate(selection(ldim(1),ldim(2),ldim(3)))
        allocate(xyz(3,npeaks), found_scores(npeaks), found_rotations(3,3,npeaks))
        call scores%get_rmat_ptr(scoremat)
        selection    = scoremat(1:ldim(1),1:ldim(2),1:ldim(3))

        neligible    = count(selection >= threshold .and. best_rotation > 0)
        allocate(eligible_scores(neligible), eligible_indices(neligible))
        neligible    = 0
        linear_index = 0
        do iz = 1, ldim(3)
            do iy = 1, ldim(2)
                do ix = 1, ldim(1)
                    linear_index = linear_index + 1
                    if(selection(ix,iy,iz) < threshold .or. best_rotation(ix,iy,iz) <= 0) cycle
                    neligible = neligible + 1
                    eligible_scores(neligible)  = selection(ix,iy,iz)
                    eligible_indices(neligible) = linear_index
                enddo
            enddo
        enddo
        write(logfhandle,'(A,I0,A,F7.3)') 'Found ', neligible, &
            ' score voxels eligible at threshold ', threshold
        call hpsort(eligible_scores, eligible_indices)

        nfound = 0
        do icandidate = neligible, 1, -1
            linear_index  = eligible_indices(icandidate)
            location(3)   = (linear_index - 1) / (ldim(1) * ldim(2)) + 1
            plane_index   = mod(linear_index - 1, ldim(1) * ldim(2))
            location(2)   = plane_index / ldim(1) + 1
            location(1)   = mod(plane_index, ldim(1)) + 1
            if(selection(location(1),location(2),location(3)) < threshold) cycle
            irot          = best_rotation(location(1),location(2),location(3))
            if(irot <= 0) cycle
            nfound        = nfound + 1
            xyz(:,nfound) = real(location - 1) * self%smpd
            found_scores(nfound) = selection(location(1),location(2),location(3))
            found_rotations(:,:,nfound) = rotations(:,:,irot)
            call suppress_neighborhood(selection, location, 2.0 / self%smpd)
            if(nfound == npeaks) exit
        enddo
        call write_candidate_files(pdbout, xyz, found_scores, found_rotations, nfound)
        write(logfhandle,'(A,I0,A,F7.3)') 'Wrote ', nfound, &
            ' non-overlapping C-alpha candidates at score >= ', threshold
        if(nfound == npeaks) write(logfhandle,'(A,I0)') &
            'Candidate output reached the requested npeaks cap of ', npeaks
        deallocate(selection, eligible_scores, eligible_indices, xyz, found_scores, found_rotations)
    end subroutine write_candidates

    subroutine suppress_neighborhood( values, center, radius_voxels )
        real,    intent(inout) :: values(:,:,:)
        integer, intent(in)    :: center(3)
        real,    intent(in)    :: radius_voxels
        integer :: ix, iy, iz, reach
        reach = ceiling(radius_voxels)
        do iz = max(1,center(3)-reach), min(size(values,3),center(3)+reach)
            do iy = max(1,center(2)-reach), min(size(values,2),center(2)+reach)
                do ix = max(1,center(1)-reach), min(size(values,1),center(1)+reach)
                    if(sum(real([ix,iy,iz] - center)**2) <= radius_voxels**2) &
                        values(ix,iy,iz) = -huge(1.)
                enddo
            enddo
        enddo
    end subroutine suppress_neighborhood

    subroutine write_candidate_files( pdbout, xyz, scores, rotations, nfound )
        class(string), intent(in) :: pdbout
        real,          intent(in) :: xyz(:,:), scores(:), rotations(:,:,:)
        integer,       intent(in) :: nfound
        type(atoms)  :: candidates
        type(string) :: csvfile
        integer :: ipeak, funit, io_stat

        if(nfound > 0)then
            call candidates%new(nfound, .true.)
            do ipeak = 1, nfound
                call candidates%set_name(ipeak, 'CA  ')
                call candidates%set_element(ipeak, 'C ')
                call candidates%set_resname(ipeak, 'UNK')
                call candidates%set_chain(ipeak, 'A')
                call candidates%set_num(ipeak, ipeak)
                call candidates%set_resnum(ipeak, ipeak)
                call candidates%set_coord(ipeak, xyz(:,ipeak))
                call candidates%set_occupancy(ipeak, 1.)
                call candidates%set_beta(ipeak, scores(ipeak))
            enddo
            call candidates%writepdb(pdbout)
            call candidates%kill()
        else
            call fopen(funit, status='REPLACE', action='WRITE', file=pdbout, iostat=io_stat)
            call fileiochk('write_candidate_files; opening empty PDB', io_stat)
            call fclose(funit)
        endif

        csvfile = get_fbody(pdbout, 'pdb') // '.csv'
        call fopen(funit, status='REPLACE', action='WRITE', file=csvfile, iostat=io_stat)
        call fileiochk('write_candidate_files; opening CSV', io_stat)
        write(funit,'(A)') 'x,y,z,score,r11,r12,r13,r21,r22,r23,r31,r32,r33'
        do ipeak = 1, nfound
            write(funit,'(3(F12.5,","),F10.6,9(",",F12.7))') xyz(:,ipeak), scores(ipeak), &
                rotations(1,1,ipeak), rotations(1,2,ipeak), rotations(1,3,ipeak), &
                rotations(2,1,ipeak), rotations(2,2,ipeak), rotations(2,3,ipeak), &
                rotations(3,1,ipeak), rotations(3,2,ipeak), rotations(3,3,ipeak)
        enddo
        call fclose(funit)
    end subroutine write_candidate_files

    subroutine kill( self )
        class(calpha_finder), intent(inout) :: self
        if(allocated(self%target_mean)) deallocate(self%target_mean)
        if(allocated(self%target_weight)) deallocate(self%target_weight)
        self%smpd      = 0.
        self%radius    = 0.
        self%halfwidth = 0
        self%nside     = 0
        self%existence = .false.
    end subroutine kill

end module simple_calpha_finder
