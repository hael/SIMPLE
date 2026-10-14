!@descr: for distributed sigma2 calculations in objfun=euclid 2D and 3D refinement
module simple_commanders_euclid_distr
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use, intrinsic :: iso_fortran_env, only: int64, real32
use simple_commanders_api
use simple_ptcl_layout, only: ptcl_layout_digest
use simple_sigma2_pspec_sums, only: pspec_part_sums_fname, read_pspec_part_sums, bootstrap_seed_spectra, sum_valid_records
use simple_sigma2_state, only: sigma2_state_candidate_path
use simple_sigma2_state_file, only: sigma2_state_header, sigma2_state_create_candidate, &
    &sigma2_state_init_header, sigma2_state_read_header, sigma2_state_read_particles, &
    &sigma2_state_write_particles, sigma2_state_write_global, sigma2_state_validate_file, &
    &sigma2_state_publish, SIGMA2_PROV_PSPEC, SIGMA2_STATE_FNAME
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_calc_pspec_assemble
  contains
    procedure :: execute      => exec_calc_pspec_assemble
end type commander_calc_pspec_assemble

contains

    subroutine exec_calc_pspec_assemble( self, cline )
        class(commander_calc_pspec_assemble), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(parameters)                 :: params
        type(image)                      :: avg_img
        type(builder)                    :: build
        type(sigma2_state_header)        :: state_header, previous_header
        type(string)                     :: part_fname
        type(string)                     :: state_path,candidate_path
        integer(int64)                   :: generation,layout_digest,prefix_digest
        integer(int64)                   :: global_counts(2), ninvalid
        integer                          :: iptcl,ipart,nptcls,nptcls_sel,eo,nyq,pspec_l,pspec_u
        integer                          :: first_row,last_row,nprefix,nshell,half,i
        integer                          :: state_status, nactive(2), part_nimgs, nimgs
        real(dp),          allocatable   :: seeds(:,:), half_sums(:,:), part_sums(:,:), global_sums(:,:)
        real,              allocatable   :: pspec_ave(:)
        real(real32),      allocatable   :: block(:,:), seed_records(:,:), global_means(:,:)
        integer,           allocatable   :: eo_ids(:)
        logical,           allocatable   :: pspec_covered(:), active(:)
        logical                          :: state_path_found, preserve_prefix
        character(len=STDLEN)            :: state_message
        call cline%set('mkdir', 'no')
        call cline%set('stream','no')
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.false.)
        ! set Fourier index range
        params%kfromto(1) = 1
        params%kfromto(2) = calc_fourier_index(2.*params%smpd, params%box, params%smpd)
        ! the partitions' Fourier sums of the sampled images, for the average image
        nptcls     = build%spproj_field%get_noris(consider_state=.false.)
        nptcls_sel = build%spproj_field%count_state_gt_zero()
        if( nptcls_sel < 1 ) THROW_HARD('No active particles found; exec_calc_pspec_assemble')
        call avg_img%new([params%box,params%box,1], params%smpd)
        call avg_img%zero_and_flag_ft
        do ipart = 1,params%nparts
            call build%img%zero_and_flag_ft
            part_fname = 'sum_img_part'//int2str_pad(ipart,params%numlen)//params%ext%to_char()
            call build%img%read(part_fname)
            call avg_img%add(build%img)
            call del_file(part_fname)
        enddo
        nyq = avg_img%get_nyq()
        ! the partitions' per-half sums of the sampled spectra, over the particle ranges they cover, and the number
        ! of images their Fourier sums hold
        allocate(half_sums(nyq,2), source=0.d0)
        allocate(pspec_covered(params%nptcls), source=.false.)
        nimgs = 0
        do ipart = 1,params%nparts
            part_fname = pspec_part_sums_fname(ipart)
            call read_pspec_part_sums(part_fname%to_char(), nyq, pspec_l, pspec_u, part_nimgs, part_sums, &
                &state_status, state_message)
            if( state_status /= 0 ) THROW_HARD(trim(state_message))
            if( (pspec_l<1).or.(pspec_u>params%nptcls) )then
                write(logfhandle,*) 'file/ptcl range/nptcls: ', part_fname%to_char(), pspec_l, pspec_u, params%nptcls
                THROW_HARD('commander_euclid; exec_calc_pspec_assemble; invalid particle spectra range')
            end if
            if( any(pspec_covered(pspec_l:pspec_u)) )then
                THROW_HARD('commander_euclid; exec_calc_pspec_assemble; overlapping particle spectra ranges')
            endif
            pspec_covered(pspec_l:pspec_u) = .true.
            if( .not. all(ieee_is_finite(part_sums)) ) THROW_HARD('a partition sum of sampled spectra is not finite; exec_calc_pspec_assemble')
            half_sums = half_sums + part_sums
            nimgs     = nimgs + part_nimgs
            deallocate(part_sums)
        end do
        ! the average image is the mean of the images summed (sampled, with a valid spectrum); with none, the
        ! sum stays zero and so does its power
        if( nimgs > 0 ) call avg_img%div(real(nimgs))
        call avg_img%spectrum('power',pspec_ave,norm=.true.)
        pspec_ave = pspec_ave / 2.0
        call avg_img%kill
        ! every active particle is covered; the halves' sizes are those of the project's active set
        nactive = 0
        do iptcl = 1,nptcls
            if( build%spproj_field%get_state(iptcl) <= 0 ) cycle
            if( .not.pspec_covered(iptcl) )then
                write(logfhandle,*) 'iptcl: ', iptcl
                THROW_HARD('active particle sigma spectrum was not covered; exec_calc_pspec_assemble')
            endif
            eo = build%spproj_field%get_eo(iptcl)
            if( eo < 0 .or. eo > 1 ) THROW_HARD('active particle has invalid even/odd assignment; exec_calc_pspec_assemble')
            nactive(eo+1) = nactive(eo+1) + 1
        enddo
        allocate(seeds(nyq,2))
        call bootstrap_seed_spectra(half_sums, nactive, pspec_ave, seeds)
        if( trim(params%oritype) /= 'ptcl2D' .and. trim(params%oritype) /= 'ptcl3D' ) &
            &THROW_HARD('canonical sigma2 bootstrap requires a particle orientation field')
        call build%spproj%get_sigma2_state_path(state_path, state_path_found)
        if( .not. state_path_found )then
            ! execution-local state: registered by name, it lives next to
            ! the project file that owns it
            call build%spproj%set_sigma2_state_path(string(SIGMA2_STATE_FNAME))
            call build%spproj%get_sigma2_state_path(state_path, state_path_found)
        endif
        generation = 1_int64
        preserve_prefix = .false.
        if( file_exists(state_path) )then
            call sigma2_state_read_header(state_path%to_char(), previous_header, state_status, state_message)
            if( state_status == 0 )then
                generation = previous_header%generation + 1_int64
                if( previous_header%nptcls > 0 .and. previous_header%nptcls < nptcls )then
                    prefix_digest = ptcl_layout_digest(build%spproj, build%spproj_field, &
                        &int(previous_header%nptcls))
                    preserve_prefix = prefix_digest == previous_header%layout_digest .and. &
                        &previous_header%box == params%box .and. &
                        &previous_header%kfrom == params%kfromto(1) .and. &
                        &previous_header%kto == params%kfromto(2) .and. &
                        &abs(real(previous_header%smpd)-params%smpd) <= 1.e-5*max(1.,params%smpd)
                    if( preserve_prefix )then
                        call sigma2_state_validate_file(state_path%to_char(), state_status, state_message)
                        preserve_prefix = state_status == 0
                    endif
                endif
            endif
        endif
        layout_digest = ptcl_layout_digest(build%spproj, build%spproj_field)
        if( layout_digest == 0_int64 ) THROW_HARD('cannot derive canonical sigma2 particle layout identity')
        ! one global spectrum per half, bootstrapped from a global even/odd sample
        call sigma2_state_init_header(state_header, params%kfromto(1), params%kfromto(2), &
            &nptcls, params%box, params%smpd, generation, &
            &layout_digest, SIGMA2_PROV_PSPEC)
        candidate_path = sigma2_state_candidate_path(state_path%to_char(), generation)
        call sigma2_state_create_candidate(candidate_path%to_char(), state_header, state_status, state_message)
        if( state_status /= 0 ) THROW_HARD(trim(state_message))
        ! the records in blocks: the retained committed prefix, then every active row of a half with its seed (the
        ! inactive rows stay zero); the global spectra are the means of the valid active records, summed while
        ! written, in row order as a reduction would read them
        nshell = params%kfromto(2) - params%kfromto(1) + 1
        allocate(seed_records(nshell,2))
        do half = 1, 2
            seed_records(:,half) = real(seeds(params%kfromto(1):params%kfromto(2),half), real32)
            if( .not.(all(ieee_is_finite(seed_records(:,half))) .and. all(seed_records(:,half) > 0.0_real32)) ) &
                &THROW_HARD('the bootstrap seed spectrum is not finite and positive; exec_calc_pspec_assemble')
        enddo
        allocate(active(nptcls), eo_ids(nptcls))
        do iptcl = 1, nptcls
            active(iptcl) = build%spproj_field%get_state(iptcl) > 0
            eo_ids(iptcl) = build%spproj_field%get_eo(iptcl)
        enddo
        allocate(global_sums(nshell,2), source=0.d0)
        global_counts = 0_int64
        ninvalid      = 0_int64
        nprefix       = 0
        if( preserve_prefix ) nprefix = int(previous_header%nptcls)
        do first_row = 1, nprefix, 4096
            last_row = min(first_row+4095, nprefix)
            call sigma2_state_read_particles(state_path%to_char(), first_row, last_row, &
                &block, state_status, state_message)
            if( state_status /= 0 ) THROW_HARD(trim(state_message))
            call sum_valid_records(block, first_row, active, eo_ids, global_sums, global_counts, ninvalid)
            call sigma2_state_write_particles(candidate_path%to_char(), first_row, block, state_status, state_message)
            if( state_status /= 0 ) THROW_HARD(trim(state_message))
            deallocate(block)
        enddo
        if( nprefix > 0 ) write(logfhandle,'(A,I0,A)') '>>> SIGMA2 APPEND: retained ', nprefix, &
            &' committed particle spectra and bootstrapped the appended suffix'
        do first_row = nprefix + 1, nptcls, 4096
            last_row = min(first_row+4095, nptcls)
            allocate(block(nshell,last_row-first_row+1), source=0.0_real32)
            do i = 1, size(block,2)
                iptcl = first_row + i - 1
                if( active(iptcl) ) block(:,i) = seed_records(:,eo_ids(iptcl)+1)
            enddo
            call sum_valid_records(block, first_row, active, eo_ids, global_sums, global_counts, ninvalid)
            call sigma2_state_write_particles(candidate_path%to_char(), first_row, block, state_status, state_message)
            if( state_status /= 0 ) THROW_HARD(trim(state_message))
            deallocate(block)
        enddo
        if( ninvalid > 0_int64 ) write(logfhandle,'(A,I0,A)') &
            &'>>> WARNING: canonical sigma2 reduction skipped ', ninvalid, &
            &' active particle record(s) with non-finite or non-positive values'
        if( any(global_counts == 0_int64) ) THROW_HARD('canonical sigma2 has an empty even/odd half; exec_calc_pspec_assemble')
        allocate(global_means(nshell,2))
        do half = 1, 2
            global_means(:,half) = real(global_sums(:,half)/real(global_counts(half),dp), real32)
        enddo
        call sigma2_state_write_global(candidate_path%to_char(), global_means, state_status, state_message)
        if( state_status /= 0 ) THROW_HARD(trim(state_message))
        call sigma2_state_publish(candidate_path%to_char(), state_path%to_char(), state_status, state_message)
        if( state_status /= 0 ) THROW_HARD(trim(state_message))
        ! the partitions' sums stay until the seed is published, so a failed assembly can be rerun
        do ipart = 1,params%nparts
            call del_file(pspec_part_sums_fname(ipart))
        end do
        call build%spproj%set_sigma2_state_path(state_path)
        call build%spproj%write_segment_inside('projinfo', params%projfile)
        ! end gracefully
        deallocate(seeds, half_sums, pspec_covered, seed_records, global_sums, global_means, active, eo_ids)
        if( allocated(pspec_ave) ) deallocate(pspec_ave)
        call state_path%kill
        call candidate_path%kill
        call build%kill_general_tbox
        call simple_touch('CALC_PSPEC_FINISHED')
        call simple_end('**** SIMPLE_CALC_PSPEC_ASSEMBLE NORMAL STOP ****', print_simple=.false.)

    end subroutine exec_calc_pspec_assemble

end module simple_commanders_euclid_distr
