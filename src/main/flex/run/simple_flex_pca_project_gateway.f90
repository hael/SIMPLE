!@descr: flex_pca project gateway: the one owner of what the run reads from and writes into the project
!!
!! Selection and validation of the run's particles (a worker's partition list, the master's state>0
!! rows), the even/odd repair the master persists for its workers, the canonical sigma2 load with
!! the master/worker pin, and the two deliveries that mutate the run's project copy: the per-state
!! weight store registered in the out segment, and the hard state labels written into ptcl3D.
module simple_flex_pca_project_gateway
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_builder,            only: builder
use simple_cmdline,            only: cmdline
use simple_parameters,         only: parameters
use simple_sp_project,         only: sp_project
use simple_sigma2_files,       only: load_sigma2_groups
use simple_estimate_ssnr,      only: fsc2optlp_sub
use simple_flex_pca_rounds,    only: flex_pca_rounds
use simple_flex_pca_run_types, only: flex_run_settings
use simple_flex_weights_state, only: flex_weights_deliver, flex_weights_state_fname, FLEX_WEIGHTS_STALE_SCAN
use simple_flex_weights_file,  only: FLEX_WEIGHTS_PROV_FLEX_PCA, FLEX_WEIGHTS_PROV_MERGED
implicit none
private
#include "simple_local_flags.inc"

public :: validate_covariance_inputs, load_and_validate_sigma
public :: write_flex_weights_store, write_discrete_state_project
public :: prepare_project_fsc_lowpass_filters, publish_state_volume, write_out_segment

character(len=*), parameter :: SIGMA_STATE_FNAME = 'flex_pca_sigma_state.txt'

contains

    !> Read a part's particle-index list (one integer per line, as written by the master's
    !! flex_pca_plan_partitions through arr2txtfile).
    subroutine read_pind_list( fname, pinds, nptcls )
        character(len=*),     intent(in)  :: fname
        integer, allocatable, intent(out) :: pinds(:)
        integer,              intent(out) :: nptcls
        type(string) :: fn
        integer :: funit, io_stat, i, n
        fn = trim(fname)
        if( .not. file_exists(fn) ) THROW_HARD('flex_pca worker: particle list not found: '//trim(fname))
        n = nlines(fn)
        if( n < 1 ) THROW_HARD('flex_pca worker: empty particle list: '//trim(fname))
        allocate(pinds(n))
        call fopen(funit, file=fn, action='READ', status='OLD', iostat=io_stat)
        call fileiochk('read_pind_list', io_stat)
        do i = 1, n
            read(funit,*) pinds(i)
        end do
        call fclose(funit)
        call fn%kill
        nptcls = n
    end subroutine read_pind_list

    subroutine validate_covariance_inputs( params, cfg, build, pinds, nptcls , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        type(parameters), intent(inout) :: params
        type(flex_run_settings), intent(in) :: cfg
        type(builder),    intent(inout) :: build
        integer, allocatable, intent(out) :: pinds(:)
        integer, intent(out) :: nptcls
        integer :: q, i, cnt
        integer, allocatable :: sel(:)
        if( trim(params%oritype) /= 'ptcl3D' ) THROW_HARD('flex_pca requires oritype=ptcl3D')
        if( .not. cfg%l_vol1_explicit )then
            THROW_HARD('flex_pca requires a consensus mean map: pass vol1 or register one in the project out segment')
        endif
        ! a run writes its state labels into its project copy, so a rerun must start from the original
        ! project (a subset is selected with pindfile=, never with the labels)
        if( build%spproj_field%get_n('state') /= 1 ) THROW_HARD('flex_pca works on one population but the project carries several state labels (a delivered copy): rerun from the original project, selecting particles with pindfile= if needed')
        if( cfg%l_pindfile )then
            ! a distributed worker takes the master's partition as it was planned: its own
            ! particle-index list, never a re-derivation of the selection from fromp/top
            call read_pind_list(params%pindfile%to_char(), pinds, nptcls)
            write(logfhandle,'(A,I0,A,A)') '>>> FLEX_PCA (WORKER) ',nptcls,' particles from ', &
                &params%pindfile%to_char()
            call flush(logfhandle)
        else
            ! not sample4rec: its updatecnt > 0 condition belongs to trailing reconstruction
            allocate(sel(max(0, params%top - params%fromp + 1)))
            cnt = 0
            do i = params%fromp, params%top
                if( build%spproj_field%get_state(i) > 0 )then
                    cnt = cnt + 1; sel(cnt) = i
                endif
            end do
            nptcls = cnt
            if( allocated(pinds) ) deallocate(pinds)
            allocate(pinds(nptcls), source=sel(1:nptcls))
            deallocate(sel)
        endif
        if( nptcls < 100 ) THROW_HARD('flex_pca requires at least 100 active particles')
        if( build%spproj%os_ptcl2D%get_noris() > 0 .and. &
            &build%spproj%os_ptcl2D%get_noris() /= build%spproj%os_ptcl3D%get_noris() ) &
            &THROW_HARD('flex_pca requires matching ptcl2D and ptcl3D rows')
        ! eo split by index parity when the project's is degenerate. MASTER-ONLY and must be persisted: a
        ! worker counts over its own range, so it would silently build a different -- equally valid -- split.
        if( rounds%is_worker() )then
            if( count([(build%spproj_field%get_eo(pinds(q))==0,q=1,nptcls)]) < 1 .or. &
                &count([(build%spproj_field%get_eo(pinds(q))==1,q=1,nptcls)]) < 1 )then
                THROW_HARD('flex_pca worker sees a degenerate eo split; the master did not persist its repair')
            endif
        else if( count([(build%spproj_field%get_eo(pinds(q))==0,q=1,nptcls)]) < 20 .or. &
            &count([(build%spproj_field%get_eo(pinds(q))==1,q=1,nptcls)]) < 20 )then
            write(logfhandle,'(A)') '>>> FLEX_PCA assigning alternating even/odd halfsets (project eo was degenerate)'
            call build%spproj_field%partition_eo
            if( rounds%is_master() )then
                call build%spproj%write_segment_inside(params%oritype, params%projfile)
                write(logfhandle,'(A)') '>>> FLEX_PCA persisted the repaired eo split for the workers'
            endif
        endif
        if( count([(build%spproj_field%get_eo(pinds(q))==0,q=1,nptcls)]) < 20 .or. &
            &count([(build%spproj_field%get_eo(pinds(q))==1,q=1,nptcls)]) < 20 ) &
            &THROW_HARD('flex_pca requires populated even and odd halfsets')
    end subroutine validate_covariance_inputs

    subroutine load_and_validate_sigma( params, build, cline, pinds, loaded , rounds)
        class(flex_pca_rounds), intent(inout) :: rounds
        type(parameters), intent(inout) :: params
        type(builder),    intent(inout) :: build
        class(cmdline),   intent(inout) :: cline
        integer,          intent(in)    :: pinds(:)
        logical,          intent(out)   :: loaded
        integer :: i, k, iptcl, noris, fromp_save, top_save
        ! The sigma table is allocated on params%fromp:top. A distributed worker's particle list is the
        ! master's partition (pindfile), not the fromp/top split of the project, so widen the range to
        ! every particle for the load (the canonical reader fills rows for state>0 within the range) and
        ! restore after. The unit fallback below has always done the same.
        noris      = build%spproj_field%get_noris()
        fromp_save = params%fromp
        top_save   = params%top
        params%fromp = 1
        params%top   = noris
        call load_sigma2_groups(params, build%pftc, build%esig, build%spproj, build%spproj_field, loaded)
        params%fromp = fromp_save
        params%top   = top_save
        if( .not. loaded )then
            ! gen_fplane4rec follows norm_noise_taper_edge_pad_fft, so a unit spectrum is the correct fallback
            ! Construct the sigma object through its own API (the legacy loader's idiom: widen
            ! fromp/top to every particle so the table covers the whole project, restore after):
            ! esig%new also registers the table with the polar calculator; a bare allocate of the
            ! component inside a never-constructed object leaves the consumers with an
            ! inconsistent object
            params%fromp = 1
            params%top   = noris
            call build%esig%new(params, build%pftc, string('flex_pca_unit_sigma2.dat'), params%box)
            params%fromp = fromp_save
            params%top   = top_save
            build%esig%sigma2_noise = 1.0
            loaded = .true.
            write(logfhandle,'(A)') '>>> FLEX_PCA SIGMA: using unit white-noise spectrum after background normalization'
            write(logfhandle,'(A)') '>>> FLEX_PCA SIGMA: provide sigma2_group*.star for coloured experimental noise'
        endif
        do i = 1, size(pinds)
            iptcl = pinds(i)
            do k = lbound(build%esig%sigma2_noise,1), ubound(build%esig%sigma2_noise,1)
                if( .not. ieee_is_finite(build%esig%sigma2_noise(k,iptcl)) .or. &
                    &build%esig%sigma2_noise(k,iptcl) <= TINY )then
                    THROW_HARD('flex_pca found a nonpositive or nonfinite sigma2 value')
                endif
            end do
        end do
        params%ml_reg   = 'yes'
        params%l_ml_reg = .true.
        call cline%set('ml_reg','yes')
        ! Pin the sigma2 decision across the master/worker boundary (see write_sigma_state).
        if( rounds%is_master() )then
            call write_sigma_state(loaded)
        else if( rounds%is_worker() )then
            call check_sigma_state(loaded)
        endif
    end subroutine load_and_validate_sigma

    !> Whether the master resolved real sigma2 spectra or fell back to the unit white-noise spectrum.
    !! Discovery is re-run independently in every worker, and a worker that resolves it differently
    !! whitens against a different noise model -- which changes every inner product without any error
    !! being raised, because both outcomes are individually legitimate. Pin it and make workers assert.
    subroutine write_sigma_state( loaded )
        logical, intent(in) :: loaded
        integer :: funit, io_stat
        call fopen(funit, file=string(SIGMA_STATE_FNAME), action='WRITE', status='REPLACE', iostat=io_stat)
        call fileiochk('write_sigma_state; open', io_stat)
        write(funit,'(L1)') loaded
        call fclose(funit)
    end subroutine write_sigma_state

    subroutine check_sigma_state( loaded )
        logical, intent(in) :: loaded
        integer :: funit, io_stat
        logical :: loaded_master
        if( .not. file_exists(string(SIGMA_STATE_FNAME)) )then
            THROW_HARD('flex_pca worker found no flex_pca_sigma_state.txt from the master')
        endif
        call fopen(funit, file=string(SIGMA_STATE_FNAME), action='READ', status='OLD', iostat=io_stat)
        call fileiochk('check_sigma_state; open', io_stat)
        read(funit,'(L1)') loaded_master
        call fclose(funit)
        if( loaded_master .neqv. loaded )then
            THROW_HARD('flex_pca worker resolved sigma2 differently from the master; the parts would be whitened against different noise models')
        endif
    end subroutine check_sigma_state

    !>  The delivered weight table into one file per state (flex_weights_state_NNN.bin: every physical
    !!  row, rows outside the selection zero, that state's scalars alongside), each registered in the
    !!  out segment of the run's own project copy as imgkind flex_weights, state NNN, beside vol_flex
    !!  state NNN. The files validate their rows against the field's `state` as the run saw it, so
    !!  this runs BEFORE write_discrete_state_project overwrites those labels.
    subroutine write_flex_weights_store( params, build, pinds, weights, labels, targets, bandwidths, l_merged )
        type(parameters), intent(in)    :: params
        type(builder),    intent(inout) :: build
        integer,          intent(in)    :: pinds(:), labels(:)
        real,             intent(in)    :: weights(:,:), targets(:,:), bandwidths(:)
        logical,          intent(in)    :: l_merged
        character(len=STDLEN) :: message
        integer :: status, s, nstates
        nstates = size(weights,2)
        call flex_weights_deliver(build%spproj, build%spproj_field, params%box, params%smpd, &
            &params%box_crop, params%smpd_crop, pinds, weights, labels, targets, bandwidths, &
            &merge(FLEX_WEIGHTS_PROV_MERGED, FLEX_WEIGHTS_PROV_FLEX_PCA, l_merged), status, message)
        if( status /= 0 ) THROW_HARD('flex_pca could not deliver the state weights: '//trim(message))
        do s = 1, nstates
            call build%spproj%add_flex_weights2os_out(flex_weights_state_fname(s), s, params%box, params%smpd)
        end do
        ! entries of a previous delivery with more states would point at files the delivery removed
        do s = nstates+1, nstates+FLEX_WEIGHTS_STALE_SCAN
            if( build%spproj%isthere_in_osout('flex_weights', s) ) call build%spproj%remove_entry_from_osout('flex_weights', s)
        end do
        call build%spproj%write_segment_inside('out', params%projfile)
        write(logfhandle,'(A,I0,A,I0,A)') '>>> FLEX_PCA STATE WEIGHTS WRITTEN: flex_weights_state_001..'// &
            &int2str_pad(nstates,3)//'.bin (', size(weights,1), ' particles x ', nstates, &
            &' states; registered in the out segment as flex_weights per state)'
        call flush(logfhandle)
    end subroutine write_flex_weights_store

    !>  Write the hard state assignment INTO the run's own project: ptcl3D/state carries each embedded
    !!  particle's label, 0 elsewhere. Judge the clusters independently of the kernel-weighted backend with
    !!      simple_exec prg=reconstruct3D projfile=<projfile> nstates=<nstates>
    !!  mkdir=yes already gave the master a private copy of the project, so this rewrites that copy and
    !!  never the project the user pointed at. MUTATES the live field, so it must run after every stage
    !!  that reads the input particle selection, and on the master only.
    subroutine write_discrete_state_project( spproj, pinds, labels, nstates, projfile )
        type(sp_project), intent(inout) :: spproj
        integer,          intent(in)    :: pinds(:), labels(:), nstates
        type(string),     intent(in)    :: projfile
        logical, allocatable :: assigned(:)
        integer :: i, iptcl, state, nptcls, nexcluded
        if( size(pinds) < 1 .or. size(labels) /= size(pinds) .or. nstates < 2 ) &
            &THROW_HARD('invalid flex_pca discrete-state assignment')
        if( len_trim(projfile%to_char()) == 0 ) THROW_HARD('flex_pca discrete-state project file is empty')
        nptcls = spproj%os_ptcl3D%get_noris()
        ! validate BEFORE mutating: this overwrites the live project field rather than a private copy,
        ! so a mid-loop abort would leave the input selection half-replaced
        allocate(assigned(nptcls), source=.false.)
        nexcluded = 0
        do i = 1,size(pinds)
            iptcl = pinds(i)
            state = labels(i)
            if( iptcl < 1 .or. iptcl > nptcls ) THROW_HARD('flex_pca discrete-state particle index outside project')
            if( assigned(iptcl) ) THROW_HARD('duplicate particle in flex_pca discrete-state assignment')
            if( state > nstates ) THROW_HARD('flex_pca discrete-state label outside state range')
            if( state < 1 ) nexcluded = nexcluded + 1
            assigned(iptcl) = .true.
        end do
        do iptcl = 1,nptcls
            call spproj%os_ptcl3D%set_state(iptcl,0)
        end do
        do i = 1,size(pinds)
            if( labels(i) >= 1 ) call spproj%os_ptcl3D%set_state(pinds(i),labels(i))
        end do
        ! ptcl3D only, as refine3D does for nstates>1: a state INDEX carries no ptcl2D meaning, and only
        ! selection (0/1) is mirrored across the two segments
        call spproj%write_segment_inside('ptcl3D', projfile)
        write(logfhandle,'(A,A)') '>>> FLEX_PCA HARD STATES WRITTEN TO: ',projfile%to_char()
        do state = 1,nstates
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX_PCA DISCRETE-STATE state=',state, &
                &' population=',count(labels==state)
        end do
        if( nexcluded > 0 ) write(logfhandle,'(A,I0)') &
            &'>>> FLEX_PCA DISCRETE-STATE unassigned particles left at state=0: ',nexcluded
        write(logfhandle,'(A,A,A,I0)') '>>> RECONSTRUCT WITH: simple_exec prg=reconstruct3D projfile=', &
            &projfile%to_char(),' nstates=',nstates
        call flush(logfhandle)
        deallocate(assigned)
    end subroutine write_discrete_state_project

    !> The project's state-1 FSC as a per-state low-pass filter sized to the delivered map
    !! (filtsz = fdim(box_rec)-1); has_filter false wherever the project cannot provide one.
    subroutine prepare_project_fsc_lowpass_filters( params, filtsz, nstates, lowpass_filters, has_filter, source_state )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: filtsz, nstates
        real, allocatable, intent(out) :: lowpass_filters(:,:)
        logical, allocatable, intent(out) :: has_filter(:)
        integer, allocatable, intent(out) :: source_state(:)
        type(sp_project) :: spproj
        type(string) :: fsc_fname, imgkind_here, proj_for_fsc
        real, allocatable :: fsc(:)
        integer :: state, fsc_box, i, state1_fsc_count
        logical :: out_loaded
        allocate(lowpass_filters(filtsz,nstates),has_filter(nstates),source_state(nstates))
        lowpass_filters=0.
        has_filter=.false.
        source_state=0
        if( filtsz<1 ) return
        proj_for_fsc=params%projfile
        if( .not.file_exists(proj_for_fsc) )then
            write(logfhandle,'(A)') '>>> FLEX PRE-IMAGE project-FSC low-pass unavailable (projfile not found); states rely on their own eo-FSC'
            has_filter = .false.   ! the state's own eo-FSC decides
            call proj_for_fsc%kill
            return
        endif
        call spproj%read_segment('out',proj_for_fsc)
        out_loaded=spproj%os_out%get_noris()>0
        if( .not.out_loaded )then
            write(logfhandle,'(A)') '>>> FLEX PRE-IMAGE project-FSC low-pass unavailable (empty out segment); states rely on their own eo-FSC'
            has_filter = .false.   ! the state's own eo-FSC decides
            call proj_for_fsc%kill
            call spproj%kill
            return
        endif
        state1_fsc_count=0
        do i=1,spproj%os_out%get_noris()
            if( .not.spproj%os_out%isthere(i,'imgkind') ) cycle
            call spproj%os_out%getter(i,'imgkind',imgkind_here)
            if( imgkind_here%to_char()/='fsc' ) cycle
            if( spproj%os_out%get_state(i)==1 ) state1_fsc_count=state1_fsc_count+1
        end do
        if( state1_fsc_count/=1 )then
            write(logfhandle,'(A,I0)') '>>> FLEX PRE-IMAGE project-FSC low-pass unavailable (state=1 FSC count=', &
                &state1_fsc_count
            write(logfhandle,'(A)') '>>>   ); states rely on their own eo-FSC'
            has_filter = .false.   ! the state's own eo-FSC decides
            call imgkind_here%kill
            call proj_for_fsc%kill
            call spproj%kill
            return
        endif
        call spproj%get_fsc(1,fsc_fname,fsc_box)
        if( .not.file_exists(fsc_fname) )then
            write(logfhandle,'(A)') '>>> FLEX PRE-IMAGE project-FSC low-pass unavailable (state=1 FSC file missing); states rely on their own eo-FSC'
            has_filter = .false.   ! the state's own eo-FSC decides
            call imgkind_here%kill
            call proj_for_fsc%kill
            call spproj%kill
            return
        endif
        fsc=file2rarr(fsc_fname)
        if( size(fsc)/=filtsz )then
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX PRE-IMAGE project-FSC low-pass unavailable (state=1 FSC size mismatch; fsc_nyq=', &
                &size(fsc),' model_nyq=',filtsz
            write(logfhandle,'(A)') '>>>   ); states rely on their own eo-FSC'
            has_filter = .false.   ! the state's own eo-FSC decides
            deallocate(fsc)
            call fsc_fname%kill
            call imgkind_here%kill
            call proj_for_fsc%kill
            call spproj%kill
            return
        endif
        do state=1,nstates
            call fsc2optlp_sub(filtsz,fsc,lowpass_filters(:,state),merged=.false.)
            has_filter(state)=any(lowpass_filters(:,state)>0.)
            source_state(state)=1
        end do
        deallocate(fsc)
        call fsc_fname%kill
        call imgkind_here%kill
        call proj_for_fsc%kill
        call spproj%kill
    end subroutine prepare_project_fsc_lowpass_filters

    !> Register one delivered state map in the project's out segment (imgkind vol_flex).
    subroutine publish_state_volume( build, vol_fname, smpd, state, box )
        class(builder), intent(inout) :: build
        class(string),  intent(in)    :: vol_fname
        real,           intent(in)    :: smpd
        integer,        intent(in)    :: state, box
        call build%spproj%add_vol2os_out(vol_fname, smpd, state, 'vol_flex', box=box)
    end subroutine publish_state_volume

    !> Write the out segment back to the project file.
    subroutine write_out_segment( build, projfile )
        class(builder), intent(inout) :: build
        class(string),  intent(in)    :: projfile
        call build%spproj%write_segment_inside('out', projfile)
    end subroutine write_out_segment

end module simple_flex_pca_project_gateway
