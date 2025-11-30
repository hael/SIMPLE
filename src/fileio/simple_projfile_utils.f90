module simple_projfile_utils
include 'simple_lib.f08'
use simple_image,      only: image
use simple_parameters, only: params_glob
use simple_sp_project, only: sp_project
use simple_class_frcs
implicit none
#include "simple_local_flags.inc"

contains

    subroutine merge_chunk_projfiles( chunk_fnames, folder, merged_proj, projname_out, write_proj )
        class(string),           intent(in)    :: chunk_fnames(:) ! List of project files
        class(string),           intent(in)    :: folder          ! output folder
        class(sp_project),       intent(inout) :: merged_proj     ! output project, assumed to have compuational env info
        class(string), optional, intent(in)    :: projname_out    ! name for output project file
        logical,       optional, intent(in)    :: write_proj      ! write project file
        type(sp_project), allocatable :: chunks(:)
        real,             allocatable :: states(:)
        integer,          allocatable :: clsmap(:)
        type(class_frcs) :: frcs, frcs_chunk
        type(image)      :: img
        type(string)     :: projname, stkname, evenname, oddname, frc_fname, projfile_out, dir, cavgs
        real             :: smpd
        integer          :: ldim(3), i, ic, icls, ncls, nchunks, nallmics, nallstks, nallptcls, ncls_tot, box4frc
        integer          :: fromp, fromp_glob, top, top_glob, j, iptcl_glob, nstks, nmics, nptcls, istk
        logical          :: l_write_proj
        l_write_proj = .true.
        if(present(write_proj)) l_write_proj = write_proj
        nchunks = size(chunk_fnames)
        allocate(chunks(nchunks))
        dir = folder%to_char()//'/'
        if( present(projname_out) )then
            projfile_out = dir%to_char()//projname_out%to_char()//trim(METADATA_EXT)
        else
            projfile_out = dir%to_char()//'set'//METADATA_EXT
        endif
        call merged_proj%os_mic%kill
        call merged_proj%os_stk%kill
        call merged_proj%os_ptcl2D%kill
        call merged_proj%os_ptcl3D%kill
        call merged_proj%os_cls2D%kill
        call merged_proj%os_cls3D%kill
        call merged_proj%os_out%kill
        cavgs     = dir%to_char()//'cavgs.mrc'
        nallptcls = 0
        nallstks  = 0
        nallmics  = 0
        icls      = 0
        do ic = 1,nchunks
            projname = chunk_fnames(ic)
            call chunks(ic)%read_data_info(projname, nmics, nstks, nptcls)
            nallmics  = nallmics  + nmics
            nallstks  = nallstks  + nstks
            nallptcls = nallptcls + nptcls
            call chunks(ic)%read_segment('out',  projname)
            call chunks(ic)%read_segment('cls2D',  projname)
            call chunks(ic)%get_cavgs_stk(stkname, ncls, smpd, imgkind='cavg')
            call find_ldim_nptcls(stkname, ldim, ncls)
            ldim(3) = 1
            call img%new(ldim, smpd)
            evenname = add2fbody(stkname, params_glob%ext, '_even')
            oddname  = add2fbody(stkname, params_glob%ext, '_odd')
            do i = 1,ncls
                if( chunks(ic)%os_cls2D%get_state(i) == 0 ) cycle
                icls = icls+1
                call img%read(stkname,i)
                call img%write(cavgs,icls)
                call img%read(evenname,i)
                call img%write(dir//'cavgs_even.mrc',icls)
                call img%read(oddname,i)
                call img%write(dir//'cavgs_odd.mrc',icls)
            enddo
        enddo
        call img%kill
        ncls_tot = icls
        ! micrographs
        if( nallmics > 0 )then
            call merged_proj%os_mic%new(nallmics,.false.)
            j = 0
            do ic = 1,nchunks
                projname = chunk_fnames(ic)
                call chunks(ic)%read_segment('mic', projname)
                do i = 1,chunks(ic)%os_mic%get_noris()
                    j = j+1
                    call merged_proj%os_mic%transfer_ori(j, chunks(ic)%os_mic, i)
                enddo
                call chunks(ic)%os_mic%kill
            enddo
        endif
        ! particles, stacks and classes frcs & metadata
        call merged_proj%os_cls2D%new(ncls_tot,  .false.)
        call merged_proj%os_ptcl2D%new(nallptcls,.true.)
        call merged_proj%os_stk%new(nallstks,    .false.)
        icls       = 0
        istk       = 0
        iptcl_glob = 0
        fromp_glob = 1
        do ic = 1,nchunks
            projname = chunk_fnames(ic)
            call chunks(ic)%read_segment('stk', projname)
            call chunks(ic)%read_segment('ptcl2D',projname)
            ! classes frcs & info
            call chunks(ic)%get_frcs(frc_fname, 'frc2D')
            ncls = chunks(ic)%os_cls2D%get_noris()
            call frcs_chunk%read(frc_fname)
            if( ic == 1 )then
                box4frc = frcs_chunk%get_box()
                call frcs%new(ncls_tot, box4frc, smpd)
            endif
            allocate(clsmap(ncls),source=0)
            do i = 1,ncls
                if( chunks(ic)%os_cls2D%get_state(i) == 0 ) cycle
                icls      = icls+1
                clsmap(i) = icls
                call merged_proj%os_cls2D%transfer_ori(icls,   chunks(ic)%os_cls2D, i)
                call merged_proj%os_cls2D%set_class(icls, icls)
                call merged_proj%os_cls2D%set(icls,'origclass',i)
                call merged_proj%os_cls2D%set(icls,'chunk',    projname)
                call frcs%set_frc(icls, frcs_chunk%get_frc(i,  box4frc))
            enddo
            ! particles and stacks
            nstks  = chunks(ic)%os_stk%get_noris()
            do i = 1,nstks
                istk  = istk + 1
                fromp = chunks(ic)%os_stk%get_fromp(i)
                top   = chunks(ic)%os_stk%get_top(i)
                do j = fromp,top
                    iptcl_glob = iptcl_glob + 1
                    if( chunks(ic)%os_ptcl2D%get_state(j) > 0 )then
                        call chunks(ic)%os_ptcl2D%set_class(j, clsmap(chunks(ic)%os_ptcl2D%get_class(j)))
                    endif
                    call chunks(ic)%os_ptcl2D%set_stkind(j, istk)
                    call merged_proj%os_ptcl2D%transfer_ori(iptcl_glob, chunks(ic)%os_ptcl2D, j)
                enddo
                top_glob = fromp_glob + top - fromp
                call chunks(ic)%os_stk%set(i, 'fromp', fromp_glob)
                call chunks(ic)%os_stk%set(i, 'top',   top_glob)
                fromp_glob = top_glob+1
                call merged_proj%os_stk%transfer_ori(istk, chunks(ic)%os_stk, i)
            enddo
            deallocate(clsmap)
            ! making sure the compenv is informed
            if( (ic == 1) .and. (merged_proj%compenv%get_noris() == 0) )then
                call chunks(1)%read_non_data_segments(chunk_fnames(1))
                merged_proj%compenv = chunks(1)%compenv
                call merged_proj%update_projinfo(projfile_out)
            endif
            ! cleanup
            call chunks(ic)%kill
        enddo
        deallocate(chunks)
        ! add classes, frcs
        call frcs%write(dir//trim(FRCS_FILE))
        call merged_proj%add_frcs2os_out(dir//trim(FRCS_FILE), 'frc2D')
        call merged_proj%add_cavgs2os_out(cavgs, smpd, imgkind='cavg')
        states = merged_proj%os_cls2D%get_all('state')
        call merged_proj%os_cls3D%new(ncls_tot, .false.)
        call merged_proj%os_cls3D%set_all('state', states)
        merged_proj%os_ptcl3D = merged_proj%os_ptcl2D
        call merged_proj%os_ptcl3D%delete_2Dclustering
        ! write
        if(l_write_proj) call merged_proj%write(projfile_out)
        ! cleanup
        call frcs%kill
        call frcs_chunk%kill
    end subroutine merge_chunk_projfiles

end module simple_projfile_utils