submodule(simple_sp_project) simple_sp_project_cls
implicit none
#include "simple_local_flags.inc"
contains 

    module subroutine shape_ranked_cavgs2jpg( self, cavg_inds, jpgname, xtiles, ytiles, mskdiam_px )
        class(sp_project),    intent(inout) :: self
        integer, allocatable, intent(inout) :: cavg_inds(:)
        class(string),        intent(in)    :: jpgname
        integer,              intent(out)   :: xtiles, ytiles
        integer, optional,    intent(in)    :: mskdiam_px
        integer,          allocatable :: shape_ranks(:)
        type(string)   :: cavgsstk
        type(image)    :: img, jpegimg
        type(stack_io) :: stkio_r
        integer        :: icls, ncls, ncls_stk, ldim_read(3), ncls_sel, ix, iy, ntiles
        real           :: smpd
        ncls = self%os_cls2D%get_noris()
        if( ncls == 0 ) return
        if( .not. self%os_cls2D%isthere('shape_rank') ) THROW_HARD('No shape ranking available!')
        if( allocated(cavg_inds) ) deallocate(cavg_inds)
        allocate(shape_ranks(ncls), cavg_inds(ncls))
        do icls = 1, ncls
            shape_ranks(icls) = self%os_cls2D%get_int(icls,'shape_rank')
            cavg_inds(icls)   = icls
        end do
        cavg_inds   = pack(cavg_inds,   mask=shape_ranks > 0)
        shape_ranks = pack(shape_ranks, mask=shape_ranks > 0)
        call hpsort(shape_ranks, cavg_inds)
        call self%get_cavgs_stk(cavgsstk, ncls_stk, smpd, imgkind='cavg')
        if(.not. file_exists(cavgsstk)) THROW_HARD('cavgs stk does not exist')
        if( ncls /= ncls_stk ) THROW_HARD('Inconsistent # cavgs in spproj and stack file')
        call stkio_r%open(cavgsstk, smpd, 'read', bufsz=ncls)
        ldim_read    = stkio_r%get_ldim()
        ldim_read(3) = 1
        call stkio_r%read_whole
        ncls_sel = size(cavg_inds)
        xtiles   = floor(sqrt(real(ncls_sel)))
        ytiles   = ceiling(real(ncls_sel) / real(xtiles))
        call jpegimg%new([xtiles * JPEG_DIM, ytiles * JPEG_DIM, 1], smpd)
        ix     = 1
        iy     = 1
        ntiles = 0
        do icls = 1,ncls_sel
            call img%new(ldim_read, smpd)
            call stkio_r%get_image(cavg_inds(icls), img)
            if(present(mskdiam_px)) call img%mask(mskdiam_px / 2.0, 'softavg')
            call img%fft
            if(ldim_read(1) > JPEG_DIM) then
                call img%clip_inplace([JPEG_DIM,JPEG_DIM,1])
            else
                call img%pad_inplace([JPEG_DIM,JPEG_DIM,1], backgr=0., antialiasing=.false.)
            end if
            call img%ifft
            call jpegimg%tile(img, ix, iy)
            ntiles = ntiles + 1
            ix = ix + 1
            if(ix > xtiles) then
                ix = 1
                iy = iy + 1
            end if
        enddo
        call jpegimg%write_jpg(jpgname)
        call stkio_r%close
        call img%kill
        call jpegimg%kill
    end subroutine shape_ranked_cavgs2jpg

    module subroutine cavgs2jpg( self, cavg_inds, jpgname, xtiles, ytiles )
        class(sp_project),    intent(inout) :: self
        integer, allocatable, intent(inout) :: cavg_inds(:)
        class(string),        intent(in)    :: jpgname
        integer,              intent(out)   :: xtiles, ytiles
        logical, allocatable       :: cls_mask(:)
        type(string)   :: cavgsstk
        type(image)    :: img, jpegimg
        type(stack_io) :: stkio_r
        integer        :: icls, isel, ncls, ncls_stk, ldim_read(3), ncls_sel, ix, iy, ntiles
        real           :: smpd
        ncls = self%os_cls2D%get_noris()
        if( ncls == 0 ) return
        if( allocated(cavg_inds) ) deallocate(cavg_inds)
        allocate(cls_mask(ncls), cavg_inds(ncls))
        cls_mask = .true.
        do icls = 1, ncls
            if( self%os_cls2D%get_int(icls,'pop') == 0 ) cls_mask(icls) = .false.
            if( self%os_cls2D%get_state(icls)     == 0 ) cls_mask(icls) = .false.
            cavg_inds(icls) = icls
        enddo
        ncls_sel = count(cls_mask)
        if( ncls_sel == 0 ) return
        call self%get_cavgs_stk(cavgsstk, ncls_stk, smpd, imgkind='cavg')
        if(.not. file_exists(cavgsstk)) THROW_HARD('cavgs stk does not exist')
        if( ncls /= ncls_stk ) THROW_HARD('Inconsistent # cavgs in spproj and stack file')
        call stkio_r%open(cavgsstk, smpd, 'read', bufsz=ncls)
        ldim_read    = stkio_r%get_ldim()
        ldim_read(3) = 1
        call stkio_r%read_whole
        xtiles = floor(sqrt(real(ncls_sel)))
        ytiles = ceiling(real(ncls_sel) / real(xtiles))
        call jpegimg%new([xtiles * JPEG_DIM, ytiles * JPEG_DIM, 1], smpd)
        isel   = 0
        ix     = 1
        iy     = 1
        ntiles = 0
        do icls = 1,ncls
            if( cls_mask(icls) ) then
                isel = isel + 1
                call img%new(ldim_read, smpd)
                call stkio_r%get_image(icls, img)
                call img%fft
                if(ldim_read(1) > JPEG_DIM) then
                    call img%clip_inplace([JPEG_DIM,JPEG_DIM,1])
                else
                    call img%pad_inplace([JPEG_DIM,JPEG_DIM,1], backgr=0., antialiasing=.false.)
                end if
                call img%ifft
                call jpegimg%tile(img, ix, iy)
                ntiles = ntiles + 1
                ix = ix + 1
                if(ix > xtiles) then
                    ix = 1
                    iy = iy + 1
                end if
            endif
        enddo
        call jpegimg%write_jpg(jpgname)
        call stkio_r%close
        call img%kill
        call jpegimg%kill
    end subroutine cavgs2jpg

    module subroutine cavgs2mrc( self )
        class(sp_project), intent(inout) :: self
        logical, allocatable       :: cls_mask(:)
        type(string)   :: cavgsstk
        type(image)    :: img
        type(stack_io) :: stkio_r, stkio_w
        integer        :: icls, ncls, ncls_stk, isel, ldim_read(3), ncls_sel
        real           :: smpd
        ncls = self%os_cls2D%get_noris()
        if( ncls == 0 ) return
        allocate(cls_mask(ncls))
        cls_mask = .true.
        do icls = 1, ncls
            if( self%os_cls2D%get_state(icls)     == 0 ) cls_mask(icls) = .false.
        enddo
        ncls_sel = count(cls_mask)
        if( ncls_sel == 0 ) return
        call self%get_cavgs_stk(cavgsstk, ncls_stk, smpd, imgkind='cavg')
        if(.not. file_exists(cavgsstk)) THROW_HARD('cavgs stk does not exist')
        if( ncls /= ncls_stk ) THROW_HARD('Inconsistent # cavgs in spproj and stack file')
        call stkio_r%open(cavgsstk, smpd, 'read', bufsz=ncls)
        ldim_read    = stkio_r%get_ldim()
        ldim_read(3) = 1
        call stkio_r%read_whole
        call stkio_w%open(string("cls2D_selected.mrcs"), smpd, 'write', box=ldim_read(1), bufsz=ncls_sel)
        isel = 0
        do icls = 1,ncls
            if( cls_mask(icls) ) then
                isel = isel + 1
                call img%new(ldim_read, smpd)
                call stkio_r%get_image(icls, img)
                call stkio_w%write(isel, img)
            endif
        enddo
        call stkio_r%close
        call stkio_w%close
        call img%kill
        if(allocated(cls_mask)) deallocate(cls_mask)
    end subroutine cavgs2mrc

    ! Getters

    module function get_selected_clsinds( self ) result( clsinds )
        class(sp_project),    intent(inout) :: self
        integer, allocatable :: tmpinds(:), states_cavgs(:), clsinds(:)
        integer :: ncls, icls
        ncls         = self%os_cls2D%get_noris()
        if( ncls == 0 ) THROW_HARD('no entries in cls2D field of project')
        tmpinds      = (/(icls,icls=1,ncls)/)
        states_cavgs = self%os_cls2D%get_all_asint('state')
        if( any(states_cavgs == 0 ) )then
            clsinds = pack(tmpinds, mask=states_cavgs > 0)
        else
            THROW_WARN('no deletions in cls2D field, returning contiguous array of class indices')
            clsinds = tmpinds
        endif
        deallocate(tmpinds, states_cavgs)
    end function get_selected_clsinds

    module subroutine get_cavgs_stk( self, stkname, ncls, smpd, imgkind, fail, out_ind, box )
        class(sp_project),          intent(inout) :: self
        class(string),              intent(inout) :: stkname
        integer,                    intent(out)   :: ncls
        real,                       intent(out)   :: smpd
        character(len=*), optional, intent(in)    :: imgkind
        logical,          optional, intent(in)    :: fail
        integer,          optional, intent(inout) :: out_ind
        real,             optional, intent(inout) :: box
        type(string) :: ikind, iimgkind, stkpath
        integer :: n_os_out, ind, i, cnt
        logical :: fail_here
        if( present(imgkind) )then
            iimgkind = trim(imgkind)
        else
            iimgkind ='cavg'
        endif
        fail_here = .true.
        if( present(fail) )fail_here = fail
        ! check if field is empty
        n_os_out = self%os_out%get_noris()
        if( n_os_out == 0 )then
            if( fail_here )then
                THROW_HARD('trying to fetch from empty os_out field; get_cavgs_stk')
            else
                stkname = NIL
                ncls    = 0
                smpd    = 0.
                return
            endif
        endif
        ! look for cavgs
        ind = 0
        cnt = 0
        do i=1,n_os_out
            if( self%os_out%isthere(i,'imgkind') )then
                ikind = self%os_out%get_str(i,'imgkind')
                if(ikind .eq. iimgkind)then
                    ind = i
                    cnt = cnt + 1
                endif
            endif
        end do
        if( fail_here )then
            if( cnt > 1 )  THROW_HARD('multiple os_out entries with imgkind='//iimgkind%to_char()//', aborting... get_cavgs_stk')
            if( cnt == 0 ) THROW_HARD('no os_out entry with imgkind='//iimgkind%to_char()//' identified, aborting... get_cavgs_stk')
        else if( cnt > 1 .or. cnt == 0 )then
            stkname = NIL
            ncls    = 0
            smpd    = 0.
            return
        endif
        ! set return values
        stkname = self%os_out%get_str(ind,'stk')
        if( .not. file_exists(stkname) )then
            if(self%os_out%isthere(ind, 'stkpath'))then
                stkpath = self%os_out%get_str(ind,'stkpath')
                stkname = stkpath%to_char()//'/'//stkname%to_char()
            endif
        endif
        ncls    = self%os_out%get_int(ind, 'nptcls')
        smpd    = self%os_out%get(ind, 'smpd')
        if(present(out_ind)) out_ind = ind
        if(present(box)    ) box     = self%os_out%get(ind, 'box')
    end subroutine get_cavgs_stk

    module subroutine set_cavgs_thumb( self, projfile )
        class(sp_project),  intent(inout) :: self
        class(string),      intent(in)    :: projfile
        integer,            allocatable   :: clsinds(:)
        logical,            allocatable   :: clsmsk(:)
        type(string)                      :: cavgsstk
        integer                           :: ncls, n_thumbnails, out_ind, iori, ithumb, nxtiles, nytiles
        real                              :: smpd, thumbnail_scale
        call self%get_cavgs_stk(cavgsstk, ncls, smpd, out_ind=out_ind)
        call self%os_cls2D%mask_from_state( 1, clsmsk, clsinds )
        call mrc2jpeg_tiled(cavgsstk, stemname(projfile) //"/thumb2D.jpeg", scale=thumbnail_scale, ntiles=n_thumbnails, msk=clsmsk, n_xtiles=nxtiles, n_ytiles=nytiles)
        ithumb = 1
        do iori=1, self%os_cls2D%get_noris()
            call self%os_cls2D%set(iori, "thumb",    stemname(projfile) //"/thumb2D.jpeg")
            call self%os_cls2D%set(iori, "thumbn",   count(clsmsk))
            call self%os_cls2D%set(iori, "thumbdim", JPEG_DIM)
            call self%os_cls2D%set(iori, "thumbnx",  nxtiles)
            call self%os_cls2D%set(iori, "thumbny",  nytiles)
            if(clsmsk(iori)) then
                call self%os_cls2D%set(iori, "thumbidx", ithumb)
                ithumb = ithumb + 1
            end if
        end do
        call cavgsstk%kill
    end subroutine set_cavgs_thumb

end submodule simple_sp_project_cls