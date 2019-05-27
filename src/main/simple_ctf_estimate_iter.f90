! ctf_estimate iterator
module simple_ctf_estimate_iter
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_image,      only: image
implicit none

public :: ctf_estimate_iter
private
#include "simple_local_flags.inc"

logical,       parameter :: l_patch = .false.

type :: ctf_estimate_iter
    type(image) :: micrograph, pspec, thumbnail, img_jpg, pspec4viz
    type(image), allocatable :: windows(:,:)
  contains
      procedure          :: iterate
      procedure, private :: gen_windows
      procedure, private :: mic2spec_patch
      procedure          :: cleanup
end type ctf_estimate_iter

contains

    subroutine iterate( self, ctfvars, moviename_forctf, orientation, dir_out, l_gen_thumb )
        use simple_ori, only: ori
        use simple_ctf_estimate_patch
        class(ctf_estimate_iter), intent(inout) :: self
        class(ctfparams),         intent(inout) :: ctfvars
        character(len=*),         intent(in)    :: moviename_forctf
        class(ori),               intent(inout) :: orientation
        character(len=*),         intent(in)    :: dir_out
        logical,                  intent(in)    :: l_gen_thumb
        type(ctfparams)               :: ctfvars_glob
        type(ctf_estimate_patch)      :: ctffit
        character(len=:), allocatable :: fname_diag
        character(len=LONGSTRLEN)     :: moviename_thumb, rel_moviename_thumb, rel_ctfjpg
        real                          :: cc, ctfscore, cc90, scale
        integer                       :: nframes, ldim(3), ldim_thumb(3), i,j, start(2), incr(2), center(2)
        if( .not. file_exists(moviename_forctf) )&
        & write(logfhandle,*) 'inputted micrograph does not exist: ', trim(adjustl(moviename_forctf))
        call find_ldim_nptcls(trim(adjustl(moviename_forctf)), ldim, nframes)
        if( nframes /= 1 )then
            write(logfhandle,*) 'nframes: ', nframes
            THROW_HARD('single frame input to ctf_estimate assumed; iterate')
        endif
        ldim(3) = 1
        call self%micrograph%new(ldim, ctfvars%smpd)
        call self%micrograph%read(trim(adjustl(moviename_forctf)), 1)
        ! filter out frequencies lower than the box can express to avoid aliasing
        call self%micrograph%bp(real(params_glob%pspecsz) * ctfvars%smpd, 0.)
        ! extract powerspectra
         call self%micrograph%mic2spec(params_glob%pspecsz, 'sqrt', params_glob%hp, self%pspec)
        if( l_gen_thumb )then
            ! generate thumbnail
            scale         = real(params_glob%pspecsz)/real(ldim(1))
            ldim_thumb(1) = round2even(real(ldim(1))*scale)
            ldim_thumb(2) = round2even(real(ldim(2))*scale)
            ldim_thumb(3) = 1
            call self%thumbnail%new(ldim_thumb, ctfvars%smpd)
            call self%micrograph%fft()
            call self%micrograph%clip(self%thumbnail)
            call self%thumbnail%ifft()
            ! jpeg output
            call self%pspec4viz%copy(self%pspec)
            call self%pspec4viz%scale_pspec4viz
            call self%pspec4viz%collage(self%thumbnail, self%img_jpg)
            moviename_thumb = trim(get_fbody(basename(trim(moviename_forctf)), params_glob%ext, separator=.false.))
            moviename_thumb = swap_suffix(moviename_thumb, THUMBNAIL_SUFFIX, INTGMOV_SUFFIX)
            moviename_thumb = trim(dir_out)//trim(adjustl(moviename_thumb))//trim(JPG_EXT)
            call self%img_jpg%write_jpg(moviename_thumb, norm=.true., quality=90)
            call make_relativepath(CWD_GLOB,moviename_thumb, rel_moviename_thumb)
            call orientation%set('thumb', trim(rel_moviename_thumb))
        endif
        ! deal with output
        fname_diag = trim(get_fbody(basename(trim(moviename_forctf)), params_glob%ext, separator=.false.))
        if( str_has_substr(fname_diag, FORCTF_SUFFIX) )then
            fname_diag = swap_suffix(fname_diag, '_ctf_estimate_diag', FORCTF_SUFFIX)
        else if( str_has_substr(fname_diag, INTGMOV_SUFFIX) )then
            fname_diag = swap_suffix(fname_diag, '_ctf_estimate_diag', INTGMOV_SUFFIX)
        endif
        fname_diag = filepath(trim(dir_out),trim(fname_diag)//trim(JPG_EXT))
        ! fitting
        call ctffit%new(self%pspec, ctfvars, [params_glob%dfmin,params_glob%dfmax],&
            &[params_glob%hp,params_glob%lp], params_glob%astigtol, .false.)
        call ctffit%fit( ctfvars, fname_diag )
        cc       = ctffit%get_ccfit()
        cc90     = ctffit%get_cc90()
        ctfscore = ctffit%get_ctfscore()
        call ctffit%kill
        call make_relativepath(CWD_GLOB,fname_diag, rel_ctfjpg)
        if( l_patch )then
            call self%micrograph%ifft
            call self%gen_windows(params_glob%pspecsz,'sqrt')
            print *,ctfvars%dfx,ctfvars%dfy,ctfvars%angast,ctffit%get_ccfit()
            ctfvars_glob%dfx    = ctfvars%dfx
            ctfvars_glob%dfy    = ctfvars%dfy
            ctfvars_glob%angast = ctfvars%angast
            incr   = floor(real(ldim(1:2))/3.)
            start  = floor(real(incr)/2.)
            center = start
            do i=1,3
                center(1) = start(1) +(i-1)*incr(1)
                do j=1,3
                    center(2) = start(2) +(j-1)*incr(2)
                    ctfvars%dfx    = ctfvars_glob%dfx
                    ctfvars%dfy    = ctfvars_glob%dfy
                    ctfvars%angast = ctfvars_glob%angast
                    call self%mic2spec_patch(center,params_glob%hp)
                    fname_diag = 'pspec_'//int2str_pad(i,2)//'_'//int2str_pad(j,2)//'.jpg'
                    call ctffit%new(self%pspec, ctfvars, [params_glob%dfmin,params_glob%dfmax],&
                        &[params_glob%hp,params_glob%lp], params_glob%astigtol, .true.)
                    call ctffit%fit(ctfvars, fname_diag)
                    print *,ctfvars%dfx,ctfvars%dfy,ctfvars%angast,ctffit%get_ccfit()
                enddo
            enddo
        endif
        ! reporting
        call orientation%set('dfx',      ctfvars%dfx)
        call orientation%set('dfy',      ctfvars%dfy)
        call orientation%set('angast',   ctfvars%angast)
        call orientation%set('phshift',  ctfvars%phshift)
        call orientation%set('ctf_estimatecc', cc)
        call orientation%set('ctfscore', ctfscore)
        call orientation%set('cc90',     cc90)
        call orientation%set('ctfjpg',   rel_ctfjpg)
    end subroutine iterate

    subroutine gen_windows(self, box, speckind)
        class(ctf_estimate_iter), intent(inout) :: self
        integer,                  intent(in)    :: box
        character(len=*),         intent(in)    :: speckind
        type(image) :: tmp
        real        :: smpd
        integer     :: xind,yind, ldim(3), i,j, nx,ny, ldim2(3)
        logical     :: outside
        smpd = self%micrograph%get_smpd()
        if( allocated(self%windows) )then
                ldim2 = self%windows(1,1)%get_ldim()
                if( .not.is_equal(smpd,self%windows(1,1)%get_smpd()) .or. ldim2(1).ne.box )then
                do xind=1,size(self%windows,dim=1)
                do yind=1,size(self%windows,dim=2)
                    call self%windows(xind,yind)%kill
                enddo
                enddo
                deallocate(self%windows)
            endif
        endif
        smpd = self%micrograph%get_smpd()
        ldim = self%micrograph%get_ldim()
        call tmp%new([box,box,1], smpd)
        nx = 0
        do xind=0,ldim(1)-box,box/2
            nx = nx+1
        end do
        ny = 0
        do yind=0,ldim(2)-box,box/2
            ny = ny+1
        end do
        allocate(self%windows(nx,ny))
        i = 0
        do xind=0,ldim(1)-box,box/2
            i = i+1
            j = 0
            do yind=0,ldim(2)-box,box/2
                j = j+1
                call self%windows(i,j)%new([box,box,1], smpd)
                call self%micrograph%window_slim([xind,yind],box,tmp,outside)
                call tmp%norm
                call tmp%zero_edgeavg
                call tmp%fft
                call tmp%ft2img(speckind, self%windows(i,j))
                call tmp%zero_and_unflag_ft
                call self%windows(i,j)%dampen_pspec_central_cross
                call self%windows(i,j)%write_jpg(int2str_pad(i,2)//'_'//int2str_pad(j,2)//'.jpg')
            end do
        end do
        call tmp%kill
    end subroutine gen_windows

    !> mic2spec calculates the average powerspectrum over a micrograph
    !!          the resulting spectrum has dampened central cross and subtracted background
    subroutine mic2spec_patch(self, center, lp_backgr_subtr )
        class(ctf_estimate_iter), intent(inout) :: self
        integer,          intent(in)    :: center(2)
        real,             intent(in)    :: lp_backgr_subtr
        real, allocatable :: dists(:)
        real        :: smpd, dist, dist_thresh
        integer     :: ldim_win(3),xind, yind, cnt, ldim(3), center_win(2), box,i,j
        smpd     = self%micrograph%get_smpd()
        ldim     = self%micrograph%get_ldim()
        ldim_win = self%windows(1,1)%get_ldim()
        box = ldim_win(1)
        call self%pspec%new([box,box,1], smpd)
        cnt = 0
        do xind=0,ldim(1)-box,box/2
            do yind=0,ldim(2)-box,box/2
                cnt = cnt+1
            end do
        end do
        allocate(dists(cnt))
        cnt = 0
        do xind=0,ldim(1)-box,box/2
            do yind=0,ldim(2)-box,box/2
                cnt = cnt+1
                center_win = [xind, yind] + box/2
                dists(cnt) = sqrt(real(sum((center-center_win)**2)))
            end do
        end do
        call hpsort(dists)
        dist_thresh = dists(nint(real(cnt)*0.25))
        cnt = 0
        i = 0
        do xind=0,ldim(1)-box,box/2
            i = i+1
            j = 0
            do yind=0,ldim(2)-box,box/2
                j = j+1
                center_win = [xind, yind] + box/2
                dist = sqrt(real(sum((center-center_win)**2)))
                if( dist > dist_thresh ) cycle
                cnt = cnt+1
                call self%pspec%add(self%windows(i,j))
            end do
        end do
        call self%pspec%div(real(cnt))
        call self%pspec%dampen_pspec_central_cross
        call self%pspec%subtr_backgr(lp_backgr_subtr)
    end subroutine mic2spec_patch

    subroutine cleanup( self )
        class(ctf_estimate_iter), intent(inout) :: self
        integer :: i,j
        call self%micrograph%kill
        call self%pspec%kill
        call self%thumbnail%kill
        call self%img_jpg%kill
        call self%pspec4viz%kill
        if( allocated(self%windows) )then
            do i=1,size(self%windows,dim=1)
            do j=1,size(self%windows,dim=2)
                call self%windows(i,j)%kill
            enddo
            enddo
            deallocate(self%windows)
        endif
    end subroutine cleanup

end module simple_ctf_estimate_iter
