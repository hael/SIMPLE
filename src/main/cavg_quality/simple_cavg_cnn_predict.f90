!@descr: the body of program cavg_cnn_predict: the class-average CNN's group (junk, secondary structure, likely protein) of every class average of a project, written to a table
! Read-only: the project is read, never written; doc/how2s/how2train_cavg_cnn.md. cavg_cnn_class_probs gives the
! probabilities of a stack of class averages to other callers (particle_sieving_2).
module simple_cavg_cnn_predict
use simple_defs,           only: sp, dp, logfhandle
use simple_error,          only: simple_exception
use simple_string,         only: string
use simple_parameters,     only: parameters
use simple_sp_project,     only: sp_project
use simple_image,          only: image
use simple_imgarr_utils,   only: read_stk_into_imgarr, dealloc_imgarr
use simple_fileio,         only: add2fbody, fname2ext, fopen, fclose, fileiochk
use simple_syslib,         only: file_exists
use simple_cavg_cnn,       only: cavg_cnn
use simple_cavg_cnn_input, only: cavg_cnn_input
use simple_defs_cavg_cnn,  only: CAVG_CNN_PRED_TABLE, CAVG_CNN_GROUP_NAMES, CAVG_CNN_GROUP_EMPTY
implicit none

public :: cavg_cnn_predict_exec, cavg_cnn_class_probs
private
#include "simple_local_flags.inc"

contains

    !> project: reads the out and cls2D segments of params%projfile, writes none; the table goes to the working folder
    subroutine cavg_cnn_predict_exec( params )
        class(parameters), intent(inout) :: params
        type(sp_project)         :: spproj
        type(string)             :: stk
        real(dp),    allocatable :: pclass(:,:)
        integer,     allocatable :: group(:)
        logical,     allocatable :: empty(:)
        real                     :: smpd, smpd_net, pop, res
        integer                  :: ncls, icls, funit, ios, g
        if( .not. params%infile%is_allocated() ) THROW_HARD('infile, the weights file, is required; cavg_cnn_predict')
        if( params%infile%strlen_trim() == 0 )   THROW_HARD('infile, the weights file, is required; cavg_cnn_predict')
        call spproj%read_segment('out',   params%projfile)
        call spproj%read_segment('cls2D', params%projfile)
        call spproj%get_cavgs_stk(stk, ncls, smpd)
        call cavg_cnn_class_probs(params%infile, stk, ncls, smpd, pclass, empty, smpd_net)
        allocate(group(ncls))
        group = CAVG_CNN_GROUP_EMPTY
        do icls = 1,ncls
            if( .not. empty(icls) ) group(icls) = maxloc(pclass(:,icls), dim=1)
        enddo
        ! the table
        call fopen(funit, string(CAVG_CNN_PRED_TABLE), 'replace', 'write', iostat=ios)
        call fileiochk('cavg_cnn_predict; open '//CAVG_CNN_PRED_TABLE, ios)
        do icls = 1,ncls
            pop = 0.
            res = 0.
            if( icls <= spproj%os_cls2D%get_noris() )then
                if( spproj%os_cls2D%isthere(icls, 'pop') ) pop = spproj%os_cls2D%get(icls, 'pop')
                if( spproj%os_cls2D%isthere(icls, 'res') ) res = spproj%os_cls2D%get(icls, 'res')
            endif
            write(funit,'(A,I0,A,I0,A,A,A,A,A,A,A,A,A,A)') 'class=', icls, ' pop=', nint(pop), ' res=', fmt_real(real(res, dp), '(F10.2)'), &
                &' p_junk=', fmt_real(pclass(1,icls), '(ES16.9)'), ' p_ss=', fmt_real(pclass(2,icls), '(ES16.9)'), &
                &' p_likely=', fmt_real(pclass(3,icls), '(ES16.9)'), ' group=', trim(CAVG_CNN_GROUP_NAMES(group(icls)))
        enddo
        call fclose(funit)
        write(logfhandle,'(A,I0,A,F7.3,A,F7.3,A)') '>>> CLASS-AVERAGE CNN: ', ncls, ' classes at ', smpd, &
            &' A (network ', smpd_net, ' A)'
        do g = 1,size(CAVG_CNN_GROUP_NAMES)
            write(logfhandle,'(A,A,A,I0)') '>>> ', CAVG_CNN_GROUP_NAMES(g), ': ', count(group == g)
        enddo
        write(logfhandle,'(A)') '>>> TABLE: '//CAVG_CNN_PRED_TABLE
        call spproj%kill
        deallocate(empty, pclass, group)

      contains

        !> a number without blanks, for a key=value field
        function fmt_real( v, f ) result( str )
            real(dp),         intent(in)  :: v
            character(len=*), intent(in)  :: f
            character(len=:), allocatable :: str
            character(len=32)             :: buf
            write(buf, f) v
            str = trim(adjustl(buf))
        end function fmt_real

    end subroutine cavg_cnn_predict_exec

    !> The three probabilities (junk, secondary structure, likely protein) of each of the ncls class averages in stk at
    !! smpd under the weights file, and whether it is empty (no particles: probabilities 0); smpd_net is the network's
    subroutine cavg_cnn_class_probs( weights, stk, ncls, smpd, pclass, empty, smpd_net )
        class(string),         intent(in)               :: weights, stk
        integer,               intent(in)               :: ncls
        real,                  intent(in)               :: smpd
        real(dp),              allocatable, intent(out) :: pclass(:,:)
        logical,               allocatable, intent(out) :: empty(:)
        real,                  intent(out)              :: smpd_net
        type(cavg_cnn)           :: net
        type(cavg_cnn_input)     :: inp
        type(image), allocatable :: avgs(:), evens(:), odds(:)
        type(string)             :: ext
        real(sp),    allocatable :: raw(:,:,:,:), x(:,:,:,:), r(:,:,:)
        real(dp),    allocatable :: probs(:,:)
        integer,     allocatable :: keep(:)
        integer                  :: nch, n, box, icls, nkeep
        if( .not. file_exists(stk) ) THROW_HARD('class averages not found: '//stk%to_char())
        call net%new(weights)
        nch      = net%get_nchannels()
        box      = net%get_box()
        smpd_net = net%get_smpd()
        ! the class averages, and their halves for a network that reads them
        avgs = read_stk_into_imgarr(stk)
        if( size(avgs) /= ncls ) THROW_HARD('the class-average stack and the project disagree on the class count')
        ext = fname2ext(stk)
        ext = '.'//ext%to_char()
        if( nch == 3 )then
            evens = read_stk_into_imgarr(add2fbody(stk, ext, string('_even')))
            odds  = read_stk_into_imgarr(add2fbody(stk, ext, string('_odd')))
            if( size(evens) /= ncls .or. size(odds) /= ncls ) THROW_HARD('the halves hold another class count')
        endif
        r = avgs(1)%get_rmat()
        n = size(r,1)
        allocate(raw(n, n, nch, ncls), x(box, box, nch, ncls), empty(ncls))
        do icls = 1,ncls
            r = avgs(icls)%get_rmat()
            raw(:,:,1,icls) = r(:,:,1)
            if( nch == 3 )then
                r = evens(icls)%get_rmat()
                raw(:,:,2,icls) = r(:,:,1)
                r = odds(icls)%get_rmat()
                raw(:,:,3,icls) = r(:,:,1)
            endif
        enddo
        call dealloc_imgarr(avgs)
        if( nch == 3 )then
            call dealloc_imgarr(evens)
            call dealloc_imgarr(odds)
        endif
        call inp%new(n, smpd, box, net%get_smpd())
        !$omp parallel do default(shared) private(icls) schedule(dynamic,1)
        do icls = 1,ncls
            call inp%prepare(raw(:,:,:,icls), x(:,:,:,icls), empty(icls))
        enddo
        !$omp end parallel do
        deallocate(raw)
        ! the network on the classes with particles
        nkeep = count(.not. empty)
        keep  = pack([(icls, icls = 1,ncls)], .not. empty)
        allocate(probs(3, nkeep), pclass(3, ncls))
        if( nkeep > 0 ) call net%predict(x(:,:,:,keep), probs)
        pclass = 0.d0
        do icls = 1,nkeep
            pclass(:,keep(icls)) = probs(:,icls)
        enddo
        call net%kill
        call inp%kill
        deallocate(x, probs, keep)
    end subroutine cavg_cnn_class_probs

end module simple_cavg_cnn_predict
