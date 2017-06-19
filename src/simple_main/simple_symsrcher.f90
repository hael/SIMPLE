module simple_symsrcher
use simple_defs
use simple_math
use simple_cmdline, only: cmdline
use simple_params,  only: params
use simple_oris,    only: oris
use simple_image,   only: image
use simple_sym,     only: sym
implicit none

public :: dsym_cylinder
private

contains

    subroutine dsym_cylinder(p, dsym_os, cylinder)
        class(params), intent(in)    :: p
        class(oris),   intent(inout) :: dsym_os
        class(image),  intent(out)   :: cylinder
        type(image)          :: read_img, img_msk, dist_img, roavg_img, topview
        type(sym)            :: se
        real,    allocatable :: corrs(:), forsort(:), e2(:), radii(:)
        integer, allocatable :: labels(:)
        integer  :: halfnoris, cnt1, cnt2, i, l, noris
        real     :: cen1, cen2, sum1, sum2, sumvals
        real     :: minmax(2), width, height, sh1(3), ang
        if(p%pgrp(1:1).ne.'d' .and. p%pgrp(1:1).ne.'D')&
        &stop 'only intended for symmetry of order D; simple_symsrcher%dsym_cylinder'
        ! init
        noris = dsym_os%get_noris()
        allocate(corrs(noris),  source=-1.)
        allocate(radii(noris),  source=0.)
        allocate(e2(noris),     source=0.)
        allocate(labels(noris), source=0)
        call se%new(p%pgrp)
        call cylinder%new( [p%box, p%box, p%box], p%smpd)
        call read_img%new( [p%box, p%box, 1], p%smpd)
        call topview%new( [p%box, p%box, 1], p%smpd)
        call roavg_img%new([p%box, p%box, 1], p%smpd)
        call dist_img%new( [p%box, p%box, 1], p%smpd)
        call dist_img%cendist
        ang = 360. / real(se%get_nsym()/2)
        ! prep mask
        call img_msk%new([p%box,p%box,1],p%smpd)
        img_msk = 1.
        call img_msk%mask(p%msk, 'hard')
        ! centers, calculates self to rotational averages images & radii
        do i = 1, noris
            call read_img%read(p%stk,i)
            ! center image
            sh1 = read_img%center(p%cenlp, 'no', p%msk)
            call dsym_os%set(i, 'x', -sh1(1))
            call dsym_os%set(i, 'y', -sh1(2))
            ! rotational image
            call read_img%roavg(ang, roavg_img)
            call roavg_img%write('roavg.mrc',i)
            corrs(i) = read_img%real_corr(roavg_img, img_msk)
            ! radii
            call read_img%bp(0., p%cenlp)
            call read_img%mask(p%msk, 'soft')
            call read_img%bin('msk', p%msk)
            call read_img%write('bin.mrc',i)
            call read_img%mul(dist_img)
            minmax = read_img%minmax()
            radii(i) = min(p%msk, minmax(2))
        enddo
        call dsym_os%set_all('corr', corrs)
        ! indentify views along z and x-/y- axes
        forsort = corrs
        call hpsort(noris, forsort)
        halfnoris = nint(real(noris)/2.)
        cen1      = sum(forsort(1:halfnoris)) / real(halfnoris)
        cen2      = sum(forsort(halfnoris+1:noris)) / real(noris-halfnoris)
        sumvals   = sum(forsort)
        ! do 100 iterations of k-means
        do l = 1, 100
            sum1 = 0.
            cnt1 = 0
            do i = 1, noris
                if( (cen1-forsort(i))**2. < (cen2-forsort(i))**2. )then
                    cnt1 = cnt1 + 1
                    sum1 = sum1 + forsort(i)
                endif
            end do
            cnt2 = noris - cnt1
            sum2 = sumvals - sum1
            cen1 = sum1 / real(cnt1)
            cen2 = sum2 / real(cnt2)
        end do
        ! label views
        if(cen1 > cen2)then
            ! cen1: along z-axis
            labels = 2
            where( (cen1-corrs)**2. < (cen2-corrs)**2. )labels = 1
        else
            ! cen2: along z-axis
            labels = 1
            where( (cen1-corrs)**2. < (cen2-corrs)**2. )labels = 2
        endif
        e2 = 90.
        where(labels == 1) e2 = 0.
        call dsym_os%set_all('e2', e2)
        call dsym_os%set_all('class', real(labels))
        ! rotates input oris to asymmetric unit
        call se%rotall_to_asym(dsym_os)
        ! estimate height and cylinder radius (0.9 to account for overestimation)
        width  = 0.9 * (sum(radii, mask=(labels==1)) / real(count(labels==1)))
        height = 0.9 * (2. * sum(radii, mask=(labels==2)) / real(count(labels==2)))
        !call cylinder%bin_cylinder(width, height)
        ! dummy top view
        topview = 0.
        do i=1,noris
            if(labels(i)==1)then
                call read_img%read(p%stk,i)
                call read_img%rtsq(0., -dsym_os%get(i,'x'), -dsym_os%get(i,'y'))
                call topview%add(read_img)
            endif
        enddo
        call topview%div( real(count(labels==1)) )
        call topview%roavg(ang, roavg_img)
        topview = roavg_img
        call topview%norm
        call topview%mask(p%msk, 'soft')
        cylinder = 0.     
        do i=1,p%box
            if( abs( real(i-1)-real(p%box)/2. ) < height/2.)then
                call cylinder%set_slice( i, topview )
            endif
        enddo
        ! cleanup
        deallocate(corrs, e2, radii, labels)
        call dist_img%kill
        call img_msk%kill
        call read_img%kill
        call se%kill
        call roavg_img%kill
        call topview%kill
    end subroutine dsym_cylinder

end module simple_symsrcher

! public :: symsrch_master
! private

! type(ori)                     :: orientation
! type(oris)                    :: subgrps_oris
! type(sym)                     :: symobj
! integer                       :: i, s, n_subgrps, bestind, iptcl, nptcls, lfoo(3)
! character(len=3), allocatable :: subgrps(:)
! real                          :: shvec(3)
! character(len=32), parameter  :: SYMSHTAB = 'sym_3dshift.txt'

! contains

!     subroutine symsrch_master( cline, p, b, os )
!         use simple_filehandling,   only: nlines
!         use simple_projector_hlev, only: projvol
!         class(cmdline), intent(inout) :: cline
!         class(params),  intent(inout) :: p
!         class(build),   intent(inout) :: b
!         class(oris),    intent(out)   :: os
!         type(oris) :: oshift
!         type(ori)  :: o         
!         ! center volume
!         call b%vol%read(p%vols(1))
!         shvec = b%vol%center(p%cenlp,'no',p%msk)
!         if( p%l_distr_exec .and. p%part.eq.1 )then
!             ! writes shifts for distributed execution
!             call oshift%new(1)
!             call oshift%set(1,'x',shvec(1))
!             call oshift%set(1,'y',shvec(2))
!             call oshift%set(1,'z',shvec(3))
!             call oshift%write(trim(SYMSHTAB))
!             call oshift%kill
!         endif
!         ! main fork
!         if( p%compare .eq. 'no' )then
!             ! single symmetry search
!             ! generate projections
!             call b%vol%mask(p%msk, 'soft')
!             b%ref_imgs(1,:) = projvol(b%vol, b%e, p)
!             ! do the symmetry search
!             call single_symsrch( b, p, orientation )
!         else
            ! ! generate projections
            ! call b%vol%mask(p%msk, 'soft')
            ! b%ref_imgs(1,:) = projvol(b%vol, b%e, p)




            ! symmetry subgroup scan
            ! get subgroups
            ! subgrps      = b%se%get_all_subgrps()
            ! n_subgrps    = size(subgrps)
            ! subgrps_oris = oris(n_subgrps)
            ! ! loop over subgroups
            ! do i=1,n_subgrps
            !     symobj = sym( subgrps(i) )         ! current symmetry object
            !     call updates_tboxs( b, p, symobj ) ! updates symmetry variables in builder
            !     write(*,'(A,A3)') '>>> PERFORMING SYMMETRY SEARCH FOR POINT GROUP ', p%pgrp
            !     ! do the search
            !     call single_symsrch( b, p, orientation )
            !     call subgrps_oris%set_ori( i, orientation ) ! stashes result
            !     write(*,'(A)') '>>> '
            ! enddo
            ! ! aggregates correlations & probabilities
            ! call sym_probs( bestind )
            ! orientation = subgrps_oris%get_ori(bestind)

        ! endif
        ! if( p%l_distr_exec )then
        !     ! alles klar
        ! else ! we have 3D shifts to deal with
        !     ! rotate the orientations & transfer the 3d shifts to 2d
        !     shvec = -1.*shvec
        !     os    = oris(nlines(p%oritab))
        !     call os%read(p%oritab)
        !     if( cline%defined('state') )then
        !         do i=1,os%get_noris()
        !             s = nint(os%get(i, 'state'))
        !             if(s .ne. p%state)cycle
        !             call os%map3dshift22d(i, shvec)
        !             call os%rot(i,orientation)
        !             o = os%get_ori(i)
        !             call b%se%rot_to_asym(o)
        !             call os%set_ori(i, o) 
        !         end do
        !     else
        !         call os%map3dshift22d( shvec ) 
        !         call os%rot(orientation)
        !         call b%se%rotall_to_asym(os)
        !     endif
        ! endif
    ! end subroutine symsrch_master

    !>  \brief  Updates general and common lines toolboxes with current symmetry
    ! subroutine updates_tboxs( b, p, so )
    !     use simple_comlin, only: comlin
    !     type(build)  :: b
    !     type(params) :: p
    !     type(sym)    :: so
    !     integer      :: i, k, alloc_stat
    !     ! Update general toolbox sym functionality
    !     p%pgrp    = so%get_pgrp()
    !     call b%se%new(p%pgrp)
    !     p%nsym    = b%se%get_nsym()
    !     p%eullims = b%se%srchrange()
    !     call b%a%new(p%nptcls)
    !     call b%a%spiral(p%nsym, p%eullims)
    !     call b%e%new(p%nspace)
    !     call b%e%spiral(p%nsym, p%eullims)
    !     ! Update common lines toolbox
    !     if( p%pgrp /= 'c1' )then
    !         call b%a%symmetrize(p%nsym)
    !         do i=1,size(b%imgs_sym)
    !             call b%imgs_sym(i)%kill
    !         enddo
    !         do k=1,p%nstates
    !             do i=1,p%nspace
    !                 call b%ref_imgs(k,i)%kill
    !             enddo
    !         enddo
    !         deallocate(b%imgs_sym, b%ref_imgs)
    !         allocate( b%imgs_sym(1:p%nsym*p%nptcls), b%ref_imgs(p%nstates,p%nspace), stat=alloc_stat )
    !         call alloc_err( 'build_update_tboxs; simple_symsrch, 1', alloc_stat )            
    !         do i=1,p%nptcls*p%nsym
    !             call b%imgs_sym(i)%new([p%box,p%box,1], p%smpd)
    !         end do 
    !         b%clins = comlin(b%a, b%imgs_sym, p%lp)
    !     endif
    ! end subroutine updates_tboxs
    
    !>  \brief does common lines symmetry search over one symmetry pointgroup
    ! subroutine single_symsrch( b, p, best_o, fromto )
    !     use simple_comlin_srch ! use all in there
    !     type(build)       :: b
    !     type(params)      :: p
    !     type(ori)         :: best_o
    !     integer, optional :: fromto(2)
    !     integer           :: i, j, cnt
    !     ! expand over symmetry group
    !     cnt = 0
    !     do i=1,p%nptcls
    !         do j=1,p%nsym
    !             cnt = cnt+1
    !             b%imgs_sym(cnt) = b%ref_imgs(1,i)
    !             call b%imgs_sym(cnt)%fwd_ft
    !         end do
    !     end do
    !     ! search for the axis
    !     call comlin_srch_init( b, p, 'simplex', 'sym')
    !     call comlin_srch_symaxis( best_o, fromto )
    ! end subroutine single_symsrch
    
    !>  \brief  calculates probabilities and returns best axis index
    ! subroutine sym_probs( bestind )
    !     integer           :: bestind
    !     type(sym)         :: se, subse
    !     real, allocatable :: probs(:)
    !     character(len=3)  :: pgrp
    !     real              :: bestp
    !     integer           :: i, j, k, alloc_stat
    !     allocate( probs(n_subgrps), stat=alloc_stat)
    !     call alloc_err( 'sym_probs; simple_sym', alloc_stat )            
    !     probs = 0.
    !     do i=1,n_subgrps
    !         pgrp = subgrps(i)
    !         se   = sym(pgrp)
    !         do j=1,se%get_nsubgrp()
    !             subse = se%get_subgrp(j)
    !             do k=1,n_subgrps
    !                 if( subse%get_pgrp().eq.subgrps(k) )then
    !                     probs(i) = probs(i) + subgrps_oris%get(k,'corr') 
    !                 endif
    !             enddo
    !         enddo
    !         probs(i) = probs(i)/real(se%get_nsubgrp())
    !     enddo
    !     probs   = probs/sum(probs) 
    !     bestp   = -1.
    !     bestind = 0
    !     do i=1,n_subgrps
    !         write(*,'(A,A3,A,F4.2)')'>>> PROBABILITY FOR POINT-GROUP ', subgrps(i),' : ',probs(i)
    !         if( probs(i).gt.bestp )then
    !             bestp   = probs(i)
    !             bestind = i
    !         endif
    !     enddo
    !     deallocate(probs)
    ! end subroutine sym_probs
    