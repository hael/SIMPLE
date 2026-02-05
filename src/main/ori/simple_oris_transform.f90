submodule (simple_oris) simple_oris_transform
use simple_ori_api
use simple_ori, only: ori
implicit none
#include "simple_local_flags.inc"

contains

    module subroutine round_shifts( self )
        class(oris), intent(inout) :: self
        integer :: i
        do i=1,self%n
            call self%o(i)%round_shifts
        end do
    end subroutine round_shifts

    module subroutine introd_alig_err( self, angerr, sherr )
        class(oris), intent(inout) :: self
        real,        intent(in)    :: angerr, sherr
        real    :: x, y, e1, e2, e3
        integer :: i
        do i=1,self%n
            e1 = self%e1get(i)+ran3()*angerr-angerr/2.
            if( e1 > 360. ) e1 = e1-360.
            if( e1 < 0. ) e1 = e1+360.
            call self%o(i)%e1set(e1)
            e2 = self%e2get(i)+ran3()*angerr-angerr/2.
            if( e2 > 180. ) e2 = e2-180.
            if( e2 < 0. ) e2 = e2+180.
            call self%o(i)%e2set(e2)
            e3 = self%e3get(i)+ran3()*angerr-angerr/2.
            if( e3 > 360. ) e3 = e3-360.
            if( e3 < 0. ) e3 = e3+360.
            call self%o(i)%e3set(e3)
            x = self%o(i)%get('x')
            y = self%o(i)%get('y')
            x = x+ran3()*sherr-sherr/2.
            y = y+ran3()*sherr-sherr/2.
            call self%o(i)%set('x', x)
            call self%o(i)%set('y', y)
        end do
    end subroutine introd_alig_err

    module subroutine introd_ctf_err( self, dferr )
        class(oris), intent(inout) :: self
        real,        intent(in)    :: dferr
        real    :: dfx, dfy
        integer :: i
        do i=1,self%n
            if( self%o(i)%isthere('dfx') )then
                do
                    dfx = self%o(i)%get_dfx()+ran3()*dferr-dferr/2.
                    if( dfx > 0. ) exit
                end do
                call self%o(i)%set_dfx(dfx)
            endif
            if( self%o(i)%isthere('dfy') )then
                do
                    dfy = self%o(i)%get_dfy()+ran3()*dferr-dferr/2.
                    if( dfy > 0. ) exit
                end do
                call self%o(i)%set_dfy(dfy)
            endif
        end do
    end subroutine introd_ctf_err

    module subroutine rot_1( self, e )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: e
        type(ori) :: o_tmp
        integer   :: i
        do i=1,self%n
             call e%compose(self%o(i), o_tmp)
             call self%o(i)%set_euler(o_tmp%get_euler())
         end do
         call o_tmp%kill
    end subroutine rot_1

    module subroutine rot_2( self, i, e )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        class(ori),  intent(in)    :: e
        type(ori)                  :: o_tmp
        call e%compose(self%o(i), o_tmp)
        call self%o(i)%set_euler(o_tmp%get_euler())
        call o_tmp%kill
    end subroutine rot_2

    module subroutine rot_transp_1( self, e )
        class(oris), intent(inout) :: self
        class(ori),  intent(in)    :: e
        type(ori) :: o_tmp, e_transp
        integer   :: i
        e_transp = e
        call e_transp%transp
        do i=1,self%n
             call e_transp%compose(self%o(i), o_tmp)
             call self%o(i)%set_euler(o_tmp%get_euler())
         end do
         call o_tmp%kill
         call e_transp%kill
    end subroutine rot_transp_1

    module subroutine rot_transp_2( self, i, e )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: i
        class(ori),  intent(in)    :: e
        type(ori) :: o_tmp, e_transp
        e_transp = e
        call e_transp%transp
        call e_transp%compose(self%o(i), o_tmp)
        call self%o(i)%set_euler(o_tmp%get_euler())
        call o_tmp%kill
        call e_transp%kill
    end subroutine rot_transp_2

    module subroutine map3dshift22d_1( self, sh3d, state )
        class(oris),       intent(inout) :: self
        real,              intent(in)    :: sh3d(3)
        integer, optional, intent(in)    :: state
        integer :: i
        do i=1,self%n
            call self%map3dshift22d_2(i, sh3d, state)
        end do
    end subroutine map3dshift22d_1

    module subroutine map3dshift22d_2( self, i, sh3d, state )
        class(oris),       intent(inout) :: self
        integer,           intent(in)    :: i
        real,              intent(in)    :: sh3d(3)
        integer, optional, intent(in)    :: state
        integer :: mystate
        if( present(state) )then
            mystate = self%o(i)%get_state()
            if( mystate == state ) call self%o(i)%map3dshift22d(sh3d)
        else
            call self%o(i)%map3dshift22d(sh3d)
        endif
    end subroutine map3dshift22d_2

    module subroutine calc_avg_offset2D( self, class, offset_avg )
        class(oris),    intent(in)    :: self
        integer,        intent(in)    :: class
        real,           intent(inout) :: offset_avg(2)
        real(dp) :: avg(2)
        real     :: sh(2)
        integer  :: i, n
        n   = 0
        avg = 0.d0
        do i=1,self%n
            if( self%o(i)%isstatezero() ) cycle
            if( self%o(i)%get_class() == class )then
                call self%o(i)%calc_offset2D(sh)
                n   = n+1
                avg = avg + real(sh,dp)
            endif
        end do
        if( n > 0 )then
            avg        = avg / real(n,dp)
            offset_avg = real(avg)
        else
            offset_avg = 0.0
        endif
    end subroutine calc_avg_offset2D

    module subroutine calc_avg_offset3D( self, offset_avg, state )
        class(oris),       intent(inout)    :: self
        real,              intent(inout) :: offset_avg(3)
        integer, optional, intent(in)    :: state
        integer,    parameter :: N = 300
        integer,    parameter :: PROJDIRMINPOP = 10
        type(oris)            :: spiral
        type(ori)             :: o, oxy, oxz, oyz
        integer,  allocatable :: closest_proj(:)
        real(dp) :: avg(3), offset(3), w(3), sumw(3)
        real     :: sh3d(3)
        integer  :: i,j, istate, pop, npop
        istate = 1
        if( present(state) ) istate = state
        call spiral%new(N,.true.)
        call spiral%spiral
        call spiral%set_all2single('state', 1.)
        allocate(closest_proj(self%n),source=0)
        call self%remap_projs(spiral, closest_proj)
        avg  = 0.d0
        npop = 0
        o    = ori(.true.)
        oxy  = ori(.true.)
        oxz  = ori(.true.)
        oyz  = ori(.true.)
        call oxy%set_euler([ 0., 0.,0.])
        call oxz%set_euler([90.,90.,0.])
        call oyz%set_euler([ 0.,90.,0.])
        !$omp parallel do default(shared) private(i,j,pop,offset,o,sh3d,w,sumw) &
        !$omp schedule(static) proc_bind(close) reduction(+:npop,avg)
        do i = 1,N
            pop    = 0
            offset = 0.d0
            sumw   = 0.d0
            call spiral%get_ori(i, o)
            do j = 1,self%n
                if( self%o(j)%get_state() /= istate ) cycle
                if( closest_proj(j) /= i ) cycle
                call self%o(j)%compose2dshift3d(sh3d)
                w(1)   = abs(sin(oyz.euldist.self%o(j)))
                w(2)   = abs(sin(oxz.euldist.self%o(j)))
                w(3)   = abs(sin(oxy.euldist.self%o(j)))
                offset = offset + w*real(sh3d,dp)
                sumw   = sumw   + w**2
                pop    = pop    + 1
            enddo
            if( pop < PROJDIRMINPOP ) cycle
            where( sumw < 1.d-6 )
                offset = 0.d0
            else where
                offset = offset / sumw
            end where
            avg  = avg  + offset
            npop = npop + 1
        enddo
        !$omp end parallel do
        avg        = avg / real(npop,dp)
        offset_avg = real(avg)
        call spiral%kill
        call o%kill
        call oxy%kill
        call oxz%kill
        call oyz%kill
        deallocate(closest_proj)
    end subroutine calc_avg_offset3D

    module subroutine mirror2d( self )
        class(oris), intent(inout) :: self
        integer :: i
        do i = 1, self%n
            call self%o(i)%mirror2d()
        enddo
    end subroutine mirror2d

    module subroutine mirror3d( self )
        class(oris), intent(inout) :: self
        integer :: i
        do i = 1, self%n
            call self%o(i)%mirror3d()
        enddo
    end subroutine mirror3d

    module subroutine add_shift2class( self, class, sh2d )
        class(oris), intent(inout) :: self
        integer,     intent(in)    :: class
        real,        intent(in)    :: sh2d(2)
        integer :: i
        real    :: sh3d(3)
        sh3d(1:2) = sh2d
        sh3d(3)   = 0.
        do i=1,self%n
            if( self%o(i)%get_class() == class )then
                call self%o(i)%map3dshift22d(sh3d)
            endif
        end do
    end subroutine add_shift2class

end submodule simple_oris_transform
