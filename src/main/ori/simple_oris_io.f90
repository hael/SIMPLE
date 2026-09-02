!@descr: input/output routines for oris object
submodule (simple_oris) simple_oris_io
use simple_ori_api
use simple_ori, only: ori
use simple_jpg, only: jpg_img
implicit none
#include "simple_local_flags.inc"

contains

    module subroutine print( self, i )
        class(oris), intent(in) :: self
        integer,     intent(in) :: i
        call self%o(i)%print_ori()
    end subroutine print

    module subroutine print_matrices( self )
        class(oris), intent(inout) :: self
        integer :: i
        write(logfhandle,*) 'ORDER OF ROTATION MATRIX ELEMENTS: (1,1) (1,2) (1,3) (2,1) (2,2) (2,3) (3,1) (3,2) (3,3)'
        do i=1,self%n
            call self%o(i)%print_mat()
        end do
    end subroutine print_matrices

    module subroutine read( self, orifile, fromto, nst )
        class(oris),       intent(inout) :: self
        class(string),     intent(in)    :: orifile
        integer, optional, intent(in)    :: fromto(2)
        integer, optional, intent(out)   :: nst
        character(len=100) :: io_message
        integer :: file_stat, i, fnr, state, istart, iend
        if( .not. file_exists(orifile) )then
            THROW_HARD("the file you are trying to read: "//orifile%to_char()//' does not exist in cwd' )
        endif
        if( fname2ext(orifile) == 'bin' )then
            THROW_HARD('this method does not support binary files; read')
        endif
        io_message='No error'
        call fopen(fnr, FILE=orifile, STATUS='OLD', action='READ', iostat=file_stat,iomsg=io_message)
        call fileiochk("oris ; read ,Error when opening file for reading: "//orifile%to_char()//':'//trim(io_message), file_stat)
        if( present(nst) ) nst = 0
        if( present(fromto) )then
            istart = fromto(1)
            iend   = fromto(2)
            if(istart <      1) THROW_HARD('Invalid index; read')
            if(iend   > self%n) THROW_HARD('Invalid index; read')
        else
            istart = 1
            iend   = self%n
        endif
        do i = istart, iend
            call self%o(i)%read(fnr)
            if( present(nst) )then
                state = self%o(i)%get_state()
                nst   = max(1,max(state,nst))
            endif
        end do
        call fclose(fnr)
    end subroutine read

    module subroutine read_ctfparams_state_eo( self, ctfparamfile )
        class(oris),    intent(inout) :: self
        class(string),  intent(in)    :: ctfparamfile
        logical    :: params_are_there(10)
        integer    :: i
        type(oris) :: os_tmp
        if( .not. file_exists(ctfparamfile) )then
            THROW_HARD ("read_ctfparams_state_eo; The file you are trying to read: "//ctfparamfile%to_char()//' does not exist')
        endif
        if( ctfparamfile%has_substr('.bin') )then
            THROW_HARD('this method does not support binary files; read_ctfparams_state_eo')
        endif
        ! Parse into a hash-backed container regardless of the destination typing.
        ! Particle oris keep these keys in fixed slots that always report present, so
        ! reading into one would make every presence test below unconditionally true
        ! and a file lacking a column would overwrite the stored value with zero.
        call os_tmp%new(self%n, is_ptcl=.false.)
        call os_tmp%read(ctfparamfile)
        params_are_there(1)  = os_tmp%isthere('smpd')
        params_are_there(2)  = os_tmp%isthere('kv')
        params_are_there(3)  = os_tmp%isthere('cs')
        params_are_there(4)  = os_tmp%isthere('fraca')
        params_are_there(5)  = os_tmp%isthere('phshift')
        params_are_there(6)  = os_tmp%isthere('dfx')
        params_are_there(7)  = os_tmp%isthere('dfy')
        params_are_there(8)  = os_tmp%isthere('angast')
        params_are_there(9)  = os_tmp%isthere('state')
        params_are_there(10) = os_tmp%isthere('eo')
        do i=1,self%n
            if( params_are_there(1)  ) call self%set(i, 'smpd',    os_tmp%get(i, 'smpd')   )
            if( params_are_there(2)  ) call self%set(i, 'kv',      os_tmp%get(i, 'kv')     )
            if( params_are_there(3)  ) call self%set(i, 'cs',      os_tmp%get(i, 'cs')     )
            if( params_are_there(4)  ) call self%set(i, 'fraca',   os_tmp%get(i, 'fraca')  )
            if( params_are_there(5)  ) call self%set(i, 'phshift', os_tmp%get(i, 'phshift'))
            if( params_are_there(6)  ) call self%set_dfx(i,        os_tmp%get_dfx(i)       )
            if( params_are_there(7)  ) call self%set_dfy(i,        os_tmp%get_dfy(i)       )
            if( params_are_there(8)  ) call self%set(i, 'angast',  os_tmp%get(i, 'angast') )
            if( params_are_there(9)  ) call self%set(i, 'state',   os_tmp%get(i, 'state')  )
            if( params_are_there(10) ) call self%set(i, 'eo',      os_tmp%get(i, 'eo')     )
        end do
        call os_tmp%kill
    end subroutine read_ctfparams_state_eo

    module subroutine write_1( self, orifile, fromto )
        class(oris),       intent(in) :: self
        class(string),     intent(in) :: orifile
        integer, optional, intent(in) :: fromto(2)
        character(len=100) :: io_message
        integer            :: file_stat, fnr, i, ffromto(2), cnt
        ffromto(1) = 1
        ffromto(2) = self%n
        if( present(fromto) ) ffromto = fromto
        call fopen(fnr, orifile, status='REPLACE', action='WRITE', iostat=file_stat, iomsg=io_message)
        call fileiochk(' Error opening file for writing: '//orifile%to_char()//' ; '//trim(io_message), file_stat)
        cnt = 0
        do i=ffromto(1),ffromto(2)
            cnt = cnt + 1
            call self%o(i)%write(fnr)
        end do
        call fclose(fnr)
    end subroutine write_1

    module subroutine write_2( self, i, orifile )
        class(oris),   intent(inout) :: self
        class(string), intent(in)    :: orifile
        integer,       intent(in)    :: i
        integer :: fnr, file_stat
        call fopen(fnr, orifile, status='UNKNOWN', action='WRITE', position='APPEND', iostat=file_stat)
        call fileiochk( 'In: write_2, module: simple_oris.f90  opening '//orifile%to_char(), file_stat )
        call self%o(i)%write(fnr)
        call fclose(fnr)
    end subroutine write_2

    module subroutine write2bild( self, file )
        class(oris),   intent(inout) :: self
        class(string), intent(in)    :: file
        integer :: i,funit, file_stat
        call fopen(funit, file, status='REPLACE', action='WRITE',iostat=file_stat)
        call fileiochk( 'In: write2bild, module: simple_oris.f90  opening '//file%to_char(), file_stat )
        ! header
        write(funit,'(A)')".translate 0.0 0.0 0.0"
        write(funit,'(A)')".scale 10"
        write(funit,'(A)')".comment -- unit sphere --"
        write(funit,'(A)')".color 0.8 0.8 0.8"
        write(funit,'(A)')".sphere 0 0 0 1.0"
        write(funit,'(A)')".comment -- planes --"
        write(funit,'(A)')".color 0.3 0.3 0.3"
        write(funit,'(A)')".cylinder -0.02 0 0 0.02 0 0 1.02"
        write(funit,'(A)')".cylinder 0 -0.02 0 0 0.02 0 1.02"
        write(funit,'(A)')".cylinder 0 0 -0.02 0 0 0.02 1.02"
        write(funit,'(A)')".comment -- x-axis --"
        write(funit,'(A)')".color 1 0 0"
        write(funit,'(A)')".cylinder -1.5 0 0 1.5 0 0 0.02"
        write(funit,'(A)')".comment -- y-axis --"
        write(funit,'(A)')".color 0 1 0"
        write(funit,'(A)')".cylinder 0 -1.5 0 0 1.5 0 0.02"
        write(funit,'(A)')".comment -- z-axis --"
        write(funit,'(A)')".color 0 0 1"
        write(funit,'(A)')".cylinder 0 0 -1.5 0 0 1.5 0.02"
        write(funit,'(A)')".comment -- north pole --"
        write(funit,'(A)')".color 0 0 1"
        write(funit,'(A)')".sphere 0 0 1.5 0.1"
        ! body
        write(funit,'(A)')".color 0.4 0.4 0.4"
        do i=1,self%n
            call self%o(i)%write2bild(funit)
        enddo
        call fclose(funit)
    end subroutine write2bild

    module subroutine write_projdir_heatmap( self, state, nspace, fname, legend )
        class(oris),            intent(in) :: self
        integer,                intent(in) :: state, nspace
        class(string),          intent(in) :: fname
        type(string), optional, intent(in) :: legend
        integer, parameter :: JPEG_QUALITY      = 96
        integer, parameter :: ANGLE_LABEL_WIDTH = 24, ANGLE_LABEL_HEIGHT = 24
        integer, parameter :: TICK_LABEL_WIDTH  = 10, TICK_LABEL_HEIGHT  = 14, TICK_LABEL_ADVANCE = 12
        integer, parameter :: PLOT_WIDTH        = 540, PLOT_HEIGHT       = 270
        integer, parameter :: LEFT_MARGIN       = 72,  RIGHT_MARGIN      = 8
        integer, parameter :: TOP_MARGIN        = 32,  BOTTOM_MARGIN     = 72
        integer, parameter :: COLORBAR_GAP      = 24,  COLORBAR_WIDTH    = 16
        integer, parameter :: SCALE_LABEL_WIDTH = 80
        type(jpg_img)     :: jpg
        real, allocatable :: canvas(:,:), heatmap(:,:)
        character(len=64) :: legend_here
        character(len=16) :: min_label, max_label
        character(len=3)  :: pi_label
        character(len=4)  :: two_pi_label
        real    :: val, vmin, vmax, range_scale, scaled, color, black, white
        integer :: canvas_width, canvas_height, plot_x0, plot_x1, plot_y0, plot_y1
        integer :: colorbar_x0, colorbar_x1, x, y, status, red, green, blue
        logical :: error
        if( self%n < 1 ) THROW_HARD('Projection-direction heat-map requires a nonempty oris object')
        if( .not. self%isthere('proj') ) THROW_HARD('Projection-direction heat-map requires PROJ keys')
        if( nspace < 1 ) THROW_HARD('Projection-direction heat-map requires a positive nspace')
        if( state < 1 ) THROW_HARD('Projection-direction heat-map requires a positive state')
        legend_here = ''
        if( present(legend) ) legend_here = trim(legend%to_char())
        pi_label     = achar(92)//'PI'
        two_pi_label = '2'//pi_label
        allocate(heatmap(PLOT_WIDTH,PLOT_HEIGHT), source=0.)
        call calc_heatmap(heatmap, error)
        if( error ) return
        vmin = minval(heatmap)
        vmax = maxval(heatmap)
        if( vmax <= 0. ) THROW_HARD('Projection-direction heat-map contains no positive population')
        canvas_width  = LEFT_MARGIN + PLOT_WIDTH + COLORBAR_GAP + COLORBAR_WIDTH + SCALE_LABEL_WIDTH + RIGHT_MARGIN
        canvas_height = TOP_MARGIN + PLOT_HEIGHT + BOTTOM_MARGIN
        plot_x0       = LEFT_MARGIN + 1
        plot_x1       = plot_x0 + PLOT_WIDTH - 1
        plot_y0       = TOP_MARGIN + 1
        plot_y1       = plot_y0 + PLOT_HEIGHT - 1
        colorbar_x0   = plot_x1 + COLORBAR_GAP
        colorbar_x1   = colorbar_x0 + COLORBAR_WIDTH - 1
        black = rgb_code(0, 0, 0)
        white = rgb_code(255, 255, 255)
        allocate(canvas(canvas_width,canvas_height), source=white)
        if( vmax > vmin )then
            range_scale = 1. / (vmax - vmin)
            ! Heat-map theta index 1 is zero; the canvas has pi at its top edge.
            do y = plot_y0, plot_y1
                do x = plot_x0, plot_x1
                    val = heatmap(x - plot_x0 + 1, plot_y1 - y + 1)
                    scaled = max(0., min(1., (val - vmin) * range_scale))
                    call heat_color(scaled, red, green, blue)
                    canvas(x,y) = rgb_code(red, green, blue)
                enddo
            enddo
        else
            call heat_color(0.5, red, green, blue)
            canvas(plot_x0:plot_x1,plot_y0:plot_y1) = rgb_code(red, green, blue)
        endif
        do y = plot_y0, plot_y1
            scaled = real(plot_y1 - y) / real(max(1, PLOT_HEIGHT - 1))
            call heat_color(scaled, red, green, blue)
            color = rgb_code(red, green, blue)
            canvas(colorbar_x0:colorbar_x1,y) = color
        enddo
        call draw_hline(canvas, plot_x0 - 1, plot_x1 + 1, plot_y0 - 1, black)
        call draw_hline(canvas, plot_x0 - 1, plot_x1 + 1, plot_y1 + 1, black)
        call draw_vline(canvas, plot_y0 - 1, plot_y1 + 1, plot_x0 - 1, black)
        call draw_vline(canvas, plot_y0 - 1, plot_y1 + 1, plot_x1 + 1, black)
        call draw_hline(canvas, colorbar_x0 - 1, colorbar_x1 + 1, plot_y0 - 1, black)
        call draw_hline(canvas, colorbar_x0 - 1, colorbar_x1 + 1, plot_y1 + 1, black)
        call draw_vline(canvas, plot_y0 - 1, plot_y1 + 1, colorbar_x0 - 1, black)
        call draw_vline(canvas, plot_y0 - 1, plot_y1 + 1, colorbar_x1 + 1, black)
        call draw_vline(canvas, plot_y1 + 2, plot_y1 + 7, plot_x0, black)
        call draw_vline(canvas, plot_y1 + 2, plot_y1 + 7, plot_x1, black)
        call draw_hline(canvas, plot_x0 - 7, plot_x0 - 2, plot_y0, black)
        call draw_hline(canvas, plot_x0 - 7, plot_x0 - 2, plot_y1, black)
        call draw_smooth_tick_text(canvas, plot_x0, plot_y1 + 12, '0', black)
        call draw_smooth_tick_text(canvas, plot_x1 - smooth_tick_text_width(two_pi_label) + 1, &
            &plot_y1 + 12, two_pi_label, black)
        call draw_smooth_tick_text(canvas, plot_x0 - smooth_tick_text_width(pi_label) - 12, &
            &plot_y0 - 6, pi_label, black)
        call draw_smooth_tick_text(canvas, plot_x0 - smooth_tick_text_width('0') - 12, &
            &plot_y1 - 6, '0', black)
        call draw_smooth_angle_glyph(canvas, plot_x0 + (PLOT_WIDTH - ANGLE_LABEL_WIDTH) / 2, &
            &plot_y1 + 42, 'PHI', black)
        call draw_smooth_angle_glyph(canvas, 12, &
            &plot_y0 + (PLOT_HEIGHT - ANGLE_LABEL_HEIGHT) / 2, 'THETA', black)
        write(min_label,'(ES9.2E1)') vmin
        write(max_label,'(ES9.2E1)') vmax
        if( present(legend) )then
            call draw_text(canvas, plot_x0 + (PLOT_WIDTH - text_width(legend_here, 2)) / 2, &
                &12, trim(legend_here), 2, black)
        endif
        call draw_text(canvas, colorbar_x1 + 8, plot_y0 - 2, adjustl(max_label), 1, black)
        call draw_text(canvas, colorbar_x1 + 8, plot_y1 - 6, adjustl(min_label), 1, black)
        ! Write the canvas to a JPEG file
        status = jpg%writejpg(fname%to_char(), canvas, quality=JPEG_QUALITY, colorspec=3)
        if( status /= 0 ) THROW_WARN('Failed to write projection-direction heat-map JPEG')
        deallocate(canvas, heatmap)
    contains

        subroutine calc_heatmap( map, error )
            real,    intent(inout) :: map(:,:)
            logical, intent(out)   :: error
            real,    allocatable :: proj_pops(:), angles(:,:), stencil_weights(:)
            integer, allocatable :: stencil_dx(:), stencil_dy(:)
            logical, allocatable :: proj_seen(:)
            real    :: phi, theta, population, weight
            real    :: phi_step, theta_step, delta, sig2, cutoff2, dsq
            integer :: i, p, center_x, center_y, dx, dy, target_x, target_y, s, stencil_index
            integer :: halfwin_x, halfwin_y, stencil_capacity, stencil_size
            error = .false.
            allocate(proj_pops(nspace), angles(2,nspace), source=0.)
            allocate(proj_seen(nspace), source=.false.)
            do i = 1,self%n
                s = self%get_state(i)
                if( s /= state ) cycle
                p = self%get_proj(i)
                if( p == 0 )then
                    if( .not. self%isthere(i, 'proj') ) THROW_HARD('Active orientation is missing a PROJ key')
                    cycle
                endif
                if( p < 0 .or. p > nspace ) THROW_HARD('PROJ index is outside [1,nspace]')
                if( .not. proj_seen(p) )then
                    angles(1,p) = modulo(self%e1get(i), 360.)
                    angles(2,p) = max(0., min(180., self%e2get(i)))
                    proj_seen(p) = .true.
                endif
                weight = 1.
                if( self%isthere(i, 'w') ) weight = self%get(i, 'w')
                proj_pops(p) = proj_pops(p) + max(0., weight)
            enddo
            if( sum(proj_pops) <= 0.001 )then
                THROW_WARN('Projection-direction heat-map is empty: state '//int2str(state))
                error = .true.
                return
            endif
            ! Expected angular separation for nspace points distributed on a sphere.
            delta      = rad2deg(3.809 / sqrt(real(nspace)))
            cutoff2    = 6.0 * delta**2
            sig2       = 0.8 * delta**2
            phi_step   = 360. / real(PLOT_WIDTH)
            theta_step = 180. / real(PLOT_HEIGHT - 1)
            halfwin_x  = ceiling(sqrt(cutoff2) / phi_step)
            halfwin_y  = ceiling(sqrt(cutoff2) / theta_step)
            stencil_capacity = (2 * halfwin_x + 1) * (2 * halfwin_y + 1)
            allocate(stencil_dx(stencil_capacity), stencil_dy(stencil_capacity), stencil_weights(stencil_capacity))
            stencil_size = 0
            do dy = -halfwin_y, halfwin_y
                do dx = -halfwin_x, halfwin_x
                    dsq = (real(dx) * phi_step)**2 + (real(dy) * theta_step)**2
                    if( dsq > cutoff2 ) cycle
                    stencil_size = stencil_size + 1
                    stencil_dx(stencil_size) = dx
                    stencil_dy(stencil_size) = dy
                    stencil_weights(stencil_size) = exp(-0.5 * dsq / sig2)
                enddo
            enddo
            do p = 1,nspace
                population = proj_pops(p)
                if( .not. proj_seen(p) .or. population <= 0. ) cycle
                phi   = angles(1,p)
                theta = angles(2,p)
                center_x = 1 + modulo(nint(phi / phi_step), PLOT_WIDTH)
                center_y = 1 + nint(theta / theta_step)
                center_y = max(1, min(PLOT_HEIGHT, center_y))
                do stencil_index = 1, stencil_size
                    target_y = center_y + stencil_dy(stencil_index)
                    if( target_y < 1 .or. target_y > PLOT_HEIGHT ) cycle
                    target_x = 1 + modulo(center_x - 1 + stencil_dx(stencil_index), PLOT_WIDTH)
                    map(target_x,target_y) = map(target_x,target_y) + &
                        &population * stencil_weights(stencil_index)
                enddo
            enddo
            map = map * (maxval(proj_pops) / maxval(map))
        end subroutine calc_heatmap

        pure real function rgb_code( r, g, b )
            integer, intent(in) :: r, g, b
            integer :: packed
            packed = ishft(max(0, min(255, r)), 16) + ishft(max(0, min(255, g)), 8) + max(0, min(255, b))
            rgb_code = real(packed) / 16777215.0
        end function rgb_code

        pure subroutine heat_color( fraction, r, g, b )
            integer, parameter :: N_COLORS = 5
            integer, parameter :: HEAT_PALETTE(3,N_COLORS) = reshape([ &
                &239, 243, 255, & ! #EFF3FF, low
                &189, 215, 231, & ! #BDD7E7
                &107, 174, 214, & ! #6BAED6
                & 49, 130, 189, & ! #3182BD
                &  8,  81, 156  & ! #08519C, high
                &], [3,N_COLORS])
            real,    intent(in)  :: fraction
            integer, intent(out) :: r, g, b
            real :: local_fraction
            integer :: segment
            integer :: rgb0(3), rgb1(3)
            real :: fraction_here
            fraction_here = max(0., min(1., fraction))
            segment = min(N_COLORS - 1, int(real(N_COLORS - 1) * fraction_here) + 1)
            local_fraction = real(N_COLORS - 1) * fraction_here - real(segment - 1)
            rgb0 = HEAT_PALETTE(:,segment)
            rgb1 = HEAT_PALETTE(:,segment+1)
            r = nint(real(rgb0(1)) + local_fraction * real(rgb1(1) - rgb0(1)))
            g = nint(real(rgb0(2)) + local_fraction * real(rgb1(2) - rgb0(2)))
            b = nint(real(rgb0(3)) + local_fraction * real(rgb1(3) - rgb0(3)))
        end subroutine heat_color

        pure subroutine draw_hline( bitmap, x0, x1, ypos, ink )
            real,    intent(inout) :: bitmap(:,:)
            integer, intent(in)    :: x0, x1, ypos
            real,    intent(in)    :: ink
            integer :: first, last
            if( ypos < 1 .or. ypos > size(bitmap,2) ) return
            first = max(1, min(x0, x1))
            last  = min(size(bitmap,1), max(x0, x1))
            if( first <= last ) bitmap(first:last,ypos) = ink
        end subroutine draw_hline

        pure subroutine draw_vline( bitmap, y0, y1, xpos, ink )
            real,    intent(inout) :: bitmap(:,:)
            integer, intent(in)    :: y0, y1, xpos
            real,    intent(in)    :: ink
            integer :: first, last
            if( xpos < 1 .or. xpos > size(bitmap,1) ) return
            first = max(1, min(y0, y1))
            last  = min(size(bitmap,2), max(y0, y1))
            if( first <= last ) bitmap(xpos,first:last) = ink
        end subroutine draw_vline

        pure subroutine fill_rect( bitmap, x0, y0, x1, y1, ink )
            real,    intent(inout) :: bitmap(:,:)
            integer, intent(in)    :: x0, y0, x1, y1
            real,    intent(in)    :: ink
            integer :: xa, xb, ya, yb
            xa = max(1, min(x0, x1))
            xb = min(size(bitmap,1), max(x0, x1))
            ya = max(1, min(y0, y1))
            yb = min(size(bitmap,2), max(y0, y1))
            if( xa <= xb .and. ya <= yb ) bitmap(xa:xb,ya:yb) = ink
        end subroutine fill_rect

        pure subroutine draw_smooth_angle_glyph( bitmap, xpos, ypos, glyph, ink )
            integer, parameter :: SUPERSAMPLE = 4
            real, parameter :: INNER_RADIUS_SQ = 0.88**2, OUTER_RADIUS_SQ = 1.12**2
            real,             intent(inout) :: bitmap(:,:)
            integer,          intent(in)    :: xpos, ypos
            character(len=*), intent(in)    :: glyph
            real,             intent(in)    :: ink
            real :: coverage, ellipse_radius_sq, sample_x, sample_y
            integer :: ix, iy, sample_ix, sample_iy, covered_samples
            logical :: is_glyph_sample
            do iy = 0, ANGLE_LABEL_HEIGHT - 1
                if( ypos + iy < 1 .or. ypos + iy > size(bitmap,2) ) cycle
                do ix = 0, ANGLE_LABEL_WIDTH - 1
                    if( xpos + ix < 1 .or. xpos + ix > size(bitmap,1) ) cycle
                    covered_samples = 0
                    do sample_iy = 1, SUPERSAMPLE
                        sample_y = (real(iy) + (real(sample_iy) - 0.5) / real(SUPERSAMPLE)) / &
                            &real(ANGLE_LABEL_HEIGHT)
                        do sample_ix = 1, SUPERSAMPLE
                            sample_x = (real(ix) + (real(sample_ix) - 0.5) / real(SUPERSAMPLE)) / &
                                &real(ANGLE_LABEL_WIDTH)
                            ellipse_radius_sq = ((sample_x - 0.5) / 0.35)**2 + &
                                &((sample_y - 0.5) / 0.35)**2
                            is_glyph_sample = ellipse_radius_sq >= INNER_RADIUS_SQ .and. &
                                &ellipse_radius_sq <= OUTER_RADIUS_SQ
                            select case(glyph)
                            case('PHI')
                                is_glyph_sample = is_glyph_sample .or. abs(sample_x - 0.5) <= 0.06
                            case('THETA')
                                is_glyph_sample = is_glyph_sample .or. &
                                    &(abs(sample_y - 0.5) <= 0.05 .and. abs(sample_x - 0.5) <= 0.35)
                            end select
                            if( is_glyph_sample ) covered_samples = covered_samples + 1
                        enddo
                    enddo
                    coverage = real(covered_samples) / real(SUPERSAMPLE * SUPERSAMPLE)
                    if( coverage > 0.0 ) bitmap(xpos+ix,ypos+iy) = &
                        &coverage * ink + (1.0 - coverage) * bitmap(xpos+ix,ypos+iy)
                enddo
            enddo
        end subroutine draw_smooth_angle_glyph

        pure integer function smooth_tick_text_width( text )
            character(len=*), intent(in) :: text
            integer :: icharacter, nglyphs, token_length
            icharacter = 1
            nglyphs = 0
            do while( icharacter <= len_trim(text) )
                token_length = greek_token_length(text, icharacter)
                if( token_length > 0 )then
                    icharacter = icharacter + token_length
                else
                    icharacter = icharacter + 1
                endif
                nglyphs = nglyphs + 1
            enddo
            smooth_tick_text_width = max(0, nglyphs * TICK_LABEL_ADVANCE - &
                &(TICK_LABEL_ADVANCE - TICK_LABEL_WIDTH))
        end function smooth_tick_text_width

        pure subroutine draw_smooth_tick_text( bitmap, xpos, ypos, text, ink )
            real,             intent(inout) :: bitmap(:,:)
            integer,          intent(in)    :: xpos, ypos
            character(len=*), intent(in)    :: text
            real,             intent(in)    :: ink
            integer :: icharacter, token_length, xorigin
            xorigin = xpos
            icharacter = 1
            do while( icharacter <= len_trim(text) )
                token_length = greek_token_length(text, icharacter)
                if( token_length == 3 )then
                    call draw_smooth_tick_glyph(bitmap, xorigin, ypos, 'PI', ink)
                    icharacter = icharacter + token_length
                else
                    call draw_smooth_tick_glyph(bitmap, xorigin, ypos, text(icharacter:icharacter), ink)
                    icharacter = icharacter + 1
                endif
                xorigin = xorigin + TICK_LABEL_ADVANCE
            enddo
        end subroutine draw_smooth_tick_text

        pure subroutine draw_smooth_tick_glyph( bitmap, xpos, ypos, glyph, ink )
            integer, parameter :: SUPERSAMPLE = 4
            real, parameter :: ZERO_INNER_SQ = 0.50, ZERO_OUTER_SQ = 1.00
            real, parameter :: STROKE_RADIUS_SQ = 0.75**2
            real,             intent(inout) :: bitmap(:,:)
            integer,          intent(in)    :: xpos, ypos
            character(len=*), intent(in)    :: glyph
            real,             intent(in)    :: ink
            real :: coverage, ellipse_radius_sq, sample_x, sample_y
            integer :: ix, iy, sample_ix, sample_iy, covered_samples
            logical :: is_glyph_sample
            do iy = 0, TICK_LABEL_HEIGHT - 1
                if( ypos + iy < 1 .or. ypos + iy > size(bitmap,2) ) cycle
                do ix = 0, TICK_LABEL_WIDTH - 1
                    if( xpos + ix < 1 .or. xpos + ix > size(bitmap,1) ) cycle
                    covered_samples = 0
                    do sample_iy = 1, SUPERSAMPLE
                        sample_y = real(iy) + (real(sample_iy) - 0.5) / real(SUPERSAMPLE)
                        do sample_ix = 1, SUPERSAMPLE
                            sample_x = real(ix) + (real(sample_ix) - 0.5) / real(SUPERSAMPLE)
                            select case(glyph)
                            case('0')
                                ellipse_radius_sq = ((sample_x - 5.0) / 3.6)**2 + &
                                    &((sample_y - 7.0) / 5.5)**2
                                is_glyph_sample = ellipse_radius_sq >= ZERO_INNER_SQ .and. &
                                    &ellipse_radius_sq <= ZERO_OUTER_SQ
                            case('2')
                                is_glyph_sample = &
                                    &segment_distance_sq(sample_x, sample_y, 1.4, 3.4, 2.7, 1.8) <= STROKE_RADIUS_SQ .or. &
                                    &segment_distance_sq(sample_x, sample_y, 2.7, 1.8, 6.4, 1.8) <= STROKE_RADIUS_SQ .or. &
                                    &segment_distance_sq(sample_x, sample_y, 6.4, 1.8, 8.4, 3.2) <= STROKE_RADIUS_SQ .or. &
                                    &segment_distance_sq(sample_x, sample_y, 8.4, 3.2, 7.8, 5.4) <= STROKE_RADIUS_SQ .or. &
                                    &segment_distance_sq(sample_x, sample_y, 7.8, 5.4, 1.7, 11.7) <= STROKE_RADIUS_SQ .or. &
                                    &segment_distance_sq(sample_x, sample_y, 1.7, 11.7, 8.4, 11.7) <= STROKE_RADIUS_SQ
                            case('PI')
                                is_glyph_sample = &
                                    &segment_distance_sq(sample_x, sample_y, 1.2, 2.3, 8.8, 2.3) <= STROKE_RADIUS_SQ .or. &
                                    &segment_distance_sq(sample_x, sample_y, 3.0, 2.3, 2.7, 12.0) <= STROKE_RADIUS_SQ .or. &
                                    &segment_distance_sq(sample_x, sample_y, 7.0, 2.3, 7.3, 12.0) <= STROKE_RADIUS_SQ
                            case default
                                is_glyph_sample = .false.
                            end select
                            if( is_glyph_sample ) covered_samples = covered_samples + 1
                        enddo
                    enddo
                    coverage = real(covered_samples) / real(SUPERSAMPLE * SUPERSAMPLE)
                    if( coverage > 0.0 ) bitmap(xpos+ix,ypos+iy) = &
                        &coverage * ink + (1.0 - coverage) * bitmap(xpos+ix,ypos+iy)
                enddo
            enddo
        end subroutine draw_smooth_tick_glyph

        pure real function segment_distance_sq( x, y, x0, y0, x1, y1 )
            real, intent(in) :: x, y, x0, y0, x1, y1
            real :: dx, dy, segment_length_sq, position, nearest_x, nearest_y
            dx = x1 - x0
            dy = y1 - y0
            segment_length_sq = dx**2 + dy**2
            if( segment_length_sq <= TINY )then
                segment_distance_sq = (x - x0)**2 + (y - y0)**2
                return
            endif
            position = max(0.0, min(1.0, ((x - x0) * dx + (y - y0) * dy) / segment_length_sq))
            nearest_x = x0 + position * dx
            nearest_y = y0 + position * dy
            segment_distance_sq = (x - nearest_x)**2 + (y - nearest_y)**2
        end function segment_distance_sq

        pure integer function text_width( text, scale )
            character(len=*), intent(in) :: text
            integer,          intent(in) :: scale
            integer :: icharacter, nglyphs, token_length
            icharacter = 1
            nglyphs = 0
            do while( icharacter <= len_trim(text) )
                token_length = greek_token_length(text, icharacter)
                if( token_length > 0 )then
                    icharacter = icharacter + token_length
                else
                    icharacter = icharacter + 1
                endif
                nglyphs = nglyphs + 1
            enddo
            text_width = max(0, nglyphs * 6 * scale - scale)
        end function text_width

        pure subroutine draw_text( bitmap, xpos, ypos, text, scale, ink )
            real,             intent(inout) :: bitmap(:,:)
            integer,          intent(in)    :: xpos, ypos, scale
            character(len=*), intent(in)    :: text
            real,             intent(in)    :: ink
            integer :: icharacter, row, column, token_length, xorigin
            integer :: rows(7)
            xorigin = xpos
            icharacter = 1
            do while( icharacter <= len_trim(text) )
                token_length = greek_token_length(text, icharacter)
                if( token_length > 0 )then
                    rows = glyph_rows(greek_token_glyph(text, icharacter))
                    icharacter = icharacter + token_length
                else
                    rows = glyph_rows(upper_character(text(icharacter:icharacter)))
                    icharacter = icharacter + 1
                endif
                do row = 1, 7
                    do column = 0, 4
                        if( btest(rows(row), 4 - column) )then
                            call fill_rect(bitmap, xorigin + column * scale, ypos + (row - 1) * scale, &
                                &xorigin + (column + 1) * scale - 1, ypos + row * scale - 1, ink)
                        endif
                    enddo
                enddo
                xorigin = xorigin + 6 * scale
            enddo
        end subroutine draw_text

        pure integer function greek_token_length( text, position )
            character(len=*), intent(in) :: text
            integer,          intent(in) :: position
            greek_token_length = 0
            if( token_matches(text, position, 'THETA') )then
                greek_token_length = 6
            elseif( token_matches(text, position, 'PHI') )then
                greek_token_length = 4
            elseif( token_matches(text, position, 'PI') )then
                greek_token_length = 3
            endif
        end function greek_token_length

        pure character function greek_token_glyph( text, position )
            character(len=*), intent(in) :: text
            integer,          intent(in) :: position
            select case(greek_token_length(text, position))
            case(3)
                greek_token_glyph = '~'
            case(4)
                greek_token_glyph = '{'
            case(6)
                greek_token_glyph = '}'
            case default
                greek_token_glyph = ' '
            end select
        end function greek_token_glyph

        pure logical function token_matches( text, position, name )
            character(len=*), intent(in) :: text, name
            integer,          intent(in) :: position
            integer :: icharacter
            token_matches = .false.
            if( position < 1 .or. position + len(name) > len_trim(text) ) return
            if( iachar(text(position:position)) /= 92 ) return
            do icharacter = 1, len(name)
                if( upper_character(text(position+icharacter:position+icharacter)) /= &
                    &name(icharacter:icharacter) ) return
            enddo
            token_matches = .true.
        end function token_matches

        pure character function upper_character( character_in )
            character, intent(in) :: character_in
            integer :: code
            code = iachar(character_in)
            if( code >= iachar('a') .and. code <= iachar('z') )then
                upper_character = achar(code - iachar('a') + iachar('A'))
            else
                upper_character = character_in
            endif
        end function upper_character

        pure function glyph_rows( glyph ) result( rows )
            character, intent(in) :: glyph
            integer :: rows(7)
            rows = 0
            select case(glyph)
            case('0'); rows = [14,17,19,21,25,17,14]
            case('1'); rows = [4,12,4,4,4,4,14]
            case('2'); rows = [14,17,1,2,4,8,31]
            case('3'); rows = [30,1,1,14,1,1,30]
            case('4'); rows = [2,6,10,18,31,2,2]
            case('5'); rows = [31,16,16,30,1,1,30]
            case('6'); rows = [14,16,16,30,17,17,14]
            case('7'); rows = [31,1,2,4,8,8,8]
            case('8'); rows = [14,17,17,14,17,17,14]
            case('9'); rows = [14,17,17,15,1,1,14]
            case('A'); rows = [14,17,17,31,17,17,17]
            case('B'); rows = [30,17,17,30,17,17,30]
            case('C'); rows = [14,17,16,16,16,17,14]
            case('D'); rows = [30,17,17,17,17,17,30]
            case('E'); rows = [31,16,16,30,16,16,31]
            case('F'); rows = [31,16,16,30,16,16,16]
            case('G'); rows = [14,17,16,23,17,17,15]
            case('H'); rows = [17,17,17,31,17,17,17]
            case('I'); rows = [14,4,4,4,4,4,14]
            case('J'); rows = [7,2,2,2,2,18,12]
            case('K'); rows = [17,18,20,24,20,18,17]
            case('L'); rows = [16,16,16,16,16,16,31]
            case('M'); rows = [17,27,21,21,17,17,17]
            case('N'); rows = [17,25,21,19,17,17,17]
            case('O'); rows = [14,17,17,17,17,17,14]
            case('P'); rows = [30,17,17,30,16,16,16]
            case('Q'); rows = [14,17,17,17,21,18,13]
            case('R'); rows = [30,17,17,30,20,18,17]
            case('S'); rows = [15,16,16,14,1,1,30]
            case('T'); rows = [31,4,4,4,4,4,4]
            case('U'); rows = [17,17,17,17,17,17,14]
            case('V'); rows = [17,17,17,17,17,10,4]
            case('W'); rows = [17,17,17,21,21,21,10]
            case('X'); rows = [17,17,10,4,10,17,17]
            case('Y'); rows = [17,17,10,4,4,4,4]
            case('Z'); rows = [31,1,2,4,8,16,31]
            case('-'); rows = [0,0,0,31,0,0,0]
            case('+'); rows = [0,4,4,31,4,4,0]
            case('.'); rows = [0,0,0,0,0,12,12]
            case(':'); rows = [0,12,12,0,12,12,0]
            case('('); rows = [2,4,8,8,8,4,2]
            case(')'); rows = [8,4,2,2,2,4,8]
            case('/'); rows = [1,2,2,4,8,8,16]
            case('~'); rows = [31,10,10,10,10,10,9] ! \PI token
            case('{'); rows = [4,14,21,21,14,4,4]    ! Greek phi token
            case('}'); rows = [4,14,17,31,17,14,4]  ! Greek theta token
            end select
        end function glyph_rows

    end subroutine write_projdir_heatmap

end submodule simple_oris_io
