! various mathematical subroutines and functions
module simple_math
use simple_defs
use simple_error, only: allocchk, simple_exception
implicit none

interface is_a_number
    module procedure is_a_number_1
    module procedure is_a_number_2
end interface

interface is_zero
    module procedure is_zero_0
    module procedure is_zero_1
    module procedure is_zero_2
end interface is_zero

interface is_gt_zero
    module procedure is_gt_zero_0
    module procedure is_gt_zero_1
    module procedure is_gt_zero_2
end interface is_gt_zero

interface is_equal
    module procedure is_equal_1
    module procedure is_equal_2
end interface is_equal

interface check4nans3D
    module procedure check4nans3D_1
    module procedure check4nans3D_2
end interface

interface check4nans2D
    module procedure check4nans2D_1
    module procedure check4nans2D_2
end interface

interface check4nans
    module procedure check4nans_1
    module procedure check4nans_2
end interface

interface is_nan
    module procedure is_a_number_1
    module procedure is_a_number_2
    module procedure is_a_number_3
end interface

interface is_even
    module procedure is_even_1
    module procedure is_even_2
end interface

interface cosedge
    module procedure cosedge_1
    module procedure cosedge_2
end interface

interface cosedge_inner
    module procedure cosedge_inner_1
    module procedure cosedge_inner_2
end interface

interface hardedge
    module procedure hardedge_1
    module procedure hardedge_2
    module procedure hardedge_3
    module procedure hardedge_4
end interface

interface hardedge_inner
    module procedure hardedge_inner_1
    module procedure hardedge_inner_2
    module procedure hardedge_inner_3
    module procedure hardedge_inner_4
end interface

interface find
    module procedure find_1
    module procedure find_2
end interface

interface locate
   module procedure locate_1
   module procedure locate_2
end interface

interface hpsort
    module procedure hpsort_1
    module procedure hpsort_2
    module procedure hpsort_3
    module procedure hpsort_4
    module procedure hpsort_5
end interface

interface reverse
    module procedure reverse_iarr
    module procedure reverse_rarr
end interface

interface nvoxfind
    module procedure nvoxfind_1
    module procedure nvoxfind_2
end interface

interface rad2deg
    module procedure rad2deg_1
    module procedure rad2deg_2
end interface

interface csq
    module procedure csq_1
    module procedure csq_2
end interface

interface sqwin_1d
    module procedure sqwin_1d_1
    module procedure sqwin_1d_2
end interface

interface sqwin_2d
    module procedure sqwin_2d_1
    module procedure sqwin_2d_2
end interface

interface sqwin_3d
    module procedure sqwin_3d_1
    module procedure sqwin_3d_2
end interface

logical, parameter,private :: warn=.false.

interface norm_2
    module procedure norm_2_sp
    module procedure norm_2_dp
end interface norm_2

interface deg2rad
    module procedure deg2rad_sp
    module procedure deg2rad_dp
end interface deg2rad

interface myacos
    module procedure myacos_sp
    module procedure myacos_dp
end interface myacos

!!!!!!!!ADDED BY CHIARA!!!!!!!!!!
interface vis_mat
    module procedure vis_2Dreal_mat
    module procedure vis_2Dinteger_mat
    module procedure vis_3Dreal_mat
    module procedure vis_3Dinteger_mat
end interface vis_mat

!!!!!!ADDED BY CHIARA!!!!!
interface mode
    module procedure mode_1
    module procedure mode_2
end interface mode

interface otsu
    module procedure otsu_1
    module procedure otsu_2
end interface otsu

contains

    !> \brief nvoxfind_1  to find the volume in number of voxels, given molecular weight
    !! \details SI units \f$ \si{\kilo\dalton}= \SI{1e-10}{\metre} \f$, one dalton is defined as 1/12 of the mass of an atom of Carbon 12, or 1 amu \f$ \SI{1.66053892173e-27}{\kilo\gram} \f$
    !!
    !!  Protein density \f$ \rho_{\mathrm{prot}} \f$ is defined as \f$ \si{1.43}{\gram\per\centimetre\cubed} \f$
    !!   see  Quillin and Matthews, Accurate calculation of the density of proteins, DOI 10.1107/S090744490000679X
    !! protein density in \f$ \si{\gram\per\angstrom\cubed},\ \rho_{\mathrm{prot}\si{\angstrom}} = \rho_{\mathrm{prot}} \num{1e-24} \f$
    !! unit voxel volume v in \f$ \si{\per\angstrom\cubed},\ v  = \mathrm{smpd}^3 \f$
    !! mass of protein in Da, \f$ M_\mathrm{prot Da} = \mathrm{mwkda}\times\num{1e3} \f$
    !! mass of protein in kg \f$ M_\mathrm{kg prot} = M_\mathrm{prot Da}*M_{Hydrogen (\mathrm{kg/Da}) \f$
    !! mass of protein in g \f$ M_\mathrm{g prot} = ((mwkda*1e3)*one_da)*1e3 \f$
    !! therefore number of voxels in protein is:
    !! \f[ \begin{align*} N_{\mathrm{vox}} &=& \frac{ M_{\mathrm{molecule}} }{ \rho_{\mathrm{prot}}
    !! \times \frac{1}{v^3}\\ & \sim& \nint {frac{\mathrm{mwkda}\times\num{1e3}\times\num{1.66e-27)}
    !! {\num{1.43}\times\num{1e-24}} \times \frac{1}{\mathrm{smpd}^3}}} \end{align*}
    !! \f]
    !! we must also assume the number of voxels are discrete, hence we must round to the nearest integer
    !! \param smpd sampling distance in angstroms (SI unit \f$ \si{\angstrom}= \si{1e-10}{\metre} \f$)
    !! \param mwkda molecular weight \f$ M_{\mathrm{molecule}} \f$ in kDa
    pure function nvoxfind_1( smpd, mwkda ) result( nvox )
        real, intent(in) :: smpd             !< sampling distance
        real, intent(in) :: mwkda            !< molecular weight
        integer          :: nvox             !< nr of voxels
        double precision , parameter :: prot_d = 1.43d0            ! g/cm**3
        double precision , parameter :: one_da = 1.66053892173e-27 ! kg/Da
        ! nvox = nint( ( ( (mwkda * 1e3) * one_da ) * 1e3 ) / ( prot_d * 1e-24 * (smpd**3) )  )
        nvox = nint((mwkda*one_da*1e30) / (prot_d * (smpd**3.)))
    end function nvoxfind_1

    !>   to find the volume in number of voxels, given molecular weight
    pure function nvoxfind_2( smpd, mwkda, dens ) result( nvox )
        real, intent(in) :: smpd             !< sampling distance
        real, intent(in) :: mwkda            !< molecular weight
        real, intent(in) :: dens             !< density in Da/A3  \f$ \si{\dalton\per\angstrom\cubed} \f$
        real             :: vol_per_pix, vol !< volume per pixel & volume
        integer          :: nvox             !< nr of voxels
        vol_per_pix = smpd**3.
        vol = (mwkda*1000.)/dens
        nvox = nint(vol/vol_per_pix)
    end function nvoxfind_2

    !>   converts between radians and degrees
    elemental function deg2rad_sp( deg ) result( rad )
        real, intent(in) :: deg  !< angle (degrees)
        real             :: rad  !< angle (radians)
        rad = (deg/180.)*pi
    end function deg2rad_sp

    !>   converts between radians and degrees
    elemental function deg2rad_dp( deg ) result( rad )
        real(dp), intent(in) :: deg  !< angle (degrees)
        real(dp)             :: rad  !< angle (radians)
        rad = (deg/180._dp)*dpi
    end function deg2rad_dp

    !>   converts from radians to degrees
    elemental function rad2deg_1( rad ) result( deg )
        real(sp), intent(in) :: rad  !< angle (radians)
        real(sp)             :: deg  !< angle (degrees)
        deg = (rad/PI)*180.
    end function rad2deg_1

    !>   converts from radians to degrees
    elemental function rad2deg_2( rad ) result( deg )
        real(dp), intent(in) :: rad  !< angle (radians)
        real(dp)             :: deg  !< angle (degrees)
        deg = (rad/DPI)*180.d0
    end function rad2deg_2

    !>   to check if val is even
    elemental logical function is_even_1( val )
        integer, intent(in) :: val  !< query val
        is_even_1 = mod(val,2) == 0
    end function is_even_1

    !>   to check if all vals in array are even
    pure function is_even_2( arr ) result( yep )
        integer, intent(in) :: arr(:)     !< query vector
        logical :: yep
        logical :: test(size(arr))
        integer :: i
        test = .false.
        do i=1,size(arr)
            test(i) = is_even_1(arr(i))
        end do
        yep = all(test)
    end function is_even_2

    !>    is for checking the numerical soundness of an vector
    subroutine check4nans3D_1( arr )
        real, intent(in)  :: arr(:,:,:)   !< query vector
        real, allocatable :: arr1d(:)
        arr1d = pack(arr, .true.)
        call check4nans_1(arr1d)
        deallocate(arr1d)
    end subroutine check4nans3D_1

    !>    is for checking the numerical soundness of an vector
    subroutine check4nans3D_2( arr )
        complex, intent(in)  :: arr(:,:,:)    !< query vector
        complex, allocatable :: arr1d(:)
        arr1d = pack(arr, .true.)
        call check4nans_2(arr1d)
        deallocate(arr1d)
    end subroutine check4nans3D_2

    !>    is for checking the numerical soundness of an vector
    subroutine check4nans2D_1( arr )
        real, intent(in)  :: arr(:,:)   !< query vector
        real, allocatable :: arr1d(:)
        arr1d = pack(arr, .true.)
        call check4nans_1(arr1d)
        deallocate(arr1d)
    end subroutine check4nans2D_1

    !>    is for checking the numerical soundness of an vector
    subroutine check4nans2D_2( arr )
        complex, intent(in)  :: arr(:,:)   !< query vector
        complex, allocatable :: arr1d(:)
        arr1d = pack(arr, .true.)
        call check4nans_2(arr1d)
        deallocate(arr1d)
    end subroutine check4nans2D_2

    !>    is for checking the numerical soundness of an vector
    subroutine check4nans_1( arr )
        real, intent(in) :: arr(:)   !< query vector
        integer :: i, n_nans
        n_nans = 0
        do i=1,size(arr)
            if( is_a_number(arr(i)) )then
                ! alles gut
            else
                n_nans = n_nans+1
            endif
        end do
        if( n_nans > 0 )then
            write(logfhandle,*) 'found NaNs in inputted vector; simple_math::check4nans_1', n_nans
        endif
    end subroutine check4nans_1

    !>    is for checking the numerical soundness of an vector
    subroutine check4nans_2( arr )
        complex, intent(in) :: arr(:)   !< query vector
        integer :: i, n_nans
        n_nans = 0
        do i=1,size(arr)
            if( is_a_number_2(arr(i)) )then
                ! alles gut
            else
                n_nans = n_nans+1
            endif
        end do
        if( n_nans > 0 )then
            write(logfhandle,*) 'found NaNs in inputted vector; simple_math::check4nans_2', n_nans
        endif
    end subroutine check4nans_2

    !>  returns true if the argument is odd
    elemental logical function is_odd(i)
        integer, intent(in) :: i    !< query value
        is_odd = btest(i,0)
    end function is_odd

    !>   for rounding to closest even
    elemental function round2even( val ) result( ev )
        real, intent(in) :: val            !< query value
        integer :: ev, rounded, remainer
        rounded = nint(val)
        remainer = nint(val-real(rounded))
        if( mod(rounded,2) == 0 )then
            ev = rounded
        else
            if( abs(real(rounded-1)-val) <= abs(real(rounded+1)-val) )then
                ev = rounded-1
            else
                ev = rounded+1
            endif
        endif
    end function round2even

    !>   for rounding to closest even
    elemental function round2odd( val ) result( ev )
        real, intent(in) :: val          !< query value
        integer :: ev, rounded, remainer
        rounded = nint(val)
        remainer = nint(val-real(rounded))
        if( mod(rounded,2) /= 0 )then
            ev = rounded
        else
            if( abs(real(rounded-1)-val) <= abs(real(rounded+1)-val) )then
                ev = rounded-1
            else
                ev = rounded+1
            endif
        endif
    end function round2odd

    !>   get the resolution in angstrom, given angle and diameter
    !! \param ang,diam angular resolution (degrees) and diameter (\f$\si{\angstrom}\f$)
    pure function resang( ang, diam ) result( res )
        real, intent(in)  :: ang, diam
        real :: res                      !< spatial resolution (\f$\si{\per\angstrom}\f$)
        res = (ang/360.)*(pi*diam)
    end function resang

     !>   checking for is_a_number
    pure logical function is_a_number_1( number )
        real, intent(in) :: number  !< input variable for checking

        is_a_number_1 = .true.
        if( number > 0. )then
        else if( number <= 0. )then
        else
            is_a_number_1 = .false.
        endif
    end function is_a_number_1

    !>   validity check of complex number (so that it is not nan)
    pure logical function is_a_number_2( complex_number )
        complex, intent(in) :: complex_number !< input variable for checking

        is_a_number_2 = is_a_number_1(real(complex_number)) .and. is_a_number_1(aimag(complex_number))
    end function is_a_number_2

     !>   checking for is_a_number
    pure logical function is_a_number_3( number )
        real(dp), intent(in) :: number  !< input variable for checking

        is_a_number_3 = .true.
        if( number > 0. )then
        else if( number <= 0. )then
        else
            is_a_number_3 = .false.
        endif
    end function is_a_number_3

    !>   validity check of complex number (so that it is not nan)
    pure logical function is_a_number_4( complex_number )
        complex(dp), intent(in) :: complex_number !< input variable for checking

        is_a_number_4 = is_a_number_3(real(complex_number)) .and. is_a_number_3(aimag(complex_number))
    end function is_a_number_4

    !>   to check if val is zero
    elemental logical function is_zero_0( val )
        integer, intent(in) :: val  !< query val
        is_zero_0 = abs(val) < TINY
    end function is_zero_0

    !>   to check if val is zero
    elemental logical function is_zero_1( val )
        real, intent(in) :: val  !< query val
        is_zero_1 = abs(val) == 0
    end function is_zero_1

        !>   to check if val is zero
    elemental logical function is_zero_2( val )
        real(8), intent(in) :: val  !< query val
        is_zero_2 = abs(val) < DTINY
    end function is_zero_2

    !>   to check if val is zero
    elemental logical function is_gt_zero_0( val )
        integer, intent(in) :: val  !< query val
        is_gt_zero_0 = val > 0
    end function

    !>   to check if val is zero
    elemental logical function is_gt_zero_1( val )
        real, intent(in) :: val  !< query val
        is_gt_zero_1 = val > TINY
    end function is_gt_zero_1

    !>   to check if val is zero
    elemental logical function is_gt_zero_2( val )
        real(8), intent(in) :: val  !< query val
        is_gt_zero_2 = val > DTINY
    end function is_gt_zero_2

    !>   to check if val is zero
    elemental logical function is_equal_1( val1 , val2)
        real, intent(in) :: val1, val2  !< query val
        is_equal_1 = abs(val1-val2) < TINY
    end function

        !>   to check if val is zero
    elemental logical function is_equal_2( val1, val2 )
        real(8), intent(in) :: val1, val2  !< query val
        is_equal_2 = abs(val1-val2) < DTINY
    end function is_equal_2

    !>   to put the which element (if it exists) last in the array,
    !!         swapping it with its present position
    subroutine put_last( which, arr )
        integer, intent(in)    :: which  !< index
        integer, intent(inout) :: arr(:) !< array for modification
        integer :: pos, i, sz
        sz = size(arr)
        if( sz > 1 )then
            pos = 0
            do i=1,sz
                if( arr(i) == which ) pos = i
            end do
            if( pos > 0 )then
                if( pos == sz ) return ! it is last already
                ! swap
                arr(pos) = arr(sz)
                arr(sz)  = which
            endif
        endif
    end subroutine put_last

    !>   reverses an integer array
    subroutine reverse_iarr( iarr )
        integer, intent(inout) :: iarr(:) !< array for modification
        integer                :: i, j, iswap, sz, en
        sz = size(iarr,1)
        if( sz < 2 )then
            return
        endif
        if( mod(sz,2) == 0 )then
            en = sz/2+1
        else
            en = (sz+1)/2+1
        endif
        j = 0
        do i = sz,en,-1
            j = j+1
            iswap   = iarr(j)
            iarr(j) = iarr(i)
            iarr(i) = iswap
        end do
    end subroutine reverse_iarr

    !>   reverses a real array
    subroutine reverse_rarr( rarr )
        real, intent(inout) :: rarr(:) !< array for modification
        integer             :: i, j, sz, en
        real                :: rswap
        sz = size(rarr,1)
        if( sz < 2 )then
            return
        endif
        if( mod(sz,2) == 0 )then
            en = sz/2+1
        else
            en = (sz+1)/2+1
        endif
        j = 0
        do i = sz,en,-1
            j = j+1
            rswap   = rarr(j)
            rarr(j) = rarr(i)
            rarr(i) = rswap
        end do
    end subroutine reverse_rarr

    !>   implements the sortmeans algorithm
    subroutine sortmeans( dat, maxits, means, labels )
        real,                 intent(in)  :: dat(:)    !< array for input
        integer,              intent(in)  :: maxits    !< limit sort
        real,                 intent(out) :: means(:)  !< array for output
        integer, allocatable, intent(out) :: labels(:) !< labels for output
        logical, allocatable :: mask(:)
        integer :: ncls, ndat, clssz, i, j, cnt_means, loc(1), changes
        real, allocatable :: dat_sorted(:)
        ncls = size(means)
        ndat = size(dat)
        if( allocated(labels) ) deallocate(labels)
        allocate(  mask(ndat), labels(ndat), stat=alloc_stat )
        if(alloc_stat /= 0) call allocchk("sortmeans; simple_math", alloc_stat)
        ! initialization by sorting
        dat_sorted = dat
        call hpsort(dat_sorted)
        clssz = int(real(ndat)/real(ncls))
        cnt_means = 0
        do i=1,ndat,clssz
            cnt_means = cnt_means + 1
            means(cnt_means) = dat_sorted(i)
            if( cnt_means == ncls ) exit
        end do
        ! the kmeans step
        labels = 1
        do j=1,maxits
            changes = 0
            ! find closest
            do i=1,ndat
                loc = minloc((means-dat(i))**2.0)
                if( labels(i) /= loc(1) ) changes = changes + 1
                labels(i) = loc(1)
            end do
            ! update means
            do i=1,ncls
                where( labels == i )
                    mask = .true.
                else where
                    mask = .false.
                end where
                means(i) = sum(dat, mask=mask)/real(count(mask))
            end do
            if( changes == 0 ) exit
        end do
        deallocate(dat_sorted, mask)
    end subroutine sortmeans

    !>   one-dimensional cyclic index generation
    !! \param lims limits
    pure function cyci_1d( lims, i ) result( ind )
        integer, intent(in) :: lims(2), i !< input var
        integer :: ind, del
        ind = i
        if( ind > lims(2) )then
            del = ind-lims(2)
            ind = lims(1)+del-1
        else if( ind < lims(1) )then
            del = lims(1)-ind
            ind = lims(2)-del+1
        endif
    end function cyci_1d

    !>  as per cyci_1d, assumes lower limit equals 1
    pure function cyci_1d_static( ulim, i ) result( ind )
        integer, intent(in) :: ulim, i !< input vars upper limit and index
        integer :: ind ! return index
        ind = merge(ulim+i, i , i<1)
        ind = merge(i-ulim,ind ,i>ulim)
    end function cyci_1d_static

    !>   4 shifting variables
    subroutine shft(a,b,c,d)
        real, intent(out)   :: a   !< new pos 1
        real, intent(inout) :: b,c !< pos 2 and 3
        real, intent(in)    :: d   !< pos4
        a = b
        b = c
        c = d
    end subroutine shft

    !> one-dimensional symmetric hard window
    pure subroutine sqwin_1d_1( x, winsz, win )
        real,    intent(in)  :: x      !< input point
        real,    intent(in)  :: winsz  !< window size
        integer, intent(out) :: win(2) !< window
        integer :: iwinsz
        win(:) = nint(x)
        iwinsz = ceiling(winsz - 0.5)
        win(1) = win(1)-iwinsz
        win(2) = win(2)+iwinsz
    end subroutine sqwin_1d_1

    !> one-dimensional symmetric hard window with limits
    pure subroutine sqwin_1d_2( x, winsz, lims, win )
        real,    intent(in)  :: x       !< input point
        real,    intent(in)  :: winsz   !< window size
        integer, intent(in)  :: lims(2) !< bounds
        integer, intent(out) :: win(2)  !< window
        integer :: iwinsz
        win(:) = nint(x)
        iwinsz = ceiling(winsz - 0.5)
        win(1) = max(lims(1), win(1) - iwinsz)
        win(2) = min(lims(2), win(2) + iwinsz)
    end subroutine sqwin_1d_2

    !> two-dimensional symmetric hard window
    pure subroutine sqwin_2d_1( x, y, winsz, win )
        real,    intent(in)  :: x,y      !< input point
        real,    intent(in)  :: winsz    !< window size
        integer, intent(out) :: win(2,2) !< window
        call sqwin_1d_1(x, winsz, win(1,:))
        call sqwin_1d_1(y, winsz, win(2,:))
    end subroutine sqwin_2d_1

    !> two-dimensional symmetric hard window with limits
    pure subroutine sqwin_2d_2( x, y, winsz, lims, win )
        real,    intent(in)  :: x,y       !< input point
        real,    intent(in)  :: winsz     !< window size
        integer, intent(in)  :: lims(2,2) !< bounds
        integer, intent(out) :: win(2,2)  !< window
        call sqwin_1d_2(x,winsz,lims(1,:), win(1,:))
        call sqwin_1d_2(y,winsz,lims(2,:), win(2,:))
    end subroutine sqwin_2d_2

    !> three-dimensional symmetric hard window
    pure subroutine sqwin_3d_1( x, y, z, winsz, win )
        real,    intent(in)  :: x,y,z    !< input point
        real,    intent(in)  :: winsz    !< window size
        integer, intent(out) :: win(3,2) !< window
        call sqwin_1d_1(x, winsz, win(1,:))
        call sqwin_1d_1(y, winsz, win(2,:))
        call sqwin_1d_1(z, winsz, win(3,:))
    end subroutine sqwin_3d_1

    !> three-dimensional symmetric hard window with limits
    pure subroutine sqwin_3d_2( x, y, z, winsz, lims, win )
        real,    intent(in)  :: x,y,z     !< input point
        real,    intent(in)  :: winsz     !< window size
        integer, intent(in)  :: lims(3,2) !< bounds
        integer, intent(out) :: win(3,2)  !< window
        call sqwin_1d_2(x,winsz,[lims(1,1), lims(1,2)], win(1,:))
        call sqwin_1d_2(y,winsz,[lims(2,1), lims(2,2)], win(2,:))
        call sqwin_1d_2(z,winsz,[lims(3,1), lims(3,2)], win(3,:))
    end subroutine sqwin_3d_2

    ! USEFUL MATHEMATICAL FUNCTIONS

    !>   returns acos with the argument's absolute value limited to 1.
    !!         this seems to be necessary due to small numerical inaccuracies.
    pure function myacos_sp( arg ) result( r )
        real, intent(in) :: arg     !< input (radians)
        real             :: r, x, y
        x = min(1.,abs(arg))
        y = sign(x,arg)
        r = acos(y)
    end function myacos_sp

    !>   returns acos with the argument's absolute value limited to 1.
    !!         this seems to be necessary due to small numerical inaccuracies.
    pure function myacos_dp( arg ) result( r )
        real(dp), intent(in) :: arg     !< input (radians)
        real(dp)             :: r, x, y
        x = min(1._dp,abs(arg))
        y = sign(x,arg)
        r = acos(y)
    end function myacos_dp

    !>   sinc function
    function sinc( x ) result( r )
        real, intent(in) :: x       !< input (radians)
        real             :: r, arg
        if( abs(x) < 1e-8 ) then
            r = 1.
        else
            arg = pi*x
            r = sin(arg)/(arg)
        endif
    end function sinc

    !>   is a truncated Gaussian window function
    function gauwfun( x, alpha ) result( w )
        real, intent(in) :: x, alpha
        real             :: w, var
        var = alpha*2.
        w = 2.**(-(x/var)**2.)
    end function gauwfun

    real function gaussian2D( center_coords, x, y, xsigma, ysigma )
        real, intent(in) :: center_coords(2), x, y, xsigma, ysigma
        gaussian2D = exp( -((x - center_coords(1))**2.0 / (2.0 * xsigma * xsigma) +&
                            (y - center_coords(2))**2.0 / (2.0 * ysigma * ysigma)) )
    end function gaussian2D

    ! COMPLEX STUFF

    !>   is for calculating the phase angle of a Fourier component
    !! \return phase phase angle only meaningful when cabs(comp) is well above the noise level
    elemental function phase_angle( comp ) result( phase )
        complex, intent(in) :: comp !< Fourier component
        real :: nom, denom, phase
        nom = aimag(comp)
        denom = real(comp)
        if( is_zero(denom) )then
           if( is_zero(nom) )then
                phase = 0.
            else if( nom > 0. )then
                phase = pi/2.
            else
                phase = -pi/2.
            endif
            return
        endif
        phase = atan(nom/denom)
    end function phase_angle

    !>   is for complex squaring
    elemental function csq_1( a ) result( sq )
        complex(sp), intent(in) :: a !< complx component
        real(sp) :: sq, x, y, frac
        x = abs(real(a))
        y = abs(aimag(a))
        if( is_zero(x)) then
            sq = y*y
        else if( is_zero(y) ) then
           sq = x*x
        else if( x > y ) then
            frac = y/x
            sq = x*x*(1.+frac*frac)
        else
            frac = x/y
            sq = y*y*(1.+frac*frac)
        endif
    end function csq_1

    !>  is for double complex squaring
    elemental function csq_2( a ) result( sq )
        complex(dp), intent(in) :: a !< complx component
        real(dp) :: sq, x, y, frac
        x = abs(real(a))
        y = abs(aimag(a))
        if( is_zero(x)) then
            sq = y*y
        else if( is_zero(y) ) then
            sq = x*x
        else if( x > y ) then
            frac = y/x
            sq = x*x*(1.+frac*frac)
        else
            frac = x/y
            sq = y*y*(1.+frac*frac)
        endif
    end function csq_2

    !>   is for calculating complex arg/abs/modulus, from numerical recipes
    elemental function mycabs( a ) result( myabs )
        complex, intent(in) :: a      !< complx component
        real                :: myabs, x, y, frac
        x = abs(real(a))
        y = abs(aimag(a))
        if( is_zero(x) ) then
            myabs = y
        else if( is_zero(y)  ) then
           myabs = x
        else if( x > y ) then
            frac = y/x
            myabs = x*sqrt(1.+frac*frac)
        else
            frac = x/y
            myabs = y*sqrt(1.+frac*frac)
        endif
    end function mycabs

    pure function nextPow2(v) result(w)
        integer, intent(in) :: v
        integer :: w
        w = v - 1
        w = IOR(w, ISHFT(w,-1))
        w = IOR(w, ISHFT(w,-2))
        w = IOR(w, ISHFT(w,-4))
        w = IOR(w, ISHFT(w,-8))
        w = IOR(w, ISHFT(w,-16))
        w = w + 1
    end function nextPow2

    ! edge functions

    !>   two-dimensional hard edge
    !! \f$r^2 < x^2+y^2\f$.
    !! \param x x position
    !! \param y y position
    !! \param mskrad masking radius
    !! \return w on or off
    !!
    pure function hardedge_1( x, y, mskrad ) result( w )
        real,intent(in) :: x, y, mskrad
        real :: w
        w = 1.
        if( x * x + y * y > mskrad * mskrad ) w = 0.
    end function hardedge_1

    !>   three-dimensional hard edge
    !! \f$r^2 < x^2+y^2+z^2\f$.
    !! \param x x position
    !! \param y y position
    !! \param mskrad masking radius
    !! \return w on or off
    !!
    pure function hardedge_2( x, y, z, mskrad ) result( w )
        real,intent(in) :: x, y, z, mskrad
        real :: w
        w = 1.
        if( x * x + y * y + z * z > mskrad * mskrad ) w = 0.
    end function hardedge_2

   pure function hardedge_3( x, y, mskrad ) result( w )
        integer,intent(in) :: x, y
        real, intent(in)   :: mskrad
        real :: w
        w = 1.
        if( real(x * x + y * y) > mskrad * mskrad ) w = 0.
    end function hardedge_3
    !!
    pure function hardedge_4( x, y, z, mskrad ) result( w )
        integer,intent(in) :: x, y, z
        real, intent(in)   :: mskrad
        real :: w
        w = 1.
        if( real(x * x + y * y + z * z) > mskrad * mskrad ) w = 0.
    end function hardedge_4

    !>   two-dimensional hard edge
    !! \f$r < \sqrt{x^2+y^2}\f$.
    !! \return w on or off
    !!
    pure function hardedge_inner_1( x, y, mskrad ) result( w )
        real,intent(in) :: x, y, mskrad
        real :: w
        w = 0.
        if( x * x + y * y > mskrad * mskrad ) w = 1.
    end function hardedge_inner_1

    pure function hardedge_inner_2( x, y, mskrad ) result( w )
        integer,intent(in) :: x, y
        real, intent(in) ::mskrad
        real :: w
        w = 0.
        if( real(x * x + y * y) > mskrad * mskrad ) w = 1.
    end function hardedge_inner_2

    !>   three-dimensional hard edge
    pure function hardedge_inner_3( x, y, z, mskrad ) result( w )
        real,intent(in) :: x, y, z, mskrad
        real :: w
        w = 0.
        if( x * x + y * y + z * z > mskrad * mskrad ) w = 1.
    end function hardedge_inner_3

    pure function hardedge_inner_4( x, y, z, mskrad ) result( w )
        integer,intent(in) :: x, y,z
        real, intent(in) ::mskrad
        real :: w
        w = 0.
        if( real(x * x + y * y + z * z) > mskrad * mskrad ) w = 1.
    end function hardedge_inner_4

    !>   two-dimensional gaussian edge
    !! \param x x position
    !! \param y y position
   pure function cosedge_1( x, y, box, mskrad ) result( w )
       real, intent(in)    :: x, y     !< input points
       integer, intent(in) :: box      !< window size
       real, intent(in)    :: mskrad   !< mask radius
        real                :: w, rad, width, maxrad
        maxrad = real(box/2)
        rad    = sqrt(x**2.+y**2.)
        width  = 2.*(maxrad-mskrad)
        w      = 1.
        if( rad .ge. maxrad )then
            w = 0.
        else if( rad .ge. (maxrad-width) )then
            w = (cos(((rad-(maxrad-width))/width)*pi)+1.)/2.
        endif
    end function cosedge_1

    !>   three-dimensional gaussian edge
    !! \f$r = \cos{(1+{(\pi{r - (d+2m)/(d-2m)})})}\f$.
    !! \param x x position
    !! \param y y position
    pure function cosedge_2( x, y, z, box, mskrad ) result( w )
        real, intent(in)    :: x, y, z   !< input points
        integer, intent(in) :: box       !< window size
        real, intent(in)    :: mskrad    !< mask radius
        real                :: w, rad, maxrad, width
        maxrad = real(box/2)
        rad    = sqrt(x**2.+y**2.+z**2.)
        width  = 2.*(maxrad-mskrad)
        w      = 1.
        if( rad .ge. maxrad )then
            w = 0.
        else if( rad .ge. (maxrad-width) )then
            w = (cos(((rad-(maxrad-width))/width)*pi)+1.)/2.
        endif
    end function cosedge_2

    !> \brief  two-dimensional gaussian edge
    !! \param x x position
    !! \param y y position
    pure function cosedge_inner_1( x, y, width, mskrad ) result( w )
        real, intent(in) :: x, y, width  !< input points and width
        real, intent(in) :: mskrad       !< mask radius
        real             :: w, rad
        rad = sqrt(x**2.+y**2.)
        if( rad .lt. mskrad-width )then
            w = 0.
        else if( rad .gt. mskrad )then
            w = 1.
        else
            w = (1.+cos(pi*(mskrad-rad)/width))/2.0
        endif
    end function cosedge_inner_1

    !>   two-dimensional gaussian edge
    !! \param x x position
    !! \param y y position
    !! \param z z position
    pure function cosedge_inner_2( x, y, z, width, mskrad ) result( w )
        real, intent(in) :: x, y, z, width !< inner mask radius
        real, intent(in) :: mskrad !< mask radius
        real             :: w, rad
        rad = sqrt(x**2.+y**2.+z**2.)
        if( rad .lt. mskrad-width )then
            w = 0.
        else if( rad .gt. mskrad )then
            w = 1.
        else
            w = (1.+cos(pi*(mskrad-rad)/width))/2.0
        endif
    end function cosedge_inner_2

    ! FOURIER STUFF

    !>   is for working out the fourier dimension
    pure function fdim( d ) result( fd )
        integer, intent(in) :: d !< dimension
        integer :: fd
        if(is_even(d))then
            fd = d/2+1
        else
            fd = (d-1)/2+1
        endif
    end function fdim

    !>   calculates the resolution values given corrs and res params
    !! \param corrs Fourier shell correlations
    !! \param res resolution value
    subroutine get_resolution( corrs, res, fsc05, fsc0143 )
        real, intent(in)  :: corrs(:), res(:) !<  corrs Fourier shell correlation
        real, intent(out) :: fsc05, fsc0143   !<  fsc05 resolution at FSC=0.5,  fsc0143 resolution at FSC=0.143
        integer           :: n, ires0143, ires05
        n = size(corrs)
        ires0143 = 1
        do while( ires0143 <= n )
            if( corrs(ires0143) >= 0.143 )then
                ires0143 = ires0143 + 1
                cycle
            else
                exit
            endif
        end do
        ires0143 = ires0143-1
        if( ires0143 == 0 )then
            fsc0143 = 0.
        else
            fsc0143 = res(ires0143)
        endif
        ires05 = 1
        do while( ires05 <= n )
            if( corrs(ires05) >= 0.5 )then
                ires05 = ires05+1
                cycle
            else
                exit
            endif
        end do
        ires05 = ires05-1
        if( ires05 == 0 )then
            fsc05 = 0.
        else
            fsc05 = res(ires05)
        endif
    end subroutine get_resolution

    subroutine phaseplate_correct_fsc( fsc, find_plate )
        real,    intent(inout) :: fsc(:)
        integer, intent(out)   :: find_plate
        logical, allocatable :: peakpos(:)
        integer :: k, n
        real    :: peakavg
        n = size(fsc)
        ! find FSC peaks
        peakpos = peakfinder_2(fsc)
        ! filter out all peaks FSC < 0.5
        where(fsc < 0.5) peakpos = .false.
        ! calculate peak average
        peakavg = sum(fsc, mask=peakpos)/real(count(peakpos))
        ! identify peak with highest frequency
        do k=n,1,-1
            if( peakpos(k) )then
                find_plate = k
                exit
            endif
        end do
        ! replace with average FSC peak value up to last peak (find_plate)
        do k=1,find_plate
            fsc(k) = peakavg
        end do
    end subroutine phaseplate_correct_fsc

    !>   returns the Fourier index of the resolution limit at corr
    function get_lplim_at_corr( fsc, corr ) result( k )
        real, intent(in) :: fsc(:), corr
        integer          :: n, k, h
        n = size(fsc)
        if( n < 3 )then
            call simple_exception('nonconforming size of fsc array; get_lplim_at_corr', __FILENAME__ , __LINE__)
        endif
        k = n-1
        do h=3,n-1
            if( fsc(h) >= corr )then
                cycle
            else
                k = h
                exit
            endif
        end do
    end function get_lplim_at_corr

    !>   returns the Fourier index of resolution \f$ (\si{\per\angstrom}) \f$
    !< \param smpd pixel size
    !! \param res resolution \f$ (\si{\angstrom}) \f$
    integer pure function calc_fourier_index( res, box, smpd )
        real, intent(in)    :: res, smpd
        integer, intent(in) :: box       !< box size
        calc_fourier_index = nint((real(box-1)*smpd)/res)
    end function calc_fourier_index

    !>   returns the Fourier index of res
    real pure function calc_lowpass_lim( find, box, smpd )
        integer, intent(in) :: find, box !< box size
        real, intent(in)    :: smpd      !< smpd pixel size \f$ (\si{\angstrom}) \f$
        calc_lowpass_lim = (real(box-1)*smpd)/real(find)
    end function calc_lowpass_lim

    !>  \brief get array of resolution steps
    !> get_res
    !! \return  res
    !!
    function get_resarr( box, smpd ) result( res )
        integer, intent(in) :: box
        real,    intent(in) :: smpd
        real, allocatable   :: res(:)
        integer :: n, k, alloc_stat
        n = fdim(box) - 1
        allocate( res(n), stat=alloc_stat )
        if(alloc_stat .ne. 0)call allocchk('In: get_res, module: simple_math', alloc_stat)
        do k=1,n
            res(k) = calc_lowpass_lim(k, box, smpd)
        end do
    end function get_resarr

    ! OLD
    ! subroutine calc_norm_bfac_weights( box, bfac, smpd, kweights )
    !     integer, intent(in)  :: box
    !     real,    intent(in)  :: bfac, smpd
    !     real,    intent(out) :: kweights(0:box)
    !     real    :: bfac_sc, res
    !     integer :: k
    !     bfac_sc = bfac / 4.
    !     do k=1,box
    !         res         = real(k) / (real(box) * smpd)
    !         kweights(k) = max(0., exp(-bfac_sc * res * res))
    !     end do
    !     ! normalize
    !     kweights(1:box) = kweights(1:box) / sum(kweights(1:box))
    !     kweights(0)     = 1.
    ! end subroutine calc_norm_bfac_weights

    subroutine calc_norm_bfac_weights( box, bfac, smpd, kweights )
        integer, intent(in)  :: box
        real,    intent(in)  :: bfac, smpd
        real,    intent(out) :: kweights(0:box)
        real    :: bfac_sc, res
        integer :: k, nyq
        nyq     = fdim(box)
        bfac_sc = bfac / 4.
        kweights(0) = 1.
        do k=1,box
            res         = real(k) / (real(box) * smpd)
            kweights(k) = max(0., exp(-bfac_sc * res * res))
        end do
        ! normalize
        kweights(0:box) = kweights(0:box) / sum(kweights(0:nyq))
    end subroutine calc_norm_bfac_weights

    ! LINEAR ALGEBRA STUFF

    !>   subroutine to find the inverse of a square matrix
    !!         author : louisda16th a.k.a ashwith j. rego
    !!         reference : algorithm explained at:
    !!         http://math.uww.edu/~mcfarlat/inverse.htm
    !!         http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
    subroutine matinv(matrix, inverse, n, errflg)
        integer, intent(in)  :: n       !< size of square matrix
        integer, intent(out) :: errflg  !< return error status. -1 for error, 0 for normal
        real, intent(in), dimension(n,n)  :: matrix  !< input matrix
        real, intent(out), dimension(n,n) :: inverse !< inverted matrix
        logical :: flag = .true.
        integer :: i, j, k
        real :: m
        real, dimension(n,2*n) :: augmatrix !< augmented matrix
        ! augment input matrix with an identity matrix
        do i=1,n
            do j=1,2*n
                if(j <= n )then
                    augmatrix(i,j) = matrix(i,j)
                else if((i+n) == j)then
                    augmatrix(i,j) = 1.
                else
                    augmatrix(i,j) = 0.
                endif
            end do
        end do
        ! reduce augmented matrix to upper traingular form
        do k=1,n-1
              if( is_zero(augmatrix(k,k)) )then
                 flag = .false.
                do i=k+1,n
                    if( augmatrix(i,k) > 0. )then
                        do j=1,2*n
                            augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
                        end do
                        flag = .true.
                        exit
                    endif
                    if(flag .eqv. .false.)then
                        inverse = 0.
                        errflg = -1
                        return
                    endif
                end do
            endif
            do j=k+1, n
                m = augmatrix(j,k)/augmatrix(k,k)
                do i=k,2*n
                    augmatrix(j,i) = augmatrix(j,i)-m*augmatrix(k,i)
                end do
            end do
        end do
        ! test for invertibility
        do i=1,n
            if( is_zero(augmatrix(i,i)) )then
                inverse = 0
                errflg = -1
                return
            endif
        end do
        ! make diagonal elements as 1
        do i=1,n
            m = augmatrix(i,i)
            do j=i,2*n
                augmatrix(i,j) = augmatrix(i,j)/m
            end do
        end do
        ! reduced right side half of augmented matrix to identity matrix
        do k=n-1,1,-1
            do i=1,k
                m = augmatrix(i,k+1)
                do j = k,2*n
                    augmatrix(i,j) = augmatrix(i,j)-augmatrix(k+1,j)*m
                end do
            end do
        end do
        ! store answer
        do i=1,n
            do j=1,n
                inverse(i,j) = augmatrix(i,j+n)
            end do
        end do
        errflg = 0
    end subroutine matinv

    !>   least squares straight line fit, from Sjors, he took it from
    !!         http://mathworld.wolfram.com/LeastSquaresFitting.html
    !!         ss_xx = sum_i x_i^2 - n * ave_x^2
    !!         ss_yy = sum_i y_i^2 - n * ave_y^2
    !!         ss_xy = sum_i x_i * y_i - n * ave_x * n ave_y
    !! \param  slope  Linear gradient: slope = xx_xy / ss_xx
    !!  \param intercept y-intercept:  intercept = ave_y - slope * ave_x
    !!  \param corr Correlation coefficient: corr_coeff = ss_xy^2 / (ss_xx * ss_yy)
    subroutine fit_straight_line( n, datavec, slope, intercept, corr )
        integer, intent(in) :: n                                            !< size of vector
        real, intent(in)    :: datavec(n,2)                                 !< input vector
        real, intent(out)   :: slope, intercept, corr
        double precision    :: ave_x, ave_y, ss_xx, ss_yy, ss_xy, x, y, dn
        integer             :: i
        ave_x = 0.d0
        ave_y = 0.d0
        ss_xx = 0.d0
        ss_yy = 0.d0
        ss_xy = 0.d0
        do i=1,n
            x = datavec(i,1)
            y = datavec(i,2)
            ave_x = ave_x+x
            ave_y = ave_y+y
            ss_xx = ss_xx+x*x
            ss_yy = ss_yy+y*y
            ss_xy = ss_xy+x*y
        end do
        dn = dble(n)
        ave_x = ave_x/dn
        ave_y = ave_y/dn
        ss_xx = ss_xx-dn*ave_x*ave_x
        ss_yy = ss_yy-dn*ave_y*ave_y
        ss_xy = ss_xy-dn*ave_x*ave_y
        slope     = real(ss_xy/ss_xx)
        intercept = real(ave_y-slope*ave_x)
        corr      = real((ss_xy*ss_xy)/(ss_xx*ss_yy))
    end subroutine fit_straight_line

    !>    is for calculating the radius
    pure function hyp( x1, x2, x3 ) result( h )
        real, intent(in) :: x1, x2
        real, intent(in), optional :: x3
        real :: h
        if( present(x3) )then
            h = sqrt(x1**2.+x2**2.+x3**2.)
        else
            h = sqrt(x1**2.+x2**2.)
        endif
    end function hyp

    !>   calculates the euclidean distance between two vectors of dimension _n_
    pure function euclid( vec1, vec2 ) result( dist )
        real, intent(in)    :: vec1(:), vec2(:)
        real                :: dist
        dist = sqrt(sum((vec1-vec2)**2.))
    end function euclid

    !>   calculates the argument of a vector
    pure function arg( vec ) result( length )
        real, intent(in) :: vec(:)
        real :: length
        length = sqrt(sum(vec*vec))
    end function arg

    !>   projects a 3d vector in the _z_-direction
    subroutine projz( vec3, vec2 )
        real, intent(in)  :: vec3(3)
        real, intent(out) :: vec2(2)
        vec2(1) = vec3(1)
        vec2(2) = vec3(2)
    end subroutine projz

    !>   is for 2d rotation matrix generation
    pure subroutine rotmat2d( ang, mat )
        real, intent(in)  :: ang ! in degrees
        real, intent(out) :: mat(2,2)
        real :: ang_in_rad
        ang_in_rad = ang*pi/180.
        mat(1,1) = cos(ang_in_rad)
        mat(1,2) = sin(ang_in_rad)
        mat(2,1) = -mat(1,2)
        mat(2,2) = mat(1,1)
    end subroutine rotmat2d

    !>  in-plane parameters to 3x3 transformation matrix
    function make_transfmat( psi, tx, ty )result( R )
        real,intent(in) :: psi,tx,ty
        real            :: R(3,3),radpsi,cospsi,sinpsi
        radpsi = deg2rad( psi )
        cospsi = cos( radpsi )
        sinpsi = sin( radpsi )
        R      = 0.
        R(1,1) = cospsi
        R(2,2) = cospsi
        R(1,2) = sinpsi
        R(2,1) = -sinpsi
        R(1,3) = tx
        R(2,3) = ty
        R(3,3) = 1.
    end function make_transfmat

    !>    quadratic interpolation in 2D, from spider
    function quadri(xx, yy, fdata, nx, ny) result(q)
        integer, intent(in) :: nx, ny
        real, intent(in)    :: xx, yy, fdata(nx,ny)
        real    :: x,y,dx0,dy0,f0,c1,c2,c3,c4,c5,dxb,dyb,q
        integer :: i,j,ip1,im1,jp1,jm1,hxc,hyc,ic,jc
        x = xx
        y = yy
        ! circular closure
        if(x.lt.1.0) x = x+(1-ifix(x)/nx)*nx
        if(x.gt.float(nx)+0.5) x = amod(x-1.0,float(nx))+1.0
        if(y.lt.1.0) y = y+(1-ifix(y)/ny)*ny
        if(y.gt.float(ny)+0.5) y = amod(y-1.0,float(ny))+1.0
        i = ifix(x)
        j = ifix(y)
        dx0 = x-i
        dy0 = y-j
        ip1 = i+1
        im1 = i-1
        jp1 = j+1
        jm1 = j-1
        if(ip1 .gt. nx) ip1 = ip1-nx
        if(im1 .lt. 1) im1 = im1+nx
        if(jp1 .gt. ny) jp1 = jp1-ny
        if(jm1 .lt. 1) jm1 = jm1+ny
        f0 = fdata(i,j)
        c1 = fdata(ip1,j)-f0
        c2 = (c1-f0+fdata(im1,j))*0.5
        c3 = fdata(i,jp1)-f0
        c4 = (c3-f0+fdata(i,jm1))*0.5
        dxb = (dx0-1)
        dyb = (dy0-1)
        ! hxc & hyc are either 1 or -1
        hxc = int(sign(1.0,dx0))
        hyc = int(sign(1.0,dy0))
        ic = i+hxc
        jc = j+hyc
        if(ic .gt. nx)then
           ic = ic-nx
        elseif (ic .lt. 1)then
           ic = ic+nx
        endif
        if(jc .gt. ny)then
           jc = jc-ny
        elseif (jc .lt. 1)then
           jc = jc+ny
        endif
        c5 = ((fdata(ic,jc)-f0-hxc*c1-(hxc*(hxc-1.0))&
        *c2-hyc*c3-(hyc*(hyc-1.0))*c4)*(hxc*hyc))
        q = f0+dx0*(c1+dxb*c2+dy0*c5)+dy0*(c3+dyb*c4)
    end function quadri

    ! SORTING

    function peakfinder_2( vals ) result( peakpos )
        real,    intent(in)  :: vals(:)
        logical, allocatable :: peakpos(:)
        integer :: n, i
        n = size(vals)
        allocate(peakpos(n), source=.false.)
        if( vals(1) > vals(2) )  peakpos(1) = .true.
        do i=2,n-1
            if( vals(i) >= vals(i-1) .and. vals(i) >= vals(i+1) ) peakpos(i) = .true.
        end do
        if( vals(n) > vals(n-1) ) peakpos(n) = .true.
    end function peakfinder_2

    subroutine peakfinder_inplace( vals, peakpos )
        real,    intent(in)  :: vals(:)
        logical, intent(out) :: peakpos(:)
        integer :: n, i
        n = size(vals)
        peakpos = .false.
        if( vals(1) > vals(2) )  peakpos(1) = .true.
        do i=2,n-1
            if( vals(i) >= vals(i-1) .and. vals(i) >= vals(i+1) ) peakpos(i) = .true.
        end do
        if( vals(n) > vals(n-1) ) peakpos(n) = .true.
    end subroutine peakfinder_inplace

    !>   rheapsort from numerical recipes (largest last)
    subroutine hpsort_1( rarr, iarr )
        real,    intent(inout) :: rarr(:)
        integer, intent(inout) :: iarr(:)
        integer :: i, ir, j, l, ia, n
        real    :: ra
        n = size(rarr)
        if( n /= size(iarr) )&
        &call simple_exception('nonconforming array sizes; hpsort_1', __FILENAME__ , __LINE__)
        if( n < 2 ) return
        l  = n/2+1
        ir = n
        do
            if(l > 1)then
                l  = l-1
                ra = rarr(l)
                ia = iarr(l)
            else
                ra = rarr(ir)
                ia = iarr(ir)
                rarr(ir) = rarr(1)
                iarr(ir) = iarr(1)
                ir = ir-1
                if(ir == 1)then
                    rarr(1) = ra
                    iarr(1) = ia
                    return
                endif
            endif
            i = l
            j = l+l
            do while(j <= ir)
                if(j < ir) then
                    if(rarr(j) < rarr(j+1)) j = j+1
                endif
                if(ra < rarr(j))then
                    rarr(i) = rarr(j)
                    iarr(i) = iarr(j)
                    i = j
                    j = j+j
                else
                    j = ir+1
                endif
                rarr(i) = ra
                iarr(i) = ia
            end do
        end do
    end subroutine hpsort_1

    !>   rheapsort from numerical recepies (largest last)
    subroutine hpsort_2( iarr )
        integer, intent(inout) :: iarr(:)
        integer :: i, ir, j, l, ra, n
        n = size(iarr)
        if( n < 2) return
        l  = n/2+1
        ir = n
        do
            if(l > 1)then
                l  = l-1
                ra = iarr(l)
            else
                ra = iarr(ir)
                iarr(ir) = iarr(1)
                ir = ir-1
                if(ir == 1)then
                    iarr(1) = ra
                    return
                endif
            endif
            i = l
            j = l+l
            do while(j <= ir)
                if(j < ir) then
                    if(iarr(j) < iarr(j+1)) j = j+1
                endif
                if(ra < iarr(j))then
                    iarr(i) = iarr(j)
                    i = j
                    j = j+j
                else
                    j = ir+1
                endif
                iarr(i) = ra
            end do
        end do
    end subroutine hpsort_2

    !>   rheapsort from numerical recepies (largest last)
    subroutine hpsort_3( iarr, p1_lt_p2 )
        integer, intent(inout) :: iarr(:)
        interface
            function p1_lt_p2( p1, p2 ) result( val )
                integer, intent(in) :: p1, p2
                logical :: val
            end function p1_lt_p2
        end interface
        integer :: i, ir, j, l, ra, n
        n = size(iarr)
        if( n < 2) return
        l  = n/2+1
        ir = n
        do
            if(l > 1)then
                l  = l-1
                ra = iarr(l)
            else
                ra = iarr(ir)
                iarr(ir) = iarr(1)
                ir = ir-1
                if(ir == 1)then
                    iarr(1) = ra
                    return
                endif
            endif
            i = l
            j = l+l
            do while(j <= ir)
                if(j < ir) then
                    if(p1_lt_p2(iarr(j),iarr(j+1))) j = j+1
                endif
                if(p1_lt_p2(ra,iarr(j)))then
                    iarr(i) = iarr(j)
                    i = j
                    j = j+j
                else
                    j = ir+1
                endif
                iarr(i) = ra
            end do
        end do
    end subroutine hpsort_3

    !>   rheapsort from numerical recepies (largest last)
    subroutine hpsort_4( rarr )
        real, intent(inout) :: rarr(:)
        integer :: i, ir, j, l, n
        real    :: ra
        n = size(rarr)
        if( n < 2) return
        l  = n/2+1
        ir = n
        do
            if(l > 1)then
                l  = l-1
                ra = rarr(l)
            else
                ra = rarr(ir)
                rarr(ir) = rarr(1)
                ir = ir-1
                if(ir == 1)then
                    rarr(1) = ra
                    return
                endif
            endif
            i = l
            j = l+l
            do while(j <= ir)
                if(j < ir) then
                    if(rarr(j) < rarr(j+1)) j = j+1
                endif
                if(ra < rarr(j))then
                    rarr(i) = rarr(j)
                    i = j
                    j = j+j
                else
                    j = ir+1
                endif
                rarr(i) = ra
            end do
        end do
    end subroutine hpsort_4

    !>   rheapsort from numerical recepies (largest last)
    subroutine hpsort_5( rarr, rarr2 )
        real, intent(inout) :: rarr(:), rarr2(:)
        integer :: i, ir, j, l, n
        real    :: ra, ra2
        n = size(rarr)
        if( n /= size(rarr2) )&
        &call simple_exception('nonconforming array sizes; hpsort_5', __FILENAME__ , __LINE__)
        if( n < 2) return
        l  = n/2+1
        ir = n
        do
            if(l > 1)then
                l  = l-1
                ra = rarr(l)
                ra2 = rarr2(l)
            else
                ra = rarr(ir)
                ra2 = rarr2(ir)
                rarr(ir) = rarr(1)
                rarr2(ir) = rarr2(1)
                ir = ir-1
                if(ir == 1)then
                    rarr(1) = ra
                    rarr2(1) = ra2
                    return
                endif
            endif
            i = l
            j = l+l
            do while(j <= ir)
                if(j < ir) then
                    if(rarr(j) < rarr(j+1)) j = j+1
                endif
                if(ra < rarr(j))then
                    rarr(i) = rarr(j)
                    rarr2(i) = rarr2(j)
                    i = j
                    j = j+j
                else
                    j = ir+1
                endif
                rarr(i) = ra
                rarr2(i) = ra2
            end do
        end do
    end subroutine hpsort_5

    !> Vector angle between two non-zero vectors in R^n
    !! magnitude of v and w must equal 1
    !! \theta = \arccos \frac {\mathbf v \cdot \mathbf w}{\left\Vert{\mathbf v}\right\Vert \left\Vert{\mathbf w}\right\Vert}
    !! by definition of the cosine formula for dot product.
    !! This function assumes the input vectors are normals!!
    pure function vector_angle_norm(v, w) result(eulerdist)
        real, dimension(3), intent(in) :: v,w
        real :: eulerdist
        eulerdist = acos( dot_product(v, w) )
    end function vector_angle_norm

    !> sort and return the 3 lowest
    function min3( rarr ) result(min_3)
        real, intent(in) :: rarr(:)
        real :: min_3(3)
        integer :: j, n
        real    :: ra
        n = size(rarr)
        if( n < 4)then
            min_3(:n) = rarr
            return
        end if
        min_3 = rarr(:3)
        call hpsort(min_3)
        do j=4,n
            if(rarr(j) < min_3(3))then
                ra = rarr(j)
                if(ra < min_3(2))then
                    if(ra < min_3(1))then
                        min_3 = (/ ra,  min_3(1), min_3(2) /)
                    else
                        min_3 = (/ min_3(1), ra, min_3(2) /)
                    end if
                else
                    min_3 = (/ min_3(1), min_3(2), ra /)
                end if
            end if
        end do
    end function min3

    function maxnloc( rarr, n ) result( loc )
        real,    intent(in) :: rarr(:)
        integer, intent(in) :: n
        real    :: arr(n), val
        integer :: loc(n), i, j, sz
        logical :: val_lt(n)
        sz = size(rarr)
        if( sz < n )&
        &call simple_exception('cannot identify more maxima than elements in the array; maxnloc_2', __FILENAME__ , __LINE__)
        loc = (/(i,i=1,n)/)
        arr = rarr(:n)
        call hpsort(arr, loc)
        call reverse(loc)
        if( sz == n ) return
        do i=n+1,sz
            val    = rarr(i)
            val_lt = val > rarr(loc)
            if( any(val_lt) )then
                do j=1,n
                    if( val_lt(j) )then
                        if( j == 1 )then
                            loc = [i,loc(1:n-1)]
                            exit
                        else if( j == n )then
                            loc(n) = i
                            exit
                        else
                            ! 1:j-1 untouched
                            ! insert i @ location j
                            ! j+1:n-1 forms the second part of the index array
                            loc = [loc(1:j-1),i,loc(j:n-1)]
                            exit
                        endif
                    endif
                end do
            endif
        enddo
    end function maxnloc

    ! SEARCHING

    !>   given an array arr(1:n), and given value x, locate returns a value j such that x is
    !!         between arr(j) and arr(j+1). arr(1:n) must be monotonic, either increasing or decreasing.
    !!         j=0 or j=n is returned to indicate that x is out of range, from numerical recepies
    pure function locate_1( arr, n, x ) result( j )
        integer, intent(in) :: n
        real, intent(in)    :: arr(n), x
        integer             :: jl, jm, ju, j
        jl = 0                ! initialize lower
        ju = n+1              ! and upper limits
        do while( ju-jl > 1 ) ! if we are not yet done
            jm = (ju+jl)/2    ! compute a midpoint
            if((arr(n) >= arr(1)).eqv.(x >= arr(jm))) then
                jl = jm       ! and replace either the lower limit
            else
                ju = jm       ! or the upper limit, as appropriate
            endif
        end do
        if( abs(x-arr(1)) < TINY )then
            j = 1
        else if( abs(x-arr(n)) < TINY )then
            j = n-1
        else
            j = jl
        endif
    end function locate_1

    !>   given an array arr(1:n), and given value x, locate returns a value j such that x is
    ! between arr(j) and arr(j+1). arr(1:n) must be monotonic, either increasing or decreasing.
    ! j=0 or j=n is returned to indicate that x is out of range, from numerical recepies
    pure function locate_2( arr, n, x ) result( j )
        integer, intent(in) :: n             !< size of list
        integer, intent(in) :: arr(n), x     !< list and search value
        integer             :: jl, jm, ju, j
        jl = 0                ! initialize lower
        ju = n+1              ! and upper limits
        do while( ju-jl > 1 ) ! if we are not yet done
            jm = (ju+jl)/2    ! compute a midpoint
            if((arr(n) >= arr(1)).eqv.(x >= arr(jm))) then
                jl = jm       ! and replace either the lower limit
            else
                ju = jm       ! or the upper limit, as appropriate
            endif
        end do
        if( x == arr(1) )then
            j = 1
        else if( x == arr(n) )then
            j = n-1
        else
            j = jl
        endif
    end function locate_2

    !>   for finding closest element in an ordered list
    subroutine find_1( arr, n, x, j, dist1 )
        integer, intent(in)  :: n           !< size of list
        real,    intent(in)  :: arr(n), x   !< list and search value
        real,    intent(out) :: dist1
        integer, intent(out) :: j
        real                 :: dist2
        j = max(1,locate_1( arr, n, x ))
        dist1 = arr(j)-x
        if( j < n ) then
            dist2 = arr(j+1)-x
            if( abs(dist1) >= abs(dist2) )then
                j = j+1
                dist1 = dist2
            endif
        endif
    end subroutine find_1

    !>   for finding closest element in an ordered list
    subroutine find_2( arr, n, x, j, dist1 )
        integer, intent(in)  :: n         !< size of list
        integer, intent(in)  :: arr(n), x !< list and search value
        integer, intent(out) :: dist1
        integer, intent(out) :: j
        integer              :: dist2
        j = max(1,locate_2( arr, n, x ))
        dist1 = arr(j)-x
        if( j < n ) then
            dist2 = arr(j+1)-x
            if( abs(dist1) >= abs(dist2) )then
                j = j+1
                dist1 = dist2
            endif
        endif
    end subroutine find_2

    !>   for calculating the median
    function median( arr ) result( val )
        real, intent(in)  :: arr(:)
        real              :: copy(size(arr))
        real    :: val, val1, val2
        integer :: n, pos1, pos2
        n = size(arr)
        if( is_even(n) )then
            pos1 = n/2
            pos2 = pos1+1
        else
            pos1 = nint(real(n)/2.)
            pos2 = pos1
        endif
        copy = arr
        if( pos1 == pos2 )then
            val  = selec(pos1,n,copy)
        else
            val1 = selec(pos1,n,copy)
            val2 = selec(pos2,n,copy)
            val  = (val1+val2)/2.
        endif
    end function median

    !>   for calculating the median
    function median_nocopy( arr ) result( val )
        real, intent(inout) :: arr(:)
        real    :: val, val1, val2
        integer :: n, pos1, pos2
        n = size(arr)
        if( is_even(n) )then
            pos1 = n/2
            pos2 = pos1+1
        else
            pos1 = nint(real(n)/2.)
            pos2 = pos1
        endif
        if( pos1 == pos2 )then
            val  = selec(pos1,n,arr)
        else
            val1 = selec(pos1,n,arr)
            val2 = selec(pos2,n,arr)
            val  = (val1+val2)/2.
        endif
    end function median_nocopy

    function stdev (a)
        real, intent(in) :: a(:)
        real    :: stdev, avg, SumSQR, var, n
        n= real(size(a))
        if(n < 2) return
        avg = sum(a)/n
        SumSQR = sum(sqrt(a))/n
        var = (SumSQR - avg*avg)/(n-1) !CHIARA
        ! var = (SumSQR - avg*avg/n)/(n-1)
        stdev = sqrt(var)
    end function stdev

    ! from Numerical Recipes in Fortran 77.
    real function quickselect(arr, k)
        real,    intent(inout) :: arr(:)
        integer, intent(in)    :: k
        integer :: i,ir,j,l,mid,n
        real    :: a,temp
        n  = size(arr)
        l  = 1
        ir = n
        do while (ir-l.gt.1)
            mid = (l+ir)/2
            temp = arr(mid)
            arr(mid) = arr(l+1)
            arr(l+1) = temp
            if (arr(l).gt.arr(ir)) then
                temp = arr(l)
                arr(l) = arr(ir)
                arr(ir) = temp
            endif
            if (arr(l+1).gt.arr(ir)) then
                temp = arr(l+1)
                arr(l+1) = arr(ir)
                arr(ir) = temp
            endif
            if (arr(l).gt.arr(l+1)) then
                temp = arr(l)
                arr(l) = arr(l+1)
                arr(l+1) = temp
            endif
            i = l+1
            j = ir
            a = arr(l+1)
            do
                i = i+1
                if (arr(i).lt.a) cycle
                j = j-1
                do while (arr(j).gt.a)
                    j=j-1
                end do
                if (j.lt.i) exit
                temp = arr(i)
                arr(i) = arr(j)
                arr(j) = temp
            end do
            arr(l+1) = arr(j)
            arr(j) = a
            if (j.ge.k) ir = j-1
            if (j.le.k) l = i
        end do
        if (ir-1.eq.1) then
            if (arr(ir).lt.arr(l)) then
                temp = arr(l)
                arr(l) = arr(ir)
                arr(ir) = temp
            endif
        endif
        quickselect = arr(k)
    end function quickselect

    !>   for selecting kth largest, array is modified
    real function selec(k,n,arr)
        integer, intent(in)    :: k,n
        real,    intent(inout) :: arr(n)
        integer :: i,ir,j,l,mid
        real    :: a,temp
        l = 1
        ir = n
    22  if (ir-l.le.1) then
            if (ir-1.eq.1) then
                if (arr(ir).lt.arr(l)) then
                    temp = arr(l)
                    arr(l) = arr(ir)
                    arr(ir) = temp
                endif
            endif
            selec = arr(k)
            return
        else
            mid = (l+ir)/2
            temp = arr(mid)
            arr(mid) = arr(l+1)
            arr(l+1) = temp
            if (arr(l).gt.arr(ir)) then
                temp = arr(l)
                arr(l) = arr(ir)
                arr(ir) = temp
            endif
            if (arr(l+1).gt.arr(ir)) then
                temp = arr(l+1)
                arr(l+1) = arr(ir)
                arr(ir) = temp
            endif
            if (arr(l).gt.arr(l+1)) then
                temp = arr(l)
                arr(l) = arr(l+1)
                arr(l+1) = temp
            endif
            i = l+1
            j = ir
            a = arr(l+1)
    23      continue
            i = i+1
            if (arr(i).lt.a) goto 23
    24      continue
            j = j-1
            if (arr(j).gt.a) goto 24
            if (j.lt.i)      goto 25
            temp = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
            goto 23
    25      arr(l+1) = arr(j)
            arr(j) = a
            if (j.ge.k) ir = j-1
            if (j.le.k) l = i
        endif
        goto 22
    end function selec

    pure function norm_2_sp(v) result(r)
        real, intent(in) :: v(:)
        real             :: r
        r = sqrt(dot_product(v,v))
    end function norm_2_sp

    pure function norm_2_dp(v) result(r)
        real(dp), intent(in) :: v(:)
        real(dp)             :: r
        r = sqrt(dot_product(v,v))
    end function norm_2_dp

    ! GCD computes greatest common divisors
    pure function gcd (p,q)
        integer :: gcd
        integer, intent(in) :: p,q
        integer x,y
        x = abs ( p )
        y = abs ( q )
        do
            if ( x .gt. y ) then
                x = mod ( x , y )
                if ( x .eq. 0 ) exit
            else if ( x .lt. y ) then
                y = mod ( y , x )
                if ( y .eq. 0 ) exit
            end if
        end do
        gcd = max ( x , y )
    end function gcd

    ! LCM computes the least common multiple of integers P,Q
    pure integer function lcm( i, j )
        integer, intent( in ) :: i, j
        lcm = ( i * j ) / gcd( i, j )
    end function lcm

    ! TRACE calculates the trace of a real 2D matrix
    pure function trace(mat) result (tr)
        real, intent(in) :: mat(:,:)
        real             :: tr
        integer          :: i
        tr = 0.
        do i = 1, size(mat, dim = 1)
            tr = tr + mat(i,i)
        enddo
    end function trace

    !>  \brief  generates a binary mask from a logical one
    function logical2bin( mask ) result( matrix )
        logical, intent(in) :: mask(:,:,:)
        real, allocatable   :: matrix(:,:,:)
        integer      :: s(3), i, j
        s = shape(mask)
        if(allocated(matrix)) deallocate(matrix)
        allocate(matrix(s(1),s(2),1), source = 0.)
        do i = 1, s(1)
            do j = 1, s(2)
                if(mask(i,j,1)) matrix(i,j,1)=1.
            enddo
        enddo
    end function logical2bin

    ! This function takes in input arraya and gives as output arrayb
    ! which is the same as array, but with NO repetition in the elements.
    subroutine elim_dup(arraya, arrayb)
        real,              intent(in)  :: arraya(:)
        real, allocatable, intent(out) :: arrayb(:)
        integer              :: ix
        logical, allocatable :: mask(:)
        allocate(mask(size(arraya)))
        mask = .true.
        do ix = size(arraya),2,-1
            mask(ix) = .not.(any(arraya(:ix-1)==arraya(ix)))
        enddo
        arrayb = pack(arraya, mask)
        deallocate(mask)
    end subroutine elim_dup

    !!!!!!ADDED BY CHIARA!!!!!
    ! All the 4 following routines have testing purposes and allow
    ! the user to print a matrix so that the output is similar to
    ! the usual representation of a matrix.
    subroutine vis_2Dreal_mat(mat)
        real, intent(in) :: mat(:,:)
        integer :: j, s(2)
        s = shape(mat)
        do j = 1, s(1)
            write(logfhandle,*) mat(j,:)
        enddo
    end subroutine vis_2Dreal_mat

    subroutine vis_2Dinteger_mat(mat)
        integer, intent(in) :: mat(:,:)
        integer :: j, s(2)
        s = shape(mat)
        do j = 1, s(1)
            write(logfhandle,*) mat(j,:)
        enddo
    end subroutine vis_2Dinteger_mat

    subroutine vis_3Dreal_mat(mat)
        real, intent(in) :: mat(:,:,:)
        integer :: j, s(3)
        s = shape(mat)
        do j = 1, s(1)
            write(logfhandle,*) mat(j,:,1)
        enddo
    end subroutine vis_3Dreal_mat

    subroutine vis_3Dinteger_mat(mat)
        integer, intent(in) :: mat(:,:,:)
        integer :: j, s(3)
        s = shape(mat)
        do j = 1, s(1)
            write(logfhandle,*) mat(j,:,1)
        enddo
    end subroutine vis_3Dinteger_mat

    ! calculates the mode of an array
    subroutine mode_1( arr, m, npxls_at_mode )
        integer,           intent(in)  :: arr(:)
        integer,           intent(out) :: m ! mode
        integer, optional, intent(out) :: npxls_at_mode
        integer, allocatable :: counts(:)
        integer :: i, loc(1)
        ! initialise array to count occurrences of each value
        allocate(counts(minval(arr):maxval(arr)), source=0)
        ! loop over inputted array, counting occurrence of each value
        do i=1,size(arr)
            counts(arr(i)) = counts(arr(i)) + 1
        end do
        ! Find the mode
        loc = maxloc(counts)
        m   = loc(1)
        if(present(npxls_at_mode)) npxls_at_mode = maxval(counts)
    end subroutine mode_1

    subroutine mode_2(arr, m, npxls_at_mode)  !REAL VECTORS
        real, intent(in)  :: arr(:)
        real, intent(out) :: m(1)
        integer, optional, intent(out) :: npxls_at_mode
        integer              :: n !the number of intervals
        real,    allocatable :: xhist(:)
        integer, allocatable :: yhist(:)
        integer              :: i, j
        real                 :: xmin, xmax, dx
        integer, dimension(:), allocatable :: counts
        integer :: astat
        character(len=128) :: error_str
        xmin=minval(arr)
        xmax=maxval(arr)
        n = 2*(nint(xmax-xmin)+1) !this value influence the result because it determins how to approximate the steps
        allocate(xhist(0:n-1))
        allocate(yhist(n), source = 0 )
        dx=(xmax+1-xmin)/n
        do i=0,n-1
            xhist(i)=xmin+i*dx
        end do
        do i=1,size(arr, dim =1)
            j = nint((arr(i)-xmin)/dx)+1
            yhist(j)=yhist(j)+1
        end do
        m(:) = xhist(maxloc(yhist)-1)
        if(present(npxls_at_mode)) npxls_at_mode = maxval(yhist)
    end subroutine mode_2

    ! takes in input a vector x and an integer n
    ! stores in xhist the discretisation of the values in x
    ! returns the number of occurences in yhist
    subroutine create_hist_vector(x,n,xhist,yhist)
        real,                 intent(in)  :: x(:)     !data
        integer,              intent(in)  :: n        !number of intervals
        real,    allocatable, intent(out) :: xhist(:) !discretization of the values
        integer, allocatable, intent(out) :: yhist(:) !number of occurences
        integer   :: i, j
        real      :: xmin, xmax, dx
        if(allocated(xhist)) deallocate(xhist)
        if(allocated(yhist)) deallocate(yhist)
        allocate(xhist(0:n-1), source = 0.)
        allocate(yhist(n), source = 0 )
        xmin=minval(x)
        xmax=maxval(x)
        dx=(xmax+1-xmin)/n
        do i=0,n-1
            xhist(i)=xmin+i*dx
        end do
        do i=1,size(x)
            j=nint((x(i)-xmin)/dx)+1
            yhist(j)=yhist(j)+1
        end do
    end subroutine create_hist_vector

    !!!!!!!!!!!!!!!!!!!ADDED BY CHIARAAAA!!!!!!!!!!!!!!
    !This function is meant to be a support for performing histogram stretching.
    !It takes in input the histogram vectors xhist and yhist (see function create_hist_vector
    !in the simple_image file)
    !and gives as an output the mode m of the histogram, the corresponding number of pixels
    !at mode npxls_at_mode and calculates the min and max gray-level values as limits
    !for histogram stretching (stretch_lim). By default, the limits specify the bottom
    !1% and the top 1% of all pixel values.
    subroutine find_stretch_minmax(xhist,yhist,m,npxls_at_mode,stretch_lim)
        real,    intent(in)  :: xhist(:)
        integer, intent(in)  :: yhist(:)
        real,    intent(out) :: m(1)
        integer, intent(out) :: npxls_at_mode
        real,    intent(out) :: stretch_lim(2)
        real, parameter :: OFFSET = 0. !to change??
        integer         ::  i, ilim(1)
        m(:)    = xhist(maxloc(yhist)+1) !+1 cuz of different index
        ilim(:) =       maxloc(yhist)
        npxls_at_mode = maxval(yhist)
        do i = 1, ilim(1)-1
            if(yhist(i) > real(npxls_at_mode/100) - OFFSET) goto 128  !bottom 1%
        enddo
        128 continue
        stretch_lim(1) = xhist(i-1)
        do i = ilim(1)+1, size(yhist, dim = 1)
            if(yhist(i) < real(npxls_at_mode/100) + OFFSET) goto 132  !top 1%
        enddo
        132 continue
        stretch_lim(2) = xhist(i)
    end subroutine find_stretch_minmax

    ! rescales the vector to a new input range
    subroutine scale_vector( x, new_range )
        real, intent(inout) :: x(:)
        real, intent(in)    :: new_range(2)
        real :: old_range(2), sc
        old_range(1) = minval(x)
        old_range(2) = maxval(x)
        sc = (new_range(2) - new_range(1))/(old_range(2) - old_range(1))
        x(:) = sc*x(:)+new_range(1)-sc*old_range(1)
    end subroutine scale_vector


    ! Otsu's method
    ! the algorithm assumes that x contains two classes of numbers following bi-modal histogram (foreground and background)
    ! it then calculates the optimum threshold separating the two classes so that their combined spread (intra-class variance) is minimal
    subroutine otsu_1( x, thresh )
        real, intent(inout) :: x(:)
        real, intent(out)   :: thresh ! threshold for binarisation
        integer, parameter  :: MIN_VAL = 0, MAX_VAL = 255
        integer :: i, T
        real    :: q1, q2, m1, m2, sigma1, sigma2, sigma, sigma_next, sum1, sum2
        real,    allocatable :: p(:), xhist(:)
        integer, allocatable :: yhist(:)
        call scale_vector(x,real([MIN_VAL,MAX_VAL]))
        allocate(p(MIN_VAL:MAX_VAL), source=0.)
        call create_hist_vector(x, MAX_VAL-MIN_VAL+1, xhist, yhist)
        do i = MIN_VAL, MAX_VAL
            p(i) = yhist(i+1)
        enddo
        deallocate(xhist,yhist)
        p    = p / size(x, dim = 1) ! normalise, it's a probability
        q1   = 0.
        q2   = sum(p)
        sum1 = 0.
        sum2 = 0.
        do i = MIN_VAL, MAX_VAL
            sum2 = sum2 + i * p(i)
        enddo
        sigma = HUGE(sigma) ! initialisation, to get into the loop
        do T = MIN_VAL, MAX_VAL - 1
          q1   = q1 + p(T)
          q2   = q2 - p(T)
          sum1 = sum1 + T*p(T)
          sum2 = sum2 - T*p(T)
          m1   = sum1/q1
          m2   = sum2/q2
          sigma1 = 0.
          if( T > MIN_VAL )then ! do not know if it is necessary
              do i = MIN_VAL, T
                  sigma1 = sigma1 + (i - m1)**2 * p(i)
              enddo
          else
              sigma1 = 0.
          endif
          sigma1 = sigma1 / q1
          sigma2 = 0.
          if( T < MAX_VAL - 1 )then
              do i = T + 1, MAX_VAL
                  sigma2 = sigma2 + (i - m2)**2 * p(i)
              enddo
          else
              sigma2 = 0.
          endif
          sigma2 = sigma2 / q2
          sigma_next = q1 * sigma1 + q2 * sigma2
          if( sigma_next < sigma .and. T > MIN_VAL )then
              thresh = T ! keep the minimum
              sigma  = sigma_next
          elseif( T == MIN_VAL )then
              sigma = sigma_next
          endif
        enddo
    end subroutine otsu_1

    ! Otsu's method, see above
    subroutine otsu_2(x, x_out, thresh)
        real,              intent(inout) :: x(:)
        real, allocatable, intent(out)   :: x_out(:)
        real, optional,    intent(out)   :: thresh  !selected threshold for binarisation
        real :: threshold ! threshold for binarisation
        call otsu_1(x, threshold)
        if( allocated(x_out) ) deallocate(x_out)
        allocate(x_out(size(x)), source=0.)
        where(x > threshold)
            x_out = 1.
        elsewhere
            x_out = 0.
        endwhere
        if(present(thresh)) thresh = threshold
    end subroutine otsu_2

end module simple_math
