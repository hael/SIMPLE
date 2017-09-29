! various mathematical subroutines and functions

module simple_math
use simple_defs
use simple_syslib, only: alloc_errchk 
implicit none

interface is_a_number
    module procedure is_a_number_1
    module procedure is_a_number_2
end interface

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

interface is_even
    module procedure is_even_1
    module procedure is_even_2
end interface

interface cosedge
    module procedure cosedge_1
    module procedure cosedge_2
    module procedure cosedge_3
end interface

interface cosedge_inner
    module procedure cosedge_inner_1
    module procedure cosedge_inner_2
end interface

interface hardedge
    module procedure hardedge_1
    module procedure hardedge_2
end interface

interface hardedge_inner
    module procedure hardedge_inner_1
    module procedure hardedge_inner_2
end interface

interface find
    module procedure find_1
    module procedure find_2
end interface

interface locate
    module procedure locate_1
    module procedure locate_2
end interface

interface selec
    module procedure selec_1
    module procedure selec_2
end interface

interface hpsort
    module procedure hpsort_1
    module procedure hpsort_2
    module procedure hpsort_3
    module procedure hpsort_4
    module procedure hpsort_5
end interface

interface hpsel
    module procedure hpsel_1
    module procedure hpsel_2
end interface

interface reverse
    module procedure reverse_iarr
    module procedure reverse_rarr
    module procedure reverse_carr
end interface

interface nvoxfind
    module procedure nvoxfind_1
    module procedure nvoxfind_2
end interface

interface gaussian
    module procedure gaussian_1
    module procedure gaussian_2
end interface

interface zeros
    module procedure zeros_1
    module procedure zeros_2
end interface

interface rad2deg
    module procedure rad2deg_1
    module procedure rad2deg_2
end interface

interface csq
    module procedure csq_1
    module procedure csq_2
end interface

logical, parameter,private :: warn=.false.

contains

    ! JIFFYS

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
    elemental function deg2rad( deg ) result( rad )
        real, intent(in) :: deg  !< angle (degrees)
        real             :: rad  !< angle (radians)
        rad = (deg/180.)*pi
    end function deg2rad

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

    !>   converts from correlation to euclidean distance
    ! pure function corr2dist( corr ) result( dist )
    !     real, intent(in) :: corr   !< query correlation
    !     real :: dist
    !     dist = 1. - corr
    ! end function corr2dist

    ! !>   converts from euclidean distance to correlation
    ! pure function dist2corr( dist ) result( corr )
    !     real, intent(in) :: dist  !< query distance
    !     real :: corr
    !     corr = dist + 1
    ! end function

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
            write(*,*) 'found NaNs in inputted vector; simple_math::check4nans_1', n_nans
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
            write(*,*) 'found NaNs in inputted vector; simple_math::check4nans_2', n_nans
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

    !>   get the angular resultion in degrees, given diameter and resolution
    !! \param res,diam spatial resolution (\f$\si{\per\angstrom}\f$)
    !! \param diam diameter (\f$\si{\angstrom}\f$)
    pure function angres( res, diam ) result( ang )
        real, intent(in)  :: res, diam
        real :: ang                     !< angular resolution in degrees
        ang = (res/(pi*diam))*360.
    end function angres

    !>   get the resolution in angstrom, given angle and diameter
    !! \param ang,diam angular resolution (degrees) and diameter (\f$\si{\angstrom}\f$)
    pure function resang( ang, diam ) result( res )
        real, intent(in)  :: ang, diam
        real :: res                      !< spatial resolution (\f$\si{\per\angstrom}\f$)
        res = (ang/360.)*(pi*diam)
    end function resang

     !>   checking for is_a_number
    pure function is_a_number_1( number ) result( is )
        real, intent(in) :: number  !< input variable for checking
        logical :: is
        is = .true.
        if( number > 0. )then
        else if( number <= 0. )then
        else
            is = .false.
        endif
    end function is_a_number_1

    !>   validity check of complex number (so that it is not nan)
    pure function is_a_number_2( complex_number ) result( is )
        complex, intent(in) :: complex_number !< input variable for checking
        logical             :: is
        is = is_a_number_1(real(complex_number)) .and. is_a_number_1(aimag(complex_number))
    end function is_a_number_2

    !>   converts string descriptors of c and d pointgroups to euler angle limits
    !! \param  t1,t2,p1,p2  euler angle limits
    subroutine pgroup_to_lim(pgroup, p1, p2, t1, t2, csym )
        character(len=*), intent(in) :: pgroup !< pointgroup
        integer, intent(out)         :: csym   !< symmetry token
        real, intent(out)            :: t1, t2, p1, p2
        if( pgroup(1:1) .eq. 'c' )then
            t1     = 0.
            t2     = 180.
        else if( pgroup(1:1) .eq. 'd' )then
            t1     = 0.
            t2     = 90.
        else
            write(*,*) 'error, presently only c and d symmetries are supported!'
            write(*,*) 'in: pgroup_to_lim, module: simple_math.f90'
            stop
        endif
        read(pgroup(2:),'(i2)') csym
        p1 = 0.
        p2 = 359.9999/real(csym)
    end subroutine pgroup_to_lim

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

    !>   reverses a complex array
    subroutine reverse_carr( carr )
        complex, intent(inout) :: carr(:) !< array for modification
        integer                :: i, j, sz, en
        complex                :: cswap
        sz = size(carr,1)
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
            cswap   = carr(j)
            carr(j) = carr(i)
            carr(i) = cswap
        end do 
    end subroutine reverse_carr

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
        if(alloc_stat /= 0) call alloc_errchk("sortmeans; simple_math", alloc_stat)
        ! initialization by sorting
        dat_sorted = dat ! reallocation  dat_sorted(ndat),
        call hpsort(ndat, dat_sorted)
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

    !>   calculates the number of common integers in two arrays
    function common_ints( arr1, arr2 ) result( n )
        integer, intent(in) :: arr1(:), arr2(:) !< arrays for comparison
        integer :: i, j, n
        n = 0
        do i=1,size(arr1)
            do j=1,size(arr2)
                if( arr1(i) == arr2(j) ) n = n + 1
            end do
        end do
    end function common_ints

    !>   to enforce cyclic limit
    pure subroutine enforce_cyclic_limit( x, lim )
        real, intent(inout) :: x   !< input var
        real, intent(in)    :: lim !< limit
        do while( x > lim )
            x = x-lim
        end do
        do while( x < 0. )
            x = x+lim
        end do
    end subroutine enforce_cyclic_limit

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

    !>   4 shifting variables
    subroutine shft(a,b,c,d)
        real, intent(out)   :: a   !< new pos 1
        real, intent(inout) :: b,c !< pos 2 and 3
        real, intent(in)    :: d   !< pos4
        a = b
        b = c
        c = d
    end subroutine shft

    ! !>    one-dimensional symmetric hard window
    pure function sqwin_1d( x, winsz ) result( win )
        real, intent(in) :: x       !< input point
        real, intent(in) :: winsz   !< window size
        integer          :: iwinsz, win(2) !< starts & stops
        win(:) = nint(x)
        iwinsz = ceiling(winsz)
        win(1) = win(1)-iwinsz
        win(2) = win(2)+iwinsz
    end function sqwin_1d

    !>    two-dimensional symmetric hard window
    !! \param x,y      input points
    pure function sqwin_2d( x, y, winsz ) result( win )
        real, intent(in) :: x,y      !< input point
        real, intent(in) :: winsz    !< window size
        integer          :: win(2,2) !< starts & stops
        win(1,:) = sqwin_1d(x,winsz)
        win(2,:) = sqwin_1d(y,winsz)
    end function sqwin_2d

    !>    three-dimensional symmetric hard window
    !! \param x,y,z      input points
    pure function sqwin_3d( x, y, z, winsz ) result( win )
        real, intent(in) :: x,y,z    !< input point
        real, intent(in) :: winsz    !< window size
        integer          :: win(3,2) !< starts & stops
        win(1,:) = sqwin_1d(x,winsz)
        win(2,:) = sqwin_1d(y,winsz)
        win(3,:) = sqwin_1d(z,winsz)
    end function sqwin_3d

    !>    one-dimensional hard window
    pure function recwin_1d( x, winsz ) result( win )
        real, intent(in) :: x       !< input point
        real, intent(in) :: winsz   !< window size
        integer          :: win(2) !< starts & stops
        win(1) = floor(x-real(winsz))
        win(2) = ceiling(x+real(winsz))
    end function recwin_1d

    !>    two-dimensional hard window
    !! \param x,y      input points
    pure function recwin_2d( x, y, winsz ) result( win )
        real, intent(in) :: x, y      !< input point
        real, intent(in) :: winsz     !< window size
        integer          :: win(2,2) !< starts & stops
        win(1,:) = recwin_1d(x,winsz)
        win(2,:) = recwin_1d(y,winsz)
    end function recwin_2d

    !>    three-dimensional hard window
    !! \param x,y,z      input points
    pure function recwin_3d( x, y, z, winsz ) result( win )
        real, intent(in) :: x, y, z   !< input point
        real, intent(in) :: winsz     !< window size
        integer          :: win(3,2) !< starts & stops
        win(1,:) = recwin_1d(x,winsz)
        win(2,:) = recwin_1d(y,winsz)
        win(3,:) = recwin_1d(z,winsz)
    end function recwin_3d

    ! USEFUL MATHEMATICAL FUNCTIONS

    !>   takes the logarithm of the positive elements of an array
    subroutine logarr(arr)
        real, intent(inout) :: arr(:)  !< input array
        integer :: i
        do i=1,size(arr)
            if( arr(i) > 0. )then
                arr(i) = log(arr(i))
            endif
        end do
    end subroutine logarr

    !>   returns acos with the argument's absolute value limited to 1.
    !!         this seems to be necessary due to small numerical inaccuracies.
    pure function myacos( arg ) result( r )
        real, intent(in) :: arg     !< input (radians)
        real             :: r, x, y
        x = min(1.,abs(arg))
        y = sign(x,arg)
        r = acos(y)
    end function myacos

    !>   sinc function
    function sinc( x ) result( r )
        real, intent(in) :: x       !< input (radians)
        real             :: r, arg
        if( abs(x) < 0.00000001 ) then
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

    !>   is a Gaussian function (for calculating the probability given
    !! \param x input value
    !! \param mean Gaussian mean
    !! \param sigma Gaussian std. dev.
    function gaussian_1( x, mean, sigma ) result( r )
        real, intent(in) :: x, mean, sigma
        real :: r
        r = (1./(sigma*sqrt(twopi)))*exp(-(x-mean)**2./(2.*(sigma*sigma)))
    end function gaussian_1

    !>   is a Gaussian function (for calculating the probability given dist & sigma)
    function gaussian_2( dist, sigma ) result( r )
        real, intent(in) :: dist, sigma
        real :: r
        r = (1./(sigma*sqrt(twopi)))*exp(-dist/(2.*(sigma*sigma)))
    end function gaussian_2

    ! COMPLEX STUFF

    !>   is for calculating the phase angle of a Fourier component
    !! \return phase phase angle only meaningful when cabs(comp) is well above the noise level
    elemental function phase_angle( comp ) result( phase )
        complex, intent(in) :: comp !< Fourier component
        real :: nom, denom, phase
        nom = aimag(comp)
        denom = real(comp)
#ifdef USETINY
        if(  abs(denom) < TINY )then
           if( abs(nom) < TINY )then
#else
        if( denom == 0. )then
           if( nom == 0. )then
#endif
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
#ifdef USETINY
        if( x < TINY ) then
            sq = y*y
        else if( y < TINY ) then
#else
        if( x == 0.) then
            sq = y*y
        else if( y == 0. ) then
#endif
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
#ifdef USETINY
        if( x < TINY ) then
            sq = y*y
        else if( y < TINY ) then
#else
        if( x == 0.) then
            sq = y*y
        else if( y == 0. ) then
#endif  
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
#ifdef USETINY
        if( x < TINY ) then
            myabs = y
        else if( y < TINY ) then
#else
        if( x == 0. ) then
            myabs = y
        else if( y == 0. ) then
#endif
           myabs = x
        else if( x > y ) then
            frac = y/x
            myabs = x*sqrt(1.+frac*frac)
        else
            frac = x/y
            myabs = y*sqrt(1.+frac*frac)
        endif
    end function mycabs

    !>   normalized correlation coefficient between two complex numbers
    function ccorr( c1, c2 ) result( corr )
        complex, intent(in) :: c1, c2 !< complx components
        real    :: corr, c1sq, c2sq
        corr = 0.
        c1sq = csq(c1)
        c2sq = csq(c2)
#ifdef USETINY
        if( c1sq > TINY .and. c2sq > TINY )then
#else
        if( c1sq /= 0. .and. c2sq /= 0. )then
#endif
            corr = calc_corr(real(c1*conjg(c2)),sqrt(c1sq*c2sq))
        endif
   end function ccorr

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

    !>   two-dimensional hard edge
    !! \f$r < \sqrt{x^2+y^2}\f$.
    !! \param x x position
    !! \param y y position
    !! \param mskrad masking radius
    !! \return w on or off
    !!
    pure function hardedge_inner_1( x, y, mskrad ) result( w )
        real,intent(in) :: x, y, mskrad
        real :: w
        w = 0.
        if( x * x + y * y > mskrad * mskrad ) w = 1.
    end function hardedge_inner_1

    !>   three-dimensional hard edge
    !! \f$r < \sqrt{ x^2+y^2+z^2 }\f$.
    !! \param x x position
    !! \param y y position
    !! \param z z position
    !! \param mskrad masking radius
    !! \return w on or off
    !!
    pure function hardedge_inner_2( x, y, z, mskrad ) result( w )
        real,intent(in) :: x, y, z, mskrad
        real :: w
        w = 0.
        if( x * x + y * y + z * z > mskrad * mskrad ) w = 1.
    end function

    !>   low-end one-dimensional gaussian edge
    !! \param d distance
    !! \param mskrad masking radius
    !! \return w off or graded edge
    !!
    pure function cosedge_1( d,rad ) result( w )
        real, intent(in) :: d,rad
        real             :: w,halfrad
        halfrad = rad/2.
        if( d.lt.halfrad )then
            w = (cos(pi*d/halfrad)+1.)/2.
        else
            w = 0.
        endif
    end function cosedge_1

    !>   two-dimensional gaussian edge
    !! \param x x position
    !! \param y y position
   pure function cosedge_2( x, y, box, mskrad ) result( w )
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
    end function cosedge_2

    !>   three-dimensional gaussian edge
    !! \f$r = \cos{(1+{(\pi{r - (d+2m)/(d-2m)})})}\f$.
    !! \param x x position
    !! \param y y position
    pure function cosedge_3( x, y, z, box, mskrad ) result( w )
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
    end function cosedge_3

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

    !>   returns the shell to which the voxel with address h,k,l belongs
    !! \f$s = \nint{\sqrt{h^2+k^2+l^2}}\f$.
    !! \param h,k,l voxel indices
    pure function shell( h, k, l ) result( sh )
        integer, intent(in) :: h, k, l
        integer :: sh
        sh = nint(sqrt(real(h**2+k**2+l**2)))
    end function shell

    !>   calculates the resolution values given corrs and res params
    !! \param corrs Fourier shell correlations
    !! \param res resolution value
    subroutine get_resolution( corrs, res, fsc05, fsc0143 )
        real, intent(in)  :: corrs(:), res(:) !<  corrs Fourier shell correlation
        real, intent(out) :: fsc05, fsc0143 !<  fsc05 resolution at FSC=0.5,  fsc0143 resolution at FSC=0.143
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
        if( ires0143 == n .or. ires0143 == 0 )then
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
        if( ires05 == n .or. ires05 == 0 )then
            fsc05 = 0.
        else
            fsc05 = res(ires05)
        endif
    end subroutine get_resolution

    !>   returns the Fourier index of the resolution limit
    function get_lplim( fsc ) result( k )
        real, intent(in) :: fsc(:) !< Fourier shell correlation array
        integer :: n, k, h
        n = size(fsc)
        if( n < 3 )then
            stop 'nonconforming size of fsc array; get_lplim; simple_math'
        endif
        k = n-1
        do h=3,n-1
            if( fsc(h) >= 0.143 )then
                cycle
            else
                k = h
                exit
            endif
        end do
    end function get_lplim

    !>   compares two Fourier Shell Correlation (FSC) functions
    logical function fsc1_ge_fsc2( fsc1, fsc2 )
        real, intent(in) :: fsc1(:), fsc2(:) !< Fourier shell correlation arrays
        integer          :: n, hplim, lplim, nfreq, n1, n2, h
        real, parameter  :: FSC_MAX = 0.98
        n = size(fsc1)
        if( n /= size(fsc2) ) stop 'nonconforming fsc sizes; simple_math :: fsc1_ge_fsc2'
        lplim = max(get_lplim(fsc1),get_lplim(fsc2))
        hplim = 2
        do h=3,n
            if( fsc1(h) >= FSC_MAX .and. fsc2(h) >= FSC_MAX )then
                cycle
            else
                hplim = h
                exit
            endif
        end do
        nfreq = lplim - hplim + 1
        n1    = count(fsc1(hplim:lplim) >= fsc2(hplim:lplim))
        n2    = nfreq - n1
        fsc1_ge_fsc2 = (n1 >= n2)
    end function fsc1_ge_fsc2

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
        call alloc_errchk('In: get_res, module: simple_math', alloc_stat)
        do k=1,n
            res(k) = calc_lowpass_lim(k, box, smpd)
        end do
    end function get_resarr

    !>   calculates a corr coeff based on sum cross prod and denominator
    !! \param sxy cross prod sum
    !! \param den denominator
    function calc_corr( sxy, den ) result( corr )
        real, intent(in) :: sxy, den
        real :: corr                 !< output corr coeff
        if( den > 0. )then
            corr = sxy/sqrt(den)
            if( is_a_number(corr) )then
                if( corr > 1. .and. warn )then
                    write(*,*) 'WARNING! corr > 1, numerical errors', corr
                else if( corr < -1. .and. warn )then
                    write(*,*) 'WARNING! corr < -1, numerical errors', corr
                endif
                corr = min(1.,max(-1.,corr))
            else
                corr = 0.
            endif
        else
            corr = 0.
        endif
    end function calc_corr

    ! NUMERICAL STUFF

    !>    This routine computes the nth stage of refinement of an extended trapezoidal
    !!          rule. When
    !!          called with n=1, the routine returns as s the crudest estimate of the integral.
    !!          Subsequent calls with n=2,3,... (in that sequential order) will improve the
    !!          accuracy of s by adding 2**(n-2) additional interior points. s should not be
    !!          modified between sequential calls (from Numerical Recepies)
    !! \param  func is the function to be integrated between limits a and b.
    !! \param  a lower limits of func
    !! \param  b upper limits of func
    !!
    subroutine trapzd( func, a, b, s, n )
        interface
            function func( point ) result( val )
                real, intent(in) :: point
                real :: val
            end function
        end interface
        real, intent(in)    :: a, b
        real, intent(inout) :: s      !< integral output
        integer, intent(in) :: n      !< query index
        integer             :: it, j
        real                :: del, sum, tnm, x
        if( n .eq. 1 )then
            s = 0.5*(b-a)*(func(a)+func(b))
        else
            it  = 2**(n-2)
            tnm = it
            del = (b-a)/tnm ! this is the spacing of the points to be added
            x   = a+0.5*del
            sum = 0.
            do j=1,it
                sum = sum+func(x)
                x = x+del
            end do
            s = 0.5*(s+(b-a)*sum/tnm) ! replaces s by its refined value
        endif
    end subroutine trapzd

    !>    Returns as s the integral of the function func from a to b. The parameters EPS
    !!          can be set to the desired fractional accuracy and JMAX so that 2 to the power
    !!          JMAX-1 is the maximum allowed number of steps. Integration is performed by Simpson's
    !!          rule (from Numerical Recepies)
    !! \param  func is the function to be integrated between limits a and b.
    !! \param  a limits of func
    function qsimp( func, a, b ) result( s )
        interface
            function func( point ) result( val )
                real, intent(in) :: point
                real :: val
            end function
        end interface
        real, intent(in)   :: a, b
        real               :: s, os, ost, st
        real, parameter    :: EPS=1.e-5
        integer, parameter :: JMAX=20
        integer            :: j
        ost = -1.e30
        os  = -1.e30
        do j=1,JMAX
            call trapzd(func,a,b,st,j)
            s = (4.*st-ost)/3.
            if( abs(s-os) .lt. EPS*abs(os) ) return
            if( s == 0. .and. os == 0. .and. j.gt.6 ) return
            os = s
            ost = st
        end do
        !write(*,'(a)') 'too many steps in qsimp; simple_math'
    end function qsimp

    !>    Returns the derivative of a function func at a point x by Ridders'
    !!          method of polynomial extrapolation. The value h is input as an estimated
    !!          initial stepsize; it need not to be small but rather should be an
    !!          increment in x over which func changes substantially. An estimate of
    !!          the error in the derivative is returned in err
    !! \param func is the function to be derived at point x.
    !! \param x point query
    !! \param err  estimated error 
    subroutine numderiv( func, x, h, err, df )
    ! Parameters: Stepsize is decreased by CON at each iteration. Max size of tableu is
    ! set by NTAB. Return when error is SAFE worse than the best so far. The number of
    ! function evaluations are typically 6 to 12 byt is allowed to be as many as 2*NTAB.
    ! You should therefore select a fairly large value for h but monitor the returned
    ! err, decreasing h if it is not small. For functions whose characteristic x scale
    ! is of order unity-take h to be a few tenths (5 degrees for correlation search)
        interface
            function func( point ) result( val )
                real, intent(in) :: point
                real :: val
            end function
        end interface
        integer, parameter :: NTAB=10                                   !< max size of tableu
        real, parameter    :: CON=1.4, CON2=CON*CON, BIG=1.E30, SAFE=2. !< real params
        real, intent(in)   :: x, h                                      !< point & step
        real, intent(out)  :: err, df                                   !< gradient
        real               :: errt, fac, hh, a(NTAB,NTAB)
        integer            :: i, j
        if(abs(h)<TINY) stop 'h must be nonzero; numderiv; simple_math'
        hh = h
        a(1,1) = (func(x+hh)-func(x-hh))/(2.*hh)
        err = BIG
        do i=2,NTAB          ! successive columns in the Neville tableu will go to
            hh     = hh/CON  ! smaller stepsizes and higher orders of extrapolation
            a(1,i) = (func(x+hh)-func(x-hh))/(2.*hh) ! try new smaller stepsize
            fac    = CON2
            do j=2,NTAB      ! compute extrapolations of various orders, no new evals required
                a(j,i) = (a(j-1,i)*fac-a(j-1,i-1))/(fac-1.)
                fac    = CON2*fac
                errt   = max(abs(a(j,i)-a(j-1,i)), abs(a(j,i)-a(j-1,i-1)))
                ! The error strategy is to compare each new extrapolation to one order lower,
                ! both at the present stepsize and the previous one
                if( errt .le. err )then ! If error decreased save the better answer
                    err = errt
                    df  = a(j,i)
                endif
            end do
            ! If higher order is worse by a significant factor SAFE, then quit early
            if( abs(a(i,i)-a(i-1,i-1)) .ge. SAFE*err ) return
        end do
    end subroutine numderiv

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
#ifdef USETINY
           if( abs(augmatrix(k,k)) < TINY )then
#else
              if( augmatrix(k,k) == 0. )then
#endif
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
        !test for invertibility
        do i=1,n
            if( augmatrix(i,i) == 0. )then
                inverse = 0
                errflg = -1
                return
            endif
        end do
        !make diagonal elements as 1
        do i=1,n
            m = augmatrix(i,i)
            do j=i,2*n
                augmatrix(i,j) = augmatrix(i,j)/m
            end do
        end do
        !reduced right side half of augmented matrix to identity matrix
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

    !>   subroutine to find the inverse of a square matrix
    !!         author : louisda16th a.k.a ashwith j. rego
    !!         reference : algorithm explained at:
    !!         http://math.uww.edu/~mcfarlat/inverse.htm
    !!         http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
    subroutine matinv_D(matrix, inverse, n, errflg)
        integer, intent(in)  :: n                          !< size of square matrix
        integer, intent(out) :: errflg                   !< return error status. -1 for error, 0 for normal
        real(dp), intent(in), dimension(n,n)  :: matrix  !< input matrix
        real(dp), intent(out), dimension(n,n) :: inverse !< inverted matrix
        logical :: flag = .true.
        integer :: i, j, k
        real(dp) :: m
        real(dp), dimension(n,2*n) :: augmatrix !augmented matrix
        ! augment input matrix with an identity matrix
        do i=1,n
            do j=1,2*n
                if(j <= n )then
                    augmatrix(i,j) = matrix(i,j)
                else if((i+n) == j)then
                    augmatrix(i,j) = 1.0d0
                else
                    augmatrix(i,j) = 0.0d0
                endif
            end do
        end do
        ! reduce augmented matrix to upper traingular form
        do k=1,n-1
            if( augmatrix(k,k) == 0.d0 )then
                flag = .false.
                do i=k+1,n
                    if( augmatrix(i,k) > 0.d0 )then
                        do j=1,2*n
                            augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
                        end do
                        flag = .true.
                        exit
                    endif
                    if(flag .eqv. .false.)then
                        inverse = 0.d0
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
        !test for invertibility
        do i=1,n
            if( augmatrix(i,i) == 0.d0 )then
                inverse = 0.d0
                errflg = -1
                return
            endif
        end do
        !make diagonal elements as 1
        do i=1,n
            m = augmatrix(i,i)
            do j=i,2*n
                augmatrix(i,j) = augmatrix(i,j)/m
            end do
        end do
        !reduced right side half of augmented matrix to identity matrix
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
    end subroutine matinv_D

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

    !>   trace of a square matrix
    pure function tr( matrix, n ) result( sum )
        integer, intent(in) :: n          !< matrix size
        real, intent(in), dimension(n,n) :: matrix !< input matrix
        real :: sum
        integer :: i
        sum = 0.
        do i=1,n
            sum = sum+matrix(i,i)
        end do
    end function tr

    !>   determinant of a square matrix, matrix is modified
    function det( matrix, n ) result( d )
        integer, intent(in)   :: n          !< matrix size
        real, dimension(n,n)  :: matrix     !< input matrix
        integer, dimension(n) :: indx(n)
        real                  :: d
        logical               :: err
        integer               :: j
        call ludcmp(matrix,n,n,indx,d,err)
        if( err )then
            d = 0.
            return
        endif
        do j=1,n
            d = d*matrix(j,j) ! this returns d as +-1
        end do
    end function det

    !>   determinant of a square matrix, matrix is modified
    function det_D( matrix, n ) result( d )
        integer, intent(in)   :: n           !< matrix size
        real(dp), dimension(n,n)  :: matrix
        integer, dimension(n) :: indx(n)
        real(dp)                  :: d
        logical               :: err
        integer               :: j
        call ludcmp_D(matrix,n,n,indx,d,err)
        if( err )then
            d = 0.d0
            return
        endif
        do j=1,n
            d = d*matrix(j,j) ! this returns d as +-1
        end do
    end function det_D

    !>   diagonal of a square matrix
    pure function diag( diagvals, n ) result( mat )
        integer, intent(in) :: n          !< matrix size
        real,    intent(in) :: diagvals(n) !< input vector
        real :: mat(n,n)
        integer :: i
        mat = 0.
        do i=1,n
            mat(i,i) = diagvals(i)
        end do
    end function diag

    !>   to allocate a one-dim array with zeros
    function zeros_1( n ) result( a )
        integer, intent(in) :: n          !< matrix size
        real, allocatable   :: a(:)
        allocate( a(n), stat=alloc_stat )
        if(alloc_stat /= 0) call alloc_errchk("In: zeros_1; simple_math", alloc_stat)
        a = 0.
    end function zeros_1

    !>   to allocate a two-dim array with zeros
    function zeros_2( n1, n2 ) result( a )
        integer, intent(in) :: n1, n2    !< matrix size
        real, allocatable   :: a(:,:)
        allocate( a(n1,n2), stat=alloc_stat )
        if(alloc_stat /= 0) call alloc_errchk("In: zeros_2; simple_math", alloc_stat)
        a = 0.
    end function zeros_2

    !>   lu decomposition, nr
    !! \see Numerical recipes

    subroutine ludcmp(a,n,np,indx,d,err)
        integer :: n,np,indx(n),nmax
        real    :: d,a(np,np), tiny
        parameter (nmax=100,tiny=1.0e-20)
        integer :: i,imax=0,j,k
        real    :: aamax,dum,sum,vv(nmax)
        logical :: err
        err = .false.
        d=1.
        do i=1,n
          aamax=0.
          do j=1,n
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
          end do
          if( aamax == 0. )then
              err = .true.
              return
          endif
          vv(i)=1./aamax
        end do
        do j=1,n
          if (j.gt.1) then
            do i=1,j-1
              sum=a(i,j)
              if (i.gt.1)then
                do k=1,i-1
                  sum=sum-a(i,k)*a(k,j)
                end do
                a(i,j)=sum
              endif
            end do
          endif
          aamax=0.
          do i=j,n
            sum=a(i,j)
            if (j.gt.1)then
              do k=1,j-1
                sum=sum-a(i,k)*a(k,j)
              end do
              a(i,j)=sum
            endif
            dum=vv(i)*abs(sum)
            if (dum.ge.aamax) then
              imax=i
              aamax=dum
            endif
          end do
          if (j.ne.imax)then
            do k=1,n
              dum=a(imax,k)
              a(imax,k)=a(j,k)
              a(j,k)=dum
            end do
            d=-d
            vv(imax)=vv(j)
          endif
          indx(j)=imax
          if(j.ne.n)then
            if(abs(a(j,j))<tiny) a(j,j)=tiny
            dum=1./a(j,j)
            do i=j+1,n
              a(i,j)=a(i,j)*dum
          end do
          endif
        end do
        if(abs(a(n,n))<tiny)a(n,n)=tiny
    end subroutine ludcmp

    !>   lu decomposition, nr
    subroutine ludcmp_D(a,n,np,indx,d,err)
        integer :: n,np,indx(n),nmax
        real(dp)    :: d,a(np,np), tiny
        parameter (nmax=100,tiny=1.0e-20)
        integer :: i,imax=0,j,k
        real(dp)    :: aamax,dum,sum,vv(nmax)
        logical :: err
        err = .false.
        d=1.
        do i=1,n
          aamax=0.d0
          do j=1,n
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
          end do
          if( aamax == 0.d0 )then
              err = .true.
              return
          endif
          vv(i)=1./aamax
        end do
        do j=1,n
          if (j.gt.1) then
            do i=1,j-1
              sum=a(i,j)
              if (i.gt.1)then
                do k=1,i-1
                  sum=sum-a(i,k)*a(k,j)
                end do
                a(i,j)=sum
              endif
            end do
          endif
          aamax=0.d0
          do i=j,n
            sum=a(i,j)
            if (j.gt.1)then
              do k=1,j-1
                sum=sum-a(i,k)*a(k,j)
              end do
              a(i,j)=sum
            endif
            dum=vv(i)*abs(sum)
            if (dum.ge.aamax) then
              imax=i
              aamax=dum
            endif
          end do
          if (j.ne.imax)then
            do k=1,n
              dum=a(imax,k)
              a(imax,k)=a(j,k)
              a(j,k)=dum
            end do
            d=-d
            vv(imax)=vv(j)
          endif
          indx(j)=imax
          if(j.ne.n)then
            if(abs(a(j,j))<tiny) a(j,j)=tiny
            dum=1./a(j,j)
            do i=j+1,n
              a(i,j)=a(i,j)*dum
          end do
          endif
        end do
        if(abs(a(n,n))<tiny)a(n,n)=tiny
    end subroutine ludcmp_D

    !>  jacobi SVD, NR
    !! \see Numerical recipes

    subroutine jacobi( a, n, np, d, v, nrot)
        integer, intent(in)    :: n,np
        real,    intent(inout) :: a(np,np), v(np,np), d(np)
        integer, intent(inout) :: nrot
        real                   :: c,g,h,s,sm,t,tau,theta,tresh,b(n), z(n)
        integer                :: i,j,ip,iq
        v = 0.
        do i=1,n
            v(i,i) = 1.
            b(i)   = a(i,i)
        enddo
        d    = b
        z    = 0.
        nrot = 0
        do i=1,50
            sm=0.
            do ip=1,n-1
                do iq=ip+1,n
                    sm=sm+abs(a(ip,iq))
                enddo
            enddo
            if(sm.eq.0.)return
            if(i.lt.4)then
                tresh=0.2*sm/n**2
            else
                tresh=0.
            endif
            do ip=1,n-1
                do iq=ip+1,n
                    g=100.*abs(a(ip,iq))
                    if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))) &
                        .and.(abs(d(iq))+g.eq.abs(d(iq))))then
                        a(ip,iq)=0.
                    else if(abs(a(ip,iq)).gt.tresh)then
                        h=d(iq)-d(ip)
                        if(abs(h)+g.eq.abs(h))then
                            t=a(ip,iq)/h
                        else
                            theta=0.5*h/a(ip,iq)
                            t=1./(abs(theta)+sqrt(1.+theta**2))
                            if(theta.lt.0.)t=-t
                        endif
                        c=1./sqrt(1+t**2)
                        s=t*c
                        tau=s/(1.+c)
                        h=t*a(ip,iq)
                        z(ip)=z(ip)-h
                        z(iq)=z(iq)+h
                        d(ip)=d(ip)-h
                        d(iq)=d(iq)+h
                        a(ip,iq)=0.
                        do j=1,ip-1
                            g=a(j,ip)
                            h=a(j,iq)
                            a(j,ip)=g-s*(h+g*tau)
                            a(j,iq)=h+s*(g-h*tau)
                        enddo
                        do j=ip+1,iq-1
                            g=a(ip,j)
                            h=a(j,iq)
                            a(ip,j)=g-s*(h+g*tau)
                            a(j,iq)=h+s*(g-h*tau)
                        enddo
                        do j=iq+1,n
                            g=a(ip,j)
                            h=a(iq,j)
                            a(ip,j)=g-s*(h+g*tau)
                            a(iq,j)=h+s*(g-h*tau)
                        enddo
                        do j=1,n
                            g=v(j,ip)
                            h=v(j,iq)
                            v(j,ip)=g-s*(h+g*tau)
                            v(j,iq)=h+s*(g-h*tau)
                        enddo
                        nrot=nrot+1
                    endif
                enddo
            enddo
            do ip=1,n
                b(ip)=b(ip)+z(ip)
                d(ip)=b(ip)
                z(ip)=0.
            enddo
        enddo
        write(*,*)' Too many iterations in Jacobi'
    end subroutine jacobi

    !>   sorts eigenvalues and eigenvectors from jacobi routine in descending order
    subroutine eigsrt(d,v,n,np)
        integer, intent(in)    :: n,np
        real,    intent(inout) :: d(np),v(np,np)
        ! Given the eigenvalues d and eigenvectors v as output from jacobi (section 11.1) or tqli (section 11.3),
        ! this routine sorts the eigenvalues into descending order,
        ! and rearranges the columns of v correspondingly. The method is straight insertion.
        integer             :: i,j,k
        real                :: p
        do i=1,n-1
            k = i
            p = d(i)
            do j=i+1,n
                if(d(j)>p)then
                    k = j
                    p = d(j)
                endif
            enddo
            if(k.ne.i)then
                d(k) = d(i)
                d(i) = p
                do j=1,n
                    p      = v(j,i)
                    v(j,i) = v(j,k)
                    v(j,k) = p
                enddo
            endif
        enddo
    end subroutine eigsrt

    ! Given a matrix a[1..m][1..n], this routine computes its singular value decomposition, A =
    ! U ·W ·V T . The matrix U replaces a on output. The diagonal matrix of singular values W is output
    ! as a vector w[1..n]. The matrix V (not the transpose V T ) is output as v[1..n][1..n].
    subroutine svdcmp(a,m,n,mp,np,w,v)
      integer :: m,mp,n,np!,nmax
      real, intent(inout) :: a(mp,np),v(np,np),w(np)
      integer ::i,its,j,k,l,nm !,jj
      real :: anorm,c,f,g,h,s,scale,x,y,z,rv1(n)!,pythag
      g=0.0
      scale=0.0
      anorm=0.0
      do 25 i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.0
        s=0.0
        scale=0.0
        if (i.le.m) then
          do 11 k=i,m
            scale=scale+abs(a(k,i))
11        continue
          if (scale.ne.0.) then
            do 12 k=i,m
              a(k,i)=a(k,i)/scale
              s=s+a(k,i)*a(k,i)
12          continue
            f=a(i,i)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            if (i.ne.n) then
              do 15 j=l,n
                s=0.0
                do 13 k=i,m
                  s=s+a(k,i)*a(k,j)
13              continue
                f=s/h
                do 14 k=i,m
                  a(k,j)=a(k,j)+f*a(k,i)
14              continue
15            continue
           endif
            do 16 k= i,m
              a(k,i)=scale*a(k,i)
16          continue
          endif
        endif
        w(i)=scale *g
        g=0.0
        s=0.0
        scale=0.0
        if ((i.le.m).and.(i.ne.n)) then
          do 17 k=l,n
            scale=scale+abs(a(i,k))
17        continue
          if (scale.ne.0.) then
            do 18 k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
18          continue
            f=a(i,l)
            g=-sign(sqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do 19 k=l,n
              rv1(k)=a(i,k)/h
19          continue
            if (i.ne.m) then
              do 23 j=l,m
                s=0.0
                do 21 k=l,n
                  s=s+a(j,k)*a(i,k)
21              continue
                do 22 k=l,n
                  a(j,k)=a(j,k)+s*rv1(k)
22              continue
23            continue
            endif
            do 24 k=l,n
              a(i,k)=scale*a(i,k)
24          continue
          endif
        endif
        anorm=max(anorm,(abs(w(i))+abs(rv1(i))))
25    continue
      do 32 i=n,1,-1
        if (i.lt.n) then
          if (g.ne.0.) then
            do 26 j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
26          continue
            do 29 j=l,n
              s=0.0
              do 27 k=l,n
                s=s+a(i,k)*v(k,j)
27            continue
              do 28 k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
28            continue
29          continue
          endif
          do 31 j=l,n
            v(i,j)=0.0
            v(j,i)=0.0
31        continue
        endif
        v(i,i)=1.0
        g=rv1(i)
        l=i
32    continue
      do 39 i=n,1,-1
        l=i+1
        g=w(i)
        if (i.lt.n) then
          do 33 j=l,n
            a(i,j)=0.0
33        continue
        endif
        if (g.ne.0.) then
          g=1.0/g
          if (i.ne.n) then
            do 36 j=l,n
              s=0.0
              do 34 k=l,m
                s=s+a(k,i)*a(k,j)
34            continue
              f=(s/a(i,i))*g
              do 35 k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
35            continue
36          continue
          endif
          do 37 j=i,m
            a(j,i)=a(j,i)*g
37        continue
        else
          do 38 j= i,m
            a(j,i)=0.0
38        continue
        endif
        a(i,i)=a(i,i)+1.0
39    continue
      do 49 k=n,1,-1
        do 48 its=1,30
          do 41 l=k,1,-1
            nm=l-1
            if ((abs(rv1(l))+anorm).eq.anorm)  go to 2
            if ((abs(w(nm))+anorm).eq.anorm)  go to 1
41        continue
1         c=0.0
          s=1.0
          do 43 i=l,k
            f=s*rv1(i)
            if ((abs(f)+anorm).ne.anorm) then
              g=w(i)
              h=sqrt(f*f+g*g)
              w(i)=h
              h=1.0/h
              c= (g*h)
              s=-(f*h)
              do 42 j=1,m
                y=a(j,nm)
                z=a(j,i)
                a(j,nm)=(y*c)+(z*s)
                a(j,i)=-(y*s)+(z*c)
42            continue
            endif
43        continue
2         z=w(k)
          if (l.eq.k) then
            if (z.lt.0.0) then
              w(k)=-z
              do 44 j=1,n
                v(j,k)=-v(j,k)
44            continue
            endif
            go to 3
          endif
          if (its.eq.30) stop 'no convergence in 30 iterations'
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
          g=sqrt(f*f+1.0)
          f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
          c=1.0
          s=1.0
          do 47 j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=sqrt(f*f+h*h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do 45 nm=1,n
              x=v(nm,j)
              z=v(nm,i)
              v(nm,j)= (x*c)+(z*s)
              v(nm,i)=-(x*s)+(z*c)
45          continue
            z=sqrt(f*f+h*h)
            w(j)=z
            if (z.ne.0.0) then
              z=1.0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do 46 nm=1,m
              y=a(nm,j)
              z=a(nm,i)
              a(nm,j)= (y*c)+(z*s)
              a(nm,i)=-(y*s)+(z*c)
46          continue
47        continue
          rv1(l)=0.0
          rv1(k)=f
          w(k)=x
48      continue
3       continue
49    continue
      return
    end subroutine svdcmp

    !>   for generating a spiraling path in 2d
    subroutine spiral_2d( np, trs, vecs )
        integer, intent(in) :: np         !< total num of points
        real, intent(in)    :: trs        !< radial step size
        real, intent(out)   :: vecs(np,2) !< output vector of spiral coords
        real                :: angstep, ang, normstep, norm
        integer             :: i
        angstep  = 360./(real(np)/trs)
        normstep = trs/real(np)
        ang      = 0.
        norm     = 0.
        do i=1,np
            call get_radial_line( ang, vecs(i,:) )
            vecs(i,:) = norm*vecs(i,:)
            norm = norm+normstep
            ang  = ang+angstep
        end do
    end subroutine spiral_2d

    !>    is for calculating the radius
    function hyp( x1, x2, x3 ) result( h )
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

    !>   calculates the summmed difference between two vectors
    pure function sum_diff( vec1, vec2 ) result( diff_sum )
        real, intent(in)    :: vec1(:), vec2(:)
        real                :: diff_sum
        diff_sum = sum(vec1-vec2)
    end function sum_diff

    !>   calculates the summmed relative difference between two vectors
    pure function rel_sum_diff( vec1, vec2 ) result( diff_sum_rel )
        real, intent(in)    :: vec1(:), vec2(:)
        real                :: diff_sum_rel
        diff_sum_rel = sum((vec1-vec2)/vec2)
    end function rel_sum_diff

    !>   calculates the argument of a vector
    pure function arg( vec ) result( length )
        real, intent(in) :: vec(:)
        real :: length
        length = sqrt(sum(vec**2.))
    end function arg

    !>   normalizes the length of a vector to 1
    subroutine normvec( vec )
        real, intent(inout) :: vec(:)
        vec = vec/arg(vec)
    end subroutine normvec

    !>   projects a 3d vector in the _z_-direction
    subroutine projz( vec3, vec2 )
        real, intent(in)  :: vec3(3)
        real, intent(out) :: vec2(2)
        vec2(1) = vec3(1)
        vec2(2) = vec3(2)
    end subroutine projz

    !>   calculates a 2d vector in the _xy_ plane rotated _ang_ degrees.
    subroutine get_radial_line( ang, lin )
        real, intent(in) :: ang
        real, intent(out), dimension(2) :: lin
        real :: mat(2,2), u(2)
        mat = rotmat2d( ang )
        u(1) = 0.
        u(2) = 1.
        lin = matmul(u,mat)
    end subroutine get_radial_line

    !>   is for 2d rotation matrix generation
    pure function rotmat2d( ang ) result( mat )
        real, intent(in) :: ang ! in degrees
        real :: mat(2,2), ang_in_rad
        ang_in_rad = ang*pi/180.
        mat(1,1) = cos(ang_in_rad)
        mat(1,2) = sin(ang_in_rad)
        mat(2,1) = -mat(1,2)
        mat(2,2) = mat(1,1)
    end function rotmat2d

    !>  extracts in-plane parameters from transformation matrix
    subroutine transfmat2inpls( R, psi, tx, ty )
        real,intent(inout) :: psi,tx,ty
        real,intent(in)    :: R(3,3)
        psi = rad2deg( myacos( R(1,1) ))
        if( R(1,2)<0. )psi=360.-psi
        tx  = R(1,3)
        ty  = R(2,3)
    end subroutine transfmat2inpls

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

    !>   is for converting a rotation matrix to axis-angle representation
    subroutine rotmat2axis( R, axis )
        real, intent(in)  :: R(3,3)
        real, intent(out) :: axis(3)
        axis(1) = R(3,2) - R(2,3)
        axis(2) = R(1,3) - R(3,1)
        axis(3) = R(2,1) - R(1,2)
        ! normalization
        axis = axis / sqrt(dot_product(axis,axis))
    end subroutine rotmat2axis

    !>   generates polar coordinates
    subroutine gen_polar_coords( kfromto, ring2, coords, angtab )
        integer, intent(in)            :: kfromto(2), ring2
        real, allocatable, intent(out) :: coords(:,:,:), angtab(:)
        integer                        :: nradial_lines
        real                           :: dang
        integer                        :: i, j
        nradial_lines = round2even(twopi*real(ring2))
        if( allocated(coords) )then
            deallocate(coords)
        end if
        if( allocated(angtab) ) then
            deallocate(angtab)
        end if
        allocate( coords(nradial_lines,kfromto(1):kfromto(2),2), angtab(nradial_lines), stat=alloc_stat )
        if(alloc_stat /= 0) call alloc_errchk("In: gen_polar_coords_1; simple_math coords/angtab alloc", alloc_stat)
        dang = twopi/real(nradial_lines)
        do i=1,nradial_lines
            angtab(i) = real(i-1)*dang
            do j=kfromto(1),kfromto(2)
                coords(i,j,1) = cos(angtab(i))*real(j)
                coords(i,j,2) = sin(angtab(i))*real(j)
            end do
            angtab(i) = rad2deg(angtab(i))
        end do
    end subroutine gen_polar_coords

    ! INTERPOLATION

    !>    rational function interpolation & extrapolation, from NR
    !!          Given arrays xa(:) and ya(:) and a value of x, this routine
    !!          returns a value of y and an accuracy estimate dy
    !! \param xa,ya  input line
    subroutine ratint(xa, ya, x, y, dy)
        real, intent(in)  :: xa(:), ya(:), x  !< x position
        real, intent(out) :: y, dy            !< est value and accuracy
        real, parameter   :: TINY = 1.e-25
        real, allocatable :: c(:), d(:)
        integer :: n, i, m, ns
        real    :: dd, h, hh, t, w
        n = size(xa)
        if( n /= size(ya) ) stop 'incompatible array sizes; ratint; simple_math'
        allocate(c(n), d(n), stat=alloc_stat)
        if(alloc_stat /= 0) call alloc_errchk("In: ratint; simple_math", alloc_stat)
        ns = 1
        hh = abs(x-xa(1))
        do i=1,n
            h = abs(x-xa(i))
            if( h == 0. )then
                y  = ya(i)
                dy = 0.
                return
            else if( h.lt.hh )then
                ns = i
                hh = h
            endif
            c(i) = ya(i)
            d(i) = ya(i)+TINY ! TINY needed to prevent rare zero-over-zero
        end do
        y  = ya(ns)
        ns = ns-1
        do m=1,n-1
            do i=1,n-m
                w = c(i+1)-d(i)
                h = xa(i+m)-x
                t = (xa(i)-x)*d(i)/h ! h never zero, since tested in initialization loop
                dd = t-c(i+1)
                if( dd == 0. ) stop 'failure in ratint; simple_math'
                dd = w/dd
                d(i) = c(i+1)*dd
                c(i) = t*dd
            end do
            if( 2*ns < n-m )then
                dy = c(ns+1)
            else
                dy = d(ns)
                ns = ns-1
            endif
            y = y+dy
        end do
        deallocate(c,d,stat=alloc_stat)
        call alloc_errchk("In: ratint; simple_math", alloc_stat)
    end subroutine ratint

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

    subroutine parabl(z,xsh,ysh,peakv)
        real, intent(in) :: z(3,3)
        real(dp) :: dz(3,3)
        real(dp) :: c1,c2,c3,c4,c5,c6,denom,dpeakv
        real,intent(out) :: xsh, ysh, peakv
        dz = dble(z)
          c1 = (26.*dz(1,1)-dz(1,2)+2*dz(1,3)-dz(2,1)-19.*dz(2,2)-7.*dz(2,3)+2.*dz(3,1)-7.*dz(3,2)+14.*dz(3,3))/9.
          c2 = (8.* dz(1,1)-8.*dz(1,2)+5.*dz(2,1)-8.*dz(2,2)+3.*dz(2,3)+2.*dz(3,1)-8.*dz(3,2)+6.*dz(3,3))/(-6.)
          c3 = (dz(1,1)-2.*dz(1,2)+dz(1,3)+dz(2,1)-2.*dz(2,2)+dz(2,3)+dz(3,1)-2.*dz(3,2)+dz(3,3))/6.
          c4 = (8.*dz(1,1)+5.*dz(1,2)+2.*dz(1,3)-8.*dz(2,1)-8.*dz(2,2)-8.*dz(2,3)+3.*dz(3,2)+6.*dz(3,3))/(-6.)
          c5 = (dz(1,1)-dz(1,3)-dz(3,1)+dz(3,3))/4.
          c6 = (dz(1,1)+dz(1,2)+dz(1,3)-2.*dz(2,1)-2.*dz(2,2)-2.*dz(2,3)+dz(3,1)+dz(3,2)+dz(3,3))/6.
        ! the peak coordinates of the paraboloid can now be evaluated as:
        ysh = 0.
        xsh = 0.
        denom = 4.*c3*c6-c5*c5
        if(abs(denom)<TINY) return
        ysh   = real((c4*c5-2.d0*c2*c6)/denom-2.d0)
        xsh   = real((c2*c5-2.*c4*c3)/denom-2.d0)
        dpeakv = 4.*c1*c3*c6-c1*c5*c5-c2*c2*c6+c2*c4*c5-c4*c4*c3
        peakv = real(dpeakv/denom)
        ! limit interplation to +/- 1. range
        if(ysh .lt. -1.) ysh = -1.
        if(ysh .gt.  1.) ysh = 1.
        if(xsh .lt. -1.) xsh = -1.
        if(xsh .gt.  1.) xsh = 1.
    end subroutine parabl

    ! SORTING

    function peakfinder( vals, npeaks ) result( peakpos )
        real,    intent(in)  :: vals(:)
        integer, intent(in)  :: npeaks
        integer, allocatable :: peakpos(:)
        logical, allocatable :: mask(:)
        integer :: n, ipeak, loc(1)
        n = size(vals)
        allocate(peakpos(npeaks), mask(n),stat=alloc_stat)
        if(alloc_stat /= 0) call alloc_errchk("In: peakfinder; simple_math", alloc_stat)
        mask = .true.
        do ipeak=1,npeaks
            loc = maxloc(vals, mask=mask)
            peakpos(ipeak) = loc(1)
            mask(loc(1)) = .false.
        end do
    end function peakfinder

    !>   rheapsort from numerical recepies (largest last)
    subroutine hpsort_1( n, rarr, iarr )
        integer, intent(in)    :: n
        real, intent(inout)    :: rarr(n)
        integer, intent(inout) :: iarr(n)
        integer                :: i, ir, j, l, ia
        real                   :: ra
        if( n < 2) return
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
    subroutine hpsort_2( n, iarr )
        integer, intent(in)    :: n
        integer, intent(inout) :: iarr(n)
        integer                :: i, ir, j, l, ra
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
    subroutine hpsort_3( n, iarr, p1_lt_p2 )
        integer, intent(in)    :: n
        integer, intent(inout) :: iarr(n)
        interface
            function p1_lt_p2( p1, p2 ) result( val )
                integer, intent(in) :: p1, p2
                logical :: val
            end function p1_lt_p2
        end interface
        integer                :: i, ir, j, l, ra
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
    subroutine hpsort_4( n, rarr )
        integer, intent(in) :: n
        real, intent(inout) :: rarr(n)
        integer             :: i, ir, j, l
        real                :: ra
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
    subroutine hpsort_5( n, rarr, rarr2 )
        integer, intent(in)    :: n
        real, intent(inout)    :: rarr(n), rarr2(n)
        integer                :: i, ir, j, l
        real                   :: ra, ra2
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

    !>   selecting the size(rheap) largest
    subroutine hpsel_1( rarr, rheap )
        real, intent(inout) :: rarr(:), rheap(:)
        integer :: i, j, k, m, n
        real    :: swap
        m = size(rheap)
        n = size(rarr)
        if(m.gt.n/2.or.m.lt.1) write(*,*) 'WARNING! probable misuse of hspel_1; simple_math'
        do i=1,m
            rheap(i) = rarr(i)
        end do
        call hpsort_4(m,rheap)
        do i=m+1,n
            if(rarr(i).gt.rheap(1))then
                rheap(1) = rarr(i)
                j=1
                do
                    k=2*j
                    if(k.gt.m) exit
                    if(k.ne.m)then
                        if(rheap(k).gt.rheap(k+1)) k=k+1
                    endif
                    if(rheap(j).le.rheap(k)) exit
                    swap = rheap(k)
                    rheap(k) = rheap(j)
                    rheap(j) = swap
                    j=k
                end do
            endif
        end do 
    end subroutine hpsel_1

    !>   selecting the size(rheap) largest
    subroutine hpsel_2( rarr, rheap, iheap )
        real, intent(inout)    :: rarr(:), rheap(:)
        integer, intent(inout) :: iheap(:)
        integer :: i, j, k, m, n, mi, iswap
        real    :: rswap
        m  = size(rheap)
        n  = size(rarr)
        mi = size(iheap)
        if( m /= mi ) stop 'nonconformable sizes of rheap & iheap; hpsel_2; simple_math'
        if(m.gt.n/2.or.m.lt.1) write(*,*) 'WARNING! probable misuse of hpsel_2; simple_math'
        do i=1,m
            rheap(i) = rarr(i)
            iheap(i) = i
        end do
        call hpsort_1(m,rheap,iheap)
        do i=m+1,n
            if(rarr(i).gt.rheap(1))then
                rheap(1) = rarr(i)
                iheap(1) = i
                j=1
                do
                    k=2*j
                    if(k.gt.m) exit
                    if(k.ne.m)then
                        if(rheap(k).gt.rheap(k+1)) k=k+1
                    endif
                    if(rheap(j).le.rheap(k)) exit
                    ! swap reals
                    rswap    = rheap(k)
                    rheap(k) = rheap(j)
                    rheap(j) = rswap
                    ! swap ints
                    iswap    = iheap(k)
                    iheap(k) = iheap(j)
                    iheap(j) = iswap
                    j=k
                end do
            endif
        end do
    end subroutine hpsel_2

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
        real, allocatable :: copy(:)
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
        !allocate( copy(n), stat=alloc_stat )
        !if(alloc_stat /= 0) call alloc_errchk('In: median, module: simple_math', alloc_stat)
        copy = arr  ! compiler will realloc here anyway
        if( pos1 == pos2 )then
            val  = selec(pos1,n,copy)
        else
            val1 = selec(pos1,n,copy)
            val2 = selec(pos2,n,copy)
            val  = (val1+val2)/2.
        endif
        deallocate(copy)
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

    !>   for selecting kth largest, array is modified
    real function selec_1(k,n,arr)
      integer k,n
      real arr(n)
      integer i,ir,j,l,mid
      real a,temp
      l = 1
      ir = n
 2    if (ir-l.le.1) then
         if (ir-1.eq.1) then
            if (arr(ir).lt.arr(l)) then
               temp = arr(l)
               arr(l) = arr(ir)
               arr(ir) = temp
            endif
         endif
         selec_1 = arr(k)
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
 3       continue
         i = i+1
         if (arr(i).lt.a) goto 3
 4       continue
         j = j-1
         if (arr(j).gt.a) goto 4
         if (j.lt.i) goto 5
         temp = arr(i)
         arr(i) = arr(j)
         arr(j) = temp
         goto 3
 5       arr(l+1) = arr(j)
         arr(j) = a
         if (j.ge.k) ir = j-1
         if (j.le.k) l = i
      endif
      goto 2
    end function selec_1

    !>   selecting kth largest, array is modified
    integer function selec_2(k,n,arr)
      integer :: k,n
      integer :: arr(n)
      integer :: i,ir,j,l,mid
      integer ::  a,temp
      l = 1
      ir = n
 2    if (ir-l.le.1) then
         if (ir-1.eq.1) then
            if (arr(ir).lt.arr(l)) then
               temp = arr(l)
               arr(l) = arr(ir)
               arr(ir) = temp
            endif
         endif
         selec_2 = arr(k)
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
 3       continue
         i = i+1
         if (arr(i).lt.a) goto 3
 4       continue
         j = j-1
         if (arr(j).gt.a) goto 4
         if (j.lt.i) goto 5
         temp = arr(i)
         arr(i) = arr(j)
         arr(j) = temp
         goto 3
 5       arr(l+1) = arr(j)
         arr(j) = a
         if (j.ge.k) ir = j-1
         if (j.le.k) l = i
      endif
      goto 2
    end function selec_2

end module simple_math
