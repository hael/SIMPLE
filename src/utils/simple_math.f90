! various mathematical subroutines and functions
module simple_math
use simple_defs
use simple_error, only: simple_exception
use simple_srch_sort_loc
use simple_is_check_assert
use simple_linalg
implicit none

interface otsu
    module procedure otsu_1, otsu_2, otsu_3
end interface otsu

interface pixels_dist
   module procedure pixels_dist_1, pixels_dist_2
end interface

interface nvoxfind
    module procedure nvoxfind_1, nvoxfind_2
end interface

interface vis_mat
    module procedure vis_2Dreal_mat, vis_2Dinteger_mat, vis_3Dreal_mat, vis_3Dinteger_mat
end interface vis_mat

interface mode
    module procedure mode_1, mode_2
end interface mode

logical, parameter, private :: warn=.false.

contains

    ! array operations

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

    pure subroutine bounds_from_mask3D( l_mask, lb, ub )
        logical, intent(in)    :: l_mask(:,:,:)
        integer, intent(inout) :: lb(3), ub(3)
        integer :: ldim(3)
        ldim(1) = size(l_mask, dim=1)
        ldim(2) = size(l_mask, dim=2)
        ldim(3) = size(l_mask, dim=3)
        lb(1) = 1
        do while( lb(1) <= ldim(1) / 2 )
            if( any(l_mask(lb(1),:,:)) ) exit
            lb(1) = lb(1) + 1
        end do
        lb(2) = 1
        do while( lb(2) <= ldim(2) / 2 )
            if( any(l_mask(:,lb(2),:)) ) exit
            lb(2) = lb(2) + 1
        end do
        lb(3) = 1
        do while( lb(3) <= ldim(3) / 2 )
            if( any(l_mask(:,:,lb(3))) ) exit
            lb(3) = lb(3) + 1
        end do
        ub(1) = ldim(1)
        do while( ub(1) >= ldim(1) / 2 )
            if( any(l_mask(ub(1),:,:)) ) exit
            ub(1) = ub(1) - 1
        end do
        ub(2) = ldim(2)
        do while( ub(2) >= ldim(2) / 2 )
            if( any(l_mask(:,ub(2),:)) ) exit
            ub(2) = ub(2) - 1
        end do
        ub(3) = ldim(3)
        do while( ub(3) >= ldim(3) / 2 )
            if( any(l_mask(:,:,ub(3))) ) exit
            ub(3) = ub(3) - 1
        end do
    end subroutine bounds_from_mask3D

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

    ! vector quantization

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
        allocate( mask(ndat), labels(ndat) )
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
                loc = minloc(abs(means-dat(i)))
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

    subroutine detect_peak_thres( n, n_ub, x, t )
        integer, intent(in)    :: n, n_ub
        real,    intent(in)    :: x(n)
        real,    intent(inout) :: t
        real, allocatable  :: arr(:)
        real    :: ts(2), smd, y
        integer :: narr, i
        ts(1) = -huge(y)
        do
            narr = count(x >= ts(1))
            if( narr < n_ub )then
                t = ts(1)
                return
            endif
            arr  = pack(x, mask=x >= ts(1))
            call otsu(narr, arr, ts(2))
            ts(1) = ts(2)
            deallocate(arr)
        end do
    end subroutine detect_peak_thres

    ! Source https://www.mathworks.com/help/stats/hierarchical-clustering.html#bq_679x-10
    subroutine hac_1d( vec, thresh, labels, centroids, populations )
        real,    intent(in)  :: vec(:)       ! input vector
        real,    intent(in)  :: thresh       ! threshold for class merging
        integer, intent(out) :: labels(:)    ! labels of the elements in vec
        real,    allocatable, intent(out) :: centroids(:)   ! centroids of the classes
        integer, allocatable, intent(out) :: populations(:) ! number of elements belonging to each class
        real,    allocatable :: mat(:,:)     ! matrix that contains the distances
        logical, allocatable :: mask(:), outliers(:)
        real    :: d
        integer :: N, i, j, index(1), loc1(1), loc2(1), cnt, ncls
        N = size(vec) ! number of data points
        ! 1) calc all the couples of distances, using Manhattan dist
        allocate(mat(N,N), source = 0.)
        do i = 1, N-1
            do j = i+1, N
                mat(i,j) = abs(vec(i)-vec(j))
                mat(j,i) = mat(i,j)
            enddo
        enddo
        ! 2) Generate binary clusters
        allocate(mask(N),source = .true. )
        allocate(outliers(N), source = .false.)
        ncls = 0
        do i = 1, N ! find the class for vec(i)
            if(mask(i)) then ! if it's not already been clustered
                mask(i) = .false.
                ! find the index of the couple
                d = minval(mat(i,:), mask)
                index(:) = minloc(mat(i,:), mask)
                ncls = ncls + 1
                ! assign labels
                labels(i) = ncls
                if(d <= thresh) then ! if it's not an outlier (it has a couple)
                    labels(index(1)) = labels(i)
                    mask(index(1)) = .false. ! index(1) has already been clustered
                else
                    outliers(i) = .true.
                endif
            endif
        enddo
        ! 3) Calculate centroids
        allocate(centroids(ncls), source = 0.)
        mask = .true. ! reset
        do i = 1, ncls
            ! find the first member of the class
            loc1(:) = minloc(abs(labels-i))
            if(.not. outliers(loc1(1))) then
                mask(loc1(1)) = .false.
                ! find the second member of the class
                loc2(:) = minloc(abs(labels-i), mask)
                mask(loc2(1)) = .false.
                centroids(i) = (vec(loc1(1)) + vec(loc2(1)))/2.
            else ! the class has just one element
                loc1(:) = minloc(abs(labels-i))
                centroids(i) = vec(loc1(1))
                mask(loc1(1)) = .false.
            endif
        enddo
        mask  = .true. ! reset
        ! 4) merge classes
        do i = 1, ncls-1
            do j = i+1, ncls
                if(abs(centroids(i)-centroids(j)) <= thresh) then ! merge classes
                    ! change label to class j
                    where(labels == j) labels = i
                endif
            enddo
        enddo
        ! 5) Reorder labels
        cnt = 0
        do i = 1, ncls
            if(any(labels== i)) then ! there is a class labelled i
                cnt = cnt + 1        ! increasing order
                where(labels == i) labels = cnt
            endif
        enddo
        ! 6) recalculate centroids
        deallocate(centroids)
        ncls = maxval(labels) ! the nr of classes is maxval(labels)
        allocate(centroids(ncls),   source = 0.)
        allocate(populations(ncls), source = 0 )
        mask = .true. ! reset
        do i = 1, ncls
            populations(i) = count(labels == i) ! counts the nr of elements in the class
            ! find all cnt member of the class and update the centroids
            do j = 1, populations(i)
                loc1(:) = minloc(abs(labels-i), mask)
                mask(loc1(1)) = .false. ! do not consider this member of the class anymore
                centroids(i) = centroids(i)+ vec(loc1(1))
            enddo
            centroids(i) = centroids(i)/real(populations(i))
        enddo
    end subroutine hac_1d

    ! image processing

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

    subroutine create_hist_vector( x, n, xhist, yhist )
        real,    intent(in)    :: x(:)         ! data
        integer, intent(in)    :: n            ! number of intervals
        real,    intent(inout) :: xhist(0:n-1) ! discretization of the values
        integer, intent(inout) :: yhist(n)     ! number of occurences
        integer :: i, j
        real    :: xmin, xmax, dx
        xmin  = minval(x)
        xmax  = maxval(x)
        xhist = 0.
        yhist = 0
        dx    = ( xmax + 1 - xmin )/n
        do i = 0, n-1
            xhist(i) = xmin + i*dx
        end do
        do i = 1, size(x)
            j = nint( (x(i) - xmin)/dx )+1
            if( j <= size(yhist) )then
                yhist(j) = yhist(j) + 1
            else
                yhist(size(yhist)) = yhist(size(yhist)) + 1
            endif
        end do
    end subroutine create_hist_vector

    ! Otsu's method
    ! the algorithm assumes that x contains two classes of numbers following bi-modal histogram (foreground and background)
    ! it then calculates the optimum threshold separating the two classes so that their combined spread (intra-class variance) is minimal
    subroutine otsu_1( n, x, thresh )
        integer, intent(in)    :: n
        real,    intent(inout) :: x(:)
        real,    intent(inout) :: thresh ! threshold for binarisation, in the old value range
        integer, parameter  :: MIN_VAL = 0, MAX_VAL = 255, NHIST = MAX_VAL-MIN_VAL+1
        integer :: i, T, yhist(NHIST)
        real    :: q1, q2, m1, m2, sigma1, sigma2, sigma, sigma_next, sum1, sum2
        real    :: p(MIN_VAL:MAX_VAL), xhist(0:NHIST-1), x_copy(n), sc, old_range(2) ! scale factor and old pixel range
        x_copy = x
        ! scale vecor
        old_range(1) = minval(x_copy)
        old_range(2) = maxval(x_copy)
        sc = (MAX_VAL - MIN_VAL)/(old_range(2) - old_range(1))
        x_copy(:) = sc * x_copy(:) + MIN_VAL - sc * old_range(1)
        p = 0.
        call create_hist_vector(x_copy, NHIST, xhist, yhist)
        do i = MIN_VAL, MAX_VAL
            p(i) = yhist(i+1)
        enddo
        p    = p / size(x_copy, dim = 1) ! normalise, it's a probability
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
        ! rescale in the old range
        thresh = thresh / sc + old_range(1)
    end subroutine otsu_1

    ! Otsu's method, see above
    subroutine otsu_2(n, x, x_out, thresh)
        integer,        intent(in)    :: n
        real,           intent(inout) :: x(n)
        real,           intent(inout) :: x_out(n)
        real, optional, intent(inout) :: thresh ! selected threshold for binarisation
        real :: threshold ! threshold for binarisation
        call otsu_1(n, x, threshold)
        where(x > threshold)
            x_out = 1.
        elsewhere
            x_out = 0.
        endwhere
        if(present(thresh)) thresh = threshold
    end subroutine otsu_2

    ! Otsu's method, see above
    subroutine otsu_3(n, x, mask)
        integer, intent(in)    :: n
        real,    intent(inout) :: x(n)
        logical, intent(inout) :: mask(n)
        real :: threshold ! threshold for binarisation
        call otsu_1(n, x, threshold)
        where(x > threshold)
            mask = .true.
        elsewhere
            mask = .false.
        endwhere
    end subroutine otsu_3

    function calc_score_thres( n, scores, npeaks ) result( t )
        integer, intent(in) :: n, npeaks
        real,    intent(in) :: scores(n)
        real    :: t, scores_sorted(n)
        integer :: i, cnt
        scores_sorted = scores
        call hpsort(scores_sorted) ! largest last
        if( npeaks >= n )then
            t = scores_sorted(n)
        else 
            cnt = 0
            do i = n, 1, -1
                t   = scores_sorted(i)
                cnt = cnt + 1
                if( cnt == npeaks ) exit
            end do
        endif
    end function calc_score_thres

    !>   calculates the euclidean distance between one pixel and a list of other pixels.
    ! if which == 'max' then distance is the maximum value of the distance between
    !              the selected pixel and all the others
    ! if which == 'min' then distance is the minimum value of the distance between
    !              the selected pixel and all the others
    ! if which == 'sum' then distance is the sum of the distances between the
    !              selected pixel and all the others.
    function pixels_dist_1( px, vec, which, mask, location) result( dist )
        integer,           intent(in)     :: px(3)
        integer,           intent(in)     :: vec(:,:)
        character(len=*),  intent(in)     :: which
        logical,           intent(inout)  :: mask(:)
        integer, optional, intent(out)    :: location(1)
        real    :: dist
        integer :: i
        if(size(mask,1) .ne. size(vec, dim = 2)) write(logfhandle,*)'Error! Incompatible sizes mask and input vector; pixels_dist_1'
        if((any(mask .eqv. .false.)) .and. which .eq. 'sum') write(logfhandle,*) 'Attention! Not considering mask for sum; pixels_dist_1'
        select case(which)
        case('max')
            dist =  maxval(sqrt((real(px(1)-vec(1,:)))**2+(real(px(2)-vec(2,:)))**2+(real(px(3)-vec(3,:)))**2), mask)
            if(present(location)) location =  maxloc(sqrt((real(px(1)-vec(1,:)))**2+(real(px(2)-vec(2,:)))**2+(real(px(3)-vec(3,:)))**2), mask)
        case('min')
            !to calculation of the 'min' excluding the pixel itself, otherwise it d always be 0
            do i = 1, size(vec, dim = 2)
                if( px(1)==vec(1,i) .and. px(2)==vec(2,i) .and. px(3)==vec(3,i) )then
                    if(which .ne. 'sum') mask(i) = .false.
                endif
            enddo
            dist =  minval(sqrt((real(px(1)-vec(1,:)))**2+(real(px(2)-vec(2,:)))**2+(real(px(3)-vec(3,:)))**2), mask)
            if(present(location)) location =  minloc(sqrt((real(px(1)-vec(1,:)))**2+(real(px(2)-vec(2,:)))**2+(real(px(3)-vec(3,:)))**2), mask)
        case('sum')
            if(present(location))   write(logfhandle,*)'Error! Unsupported location parameter with sum mode; pixels_dist_1'
            dist =  sum   (sqrt((real(px(1)-vec(1,:)))**2+(real(px(2)-vec(2,:)))**2+(real(px(3)-vec(3,:)))**2))
        case DEFAULT
            write(logfhandle,*) 'Pixels_dist kind: ', trim(which)
            write(logfhandle,*)'Error! Unsupported pixels_dist kind; pixels_dist_1'
        end select
    end function pixels_dist_1

    function pixels_dist_2( px, vec, which, mask, location, keep_zero) result( dist )
        real,              intent(in)     :: px(3)
        real,              intent(in)     :: vec(:,:)
        character(len=*),  intent(in)     :: which
        logical,           intent(inout)  :: mask(:)
        integer, optional, intent(out)    :: location(1)
        logical, optional, intent(in)     :: keep_zero   ! when calculating min dist, let 0 be the output
        real    :: dist
        integer :: i
        logical :: kkeep_zero
        if( size(mask) /= size(vec, dim=2) )then
            write(logfhandle,*)'Error! Incompatible sizes mask and input vector; pixels_dist_2'
            print *, 'size(mask)      : ', size(mask)
            print *, 'size(vec, dim=2): ', size(vec, dim=2)
            stop
        endif
        if( any(mask .eqv. .false.) .and. which .eq. 'sum' ) write(logfhandle,*) 'Attention! Not considering mask for sum; pixels_dist_2'
        select case(which)
            case('max')
                dist =  maxval(sqrt((px(1)-vec(1,:))**2+(px(2)-vec(2,:))**2+(px(3)-vec(3,:))**2), mask)
                if(present(location)) location = maxloc(sqrt((px(1)-vec(1,:))**2+(px(2)-vec(2,:))**2+(px(3)-vec(3,:))**2), mask)
            case('min')
                kkeep_zero = .false. ! default
                if(present(keep_zero)) kkeep_zero = keep_zero
                if(.not. kkeep_zero ) then
                    ! to calculation of the 'min' excluding the pixel itself, otherwise it d always be 0
                    do i = 1, size(vec, dim = 2)
                        if(      abs(px(1)-vec(1,i)) < TINY .and. abs(px(2)-vec(2,i)) < TINY  &
                        &  .and. abs(px(3)-vec(3,i)) < TINY )then
                            mask(i) = .false.
                            exit
                        endif
                    enddo
                endif
                dist =  minval(sqrt((px(1)-vec(1,:))**2+(px(2)-vec(2,:))**2+(px(3)-vec(3,:))**2), mask)
                if(present(location)) location = minloc(sqrt((px(1)-vec(1,:))**2+(px(2)-vec(2,:))**2+(px(3)-vec(3,:))**2), mask)
            case('sum')
                if(present(location))   write(logfhandle,*)'Error! Unsupported location parameter with sum mode; pixels_dist_2'
                dist =  sum(sqrt((px(1)-vec(1,:))**2+(px(2)-vec(2,:))**2+(px(3)-vec(3,:))**2))
            case DEFAULT
                write(logfhandle,*) 'Pixels_dist kind: ', trim(which)
                write(logfhandle,*)'Error! Unsupported pixels_dist kind; pixels_dist_2'
        end select
    end function pixels_dist_2

    ! This subroutine stores in pos the indexes corresponding to
    ! the pixels with value > 0 in the binary matrix imat.
    subroutine get_pixel_pos(imat, pos)
        integer,              intent(in)  :: imat(:,:,:)
        integer, allocatable, intent(out) :: pos(:,:)
        integer ::  i, j, k, cnt
        integer :: s(3)
        if(allocated(pos)) deallocate(pos)
        s = shape(imat)
        allocate(pos(3,count(imat(:s(1),:s(2),:s(3)) > 0.5)), source = 0)
        cnt = 0
        do i = 1, s(1)
              do j = 1, s(2)
                  do k = 1, s(3)
                      if(imat(i,j,k) > 0.5) then
                          cnt = cnt + 1
                          pos(:3,cnt) = [i,j,k]
                      endif
                  enddo
              enddo
          enddo
      end subroutine get_pixel_pos

    ! helpers

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

    !>   for rounding to closest odd
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

    !>   4 shifting variables
    subroutine shft(a,b,c,d)
        real, intent(out)   :: a   !< new pos 1
        real, intent(inout) :: b,c !< pos 2 and 3
        real, intent(in)    :: d   !< pos4
        a = b
        b = c
        c = d
    end subroutine shft

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

    pure function nvoxfind_1( smpd, mwkda ) result( nvox )
        real, intent(in) :: smpd             !< sampling distance
        real, intent(in) :: mwkda            !< molecular weight
        integer          :: nvox             !< nr of voxels
        double precision , parameter :: prot_d = 1.43d0            ! g/cm**3
        double precision , parameter :: one_da = 1.66053892173e-27 ! kg/Da
        nvox = nint((mwkda*one_da*1e30) / (prot_d * (smpd**3.)))
    end function nvoxfind_1

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

    pure function mwkdafind( smpd, nvox ) result( mwkda )
        real,    intent(in) :: smpd !< sampling distance
        integer, intent(in) :: nvox !< # voxels
        double precision , parameter :: prot_d = 1.43d0            ! g/cm**3
        double precision , parameter :: one_da = 1.66053892173e-27 ! kg/Da
        real :: pixv, mwkda
        pixv = smpd * smpd * smpd
        mwkda = (prot_d * pixv * real(nvox)) / (one_da * 1e30)
    end function mwkdafind

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

    ! mathematical functions

    !>   Compute the cross product of 2 3D real vectors
    function cross( a, b ) result( c )
        real, intent(in) :: a(3),b(3)
        real :: c(3)
        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
    end function cross

    !>   given two points on the sphere, find N sample points on the great circle formed by two points
    function great_circle_samples( p1, p2, N ) result( points )
        real,    intent(in) :: p1(3), p2(3)
        integer, intent(in) :: N
        real    :: points(3, N), dt, t, u(3), w(3), v(3)
        integer :: k
        u  = p1/norm2(p1)
        w  = cross(p1, p2)
        v  = cross(u, w/norm2(w))
        dt = 2.*pi/N
        do k = 1, N
            t = (k - 1)*dt
            points(:, k) = u * cos(t) + v * sin(t)
        enddo
    end function great_circle_samples

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

    real function gaussian1D( x, avg, sigma_sq )
        real, intent(in) ::   x, avg, sigma_sq
        gaussian1D = (1./(sqrt(2.*PI*sigma_sq)))*exp((1./(2.*sigma_sq))*(x-avg)**2)
    end function gaussian1D

    real function gaussian2D( center_coords, x, y, xsigma, ysigma )
        real, intent(in) :: center_coords(2), x, y, xsigma, ysigma
        gaussian2D = exp( -((x - center_coords(1))**2.0 / (2.0 * xsigma * xsigma) +&
                            (y - center_coords(2))**2.0 / (2.0 * ysigma * ysigma)) )
    end function gaussian2D

    real function gaussian3D( center_coords, x, y, z, xsigma, ysigma, zsigma )
        real, intent(in) :: center_coords(3), x, y, z, xsigma, ysigma, zsigma
        gaussian3D = exp( -((x - center_coords(1))**2.0 / (2.0 * xsigma * xsigma) +&
                           &(y - center_coords(2))**2.0 / (2.0 * ysigma * ysigma) +&
                           &(z - center_coords(3))**2.0 / (2.0 * zsigma * zsigma)) )
    end function gaussian3D

    ! Savitzky-Golay filter

    ! From numerical recipes ch14.8
    ! y: input vector (size n) is returned filtered.
    ! "To summarize: Within limits, Savitzky-Golay filtering does manage to provide smoothing
    ! without loss of resolution. It does this by assuming that relatively distant data points have
    ! some significant redundancy that can be used to reduce the level of noise. The specific nature
    ! of the assumed redundancy is that the underlying function should be locally well-fitted by a
    ! polynomial. When this is true, as it is for smooth line profiles not too much narrower than
    ! the filter width, then the performance of Savitzky-Golay filters can be spectacular. When it
    ! is not true, then these filters have no compelling advantage over other classes of smoothing
    ! filter coefficients."
    subroutine SavitzkyGolay_filter( n, y )
        integer, intent(in)    :: n
        real,    intent(inout) :: y(n)
        integer, parameter :: nl = 3  ! number of points to the left
        integer, parameter :: nr = 3  ! number of points to the right
        integer, parameter :: NP = 7 ! NP = nl+nr+1
        integer, parameter :: m  = 4  ! smoothing polynomial order
        real    :: ysave(n), c(NP)
        integer :: i,j, index(NP)
        ysave = y
        ! seek shift index for given case nl, nr, m (see savgol).
        index(1)=0
        j=3
        do i=2, nl+1
            index(i)=i-j
            j=j+2
        end do
        j=2
        do i=nl+2, nl+nr+1
            index(i)=i-j
            j=j+2
        end do
        ! calculate Savitzky-Golay filter coefficients
        call savgol(nl+nr+1,0)
        ! Apply filter to input data
        do i=1, n-nr
            y(i)=0.
            do j=1, nl+nr+1
                !skip left points that do not exist
                if (i+index(j).gt.0) y(i) = y(i)+c(j)*ysave(i+index(j))
            end do
        enddo

        contains

            subroutine savgol(np_here,ld)
                integer, intent(in)  :: np_here, ld
                integer, parameter   :: MMAX = 6
                integer :: d,icode,imj,ipj,j,k,kk,mm,indx(MMAX+1)
                real    :: fac,sum,a(MMAX+1,MMAX+1),b(MMAX+1)
                if(np_here.lt.nl+nr+1.or.nl.lt.0.or.nr.lt.0.or.ld.gt.m.or.m.gt.MMAX  &
                .or.nl+nr.lt.m) stop ' Bad args in savgol'
                do ipj=0,2*m        !Set up the normal equations of the desired leastsquares fit.
                    sum=0.
                    if(ipj.eq.0) sum=1.
                    do k=1,nr
                        sum=sum+float(k)**ipj
                    end do
                    do k=1,nl
                        sum=sum+float(-k)**ipj
                    end do
                    mm=min(ipj,2*m-ipj)
                    do imj=-mm,mm,2
                        a(1+(ipj+imj)/2,1+(ipj-imj)/2)=sum
                    end do
                end do
                call ludcmp(a,m+1,MMAX+1,indx,d,icode)    !Solve them: LU decomposition.
                do j=1,m+1
                    b(j)=0.
                end do
                b(ld+1)=1.      !Right-hand side vector is unit vector, depending on which derivative we want.
                call lubksb(a,m+1,MMAX+1,indx,b)   !Backsubstitute, giving one row of the inverse matrix.
                do kk=1,np_here                    !Zero the output array (it may be bigger than the number
                    c(kk)=0.                         !of coefficients).
                end do
                do k=-nl,nr                        !Each Savitzky-Golay coefficient is the dot product
                    sum=b(1)                         !of powers of an integer with the inverse matrix row.
                    fac=1.
                    do mm=1,m
                        fac=fac*k
                        sum=sum+b(mm+1)*fac
                    end do
                    kk=mod(np_here-k,np_here)+1                !Store in wrap-around order.
                    c(kk)=sum
                end do
            end subroutine savgol

    end subroutine SavitzkyGolay_filter

end module simple_math
