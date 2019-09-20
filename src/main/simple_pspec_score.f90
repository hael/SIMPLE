module simple_pspec_score
include 'simple_lib.f08'
use simple_parameters, only: params_glob
use simple_image,      only: image
use simple_math,       only: deg2rad
use simple_stat
implicit none

private
public :: pspec_score
#include "simple_local_flags.inc"

type pspec_score
    type(image)         :: micrograph, pspec
    real                :: icescore
    real                :: graphenescore
    integer             :: grapheneindex
    real                :: averageheight
    real                :: avgofheights
    real, allocatable   :: rmat(:,:,:) 
    real, allocatable   :: pspecpolar(:,:) ! r, theta
    real, allocatable   :: plot1d(:)
    
contains
    procedure :: pspec2polar
    procedure :: score
    procedure :: ice_score
    procedure :: graphene_index
    procedure :: polar1d
    procedure :: kill
end type pspec_score

contains

    subroutine pspec2polar(self, ctfvars)
        class(pspec_score),     intent(inout)   :: self
        type(ctfparams),        intent(in)      :: ctfvars
        integer                                 :: x, y, xcen, ycen, radius
        real                                    :: angle, df
        integer                                 :: r, theta
        
        if(.NOT. allocated(self%pspecpolar)) then
           allocate(self%pspecpolar(SIZE(self%rmat,1) / 2, 720))
        end if
        
        do r = 1, SIZE(self%pspecpolar, 1)
            do theta = 1, SIZE(self%pspecpolar, 2)
            df = (ctfvars%dfx + ctfvars%dfy + cos(2.0 * deg2rad_sp(real(theta)/4 - ctfvars%angast)) * (ctfvars%dfx - ctfvars%dfy))
            x = floor(r*df*sin(deg2rad_sp(real(theta)/4)) / (ctfvars%dfx + ctfvars%dfy))
            y = floor(r*df*cos(deg2rad_sp(real(theta)/4)) / (ctfvars%dfx + ctfvars%dfy))
                self%pspecpolar(r, theta) = self%rmat(((SIZE(self%rmat,1) / 2) + x), ((SIZE(self%rmat,1) / 2) - y), 1)
            end do
        end do
        
    end subroutine pspec2polar
    
    subroutine graphene_index(self, ctfvars)
        class(pspec_score),     intent(inout)   :: self
        type(ctfparams),        intent(in)      :: ctfvars
        real, allocatable                       :: radial(:)
        real, allocatable                       :: radial60(:)
        real, allocatable                       :: zscore60(:)
        real, allocatable                       :: smooth(:)
        integer                                 :: r, theta
        real                                    :: freq, minvalue, maxvalue
        real                                    :: area,maxarea
        integer                                 :: addtheta, peakwidth
        real                                    :: addterm      ! for radial
        real                                    :: totalarea
        real                                    :: heightsum
        
        
        if(.NOT. allocated(radial)) then
            allocate(radial(SIZE(self%pspecpolar,2)))
        end if
        
        do theta = 1, SIZE(radial,1)
            radial(theta) = 0
        end do
        
        do r = 1, SIZE(self%pspecpolar,1)
            freq = r/(ctfvars%smpd*SIZE(self%rmat,1))
            if ( (freq .GT. 0.4 ) .AND. (freq .LT. 0.54 )) then !(theoretical peak @ 2.13A (0.47) did have between 0.46 and 0.48 but increased due to astig
                do theta = 1, SIZE(radial,1)
                
                    ! average of 3x3 square
                    addterm = self%pspecpolar(r, theta) + self%pspecpolar(r-1, theta) + self%pspecpolar(r+1, theta)
                    if (theta .eq. 1) then
                        addterm = addterm + self%pspecpolar(r, SIZE(self%pspecpolar,2)) + self%pspecpolar(r-1, SIZE(self%pspecpolar,2)) + self%pspecpolar(r+1, SIZE(self%pspecpolar,2))
                    else
                        addterm = self%pspecpolar(r, theta-1) + self%pspecpolar(r-1, theta-1) + self%pspecpolar(r+1, theta-1)
                    end if
                    if (theta .eq. SIZE(self%pspecpolar,2)) then
                        addterm = addterm + self%pspecpolar(r, 1) + self%pspecpolar(r-1, 1) + self%pspecpolar(r+1, 1)
                    else
                        addterm = addterm + self%pspecpolar(r, theta+1) + self%pspecpolar(r-1, theta+1) + self%pspecpolar(r+1, theta+1)
                    end if
                    addterm = addterm/9
                    
                    
                    ! highest value
                    if(addterm .gt. radial(theta)) then
                        radial(theta) = addterm
                    end if
                    
                end do 
            end if
        end do
        
        
        if(.NOT. allocated(radial60)) then
            allocate(radial60(SIZE(self%pspecpolar,2) / 3))
        end if
        
        if(.NOT. allocated(zscore60)) then
            allocate(zscore60(SIZE(self%pspecpolar,2) / 3))
        end if
        
        if(.NOT. allocated(smooth)) then
            allocate(smooth(SIZE(self%pspecpolar,2) / 3))
        end if
        
        
        do theta = 1, SIZE(radial60,1)
            radial60(theta) = 0
        end do
        
        ! sum to 60 degrees
        do theta = 1, SIZE(radial,1)
            radial60(MODULO(theta-1,SIZE(radial60,1)) + 1) = radial60(MODULO(theta-1,SIZE(radial60,1)) + 1) + radial(theta)
        end do
        
        !zscores of sum
        zscore60 = z_scores(radial60)
        
        minvalue = 0
        maxvalue = 0
        
        do theta = 1, SIZE(zscore60,1)
            ! sum over box width 3
            smooth(theta) = (zscore60(theta) + zscore60(MODULO(theta,SIZE(zscore60,1)) + 1) + zscore60(MODULO(theta+1,SIZE(zscore60,1)) + 1))/3
            if (smooth(theta) .lt. minvalue) then
                minvalue = smooth(theta)
            end if
            if (smooth(theta) .gt. maxvalue) then
                maxvalue = smooth(theta)
            end if
        end do
        
        if (maxvalue .gt. (-2)*minvalue) then
            print*, "GRAPHENE!"
        else
            print*, "NO GRAPHENE!"
            if(allocated(radial))deallocate(radial)
            if(allocated(radial60))deallocate(radial60)
            if(allocated(zscore60))deallocate(zscore60)
            if(allocated(smooth))deallocate(smooth)
            self%grapheneindex = 0
            self%graphenescore = 100
            return      ! exit subroutine
        end if
        
        addtheta = 0
        maxarea = 0
        
        do while (smooth(1 + addtheta) .gt. 0)
            ! don't start in the middle of a peak
            addtheta = addtheta + 1
        end do
        
        ! find biggest peak
        do theta = 1, SIZE(smooth)
            if (smooth(MODULO(theta-1+addtheta,SIZE(zscore60,1)) + 1) .lt. 0) then
                area = 0
            else
                area = area + smooth(MODULO(theta-1+addtheta,SIZE(zscore60,1)) + 1)
                if ( (area .gt. maxarea) ) then
                    maxarea = area
                end if
            end if
        end do
        
        self%grapheneindex = 0      ! peak count
        self%graphenescore = 0
        heightsum = 0
        totalarea = 0
        
        ! count peaks
        do theta = 1, SIZE(smooth)
            if (smooth(MODULO(theta-1+addtheta,SIZE(zscore60,1)) + 1) .lt. 0) then
                area = 0
                peakwidth = 0
            else
                area = area + smooth(MODULO(theta-1+addtheta,SIZE(zscore60,1)) + 1)
                peakwidth = peakwidth + 1
                if ( (smooth(MODULO(theta+addtheta,SIZE(zscore60,1)) + 1) .lt. 0 ) .and. (area .gt. maxarea/10) ) then
                    self%grapheneindex = self%grapheneindex + 1
                    self%graphenescore = self%graphenescore + peakwidth
                    totalarea = totalarea + area
                    heightsum = heightsum + (area/peakwidth)
                end if
            end if
        end do
        
        if( self%graphenescore .gt. 100) then
            self%graphenescore = 100
        end if
        
        self%averageheight = totalarea/self%graphenescore
        self%avgofheights = heightsum/self%grapheneindex  
        
        if(allocated(radial))deallocate(radial)
        if(allocated(radial60))deallocate(radial60)
        if(allocated(zscore60))deallocate(zscore60)
        if(allocated(smooth))deallocate(smooth)
        
    end subroutine graphene_index
    
    subroutine polar1d(self)
        class(pspec_score),     intent(inout)   :: self
        integer                                 :: r, theta
        
        if(.NOT. allocated(self%plot1d)) then
            allocate(self%plot1d(SIZE(self%pspecpolar,1)))
        end if
        
        do r = 1, SIZE(self%plot1d,1)
            self%plot1d(r) = 0 
            do theta = 1, SIZE(self%pspecpolar, 2)
                self%plot1d(r) = self%plot1d(r) + self%pspecpolar(r, theta)
            end do
            self%plot1d(r) = self%plot1d(r) / SIZE(self%pspecpolar, 2)
        end do
        
        
    end subroutine polar1d
    
    subroutine ice_score(self, ctfvars)
        class(pspec_score),     intent(inout)   :: self
        type(ctfparams),        intent(in)      :: ctfvars
        real                                    :: freq
        integer                                 :: r
        
        self%icescore = 0
        
        do r = 1, SIZE(self%plot1d,1)
            freq = r/(ctfvars%smpd*SIZE(self%rmat,1))
            if ( (freq .GT. 0.26 ) .AND. (freq .LT. 0.28 )) then
                self%icescore = self%icescore + self%plot1d(r)
            end if
        end do
    end subroutine ice_score
     
    subroutine score( self, ctfvars, moviename_forctf )
        class(pspec_score),     intent(inout)   :: self
        type(ctfparams),        intent(in)      :: ctfvars
        character(len=*),       intent(in)      :: moviename_forctf
        integer                                 :: ldim(3)
        integer                                 :: nframes
        
        call find_ldim_nptcls(trim(adjustl(moviename_forctf)), ldim, nframes)
        call self%micrograph%new(ldim, ctfvars%smpd)
        call self%micrograph%read(trim(adjustl(moviename_forctf)), 1)
        call self%micrograph%bp(real(params_glob%pspecsz) * ctfvars%smpd, 0.0)
        call self%pspec%new([params_glob%pspecsz, params_glob%pspecsz, 1], ctfvars%smpd)
        call self%micrograph%mic2spec(params_glob%pspecsz, 'sqrt', params_glob%hp, self%pspec)
        call self%pspec%write('pspec.mrc', 1, .true.)
        self%rmat = self%pspec%get_rmat()
        
        call self%pspec2polar(ctfvars)
        call self%polar1d()
        call self%ice_score(ctfvars)
        call self%graphene_index(ctfvars)
        
        call self%micrograph%kill
        call self%pspec%kill
        if(allocated(self%rmat))deallocate(self%rmat)
        if(allocated(self%pspecpolar))deallocate(self%pspecpolar)
        if(allocated(self%plot1d))deallocate(self%plot1d)
    end subroutine score
    
    subroutine kill( self )
        class(pspec_score), intent(inout) :: self
        call self%micrograph%kill
        call self%pspec%kill
        if(allocated(self%rmat))deallocate(self%rmat)
        if(allocated(self%pspecpolar))deallocate(self%pspecpolar)
        if(allocated(self%plot1d))deallocate(self%plot1d)
    end subroutine kill
    
end module simple_pspec_score
