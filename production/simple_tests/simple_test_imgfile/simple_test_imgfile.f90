program simple_test_imgfile
    use simple_image,   only: image
    use simple_imgfile, only: imgfile
    use simple_imghead, only: imghead, SpiImgHead, test_imghead
    use simple_defs
    integer       :: ldim(3), i, j, cnt
    real          :: smpd=2., corr, corrs(20)
    type(image)   :: img, img_2
    type(image)   :: imgs(20)
    logical       :: ft=.true.
    
    ! SELF-CONSISTENCY TESTS
    allocate( endconv, source='NATIVE' )
    
    ! create a square
    ldim = [120,120,1]
    call img%new(ldim, smpd)
    call img%square(20)
    ! write stacks of 5 squares 
    do i=1,5
        if( ft ) call img%fwd_ft
        call img%write('squares_spider.spi',i)
        call img%write('squares_mrc.mrc',i) 
    end do
    ! convert the squares from SPIDER to MRC & vice versa
    do i=1,5
        call img%read('squares_spider.spi',i)
        if( ft ) call img%bwd_ft
        call img%write('squares_spider_converted.mrc',i)
        call img%read('squares_mrc.mrc',i)
        if( ft ) call img%bwd_ft
        call img%write('squares_mrc_converted.spi',i)
    end do
    ! test SPIDER vs. MRC & converted vs. nonconverted
    do i=1,20
        call imgs(i)%new(ldim, smpd)
    end do
    cnt = 0
    do i=1,5
        cnt = cnt+1
        call imgs(cnt)%read('squares_spider.spi',i)
        if( ft ) call imgs(cnt)%bwd_ft
    end do
    do i=1,5
        cnt = cnt+1
        call imgs(cnt)%read('squares_spider_converted.mrc',i)
        if( ft ) call imgs(cnt)%bwd_ft
    end do
    do i=1,5
        cnt = cnt+1
        call imgs(cnt)%read('squares_mrc.mrc',i)
        if( ft ) call imgs(cnt)%bwd_ft
    end do
    do i=1,5
        cnt = cnt+1
        call imgs(cnt)%read('squares_mrc_converted.spi',i)
        if( ft ) call imgs(cnt)%bwd_ft
    end do
    do i=1,19
        do j=i+1,20
            corr = imgs(i)%corr(imgs(j))
            if( corr < 0.99999 )then
                stop 'SPIDER vs. MRC & converted vs. nonconverted test failed'
            endif
        end do
    end do
    
    ! create a cube
    ldim = [120,120,120]
    call img%new(ldim, smpd)
    call img%square(20)
    ! write volume files
    do i=1,5
        if( ft ) call img%fwd_ft
        call img%write('cube_spider.spi')
        call img%write('cube_mrc.mrc') 
    end do
    ! convert the cubes from SPIDER to MRC & vice versa
    do i=1,5
        call img%read('cube_spider.spi')
        if( ft ) call img%bwd_ft
        call img%write('cube_spider_converted.mrc')
        call img%read('cube_mrc.mrc')
        if( ft ) call img%bwd_ft
        call img%write('cube_mrc_converted.spi')
    end do
    ! test SPIDER vs. MRC & converted vs. nonconverted
    do i=1,4
        call imgs(i)%new(ldim, smpd)
        call imgs(i)%read('cube_spider.spi')
        call imgs(i)%read('cube_spider_converted.mrc')
        call imgs(i)%read('cube_mrc.mrc')
        call imgs(i)%read('cube_mrc_converted.spi')
        if( ft ) call imgs(i)%bwd_ft
    end do
    do i=1,3
        do j=i+1,4
            corr = imgs(i)%corr(imgs(j))
            if( corr < 0.99999 )then
                stop 'SPIDER vs. MRC & converted vs. nonconverted test failed'
            endif
        end do
    end do
    stop
    
    ! EXTERNAL TESTS
!****(1) downloaded my favourite molecule (RNA olymerase II, PDB: 1WCM)
!****(2) converted to MRC using EMAN 1.9
!          >> pdb2mrc 1WCM.pdb 1WCM.mrc apix=1.5 center
!****(3) converted to SPIDER using EMAN 1.9
!         >> proc3d 1WCM.mrc 1WCM_from_eman.spi spidersingle
!****(4) coverted pdb2spider using SPIDER:
    
!     \__`O O'__/        SPIDER -- COPYRIGHT
!     ,__xXXXx___        HEALTH RESEARCH INC., ALBANY, NY.
!      __xXXXx__
!     /  /xxx\  \        VERSION:  UNIX  21.10 ISSUED: 11/12/2013
!       /     \          DATE:     11-MAY-2015    AT  03:40:48
    
!     .OPERATION: cp from pdb
!      cp from pdb
!      .PDB INPUT FILE: 1WCM.pdb
!      1WCM.pdb
!       OPENED (SF): 1WCM.pdb
!       Cell size:   146.8,  140.9,  156.8
!      .VOXEL SIZE [A]: 1.5
!           1.50
!      .CENTER? (Y/N): y
!          Y
!      .ATOMS OR TEMPERATURE? (A/T): A
!          A
!       Minimum size needed for this volume:    123   138   123
!      .SPIDER VOLUME SIZE; NX, NY & NZ: 144,144,144
!            144     144     144
!      .SPIDER OUTPUT FILE: 1WCM
! 
!****(5) generated an angle doc using VO EA
!     .OPERATION: vo ea
!      vo ea
!      .DELTA THETA: 30
!           30.0
!      .RANGE OF THETA (0,90):
!      .RANGE OF PHI (0,359.9):
!      .DOCUMENT FILE: angs
!      angs
!       11-MAY-2015 AT 03:48:37    OPENED NEW DOC FILE: angs
!       Total number of directions:  20
!     (6) Projected the two spider volumes (1WCM.spi & 1WCM_from_eman.spi) using PJ 3G
!     .OPERATION: pj 3g
!      pj 3g
!      .3-D INPUT FILE: 1WCM
!      1WCM
!       1WCM
!       (R3)   144   144   144 CREATED 11-MAY-2015 AT 03:41:28  O HEADER BYTES:   1152
!      .FILE NUMBERS OR SELECTION DOC. FILE NAME: 1-20
!          1-20
!      .ANGLES DOC FILE: angs
!      angs
!       OPENED EXISTING DOC FILE: angs
!       Number of keys recovered:      20
!       New Reverse FFTW3 Plan:     4341136304 (   288,   288,     1) Threads:  1
!       New Forward FFTW3 Plan:     4343201904 (   288,   288,   288) Threads:  1
!      .ENTER TEMPLATE FOR 2-D PROJECTION: projs@***
! Note: the volumes generated by SPIDER & EMAN not docked so cross-comparision not possible
    
    ! read the volumes, convert them, read back the converted files and compare with originals
    ldim = [144,144,144]
    call img%new(ldim, 1.5)
    do i=1,6
        imgs(i) = img
    end do
    call imgs(1)%read( '1WCM.mrc'                    )
    call imgs(2)%read( '1WCM.spi'                    )
    call imgs(3)%read( '1WCM_from_eman.spi'          )
    call imgs(1)%write('1WCM_converted.spi'          )
    call imgs(2)%write('1WCM_converted.mrc'          )
    call imgs(3)%write('1WCM_from_eman_converted.mrc')
    call imgs(4)%read( '1WCM_converted.spi'          )
    call imgs(5)%read( '1WCM_converted.mrc'          )
    call imgs(6)%read( '1WCM_from_eman_converted.mrc')
    corrs(1) = imgs(1)%corr(imgs(4))
    corrs(2) = imgs(2)%corr(imgs(5))
    corrs(3) = imgs(3)%corr(imgs(6))
    if( any(corrs(1:3) < 0.99999) )then
        stop 'External SPIDER vs. MRC & converted vs. nonconverted test failed'
    endif

    ! read the projections, convert them, read back the converted files and compare with originals
    ldim = [144,144,1]
    call img%new(ldim, 1.5)
    img_2 = img
    do i=1,20
        call img%read('projs.spi',i)
        call img%write('projs_converted.mrc',i)
        call img%read('projs_from_eman.spi',i)
        call img%write('projs_from_eman_converted.mrc',i)
    end do
    do i=1,20
        call img%read('projs.spi',i)
        call img_2%read('projs_converted.mrc',i)
        if( img%corr(img_2) < 0.99999 )then
            stop 'SPIDER vs. MRC & converted vs. nonconverted projection test failed'
        endif
        call img%read('projs_from_eman.spi',i)
        call img_2%read('projs_from_eman_converted.mrc',i)
        if( img%corr(img_2) < 0.99999 )then
            stop 'SPIDER vs. MRC & converted vs. nonconverted projection test failed'
        endif
    end do
    
    ! OTHER TESTS
    ! * checked that UCSF Chimera can read the volume files
    ! * checked that EMAN 1.9 can read all files produced
    ! FAILURE: SPIDER CANNOT READ THE STACKS, SEE BELOW
    
!     .OPERATION: fi
!      fi
!      .INPUT FILE: squares_spider@001
!      squares_spider@001
!       NON-NATIVE BYTE ORDERED SPIDER FILE
!       squares_spider@001
!       (S2)   120   120 (..      5) CREATED 19 AT   O HEADER BYTES:   1440
!      *** ERROR -- STACK LACKS IMAGE:        1
    
! SPIDER does, however, succeed in reading all the volumes & converting 
! all the MRC volumes to EMAN-readable SPIDER volumes using the command "cp from mrc"

! Trying to just copy the stack and see if SPIDER likes that

! do i=1,20
!     call img%read('projs.spi',i)
!     call img%write('projs_copy.spi',i)
! end do

! Nope, get same error *** ERROR -- STACK LACKS IMAGE:        1
        
end program
