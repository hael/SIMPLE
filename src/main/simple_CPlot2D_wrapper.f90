module CPlot2D_wrapper_module
    use, intrinsic :: ISO_C_Binding, only: C_int, C_ptr, C_NULL_ptr, C_char, C_double, C_bool
    implicit none

    type :: CDataPoint_type
        private
        type(C_ptr) :: object = C_NULL_ptr
    end type CDataPoint_type
    type :: CDataSet_type
        private
        type(C_ptr) :: object = C_NULL_ptr
    end type CDataSet_type
    type :: CPlot2D_type
        private
        type(C_ptr) :: object = C_NULL_ptr
    end type CPlot2D_type
    interface
        function C_CPlot2D__new(title) result(this) bind(C,name="CPlot2D__new")
            import
            type(C_ptr) :: this
            character(kind=c_char), dimension(*) :: title
        end function C_CPlot2D__new
        subroutine C_CPlot2D__SetXAxisSize(this, val) bind(C,name="CPlot2D__SetXAxisSize")
            import
            type(C_ptr), value :: this
            real(C_double), value :: val
        end subroutine C_CPlot2D__SetXAxisSize
        subroutine C_CPlot2D__SetYAxisSize(this, val) bind(C,name="CPlot2D__SetYAxisSize")
            import
            type(C_ptr), value :: this
            real(C_double), value :: val
        end subroutine C_CPlot2D__SetYAxisSize
        subroutine C_CPlot2D__SetDrawLegend(this, flag) bind(C,name="CPlot2D__SetDrawLegend")
            import
            type(C_ptr), value :: this
            logical(c_bool), value :: flag
        end subroutine C_CPlot2D__SetDrawLegend
        subroutine C_CPlot2D__SetFlipY(this, flag) bind(C,name="CPlot2D__SetFlipY")
            import
            type(C_ptr), value :: this
            logical(C_bool), value :: flag
        end subroutine C_CPlot2D__SetFlipY
        subroutine C_CPlot2D__AddDataSet(this, dataSet) bind(C,name="CPlot2D__AddDataSet")
            import
            type(C_ptr), value :: this
            type(C_ptr), value :: dataSet
        end subroutine C_CPlot2D__AddDataSet
        subroutine C_CPlot2D__SetXAxisTitle(this, font) bind(C,name="CPlot2D__SetXAxisTitle")
            import
            type(C_ptr), value :: this
            character(kind=c_char), dimension(*) :: font
        end subroutine C_CPlot2D__SetXAxisTitle
        subroutine C_CPlot2D__SetYAxisTitle(this, font) bind(C,name="CPlot2D__SetYAxisTitle")
            import
            type(C_ptr), value :: this
            character(kind=c_char), dimension(*) :: font
        end subroutine C_CPlot2D__SetYAxisTitle
        subroutine C_CPlot2D__OutputPostScriptPlot(this, fileName) bind(C,name="CPlot2D__OutputPostScriptPlot")
            import
            type(C_ptr), value :: this
            character(kind=c_char), dimension(*) :: fileName
        end subroutine C_CPlot2D__OutputPostScriptPlot
        subroutine C_CPlot2D__delete(this) bind(C,name="CPlot2D__delete")
            import
            type(C_ptr), value :: this
        end subroutine C_CPlot2D__delete
        function C_CDataSet__new() result(this) bind(C,name="CDataSet__new")
            import
            type(C_ptr) :: this
        end function C_CDataSet__new
        subroutine C_CPlot2D__SetDrawXAxisGridLines(this, flag) bind(C,name="CPlot2D__SetDrawXAxisGridLines")
            import
            type(C_ptr), value :: this
            logical(C_bool), value :: flag
        end subroutine C_CPlot2D__SetDrawXAxisGridLines
        subroutine C_CPlot2D__SetDrawYAxisGridLines(this, flag) bind(C,name="CPlot2D__SetDrawYAxisGridLines")
            import
            type(C_ptr), value :: this
            logical(C_bool), value :: flag
        end subroutine C_CPlot2D__SetDrawYAxisGridLines
        subroutine C_CDataSet__SetDrawMarker(this, flag) bind(C,name="CDataSet__SetDrawMarker")
            import
            type(C_ptr), value :: this
            logical(C_bool), value :: flag
        end subroutine C_CDataSet__SetDrawMarker
        subroutine C_CDataSet__SetDrawLine(this, flag) bind(C,name="CDataSet__SetDrawLine")
            import
            type(C_ptr), value :: this
            logical(C_bool), value :: flag
        end subroutine C_CDataSet__SetDrawLine
        subroutine C_CDataSet__SetMarkerSize(this, size) bind(C,name="CDataSet__SetMarkerSize")
            import
            type(C_ptr), value :: this
            real(C_double), value :: size
        end subroutine C_CDataSet__SetMarkerSize
        subroutine C_CDataSet__SetDatasetColor(this, r, g, b) bind(C,name="CDataSet__SetDatasetColor")
            import
            type(C_ptr), value :: this
            real(C_double), value :: r, g, b
        end subroutine C_CDataSet__SetDatasetColor
        function C_CDataPoint__new2(x, y) result(this) bind(C,name="CDataPoint__new2")
            import
            type(C_ptr) :: this
            real(C_double), value :: x, y
        end function C_CDataPoint__new2
        subroutine C_CDataPoint__delete(this) bind(C,name="CDataPoint__delete")
            import
            type(C_ptr), value :: this
        end subroutine C_CDataPoint__delete
        subroutine C_CDataSet__AddDataPoint(this, point) bind(C,name="CDataSet__AddDataPoint")
            import
            type(C_ptr), value :: this
            type(C_ptr), value :: point
        end subroutine C_CDataSet__AddDataPoint
        subroutine C_CDataSet__delete(this) bind(C,name="CDataSet__delete")
            import
            type(C_ptr), value :: this
        end subroutine C_CDataSet__delete
    end interface
contains
    subroutine CPlot2D__new(this, title)
        type(CPlot2D_type), intent(out) :: this
        character(kind=c_char), dimension(*), intent(in) :: title
        this%object = C_CPlot2D__new(title)
    end subroutine CPlot2D__new
    subroutine CPlot2D__SetXAxisSize(this, val)
        type(CPlot2D_type), intent(inout) :: this
        real(C_double), intent(in) :: val
        call C_CPlot2D__SetXAxisSize(this%object, val)
    end subroutine CPlot2D__SetXAxisSize
    subroutine CPlot2D__SetYAxisSize(this, val)
        type(CPlot2D_type), intent(inout) :: this
        real(C_double), intent(in) :: val
        call C_CPlot2D__SetYAxisSize(this%object, val)
    end subroutine CPlot2D__SetYAxisSize
    subroutine CPlot2D__SetDrawLegend(this, flag)
        type(CPlot2D_type), intent(inout) :: this
        logical(c_bool), intent(in) :: flag
        call C_CPlot2D__SetDrawLegend(this%object, flag)
    end subroutine CPlot2D__SetDrawLegend
    subroutine CPlot2D__SetFlipY(this, flag)
        type(CPlot2D_type), intent(inout) :: this
        logical(C_bool), intent(in) :: flag
        call C_CPlot2D__SetFlipY(this%object, flag)
    end subroutine CPlot2D__SetFlipY
    subroutine CPlot2D__AddDataSet(this, dataSet)
        type(CPlot2D_type), intent(inout) :: this
        type(CDataSet_type), intent(inout) :: dataSet
        call C_CPlot2D__AddDataSet(this%object, dataSet%object)
    end subroutine CPlot2D__AddDataSet
    subroutine CPlot2D__SetXAxisTitle(this, font)
        type(CPlot2D_type), intent(inout) :: this
        character(kind=c_char), dimension(*), intent(in) :: font
        call C_CPlot2D__SetXAxisTitle(this%object, font)
    end subroutine CPlot2D__SetXAxisTitle
    subroutine CPlot2D__SetYAxisTitle(this, font)
        type(CPlot2D_type), intent(inout) :: this
        character(kind=c_char), dimension(*), intent(in) :: font
        call C_CPlot2D__SetYAxisTitle(this%object, font)
    end subroutine CPlot2D__SetYAxisTitle
    subroutine CPlot2D__OutputPostScriptPlot(this, fileName)
        type(CPlot2D_type), intent(inout) :: this
        character(kind=c_char), dimension(*), intent(in) :: fileName
        call C_CPlot2D__OutputPostScriptPlot(this%object, fileName)
    end subroutine CPlot2D__OutputPostScriptPlot
    subroutine CPlot2D__delete(this)
        type(CPlot2D_type), intent(inout) :: this
        call C_CPlot2D__delete(this%object)
        this%object = C_NULL_PTR
    end subroutine CPlot2D__delete
    subroutine CDataSet__new(this)
        type(CDataSet_type), intent(out) :: this
        this%object = C_CDataSet__new()
    end subroutine CDataSet__new
    subroutine CPlot2D__SetDrawYAxisGridLines(this, flag)
        type(CPlot2D_type), intent(inout) :: this
        logical(C_bool), intent(in) :: flag
        call C_CPlot2D__SetDrawYAxisGridLines(this%object, flag)
    end subroutine CPlot2D__SetDrawYAxisGridLines
    subroutine CPlot2D__SetDrawXAxisGridLines(this, flag)
        type(CPlot2D_type), intent(inout) :: this
        logical(C_bool), intent(in) :: flag
        call C_CPlot2D__SetDrawXAxisGridLines(this%object, flag)
    end subroutine CPlot2D__SetDrawXAxisGridLines
    subroutine CDataSet__SetDrawMarker(this, flag)
        type(CDataSet_type), intent(inout) :: this
        logical(C_bool), intent(in) :: flag
        call C_CDataSet__SetDrawMarker(this%object, flag)
    end subroutine CDataSet__SetDrawMarker
    subroutine CDataSet__SetDrawLine(this, flag)
        type(CDataSet_type), intent(inout) :: this
        logical(C_bool), intent(in) :: flag
        call C_CDataSet__SetDrawLine(this%object, flag)
    end subroutine CDataSet__SetDrawLine
    subroutine CDataSet__SetMarkerSize(this, size)
        type(CDataSet_type), intent(inout) :: this
        real(C_double), intent(in) :: size
        call C_CDataSet__SetMarkerSize(this%object, size)
    end subroutine CDataSet__SetMarkerSize
    subroutine CDataSet__SetDatasetColor(this, r, g, b)
        type(CDataSet_type), intent(inout) :: this
        real(C_double), intent(in) :: r, g, b
        call C_CDataSet__SetDatasetColor(this%object, r, g, b)
    end subroutine CDataSet__SetDatasetColor
    subroutine CDataSet__delete(this)
        type(CDataSet_type), intent(inout) :: this
        call C_CDataSet__delete(this%object)
        this%object = C_NULL_PTR
    end subroutine CDataSet__delete
    subroutine CDataPoint__new2(x, y, this)
        type(CDataPoint_type), intent(out) :: this
        real(C_double), intent(in) :: x, y
        this%object = C_CDataPoint__new2(x, y)
    end subroutine CDataPoint__new2
    subroutine CDataSet__AddDataPoint(this, point)
        type(CDataSet_type), intent(inout) :: this
        type(CDataPoint_type), intent(inout) :: point
        call C_CDataSet__AddDataPoint(this%object, point%object)
    end subroutine CDataSet__AddDataPoint
    subroutine CDataPoint__delete(this)
        type(CDataPoint_type), intent(inout) :: this
        call C_CDataPoint__delete(this%object)
        this%object = C_NULL_PTR
    end subroutine CDataPoint__delete

    subroutine CDataSet_addpoint(this, x, y)
        type(CDataSet_type), intent(inout) :: this
        real,                  intent(in)    :: x,y
        type(CDataPoint_type) :: point
        call CDataPoint__new2( real(x, c_double), real(y, c_double), point)
        call CDataSet__AddDataPoint(this, point)
        call CDataPoint__delete(point)
    end subroutine CDataSet_addpoint

    subroutine test_CPlot2D
        include 'simple_lib.f08'
        ! use CPlot2D_wrapper_module
        type(str4arr)             :: title
        type(CPlot2D_type)        :: plot2D
        type(CDataSet_type)       :: dataSet
        character(len=LONGSTRLEN) :: ps2pdf_cmd, fname_pdf, ps2jpeg_cmd, fname_jpeg, fname_eps
        real    :: vals1(50,2), vals2(50,2), left, right, up, down, xdim,ydim
        integer :: l, i, iostat
        call seed_rnd
        do i =1,50
            vals1(i,:) = [ran3(), ran3()]
            vals2(i,:) = [ran3(), ran3()]
        enddo
        call CPlot2D__new(plot2D, 'some name'//C_NULL_CHAR)
        call CPlot2D__SetXAxisSize(plot2D, 400.0_c_double)
        call CPlot2D__SetYAxisSize(plot2D, 400.0_c_double)
        call CPlot2D__SetDrawXAxisGridLines(plot2D, C_TRUE)
        call CPlot2D__SetDrawYAxisGridLines(plot2D, C_TRUE)
        call CPlot2D__SetDrawLegend(plot2D, C_FALSE)
        call CPlot2D__SetFlipY(plot2D, C_TRUE)
        ! #0 dummy points for padding, arbitrary
        left  = min(minval(vals1(:,1)),minval(vals2(:,1)))
        right = max(maxval(vals1(:,1)),maxval(vals2(:,1)))
        down  = min(minval(vals1(:,2)),minval(vals2(:,2)))
        up    = max(maxval(vals1(:,2)),maxval(vals2(:,2)))
        xdim  = right-left
        ydim  = up-down
        left  = left  - 0.1*xdim; right = right + 0.1*xdim
        up    = up    + 0.1*ydim; down  = down  - 0.1*ydim
        call CDataSet__new(dataSet)
        call CDataSet__SetDrawMarker(dataSet, C_FALSE)
        call CDataSet__SetDatasetColor(dataSet, 1.0_c_double, 1.0_c_double, 1.0_c_double)
        call CDataSet_addpoint(dataSet, left, down)
        call CDataSet_addpoint(dataSet, right,down)
        call CDataSet_addpoint(dataSet, right,up)
        call CDataSet_addpoint(dataSet, left, up)
        call CPlot2D__AddDataSet(plot2D, dataset)
        call CDataSet__delete(dataset)
        ! #1, line, no marker, one dataset
        call CDataSet__new(dataSet)
        call CDataSet__SetDrawMarker(dataSet, C_FALSE) ! turn to C_TRUE for markers
        call CDataSet__SetDatasetColor(dataSet, 0.0_c_double, 0.0_c_double, 1.0_c_double)
        do i = 1,50
            call CDataSet_addpoint(dataSet, vals1(i,1), vals1(i,2))
        enddo
        call CPlot2D__AddDataSet(plot2D, dataset)
        call CDataSet__delete(dataset)
        ! #2, markers, no line, each point is a dataset
        do i = 1,50
            call CDataSet__new(dataSet)
            call CDataSet__SetDrawMarker(dataSet, C_TRUE)
            call CDataSet__SetMarkerSize(dataSet, real(max(5.,real(i)/2.),c_double))
            call CDataSet__SetDatasetColor(dataSet, 1.0_c_double, real(real(i)/100.+0.5,c_double), 0.0_c_double)
            call CDataSet_addpoint(dataSet, vals2(i,1), vals2(i,2))
            call CPlot2D__AddDataSet(plot2D, dataset)
            call CDataSet__delete(dataset)
        enddo
        ! legend
        title%str = 'X axis legend'//C_NULL_CHAR
        call CPlot2D__SetXAxisTitle(plot2D, title%str)
        title%str = 'Y axis legend'//C_NULL_CHAR
        call CPlot2D__SetYAxisTitle(plot2D, title%str)
        ! output
        fname_eps = 'test.eps'//C_NULL_CHAR
        call CPlot2D__OutputPostScriptPlot(plot2D, fname_eps)
        call CPlot2D__delete(plot2D)
        l = len_trim(fname_eps)
        fname_eps = fname_eps(:l-1) ! removing trailing C NULL character
        ! conversion to JPEG
        fname_jpeg  = trim(get_fbody(fname_eps,'eps'))//'.jpeg'
        ps2jpeg_cmd = 'gs -q -sDEVICE=jpeg -dJPEGQ=96 -dNOPAUSE -dBATCH -dSAFER -dDEVICEWIDTHPOINTS=512 -dDEVICEHEIGHTPOINTS=512 -sOutputFile='&
            //trim(fname_jpeg)//' '//trim(fname_eps)
        call exec_cmdline(trim(adjustl(ps2jpeg_cmd)), suppress_errors=.true., exitstat=iostat)
        ! conversion to PDF
        fname_pdf  = trim(get_fbody(fname_eps,'eps'))//'.pdf'
        ps2pdf_cmd = 'gs -q -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dDEVICEWIDTHPOINTS=512 -dDEVICEHEIGHTPOINTS=512 -sOutputFile='&
            //trim(fname_pdf)//' '//trim(fname_eps)
        call exec_cmdline(trim(adjustl(ps2pdf_cmd)), suppress_errors=.true., exitstat=iostat)
        ! remove eps
        if( iostat == 0 ) call del_file(fname_eps)
    end subroutine test_CPlot2D

end module CPlot2D_wrapper_module
