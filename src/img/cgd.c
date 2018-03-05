#define HAVE_LIBPNG 1
#define HAVE_LIBJPEG 1
#define CGD_IMAGE_SET_AA                 cgd_image_set_aa_
#define CGD_IMAGE_SET_AA_NB              cgd_image_set_aa_nb_
#define CGD_IMAGE_SET_CLIP               cgd_image_set_clip_
#define CGD_IMAGE_GET_CLIP               cgd_image_get_clip_
#define CGD_IMAGE_CREATE_TRUECOLOR       cgd_image_create_truecolor_
#define CGD_IMAGE_CREATE                 cgd_image_create_
#define CGD_IMAGE_COLOR_ALLOCATE         cgd_image_color_allocate_
#define CGD_IMAGE_COLOR_ALLOCATE_ALPHA   cgd_image_color_allocate_alpha_
#define CGD_IMAGE_COLOR_DEALLOCATE       cgd_image_color_deallocate_
#define CGD_IMAGE_SET_PIXEL              cgd_image_set_pixel_
#define CGD_IMAGE_GET_PIXEL              cgd_image_get_pixel_
#define CGD_IMAGE_SET_LINE_THICKNESS     cgd_image_set_line_thickness_
#define CGD_IMAGE_LINE                   cgd_image_line_
#define CGD_IMAGE_RECTANGLE              cgd_image_rectangle_
#define CGD_IMAGE_FILLED_RECTANGLE       cgd_image_filled_rectangle_
#define CGD_IMAGE_ELLIPSE                cgd_image_ellipse_
#define CGD_IMAGE_FILLED_ELLIPSE         cgd_image_filled_ellipse_
#define CGD_IMAGE_ARC                    cgd_image_arc_
#define CGD_IMAGE_FILLED_ARC             cgd_image_filled_arc_
#define CGD_IMAGE_DESTROY                cgd_image_destroy_
#define CGD_IMAGE_GIF                    cgd_image_gif_
#define CGD_IMAGE_PNG                    cgd_image_png_
#define CGD_IMAGE_GD                     cgd_image_gd_
#define CGD_IMAGE_JPEG                   cgd_image_jpeg_
#define CGD_IMAGE_XPM                    cgd_image_xpm_
#define CGD_FONT_TINY                    cgd_font_tiny_
#define CGD_FONT_SMALL                   cgd_font_small_
#define CGD_FONT_MEDIUM_BOLD             cgd_font_medium_bold_
#define CGD_FONT_LARGE                   cgd_font_large_
#define CGD_FONT_GIANT                   cgd_font_giant_
#define CGD_FONT_WIDTH                   cgd_font_width_
#define CGD_FONT_HEIGHT                  cgd_font_height_
#define CGD_IMAGE_STRING                 cgd_image_string_
#define CGD_IMAGE_STRING_UP              cgd_image_string_up_
#define CGD_IMAGE_CHAR                   cgd_image_char_
#define CGD_IMAGE_CHAR_UP                cgd_image_char_up_
#define CGD_IMAGE_POLYGON                cgd_image_polygon_
#define CGD_IMAGE_FILLED_POLYGON         cgd_image_filled_polygon_
#define CGD_IMAGE_FILL                   cgd_image_fill_
#define CGD_IMAGE_FILL_TO_BORDER         cgd_image_fill_to_border_
#define CGD_WIDTH                        cgd_width_
#define CGD_HEIGHT                       cgd_height_
#define CGD_IMAGE_CREATE_FROM_GIF        cgd_image_create_from_gif_
#define CGD_IMAGE_CREATE_FROM_GD         cgd_image_create_from_gd_
#define CGD_IMAGE_CREATE_FROM_PNG        cgd_image_create_from_png_
#define CGD_IMAGE_CREATE_FROM_JPEG       cgd_image_create_from_jpeg_
#define CGD_IMAGE_CREATE_FROM_XBM        cgd_image_create_from_xbm_
#define CGD_IMAGE_CREATE_FROM_XPM        cgd_image_create_from_xpm_
#define CGD_IMAGE_SET_BRUSH              cgd_image_set_brush_
#define CGD_IMAGE_SET_TILE               cgd_image_set_tile_
#define CGD_IMAGE_SET_STYLE              cgd_image_set_style_
#define CGD_RED                          cgd_red_
#define CGD_GREEN                        cgd_green_
#define CGD_BLUE                         cgd_blue_
#define CGD_ALPHA                        cgd_alpha_
#define CGD_GREYSCALE                    cgd_greyscale_
#define CGD_SCALE                        cgd_scale_
#define CGD_PIXELCOLOR                   cgd_pixel_color_
#define CGD_IMAGE_TRUECOLOR_TO_PALETTE   cgd_image_truecolor2palette_
#define CGD_IMAGE_BOUNDS_SAFE            cgd_image_bounds_safe_
#define CGD_IMAGE_COLOR_CLOSEST          cgd_image_color_closest_
#define CGD_IMAGE_COLOR_CLOSEST_HWB      cgd_image_color_closest_hwb_
#define CGD_IMAGE_COLOR_CLOSEST_ALPHA    cgd_image_color_closest_alpha_
#define CGD_IMAGE_COLOR_RESOLVE          cgd_image_color_resolve_
#define CGD_IMAGE_COLOR_RESOLVE_ALPHA    cgd_image_color_resolve_alpha_
#define CGD_IMAGE_COLOR_EXACT            cgd_image_color_exact_
#define CGD_IMAGE_COLORS_TOTAL           cgd_image_colors_total_
#define CGD_IMAGE_GET_INTERLACED         cgd_image_get_interlaced_
#define CGD_IMAGE_GET_TRANSPARENT        cgd_image_get_transparent_
#define CGD_IMAGE_TRANSPARENT            cgd_image_transparent_
#define CGD_IMAGE_INTERLACE              cgd_image_interlace_
#define CGD_IMAGE_SET_ALPHA_BLENDING     cgd_image_set_alpha_blending_
#define CGD_IMAGE_SAVE_ALPHA             cgd_image_save_alpha_
#define CGD_IMAGE_STRING_FT              cgd_image_string_ft_
#define CGD_BRUSHED                      cgd_brushed_
#define CGD_STYLED                       cgd_styled_
#define CGD_STYLED_BRUSHED               cgd_styled_brushed_
#define CGD_TILED                        cgd_tiled_
#define CGD_TRANSPARENT                  cgd_transparent_
#define CGD_PIE                          cgd_pie_
#define CGD_CHORD                        cgd_chord_
#define CGD_NOFILL                       cgd_nofill_
#define CGD_EDGED                        cgd_edged_
#define CGD_ANTI_ALIASED                 cgd_anti_aliased_
#define CGD_IMAGE_COMPARE                cgd_image_compare_
#define CGD_IMAGE_JPEG_BUFFER_PUT        cgd_image_jpeg_buffer_put_
#define CGD_IMAGE_JPEG_BUFFER_GET        cgd_image_jpeg_buffer_get_
#include <sys/types.h>
#include <sys/stat.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <gd.h>
#include <gd_io.h>

/******************************************************************************/
void CGD_IMAGE_CREATE_TRUECOLOR(int *x, int *y, gdImagePtr *ptr)
{
    *ptr = gdImageCreateTrueColor(*x, *y);
    return;
}
/******************************************************************************/
void CGD_IMAGE_TRUECOLOR_TO_PALETTE(gdImagePtr *ptr, int *c1, int *c2)
{
    gdImageTrueColorToPalette(*ptr, *c1, *c2);
    return;
}
/******************************************************************************/
void CGD_IMAGE_CREATE(int *x, int *y, gdImagePtr *ptr)
{
    *ptr = gdImageCreate(*x, *y);
    return;
}
/******************************************************************************/
void CGD_IMAGE_COLOR_ALLOCATE(gdImagePtr *ptr,
                              int *r, int *g, int *b,
                              int *color)
{
    *color = gdImageColorAllocate(*ptr, *r, *g, *b);
    return;
}
/******************************************************************************/
void CGD_IMAGE_COLOR_ALLOCATE_ALPHA(gdImagePtr *ptr,
                                    int *r, int *g, int *b, int *a,
                                    int *color)
{
    *color = gdImageColorAllocateAlpha(*ptr, *r, *g, *b, *a);
    return;
}
/******************************************************************************/
void CGD_IMAGE_COLOR_DEALLOCATE(gdImagePtr *ptr, int *color)
{
    gdImageColorDeallocate(*ptr, *color);
    return;
}

/******************************************************************************/
void CGD_IMAGE_SET_PIXEL(gdImagePtr *ptr, int *x, int *y, int *color)
{
    gdImageSetPixel(*ptr, *x, *y, *color);
    return;
}
/******************************************************************************/
void CGD_IMAGE_GET_PIXEL(gdImagePtr *ptr, int *x, int *y, int *color)
{
    *color = gdImageGetPixel(*ptr, *x, *y);
    return;
}
/******************************************************************************/
void CGD_IMAGE_SET_LINE_THICKNESS(gdImagePtr *ptr, int *thickness)
{
    gdImageSetThickness(*ptr, *thickness);
    return;
}
/******************************************************************************/
void CGD_IMAGE_LINE(gdImagePtr *ptr,
                    int *x1, int *y1, int *x2, int *y2,
                    int *color)
{
    gdImageLine(*ptr, *x1, *y1, *x2, *y2, *color);
    return;
}
/******************************************************************************/
void CGD_IMAGE_RECTANGLE(gdImagePtr *ptr,
                         int *x1, int *y1, int *x2, int *y2,
                         int *color)
{
    gdImageRectangle(*ptr, *x1, *y1, *x2, *y2, *color);
    return;
}
/******************************************************************************/
void CGD_IMAGE_FILLED_RECTANGLE(gdImagePtr *ptr,
                                int *x1, int *y1, int *x2, int *y2,
                                int *color)
{
    gdImageFilledRectangle(*ptr, *x1, *y1, *x2, *y2, *color);
    return;
}

/******************************************************************************/
void CGD_IMAGE_FILLED_ELLIPSE(gdImagePtr *ptr,
                              int *x, int *y, int *w, int *h,
                              int *color)
{
    gdImageFilledEllipse(*ptr, *x, *y, *w, *h, *color);
    return;
}
/******************************************************************************/
void CGD_IMAGE_ARC(gdImagePtr *ptr,
                   int *x, int *y, int *w, int *h, int *s, int *e,
                   int *color)
{
    gdImageArc(*ptr, *x, *y, *w, *h, *s, *e, *color);
    return;
}
/******************************************************************************/
void CGD_IMAGE_FILLED_ARC(gdImagePtr *ptr,
                          int *x, int *y, int *w, int *h, int *s, int *e,
                          int *color, int *style)
{
    gdImageFilledArc(*ptr, *x, *y, *w, *h, *s, *e, *color, *style);
    return;
}
/******************************************************************************/
void CGD_IMAGE_DESTROY(gdImagePtr *ptr)
{
    gdImageDestroy(*ptr);
    return;
}
/******************************************************************************/
void CGD_IMAGE_FILE(gdImagePtr *ptr, const char *file, int *status)
{
    if(strlen(file) > 0 && strcmp(file, "-") != 0) {
        *status = gdImageFile(*ptr, file);
    } else  *status = 2;
    return;
}
/******************************************************************************/
void CGD_IMAGE_PNG(gdImagePtr *ptr, const char *file)
{
    if(strcmp(file, "-") == 0) {
        gdImagePng(*ptr, stdout);
    } else {
        FILE *f;
        f = fopen(file, "wb");
        if(f == NULL) {
            *ptr = 0;
        } else {
            gdImagePng(*ptr, f);
            (void)fclose(f);
        }
    }
    return;
}
/******************************************************************************/
void CGD_IMAGE_GD(gdImagePtr *ptr, const char *file)
{
    if(strcmp(file, "-") == 0) {
        gdImageGd(*ptr, stdout);
    } else {
        FILE *f;
        f = fopen(file, "wb");
        if(f == NULL) {
            *ptr = 0;
        } else {
            gdImageGd(*ptr, f);
            (void)fclose(f);
        }
    }
    return;
}
void CGD_IMAGE_POLYGON(gdImagePtr *ptr,
                       int *n, int x[], int y[],
                       int *color)
{
    /*
        gdPoint points[*n];
        int i;

        for(i=0;i<*n;i++)
        {
            points[i].x = x[i];
            points[i].y = y[i];
        }
        gdImagePolygon(*ptr,points,*n,*color);
        return;
    */

    int i;
    gdPoint *points;

    points = malloc((*n) * sizeof(gdPoint));
    for(i = 0; i < *n; i++) {
        points[i].x = x[i];
        points[i].y = y[i];
    }
    gdImagePolygon(*ptr, points, *n, *color);
    free(points);
    return;
}
/******************************************************************************/
void CGD_IMAGE_FILLED_POLYGON(gdImagePtr *ptr,
                              int *n, int x[], int y[],
                              int *color)
{
    int i;
    gdPoint *points;

    points = malloc((*n) * sizeof(gdPoint));
    for(i = 0; i < *n; i++) {
        points[i].x = x[i];
        points[i].y = y[i];
    }
    gdImageFilledPolygon(*ptr, points, *n, *color);
    free(points);
    return;
    /*
        gdPoint points[*n];
        int i;

        for(i=0;i<*n;i++)
        {
            points[i].x = x[i];
            points[i].y = y[i];
        }
        gdImageFilledPolygon(*ptr,points,*n,*color);
        return;
    */
}
/******************************************************************************/
void CGD_IMAGE_FILL(gdImagePtr *ptr, int *x, int *y, int *color)
{
    gdImageFill(*ptr, *x, *y, *color);
    return;
}
/******************************************************************************/
void CGD_IMAGE_FILL_TO_BORDER(gdImagePtr *ptr, int *x, int *y, int *bc, int *c)
{
    gdImageFillToBorder(*ptr, *x, *y, *bc, *c);
    return;
}
/******************************************************************************/
void CGD_WIDTH(gdImagePtr *ptr, int *width)
{
    *width = gdImageSX(*ptr);
    return;
}
/******************************************************************************/
void CGD_HEIGHT(gdImagePtr *ptr, int *height)
{
    *height = gdImageSY(*ptr);
    return;
}
/******************************************************************************/
void CGD_IMAGE_CREATE_FROM_FILE(const char *file, gdImagePtr *ptr)
{

    if(0 == strcmp(file, "")) {
        *ptr = 0;
    } else {
        *ptr = gdImageCreateFromFile(file);
    }
    return;
}
/******************************************************************************/
void CGD_IMAGE_CREATE_FROM_GD(const char *file, gdImagePtr *ptr)
{
    FILE *f;
    f = fopen(file, "rb");
    if(f == NULL) {
        *ptr = 0;
    } else {
        *ptr = gdImageCreateFromGd(f);
        (void)fclose(f);
    }
    return;
}
/******************************************************************************/
void CGD_IMAGE_CREATE_FROM_PNG(const char *file, gdImagePtr *ptr)
{
    FILE *f;
    f = fopen(file, "rb");
    if(f == NULL) {
        *ptr = 0;
    } else {
        *ptr = gdImageCreateFromPng(f);
        (void)fclose(f);
    }
    return;
}
void CGD_IMAGE_JPEG(gdImagePtr *ptr, const char *file, int *quality)
{
    if(0 == strcmp(file, "-")) {
        gdImageJpeg(*ptr, stdout, *quality);
    } else {
        FILE *f;
        f = fopen(file, "wb");
        if(f == NULL) {
            *ptr = 0;
        } else {
            gdImageJpeg(*ptr, f, *quality);
            (void)fclose(f);
        }
    }
#ifdef _DEBUG
    	FILE*in = fopen (file, "rb");
      if (!in) {
        fprintf(stderr, "Can't open file test/gdtest.jpg.\n");
        exit (1);
      }
      gdImagePtr im2 = gdImageCreateFromJpeg (in);
      fclose (in);
      if (!im2) {
        fprintf(stderr, "gdImageCreateFromJpeg failed.\n");
        exit (1);
      }
      gdImageDestroy (im2);
#endif
    return;
}

void CGD_IMAGE_CREATE_FROM_JPEG(const char *file, gdImagePtr *ptr)
{
    FILE *f;
    f = fopen(file, "rb");
    if(f == NULL) {
        *ptr = 0;
    } else {
        *ptr = gdImageCreateFromJpeg(f);
        if (!*ptr) {
          fprintf(stderr, "gdImageCreateFromJpeg failed. Attempted to open %s\n", file);
          exit (1);
        }
        (void)fclose(f);
    }
    return;
}

/******************************************************************************/
void CGD_IMAGE_CREATE_FROM_XBM(const char *file, gdImagePtr *ptr)
{
    FILE *f;
    f = fopen(file, "rb");
    if(f == NULL) {
        *ptr = 0;
    } else {
        *ptr = gdImageCreateFromXbm(f);
        (void)fclose(f);
    }
    return;
}
/******************************************************************************/
#ifdef WITH_XPM
void CGD_IMAGE_CREATE_FROM_XPM(char *file, gdImagePtr *ptr)
{
    *ptr = gdImageCreateFromXpm(file);
    return;
}

void CGD_IMAGE_XPM(gdImagePtr *ptr, const char *file)
{
    int status;
    if(0 != strcmp(file, "")) {
        status = gdImageFile(*ptr, file);
    } else {
        *ptr = 0;

    }
    return;
}

#endif
/******************************************************************************/
void CGD_IMAGE_SET_BRUSH(gdImagePtr *ptr, gdImagePtr *brush)
{
    gdImageSetBrush(*ptr, *brush);
    return;
}
/******************************************************************************/
void CGD_IMAGE_SET_TILE(gdImagePtr *ptr, gdImagePtr *tile)
{
    gdImageSetTile(*ptr, *tile);
    return;
}
/******************************************************************************/
void CGD_IMAGE_SET_STYLE(gdImagePtr *ptr, int *style, int* length)
{
    gdImageSetStyle(*ptr, style, *length);
    return;
}
/******************************************************************************/
void CGD_RED(gdImagePtr *ptr, int *color, int *red)
{
    *red = gdImageRed(*ptr, *color);
    return;
}
/******************************************************************************/
void CGD_GREEN(gdImagePtr *ptr, int *color, int *green)
{
    *green = gdImageGreen(*ptr, *color);
    return;
}
/******************************************************************************/
void CGD_BLUE(gdImagePtr *ptr, int *color, int *blue)
{
    *blue = gdImageBlue(*ptr, *color);
    return;
}
/******************************************************************************/
void CGD_ALPHA(gdImagePtr *ptr, int *color, int *alpha)
{
    *alpha = gdImageAlpha(*ptr, *color);
    return;
}
/******************************************************************************/
void CGD_PIXELCOLOR(gdImagePtr *ptr, int *x, int *y, int *color)
{
    *color = gdImagePalettePixel(*ptr, *x, *y);
    return;
}
/******************************************************************************/
void CGD_SCALE(gdImagePtr *ptr, int *c1, int *c2, gdImagePtr *ptr2)
{
    *ptr2 = gdImageScale(*ptr, *c1, *c2);
    return;
}/******************************************************************************/
void CGD_GREYSCALE(gdImagePtr *ptr, int *grey)
{
    *grey = gdImageGrayScale(*ptr);
    return;
}
/******************************************************************************/
void CGD_IMAGE_BOUNDS_SAFE(gdImagePtr *ptr, int *x, int *y, int *safe)
{
    *safe = gdImageBoundsSafe(*ptr, *x, *y);
    return;
}
/******************************************************************************/
void CGD_IMAGE_COLOR_CLOSEST(gdImagePtr *ptr,
                             int *r, int *g, int *b,
                             int *color)
{
    *color = gdImageColorClosest(*ptr, *r, *g, *b);
    return;
}
/******************************************************************************/
void CGD_IMAGE_COLOR_CLOSEST_HWB(gdImagePtr *ptr,
                                 int *r, int *g, int *b,
                                 int *color)
{
    *color = gdImageColorClosestHWB(*ptr, *r, *g, *b);
    return;
}
/******************************************************************************/
void CGD_IMAGE_COLOR_RESOLVE_ALPHA(gdImagePtr *ptr,
                                   int *r, int *g, int *b, int *a,
                                   int *color)
{
    *color = gdImageColorResolveAlpha(*ptr, *r, *g, *b, *a);
    return;
}
/******************************************************************************/
void CGD_IMAGE_COLOR_RESOLVE(gdImagePtr *ptr,
                             int *r, int *g, int *b,
                             int *color)
{
    *color = gdImageColorResolve(*ptr, *r, *g, *b);
    return;
}
/******************************************************************************/

void CGD_IMAGE_COLOR_EXACT(gdImagePtr *ptr, int *r, int *g, int *b, int *color)
{
    *color = gdImageColorExact(*ptr, *r, *g, *b);
    return;
}
/******************************************************************************/
void CGD_IMAGE_COLORS_TOTAL(gdImagePtr *ptr, int *n)
{
    *n = gdImageColorsTotal(*ptr);
    return;
}
/******************************************************************************/
void CGD_IMAGE_GET_INTERLACED(gdImagePtr *ptr, int *n)
{
    *n = gdImageGetInterlaced(*ptr);
    return;
}
/******************************************************************************/
void CGD_IMAGE_GET_TRANSPARENT(gdImagePtr *ptr, int *n)
{
    *n = gdImageGetTransparent(*ptr);
    return;
}
/******************************************************************************/
void CGD_IMAGE_TRANSPARENT(gdImagePtr *ptr, int *color)
{
    gdImageColorTransparent(*ptr, *color);
    return;
}
/******************************************************************************/
void CGD_IMAGE_INTERLACE(gdImagePtr *ptr, int *interlace)
{
    gdImageInterlace(*ptr, *interlace);
    return;
}
/******************************************************************************/
void CGD_IMAGE_SET_ALPHA_BLENDING(gdImagePtr *ptr, int *blending)
{
    gdImageAlphaBlending(*ptr, *blending);
    return;
}
/******************************************************************************/
void CGD_IMAGE_SAVE_ALPHA(gdImagePtr *ptr, int *save)
{
    gdImageSaveAlpha(*ptr, *save);
    return;
}

/******************************************************************************/
void CGD_IMAGE_COLOR_CLOSEST_ALPHA(gdImagePtr *ptr,
                                   int *r, int *g, int *b, int *a,
                                   int *color)
{
    *color = gdImageColorClosestAlpha(*ptr, *r, *g, *b, *a);
    return;
}
/******************************************************************************/
void CGD_BRUSHED(int *n)
{
    *n = gdBrushed;
    return;
}
/******************************************************************************/
void CGD_STYLED(int *n)
{
    *n = gdStyled;
    return;
}
/******************************************************************************/
void CGD_STYLED_BRUSHED(int *n)
{
    *n = gdStyledBrushed;
    return;
}
/******************************************************************************/
void CGD_TILED(int *n)
{
    *n = gdTiled;
    return;
}
/******************************************************************************/
void CGD_TRANSPARENT(int *n)
{
    *n = gdTransparent;
    return;
}
/******************************************************************************/
void CGD_PIE(int *n)
{
    *n = gdPie;
    return;
}
/******************************************************************************/
void CGD_CHORD(int *n)
{
    *n = gdChord;
    return;
}
/******************************************************************************/
void CGD_NOFILL(int *n)
{
    *n = gdNoFill;
    return;
}
/******************************************************************************/
void CGD_EDGED(int *n)
{
    *n = gdEdged;
    return;
}
/******************************************************************************/
void CGD_ANTI_ALIASED(int *n)
{
    *n = gdAntiAliased;
    return;
}
/******************************************************************************/
void CGD_IMAGE_SET_CLIP(gdImagePtr *ptr, int *x1, int *y1, int *x2, int *y2)
{
    gdImageSetClip(*ptr, *x1, *y1, *x2, *y2);
    return;
}
/******************************************************************************/
void CGD_IMAGE_GET_CLIP(gdImagePtr *ptr, int *x1, int *y1, int *x2, int *y2)
{
    gdImageGetClip(*ptr, x1, y1, x2, y2);
    return;
}
/******************************************************************************/
void CGD_IMAGE_SET_AA(gdImagePtr *ptr, int *c)
{
    gdImageSetAntiAliased(*ptr, *c);
    return;
}
/******************************************************************************/
void CGD_IMAGE_SET_AA_NB(gdImagePtr *ptr, int *c, int *d)
{
    gdImageSetAntiAliasedDontBlend(*ptr, *c, *d);
    return;
}
/******************************************************************************/
void CGD_IMAGE_COMPARE(gdImagePtr *ptr1, gdImagePtr *ptr2, int*c)
{
    *c = gdImageCompare(*ptr1, *ptr2);
    return;
}


void CGD_IMAGE_PNG_BUFFER_PUT(gdImagePtr *ptr, int  *size,  int **array)
{
    struct gdIOCtx *ctx; int status; const void * tmp =  array;
    gdImagePngCtx(*ptr, ctx);
    status = gdPutBuf(tmp, *size, ctx);
}
void CGD_IMAGE_PNG_BUFFER_GET(gdImagePtr *ptr, int  *size,  int **array)
{
    gdIOCtx *ctx; void ** tmp;
    int status;
    gdImagePngCtx(*ptr, ctx);
    status = gdGetBuf(tmp, *size, ctx);
    *array = *((int *)tmp);
}

void CGD_IMAGE_JPEG_BUFFER_PUT(gdImagePtr *ptr, int  *size,  int **array)
{
    gdIOCtx *ctx; const void * tmp =  *array;
    gdImageJpegCtx(*ptr, ctx, *size);
    gdPutBuf(tmp, *size, ctx);
}

void CGD_IMAGE_JPEG_BUFFER_GET(gdImagePtr *ptr, int  *size,  int **array)
{
    gdIOCtx *ctx; void * tmp;
    int status;
    gdImageJpegCtx(*ptr, ctx, *size);
    status = gdGetBuf(tmp, *size, ctx);
    *array = (int *)(tmp);
    return;
}
