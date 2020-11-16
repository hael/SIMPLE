#ifndef __STARFILE_GENERIC_H_
#define __STARFILE_GENERIC_H_

typedef float RFLOAT;
#ifndef ABS
#define ABS(x) (((x) >= 0) ? (x) : (-(x)))
#endif
#define REPORT_ERROR(ErrormMsg) throw SimpleError(ErrormMsg, __FILE__, __LINE__)
#define FLOOR(x) (((x) == (int)(x)) ? (int)(x):(((x) > 0) ? (int)(x) : (int)((x) - 1)))
#ifndef SGN
#define SGN(x) (((x) >= 0) ? 1 : -1)
#endif
#ifndef ROUND
#define ROUND(x) (((x) > 0) ? (int)((x) + 0.5) : (int)((x) - 0.5))
#endif


#endif
