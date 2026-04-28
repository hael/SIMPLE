#ifdef _WIN32
#include <errno.h>
// Windows stub for symlink: always fails
int symlink(const char *target, const char *linkpath) {
    errno = ENOSYS; // Function not implemented
    return -1;
}
#endif
