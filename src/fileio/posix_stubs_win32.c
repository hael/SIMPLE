#ifdef _WIN32
#include <errno.h>
#include <stddef.h>
#include <sys/types.h>

// File sync
int fsync(int fd) { errno = ENOSYS; return -1; }

// Socket functions
int connect(int sockfd, const void *addr, size_t addrlen) { errno = ENOSYS; return -1; }
ssize_t send(int sockfd, const void *buf, size_t len, int flags) { errno = ENOSYS; return -1; }
int inet_pton(int af, const char *src, void *dst) { errno = ENOSYS; return -1; }
unsigned short htons(unsigned short hostshort) { errno = ENOSYS; return 0; }
int socket(int domain, int type, int protocol) { errno = ENOSYS; return -1; }
int accept(int sockfd, void *addr, void *addrlen) { errno = ENOSYS; return -1; }
int listen(int sockfd, int backlog) { errno = ENOSYS; return -1; }
int bind(int sockfd, const void *addr, size_t addrlen) { errno = ENOSYS; return -1; }
int setsockopt(int sockfd, int level, int optname, const void *optval, size_t optlen) { errno = ENOSYS; return -1; }

// Process control
pid_t waitpid(pid_t pid, int *status, int options) { errno = ENOSYS; return -1; }
int pipe(int pipefd[2]) { errno = ENOSYS; return -1; }
int fcntl(int fd, int cmd, ...) { errno = ENOSYS; return -1; }
int c_fcntl(int fd, int cmd, int arg) { errno = ENOSYS; return -1; } /* unix_fcntl's fixed-argument binding; unix_macro.c on POSIX */
pid_t fork(void) { errno = ENOSYS; return -1; }
int kill(pid_t pid, int sig) { errno = ENOSYS; return -1; }

// poll (POSIX I/O multiplexing - not available in Windows UCRT)
int poll(void *fds, unsigned long nfds, int timeout) { errno = ENOSYS; return -1; }

#endif // _WIN32
