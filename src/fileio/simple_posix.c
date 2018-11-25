/**
 *   simple_posix.c
 *   \brief POSIX system calls for SIMPLE
 *
 *   OS file utilities to allow simple to interface to POSIX functions without using system calls
 *
 *   Michael Eager   2018
 */
#define  _POSIX_C_SOURCE 200809L
#define _THREAD_SAFE
#define _SVID_SOURCE
#define _GNU_SOURCE             /* See feature_test_macros(7) */
#ifdef __APPLE__
#define _DARWIN_C_SOURCE
#include <MacTypes.h>
#define _BSD_SOURCE
#define __DARWIN_C_SOURCE
#endif
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fts.h>               /* file traversal */
#include <fcntl.h>           /* Definition of AT_* constants */
#include <sys/wait.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>          /*     extern int errno;   */
#include <pwd.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>       /* getpid()  */
#include <dirent.h>      /* DIR and scandir,DT_DIR */
#include <time.h>
#include <glob.h>
#include <regex.h>
#include <limits.h>      /* PATH_MAX */
#ifdef __linux__
#include<linux/limits.h>
#endif
#include <ftw.h>
/* By default, print all messages of severity info and above.  */
#ifdef _DEBUG
static int global_debug = 3;
#else
static int global_debug = 1;
#endif
/* Name of this program */
static char *global_progname = "simple_posix";
#ifdef _DEBUG
#define dgprintf if (global_debug >= 3) \
        fprintf (stderr, "%s: debug: (%d) ", global_progname, __LINE__), \
        fprintf
#else
void noprint(FILE *f,...){}
#define dgprintf noprint
#endif
#define eprintf if (global_debug >= 0) \
        fprintf (stderr, "%s: error: (%d) ", global_progname, __LINE__), \
        fprintf
#define LONGSTRLEN 1024
#define MAKEDIR makedir_
#define REMOVEDIR removedir_
#define GET_FILE_LIST get_file_list_
#define LIST_DIRS list_dirs_
#define GLOB_FILE_LIST glob_file_list_
#define GLOB_RM_ALL glob_rm_all_
#define GET_ABSOLUTE_PATHNAME get_absolute_pathname_
#define SUBPROCESS subprocess_
#define WAIT_PID wait_pid_
#define TOUCH touch_
#define GET_SYSINFO get_sysinfo_
#define RMDIRECTORIES rmDirectories_
#define REMOVE_DIRECTORY remove_directory_

// Protect fortran strings by appending '\0' to end
char *F90toCstring(char *str, int len)
{
    char *res; /* C arrays are from 0:len-1 */
    if(len > LONGSTRLEN) {
        fprintf(stderr, "Warning: simple_posix.c: possible strlen error in converting f90 string to c-string\n");
    }
    if((res = (char*)malloc(len + 1))) {
        strncpy(res, str, len);
        res[len] = '\0';
    }
    return res;
}

int fstrlen(char *f, int max_f)
{
    for(; max_f > 0 && (isspace(f[max_f - 1]) || f[max_f - 1] == '\0'); max_f--);
    return max_f;
}

static char *c2fstr(char* cstr, char *fstr, int elem_len, int sizeofcstr)
{
    int i, j;
    /* elem_len includes \0 for C strings. Fortran strings don't have term. \0.
       Useful size of string must be the same in both languages. */
    for(i = 0; i < sizeofcstr / elem_len; i++) {
        for(j = 1; j < elem_len && *cstr; j++) *fstr++ = *cstr++;
        cstr += 1 + elem_len - j;
        for(; j < elem_len; j++) *fstr++ = ' ';
    } /* f951 - Seems to be returning the original fstr. */
    return fstr - sizeofcstr + sizeofcstr / elem_len;
}

void f2cstr(char* f_str, char* c_str, unsigned length)
{
    /* fprintf (stderr, "f2cstr in '%.*s'\n", length, f_str); */

    unsigned i = length - 1;
    strncpy(c_str, f_str, length);

    while(f_str[i] == ' ') {
        c_str[i] = '\0';
        if(i == 0)
            break;
        i--;
    }
    /* fprintf (stderr, "f2cstr out '%s'\n", c_str); */
}

int isdir(char* pathname, int* len, size_t ivf_pathname)
{
    if(pathname == NULL  || *len < 1) {
        perror("Failed : simple_posix.c::isdir inputs poorly defined");
        return 0;
    }
    //printf("isdir pathname :%s: len %d\n",  pathname, *len);
    char* cpathname = F90toCstring(pathname, *len);
    if(!cpathname) {
        printf("%d %s\nisdir failed to create str %s\n", errno, strerror(errno), pathname);
        perror("Failed : simple_posix.c::isdir ");
        return 0;
    }
    dgprintf(stderr, "DEBUG: In simple_posix.c::isdir pathname :%s:\n",  cpathname);

    struct stat ibuf;
    int i = stat(cpathname, &ibuf);
    free(cpathname);
    if(i != 0) {

        dgprintf(stderr, "DEBUG: In simple_posix.c::isdir pathname :%s: does not exist\n",  cpathname);
        //printf("%d %s\nisdir failed to create stat struct %s\n", errno, strerror(errno), cpathname);
        //perror("Failed : simple_posix.c::isdir ");
        return 0;
    }
    dgprintf(stderr, "DEBUG: In isdir stat mode %ud\t isdir %d\n", ibuf.st_mode, S_ISDIR(ibuf.st_mode));
    dgprintf(stderr, "DEBUG: In isdir stat mode s_ifdir %d\n" , (ibuf.st_mode & S_IFMT) == S_IFDIR);

    return S_ISDIR(ibuf.st_mode) ? 1 : 0;
}

// Recursive make directory
int makedir(char *path,
    int* charLen,
    size_t ivf_path)
{
    /* https://gist.github.com/JonathonReinhart/8c0d90191c38af2dcadb102c4e202950 */
    /* Adapted from http://stackoverflow.com/a/2336245/119527 */

    /*
      EACCES The parent directory does not allow write permission to the process,
      or one of the directories in pathname did not allow search permission. (See
      also path_resolution(7).)

      EDQUOT The user's quota of disk blocks or inodes on the filesystem has been
      exhausted.

      EEXIST pathname already exists (not necessarily as a directory). This
      includes the case where pathname is a symbolic link, dangling or not.

      EFAULT pathname points outside your accessible address space.

      ELOOP  Too many symbolic links were encountered in resolving pathname.

      EMLINK The number of links to the parent directory would exceed LINK_MAX.

      ENAMETOOLONG
      pathname was too long.

      ENOENT A directory component in pathname does not exist or is a dangling
      symbolic link.

      ENOMEM Insufficient kernel memory was available.

      ENOSPC The device containing pathname has no room for the new directory.

      ENOSPC The new directory cannot be created because the user's disk quota is
      exhausted.

      ENOTDIR
      A component used as a directory in pathname is not, in fact, a directory.

      EPERM The filesystem containing pathname does not support the creation of
      directories.

      EROFS  pathname refers to a file on a read-only filesystem.

      The following additional errors can occur for mkdirat():

      EBADF  dirfd is not a valid file descriptor.

      ENOTDIR pathname is relative and dirfd is a file descriptor referring to a
      file other than a directory.
    */
    char *cpath = F90toCstring(path, *charLen);
    // fprintf(stderr, "makedir  %d  %d %s  %s\n",  *charLen, strlen(path), path, cpath);
    if(cpath == NULL) {
        printf("%d %s\n makedir failed to convert string (unprotected) %s\n", errno, strerror(errno), path);
        perror("Failed : simple_posix.c::remove_dir ");
        return -2;
    }
    char _path[LONGSTRLEN];
    char *p;
    extern int errno;
    errno = 0;

    /* Copy string so its mutable */
    if(*charLen > sizeof(_path) - 1) {
        errno = ENAMETOOLONG;
        free(cpath);
        return -1;
    }

    strcpy(_path, cpath);
    // fprintf(stderr, "makedir %d [%s]\n" , strlen(_path), _path);
    free(cpath);
    /* Iterate the string */
    for(p = _path + 1; *p; p++) {
        if(*p == '/') {
            /* Temporarily truncate */
            *p = '\0';

            if(mkdir(_path, S_IRWXU|S_IRWXG) != 0) {
                if(errno != EEXIST) {
                    fprintf(stderr, "makedir %s\nerrno:%d msg:%s\n", _path, errno, strerror(errno));
                    perror("Failed : mkdir in simple_posix::makedir");
                    return -1;
                }
            }
            *p = '/';
        }
    }
    if(mkdir(_path, S_IRWXU|S_IRWXG) != 0) {
        if(errno != EEXIST) {
            fprintf(stderr, "makedir %s\nerrno:%d msg:%s\n", _path, errno, strerror(errno));
            perror("Failed : mkdir in simple_posix::makedir");
            return -1;
        }
    }
    return 0;
}

// Recursive remove directory
int removedir(char *path, int* len, int* count, size_t ivf_path)
{
    /* https://gist.github.com/JonathonReinhart/8c0d90191c38af2dcadb102c4e202950 */
    /* Adapted from http://stackoverflow.com/a/2336245/119527 */

    /*
      EACCES Write access to the directory containing pathname was not allowed, or
      one of the directories in the path prefix of pathname did not allow search
      permission. (See also path_resolution(7).

      EBUSY pathname is currently in use by the system or some process that
      prevents its removal. On Linux this means pathname is currently used as a
      mount point or is the root directory of the calling process.

      EFAULT pathname points outside your accessible address space.

      EINVAL pathname has .  as last component.

      ELOOP  Too many symbolic links were encountered in resolving pathname.

      ENAMETOOLONG
      pathname was too long.

      ENOENT A directory component in pathname does not exist or is a dangling
      symbolic link.

      ENOMEM Insufficient kernel memory was available.

      ENOTDIR pathname, or a component used as a directory in pathname, is not, in
      fact, a directory.

      ENOTEMPTY pathname contains entries other than . and .. ; or, pathname has
      .. as its final component. POSIX.1 also allows EEXIST for this condition.

      EPERM The directory containing pathname has the sticky bit (S_ISVTX) set and
      the process's effective user ID is neither the user ID of the file to be
      deleted nor that of the directory containing it, and the process is not
      privileged (Linux: does not have the CAP_FOWNER capability).

      EPERM The filesystem containing pathname does not support the removal of
      directories.

      EROFS  pathname refers to a directory on a read-only filesystem.
    */
    char *cpath = F90toCstring(path, *len);
    if(cpath == NULL) {
        printf("%d %s\n removedir failed to convert string (unprotected) %s\n", errno, strerror(errno), path);
        perror("Failed : simple_posix.c::remove_dir ");
        return -2;
    }
    //  const size_t len = strlen(cpath);
    char _path[LONGSTRLEN];
    char *p;
    errno = 0;

    /* Copy string so its mutable */
    if(*len > (int)sizeof(_path) - 1) {
        errno = ENAMETOOLONG;
        free(cpath);
        return -1;
    }
    strncpy(_path, cpath, *len);
    free(cpath);
    printf("DEBUG:in removedir  rmdir %s\n", _path);
    if(rmdir(_path) != 0) {
        if(errno != ENOTEMPTY)
            return -1;
    }
    /* Iterate the string */
    for(p = _path + *len - 1; *p; p--) {
        if(*p == '/') {
            /* Temporarily truncate */
            *p = '\0';
            printf("DEBUG:  rmdir %s\n", _path);
            if(rmdir(_path) != 0) {
                if(errno != ENOTEMPTY)
                    printf("rmdir %s\nerrno:%d msg:%s\n", _path, errno, strerror(errno));
                perror("Failed : rmdir in simple_posix::removedir");
                return -1;
            }

            *p = '/';
        }
    }
    if(rmdir(_path) != 0) {
        if(errno != ENOTEMPTY) {
            printf("rmdir %s\nerrno:%d msg:%s\n", _path, errno, strerror(errno));
            perror("Failed : rmdir in simple_posix::removedir");
            return -1;
        }
    }
    *count = 1;
    printf("DEBUG:in removedir  completed rmdir %s\n", _path);
    return 0;
}

// get_file_list emulates terminal command 'ls '
int get_file_list(char * path, int*len, char * ext, int *count, size_t ivf_path)
{
    char *cpath = F90toCstring(path, *len);
    if(cpath == NULL) {
        printf("%d %s\nget_file_list failed to convert string (unprotected) %s\n", errno, strerror(errno), path);
        perror("Failed : simple_posix.c::get_file_list ");
        return -1;
    }
    dgprintf(stderr, "DEBUG: In get_file_list, %15s:%s\n", "path", cpath);
    dgprintf(stderr, "DEBUG: In get_file_list, %15s:%s\n", "ext", ext);
    DIR *d;
    d = opendir(cpath); free(cpath);
    if(d) {
        struct dirent *elem;
        int fcount = 0;
        while((elem = readdir(d)) != NULL) {
            if(elem-> d_type != DT_DIR) {
                fcount++;
            }
        }
        rewinddir(d);
        dgprintf(stderr, " get_file_list size %d\n", fcount);
        FILE* f = fopen("__simple_filelist__", "w");
        fcount = 0;
        char *cext = F90toCstring(ext, 3);
        if(strcmp(ext, "")) {
            while((elem = readdir(d)) != NULL) {
                if(elem->d_type != DT_DIR) {
                    fprintf(f, "%s\n", elem->d_name);
                    fcount++;
                }
            }
        } else {
            while((elem = readdir(d)) != NULL) {
                if(elem->d_type != DT_DIR) {
                    if(strstr(elem->d_name, ext) != NULL) {
                        fprintf(f, "%s\n", elem->d_name);
                        fcount++;
                    }
                }
            }
        }
        *count = fcount;
        free(cext);
        fclose(f);
        closedir(d);
    } else {
        printf("%d %s\nget_file_list opendir failed to open %s\n", errno, strerror(errno), path);
        perror("Failed : simple_posix.c::get_file_list ");
        return -1;
    }
    return 0;
}

// list directories in directory 'path'
int list_dirs(char * path, int*len, char * fout, int*len_fout, int* count, size_t ivf_path)
{
    char *cpath = F90toCstring(path, *len);
    char *cfout = F90toCstring(fout, *len_fout);
    if(cpath == NULL) {
        printf("%d %s\nlist_dirs failed to open convert string (unprotected) %s\n", errno, strerror(errno), path);
        perror("Failed : simple_posix.c::list_dirs ");
        return -1;
    }
    if(cfout == NULL) {
        printf("%d %s\nlist_dirs failed to convert string (unprotected) %s\n", errno, strerror(errno), fout);
        perror("Failed : simple_posix.c::list_dirs ");
        return -1;
    }
    extern int errno;
    DIR *d;
    int fcount = 0;
    d = opendir(cpath);
    free(cpath);
    if(d) {
        FILE* f = fopen(cfout, "w");
        if(f) {
            fcount = 0;
            struct dirent *dir;
            while((dir = readdir(d)) != NULL) {
                if(dir-> d_type == DT_DIR && strcmp(dir->d_name, ".") != 0 && strcmp(dir->d_name, "..") != 0) {
                    fprintf(f, "%s\n", dir->d_name);
                    fcount++;
                }
            }
            fclose(f);
        } else {
            printf("%d %s\nlist_dirs opendir failed to open temp file\n", errno, strerror(errno));
            perror("Failed : simple_posix.c::list_dirs ");
            return -1;
        }
        closedir(d);
    } else {
        printf("%d %s\nlist_dirs opendir failed to open %s\n", errno, strerror(errno), path);
        perror("Failed : simple_posix.c::list_dirs ");
        return -1;
    }
    *count = fcount;
    return 0;
}

int wait_pid(int* cpid)
{
    int status;
    pid_t w;
    if(*cpid > 0) {
        do {
          pid_t c = (pid_t) * cpid;
            w = waitpid(c , &status, WNOHANG | WUNTRACED | WCONTINUED);
            if(w == -1) {
                perror("waitpid");
                exit(EXIT_FAILURE);
            }

            if(WIFEXITED(status)) {
                printf("exited, status=%d\n", WEXITSTATUS(status));
            } else if(WIFSIGNALED(status)) {
                printf("killed by signal %d\n", WTERMSIG(status));
            } else if(WIFSTOPPED(status)) {
                printf("stopped by signal %d\n", WSTOPSIG(status));
            } else if(WIFCONTINUED(status)) {
                printf("continued\n");
            }
        } while(!WIFEXITED(status) && !WIFSIGNALED(status));
    } else return -1;

    return (int) w;
}

char * lrealpath(const char *filename)
{
    /* Method 3: Now we're getting desperate!  The system doesn't have a
       compile time buffer size and no alternative function.  Query the
       OS, using pathconf(), for the buffer limit.  Care is needed
       though, some systems do not limit PATH_MAX (return -1 for
       pathconf()) making it impossible to pass a correctly sized buffer
       to realpath() (it could always overflow).  On those systems, we
       skip this.  */
    /* Find out the max path size.  */
    long path_max = pathconf("/", _PC_PATH_MAX);
    if(path_max > 0) {
        /* PATH_MAX is bounded.  */
        char *buf, *rp, *ret;
        buf = (char *) malloc(path_max);
        if(buf == NULL)
            return NULL;
        rp = realpath(filename, buf);
        ret = strdup(rp ? rp : filename);
        free(buf);
        return ret;
    }
    return NULL;
}

// returns a null-terminated string containing the canonicalized absolute pathname corresponding to path.
int  get_absolute_pathname(char* in, int* inlen, char* out, int* outlen)
{
    // realpath(in, resolved_path);
    //      printf("\n%s\n",resolved_path);
    //      return 0;
    char *filein = F90toCstring(in, *inlen);
#if defined(HAVE_CANONICALIZE_FILE_NAME)
    char *resolved = canonicalize_file_name(filein);
#else
    char *resolved = lrealpath(filein);
#endif

    if(resolved == NULL) {
      fprintf(stderr, "%d %s\nget_absolute_path failed to canonicalize  file %s\n", errno, strerror(errno), filein);
        perror("Failed : simple_posix.c::get_absolute_pathname ");
        return 1;
    }
    *outlen = strlen(resolved) ;
    if (*outlen > LONGSTRLEN){
      fprintf(stderr, "get_absolute_path: lrealpath returned string longer than str max (%d): \nstrlen %d  \npath:%s \n",LONGSTRLEN, *outlen, resolved);
    }
    strncpy(out, resolved, *outlen);
    out[*outlen] = '\0';

    c2fstr(resolved, out, *outlen, sizeof(resolved));
    out[0] = '/';

    free(filein);
    free(resolved);

    return 0;
}

#ifdef __APPLE__
#include <mach/mach.h>
#include <mach/vm_statistics.h>
#include <mach/mach_types.h>
#include <mach/mach_init.h>
#include <mach/mach_host.h>
#include <sys/sysctl.h>
#include <libproc.h>
#include <mach/vm_types.h>
#include <mach/vm_task.h>
#include <mach/task.h>
#else
#include <sys/sysinfo.h>
#endif
int get_sysinfo(long* HWMusage, long*totalram, long* sharedram, long* bufferram, long* totalhigh)
{
    *sharedram = 0L;       /* Amount of shared memory */
    *bufferram = 0L;       /* Memory used by buffers */
    *totalhigh = 0L;       /* Total high memory size */
    *HWMusage = 0L;        /* high water mark -- last peak memory used */

#ifdef __APPLE__
//Total RAM used
    int mib[2];
    int64_t physical_memory;
    mib[0] = CTL_HW;
    mib[1] = HW_MEMSIZE;
    size_t length = sizeof(int64_t);
    sysctl(mib, 2, &physical_memory, &length, NULL, 0);

    struct proc_taskinfo p_info;
    int result = proc_pidinfo( getpid(), PROC_PIDTASKINFO, 0,  & p_info, sizeof( p_info ) ) ;

    //Virtual Memory Currently Used by my Process
    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

    if(KERN_SUCCESS != task_info(mach_task_self(),
                                 TASK_BASIC_INFO, (task_info_t)&t_info,
                                 &t_info_count)) {
        perror("unable to get virtual memswap usage by calling taskinfo");
        return -1;
    }
// resident size is in t_info.resident_size;
// virtual size is in t_info.virtual_size;
    // RAM currently used
    vm_size_t page_size;
    mach_port_t mach_port;
    mach_msg_type_number_t count;
    vm_statistics64_data_t vm_stats;
    long long free_memory = 0;
    long long used_memory = 0;
    mach_port = mach_host_self();
    count = sizeof(vm_stats) / sizeof(natural_t);
    if(KERN_SUCCESS == host_page_size(mach_port, &page_size) &&
       KERN_SUCCESS == host_statistics64(mach_port, HOST_VM_INFO,
                                         (host_info64_t)&vm_stats, &count)) {
        free_memory = (int64_t)vm_stats.free_count * (int64_t)page_size;

        used_memory = ((int64_t)vm_stats.active_count +
                                 (int64_t)vm_stats.inactive_count +
                                 (int64_t)vm_stats.wire_count) * (int64_t)page_size;
        printf("free memory: %lld\nused memory: %lld\n", free_memory, used_memory);
    }


    *totalram =  physical_memory;           /* Total usable main memory size (bytes)*/
    *sharedram = t_info.virtual_size;       /* Amount of shared memory */
    *bufferram = t_info.resident_size;      /* Memory used by buffers */
    *totalhigh = free_memory + used_memory; /* Total high water mark memory size */
    *HWMusage =  used_memory;               /* high memory size used */

#else
    struct sysinfo s;
    if( sysinfo(&s) ) {
      perror("simple_posix.c::get_sysinfo unable to get mem usage by calling sysinfo");
      return -1;
    } else {
      *totalram = s.totalram;           /* Total usable main memory size */
      //*freeram=s.freeram;               /* Available memory size */
      *sharedram = s.sharedram;       /* Amount of shared memory */
      *bufferram = s.bufferram;              /* Memory used by buffers */
      *totalhigh = s.totalhigh;              /* Total high memory size */
      *HWMusage = s.totalhigh - s.freehigh; /* high memory size used */
    }
#endif
    return 0;
}

int unlink_cb(const char *fpath, const struct stat *sb, int typeflag, struct FTW *ftwbuf)
{
    int rv = remove(fpath);
    if(rv)
        perror(fpath);
    return rv;
}

static int rmFiles(const char *pathname, const struct stat *sbuf, int type, struct FTW *ftwb)
{
    if(remove(pathname) < 0) {
        perror("ERROR: remove");
        return -1;
    }
    return 0;
}

int remove_directory_recursive(char *path, int *len, int*count)
{
    char *cpath = F90toCstring(path, *len);
    DIR *d = opendir(cpath);
    size_t path_len = strlen(cpath);
    int r = -1;
    free(cpath);
    if(d) {
        struct dirent *p;
        r = 0;
        while(!r && (p = readdir(d))) {
            int r2 = -1;
            char *buf;
            size_t len;
            /* Skip the names "." and ".." as we don't want to recurse on them. */
            if(!strcmp(p->d_name, ".") || !strcmp(p->d_name, ".."))
                continue;
            len = path_len + strlen(p->d_name) + 2;
            buf = malloc(len);
            if(buf) {
                struct stat statbuf;
                snprintf(buf, len, "%s/%s", path, p->d_name);
                if(!stat(buf, &statbuf)) {
                    if(S_ISDIR(statbuf.st_mode)) {
                        int len2 = (int)strlen(buf);
                        r2 = remove_directory_recursive(buf, &len2 , count);

                    } else {
                        r2 = unlink(buf);
                        *count = *count + 1;
                        FILE* f = fopen("__simple_filelist__", "a");
                        fprintf(f, "%s\n", buf);
                        fclose(f);
                    }
                }
                free(buf);
            }
            r = r2;
        }
        closedir(d);
    }
    if(!r) {
        r = rmdir(path);
        *count = *count + 1;
        FILE* f = fopen("__simple_filelist__", "a");
        fprintf(f, "%s\n", path);
        fclose(f);
    }
    return r;
}

int remove_directory(char *path, int *len, int*count)
{
    FILE* f = fopen("__simple_filelist__", "w");
    fclose(f);
    *count = 0;
    return remove_directory_recursive(path, len, count);
}

int touch(char* path, int*len)
{
    char*cpath = F90toCstring(path, *len);
    int fd2 = open(cpath, O_RDWR | O_CREAT, 0777);
    free(cpath);
    if(fd2 != -1) {
        // use file descriptor
        close(fd2);
    }
}

// POSIX regular expression
int regexp_match(char* srcstr,  int*srclen, char* regexString , int*rgxlen)
{
  printf("Regexp_match:\nsrcstr:%s\nregex:%s\n", srcstr, regexString);
  regex_t regex;
  int reti;
  char msgbuf[100];
  char*rgxstr = F90toCstring(regexString, *rgxlen);
  /* Compile regular expression */
  reti = regcomp(&regex, rgxstr, REG_EXTENDED|REG_ICASE|REG_NOSUB);
  if (reti) {
    fprintf(stderr, "Could not compile regex\n");
    free(rgxstr);
    return -1;
  }

  char*instr = F90toCstring(srcstr, *srclen);
  /* Execute regular expression */
  reti = regexec(&regex, instr, 0, NULL, 0);
  if (!reti) {
    puts("Match");
  }
  else if (reti == REG_NOMATCH) {
    puts("No match");
  }
  else {
    regerror(reti, &regex, msgbuf, sizeof(msgbuf));
    fprintf(stderr, "Regex match failed: %s\n", msgbuf);
    free(rgxstr);
    free(instr);
    return -1;
  }
  /* Free memory allocated to the pattern buffer by regcomp() */
  regfree(&regex);
  free(rgxstr);
  free(instr);
  return reti;
}
