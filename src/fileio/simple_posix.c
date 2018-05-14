/**
 *   simple_posix.c
 *   \brief POSIX system calls for SIMPLE
 *
 *   OS file utilities to allow simple to interface to POSIX functions without using system calls
 *
 *   Michael Eager   2018
 */
#define _DEBUG 1
#define  _POSIX_C_SOURCE 200809L
#define _THREAD_SAFE
//#define _FORTIFY_SOURCE
#define _SVID_SOURCE
#define _GNU_SOURCE             /* See feature_test_macros(7) */
#ifdef __APPLE__
#define _DARWIN_C_SOURCE
#include <MacTypes.h>
#define _BSD_SOURCE
#define __DARWIN_C_SOURCE
//#define DT_DIR 4
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
#include <limits.h>      /* PATH_MAX */
#ifdef __linux__
#include<linux/limits.h>
#endif
#include <ftw.h>
/* By default, print all messages of severity info and above.  */
#ifdef _DEBUG
static int global_debug = 3;
#else
static int global_debug = 2;
#endif
/* Name of this program */
static char *global_progname = "simple_posix";

#define dprintf if (global_debug >= 3) \
        fprintf (stdout, "%s: debug: (%d) ", global_progname, __LINE__), \
        fprintf

#define eprintf if (global_debug >= 0) \
        fprintf (stderr, "%s: error: (%d) ", global_progname, __LINE__), \
        fprintf


// testing for non-GNU systems #undef versionsort

#define MAX_CHAR_FILENAME 256
#define MAKEDIR makedir_
#define REMOVEDIR removedir_
#define GET_FILE_LIST get_file_list_
#define GET_FILE_LIST_MODIFIED get_file_list_modified_
#define LIST_DIRS list_dirs_
#define GLOB_FILE_LIST glob_file_list_
#define GLOB_RM_ALL glob_rm_all_
#define GET_ABSOLUTE_PATHNAME get_absolute_pathname_
#define FCOPY fcopy_
#define FCOPY_MMAP fcopy_mmap
#define FCOPY_SENDFILE fcopy_sendfile
#define FREE_FILE_LIST free_file_list_
#define SUBPROCESS subprocess_
#define WAIT_PID wait_pid_
#define TOUCH touch_
#define GET_SYSINFO get_sysinfo_
#define RMDIRECTORIES rmDirectories_
#define REMOVE_DIRECTORY remove_directory_
#define RECURSIVE_DELETE recursive_delete_

//#define __POSIX_PURE__
/**
 *  \brief Protect fortran strings by appending '\0' to end
 *
 *  \param str Fortran string
 *  \return C-style string
 */
char *F90toCstring(char *str, int len)
{
    char *res; /* C arrays are from 0:len-1 */
    if(len > PATH_MAX) {
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
    dprintf(stderr, "DEBUG: In simple_posix.c::isdir pathname :%s:\n",  cpathname);

    struct stat ibuf;
    int i = stat(cpathname, &ibuf);
    free(cpathname);
    if(i != 0) {

        dprintf(stderr, "DEBUG: In simple_posix.c::isdir pathname :%s: does not exist\n",  cpathname);
        //printf("%d %s\nisdir failed to create stat struct %s\n", errno, strerror(errno), cpathname);
        //perror("Failed : simple_posix.c::isdir ");
        return 0;
    }
    dprintf(stderr, "DEBUG: In isdir stat mode %ud\t isdir %d\n", ibuf.st_mode, S_ISDIR(ibuf.st_mode));
    dprintf(stderr, "DEBUG: In isdir stat mode s_ifdir %d\n" , (ibuf.st_mode & S_IFMT) == S_IFDIR);

    return S_ISDIR(ibuf.st_mode) ? 1 : 0;
}

// Recursive make directory
int makedir(char *path, size_t ivf_path)
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

    const size_t len = strlen(path);
    char _path[PATH_MAX];
    char *p;

    errno = 0;

    /* Copy string so its mutable */
    if(len > sizeof(_path) - 1) {
        errno = ENAMETOOLONG;
        return -1;
    }
    strcpy(_path, path);

    /* Iterate the string */
    for(p = _path + 1; *p; p++) {
        if(*p == '/') {
            /* Temporarily truncate */
            *p = '\0';

            if(mkdir(_path, S_IRWXU) != 0) {
                if(errno != EEXIST) {
                    fprintf(stderr, "mkdir %s\nerrno:%d msg:%s\n", _path, errno, strerror(errno));
                    perror("Failed : rmdir in simple_posix::makedir");
                    return -1;
                }
            }

            *p = '/';
        }
    }
    if(mkdir(_path, S_IRWXU) != 0) {
        if(errno != EEXIST) {
            fprintf(stderr, "mkdir %s\nerrno:%d msg:%s\n", _path, errno, strerror(errno));
            perror("Failed : rmdir in simple_posix::makedir");
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
    char _path[PATH_MAX];
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

/* let us make a recursive function to print the content of a given folder */
void show_dir_content_recursive(char * path, size_t ivf_path)
{
    char *cpath = F90toCstring(path, strlen(path));
    if(cpath == NULL) {
        printf("%d %s\n show_dir_content_recursive failed to convert string (unprotected) %s\n", errno, strerror(errno), path);
        perror("Failed : simple_posix.c:: show_dir_content_recursive ");
        return;
    }
    DIR * d = opendir(cpath); // open the path
    free(cpath);
    if(d == NULL) return; // if was not able return
    struct dirent * dir; // for the directory entries
    while((dir = readdir(d)) != NULL) { // if we were able to read somehting from the directory
        if(dir-> d_type != DT_DIR) // if the type is not directory just print it with blue
            printf("%s%s\n", "\x1B[34m", dir->d_name);
        else if(dir -> d_type == DT_DIR && strcmp(dir->d_name, ".") != 0 &&
                strcmp(dir->d_name, "..") != 0) {
            // if it is a directory
            printf("%s%s\n", "\x1B[32m", dir->d_name); // print its name in green
            char d_path[255]; // here I am using sprintf which is safer than strcat
            sprintf(d_path, "%s/%s", path, dir->d_name);
            show_dir_content_recursive(d_path, ivf_path); // recall with the new path
        }
    }
    closedir(d); // finally close the directory
}


/**
 *  \brief get_file_list emulates terminal command 'ls '
 *
 *
 *  \param path Directory string
 *  \return success status
 */
int get_file_list(char * path, char * ext, int *count, size_t ivf_path)
{
    char *cpath = F90toCstring(path, strlen(path));
    if(cpath == NULL) {
        printf("%d %s\nget_file_list failed to convert string (unprotected) %s\n", errno, strerror(errno), path);
        perror("Failed : simple_posix.c::get_file_list ");
        return -1;
    }

    dprintf(stderr, "DEBUG: In get_file_list, %15s:%s\n", "path", cpath);
    dprintf(stderr, "DEBUG: In get_file_list, %15s:%s\n", "ext", ext);


    DIR           *d;

    //  char ** results,

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
        /*option 1*/
        //results = malloc ( sizeof(char*) * fcount );
        dprintf(stderr, " get_file_list size %d\n", fcount);
        /*option 2*/
        FILE* f = fopen("__simple_filelist__", "w");
        fcount = 0;
        char *cext = F90toCstring(ext, 3);
        if(strcmp(ext, "")) {
            while((elem = readdir(d)) != NULL) {
                if(elem->d_type != DT_DIR) {
                    /*option 1*/
                    //        results[fcount] = malloc ( sizeof(char) * (elem->d_reclen+1));
                    //        strncpy(results[fcount], elem->d_name, elem->d_reclen );
                    /*option 2*/
                    fprintf(f, "%s\n", elem->d_name);
                    // printf(" get_file_list elem %d %s\n", fcount, elem->d_name);
                    fcount++;
                }
            }
        } else {
            while((elem = readdir(d)) != NULL) {
                if(elem->d_type != DT_DIR) {
                    if(strstr(elem->d_name, ext) != NULL) {
                        fprintf(f, "%s\n", elem->d_name);
                        //  printf(" get_file_list elem %d %s\n", fcount, elem->d_name);
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

/**
 *  \brief list files in directory 'path' , emulates 'ls -tr' if flags are set
 *
 *  \param path Fortran string input path
 *  \param ext filename extension
 *  \param count return ptr for number of files found
 *  \param flag dictates [reverse order] | [last modified time or alphanumeric sort]
 *  \return return status success=0
 */
int get_file_list_modified(char * path, char* ext, int* count, int flag, size_t ivf_path)
{
    /*equivalent to 'ls -tr path' */
    struct dirent **namelist;
    int n, nscandir, reverse = 0, mtime_alpha = 0;
    *count = 0;
    n=0;
    mtime_alpha = (flag >> 1) & 1;
    reverse = (flag) & 1;
    ssize_t ri = strlen(path);
    dprintf(stderr, "DEBUG:In get_file_list_modified %15s:%zd\n", "strlen(path)", ri);
    if(ri > 256) {
        ri = index(path, ' ') - path;
    }
    char *cpath = F90toCstring(path, ri);
    if(cpath == NULL) {
        printf("%d %s\nget_file_list_modified failed to convert string (unprotected) %s\n", errno, strerror(errno), path);
        perror("Failed : simple_posix.c::get_file_list_reverse_modified ");
        return -1;
    }

    dprintf(stderr, "DEBUG: In get_file_list_modified %15s:%s\n", "path", path);
    dprintf(stderr, "DEBUG: In get_file_list_modified %15s:%zd\n", " strlen(path)", strlen(path));
    dprintf(stderr, "DEBUG: In get_file_list_modified %15s:%d  %15s:%3s\n", "flag", flag, "ext", ext);
    dprintf(stderr, "DEBUG: In get_file_list_modified %15s:%zd\n", "strlen ext", strlen(ext));
    // dprintf(stderr, "DEBUG: In get_file_list_modified %15s:%12lu\n", "ivf_pathLen", (unsigned long) ivf_path);

    // Remove temp file if it exists
    if(access("./__simple_filelist__", F_OK) != -1) {
        if(remove("./__simple_filelist__") != 0) {
            printf("%d %s\nget_file_list_modified failed to remove __simple_filelist__ in path %s\n", errno, strerror(errno), cpath);
            perror("Failed : simple_posix.c::get_file_list_reverse_modified ");
            return -1;
        }
    }
    if(mtime_alpha) { /* modified time sort */
#if defined(versionsort)
        nscandir = scandir(cpath, &namelist, NULL, versionsort);
#else
        nscandir = scandir(cpath, &namelist, NULL, alphasort);
#endif
    } else { /* alphanumeric sort */
        nscandir = scandir(cpath, &namelist, NULL, alphasort);
    }

    free(cpath);
    FILE* f = fopen("__simple_filelist__", "w");
    if(!f) {
        printf("%d %s\nget_file_list_modified failed to open __simple_filelist__ in %s\n", errno, strerror(errno), cpath);
        perror("Failed : simple_posix.c::get_file_list_modified ");
        return -1;
    }
    if(nscandir < 0)
        perror("scandir failed to find files ; get_file_list_modified");
    else {

        int i = nscandir; // for reverse cases
        n = 0; // start count
        if(strcmp(ext, "") == 0) { /* Ignoring extension */
            if(reverse) {
                while(i--) {
                    if(namelist[i]->d_type != DT_DIR) {
                        fprintf(f, "%s\n", namelist[i]->d_name);
                        n++;
                    }
                    free(namelist[i]);
                }
            } else {
                for(i = 0; i < nscandir; i++) {
                    if(namelist[i]->d_type != DT_DIR) {
                        fprintf(f, "%s\n", namelist[i]->d_name);
                        n++;
                    }
                    free(namelist[i]);
                }
            }

        } else {    /*  Compare with extension */
            char *cext = F90toCstring(ext, 3);
            if(reverse) {
                while(i--) {
                    if((namelist[i]->d_type != DT_DIR) && (strstr(namelist[i]->d_name, cext) != NULL)) {
                        fprintf(f, "%s\n", namelist[i]->d_name);
                        n++;
                    }
                    free(namelist[i]);
                }
            } else {
                for(i = 0; i < nscandir; i++) {
                    if((namelist[i]->d_type != DT_DIR) && (strstr(namelist[i]->d_name, cext) != NULL)) {
                        fprintf(f, "%s\n", namelist[i]->d_name);
                        n++;
                    }
                    free(namelist[i]);
                }
            }
            free(cext);
        }

        free(namelist);
    }
    *count = n;
    fclose(f);
    return 0;
}

/**
 *  \brief list directories in directory 'path'
 *
 *
 *  \param path Fortran string input path
 *  \param count return ptr for number of dirss found
 *  \return return status success=0
 */
int list_dirs(char * path, int* count, size_t ivf_path)
{
    char *cpath = F90toCstring(path, strlen(path));
    if(cpath == NULL) {
        printf("%d %s\nlist_dirs failed to open convert string (unprotected) %s\n", errno, strerror(errno), path);
        perror("Failed : simple_posix.c::list_dirs ");
        return -1;
    }
    extern int errno;
    /*option 2*/
    //  char ** results;
    DIR           *d;
    int fcount = 0;
    d = opendir(cpath);
    free(cpath);
    if(d) {
        // // Prescan
        // struct dirent *dir;
        //    while((dir = readdir(d)) != NULL) {
        //  if(dir-> d_type == DT_DIR && strcmp(dir->d_name, ".") != 0 && strcmp(dir->d_name, "..") != 0) {
        //    fcount++;
        //  }
        //}
        // rewinddir(d);
        // //Allocate string array for Fortran
        //  results = malloc ( sizeof(char*) * fcount );
        //printf(" list_dirs size %d\n", fcount);
        //if (fcount > 0){
        FILE* f = fopen("__simple_filelist__", "w");
        if(f) {
            fcount = 0;
            struct dirent *dir;
            while((dir = readdir(d)) != NULL) {
                if(dir-> d_type == DT_DIR && strcmp(dir->d_name, ".") != 0 && strcmp(dir->d_name, "..") != 0) {
                    /*option 1*/
                    // results[fcount] = malloc ( sizeof(char) * (dir->d_reclen+1));
                    // strncpy(results[fcount], dir->d_name, (dir->d_reclen+1) );
                    /*option 2*/
                    fprintf(f, "%s\n", dir->d_name);
                    //              printf(" list_dirs elem %d %s\n", fcount, dir->d_name);
                    fcount++;
                }
            }
            fclose(f);
        } else {
            printf("%d %s\nlist_dirs opendir failed to open temp file\n", errno, strerror(errno));
            perror("Failed : simple_posix.c::list_dirs ");
            return -1;
        }
        //}
        closedir(d);
    } else {
        printf("%d %s\nlist_dirs opendir failed to open %s\n", errno, strerror(errno), path);
        perror("Failed : simple_posix.c::list_dirs ");
        return -1;
    }
    *count = fcount;
    return 0;
}

int qs_struct(struct dirent **namelist, int left, int right)
{

    register int i, j;

    struct dirent *temp;
    struct stat ibuf;
    struct stat jbuf;
    struct stat xbuf;

    i = left; j = right;
    stat(namelist[i]->d_name, &ibuf);
    stat(namelist[j]->d_name, &jbuf);
    stat(namelist[(left + right) / 2]->d_name, &xbuf);

    do {
        while((ibuf.st_mtime < xbuf.st_mtime) && (i < right)) {
            i++;
            stat(namelist[i]->d_name, &ibuf);
        }
        while((jbuf.st_mtime > xbuf.st_mtime) && (j > left))  {
            j--;
            stat(namelist[j]->d_name, &jbuf);
        }
        if(i <= j) {
            temp = namelist[i];
            namelist[i] = namelist[j];
            namelist[j] = temp;
            i++; j--;
        }

    } while(i <= j);

    if(left < j) qs_struct(namelist, left, j);
    if(i < right) qs_struct(namelist, i, right);

    return 0;
}

int qs_glob(glob_t*globlist, int left, int right)
{
    register int i, j;
    char *temp;
    struct stat ibuf;
    struct stat jbuf;
    struct stat xbuf;

    i = left; j = right;
    stat(globlist->gl_pathv[i], &ibuf);
    stat(globlist->gl_pathv[j], &jbuf);
    stat(globlist->gl_pathv[(left + right) / 2], &xbuf);
    do {
        while((ibuf.st_mtime < xbuf.st_mtime) && (i < right)) {
            i++;
            stat(globlist->gl_pathv[i], &ibuf);
        }
        while((jbuf.st_mtime > xbuf.st_mtime) && (j > left))  {
            j--;
            stat(globlist->gl_pathv[j], &jbuf);
        }
        if(i <= j) {
            temp = malloc(strlen(globlist->gl_pathv[i]) + 1);
            strcpy(globlist->gl_pathv[i], temp);
            globlist->gl_pathv[i] = realloc(globlist->gl_pathv[i], strlen(globlist->gl_pathv[j]) + 1);
            strcpy(globlist->gl_pathv[j], globlist->gl_pathv[i]);
            globlist->gl_pathv[j] = realloc(globlist->gl_pathv[j], strlen(temp) + 1);
            strcpy(temp                 , globlist->gl_pathv[j]);
            i++; j--;
            free(temp);
        }
    } while(i <= j);

    if(left < j) qs_glob(globlist, left, j);
    if(i < right) qs_glob(globlist, i, right);
    return 0;
}

/* A Quicksort for structures of type glob. */
void quicksort_glob(glob_t*globlist, int start, int count)
{
    qs_glob(globlist, start, count);
}

/**
 *  \brief list files in relative directory using GLOB
 *
 *  \param match Fortran string input glob
 *  \param count return ptr for number of files found
 *  \param sort_by_time flag that dictates last modified time or alphanumeric sorting
 *  \return return status success=0
 */
int glob_file_list(char *match,  int*count, int* sort_by_time, size_t ivf_match)
{
    *count=0;
    /* Check input char pointer */
    if(match == NULL || strlen(match) <= 0) {
        fprintf(stderr, "simple_posix.c::glob_file_list match string null .. Failing\n");
        return -1;
    }
    extern int errno;
    char *cmatch = F90toCstring(match, strlen(match));
    if(cmatch == NULL) {
        printf("%d %s\n glob_file_list failed to convert string (unprotected) %s\n", errno, strerror(errno), match);
        perror("Failed : simple_posix.c::glob_file_list 1");
        return -1;
    }
    glob_t        globlist;
    int err,n;
    err = 0;n=0;

    dprintf(stderr, "DEBUG: In glob_file_list size cmatch:%zd size match:%zd\n", strlen(cmatch), strlen(match));
    dprintf(stderr, "DEBUG: In glob_file_list cmatch:%s\n", cmatch);
    dprintf(stderr, "DEBUG: In glob_file_list flag:%d \n", *sort_by_time);


    int glob_val = glob(cmatch, 0, NULL, &globlist);
    free(cmatch);
    // Remove temp file if it exists
    if(access("./__simple_filelist__", F_OK) != -1) {
        if(remove("./__simple_filelist__") != 0) {
            printf("%d %s\nglob_file_list failed to remove __simple_filelist__ in glob %s\n", errno, strerror(errno), cmatch);
            perror("Failed : simple_posix.c::glob_file_list ");
            return -1;
        }
    }

    FILE* f = fopen("__simple_filelist__", "w");
    if(!f) {
        printf("%d %s\nglob_file_list failed to open __simple_filelist__\n", errno, strerror(errno));
        perror("Failed : simple_posix.c::glob_file_list 2");
        return -1;
    }

    if(glob_val == GLOB_NOSPACE) {
        fprintf(stderr, "simple_posix.c::glob_file_list NOSPACE no memory left \n");
        err = -2;
    } else if(glob_val == GLOB_NOMATCH) {
        fprintf(stderr, "simple_posix.c::glob_file_list NOMATCH  glob:%s\n", cmatch);
        err = 0; // Do nothing
    } else if(glob_val == GLOB_ABORTED) {
        fprintf(stderr, "simple_posix.c::glob_file_list GLOB_ABORTED\n");
        err = -1;
    } else {
        n = (int) globlist.gl_pathc; //account for __simple_filelist__
        printf(" glob_file_list size before trimming %d \n", *count);
        if(n > 0) {
            if(*sort_by_time) quicksort_glob(&globlist, 0, n);
            n = 0;
            int i = 0;
            while(globlist.gl_pathv[i]) {
                //          fprintf(stderr, "GLOB>> %s\n", globlist.gl_pathv[i]);
                if(strstr(globlist.gl_pathv[i], "__simple_filelist__") == NULL) {

                    fprintf(f, "%s\n", globlist.gl_pathv[i]);
                    n = n + 1;
                }
                i++;
            }
            /*
            #if 0
            (*results)[*count] = malloc(sizeof(char) * (*count));
            for(i = 0; i < *count; i++) {
                int lenstr = strlen(globlist.gl_pathv[i]) + 1;
                results[i] = malloc(sizeof(char) * lenstr);
                strncpy(results[i], globlist.gl_pathv[i], lenstr - 1);
            }
            printf(" glob_file_list allocated \n");

            i = 0;
            while(globlist.gl_pathv[i]) {
                printf("%s\n", globlist.gl_pathv[i]);
                i++;
            }
            printf(" glob_file_list freeing globlist struct \n");
            #endif
            */
        }
        dprintf(stderr, " glob_file_list size after trimming %d \n", *count);
    }
    globfree(&globlist);
    fclose(f);
    *count=n;
    return err;
}

void free_file_list(char**ptr, int n)
{
    for(int i = 0; i < n; i++) free(ptr[i]);
    free(ptr);
}

int subprocess(char* args, int* arglen, size_t ivf_cmd)
{
    char * cargs = F90toCstring(args, *arglen);

    char * const argsC[] = { "-l", "-c", cargs, "\r", NULL };

    dprintf(stderr, "DEBUG:In subprocess %30s:%s:\n", " args ", argsC[2]);
    dprintf(stderr, "DEBUG:In subprocess %30s:%zd  %d\n", "strlen args", strlen(args), *arglen);

    int pid = fork();
    if(pid == 0) {
        execvp("/bin/bash", argsC);
    }
    free(cargs);
    return pid;
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




//#ifdef _POSIX_MAPPED_FILES
int fcopy_mmap(char* file1,  int*len1, char* file2, int*len2)
{

    int sfd, dfd;
    dprintf(stderr, "MMAP FCOPY In fcopy mmap/memcpy\n");
    dprintf(stderr, "DEBUG: In fcopy mmap file1:%s \n size :%zu\t%d\n", file1, strlen(file1), *len1);
    dprintf(stderr, "DEBUG: In fcopy mmap file2:%s \n size :%zu\t%d\n", file2, strlen(file2), *len2);
    size_t filesize;
    extern int errno;
    char *buffer1, *buffer2;
    char *src, *dest;
    int flag = 0;
    if(!(buffer1 = F90toCstring(file1, *len1))) {
        printf("%d %s\nfcopy failed to interpret file str %30s\n", errno, strerror(errno), file1);
        perror("Failed : F90toCstring (buf) in simple_posix.c:fcopy");
        free(buffer1);
        flag = 1;
    } else if(!(buffer2 = F90toCstring(file2, *len2))) {
        printf("%d %s\nfcopy failed to interpret file str %30s\n", errno, strerror(errno), file2);
        perror("Failed : F90toCstring (buf) in simple_posix::fcopy");
        free(buffer2);
        flag = 1;
    } else {
        dprintf(stderr, "In fcopy opening files\n");
        /* SOURCE */
        sfd = open(buffer1, O_RDONLY);
        if(sfd) {
            filesize = lseek(sfd, 0, SEEK_END);
            src = mmap(NULL, filesize, PROT_READ, MAP_PRIVATE, sfd, 0);
            /* DESTINATION */
            dfd = open(buffer2, O_RDWR | O_CREAT, 0666);
            if(dfd) {
                if(ftruncate(dfd, filesize) != 0) {
                    printf("%d %s\nfcopy failed to truncate  file %s\n", errno, strerror(errno), buffer2);
                    perror("Failed : simple_posix.c::fcopy ftruncate\n");
                    return -1;

                }

                dest = mmap(NULL, filesize, PROT_READ | PROT_WRITE, MAP_SHARED, dfd, 0);

                /* COPY */

                memcpy(dest, src, filesize);

                munmap(src, filesize);
                munmap(dest, filesize);

                close(dfd);
                close(sfd);
            } else {
                printf("%d %s\nfcopy failed to open  file %s\n", errno, strerror(errno), buffer2);
                perror("Failed : simple_posix.c::fcopy opening output file\n");
                flag = 1;
                close(dfd);
            }
        } else {
            printf("%d %s\nfcopy failed to open file %s\n", errno, strerror(errno), buffer1);
            perror("Failed : simple_posix.c::fcopy ");
            flag = 1;
            close(sfd);
        }
        free(buffer1);
        free(buffer2);
    }
    return flag;
}


// Copy file in chunks of 1k bytes
int fcopy(char* file1,  int*len1, char* file2, int*len2, size_t ivf_file1, size_t ivf_file2)
{
    extern int errno;
    fprintf(stderr, "DEBUG: In fcopy file1:%s \n size :%zu\t%d\t%zu\n", file1, strlen(file1), *len1, ivf_file1);
    fprintf(stderr, "DEBUG: In fcopy file2:%s \n size :%zu\t%d\t%zu\n", file2, strlen(file2), *len2, ivf_file2);
    char *buffer1, *buffer2;
    int flag = 0;
    if(!(buffer1 = F90toCstring(file1, *len1))) {
        printf("%d %s\nfcopy failed to interpret file str %30s\n", errno, strerror(errno), file1);
        perror("Failed : F90toCstring (buf) in simple_posix.c:fcopy");
        free(buffer1);
        flag = 1;
    } else if(!(buffer2 = F90toCstring(file2, *len2))) {
        printf("%d %s\nfcopy failed to interpret file str %30s\n", errno, strerror(errno), file2);
        perror("Failed : F90toCstring (buf) in simple_posix::fcopy");
        free(buffer2);
        flag = 1;
    } else {
        dprintf(stderr, "In fcopy opening files\n");
        FILE *f1, *f2;
        f1 = fopen(buffer1, "r");
        if(f1) {
            f2 = fopen(buffer2, "w");
            if(f2) {
                char            buffer[1024];
                size_t          n;
                while((n = fread(buffer, sizeof(char), sizeof(buffer), f1)) > 0) {
                    if(fwrite(buffer, sizeof(char), n, f2) != n) {
                        flag = 1;
                        break;
                    }
                }
                fclose(f2);
            } else {
                printf("%d %s\nfcopy failed to open  file %s\n", errno, strerror(errno), buffer2);
                perror("Failed : simple_posix.c::fcopy opening output file\n");
                flag = 1;
            }
            fclose(f1);
        } else {
            printf("%d %s\nfcopy failed to open file %s\n", errno, strerror(errno), buffer1);
            perror("Failed : simple_posix.c::fcopy ");
            flag = 1;
        }
        free(buffer1); free(buffer2);
    }
    return flag;
}


#ifdef __APPLE__
#include <sys/socket.h>
#include <sys/uio.h>
int fcopy_sendfile(char* file1,  int*len1, char* file2, int*len2, size_t ivf_file1, size_t ivf_file2)
{ return fcopy(file1,len1, file2,len2, ivf_file1, ivf_file2);
}
#else
#include <sys/sendfile.h>


int fcopy_sendfile(char* file1,  int*len1, char* file2, int*len2, size_t ivf_file1, size_t ivf_file2)
{
    dprintf(stderr, "In fcopy2\n");
    dprintf(stderr, "DEBUG: In fcopy2 file1:%s \n size :%zu\t%d\t%zu\n", file1, strlen(file1), *len1, ivf_file1);
    dprintf(stderr, "DEBUG: In fcopy2 file2:%s \n size :%zu\t%d\t%zu\n", file2, strlen(file2), *len2, ivf_file2);
    extern int errno;
    ssize_t sent = 0;
    char *buffer1, *buffer2;
    int flag = 0;
    if(!(buffer1 = F90toCstring(file1, *len1))) {
        printf("%d %s\nfcopy2 failed to interpret file str %30s\n", errno, strerror(errno), file1);
        perror("Failed : F90toCstring (buf) in simple_posix.c:fcopy2");
        free(buffer1);
        flag = 1;
    } else if(!(buffer2 = F90toCstring(file2, *len2))) {
        printf("%d %s\nfcopy2 failed to interpret file str %30s\n", errno, strerror(errno), file2);
        perror("Failed : F90toCstring (buf) in simple_posix::fcopy2");
        free(buffer2);
        flag = 1;
    } else {
        struct stat st;
        if(0 == stat(buffer1, &st)) {
            int INPUTsize = st.st_size;
	    //            off_t OSX_INPUTsize = (off_t) INPUTsize;
            dprintf(stderr, "In fcopy2 opening files\n");
            FILE *fin = fopen(buffer1, "r");
            if(fin) {
                FILE *fout = fopen(buffer2, "w");
                if(fout) {
                    fflush(fout);
                    off_t offset;

                    //        fseek(fin, 0L, SEEK_END);
                    //size_t toCopy = ftell(fin);
                    //rewind(fin); //fseek(fin, 0L, SEEK_SET);

                    offset = ftello(fin);
                    sent = sendfile(fileno(fout), fileno(fin), &offset, INPUTsize);
                    fseeko(fin, offset, SEEK_SET);
                    if(sent != INPUTsize) {
                        printf("%d %s\nfcopy2 failed to copy file %s to %s\n", errno, strerror(errno), buffer1, buffer2);
                        perror("Failed : simple_posix.c::fcopy2 copying file\n");
                        flag = 1;
                    }
                    fclose(fout);
                } else {
                    printf("%d %s\nfcopy2 failed to open  file %s\n", errno, strerror(errno), buffer2);
                    perror("Failed : simple_posix.c::fcopy2 opening output file\n");
                    flag = 1;
                }
                fclose(fin);
            } else {
                printf("%d %s\nfcopy2 failed to open file %s\n", errno, strerror(errno), buffer1);
                perror("Failed : simple_posix.c::fcopy2 ");
                flag = 1;
            }
        } else {
            printf("%d %s\nfcopy2 failed to stat file %s\n", errno, strerror(errno), buffer1);
            perror("Failed : simple_posix.c::fcopy2 ");
            flag = 1;
        }
        free(buffer1);
        free(buffer2);
    }
    return flag;
}
#endif

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


/**
 * returns a null-terminated string containing the canonicalized absolute pathname corresponding to path.
 */
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
    dprintf(stderr, "DEBUG:In get_absolute_pathname %30s:%s:\n", "input", filein);
    dprintf(stderr, "DEBUG:In get_absolute_pathname %30s:%zd\n", "strlen(input path)", strlen(filein));
    dprintf(stderr, "DEBUG:%30s: resolved:%s\n", "DEBUG: In  get_absolute_pathname", resolved);


    if(resolved == NULL) {
        printf("%d %s\nget_absolute_path failed to canonicalize  file %s\n", errno, strerror(errno), filein);
        perror("Failed : simple_posix.c::get_absolute_pathname ");
        return 1;
    }
    *outlen = strlen(resolved) ;
    strncpy(out, resolved, *outlen); //for(int i = *outlen; i < MAX_CHAR_FILENAME; i++) out[i] = '\0';
    out[*outlen] = '\0';
    dprintf(stderr, "DEBUG:In get_absolute_pathname %30s:%s:\n", " out path", out);
    dprintf(stderr, "DEBUG:%30s: strlen out path:%zd\n", "DEBUG: In  get_absolute_pathname", strlen(out));

    c2fstr(resolved, out, *outlen, sizeof(resolved));
    out[0] = '/';

    dprintf(stderr, "DEBUG: In get_absolute_pathname c2fstr: %s\n", out);
    dprintf(stderr, "DEBUG: In get_absolute_pathname %30s:%d\n", " out path addr", *out);
    dprintf(stderr, "DEBUG: In get_absolute_pathname strlen(out):%zd\n", strlen(out));

    free(filein);
    free(resolved);
    //out = resolved;
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
    long long free_memory;
    long long used_memory;
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


    *totalram =  physical_memory;          /* Total usable main memory size (bytes)*/
    *sharedram = t_info.virtual_size;     /* Amount of shared memory */
    *bufferram = t_info.resident_size;    /* Memory used by buffers */
    *totalhigh = free_memory + used_memory; /* Total high water mark memory size */
    *HWMusage =  used_memory;             /* high memory size used */

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

// int rmrf(char *path)
// {
//     return nftw(path, unlink_cb, 64, FTW_DEPTH | FTW_PHYS);
// }


static int rmFiles(const char *pathname, const struct stat *sbuf, int type, struct FTW *ftwb)
{
    if(remove(pathname) < 0) {
        perror("ERROR: remove");
        return -1;
    }
    return 0;
}

// int rmDirectories(char* path, int *len, int*count)
// {
//     char *cpath = F90toCstring(path, *len);
//     *count = 0;
//     if(nftw(cpath, rmFiles, 10, FTW_DEPTH | FTW_MOUNT | FTW_PHYS) < 0) {
//         perror("ERROR: ntfw");
//         exit(1);
//     }
//     *count = *count + 1;
//     FILE* f = fopen("__simple_filelist__", "a");
//     fprintf(f, "%s\n", cpath);
//     fclose(f);
//     free(cpath);
//     return 0;
// }


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

int recursive_delete(char *dir, int*len, int*count)
{
    int ret = 0;
    FTS *ftsp = NULL;
    extern int errno;
    FTSENT *curr;
    // Cast needed (in C) because fts_open() takes a "char * const *", instead
    // of a "const char * const *", which is only allowed in C++. fts_open()
    // does not modify the argument.
    char *files[] = { (char *) F90toCstring(dir, *len), NULL
                    };
    // FTS_NOCHDIR  - Avoid changing cwd, which could cause unexpected behavior
    //                in multithreaded programs
    // FTS_PHYSICAL - Don't follow symlinks. Prevents deletion of files outside
    //                of the specified directory
    // FTS_XDEV     - Don't cross filesystem boundaries
    ftsp = fts_open(files, FTS_NOCHDIR | FTS_PHYSICAL | FTS_XDEV, NULL);
    free(files[0]);
    if(!ftsp) {
        fprintf(stderr, "simple_posix.c:recursive_delete()  fts_open failed :%s\nERROR %d %s\n", dir, errno, strerror(errno));
        perror("simple_posix.c::recursive_delete fts_open failed");

        /* fts_open errors
           [EACCES] Search permission is denied for any component of
         path or read permission is denied for path, or fn returns -1 and does not reset
         errno.

         [ELOOP] A loop exists in symbolic links encountered during resolution of the
         path argument.

         [ENAMETOOLONG] The length of a component of a pathname is longer than
         {NAME_MAX}.

         [ENOENT] A component of path does not name an existing file or path is an empty
         string.

         [ENOTDIR] A component of path names an existing file that is neither a
         directory nor a symbolic link to a directory.

         [EOVERFLOW] A field in the stat structure cannot be represented correctly in
         the current programming environment for one or more files found in the file
         hierarchy. */
        ret = -1;
        return ret;
    }


    while((curr = fts_read(ftsp))) {
        switch(curr->fts_info) {
        case FTS_NS:
        case FTS_DNR:
        case FTS_ERR:
            fprintf(stderr, "%s: fts_read error: %s\n",
                    curr->fts_accpath, strerror(curr->fts_errno));
            break;
        case FTS_DC:
        case FTS_DOT:
        case FTS_NSOK:
            // Not reached unless FTS_LOGICAL, FTS_SEEDOT, or FTS_NOSTAT were
            // passed to fts_open()
            break;
        case FTS_D:
            // Do nothing. Need depth-first search, so directories are deleted
            // in FTS_DP
            break;
        case FTS_DP:
        case FTS_F:
        case FTS_SL:
        case FTS_SLNONE:
        case FTS_DEFAULT:
            if(remove(curr->fts_accpath) < 0) {
                fprintf(stderr, "%s: Failed to remove: %s\n",
                        curr->fts_path, strerror(errno));
                ret = -1;
            } else {
                *count = *count + 1;
                FILE* f = fopen("__simple_filelist__", "a");
                fprintf(f, "%s\n", curr->fts_accpath);
                fclose(f);
            }
            break;
        }
    }
    if(ftsp) {
        fts_close(ftsp);
    }
    return ret;
}

/**
 *  \brief glob_rm_all removes dirs/files in relative directory using GLOB
 *
 *  \param match Fortran string input glob
 *  \param count return ptr for number of files/dirs found
 *  \return return status success=0
 */
int glob_rm_all(char *match,  int*count, size_t ivf_match)
{
    /* Check input char pointer */
    if(match == NULL || strlen(match) <= 0) {
        fprintf(stderr, "simple_posix.c::glob_rm_all match string null .. Failing\n");
        return -1;
    }
    extern int errno;
    char *cmatch = F90toCstring(match, strlen(match));
    if(cmatch == NULL) {
        printf("%d %s\n glob_rm_all failed to convert string (unprotected) %s\n", errno, strerror(errno), match);
        perror("Failed : simple_posix.c::glob_rm_all 1");
        return -1;
    }
    glob_t        globlist;
    FILE* f;
    int err;
    err = 0;
    dprintf(stderr, "DEBUG: In glob_rm_all size cmatch:%zd size match:%zd\n", strlen(cmatch), strlen(match));
    dprintf(stderr, "DEBUG: In glob_rm_all cmatch:%s\n", cmatch);

    int glob_val = glob(cmatch, 0, NULL, &globlist);
    free(cmatch);
    f = fopen("__simple_filelist__", "w");
    if(!f) {
        printf("%d %s\nglob_rm_all failed to open __simple_filelist__\n", errno, strerror(errno));
        perror("Failed : simple_posix.c::glob_rm_all 2");
        return -1;
    }
    fclose(f);


    if(glob_val == GLOB_NOSPACE) {
        fprintf(stderr, "simple_posix.c::glob_rm_all NOSPACE no memory left \n");
        err = -2;
    } else if(glob_val == GLOB_NOMATCH) {
        fprintf(stderr, "simple_posix.c::glob_rm_all NOMATCH  glob:%s\n", cmatch);
        err = 0; // Do nothing
    } else if(glob_val == GLOB_ABORTED) {
        fprintf(stderr, "simple_posix.c::glob_rm_all GLOB_ABORTED\n");
        err = -1;
    } else {
        *count = (int) globlist.gl_pathc; //account for __simple_filelist__
        dprintf(stderr, " glob_rm_all size before trimming %d \n", *count);

        *count = 0;
        int i = 0;
        while(globlist.gl_pathv[i]) {

            dprintf(stderr, "GLOB>> %s ... ", globlist.gl_pathv[i]);
            if(strstr(globlist.gl_pathv[i], "__simple_filelist__") == NULL) {
                if(strcmp(globlist.gl_pathv[i], ".") != 0 || strcmp(globlist.gl_pathv[i], "..") != 0) {
                    dprintf(stderr, " for deletion\n");

                    //fprintf(f, "%s\n", globlist.gl_pathv[i]);
                    *count = *count + 1;
                } else {
                    dprintf(stderr, "dot dir found, CANNOT DELETE CWD, stepping over\n");

                }
            } else {

                dprintf(stderr, "temp file found, stepping over\n");

            }
            i++;
        }

        dprintf(stderr, "Total glob items for deletion:  %d\n", *count);
        *count = 0;
        i = 0;
        while(globlist.gl_pathv[i]) {
            if(strstr(globlist.gl_pathv[i], "__simple_filelist__") == NULL) {
                if(strcmp(globlist.gl_pathv[i], ".") != 0 || strcmp(globlist.gl_pathv[i], "..") != 0) {
                    fprintf(stderr, "GLOB>> %s ... ", globlist.gl_pathv[i]);
                    struct stat s;
                    int staterr = stat(globlist.gl_pathv[i], &s);
                    if(staterr == 0) {
#if _DEBUG
                        if(S_ISDIR(s.st_mode)) {
                            fprintf(stderr, " dir ");
                        } else if(S_ISREG(s.st_mode)) {
                            fprintf(stderr, " file ");
                        } else if(S_ISLNK(s.st_mode)) {
                            fprintf(stderr, " link ");
                        } else {
                            fprintf(stderr, " unsupported mode ");
                        }
#endif
                        // remove is a stdlib function that uses unlink for files and rmdir for directories
                        if(remove(globlist.gl_pathv[i]) == 0) {

                            dprintf(stderr, " deleted\n");

                            *count = *count + 1;
                            FILE* f = fopen("__simple_filelist__", "a");
                            fprintf(f, "%s\n", globlist.gl_pathv[i]);
                            fclose(f);
                        } else {
                            if(errno == ENOTEMPTY) {
                                fprintf(stderr, " calling recursive delete\n");
                                int len2 = (int)strlen(globlist.gl_pathv[i]);
                                err = remove_directory_recursive(globlist.gl_pathv[i], &len2, count);
                            }
                            if(err != 0) {
                                fprintf(stderr, "\n%d %s\nglob_rm_all failed to remove %s\n", errno, strerror(errno), globlist.gl_pathv[i]);
                                perror("Failed : simple_posix.c::glob_rm_all ");
                            }

                        }
                    } else {
                        fprintf(stderr, "\n%d %s\nglob_rm_all failed to stat %s\n", errno, strerror(errno), globlist.gl_pathv[i]);
                        perror("Failed : simple_posix.c::glob_rm_all ");
                    }

                }
            }
            i++;
        }

        dprintf(stderr, " glob_rm_all size after trimming %d \n", *count);

    }
    globfree(&globlist);

    return err;
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
