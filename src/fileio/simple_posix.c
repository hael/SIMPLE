/**
 *   simple_posix.c
 *   \brief POSIX system calls for SIMPLE
 *
 *   OS file utilities to allow simple to interface to POSIX functions without using system calls
 *
 *   Michael Eager   2018
 */

#define _SVID_SOURCE
#define _GNU_SOURCE             /* See feature_test_macros(7) */
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>           /* Definition of AT_* constants */
#include <sys/wait.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <pwd.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <dirent.h>
#include <time.h>
#include <glob.h>

#ifdef __linux__
#include<linux/limits.h>
#endif

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
#define FREE_FILE_LIST free_file_list_
#define SUBPROCESS subprocess_

/**
 *  \brief Protect fortran strings by appending '\0' to end
 *
 *  \param str Fortran string
 *  \return C-style string
 */
char *F90toCstring(char *str, int len)
{
    char *res; /* C arrays are from 0:len-1 */
    if(len > MAX_CHAR_FILENAME) {
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
    for (; max_f > 0 && (isspace(f[max_f-1]) || f[max_f-1]=='\0'); max_f--);
    return max_f;
}
#if 0
void f2cstr(char *f, int max_f, char *c, int max_c)
{
    int i;
    i = min(fstrlen(f,max_f),max_c);
    strncpy(c,f,i);
    c[i]='\0';
}

void c2fstr(char *c, int max_c, char *f, int max_f)
{
    int i;
    i = min(strlen(c),max_f);
    strncpy(f,c,i);
    for( ; i<max_f; i++) f[i]=' ';
}
#endif

static char *c2fstr(char* cstr, char *fstr, int elem_len, int sizeofcstr)
{
  int i,j;
  /* elem_len includes \0 for C strings. Fortran strings don't have term. \0.
     Useful size of string must be the same in both languages. */
  for (i=0; i<sizeofcstr/elem_len; i++) {
    for (j=1; j<elem_len && *cstr; j++) *fstr++ = *cstr++;
    cstr += 1+elem_len-j;
    for (; j<elem_len; j++) *fstr++ = ' ';
  } /* f951 - Seems to be returning the original fstr. */
  return fstr-sizeofcstr+sizeofcstr/elem_len;
}

void f2cstr (const char* f_str, char* c_str, unsigned length)
{
  /* fprintf (stderr, "f2cstr in '%.*s'\n", length, f_str); */

  unsigned i = length-1;
  strncpy (c_str, f_str, length);

  while (f_str[i] == ' ') {
    c_str[i] = '\0';
    if (i == 0)
      break;
    i--;
  }
  /* fprintf (stderr, "f2cstr out '%s'\n", c_str); */
}



int makedir(const char *pathname, size_t ivf_pathname)
{
    mode_t mask = umask(0);
    umask(mask);
    return mkdir(pathname, mask);
}

/* let us make a recursive function to print the content of a given folder */
void show_dir_content_recursive(const char * path, size_t ivf_path)
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
        else if(dir -> d_type == DT_DIR && strcmp(dir->d_name, ".") != 0 && strcmp(dir->d_name, "..") != 0) { // if it is a directory
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
 *  Detailed description
 *
 *  \param param
 *  \return return type
 */
int get_file_list(const char * path,  const char * ext, int *count, size_t ivf_path)
{
    char *cpath = F90toCstring(path, strlen(path));
    if(cpath == NULL) {
        printf("%d %s\nget_file_list failed to convert string (unprotected) %s\n", errno, strerror(errno), path);
        perror("Failed : simple_posix.c::get_file_list ");
        return -1;
    }
    printf("DEBUG:  In get_file_list %s\n", cpath);
    extern int errno;
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
        printf(" get_file_list size %d\n", fcount);
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
 *  \param ext [optional] filename extension
 *  \param count return ptr for number of files found
 *  \param flag dictates [reverse order] | [last modified time or alphanumeric sort]
 *  \return return status success=0
 */
int get_file_list_modified(const char * path, const char* ext, int* count, int flag, size_t ivf_path)
{
    /*equivalent to 'ls -tr path' */
    struct dirent **namelist;
    int   n, reverse = 0, mtime_alpha = 0;
    *count = 0;
    mtime_alpha = (flag >> 1) & 1;
    reverse = (flag) & 1;
    size_t ri = strlen(path);
    printf("%30s: path:%zd\n", "DEBUG: In get_file_list_modified", ri);
    if(ri > 256) {
        ri = index(path, ' ') - path;
    }
    char *cpath = F90toCstring(path, ri);
    if(cpath == NULL) {
        printf("%d %s\nget_file_list_modified failed to convert string (unprotected) %s\n", errno, strerror(errno), path);
        perror("Failed : simple_posix.c::get_file_list_reverse_modified ");
        return -1;
    }
    printf("%30s: path:%s\n", "DEBUG: In get_file_list_modified", path);
    printf("%30s: strlen path:%zd\n", "DEBUG: In get_file_list_modified", strlen(path));
    printf("%30s:  flag:%d  ext:%3s\n", "DEBUG: In get_file_list_modified", flag, ext);
    printf("%30s: strlen ext:%zd\n", "DEBUG: In get_file_list_modified", strlen(ext));
    printf("%30s: ivf_pathLen:%12u\n", "DEBUG: In get_file_list_modified", (long) ivf_path);

    // Remove temp file if it exists
    if(access( "./__simple_filelist__", F_OK ) != -1){
        if(remove("./__simple_filelist__")!=0) {
          printf("%d %s\nget_file_list_modified failed to remove __simple_filelist__ in path %s\n", errno, strerror(errno), cpath);
          perror("Failed : simple_posix.c::get_file_list_reverse_modified ");
          return -1;
        }
    }
    if(mtime_alpha) { /* modified time sort */
        n = scandir(cpath, &namelist, NULL, versionsort);
    } else { /* alphanumeric sort */
        n = scandir(cpath, &namelist, NULL, alphasort);
    }

    free(cpath);
    FILE* f = fopen("__simple_filelist__", "w");
    if(!f) {
      printf("%d %s\nget_file_list_modified failed to open __simple_filelist__ in %s\n", errno, strerror(errno), cpath);
      perror("Failed : simple_posix.c::get_file_list_modified ");
      return -1;
    }
    if(n < 0)
        perror("scandir failed to find files ; get_file_list_modified");
    else {
        *count = n; // count is a temp here since not all elements in scandir are directories
        int i = n; // for reverse cases
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
                for(i = 0; i < *count; i++) {
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
                for(i = 0; i < *count; i++) {
                    if((namelist[i]->d_type != DT_DIR) && (strstr(namelist[i]->d_name, cext) != NULL)) {
                        fprintf(f, "%s\n", namelist[i]->d_name);
                        n++;
                    }
                    free(namelist[i]);
                }
            }
            free(cext);
        }
        *count = n;
        free(namelist);
    }
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
int list_dirs(const char * path, int* count, size_t ivf_path)
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


int removedir(char *buf, size_t ivf_buf)
{
    char *buffer;
    extern int errno;
    int status = -1;

    if(!(buffer = F90toCstring(buf, strlen(buf)))) {
        perror("Failed : appendtostring (buf) in posixwrapper:fortranrmdir");
        exit(1);
    }
    errno = 0;
    status = rmdir(buffer);
    if(status == -1) {
        printf("%d %s\n", errno, strerror(errno));
        perror("Failed : rmdir in posixwrapper:fortranrmdir");
        exit(1);
    }
    free(buffer);
    return status;
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
void quick_glob(glob_t*globlist, int start, int count)
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
int glob_file_list(const char *match,  int*count, int* sort_by_time, size_t ivf_match) //void** results, int *count)
{
    /* Check input char pointer */
    if(match == NULL || strlen(match) <= 0) {
        fprintf(stderr, "simple_posix.c::glob_file_list match string null .. Failing\n");
        return -1;;
    }
    char *cmatch = F90toCstring(match, strlen(match));
    if(cmatch == NULL) {
        printf("%d %s\n glob_file_list failed to convert string (unprotected) %s\n", errno, strerror(errno), match);
        perror("Failed : simple_posix.c::glob_file_list 1");
        return -1;
    }
    glob_t        globlist;
    FILE* f;
    int err;
    err = 0;
    printf("DEBUG: In glob_file_list size cmatch:%zd size match:%zd\n", strlen(cmatch), strlen(match));
    printf("DEBUG: In glob_file_list cmatch:%s\n", cmatch);
    printf("DEBUG: In glob_file_list flag:%d \n", *sort_by_time);

    int glob_val = glob(cmatch, GLOB_PERIOD, NULL, &globlist);
    free(cmatch);
    f = fopen("__simple_filelist__", "w");
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
        *count = (int) globlist.gl_pathc; //account for __simple_filelist__
        printf(" glob_file_list size before trimming %d \n", *count);
        if(*count > 2) {
            if(*sort_by_time) quick_glob(&globlist, 2, *count);
            *count = 0;
            int i = 2;
            while(globlist.gl_pathv[i]) {
                //          fprintf(stderr, "GLOB>> %s\n", globlist.gl_pathv[i]);
                if(strstr(globlist.gl_pathv[i], "__simple_filelist__") == NULL) {

                    fprintf(f, "%s\n", globlist.gl_pathv[i]);
                    *count = *count + 1;
                }
                i++;
            }
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
        }
        fprintf(stderr, " glob_file_list size after trimming %d \n", *count);
    }
    globfree(&globlist);
    fclose(f);

    return err;
}

/**
 *  \brief glob_rm_all removes dirs/files in relative directory using GLOB
 *
 *  \param match Fortran string input glob
 *  \param count return ptr for number of files/dirs found
 *  \return return status success=0
 */
int glob_rm_all(const char *match,  int*count, size_t ivf_match)
{
    /* Check input char pointer */
    if(match == NULL || strlen(match) <= 0) {
        fprintf(stderr, "simple_posix.c::glob_rm_all match string null .. Failing\n");
        return -1;;
    }
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
    printf("DEBUG: In glob_rm_all size cmatch:%zd size match:%zd\n", strlen(cmatch), strlen(match));
    printf("DEBUG: In glob_rm_all cmatch:%s\n", cmatch);

    int glob_val = glob(cmatch, GLOB_PERIOD, NULL, &globlist);
    free(cmatch);
    f = fopen("__simple_filelist__", "w");
    if(!f) {
        printf("%d %s\nglob_rm_all failed to open __simple_filelist__\n", errno, strerror(errno));
        perror("Failed : simple_posix.c::glob_rm_all 2");
        return -1;
    }

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
        printf(" glob_rm_all size before trimming %d \n", *count);
        if(*count > 2) {
            *count = 0;
            int i = 0;
            while(globlist.gl_pathv[i]) {
                fprintf(stderr, "GLOB>> %s ... ", globlist.gl_pathv[i]);
                if(strstr(globlist.gl_pathv[i], "__simple_filelist__") == NULL) {
                    if(strcmp(globlist.gl_pathv[i], ".") != 0 || strcmp(globlist.gl_pathv[i], "..") != 0) {
                        fprintf(stderr, " for deletion\n");
                        fprintf(f, "%s\n", globlist.gl_pathv[i]);
                        *count = *count + 1;
                    } else fprintf(stderr, "dot dir found, CANNOT DELETE CWD, stepping over\n");
                } else fprintf(stderr, "temp file found, stepping over\n");
                i++;
            }
            i = 0;
            while(globlist.gl_pathv[i]) {
                if(strstr(globlist.gl_pathv[i], "__simple_filelist__") == NULL) {
                    if(strcmp(globlist.gl_pathv[i], ".") != 0 || strcmp(globlist.gl_pathv[i], "..") != 0) {
                        fprintf(stderr, "GLOB>> %s ... ", globlist.gl_pathv[i]);
                        struct stat s;
                        int staterr = stat(globlist.gl_pathv[i], &s);
                        if(S_ISDIR(s.st_mode)) {
                            fprintf(stderr, " dir ");
                        } else if(S_ISREG(s.st_mode)) {
                            fprintf(stderr, " file ");
                        } else if(S_ISLNK(s.st_mode)) {
                            fprintf(stderr, " link ");
                        } else {
                            fprintf(stderr, " unsupported ");
                        }
                        // remove is a stdlib function that uses unlink for files and rmdir for directories
                        if(remove(globlist.gl_pathv[i]) == 0) {
                          fprintf(stderr, "deleted\n");
                        }
                        else{
                          fprintf(stderr, "\n%d %s\nglob_rm_all failed to remove %s\n", errno, strerror(errno), globlist.gl_pathv[i]);
                          perror("Failed : simple_posix.c::glob_rm_all ");
                        }
                    }
                }
                i++;
            }
        }
        fprintf(stderr, " glob_rm_all size after trimming %d \n", *count);
    }
    globfree(&globlist);
    fclose(f);

    return err;
}



void free_file_list(char**ptr, int n)
{
    for(int i = 0; i < n; i++) free(ptr[i]);
    free(ptr);
}

int subprocess(const char* cmd, char* const args, size_t ivf_cmd)
{
    int pid = fork();
    if(pid == 0) {
        execvp(cmd , &args);
    }
    return pid;

}

/**
 * Copy file in chunks of MAX_CHAR_FILENAME bytes
 *
 */
int fcopy(char* file1,  char* file2, size_t ivf_file1)
{
    char *buffer1, *buffer2;
    int flag = 0;
    if(!(buffer1 = F90toCstring(file1, strlen(file1)))) {
        perror("Failed : appendtostring (buf) in posixwrapper:fortranrmdir");
        free(buffer1);
        return 1;
    }
    if(!(buffer2 = F90toCstring(file2, strlen(file2)))) {
        perror("Failed : appendtostring (buf) in posixwrapper:fortranrmdir");
        free(buffer2);
        return 1;
    }
    FILE *f1, *f2;
    f1 = fopen(buffer1, "r");
    f2 = fopen(buffer2, "w");
    free(buffer1); free(buffer2);
    if(f1 && f2) {
        char            buffer[MAX_CHAR_FILENAME];
        size_t          n;

        while((n = fread(buffer, sizeof(char), sizeof(buffer), f1)) > 0) {
            if(fwrite(buffer, sizeof(char), n, f2) != n) {
                flag = 1;
                break;
            }
        }
        fclose(f1); fclose(f2);

    } else
        flag = 1;

    if(flag == 0) return 0;
    else return -1;
}

/**
 * returns a null-terminated string containing the canonicalized absolute pathname corresponding to path.
 */
int  get_absolute_pathname(const char* in, char* out, int* outlen, size_t ivf_in)
{
    // realpath(in, resolved_path);
    //      printf("\n%s\n",resolved_path);
    //      return 0;
    char *filein = F90toCstring(in, strlen(in));
    char *resolved = canonicalize_file_name(filein); free(filein);
    if(resolved == NULL) {
        printf("%d %s\nget_absolute_path failed to canonicalize  file %s\n", errno, strerror(errno), filein);
        perror("Failed : simple_posix.c::glob_file_list ");
        return 1;
    }
    *outlen = strlen(resolved);
    f2cstr(out , resolved, *outlen);
    //strncpy(out, resolved, *outlen);
    free(resolved);
    //out = resolved;
    return 0;
}

#include <sys/sysinfo.h>
int get_sysinfo(long* HWMusage, long*totalram, long* sharedram, long* bufferram, long* totalhigh)
{
    struct sysinfo s;
    int i = sysinfo(&s);
    *totalram = s.totalram;           /* Total usable main memory size */
    //*freeram=s.freeram;               /* Available memory size */
    *sharedram = s.sharedram;       /* Amount of shared memory */
    *bufferram = s.bufferram;              /* Memory used by buffers */
    *totalhigh = s.totalhigh;              /* Total high memory size */
    *HWMusage = s.totalhigh - s.freehigh; /* high memory size used */
    return i;
}
