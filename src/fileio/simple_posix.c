/*****************************************************************************/
/*                       POSIX system calls for SIMPLE                       */
/* Author:  Michael Eager                                                    */
/*                                                                           */
/*****************************************************************************/
#define _SVID_SOURCE
#define _GNU_SOURCE             /* See feature_test_macros(7) */
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <sys/wait.h>
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
#define NORMAL_COLOR  "\x1B[0m"
#define GREEN  "\x1B[32m"
#define BLUE  "\x1B[34m"


char *F90toCstring(char *str, int len)
{
  char *res; /* C arrays are from 0:len-1 */
  if (len > MAX_CHAR_FILENAME){
    printf( "simple_posix.c possible strlen error in converting f90 string to c-string\n");
  }
  if((res = (char*)malloc(len + 1))) {
    strncpy(res, str, len);
    res[len] = '\0';
  }
  return res;
}

#define MAKEDIR makedir_
int makedir(const char *pathname)
{
  mode_t mask = umask(0);
  umask(mask);
  return mkdir(pathname, mask);
}

/* let us make a recursive function to print the content of a given folder */
#define SHOW_DIR_CONTENT_RECURSIVE show_dir_content_recursive_
void show_dir_content_recursive(const char * path)
{
  DIR * d = opendir(path); // open the path
  if(d == NULL) return; // if was not able return
  struct dirent * dir; // for the directory entries
  while((dir = readdir(d)) != NULL) { // if we were able to read somehting from the directory
    if(dir-> d_type != DT_DIR) // if the type is not directory just print it with blue
      printf("%s%s\n", BLUE, dir->d_name);
    else if(dir -> d_type == DT_DIR && strcmp(dir->d_name, ".") != 0 && strcmp(dir->d_name, "..") != 0) { // if it is a directory
      printf("%s%s\n", GREEN, dir->d_name); // print its name in green
      char d_path[255]; // here I am using sprintf which is safer than strcat
      sprintf(d_path, "%s/%s", path, dir->d_name);
      show_dir_content_recursive(d_path); // recall with the new path
    }
  }
  closedir(d); // finally close the directory
}


#define GET_FILE_LIST get_file_list_
int get_file_list(const char * path,  const char * ext,int *count)
{
  char *cpath = F90toCstring(path, strlen(path));
  printf(" get_file_list %s\n", path);
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
    if(strcmp(ext,"")){
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
          if (strstr(elem->d_name, ext)!=NULL){
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


/* print files in current directory in reverse order */
#define GET_FILE_LIST_REVERSE_MODIFIED get_file_list_reverse_modified_
int get_file_list_reverse_modified(const char * path, const char* ext, int* count)
{
  /*equivalent to 'ls -tr path' */
  struct dirent **namelist;
  int n;

  *count=0;
  FILE* f = fopen("__simple_filelist__", "w");
  if(!f){
    printf("%d %s\nget_file_list_reverse_modified failed to open __simple_filelist__ in %s\n", errno, strerror(errno), path);
    perror("Failed : simple_posix.c::get_file_list ");
    return -1;
  }
  char *cpath = F90toCstring(path, strlen(path));
  n = scandir(cpath, &namelist, NULL, versionsort); free(cpath);
  if (n < 0)
    perror("scandir");
  else {
    *count = n;
    char *cext = F90toCstring(ext, 3);
    if(strcmp(ext,"")){
      while (n--) {
        if(namelist[n]->d_type != DT_DIR) fprintf(f, "%s\n", namelist[n]->d_name);
      }
      free(namelist[n]);
    } else {
      while (n--) {
        if( (namelist[n]->d_type != DT_DIR) && (strstr(namelist[n]->d_name, ext)!=NULL) )
          fprintf(f, "%s\n", namelist[n]->d_name);
      }
      free(namelist[n]);
    }
    *count = n;
    free(cext);
    free(namelist);
  }
  fclose(f);
  return 0;
}


#define LIST_DIRS list_dirs_
int list_dirs(const char * path, int* count)
{
  char *cpath = F90toCstring(path, strlen(path));
  //  printf(" simple_posix.c:: list_dirs %s\n", cpath);
  extern int errno;
  /*option 2*/
  //  char ** results;
  DIR           *d;
  int fcount = 0;
  d = opendir(cpath);
  free(cpath);
  if(d) {

    struct dirent *dir;
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
    if(f){
    fcount = 0;
    while((dir = readdir(d)) != NULL) {
      if(dir-> d_type == DT_DIR && strcmp(dir->d_name, ".") != 0 && strcmp(dir->d_name, "..") != 0) {
          /*option 1*/
          // results[fcount] = malloc ( sizeof(char) * (dir->d_reclen+1));
          // strncpy(results[fcount], dir->d_name, (dir->d_reclen+1) );
          /*option 2*/
        fprintf(f,"%s\n", dir->d_name);
          //              printf(" list_dirs elem %d %s\n", fcount, dir->d_name);
        fcount++;
      }
    }
    fclose(f);
    }else{
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


void c_rmdir_(char *buf, int *status, int buflength)
{
  char *buffer;
  extern int errno;

  *status = -1;

  if(!(buffer = F90toCstring(buf, buflength))) {
    perror("Failed : appendtostring (buf) in posixwrapper:fortranrmdir");
    exit(1);
  }
  errno = 0;
  *status = rmdir(buffer);

  if(*status == -1) {
    printf("%d %s\n", errno, strerror(errno));
    perror("Failed : rmdir in posixwrapper:fortranrmdir");
    exit(1);
  }

  *status = 0;
  free(buffer);
}



int qs_struct(struct dirent **namelist, int left, int right) {

  register int i, j;

  struct dirent *temp;
  struct stat ibuf;
  struct stat jbuf;
  struct stat xbuf;

  i = left; j = right;

  stat(namelist[i]->d_name, &ibuf);
  stat(namelist[j]->d_name, &jbuf);
  stat(namelist[(left+right)/2]->d_name, &xbuf);


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
  stat(globlist->gl_pathv[(left+right)/2], &xbuf);


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
      temp = malloc(strlen(globlist->gl_pathv[i]) +1);
      strcpy(globlist->gl_pathv[i],temp                  );
      globlist->gl_pathv[i] = realloc(globlist->gl_pathv[i], strlen(globlist->gl_pathv[j]) +1 );
      strcpy(globlist->gl_pathv[j],globlist->gl_pathv[i] );
      globlist->gl_pathv[j] = realloc( globlist->gl_pathv[j], strlen(temp) + 1 );
      strcpy(temp                 ,globlist->gl_pathv[j] );
      i++; j--;
      free(temp);
    }

  } while(i <= j);

  if(left < j) qs_glob(globlist, left, j);
  if(i < right) qs_glob(globlist, i, right);
  return 0;
}


/* A Quicksort for structures of type glob. */
void quick_glob(glob_t*globlist, int start, int count) {
  qs_glob(globlist,start,count);
}


#define GLOB_FILE_LIST glob_file_list_
int glob_file_list(const char *match, int*len_match, int*count, int* sort_by_time) //void** results, int *count)
{
  /* Check input char pointer */
  if(match == NULL || *len_match <= 0 ) {
    fprintf(stderr, "simple_posix.c::glob_file_list match string null .. Failing\n");
    exit(1);
  }

  glob_t        globlist;
  FILE* f;
  /* char * base; */
  /* if(rindex(match, '/')==NULL){ */
  /*   base=""; */
  /* } */

  //  results=NULL;
  f = fopen("__simple_filelist__", "w");
  if(!f){
    printf("%d %s\nglob_file_list failed to open __simple_filelist__\n", errno, strerror(errno));
    perror("Failed : simple_posix.c::glob_file_list ");
    return -1;
  }
  int err;
  err = 0;
  if(!f) {
    fprintf(stderr, "simple_posix.c::glob_file_list cannot open tmp file \n");
    err = -3;
  } else {
    char *cmatch = F90toCstring(match, strlen(match));
    int glob_val= glob(match, GLOB_PERIOD, NULL, &globlist);free(cmatch);
    if(glob_val == GLOB_NOSPACE){
      fprintf(stderr, "simple_posix.c::glob_file_list NOSPACE no memory left \n");
      err = -2;
    }else if(glob_val == GLOB_NOMATCH) {
      fprintf(stderr, "simple_posix.c::glob_file_list NOMATCH  glob:%s\n", match);
      err = 0; // Do nothing
    } else if( glob_val== GLOB_ABORTED) {
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
            /* if( strcmp(globlist.gl_pathv[i],".")!=0 && strcmp(globlist.gl_pathv[i],"..")!=0 ) */
            /* { */
            /*   *count --; */
            /* }else{ */
            fprintf(f, "%s\n", globlist.gl_pathv[i]);
            *count = *count + 1;
            //  }
          }
          i++;
        }
        /* (*results)[*count] = malloc( sizeof(char)*(*count)); */
        /* for(i=0;i<*count;i++) { */
        /*   int lenstr = strlen(globlist.gl_pathv[i])+1; */
        /*   results[i] = malloc ( sizeof(char) * lenstr); */
        /*   strncpy( results[i], globlist.gl_pathv[i], lenstr-1); */
        /* } */
        /* printf(" glob_file_list allocated \n"); */
        /* //#ifdef _DEBUG */
        /* i=0; */
        /* while (globlist.gl_pathv[i]) */
        /* { */
        /*   printf("%s\n", globlist.gl_pathv[i]); */
        /*   i++; */
        /* } */
        /* printf(" glob_file_list freeing globlist struct \n"); */
        // #endif
      }
      fprintf(stderr," glob_file_list size after trimming %d \n", *count);
    }
    globfree(&globlist);
    fclose(f);
  }
  return err;
}

#define FREE_FILE_LIST free_file_list_
void free_file_list(char**ptr, int n)
{
  for(int i = 0; i < n; i++) free(ptr[i]);
  free(ptr);
}


#define SUBPROCESS subprocess_
int subprocess(const char* cmd, char* const args)
{
  int pid = fork();

  if(pid == 0) {
    execvp(cmd , &args);
  }
  return pid;

}

#define FCOPY fcopy_
int fcopy(char* file1,  char* file2)
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



#define GET_ABSOLUTE_PATHNAME get_absolute_pathname_
int  get_absolute_pathname(const char* in, char out[], int* outlen)
{
  // realpath(in, resolved_path);
  //      printf("\n%s\n",resolved_path);
  //      return 0;
  char *filein = F90toCstring(in, strlen(in));
  char *resolved = canonicalize_file_name(filein);free(filein);
  if(resolved ==NULL){
    printf("%d %s\nget_absolute_path failed to canonicalize  file %s\n", errno, strerror(errno), filein);
    perror("Failed : simple_posix.c::glob_file_list ");
    return 1;
  }
  *outlen = strlen(resolved);
  strncpy(out, resolved, *outlen);
  free(resolved);
  //out = resolved;
  return 0;
}
