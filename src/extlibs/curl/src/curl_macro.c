/* curl_macro.c */
#include <curl/curl.h>

#ifdef __cplusplus
extern "C" {
#endif

int curl_easy_setopt_c_char(CURL *, CURLoption, char *);
int curl_easy_setopt_c_funptr(CURL *, CURLoption, void *);
int curl_easy_setopt_c_int(CURL *, CURLoption, int);
int curl_easy_setopt_c_long(CURL *, CURLoption, long);
int curl_easy_setopt_c_ptr(CURL *, CURLoption, void *);
int curl_easy_getinfo_c_double(CURL *, CURLINFO, double *);
int curl_easy_getinfo_c_long(CURL *, CURLINFO, long *);
int curl_easy_getinfo_c_ptr(CURL *, CURLINFO, void **);
int curl_version_now(void);

/* Non-variadic wrappers to `curl_easy_setopt()`. */
int curl_easy_setopt_c_char(CURL *curl, CURLoption option, char *value)
{
    return curl_easy_setopt(curl, option, value);
}

int curl_easy_setopt_c_funptr(CURL *curl, CURLoption option, void *value)
{
    return curl_easy_setopt(curl, option, value);
}

int curl_easy_setopt_c_int(CURL *curl, CURLoption option, int value)
{
    return curl_easy_setopt(curl, option, value);
}

int curl_easy_setopt_c_long(CURL *curl, CURLoption option, long value)
{
    return curl_easy_setopt(curl, option, value);
}

int curl_easy_setopt_c_ptr(CURL *curl, CURLoption option, void *value)
{
    return curl_easy_setopt(curl, option, value);
}

/* Non-variadic wrappers to `curl_easy_getinfo()`. */
int curl_easy_getinfo_c_double(CURL *curl, CURLINFO info, double *value)
{
    return curl_easy_getinfo(curl, info, value);
}

int curl_easy_getinfo_c_long(CURL *curl, CURLINFO info, long *value)
{
    return curl_easy_getinfo(curl, info, value);
}

int curl_easy_getinfo_c_ptr(CURL *curl, CURLINFO info, void **value)
{
    return curl_easy_getinfo(curl, info, value);
}

/* Wrapper function that provides access to the constant `CURLVERSION_NOW`. */
int curl_version_now(void)
{
    return CURLVERSION_NOW;
}

#ifdef __cplusplus
}
#endif
