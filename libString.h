/*
   xyza2pipe - a cross conversion environment of NMR spectra
   Copyright 2017 Masashi Yokochi

   https://github.com/yokochi47/xyza2pipe
    forked from http://fermi.pharm.hokudai.ac.jp/olivia/

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
 */

#ifndef __LIB_STRING__
#define __LIB_STRING__

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <ctype.h>
#include <time.h>
#include <assert.h>

#define MAXVARS		256
#define MAXCHAR		1024

void path2uscore(char *path);

int column_of(const char *str);
int row_of(const char *str);
int quote_of(char *str);

int is_integer(const char *str);

int strsearch(const char *haystack, const char *needle);
int strcasesearch(const char *haystack, const char *needle);

int strmatchcount(const char *haystack, const char *needle);
int strcasematchcount(const char *haystack, const char *needle);

char *struppercpy(char *dst, const char *src);
char *strlowercpy(char *dst, const char *src);
char *strclasscpy(char *dst, const char *src);
char *strunquotecpy(char *dst, const char *src);
char *strseqoffcpy(char *dst, const char *src);

char *dirnamecpy(char *dir, const char *path);
char *basenamecpy(char *baes, const char *path);

char *strspacefillcpy(char *dst, const char *src, size_t len);
char *strspacefill(char *src, size_t len);

char *strreplacecpy(char *dst, const char *src, const char *key, const char *trl);
char *strreplace(char *src, const char *key, const char *trl);

char *strsplit(const char *src, const char *key);

char *strcenteringcpy(char *dst, const char *src, size_t len);
char *strquarteringcpy(char *dst, const char *src, size_t len);

void swrite2bin(FILE * fp, int16_t a);
void iwrite2bin(FILE * fp, int32_t a);
void lwrite2bin(FILE * fp, long a);
void fwrite2bin(FILE * fp, float a);

void spwrite2bin(FILE * fp, int16_t * a, size_t n);
void ipwrite2bin(FILE * fp, int32_t * a, size_t n);
void lpwrite2bin(FILE * fp, long *a, size_t n);
void fpwrite2bin(FILE * fp, float *a, size_t n);

void swrite2mem(void *mem, int16_t a);
void iwrite2mem(void *mem, int32_t a);
void lwrite2mem(void *mem, long a);
void fwrite2mem(void *mem, float a);

void swrite2bin_swap(FILE * fp, int16_t a, int swap);
void iwrite2bin_swap(FILE * fp, int32_t a, int swap);
void lwrite2bin_swap(FILE * fp, long a, int swap);
void fwrite2bin_swap(FILE * fp, float a, int swap);

void spwrite2bin_swap(FILE * fp, int16_t * a, size_t n, int swap);
void ipwrite2bin_swap(FILE * fp, int32_t * a, size_t n, int swap);
void lpwrite2bin_swap(FILE * fp, long *a, size_t n, int swap);
void fpwrite2bin_swap(FILE * fp, float *a, size_t n, int swap);

void swrite2mem_swap(void *mem, int16_t a, int swap);
void iwrite2mem_swap(void *mem, int32_t a, int swap);
void lwrite2mem_swap(void *mem, long a, int swap);
void fwrite2mem_swap(void *mem, float a, int swap);

int is_big_endian();
int is_little_endian();

int is_big_endian_float(char *mem, float real);
int is_little_endian_float(char *mem, float real);

void swapbyte(unsigned int unit_size, unsigned int mem_size, char *mem);

char *parse2str(char *p_char, const char *str, unsigned int start, unsigned int stop);
int *parse2int(int *p_int, const char *str, unsigned int start, unsigned int stop);
float *parse2float(float *p_float, const char *str, unsigned int start, unsigned int stop);

int line2arg(const char *str, unsigned char fs_code, char *argv[]);
int line4arg(const char *str, unsigned char fs_code, char *argv[], unsigned char start_code, unsigned char stop_code);

/* NOTE THAT FOLLOWING FUNCTIONS ARE THREAD-UNSAFE FUNCTIONS */

void alloc_arg(unsigned int thrd_id, int size, char *argv[]);
void free_arg(unsigned int thrd_id);

char *_dirname_(const char *path);
char *_basename_(const char *path);

#endif
