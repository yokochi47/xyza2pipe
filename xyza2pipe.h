/*
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

/* xyza2pipe: Cross conversion environment of NMR spectra
   https://github.com/yokochi47/xyza2pipe forked from http://fermi.pharm.hokudai.ac.jp/olivia/
 */

#ifndef __XYZA2PIPE__
#define __XYZA2PIPE__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <ctype.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>
#include <getopt.h>
#include <dirent.h>

#include "libString.h"
#include "libMath.h"
#include "libMatrix.h"

#define MAXASSNAME		8
#define MAXFILENAME		64
#define MAXCHAR			1024
#define MAXLONGNAME		8192

#define PIPE_HEADER_0		0.0
#define PIPE_HEADER_1		-286331168.0
#define PIPE_HEADER_2		2.345

#define PIPE_HEADERSIZE		2048
#define NV_HEADERSIZE		2048

#define UCSF_MAXBLOCKSIZE	8192	/* Max Block Size: 32KB */
#define NV_MAXBLOCKSIZE		8192	/* Max Block Size: 32KB */
#define XEASY_MAXBLOCKSIZE	8192	/* Max Block Size: 16KB */

#define NULLPPM			-100.0

/* APODIZATION */

typedef enum { NO_APODIZATION /* NA */ ,
	SINE_BELL /* SP */ ,
	EXPONENTIAL /* EM */ ,
	LORENTZ_GAUSS /* GM */ ,
	TRAPEZOID /* TM */ ,
	ZERO_RANGE /* ZE */ ,
	TRIANGLE /* TRI */ ,
	EXPONENTIAL_GAUSS		/* GMB */
} enum_apod_code;

/* COMBINE OPERATOR */

typedef enum { COMBINE_ADD, COMBINE_SUB, COMBINE_MUL, COMBINE_DIV } enum_combine_opr;

extern char header[PIPE_HEADERSIZE];

extern char axislabel[4][MAXASSNAME + 1];
extern char axisname[4][MAXASSNAME + 1];

extern char dimension, _dimension, byteswap;

extern char swapdata, swappar, usrlabel, usrshift;
extern char leftcar, extleft, adjcar, adjh2o, relyof;
extern char usrphase[4];

extern short headersize;
extern short blocksize[4];

extern int datasize[4], datasize_orig[4];
extern int unitsize[4];

extern float obsfreq[4];
extern float spcenter[4];
extern float origfreq[4];
extern float spwidth[4];
extern float usrcenter[4];
extern float phase[4][2];

extern const float PIPE_HEADER[3];

void checklabel(char string[]);

void cnvhdr(const char axis_option, const char fb);

int checkpipe();
void make_dir(char filename[]);
void set_clean_string(char clean_string[]);

int openpipe2d(float **mat2d);

int checkxyza(char filename[]);

int pushxyza2d(char spectra2d[], const char axis_option);
int pushxyza3d(char spectra3d[], const char axis_option);
int pushxyza4d(char spectra4d[], const char axis_option);

int openxyza2d(char spectra2d[], float **mat2d);
int openxyza3d(char spectra3d[], const int z, float **mat2d);
int openxyza4d(char spectra4d[], const int z, const int a, float **mat2d);

int pullxyza2d(char spectra2d[], const char axis_option);
int pullxyza3d(char spectra3d[], const char axis_option);
int pullxyza4d(char spectra4d[], const char axis_option);

int checkucsf(char filename[]);

int pushucsf2d(char spectra2d[], const char axis_option);
int pushucsf3d(char spectra3d[], const char axis_option);
int pushucsf4d(char spectra4d[], const char axis_option);

int openucsf2d(char spectra2d[], float **mat2d);
int openucsf3d(char spectra3d[], float ***mat3d);
int openucsf4d(char spectra4d[], float ****mat4d);

int pullucsf2d(char spectra2d[]);
int pullucsf3d(char spectra3d[]);
int pullucsf4d(char spectra4d[]);

int checknv(char filename[]);

int pushnv2d(char spectra2d[], const char axis_option);
int pushnv3d(char spectra3d[], const char axis_option);
int pushnv4d(char spectra4d[], const char axis_option);

int opennv2d(char spectra2d[], float **mat2d);
int opennv3d(char spectra3d[], float ***mat3d);
int opennv4d(char spectra4d[], float ****mat4d);

int pullnv2d(char spectra2d[]);
int pullnv3d(char spectra3d[]);
int pullnv4d(char spectra4d[]);

int checkxeasy(char filename[]);

int pushxeasy2d(char spectra2d[], const char axis_option);
int pushxeasy3d(char spectra3d[], const char axis_option);
int pushxeasy4d(char spectra4d[], const char axis_option);

int openxeasy2d(char spectra2d[], float **mat2d);
int openxeasy3d(char spectra3d[], float ***mat3d);
int openxeasy4d(char spectra4d[], float ****mat4d);

int pullxeasy2d(char spectra2d[]);
int pullxeasy3d(char spectra3d[]);
int pullxeasy4d(char spectra4d[]);

int checkazara(char filename[]);

int pushazara2d(char spectra2d[], const char axis_option);
int pushazara3d(char spectra3d[], const char axis_option);
int pushazara4d(char spectra4d[], const char axis_option);

int openazara2d(char spectra2d[], float **mat2d);
int openazara3d(char spectra3d[], float ***mat3d);
int openazara4d(char spectra4d[], float ****mat4d);

int pullazara2d(char spectra2d[]);
int pullazara3d(char spectra3d[]);
int pullazara4d(char spectra4d[]);

int checkvnmr(char filename[], char pardir[], char monofile[]);

int pushvnmr2d(char monofile[], char pardir[], const char axis_option);
int pushvnmr3d(char monofile[], char pardir[], const char axis_option);
int pushvnmr4d(char monofile[], char pardir[], const char axis_option);

int openvnmr2d(char monofile[], char pardir[], float **mat2d);
int openvnmr3d(char monofile[], char pardir[], float ***mat3d);
int openvnmr4d(char monofile[], char pardir[], float ****mat4d);

int checkxwnmr(char filename[]);

int pushxwnmr2d(char spectra2d[], const char axis_option);
int pushxwnmr3d(char spectra3d[], const char axis_option);
int pushxwnmr4d(char spectra4d[], const char axis_option);

int openxwnmr2d(char spectra2d[], float **mat2d);
int openxwnmr3d(char spectra3d[], float ***mat3d);
int openxwnmr4d(char spectra4d[], float ****mat4d);

int pullproj2d(char spectra2d[], const int abs_mode);
int pullproj3d(char spectra3d[], const int abs_mode);
int pullproj4d(char spectra4d[], const int abs_mode);

int pushadd2d(char spectra2d1[], char spectra2d2[], const float c1, const float c2, enum_combine_opr opr_code);
int pushadd3d(char spectra3d1[], char spectra3d2[], const float c1, const float c2, enum_combine_opr opr_code);
int pushadd4d(char spectra4d1[], char spectra4d2[], const float c1, const float c2, enum_combine_opr opr_code);

int pushadducsf2d(char spectra2d1[], char spectra2d2[], const float c1, const float c2, enum_combine_opr opr_code);
int pushadducsf3d(char spectra3d1[], char spectra3d2[], const float c1, const float c2, enum_combine_opr opr_code);
int pushadducsf4d(char spectra4d1[], char spectra4d2[], const float c1, const float c2, enum_combine_opr opr_code);

int pushaddnv2d(char spectra2d1[], char spectra2d2[], const float c1, const float c2, enum_combine_opr opr_code);
int pushaddnv3d(char spectra3d1[], char spectra3d2[], const float c1, const float c2, enum_combine_opr opr_code);
int pushaddnv4d(char spectra4d1[], char spectra4d2[], const float c1, const float c2, enum_combine_opr opr_code);

int pushaddxeasy2d(char spectra2d1[], char spectra2d2[], const float c1, const float c2, enum_combine_opr opr_code);
int pushaddxeasy3d(char spectra3d1[], char spectra3d2[], const float c1, const float c2, enum_combine_opr opr_code);
int pushaddxeasy4d(char spectra4d1[], char spectra4d2[], const float c1, const float c2, enum_combine_opr opr_code);

int pushaddazara2d(char spectra2d1[], char spectra2d2[], const float c1, const float c2, enum_combine_opr opr_code);
int pushaddazara3d(char spectra3d1[], char spectra3d2[], const float c1, const float c2, enum_combine_opr opr_code);
int pushaddazara4d(char spectra4d1[], char spectra4d2[], const float c1, const float c2, enum_combine_opr opr_code);

int pushaddvnmr2d(char monofile1[], char monofile2[], char pardir1[], char pardir2[], const float c1, const float c2,
		enum_combine_opr opr_code);
int pushaddvnmr3d(char monofile1[], char monofile2[], char pardir1[], char pardir2[], const float c1, const float c2,
		enum_combine_opr opr_code);
int pushaddvnmr4d(char monofile1[], char monofile2[], char pardir1[], char pardir2[], const float c1, const float c2,
		enum_combine_opr opr_code);

int pushaddxwnmr2d(char spectra2d1[], char spectra2d2[], const float c1, const float c2, enum_combine_opr opr_code);
int pushaddxwnmr3d(char spectra3d1[], char spectra3d2[], const float c1, const float c2, enum_combine_opr opr_code);
int pushaddxwnmr4d(char spectra4d1[], char spectra4d2[], const float c1, const float c2, enum_combine_opr opr_code);

int checkdefl(char filename[]);

#endif
