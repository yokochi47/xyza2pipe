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

#ifndef __VENDOR_FID__
#define __VENDOR_FID__

#define M_CONST		-0.009552
#define B_CONST		5.011718

extern const float gH1;
extern const float gC13;
extern const float gN15;
extern const float gP31;

/* VARIAN BINARY FORMAT */

#define S_DATA		0x1
#define S_SPEC		0x2
#define S_32		0x4
#define S_FLOAT		0x8
#define S_COMPLEX	0x10
#define S_HYPERCOMPLEX	0x20

#define S_ACQPAR	0x80
#define S_SECND		0x100
#define S_TRANSF	0x200
#define S_NP		0x800
#define S_NF		0x1000
#define S_NI		0x2000
#define S_NI2		0x4000
#define S_NI3		0x8000

#define MORE_BLOCKS	0x80
#define NP_CMPLX	0x100
#define NF_CMPLX	0x200
#define NI_CMPLX	0x400
#define NI2_CMPLX	0x800
#define NI3_CMPLX	0x1000

#define NP_PHMODE	0x1
#define NP_AVMODE	0x2
#define NP_PWRMODE	0x4
#define NF_PHMODE	0x10
#define NF_AVMODE	0x20
#define NF_PWRMODE	0x40
#define NI_PHMODE	0x100
#define NI_AVMODE	0x200
#define NI_PWRMODE	0x400
#define NI2_PHMODE	0x1000
#define NI2_AVMODE	0x2000
#define NI2_PWRMODE	0x4000
#define NI3_PHMODE	0x8
#define NI3_AVMODE	0x80
#define NI3_PWRMODE	0x800

#define U_HYPERCOMPLEX	0x2

/* VARIAN BINARY HEADER STRUCTURE */

typedef struct {
	int32_t nblocks, ntraces, np, ebytes, tbytes, bbytes;
	short vers_id, status;
	int32_t nbheaders;
} datafilehead;

typedef struct {
	short scale, status, index, mode;
	int32_t ctcount;
	float lpval, rpval, lvl, tlt;
} datablockhead;

typedef struct {
	short s_spare1, status, s_spare2, s_spare_3;
	int32_t l_spare1;
	float lpval1, rpval1, f_spare1, f_spare2;
} hypercmplxbhead;

int get_varian_parsize(unsigned int dim, char pardir[], char parname[]);
int get_varian_parameter(unsigned int dim, char pardir[], char parname[], unsigned int col, char string[]);
int get_varian_dimension_from_file(unsigned int dim, char pardir[]);
int guess_varian_dimension_from_file(char pardir[]);

int get_bruker_proc_parsize(unsigned int dim, char procdir[], unsigned int axis, char parname[]);
int get_bruker_proc_parameter(unsigned int dim, char procdir[], unsigned int axis, char parname[], unsigned int col,
		char string[]);
int get_bruker_acq_parsize(unsigned int dim, char acqdir[], unsigned int axis, char parname[]);
int get_bruker_acq_parameter(unsigned int dim, char acqdir[], unsigned int axis, char parname[], unsigned int col,
		char string[]);
int get_bruker_dimension_from_acq_file(unsigned int dim, char acqdir[]);
int guess_bruker_dimension_from_acq_file(char acqdir[]);
int guess_bruker_dimension_from_proc_file(char procdir[]);

#endif
