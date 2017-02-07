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

#ifndef __LIB_MATH__
#define __LIB_MATH__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>

typedef struct {
	float r, i;
} complex;

typedef struct {
	double r, i;
} doublecomplex;

typedef struct {
	complex x, y;
} hypercomplex;

typedef struct {
	doublecomplex x, y;
} doublehypercomplex;

/* complex */

complex Cadd(complex a, complex b);
complex Csub(complex a, complex b);
complex Cmul(complex a, complex b);
complex Conj(complex z);
complex Cdiv(complex a, complex b);
float Cabs(complex z);
complex Csqrt(complex z);
complex RCmul(float x, complex a);
complex Cexp(complex z);

/* doublecomplex */

doublecomplex DCadd(doublecomplex a, doublecomplex b);
doublecomplex DCsub(doublecomplex a, doublecomplex b);
doublecomplex DCmul(doublecomplex a, doublecomplex b);
doublecomplex DConj(doublecomplex z);
doublecomplex DCdiv(doublecomplex a, doublecomplex b);
double DCabs(doublecomplex z);
doublecomplex DCsqrt(doublecomplex z);
doublecomplex RDCmul(double x, doublecomplex a);
doublecomplex DCexp(doublecomplex z);

/* hypercomplex */

hypercomplex HCadd(hypercomplex a, hypercomplex b);
hypercomplex HCsub(hypercomplex a, hypercomplex b);
hypercomplex HCmul(hypercomplex a, hypercomplex b);
hypercomplex HConj(hypercomplex z);
hypercomplex HCtp(hypercomplex z);

#endif
