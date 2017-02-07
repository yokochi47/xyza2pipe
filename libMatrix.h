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

#ifndef __LIB_MATRIX__
#define __LIB_MATRIX__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <assert.h>

float *fmalloc1d(const size_t size1);
float **fmalloc2d(const size_t size1, const size_t size2);
float ***fmalloc3d(const size_t size1, const size_t size2, const size_t size3);
float ****fmalloc4d(const size_t size1, const size_t size2, const size_t size3, const size_t size4);

void free_fmatrix1d(float *p);
void free_fmatrix2d(float **p);
void free_fmatrix3d(float ***p);
void free_fmatrix4d(float ****p);

#endif
