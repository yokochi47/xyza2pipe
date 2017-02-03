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

#include "xyza2pipe.h"

char header[PIPE_HEADERSIZE];
char axislabel[4][MAXASSNAME + 1] = { {""} };

char axisname[4][MAXASSNAME + 1];

char dimension, _dimension, byteswap;

char swapdata = 0, swappar = 0, usrlabel = 0, usrshift = 0;
char leftcar = 0, extleft = 0, adjcar = 0, adjh2o = 0, relyof = 0;
char usrphase[4] = { 0 };

short headersize;
short blocksize[4];

int datasize[4], datasize_orig[4];
int unitsize[4];

float obsfreq[4];
float spcenter[4];
float origfreq[4];
float spwidth[4];
float usrcenter[4] = { NULLPPM };

float phase[4][2];

const float PIPE_HEADER[3] = { PIPE_HEADER_0, PIPE_HEADER_1, PIPE_HEADER_2 };
