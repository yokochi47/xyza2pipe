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

#include "libMatrix.h"

float *fmalloc1d(const size_t size1)
{
	float *p;

	assert(p = (float *) malloc(size1 * sizeof(float)));

	return p;
}

float **fmalloc2d(const size_t size1, const size_t size2)
{
	int i;
	float **p;

	assert(p = (float **) malloc(size1 * sizeof(float *)));

	assert(p[0] = (float *) malloc(size1 * size2 * sizeof(float)));

	for (i = 1; i < size1; i++)
		p[i] = p[i - 1] + size2;

	return p;
}

float ***fmalloc3d(const size_t size1, const size_t size2, const size_t size3)
{
	int i, j;
	float ***p;

	assert(p = (float ***) malloc(size1 * sizeof(float **)));

	assert(p[0] = (float **) malloc(size1 * size2 * sizeof(float *)));

	assert(p[0][0] = (float *) malloc(size1 * size2 * size3 * sizeof(float)));

	for (j = 1; j < size2; j++)
		p[0][j] = p[0][j - 1] + size3;

	for (i = 1; i < size1; i++) {
		p[i] = p[i - 1] + size2;
		p[i][0] = p[i - 1][0] + size2 * size3;
		for (j = 1; j < size2; j++)
			p[i][j] = p[i][j - 1] + size3;
	}

	return p;
}

float ****fmalloc4d(const size_t size1, const size_t size2, const size_t size3, const size_t size4)
{
	int i, j, k;
	float ****p;

	assert(p = (float ****) malloc(size1 * sizeof(float ***)));

	assert(p[0] = (float ***) malloc(size1 * size2 * sizeof(float **)));

	assert(p[0][0] = (float **) malloc(size1 * size2 * size3 * sizeof(float *)));

	assert(p[0][0][0] = (float *) malloc(size1 * size2 * size3 * size4 * sizeof(float)));

	for (k = 1; k < size3; k++)
		p[0][0][k] = p[0][0][k - 1] + size4;

	for (j = 1; j < size2; j++) {
		p[0][j] = p[0][j - 1] + size3;
		p[0][j][0] = p[0][j - 1][0] + size3 * size4;
		for (k = 1; k < size3; k++)
			p[0][j][k] = p[0][j][k - 1] + size4;
	}

	for (i = 1; i < size1; i++) {
		p[i] = p[i - 1] + size2;
		p[i][0] = p[i - 1][0] + size2 * size3;
		p[i][0][0] = p[i - 1][0][0] + size2 * size3 * size4;

		for (k = 1; k < size3; k++)
			p[i][0][k] = p[i][0][k - 1] + size4;

		for (j = 1; j < size2; j++) {
			p[i][j] = p[i][j - 1] + size3;
			p[i][j][0] = p[i][j - 1][0] + size3 * size4;
		}
	}

	for (i = 1; i < size1; i++) {
		for (k = 1; k < size3; k++)
			p[i][0][k] = p[i - 1][0][k - 1] + size2 * size3 * size4 + size4;

		for (j = 1; j < size2; j++) {
			p[i][j] = p[i - 1][j - 1] + size2 * size3 + size3;
			p[i][j][0] = p[i - 1][j - 1][0] + size2 * size3 * size4 + size3 * size4;
			for (k = 1; k < size3; k++)
				p[i][j][k] = p[i - 1][j][k - 1] + size2 * size3 * size4 + size4;
		}
	}

	return p;
}

void free_fmatrix1d(float *p)
{
	if (p == NULL)
		return;

	free((char *) (p));

	p = NULL;
}

void free_fmatrix2d(float **p)
{
	if (p == NULL)
		return;

	free((char *) (p[0]));
	free((char *) (p));

	p = NULL;
}

void free_fmatrix3d(float ***p)
{
	if (p == NULL)
		return;

	free((char *) (p[0][0]));
	free((char *) (p[0]));
	free((char *) (p));

	p = NULL;
}

void free_fmatrix4d(float ****p)
{
	if (p == NULL)
		return;

	free((char *) (p[0][0][0]));
	free((char *) (p[0][0]));
	free((char *) (p[0]));
	free((char *) (p));

	p = NULL;
}
