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

#include "xyza2pipe.h"

int pushxyza2d(char spectra2d[], const char axis_option)
{
	int i, j;
	float **matrix2d;

	/* HEADER */

	fwrite(header, sizeof(char), PIPE_HEADERSIZE, stdout);

	fflush(stdout);

	matrix2d = fmalloc2d(datasize[1], datasize[0]);

	switch (axis_option) {
	case 'x':
		if (openxyza2d(spectra2d, matrix2d) != 0)
			goto escape;

		fpwrite2bin(stdout, &(matrix2d[0][0]), datasize[1] * datasize[0]);

		fflush(stdout);

		break;

	case 'y':
		if (openxyza2d(spectra2d, matrix2d) != 0)
			goto escape;

		for (i = 0; i < datasize[0]; i++) {
			for (j = 0; j < datasize[1]; j++) {
				fpwrite2bin(stdout, &(matrix2d[j][i]), 1);
			}
		}

		fflush(stdout);

		break;
	}

	free_fmatrix2d(matrix2d);

	return 0;

	escape:free_fmatrix2d(matrix2d);

	return 1;
}

int pushxyza3d(char spectra3d[], const char axis_option)
{
	int i, j, k;
	float **matrix2d, ***matrix3d_;

	/* HEADER */

	fwrite(header, sizeof(char), PIPE_HEADERSIZE, stdout);

	fflush(stdout);

	matrix2d = fmalloc2d(datasize[1], datasize[0]);

	switch (axis_option) {
	case 'x':
		for (k = 0; k < datasize[2]; k++) {

			if (openxyza3d(spectra3d, k, matrix2d) != 0)
				goto escape;

			fpwrite2bin(stdout, &(matrix2d[0][0]), datasize[1] * datasize[0]);

			fflush(stdout);
		}

		break;

	case 'y':
		for (k = 0; k < datasize[2]; k++) {

			if (openxyza3d(spectra3d, k, matrix2d) != 0)
				goto escape;

			for (i = 0; i < datasize[0]; i++) {
				for (j = 0; j < datasize[1]; j++) {
					fpwrite2bin(stdout, &(matrix2d[j][i]), 1);
				}
			}

			fflush(stdout);
		}

		break;

	case 'z':
		matrix3d_ = fmalloc3d(datasize[1], datasize[0], datasize[2]);

		for (k = 0; k < datasize[2]; k++) {

			if (openxyza3d(spectra3d, k, matrix2d) != 0)
				goto escape2;

			for (j = 0; j < datasize[1]; j++) {
				for (i = 0; i < datasize[0]; i++) {
					matrix3d_[j][i][k] = matrix2d[j][i];
				}
			}
		}

		for (j = 0; j < datasize[1]; j++) {
			fpwrite2bin(stdout, &(matrix3d_[j][0][0]), datasize[0] * datasize[2]);

			fflush(stdout);
		}

		free_fmatrix3d(matrix3d_);

		break;
	}

	free_fmatrix2d(matrix2d);

	return 0;

	escape2:free_fmatrix3d(matrix3d_);
	escape:free_fmatrix2d(matrix2d);

	return 1;
}

int pushxyza4d(char spectra4d[], const char axis_option)
{
	int i, j, k, l;
	float **matrix2d, ***matrix3d_;

	/* HEADER */

	fwrite(header, sizeof(char), PIPE_HEADERSIZE, stdout);

	fflush(stdout);

	matrix2d = fmalloc2d(datasize[1], datasize[0]);

	switch (axis_option) {
	case 'x':
		for (l = 0; l < datasize[3]; l++) {

			for (k = 0; k < datasize[2]; k++) {

				if (openxyza4d(spectra4d, k, l, matrix2d) != 0)
					goto escape;

				fpwrite2bin(stdout, &(matrix2d[0][0]), datasize[1] * datasize[0]);

				fflush(stdout);
			}
		}

		break;

	case 'y':
		for (l = 0; l < datasize[3]; l++) {

			for (k = 0; k < datasize[2]; k++) {

				if (openxyza4d(spectra4d, k, l, matrix2d) != 0)
					goto escape;

				for (i = 0; i < datasize[0]; i++) {
					for (j = 0; j < datasize[1]; j++) {
						fpwrite2bin(stdout, &(matrix2d[j][i]), 1);
					}
				}

				fflush(stdout);
			}
		}

		break;

	case 'z':
		matrix3d_ = fmalloc3d(datasize[1], datasize[0], datasize[2]);

		for (l = 0; l < datasize[3]; l++) {

			for (k = 0; k < datasize[2]; k++) {

				if (openxyza4d(spectra4d, k, l, matrix2d) != 0)
					goto escape2;

				for (j = 0; j < datasize[1]; j++) {
					for (i = 0; i < datasize[0]; i++) {
						matrix3d_[j][i][k] = matrix2d[j][i];
					}
				}
			}

			for (j = 0; j < datasize[1]; j++) {
				fpwrite2bin(stdout, &(matrix3d_[j][0][0]), datasize[0] * datasize[2]);

				fflush(stdout);
			}
		}

		free_fmatrix3d(matrix3d_);

		break;

	case 'a':
		matrix3d_ = fmalloc3d(datasize[1], datasize[0], datasize[3]);

		for (k = 0; k < datasize[2]; k++) {

			for (l = 0; l < datasize[3]; l++) {

				if (openxyza4d(spectra4d, k, l, matrix2d) != 0)
					goto escape2;

				for (j = 0; j < datasize[1]; j++) {
					for (i = 0; i < datasize[0]; i++) {
						matrix3d_[j][i][l] = matrix2d[j][i];
					}
				}
			}

			for (j = 0; j < datasize[1]; j++) {
				fpwrite2bin(stdout, &(matrix3d_[j][0][0]), datasize[0] * datasize[3]);

				fflush(stdout);
			}
		}

		free_fmatrix3d(matrix3d_);

		break;
	}

	free_fmatrix2d(matrix2d);

	return 0;

	escape2:free_fmatrix3d(matrix3d_);
	escape:free_fmatrix2d(matrix2d);

	return 1;
}
