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

int pushaddxeasy2d(char spectra2d1[], char spectra2d2[], const float c1, const float c2, enum_combine_opr opr_code)
{
	int i, j;
	float **matrix2d, **matrix2d_;

	/* HEADER */
	fwrite(header, sizeof(char), PIPE_HEADERSIZE, stdout);

	fflush(stdout);

	matrix2d = fmalloc2d(datasize[1], datasize[0]);
	matrix2d_ = fmalloc2d(datasize[1], datasize[0]);

	if (openxeasy2d(spectra2d1, matrix2d) != 0)
		goto escape;

	for (j = 0; j < datasize[1]; j++) {
		for (i = 0; i < datasize[0]; i++) {
			matrix2d_[j][i] = c1 * matrix2d[j][i];
		}
	}

	if (openxeasy2d(spectra2d2, matrix2d) != 0)
		goto escape;

	switch (opr_code) {
	case COMBINE_ADD:
		for (j = 0; j < datasize[1]; j++) {
			for (i = 0; i < datasize[0]; i++) {
				matrix2d_[j][i] += c2 * matrix2d[j][i];
			}
		}

		break;

	case COMBINE_SUB:
		for (j = 0; j < datasize[1]; j++) {
			for (i = 0; i < datasize[0]; i++) {
				matrix2d_[j][i] -= c2 * matrix2d[j][i];
			}
		}

		break;

	case COMBINE_MUL:
		for (j = 0; j < datasize[1]; j++) {
			for (i = 0; i < datasize[0]; i++) {
				matrix2d_[j][i] *= c2 * matrix2d[j][i];
			}
		}

		break;
	}

	fpwrite2bin(stdout, &(matrix2d_[0][0]), get_data_plane());

	fflush(stdout);

	free_fmatrix2d(matrix2d);
	free_fmatrix2d(matrix2d_);

	return 0;

	escape:free_fmatrix2d(matrix2d);
	free_fmatrix2d(matrix2d_);

	return 1;
}

int pushaddxeasy3d(char spectra3d1[], char spectra3d2[], const float c1, const float c2, enum_combine_opr opr_code)
{
	int i, j, k;
	float **matrix2d_, ***matrix3d, ***matrix3d_;

	/* HEADER */
	fwrite(header, sizeof(char), PIPE_HEADERSIZE, stdout);

	fflush(stdout);

	matrix2d_ = fmalloc2d(datasize[1], datasize[0]);
	matrix3d = fmalloc3d(datasize[2], datasize[1], datasize[0]);
	matrix3d_ = fmalloc3d(datasize[2], datasize[1], datasize[0]);

	if (openxeasy3d(spectra3d1, matrix3d) != 0)
		goto escape;

	if (openxeasy3d(spectra3d2, matrix3d_) != 0)
		goto escape;

	int data_plane = get_data_plane();

	for (k = 0; k < datasize[2]; k++) {

		for (j = 0; j < datasize[1]; j++) {
			for (i = 0; i < datasize[0]; i++) {
				matrix2d_[j][i] = c1 * matrix3d[k][j][i];
			}
		}

		switch (opr_code) {
		case COMBINE_ADD:
			for (j = 0; j < datasize[1]; j++) {
				for (i = 0; i < datasize[0]; i++) {
					matrix2d_[j][i] += c2 * matrix3d_[k][j][i];
				}
			}

			break;

		case COMBINE_SUB:
			for (j = 0; j < datasize[1]; j++) {
				for (i = 0; i < datasize[0]; i++) {
					matrix2d_[j][i] -= c2 * matrix3d_[k][j][i];
				}
			}

			break;

		case COMBINE_MUL:
			for (j = 0; j < datasize[1]; j++) {
				for (i = 0; i < datasize[0]; i++) {
					matrix2d_[j][i] *= c2 * matrix3d_[k][j][i];
				}
			}

			break;
		}

		fpwrite2bin(stdout, &(matrix2d_[0][0]), data_plane);

		fflush(stdout);
	}

	free_fmatrix2d(matrix2d_);
	free_fmatrix3d(matrix3d);
	free_fmatrix3d(matrix3d_);

	return 0;

	escape:free_fmatrix2d(matrix2d_);
	free_fmatrix3d(matrix3d);
	free_fmatrix3d(matrix3d_);

	return 1;
}

int pushaddxeasy4d(char spectra4d1[], char spectra4d2[], const float c1, const float c2, enum_combine_opr opr_code)
{
	int i, j, k, l;
	float **matrix2d_, ****matrix4d, ****matrix4d_;

	/* HEADER */
	fwrite(header, sizeof(char), PIPE_HEADERSIZE, stdout);

	fflush(stdout);

	matrix2d_ = fmalloc2d(datasize[1], datasize[0]);
	matrix4d = fmalloc4d(datasize[3], datasize[2], datasize[1], datasize[0]);
	matrix4d_ = fmalloc4d(datasize[3], datasize[2], datasize[1], datasize[0]);

	if (openxeasy4d(spectra4d1, matrix4d) != 0)
		goto escape;

	if (openxeasy4d(spectra4d2, matrix4d_) != 0)
		goto escape;

	int data_plane = get_data_plane();

	for (l = 0; l < datasize[3]; l++) {

		for (k = 0; k < datasize[2]; k++) {

			for (j = 0; j < datasize[1]; j++) {
				for (i = 0; i < datasize[0]; i++) {
					matrix2d_[j][i] = c1 * matrix4d[l][k][j][i];
				}
			}

			switch (opr_code) {
			case COMBINE_ADD:
				for (j = 0; j < datasize[1]; j++) {
					for (i = 0; i < datasize[0]; i++) {
						matrix2d_[j][i] += c2 * matrix4d_[l][k][j][i];
					}
				}

				break;

			case COMBINE_SUB:
				for (j = 0; j < datasize[1]; j++) {
					for (i = 0; i < datasize[0]; i++) {
						matrix2d_[j][i] -= c2 * matrix4d_[l][k][j][i];
					}
				}

				break;

			case COMBINE_MUL:
				for (j = 0; j < datasize[1]; j++) {
					for (i = 0; i < datasize[0]; i++) {
						matrix2d_[j][i] *= c2 * matrix4d_[l][k][j][i];
					}
				}

				break;
			}

			fpwrite2bin(stdout, &(matrix2d_[0][0]), data_plane);

			fflush(stdout);
		}
	}

	free_fmatrix2d(matrix2d_);
	free_fmatrix4d(matrix4d);
	free_fmatrix4d(matrix4d_);

	return 0;

	escape:free_fmatrix2d(matrix2d_);
	free_fmatrix4d(matrix4d);
	free_fmatrix4d(matrix4d_);

	return 1;
}
