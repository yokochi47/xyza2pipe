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

int pushaddvnmr2d(char monofile1[], char monofile2[], char pardir1[], char pardir2[], const float c1, const float c2,
		enum_combine_opr opr_code)
{
	char string[MAXCHAR];
	int i, j;
	float **matrix2d, **matrix2d_;

	/* HEADER */
	fwrite(header, sizeof(char), PIPE_HEADERSIZE, stdout);

	fflush(stdout);

	matrix2d = fmalloc2d(datasize[1], datasize[0]);
	matrix2d_ = fmalloc2d(datasize[1], datasize[0]);

	if (openvnmr2d(monofile1, pardir1, matrix2d) != 0)
		goto escape;

	sprintf(string, "rm -f %s", monofile1);
	system(string);

	for (j = 0; j < datasize[1]; j++) {
		for (i = 0; i < datasize[0]; i++) {
			matrix2d_[j][i] = c1 * matrix2d[j][i];
		}
	}

	if (openvnmr2d(monofile2, pardir2, matrix2d) != 0)
		goto escape;

	sprintf(string, "rm -f %s", monofile2);
	system(string);

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

	case COMBINE_DIV:
		for (j = 0; j < datasize[1]; j++) {
			for (i = 0; i < datasize[0]; i++) {
				matrix2d_[j][i] /= c2 * matrix2d[j][i];
			}
		}

		break;
	}

	fpwrite2bin(stdout, &(matrix2d_[0][0]), datasize[1] * datasize[0]);

	fflush(stdout);

	free_fmatrix2d(matrix2d);
	free_fmatrix2d(matrix2d_);

	return 0;

	escape:free_fmatrix2d(matrix2d);
	free_fmatrix2d(matrix2d_);

	return 1;
}

int pushaddvnmr3d(char monofile1[], char monofile2[], char pardir1[], char pardir2[], const float c1, const float c2,
		enum_combine_opr opr_code)
{
	char string[MAXCHAR];
	int i, j, k;
	float **matrix2d_, ***matrix3d, ***matrix3d_;

	/* HEADER */
	fwrite(header, sizeof(char), PIPE_HEADERSIZE, stdout);

	fflush(stdout);

	matrix2d_ = fmalloc2d(datasize[1], datasize[0]);
	matrix3d = fmalloc3d(datasize[2], datasize[1], datasize[0]);
	matrix3d_ = fmalloc3d(datasize[2], datasize[1], datasize[0]);

	if (openvnmr3d(monofile1, pardir1, matrix3d) != 0)
		goto escape;

	sprintf(string, "rm -f %s", monofile1);
	system(string);

	if (openvnmr3d(monofile2, pardir2, matrix3d_) != 0)
		goto escape;

	sprintf(string, "rm -f %s", monofile1);
	system(string);

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

		case COMBINE_DIV:
			for (j = 0; j < datasize[1]; j++) {
				for (i = 0; i < datasize[0]; i++) {
					matrix2d_[j][i] /= c2 * matrix3d_[k][j][i];
				}
			}

			break;
		}

		fpwrite2bin(stdout, &(matrix2d_[0][0]), datasize[1] * datasize[0]);

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

int pushaddvnmr4d(char monofile1[], char monofile2[], char pardir1[], char pardir2[], const float c1, const float c2,
		enum_combine_opr opr_code)
{
	char string[MAXCHAR];
	int i, j, k, l;
	float **matrix2d_, ****matrix4d, ****matrix4d_;

	/* HEADER */
	fwrite(header, sizeof(char), PIPE_HEADERSIZE, stdout);

	fflush(stdout);

	matrix2d_ = fmalloc2d(datasize[1], datasize[0]);
	matrix4d = fmalloc4d(datasize[3], datasize[2], datasize[1], datasize[0]);
	matrix4d_ = fmalloc4d(datasize[3], datasize[2], datasize[1], datasize[0]);

	if (openvnmr4d(monofile1, pardir1, matrix4d) != 0)
		goto escape;

	sprintf(string, "rm -f %s", monofile1);
	system(string);

	if (openvnmr4d(monofile2, pardir2, matrix4d_) != 0)
		goto escape;

	sprintf(string, "rm -f %s", monofile1);
	system(string);

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

			case COMBINE_DIV:
				for (j = 0; j < datasize[1]; j++) {
					for (i = 0; i < datasize[0]; i++) {
						matrix2d_[j][i] /= c2 * matrix4d_[l][k][j][i];
					}
				}

				break;
			}

			fpwrite2bin(stdout, &(matrix2d_[0][0]), datasize[1] * datasize[0]);

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
