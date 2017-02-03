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

int pushadd2d(char spectra2d1[], char spectra2d2[], const float c1, const float c2, enum_combine_opr opr_code)
{
	int i, j;
	float **matrix2d, **matrix2d_;

	/* HEADER */
	fwrite(header, sizeof(char), PIPE_HEADERSIZE, stdout);

	fflush(stdout);

	matrix2d = fmalloc2d(datasize[1], datasize[0]);
	matrix2d_ = fmalloc2d(datasize[1], datasize[0]);

	if (openxyza2d(spectra2d1, matrix2d) != 0)
		goto escape;

	for (j = 0; j < datasize[1]; j++) {
		for (i = 0; i < datasize[0]; i++) {
			matrix2d_[j][i] = c1 * matrix2d[j][i];
		}
	}

	if (openxyza2d(spectra2d2, matrix2d) != 0)
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

int pushadd3d(char spectra3d1[], char spectra3d2[], const float c1, const float c2, enum_combine_opr opr_code)
{
	int i, j, k;
	float **matrix2d, **matrix2d_;

	/* HEADER */
	fwrite(header, sizeof(char), PIPE_HEADERSIZE, stdout);

	fflush(stdout);

	matrix2d = fmalloc2d(datasize[1], datasize[0]);
	matrix2d_ = fmalloc2d(datasize[1], datasize[0]);

	for (k = 0; k < datasize[2]; k++) {

		if (openxyza3d(spectra3d1, k, matrix2d) != 0)
			goto escape;

		for (j = 0; j < datasize[1]; j++) {
			for (i = 0; i < datasize[0]; i++) {
				matrix2d_[j][i] = c1 * matrix2d[j][i];
			}
		}

		if (openxyza3d(spectra3d2, k, matrix2d) != 0)
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
	}

	free_fmatrix2d(matrix2d);
	free_fmatrix2d(matrix2d_);

	return 0;

	escape:free_fmatrix2d(matrix2d);
	free_fmatrix2d(matrix2d_);

	return 1;
}

int pushadd4d(char spectra4d1[], char spectra4d2[], const float c1, const float c2, enum_combine_opr opr_code)
{
	int i, j, k, l;
	float **matrix2d, **matrix2d_;

	/* HEADER */
	fwrite(header, sizeof(char), PIPE_HEADERSIZE, stdout);

	fflush(stdout);

	matrix2d = fmalloc2d(datasize[1], datasize[0]);
	matrix2d_ = fmalloc2d(datasize[1], datasize[0]);

	for (l = 0; l < datasize[3]; l++) {

		for (k = 0; k < datasize[2]; k++) {

			if (openxyza4d(spectra4d1, k, l, matrix2d) != 0)
				goto escape;

			for (j = 0; j < datasize[1]; j++) {
				for (i = 0; i < datasize[0]; i++) {
					matrix2d_[j][i] = c1 * matrix2d[j][i];
				}
			}

			if (openxyza4d(spectra4d2, k, l, matrix2d) != 0)
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
		}
	}

	free_fmatrix2d(matrix2d);
	free_fmatrix2d(matrix2d_);

	return 0;

	escape:free_fmatrix2d(matrix2d);
	free_fmatrix2d(matrix2d_);

	return 1;
}
