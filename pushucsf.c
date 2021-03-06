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

int pushucsf2d(char spectra2d[], const char axis_option)
{
	int i, j;
	float **matrix2d;

	if (datasize[0] > datasize_orig[0] || datasize[1] > datasize_orig[1])
		goto resize;

	/* HEADER */
	fwrite(header, sizeof(char), PIPE_HEADERSIZE, stdout);

	fflush(stdout);

	matrix2d = fmalloc2d(datasize[1], datasize[0]);

	if (openucsf2d(spectra2d, matrix2d) != 0)
		goto escape;

	switch (axis_option) {
	case 'x':
		fpwrite2bin(stdout, &(matrix2d[0][0]), get_data_plane());

		fflush(stdout);

		break;

	case 'y':
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

	return EXIT_FAILURE;

	resize:cnvhdr(axis_option, 'b');

	fwrite2mem(header + 396, (float) (datasize_orig[0]));
	fwrite2mem(header + 876, (float) (datasize_orig[1]));

	cnvhdr(axis_option, 'f');

	/* HEADER */
	fwrite(header, sizeof(char), PIPE_HEADERSIZE, stdout);

	fflush(stdout);

	matrix2d = fmalloc2d(datasize[1], datasize[0]);

	if (openucsf2d(spectra2d, matrix2d) != 0)
		goto escape;

	switch (axis_option) {
	case 'x':
		for (j = 0; j < datasize_orig[1]; j++) {
			for (i = 0; i < datasize_orig[0]; i++) {
				fpwrite2bin(stdout, &(matrix2d[j][i]), 1);
			}
		}

		fflush(stdout);

		break;

	case 'y':
		for (i = 0; i < datasize_orig[0]; i++) {
			for (j = 0; j < datasize_orig[1]; j++) {
				fpwrite2bin(stdout, &(matrix2d[j][i]), 1);
			}
		}

		fflush(stdout);

		break;
	}

	free_fmatrix2d(matrix2d);

	return 0;
}

int pushucsf3d(char spectra3d[], const char axis_option)
{
	int i, j, k;
	int data_plane = get_data_plane();
	float ***matrix3d;

	if (datasize[0] > datasize_orig[0] || datasize[1] > datasize_orig[1] || datasize[2] > datasize_orig[2])
		goto resize;

	/* HEADER */
	fwrite(header, sizeof(char), PIPE_HEADERSIZE, stdout);

	fflush(stdout);

	matrix3d = fmalloc3d(datasize[2], datasize[1], datasize[0]);

	if (openucsf3d(spectra3d, matrix3d) != 0)
		goto escape;

	switch (axis_option) {
	case 'x':
		for (k = 0; k < datasize[2]; k++) {
			fpwrite2bin(stdout, &(matrix3d[k][0][0]), data_plane);

			fflush(stdout);
		}

		break;

	case 'y':
		for (k = 0; k < datasize[2]; k++) {
			for (i = 0; i < datasize[0]; i++) {
				for (j = 0; j < datasize[1]; j++) {
					fpwrite2bin(stdout, &(matrix3d[k][j][i]), 1);
				}
			}

			fflush(stdout);
		}

		break;

	case 'z':
		for (j = 0; j < datasize[1]; j++) {
			for (i = 0; i < datasize[0]; i++) {
				for (k = 0; k < datasize[2]; k++) {
					fpwrite2bin(stdout, &(matrix3d[k][j][i]), 1);
				}
			}

			fflush(stdout);
		}

		break;
	}

	free_fmatrix3d(matrix3d);

	return 0;

	escape:free_fmatrix3d(matrix3d);

	return EXIT_FAILURE;

	resize:cnvhdr(axis_option, 'b');

	fwrite2mem(header + 396, (float) (datasize_orig[0]));
	fwrite2mem(header + 876, (float) (datasize_orig[1]));
	fwrite2mem(header + 60, (float) (datasize_orig[2]));

	fwrite2mem(header + 1768, (float) get_orig_indirect_planes());

	cnvhdr(axis_option, 'f');

	/* HEADER */
	fwrite(header, sizeof(char), PIPE_HEADERSIZE, stdout);

	fflush(stdout);

	matrix3d = fmalloc3d(datasize[2], datasize[1], datasize[0]);

	if (openucsf3d(spectra3d, matrix3d) != 0)
		goto escape;

	switch (axis_option) {
	case 'x':
		for (k = 0; k < datasize_orig[2]; k++) {
			for (j = 0; j < datasize_orig[1]; j++) {
				for (i = 0; i < datasize_orig[0]; i++) {
					fpwrite2bin(stdout, &(matrix3d[k][j][i]), 1);
				}
			}

			fflush(stdout);
		}

		break;

	case 'y':
		for (k = 0; k < datasize_orig[2]; k++) {
			for (i = 0; i < datasize_orig[0]; i++) {
				for (j = 0; j < datasize_orig[1]; j++) {
					fpwrite2bin(stdout, &(matrix3d[k][j][i]), 1);
				}
			}

			fflush(stdout);
		}

		break;

	case 'z':
		for (j = 0; j < datasize_orig[1]; j++) {
			for (i = 0; i < datasize_orig[0]; i++) {
				for (k = 0; k < datasize_orig[2]; k++) {
					fpwrite2bin(stdout, &(matrix3d[k][j][i]), 1);
				}
			}

			fflush(stdout);
		}

		break;
	}

	free_fmatrix3d(matrix3d);

	return 0;
}

int pushucsf4d(char spectra4d[], const char axis_option)
{
	int i, j, k, l;
	int data_plane = get_data_plane();
	float ****matrix4d;

	if (datasize[0] > datasize_orig[0] || datasize[1] > datasize_orig[1] || datasize[2] > datasize_orig[2] || datasize[3] > datasize_orig[3])
		goto resize;

	/* HEADER */
	fwrite(header, sizeof(char), PIPE_HEADERSIZE, stdout);

	fflush(stdout);

	matrix4d = fmalloc4d(datasize[3], datasize[2], datasize[1], datasize[0]);

	if (openucsf4d(spectra4d, matrix4d) != 0)
		goto escape;

	switch (axis_option) {
	case 'x':
		for (l = 0; l < datasize[3]; l++) {
			for (k = 0; k < datasize[2]; k++) {
				fpwrite2bin(stdout, &(matrix4d[l][k][0][0]), data_plane);

				fflush(stdout);
			}
		}

		break;

	case 'y':
		for (l = 0; l < datasize[3]; l++) {
			for (k = 0; k < datasize[2]; k++) {
				for (i = 0; i < datasize[0]; i++) {
					for (j = 0; j < datasize[1]; j++) {
						fpwrite2bin(stdout, &(matrix4d[l][k][j][i]), 1);
					}
				}

				fflush(stdout);
			}
		}

		break;

	case 'z':
		for (l = 0; l < datasize[3]; l++) {
			for (j = 0; j < datasize[1]; j++) {
				for (i = 0; i < datasize[0]; i++) {
					for (k = 0; k < datasize[2]; k++) {
						fpwrite2bin(stdout, &(matrix4d[l][k][j][i]), 1);
					}
				}

				fflush(stdout);
			}
		}

		break;

	case 'a':
		for (k = 0; k < datasize[2]; k++) {
			for (j = 0; j < datasize[1]; j++) {
				for (i = 0; i < datasize[0]; i++) {
					for (l = 0; l < datasize[3]; l++) {
						fpwrite2bin(stdout, &(matrix4d[l][k][j][i]), 1);
					}
				}

				fflush(stdout);
			}
		}

		break;
	}

	free_fmatrix4d(matrix4d);

	return 0;

	escape:free_fmatrix4d(matrix4d);

	return EXIT_FAILURE;

	resize:cnvhdr(axis_option, 'b');

	fwrite2mem(header + 396, (float) (datasize_orig[0]));
	fwrite2mem(header + 876, (float) (datasize_orig[1]));
	fwrite2mem(header + 60, (float) (datasize_orig[2]));
	fwrite2mem(header + 128, (float) (datasize_orig[3]));

	fwrite2mem(header + 1768, (float) get_orig_indirect_planes());

	cnvhdr(axis_option, 'f');

	/* HEADER */
	fwrite(header, sizeof(char), PIPE_HEADERSIZE, stdout);

	fflush(stdout);

	/* HEADER */
	fwrite(header, sizeof(char), PIPE_HEADERSIZE, stdout);

	fflush(stdout);

	matrix4d = fmalloc4d(datasize[3], datasize[2], datasize[1], datasize[0]);

	if (openucsf4d(spectra4d, matrix4d) != 0)
		goto escape;

	switch (axis_option) {
	case 'x':
		for (l = 0; l < datasize_orig[3]; l++) {
			for (k = 0; k < datasize_orig[2]; k++) {
				for (j = 0; j < datasize_orig[1]; j++) {
					for (i = 0; i < datasize_orig[0]; i++) {
						fpwrite2bin(stdout, &(matrix4d[l][k][j][i]), 1);
					}
				}

				fflush(stdout);
			}
		}

		break;

	case 'y':
		for (l = 0; l < datasize_orig[3]; l++) {
			for (k = 0; k < datasize_orig[2]; k++) {
				for (i = 0; i < datasize_orig[0]; i++) {
					for (j = 0; j < datasize_orig[1]; j++) {
						fpwrite2bin(stdout, &(matrix4d[l][k][j][i]), 1);
					}
				}

				fflush(stdout);
			}
		}

		break;

	case 'z':
		for (l = 0; l < datasize_orig[3]; l++) {
			for (j = 0; j < datasize_orig[1]; j++) {
				for (i = 0; i < datasize_orig[0]; i++) {
					for (k = 0; k < datasize_orig[2]; k++) {
						fpwrite2bin(stdout, &(matrix4d[l][k][j][i]), 1);
					}
				}

				fflush(stdout);
			}
		}

		break;

	case 'a':
		for (k = 0; k < datasize_orig[2]; k++) {
			for (j = 0; j < datasize_orig[1]; j++) {
				for (i = 0; i < datasize_orig[0]; i++) {
					for (l = 0; l < datasize_orig[3]; l++) {
						fpwrite2bin(stdout, &(matrix4d[l][k][j][i]), 1);
					}
				}

				fflush(stdout);
			}
		}

		break;
	}

	free_fmatrix4d(matrix4d);

	return 0;
}
