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

int pushxwnmr2d(char spectra2d[], const char axis_option)
{
	int i, j;
	float **matrix2d;

	matrix2d = fmalloc2d(datasize[1], datasize[0]);

	if (openxwnmr2d(spectra2d, matrix2d) != 0)
		goto escape;

	if (extleft && fabs(spcenter[0] - 4.75) > 0.5) {
		fprintf(stderr, "'--extLeft' option was overriden because of spectral center (%.2f [ppm]).\n", spcenter[0]);
		extleft = 0;
	}

	if (extleft) {
		datasize[0] /= 2;
		spwidth[0] /= 2.0;
		origfreq[0] += spwidth[0];
		spcenter[0] = (origfreq[0] + spwidth[0] / 2.0) / obsfreq[0];

		fwrite2mem(header + 400, spwidth[0]);
		fwrite2mem(header + 264, spcenter[0]);
		fwrite2mem(header + 404, origfreq[0]);
		fwrite2mem(header + 396, (float) (datasize[0]));
		fwrite2mem(header + 876, (float) (datasize[1]));

		cnvhdr(axis_option, 'f');
	}

	fwrite(header, sizeof(char), PIPE_HEADERSIZE, stdout);

	fflush(stdout);

	switch (axis_option) {
	case 'x':
		if (extleft) {
			for (j = 0; j < datasize[1]; j++) {
				for (i = 0; i < datasize[0]; i++) {
					fpwrite2bin(stdout, &(matrix2d[j][i]), 1);
				}
			}
		}

		else
			fpwrite2bin(stdout, &(matrix2d[0][0]), datasize[1] * datasize[0]);

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

	return 1;
}

int pushxwnmr3d(char spectra3d[], const char axis_option)
{
	int i, j, k;
	float ***matrix3d;

	matrix3d = fmalloc3d(datasize[2], datasize[1], datasize[0]);

	if (openxwnmr3d(spectra3d, matrix3d) != 0)
		goto escape;

	if (extleft && fabs(spcenter[0] - 4.75) > 0.5) {
		fprintf(stderr, "'--extLeft' option was overriden because of spectral center (%.2f [ppm]).\n", spcenter[0]);
		extleft = 0;
	}

	if (extleft) {
		datasize[0] /= 2;
		spwidth[0] /= 2.0;
		origfreq[0] += spwidth[0];
		spcenter[0] = (origfreq[0] + spwidth[0] / 2.0) / obsfreq[0];

		fwrite2mem(header + 400, spwidth[0]);
		fwrite2mem(header + 264, spcenter[0]);
		fwrite2mem(header + 404, origfreq[0]);
		fwrite2mem(header + 396, (float) (datasize[0]));
		fwrite2mem(header + 876, (float) (datasize[1]));
		fwrite2mem(header + 60, (float) (datasize[2]));

		cnvhdr(axis_option, 'f');
	}

	fwrite(header, sizeof(char), PIPE_HEADERSIZE, stdout);

	fflush(stdout);

	switch (axis_option) {
	case 'x':
		if (extleft) {
			for (k = 0; k < datasize[2]; k++) {
				for (j = 0; j < datasize[1]; j++) {
					for (i = 0; i < datasize[0]; i++) {
						fpwrite2bin(stdout, &(matrix3d[k][j][i]), 1);
					}
				}

				fflush(stdout);
			}
		}

		else {
			for (k = 0; k < datasize[2]; k++) {
				fpwrite2bin(stdout, &(matrix3d[k][0][0]), datasize[1] * datasize[0]);

				fflush(stdout);
			}
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

	return 1;
}

int pushxwnmr4d(char spectra4d[], const char axis_option)
{
	int i, j, k, l;
	float ****matrix4d;

	matrix4d = fmalloc4d(datasize[3], datasize[2], datasize[1], datasize[0]);

	if (openxwnmr4d(spectra4d, matrix4d) != 0)
		goto escape;

	if (extleft && fabs(spcenter[0] - 4.75) > 0.5) {
		fprintf(stderr, "'--extLeft' option was overriden because of spectral center (%.2f [ppm]).\n", spcenter[0]);
		extleft = 0;
	}

	if (extleft) {
		datasize[0] /= 2;
		spwidth[0] /= 2.0;
		origfreq[0] += spwidth[0];
		spcenter[0] = (origfreq[0] + spwidth[0] / 2.0) / obsfreq[0];

		fwrite2mem(header + 400, spwidth[0]);
		fwrite2mem(header + 264, spcenter[0]);
		fwrite2mem(header + 404, origfreq[0]);
		fwrite2mem(header + 396, (float) (datasize[0]));
		fwrite2mem(header + 876, (float) (datasize[1]));
		fwrite2mem(header + 60, (float) (datasize[2]));
		fwrite2mem(header + 128, (float) (datasize[3]));

		cnvhdr(axis_option, 'f');
	}

	fwrite(header, sizeof(char), PIPE_HEADERSIZE, stdout);

	fflush(stdout);

	switch (axis_option) {
	case 'x':
		if (extleft) {
			for (l = 0; l < datasize[3]; l++) {
				for (k = 0; k < datasize[2]; k++) {
					for (j = 0; j < datasize[1]; j++) {
						for (i = 0; i < datasize[0]; i++) {
							fpwrite2bin(stdout, &(matrix4d[l][k][j][i]), 1);
						}
					}

					fflush(stdout);
				}
			}
		}

		else {
			for (l = 0; l < datasize[3]; l++) {
				for (k = 0; k < datasize[2]; k++) {
					fpwrite2bin(stdout, &(matrix4d[l][k][0][0]), datasize[1] * datasize[0]);

					fflush(stdout);
				}
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

	return 1;
}
