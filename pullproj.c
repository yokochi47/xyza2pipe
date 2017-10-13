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

int pullproj2d(char filename[], const int abs_mode)
{
	FILE *fp = (filename[0] == 0 ? stdout : NULL);
	int i, j, k, count = 0, axisorder[4] = { 0 }, orderaxis[4];
	int data_plane = get_data_plane();
	float *vp;
	float **matrix2d;

	if (fp != stdout)
		make_dir(filename);

	set_clean_string(clean_string);

	vp = (float *) ((void *) (&header[96]));
	if (*vp == 2.0)
		axisorder[0] = 0;
	if (*vp == 1.0)
		axisorder[0] = 1;
	if (*vp == 3.0)
		axisorder[0] = 2;
	if (*vp == 4.0)
		axisorder[0] = 3;

	vp = (float *) ((void *) (&header[100]));
	if (*vp == 2.0)
		axisorder[1] = 0;
	if (*vp == 1.0)
		axisorder[1] = 1;
	if (*vp == 3.0)
		axisorder[1] = 2;
	if (*vp == 4.0)
		axisorder[1] = 3;

	vp = (float *) ((void *) (&header[104]));
	if (*vp == 2.0)
		axisorder[2] = 0;
	if (*vp == 1.0)
		axisorder[2] = 1;
	if (*vp == 3.0)
		axisorder[2] = 2;
	if (*vp == 4.0)
		axisorder[2] = 3;

	vp = (float *) ((void *) (&header[108]));
	if (*vp == 2.0)
		axisorder[3] = 0;
	if (*vp == 1.0)
		axisorder[3] = 1;
	if (*vp == 3.0)
		axisorder[3] = 2;
	if (*vp == 4.0)
		axisorder[3] = 3;

	for (j = 0; j < 4; j++) {
		for (k = 0; k < 4; k++) {
			if (axisorder[k] == j)
				orderaxis[j] = k;
		}
	}

	if (leftcar) {
		if (orderaxis[0] == 0)
			fwrite2mem(header + 264, spcenter[0] + leftcar * spwidth[0] / 2.0 / obsfreq[0]);
		else if (orderaxis[1] == 0)
			fwrite2mem(header + 268, spcenter[0] + leftcar * spwidth[0] / 2.0 / obsfreq[0]);
	}

	matrix2d = fmalloc2d(datasize[1], datasize[0]);

	fprintf(stderr, "%s", clean_string);
	fprintf(stderr, "Receiving %d of %d ( -> pipe2proj) ...\r", ++count, 1);

	if (openpipe2d(matrix2d) != 0)
		goto escape;

	fprintf(stderr, "%s", clean_string);
	if (fp != stdout)
		fprintf(stderr, "Writing %s <- pipe2proj ...\r", filename);

	if (swapdata != 0)
		swapbyte(sizeof(float), data_plane * sizeof(float), (char *) (&(matrix2d[0][0])));

	if (abs_mode) {
		for (j = 0; j < datasize[1]; j++) {
			for (i = 0; i < datasize[0]; i++) {
				matrix2d[j][i] = fabs(matrix2d[j][i]);
			}
		}
	}

	if (fp != stdout && (fp = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "Permission denied to write %s.\n", filename);
		exit(EXIT_FAILURE);
	}

	fwrite2mem(header + 36, 2.0);

	fwrite(header, sizeof(char), PIPE_HEADERSIZE, fp);

	fpwrite2bin(fp, &(matrix2d[0][0]), data_plane);

	if (fp != stdout)
		fclose(fp);
	else
		fflush(stdout);

	free_fmatrix2d(matrix2d);

	fputc('\n', stderr);

	return 0;

	escape:free_fmatrix2d(matrix2d);

	return EXIT_FAILURE;
}

int pullproj3d(char filename[], const int abs_mode)
{
	FILE *fp = (filename[0] == 0 ? stdout : NULL);
	int i, j, k, count = 0, axisorder[4] = { 0 }, orderaxis[4];
	int data_plane = get_data_plane(), rcv_planes;
	float *vp, **matrix2d, **matrix2d_;

	if (fp != stdout)
		make_dir(filename);

	set_clean_string(clean_string);

	vp = (float *) ((void *) (&header[96]));
	if (*vp == 2.0)
		axisorder[0] = 0;
	if (*vp == 1.0)
		axisorder[0] = 1;
	if (*vp == 3.0)
		axisorder[0] = 2;
	if (*vp == 4.0)
		axisorder[0] = 3;

	vp = (float *) ((void *) (&header[100]));
	if (*vp == 2.0)
		axisorder[1] = 0;
	if (*vp == 1.0)
		axisorder[1] = 1;
	if (*vp == 3.0)
		axisorder[1] = 2;
	if (*vp == 4.0)
		axisorder[1] = 3;

	vp = (float *) ((void *) (&header[104]));
	if (*vp == 2.0)
		axisorder[2] = 0;
	if (*vp == 1.0)
		axisorder[2] = 1;
	if (*vp == 3.0)
		axisorder[2] = 2;
	if (*vp == 4.0)
		axisorder[2] = 3;

	vp = (float *) ((void *) (&header[108]));
	if (*vp == 2.0)
		axisorder[3] = 0;
	if (*vp == 1.0)
		axisorder[3] = 1;
	if (*vp == 3.0)
		axisorder[3] = 2;
	if (*vp == 4.0)
		axisorder[3] = 3;

	for (j = 0; j < 4; j++) {
		for (k = 0; k < 4; k++) {
			if (axisorder[k] == j)
				orderaxis[j] = k;
		}
	}

	if (leftcar) {
		if (orderaxis[0] == 0)
			fwrite2mem(header + 264, spcenter[0] + leftcar * spwidth[0] / 2.0 / obsfreq[0]);
		else if (orderaxis[1] == 0)
			fwrite2mem(header + 268, spcenter[0] + leftcar * spwidth[0] / 2.0 / obsfreq[0]);
		else if (orderaxis[2] == 0)
			fwrite2mem(header + 272, spcenter[0] + leftcar * spwidth[0] / 2.0 / obsfreq[0]);
	}

	matrix2d = fmalloc2d(datasize[1], datasize[0]);
	matrix2d_ = fmalloc2d(datasize[1], datasize[0]);

	for (j = 0; j < datasize[1]; j++) {
		for (i = 0; i < datasize[0]; i++) {
			matrix2d_[j][i] = 0.0;
		}
	}

	rcv_planes = datasize[2];

	for (k = 0; k < datasize[2]; k++) {
		fprintf(stderr, "%s", clean_string);
		fprintf(stderr, "Receiving %d of %d ( -> pipe2proj) ...\r", ++count, rcv_planes);

		if (openpipe2d(matrix2d) != 0)
			goto escape;

		fprintf(stderr, "%s", clean_string);
		if (fp != stdout)
			fprintf(stderr, "Writing %s <- pipe2proj ...\r", filename);

		if (swapdata != 0)
			swapbyte(sizeof(float), data_plane * sizeof(float), (char *) (&(matrix2d[0][0])));

		if (abs_mode) {
			for (j = 0; j < datasize[1]; j++) {
				for (i = 0; i < datasize[0]; i++) {
					matrix2d_[j][i] += fabs(matrix2d[j][i]);
				}
			}
		}

		else {
			for (j = 0; j < datasize[1]; j++) {
				for (i = 0; i < datasize[0]; i++) {
					matrix2d_[j][i] += matrix2d[j][i];
				}
			}
		}
	}

	if (fp != stdout && (fp = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "Permission denied to write %s.\n", filename);
		exit(EXIT_FAILURE);
	}

	fwrite2mem(header + 36, 2.0);

	fwrite(header, sizeof(char), PIPE_HEADERSIZE, fp);

	fpwrite2bin(fp, &(matrix2d_[0][0]), data_plane);

	if (fp != stdout)
		fclose(fp);
	else
		fflush(stdout);

	free_fmatrix2d(matrix2d);
	free_fmatrix2d(matrix2d_);

	fputc('\n', stderr);

	return 0;

	escape:free_fmatrix2d(matrix2d);
	free_fmatrix2d(matrix2d_);

	return EXIT_FAILURE;
}

int pullproj4d(char filename[], const int abs_mode)
{
	FILE *fp = (filename[0] == 0 ? stdout : NULL);
	int i, j, k, l, count = 0, axisorder[4] = { 0 }, orderaxis[4];
	int data_plane = get_data_plane(), rcv_planes;
	float *vp, **matrix2d, **matrix2d_;

	if (fp != stdout)
		make_dir(filename);

	set_clean_string(clean_string);

	vp = (float *) ((void *) (&header[96]));
	if (*vp == 2.0)
		axisorder[0] = 0;
	if (*vp == 1.0)
		axisorder[0] = 1;
	if (*vp == 3.0)
		axisorder[0] = 2;
	if (*vp == 4.0)
		axisorder[0] = 3;

	vp = (float *) ((void *) (&header[100]));
	if (*vp == 2.0)
		axisorder[1] = 0;
	if (*vp == 1.0)
		axisorder[1] = 1;
	if (*vp == 3.0)
		axisorder[1] = 2;
	if (*vp == 4.0)
		axisorder[1] = 3;

	vp = (float *) ((void *) (&header[104]));
	if (*vp == 2.0)
		axisorder[2] = 0;
	if (*vp == 1.0)
		axisorder[2] = 1;
	if (*vp == 3.0)
		axisorder[2] = 2;
	if (*vp == 4.0)
		axisorder[2] = 3;

	vp = (float *) ((void *) (&header[108]));
	if (*vp == 2.0)
		axisorder[3] = 0;
	if (*vp == 1.0)
		axisorder[3] = 1;
	if (*vp == 3.0)
		axisorder[3] = 2;
	if (*vp == 4.0)
		axisorder[3] = 3;

	for (j = 0; j < 4; j++) {
		for (k = 0; k < 4; k++) {
			if (axisorder[k] == j)
				orderaxis[j] = k;
		}
	}

	if (leftcar) {
		if (orderaxis[0] == 0)
			fwrite2mem(header + 264, spcenter[0] + leftcar * spwidth[0] / 2.0 / obsfreq[0]);
		else if (orderaxis[1] == 0)
			fwrite2mem(header + 268, spcenter[0] + leftcar * spwidth[0] / 2.0 / obsfreq[0]);
		else if (orderaxis[2] == 0)
			fwrite2mem(header + 272, spcenter[0] + leftcar * spwidth[0] / 2.0 / obsfreq[0]);
		else if (orderaxis[3] == 0)
			fwrite2mem(header + 276, spcenter[0] + leftcar * spwidth[0] / 2.0 / obsfreq[0]);
	}

	matrix2d = fmalloc2d(datasize[1], datasize[0]);
	matrix2d_ = fmalloc2d(datasize[1], datasize[0]);

	for (j = 0; j < datasize[1]; j++) {
		for (i = 0; i < datasize[0]; i++) {
			matrix2d_[j][i] = 0.0;
		}
	}

	rcv_planes = datasize[2] * datasize[3];

	for (l = 0; l < datasize[3]; l++) {
		for (k = 0; k < datasize[2]; k++) {
			fprintf(stderr, "%s", clean_string);
			fprintf(stderr, "Receiving %d of %d ( -> pipe2proj) ...\r", ++count, rcv_planes);

			if (openpipe2d(matrix2d) != 0)
				goto escape;

			fprintf(stderr, "%s", clean_string);
			if (fp != stdout)
				fprintf(stderr, "Writing %s <- pipe2proj ...\r", filename);

			if (swapdata != 0)
				swapbyte(sizeof(float), data_plane * sizeof(float), (char *) (&(matrix2d[0][0])));

			if (abs_mode) {
				for (j = 0; j < datasize[1]; j++) {
					for (i = 0; i < datasize[0]; i++) {
						matrix2d_[j][i] += fabs(matrix2d[j][i]);
					}
				}
			}

			else {
				for (j = 0; j < datasize[1]; j++) {
					for (i = 0; i < datasize[0]; i++) {
						matrix2d_[j][i] += matrix2d[j][i];
					}
				}
			}
		}
	}

	if (fp != stdout && (fp = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "Permission denied to write %s.\n", filename);
		exit(EXIT_FAILURE);
	}

	fwrite2mem(header + 36, 2.0);

	fwrite(header, sizeof(char), PIPE_HEADERSIZE, fp);

	fpwrite2bin(fp, &(matrix2d_[0][0]), data_plane);

	if (fp != stdout)
		fclose(fp);
	else
		fflush(stdout);

	free_fmatrix2d(matrix2d);
	free_fmatrix2d(matrix2d_);

	fputc('\n', stderr);

	return 0;

	escape:free_fmatrix2d(matrix2d);
	free_fmatrix2d(matrix2d_);

	return EXIT_FAILURE;
}
