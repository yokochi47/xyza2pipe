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

int pullxyza2d(char spectra2d[], const char axis_option)
{
	FILE *fp = (spectra2d[0] == 0 ? stdout : NULL);
	int i, j, k, count = 0, axisorder[4] = { 0 }, orderaxis[4];
	float *vp, **matrix2d;

	if (fp != stdout)
		make_dir(spectra2d);

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

	switch (axis_option) {
	case 'x':
		fprintf(stderr, "%s", clean_string);
		fprintf(stderr, "Receiving %d of %d ( -> pipe2xyza) ...\r", ++count, 1);

		if (openpipe2d(matrix2d) != 0)
			goto escape;

		if (fp != stdout && (fp = fopen(spectra2d, "w")) == NULL) {
			fprintf(stderr, "Permission denied to write %s.\n", spectra2d);
			exit(EXIT_FAILURE);
		}

		fprintf(stderr, "%s", clean_string);
		if (fp != stdout)
			fprintf(stderr, "Writing %s <- pipe2xyza ...\r", spectra2d);

		fwrite(header, sizeof(char), PIPE_HEADERSIZE, fp);

		fpwrite2bin_swap(fp, &(matrix2d[0][0]), get_data_plane(), swapdata);

		if (fp != stdout)
			fclose(fp);
		else
			fflush(stdout);

		break;

	case 'y':
		fprintf(stderr, "%s", clean_string);
		fprintf(stderr, "Receiving %d of %d ( -> pipe2xyza) ...\r", ++count, 1);

		if (openpipe2d(matrix2d) != 0)
			goto escape;

		if (fp != stdout && (fp = fopen(spectra2d, "w")) == NULL) {
			fprintf(stderr, "Permission denied to write %s.\n", spectra2d);
			exit(EXIT_FAILURE);
		}

		fprintf(stderr, "%s", clean_string);
		if (fp != stdout)
			fprintf(stderr, "Writing %s <- pipe2xyza ...\r", spectra2d);

		fwrite(header, sizeof(char), PIPE_HEADERSIZE, fp);

		for (i = 0; i < datasize[0]; i++) {
			for (j = 0; j < datasize[1]; j++) {
				fpwrite2bin_swap(fp, &(matrix2d[j][i]), 1, swapdata);
			}
		}

		if (fp != stdout)
			fclose(fp);
		else
			fflush(stdout);

		break;
	}

	free_fmatrix2d(matrix2d);

	fputc('\n', stderr);

	return 0;

	escape:free_fmatrix2d(matrix2d);

	return EXIT_FAILURE;
}

int pullxyza3d(char spectra3d[], const char axis_option)
{
	FILE *fp = (spectra3d[0] == 0 ? stdout : NULL);
	char filename[MAXLONGNAME];
	int i, j, k, count = 0, axisorder[4] = { 0 }, orderaxis[4];
	int data_plane = get_data_plane(), rcv_planes, dataw_plane;
	float *vp, **matrix2d, ***matrix3d_;

	if (fp != stdout)
		make_dir(spectra3d);

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

	rcv_planes = datasize[2];

	switch (axis_option) {
	case 'x':
	case 'a':
		for (k = 0; k < datasize[2]; k++) {
			fprintf(stderr, "%s", clean_string);
			fprintf(stderr, "Receiving %d of %d ( -> pipe2xyza) ...\r", ++count, rcv_planes);

			if (openpipe2d(matrix2d) != 0)
				goto escape;

			if (fp != stdout) {
				sprintf(filename, spectra3d, k + 1);

				if ((fp = fopen(filename, "w")) == NULL) {
					fprintf(stderr, "Permission denied to write %s.\n", filename);
					exit(EXIT_FAILURE);
				}
			}

			fprintf(stderr, "%s", clean_string);
			if (fp != stdout)
				fprintf(stderr, "Writing %s <- pipe2xyza ...\r", filename);

			if (k == 0 || fp != stdout)
				fwrite(header, sizeof(char), PIPE_HEADERSIZE, fp);

			fpwrite2bin_swap(fp, &(matrix2d[0][0]), data_plane, swapdata);

			if (fp != stdout)
				fclose(fp);
			else
				fflush(stdout);
		}

		break;

	case 'y':
		for (k = 0; k < datasize[2]; k++) {
			fprintf(stderr, "%s", clean_string);
			fprintf(stderr, "Receiving %d of %d ( -> pipe2xyza) ...\r", ++count, rcv_planes);

			if (openpipe2d(matrix2d) != 0)
				goto escape;

			if (fp != stdout) {
				sprintf(filename, spectra3d, k + 1);

				if ((fp = fopen(filename, "w")) == NULL) {
					fprintf(stderr, "Permission denied to write %s.\n", filename);
					exit(EXIT_FAILURE);
				}
			}

			fprintf(stderr, "%s", clean_string);
			if (fp != stdout)
				fprintf(stderr, "Writing %s <- pipe2xyza ...\r", filename);

			if (k == 0 || fp != stdout)
				fwrite(header, sizeof(char), PIPE_HEADERSIZE, fp);

			for (i = 0; i < datasize[0]; i++) {
				for (j = 0; j < datasize[1]; j++) {
					fpwrite2bin_swap(fp, &(matrix2d[j][i]), 1, swapdata);
				}
			}

			if (fp != stdout)
				fclose(fp);
			else
				fflush(stdout);
		}

		break;

	case 'z':
		matrix3d_ = fmalloc3d(datasize[0], datasize[2], datasize[1]);

		for (k = 0; k < datasize[2]; k++) {
			fprintf(stderr, "%s", clean_string);
			fprintf(stderr, "Receiving %d of %d ( -> pipe2xyza) ...\r", ++count, rcv_planes);

			if (openpipe2d(matrix2d) != 0)
				goto escape2;

			for (i = 0; i < datasize[0]; i++) {
				for (j = 0; j < datasize[1]; j++) {
					matrix3d_[i][k][j] = matrix2d[j][i];
				}
			}
		}

		dataw_plane = datasize[2] * datasize[1];

		for (i = 0; i < datasize[0]; i++) {

			if (fp != stdout) {
				sprintf(filename, spectra3d, i + 1);

				if ((fp = fopen(filename, "w")) == NULL) {
					fprintf(stderr, "Permission denied to write %s.\n", filename);
					exit(EXIT_FAILURE);
				}
			}

			fprintf(stderr, "%s", clean_string);
			if (fp != stdout)
				fprintf(stderr, "Writing %s <- pipe2xyza ...\r", filename);

			if (i == 0 || fp != stdout)
				fwrite(header, sizeof(char), PIPE_HEADERSIZE, fp);

			fpwrite2bin_swap(fp, &(matrix3d_[i][0][0]), dataw_plane, swapdata);

			if (fp != stdout)
				fclose(fp);
			else
				fflush(stdout);
		}

		free_fmatrix3d(matrix3d_);

		break;
	}

	free_fmatrix2d(matrix2d);

	fputc('\n', stderr);

	return 0;

	escape2:free_fmatrix3d(matrix3d_);

	escape:free_fmatrix2d(matrix2d);

	return EXIT_FAILURE;
}

int pullxyza4d(char spectra4d[], const char axis_option)
{
	FILE *fp = (spectra4d[0] == 0 ? stdout : NULL);
	char filename[MAXLONGNAME];
	int i, j, k, l, count = 0, axisorder[4] = { 0 }, orderaxis[4];
	int data_plane = get_data_plane(), rcv_planes, dataw_plane;
	float *vp, **matrix2d, ***matrix3d_;

	if (fp != stdout)
		make_dir(spectra4d);

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

	rcv_planes = datasize[2] * datasize[3];
	dataw_plane = datasize[2] * datasize[1];

	switch (axis_option) {
	case 'x':
		for (l = 0; l < datasize[3]; l++) {
			for (k = 0; k < datasize[2]; k++) {
				fprintf(stderr, "%s", clean_string);
				fprintf(stderr, "Receiving %d of %d ( -> pipe2xyza) ...\r", ++count, rcv_planes);

				if (openpipe2d(matrix2d) != 0)
					goto escape;

				if (fp != stdout) {
					sprintf(filename, spectra4d, l + 1, k + 1);

					if ((fp = fopen(filename, "w")) == NULL) {
						fprintf(stderr, "Permission denied to write %s.\n", filename);
						exit(EXIT_FAILURE);
					}
				}

				fprintf(stderr, "%s", clean_string);
				if (fp != stdout)
					fprintf(stderr, "Writing %s <- pipe2xyza ...\r", filename);

				if (l + k == 0 || fp != stdout)
					fwrite(header, sizeof(char), PIPE_HEADERSIZE, fp);

				fpwrite2bin_swap(fp, &(matrix2d[0][0]), data_plane, swapdata);

				if (fp != stdout)
					fclose(fp);
				else
					fflush(stdout);
			}
		}

		break;

	case 'y':
		for (l = 0; l < datasize[3]; l++) {
			for (k = 0; k < datasize[2]; k++) {
				fprintf(stderr, "%s", clean_string);
				fprintf(stderr, "Receiving %d of %d ( -> pipe2xyza) ...\r", ++count, rcv_planes);

				if (openpipe2d(matrix2d) != 0)
					goto escape;

				if (fp != stdout) {
					sprintf(filename, spectra4d, l + 1, k + 1);

					if ((fp = fopen(filename, "w")) == NULL) {
						fprintf(stderr, "Permission denied to write %s.\n", filename);
						exit(EXIT_FAILURE);
					}
				}

				fprintf(stderr, "%s", clean_string);
				if (fp != stdout)
					fprintf(stderr, "Writing %s <- pipe2xyza ...\r", filename);

				if (l + k == 0 || fp != stdout)
					fwrite(header, sizeof(char), PIPE_HEADERSIZE, fp);

				for (i = 0; i < datasize[0]; i++) {
					for (j = 0; j < datasize[1]; j++) {
						fpwrite2bin_swap(fp, &(matrix2d[j][i]), 1, swapdata);
					}
				}

				if (fp != stdout)
					fclose(fp);
				else
					fflush(stdout);
			}
		}

		break;

	case 'z':
		matrix3d_ = fmalloc3d(datasize[0], datasize[2], datasize[1]);

		for (l = 0; l < datasize[3]; l++) {

			for (k = 0; k < datasize[2]; k++) {
				fprintf(stderr, "%s", clean_string);
				fprintf(stderr, "Receiving %d of %d ( -> pipe2xyza) ...\r", ++count, rcv_planes);

				if (openpipe2d(matrix2d) != 0)
					goto escape2;

				for (i = 0; i < datasize[0]; i++) {
					for (j = 0; j < datasize[1]; j++) {
						matrix3d_[i][k][j] = matrix2d[j][i];
					}
				}
			}

			for (i = 0; i < datasize[0]; i++) {

				if (fp != stdout) {
					sprintf(filename, spectra4d, l + 1, i + 1);

					if ((fp = fopen(filename, "w")) == NULL) {
						fprintf(stderr, "Permission denied to write %s.\n", filename);
						exit(EXIT_FAILURE);
					}
				}

				fprintf(stderr, "%s", clean_string);
				if (fp != stdout)
					fprintf(stderr, "Writing %s <- pipe2xyza ...\r", filename);

				if (l + k == 0 || fp != stdout)
					fwrite(header, sizeof(char), PIPE_HEADERSIZE, fp);

				fpwrite2bin_swap(fp, &(matrix3d_[i][0][0]), dataw_plane, swapdata);

				if (fp != stdout)
					fclose(fp);
				else
					fflush(stdout);
			}
		}

		free_fmatrix3d(matrix3d_);

		break;

	case 'a':
		matrix3d_ = fmalloc3d(datasize[0], datasize[2], datasize[1]);

		for (l = 0; l < datasize[3]; l++) {

			for (k = 0; k < datasize[2]; k++) {
				fprintf(stderr, "%s", clean_string);
				fprintf(stderr, "Receiving %d of %d ( -> pipe2xyza) ...\r", ++count, rcv_planes);

				if (openpipe2d(matrix2d) != 0)
					goto escape2;

				for (i = 0; i < datasize[0]; i++) {
					for (j = 0; j < datasize[1]; j++) {
						matrix3d_[i][k][j] = matrix2d[j][i];
					}
				}
			}

			for (i = 0; i < datasize[0]; i++) {

				if (fp != stdout) {
					sprintf(filename, spectra4d, i + 1, l + 1);

					if ((fp = fopen(filename, "w")) == NULL) {
						fprintf(stderr, "Permission denied to write %s.\n", filename);
						exit(EXIT_FAILURE);
					}
				}

				fprintf(stderr, "%s", clean_string);
				if (fp != stdout)
					fprintf(stderr, "Writing %s <- pipe2xyza ...\r", filename);

				if (l + k == 0 || fp != stdout)
					fwrite(header, sizeof(char), PIPE_HEADERSIZE, fp);

				fpwrite2bin_swap(fp, &(matrix3d_[i][0][0]), dataw_plane, swapdata);

				if (fp != stdout)
					fclose(fp);
				else
					fflush(stdout);
			}
		}

		free_fmatrix3d(matrix3d_);

		break;
	}

	free_fmatrix2d(matrix2d);

	fputc('\n', stderr);

	return 0;

	escape2:free_fmatrix3d(matrix3d_);

	escape:free_fmatrix2d(matrix2d);

	return EXIT_FAILURE;
}
