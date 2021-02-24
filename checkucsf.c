/*
    xyza2pipe - a cross conversion environment of NMR spectra
    Copyright 2017-2021 Masashi Yokochi

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

int checkucsf(char filename[])
{
	struct stat _stat;
	FILE *fp;
	char buffer[MAXCHAR];
	int fidsize[4] = { 0 };
	short apod_code[4] = { SINE_BELL };
	int j, k, _dimension = dimension, shift[4] = { 0 };
	int data_volume;
	long size = 0;
	float apod_par[4][3] = { {0.0} };
	float *vp;
	short *vs;
	unsigned short *uvs;

	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Spectra file %s: Couldn't open.\n", filename);
		return EXIT_FAILURE;
	}

	stat(filename, &_stat);
	size = _stat.st_size;

	fprintf(stderr, "Reading %s ...\n", filename);

	fread(buffer, sizeof(char), size > MAXCHAR ? MAXCHAR : size, fp);

	fclose(fp);

	/* HEADER SIZE */

	dimension = 0;

	j = 134;
	vs = (short *) ((void *) (&buffer[j]));
	headersize = *vs;

	byteswap = 0;

	switch (headersize) {
	case 436:
		dimension = 2;
		break;
	case 564:
		dimension = 3;
		break;
	case 692:
		dimension = 4;
		break;
	default:
		swapbyte(sizeof(short), sizeof(short), buffer + j);
		headersize = *vs;

		switch (headersize) {
		case 436:
			dimension = 2;
			byteswap = 1;
			break;
		case 564:
			dimension = 3;
			byteswap = 1;
			break;
		case 692:
			dimension = 4;
			byteswap = 1;
			break;
		}
	}

	/* FILE SEARCH */
	if (strncmp(buffer, "UCSF", 4) != 0) {
		fprintf(stderr, "Spectra file: Not UCSF file.\n");
		return EXIT_FAILURE;
	}

	if (dimension < 2 || dimension > 4) {

		if (_dimension >= 2 && _dimension <= 4)
			dimension = _dimension;
		else {
			fprintf(stderr, "Spectra file: Not 2D/3D/4D experiment.\n");
			dimension = 2;
		}

		fprintf(stderr, "Try to read %s as %dD spectra. (headersize=%d)\n", filename, dimension, headersize);
		switch (dimension) {
		case 2:
			headersize = 436;
			break;
		case 3:
			headersize = 564;
			break;
		case 4:
			headersize = 692;
			break;
		}
		byteswap = 0;
	}

	/* READ AXIS NAME */
	j = 180;

	for (k = 0; k < dimension; k++) {
		memcpy(axisname[k], buffer + j + 128 * k, MAXAXISNAME);

		if (*axislabel[k] != 0)
			memcpy(axisname[k], axislabel[k], MAXAXISNAME);

		if (usrlabel == 0)
			checklabel(axisname[k]);
	}

	switch (dimension) {
	case 2:
		fprintf(stderr, "Axis Label | %8s  %8s\n", axisname[0], axisname[1]);

		/* READ DATA SIZE */
		for (k = 0; k < dimension; k++) {
			j = 190 + 128 * k;

			if (byteswap != 0)
				swapbyte(sizeof(short), sizeof(short), buffer + j);

			uvs = (unsigned short *) ((void *) (&buffer[j]));
			datasize[k] = (int) *uvs;

			fidsize[k] = datasize[k] / 2;
		}

		memcpy(&(datasize_orig[0]), &(datasize[0]), dimension * sizeof(int));

		fprintf(stderr, "Data Size  | %8d  %8d\n", datasize[0], datasize[1]);

		/* READ BLOCK SIZE */
		for (k = 0; k < dimension; k++) {
			j = 198 + 128 * k;

			if (byteswap != 0)
				swapbyte(sizeof(short), sizeof(short), buffer + j);

			vs = (short *) ((void *) (&buffer[j]));
			blocksize[k] = *vs;
		}

		fprintf(stderr, "Block Size | %8d  %8d\n", blocksize[0], blocksize[1]);

		for (k = 0; k < dimension; k++) {
			if (blocksize[k] > datasize[k] || blocksize[k] > 1024)
				break;
		}

		if (k < dimension) {
			byteswap = (byteswap + 1) % 2;

			/* READ DATA SIZE */
			for (k = 0; k < dimension; k++) {
				j = 190 + 128 * k;

				swapbyte(sizeof(short), sizeof(short), buffer + j);

				uvs = (unsigned short *) ((void *) (&buffer[j]));
				datasize[k] = (int) *uvs;

				fidsize[k] = datasize[k] / 2;
			}

			memcpy(&(datasize_orig[0]), &(datasize[0]), dimension * sizeof(int));

			fprintf(stderr, "Data Size  | %8d  %8d (fixed)\n", datasize[0], datasize[1]);

			/* READ BLOCK SIZE */
			for (k = 0; k < dimension; k++) {
				j = 198 + 128 * k;

				swapbyte(sizeof(short), sizeof(short), buffer + j);

				vs = (short *) ((void *) (&buffer[j]));
				blocksize[k] = *vs;
			}

			fprintf(stderr, "Block Size | %8d  %8d (fixed)\n", blocksize[0], blocksize[1]);
		}

		data_volume = get_data_volume();

		if (byteswap != 0)
			swapbyte(sizeof(float), headersize, buffer);

		/* READ OBS. FREQ. */
		j = 200;
		for (k = 0; k < dimension; k++) {
			vp = (float *) ((void *) (&buffer[j + 128 * k]));
			obsfreq[k] = *vp;
		}

		/* READ SPECTRAL WIDTH */
		j = 204;
		for (k = 0; k < dimension; k++) {
			vp = (float *) ((void *) (&buffer[j + 128 * k]));
			spwidth[k] = *vp;
		}

		/* READ SPECTRAL CENTER */
		j = 208;
		for (k = 0; k < dimension; k++) {
			vp = (float *) ((void *) (&buffer[j + 128 * k]));
			spcenter[k] = *vp;
		}

		if (obsfreq[0] < 1.0e+1 || obsfreq[1] < 1.0e+1 || obsfreq[0] > 1.0e+5 || obsfreq[1] > 1.0e+5 || spcenter[0] < -1.0e+2
				|| spcenter[1] < -1.0e+2 || spcenter[0] > 1.0e+3 || spcenter[1] > 1.0e+3 || spwidth[0] < 1.0e+2
				|| spwidth[1] < 1.0e+2 || spwidth[0] > 1.0e+6 || spwidth[1] > 1.0e+6) {
			swapbyte(sizeof(float), headersize, buffer);

			byteswap = (byteswap + 1) % 2;

			/* READ OBS. FREQ. */
			j = 200;
			for (k = 0; k < dimension; k++) {
				vp = (float *) ((void *) (&buffer[j + 128 * k]));
				obsfreq[k] = *vp;
			}

			/* READ SPECTRAL WIDTH */
			j = 204;
			for (k = 0; k < dimension; k++) {
				vp = (float *) ((void *) (&buffer[j + 128 * k]));
				spwidth[k] = *vp;
			}

			/* READ SPECTRAL CENTER */
			j = 208;
			for (k = 0; k < dimension; k++) {
				vp = (float *) ((void *) (&buffer[j + 128 * k]));
				spcenter[k] = *vp;
			}
		}

		fprintf(stderr, "Obs. Freq. | %8.3f  %8.3f [MHz]\n", obsfreq[0], obsfreq[1]);
		fprintf(stderr, "Spec.Width | %8.2f  %8.2f [Hz]\n", spwidth[0], spwidth[1]);

		/* READ ORIG. FREQ. */

		for (j = 0; j < dimension; j++) {
			origfreq[j] = spcenter[j] * obsfreq[j] - spwidth[j] / 2.0 + spwidth[j] / (float) (datasize[j]);

			if (*axisname[j] == 'H' && strncmp(axisname[j], "HA", MAXAXISNAME) != 0 && *axisname[(j + 1) % dimension] == 'N'
					&& spcenter[j] >= 4.5 && spcenter[j] <= 5.0 && spwidth[j] / obsfreq[j] < 8.5)
				origfreq[j] += spwidth[j] / 2.0;
		}

		fprintf(stderr, "Orig.Freq. | %8.2f  %8.2f [Hz]\n", origfreq[0], origfreq[1]);

		/* CALCULATE SHIFT */

		for (j = 0; j < dimension; j++) {
			if ((int) (fabs(spcenter[j] * obsfreq[j] - origfreq[j]) / (spwidth[j] * 2.0) * (float) (datasize[j])) != 0)
				shift[j] =
						(int) ((spcenter[j] * obsfreq[j] - spwidth[j] / 2.0 * (datasize[j] - 2) / (float) (datasize[j]) -
								origfreq[j]) / (spwidth[j] / (float) (datasize[j])) + 0.5);
		}

		for (j = 0; j < dimension; j++)
			spcenter[j] -= shift[j] * spwidth[j] / (float) (datasize[j]);

		fprintf(stderr, "Spec.Center| %8.2f  %8.2f [ppm]\n", spcenter[0], spcenter[1]);

		if (usrshift != 0) {

			for (j = 0; j < dimension; j++) {

				if (usrcenter[j] == NULLPPM)
					continue;

				origfreq[j] += (usrcenter[j] - spcenter[j]) * obsfreq[j];
				spcenter[j] = usrcenter[j];
			}

			fprintf(stderr, "Spec.Center| %8.2f  %8.2f [ppm] (USER)\n", spcenter[0], spcenter[1]);
			fprintf(stderr, "Orig.Freq. | %8.2f  %8.2f [Hz]  (USER)\n", origfreq[0], origfreq[1]);
		}

		fputc('\n', stderr);

		/* CHECK FILE SIZE */
		if (size != sizeof(float) * data_volume + headersize) {
			fprintf(stderr, "Spectra file %s: Partially broken. (Actual=%d Expected=%d)\n", filename, (int) (size),
					(int) (sizeof(float)) * data_volume + headersize);

			fprintf(stderr, "Try to fix data points from file size...\n");

			for (j = 0; j < dimension; j++) {

				if (datasize[j] % blocksize[j] == 0)
					continue;

				datasize[j] = (size - headersize) / sizeof(float) / datasize[(j + 1) % dimension];
			}

			for (j = 0; j < dimension; j++) {

				if (datasize[j] % blocksize[j] == 0)
					continue;

				memcpy(&(datasize[0]), &(datasize_orig[0]), dimension * sizeof(int));

				data_volume = reset_datasize_by_blocksize();

				if (size != sizeof(float) * data_volume + headersize) {

					fprintf(stderr, "Failed.\n");

					return EXIT_FAILURE;
				}

				memcpy(&(datasize_orig[0]), &(datasize[0]), dimension * sizeof(int));

				break;
			}

			fprintf(stderr, "Data Size  | %8d  %8d (fixed)\n", datasize[0], datasize[1]);
		}
		break;

	case 3:
		fprintf(stderr, "Axis Label | %8s  %8s  %8s\n", axisname[0], axisname[1], axisname[2]);

		/* READ DATA SIZE */
		for (k = 0; k < dimension; k++) {
			j = 190 + 128 * k;

			if (byteswap != 0)
				swapbyte(sizeof(short), sizeof(short), buffer + j);

			uvs = (unsigned short *) ((void *) (&buffer[j]));
			datasize[k] = (int) *uvs;

			fidsize[k] = datasize[k] / 2;
		}

		memcpy(&(datasize_orig[0]), &(datasize[0]), dimension * sizeof(int));

		fprintf(stderr, "Data Size  | %8d  %8d  %8d\n", datasize[0], datasize[1], datasize[2]);

		/* READ BLOCK SIZE */
		for (k = 0; k < dimension; k++) {
			j = 198 + 128 * k;
			vs = (short *) ((void *) (&buffer[j]));

			if (byteswap != 0)
				swapbyte(sizeof(short), sizeof(short), buffer + j);

			blocksize[k] = *vs;
		}

		fprintf(stderr, "Block Size | %8d  %8d  %8d\n", blocksize[0], blocksize[1], blocksize[2]);

		for (k = 0; k < dimension; k++) {
			if (blocksize[k] > datasize[k] || blocksize[k] > 512)
				break;
		}

		if (k < dimension) {
			byteswap = (byteswap + 1) % 2;

			/* READ DATA SIZE */
			for (k = 0; k < dimension; k++) {
				j = 190 + 128 * k;

				swapbyte(sizeof(short), sizeof(short), buffer + j);

				uvs = (unsigned short *) ((void *) (&buffer[j]));
				datasize[k] = (int) *uvs;

				fidsize[k] = datasize[k] / 2;
			}

			memcpy(&(datasize_orig[0]), &(datasize[0]), dimension * sizeof(int));

			fprintf(stderr, "Data Size  | %8d  %8d  %8d (fixed)\n", datasize[0], datasize[1], datasize[2]);

			/* READ BLOCK SIZE */
			for (k = 0; k < dimension; k++) {
				j = 198 + 128 * k;
				vs = (short *) ((void *) (&buffer[j]));

				swapbyte(sizeof(short), sizeof(short), buffer + j);

				blocksize[k] = *vs;
			}

			fprintf(stderr, "Block Size | %8d  %8d  %8d (fixed)\n", blocksize[0], blocksize[1], blocksize[2]);
		}

		data_volume = get_data_volume();

		if (byteswap != 0)
			swapbyte(sizeof(float), headersize, buffer);

		/* READ OBS. FREQ. */
		j = 200;
		for (k = 0; k < dimension; k++) {
			vp = (float *) ((void *) (&buffer[j + 128 * k]));
			obsfreq[k] = *vp;
		}

		/* READ SPECTRAL WIDTH */
		j = 204;
		for (k = 0; k < dimension; k++) {
			vp = (float *) ((void *) (&buffer[j + 128 * k]));
			spwidth[k] = *vp;
		}

		/* READ SPECTRAL CENTER */
		j = 208;
		for (k = 0; k < dimension; k++) {
			vp = (float *) ((void *) (&buffer[j + 128 * k]));
			spcenter[k] = *vp;
		}

		if (obsfreq[0] < 1.0e+1 || obsfreq[1] < 1.0e+1 || obsfreq[2] < 1.0e+1 || obsfreq[0] > 1.0e+5 || obsfreq[1] > 1.0e+5
				|| obsfreq[2] > 1.0e+5 || spcenter[0] < -1.0e+2 || spcenter[1] < -1.0e+2 || spcenter[2] < -1.0e+2
				|| spcenter[0] > 1.0e+3 || spcenter[1] > 1.0e+3 || spcenter[2] > 1.0e+3 || spwidth[0] < 1.0e+2
				|| spwidth[1] < 1.0e+2 || spwidth[2] < 1.0e+2 || spwidth[0] > 1.0e+6 || spwidth[1] > 1.0e+6
				|| spwidth[2] > 1.0e+6) {
			swapbyte(sizeof(float), headersize, buffer);

			byteswap = (byteswap + 1) % 2;

			/* READ OBS. FREQ. */
			j = 200;
			for (k = 0; k < dimension; k++) {
				vp = (float *) ((void *) (&buffer[j + 128 * k]));
				obsfreq[k] = *vp;
			}

			/* READ SPECTRAL WIDTH */
			j = 204;
			for (k = 0; k < dimension; k++) {
				vp = (float *) ((void *) (&buffer[j + 128 * k]));
				spwidth[k] = *vp;
			}

			/* READ SPECTRAL CENTER */
			j = 208;
			for (k = 0; k < dimension; k++) {
				vp = (float *) ((void *) (&buffer[j + 128 * k]));
				spcenter[k] = *vp;
			}
		}

		fprintf(stderr, "Obs. Freq. | %8.3f  %8.3f  %8.3f [MHz]\n", obsfreq[0], obsfreq[1], obsfreq[2]);
		fprintf(stderr, "Spec.Width | %8.2f  %8.2f  %8.2f [Hz]\n", spwidth[0], spwidth[1], spwidth[2]);

		/* READ ORIG. FREQ. */

		for (j = 0; j < dimension; j++) {
			origfreq[j] = spcenter[j] * obsfreq[j] - spwidth[j] / 2.0 + spwidth[j] / (float) (datasize[j]);

			if (*axisname[j] == 'H' && strncmp(axisname[j], "HA", MAXAXISNAME) != 0 && (*axisname[(j + 1) % dimension] == 'N'
					|| *axisname[(j + 2) % dimension] == 'N') && spcenter[j] >= 4.5 && spcenter[j] <= 5.0
					&& spwidth[j] / obsfreq[j] < 8.5)
				origfreq[j] += spwidth[j] / 2.0;
		}

		fprintf(stderr, "Orig.Freq. | %8.2f  %8.2f  %8.2f [Hz]\n", origfreq[0], origfreq[1], origfreq[2]);

		/* CALCULATE SHIFT */

		for (j = 0; j < dimension; j++) {
			if ((int) (fabs(spcenter[j] * obsfreq[j] - origfreq[j]) / (spwidth[j] * 2.0) * (float) (datasize[j])) != 0)
				shift[j] =
						(int) ((spcenter[j] * obsfreq[j] - spwidth[j] / 2.0 * (datasize[j] - 2) / (float) (datasize[j]) -
								origfreq[j]) / (spwidth[j] / (float) (datasize[j])) + 0.5);
		}

		for (j = 0; j < dimension; j++)
			spcenter[j] -= shift[j] * spwidth[j] / (float) (datasize[j]);

		fprintf(stderr, "Spec.Center| %8.2f  %8.2f  %8.2f [ppm]\n", spcenter[0], spcenter[1], spcenter[2]);

		if (usrshift != 0) {

			for (j = 0; j < dimension; j++) {

				if (usrcenter[j] == NULLPPM)
					continue;

				origfreq[j] += (usrcenter[j] - spcenter[j]) * obsfreq[j];
				spcenter[j] = usrcenter[j];
			}

			fprintf(stderr, "Spec.Center| %8.2f  %8.2f  %8.2f [ppm] (USER)\n", spcenter[0], spcenter[1], spcenter[2]);
			fprintf(stderr, "Orig.Freq. | %8.2f  %8.2f  %8.2f [Hz]  (USER)\n", origfreq[0], origfreq[1], origfreq[2]);
		}

		fputc('\n', stderr);

		/* CHECK FILE SIZE */
		if (size != sizeof(float) * data_volume + headersize) {
			fprintf(stderr, "Spectra file %s: Partially broken. (Actual=%d Expected=%d)\n", filename, (int) (size),
					(int) (sizeof(float)) * data_volume + headersize);

			fprintf(stderr, "Try to fix data points from file size...\n");

			for (j = 0; j < dimension; j++) {

				if (datasize[j] % blocksize[j] == 0)
					continue;

				datasize[j] = (size - headersize) / sizeof(float) / datasize[(j + 1) % dimension] / datasize[(j + 2) % dimension];
			}

			for (j = 0; j < dimension; j++) {

				if (datasize[j] % blocksize[j] == 0)
					continue;

				memcpy(&(datasize[0]), &(datasize_orig[0]), dimension * sizeof(int));

				data_volume = reset_datasize_by_blocksize();

				if (size != sizeof(float) * data_volume + headersize) {

					fprintf(stderr, "Failed.\n");

					return EXIT_FAILURE;
				}

				memcpy(&(datasize_orig[0]), &(datasize[0]), dimension * sizeof(int));

				break;
			}

			fprintf(stderr, "Data Size  | %8d  %8d  %8d (fixed)\n", datasize[0], datasize[1], datasize[2]);
		}
		break;

	case 4:
		fprintf(stderr, "Axis Label | %8s  %8s  %8s  %8s\n", axisname[0], axisname[1], axisname[2], axisname[3]);

		/* READ DATA SIZE */
		for (k = 0; k < dimension; k++) {
			j = 190 + 128 * k;

			if (byteswap != 0)
				swapbyte(sizeof(short), sizeof(short), buffer + j);

			uvs = (unsigned short *) ((void *) (&buffer[j]));
			datasize[k] = (int) *uvs;

			fidsize[k] = datasize[k] / 2;
		}

		memcpy(&(datasize_orig[0]), &(datasize[0]), dimension * sizeof(int));

		fprintf(stderr, "Data Size  | %8d  %8d  %8d  %8d\n", datasize[0], datasize[1], datasize[2], datasize[3]);

		/* READ BLOCK SIZE */
		for (k = 0; k < dimension; k++) {
			j = 198 + 128 * k;

			if (byteswap != 0)
				swapbyte(sizeof(short), sizeof(short), buffer + j);

			vs = (short *) ((void *) (&buffer[j]));
			blocksize[k] = *vs;
		}

		fprintf(stderr, "Block Size | %8d  %8d  %8d  %8d\n", blocksize[0], blocksize[1], blocksize[2], blocksize[3]);

		for (k = 0; k < dimension; k++) {
			if (blocksize[k] > datasize[k] || blocksize[k] > 256)
				break;
		}

		if (k < dimension) {
			byteswap = (byteswap + 1) % 2;

			/* READ DATA SIZE */
			for (k = 0; k < dimension; k++) {
				j = 190 + 128 * k;

				swapbyte(sizeof(short), sizeof(short), buffer + j);

				uvs = (unsigned short *) ((void *) (&buffer[j]));
				datasize[k] = (int) *uvs;

				fidsize[k] = datasize[k] / 2;
			}

			memcpy(&(datasize_orig[0]), &(datasize[0]), dimension * sizeof(int));

			fprintf(stderr, "Data Size  | %8d  %8d  %8d  %8d (fixed)\n", datasize[0], datasize[1], datasize[2], datasize[3]);

			/* READ BLOCK SIZE */
			for (k = 0; k < dimension; k++) {
				j = 198 + 128 * k;

				swapbyte(sizeof(short), sizeof(short), buffer + j);

				vs = (short *) ((void *) (&buffer[j]));
				blocksize[k] = *vs;
			}

			fprintf(stderr, "Block Size | %8d  %8d  %8d  %8d (fixed)\n", blocksize[0], blocksize[1], blocksize[2], blocksize[3]);
		}

		data_volume = get_data_volume();

		if (byteswap != 0)
			swapbyte(sizeof(float), headersize, buffer);

		/* READ OBS. FREQ. */
		j = 200;
		for (k = 0; k < dimension; k++) {
			vp = (float *) ((void *) (&buffer[j + 128 * k]));
			obsfreq[k] = *vp;
		}

		/* READ SPECTRAL WIDTH */
		j = 204;
		for (k = 0; k < dimension; k++) {
			vp = (float *) ((void *) (&buffer[j + 128 * k]));
			spwidth[k] = *vp;
		}

		/* READ SPECTRAL CENTER */
		j = 208;
		for (k = 0; k < dimension; k++) {
			vp = (float *) ((void *) (&buffer[j + 128 * k]));
			spcenter[k] = *vp;
		}

		if (obsfreq[0] < 1.0e+1 || obsfreq[1] < 1.0e+1 || obsfreq[2] < 1.0e+1 || obsfreq[3] < 1.0e+1 || obsfreq[0] > 1.0e+5
				|| obsfreq[1] > 1.0e+5 || obsfreq[2] > 1.0e+5 || obsfreq[3] > 1.0e+5 || spcenter[0] < -1.0e+2
				|| spcenter[1] < -1.0e+2 || spcenter[2] < -1.0e+2 || spcenter[3] < -1.0e+2 || spcenter[0] > 1.0e+3
				|| spcenter[1] > 1.0e+3 || spcenter[2] > 1.0e+3 || spcenter[3] > 1.0e+3 || spwidth[0] < 1.0e+2
				|| spwidth[1] < 1.0e+2 || spwidth[2] < 1.0e+2 || spwidth[3] < 1.0e+2 || spwidth[0] > 1.0e+6 || spwidth[1] > 1.0e+6
				|| spwidth[2] > 1.0e+6 || spwidth[3] > 1.0e+6) {
			swapbyte(sizeof(float), headersize, buffer);

			byteswap = (byteswap + 1) % 2;

			/* READ OBS. FREQ. */
			j = 200;
			for (k = 0; k < dimension; k++) {
				vp = (float *) ((void *) (&buffer[j + 128 * k]));
				obsfreq[k] = *vp;
			}

			/* READ SPECTRAL WIDTH */
			j = 204;
			for (k = 0; k < dimension; k++) {
				vp = (float *) ((void *) (&buffer[j + 128 * k]));
				spwidth[k] = *vp;
			}

			/* READ SPECTRAL CENTER */
			j = 208;
			for (k = 0; k < dimension; k++) {
				vp = (float *) ((void *) (&buffer[j + 128 * k]));
				spcenter[k] = *vp;
			}
		}

		fprintf(stderr, "Obs. Freq. | %8.3f  %8.3f  %8.3f  %8.3f [MHz]\n", obsfreq[0], obsfreq[1], obsfreq[2], obsfreq[3]);
		fprintf(stderr, "Spec.Width | %8.2f  %8.2f  %8.2f  %8.2f [Hz]\n", spwidth[0], spwidth[1], spwidth[2], spwidth[3]);

		/* READ ORIG. FREQ. */

		for (j = 0; j < dimension; j++) {
			origfreq[j] = spcenter[j] * obsfreq[j] - spwidth[j] / 2.0 + spwidth[j] / (float) (datasize[j]);

			if (*axisname[j] == 'H' && strncmp(axisname[j], "HA", MAXAXISNAME) != 0 && (*axisname[(j + 1) % dimension] == 'N'
					|| *axisname[(j + 2) % dimension] == 'N' || *axisname[(j + 3) % dimension] == 'N') && spcenter[j] >= 4.5
					&& spcenter[j] <= 5.0 && spwidth[j] / obsfreq[j] < 8.5)
				origfreq[j] += spwidth[j] / 2.0;
		}

		fprintf(stderr, "Orig.Freq. | %8.2f  %8.2f  %8.2f  %8.2f [Hz]\n", origfreq[0], origfreq[1], origfreq[2], origfreq[3]);

		/* CALCULATE SHIFT */

		for (j = 0; j < dimension; j++) {
			if ((int) (fabs(spcenter[j] * obsfreq[j] - origfreq[j]) / (spwidth[j] * 2.0) * (float) (datasize[j])) != 0)
				shift[j] =
						(int) ((spcenter[j] * obsfreq[j] - spwidth[j] / 2.0 * (datasize[j] - 2) / (float) (datasize[j]) -
								origfreq[j]) / (spwidth[j] / (float) (datasize[j])) + 0.5);
		}

		for (j = 0; j < dimension; j++)
			spcenter[j] -= shift[j] * spwidth[j] / (float) (datasize[j]);

		fprintf(stderr, "Spec.Center| %8.2f  %8.2f  %8.2f  %8.2f [ppm]\n", spcenter[0], spcenter[1], spcenter[2],
				spcenter[3]);

		if (usrshift != 0) {

			for (j = 0; j < dimension; j++) {

				if (usrcenter[j] == NULLPPM)
					continue;

				origfreq[j] += (usrcenter[j] - spcenter[j]) * obsfreq[j];
				spcenter[j] = usrcenter[j];
			}

			fprintf(stderr, "Spec.Center| %8.2f  %8.2f  %8.2f  %8.2f [ppm] (USER)\n", spcenter[0], spcenter[1], spcenter[2],
					spcenter[3]);
			fprintf(stderr, "Orig.Freq. | %8.2f  %8.2f  %8.2f  %8.2f [Hz]  (USER)\n", origfreq[0], origfreq[1], origfreq[2],
					origfreq[3]);
		}

		fputc('\n', stderr);

		/* CHECK FILE SIZE */
		if (size != sizeof(float) * data_volume + headersize) {
			fprintf(stderr, "Spectra file %s: Partially broken. (Actual=%d Expected=%d)\n", filename, (int) (size),
					(int) (sizeof(float)) * data_volume + headersize);

			fprintf(stderr, "Try to fix data points from file size...\n");

			for (j = 0; j < dimension; j++) {

				if (datasize[j] % blocksize[j] == 0)
					continue;

				datasize[j] =
						(size - headersize) / sizeof(float) / datasize[(j + 1) % dimension] / datasize[(j +
								2) % dimension] / datasize[(j + 3) % dimension];
			}

			for (j = 0; j < dimension; j++) {

				if (datasize[j] % blocksize[j] == 0)
					continue;

				memcpy(&(datasize[0]), &(datasize_orig[0]), dimension * sizeof(int));

				data_volume = reset_datasize_by_blocksize();

				if (size != sizeof(float) * data_volume + headersize) {

					fprintf(stderr, "Failed.\n");

					return EXIT_FAILURE;
				}

				memcpy(&(datasize_orig[0]), &(datasize[0]), dimension * sizeof(int));

				break;
			}

			fprintf(stderr, "Data Size  | %8d  %8d  %8d  %8d (fixed)\n", datasize[0], datasize[1], datasize[2], datasize[3]);
		}
		break;
	}

	/* APODIZATION (assumed) */

	for (j = 0; j < dimension; j++) {
		apod_code[j] = SINE_BELL;
		apod_par[j][0] = 0.5;
		apod_par[j][1] = 0.98;
		apod_par[j][2] = (j == 0 ? 2.0 : 1.0);
	}

	memset(&(header[0]), 0, PIPE_HEADERSIZE * sizeof(char));

	fwrite2mem(header, PIPE_HEADER[0]);
	fwrite2mem(header + 4, PIPE_HEADER[1]);
	fwrite2mem(header + 8, PIPE_HEADER[2]);

	fwrite2mem(header + 36, (float) (dimension));

	fwrite2mem(header + 880, 1.0);
	fwrite2mem(header + 888, 1.0);

	if (dimension >= 3)
		fwrite2mem(header + 52, 1.0);
	if (dimension >= 4)
		fwrite2mem(header + 124, 1.0);

	fwrite2mem(header + 220, 1.0);
	fwrite2mem(header + 224, 1.0);

	if (dimension >= 3)
		fwrite2mem(header + 204, 1.0);
	if (dimension >= 4)
		fwrite2mem(header + 216, 1.0);

	fwrite2mem(header + 884, 0.0);

	fwrite2mem(header + 96, 2.0);
	fwrite2mem(header + 100, 1.0);
	fwrite2mem(header + 104, 3.0);
	fwrite2mem(header + 108, 4.0);

	j = 64;
	for (k = 0; k < dimension; k++) {
		memcpy(header + j + k * MAXAXISNAME, axisname[k], MAXAXISNAME);
		unitsize[k] = datasize[k] / blocksize[k];
	}

	switch (dimension) {
	case 4:
		fwrite2mem(header + 112, obsfreq[3]);

		fwrite2mem(header + 276, spcenter[3]);

		fwrite2mem(header + 120, origfreq[3]);

		fwrite2mem(header + 116, spwidth[3]);

		fwrite2mem(header + 128, (float) (datasize[3]));

		fwrite2mem(header + 1556, (float) (fidsize[3]));

		fwrite2mem(header + 212, (float) (fidsize[3]));

		fwrite2mem(header + 1620, (float) (apod_code[3]));

		for (j = 0; j < 3; j++)
			fwrite2mem(header + 1624 + j * sizeof(float), apod_par[3][j]);

		fwrite2mem(header + 248, phase[3][0]);
		fwrite2mem(header + 252, phase[3][1]);
		/* no break */

	case 3:
		fwrite2mem(header + 40, obsfreq[2]);

		fwrite2mem(header + 272, spcenter[2]);

		fwrite2mem(header + 48, origfreq[2]);

		fwrite2mem(header + 44, spwidth[2]);

		fwrite2mem(header + 60, (float) (datasize[2]));

		fwrite2mem(header + 1552, (float) (fidsize[2]));

		fwrite2mem(header + 200, (float) (fidsize[2]));

		fwrite2mem(header + 1600, (float) (apod_code[2]));

		for (j = 0; j < 3; j++)
			fwrite2mem(header + 1604 + j * sizeof(float), apod_par[2][j]);

		fwrite2mem(header + 240, phase[2][0]);
		fwrite2mem(header + 244, phase[2][1]);
		/* no break */

	case 2:
		fwrite2mem(header + 476, obsfreq[0]);
		fwrite2mem(header + 872, obsfreq[1]);

		fwrite2mem(header + 264, spcenter[0]);
		fwrite2mem(header + 268, spcenter[1]);

		fwrite2mem(header + 404, origfreq[0]);
		fwrite2mem(header + 996, origfreq[1]);

		fwrite2mem(header + 400, spwidth[0]);
		fwrite2mem(header + 916, spwidth[1]);

		fwrite2mem(header + 396, (float) (datasize[0]));
		fwrite2mem(header + 876, (float) (datasize[1]));

		fwrite2mem(header + 1544, (float) (fidsize[0]));
		fwrite2mem(header + 1548, (float) (fidsize[1]));

		fwrite2mem(header + 380, (float) (fidsize[0]));
		fwrite2mem(header + 1712, (float) (fidsize[1]));

		fwrite2mem(header + 1652, (float) (apod_code[0]));
		fwrite2mem(header + 1656, (float) (apod_code[1]));

		for (j = 0; j < 3; j++) {
			fwrite2mem(header + 1660 + j * sizeof(float), apod_par[0][j]);
			fwrite2mem(header + 1680 + j * sizeof(float), apod_par[1][j]);
		}

		fwrite2mem(header + 436, phase[0][0]);
		fwrite2mem(header + 440, phase[0][1]);
		fwrite2mem(header + 980, phase[1][0]);
		fwrite2mem(header + 984, phase[1][1]);
		break;
	}

	/* INDIRECT PLANES */
	fwrite2mem(header + 1768, (float) get_indirect_planes());

	/* QUADFLAG */
	fwrite2mem(header + 424, 1.0);

	/* STATES */
	fwrite2mem(header + 1024, 2.0);

	fwrite2mem(header + 1028, 1.0);
	fwrite2mem(header + 1036, 1.0);
	fwrite2mem(header + 1044, 1.0);
	fwrite2mem(header + 1052, 1.0);

	return 0;
}
