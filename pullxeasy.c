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

extern void float2xeasy(float x, unsigned char *lo16);

int pullxeasy2d(char spectra2d[])
{
	FILE *fp;
	char parfile[MAXLONGNAME];
	unsigned char x16[2];
	int i, j, count = 0;
	int block_volume, block_i, block_j, block_id;
	int offset_i, offset_j, offset;
	float **matrix2d;

	make_dir(spectra2d);

	set_clean_string(clean_string);

	if (strstr(spectra2d, ".16"))
		strreplacecpy(parfile, spectra2d, ".16", ".param");
	else
		sprintf(parfile, "%s.param", spectra2d);

	if ((fp = fopen(parfile, "w")) == NULL) {
		fprintf(stderr, "Permission denied to write %s.\n", parfile);
		exit(EXIT_FAILURE);
	}

	block_volume = set_block_volume();

	fprintf(fp, "Version ....................... 1\n");
	fprintf(fp, "Number of dimensions .......... 2\n");
	fprintf(fp, "16 or 8 bit file type ......... 16\n");
	for (j = 0; j < dimension; j++)
		fprintf(fp, "Spectrometer frequency in w%d .. %f\n", j + 1, obsfreq[j]);
	for (j = 0; j < dimension; j++)
		fprintf(fp, "Spectral sweep width in w%d .... %f\n", j + 1, spwidth[j] / obsfreq[j]);
	for (j = 0; j < dimension; j++)
		fprintf(fp, "Maximum chemical shift in w%d .. %f\n", j + 1,
				(relyof ==
						0 ? spcenter[j] : (origfreq[j] + spwidth[j] / 2.0 - spwidth[j] / (float) (datasize[j])) / obsfreq[j]) + (j ==
								0 && leftcar ? leftcar * spwidth[0] / 2.0 / obsfreq[0] : 0.0) + spwidth[j] / obsfreq[j] / 2.0);
	for (j = 0; j < dimension; j++)
		fprintf(fp, "Size of spectrum in w%d ........ %d\n", j + 1, datasize[j]);
	for (j = 0; j < dimension; j++)
		fprintf(fp, "Submatrix size in w%d .......... %d\n", j + 1, blocksize[j]);
	for (j = 0; j < dimension; j++)
		fprintf(fp, "Permutation for w%d ............ %d\n", j + 1, j + 1);
	for (j = 0; j < dimension; j++)
		fprintf(fp, "Folding in w%d ................. NO\n", j + 1);
	fprintf(fp, "Type of spectrum .............. -\n");
	for (j = 0; j < dimension; j++)
		fprintf(fp, "Identifier for dimension w%d ... %s\n", j + 1, axisname[j]);

	fclose(fp);

	if ((fp = fopen(spectra2d, "w")) == NULL) {
		fprintf(stderr, "Permission denied to write %s.\n", spectra2d);
		exit(EXIT_FAILURE);
	}

	fprintf(stderr, "%s", clean_string);
	fprintf(stderr, "Writing %s <- pipe2xeasy ...\r", spectra2d);

	matrix2d = fmalloc2d(datasize[1], datasize[0]);

	fprintf(stderr, "%s", clean_string);
	fprintf(stderr, "Receiving %d of %d ( -> pipe2xeasy) ...\r", ++count, 1);

	if (openpipe2d(matrix2d) != 0)
		goto escape;

	for (j = 0; j < datasize[1]; j++) {
		block_j = (int) (j / blocksize[1]);
		offset_j = j - block_j * blocksize[1];

		for (i = 0; i < datasize[0]; i++) {
			block_i = (int) (i / blocksize[0]);
			offset_i = i - block_i * blocksize[0];

			block_id = block_i + block_j * unitsize[0];

			offset = offset_i + offset_j * blocksize[0];

			fseek(fp, (long) ((offset + block_id * block_volume) * sizeof(short)), SEEK_SET);

			float2xeasy(matrix2d[j][i], x16);

			spwrite2bin_swap(fp, (int16_t *) (x16), 1, swapdata);
		}
	}

	free_fmatrix2d(matrix2d);

	fclose(fp);

	fputc('\n', stderr);

	return 0;

	escape:free_fmatrix2d(matrix2d);

	fclose(fp);

	return EXIT_FAILURE;
}

int pullxeasy3d(char spectra3d[])
{
	FILE *fp;
	char parfile[MAXLONGNAME];
	unsigned char x16[2];
	int i, j, k, count = 0;
	int block_volume, block_i, block_j, block_k, block_id;
	int offset_i, offset_j, offset_k, offset;
	int rcv_planes;
	float **matrix2d;

	make_dir(spectra3d);

	set_clean_string(clean_string);

	if (strstr(spectra3d, ".16"))
		strreplacecpy(parfile, spectra3d, ".16", ".param");
	else
		sprintf(parfile, "%s.param", spectra3d);

	if ((fp = fopen(parfile, "w")) == NULL) {
		fprintf(stderr, "Permission denied to write %s.\n", parfile);
		exit(EXIT_FAILURE);
	}

	block_volume = set_block_volume();

	fprintf(fp, "Version ....................... 1\n");
	fprintf(fp, "Number of dimensions .......... 3\n");
	fprintf(fp, "16 or 8 bit file type ......... 16\n");
	for (j = 0; j < dimension; j++)
		fprintf(fp, "Spectrometer frequency in w%d .. %f\n", j + 1, obsfreq[j]);
	for (j = 0; j < dimension; j++)
		fprintf(fp, "Spectral sweep width in w%d .... %f\n", j + 1, spwidth[j] / obsfreq[j]);
	for (j = 0; j < dimension; j++)
		fprintf(fp, "Maximum chemical shift in w%d .. %f\n", j + 1,
				(relyof ==
						0 ? spcenter[j] : (origfreq[j] + spwidth[j] / 2.0 - spwidth[j] / (float) (datasize[j])) / obsfreq[j]) + (j ==
								0 && leftcar ? leftcar * spwidth[0] / 2.0 / obsfreq[0] : 0.0) + spwidth[j] / obsfreq[j] / 2.0);
	for (j = 0; j < dimension; j++)
		fprintf(fp, "Size of spectrum in w%d ........ %d\n", j + 1, datasize[j]);
	for (j = 0; j < dimension; j++)
		fprintf(fp, "Submatrix size in w%d .......... %d\n", j + 1, blocksize[j]);
	for (j = 0; j < dimension; j++)
		fprintf(fp, "Permutation for w%d ............ %d\n", j + 1, j + 1);
	for (j = 0; j < dimension; j++)
		fprintf(fp, "Folding in w%d ................. NO\n", j + 1);
	fprintf(fp, "Type of spectrum .............. -\n");
	for (j = 0; j < dimension; j++)
		fprintf(fp, "Identifier for dimension w%d ... %s\n", j + 1, axisname[j]);

	fclose(fp);

	if ((fp = fopen(spectra3d, "w")) == NULL) {
		fprintf(stderr, "Permission denied to write %s.\n", spectra3d);
		exit(EXIT_FAILURE);
	}

	fprintf(stderr, "%s", clean_string);
	fprintf(stderr, "Writing %s <- pipe2xeasy ...\r", spectra3d);

	matrix2d = fmalloc2d(datasize[1], datasize[0]);

	rcv_planes = datasize[2];

	for (k = 0; k < datasize[2]; k++) {
		fprintf(stderr, "%s", clean_string);
		fprintf(stderr, "Receiving %d of %d ( -> pipe2xeasy) ...\r", ++count, rcv_planes);

		if (openpipe2d(matrix2d) != 0)
			goto escape;

		block_k = (int) (k / blocksize[2]);
		offset_k = k - block_k * blocksize[2];

		for (j = 0; j < datasize[1]; j++) {
			block_j = (int) (j / blocksize[1]);
			offset_j = j - block_j * blocksize[1];

			for (i = 0; i < datasize[0]; i++) {
				block_i = (int) (i / blocksize[0]);
				offset_i = i - block_i * blocksize[0];

				block_id = block_i + (block_j + block_k * unitsize[1]) * unitsize[0];

				offset = offset_i + (offset_j + offset_k * blocksize[1]) * blocksize[0];

				fseek(fp, (long) ((offset + block_id * block_volume) * sizeof(short)), SEEK_SET);

				float2xeasy(matrix2d[j][i], x16);

				spwrite2bin_swap(fp, (int16_t *) (x16), 1, swapdata);
			}
		}
	}

	free_fmatrix2d(matrix2d);

	fclose(fp);

	fputc('\n', stderr);

	return 0;

	escape:free_fmatrix2d(matrix2d);

	fclose(fp);

	return EXIT_FAILURE;
}

int pullxeasy4d(char spectra4d[])
{
	FILE *fp;
	char parfile[MAXLONGNAME];
	unsigned char x16[2];
	int i, j, k, l, count = 0;
	int block_volume, block_i, block_j, block_k, block_l, block_id;
	int offset_i, offset_j, offset_k, offset_l, offset;
	int rcv_planes;
	float **matrix2d;

	make_dir(spectra4d);

	set_clean_string(clean_string);

	if (strstr(spectra4d, ".16"))
		strreplacecpy(parfile, spectra4d, ".16", ".param");
	else
		sprintf(parfile, "%s.param", spectra4d);

	if ((fp = fopen(parfile, "w")) == NULL) {
		fprintf(stderr, "Permission denied to write %s.\n", parfile);
		exit(EXIT_FAILURE);
	}

	block_volume = set_block_volume();

	fprintf(fp, "Version ....................... 1\n");
	fprintf(fp, "Number of dimensions .......... 4\n");
	fprintf(fp, "16 or 8 bit file type ......... 16\n");
	for (j = 0; j < dimension; j++)
		fprintf(fp, "Spectrometer frequency in w%d .. %f\n", j + 1, obsfreq[j]);
	for (j = 0; j < dimension; j++)
		fprintf(fp, "Spectral sweep width in w%d .... %f\n", j + 1, spwidth[j] / obsfreq[j]);
	for (j = 0; j < dimension; j++)
		fprintf(fp, "Maximum chemical shift in w%d .. %f\n", j + 1,
				(relyof ==
						0 ? spcenter[j] : (origfreq[j] + spwidth[j] / 2.0 - spwidth[j] / (float) (datasize[j])) / obsfreq[j]) + (j ==
								0 && leftcar ? leftcar * spwidth[0] / 2.0 / obsfreq[0] : 0.0) + spwidth[j] / obsfreq[j] / 2.0);
	for (j = 0; j < dimension; j++)
		fprintf(fp, "Size of spectrum in w%d ........ %d\n", j + 1, datasize[j]);
	for (j = 0; j < dimension; j++)
		fprintf(fp, "Submatrix size in w%d .......... %d\n", j + 1, blocksize[j]);
	for (j = 0; j < dimension; j++)
		fprintf(fp, "Permutation for w%d ............ %d\n", j + 1, j + 1);
	for (j = 0; j < dimension; j++)
		fprintf(fp, "Folding in w%d ................. NO\n", j + 1);
	fprintf(fp, "Type of spectrum .............. -\n");
	for (j = 0; j < dimension; j++)
		fprintf(fp, "Identifier for dimension w%d ... %s\n", j + 1, axisname[j]);

	fclose(fp);

	if ((fp = fopen(spectra4d, "w")) == NULL) {
		fprintf(stderr, "Permission denied to write %s.\n", spectra4d);
		exit(EXIT_FAILURE);
	}

	fprintf(stderr, "%s", clean_string);
	fprintf(stderr, "Writing %s <- pipe2xeasy ...\r", spectra4d);

	matrix2d = fmalloc2d(datasize[1], datasize[0]);

	rcv_planes = datasize[2] * datasize[3];

	for (l = 0; l < datasize[3]; l++) {
		block_l = (int) (l / blocksize[3]);
		offset_l = l - block_l * blocksize[3];

		for (k = 0; k < datasize[2]; k++) {
			fprintf(stderr, "%s", clean_string);
			fprintf(stderr, "Receiving %d of %d ( -> pipe2xeasy) ...\r", ++count, rcv_planes);

			if (openpipe2d(matrix2d) != 0)
				goto escape;

			block_k = (int) (k / blocksize[2]);
			offset_k = k - block_k * blocksize[2];

			for (j = 0; j < datasize[1]; j++) {
				block_j = (int) (j / blocksize[1]);
				offset_j = j - block_j * blocksize[1];

				for (i = 0; i < datasize[0]; i++) {
					block_i = (int) (i / blocksize[0]);
					offset_i = i - block_i * blocksize[0];

					block_id = block_i + (block_j + (block_k + block_l * unitsize[2]) * unitsize[1]) * unitsize[0];

					offset = offset_i + (offset_j + (offset_k + offset_l * blocksize[2]) * blocksize[1]) * blocksize[0];

					fseek(fp, (long) ((offset + block_id * block_volume) * sizeof(short)), SEEK_SET);

					float2xeasy(matrix2d[j][i], x16);

					spwrite2bin_swap(fp, (int16_t *) (x16), 1, swapdata);
				}
			}
		}
	}

	free_fmatrix2d(matrix2d);

	fclose(fp);

	fputc('\n', stderr);

	return 0;

	escape:free_fmatrix2d(matrix2d);

	fclose(fp);

	return EXIT_FAILURE;
}
