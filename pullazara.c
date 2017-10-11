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

int pullazara2d(char spectra2d[])
{
	FILE *fp;
	char parfile[MAXLONGNAME];
	int i, j, count = 0;
	int block_volume, block_i, block_j, block_id;
	int offset_i, offset_j, offset;
	short swapped = 0;
	float **matrix2d;

	if ((is_big_endian() && is_little_endian_float(header + sizeof(float) * 2, PIPE_HEADER[2])) || (is_little_endian()
			&& is_big_endian_float(header + sizeof(float) * 2, PIPE_HEADER[2])))
		swapped = 1;

	swapdata += swapped;

	make_dir(spectra2d);

	set_clean_string(clean_string);

	sprintf(parfile, "%s.par", spectra2d);

	if ((fp = fopen(parfile, "w")) == NULL) {
		fprintf(stderr, "Permission denied to write %s.\n", parfile);
		exit(EXIT_FAILURE);
	}

	fprintf(fp, "! this file = %s\n\n", parfile);

	fprintf(fp, "ndim 2\n");
	fprintf(fp, "file %s\n", spectra2d);

	block_volume = set_block_volume();

	for (j = 0; j < dimension; j++) {
		fprintf(fp, "\ndim %d\n", j + 1);
		fprintf(fp, "npts %d\n", datasize[j]);
		fprintf(fp, "block %d\n", blocksize[j]);
		fprintf(fp, "sw %f\n", spwidth[j]);
		fprintf(fp, "sf %f\n", obsfreq[j]);
		fprintf(fp, "refppm %f\n", (relyof == 0 ? spcenter[j] : (origfreq[j] + spwidth[j] / 2.0) / obsfreq[j]) + (j == 0
				&& leftcar ? leftcar * spwidth[0] / 2.0 / obsfreq[0] : 0.0));
		fprintf(fp, "refpt %f\n", datasize[j] + 0.5);
		fprintf(fp, "nuc %s\n", axisname[j]);
	}

	fclose(fp);

	if ((fp = fopen(spectra2d, "w")) == NULL) {
		fprintf(stderr, "Permission denied to write %s.\n", spectra2d);
		exit(EXIT_FAILURE);
	}

	fprintf(stderr, "%s", clean_string);
	fprintf(stderr, "Writing %s <- pipe2azara ...\r", spectra2d);

	matrix2d = fmalloc2d(datasize[1], datasize[0]);

	fprintf(stderr, "%s", clean_string);
	fprintf(stderr, "Receiving %d of %d ( -> pipe2azara) ...\r", ++count, 1);

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

			fseek(fp, (long) (headersize + (offset + block_id * block_volume) * sizeof(float)), SEEK_SET);

			fpwrite2bin_swap(fp, &(matrix2d[j][i]), 1, swapdata);
		}
	}

	free_fmatrix2d(matrix2d);

	fclose(fp);

	fputc('\n', stderr);

	return 0;

	escape:free_fmatrix2d(matrix2d);

	fclose(fp);

	return 1;
}

int pullazara3d(char spectra3d[])
{
	FILE *fp;
	char parfile[MAXLONGNAME];
	int i, j, k, count = 0;
	int block_volume, block_i, block_j, block_k, block_id;
	int offset_i, offset_j, offset_k, offset;
	short swapped = 0;
	float **matrix2d;

	if ((is_big_endian() && is_little_endian_float(header + sizeof(float) * 2, PIPE_HEADER[2])) || (is_little_endian()
			&& is_big_endian_float(header + sizeof(float) * 2, PIPE_HEADER[2])))
		swapped = 1;

	swapdata += swapped;

	make_dir(spectra3d);

	set_clean_string(clean_string);

	sprintf(parfile, "%s.par", spectra3d);

	if ((fp = fopen(parfile, "w")) == NULL) {
		fprintf(stderr, "Permission denied to write %s.\n", parfile);
		exit(EXIT_FAILURE);
	}

	fprintf(fp, "! this file = %s\n\n", parfile);

	fprintf(fp, "ndim 3\n");
	fprintf(fp, "file %s\n", spectra3d);

	block_volume = set_block_volume();

	for (j = 0; j < dimension; j++) {
		fprintf(fp, "\ndim %d\n", j + 1);
		fprintf(fp, "npts %d\n", datasize[j]);
		fprintf(fp, "block %d\n", blocksize[j]);
		fprintf(fp, "sw %f\n", spwidth[j]);
		fprintf(fp, "sf %f\n", obsfreq[j]);
		fprintf(fp, "refppm %f\n", (relyof == 0 ? spcenter[j] : (origfreq[j] + spwidth[j] / 2.0) / obsfreq[j]) + (j == 0
				&& leftcar ? leftcar * spwidth[0] / 2.0 / obsfreq[0] : 0.0));
		fprintf(fp, "refpt %f\n", datasize[j] + 0.5);
		fprintf(fp, "nuc %s\n", axisname[j]);
	}

	fclose(fp);

	if ((fp = fopen(spectra3d, "w")) == NULL) {
		fprintf(stderr, "Permission denied to write %s.\n", spectra3d);
		exit(EXIT_FAILURE);
	}

	fprintf(stderr, "%s", clean_string);
	fprintf(stderr, "Writing %s <- pipe2azara ...\r", spectra3d);

	matrix2d = fmalloc2d(datasize[1], datasize[0]);

	for (k = 0; k < datasize[2]; k++) {
		fprintf(stderr, "%s", clean_string);
		fprintf(stderr, "Receiving %d of %d ( -> pipe2azara) ...\r", ++count, datasize[2]);

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

				fseek(fp, (long) (headersize + (offset + block_id * block_volume) * sizeof(float)), SEEK_SET);

				fpwrite2bin_swap(fp, &(matrix2d[j][i]), 1, swapdata);
			}
		}
	}

	free_fmatrix2d(matrix2d);

	fclose(fp);

	fputc('\n', stderr);

	return 0;

	escape:free_fmatrix2d(matrix2d);

	fclose(fp);

	return 1;
}

int pullazara4d(char spectra4d[])
{
	FILE *fp;
	char parfile[MAXLONGNAME];
	int i, j, k, l, count = 0;
	int block_volume, block_i, block_j, block_k, block_l, block_id;
	int offset_i, offset_j, offset_k, offset_l, offset;
	short swapped = 0;
	float **matrix2d;

	if ((is_big_endian() && is_little_endian_float(header + sizeof(float) * 2, PIPE_HEADER[2])) || (is_little_endian()
			&& is_big_endian_float(header + sizeof(float) * 2, PIPE_HEADER[2])))
		swapped = 1;

	swapdata += swapped;

	make_dir(spectra4d);

	set_clean_string(clean_string);

	sprintf(parfile, "%s.par", spectra4d);

	if ((fp = fopen(parfile, "w")) == NULL) {
		fprintf(stderr, "Permission denied to write %s.\n", parfile);
		exit(EXIT_FAILURE);
	}

	fprintf(fp, "! this file = %s\n\n", parfile);

	fprintf(fp, "ndim 4\n");
	fprintf(fp, "file %s\n", spectra4d);

	block_volume = set_block_volume();

	for (j = 0; j < dimension; j++) {
		fprintf(fp, "\ndim %d\n", j + 1);
		fprintf(fp, "npts %d\n", datasize[j]);
		fprintf(fp, "block %d\n", blocksize[j]);
		fprintf(fp, "sw %f\n", spwidth[j]);
		fprintf(fp, "sf %f\n", obsfreq[j]);
		fprintf(fp, "refppm %f\n", (relyof == 0 ? spcenter[j] : (origfreq[j] + spwidth[j] / 2.0) / obsfreq[j]) + (j == 0
				&& leftcar ? leftcar * spwidth[0] / 2.0 / obsfreq[0] : 0.0));
		fprintf(fp, "refpt %f\n", datasize[j] + 0.5);
		fprintf(fp, "nuc %s\n", axisname[j]);
	}

	fclose(fp);

	if ((fp = fopen(spectra4d, "w")) == NULL) {
		fprintf(stderr, "Permission denied to write %s.\n", spectra4d);
		exit(EXIT_FAILURE);
	}

	fprintf(stderr, "%s", clean_string);
	fprintf(stderr, "Writing %s <- pipe2azara ...\r", spectra4d);

	matrix2d = fmalloc2d(datasize[1], datasize[0]);

	for (l = 0; l < datasize[3]; l++) {
		block_l = (int) (l / blocksize[3]);
		offset_l = l - block_l * blocksize[3];

		for (k = 0; k < datasize[2]; k++) {
			fprintf(stderr, "%s", clean_string);
			fprintf(stderr, "Receiving %d of %d ( -> pipe2azara) ...\r", ++count, datasize[2] * datasize[3]);

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

					fseek(fp, (long) (headersize + (offset + block_id * block_volume) * sizeof(float)), SEEK_SET);

					fpwrite2bin_swap(fp, &(matrix2d[j][i]), 1, swapdata);
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

	return 1;
}
