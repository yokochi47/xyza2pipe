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

static char clean_string[MAXCHAR] =
{ "\r                                                                               \r" };

int pullnv2d(char spectra2d[])
{
	FILE *fp;
	int i, j, k, l, count = 0;
	int block_size, block_i, block_j, block_id;
	int offset_i, offset_j, offset;
	short swapped = 0;
	float **matrix2d;

	if ((is_big_endian() && is_little_endian_float(header + sizeof(float) * 2, PIPE_HEADER[2])) || (is_little_endian()
			&& is_big_endian_float(header + sizeof(float) * 2, PIPE_HEADER[2])))
		swapped = 1;

	swappar += swapped;
	swapdata += swapped;

	make_dir(spectra2d);

	set_clean_string(clean_string);

	if ((fp = fopen(spectra2d, "w")) == NULL) {
		fprintf(stderr, "Permission denied to write %s.\n", spectra2d);
		exit(EXIT_FAILURE);
	}

	headersize = NV_HEADERSIZE;

	memset(&(header[0]), 0, headersize * sizeof(char));

	j = 12;
	iwrite2mem_swap(header + j, headersize, swappar);

	j = 24;
	iwrite2mem_swap(header + j, dimension, swappar);

	sprintf(header, "%c%c%c%c", 0xcd, 0xab, 0x18, 0x34);

	j = 1076;
	for (i = 0; i < dimension; i++)
		memcpy(header + j + 128 * i, axisname[i], MAXASSNAME);

	for (i = 0; i < dimension; i++)
		blocksize[i] = datasize[i];

	while ((block_size = blocksize[0] * blocksize[1]) > NV_MAXBLOCKSIZE) {
		for (i = k = 0; i < dimension; i++) {
			for (j = 2; j < datasize[i] / 8; j++) {
				if (blocksize[i] % j == 0) {
					blocksize[i] /= j;
					break;
				}
			}

			if (j == datasize[i] / 8)
				k++;

			else if ((block_size = blocksize[0] * blocksize[1]) <= NV_MAXBLOCKSIZE)
				break;
		}

		if (i < dimension || (k > 0 && (block_size = blocksize[0] * blocksize[1]) <= pow(NV_MAXBLOCKSIZE, (float) dimension / (dimension - k))))
			break;
	}

	for (i = 0; i < dimension; i++)
		unitsize[i] = datasize[i] / blocksize[i];

	j = 20;
	iwrite2mem_swap(header + j, block_size, swappar);

	j = 1024;
	for (i = 0; i < dimension; i++)
		iwrite2mem_swap(header + j + 128 * i, datasize[i], swappar);

	j = 1028;
	for (i = 0; i < dimension; i++)
		iwrite2mem_swap(header + j + 128 * i, blocksize[i], swappar);

	j = 1032;
	for (i = 0; i < dimension; i++)
		iwrite2mem_swap(header + j + 128 * i, 16, swappar);

	j = 1036;
	for (i = 0, l = 1; i < dimension; i++) {
		iwrite2mem_swap(header + j + 128 * i, l, swappar);
		l *= 16;
	}

	j = 1040;
	for (i = 0; i < dimension; i++)
		iwrite2mem_swap(header + j + 128 * i, blocksize[i] - 1, swappar);

	j = 1044;
	for (i = 0, l = 1; i < dimension; i++) {
		iwrite2mem_swap(header + j + 128 * i, l, swappar);
		l += (5 - i * 2);
	}

	j = 1048;
	for (i = 0; i < dimension; i++)
		fwrite2mem_swap(header + j + 128 * i, obsfreq[i], swappar);

	j = 1052;
	for (i = 0; i < dimension; i++)
		fwrite2mem_swap(header + j + 128 * i, spwidth[i], swappar);

	if (leftcar)
		spcenter[0] += leftcar * spwidth[0] / 2.0 / obsfreq[0];

	j = 1056;
	for (i = 0; i < dimension; i++)
		fwrite2mem_swap(header + j + 128 * i, datasize[i] / 2.0, swappar);

	j = 1060;
	for (i = 0; i < dimension; i++)
		fwrite2mem_swap(header + j + 128 * i, (relyof == 0 ? spcenter[i] : (origfreq[i] + spwidth[i] / 2.0) / obsfreq[i]),
				swappar);

	j = 1064;
	for (i = 0; i < dimension; i++)
		header[j + 128 * i] = 0x03;

	fwrite(header, sizeof(char), headersize, fp);

	fprintf(stderr, "%s", clean_string);
	fprintf(stderr, "Writing %s <- pipe2nv ...\r", spectra2d);

	matrix2d = fmalloc2d(datasize[1], datasize[0]);

	fprintf(stderr, "%s", clean_string);
	fprintf(stderr, "Receiving %d of %d ( -> pipe2nv) ...\r", ++count, 1);

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

			fseek(fp, (long) (headersize + (offset + block_id * block_size) * sizeof(float)), SEEK_SET);

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

int pullnv3d(char spectra3d[])
{
	FILE *fp;
	int i, j, k, l, count = 0;
	int block_size, block_i, block_j, block_k, block_id;
	int offset_i, offset_j, offset_k, offset;
	short swapped = 0;
	float **matrix2d;

	if ((is_big_endian() && is_little_endian_float(header + sizeof(float) * 2, PIPE_HEADER[2])) || (is_little_endian()
			&& is_big_endian_float(header + sizeof(float) * 2, PIPE_HEADER[2])))
		swapped = 1;

	swappar += swapped;
	swapdata += swapped;

	make_dir(spectra3d);

	set_clean_string(clean_string);

	if ((fp = fopen(spectra3d, "w")) == NULL) {
		fprintf(stderr, "Permission denied to write %s.\n", spectra3d);
		exit(EXIT_FAILURE);
	}

	headersize = NV_HEADERSIZE;

	memset(&(header[0]), 0, headersize * sizeof(char));

	j = 12;
	iwrite2mem_swap(header + j, headersize, swappar);

	j = 24;
	iwrite2mem_swap(header + j, dimension, swappar);

	sprintf(header, "%c%c%c%c", 0xcd, 0xab, 0x18, 0x34);

	j = 1076;
	for (i = 0; i < dimension; i++)
		memcpy(header + j + 128 * i, axisname[i], MAXASSNAME);

	for (i = 0; i < dimension; i++)
		blocksize[i] = datasize[i];

	while ((block_size = blocksize[0] * blocksize[1] * blocksize[2]) > NV_MAXBLOCKSIZE) {
		for (i = k = 0; i < dimension; i++) {
			for (j = 2; j < datasize[i] / 8; j++) {
				if (blocksize[i] % j == 0) {
					blocksize[i] /= j;
					break;
				}
			}

			if (j == datasize[i] / 8)
				k++;

			else if ((block_size = blocksize[0] * blocksize[1] * blocksize[2]) <= NV_MAXBLOCKSIZE)
				break;
		}

		if (i < dimension || (k > 0 && (block_size = blocksize[0] * blocksize[1] * blocksize[2]) <= pow(NV_MAXBLOCKSIZE, (float) dimension / (dimension - k))))
			break;
	}

	for (i = 0; i < dimension; i++)
		unitsize[i] = datasize[i] / blocksize[i];

	j = 20;
	iwrite2mem_swap(header + j, block_size, swappar);

	j = 1024;
	for (i = 0; i < dimension; i++)
		iwrite2mem_swap(header + j + 128 * i, datasize[i], swappar);

	j = 1028;
	for (i = 0; i < dimension; i++)
		iwrite2mem_swap(header + j + 128 * i, blocksize[i], swappar);

	j = 1032;
	for (i = 0; i < dimension; i++)
		iwrite2mem_swap(header + j + 128 * i, 16, swappar);

	j = 1036;
	for (i = 0, l = 1; i < dimension; i++) {
		iwrite2mem_swap(header + j + 128 * i, l, swappar);
		l *= 16;
	}

	j = 1040;
	for (i = 0; i < dimension; i++)
		iwrite2mem_swap(header + j + 128 * i, blocksize[i] - 1, swappar);

	j = 1044;
	for (i = 0, l = 1; i < dimension; i++) {
		iwrite2mem_swap(header + j + 128 * i, l, swappar);
		l += (5 - i * 2);
	}

	j = 1048;
	for (i = 0; i < dimension; i++)
		fwrite2mem_swap(header + j + 128 * i, obsfreq[i], swappar);

	j = 1052;
	for (i = 0; i < dimension; i++)
		fwrite2mem_swap(header + j + 128 * i, spwidth[i], swappar);

	if (leftcar)
		spcenter[0] += leftcar * spwidth[0] / 2.0 / obsfreq[0];

	j = 1056;
	for (i = 0; i < dimension; i++)
		fwrite2mem_swap(header + j + 128 * i, datasize[i] / 2.0, swappar);

	j = 1060;
	for (i = 0; i < dimension; i++)
		fwrite2mem_swap(header + j + 128 * i, (relyof == 0 ? spcenter[i] : (origfreq[i] + spwidth[i] / 2.0) / obsfreq[i]),
				swappar);

	j = 1064;
	for (i = 0; i < dimension; i++)
		header[j + 128 * i] = 0x03;

	fwrite(header, sizeof(char), headersize, fp);

	fprintf(stderr, "%s", clean_string);
	fprintf(stderr, "Writing %s <- pipe2nv ...\r", spectra3d);

	matrix2d = fmalloc2d(datasize[1], datasize[0]);

	for (k = 0; k < datasize[2]; k++) {
		fprintf(stderr, "%s", clean_string);
		fprintf(stderr, "Receiving %d of %d ( -> pipe2nv) ...\r", ++count, datasize[2]);

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

				fseek(fp, (long) (headersize + (offset + block_id * block_size) * sizeof(float)), SEEK_SET);

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

int pullnv4d(char spectra4d[])
{
	FILE *fp;
	int i, j, k, l, count = 0;
	int block_size, block_i, block_j, block_k, block_l, block_id;
	int offset_i, offset_j, offset_k, offset_l, offset;
	short swapped = 0;
	float **matrix2d;

	if ((is_big_endian() && is_little_endian_float(header + sizeof(float) * 2, PIPE_HEADER[2])) || (is_little_endian()
			&& is_big_endian_float(header + sizeof(float) * 2, PIPE_HEADER[2])))
		swapped = 1;

	swappar += swapped;
	swapdata += swapped;

	make_dir(spectra4d);

	set_clean_string(clean_string);

	if ((fp = fopen(spectra4d, "w")) == NULL) {
		fprintf(stderr, "Permission denied to write %s.\n", spectra4d);
		exit(EXIT_FAILURE);
	}

	headersize = NV_HEADERSIZE;

	memset(&(header[0]), 0, headersize * sizeof(char));

	j = 12;
	iwrite2mem_swap(header + j, headersize, swappar);

	j = 24;
	iwrite2mem_swap(header + j, dimension, swappar);

	sprintf(header, "%c%c%c%c", 0xcd, 0xab, 0x18, 0x34);

	j = 1076;
	for (i = 0; i < dimension; i++)
		memcpy(header + j + 128 * i, axisname[i], MAXASSNAME);

	for (i = 0; i < dimension; i++)
		blocksize[i] = datasize[i];

	while ((block_size = blocksize[0] * blocksize[1] * blocksize[2] * blocksize[3]) > NV_MAXBLOCKSIZE) {
		for (i = k = 0; i < dimension; i++) {
			for (j = 2; j < datasize[i] / 8; j++) {
				if (blocksize[i] % j == 0) {
					blocksize[i] /= j;
					break;
				}
			}

			if (j == datasize[i] / 8)
				k++;

			else if ((block_size = blocksize[0] * blocksize[1] * blocksize[2] * blocksize[3]) <= NV_MAXBLOCKSIZE)
				break;
		}

		if (i < dimension || (k > 0 && (block_size = blocksize[0] * blocksize[1] * blocksize[2] * blocksize[3]) <= pow(NV_MAXBLOCKSIZE, (float) dimension / (dimension - k))))
			break;
	}

	for (i = 0; i < dimension; i++)
		unitsize[i] = datasize[i] / blocksize[i];

	j = 20;
	iwrite2mem_swap(header + j, block_size, swappar);

	j = 1024;
	for (i = 0; i < dimension; i++)
		iwrite2mem_swap(header + j + 128 * i, datasize[i], swappar);

	j = 1028;
	for (i = 0; i < dimension; i++)
		iwrite2mem_swap(header + j + 128 * i, blocksize[i], swappar);

	j = 1032;
	for (i = 0; i < dimension; i++)
		iwrite2mem_swap(header + j + 128 * i, 16, swappar);

	j = 1036;
	for (i = 0, l = 1; i < dimension; i++) {
		iwrite2mem_swap(header + j + 128 * i, l, swappar);
		l *= 16;
	}

	j = 1040;
	for (i = 0; i < dimension; i++)
		iwrite2mem_swap(header + j + 128 * i, blocksize[i] - 1, swappar);

	j = 1044;
	for (i = 0, l = 1; i < dimension; i++) {
		iwrite2mem_swap(header + j + 128 * i, l, swappar);
		l += (5 - i * 2);
	}

	j = 1048;
	for (i = 0; i < dimension; i++)
		fwrite2mem_swap(header + j + 128 * i, obsfreq[i], swappar);

	j = 1052;
	for (i = 0; i < dimension; i++)
		fwrite2mem_swap(header + j + 128 * i, spwidth[i], swappar);

	if (leftcar)
		spcenter[0] += leftcar * spwidth[0] / 2.0 / obsfreq[0];

	j = 1056;
	for (i = 0; i < dimension; i++)
		fwrite2mem_swap(header + j + 128 * i, datasize[i] / 2.0, swappar);

	j = 1060;
	for (i = 0; i < dimension; i++)
		fwrite2mem_swap(header + j + 128 * i, (relyof == 0 ? spcenter[i] : (origfreq[i] + spwidth[i] / 2.0) / obsfreq[i]),
				swappar);

	j = 1064;
	for (i = 0; i < dimension; i++)
		header[j + 128 * i] = 0x03;

	fwrite(header, sizeof(char), headersize, fp);

	fprintf(stderr, "%s", clean_string);
	fprintf(stderr, "Writing %s <- pipe2nv ...\r", spectra4d);

	matrix2d = fmalloc2d(datasize[1], datasize[0]);

	for (l = 0; l < datasize[3]; l++) {
		block_l = (int) (l / blocksize[3]);
		offset_l = l - block_l * blocksize[3];

		for (k = 0; k < datasize[2]; k++) {
			fprintf(stderr, "%s", clean_string);
			fprintf(stderr, "Receiving %d of %d ( -> pipe2nv) ...\r", ++count, datasize[2] * datasize[3]);

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

					fseek(fp, (long) (headersize + (offset + block_id * block_size) * sizeof(float)), SEEK_SET);

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
