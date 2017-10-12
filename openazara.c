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

int openazara2d(char spectra2d[], float **mat2d)
{
	struct stat _stat;
	FILE *fp;
	int i, j;
	int block_volume = get_block_volume(), block_i, block_j, block_id;
	int offset_i, offset_j, offset;
	int data_volume = get_data_volume();
	long size = 0;

	if ((fp = fopen(spectra2d, "r")) == NULL) {
		fprintf(stderr, "Spectra file %s: Couldn't open.\n", spectra2d);
		return 1;
	}

	stat(spectra2d, &_stat);
	size = _stat.st_size;

	/* CHECK FILE SIZE */
	if (size != sizeof(float) * data_volume) {
		fprintf(stderr, "Spectra file %s: Partially broken. (Actual=%d Expected=%d)\n", spectra2d, (int) (size),
				(int) (sizeof(float)) * data_volume);
		goto escape;
	}

	for (j = 0; j < datasize[1]; j++) {
		block_j = (int) (j / blocksize[1]);
		offset_j = j - block_j * blocksize[1];

		for (i = 0; i < datasize[0]; i++) {
			block_i = (int) (i / blocksize[0]);
			offset_i = i - block_i * blocksize[0];

			block_id = block_i + block_j * unitsize[0];

			offset = offset_i + offset_j * blocksize[0];

			fseek(fp, (long) ((offset + block_id * block_volume) * sizeof(float)), SEEK_SET);

			fread(&(mat2d[j][i]), sizeof(char), sizeof(float), fp);
		}
	}

	if (byteswap != 0)
		swapbyte(sizeof(float), data_volume * sizeof(float), (char *) (&(mat2d[0][0])));

	fclose(fp);

	return 0;

	escape:fclose(fp);

	return 1;
}

int openazara3d(char spectra3d[], float ***mat3d)
{
	struct stat _stat;
	FILE *fp;
	int i, j, k;
	int block_volume = get_block_volume(), block_i, block_j, block_k, block_id;
	int offset_i, offset_j, offset_k, offset;
	int data_volume = get_data_volume();
	long size = 0;

	if ((fp = fopen(spectra3d, "r")) == NULL) {
		fprintf(stderr, "Spectra file %s: Couldn't open.\n", spectra3d);
		return 1;
	}

	stat(spectra3d, &_stat);
	size = _stat.st_size;

	/* CHECK FILE SIZE */
	if (size != sizeof(float) * data_volume) {
		fprintf(stderr, "Spectra file %s: Partially broken. (Actual=%d Expected=%d)\n", spectra3d, (int) (size),
				(int) (sizeof(float)) * data_volume);
		goto escape;
	}

	for (k = 0; k < datasize[2]; k++) {
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

				fseek(fp, (long) ((offset + block_id * block_volume) * sizeof(float)), SEEK_SET);

				fread(&(mat3d[k][j][i]), sizeof(char), sizeof(float), fp);
			}
		}
	}

	if (byteswap != 0)
		swapbyte(sizeof(float), data_volume * sizeof(float), (char *) (&(mat3d[0][0][0])));

	fclose(fp);

	return 0;

	escape:fclose(fp);

	return 1;
}

int openazara4d(char spectra4d[], float ****mat4d)
{
	struct stat _stat;
	FILE *fp;
	int i, j, k, l;
	int block_volume = get_block_volume(), block_i, block_j, block_k, block_l, block_id;
	int offset_i, offset_j, offset_k, offset_l, offset;
	int data_volume = get_data_volume();
	long size = 0;

	if ((fp = fopen(spectra4d, "r")) == NULL) {
		fprintf(stderr, "Spectra file %s: Couldn't open.\n", spectra4d);
		return 1;
	}

	stat(spectra4d, &_stat);
	size = _stat.st_size;

	/* CHECK FILE SIZE */
	if (size != sizeof(float) * data_volume) {
		fprintf(stderr, "Spectra file %s: Partially broken. (Actual=%d Expected=%d)\n", spectra4d, (int) (size),
				(int) (sizeof(float)) * data_volume);
		goto escape;
	}

	for (l = 0; l < datasize[3]; l++) {
		block_l = (int) (l / blocksize[3]);
		offset_l = l - block_l * blocksize[3];

		for (k = 0; k < datasize[2]; k++) {
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

					fseek(fp, (long) ((offset + block_id * block_volume) * sizeof(float)), SEEK_SET);

					fread(&(mat4d[l][k][j][i]), sizeof(char), sizeof(float), fp);
				}
			}
		}
	}

	if (byteswap != 0)
		swapbyte(sizeof(float), data_volume * sizeof(float), (char *) (&(mat4d[0][0][0][0])));

	fclose(fp);

	return 0;

	escape:fclose(fp);

	return 1;
}
