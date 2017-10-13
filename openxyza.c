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

int openxyza2d(char spectra2d[], float **mat2d)
{
	struct stat _stat;
	FILE *fp;
	char buffer[PIPE_HEADERSIZE];
	int data_plane = get_data_plane();
	long size = 0;
	float *vp;

	if ((fp = fopen(spectra2d, "r")) == NULL) {
		fprintf(stderr, "Spectra file %s: Couldn't open.\n", spectra2d);
		return EXIT_FAILURE;
	}

	stat(spectra2d, &_stat);
	size = _stat.st_size;

	fread(buffer, sizeof(char), size > PIPE_HEADERSIZE ? PIPE_HEADERSIZE : size, fp);

	/* CHECK FILE SIZE */
	vp = (float *) ((void *) (buffer));

	byteswap = 0;

	if (*vp != PIPE_HEADER[0] || *(vp + 2) != PIPE_HEADER[2]) {
		byteswap = 1;
		swapbyte(sizeof(float), size > PIPE_HEADERSIZE ? PIPE_HEADERSIZE : size, buffer);

		if (*vp != PIPE_HEADER[0] || *(vp + 2) != PIPE_HEADER[2]) {
			fprintf(stderr, "Spectra file: Not NMRPipe file.\n");
			goto escape;
		}
	}

	if (size != sizeof(float) * data_plane + PIPE_HEADERSIZE) {
		fprintf(stderr, "Spectra file %s: Partially broken. (Actual=%d Expected=%d)\n", spectra2d, (int) (size),
				(int) (sizeof(float)) * data_plane + PIPE_HEADERSIZE);
		goto escape;
	}

	fseek(fp, PIPE_HEADERSIZE, SEEK_SET);
	fread(&(mat2d[0][0]), sizeof(float), data_plane, fp);

	if (byteswap != 0)
		swapbyte(sizeof(float), data_plane * sizeof(float), (char *) (&(mat2d[0][0])));

	fclose(fp);

	return 0;

	escape:fclose(fp);

	return EXIT_FAILURE;
}

int openxyza3d(char spectra3d[], const int z, float **mat2d)
{
	char filename[MAXLONGNAME];

	sprintf(filename, spectra3d, z + 1);

	return openxyza2d(filename, mat2d);
}

int openxyza4d(char spectra4d[], const int z, const int a, float **mat2d)
{
	char filename[MAXLONGNAME];

	sprintf(filename, spectra4d, a + 1, z + 1);

	return openxyza2d(filename, mat2d);
}
