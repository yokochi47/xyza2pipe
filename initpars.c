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

char header[PIPE_HEADERSIZE];
char axislabel[4][MAXASSNAME + 1] = { {""} };

char axisname[4][MAXASSNAME + 1];

char dimension, _dimension, byteswap;

char swapdata = 0, swappar = 0, usrlabel = 0, usrshift = 0;
char leftcar = 0, extleft = 0, adjcar = 0, adjh2o = 0, relyof = 0;
char usrphase[4] = { 0 };

short headersize;
short blocksize[4];

int datasize[4], datasize_orig[4];
int unitsize[4];

float obsfreq[4];
float spcenter[4];
float origfreq[4];
float spwidth[4];
float usrcenter[4] = { NULLPPM };

float phase[4][2];

const float PIPE_HEADER[3] = { PIPE_HEADER_0, PIPE_HEADER_1, PIPE_HEADER_2 };

char clean_string[MAXCHAR] =
{ "\r                                                                               \r" };

float get_indirect_planes()
{

	switch (dimension) {
	case 2:
		return 1.0;
	case 3:
		return (float) datasize[2];
	case 4:
		return (float) datasize[2] * datasize[3];
	default:
		return 0.0;
	}

}

float get_orig_indirect_planes()
{

	switch (dimension) {
	case 2:
		return 1.0;
	case 3:
		return (float) datasize_orig[2];
	case 4:
		return (float) datasize_orig[2] * datasize_orig[3];
	default:
		return 0.0;
	}

}

int get_data_plane()
{
	return datasize[1] * datasize[0];
}

int get_data_volume()
{
	int data_volume = get_data_plane();

	switch (dimension) {
	case 4:
		data_volume *= datasize[3];
		/* no break */

	case 3:
		data_volume *= datasize[2];
		/* no break */

	}

	return data_volume;
}

int get_block_plane()
{
	return blocksize[1] * blocksize[0];
}

int get_block_volume()
{
	int block_volume = get_block_plane();

	switch (dimension) {
	case 4:
		block_volume *= blocksize[3];
		/* no break */

	case 3:
		block_volume *= blocksize[2];
		/* no break */

	}

	return block_volume;
}

int set_block_volume()
{
	int i, j, k;

	for (i = 0; i < dimension; i++)
		blocksize[i] = datasize[i];

	while (get_block_volume() > MAXBLOCKSIZE) {
		for (i = k = 0; i < dimension; i++) {
			for (j = 2; j < datasize[i] / 8; j++) {
				if (blocksize[i] % j == 0) {
					blocksize[i] /= j;

					break;

				}
			}

			if (j == datasize[i] / 8)
				k++;

			else if (get_block_volume() <= MAXBLOCKSIZE)
				break;
		}

		if (i < dimension || (k > 0 && get_block_volume() <= pow(MAXBLOCKSIZE, (float) dimension / (dimension - k))))
			break;
	}

	for (i = 0; i < dimension; i++)
		unitsize[i] = datasize[i] / blocksize[i];

	return get_block_volume();
}
