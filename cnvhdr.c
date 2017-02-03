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

static char _header_[PIPE_HEADERSIZE];
static float filesize;

static void cnvhdr2d(const char axis_option, const char fb, const int size, const int x, const int y);
static void cnvhdr3d(const char axis_option, const char fb, const int size, const int x, const int y, const int z);
static void cnvhdr4d(const char axis_option, const char fb, const int size, const int x, const int y, const int z,
		const int a);

void cnvhdr(const char axis_option, const char fb)
{
	switch (dimension) {
	case 2:
		cnvhdr2d(axis_option, fb, sizeof(float), 96, 100);
		cnvhdr2d(axis_option, fb, sizeof(float), 396, 876);
		break;

	case 3:
		cnvhdr3d(axis_option, fb, sizeof(float), 96, 100, 104);
		cnvhdr3d(axis_option, fb, sizeof(float), 396, 876, 60);
		break;

	case 4:
		cnvhdr4d(axis_option, fb, sizeof(float), 96, 100, 104, 108);
		cnvhdr4d(axis_option, fb, sizeof(float), 396, 876, 60, 128);
		break;
	}

	fwrite2mem(header + 1768, filesize);
}

static void cnvhdr2d(const char axis_option, const char fb, const int size, const int x, const int y)
{
	filesize = 1.0;

	switch (axis_option) {

	case 'x':
		break;

	case 'y':
		memcpy(_header_, header, PIPE_HEADERSIZE);

		if (fb == 'f' || fb == 'b') {
			memcpy(header + x, _header_ + y, size);
			memcpy(header + y, _header_ + x, size);
		}

		break;
	}
}

static void cnvhdr3d(const char axis_option, const char fb, const int size, const int x, const int y, const int z)
{
	switch (axis_option) {

	case 'x':
		filesize = datasize[2];
		break;

	case 'y':
		memcpy(_header_, header, PIPE_HEADERSIZE);

		if (fb == 'f' || fb == 'b') {
			memcpy(header + x, _header_ + y, size);
			memcpy(header + y, _header_ + x, size);

			filesize = datasize[2];
		}

		break;

	case 'z':
		memcpy(_header_, header, PIPE_HEADERSIZE);

		if (fb == 'f') {
			memcpy(header + x, _header_ + z, size);
			memcpy(header + y, _header_ + x, size);
			memcpy(header + z, _header_ + y, size);

			filesize = datasize[1];
		}

		if (fb == 'b') {
			memcpy(header + x, _header_ + y, size);
			memcpy(header + y, _header_ + z, size);
			memcpy(header + z, _header_ + x, size);

			filesize = datasize[0];
		}

		break;
	}
}

static void cnvhdr4d(const char axis_option, const char fb, const int size, const int x, const int y, const int z,
		const int a)
{
	switch (axis_option) {

	case 'x':
		filesize = datasize[2] * datasize[3];
		break;

	case 'y':
		memcpy(_header_, header, PIPE_HEADERSIZE);

		if (fb == 'f' || fb == 'b') {
			memcpy(header + x, _header_ + y, size);
			memcpy(header + y, _header_ + x, size);

			filesize = datasize[2] * datasize[3];
		}

		break;

	case 'z':
		memcpy(_header_, header, PIPE_HEADERSIZE);

		if (fb == 'f') {
			memcpy(header + x, _header_ + z, size);
			memcpy(header + y, _header_ + x, size);
			memcpy(header + z, _header_ + y, size);

			filesize = datasize[1] * datasize[3];
		}

		if (fb == 'b') {
			memcpy(header + x, _header_ + y, size);
			memcpy(header + y, _header_ + z, size);
			memcpy(header + z, _header_ + x, size);

			filesize = datasize[0] * datasize[3];
		}

		break;

	case 'a':
		memcpy(_header_, header, PIPE_HEADERSIZE);

		if (fb == 'f') {
			memcpy(header + x, _header_ + a, size);
			memcpy(header + y, _header_ + x, size);
			memcpy(header + z, _header_ + y, size);
			memcpy(header + a, _header_ + z, size);

			filesize = datasize[1] * datasize[2];
		}

		if (fb == 'b') {
			memcpy(header + x, _header_ + y, size);
			memcpy(header + y, _header_ + z, size);
			memcpy(header + z, _header_ + a, size);
			memcpy(header + a, _header_ + x, size);

			filesize = datasize[0] * datasize[3];
		}

		break;
	}
}
