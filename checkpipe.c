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

static char _out_ = 0;

int checkpipe()
{
	char buffer[PIPE_HEADERSIZE];
	int fidsize[4];
	int j, k, axisorder[4] = { 0 }, orderaxis[4];
	float *vp;

	if (usrlabel != 0 || usrshift != 0)
		_out_ = 1;

	if (_out_)
		fprintf(stderr, "Reading (stdin) stream ...\n");

	fread(buffer, sizeof(char), PIPE_HEADERSIZE, stdin);

	/* FILE SEARCH */
	vp = (float *) ((void *) (buffer));

	byteswap = 0;

	if (*vp != PIPE_HEADER[0] || *(vp + 2) != PIPE_HEADER[2]) {
		byteswap = (byteswap + 1) % 2;
		swapbyte(sizeof(float), PIPE_HEADERSIZE, buffer);

		if (*vp != PIPE_HEADER[0] || *(vp + 2) != PIPE_HEADER[2]) {
			fprintf(stderr, "Spectra file: Not NMRPipe file.\n");
			return 1;
		}
	}

	vp = (float *) ((void *) (&buffer[36]));

	if (*vp < 2.0 || *vp > 4.0) {
		fprintf(stderr, "Spectra file: Not 2D/3D/4D experiment.\n");
		return 1;
	}

	dimension = 0;
	if (*vp == 2.0)
		dimension = 2;
	else if (*vp == 3.0)
		dimension = 3;
	else if (*vp == 4.0)
		dimension = 4;

	vp = (float *) ((void *) (&buffer[880]));

	if (*vp == 0.0) {
		fprintf(stderr, "Spectra file: Time domain dataset.\nExecute Fourier Transform.\n");
		return 1;
	}

	vp = (float *) ((void *) (&buffer[888]));

	if (*vp == 0.0) {
		fprintf(stderr, "Spectra file: Time domain dataset.\nExecute Fourier Transform.\n");
		return 1;
	}

	vp = (float *) ((void *) (&buffer[220]));

	if (*vp != 1.0) {
		fprintf(stderr, "Spectra file: Complex dataset.\nRemove imaginary part.\n");
		return 1;
	}

	vp = (float *) ((void *) (&buffer[224]));

	if (*vp != 1.0) {
		fprintf(stderr, "Spectra file: Complex dataset.\nRemove imaginary part.\n");
		return 1;
	}

	vp = (float *) ((void *) (&buffer[884]));

	if (*vp != 0.0) {
		fprintf(stderr, "Spectra file: Warning! Transposed dataset.\n");
		fwrite2mem(header + 884, 0.0);
	}

	/* READ AXIS ORDER */
	vp = (float *) ((void *) (&buffer[96]));
	if (*vp == 2.0)
		axisorder[0] = 0;
	if (*vp == 1.0)
		axisorder[0] = 1;
	if (*vp == 3.0)
		axisorder[0] = 2;
	if (*vp == 4.0)
		axisorder[0] = 3;

	vp = (float *) ((void *) (&buffer[100]));
	if (*vp == 2.0)
		axisorder[1] = 0;
	if (*vp == 1.0)
		axisorder[1] = 1;
	if (*vp == 3.0)
		axisorder[1] = 2;
	if (*vp == 4.0)
		axisorder[1] = 3;

	vp = (float *) ((void *) (&buffer[104]));
	if (*vp == 2.0)
		axisorder[2] = 0;
	if (*vp == 1.0)
		axisorder[2] = 1;
	if (*vp == 3.0)
		axisorder[2] = 2;
	if (*vp == 4.0)
		axisorder[2] = 3;

	vp = (float *) ((void *) (&buffer[108]));
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

	switch (dimension) {
	case 2:
		if (axisorder[0] + axisorder[1] != 1 || axisorder[0] * axisorder[1] != 0 || axisorder[0] == axisorder[1]) {
			fprintf(stderr, "Spectra file: Warning! Irregularly transposed dataset.\n");

			for (j = 0; j < 4; j++)
				axisorder[j] = orderaxis[j] = j;

			fwrite2mem(header + 96, 2.0);
			fwrite2mem(header + 100, 1.0);
			fwrite2mem(header + 104, 3.0);
			fwrite2mem(header + 108, 4.0);
		}
		break;
	case 3:
		if (axisorder[0] + axisorder[1] + axisorder[2] != 3 || axisorder[0] * axisorder[1] * axisorder[2] != 0
				|| axisorder[0] == axisorder[1] || axisorder[1] == axisorder[2] || axisorder[2] == axisorder[0]) {
			fprintf(stderr, "Spectra file: Warning! Irregularly transposed dataset.\n");

			for (j = 0; j < 4; j++)
				axisorder[j] = orderaxis[j] = j;

			fwrite2mem(header + 96, 2.0);
			fwrite2mem(header + 100, 1.0);
			fwrite2mem(header + 104, 3.0);
			fwrite2mem(header + 108, 4.0);
		}
		break;
	case 4:
		if (axisorder[0] + axisorder[1] + axisorder[2] + axisorder[3] != 6
				|| axisorder[0] * axisorder[1] * axisorder[2] * axisorder[3] != 0 || axisorder[0] == axisorder[1]
																											   || axisorder[0] == axisorder[2] || axisorder[0] == axisorder[3] || axisorder[1] == axisorder[2]
																																																			|| axisorder[1] == axisorder[3] || axisorder[2] == axisorder[3]) {
			fprintf(stderr, "Spectra file: Warning! Irregularly transposed dataset.\n");

			for (j = 0; j < 4; j++)
				axisorder[j] = orderaxis[j] = j;

			fwrite2mem(header + 96, 2.0);
			fwrite2mem(header + 100, 1.0);
			fwrite2mem(header + 104, 3.0);
			fwrite2mem(header + 108, 4.0);
		}
		break;
	}

	/* READ AXIS NAME */
	j = 64;

	if (byteswap == 1)
		swapbyte(sizeof(float), 32, buffer + j);

	for (k = 0; k < dimension; k++) {
		memcpy(axisname[k], buffer + j + MAXASSNAME * axisorder[k], MAXASSNAME);
		if (*axislabel[k] != 0)
			memcpy(axisname[k], axislabel[k], MAXASSNAME);

		if (usrlabel == 0)
			checklabel(axisname[k]);
	}

	switch (dimension) {
	case 2:
		if (_out_)
			fprintf(stderr, "Axis Label | %8s  %8s\n", axisname[0], axisname[1]);

		/* READ OBS. FREQ. */
		vp = (float *) ((void *) (&buffer[476]));
		obsfreq[orderaxis[0]] = *vp;

		vp = (float *) ((void *) (&buffer[872]));
		obsfreq[orderaxis[1]] = *vp;

		if (_out_)
			fprintf(stderr, "Obs. Freq. | %8.3f  %8.3f [MHz]\n", obsfreq[0], obsfreq[1]);

		/* READ SPECTRAL CENTER */
		vp = (float *) ((void *) (&buffer[264]));
		spcenter[orderaxis[0]] = *vp;

		vp = (float *) ((void *) (&buffer[268]));
		spcenter[orderaxis[1]] = *vp;

		if (_out_)
			fprintf(stderr, "Spec.Center| %8.2f  %8.2f [ppm]\n", spcenter[0], spcenter[1]);

		/* READ ORIG. FREQ. */
		vp = (float *) ((void *) (&buffer[404]));
		origfreq[orderaxis[0]] = *vp;

		vp = (float *) ((void *) (&buffer[996]));
		origfreq[orderaxis[1]] = *vp;

		if (_out_)
			fprintf(stderr, "Orig.Freq. | %8.2f  %8.2f [Hz]\n", origfreq[0], origfreq[1]);

		/* READ SPECTRAL WIDTH */
		vp = (float *) ((void *) (&buffer[400]));
		spwidth[orderaxis[0]] = *vp;

		vp = (float *) ((void *) (&buffer[916]));
		spwidth[orderaxis[1]] = *vp;

		if (_out_)
			fprintf(stderr, "Spec.Width | %8.2f  %8.2f [Hz]\n", spwidth[0], spwidth[1]);

		/* READ DATA SIZE */
		vp = (float *) ((void *) (&buffer[396]));
		datasize[0] = (int) (*vp);

		vp = (float *) ((void *) (&buffer[876]));
		datasize[1] = (int) (*vp);

		if (_out_)
			fprintf(stderr, "Data Size  | %8d  %8d\n", datasize[0], datasize[1]);

		/* READ FID SIZE */
		vp = (float *) ((void *) (&buffer[1544]));
		fidsize[orderaxis[0]] = (int) (*vp);

		vp = (float *) ((void *) (&buffer[1548]));
		fidsize[orderaxis[1]] = (int) (*vp);

		if (_out_)
			fprintf(stderr, "FID  Size  | %8d  %8d\n", fidsize[0], fidsize[1]);

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

		break;

	case 3:
		if (_out_)
			fprintf(stderr, "Axis Label | %8s  %8s  %8s\n", axisname[0], axisname[1], axisname[2]);

		/* READ OBS. FREQ. */
		vp = (float *) ((void *) (&buffer[476]));
		obsfreq[orderaxis[0]] = *vp;

		vp = (float *) ((void *) (&buffer[872]));
		obsfreq[orderaxis[1]] = *vp;

		vp = (float *) ((void *) (&buffer[40]));
		obsfreq[orderaxis[2]] = *vp;

		if (_out_)
			fprintf(stderr, "Obs. Freq. | %8.3f  %8.3f  %8.3f [MHz]\n", obsfreq[0], obsfreq[1], obsfreq[2]);

		/* READ SPECTRAL CENTER */
		vp = (float *) ((void *) (&buffer[264]));
		spcenter[orderaxis[0]] = *vp;

		vp = (float *) ((void *) (&buffer[268]));
		spcenter[orderaxis[1]] = *vp;

		vp = (float *) ((void *) (&buffer[272]));
		spcenter[orderaxis[2]] = *vp;

		if (_out_)
			fprintf(stderr, "Spec.Center| %8.2f  %8.2f  %8.2f [ppm]\n", spcenter[0], spcenter[1], spcenter[2]);

		/* READ ORIG. FREQ. */
		vp = (float *) ((void *) (&buffer[404]));
		origfreq[orderaxis[0]] = *vp;

		vp = (float *) ((void *) (&buffer[996]));
		origfreq[orderaxis[1]] = *vp;

		vp = (float *) ((void *) (&buffer[48]));
		origfreq[orderaxis[2]] = *vp;

		if (_out_)
			fprintf(stderr, "Orig.Freq. | %8.2f  %8.2f  %8.2f [Hz]\n", origfreq[0], origfreq[1], origfreq[2]);

		/* READ SPECTRAL WIDTH */
		vp = (float *) ((void *) (&buffer[400]));
		spwidth[orderaxis[0]] = *vp;

		vp = (float *) ((void *) (&buffer[916]));
		spwidth[orderaxis[1]] = *vp;

		vp = (float *) ((void *) (&buffer[44]));
		spwidth[orderaxis[2]] = *vp;

		if (_out_)
			fprintf(stderr, "Spec.Width | %8.2f  %8.2f  %8.2f [Hz]\n", spwidth[0], spwidth[1], spwidth[2]);

		/* READ DATA SIZE */
		vp = (float *) ((void *) (&buffer[396]));
		datasize[0] = (int) (*vp);

		vp = (float *) ((void *) (&buffer[876]));
		datasize[1] = (int) (*vp);

		vp = (float *) ((void *) (&buffer[60]));
		datasize[2] = (int) (*vp);

		if (_out_)
			fprintf(stderr, "Data Size  | %8d  %8d  %8d\n", datasize[0], datasize[1], datasize[2]);

		/* READ FID SIZE */
		vp = (float *) ((void *) (&buffer[1544]));
		fidsize[orderaxis[0]] = (int) (*vp);

		vp = (float *) ((void *) (&buffer[1548]));
		fidsize[orderaxis[1]] = (int) (*vp);

		vp = (float *) ((void *) (&buffer[1552]));
		fidsize[orderaxis[2]] = (int) (*vp);

		if (_out_)
			fprintf(stderr, "FID  Size  | %8d  %8d  %8d\n", fidsize[0], fidsize[1], fidsize[2]);

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

		break;

	case 4:
		if (_out_)
			fprintf(stderr, "Axis Label | %8s  %8s  %8s  %8s\n", axisname[0], axisname[1], axisname[2], axisname[3]);

		/* READ OBS. FREQ. */
		vp = (float *) ((void *) (&buffer[476]));
		obsfreq[orderaxis[0]] = *vp;

		vp = (float *) ((void *) (&buffer[872]));
		obsfreq[orderaxis[1]] = *vp;

		vp = (float *) ((void *) (&buffer[40]));
		obsfreq[orderaxis[2]] = *vp;

		vp = (float *) ((void *) (&buffer[112]));
		obsfreq[orderaxis[3]] = *vp;

		if (_out_)
			fprintf(stderr, "Obs. Freq. | %8.3f  %8.3f  %8.3f  %8.3f [MHz]\n", obsfreq[0], obsfreq[1], obsfreq[2], obsfreq[3]);

		/* READ SPECTRAL CENTER */
		vp = (float *) ((void *) (&buffer[264]));
		spcenter[orderaxis[0]] = *vp;

		vp = (float *) ((void *) (&buffer[268]));
		spcenter[orderaxis[1]] = *vp;

		vp = (float *) ((void *) (&buffer[272]));
		spcenter[orderaxis[2]] = *vp;

		vp = (float *) ((void *) (&buffer[276]));
		spcenter[orderaxis[3]] = *vp;

		if (_out_)
			fprintf(stderr, "Spec.Center| %8.2f  %8.2f  %8.2f  %8.2f [ppm]\n", spcenter[0], spcenter[1], spcenter[2],
					spcenter[3]);

		/* READ ORIG. FREQ. */
		vp = (float *) ((void *) (&buffer[404]));
		origfreq[orderaxis[0]] = *vp;

		vp = (float *) ((void *) (&buffer[996]));
		origfreq[orderaxis[1]] = *vp;

		vp = (float *) ((void *) (&buffer[48]));
		origfreq[orderaxis[2]] = *vp;

		vp = (float *) ((void *) (&buffer[120]));
		origfreq[orderaxis[3]] = *vp;

		if (_out_)
			fprintf(stderr, "Orig.Freq. | %8.2f  %8.2f  %8.2f  %8.2f [Hz]\n", origfreq[0], origfreq[1], origfreq[2],
					origfreq[3]);

		/* READ SPECTRAL WIDTH */
		vp = (float *) ((void *) (&buffer[400]));
		spwidth[orderaxis[0]] = *vp;

		vp = (float *) ((void *) (&buffer[916]));
		spwidth[orderaxis[1]] = *vp;

		vp = (float *) ((void *) (&buffer[44]));
		spwidth[orderaxis[2]] = *vp;

		vp = (float *) ((void *) (&buffer[116]));
		spwidth[orderaxis[3]] = *vp;

		if (_out_)
			fprintf(stderr, "Spec.Width | %8.2f  %8.2f  %8.2f  %8.2f [Hz]\n", spwidth[0], spwidth[1], spwidth[2], spwidth[3]);

		/* READ DATA SIZE */
		vp = (float *) ((void *) (&buffer[396]));
		datasize[0] = (int) (*vp);

		vp = (float *) ((void *) (&buffer[876]));
		datasize[1] = (int) (*vp);

		vp = (float *) ((void *) (&buffer[60]));
		datasize[2] = (int) (*vp);

		vp = (float *) ((void *) (&buffer[128]));
		datasize[3] = (int) (*vp);

		if (_out_)
			fprintf(stderr, "Data Size  | %8d  %8d  %8d  %8d\n", datasize[0], datasize[1], datasize[2], datasize[3]);

		/* READ FID SIZE */
		vp = (float *) ((void *) (&buffer[1544]));
		fidsize[orderaxis[0]] = (int) (*vp);

		vp = (float *) ((void *) (&buffer[1548]));
		fidsize[orderaxis[1]] = (int) (*vp);

		vp = (float *) ((void *) (&buffer[1552]));
		fidsize[orderaxis[2]] = (int) (*vp);

		vp = (float *) ((void *) (&buffer[1556]));
		fidsize[orderaxis[3]] = (int) (*vp);

		if (_out_)
			fprintf(stderr, "FID  Size  | %8d  %8d  %8d  %8d\n", fidsize[0], fidsize[1], fidsize[2], fidsize[3]);

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

		break;
	}

	/* CHECK FILE SIZE */
	memcpy(header, buffer, PIPE_HEADERSIZE);

	if (usrshift != 0) {
		switch (dimension) {
		case 4:
			fwrite2mem(header + 276, spcenter[3]);

			fwrite2mem(header + 120, origfreq[3]);
			/* no break */

		case 3:
			fwrite2mem(header + 272, spcenter[2]);

			fwrite2mem(header + 48, origfreq[2]);
			/* no break */

		case 2:
			fwrite2mem(header + 264, spcenter[0]);
			fwrite2mem(header + 268, spcenter[1]);

			fwrite2mem(header + 404, origfreq[0]);
			fwrite2mem(header + 996, origfreq[1]);

			break;
		}
	}

	return 0;
}

void make_dir(char filename[])
{
	struct stat _stat;
	char _dirname[MAXLONGNAME], command[MAXLONGNAME];

	dirnamecpy(_dirname, filename);

	if (strcmp(_dirname, ".") != 0) {

		if (stat(_dirname, &_stat) != 0 || (_stat.st_mode & S_IFMT) != S_IFDIR) {
			sprintf(command, "mkdir -p %s", _dirname);
			system(command);

			if (stat(_dirname, &_stat) != 0 || (_stat.st_mode & S_IFMT) != S_IFDIR) {
				fprintf(stderr, "Permission denied to make %s.\n", _dirname);
				exit(EXIT_FAILURE);
			}
		}
	}
}

void set_clean_string(char clean_string[])
{
	FILE *fp;
	char string[MAXCHAR];
	int cols = 0;

	if ((fp = popen("tput cols", "r")) != NULL) {
		fgets(string, MAXCHAR, fp);
		cols = atoi(string);
		pclose(fp);
	}

	if (cols > 1 && cols + 1 < MAXCHAR) {
		clean_string[0] = '\r';
		memset(&(clean_string[1]), ' ', (cols - 1) * sizeof(char));
		clean_string[cols] = '\r';
		clean_string[cols + 1] = '\0';
	}
}
