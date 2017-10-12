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

int checkazara(char filename[])
{
	struct stat _stat;
	FILE *fp, *par, *scr, *fidpar;
	char parline[MAXCHAR], *argv[MAXVARS], buffer[MAXCHAR];
	char parfile[MAXLONGNAME], scrfile[MAXLONGNAME], fidparfile[MAXLONGNAME];
	int fidsize[4] = { 0 };
	short apod_code[4] = { SINE_BELL };
	int j = 0, k, l, repeat, column_len;
	long size = 0;
	float *vp, *vp2;
	float fidwidth[4] = { 0.0 };
	float apod_par[4][3] = { {0.0} };

	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Spectra file %s: Couldn't open.\n", filename);
		return 1;
	}

	stat(filename, &_stat);
	size = _stat.st_size;

	fprintf(stderr, "Reading %s ...\n", filename);

	fread(buffer, sizeof(char), size > MAXCHAR ? MAXCHAR : size, fp);

	if (byteswap == 1)
		swapbyte(sizeof(float), size > MAXCHAR ? MAXCHAR : size, buffer);

	vp = (float *) ((void *) (buffer));

	vp2 = (float *) ((void *) (&buffer[sizeof(float)]));

	fclose(fp);

	if (fabs(*vp - *vp2) < 1.0e-6)
		byteswap = (byteswap + 1) % 2;

	sprintf(parfile, "%s.par", filename);

	if ((par = fopen(parfile, "r")) == NULL) {
		fprintf(stderr, "Parameter file %s: Couldn't open.\n", parfile);
		return 1;
	}

	/* READ PARAMETERS */
	alloc_arg(0, MAXFILENAME, argv);

	fgets(parline, MAXCHAR, par);

	if (strstr(parline, "! this file = ") == NULL) {
		fprintf(stderr, "Parameter file %s: Not Azara parameter file.\nCouldn't find description such as '! this file = '.\n",
				parfile);
		goto escape;
	}

	while (strncmp(parline, "ndim", 4) != 0) {
		if (fgets(parline, MAXCHAR, par) == NULL)
			break;
	}

	line2arg(parline, ' ', argv);

	dimension = atoi(argv[1]);

	if (dimension < 2 || dimension > 4) {
		fprintf(stderr, "Spectra file: Not 2D/3D/4D experiment.\n");
		goto escape;
	}

	for (repeat = 0; repeat < dimension; repeat++) {

		while (strncmp(parline, "dim", 3) != 0) {
			if (fgets(parline, MAXCHAR, par) == NULL)
				break;
		}

		line2arg(parline, ' ', argv);

		j = atoi(argv[1]) - 1;

		if (j < 0 || j >= dimension) {
			fprintf(stderr, "Spectra file: Invalid dimension number.\n");
			goto escape;
		}

		do {

			if (fgets(parline, MAXCHAR, par) == NULL)
				break;

			column_len = line2arg(parline, ' ', argv);

			if (strcmp(argv[0], "npts") == 0) {

				/* READ DATA SIZE */

				datasize[j] = atoi(argv[1]);
			}

			else if (strcmp(argv[0], "block") == 0) {

				/* READ BLOCK SIZE */

				blocksize[j] = atoi(argv[1]);
				if ((datasize[j] % blocksize[j]) != 0)
					datasize[j] = ((int) (datasize[j] / blocksize[j]) + 1) * blocksize[j];
			}

			else if (strcmp(argv[0], "sw") == 0) {

				/* READ SPECTRAL WIDTH */

				spwidth[j] = atof(argv[1]);
			}

			else if (strcmp(argv[0], "sf") == 0) {

				/* READ OBS. FREQ. */

				obsfreq[j] = atof(argv[1]);
			}

			else if (strcmp(argv[0], "refppm") == 0) {

				/* READ SPECTRAL CENTER */

				spcenter[j] = atof(argv[1]);
			}

			else if (strcmp(argv[0], "refpt") == 0) {
			}

			else if (strcmp(argv[0], "nuc") == 0) {

				/* READ AXIS NAME */

				if (*axislabel[j] != 0)
					strncpy(axisname[j], axislabel[j], MAXASSNAME);
				else
					strncpy(axisname[j], argv[1], MAXASSNAME);

				if (usrlabel == 0)
					checklabel(axisname[j]);
			}

		} while (column_len != 0);
	}

	/* READ SCRIPT FILE */

	rewind(par);

	while (fgets(parline, MAXCHAR, par) != NULL) {

		strreplace(parline, "\t", " ");

		line2arg(parline, ' ', argv);

		if (strcmp(argv[1], "script") != 0 || strcmp(argv[2], "file") != 0)
			continue;

		strncpy(scrfile, argv[4], MAXLONGNAME);

		if ((scr = fopen(scrfile, "r")) == NULL) {
			fprintf(stderr, "Script file %s: Couldn't open.\n", scrfile);
			goto escape;
		}

		while (fgets(parline, MAXCHAR, scr) != NULL) {
			line2arg(parline, ' ', argv);

			if (strcmp(argv[0], "input") != 0)
				continue;

			strncpy(fidparfile, argv[1], MAXLONGNAME);

			if ((fidpar = fopen(fidparfile, "r")) == NULL) {
				fclose(scr);

				fprintf(stderr, "Parameter file %s: Couldn't open.\n", fidparfile);
				goto escape;
			}

			fgets(parline, MAXCHAR, fidpar);

			while (strncmp(parline, "ndim", 4) != 0) {
				if (fgets(parline, MAXCHAR, fidpar) == NULL)
					break;
			}

			line2arg(parline, ' ', argv);

			if (dimension != atoi(argv[1])) {
				fclose(scr);
				fclose(fidpar);

				fprintf(stderr, "Parameter file %s: Unmached dimension number.\n", fidparfile);
				goto escape;
			}

			for (repeat = 0; repeat < dimension; repeat++) {

				while (strncmp(parline, "dim", 3) != 0) {
					if (fgets(parline, MAXCHAR, fidpar) == NULL)
						break;
				}

				line2arg(parline, ' ', argv);

				j = atoi(argv[1]) - 1;

				if (j < 0 || j >= dimension) {
					fprintf(stderr, "Parameter file %s: Invalid dimension number.\n", fidparfile);
					goto escape;
				}

				do {

					if (fgets(parline, MAXCHAR, fidpar) == NULL)
						break;

					column_len = line2arg(parline, ' ', argv);

					if (strcmp(argv[0], "npts") == 0) {

						/* READ DATA SIZE */

						fidsize[j] = datasize_orig[j] = atoi(argv[1]) / 2;
					}

					else if (strcmp(argv[0], "sw") == 0) {

						/* READ SPECTRAL WIDTH */

						fidwidth[j] = atof(argv[1]);
					}

				} while (column_len != 0);
			}

			fclose(fidpar);

			break;
		}

		fclose(scr);

		break;
	}

	switch (dimension) {
	case 2:
		for (j = 0; j < dimension; j++) {
			origfreq[j] = spcenter[j] * obsfreq[j] - spwidth[j] / 2.0 + spwidth[j] / (float) (datasize[j]);

			if (*axisname[j] == 'H' && strncmp(axisname[j], "HA", MAXASSNAME) != 0 && *axisname[(j + 1) % dimension] == 'N'
					&& spcenter[j] >= 4.5 && spcenter[j] <= 5.0 && spwidth[j] / obsfreq[j] < 8.5)
				origfreq[j] += spwidth[j] / 2.0;
		}

		break;
	case 3:
		for (j = 0; j < dimension; j++) {
			origfreq[j] = spcenter[j] * obsfreq[j] - spwidth[j] / 2.0 + spwidth[j] / (float) (datasize[j]);

			if (*axisname[j] == 'H' && strncmp(axisname[j], "HA", MAXASSNAME) != 0 && (*axisname[(j + 1) % dimension] == 'N'
					|| *axisname[(j + 2) % dimension] == 'N') && spcenter[j] >= 4.5 && spcenter[j] <= 5.0
					&& spwidth[j] / obsfreq[j] < 8.5)
				origfreq[j] += spwidth[j] / 2.0;
		}

		break;
	case 4:
		for (j = 0; j < dimension; j++) {
			origfreq[j] = spcenter[j] * obsfreq[j] - spwidth[j] / 2.0 + spwidth[j] / (float) (datasize[j]);

			if (*axisname[j] == 'H' && strncmp(axisname[j], "HA", MAXASSNAME) != 0 && (*axisname[(j + 1) % dimension] == 'N'
					|| *axisname[(j + 2) % dimension] == 'N' || *axisname[(j + 3) % dimension] == 'N') && spcenter[j] >= 4.5
					&& spcenter[j] <= 5.0 && spwidth[j] / obsfreq[j] < 8.5)
				origfreq[j] += spwidth[j] / 2.0;
		}

		break;
	}

	/* READ PROCESS LOG */

	rewind(par);

	while (fgets(parline, MAXCHAR, par) != NULL) {

		strreplace(parline, "\t", " ");

		line2arg(parline, ' ', argv);

		if (strcmp(argv[1], "Script") == 0) {

			if (strstr(parline, "(dimension 1)"))
				j = 0;
			else if (strstr(parline, "(dimension 2)"))
				j = 1;
			else if (strstr(parline, "(dimension 3)"))
				j = 2;
			else if (strstr(parline, "(dimension 4)"))
				j = 3;

			continue;
		}

		else if (strcmp(argv[1], "Command") == 0) {

			/* arrange */

			if (strcmp(argv[3], "lower") == 0) {
				l = blocksize[j];
				while (l < atoi(argv[4]) - 1)
					l += blocksize[j];

				if (l != atoi(argv[4]) - 1)
					l -= blocksize[j];

				if (l != atoi(argv[4]) - 1)
					origfreq[j] += (float) (l - (atoi(argv[4]) - 1)) / (float) (datasize_orig[j]) * spwidth[j];

				spcenter[j] = origfreq[j] / obsfreq[j];
			}

			else if (strcmp(argv[3], "upper") == 0) {
				l = blocksize[j];
				while (l < atoi(argv[4]))
					l += blocksize[j];

				origfreq[j] += (float) (datasize_orig[j] - l) / (float) (datasize_orig[j]) * spwidth[j];

				spcenter[j] = origfreq[j] / obsfreq[j];
			}

			else if (strcmp(argv[3], "range") == 0) {
				l = blocksize[j];
				while (l < atoi(argv[5]))
					l += blocksize[j];
			}

			else if (strcmp(argv[3], "shift") == 0) {
				l = blocksize[j];
				while (l < atoi(argv[4]))
					l += blocksize[j];
			}

			else if (strcmp(argv[3], "zerofill") == 0)
				datasize_orig[j] *= (int) (pow(2.0, atoi(argv[4])));

			/* weight */

			if (strcmp(argv[3], "decay") == 0) {
				apod_code[j] = EXPONENTIAL;
				apod_par[j][0] = -log(atof(argv[4])) * fidwidth[j] / M_PI / (fidsize[j] - 1);
				apod_par[j][1] = 0.0;
				apod_par[j][2] = 0.0;
			}

			else if (strcmp(argv[3], "decay_sw") == 0) {
				apod_code[j] = EXPONENTIAL;
				apod_par[j][0] = atof(argv[4]) / atof(argv[5]) * fidwidth[j];
				apod_par[j][1] = 0.0;
				apod_par[j][2] = 0.0;
			}

			else if (strcmp(argv[3], "gaussian") == 0) {
				apod_code[j] = LORENTZ_GAUSS;
				apod_par[j][0] = 0.0;
				apod_par[j][1] = sqrt(-log(atof(argv[5]))) * fidwidth[j] / (0.6 * M_PI * (atof(argv[4]) - 1.0) * (fidsize[j] - 1));
				apod_par[j][2] = atof(argv[4]);
			}

			else if (strcmp(argv[3], "gaussian_sw") == 0) {
				apod_code[j] = LORENTZ_GAUSS;
				apod_par[j][0] = atof(argv[4]) / atof(argv[6]) * fidwidth[j];
				apod_par[j][1] = 2.0 * apod_par[j][0];
				apod_par[j][2] = 2.0 * M_LN2 / (atof(argv[4]) * M_PI * fidsize[j] / fidwidth[j] * pow(atof(argv[5]), 2));
			}

			else if (strcmp(argv[3], "sinebell") == 0) {
				apod_code[j] = SINE_BELL;
				apod_par[j][0] = atof(argv[4]) / 180.0;
				apod_par[j][1] = 0.98;
				apod_par[j][2] = 1.0;
			}

			else if (strcmp(argv[3], "sinebell2") == 0) {
				apod_code[j] = SINE_BELL;
				apod_par[j][0] = atof(argv[4]) / 180.0;
				apod_par[j][1] = 0.98;
				apod_par[j][2] = 2.0;
			}

			/* phase */

			if (strcmp(argv[3], "phase") == 0) {
				phase[j][0] = atof(argv[4]);
				phase[j][1] = atof(argv[5]);
			}

			else if (strcmp(argv[3], "phase2") == 0) {
				phase[j][0] = atof(argv[4]) - (atoi(argv[6]) - 1) * atof(argv[5]);
				phase[j][1] = atof(argv[5]);
			}

		}

	}

	free_arg(0);

	fclose(par);

	int data_volume = get_data_volume();

	switch (dimension) {
	case 2:
		fprintf(stderr, "Axis Label | %8s  %8s\n", axisname[0], axisname[1]);
		fprintf(stderr, "Block Size | %8d  %8d\n", blocksize[0], blocksize[1]);
		fprintf(stderr, "Data Size  | %8d  %8d\n", datasize[0], datasize[1]);
		fprintf(stderr, "Obs. Freq. | %8.3f  %8.3f [MHz]\n", obsfreq[0], obsfreq[1]);
		fprintf(stderr, "Spec.Center| %8.2f  %8.2f [ppm]\n", spcenter[0], spcenter[1]);
		fprintf(stderr, "Spec.Width | %8.2f  %8.2f [Hz]\n", spwidth[0], spwidth[1]);

		/* READ ORIG. FREQ. */

		fprintf(stderr, "Orig.Freq. | %8.2f  %8.2f [Hz]\n", origfreq[0], origfreq[1]);

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
		if (size != sizeof(float) * data_volume) {
			fprintf(stderr, "Spectra file %s: Partially broken. (Actual=%d Expected=%d)\n", filename, (int) (size),
					(int) (sizeof(float)) * data_volume);

			return 1;
		}
		break;

	case 3:
		fprintf(stderr, "Axis Label | %8s  %8s  %8s\n", axisname[0], axisname[1], axisname[2]);
		fprintf(stderr, "Block Size | %8d  %8d  %8d\n", blocksize[0], blocksize[1], blocksize[2]);
		fprintf(stderr, "Data Size  | %8d  %8d  %8d\n", datasize[0], datasize[1], datasize[2]);
		fprintf(stderr, "Obs. Freq. | %8.3f  %8.3f  %8.3f [MHz]\n", obsfreq[0], obsfreq[1], obsfreq[2]);
		fprintf(stderr, "Spec.Center| %8.2f  %8.2f  %8.2f [ppm]\n", spcenter[0], spcenter[1], spcenter[2]);
		fprintf(stderr, "Spec.Width | %8.2f  %8.2f  %8.2f [Hz]\n", spwidth[0], spwidth[1], spwidth[2]);

		/* READ ORIG. FREQ. */

		fprintf(stderr, "Orig.Freq. | %8.2f  %8.2f  %8.2f [Hz]\n", origfreq[0], origfreq[1], origfreq[2]);

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
		if (size != sizeof(float) * data_volume) {
			fprintf(stderr, "Spectra file %s: Partially broken. (Actual=%d Expected=%d)\n", filename, (int) (size),
					(int) (sizeof(float)) * data_volume);

			return 1;
		}
		break;

	case 4:
		fprintf(stderr, "Axis Label | %8s  %8s  %8s  %8s\n", axisname[0], axisname[1], axisname[2], axisname[3]);
		fprintf(stderr, "Block Size | %8d  %8d  %8d  %8d\n", blocksize[0], blocksize[1], blocksize[2], blocksize[3]);
		fprintf(stderr, "Data Size  | %8d  %8d  %8d  %8d\n", datasize[0], datasize[1], datasize[2], datasize[3]);
		fprintf(stderr, "Obs. Freq. | %8.3f  %8.3f  %8.3f  %8.3f [MHz]\n", obsfreq[0], obsfreq[1], obsfreq[2], obsfreq[3]);
		fprintf(stderr, "Spec.Center| %8.2f  %8.2f  %8.2f  %8.2f [ppm]\n", spcenter[0], spcenter[1], spcenter[2],
				spcenter[3]);
		fprintf(stderr, "Spec.Width | %8.2f  %8.2f  %8.2f  %8.2f [Hz]\n", spwidth[0], spwidth[1], spwidth[2], spwidth[3]);

		/* READ ORIG. FREQ. */

		fprintf(stderr, "Orig.Freq. | %8.2f  %8.2f  %8.2f  %8.2f [Hz]\n", origfreq[0], origfreq[1], origfreq[2], origfreq[3]);

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
		if (size != sizeof(float) * data_volume) {
			fprintf(stderr, "Spectra file %s: Partially broken. (Actual=%d Expected=%d)\n", filename, (int) (size),
					(int) (sizeof(float)) * data_volume);

			return 1;
		}
		break;
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
		memcpy(header + j + k * MAXASSNAME, axisname[k], MAXASSNAME);
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
	fwrite2mem(header + 1768, get_indirect_planes());

	/* QUADFLAG */
	fwrite2mem(header + 424, 1.0);

	/* STATES */
	fwrite2mem(header + 1024, 2.0);

	fwrite2mem(header + 1028, 1.0);
	fwrite2mem(header + 1036, 1.0);
	fwrite2mem(header + 1044, 1.0);
	fwrite2mem(header + 1052, 1.0);

	return 0;

	escape:fclose(par);

	return 1;
}
