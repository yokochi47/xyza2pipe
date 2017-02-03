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

/*  written by Masashi Yokochi */
/*  modified by Alexander Eletsky - Submatrix Tiling, Permutation patches */

int checkxeasy(char filename[])
{
	struct stat _stat;
	FILE *fp, *par;
	char parline[MAXCHAR], *argv[MAXVARS];
	char parfile[MAXLONGNAME];
	int fidsize[4] = { 0 };
	short apod_code[4] = { SINE_BELL };
	int j, k, column_len;
	long size = 0;
	float apod_par[4][3] = { {0.0} };

	int permut[4] = { 0 }, itmp;
	float ftmp;
	char ctmp[MAXASSNAME];		/* Added by A.E. */

	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Spectra file %s: Couldn't open.\n", filename);
		return 1;
	}

	stat(filename, &_stat);
	size = _stat.st_size;

	fclose(fp);

	strreplacecpy(parfile, filename, ".16", ".param");

	if ((par = fopen(parfile, "r")) == NULL) {
		fprintf(stderr, "Parameter file %s: Couldn't open.\n", parfile);
		return 1;
	}

	fprintf(stderr, "Reading %s ...\n", filename);

	/* READ PARAMETERS */
	alloc_arg(0, MAXFILENAME, argv);

	fgets(parline, MAXCHAR, par);

	if (strstr(parline, "Version") == NULL) {
		fprintf(stderr, "Parameter file %s: Not XEASY parameter file.\nCouldn't find description such as 'Version'.\n",
				parfile);
		goto escape;
	}

	fgets(parline, MAXCHAR, par);

	column_len = line2arg(parline, ' ', argv);

	dimension = atoi(argv[column_len - 1]);

	if (dimension < 2 || dimension > 4) {
		fprintf(stderr, "Spectra file: Not 2D/3D/4D experiment.\n");
		goto escape;
	}

	fgets(parline, MAXCHAR, par);

	column_len = line2arg(parline, ' ', argv);

	if (atoi(argv[column_len - 1]) != 16) {
		fprintf(stderr, "Unsupported file type. (Not 16 bit format)\n");
		goto escape;
	}

	/* READ OBS. FREQ. */
	for (j = 0; j < dimension; j++) {
		fgets(parline, MAXCHAR, par);
		column_len = line2arg(parline, ' ', argv);

		obsfreq[j] = atof(argv[column_len - 1]);
	}

	/* READ SPECTRAL WIDTH */
	for (j = 0; j < dimension; j++) {
		fgets(parline, MAXCHAR, par);
		column_len = line2arg(parline, ' ', argv);

		spwidth[j] = atof(argv[column_len - 1]) * obsfreq[j];
	}

	/* READ SPECTRAL CENTER */
	for (j = 0; j < dimension; j++) {
		fgets(parline, MAXCHAR, par);
		column_len = line2arg(parline, ' ', argv);

		spcenter[j] = atof(argv[column_len - 1]);
	}

	/* READ DATA SIZE */
	for (j = 0; j < dimension; j++) {
		fgets(parline, MAXCHAR, par);
		column_len = line2arg(parline, ' ', argv);

		datasize[j] = atoi(argv[column_len - 1]);
	}

	/* CALC SPECTRAL CENTER */
	for (j = 0; j < dimension; j++) {
		if (j == 0 && fabs(spcenter[j] - spwidth[j] / obsfreq[j] - 4.75) < 0.1)
			spcenter[j] = spcenter[j] - spwidth[j] / obsfreq[j];
		else
			spcenter[j] = spcenter[j] - spwidth[j] / obsfreq[j] / 2.0;	/* Modified by A.E. */
	}

	/* SUB MATRIX */
	for (j = 0; j < dimension; j++) {
		fgets(parline, MAXCHAR, par);
		column_len = line2arg(parline, ' ', argv);

		blocksize[j] = atoi(argv[column_len - 1]);	/* Added by A.E. */
	}

	/* PERMUTATION */
	for (j = 0; j < dimension; j++) {
		fgets(parline, MAXCHAR, par);
		column_len = line2arg(parline, ' ', argv);

		permut[j] = atoi(argv[column_len - 1]) - 1;	/* Added by A.E. */
	}

	/* FOLDING */
	for (j = 0; j < dimension; j++)
		fgets(parline, MAXCHAR, par);

	/* TYPE OF SPECTRUM */
	fgets(parline, MAXCHAR, par);

	/* READ AXIS NAME */
	for (j = 0; j < dimension; j++) {
		fgets(parline, MAXCHAR, par);
		column_len = line2arg(parline, ' ', argv);

		if (*axislabel[j] != 0)
			strncpy(axisname[j], axislabel[j], MAXASSNAME);
		else
			strncpy(axisname[j], argv[column_len - 1], MAXASSNAME);

		if (usrlabel == 0)
			checklabel(axisname[j]);
	}

	fclose(par);

	/* Reshuffle parameters according to permutation - added by A.E. */

	for (j = 0; j < dimension; j++) {

		while (permut[j] != j) {
			ftmp = obsfreq[j];
			obsfreq[j] = obsfreq[permut[j]];
			obsfreq[permut[j]] = ftmp;

			ftmp = spwidth[j];
			spwidth[j] = spwidth[permut[j]];
			spwidth[permut[j]] = ftmp;

			ftmp = spcenter[j];
			spcenter[j] = spcenter[permut[j]];
			spcenter[permut[j]] = ftmp;

			itmp = datasize[j];
			datasize[j] = datasize[permut[j]];
			datasize[permut[j]] = itmp;

			itmp = blocksize[j];
			blocksize[j] = blocksize[permut[j]];
			blocksize[permut[j]] = itmp;

			strncpy(ctmp, axisname[j], MAXASSNAME);
			strncpy(axisname[j], axisname[permut[j]], MAXASSNAME);
			strncpy(axisname[permut[j]], ctmp, MAXASSNAME);

			itmp = permut[j];
			permut[j] = permut[itmp];
			permut[itmp] = itmp;
		}
	}

	switch (dimension) {
	case 2:
		fprintf(stderr, "Axis Label | %8s  %8s\n", axisname[0], axisname[1]);

		fprintf(stderr, "Data Size  | %8d  %8d\n", datasize[0], datasize[1]);

		for (j = 0; j < dimension; j++)
			fidsize[j] = datasize[j] / 2;

		fprintf(stderr, "Block Size | %8d  %8d\n", blocksize[0], blocksize[1]);
		fprintf(stderr, "Dim. Order | %8d  %8d\n", permut[0] + 1, permut[1] + 1);	/* Added by A.E. */

		fprintf(stderr, "Obs. Freq. | %8.3f  %8.3f [MHz]\n", obsfreq[0], obsfreq[1]);
		fprintf(stderr, "Spec.Center| %8.2f  %8.2f [ppm]\n", spcenter[0], spcenter[1]);
		fprintf(stderr, "Spec.Width | %8.2f  %8.2f [Hz]\n", spwidth[0], spwidth[1]);

		/* READ ORIG. FREQ. */

		for (j = 0; j < dimension; j++) {
			origfreq[j] = spcenter[j] * obsfreq[j] - spwidth[j] / 2.0 + spwidth[j] / (float) (datasize[j]);

			if (*axisname[j] == 'H' && strncmp(axisname[j], "HA", MAXASSNAME) != 0 && *axisname[(j + 1) % dimension] == 'N'
					&& spcenter[j] >= 4.5 && spcenter[j] <= 5.0 && spwidth[j] / obsfreq[j] < 8.5)
				origfreq[j] += spwidth[j] / 2.0;
		}

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
		if (size != sizeof(short) * datasize[0] * datasize[1]) {
			fprintf(stderr, "Spectra file %s: Partially broken. (Actual=%d Expected=%d)\n", filename, (int) (size),
					(int) (sizeof(short)) * datasize[0] * datasize[1]);
			return 1;
		}
		break;

	case 3:
		fprintf(stderr, "Axis Label | %8s  %8s  %8s\n", axisname[0], axisname[1], axisname[2]);

		fprintf(stderr, "Data Size  | %8d  %8d  %8d\n", datasize[0], datasize[1], datasize[2]);

		for (j = 0; j < dimension; j++)
			fidsize[j] = datasize[j] / 2;

		fprintf(stderr, "Block Size | %8d  %8d  %8d\n", blocksize[0], blocksize[1], blocksize[2]);
		fprintf(stderr, "Dim. Order | %8d  %8d  %8d\n", permut[0] + 1, permut[1] + 1, permut[2] + 1);	/* Added by A.E. */

		fprintf(stderr, "Obs. Freq. | %8.3f  %8.3f  %8.3f [MHz]\n", obsfreq[0], obsfreq[1], obsfreq[2]);
		fprintf(stderr, "Spec.Center| %8.2f  %8.2f  %8.2f [ppm]\n", spcenter[0], spcenter[1], spcenter[2]);
		fprintf(stderr, "Spec.Width | %8.2f  %8.2f  %8.2f [Hz]\n", spwidth[0], spwidth[1], spwidth[2]);

		/* READ ORIG. FREQ. */

		for (j = 0; j < dimension; j++) {
			origfreq[j] = spcenter[j] * obsfreq[j] - spwidth[j] / 2.0 + spwidth[j] / (float) (datasize[j]);

			if (*axisname[j] == 'H' && strncmp(axisname[j], "HA", MAXASSNAME) != 0 && (*axisname[(j + 1) % dimension] == 'N'
					|| *axisname[(j + 2) % dimension] == 'N') && spcenter[j] >= 4.5 && spcenter[j] <= 5.0
					&& spwidth[j] / obsfreq[j] < 8.5)
				origfreq[j] += spwidth[j] / 2.0;
		}

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
		if (size != sizeof(short) * datasize[0] * datasize[1] * datasize[2]) {
			fprintf(stderr, "Spectra file %s: Partially broken. (Actual=%d Expected=%d)\n", filename, (int) (size),
					(int) (sizeof(short)) * datasize[0] * datasize[1] * datasize[2]);
			return 1;
		}
		break;

	case 4:
		fprintf(stderr, "Axis Label | %8s  %8s  %8s  %8s\n", axisname[0], axisname[1], axisname[2], axisname[3]);

		fprintf(stderr, "Data Size  | %8d  %8d  %8d  %8d\n", datasize[0], datasize[1], datasize[2], datasize[3]);

		for (j = 0; j < dimension; j++)
			fidsize[j] = datasize[j] / 2;

		fprintf(stderr, "Block Size | %8d  %8d  %8d  %8d\n", blocksize[0], blocksize[1], blocksize[2], blocksize[3]);
		fprintf(stderr, "Dim. Order | %8d  %8d  %8d  %8d\n", permut[0] + 1, permut[1] + 1, permut[2] + 1, permut[3] + 1);	/* Added by A.E. */

		fprintf(stderr, "Obs. Freq. | %8.3f  %8.3f  %8.3f  %8.3f [MHz]\n", obsfreq[0], obsfreq[1], obsfreq[2], obsfreq[3]);
		fprintf(stderr, "Spec.Center| %8.2f  %8.2f  %8.2f  %8.2f [ppm]\n", spcenter[0], spcenter[1], spcenter[2],
				spcenter[3]);
		fprintf(stderr, "Spec.Width | %8.2f  %8.2f  %8.2f  %8.2f [Hz]\n", spwidth[0], spwidth[1], spwidth[2], spwidth[3]);

		/* READ ORIG. FREQ. */

		for (j = 0; j < dimension; j++) {
			origfreq[j] = spcenter[j] * obsfreq[j] - spwidth[j] / 2.0 + spwidth[j] / (float) (datasize[j]);

			if (*axisname[j] == 'H' && strncmp(axisname[j], "HA", MAXASSNAME) != 0 && (*axisname[(j + 1) % dimension] == 'N'
					|| *axisname[(j + 2) % dimension] == 'N' || *axisname[(j + 3) % dimension] == 'N') && spcenter[j] >= 4.5
					&& spcenter[j] <= 5.0 && spwidth[j] / obsfreq[j] < 8.5)
				origfreq[j] += spwidth[j] / 2.0;
		}

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
		if (size != sizeof(short) * datasize[0] * datasize[1] * datasize[2] * datasize[3]) {
			fprintf(stderr, "Spectra file %s: Partially broken. (Actual=%d Expected=%d)\n", filename, (int) (size),
					(int) (sizeof(short)) * datasize[0] * datasize[1] * datasize[2] * datasize[3]);
			return 1;
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

	switch (dimension) {
	case 2:
		fwrite2mem(header + 1768, 1.0);
		break;
	case 3:
		fwrite2mem(header + 1768, (float) (datasize[2]));
		break;
	case 4:
		fwrite2mem(header + 1768, (float) (datasize[2] * datasize[3]));
		break;
	}

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
