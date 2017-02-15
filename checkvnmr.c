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
#include "vendorpar.h"

#define FT3DHDR_SIZE	260

int checkvnmr(char filename[], char pardir[], char monofile[])
{
	struct stat _stat;
	FILE *fp, *par, *mono, *fp3d[32];
	char _filename[MAXLONGNAME], parfile[MAXLONGNAME], axisorder[MAXASSNAME];
	char *buffer;
	char nucleus[MAXASSNAME], _obsfreq[MAXCHAR], _spwidth[MAXCHAR], _datasize[MAXCHAR], _fidsize[MAXCHAR];
	char _rfl[MAXCHAR], _rfp[MAXCHAR], temparature[MAXCHAR];
	int fidsize[4] = { 0 };
	short apod_code[4] = { SINE_BELL };
	int i, j, k, l, flag, nfiles = 1;
	long size = 0;
	float apod_par[4][3] = { {0.0} };
	float rfl[4], rfp[4];
	float refHelz = 0.0, h2oPPM;
	datafilehead fheader, _fheader;

#ifdef __DEBUG__

	datablockhead bheader;
	hypercmplxbhead hypercmplx;

#endif

	sprintf(parfile, "%s/procpar", pardir);

	if ((par = fopen(parfile, "r")) == NULL) {

		sprintf(parfile, "%s/procpar3d", pardir);

		if ((par = fopen(parfile, "r")) == NULL) {
			sprintf(parfile, "%s/procpar4d", pardir);

			if ((par = fopen(parfile, "r")) == NULL) {
				fprintf(stderr, "Parameter file %s/procpar: Couldn't open.\n", pardir);
				return 1;
			}
		}
	}

	fclose(par);

	if (dimension < 2 || dimension > 4) {
		dimension = guess_varian_dimension_from_file(pardir);

		if (dimension < 2 || dimension > 4) {
			fprintf(stderr, "Spectra file: Not 2D/3D/4D experiment.\n");
			return 1;
		}
	}

	switch (dimension) {
	case 2:
		strncpy(_filename, filename, MAXLONGNAME);
		break;
	case 3:
	case 4:
		sprintf(_filename, filename, 1);
		break;
	}

	if ((fp = fopen(_filename, "r")) == NULL) {
		fprintf(stderr, "Spectra file %s: Couldn't open.\n", _filename);
		return 1;
	}

	stat(_filename, &_stat);
	size = _stat.st_size;

	fread(&fheader, sizeof(datafilehead), 1, fp);

	byteswap = 0;

	if (fheader.tbytes != fheader.np * fheader.ebytes) {
		byteswap = 1;
		swapbyte(sizeof(int32_t), sizeof(int32_t) * 6, (char *) &fheader);
		swapbyte(sizeof(short), sizeof(short) * 2, (char *) &fheader + sizeof(int32_t) * 6);
		swapbyte(sizeof(int32_t), sizeof(int32_t), (char *) &fheader + sizeof(int32_t) * 6 + sizeof(short) * 2);
	}
#ifdef __DEBUG__

	fprintf(stderr, "number of blocks in file = %d\n", (int) fheader.nblocks);
	fprintf(stderr, "number of traces per block = %d\n", (int) fheader.ntraces);
	fprintf(stderr, "number of elements per trace = %d\n", (int) fheader.np);
	fprintf(stderr, "number of bytes per element = %d\n", (int) fheader.ebytes);
	fprintf(stderr, "number of bytes per trace = %d\n", (int) fheader.tbytes);
	fprintf(stderr, "number of bytes per block = %d\n", (int) fheader.bbytes);
	fprintf(stderr, "version_id = %d\n", fheader.vers_id);
	fprintf(stderr, "status = %d\n", fheader.status);
	fprintf(stderr, "number of block headers = %d\n", (int) fheader.nbheaders);
	fprintf(stderr, "file size = %d\n", (int) size);

#endif

	/* file header */

	if ((fheader.status & S_DATA) != 0) {

#ifdef __DEBUG__

		fprintf(stderr, "data file\n");

#endif

	} else {
		fprintf(stderr, "Spectra file: This is not data file.\n");
		goto escape;
	}

	if ((fheader.status & S_SPEC) != 0) {

#ifdef __DEBUG__

		fprintf(stderr, "spectra\n");

#endif

	}

	if ((fheader.status & S_FLOAT) == 0) {
		if ((fheader.status & S_32) != 0)
			fprintf(stderr, "32 bit integer data is not supported.\n");
		else
			fprintf(stderr, "16 bit integer data is not supported.\n");

		goto escape;
	}
#ifdef __DEBUG__

	fprintf(stderr, "floating point\n");

#endif

#ifdef __DEBUG__

	if ((fheader.status & S_COMPLEX) != 0) {
		if ((fheader.status & S_HYPERCOMPLEX) != 0)
			fprintf(stderr, "hyper complex data\n");
		else
			fprintf(stderr, "complex data\n");
	} else
		fprintf(stderr, "real data\n");

	if ((fheader.status & S_ACQPAR) != 0)
		fprintf(stderr, "acqpar\n");

#endif

	if ((fheader.status & S_SECND) != 0) {

#ifdef __DEBUG__

		fprintf(stderr, "second FT\n");

#endif

	}

#ifdef __DEBUG__

	if ((fheader.status & S_TRANSF) != 0)
		fprintf(stderr, "transposed\n");
	else
		fprintf(stderr, "regular\n");

	if ((fheader.status & S_NP) != 0)
		fprintf(stderr, "np dimension is active\n");

	if ((fheader.status & S_NF) != 0)
		fprintf(stderr, "nf dimension is active\n");

	if ((fheader.status & S_NI) != 0)
		fprintf(stderr, "ni dimension is active\n");

	if ((fheader.status & S_NI2) != 0)
		fprintf(stderr, "ni2 dimension is active\n");

	if ((fheader.status & S_NI3) != 0)
		fprintf(stderr, "ni3 dimension is active\n");

	fread(&bheader, sizeof(datablockhead), 1, fp);

	if (byteswap != 0) {
		swapbyte(sizeof(short), sizeof(short) * 4, (char *) &bheader);
		swapbyte(sizeof(int32_t), sizeof(int32_t), (char *) &bheader + sizeof(short) * 4);
		swapbyte(sizeof(float), sizeof(float) * 4, (char *) &bheader + sizeof(short) * 4 + sizeof(int32_t));
	}

	if (fheader.nbheaders > 1) {
		fread(&hypercmplx, sizeof(hypercmplxbhead), 1, fp);

		if (byteswap != 0) {
			swapbyte(sizeof(short), sizeof(short) * 4, (char *) &hypercmplx);
			swapbyte(sizeof(int32_t), sizeof(int32_t), (char *) &hypercmplx + sizeof(short) * 4);
			swapbyte(sizeof(float), sizeof(float) * 4, (char *) &hypercmplx + sizeof(short) * 4 + sizeof(int32_t));
		}
	}
#endif

	fclose(fp);

#ifdef __DEBUG__

	if ((fheader.status & S_NP) != 0) {
		if ((bheader.status & NP_CMPLX) != 0)
			fprintf(stderr, "np complex\n");
		else
			fprintf(stderr, "np real\n");

		if ((bheader.mode & NP_PHMODE) != 0)
			fprintf(stderr, "np ph mode\n");
		else if ((bheader.mode & NP_AVMODE) != 0)
			fprintf(stderr, "np av mode\n");
		if ((bheader.mode & NP_PWRMODE) != 0)
			fprintf(stderr, "np pwr mode\n");
	}

	if ((fheader.status & S_NF) != 0) {
		if ((bheader.status & NF_CMPLX) != 0)
			fprintf(stderr, "nf complex\n");
		else
			fprintf(stderr, "nf real\n");

		if ((bheader.mode & NF_PHMODE) != 0)
			fprintf(stderr, "nf ph mode\n");
		else if ((bheader.mode & NF_AVMODE) != 0)
			fprintf(stderr, "nf av mode\n");
		if ((bheader.mode & NF_PWRMODE) != 0)
			fprintf(stderr, "nf pwr mode\n");
	}

	if ((fheader.status & S_NI) != 0) {
		if ((bheader.status & NI_CMPLX) != 0)
			fprintf(stderr, "ni complex\n");
		else
			fprintf(stderr, "ni real\n");

		if ((bheader.mode & NI_PHMODE) != 0)
			fprintf(stderr, "ni ph mode\n");
		else if ((bheader.mode & NI_AVMODE) != 0)
			fprintf(stderr, "ni av mode\n");
		if ((bheader.mode & NI_PWRMODE) != 0)
			fprintf(stderr, "ni pwr mode\n");
	}

	if ((fheader.status & S_NI2) != 0) {
		if ((bheader.status & NI2_CMPLX) != 0)
			fprintf(stderr, "ni2 complex\n");
		else
			fprintf(stderr, "ni2 real\n");

		if ((bheader.mode & NI2_PHMODE) != 0)
			fprintf(stderr, "ni2 ph mode\n");
		else if ((bheader.mode & NI2_AVMODE) != 0)
			fprintf(stderr, "ni2 av mode\n");
		if ((bheader.mode & NI2_PWRMODE) != 0)
			fprintf(stderr, "ni2 pwr mode\n");
	}

	if ((fheader.status & S_NI3) != 0) {
		if ((bheader.status & NI3_CMPLX) != 0)
			fprintf(stderr, "ni3 complex\n");
		else
			fprintf(stderr, "ni3 real\n");

		if ((bheader.mode & NI3_PHMODE) != 0)
			fprintf(stderr, "ni3 ph mode\n");
		else if ((bheader.mode & NI3_AVMODE) != 0)
			fprintf(stderr, "ni3 av mode\n");
		if ((bheader.mode & NI3_PWRMODE) != 0)
			fprintf(stderr, "ni3 pwr mode\n");
	}

	if ((bheader.status & MORE_BLOCKS) != 0)
		fprintf(stderr, "more blocks present\n");
	else
		fprintf(stderr, "more blocks absent\n");

	if (fheader.nbheaders > 1) {
		if ((hypercmplx.status & U_HYPERCOMPLEX) != 0)
			fprintf(stderr, "hypercomplex block structure\n");

		fprintf(stderr, "lpval1 = %f rpval1 = %f\n", hypercmplx.lpval1, hypercmplx.rpval1);
	}
#endif

	fprintf(stderr, "Reading %s ...\n", _filename);

	/* READ PARAMETERS */

	get_varian_parameter(dimension, pardir, "axis", 0, axisorder);
	strunquotecpy(axisorder, axisorder);

	j = 0;

	get_varian_parameter(dimension, pardir, "tn", 0, nucleus);
	strunquotecpy(nucleus, nucleus);
	get_varian_parameter(dimension, pardir, "sfrq", 0, _obsfreq);
	get_varian_parameter(dimension, pardir, "sw", 0, _spwidth);
	get_varian_parameter(dimension, pardir, "np", 0, _fidsize);
	get_varian_parameter(dimension, pardir, "fn", 0, _datasize);
	get_varian_parameter(dimension, pardir, "rfl", 0, _rfl);
	get_varian_parameter(dimension, pardir, "rfp", 0, _rfp);

	if (*axislabel[j] != 0)
		strncpy(axisname[j], axislabel[j], MAXASSNAME);
	else
		strncpy(axisname[j], nucleus, MAXASSNAME);

	if (usrlabel == 0)
		checklabel(axisname[j]);

	obsfreq[j] = atof(_obsfreq);
	spwidth[j] = atof(_spwidth);
	fidsize[j] = atoi(_fidsize) / 2;
	if ((fheader.status & S_COMPLEX) == 0 || (fheader.status & S_HYPERCOMPLEX) != 0)
		datasize[j] = atoi(_datasize) / 2;
	else
		datasize[j] = atoi(_fidsize) / 2;
	rfl[j] = atof(_rfl);
	rfp[j] = atof(_rfp);

	for (j = 1; j < dimension; j++) {

		for (l = flag = 0; l < j; l++) {

			if (axisorder[j] == axisorder[l]) {

				strncpy(axisname[j], axisname[l], MAXASSNAME);
				obsfreq[j] = obsfreq[l];
				spwidth[j] = spwidth[l];

				flag = 1;

				break;
			}
		}

		if (flag != 0)
			continue;

		switch (axisorder[j]) {
		case '1':
		case 'd':
			get_varian_parameter(dimension, pardir, "dn", 0, nucleus);
			strunquotecpy(nucleus, nucleus);

			get_varian_parameter(dimension, pardir, "dfrq", 0, _obsfreq);

			get_varian_parameter(dimension, pardir, "sw1", 0, _spwidth);

			strncpy(axisname[j], nucleus, MAXASSNAME);
			if (usrlabel == 0)
				checklabel(axisname[j]);

			obsfreq[j] = atof(_obsfreq);
			spwidth[j] = atof(_spwidth);

			break;
		case '2':
			get_varian_parameter(dimension, pardir, "dn2", 0, nucleus);
			strunquotecpy(nucleus, nucleus);

			get_varian_parameter(dimension, pardir, "dfrq2", 0, _obsfreq);

			switch (dimension) {
			case 2:
				get_varian_parameter(dimension, pardir, "sw1", 0, _spwidth);

				strncpy(axisname[j], nucleus, MAXASSNAME);
				if (usrlabel == 0)
					checklabel(axisname[j]);

				obsfreq[j] = atof(_obsfreq);
				spwidth[j] = atof(_spwidth);

				break;
			case 3:
				get_varian_parameter(dimension, pardir, "sw2", 0, _spwidth);

				strncpy(axisname[j], nucleus, MAXASSNAME);
				if (usrlabel == 0)
					checklabel(axisname[j]);

				obsfreq[j] = atof(_obsfreq);
				spwidth[j] = atof(_spwidth);

				break;
			case 4:
				get_varian_parameter(dimension, pardir, "sw3", 0, _spwidth);

				strncpy(axisname[j], nucleus, MAXASSNAME);
				if (usrlabel == 0)
					checklabel(axisname[j]);

				obsfreq[j] = atof(_obsfreq);
				spwidth[j] = atof(_spwidth);

				break;
			}
			break;
			case '3':
				get_varian_parameter(dimension, pardir, "dn3", 0, nucleus);
				strunquotecpy(nucleus, nucleus);

				get_varian_parameter(dimension, pardir, "dfrq3", 0, _obsfreq);

				switch (dimension) {
				case 2:
					get_varian_parameter(dimension, pardir, "sw1", 0, _spwidth);

					strncpy(axisname[j], nucleus, MAXASSNAME);
					if (usrlabel == 0)
						checklabel(axisname[j]);

					obsfreq[j] = atof(_obsfreq);
					spwidth[j] = atof(_spwidth);

					break;
				case 3:
					get_varian_parameter(dimension, pardir, "sw2", 0, _spwidth);

					strncpy(axisname[j], nucleus, MAXASSNAME);
					if (usrlabel == 0)
						checklabel(axisname[j]);

					obsfreq[j] = atof(_obsfreq);
					spwidth[j] = atof(_spwidth);

					break;
				case 4:
					get_varian_parameter(dimension, pardir, "sw3", 0, _spwidth);

					strncpy(axisname[j], nucleus, MAXASSNAME);
					if (usrlabel == 0)
						checklabel(axisname[j]);

					obsfreq[j] = atof(_obsfreq);
					spwidth[j] = atof(_spwidth);

					break;
				}
				break;
		}
	}

	for (j = 1; j < dimension; j++) {

		switch (j) {
		case 1:
			get_varian_parameter(dimension, pardir, "ni", 0, _fidsize);
			get_varian_parameter(dimension, pardir, "fn1", 0, _datasize);
			get_varian_parameter(dimension, pardir, "rfl1", 0, _rfl);
			get_varian_parameter(dimension, pardir, "rfp1", 0, _rfp);
			break;
		case 2:
			get_varian_parameter(dimension, pardir, "ni2", 0, _fidsize);
			get_varian_parameter(dimension, pardir, "fn2", 0, _datasize);
			get_varian_parameter(dimension, pardir, "rfl2", 0, _rfl);
			get_varian_parameter(dimension, pardir, "rfp2", 0, _rfp);
			break;
		case 3:
			get_varian_parameter(dimension, pardir, "ni3", 0, _fidsize);
			get_varian_parameter(dimension, pardir, "fn3", 0, _datasize);
			get_varian_parameter(dimension, pardir, "rfl3", 0, _rfl);
			get_varian_parameter(dimension, pardir, "rfp3", 0, _rfp);
			break;
		}

		fidsize[j] = atoi(_fidsize);
		if ((fheader.status & S_COMPLEX) == 0 || (fheader.status & S_HYPERCOMPLEX) != 0)
			datasize[j] = atoi(_datasize) / 2;
		else
			datasize[j] = atoi(_fidsize);
		rfl[j] = atof(_rfl);
		rfp[j] = atof(_rfp);
	}

	if (adjh2o != 0 || adjcar != 0) {

		j = 0;

		if (adjh2o != 0) {
			get_varian_parameter(dimension, pardir, "temp", 0, temparature);
			h2oPPM = M_CONST * atof(temparature) + B_CONST;

			spcenter[j] = h2oPPM;
			rfl[j] = spwidth[j] / 2.0 + rfp[j] - (h2oPPM * obsfreq[j]);
			refHelz = obsfreq[j] * 1.0e+6 - (spwidth[j] / 2.0 - rfl[j] + rfp[j]);
		}

		else {
			spcenter[j] = (spwidth[j] / 2.0 - rfl[j] + rfp[j]) / obsfreq[j];
			refHelz = obsfreq[j] * 1.0e+6 - (spwidth[j] / 2.0 - rfl[j] + rfp[j]);
		}

		for (j = 1; j < dimension; j++) {
			strncpy(nucleus, axisname[j], MAXASSNAME);
			checklabel(nucleus);
			switch (nucleus[0]) {
			case 'H':
				rfl[j] = (gH1 * refHelz) - (obsfreq[j] * 1.0e+6) + (spwidth[j] / 2.0);
				break;
			case 'C':
				rfl[j] = (gC13 * refHelz) - (obsfreq[j] * 1.0e+6) + (spwidth[j] / 2.0);
				break;
			case 'N':
				rfl[j] = (gN15 * refHelz) - (obsfreq[j] * 1.0e+6) + (spwidth[j] / 2.0);
				break;
			case 'P':
				rfl[j] = (gP31 * refHelz) - (obsfreq[j] * 1.0e+6) + (spwidth[j] / 2.0);
				break;
			default:
				fprintf(stderr, "Spectra file %s: Unsupported nucleus '%c'\n", _filename, nucleus[0]);
				return 1;
			}
			spcenter[j] = (spwidth[j] / 2.0 - rfl[j]) / obsfreq[j];
		}
	}

	else {
		for (j = 0; j < dimension; j++)
			spcenter[j] = (spwidth[j] / 2.0 - rfl[j] + rfp[j]) / obsfreq[j];
	}

	switch (dimension) {
	case 2:
		fprintf(stderr, "Axis Label | %8s  %8s\n", axisname[0], axisname[1]);

		fprintf(stderr, "Data Size  | %8d  %8d\n", datasize[0], datasize[1]);

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
		if (size != fheader.nblocks * fheader.bbytes + sizeof(datafilehead)) {
			fprintf(stderr, "Spectra file %s: Partially broken. (Actual=%d Expected=%d)\n", _filename, (int) (size),
					(int) (fheader.nblocks * fheader.bbytes + sizeof(datafilehead)));
			return 1;
		}

		/* monofile: skip internal headers (endian resolved) */

		if ((mono = fopen(monofile, "w")) == NULL) {
			fprintf(stderr, "Permission denied to write %s.\n", monofile);
			exit(EXIT_FAILURE);
		}

		assert(buffer = (char *) malloc(fheader.bbytes * sizeof(char)));

		if ((fp = fopen(_filename, "r")) != NULL) {
			fread(&_fheader, sizeof(datafilehead), 1, fp);

			if (byteswap != 0) {
				swapbyte(sizeof(int32_t), sizeof(int32_t) * 6, (char *) &_fheader);
				swapbyte(sizeof(short), sizeof(short) * 2, (char *) &_fheader + sizeof(int32_t) * 6);
				swapbyte(sizeof(int32_t), sizeof(int32_t), (char *) &_fheader + sizeof(int32_t) * 6 + sizeof(short) * 2);
			}

			if (fheader.bbytes != _fheader.bbytes) {
				fprintf(stderr, "Spectra file %s: File size is not same.\n", _filename);
				exit(EXIT_FAILURE);
			}

			fwrite(&_fheader, sizeof(datafilehead), 1, mono);
		}

		for (j = 0; j < fheader.nblocks; j++) {
			fread(buffer, sizeof(char), fheader.bbytes, fp);

			if (byteswap != 0)
				swapbyte(sizeof(float), fheader.tbytes * fheader.ntraces,
						buffer + fheader.bbytes - fheader.tbytes * fheader.ntraces);

			fwrite(buffer + fheader.bbytes - fheader.tbytes * fheader.ntraces, sizeof(char), fheader.tbytes * fheader.ntraces,
					mono);
		}

		free(buffer);

		fclose(fp);

		fclose(mono);

		if (stat(monofile, &_stat) != 0) {
			fprintf(stderr, "Failed in writing %s.\n", monofile);
			exit(EXIT_FAILURE);
		}

		size = _stat.st_size;

		if (size != fheader.nblocks * fheader.tbytes * fheader.ntraces + sizeof(datafilehead)) {
			fprintf(stderr, "Spectra file %s: Partially broken. (Actual=%d Expected=%d)\n", monofile, (int) (size),
					(int) (fheader.nblocks * fheader.tbytes * fheader.ntraces + sizeof(datafilehead)));
			return 1;
		}
		break;

	case 3:
		fprintf(stderr, "Axis Label | %8s  %8s  %8s\n", axisname[0], axisname[1], axisname[2]);

		fprintf(stderr, "Data Size  | %8d  %8d  %8d\n", datasize[0], datasize[1], datasize[2]);

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
		nfiles = rint(fheader.nblocks * fheader.bbytes / (size - FT3DHDR_SIZE));

		if (nfiles < 1 || nfiles > 32) {
			fprintf(stderr, "Spectra file %s: Partially broken. (Actual=%d Expected=%d)\n", _filename, (int) (size),
					(int) (fheader.nblocks * fheader.bbytes + sizeof(datafilehead)));
			return 1;
		}

		/* monofile: merge data files & skip internal headers (endian resolved) */

		if ((mono = fopen(monofile, "w")) == NULL) {
			fprintf(stderr, "Permission denied to write %s.\n", monofile);
			exit(EXIT_FAILURE);
		}

		assert(buffer = (char *) malloc(fheader.bbytes * sizeof(char)));

		for (i = 0; i < nfiles; i++) {

			sprintf(_filename, filename, i + 1);

			if ((fp3d[i] = fopen(_filename, "r")) != NULL) {
				fread(&_fheader, sizeof(datafilehead), 1, fp3d[i]);

				if (byteswap != 0) {
					swapbyte(sizeof(int32_t), sizeof(int32_t) * 6, (char *) &_fheader);
					swapbyte(sizeof(short), sizeof(short) * 2, (char *) &_fheader + sizeof(int32_t) * 6);
					swapbyte(sizeof(int32_t), sizeof(int32_t), (char *) &_fheader + sizeof(int32_t) * 6 + sizeof(short) * 2);
				}

				if (fheader.bbytes != _fheader.bbytes) {
					fprintf(stderr, "Spectra file %s: File size is not same.\n", _filename);
					exit(EXIT_FAILURE);
				}

				if (i == 0)
					fwrite(&_fheader, sizeof(datafilehead), 1, mono);
			}
		}

		for (j = 0; j < fheader.nblocks / (nfiles); j++) {

			for (i = 0; i < nfiles; i++) {
				fread(buffer, sizeof(char), fheader.bbytes, fp3d[i]);

				if (byteswap != 0)
					swapbyte(sizeof(float), fheader.tbytes * fheader.ntraces,
							buffer + fheader.bbytes - fheader.tbytes * fheader.ntraces);

				fwrite(buffer + fheader.bbytes - fheader.tbytes * fheader.ntraces, sizeof(char), fheader.tbytes * fheader.ntraces,
						mono);
			}
		}

		free(buffer);

		for (i = 0; i < nfiles; i++)
			fclose(fp3d[i]);

		fclose(mono);

		if (stat(monofile, &_stat) != 0) {
			fprintf(stderr, "Failed in writing %s.\n", monofile);
			exit(EXIT_FAILURE);
		}

		size = _stat.st_size;

		if (size != fheader.nblocks * fheader.tbytes * fheader.ntraces + sizeof(datafilehead)) {
			fprintf(stderr, "Spectra file %s: Partially broken. (Actual=%d Expected=%d)\n", monofile, (int) (size),
					(int) (fheader.nblocks * fheader.tbytes * fheader.ntraces + sizeof(datafilehead)));
			return 1;
		}
		break;

	case 4:
		fprintf(stderr, "Axis Label | %8s  %8s  %8s  %8s\n", axisname[0], axisname[1], axisname[2], axisname[3]);

		fprintf(stderr, "Data Size  | %8d  %8d  %8d  %8d\n", datasize[0], datasize[1], datasize[2], datasize[3]);

		fprintf(stderr, "Obs. Freq. | %8.3f  %8.3f  %8.3f  %8.3f [MHz]\n", obsfreq[0], obsfreq[1], obsfreq[2], obsfreq[3]);
		fprintf(stderr, "Spec.Center| %8.2f  %8.2f  %8.2f  %8.2f [ppm]\n", spcenter[0], spcenter[1], spcenter[2],
				spcenter[3]);
		fprintf(stderr, "Spec.Width | %8.2f  %8.2f  %8.2f  %8.2f [Hz]\n", spwidth[0], spwidth[1], spwidth[2], spwidth[3]);

		/* READ ORIG. FREQ. */

		for (j = 0; j < dimension; j++) {
			origfreq[j] = spcenter[j] * obsfreq[j] - spwidth[j] / 2.0 + spwidth[j] / (float) (datasize[j]);

			if (*axisname[j] == 'H' && strncmp(axisname[j], "HA", MAXASSNAME) != 0 && (*axisname[(j + 1) % dimension] == 'N'
					|| *axisname[(j + 2) % dimension] == 'N') && spcenter[j] >= 4.5 && spcenter[j] <= 5.0
					&& spwidth[j] / obsfreq[j] < 8.5)
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
		nfiles = rint(fheader.nblocks * fheader.bbytes / (size - FT3DHDR_SIZE));

		if (nfiles < 1 || nfiles > 32) {
			fprintf(stderr, "Spectra file %s: Partially broken. (Actual=%d Expected=%d)\n", _filename, (int) (size),
					(int) (fheader.nblocks * fheader.bbytes + sizeof(datafilehead)));
			return 1;
		}

		/* monofile: merge data files & skip internal headers (endian resolved) */

		if ((mono = fopen(monofile, "w")) == NULL) {
			fprintf(stderr, "Permission denied to write %s.\n", monofile);
			exit(EXIT_FAILURE);
		}

		assert(buffer = (char *) malloc(fheader.bbytes * sizeof(char)));

		for (i = 0; i < nfiles; i++) {

			sprintf(_filename, filename, i + 1);

			if ((fp3d[i] = fopen(_filename, "r")) != NULL) {
				fread(&_fheader, sizeof(datafilehead), 1, fp3d[i]);

				if (byteswap != 0) {
					swapbyte(sizeof(int32_t), sizeof(int32_t) * 6, (char *) &_fheader);
					swapbyte(sizeof(short), sizeof(short) * 2, (char *) &_fheader + sizeof(int32_t) * 6);
					swapbyte(sizeof(int32_t), sizeof(int32_t), (char *) &_fheader + sizeof(int32_t) * 6 + sizeof(short) * 2);
				}

				if (fheader.bbytes != _fheader.bbytes) {
					fprintf(stderr, "Spectra file %s: File size is not same.\n", _filename);
					exit(EXIT_FAILURE);
				}

				if (i == 0)
					fwrite(&_fheader, sizeof(datafilehead), 1, mono);
			}
		}

		for (j = 0; j < fheader.nblocks / (nfiles); j++) {

			for (i = 0; i < nfiles; i++) {
				fread(buffer, sizeof(char), fheader.bbytes, fp3d[i]);

				if (byteswap != 0)
					swapbyte(sizeof(float), fheader.tbytes * fheader.ntraces,
							buffer + fheader.bbytes - fheader.tbytes * fheader.ntraces);

				fwrite(buffer + fheader.bbytes - fheader.tbytes * fheader.ntraces, sizeof(char), fheader.tbytes * fheader.ntraces,
						mono);
			}
		}

		free(buffer);

		for (i = 0; i < nfiles; i++)
			fclose(fp3d[i]);

		fclose(mono);

		if (stat(monofile, &_stat) != 0) {
			fprintf(stderr, "Failed in writing %s.\n", monofile);
			exit(EXIT_FAILURE);
		}

		size = _stat.st_size;

		if (size != fheader.nblocks * fheader.tbytes * fheader.ntraces + sizeof(datafilehead)) {
			fprintf(stderr, "Spectra file %s: Partially broken. (Actual=%d Expected=%d)\n", monofile, (int) (size),
					(int) (fheader.nblocks * fheader.tbytes * fheader.ntraces + sizeof(datafilehead)));
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

	switch (dimension) {
	case 2:
		blocksize[0] = fheader.ntraces * 4;
		blocksize[1] = fheader.np;
		break;
	case 4:
		blocksize[3] = 1;
		/* no break */
	case 3:
		blocksize[1] = fheader.ntraces;
		blocksize[0] = fheader.np / nfiles;
		blocksize[2] = nfiles;
		break;
	}

	j = 64;
	for (k = 0; k < dimension; k++) {
		memcpy(header + j + k * MAXASSNAME, axisname[k], MAXASSNAME);
		unitsize[k] = datasize[k] / blocksize[k];
		if (unitsize[k] == 0)
			unitsize[k] = 1;
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

	escape:fclose(fp);

	return 1;
}
