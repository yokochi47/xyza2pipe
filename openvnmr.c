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
#include "vendorpar.h"

typedef enum {
	HYPERCOMPLEX_DATA, COMPLEX_DATA, REAL_DATA
} enum_data_type;

int openvnmr2d(char monofile[], char pardir[], float **mat2d)
{
	struct stat _stat;
	FILE *fp;
	char string[MAXCHAR];
	int i, j, _i, blk_inc;
	int block_size, block_i, block_j, block_id;
	int offset_i, offset_j, offset, position;
	long size = 0;
	complex mix, ret;
	enum_data_type data_type = HYPERCOMPLEX_DATA;
	hypercomplex data_raw, data_ph;
	hypercomplex phi, _phase, _phase_x;
	datafilehead fheader;

	if ((fp = fopen(monofile, "r")) == NULL) {
		fprintf(stderr, "Spectra file %s: Couldn't open.\n", monofile);
		return 1;
	}

	stat(monofile, &_stat);
	size = _stat.st_size;

	fread(&fheader, sizeof(datafilehead), 1, fp);

	if ((fheader.status & S_COMPLEX) != 0) {
		if ((fheader.status & S_HYPERCOMPLEX) != 0)
			data_type = HYPERCOMPLEX_DATA;
		else
			data_type = COMPLEX_DATA;
	} else
		data_type = REAL_DATA;

	/* CHECK FILE SIZE */
	if (size != fheader.nblocks * fheader.tbytes * fheader.ntraces + sizeof(datafilehead)) {
		fprintf(stderr, "Spectra file %s: Partially broken. (Actual=%d Expected=%d)\n", monofile, (int) (size),
				(int) (fheader.nblocks * fheader.tbytes * fheader.ntraces + sizeof(datafilehead)));
		goto escape;
	}

	if (usrphase[0] == 0 || data_type != HYPERCOMPLEX_DATA) {
		get_varian_parameter(dimension, pardir, "rp", 0, string);
		phase[0][0] = atof(string);
		get_varian_parameter(dimension, pardir, "lp", 0, string);
		phase[0][1] = atof(string);
	}

	if (usrphase[1] == 0 || data_type != HYPERCOMPLEX_DATA) {
		get_varian_parameter(dimension, pardir, "rp1", 0, string);
		phase[1][0] = atof(string);
		get_varian_parameter(dimension, pardir, "lp1", 0, string);
		phase[1][1] = atof(string);
	}

	fwrite2mem(header + 436, phase[0][0]);
	fwrite2mem(header + 440, phase[0][1]);
	fwrite2mem(header + 980, phase[1][0]);
	fwrite2mem(header + 984, phase[1][1]);

#ifdef __DEBUG__

	fprintf(stderr, "phase settings: rp= %.1f lp= %.1f [deg]\n", phase[0][0], phase[0][1]);
	fprintf(stderr, "phase settings: rp1=%.1f lp1=%.1f [deg]\n\n", phase[1][0], phase[1][1]);

#endif

	block_size = fheader.ntraces * fheader.np;

	if (data_type == HYPERCOMPLEX_DATA) {
		phi.x.r = phi.y.r = 0.0;

		for (j = 0; j < datasize[1]; j++) {
			block_j = (int) (j / blocksize[1]);
			offset_j = j - block_j * blocksize[1];

			phi.y.i = -M_PI / 180.0 * (phase[1][0] + (float) (j) / (float) (datasize[1]) * phase[1][1]);

			_phase.y = Cexp(phi.y);

			for (i = 0; i < datasize[0]; i++) {
				_i = i * 4;
				block_i = (int) (_i / blocksize[0]);
				offset_i = _i - block_i * blocksize[0];

				phi.x.i = -M_PI / 180.0 * (phase[0][0] + (float) (i) / (float) (datasize[0]) * phase[0][1]);

				_phase.x = Cexp(phi.x);

				block_id = block_i + block_j * unitsize[0];

				offset = offset_i + offset_j * blocksize[0];

				position = offset + block_id * block_size;

				fseek(fp, (long) (position * sizeof(float) + sizeof(datafilehead)), SEEK_SET);

				fread(&data_raw, sizeof(char), sizeof(hypercomplex), fp);

				_phase_x.x = _phase.x;
				_phase_x.y = _phase.x;

				data_ph = HCmul(data_raw, _phase_x);

				mix.r = data_ph.x.r;
				mix.i = -data_ph.y.r;

				ret = Cmul(mix, _phase.y);

				mat2d[j][i] = ret.r;
			}
		}
	}

	else if (data_type == COMPLEX_DATA) {
		blk_inc = fheader.nblocks * fheader.tbytes * fheader.ntraces / (datasize[0] * datasize[1] * sizeof(complex));

		for (j = 0; j < datasize[1]; j++) {
			block_j = (int) (j / blocksize[1]);
			offset_j = j - block_j * blocksize[1];

			for (i = 0; i < datasize[0]; i++) {
				_i = i * 2;
				block_i = (int) (_i / blocksize[0]);
				offset_i = _i - block_i * blocksize[0];

				block_id = block_i + block_j * unitsize[0] * blk_inc;

				offset = offset_i + offset_j * blocksize[0];

				position = offset + block_id * block_size;

				fseek(fp, (long) (position * sizeof(float) + sizeof(datafilehead)), SEEK_SET);

				fread(&(mat2d[j][i]), sizeof(char), sizeof(float), fp);
			}
		}
	}

	else {

		for (j = 0; j < datasize[1]; j++) {
			block_j = (int) (j / blocksize[1]);
			offset_j = j - block_j * blocksize[1];

			for (i = 0; i < datasize[0]; i++) {
				block_i = (int) (i / blocksize[0]);
				offset_i = i - block_i * blocksize[0];

				block_id = block_i + block_j * unitsize[0];

				offset = offset_i + offset_j * blocksize[0];

				position = offset + block_id * block_size;

				fseek(fp, (long) (position * sizeof(float) + sizeof(datafilehead)), SEEK_SET);

				fread(&(mat2d[j][i]), sizeof(char), sizeof(float), fp);
			}
		}
	}

	fclose(fp);

	return 0;

	escape:fclose(fp);

	return 1;
}

int openvnmr3d(char monofile[], char pardir[], float ***mat3d)
{
	struct stat _stat;
	FILE *fp;
	int i, j, k, _i, blk_inc;
	int block_size, block_i, block_j, block_k, block_id;
	int offset_i, offset_j, offset_k, offset, position;
	long size = 0;
	enum_data_type data_type = HYPERCOMPLEX_DATA;

#ifdef __DEBUG__

	char string[MAXCHAR];

#endif

	datafilehead fheader;

	if ((fp = fopen(monofile, "r")) == NULL) {
		fprintf(stderr, "Spectra file %s: Couldn't open.\n", monofile);
		return 1;
	}

	stat(monofile, &_stat);
	size = _stat.st_size;

	fread(&fheader, sizeof(datafilehead), 1, fp);

	/* CHECK FILE SIZE */
	if (size != fheader.nblocks * fheader.tbytes * fheader.ntraces + sizeof(datafilehead)) {
		fprintf(stderr, "Spectra file %s: Partially broken. (Actual=%d Expected=%d)\n", monofile, (int) (size),
				(int) (fheader.nblocks * fheader.tbytes * fheader.ntraces + sizeof(datafilehead)));
		goto escape;
	}

	if ((fheader.status & S_COMPLEX) != 0) {
		if ((fheader.status & S_HYPERCOMPLEX) != 0) {
			data_type = HYPERCOMPLEX_DATA;
			fprintf(stderr, "hypercomplex data is not supported.\n");
			goto escape;
		} else
			data_type = COMPLEX_DATA;
	} else
		data_type = REAL_DATA;

#ifdef __DEBUG__

	get_varian_parameter(dimension, pardir, "rp", 0, string);
	phase[0][0] = atof(string);
	get_varian_parameter(dimension, pardir, "lp", 0, string);
	phase[0][1] = atof(string);

	get_varian_parameter(dimension, pardir, "rp1", 0, string);
	phase[1][0] = atof(string);
	get_varian_parameter(dimension, pardir, "lp1", 0, string);
	phase[1][1] = atof(string);

	get_varian_parameter(dimension, pardir, "rp2", 0, string);
	phase[2][0] = atof(string);
	get_varian_parameter(dimension, pardir, "lp2", 0, string);
	phase[2][1] = atof(string);

	fwrite2mem(header + 436, phase[0][0]);
	fwrite2mem(header + 440, phase[0][1]);
	fwrite2mem(header + 980, phase[1][0]);
	fwrite2mem(header + 984, phase[1][1]);
	fwrite2mem(header + 240, phase[2][0]);
	fwrite2mem(header + 244, phase[2][1]);

	fprintf(stderr, "phase settings: rp= %.1f lp= %.1f [deg]\n", phase[0][0], phase[0][1]);
	fprintf(stderr, "phase settings: rp1=%.1f lp1=%.1f [deg]\n", phase[1][0], phase[1][1]);
	fprintf(stderr, "phase settings: rp2=%.1f lp2=%.1f [deg]\n\n", phase[2][0], phase[2][1]);

#endif

	block_size = fheader.ntraces * fheader.np;

	if (data_type == COMPLEX_DATA) {
		blk_inc =
				fheader.nblocks * fheader.tbytes * fheader.ntraces / (datasize[0] * datasize[1] * datasize[2] * sizeof(complex));

		for (k = 0; k < datasize[2]; k++) {
			block_k = (int) (k / blocksize[2]);
			offset_k = k - block_k * blocksize[2];

			for (j = 0; j < datasize[1]; j++) {
				block_j = (int) (j / blocksize[1]);
				offset_j = j - block_j * blocksize[1];

				for (i = 0; i < datasize[0]; i++) {
					_i = i * 2;
					block_i = (int) (_i / blocksize[0]);
					offset_i = _i - block_i * blocksize[0];

					block_id = block_i + (block_j + block_k * unitsize[1]) * unitsize[0] * blk_inc;

					offset = offset_i + (offset_j + offset_k * blocksize[1]) * blocksize[0];

					position = offset + block_id * block_size;

					fseek(fp, (long) (position * sizeof(float) + sizeof(datafilehead)), SEEK_SET);

					fread(&(mat3d[k][j][i]), sizeof(char), sizeof(float), fp);
				}
			}
		}
	}

	else {

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

					position = offset + block_id * block_size;

					fseek(fp, (long) (position * sizeof(float) + sizeof(datafilehead)), SEEK_SET);

					fread(&(mat3d[k][j][i]), sizeof(char), sizeof(float), fp);
				}
			}
		}
	}

	fclose(fp);

	return 0;

	escape:fclose(fp);

	return 1;
}

int openvnmr4d(char monofile[], char pardir[], float ****mat4d)
{
	struct stat _stat;
	FILE *fp;
	int i, j, k, l, _i, blk_inc;
	int block_size, block_i, block_j, block_k, block_l, block_id;
	int offset_i, offset_j, offset_k, offset_l, offset, position;
	long size = 0;
	enum_data_type data_type = HYPERCOMPLEX_DATA;

#ifdef __DEBUG__

	char string[MAXCHAR];

#endif

	datafilehead fheader;

	if ((fp = fopen(monofile, "r")) == NULL) {
		fprintf(stderr, "Spectra file %s: Couldn't open.\n", monofile);
		return 1;
	}

	stat(monofile, &_stat);
	size = _stat.st_size;

	fread(&fheader, sizeof(datafilehead), 1, fp);

	/* CHECK FILE SIZE */
	if (size != fheader.nblocks * fheader.tbytes * fheader.ntraces + sizeof(datafilehead)) {
		fprintf(stderr, "Spectra file %s: Partially broken. (Actual=%d Expected=%d)\n", monofile, (int) (size),
				(int) (fheader.nblocks * fheader.tbytes * fheader.ntraces + sizeof(datafilehead)));
		goto escape;
	}

	if ((fheader.status & S_COMPLEX) != 0) {
		if ((fheader.status & S_HYPERCOMPLEX) != 0) {
			data_type = HYPERCOMPLEX_DATA;
			fprintf(stderr, "hypercomplex data is not supported.\n");
			goto escape;
		} else
			data_type = COMPLEX_DATA;
	} else
		data_type = REAL_DATA;

#ifdef __DEBUG__

	get_varian_parameter(dimension, pardir, "rp", 0, string);
	phase[0][0] = atof(string);
	get_varian_parameter(dimension, pardir, "lp", 0, string);
	phase[0][1] = atof(string);

	get_varian_parameter(dimension, pardir, "rp1", 0, string);
	phase[1][0] = atof(string);
	get_varian_parameter(dimension, pardir, "lp1", 0, string);
	phase[1][1] = atof(string);

	get_varian_parameter(dimension, pardir, "rp2", 0, string);
	phase[2][0] = atof(string);
	get_varian_parameter(dimension, pardir, "lp2", 0, string);
	phase[2][1] = atof(string);

	get_varian_parameter(dimension, pardir, "rp3", 0, string);
	phase[3][0] = atof(string);
	get_varian_parameter(dimension, pardir, "lp3", 0, string);
	phase[3][1] = atof(string);

	fwrite2mem(header + 436, phase[0][0]);
	fwrite2mem(header + 440, phase[0][1]);
	fwrite2mem(header + 980, phase[1][0]);
	fwrite2mem(header + 984, phase[1][1]);
	fwrite2mem(header + 240, phase[2][0]);
	fwrite2mem(header + 244, phase[2][1]);
	fwrite2mem(header + 248, phase[3][0]);
	fwrite2mem(header + 252, phase[3][1]);

	fprintf(stderr, "phase settings: rp= %.1f lp= %.1f [deg]\n", phase[0][0], phase[0][1]);
	fprintf(stderr, "phase settings: rp1=%.1f lp1=%.1f [deg]\n", phase[1][0], phase[1][1]);
	fprintf(stderr, "phase settings: rp2=%.1f lp2=%.1f [deg]\n", phase[2][0], phase[2][1]);
	fprintf(stderr, "phase settings: rp3=%.1f lp3=%.1f [deg]\n\n", phase[3][0], phase[3][1]);

#endif

	block_size = fheader.ntraces * fheader.np;

	if (data_type == COMPLEX_DATA) {
		blk_inc =
				fheader.nblocks * fheader.tbytes * fheader.ntraces / (datasize[0] * datasize[1] * datasize[2] * datasize[3] *
						sizeof(complex));

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
						_i = i * 2;
						block_i = (int) (_i / blocksize[0]);
						offset_i = _i - block_i * blocksize[0];

						block_id = block_i + (block_j + (block_k + block_l * unitsize[2]) * unitsize[1]) * unitsize[0] * blk_inc;

						offset = offset_i + (offset_j + (offset_k + offset_l * blocksize[2]) * blocksize[1]) * blocksize[0];

						position = offset + block_id * block_size;

						fseek(fp, (long) (position * sizeof(float) + sizeof(datafilehead)), SEEK_SET);

						fread(&(mat4d[l][k][j][i]), sizeof(char), sizeof(float), fp);
					}
				}
			}
		}
	}

	else {

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

						position = offset + block_id * block_size;

						fseek(fp, (long) (position * sizeof(float) + sizeof(datafilehead)), SEEK_SET);

						fread(&(mat4d[l][k][j][i]), sizeof(char), sizeof(float), fp);
					}
				}
			}
		}
	}

	fclose(fp);

	return 0;

	escape:fclose(fp);

	return 1;
}
