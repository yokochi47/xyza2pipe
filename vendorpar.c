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

#define FS_TAB		x09
#define FS_COMMA	','
#define FS_SPACE	' '

const float gH1 = 1.0;
const float gC13 = 0.251449530;
const float gN15 = 0.101329118;
const float gP31 = 0.404808636;

int get_varian_parsize(unsigned int dim, char pardir[], char parname[])
{
	FILE *fp = NULL;
	char long_buffer[MAXLONGNAME], procname[MAXLONGNAME], *argv[MAXVARS];
	int ret = 0;

	sprintf(procname, "%s/procpar", pardir);

	if ((fp = fopen(procname, "r")) == NULL) {

		sprintf(procname, "%s/procpar3d", pardir);

		if ((fp = fopen(procname, "r")) == NULL)
			return 0;
	}

	alloc_arg(0, MAXLONGNAME, argv);

	while (fgets(long_buffer, MAXLONGNAME, fp) != NULL) {
		line2arg(long_buffer, FS_SPACE, argv);

		if (strcmp(argv[0], parname) == 0) {
			fgets(long_buffer, MAXLONGNAME, fp);

			line2arg(long_buffer, FS_SPACE, argv);

			ret = atoi(argv[0]);

			break;
		}
	}

	fclose(fp);

	return ret;
}

int get_varian_parameter(unsigned int dim, char pardir[], char parname[], unsigned int col, char string[])
{
	FILE *fp = NULL;
	char long_buffer[MAXLONGNAME], procname[MAXLONGNAME], *argv[MAXVARS];
	int ret = 0;

	*string = '\0';

	sprintf(procname, "%s/procpar", pardir);

	if ((fp = fopen(procname, "r")) == NULL) {

		sprintf(procname, "%s/procpar3d", pardir);

		if ((fp = fopen(procname, "r")) == NULL)
			return 0;
	}

	alloc_arg(0, MAXLONGNAME, argv);

	while (fgets(long_buffer, MAXLONGNAME, fp) != NULL) {
		line2arg(long_buffer, FS_SPACE, argv);

		if (strcmp(argv[0], parname) == 0) {
			fgets(long_buffer, MAXLONGNAME, fp);

			line2arg(long_buffer, FS_SPACE, argv);

			ret = atoi(argv[0]);

			strcpy(string, argv[col + 1]);

			break;
		}
	}

	fclose(fp);

	strreplace(string, "\n", "");

	return ret;
}

int get_varian_dimension_from_file(unsigned int dim, char pardir[])
{
	FILE *fp;
	char procname[MAXLONGNAME];

	sprintf(procname, "%s/procpar", pardir);

	if ((fp = fopen(procname, "r")) == NULL) {

		sprintf(procname, "%s/procpar3d", pardir);

		if ((fp = fopen(procname, "r")) == NULL)
			return 0;
	}

	fclose(fp);

	switch (dim) {
	case 2:
		if (get_varian_parsize(dim, pardir, "phase") <= 1)
			return 1;
		else
			return 2;

		break;

	case 3:
		if (get_varian_parsize(dim, pardir, "phase") <= 1)
			return 1;
		else if (get_varian_parsize(dim, pardir, "phase2") <= 1)
			return 2;
		else
			return 3;

		break;

	case 4:
		if (get_varian_parsize(dim, pardir, "phase") <= 1)
			return 1;
		else if (get_varian_parsize(dim, pardir, "phase2") <= 1)
			return 2;
		else if (get_varian_parsize(dim, pardir, "phase3") <= 1)
			return 3;
		else
			return 4;

		break;
	}

	return 0;
}

int guess_varian_dimension_from_file(char pardir[])
{
	char long_buffer[MAXLONGNAME], *argv[MAXVARS];
	int dim;

	get_varian_parameter(0, pardir, "array", 0, long_buffer);

	strunquotecpy(long_buffer, long_buffer);

	alloc_arg(0, MAXLONGNAME, argv);

	dim = line2arg(long_buffer, FS_COMMA, argv) + 1;

	if (dim == 1) {
		get_varian_parameter(dim, pardir, "ni", 0, long_buffer);
		fprintf(stderr, "ni %s\n", long_buffer);
		if (atoi(long_buffer) > 1)
			dim++;
		get_varian_parameter(dim, pardir, "ni2", 0, long_buffer);
		fprintf(stderr, "ni2 %s\n", long_buffer);
		if (atoi(long_buffer) > 1)
			dim++;
		get_varian_parameter(dim, pardir, "ni3", 0, long_buffer);
		fprintf(stderr, "ni3 %s\n", long_buffer);
		if (atoi(long_buffer) > 1)
			dim++;
	}

	return dim;
}

static int get_bruker_parsize(unsigned int dim, char filename[], char parname[]);

int get_bruker_proc_parsize(unsigned int dim, char procdir[], unsigned int axis, char parname[])
{
	char procname[MAXLONGNAME];
	char *axis_digit[4] = { "", "2", "3", "4" };

	if (axis >= dim)
		return 0;

	sprintf(procname, "%s/proc%ss", procdir, axis_digit[axis]);

	return get_bruker_parsize(dim, procname, parname);
}

int get_bruker_acq_parsize(unsigned int dim, char acqdir[], unsigned int axis, char parname[])
{
	char acqname[MAXLONGNAME];
	char *axis_digit[4] = { "", "2", "3", "4" };

	if (axis >= dim)
		return 0;

	sprintf(acqname, "%s/acqu%ss", acqdir, axis_digit[axis]);

	return get_bruker_parsize(dim, acqname, parname);
}

int get_bruker_parsize(unsigned int dim, char filename[], char parname[])
{
	FILE *fp = NULL;
	char long_buffer[MAXLONGNAME];
	char _parname[MAXFILENAME];
	int l, m, ret = 0;

	if ((fp = fopen(filename, "r")) == NULL)
		return 0;

	sprintf(_parname, "##$%s= ", parname);
	l = strlen(_parname);

	while (fgets(long_buffer, MAXLONGNAME, fp) != NULL) {
		if (strncmp(long_buffer, _parname, l) == 0) {

			if (long_buffer[l] == '(' && (m = strsearch(long_buffer, ")")) != -1) {
				long_buffer[m--] = '\0';

				while (m > 0 && isdigit((int) (long_buffer[m])))
					m--;

				if (m > 0)
					ret = atoi(long_buffer + (++m)) + 1;
			}

			else
				ret = 1;

			break;
		}
	}

	fclose(fp);

	return ret;
}

static int get_bruker_parameter(unsigned int dim, char filename[], char parname[], unsigned int col, char string[]);

int get_bruker_proc_parameter(unsigned int dim, char procdir[], unsigned int axis, char parname[], unsigned int col,
		char string[])
{
	char procname[MAXLONGNAME];
	char *axis_digit[4] = { "", "2", "3", "4" };

	if (axis >= dim)
		return 0;

	sprintf(procname, "%s/proc%ss", procdir, axis_digit[axis]);

	return get_bruker_parameter(dim, procname, parname, col, string);
}

int get_bruker_acq_parameter(unsigned int dim, char acqdir[], unsigned int axis, char parname[], unsigned int col,
		char string[])
{
	char acqname[MAXLONGNAME];
	char *axis_digit[4] = { "", "2", "3", "4" };

	if (axis >= dim)
		return 0;

	sprintf(acqname, "%s/acqu%ss", acqdir, axis_digit[axis]);

	return get_bruker_parameter(dim, acqname, parname, col, string);
}

int get_bruker_parameter(unsigned int dim, char filename[], char parname[], unsigned int col, char string[])
{
	FILE *fp = NULL;
	char long_buffer[MAXLONGNAME], *argv[MAXVARS];
	char _parname[MAXFILENAME];
	int l, m, n, _n, ret = 0;

	*string = '\0';

	if ((fp = fopen(filename, "r")) == NULL)
		return 0;

	sprintf(_parname, "##$%s= ", parname);
	l = strlen(_parname);

	alloc_arg(0, MAXLONGNAME, argv);

	while (fgets(long_buffer, MAXLONGNAME, fp) != NULL) {
		if (strncmp(long_buffer, _parname, l) == 0) {

			if (long_buffer[l] == '(' && (m = strsearch(long_buffer, ")")) != -1) {
				long_buffer[m--] = '\0';

				while (m > 0 && isdigit((int) (long_buffer[m])))
					m--;

				if (m > 0) {
					ret = atoi(long_buffer + (++m)) + 1;

					n = _n = 0;
					while (col > n) {
						fgets(long_buffer, MAXLONGNAME, fp);
						_n = n;
						n += line2arg(long_buffer, FS_SPACE, argv);
					}

					strcpy(string, argv[col - _n]);
				}
			}

			else {
				ret = 1;

				if (long_buffer[l] != '<')
					strcpy(string, long_buffer + l);

				else {
					strreplacecpy(string, long_buffer + l, "<", "");
					strreplace(string, ">", "");
				}

				strunquotecpy(string, string);
			}

			break;
		}
	}

	fclose(fp);

	strreplace(string, "\n", "");

	return ret;
}

int get_bruker_dimension_from_acq_file(unsigned int dim, char acqdir[])
{
	FILE *fp;
	char acqname[MAXLONGNAME], long_buffer[MAXLONGNAME];

	sprintf(acqname, "%s/acqus", acqdir);

	if ((fp = fopen(acqname, "r")) == NULL)
		return 0;

	fclose(fp);

	switch (dim) {
	case 2:
		if (get_bruker_acq_parameter(dim, acqdir, 1, "TD", 0, long_buffer) == 0)
			return 0;

		if (atoi(long_buffer) <= 1)
			return 1;
		else
			return 2;

		break;

	case 3:
		if (get_bruker_acq_parameter(dim, acqdir, 1, "TD", 0, long_buffer) == 0)
			return 0;

		if (atoi(long_buffer) <= 1)
			return 1;

		else {
			if (get_bruker_acq_parameter(dim, acqdir, 2, "TD", 0, long_buffer) == 0)
				return 2;

			if (atoi(long_buffer) <= 1)
				return 2;
			else
				return 3;
		}

		break;

	case 4:
		if (get_bruker_acq_parameter(dim, acqdir, 1, "TD", 0, long_buffer) == 0)
			return 0;

		if (atoi(long_buffer) <= 1)
			return 1;

		else {
			if (get_bruker_acq_parameter(dim, acqdir, 2, "TD", 0, long_buffer) == 0)
				return 2;

			if (atoi(long_buffer) <= 1)
				return 2;

			else {
				if (get_bruker_acq_parameter(dim, acqdir, 3, "TD", 0, long_buffer) == 0)
					return 3;

				if (atoi(long_buffer) <= 1)
					return 3;
				else
					return 4;
			}
		}

		break;
	}

	return 0;
}

int guess_bruker_dimension_from_acq_file(char acqdir[])
{
	FILE *fp = NULL;
	char acqname[MAXLONGNAME];
	char *axis_digit[4] = { "", "2", "3", "4" };
	int axis;

	for (axis = 0; axis < 4; axis++) {
		sprintf(acqname, "%s/acqu%ss", acqdir, axis_digit[axis]);

		if ((fp = fopen(acqname, "r")) == NULL)
			break;
	}

	if (fp != NULL)
		fclose(fp);

	return axis;
}

int guess_bruker_dimension_from_proc_file(char procdir[])
{
	FILE *fp = NULL;
	char procname[MAXLONGNAME];
	char *axis_digit[4] = { "", "2", "3", "4" };
	int axis;

	for (axis = 0; axis < 4; axis++) {
		sprintf(procname, "%s/proc%ss", procdir, axis_digit[axis]);

		if ((fp = fopen(procname, "r")) == NULL)
			break;
	}

	if (fp != NULL)
		fclose(fp);

	return axis;
}
