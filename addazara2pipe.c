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

static const char *usage[] = {
		"addazara2pipe -i/--in1 inTemplate1(2D/3D/4D) -j/--in2 inTemplate2(2D/3D/4D) > stdout\n", "Optional Arguments:\n",
		" -a/--add inTemplate1 + inTemplate2 (default)\n", " -s/--sub inTemplate1 - inTemplate2\n",
		" -m/--mul inTemplate1 * inTemplate2\n", " -d/--div inTemplate1 / inTemplate2\n",
		" --c1   Scale factor for inTemplate1 (default=1.0)\n", " --c2   Scale factor for inTemplate2 (default=1.0)\n",
		" --xLAB X-Axis Label\n", " --yLAB Y-Axis Label\n", " --zLAB Z-Axis Label\n", " --aLAB A-Axis Label\n",
		" --xCAR X-Axis Center [ppm]\n", " --yCAR Y-Axis Center [ppm]\n", " --zCAR Z-Axis Center [ppm]\n",
		" --aCAR A-Axis Center [ppm]\n\n", ""
};

int main(int argc, char *argv[])
{
	char axis_option = 'x';
	char filename1[MAXLONGNAME] = { "" }, filename2[MAXLONGNAME] = {
			""};
	int l, opt_idx = 0;
	float c1 = 1.0, c2 = 1.0;
	enum_combine_opr opr_code = COMBINE_ADD;

	static struct option options[] = {
			{"in1", required_argument, 0, 'i'},
			{"in2", required_argument, 0, 'j'},
			{"c1", required_argument, 0, 3},
			{"c2", required_argument, 0, 4},
			{"add", no_argument, 0, 'a'},
			{"sub", no_argument, 0, 's'},
			{"mul", no_argument, 0, 'm'},
			{"div", no_argument, 0, 'd'},
			{"xLAB", required_argument, 0, 6},
			{"yLAB", required_argument, 0, 7},
			{"zLAB", required_argument, 0, 8},
			{"aLAB", required_argument, 0, 9},
			{"xCAR", required_argument, 0, 10},
			{"yCAR", required_argument, 0, 11},
			{"zCAR", required_argument, 0, 12},
			{"aCAR", required_argument, 0, 13},
			{0}
	};

	while (1) {
		l = getopt_long(argc, argv, "i:j:asmd", options, &opt_idx);

		if (l == -1)
			break;

		switch (l) {
		case 0:
			if (options[opt_idx].flag != 0)
				break;
			fprintf(stderr, "option %s", options[opt_idx].name);
			if (optarg)
				fprintf(stderr, " with arg %s", optarg);
			fputc('\n', stderr);
			break;
		case 'i':
			if (optarg)
				strncpy(filename1, optarg, MAXLONGNAME);
			break;
		case 'j':
			if (optarg)
				strncpy(filename2, optarg, MAXLONGNAME);
			break;
		case 'a':
			opr_code = COMBINE_ADD;
			break;
		case 's':
			opr_code = COMBINE_SUB;
			break;
		case 'm':
			opr_code = COMBINE_MUL;
			break;
		case 'd':
			opr_code = COMBINE_DIV;
			break;
		case 3:
			c1 = atof(optarg);
			break;
		case 4:
			c2 = atof(optarg);
			break;
		case 6:
			if (optarg) {
				usrlabel = 1;
				strncpy(axislabel[0], optarg, MAXASSNAME);
			}
			break;
		case 7:
			if (optarg) {
				usrlabel = 1;
				strncpy(axislabel[1], optarg, MAXASSNAME);
			}
			break;
		case 8:
			if (optarg) {
				usrlabel = 1;
				strncpy(axislabel[2], optarg, MAXASSNAME);
			}
			break;
		case 9:
			if (optarg) {
				usrlabel = 1;
				strncpy(axislabel[3], optarg, MAXASSNAME);
			}
			break;
		case 10:
			usrshift = 1;
			usrcenter[0] = atof(optarg);
			break;
		case 11:
			usrshift = 1;
			usrcenter[1] = atof(optarg);
			break;
		case 12:
			usrshift = 1;
			usrcenter[2] = atof(optarg);
			break;
		case 13:
			usrshift = 1;
			usrcenter[3] = atof(optarg);
			break;
		case '?':			/* getopt_long already printed an error message. */
			break;
		}
	}

	if (argc == 3) {
		strncpy(filename1, argv[1], MAXLONGNAME);
		strncpy(filename2, argv[2], MAXLONGNAME);
	}

	else if (optind < argc) {
		fprintf(stderr, "non-option ARGV-elements: ");
		while (optind < argc)
			fprintf(stderr, "%s ", argv[optind++]);
		fputc('\n', stderr);

		return 1;
	}

	if (filename1[0] == 0 || filename2[0] == 0) {
		l = 0;
		while (strcmp(usage[l], "") != 0)
			fprintf(stderr, "%s", usage[l++]);

		return 1;
	}

	if (checkazara(filename1) != 0)
		return 1;

	_dimension = dimension;

	memcpy(&(datasize_orig[0]), &(datasize[0]), dimension * sizeof(int));

	if (checkazara(filename2) != 0)
		return 1;

	if (dimension != _dimension) {
		fprintf(stderr, "addazara2pipe error: unmatch dimension number.\n");
		return 1;
	}

	for (l = 0; l < dimension; l++) {
		if (datasize[l] != datasize_orig[l]) {
			fprintf(stderr, "addazara2pipe error: unmatch data size.\n");
			return 1;
		}
	}

	cnvhdr(axis_option, 'f');

	if (isatty(STDOUT_FILENO)) {
		fprintf(stderr, "addazara2pipe error: output to terminal.\n");
		return 1;
	}

	switch (dimension) {
	case 2:
		return pushaddazara2d(filename1, filename2, c1, c2, opr_code);
		break;

	case 3:
		return pushaddazara3d(filename1, filename2, c1, c2, opr_code);
		break;

	case 4:
		return pushaddazara4d(filename1, filename2, c1, c2, opr_code);
		break;
	}

	return 0;
}
