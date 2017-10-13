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
		"addvnmr2pipe -i/--in1 inTemplate1(2D/3D/4D) -j/--in2 inTemplate2(2D/3D/4D) > stdout\n", "Optional Arguments:\n",
		" -p/--pdir1   A directory that includes VNMR procpar file1. (default=\".\")\n",
		" -q/--pdir2   A directory that includes VNMR procpar file2. (default=\".\")\n",
		" -a/--add inTemplate1 + inTemplate2 (default)\n", " -s/--sub inTemplate1 - inTemplate2\n",
		" -m/--mul inTemplate1 * inTemplate2\n",
		" --c1   Scale factor for inTemplate1 (default=1.0)\n", " --c2   Scale factor for inTemplate2 (default=1.0)\n",
		" --xLAB X-Axis Label\n", " --yLAB Y-Axis Label\n", " --zLAB Z-Axis Label\n", " --aLAB A-Axis Label\n",
		" --xCAR X-Axis Center [ppm]\n", " --yCAR Y-Axis Center [ppm]\n", " --zCAR Z-Axis Center [ppm]\n",
		" --aCAR A-Axis Center [ppm]\n\n", ""
};

int main(int argc, char *argv[])
{
	char axis_option = 'x';
	char filename1[MAXLONGNAME] = { "" }, filename2[MAXLONGNAME] = { "" };
	char pardir1[MAXLONGNAME] = { "." }, pardir2[MAXLONGNAME] = { "." };
	char monofile1[MAXLONGNAME] = { "monofile1" }, monofile2[MAXLONGNAME] = { "monofile2" };
	int l, opt_idx = 0;
	float c1 = 1.0, c2 = 1.0;
	enum_combine_opr opr_code = COMBINE_ADD;

	static struct option options[] = {
			{"in1", required_argument, 0, 'i'},
			{"in2", required_argument, 0, 'j'},
			{"pdir1", required_argument, 0, 'p'},
			{"pdir2", required_argument, 0, 'q'},
			{"c1", required_argument, 0, 3},
			{"c2", required_argument, 0, 4},
			{"add", no_argument, 0, 'a'},
			{"sub", no_argument, 0, 's'},
			{"mul", no_argument, 0, 'm'},
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
		l = getopt_long(argc, argv, "i:j:p:q:asmd", options, &opt_idx);

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
		case 'p':
			if (optarg)
				strncpy(pardir1, optarg, MAXLONGNAME);
			break;
		case 'q':
			if (optarg)
				strncpy(pardir2, optarg, MAXLONGNAME);
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
		case 3:
			c1 = atof(optarg);
			break;
		case 4:
			c2 = atof(optarg);
			break;
		case 6:
			if (optarg) {
				usrlabel = 1;
				strncpy(axislabel[0], optarg, MAXAXISNAME);
			}
			break;
		case 7:
			if (optarg) {
				usrlabel = 1;
				strncpy(axislabel[1], optarg, MAXAXISNAME);
			}
			break;
		case 8:
			if (optarg) {
				usrlabel = 1;
				strncpy(axislabel[2], optarg, MAXAXISNAME);
			}
			break;
		case 9:
			if (optarg) {
				usrlabel = 1;
				strncpy(axislabel[3], optarg, MAXAXISNAME);
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

		return EXIT_FAILURE;
	}

	if (filename1[0] == 0 || filename2[0] == 0) {
		l = 0;
		while (strcmp(usage[l], "") != 0)
			fprintf(stderr, "%s", usage[l++]);

		return EXIT_FAILURE;
	}

	if (checkvnmr(filename1, pardir1, monofile1) != 0)
		return EXIT_FAILURE;

	_dimension = dimension;

	memcpy(&(datasize_orig[0]), &(datasize[0]), dimension * sizeof(int));

	if (checkvnmr(filename2, pardir2, monofile2) != 0)
		return EXIT_FAILURE;

	if (dimension != _dimension) {
		fprintf(stderr, "addvnmr2pipe error: unmatch dimension number.\n");
		return EXIT_FAILURE;
	}

	for (l = 0; l < dimension; l++) {
		if (datasize[l] != datasize_orig[l]) {
			fprintf(stderr, "addvnmr2pipe error: unmatch data size.\n");
			return EXIT_FAILURE;
		}
	}

	cnvhdr(axis_option, 'f');

	if (isatty(STDOUT_FILENO)) {
		fprintf(stderr, "addvnmr2pipe error: output to terminal.\n");
		return EXIT_FAILURE;
	}

	switch (dimension) {
	case 2:
		return pushaddvnmr2d(monofile1, monofile2, pardir1, pardir2, c1, c2, opr_code);
	case 3:
		return pushaddvnmr3d(monofile1, monofile2, pardir1, pardir2, c1, c2, opr_code);
	case 4:
		return pushaddvnmr4d(monofile1, monofile2, pardir1, pardir2, c1, c2, opr_code);
	default:
		fprintf(stderr, "addvnmr2pipe error: unsupported dimension number.\n");
		return EXIT_FAILURE;
	}

}
