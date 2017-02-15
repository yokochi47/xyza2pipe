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

void checklabel(char string[])
{
	char _string[MAXASSNAME + 1];
	int l;

	char *cnvtbl[][2] = {
			{"1H", "H"},
			{"1D", "DH"},
			{"H1", "H"},
			{"D1", "DH"},
			{"H-ACQ", "AH"},
			{"H-IND", "DH"},
			{"ACQ", "AH"},
			{"IND", "DH"},
			{"A", "AH"},
			{"D", "DH"},
			{"13C", "C"},
			{"C13", "C"},
			{"15N", "N"},
			{"N15", "N"},
			{"31P", "P"},
			{"P31", "P"},
			{"", ""}
	};

	strncpy(_string, string, MAXASSNAME);

	strreplace(_string, " ", "");

	struppercpy(string, _string);

	l = 0;
	while (strcmp(cnvtbl[l][0], "") != 0) {

		if (strcmp(string, cnvtbl[l][0]) == 0) {
			strncpy(string, cnvtbl[l][1], MAXASSNAME);
			break;
		}

		l++;
	}

	l = strlen(string) + 1;

	memset(&(string[l]), 0, (MAXASSNAME - l) * sizeof(char));
}
