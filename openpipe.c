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

int openpipe2d(float **mat2d)
{
	fread(&(mat2d[0][0]), sizeof(float), datasize[0] * datasize[1], stdin);

	return 0;
}
