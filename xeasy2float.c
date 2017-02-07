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

/* pipe2xeasy:some_spscan_methods.cc */
/*  written by Ralf W. Glaser */
/*  modified by Masashi Yokochi */

#include <math.h>

static char init = 0;
static float reverse[96];

static void xscale_init();
static float xexp10(float x);

static void xscale_init()
{
	int i;

	reverse[0] = 1.0;

	for (i = 1; i <= 47; i++)
		reverse[i] = (xexp10((float) (i - 1) / (float) 6.6438562) + (float) 0.5);
	for (i = 48; i <= 95; i++)
		reverse[i] = -(xexp10((float) (95 - i) / (float) 6.6438562) + (float) 0.5);

	init = 1;
}

static float xexp10(float x)
{
	return (float) (exp(x * log(10.0)));
}

float xeasy2float(unsigned char *x16)
{
	if (init == 0)
		xscale_init();

	return (float) (*(x16) + 615) * reverse[*(x16 + 1)] / 721.0;
}

void float2xeasy(float x, unsigned char *lo16)
{
	int hi, lo, sign;

	if (x < 0) {
		sign = -1;
		x = -x;
	} else
		sign = 1;

	if (x < 0.5) {
		*lo16 = (unsigned char) 0;
		*(lo16 + 1) = (unsigned char) 0;
		return;
	}

	if (init == 0)
		xscale_init();

	hi = (int) (log(x) * (float) 2.8853902 + 1.414);

	if (hi < 0)
		hi = 0;
	if (hi > 47)
		hi = 47;

	lo = (int) (712 * x / reverse[hi] - 615);

	if (lo < 0)
		lo = 0;
	if (lo > 255)
		lo = 255;

	if ((sign < 0) && hi)
		hi = 96 - hi;

	*lo16 = (unsigned char) lo;
	*(lo16 + 1) = (unsigned char) hi;
}
