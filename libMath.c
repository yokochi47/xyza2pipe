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

#include "libMath.h"

/* complex */

complex Cadd(complex a, complex b)
{
	complex c;

	c.r = a.r + b.r;
	c.i = a.i + b.i;

	return c;
}

complex Csub(complex a, complex b)
{
	complex c;

	c.r = a.r - b.r;
	c.i = a.i - b.i;

	return c;
}

complex Cmul(complex a, complex b)
{
	complex c;

	c.r = a.r * b.r - a.i * b.i;
	c.i = a.i * b.r + a.r * b.i;

	return c;
}

complex Conj(complex z)
{
	complex c;

	c.r = z.r;
	c.i = -z.i;

	return c;
}

complex Cdiv(complex a, complex b)
{
	complex c;
	float r, den;

	if (fabs(b.r) >= fabs(b.i)) {
		r = b.i / b.r;
		den = b.r + r * b.i;
		c.r = (a.r + r * a.i) / den;
		c.i = (a.i - r * a.r) / den;
	}

	else {
		r = b.r / b.i;
		den = b.i + r * b.r;
		c.r = (a.r * r + a.i) / den;
		c.i = (a.i * r - a.r) / den;
	}

	return c;
}

float Cabs(complex z)
{
	float x, y, r, ans;

	x = fabs(z.r);
	y = fabs(z.i);

	if (x == 0.0)
		ans = y;

	else if (y == 0.0)
		ans = x;

	else if (x > y) {
		r = y / x;
		ans = x * sqrt(1.0 + r * r);
	}

	else {
		r = x / y;
		ans = y * sqrt(1.0 + r * r);
	}

	return ans;
}

complex Csqrt(complex z)
{
	complex c;
	float x, y, w, r;

	if (z.r == 0.0 && z.i == 0.0) {
		c.r = 0.0;
		c.i = 0.0;
	}

	else {
		x = fabs(z.r);
		y = fabs(z.i);

		if (x >= y) {
			r = y / x;
			w = sqrt(x) * sqrt(0.5 * (1.0 + sqrt(1.0 + r * r)));
		}

		else {
			r = x / y;
			w = sqrt(y) * sqrt(0.5 * (r + sqrt(1.0 + r * r)));
		}

		if (z.r >= 0.0) {
			c.r = w;
			c.i = z.i / (2.0 * w);
		}

		else {
			c.i = (z.i >= 0.0) ? w : -w;
			c.r = z.i / (2.0 * c.i);
		}
	}

	return c;
}

complex RCmul(float x, complex a)
{
	complex c;

	c.r = x * a.r;
	c.i = x * a.i;

	return c;
}

complex Cexp(complex z)
{
	complex c;
	float r;

	r = exp(z.r);

	c.r = r * cos(z.i);
	c.i = r * sin(z.i);

	return c;
}

/* double complex */

doublecomplex DCadd(doublecomplex a, doublecomplex b)
{
	doublecomplex c;

	c.r = a.r + b.r;
	c.i = a.i + b.i;

	return c;
}

doublecomplex DCsub(doublecomplex a, doublecomplex b)
{
	doublecomplex c;

	c.r = a.r - b.r;
	c.i = a.i - b.i;

	return c;
}

doublecomplex DCmul(doublecomplex a, doublecomplex b)
{
	doublecomplex c;

	c.r = a.r * b.r - a.i * b.i;
	c.i = a.i * b.r + a.r * b.i;

	return c;
}

doublecomplex DConj(doublecomplex z)
{
	doublecomplex c;

	c.r = z.r;
	c.i = -z.i;

	return c;
}

doublecomplex DCdiv(doublecomplex a, doublecomplex b)
{
	doublecomplex c;
	double r, den;

	if (fabs(b.r) >= fabs(b.i)) {
		r = b.i / b.r;
		den = b.r + r * b.i;
		c.r = (a.r + r * a.i) / den;
		c.i = (a.i - r * a.r) / den;
	}

	else {
		r = b.r / b.i;
		den = b.i + r * b.r;
		c.r = (a.r * r + a.i) / den;
		c.i = (a.i * r - a.r) / den;
	}

	return c;
}

double DCabs(doublecomplex z)
{
	double x, y, r, ans;

	x = fabs(z.r);
	y = fabs(z.i);

	if (x == 0.0)
		ans = y;

	else if (y == 0.0)
		ans = x;

	else if (x > y) {
		r = y / x;
		ans = x * sqrt(1.0 + r * r);
	}

	else {
		r = x / y;
		ans = y * sqrt(1.0 + r * r);
	}

	return ans;
}

doublecomplex DCsqrt(doublecomplex z)
{
	doublecomplex c;
	double x, y, w, r;

	if (z.r == 0.0 && z.i == 0.0) {
		c.r = 0.0;
		c.i = 0.0;
	}

	else {
		x = fabs(z.r);
		y = fabs(z.i);

		if (x >= y) {
			r = y / x;
			w = sqrt(x) * sqrt(0.5 * (1.0 + sqrt(1.0 + r * r)));
		}

		else {
			r = x / y;
			w = sqrt(y) * sqrt(0.5 * (r + sqrt(1.0 + r * r)));
		}

		if (z.r >= 0.0) {
			c.r = w;
			c.i = z.i / (2.0 * w);
		}

		else {
			c.i = (z.i >= 0.0) ? w : -w;
			c.r = z.i / (2.0 * c.i);
		}
	}

	return c;
}

doublecomplex RDCmul(double x, doublecomplex a)
{
	doublecomplex c;

	c.r = x * a.r;
	c.i = x * a.i;

	return c;
}

doublecomplex DCexp(doublecomplex z)
{
	doublecomplex c;
	double r;

	r = exp(z.r);

	c.r = r * cos(z.i);
	c.i = r * sin(z.i);

	return c;
}

/* hyper complex */

hypercomplex HCadd(hypercomplex a, hypercomplex b)
{
	hypercomplex c;

	c.x = Cadd(a.x, b.x);
	c.y = Cadd(a.y, b.y);

	return c;
}

hypercomplex HCsub(hypercomplex a, hypercomplex b)
{
	hypercomplex c;

	c.x = Csub(a.x, b.x);
	c.y = Csub(a.y, b.y);

	return c;
}

hypercomplex HCmul(hypercomplex a, hypercomplex b)
{
	hypercomplex c;

	c.x = Cmul(a.x, b.x);
	c.y = Cmul(a.y, b.y);

	return c;
}

hypercomplex HConj(hypercomplex z)
{
	hypercomplex c;

	c.x = Conj(z.x);
	c.y = Conj(z.y);

	return c;
}

hypercomplex HCtp(hypercomplex z)
{
	hypercomplex c;

	c.x = z.y;
	c.y = z.x;

	return c;
}
