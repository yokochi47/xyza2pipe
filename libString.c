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

#include "libString.h"

#define MAXTHRDS	8

void path2uscore(char *path)
{
	const char code[] = { " /\\.:;=()[]{}\0" };
	char c;
	int j, l = 0;

	while ((c = path[l]) != 0) {
		j = 0;
		while (code[j] != '\0') {
			if (code[j] == c) {
				path[l] = '_';
				break;
			}
			j++;
		}
		l++;
	}
}

int column_of(const char *str)
{
	int l = 0, m = 0, col = 1;

	while (str[l] != 0) {
		if (str[l] == '\n' || str[l] == '\r') {
			if (m > col)
				col = m;
			m = 0;
		}
		l++;
		m++;
	}

	return col;
}

int row_of(const char *str)
{
	int l = 0, row = 1;

	while (str[l] != 0) {
		if (str[l] == '\n' || str[l] == '\r')
			row++;
		l++;
	}

	return row;
}

int quote_of(char *str)
{
	int l = 0, q1 = 0, q2 = 0;

	while (str[l] != 0) {
		if (str[l] == '\'')
			q1++;
		if (str[l] == '"')
			q2++;
		l++;
	}

	if (q1 != 0 && q2 != 0) {
		l = 0;
		while (str[l] != 0) {
			if (str[l] == '"')
				str[l] = '\'';
			l++;
		}
	}

	return q1 != 0 ? '"' : '\'';
}

int is_integer(const char *str)
{
	int l = 0, ret = 0;

	while ((str[l] == ' ' || str[l] == '"' || str[l] == '\'') && str[l] != 0)
		l++;

	while (str[l] != 0) {
		if (ret == 0) {
			if (str[l] == '+' || str[l] == '-') {
			} else if (isdigit((int) (str[l])))
				ret = 1;
			else
				return 0;
			l++;
		} else {
			if (str[l] == ' ' || isdigit((int) (str[l])))
				l++;
			else if (str[l] == '"' || str[l] == '\'' || str[l] == '\t' || str[l] == '\n' || str[l] == '\r')
				return ret;
			else
				return 0;
		}
	}

	return ret;
}

int strsearch(const char *haystack, const char *needle)
{
	char *p = strstr(haystack, needle);

	if (p != NULL)
		return (int) (p - haystack);

	return -1;
}

int strcasesearch(const char *haystack, const char *needle)
{
	char *p = strcasestr(haystack, needle);

	if (p != NULL)
		return (int) (p - haystack);

	return -1;
}

static int _strsearch(const char *haystack, const char *needle, size_t n);

int strmatchcount(const char *haystack, const char *needle)
{
	size_t l, n = 0;
	int ret = 0;

	while ((l = _strsearch(haystack, needle, n)) != -1) {
		n = l + 1;
		ret++;
	}

	return ret;
}

static int _strsearch(const char *haystack, const char *needle, size_t n)
{
	char *p = strstr(haystack + n, needle);

	if (p != NULL)
		return (int) (p - haystack);

	return -1;
}

static int _strcasesearch(const char *haystack, const char *needle, size_t n);

int strcasematchcount(const char *haystack, const char *needle)
{
	size_t l, n = 0;
	int ret = 0;

	while ((l = _strcasesearch(haystack, needle, n)) != -1) {
		n = l + 1;
		ret++;
	}

	return ret;
}

static int _strcasesearch(const char *haystack, const char *needle, size_t n)
{
	char *p = strcasestr(haystack + n, needle);

	if (p != NULL)
		return (int) (p - haystack);

	return -1;
}

char *struppercpy(char *dst, const char *src)
{
	int l = 0;

	while (src[l] != 0) {
		dst[l] = toupper(src[l]);
		l++;
	}

	dst[l] = '\0';

	return dst;
}

char *strlowercpy(char *dst, const char *src)
{
	int l = 0;

	while (src[l] != 0) {
		dst[l] = tolower(src[l]);
		l++;
	}

	dst[l] = '\0';

	return dst;
}

char *strclasscpy(char *dst, const char *src)
{
	const char code[] = { " .,:;*/|_\0" };
	char c;
	int j, l = 0;

	strcpy(dst, src);

	while ((c = dst[l]) != 0) {
		j = 0;
		while (code[j] != '\0') {
			if (code[j] == c) {
				dst[l] = '\0';
				return dst;
			}
			j++;
		}
		l++;
	}

	return dst;
}

char *strunquotecpy(char *dst, const char *src)
{
	size_t l = strlen(src);
	char q = src[0];

	if ((q != '\'' && q != '"') || l < 1 || q != src[l - 1]) {
		strcpy(dst, src);
		return dst;
	}

	strncpy(dst, src + 1, l - 2);

	dst[l - 2] = '\0';

	return dst;
}

char *strseqoffcpy(char *dst, const char *src)
{
	int l = 1;

	strcpy(dst, src);

	while (!isdigit((int) (dst[0])) && dst[0] != '-' && src[l] != 0)
		strcpy(dst, src + (l++));

	return dst;
}

char *dirnamecpy(char *dir, const char *path)
{
	char *tmp;

	if (path == NULL)
		return dir;

	assert(tmp = (char *) alloca((strlen(path) + 1) * sizeof(char)));

	strcpy(tmp, path);

	strcpy(dir, dirname(tmp));

	return dir;
}

char *basenamecpy(char *base, const char *path)
{
	char *tmp;

	if (path == NULL)
		return base;

	assert(tmp = (char *) alloca((strlen(path) + 1) * sizeof(char)));

	strcpy(tmp, path);

	strcpy(base, basename(tmp));

	return base;
}

char *strspacefillcpy(char *dst, const char *src, size_t len)
{
	size_t l = strlen(src);

	strcpy(dst, src);

	while (l < len)
		dst[l++] = ' ';

	*(dst + len) = '\0';

	return dst;
}

char *strspacefill(char *src, size_t len)
{
	size_t l = strlen(src);

	while (l < len)
		src[l++] = ' ';

	*(src + len) = '\0';

	return src;
}

char *strcenteringcpy(char *dst, const char *src, size_t len)
{
	int s = (len - strlen(src)) / 2;

	strspacefillcpy(dst, src, s > 0 ? len - s : len);

	if (s > 0) {
		memmove(dst + s, dst, len - s);
		memset(&(dst[0]), ' ', s * sizeof(char));

		*(dst + len) = '\0';
	}

	return dst;
}

char *strquarteringcpy(char *dst, const char *src, size_t len)
{
	int s = (len - strlen(src)) / 4;

	strspacefillcpy(dst, src, s > 0 ? len - s : len);

	if (s > 0) {
		memmove(dst + s, dst, len - s);
		memset(&(dst[0]), ' ', s * sizeof(char));

		*(dst + len) = '\0';
	}

	return dst;
}

char *strreplacecpy(char *dst, const char *src, const char *key, const char *trl)
{
	const size_t ini_len = strlen(src);
	const size_t key_len = strlen(key);
	char nomute = 1, *p, *q;
	char *tmp;

	strcpy(p = dst, src);

	while ((p = strstr(p, key))) {
		p += key_len;
		nomute = 0;
	}

	if (nomute)
		return dst;

	dst[0] = '\0';

	assert(tmp = (char *) alloca((ini_len + 1) * sizeof(char)));

	strcpy(tmp, src);

	for (p = tmp; (q = strsplit(p, key)) != NULL; p = q) {
		strcat(dst, p);
		strcat(dst, trl);
	}

	strcat(dst, p);

	return dst;
}

char *strreplace(char *src, const char *key, const char *trl)
{
	const size_t ini_len = strlen(src);
	const size_t key_len = strlen(key);
	char nomute = 1, *p, *q;
	char *tmp;

	assert(tmp = (char *) alloca((ini_len + 1) * sizeof(char)));

	strcpy(p = tmp, src);

	while ((p = strstr(p, key))) {
		p += key_len;
		nomute = 0;
	}

	if (nomute)
		return src;

	src[0] = '\0';

	for (p = tmp; (q = strsplit(p, key)) != NULL; p = q) {
		strcat(src, p);
		strcat(src, trl);
	}

	strcat(src, p);

	return src;
}

char *strsplit(const char *src, const char *key)
{
	char *p;

	if ((p = strstr(src, key)) == NULL)
		return NULL;

	*p = '\0';

	return p + strlen(key);
}

void swrite2bin(FILE * fp, int16_t a)
{
	fwrite((void *) &a, sizeof(int16_t), 1, fp);
}

void iwrite2bin(FILE * fp, int32_t a)
{
	fwrite((void *) &a, sizeof(int32_t), 1, fp);
}

void lwrite2bin(FILE * fp, long a)
{
	fwrite((void *) &a, sizeof(long), 1, fp);
}

void fwrite2bin(FILE * fp, float a)
{
	fwrite((void *) &a, sizeof(float), 1, fp);
}

void spwrite2bin(FILE * fp, int16_t * a, size_t n)
{
	fwrite((void *) a, sizeof(int16_t), n, fp);
}

void ipwrite2bin(FILE * fp, int32_t * a, size_t n)
{
	fwrite((void *) a, sizeof(int32_t), n, fp);
}

void lpwrite2bin(FILE * fp, long *a, size_t n)
{
	fwrite((void *) a, sizeof(long), n, fp);
}

void fpwrite2bin(FILE * fp, float *a, size_t n)
{
	fwrite((void *) a, sizeof(float), n, fp);
}

void swrite2mem(void *mem, int16_t a)
{
	memcpy(mem, (void *) &a, sizeof(int16_t));
}

void iwrite2mem(void *mem, int32_t a)
{
	memcpy(mem, (void *) &a, sizeof(int32_t));
}

void lwrite2mem(void *mem, long a)
{
	memcpy(mem, (void *) &a, sizeof(long));
}

void fwrite2mem(void *mem, float a)
{
	memcpy(mem, (void *) &a, sizeof(float));
}

void swrite2bin_swap(FILE * fp, int16_t a, int swap)
{
	int i, l;
	char *_a;

	if (swap == 0) {
		swrite2bin(fp, a);
		return;
	}

	l = sizeof(int16_t) - 1;

	for (i = 0; i < sizeof(int16_t); i++) {
		_a = (char *) (&a) + l - i;
		fwrite((void *) (_a), 1, 1, fp);
	}
}

void iwrite2bin_swap(FILE * fp, int32_t a, int swap)
{
	int i, l;
	char *_a;

	if (swap == 0) {
		iwrite2bin(fp, a);
		return;
	}

	l = sizeof(int32_t) - 1;

	for (i = 0; i < sizeof(int32_t); i++) {
		_a = (char *) (&a) + l - i;
		fwrite((void *) (_a), 1, 1, fp);
	}
}

void lwrite2bin_swap(FILE * fp, long a, int swap)
{
	int i, l;
	char *_a;

	if (swap == 0) {
		lwrite2bin(fp, a);
		return;
	}

	l = sizeof(long) - 1;

	for (i = 0; i < sizeof(long); i++) {
		_a = (char *) (&a) + l - i;
		fwrite((void *) (_a), 1, 1, fp);
	}
}

void fwrite2bin_swap(FILE * fp, float a, int swap)
{
	int i, l;
	char *_a;

	if (swap == 0) {
		fwrite2bin(fp, a);
		return;
	}

	l = sizeof(float) - 1;

	for (i = 0; i < sizeof(float); i++) {
		_a = (char *) (&a) + l - i;
		fwrite((void *) (_a), 1, 1, fp);
	}
}

void spwrite2bin_swap(FILE * fp, int16_t * a, size_t n, int swap)
{
	int i, j, l;
	char *_a;

	if (swap == 0) {
		spwrite2bin(fp, a, n);
		return;
	}

	for (j = 1; j <= n; j++) {

		l = sizeof(int16_t) * j - 1;

		for (i = 0; i < sizeof(int16_t); i++) {
			_a = (char *) (a) + l - i;
			fwrite((void *) (_a), 1, 1, fp);
		}
	}
}

void ipwrite2bin_swap(FILE * fp, int32_t * a, size_t n, int swap)
{
	int i, j, l;
	char *_a;

	if (swap == 0) {
		ipwrite2bin(fp, a, n);
		return;
	}

	for (j = 1; j <= n; j++) {

		l = sizeof(int32_t) * j - 1;

		for (i = 0; i < sizeof(int32_t); i++) {
			_a = (char *) (a) + l - i;
			fwrite((void *) (_a), 1, 1, fp);
		}
	}
}

void lpwrite2bin_swap(FILE * fp, long *a, size_t n, int swap)
{
	int i, j, l;
	char *_a;

	if (swap == 0) {
		lpwrite2bin(fp, a, n);
		return;
	}

	for (j = 1; j <= n; j++) {

		l = sizeof(long) * j - 1;

		for (i = 0; i < sizeof(long); i++) {
			_a = (char *) (a) + l - i;
			fwrite((void *) (_a), 1, 1, fp);
		}
	}
}

void fpwrite2bin_swap(FILE * fp, float *a, size_t n, int swap)
{
	int i, j, l;
	char *_a;

	if (swap == 0) {
		fpwrite2bin(fp, a, n);
		return;
	}

	for (j = 1; j <= n; j++) {

		l = sizeof(float) * j - 1;

		for (i = 0; i < sizeof(float); i++) {
			_a = (char *) (a) + l - i;
			fwrite((void *) (_a), 1, 1, fp);
		}
	}
}

void swrite2mem_swap(void *mem, int16_t a, int swap)
{
	int i, l;
	char *_mem, *_a;

	if (swap == 0) {
		swrite2mem(mem, a);
		return;
	}

	l = sizeof(int16_t) - 1;

	for (i = 0; i < sizeof(int16_t); i++) {
		_mem = (char *) mem + i;
		_a = (char *) (&a) + l - i;
		memcpy((void *) (_mem), (void *) (_a), 1);
	}
}

void iwrite2mem_swap(void *mem, int32_t a, int swap)
{
	int i, l;
	char *_mem, *_a;

	if (swap == 0) {
		iwrite2mem(mem, a);
		return;
	}

	l = sizeof(int32_t) - 1;

	for (i = 0; i < sizeof(int32_t); i++) {
		_mem = (char *) mem + i;
		_a = (char *) (&a) + l - i;
		memcpy((void *) (_mem), (void *) (_a), 1);
	}
}

void lwrite2mem_swap(void *mem, long a, int swap)
{
	int i, l;
	char *_mem, *_a;

	if (swap == 0) {
		lwrite2mem(mem, a);
		return;
	}

	l = sizeof(long) - 1;

	for (i = 0; i < sizeof(long); i++) {
		_mem = (char *) mem + i;
		_a = (char *) (&a) + l - i;
		memcpy((void *) (_mem), (void *) (_a), 1);
	}
}

void fwrite2mem_swap(void *mem, float a, int swap)
{
	int i, l;
	char *_mem, *_a;

	if (swap == 0) {
		fwrite2mem(mem, a);
		return;
	}

	l = sizeof(float) - 1;

	for (i = 0; i < sizeof(float); i++) {
		_mem = (char *) mem + i;
		_a = (char *) (&a) + l - i;
		memcpy((void *) (_mem), (void *) (_a), 1);
	}
}

int is_big_endian()
{
	union {
		char c[4];
		int32_t i32;
	} val;

	val.c[0] = 0;
	val.c[1] = 0;
	val.c[2] = 0;
	val.c[3] = 1;

	if (val.i32 == 1)
		return 1;
	else
		return 0;
}

int is_little_endian()
{
	union {
		char c[4];
		int32_t i32;
	} val;

	val.c[0] = 1;
	val.c[1] = 0;
	val.c[2] = 0;
	val.c[3] = 0;

	if (val.i32 == 1)
		return 1;
	else
		return 0;
}

int is_big_endian_float(char *mem, float real)
{
	union {
		char c[sizeof(float)];
		float f;
	} val;
	int i;

	val.f = real;

	for (i = 0; i < sizeof(float); i++) {
		if (mem[i] != val.c[i])
			break;
	}

	if (i == sizeof(float))
		return is_big_endian();

	for (i = 0; i < sizeof(float); i++) {
		if (mem[i] != val.c[sizeof(float) - 1 - i])
			break;
	}

	if (i == sizeof(float))
		return is_little_endian();

	return -1;
}

int is_little_endian_float(char *mem, float real)
{
	union {
		char c[sizeof(float)];
		float f;
	} val;
	int i;

	val.f = real;

	for (i = 0; i < sizeof(float); i++) {
		if (mem[i] != val.c[i])
			break;
	}

	if (i == sizeof(float))
		return is_little_endian();

	for (i = 0; i < sizeof(float); i++) {
		if (mem[i] != val.c[sizeof(float) - 1 - i])
			break;
	}

	if (i == sizeof(float))
		return is_big_endian();

	return -1;
}

void swapbyte(unsigned int unit_size, unsigned int mem_size, char *mem)
{
	char *tmp;
	int i, j, l = unit_size - 1;

	if (unit_size < 2 || unit_size > mem_size) {
		/*
  if (unit_size <= 0 || unit_size > 1)
   fprintf(stderr, "Invalid parameters in swapbyte().\n");
		 */
		return;
	}

	assert(tmp = (char *) alloca(unit_size * sizeof(char)));

	for (i = 0; i < mem_size; i += unit_size) {
		memcpy(tmp, mem + i, unit_size);

		for (j = 0; j < unit_size; j++)
			mem[i + j] = tmp[l - j];
	}
}

char *parse2str(char *p_char, const char *str, unsigned int start, unsigned int stop)
{
	int l = stop - start + 1;

	memset(&(p_char[0]), 0, (l + 1) * sizeof(char));

	if (l <= 0 || strlen(str) < start)
		return p_char;

	strncpy(p_char, str + start - 1, l);

	strreplace(p_char, " ", "");

	return p_char;
}

int *parse2int(int *p_int, const char *str, unsigned int start, unsigned int stop)
{
	size_t l = strlen(str);
	char *tmp;

	assert(tmp = (char *) alloca((l + 1) * sizeof(char)));

	parse2str(tmp, str, start, stop);

	*p_int = atoi(tmp);

	return p_int;
}

float *parse2float(float *p_float, const char *str, unsigned int start, unsigned int stop)
{
	size_t l = strlen(str);
	char *tmp;

	assert(tmp = (char *) alloca((l + 1) * sizeof(char)));

	parse2str(tmp, str, start, stop);

	*p_float = atof(tmp);

	return p_float;
}

static unsigned int _arg_size = 0;
static unsigned int *_argc = NULL;
static char **_argv = NULL;

int line2arg(const char *str, unsigned char fs_code, char *argv[])
{
	char *p;

	int j = 0, l, argc = 0, size = 0, quote, thrd_id;

	if (str == NULL || argv == NULL)
		return 0;

	for (thrd_id = 0; thrd_id < _arg_size; thrd_id++) {
		if (*argv == _argv[thrd_id]) {
			size = _argc[thrd_id] * MAXVARS;
			break;
		}
	}

	if (thrd_id == _arg_size)
		return 0;

	while (j < size && str[j] != 0) {

		while (str[j] == fs_code)
			j++;

		if (str[j] == '\n' || str[j] == '\r')
			break;

		p = argv[argc];

		p[l = 0] = '\0';

		if (str[j] == '\'' || str[j] == '"') {
			quote = str[j];

			p[l++] = str[j++];

			while (str[j] != quote && str[j] != 0)
				p[l++] = str[j++];

			p[l++] = str[j++];
		}

		else {
			while (str[j] != fs_code && str[j] != 0)
				p[l++] = str[j++];
		}

		if (p[l - 1] == '\n' || p[l - 1] == '\r') {
			p[--l] = '\0';

			++argc;

			break;
		}

		else {
			p[l] = '\0';

			if (++argc >= MAXVARS)
				break;
		}

		if (str[j++] == 0)
			break;
	}

	if (argc < MAXVARS)
		*argv[argc] = '\0';

	return argc;
}

int line4arg(const char *str, unsigned char fs_code, char *argv[], unsigned char start_code, unsigned char stop_code)
{
	char *p;

	int j = 0, l, argc = 0, size = 0, quote, thrd_id;

	if (str == NULL || argv == NULL)
		return 0;

	for (thrd_id = 0; thrd_id < _arg_size; thrd_id++) {
		if (*argv == _argv[thrd_id]) {
			size = _argc[thrd_id] * MAXVARS;
			break;
		}
	}

	if (thrd_id == _arg_size)
		return 0;

	while (j < size && str[j] != 0) {

		while (str[j] == fs_code)
			j++;

		if (str[j] == '\n' || str[j] == '\r')
			break;

		p = argv[argc];

		p[l = 0] = '\0';

		if (str[j] == start_code) {
			j++;

			while (str[j] != stop_code && str[j] != 0)
				p[l++] = str[j++];
		}

		else if (str[j] == '\'' || str[j] == '"') {
			quote = str[j];

			p[l++] = str[j++];

			while (str[j] != quote && str[j] != 0)
				p[l++] = str[j++];

			p[l++] = str[j++];
		}

		else {
			while (str[j] != fs_code && str[j] != 0)
				p[l++] = str[j++];
		}

		if (p[l - 1] == '\n' || p[l - 1] == '\r') {
			p[--l] = '\0';

			++argc;

			break;
		}

		else {
			p[l] = '\0';

			if (++argc >= MAXVARS)
				break;
		}

		if (str[j++] == 0)
			break;
	}

	if (argc < MAXVARS)
		*argv[argc] = '\0';

	return argc;
}

/* NOTE THAT FOLLOWING FUNCTIONS ARE THREAD-UNSAFE FUNCTIONS */

void alloc_arg(unsigned int thrd_id, int size, char *argv[])
{
	int i, _arg_size_;

	if (thrd_id >= _arg_size) {
		_arg_size_ = (int) (thrd_id / MAXTHRDS) * MAXTHRDS;
		if (thrd_id >= _arg_size_)
			_arg_size_ += MAXTHRDS;

		if (_argc == NULL) {
			assert(_argc = (unsigned int *) calloc((_arg_size = _arg_size_), sizeof(unsigned int)));
			assert(_argv = (char **) calloc(_arg_size, sizeof(char *)));
		}

		else {
			assert(_argc = (unsigned int *) realloc(_argc, _arg_size_ * sizeof(unsigned int)));
			assert(_argv = (char **) realloc(_argv, _arg_size_ * sizeof(char *)));
			memset(&(_argc[_arg_size]), 0, (_arg_size_ - _arg_size) * sizeof(unsigned int));
			memset(&(_argv[_arg_size]), 0, (_arg_size_ - _arg_size) * sizeof(char *));
			_arg_size = _arg_size_;
		}
	}

	if (size <= _argc[thrd_id] && _argv[thrd_id] != NULL) {
		size = _argc[thrd_id];

		for (i = 0; i < MAXVARS; i++)
			argv[i] = _argv[thrd_id] + i * size;

		return;
	}

	free_arg(thrd_id);

	assert(_argv[thrd_id] = (char *) malloc(MAXVARS * size * sizeof(char)));

	for (i = 0; i < MAXVARS; i++)
		argv[i] = _argv[thrd_id] + i * size;

	_argc[thrd_id] = size;
}

void free_arg(unsigned int thrd_id)
{
	if (_argv[thrd_id] != NULL)
		free(_argv[thrd_id]);

	_argv[thrd_id] = NULL;

	_argc[thrd_id] = 0;
}

static int __dirname_size__ = 0;
static char *__dirname__ = NULL;

char *_dirname_(const char *path)
{
	if (path == NULL)
		return NULL;

	if (__dirname__ == NULL)
		assert(__dirname__ = (char *) malloc((__dirname_size__ =
				(strlen(path) > MAXCHAR ? strlen(path) : MAXCHAR) + 1) * sizeof(char)));
	else if (strlen(path) + 1 > __dirname_size__)
		assert(__dirname__ = (char *) realloc(__dirname__, __dirname_size__ = strlen(path) + 1));

	strcpy(__dirname__, path);

	return dirname(__dirname__);
}

static int __basename_size__ = 0;
static char *__basename__ = NULL;

char *_basename_(const char *path)
{
	if (path == NULL)
		return NULL;

	if (__basename__ == NULL)
		assert(__basename__ = (char *) malloc((__basename_size__ =
				(strlen(path) > MAXCHAR ? strlen(path) : MAXCHAR) + 1) * sizeof(char)));
	else if (strlen(path) + 1 > __basename_size__)
		assert(__basename__ = (char *) realloc(__basename__, __basename_size__ = strlen(path) + 1));

	strcpy(__basename__, path);

	return basename(__basename__);
}
