#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "msppriv.h"

int msp_verbose = 3;

/************************
 * Lightweight snprintf *
 ************************/

#ifndef kroundup64
#define kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))
#endif

static inline void str_enlarge(kstring_t *s, int l)
{
	if (s->l + l + 1 > s->m) {
		s->m = s->l + l + 1;
		kroundup64(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
}

static inline void str_copy(kstring_t *s, const char *st, const char *en)
{
	str_enlarge(s, en - st);
	memcpy(&s->s[s->l], st, en - st);
	s->l += en - st;
}

int64_t msp_sprintf_lite(kstring_t *s, const char *fmt, ...)
{
	char buf[32]; // for integer to string conversion
	const char *p, *q;
	int64_t len = 0;
	va_list ap;
	va_start(ap, fmt);
	for (q = p = fmt; *p; ++p) {
		if (*p == '%') {
			if (p > q) {
				len += p - q;
				if (s) str_copy(s, q, p);
			}
			++p;
			if (*p == 'd') {
				int c, i, l = 0;
				unsigned int x;
				c = va_arg(ap, int);
				x = c >= 0? c : -c;
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				if (c < 0) buf[l++] = '-';
				len += l;
				if (s) {
					str_enlarge(s, l);
					for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
				}
			} else if (*p == 'l' && *(p+1) == 'd') {
				int i, l = 0;
				long c;
				unsigned long x;
				c = va_arg(ap, long);
				x = c >= 0? c : -c;
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				if (c < 0) buf[l++] = '-';
				len += l;
				if (s) {
					str_enlarge(s, l);
					for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
				}
				++p;
			} else if (*p == 'u') {
				int i, l = 0;
				uint32_t x;
				x = va_arg(ap, uint32_t);
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				len += l;
				if (s) {
					str_enlarge(s, l);
					for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
				}
			} else if (*p == 's') {
				char *r = va_arg(ap, char*);
				int l;
				l = strlen(r);
				len += l;
				if (s) str_copy(s, r, r + l);
			} else if (*p == 'c') {
				++len;
				if (s) {
					str_enlarge(s, 1);
					s->s[s->l++] = va_arg(ap, int);
				}
			} else {
				fprintf(stderr, "ERROR: unrecognized type '%%%c'\n", *p);
				abort();
			}
			q = p + 1;
		}
	}
	if (p > q) {
		len += p - q;
		if (s) str_copy(s, q, p);
	}
	va_end(ap);
	if (s) s->s[s->l] = 0;
	return len;
}

/**************************
 * Timing and peak memory *
 **************************/

#if defined(WIN32) || defined(_WIN32)
#include <windows.h>

struct timezone
{
	__int32  tz_minuteswest; /* minutes W of Greenwich */
	int  tz_dsttime;     /* type of dst correction */
};

/*
 * gettimeofday.c
 *    Win32 gettimeofday() replacement
 *    taken from PostgreSQL, according to
 *    https://stackoverflow.com/questions/1676036/what-should-i-use-to-replace-gettimeofday-on-windows
 *
 * src/port/gettimeofday.c
 *
 * Copyright (c) 2003 SRA, Inc.
 * Copyright (c) 2003 SKC, Inc.
 *
 * Permission to use, copy, modify, and distribute this software and
 * its documentation for any purpose, without fee, and without a
 * written agreement is hereby granted, provided that the above
 * copyright notice and this paragraph and the following two
 * paragraphs appear in all copies.
 *
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE TO ANY PARTY FOR DIRECT,
 * INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
 * LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS
 * DOCUMENTATION, EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE AUTHOR SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE.  THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS
 * IS" BASIS, AND THE AUTHOR HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE,
 * SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
 */

/* FILETIME of Jan 1 1970 00:00:00. */
static const unsigned __int64 epoch = ((unsigned __int64) 116444736000000000ULL);

/*
 * timezone information is stored outside the kernel so tzp isn't used anymore.
 *
 * Note: this function is not for Win32 high precision timing purpose. See
 * elapsed_time().
 */
int gettimeofday(struct timeval * tp, struct timezone *tzp)
{
	FILETIME    file_time;
	SYSTEMTIME  system_time;
	ULARGE_INTEGER ularge;

	GetSystemTime(&system_time);
	SystemTimeToFileTime(&system_time, &file_time);
	ularge.LowPart = file_time.dwLowDateTime;
	ularge.HighPart = file_time.dwHighDateTime;

	tp->tv_sec = (long) ((ularge.QuadPart - epoch) / 10000000L);
	tp->tv_usec = (long) (system_time.wMilliseconds * 1000);

	return 0;
}

// taken from https://stackoverflow.com/questions/5272470/c-get-cpu-usage-on-linux-and-windows
double msp_cputime()
{
	HANDLE hProcess = GetCurrentProcess();
	FILETIME ftCreation, ftExit, ftKernel, ftUser;
	SYSTEMTIME stKernel;
	SYSTEMTIME stUser;

	GetProcessTimes(hProcess, &ftCreation, &ftExit, &ftKernel, &ftUser);
	FileTimeToSystemTime(&ftKernel, &stKernel);
	FileTimeToSystemTime(&ftUser, &stUser);

	double kernelModeTime = ((stKernel.wHour * 60.) + stKernel.wMinute * 60.) + stKernel.wSecond * 1. + stKernel.wMilliseconds / 1000.;
	double userModeTime = ((stUser.wHour * 60.) + stUser.wMinute * 60.) + stUser.wSecond * 1. + stUser.wMilliseconds / 1000.;

	return kernelModeTime + userModeTime;
}

long msp_peakrss(void) { return 0; }
#else
#include <sys/resource.h>
#include <sys/time.h>

double msp_cputime(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

long msp_peakrss(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
	return r.ru_maxrss * 1024;
#else
	return r.ru_maxrss;
#endif
}

#endif /* WIN32 || _WIN32 */

double msp_realtime(void)
{
	static double realtime0 = -1.0;
	struct timeval tp;
	double t;
	gettimeofday(&tp, NULL);
	t = tp.tv_sec + tp.tv_usec * 1e-6;
	if (realtime0 < 0.0) realtime0 = t;
	return t - realtime0;
}

double msp_percent_cpu(void)
{
	return (msp_cputime() + 1e-6) / (msp_realtime() + 1e-6);
}
