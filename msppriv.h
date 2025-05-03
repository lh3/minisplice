#ifndef MSPPRIV_H
#define MSPPRIV_H

#include <stddef.h>
#include "minisp.h"

#define MSP_MALLOC(type, cnt)       ((type*)malloc((cnt) * sizeof(type)))
#define MSP_CALLOC(type, cnt)       ((type*)calloc((cnt), sizeof(type)))
#define MSP_REALLOC(type, ptr, cnt) ((type*)realloc((ptr), (cnt) * sizeof(type)))

#define MSP_GROW(type, ptr, __i, __m) do { \
		if ((__i) >= (__m)) { \
			(__m) = (__i) + 1; \
			(__m) += ((__m)>>1) + 16; \
			(ptr) = MSP_REALLOC(type, ptr, (__m)); \
		} \
	} while (0)

#ifndef KSTRING_T // same as kstring_t in kseq.h
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

char *msp_strdup(const char *src);
int64_t msp_sprintf_lite(kstring_t *s, const char *fmt, ...);

double msp_cputime(void);
double msp_realtime(void);
double msp_percent_cpu(void);
long msp_peakrss(void);

#ifdef __cplusplus
}
#endif

#endif
