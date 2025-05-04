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

typedef struct { uint64_t x, y, z; } msp192_t;
typedef struct { uint64_t x, y; } msp128_t;

typedef struct { int64_t n, m; msp192_t *a; } msp192a_t;
typedef struct { int64_t n, m; msp128_t *a; } msp128a_t;
typedef struct { int64_t n, m; uint64_t *a; } msp64a_t;

#ifdef __cplusplus
extern "C" {
#endif

char *msp_strdup(const char *src);
int64_t msp_sprintf_lite(kstring_t *s, const char *fmt, ...);

msp192_t *msp_bed_gen_negreg(const msp_bed_t *bed, int32_t cid, int64_t *nn);
void msp_bed_format(kstring_t *out, const msp_bed1_t *b);

double msp_cputime(void);
double msp_realtime(void);
double msp_percent_cpu(void);
long msp_peakrss(void);

#ifdef __cplusplus
}
#endif

#endif
