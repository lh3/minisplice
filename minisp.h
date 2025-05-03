#ifndef MINISP_H
#define MINISP_H

#include <stdint.h>

typedef enum { MSP_FT_FASTX, MSP_FT_BED } msp_ft_t;

typedef struct {
	msp_ft_t type;
	void *fp;
	void *buf;
} msp_file_t;

typedef struct {
	int32_t n, m;
	char **a;
	void *h;
} msp_strmap_t;

typedef struct {
	int64_t st, en;
} msp_blk1_t;

typedef struct {
	int64_t st, en;
	int64_t st2, en2; // thickStart and thickEnd
	int32_t cid, score, strand;
	int32_t n_blk;
	char *ctg, *name;
	msp_blk1_t blk[];
} msp_bed1_t;

typedef struct {
	int64_t n, m;
	msp_bed1_t **a;
	msp_strmap_t *h;
} msp_bed_t;

extern int msp_verbose;

#ifdef __cplusplus
extern "C" {
#endif

// BED reader
msp_file_t *msp_bed_open(const char *fn);
void msp_file_close(msp_file_t *f);
msp_bed_t *msp_bed_read(const char *fn);
void msp_bed_sort(msp_bed_t *bed);
int msp_bed_read1(msp_file_t *fp, msp_bed1_t **b_);
void msp_bed_destroy(msp_bed_t *bed);

// strmap
msp_strmap_t *msp_strmap_init(void);
void msp_strmap_destroy(msp_strmap_t *m);
int32_t msp_strmap_add(msp_strmap_t *m, const char *s);
int32_t msp_strmap_get(const msp_strmap_t *m, const char *s);

#ifdef __cplusplus
}
#endif

#endif
