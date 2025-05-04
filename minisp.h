#ifndef MINISP_H
#define MINISP_H

#include <stdio.h>
#include <stdint.h>

typedef enum { MSP_FT_FASTX, MSP_FT_LINE } msp_ft_t;

typedef const char *msp_cstr_t;

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
	int64_t n, off;
} msp_bedctg_t;

typedef struct {
	int64_t n, m;
	msp_bed1_t **a;
	msp_bedctg_t *c;
	msp_strmap_t *h;
} msp_bed_t;

typedef struct {
	int32_t cid;
	uint64_t x; // pos<<3 | rev<<2 | acceptor<<1 | neg
	uint8_t *seq;
} msp_tdata1_t;

typedef struct {
	int64_t len;
	int64_t n[2], m[2];
	msp_tdata1_t *a[2];
} msp_tdata_t;

extern int msp_verbose;

#ifdef __cplusplus
extern "C" {
#endif

void msp_file_close(msp_file_t *f);

msp_tdata_t *msp_gen_train(const msp_bed_t *bed, msp_file_t *fx, int32_t ext, double frac_pos);
msp_tdata_t *msp_tdata_read(const char *fn);
void msp_tdata_dump(FILE *fp, const msp_bed_t *bed, const msp_tdata_t *d);
void msp_tdata_destroy(msp_tdata_t *d);

// FASTX reader
msp_file_t *msp_fastx_open(const char *fn);
int32_t msp_fastx_read(msp_file_t *fp, msp_cstr_t *name, msp_cstr_t *seq);

// BED reader
msp_file_t *msp_bed_open(const char *fn);
msp_bed_t *msp_bed_read(const char *fn);
void msp_bed_sort(msp_bed_t *bed);
void msp_bed_idxctg(msp_bed_t *bed);
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
