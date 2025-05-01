#include <zlib.h>
#include <assert.h>
#include <stdlib.h>
#include "msppriv.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

static void msp_file_init_buf(msp_file_t *fp)
{
	kstring_t *buf;
	buf = MSP_CALLOC(kstring_t, 1);
	fp->buf = buf;
}

msp_bed1_t *msp_bed_read1(msp_file_t *fp, uint32_t *err)
{
	int32_t i, ret, ctg_len = 0, name_len = 0, tot_len, tot_cnt;
	kstream_t *ks;
	kstring_t *buf;
	msp_bed1_t t, *b = 0;
	char *p, *q, *bl = 0, *bs = 0; 

	assert(fp != 0 && fp->fp != 0 && fp->type == MSP_FT_BED);

	if (fp->buf == 0) msp_file_init_buf(fp);
	*err = 0;
	buf = (kstring_t*)fp->buf;
	ks = (kstream_t*)fp->fp;
	memset(&tmp, 0, sizeof(tmp));

	ret = ks_getuntil(ks, KS_SEP_LINE, buf, 0);
	if (ret < 0) return 0; // TODO: detect error code
	for (p = q = buf->a, i = 0;; ++p) {
		if (*p == 0 || *p == '\t') {
			int32_t c = *p; // the original character
			*p = 0;
			if (i == 0) { // ctg
				t.ctg = q, ctg_len = p - q;
			} else if (i == 1) { // start
				t.st = t.st2 = atol(q); // TODO: watch out integer overflow!
				if (t.st < 0) break;
			} else if (i == 2) { // end
				t.en = t.en2 = atol(q);
				if (t.en < 0) break;
			} else if (i == 3) {
				if (*q == '.' && p - q == 1)
					t.name = 0, name_len = 0;
				else
					t.name = q, name_len = p - q;
			} else if (i == 4) { // BED score
				t.score = *q >= '0' && *q <= '9'? atof(q) : -1.0;
			} else if (i == 5) { // strand
				t.strand = *q == '+'? 1 : *q == '-'? -1 : 0;
			} else if (i == 6) {
				t.st2 = atol(q);
			} else if (i == 7) {
				t.en2 = atol(q);
			} else if (i == 9) {
				if (*q < '0' || *q > '9') break;
				t.n_blk = atol(q);
			} else if (i == 10) {
				bl = q;
			} else if (i == 11) {
				bs = q;
			}
			++i, q = p + 1;
			if (c == 0 || i >= 12) break;
		}
	}

	if (i < 3) return 0; // need at least three fields
	if (t.st < 0 || t.en < 0 || t.st > t.en || ctg_len == 0) return 0;
	if (t.st2 < t.st || t.en2 > t.en || t.st2 > t.en2) return 0;

	tot_len = sizeof(msp_blk1_t) * t.n_blk + (ctg_len + 1) + (name_len + 1);
	tot_cnt = (tot_len + sizeof(msp_blk1_t) - 1) / sizeof(msp_blk1_t);
	b = (msp_bed1_t*)calloc(sizeof(msp_bed1_t) + sizeof(msp_blk1_t) * tot_cnt, 1);
	b->ctg = (char*)&b->blk[t.n_blk];
	b->name = t.name? b->ctg + ctg_len + 1 : 0;
	memcpy(b->ctg, t.ctg, ctg_len + 1);
	if (t.name) memcpy(b->name, t.name, name_len + 1);
	b->st = t.st, b->en = t.en, b->st2 = t.st2, b->en2 = t.en2, b->score = t.score;
	b->n_blk = t.n_blk, b->strand = t.strand;

	if (i >= 12 && bl && bs) {
		for (i = 0; i < t.n_blk; ++i) {
			msp_blk1_t *p = &b->blk[i];
			p->st = t.st  + strtol(bs, &bs, 10); ++bs;
			p->en = p->st + strtol(bl, &bl, 10); ++bl;
		}
	}
	return b;
}
