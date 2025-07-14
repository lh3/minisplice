#include <zlib.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include "msppriv.h"

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

/****************
 * General file *
 ****************/

msp_file_t *msp_file_open_by_line(const char *fn)
{
	gzFile fp;
	msp_file_t *f;
	kstream_t *ks;
	fp = fn == 0 || strcmp(fn, "-") == 0? gzdopen(0, "rb") : gzopen(fn, "rb");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	f = MSP_CALLOC(msp_file_t, 1);
	f->fp = ks;
	f->type = MSP_FT_LINE;
	return f;
}

static void msp_file_init_buf(msp_file_t *fp)
{
	kstring_t *buf;
	buf = MSP_CALLOC(kstring_t, 1);
	fp->buf = buf;
}

void msp_file_close(msp_file_t *f)
{
	kstring_t *s;
	if (f == 0) return;
	if (f->type == MSP_FT_LINE) {
		kstream_t *ks = (kstream_t*)f->fp;
		gzFile fp = ks->f;
		ks_destroy(ks);
		gzclose(fp);
	} else if (f->type == MSP_FT_FASTX) {
		kseq_t *ks = (kseq_t*)f->fp;
		gzFile fp = ks->f->f;
		kseq_destroy(ks);
		gzclose(fp);
	}
	s = (kstring_t*)f->buf;
	if (s) free(s->s);
	free(f->buf);
	free(f);
}

/****************
 * FASTX reader *
 ****************/

const uint8_t msp_nt6_table[256] = {
    5, 1, 2, 3,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

msp_file_t *msp_fastx_open(const char *fn)
{
	gzFile fp;
	msp_file_t *f;
	kseq_t *ks;
	fp = fn == 0 || strcmp(fn, "-") == 0? gzdopen(0, "rb") : gzopen(fn, "rb");
	if (fp == 0) return 0;
	ks = kseq_init(fp);
	f = MSP_CALLOC(msp_file_t, 1);
	f->fp = ks;
	f->type = MSP_FT_FASTX;
	return f;
}

int64_t msp_fastx_read(msp_file_t *fp, msp_cstr_t *name, msp_cstr_t *seq)
{
	kseq_t *ks = (kseq_t*)fp->fp;
	long ret;
	*name = 0, *seq = 0;
	ret = kseq_read(ks);
	if (ret > 0) {
		size_t i;
		for (i = 0; i < ks->seq.l; ++i)
			ks->seq.s[i] = msp_nt6_table[(uint8_t)ks->seq.s[i]] - 1;
		*seq = ks->seq.s, *name = ks->name.s;
	}
	return ret;
}

/**************
 * BED reader *
 **************/

int msp_bed_read1(msp_file_t *fp, msp_bed1_t **b_)
{
	int32_t i, ret, ctg_len = 0, name_len = 0, tot_len, tot_cnt;
	kstream_t *ks;
	kstring_t *buf;
	msp_bed1_t t, *b = 0;
	char *p, *q, *bl = 0, *bs = 0; 

	assert(fp != 0 && fp->fp != 0 && fp->type == MSP_FT_LINE);

	*b_ = 0;
	if (fp->buf == 0) msp_file_init_buf(fp);
	buf = (kstring_t*)fp->buf;
	ks = (kstream_t*)fp->fp;
	memset(&t, 0, sizeof(t));

	ret = ks_getuntil(ks, KS_SEP_LINE, buf, 0);
	if (ret < 0) return ret;
	for (p = q = buf->s, i = 0;; ++p) {
		if (*p == 0 || *p == '\t') {
			int32_t c = *p; // the original character
			*p = 0;
			if (i == 0) { // ctg
				t.ctg = q, ctg_len = p - q;
			} else if (i == 1) { // start
				t.st = atol(q); // TODO: watch out integer overflow!
				if (t.st < 0) break;
			} else if (i == 2) { // end
				t.en = atol(q);
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
				t.st2 = p - q == 0 || *q == '.'? -1 : atol(q);
			} else if (i == 7) {
				t.en2 = p - q == 0 || *q == '.'? -1 : atol(q);
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

	if (i < 3) return -2; // need at least three fields
	if (t.st < 0 || t.en < 0 || t.st > t.en || ctg_len == 0) return -3;
	if (t.st2 >= 0 && t.en2 >= 0 && (t.st2 < t.st || t.en2 > t.en || t.st2 > t.en2)) return -3;

	tot_len = sizeof(msp_blk1_t) * t.n_blk + (ctg_len + 1) + (name_len + 1);
	tot_cnt = (tot_len + sizeof(msp_blk1_t) - 1) / sizeof(msp_blk1_t);
	b = (msp_bed1_t*)calloc(sizeof(msp_bed1_t) + sizeof(msp_blk1_t) * tot_cnt, 1);
	b->cid = -1;
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
	*b_ = b;
	return 0;
}

msp_bed_t *msp_bed_read(const char *fn)
{
	msp_bed_t *bed;
	msp_bed1_t *b;
	msp_file_t *fp;
	int64_t lineno = 0;
	int rc;

	fp = msp_file_open_by_line(fn);
	if (fp == 0) return 0;
	bed = MSP_CALLOC(msp_bed_t, 1);
	bed->h = msp_strmap_init();
	while ((rc = msp_bed_read1(fp, &b)) != -1) {
		++lineno;
		if (b == 0) {
			if (msp_verbose >= 1)
				fprintf(stderr, "[E::%s] BED parsing error on line %ld\n", __func__, (long)lineno);
			continue;
		}
		b->cid = msp_strmap_add(bed->h, b->ctg);
		MSP_GROW(msp_bed1_t*, bed->a, bed->n, bed->m);
		bed->a[bed->n++] = b;
	}
	msp_file_close(fp);
	if (msp_verbose >= 3)
		fprintf(stderr, "[M::%s] read %ld BED records\n", __func__, (long)bed->n);
	return bed;
}

/****************
 * sdata reader *
 ****************/

msp_sdata_t *msp_sdata_read(const char *fn)
{
	msp_file_t *fp;
	msp_sdata_t *d;
	kstring_t str = {0,0,0}, out = {0,0,0};
	kstream_t *ks;
	int32_t ret;

	fp = msp_file_open_by_line(fn);
	if (fp == 0) return 0;
	ks = (kstream_t*)fp->fp;
	d = MSP_CALLOC(msp_sdata_t, 1);
	while ((ret = ks_getuntil(ks, KS_SEP_LINE, &str, 0)) >= 0) {
		char *p, *q, *seq = 0;
		int32_t i, j, len = -1, type = -1;
		int32_t label = -1;
		out.l = 0;
		for (p = q = str.s, i = 0;; ++p) {
			if (*p == 0 || *p == '\t' || *p == ' ') {
				int c = *p;
				*p = 0;
				if (i == 0) {
					if (*q >= '0' && *q <= '9')
						type = 0, label = atol(q); // 2-column
					else type = 1; // 6-column
				} else if (type == 1 && i == 1) {
					msp_sprintf_lite(&out, "%s", q);
				} else if (type == 1 && i >= 2 && i <= 3) {
					msp_sprintf_lite(&out, ":%s", q);
				} else if (type == 1 && i == 4) {
					msp_sprintf_lite(&out, ":%s:%s", str.s, q);
					label = atol(q);
				} else if ((type == 0 && i == 1) || (type == 1 && i == 5)) {
					seq = q, len = p - q;
				}
				++i, q = p + 1;
				if (c == 0) break;
			}
		}
		if (seq && label >= 0) {
			msp_sdata1_t *r;
			d->n_label = d->n_label > label + 1? d->n_label : label + 1;
			if (d->len <= 0) d->len = len;
			else if (d->len != len) {
				if (msp_verbose >= 2)
					fprintf(stderr, "[W::%s] training data contain sequences of different lengths: %d != %d\n", __func__, d->len, len);
				continue;
			}
			MSP_GROW(msp_sdata1_t, d->a, d->n, d->m);
			r = &d->a[d->n++];
			r->label = label;
			r->seq = MSP_CALLOC(uint8_t, d->len + out.l + 1);
			r->ctg = (char*)r->seq + d->len;
			for (j = 0; j < d->len; ++j)
				r->seq[j] = msp_nt6_table[(uint8_t)seq[j]] - 1;
			if (out.l) memcpy(r->ctg, out.s, out.l);
		}
	}
	free(out.s);
	free(str.s);
	return d;
}

void msp_sdata_destroy(msp_sdata_t *d)
{
	int32_t i;
	for (i = 0; i < d->n; ++i)
		free(d->a[i].seq);
	free(d->a); free(d);
}

/***************
 * eval reader *
 ***************/

msp_eval_t *msp_eval_read(const char *fn)
{
	msp_file_t *fp;
	msp_eval_t *e = 0;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	int32_t ret;

	fp = msp_file_open_by_line(fn);
	if (fp == 0) return 0;
	ks = (kstream_t*)fp->fp;
	while ((ret = ks_getuntil(ks, KS_SEP_LINE, &str, 0)) >= 0) {
		char *p, *q;
		int32_t i, type = 0, b = -1;
		for (p = q = str.s, i = 0;; ++p) {
			if (*p == 0 || *p == '\t' || *p == ' ') {
				int c = *p;
				*p = 0;
				if (i == 0) {
					if (strcmp(q, "ST") == 0) type = 1;
					else if (strcmp(q, "BN") == 0) type = 2;
				} else if (type == 1 && i == 1) {
					double step;
					step = atof(q);
					e = msp_eval_init(step);
				} else if (type == 2 && i == 1) {
					assert(e);
					b = atoi(q);
				} else if (type == 2 && i == 2) {
					assert(b >= 0 && b < e->n_bin);
					e->bin[b].mt = atol(q);
				} else if (type == 2 && i == 3) {
					e->bin[b].mp = atol(q);
				}
				++i, q = p + 1;
				if (c == 0) break;
			}
		}
	}
	free(str.s);
	msp_file_close(fp);
	msp_eval_update(e);
	return e;
}
