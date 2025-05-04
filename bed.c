#include <stdio.h>
#include "msppriv.h"

#include "ksort.h"
#define key_x(a) ((a).x)
KRADIX_SORT_INIT(msp192x, msp192_t, key_x, 64)
#define key_y(a) ((a).y)
KRADIX_SORT_INIT(msp192y, msp192_t, key_y, 64)
#define key_self(a) (a)
KRADIX_SORT_INIT(msp64, uint64_t, key_self, 64)

void msp_bed_destroy(msp_bed_t *bed)
{
	int64_t i;
	for (i = 0; i < bed->n; ++i)
		free(bed->a[i]);
	free(bed->a);
	free(bed->c);
	msp_strmap_destroy(bed->h);
	free(bed);
}

int msp_bed_is_sorted(const msp_bed_t *bed)
{
	int64_t i;
	for (i = 1; i < bed->n; ++i) {
		const msp_bed1_t *p = bed->a[i-1], *q = bed->a[i];
		if (!(p->cid < q->cid || (p->cid == q->cid && p->st <= q->st)))
			break;
	}
	return i == bed->n? 1 : 0;
}

void msp_bed_sort(msp_bed_t *bed)
{
	int64_t i, i0;
	msp192_t *a;
	msp_bed1_t **t;
	if (msp_bed_is_sorted(bed)) return; // no need to sort
	a = MSP_MALLOC(msp192_t, bed->n);
	t = MSP_MALLOC(msp_bed1_t*, bed->n);
	for (i = 0; i < bed->n; ++i) {
		a[i].x = bed->a[i]->cid;
		a[i].y = bed->a[i]->st;
		a[i].z = i;
		t[i] = bed->a[i];
	}
	radix_sort_msp192x(a, a + bed->n);
	for (i0 = 0, i = 1; i <= bed->n; ++i) {
		if (i == bed->n || a[i0].x != a[i].x) {
			radix_sort_msp192y(&a[i0], &a[i]);
			i0 = i;
		}
	}
	for (i = 0; i < bed->n; ++i)
		bed->a[i] = t[a[i].z];
	free(t);
	free(a);
}

void msp_bed_idxctg(msp_bed_t *bed)
{
	int64_t i0, i;
	if (!msp_bed_is_sorted(bed)) msp_bed_sort(bed);
	bed->c = MSP_CALLOC(msp_bedctg_t, bed->h->n);
	for (i0 = 0, i = 1; i <= bed->n; ++i) {
		if (i == bed->n || bed->a[i0]->cid != bed->a[i]->cid) {
			bed->c[bed->a[i0]->cid].n = i - i0;
			bed->c[bed->a[i0]->cid].off = i0;
			i0 = i;
		}
	}
}

void msp_bed_format(kstring_t *out, const msp_bed1_t *b)
{
	int32_t i;
	out->l = 0;
	msp_sprintf_lite(out, "%s\t%ld\t%ld\t", b->ctg, b->st, b->en);
	if (b->name) msp_sprintf_lite(out, "%s", b->name);
	else msp_sprintf_lite(out, ".");
	if (b->score >= 0) msp_sprintf_lite(out, "\t%d", b->score);
	else msp_sprintf_lite(out, "\t.");
	if (b->st2 != b->st || b->en2 != b->en || b->n_blk > 0) {
		msp_sprintf_lite(out, "\t%ld\t%ld", b->st2, b->en2);
		if (b->n_blk > 0) {
			msp_sprintf_lite(out, "\t.\t%d\t", b->n_blk);
			for (i = 0; i < b->n_blk; ++i)
				msp_sprintf_lite(out, "%ld,", b->blk[i].en - b->blk[i].st);
			msp_sprintf_lite(out, "\t");
			for (i = 0; i < b->n_blk; ++i)
				msp_sprintf_lite(out, "%ld,", b->blk[i].st - b->st);
		}
	}
}

/*****************************
 * Generate negative regions *
 *****************************/

static void msp_reg_merge(const msp_bed_t *bed, const msp_bedctg_t *c, msp192a_t *t, int32_t strand)
{
	int64_t st, en, i;
	t->n = t->m = 0, t->a = 0;
	if (c->n == 0) return;
	for (i = c->off, st = en = -1; i < c->off + c->n; ++i) {
		if (bed->a[i]->n_blk < 2 || bed->a[i]->strand * strand <= 0) continue;
		if (bed->a[i]->st > en) {
			if (en > 0) {
				MSP_GROW(msp192_t, t->a, t->n, t->m);
				t->a[t->n].x = st, t->a[t->n].y = en, t->a[t->n++].z = strand > 0? 0 : 1;
			}
			st = bed->a[i]->st, en = bed->a[i]->en;
		} else en = en > bed->a[i]->en? en : bed->a[i]->en;
	}
	MSP_GROW(msp192_t, t->a, t->n, t->m);
	t->a[t->n].x = st, t->a[t->n].y = en, t->a[t->n++].z = strand > 0? 0 : 1;
}

static void msp_reg_sub(msp192a_t *t, const msp192a_t *a, const msp192a_t *b)
{
	int64_t i, j;
	if (a->n < 1) return;
	for (i = j = 0; i < a->n; ++i) {
		int64_t k, st1 = a->a[i].x, en1 = a->a[i].y, x;
		while (j < b->n && b->a[j].y <= st1) ++j; // skip b-regions w/o overlaps
		for (k = j, x = st1; k < b->n; ++k) { // adapted from "bedtk sub"
			int64_t st0 = b->a[k].x, en0 = b->a[k].y;
			if (st0 >= en1) break;
			if (st0 < st1) st0 = st1;
			if (en0 > en1) en0 = en1;
			if (st0 > x) {
				MSP_GROW(msp192_t, t->a, t->n, t->m);
				t->a[t->n].x = x, t->a[t->n].y = st0, t->a[t->n++].z = a->a[i].z;
			}
			x = en0;
		}
		if (x < en1) {
			MSP_GROW(msp192_t, t->a, t->n, t->m);
			t->a[t->n].x = x, t->a[t->n].y = en1, t->a[t->n++].z = a->a[i].z;
		}
	}
}

msp192_t *msp_bed_gen_negreg(const msp_bed_t *bed, int32_t cid, int64_t *nn)
{
	int64_t i;
	msp192a_t a[2], t;
	const msp_bedctg_t *c;
	if (cid >= bed->h->n) return 0;
	c = &bed->c[cid];
	msp_reg_merge(bed, c, &a[0], 1);
	msp_reg_merge(bed, c, &a[1], -1);
	t.n = t.m = 0, t.a = 0;
	msp_reg_sub(&t, &a[0], &a[1]);
	msp_reg_sub(&t, &a[1], &a[0]);
	radix_sort_msp192x(t.a, t.a + t.n);
	for (i = 0; i < t.n; ++i) t.a[i].z ^= 1; // flip the strand
	*nn = t.n;
	return t.a;
}

/*******************************
 * Generate training sequences *
 *******************************/

static void msp_uint64_dedup(msp64a_t *t)
{
	int64_t i0, i, k;
	radix_sort_msp64(t->a, t->a + t->n);
	for (i0 = 0, i = 1, k = 0; i <= t->n; ++i) {
		if (i == t->n || t->a[i0] != t->a[i]) {
			t->a[k++] = t->a[i0];
			i0 = i;
		}
	}
	t->n = k;
}

static void msp_gen_pos(msp64a_t *td, msp64a_t *ta, const msp_bed_t *bed, int32_t cid, int64_t len, const uint8_t *seq)
{
	const msp_bedctg_t *c = &bed->c[cid];
	int64_t i, j, n_noncan = 0;
	for (i = c->off; i < c->off + c->n; ++i) {
		const msp_bed1_t *b = bed->a[i];
		if (b->n_blk < 2 || b->strand == 0) continue;
		for (j = 1; j < b->n_blk; ++j) {
			int64_t p = b->blk[j-1].en, q = b->blk[j].st;
			assert(0 < p && p < q && q < len);
			//fprintf(stderr, "%c%c-%c%c\t%c\n", "ACGTN"[seq[p]], "ACGTN"[seq[p+1]], "ACGTN"[seq[q-2]], "ACGTN"[seq[q-1]], b->strand > 0? '+' : '-');
			if (b->strand > 0) {
				if (seq[p] == 2 && seq[p+1] == 3 && seq[q-2] == 0 && seq[q-1] == 2) { // GT-AG on the forward strand
					MSP_GROW(uint64_t, td->a, td->n, td->m);
					td->a[td->n++] = p<<3 | 0<<2 | 0<<1;
					MSP_GROW(uint64_t, ta->a, ta->n, ta->m);
					ta->a[ta->n++] = q<<3 | 0<<2 | 1<<1;
				} else ++n_noncan;
			} else if (b->strand < 0) {
				if (seq[p] == 1 && seq[p+1] == 3 && seq[q-2] == 0 && seq[q-1] == 1) { // CT-AC on the reverse strand
					MSP_GROW(uint64_t, td->a, td->n, td->m);
					td->a[td->n++] = q<<3 | 1<<2 | 0<<1;
					MSP_GROW(uint64_t, ta->a, ta->n, ta->m);
					ta->a[ta->n++] = p<<3 | 1<<2 | 1<<1;
				} else ++n_noncan;
			}
		}
	}
	msp_uint64_dedup(td);
	msp_uint64_dedup(ta);
	if (msp_verbose >= 3)
		fprintf(stderr, "[M::%s] collected %ld donors and %ld acceptors from \"%s\"; dropped %ld non-canonical introns\n",
			__func__, (long)td->n, (long)ta->n, bed->h->a[cid], (long)n_noncan);
}

static inline double msp_splitmix64(uint64_t *x)
{
	union { uint64_t i; double d; } u;
	uint64_t z = ((*x) += 0x9e3779b97f4a7c15ULL);
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
	z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
	z = z ^ (z >> 31);
	u.i = 0x3FFULL << 52 | z >> 12;
	return u.d - 1.0;
}

static void msp_neg_sample(msp64a_t *t, int64_t n, uint64_t *x) // reservior sampling
{
	int64_t i, k, y;
	if (t->n < n) return; // no need to downsample
	for (i = k = 0; i < t->n; ++i) {
		y = k++ < n? k - 1 : (int64_t)(msp_splitmix64(x) * k);
		if (y < n) t->a[y] = t->a[i];
	}
	t->n = n;
}

static void msp_gen_neg(msp64a_t *td, msp64a_t *ta, int64_t len, const uint8_t *seq, int64_t n_negreg, const msp192_t *negreg)
{
	int64_t j;
	for (j = 0; j < n_negreg; ++j) {
		int64_t i, st = negreg[j].x, en = negreg[j].y, rev = negreg[j].z;
		int32_t l;
		uint8_t x = 0;
		for (i = st, l = 0; i < en; ++i) { // this is similar to the k-mer counting loop
			int c = seq[i];
			if (c < 4) {
				x = (x<<2 | c) & 0xf;
				if (++l >= 2) {
					if (!rev) {
						if (x == (2<<2|3)) { // GT on forward
							MSP_GROW(uint64_t, td->a, td->n, td->m);
							td->a[td->n++] = (i-1)<<3 | 0<<2 | 0<<1 | 1;
						} else if (x == (0<<2|2)) { // AG on forward
							MSP_GROW(uint64_t, ta->a, ta->n, ta->m);
							ta->a[ta->n++] = (i+1)<<3 | 0<<2 | 1<<1 | 1;
						}
					} else {
						if (x == (0<<2|1)) { // AC on reverse
							MSP_GROW(uint64_t, td->a, td->n, td->m);
							td->a[td->n++] = (i+1)<<3 | 1<<2 | 0<<1 | 1;
						} else if (x == (1<<2|3)) { // CT on reverse
							MSP_GROW(uint64_t, ta->a, ta->n, ta->m);
							ta->a[ta->n++] = (i-1)<<3 | 1<<2 | 1<<1 | 1;
						}
					}
				}
			} else l = 0, x = 0;
		}
	}
}

msp_tdata_t *msp_gen_train_seq(const msp_bed_t *bed, int32_t cid, int64_t len, const uint8_t *seq, int32_t ext, double frac_pos, int64_t *nn)
{
	uint64_t x = 11;
	int64_t n_negreg;
	msp192_t *negreg;
	msp64a_t pd = {0,0,0}, pa = {0,0,0};
	msp64a_t nd = {0,0,0}, na = {0,0,0};

	if (cid >= bed->h->n || frac_pos <= 0.0 || frac_pos >= 1.0) return 0;
	negreg = msp_bed_gen_negreg(bed, cid, &n_negreg);
	msp_gen_pos(&pd, &pa, bed, cid, len, seq);
	msp_gen_neg(&nd, &na, len, seq, n_negreg, negreg);
	free(negreg);
	msp_neg_sample(&nd, (int64_t)(pd.n / frac_pos * (1.0 - frac_pos) + .499), &x);
	msp_neg_sample(&na, (int64_t)(pa.n / frac_pos * (1.0 - frac_pos) + .499), &x);
	return 0;
}

msp_tdata_t *msp_gen_train(const msp_bed_t *bed, msp_file_t *fx, int32_t ext, double frac_pos)
{
	int32_t len;
	const char *name, *seq;
	while ((len = msp_fastx_read(fx, &name, &seq)) >= 0) {
		int32_t cid;
		int64_t n;
		cid = msp_strmap_get(bed->h, name);
		if (cid < 0) continue; // skip if not found in the BED file
		msp_gen_train_seq(bed, cid, len, (uint8_t*)seq, ext, frac_pos, &n);
	}
	return 0;
}
