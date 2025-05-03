#include <stdio.h>
#include "msppriv.h"

#include "ksort.h"
#define msp196x_key(a) ((a).x)
KRADIX_SORT_INIT(msp196x, msp196_t, msp196x_key, 64)
#define msp196y_key(a) ((a).y)
KRADIX_SORT_INIT(msp196y, msp196_t, msp196y_key, 64)

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
	msp196_t *a;
	msp_bed1_t **t;
	if (msp_bed_is_sorted(bed)) return; // no need to sort
	a = MSP_MALLOC(msp196_t, bed->n);
	t = MSP_MALLOC(msp_bed1_t*, bed->n);
	for (i = 0; i < bed->n; ++i) {
		a[i].x = bed->a[i]->cid;
		a[i].y = bed->a[i]->st;
		a[i].z = i;
		t[i] = bed->a[i];
	}
	radix_sort_msp196x(a, a + bed->n);
	for (i0 = 0, i = 1; i <= bed->n; ++i) {
		if (i == bed->n || a[i0].x != a[i].x) {
			radix_sort_msp196y(&a[i0], &a[i]);
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
	for (i0 = 0, i = 1; i < bed->n; ++i) {
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
