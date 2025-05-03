#include "msppriv.h"

#include "ksort.h"
#define msp196x_key(a) ((a).x)
KRADIX_SORT_INIT(msp196x, msp196_t, msp196x_key, 64)
#define msp196y_key(a) ((a).y)
KRADIX_SORT_INIT(msp196y, msp196_t, msp196y_key, 64)

void msp_bed_sort(msp_bed_t *bed)
{
	int64_t i, i0;
	msp196_t *a;
	msp_bed1_t **t;
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

msp_bedview_t *msp_bed_select_sp(const msp_bed_t *bed)
{
	msp_bedview_t *bv;
	int64_t i0, i, max_en = 0;
	if (bed->n == 0) return 0;
	bv = MSP_CALLOC(msp_bedview_t, 1);
	bv->bed = bed;
	for (i0 = 0, max_en = bed->a[0]->en, i = 1; i <= bed->n; ++i) {
		if (i == bed->n || bed->a[i0]->cid != bed->a[i]->cid || bed->a[i0]->st > max_en) {
			int64_t j, nf = 0, nr = 0;
			for (j = i0; j < i; ++j) {
				if (bed->a[j]->n_blk > 1) {
					if (bed->a[j]->strand < 0) ++nr;
					else if (bed->a[j]->strand > 0) ++nf;
				}
			}
			if (nf + nr > 0 && nf * nr == 0) {
				for (j = i0; j < i; ++j) {
					MSP_GROW(msp_bed1_p, bv->a, bv->n, bv->m);
					bv->a[bv->n++] = bed->a[j];
				}
			}
			max_en = max_en > bed->a[i]->en? max_en : bed->a[i]->en;
			i = i0;
		} else {
			max_en = max_en > bed->a[i]->en? max_en : bed->a[i]->en;
		}
	}
	return bv;
}
