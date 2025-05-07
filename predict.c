#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "kann.h"
#include "msppriv.h"
#include "ketopt.h"

void msp_predict1(msp_pdata_t *t, kann_t *ann, int64_t len, const uint8_t *seq, int32_t mb_sz, int32_t type)
{
	int32_t l, k, ext, alen, n_in, n_out, i_out;
	int64_t i, *mb2i;
	uint8_t x, *s;
	float *x1;

	for (i = 0, l = k = 0, x = 0; i < len; ++i) {
		int c = seq[i];
		uint64_t z = (uint64_t)-1;
		if (c < 4) {
			x = (x << 2 | c) & 0xf;
			if (type == 0) { // donor
				if (x == (2<<2|3)) z = (i-1)<<3 | 0<<2 | 0<<1; // GT
				else if (x == (0<<2|1)) z = (i+1)<<3 | 0<<2 | 1<<1; // AC
			} else { // acceptor
				if (x == (0<<2|2)) z = (i+1)<<3 | 1<<2 | 0<<1; // AG
				else if (x == (1<<2|3)) z = (i-1)<<3 | 1<<2 | 1<<1; // CT
			}
			if (z != (uint64_t)-1) {
				MSP_GROW(msp_pdata1_t, t->a, t->n, t->m);
				t->a[t->n++].x = z;
			}
		} else l = 0, x = 0;
	}

	i_out = kann_find(ann, KANN_F_OUT, 0);
	n_in = kann_dim_in(ann);
	n_out = kann_dim_out(ann);
	assert(n_in % 4 == 0 && n_out == 2);
	alen = n_in / 4;
	assert(alen % 2 == 0);
	ext = (alen - 2) / 2;

	s = MSP_CALLOC(uint8_t, alen);
	x1 = MSP_CALLOC(float, mb_sz * kann_dim_in(ann));
	mb2i = MSP_CALLOC(int64_t, mb_sz);
	for (i = 0; i < t->n; ++i)
		t->a[i].f = -1.0f, t->a[i].s = -1;
	for (i = 0; i < t->n;) {
		int64_t j;
		int32_t k, l;
		const float *y1;
		fprintf(stderr, "%ld\n", (long)i);
		for (j = i, k = 0; j < t->n && k < mb_sz; ++j) {
			int rc;
			rc = msp_get_seq_in_place(s, len, seq, t->a[j].x, ext);
			if (rc < 0) continue;
			msp_seq2vec(alen, s, &x1[k * n_in]);
			mb2i[k++] = j;
		}
		kann_feed_bind(ann, KANN_F_IN, 0, &x1);
		kann_eval_out(ann);
		y1 = ann->v[i_out]->x;
		for (l = 0; l < k; ++l)
			t->a[mb2i[l]].f = y1[l * n_out + 1];
		i = j;
	}
	free(mb2i); free(x1); free(s);
}

void msp_predict_print(kann_t *ann, msp_file_t *fx, int32_t mb_sz, int32_t type)
{
	int32_t len;
	const char *name, *seq;
	msp_pdata_t t = {0,0,0};
	kstring_t out = {0,0,0};
	while ((len = msp_fastx_read(fx, &name, &seq)) >= 0) {
		int64_t i;
		t.n = 0;
		msp_predict1(&t, ann, len, (uint8_t*)seq, mb_sz, type);
		for (i = 0; i < t.n; ++i) {
			out.l = 0;
			msp_sprintf_lite(&out, "%s\t%ld\t%c\t%c\t", name, (long)(t.a[i].x>>3), "+-"[t.a[i].x>>2&1], "DA"[t.a[i].x>>1&1]);
			fwrite(out.s, 1, out.l, stdout);
			printf("%.4f\n", t.a[i].f);
		}
	}
	free(out.s); free(t.a);
}
/*
void msp_eval(kann_t *ann, msp_file_t *fx, const msp_bed_t *bed, int32_t mb_size, int32_t type)
{
	int32_t len;
	const char *name, *seq;
	while ((len = msp_fastx_read(fx, &name, &seq)) >= 0) {
		int32_t cid, l, k;
		int64_t i;
		uint8_t x;
		msp128a_t td = {0,0,0}, ta = {0,0,0}, *pt;

		cid = msp_strmap_get(bed->h, name);
		if (cid < 0) continue;
		msp_gen_pos(&td, &ta, bed, cid, len, seq);
		pt = type == 0? &td : &ta;

		free(td.a); free(ta.a);
	}
}
*/
int main_predict(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int32_t c, n_thread = 1, mb_sz = 128, type = 0; // 0 for donor and 1 for acceptor
	msp_bed_t *bed = 0;
	msp_file_t *fx;
	kann_t *ann;
	char *fn_bed = 0;

	while ((c = ketopt(&o, argc, argv, 1, "adt:b:e:", 0)) >= 0) {
		if (c == 't') n_thread = atoi(o.arg);
		else if (c == 'b') mb_sz = atoi(o.arg);
		else if (c == 'a') type = 1;
		else if (c == 'd') type = 0;
		else if (c == 'e') fn_bed = o.arg;
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: minisplice predict [options] <in.kan> <in.fastx>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -a          acceptor model (donor model by default)\n");
		fprintf(stderr, "  -e FILE     annotated splice sites in BED12 []\n");
		fprintf(stderr, "  -t INT      number of threads [%d]\n", n_thread);
		fprintf(stderr, "  -b INT      minibatch size [%d]\n", mb_sz);
		return 1;
	}

	ann = kann_load(argv[o.ind]);
	if (n_thread > 1) kann_mt(ann, n_thread, mb_sz);
	fx = msp_fastx_open(argv[o.ind+1]);
	assert(ann && fx);
	if (fn_bed) {
		bed = msp_bed_read(fn_bed);
		msp_bed_idxctg(bed);
		assert(bed);
	} else {
		msp_predict_print(ann, fx, mb_sz, type);
	}

	if (bed) msp_bed_destroy(bed);
	msp_file_close(fx);
	kann_delete(ann);
	return 0;
}
